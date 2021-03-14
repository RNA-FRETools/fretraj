#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import argparse
import json
import mdtraj as md
import numba as nb
import copy
import pickle

try:
    import LabelLib as ll
except ModuleNotFoundError:
    print('\nNote: the LabelLib module is not installed. \nACV calculations will use a Python-only algorithm\n')
    _LabelLib_found = False
else:
    _LabelLib_found = True

from fretraj import export
from fretraj import fret
from fretraj import grid
from fretraj import metadata

DISTANCE_SAMPLES = 100000

package_directory = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(package_directory, 'periodic_table.json')) as f:
    _periodic_table = json.load(f)

VDW_RADIUS = dict((key, _periodic_table[key]["vdw_radius"]) for key in _periodic_table.keys())

# set None if the parameter is required
_label_dict = {'Position': {'pd_key': {'attach_id': ((int, float), None),
                                       'linker_length': ((int, float), None),
                                       'linker_width': ((int, float), None),
                                       'dye_radius1': ((int, float), None),
                                       'dye_radius2': ((int, float), None),
                                       'dye_radius3': ((int, float), None),
                                       'cv_thickness': ((int, float), 0),
                                       'grid_spacing': ((int, float), 1.0),
                                       'mol_selection': (str, 'all'),
                                       'simulation_type': (str, 'AV3'),
                                       'cv_fraction': ((int, float), 0),
                                       'state': (int, 1),
                                       'frame_mdtraj': (int, 0),
                                       'use_LabelLib': (bool, False),
                                       'contour_level_AV': ((int, float), 0),
                                       'contour_level_CV': ((int, float), 0.7),
                                       'b_factor': (int, 100),
                                       'gaussian_resolution': (int, 2),
                                       'grid_buffer': ((int, float), 2.0),
                                       'transparent_AV': (bool, True)}},
               'Distance': {'pd_key': {'R0': ((int, float), None),
                                       'n_dist': (int, 10**6)}}}

_label_dict_vals = {field: {'pd_key': {key: val[1] for key, val in _label_dict[field]['pd_key'].items()}}
                    for field in _label_dict.keys()}
_default_params = {field: {key: val[1] for key, val in _label_dict[field]['pd_key'].items() if val[1] is not None}
                   for field in _label_dict.keys()}


def parseCmd():
    """Parse the command line to get the input PDB file and the dye parameters

    Returns
    -------
    tuple of str
        tuple containing the filenames of the input PDB, the parameter file
        and the output ACV

    """
    parser = argparse.ArgumentParser(
        description='compute accessible-contact clouds for \
                     an MD trajectory or a given PDB structure')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + str(metadata['Version']))
    parser.add_argument('--path', action='version', version=f'package directory: {package_directory}',
                        help='Show package directory')
    parser.add_argument('-i', '--input', help='Input PDB structure (.pdb)', required=True)
    parser.add_argument('-p', '--parameters',
                        help='Parameter file (.json)', required=True)
    parser.add_argument('-o', '--output',
                        help='Output file of accessible contact volume (.pdb, .xyz, .dx)',
                        required=False)
    args = parser.parse_args()
    in_filePDB = args.input
    param_fileJSON = args.parameters
    out_fileACV = args.output
    return (in_filePDB, param_fileJSON, out_fileACV)


def labeling_params(param_file, verbose=True):
    """Parse the parameter file to get the configuration settings

    Parameters
    ----------
    param_file : str
        json-formatted parameter filename

    Returns
    -------
    dict
        position of labels, dye and grid parameters
    """
    with open(param_file) as f:
        labels_json = json.load(f)

    try:
        check_labels(labels_json, verbose)
    except KeyError as e:
        error_type = e.args[0]
        key = e.args[1]
        pos = e.args[2]
        field = e.args[3]
        print('{}: \'{}\' in {} {}. Exiting...'.format(error_type, key, field, pos))
    except TypeError as e:
        print('Wrong data type: {}'.format(e))
    else:
        return labels_json


def check_labels(labels, verbose=True):
    """Check integrity of parameter dictionary

    Parameters
    ----------
    labels : dict
        position of labels, dye and grid parameters
    verbose : bool, optional=True
        be verbose about missing parameter and their fall-back values
    """
    for field in _label_dict.keys():
        if field in labels.keys():
            for pos in labels[field].keys():
                if field == 'Position':
                    if 'simulation_type' not in labels[field][pos]:
                        labels[field][pos]['simulation_type'] = copy.copy(_default_params[field]['simulation_type'])
                    else:
                        if labels[field][pos]['simulation_type'] == 'AV1':
                            labels[field][pos]['dye_radius2'] = 0
                            labels[field][pos]['dye_radius3'] = 0
                    if 'state' in labels[field][pos] and 'frame_mdtraj' not in labels[field][pos]:
                        labels[field][pos]['frame_mdtraj'] = labels[field][pos]['state'] - 1
                    if 'frame_mdtraj' in labels[field][pos] and 'state' not in labels[field][pos]:
                        labels[field][pos]['state'] = labels[field][pos]['frame_mdtraj'] + 1
                elif field == 'Distance':
                    pass

                # check if all keys that are needed are defined
                for key, (t, _) in _label_dict[field]['pd_key'].items():
                    if key not in labels[field][pos]:
                        if key in _default_params[field].keys():
                            labels[field][pos][key] = copy.copy(_default_params[field][key])
                            if verbose:
                                print('Missing Key: \'{}\' in {} {}. Falling back to \"{}\"'.format(key, field, pos,
                                      _default_params[field][key]))
                        else:
                            raise KeyError('Missing Key', key, pos, field)
                    else:
                        if not isinstance(labels[field][pos][key], t):
                            raise TypeError('\'{}\' in {} {} must be of one of the following types: {}'.format(key,
                                            field, pos, t))

                # check if there are any unrecognized keys
                for key in labels[field][pos].keys():
                    if key not in _label_dict_vals[field]['pd_key'].keys():
                        raise KeyError('Unrecognized key', key, pos, field)
        else:
            labels[field] = None
            raise ValueError('Cannot read {} parameters from file: Missing field \'{}\'.'.format(field, field))


def save_obj(filename, obj):
    """
    Save a serialized object to a binary file

    Parameters
    ----------
    filename : str
        filename for pickle object
    obj : serializable object
    """
    with open(filename, 'wb') as f:
        try:
            pickle.dump(obj, f)
        except TypeError:
            try:
                obj_tmp = copy.copy(obj)
                obj_tmp.acv = copy.copy(obj.acv)
                del obj_tmp.av
                del obj_tmp.acv.ll_Grid3D
                print('Note: the LabelLib.Grid objects have been removed as they cannot be pickled.')
            except AttributeError:
                print('Error: Cannot pickle the passed object')
            else:
                pickle.dump(obj_tmp, f)


def load_obj(filename):
    """
    Load a serialized object from a binary file

    Parameters
    ----------
    filename : str
        filename of pickle obj
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)


def save_labels(filename, labels):
    """Write the ACV parameters to a .json file

    Parameters
    ----------
    filename : str
        filename for label parameters
    labels : dict
        position of labels, dye and grid parameters

    Examples
    --------

    >>> obj.save_labels('parameters.json')

    """
    with open(filename, 'w') as f:
        json.dump(labels, f, indent=2)


def printProgressBar(iteration, total, prefix='Progress:', suffix='complete', length=20, fill='█'):
    """Command line progress bar

    Parameters
    ----------
    iteration : int
        current iteration
    total : int
        total number of iterations
    prefix : str, optional='Progress:'
            string before progress bar
    suffix : str, optional='complete'
        string after progress bar
    length : int, optional=20
        length of progress bar
    fill : str, optional='█'
        fill character of progress bar

    Examples
    --------

    >>> printProgressBar(0, n)
    >>> for i in range(n):
            doSomething
            printProgressBar(i + 1, n)

    """
    percent = '{:0.0f}'.format(100 * iteration / total)
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r{} |{}| {}% {}'.format(prefix, bar, percent, suffix), end='\r')
    if iteration == total:
        print()


def pipeline_frames(structure, donor_site, acceptor_site, labels, frames, fret_pair):
    """Create a pipeline to compute multi-frame donor and acceptor ACVs and calculate a FRET trajectory

    Parameters
    ----------
    structure : mdtraj.Trajectory
                trajectory of atom coordinates loaded from a pdb, xtc or other file
    donor_site : str
                 reference identifier for the donor labeling position
    acceptor_site : str
                    reference identifier for the acceptor labeling position
    labels : dict
             dye, linker and setup parameters for the accessible volume calculation
    frames_mdtraj : list
                    list of frames on the trajectory to be used for the ACV calculation
    fret_pair : str
                Distance key specifying the donor acceptor pair

    Note
    ----
    Running a pipeline more memory-efficient than running Volume.from_frames()
    because the ACVs are stored only until the FRET efficiency is calculated.
    On the downside, no ACV trajectories(.xtc / .xyz) can be saved.
    """
    _labels = copy.copy(labels)
    acv = {}
    fret = []
    n_frames = len(frames)
    printProgressBar(0, n_frames)
    for i, frame in enumerate(frames):
        for dye, site in zip(['D', 'A'], [donor_site, acceptor_site]):
            _labels['Position'][site]['frame_mdtraj'] = frame
            _labels['Position'][site]['state'] = frame + 1
            acv[dye] = Volume(structure, site, _labels)
        fret.append(FRET(acv['D'], acv['A'], fret_pair, labels))
        printProgressBar(i + 1, n_frames)
    return fret


def save_mp_traj(filename, volume_list, units='A'):
    """Save a trajectory of dye mean positions as an xyz file

    Parameters
    ----------
    filename : str
        filename for mean position trajectory
    volume_list : array_like
        list of Volume instances
    units : {'A', 'nm'}
        distance units (Angstrom or nanometer)
    """
    mps = np.vstack([volume_list[i].acv.mp for i in range(len(volume_list)) if hasattr(volume_list[i].acv, 'mp')])
    mps = np.hstack((mps, np.ones((mps.shape[0], 1))))
    mean_mp = np.mean(mps, 0)

    xyz_str = export.xyz(mps, mean_mp, None)
    with open(filename, 'w') as fname:
        fname.write(xyz_str)


def acv_subsampling(volume_list, verbose=False):
    """Subsample the ACV to obtain identical number of points over all volumes.
    This allows to create mdtraj.Trajectory objects of the ACV which can be also save as .xtc or .xyz files

    Parameters
    ----------
    volume_list : array_like
        list of Volume instances
    verbose: bool
        be verbose about the number of subsamples
    """
    rng = np.random.default_rng()
    vols = ['AV', 'FV', 'CV']
    top = {}
    xyz_sub = {}
    for v in vols:
        if v == 'CV':
            xyz_list = list(map(lambda volume: volume.acv.cloud_xyzqt[volume.acv.cloud_xyzqt[:, 4] == 2, 0:3],
                                volume_list))
        elif v == 'FV':
            xyz_list = list(map(lambda volume: volume.acv.cloud_xyzqt[volume.acv.cloud_xyzqt[:, 4] < 2, 0:3],
                                volume_list))
        else:
            xyz_list = list(map(lambda volume: volume.acv.cloud_xyzqt[:, 0:3], volume_list))
        n_pts = min([xyz.shape[0] for xyz in xyz_list])
        top[v] = _create_volume_topology(n_pts)
        xyz_sub[v] = np.stack([rng.choice(xyz, n_pts, replace=False) for xyz in xyz_list])
        if verbose:
            print(f'{n_pts} points sampled in {v}')
    return xyz_sub, top


def create_acv_traj(volume_list, verbose=False):
    """Create a mdtraj.Trajectory object from a volume list

    Parameters
    ----------
    volume_list : array_like
        list of Volume instances
    verbose: bool
        be verbose about the number of subsamples
    """
    acv_traj = {}
    xyz_sub, top = acv_subsampling(volume_list, verbose=False)
    acv_traj['AV'] = md.Trajectory(xyz_sub['AV']/10, top['AV'])
    if any(volume_list[0].acv.tag_1d > 1):
        acv_traj['FV'] = md.Trajectory(xyz_sub['FV']/10, top['FV'])
        acv_traj['CV'] = md.Trajectory(xyz_sub['CV']/10, top['CV'])
    return acv_traj


def save_acv_traj(filename, volume_list, format='pdb', verbose=False):
    """Save a trajectory of ACVs as a full multi model PDB or a subsampled .xyz or .xtc
    (with identical number of points over all volumes)

    Parameters
    ----------
    filename : str
        filename for ACV trajectory
    volume_list : array_like
        list of Volume instances
    format : str
        file format of the ACV (for 'pdb' the full cloud is saved without subsampling;
        for 'xyz' k points in the volumes are randomly picked, where k is the number
        of points in the smallest volume)
    separate_CV : bool
        Save contact and free volume in separate files
        (choose this if you want to retain the CV information)
    """
    if format == 'pdb':
        if verbose:
            print('With format \"PDB\" no subsampling is performed.')
        file_str = ''
        for i, volume in enumerate(volume_list):
            file_str += f'MODEL {i+1}\n'
            file_str += export.pdb(volume.acv.cloud_xyzqt, volume.acv.mp)
            file_str += 'ENDMDL\n\n'
        with open(filename, 'w') as fname:
            fname.write(file_str)

    elif (format == 'xyz') or (format == 'xtc'):
        vols = ['AV']
        if any(volume_list[0].acv.tag_1d > 1):
            vols += ['FV', 'CV']
        base, suffix = os.path.splitext(filename)
        xyz_sub, top = acv_subsampling(volume_list, verbose=False)
        for v in vols:
            filename = f'{base}_{v}{suffix}'
            if format == 'xtc':
                with md.formats.XTCTrajectoryFile(filename, 'w') as f:
                    f.write(xyz_sub[v]/10)
            else:
                with md.formats.XYZTrajectoryFile(filename, 'w') as f:
                    f.write(xyz_sub[v])

            first_frame = md.Trajectory(xyz_sub[v][0]/10, top[v])
            first_frame.save_pdb(f'{base}_{v}.pdb')


def load_acv_traj(filename):
    """Load an AV, FV and CV from .xyz or .xtc file

    Parameters
    ----------
    filename : str
        filename of ACV trajectory

    Note
    ----
    Multi model PDB files cannot be loaded as a mdtraj.Trajectory because they are on purpose not subsampled
    and thus contain different number of points per volume
    """
    base, suffix = os.path.splitext(filename)
    acv_traj = {}
    for v in ['FV', 'FV', 'CV']:
        acv_traj[v] = md.load(f'{base}_{v}{suffix}', top=f'{base}_{v}.pdb')
    return acv_traj


def _create_volume_topology(n_atoms):
    """Create a topology for an accessible-contact volume

    Parameters
    ----------
    n_atoms : int
        number of points in the volume
    """
    atoms = pd.concat((pd.Series(range(n_atoms), name='serial'),
                       pd.Series(['D']*n_atoms, name='name'),
                       pd.Series(['D']*n_atoms, name='element'),
                       pd.Series([0]*n_atoms, name='resSeq'),
                       pd.Series(['X']*n_atoms, name='resName'),
                       pd.Series([0]*n_atoms, name='chainID')), axis=1)
    top = md.Topology.from_dataframe(atoms)
    return top


def save_structure_traj(filename, structure, frames, format='pdb'):
    """Save selected frames of a trajectory as multi-model PDB or XTC

    Parameters
    ----------
    filename : str
        filename for structure trajectory
    structure : mdtraj.Trajectory
        trajectory of atom coordinates loaded from a pdb, xtc or other file
    format : {'pdb', 'xtc'}
        trajectory file format
    """
    sliced_structure = structure.slice(frames)
    if format == 'xtc':
        sliced_structure.save_xtc(filename)
    else:
        sliced_structure.save_pdb(filename)


class ACV:
    """Reweighted accessible contact volume of a covalently linked fluorophore

    Parameters
    ----------
    grid_acv : LabelLib.Grid3D, optional=None
        accessible volume with added weight labels for free and contact volume
        (density label FV: 1.0, density label CV: 2.0); the two volumes are subsequently reweighted
    cv_thickness : float, optional=0
        width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
        **Note:** good first approx.: 2*min(dye_radii)
    cv_fraction : float, optional=0
        fraction of dyes that are within the contact volume
        (e.g. as determined by fluorescence anisotropy)
    cloud_xyzq : ndarray, optional=None
        array of x-,y-,z-coordinates and corresponding weights
        with a shape [n_gridpts(+), 4]
    use_LabelLib : bool, optional=True
        make use of LabelLib library to compute FRET values and distances
    """

    def __init__(self, grid_acv=None, cv_thickness=0, cv_fraction=0, cloud_xyzqt=None, use_LabelLib=True):
        if grid_acv is not None:
            self.n_gridpts = np.prod(grid_acv.shape)
            self.grid_1d, self.tag_1d = self._reweight_cv(grid_acv.grid, cv_thickness, cv_fraction)
            self.grid_3d = Volume.reshape_grid(self.grid_1d, grid_acv.shape)
            self.tag_3d = Volume.reshape_grid(self.tag_1d, grid_acv.shape)
            self.cloud_xyzqt = Volume.grid2pts(self.grid_3d, grid_acv.originXYZ, [grid_acv.discStep] * 3, self.tag_3d)
            if use_LabelLib and _LabelLib_found:
                self.ll_Grid3D = ll.Grid3D(grid_acv.shape, grid_acv.originXYZ, grid_acv.discStep)
                self.ll_Grid3D.grid = self.grid_1d
            else:
                self.ll_Grid3D = None
        else:
            self.cloud_xyzqt = cloud_xyzqt
        self.mp = Volume.mean_pos(self.cloud_xyzqt)

    def _reweight_cv(self, grid, cv_thickness, cv_fraction):
        """Reweight the accessible volume based on contact and free volume

        Parameters
        ----------
        cv_thickness : float
            width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
        cv_fraction : float
            relative fraction spent inside the contact volume (i.e. stacked to the surface)

        Returns
        -------
        grid_1d : ndarray
                  one-dimensional array of grid points of length n_gridpts

        Notes
        -----
        - A good first approximation of the CV thickness is about two times the smallest dye radius
        - The CV fraction can be calculated from the residual fluorescence anisotropy
        """
        grid_1d = np.array(grid)
        tag_1d = self._tag_volume(grid_1d)
        if cv_thickness > 0:
            weight_cv = self._weight_factor(grid_1d, cv_fraction)
            grid_1d[grid_1d > 1.0] = weight_cv
        grid_1d = np.clip(grid_1d, 0, None)
        return grid_1d, tag_1d

    @staticmethod
    def _weight_factor(grid_1d, cv_fraction):
        """Calculate the weight of contact volume grid points.

        Extended summary
        ----------------
        The weight factor corresponds to the factor by which the points of the contact volume
        (CV, dye trapped on surface) are favored over those belonging the free volume (FV, free dye diffusion).
        This accounts for the (experimentally) determined fraction of dyes populating the CV.

        Parameters
        ----------
        grid_1d : ndarray
                  one-dimensional array of grid points of length n_gridpts

        cv_fraction : float
                      fraction of dyes that are within the contact volume

        Returns
        -------
        weight_cv : float
                    weight of each grid point that belongs to the contact volume

        """
        n_CV = np.count_nonzero(grid_1d > 1.0)
        if n_CV != 0:
            n_FV = np.count_nonzero([(grid_1d > 0.0) & (grid_1d <= 1.0)])
            if n_FV != 0:
                weight_cv = n_FV / n_CV * cv_fraction / (1.0 - cv_fraction)
            else:
                print('no fraction of free dyes, all stacked')
                weight_cv = 2
        else:
            weight_cv = 1
        return weight_cv

    def _tag_volume(self, grid_1d):
        """Assign a tag to the grid values depending on their location in the cloud (1: free volume, 2: contact volume)

        Parameters
        ----------
        grid_1d : ndarray
                  one-dimensional array of grid points length n_gridpts

        Returns
        -------
        tag : ndarray
              one-dimensional array of length n_gridpts
        """
        tag_1d = np.full(self.n_gridpts, 1)
        tag_1d[grid_1d > 1.0] = 2
        return tag_1d


class FRET:
    """FRET class of distance and transfer efficiency for a pair of donor and acceptor ACVs


    Parameters
    ----------
    volume1, volume2 : instances of the Volume class
    fret_pair : str
        Distance key specifying the donor acceptor pair
    labels : dict
        dye, linker and setup parameters for the accessible volume calculation
    R_DA : ndarray
        donor acceptor distance distribution and associate weights (optional)
    verbose : bool

    Examples
    --------

    >>> ft.Molecule(volume1, volume2)

    >>> ft.Molecule(volume1, volume2, use_LabelLib=False)

    """

    def __init__(self, volume1, volume2, fret_pair, labels, R_DA=None, verbose=True):
        try:
            if volume1.acv is None or volume2.acv is None:
                raise TypeError
        except TypeError:
            if verbose:
                print('One accessible volume is empty')
        else:
            self.fret_pair = fret_pair
            self.R0 = labels['Distance'][fret_pair]["R0"]
            self.n_dist = labels['Distance'][fret_pair]["n_dist"]
            self.use_LabelLib = np.all([volume1.use_LabelLib, volume2.use_LabelLib])
            if self.use_LabelLib and _LabelLib_found:
                if R_DA is None:
                    R_DA = fret.dists_DA_ll(volume1.acv, volume2.acv, n_dist=self.n_dist, return_weights=True)
                self.mean_R_DA = fret.mean_dist_DA_ll(volume1.acv, volume2.acv, n_dist=self.n_dist)
                self.sigma_R_DA = fret.std_dist_DA(volume1.acv, volume2.acv, R_DA=R_DA)
                E_DA = fret.FRET_DA(volume1.acv, volume2.acv, R_DA=R_DA, R0=self.R0)
                self.mean_E_DA = fret.mean_FRET_DA_ll(volume1.acv, volume2.acv, R0=self.R0, n_dist=self.n_dist)
                self.sigma_E_DA = fret.std_FRET_DA(volume1.acv, volume2.acv, E_DA=E_DA)
            else:
                if R_DA is None:
                    R_DA = fret.dists_DA(volume1.acv, volume2.acv, n_dist=self.n_dist, return_weights=True)
                self.mean_R_DA = fret.mean_dist_DA(volume1.acv, volume2.acv, R_DA=R_DA)
                self.sigma_R_DA = fret.std_dist_DA(volume1.acv, volume2.acv, R_DA=R_DA)
                E_DA = fret.FRET_DA(volume1.acv, volume2.acv, R_DA=R_DA, R0=self.R0)
                self.mean_E_DA = fret.mean_FRET_DA(volume1.acv, volume2.acv, E_DA=E_DA)
                self.sigma_E_DA = fret.std_FRET_DA(volume1.acv, volume2.acv, E_DA=E_DA)
            self.mean_R_DA_E = fret.mean_dist_DA_fromFRET(volume1.acv, volume2.acv, mean_E_DA=self.mean_E_DA,
                                                          R0=self.R0)
            self.sigma_R_DA_E = fret.std_dist_DA_fromFRET(volume1.acv, volume2.acv, mean_E_DA=self.mean_E_DA,
                                                          sigma_E_DA=self.sigma_E_DA, R0=self.R0)
            self.R_attach = fret.dist_attach(volume1.attach_xyz, volume2.attach_xyz)
            self.R_mp = fret.dist_mp(volume1.acv, volume2.acv)

    @classmethod
    def from_volumes(cls, volume_list1, volume_list2, fret_pair, labels, R_DA=None):
        """Alternative constructor for the FRET class by reading in a list of donor and acceptor volumes

        Parameters
        ----------
        volume_list1, volume_list2 : array_like
            lists of Volume instances
        fret_pair : str
            identifier for donor acceptor pair
        labels : dict
            dye, linker and setup parameters for the accessible volume calculation
        R_DA : ndarray
            donor acceptor distance distribution and associate weights (optional)
        """
        n_vols1 = len(volume_list1)
        n_vols2 = len(volume_list2)
        try:
            if n_vols1 != n_vols2:
                raise ValueError
        except ValueError:
            print('The length of volume_list1 and volume_list2 is not the same')
        else:
            printProgressBar(0, n_vols1)
            fret = []
            for i in range(n_vols1):
                fret_value = FRET(volume_list1[i], volume_list2[i], fret_pair, labels, R_DA)
                if volume_list1[i].acv is None or volume_list2[i].acv is None:
                    print('Skip list entry {:d}'.format(i))
                else:
                    fret.append(fret_value)
                printProgressBar(i + 1, n_vols1)
            return fret

    def save_fret(self, filename):
        """Write the FRET calculation to a json file

        Parameters
        ----------
        filename : str
            filename for FRET prediction

        Examples
        --------

        >>> obj.save_FRET('parameters.json')
        """
        self.values.to_json(filename, indent=2)

    @property
    def values(self):
        """
        Pandas Dataframe of FRET parameters

        Returns
        -------
        pandas.DataFrame
        """
        fret_results = {'R0 (A)': (float(f'{self.R0 :0.1f}'), np.nan),
                        '<R_DA> (A)': (float(f'{self.mean_R_DA :0.1f}'), float(f'{self.sigma_R_DA :0.1f}')),
                        '<E_DA>': (float(f'{self.mean_E_DA :0.2f}'), float(f'{self.sigma_E_DA :0.2f}')),
                        '<R_DA_E> (A)': (float(f'{self.mean_R_DA_E :0.1f}'), float(f'{self.sigma_R_DA_E :0.1f}')),
                        'R_attach (A)': (float(f'{self.R_attach :0.1f}'), np.nan),
                        'R_mp (A)': (float(f'{self.R_mp :0.1f}'), np.nan)}
        fret_results = pd.DataFrame(fret_results, index=['value', 'std'])
        return fret_results


class Trajectory:
    """Trajectory class of distances and transfer efficiencies from the FRET class

    Parameters
    ----------
    fret : instance of the FRET class
    timestep : int
        time difference between two frames in picoseconds
    """
    def __init__(self, fret, timestep=None, kappasquare=None):
        n = len(fret)
        self.mean_E_DA = np.array([fret[i].mean_E_DA for i in range(n) if hasattr(fret[i], 'mean_E_DA')]).round(2)
        self.mean_R_DA = np.array([fret[i].mean_R_DA for i in range(n) if hasattr(fret[i], 'mean_R_DA')]).round(1)
        self.mean_R_DA_E = np.array([fret[i].mean_R_DA_E for i in range(n) if hasattr(fret[i], 'mean_R_DA_E')]).round(1)
        self.R_attach = np.array([fret[i].R_attach for i in range(n) if hasattr(fret[i], 'R_attach')]).round(1)
        self.R_mp = np.array([fret[i].R_mp for i in range(n) if hasattr(fret[i], 'R_mp')]).round(1)
        self.timestep = timestep
        self.kappasquare = kappasquare

    @property
    def dataframe(self):
        """Pandas Dataframe view of the Trajectory object

        Returns
        -------
        pandas.DataFrame
        """
        df = pd.DataFrame((self.mean_R_DA, self.mean_E_DA, self.mean_R_DA_E, self.R_attach, self.R_mp),
                          index=['<R_DA> (A)', '<E_DA>', '<R_DA_E> (A)', 'R_attach (A)', 'R_mp (A)']).T
        if self.timestep:
            df = pd.concat((df, pd.Series(range(df.shape[0]), name='time (ps)')*self.timestep), axis=1)
        if self.kappasquare:
            df = pd.concat((df, pd.Series(np.ones(df.shape[0]), name='kappasquare')*self.kappasquare), axis=1)
        return df

    def save_traj(self, filename, format='csv', units='A', header=True, R_kappa_only=False):
        """Save the trajectory as a .csv file

        Parameters
        ----------
        filename : str
            filename for .csv trajectory
        format : {'csv', 'txt'}
        units : {'A', 'nm'}
            distance units
        header : bool, default=True
            include header in the output file
        R_kappa_only : bool, default=False
            include only time, R_DA and kappasquare columns in the output file (as *gmx dyecoupl*, Gromacs)
        """
        df = self.dataframe
        with open(filename, 'w') as f:
            if format == 'txt':
                separator = '\t'
            else:
                separator = ','
            if units == 'nm':
                df['<R_DA> (A)'] = self.dataframe['<R_DA> (A)']/10
                df['<R_DA_E> (A)'] = self.dataframe['<R_DA_E> (A)']/10
                df['R_attach (A)'] = self.dataframe['R_attach (A)']/10
                df['R_mp (A)'] = self.dataframe['R_mp (A)']/10
                df.rename({header: header.replace('(A)', '(nm)') for header in df.columns.values},
                          axis='columns', inplace=True)
            if R_kappa_only:
                if self.timestep and self.kappasquare:
                    if units == 'nm':
                        columns = ['time (ps)', '<R_DA> (nm)', 'kappasquare']
                    else:
                        columns = ['time (ps)', '<R_DA> (A)', 'kappasquare']
                else:
                    raise KeyError('Time or kappasquare is missing in the DataFrame')
            else:
                columns = None
            f.write(df.to_csv(index=False, sep=separator, header=header, columns=columns))


class Volume:
    """Class object holding the accessible contact volume of a specific labeling position

    Parameters
    ----------
    structure : mdtraj.Trajectory
        trajectory of atom coordinates loaded from a pdb, xtc or other file
    site : str
        reference identifier for the labeling position
    labels : dict
        dye, linker and setup parameters for the accessible volume calculation
    cloud_xyzqt : ndarray
        array of x-,y-,z-coordinates with corresponding weights and tags with a shape [n_gridpts(+), 5]
    """
    def __init__(self, structure, site, labels, cloud_xyzqt=None, verbose=True):
        self.structure = structure
        self.attach_id = labels['Position'][site]['attach_id']

        if self.attach_id is not None:
            self.labeling_site = site
            self.mol_selection = labels['Position'][site]['mol_selection']
            self.linker_length = labels['Position'][site]['linker_length']
            self.linker_width = labels['Position'][site]['linker_width']
            self.simulation_type = labels['Position'][site]['simulation_type']
            self.dye_radii = np.array([labels['Position'][site]['dye_radius1'], labels['Position'][site]['dye_radius2'],
                                       labels['Position'][site]['dye_radius3']])
            self.grid_spacing = labels['Position'][site]['grid_spacing']
            self.cv_thickness = labels['Position'][site]['cv_thickness']
            self.cv_fraction = labels['Position'][site]['cv_fraction']
            self.state = labels['Position'][site]['state']
            self.frame_mdtraj = labels['Position'][site]['frame_mdtraj']
            self.use_LabelLib = labels['Position'][site]['use_LabelLib']

            try:
                self.n_atoms = self.structure.n_atoms
            except:
                print('{} is no valid mdtraj.Trajectory object'.format(self.structure))
            else:
                try:
                    if (self.frame_mdtraj != self.state - 1):
                        raise ValueError
                except ValueError:
                    print('The state {:d} and frame_mdtraj {:d} are not compatible. \
                    The frame_mdtraj should be equal to state - 1'.format(self.state, self.frame_mdtraj))
                    self.av = None
                    self.acv = None
                else:
                    try:
                        if self.state > self.structure.n_frames:
                            raise IndexError
                    except IndexError:
                        print('The state {:d} and mdtraj frame {:d} are out of range, \
                            select a state within the range 1 - {:d} \
                            or a mdtraj frame within the range 0 - {:d}'.format(self.state, self.frame_mdtraj,
                                                                                self.structure.n_frames,
                                                                                self.structure.n_frames - 1))
                        self.av = None
                        self.acv = None
                    else:
                        try:
                            if self.attach_id < 1 or self.attach_id > self.n_atoms:
                                raise IndexError
                        except IndexError:
                            print('The attachment position {:d} is out of range, \
                                   select an index within 1 - {:d}'.format(self.attach_id, self.n_atoms))
                            self.av = None
                            self.acv = None
                        else:
                            self.attach_id_mdtraj = labels['Position'][site]['attach_id'] - 1
                            self.resi_atom = self.structure.top.atom(self.attach_id_mdtraj)

                            self.av = self.calc_av(self.use_LabelLib)

                            try:
                                if not np.any(np.array(self.av.grid) > 0):
                                    raise ValueError
                            except ValueError:
                                if verbose:
                                    print('Empty Accessible volume at position {:d}. \
                                    Is your attachment point buried?'.format(self.attach_id))
                                self.acv = None
                            else:
                                self.acv = self.calc_acv(self.use_LabelLib)

        elif cloud_xyzqt is not None:
            self.acv = ACV(cloud_xyzqt=cloud_xyzqt)
        else:
            print('Attachment point is unknown')

    @classmethod
    def from_frames(cls, structure, site, labels, frames_mdtraj):
        """Alternative constructor for the Volume class by reading in one or multiple
        frames using the same dye and grid parameters

        Parameters
        ----------
        structure : mdtraj.Trajectory
            trajectory of atom coordinates loaded from a pdb, xtc or other file
        site : str
            reference identifier for the labeling position
        labels : dict
            dye, linker and setup parameters for the accessible volume calculation
        frames_mdtraj : int or list
            list of frames on the trajectory to be used for the ACV calculation

        Returns
        -------
        list
            list of Volume instances
        """
        n_fr = len(frames_mdtraj)
        printProgressBar(0, n_fr)
        multiframe_volumes = []
        _labels = copy.copy(labels)
        for i, frame in enumerate(frames_mdtraj):
            _labels['Position'][site]['frame_mdtraj'] = frame
            _labels['Position'][site]['state'] = frame + 1
            multiframe_volumes.append(cls(structure, site, _labels))
            printProgressBar(i + 1, n_fr)
        return multiframe_volumes

    @classmethod
    def from_attachID(cls, structure, site, labels, attachID):
        """Alternative constructor for the Volume class by reading in one or multiple
        attachment points using the same dye and grid parameters

        Parameters
        ----------
        structure : mdtraj.Trajectory
            trajectory of atom coordinates loaded from a pdb, xtc or other file
        site : str
            reference identifier for the labeling position
        labels : dict
            dye, linker and setup parameters for the accessible volume calculation
        attachID : int or list
            list of attachment ids on the structure to be used for the ACV calculation

        Returns
        -------
        list
            list of Volume instances

        Examples
        --------

        >>> struct = md.load_pdb('data/2m23.pdb')
        >>> labels_global = ft.labeling_params('data/labels_global.json')
        >>> from_attachID(struct, 0, labelsglobals, [2,929])

        """
        n_aID = len(attachID)
        printProgressBar(0, n_aID)
        multisite_volumes = []
        _labels = copy.copy(labels)
        for i, attach_id in enumerate(attachID):
            _labels['Position'][site]['attach_id'] = attach_id
            multisite_volumes.append(cls(structure, site, _labels))
            printProgressBar(i + 1, n_aID)
        return multisite_volumes

    @staticmethod
    def reshape_grid(grid, shape):
        """Convert a 1D-grid to a 3D-Grid

        Parameters
        ----------
        grid : array_like
            one-dimensional grid array/list of length n_gridpts
        shape : list

        Returns
        -------
        ndarray
            3-dimensional array of grid points with a shape given by n_xyz

        Examples
        --------

        >>> grid = numpy.full(27,1)
        >>> ft.cloud.Volume.reshape_grid(grid, (3,3,3))
        array([[[1., 1., 1.],
                [1., 1., 1.],
                [1., 1., 1.]],
                 ...
               [[1., 1., 1.],
                [1., 1., 1.],
                [1., 1., 1.]]])

        """
        grid_3d = np.array(grid, dtype=np.float64).reshape(shape, order='F')
        return grid_3d

    @staticmethod
    @nb.jit(forceobj=True)
    def grid2pts(grid_3d, xyz_min, d_xyz, *args):
        """Convert 3D-grid with density values to xyz coordinates with a weight (q)

        Parameters
        ----------
        grid_3d : ndarray([nx,ny,nz])
            3-dimensional array of grid points with a shape given by n_xyz
        xyz_min : list
            origin coordinates of the grid
        d_xyz : list
            grid spacing in x-,y- and z-direction

        Returns
        -------
        ndarray
            array of x-,y-,z-coordinates with corresponding weights and tags with a shape [n_gridpts(+), 5]

        Examples
        --------

        >>> grid_3d = numpy.zeros((3,3,3))
        >>> grid_3d[(1,1,1)] = 1   # FV
        >>> grid_3d[(1,2,1)] = 10  # CV
        >>> ft.cloud.Volume.grid2pts(grid_3d, [0,0,0], [1,1,1])
        array([[ 0.5,  0.5,  0.5,  1.  1],
               [ 0.5,  1. ,  0.5, 10.  1]])

        >>> tag_3d = numpy.zeros((3,3,3))
        >>> grid_3d[(1,1,1)] = 1   # FV
        >>> grid_3d[(1,1,1)] = 2   # CV
        >>> ft.cloud.Volume.grid2pts(grid_3d, [0,0,0], [1,1,1], tag_3d)
        array([[ 0.5,  0.5,  0.5,  1.  1],
               [ 0.5,  1. ,  0.5, 10.  2]])

        """
        xmin, ymin, zmin = xyz_min
        nx, ny, nz = grid_3d.shape
        dx, dy, dz = d_xyz

        if len(args) != 0:
            tag_3d = args[0]
        else:
            tag_3d = np.full(grid_3d.shape, 1)

        cloud_xyzqt = np.empty((nx * ny * nz, 5), dtype=np.float64, order='C')

        gdx = np.arange(0, nx, dtype=np.float64) * dx
        gdy = np.arange(0, ny, dtype=np.float64) * dy
        gdz = np.arange(0, nz, dtype=np.float64) * dz

        n = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    d = grid_3d[ix, iy, iz]
                    if d > 0:
                        cloud_xyzqt[n, 0] = gdx[ix] + xmin
                        cloud_xyzqt[n, 1] = gdy[iy] + ymin
                        cloud_xyzqt[n, 2] = gdz[iz] + zmin
                        cloud_xyzqt[n, 3] = d
                        cloud_xyzqt[n, 4] = tag_3d[ix, iy, iz]
                        n += 1

        return cloud_xyzqt[:n]

    @property
    def mol_xyzr(self):
        """Get coordinates and vdW radii of all selected atoms in the structure

        Returns
        -------
        ndarray
            array of x-,y-,z-coordinates and VdW radii with a shape [n_atoms, 4]

        Examples
        --------

        >>> avobj.mol_xyzr

        """
        frame_mdtraj = self.frame_mdtraj
        struct = self.structure
        try:
            radii = np.array([VDW_RADIUS[atom.element.symbol] / 100 for atom in struct.top.atoms], ndmin=2).T
        except KeyError as e:
            print(f'Atom symbol {e} in topology not recognized')
            raise
        else:
            radii[self.attach_id_mdtraj] = 0.0
            try:
                sele = struct.top.select(self.mol_selection)
            except ValueError:
                print('{} is not a valid expression, please see http://mdtraj.org/latest/atom_selection.html'.format(
                    self.mol_selection))
                print('Falling back to \"all\"')
                sele = struct.top.select('all')
            finally:
                if self.attach_id_mdtraj not in sele:
                    np.append(sele, self.attach_id_mdtraj)
                xyz = struct.xyz[frame_mdtraj] * 10
                xyzr = np.hstack((xyz[sele], radii[sele]))
                return xyzr

    @property
    def attach_xyz(self):
        """Get coordinates of the dye attachment point

        Returns
        -------
        ndarray
            one-dimensional array of x-,y-,z-coordinates of the attachment point

        Examples
        --------

        >>> avobj.attach_xyz

        """
        xyz = self.structure.xyz[self.frame_mdtraj][self.attach_id_mdtraj] * 10
        return xyz.astype(np.float64)

    def save_acv(self, filename, format='pdb', **kwargs):
        """Write accessible contact volume to file

        Parameters
        ----------
        filename : str
            filename for ACV
        format : {'pdb', 'xyz', 'openDX'}
            file format for ACV
        **kwargs : dict
            key value pairs for .xyz files are {write_weights : bool, encode_element : bool}

        Notes
        -----
        - `write_weights` includes the weights of each point as an additional column in the .xyz file
        - `encode_element` uses the element column (first col.) to encode the FV (=F) or CV (=C)

        Examples
        --------

        >>> obj.save_acv('cloud.xyz', format='xyz', write_weights=False)

        """
        if format == 'xyz':
            try:
                write_weights = kwargs['write_weights']
            except KeyError:
                write_weights = True
            try:
                encode_element = kwargs['encode_element']
            except KeyError:
                encode_element = False
            file_str = export.xyz(self.acv.cloud_xyzqt, self.acv.mp, write_weights, encode_element)
        elif format == 'open_dx':
            d_xyz = [self.grid_spacing] * 3
            xyz_min = self.av.originXYZ
            file_str = export.open_dx(self.acv.grid_3d, xyz_min, d_xyz)
        else:
            file_str = export.pdb(self.acv.cloud_xyzqt, self.acv.mp)
        with open(filename, 'w') as fname:
            fname.write(file_str)

    def calc_av(self, use_LabelLib):
        """Calculate the dye accessible volume

        Returns
        -------
        LabelLib.Grid3D
            object containing a 3-dimensional array of grid points and additional attributes (see Notes)

        Notes
        -----
        The attributes of the av LabelLib.Grid3D object are:
            - discStep : float  (the grid spacing)
            - originXYZ : list  (x-/y-/z-coordinate of the grid origin)
            - shape : list      (number of grid points in x-/y-/z-direction)
            - grid : list       (flattened list of grid point values)

        See Also
        --------
        calc_acv : Calculate accessible contact volume

        References
        ----------
        .. [3] Kalinin, S. et al. "A toolkit and benchmark study for FRET-restrained high-precision \
        structural modeling", *Nat. Methods* **9**, 1218–1225 (2012).
        .. [4] Sindbert, S. et al. "Accurate distance determination of nucleic acids via Förster \
        resonance energy transfer: implications of dye linker length and rigidity", \
        *J. Am. Chem. Soc.* **133**, 2463–2480 (2011).

        Examples
        --------

        >>> avobj.calc_av()

        """
        mol_xyzr = self.mol_xyzr
        attach_xyz = self.attach_xyz

        if use_LabelLib and _LabelLib_found:
            if self.simulation_type == 'AV1':
                av = ll.dyeDensityAV1(mol_xyzr.T, attach_xyz, self.linker_length,
                                      self.linker_width, self.dye_radii[0], self.grid_spacing)
            else:
                av = ll.dyeDensityAV3(mol_xyzr.T, attach_xyz, self.linker_length,
                                      self.linker_width, self.dye_radii, self.grid_spacing)
        else:
            av = grid.Grid3D(self.mol_xyzr, self.attach_xyz, self.linker_length, self.linker_width, self.dye_radii,
                             self.grid_spacing, self.simulation_type)
        return av

    def _dye_acc_surf(self):
        """Calculate dye accessible surface by padding the vdW radius with the thickness of the contact volume

        Returns
        -------
        ndarray
            array with dimensions (5,n_atoms) of marked coordinates and padded vdW radii
            (n_atoms = number of atoms in mdtraj.Trajectory)
        """
        cv_label = 2.0
        das_marker = np.full(self.mol_xyzr.shape[0], cv_label)
        das_xyzrm = np.vstack([self.mol_xyzr.T, das_marker])
        das_xyzrm[3] += self.cv_thickness
        return das_xyzrm

    def calc_acv(self, use_LabelLib):
        """Partition the accessible volume into a free and a contact volume

        Returns
        -------
        acv : fretraj.cloud.ACV

        Notes
        -----
        The attributes from the LabelLib.Grid3D object are:
            - discStep : float  (the grid spacing)
            - originXYZ : list  (x-/y-/z-coordinate of the grid origin)
            - shape : list      (number of grid points in x-/y-/z-direction)
            - grid : list       (flattened list of grid point values)

        Additional attributes are:
            - grid_3d : ndarray([shape])
            - n_gridpts : int
            - tag : ndarray([1,n_gridpts])

        See Also
        --------
        calc_av : Calculate accessible volume

        References
        ----------
        .. [1] Steffen, F. D., Sigel, R. K. O. & Börner, R. "An atomistic view on carbocyanine \
        photophysics in the realm of RNA", *Phys. Chem. Chem. Phys.* **18**, 29045–29055 (2016).
        .. [2] Dimura, M. et al. Quantitative FRET studies and integrative modeling unravel the structure \
        and dynamics of biomolecular systems", *Curr. Opin. Struct. Biol.* **40**, 163–185 (2016).

        Examples
        --------

        >>> avobj.calc_acv()

        """
        das_xyzrm = self._dye_acc_surf()
        if use_LabelLib and _LabelLib_found:
            grid_acv = ll.addWeights(self.av, das_xyzrm)
        else:
            self.av.grid_3d = self.av.addWeights(das_xyzrm.T)
            self.av.grid = self.av.grid_3d.flatten(order='F')
            grid_acv = self.av
        acv = ACV(grid_acv, self.cv_thickness, self.cv_fraction, use_LabelLib=use_LabelLib)
        return acv

    @staticmethod
    def mean_pos(cloud_xyzqt):
        """Calculate mean dye position of the accessible contact volume

        Parameters
        ----------
        cloud_xyzq : ndarray
            array of x-,y-,z-coordinates and corresponding weights with a shape [n_gridpts(+), 4]

        Returns
        -------
        ndarray
             mean position
        """
        x = np.dot(cloud_xyzqt[:, 0], cloud_xyzqt[:, 3])
        y = np.dot(cloud_xyzqt[:, 1], cloud_xyzqt[:, 3])
        z = np.dot(cloud_xyzqt[:, 2], cloud_xyzqt[:, 3])
        mp = np.array((x, y, z)) / cloud_xyzqt[:, 3].sum()
        return mp

    @staticmethod
    def _weighted_quantile_1D(arr_1D, weights, quantile):
        """Compute the weighted quantile of a 1D-array

        Parameters
        ----------
        arr_1D : ndarray
            one-dimensional array
        weights : ndarray
            weight of each element of the array
        quantile : float
            quantile to compute (median=0.5)

        Returns
        -------
        quantile : float
        """
        if ((quantile > 1) or (quantile < 0)):
            raise ValueError('Quantiles must be in the range [0, 1]')

        idx_sorted = np.argsort(arr_1D)
        arr_1D_sorted = arr_1D[idx_sorted]
        weights_sorted = weights[idx_sorted]
        weights_cs = np.cumsum(weights_sorted)
        xp = (weights_cs - (1-quantile)*weights_sorted) / weights_cs[-1]
        return np.interp(quantile, xp, arr_1D_sorted)

    def save_mp(self, filename, format='plain', units='A'):
        """Write mean dye position to a file

        Parameters
        ----------
        filename : str
            filename for the mean position
        format : {'plain', 'json'}
            file format (plain .txt or .json)
        units : {'A', 'nm'}
            distance units ('A': Angstroms, 'nm': nanometers)

        Examples
        --------
        >>> avobj.save_mp('mp.dat')

        """
        if units == 'nm':
            mp = self.acv.mp / 10
        else:
            mp = self.acv.mp
        with open(filename, 'w') as f:
            if format == 'json':
                mp_dict = {'x': round(mp[0], 3), 'y': round(mp[1], 3), 'z': round(mp[2], 3)}
                json.dump(mp_dict, f)
            else:
                f.write('{:0.3f}   {:0.3f}   {:0.3f}\n'.format(*mp))


def main():
    in_filePDB, param_fileJSON, out_fileACV = parseCmd()
    labels = labeling_params(param_fileJSON)
    struct = md.load_pdb(in_filePDB)
    av = Volume(struct, 0, '40', labels)

    if out_fileACV:
        _, out_fileextACV = os.path.splitext(out_fileACV)
        av.save_acv(out_fileACV, format=out_fileextACV[1:])


if __name__ == "__main__":
    main()
