#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import argparse
import json
import mdtraj as md
import numba as nb
import copy

try:
    import LabelLib as ll
except ModuleNotFoundError:
    print('\nNote: LabelLib module is not found. \nACV calculation uses a Python-only algorithm\n')
    _LabelLib_found = False
else:
    _LabelLib_found = True

from fretraj import export
from fretraj import fret
from fretraj import grid
from fretraj import metadata

DISTANCE_SAMPLES = 200000

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
                                       'use_LabelLib': (bool, True),
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
    """
    Parse the command line to get the input PDB file and the dye parameters

    Returns
    -------
    pdbFile : string
    paramFile : string
    """
    parser = argparse.ArgumentParser(
        description='compute accessible-contact clouds for \
                     an MD trajectory or a given PDB structure')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + str(metadata['Version']))
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
    return in_filePDB, param_fileJSON, out_fileACV


def labeling_params(param_file):
    """
    Parse the parameter file to get the configuration settings

    Parameters
    ----------
    param_file : str

    Returns
    -------
    labels : dict
             position of labels, dye and grid parameters
    """
    with open(param_file) as f:
        labels_json = json.load(f)

    try:
        check_labels(labels_json)
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
    """
    Check integrity of parameter dictionary

    Parameters
    ----------
    labels : dict
    verbose : bool

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


def save_labels(filename, labels):
    """
    Write the ACV parameters to a .json file

    Parameters
    ----------
    filename : str
    labels : dict
             position of labels, dye and grid parameters

    Examples
    --------

    >>> obj.save_labels('parameters.json')
    """
    with open(filename, 'w') as f:
        json.dump(labels, f, indent=2)


def printProgressBar(iteration, total, prefix='Progress:', suffix='complete', length=20, fill='█'):
    """
    Command line progress bar

    Parameters
    ----------
    iteration : int
                current iteration
    total : int
            total number of iterations
    prefix : str, optional
             string before progress bar
    suffix : str, optional
             string after progress bar
    length : int
             length of progress bar
    fill : str, optional
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


def save_mp_traj(filename, volume_list, units='A'):
    """
    Save a trajectory of dye mean positions as an xyz file

    Parameters
    ----------
    filename : str
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


def save_acv_traj(filename, volume_list, **kwargs):
    """
    Save a trajectory of ACVs as a multi model PDB

    Parameters
    ----------
    filename : str
    volume_list : array_like
                  list of Volume instances
    **kwargs
            - include_mdp : bool
    """
    try:
        include_mdp = kwargs['include_mdp']
    except KeyError:
        include_mdp = False
    file_str = ''
    for i, volume in enumerate(volume_list):
        file_str += f'MODEL {i+1}\n'
        file_str += export.pdb(volume.acv.cloud_xyzqt, volume.acv.mp, volume.acv.mdp, include_mdp=include_mdp)
        file_str += 'ENDMDL\n\n'
    with open(filename, 'w') as fname:
        fname.write(file_str)


def save_structure_traj(filename, structure, frames, format='pdb'):
    """
    Save selected frames of a trajectory

    Parameters
    ----------
    filename : str
    structure : mdtraj.Trajectory
                trajectory of atom coordinates loaded from a pdb, xtc or other file
    format : str
             trajectory file format. One of the following:
             'pdb': multi model PDB (default), 'xtc'
    """
    sliced_structure = structure.slice(frames)
    if format == 'xtc':
        sliced_structure.save_xtc(filename)
    else:
        sliced_structure.save_pdb(filename)


class ACV:
    """
    Reweighted accessible contact volume of a covalently linked fluorophore

    Parameters
    ----------
    grid_acv : LabelLib.Grid3D
             accessible volume with added weight labels for free and contact volume
             (density label FV: 1.0, density label CV: 2.0); the two volumes are subsequently reweighted
    cv_thickness : float
                   width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
                   **Note:** good first approx.: 2*min(dye_radii)
    cv_fraction : float
                  fraction of dyes that are within the contact volume
                  (e.g. as determined by fluorescence anisotropy)
    cloud_xyzq : numpy.ndarray
                 array of x-,y-,z-coordinates and corresponding weights
                 with a shape [n_gridpts(+), 4]
    use_LabelLib : bool
                   make use of LabelLib library to compute FRET values and distances

    Attributes
    ----------
    grid : list
           flattened list of grid point values of the LabelLib.Grid3D class with labeled CV points (not yet reweighted)
    shape : list
            shape of the LabelLib.Grid3D
    originXYZ : list
                origin of the LabelLib.Grid3D
    discStep : float
               spacing of the LabelLib.Grid3D

    n_gridpts : int
                total number of grid points

    grid_1d : numpy.ndarray
              one-dimensional array of grid points of length n_gridpts
    grid_3d : numpy.ndarray
              3-dimensional array of grid points with a shape given by n_xyz
    tag_3d : numpy.ndarray
             one-dimensional array of length n_gridpts
    cloud_xyzq : numpy.ndarray
                 array of x-,y-,z-coordinates and corresponding weights
                 with a shape [n_gridpts(+), 4]
    ll_Grid3D : LabelLib.Grid3D
                reweighted LabelLib.Grid3D class with attributes (shape, originXYZ, discStep, grid)
    """

    def __init__(self, grid_acv=None, cv_thickness=0, cv_fraction=0, cloud_xyzqt=None, use_LabelLib=True):
        if grid_acv is not None:
            self.grid = grid_acv.grid
            self.shape = grid_acv.shape
            self.originXYZ = grid_acv.originXYZ
            self.discStep = grid_acv.discStep

            self.n_gridpts = np.prod(self.shape)
            self.grid_1d, self.tag_1d = self._reweight_cv(cv_thickness, cv_fraction)
            self.grid_3d = Volume.reshape_grid(self.grid_1d, self.shape)
            self.tag_3d = Volume.reshape_grid(self.tag_1d, self.shape)
            self.cloud_xyzqt = Volume.grid2pts(self.grid_3d, self.originXYZ, [self.discStep] * 3, self.tag_3d)
            if use_LabelLib and _LabelLib_found:
                self.ll_Grid3D = ll.Grid3D(self.shape, self.originXYZ, self.discStep)
                self.ll_Grid3D.grid = self.grid_1d
            else:
                self.ll_Grid3D = None
        else:
            self.cloud_xyzqt = cloud_xyzqt
        self.mp = Volume.mean_pos(self.cloud_xyzqt)
        self.mdp = Volume.median_pos(self.cloud_xyzqt)

    def _reweight_cv(self, cv_thickness, cv_fraction):
        """
        Reweight the accessible volume based on contact and free volume

        Parameters
        ----------
        cv_thickness : float
                   width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
                   **Note:** good first approx.: 2*min(dye_radii)
        cv_fraction : float
                  fraction of dyes that are within the contact volume
                  (e.g. as determined by fluorescence anisotropy)

        Returns
        -------
        grid_1d : ndarray
                  one-dimensional array of grid points of length n_gridpts
        """
        grid_1d = np.array(self.grid)
        tag_1d = self._tag_volume(grid_1d)
        if cv_thickness > 0:
            weight_cv = self._weight_factor(grid_1d, cv_fraction)
            grid_1d[grid_1d > 1.0] = weight_cv
        grid_1d = np.clip(grid_1d, 0, None)
        return grid_1d, tag_1d

    @staticmethod
    def _weight_factor(grid_1d, cv_fraction):
        """
        Calculate the weight of contact volume grid points.

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
        """
        Assign a tag to the grid values depending on their location in the cloud
        (1: free volume, 2: contact volume)

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
    """
    Parameters
    ----------
    volume1 : instance of the Volume class
    volume2 : instance of the Volume class
    fret_pair : str
                Distance key specifying the donor acceptor pair
    labels : dict
             dye, linker and setup parameters for the accessible volume calculation
    R_DA : ndarray
           donor acceptor distance distribution and associate weights (optional)
    verbose : bool

    Attributes
    ----------
    volume1 : instance of the Volume class
    volume2 : instance of the Volume class
    R_DA : ndarray
           donor acceptor distance distribution in (A) and associate weights (optional)
    use_LabelLib : bool
                   make use of LabelLib library to compute FRET values and distances

    R_DA : ndarray
    mean_R_DA : float
                mean FRET
    mean_E_DA : float
    mean_R_DA_E : float


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
            self.volume1 = volume1
            self.volume2 = volume2
            self.fret_pair = fret_pair
            self.R0 = labels['Distance'][fret_pair]["R0"]
            self.n_dist = labels['Distance'][fret_pair]["n_dist"]
            self.use_LabelLib = np.all([volume1.use_LabelLib, volume2.use_LabelLib])
            if self.use_LabelLib and _LabelLib_found:
                if R_DA is None:
                    self.R_DA = fret.dists_DA_ll(volume1.acv, volume2.acv, n_dist=self.n_dist, return_weights=True)
                else:
                    self.R_DA = R_DA
                self.mean_R_DA = fret.mean_dist_DA_ll(volume1.acv, volume2.acv, n_dist=self.n_dist)
                self.sigma_R_DA = fret.std_dist_DA(volume1.acv, volume2.acv, R_DA=self.R_DA)
                self.E_DA = fret.FRET_DA(volume1.acv, volume2.acv, R_DA=self.R_DA, R0=self.R0)
                self.mean_E_DA = fret.mean_FRET_DA_ll(volume1.acv, volume2.acv, R0=self.R0, n_dist=self.n_dist)
                self.sigma_E_DA = fret.std_FRET_DA(volume1.acv, volume2.acv, E_DA=self.E_DA)
            else:
                if R_DA is None:
                    self.R_DA = fret.dists_DA(volume1.acv, volume2.acv, n_dist=self.n_dist, return_weights=True)
                else:
                    self.R_DA = R_DA
                self.mean_R_DA = fret.mean_dist_DA(volume1.acv, volume2.acv, R_DA=self.R_DA)
                self.sigma_R_DA = fret.std_dist_DA(volume1.acv, volume2.acv, R_DA=self.R_DA)
                self.E_DA = fret.FRET_DA(volume1.acv, volume2.acv, R_DA=self.R_DA, R0=self.R0)
                self.mean_E_DA = fret.mean_FRET_DA(volume1.acv, volume2.acv, E_DA=self.E_DA)
                self.sigma_E_DA = fret.std_FRET_DA(volume1.acv, volume2.acv, E_DA=self.E_DA)
            self.mean_R_DA_E = fret.mean_dist_DA_fromFRET(volume1.acv, volume2.acv, mean_E_DA=self.mean_E_DA,
                                                          R0=self.R0)
            self.sigma_R_DA_E = fret.std_dist_DA_fromFRET(volume1.acv, volume2.acv, mean_E_DA=self.mean_E_DA,
                                                          sigma_E_DA=self.sigma_E_DA, R0=self.R0)
            self.R_attach = fret.dist_attach(volume1.attach_xyz, volume2.attach_xyz)
            self.R_mp = fret.dist_mp(volume1.acv, volume2.acv)

    @classmethod
    def from_volumes(cls, volume_list1, volume_list2, fret_pair, labels, R_DA=None):
        """
        Alternative constructor for the ft.cloud.FRET class by reading in a list of donor and acceptor volumes

        Parameters
        ----------
        volume_list1 : array_like
                       list of Volume instances
        volume_list2 : array_like
                       list of Volume instances
        fret_pair : str
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
            fret_trajectory = []
            for i in range(n_vols1):
                fret_value = FRET(volume_list1[i], volume_list2[i], fret_pair, labels, R_DA)
                if volume_list1[i].acv is None or volume_list2[i].acv is None:
                    print('Skip list entry {:d}'.format(i))
                else:
                    fret_trajectory.append(fret_value)
                printProgressBar(i + 1, n_vols1)
            return fret_trajectory

    def save_fret(self, filename):
        """
        Write the FRET calculation to a json file

        Parameters
        ----------
        filename : str

        Examples
        --------

        >>> obj.save_FRET('parameters.json')
        """
        fret_results = self.values
        with open(filename, 'w') as f:
            json.dump(fret_results, f, indent=2)

    @property
    def values(self):
        """
        Returns a dictionary of FRET parameters
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
    """
    Parameters
    ----------
    fret : instance of the FRET class
    timestep : int
               time difference between two frames in picoseconds
    """
    def __init__(self, fret, timestep=None):
        n = len(fret)
        self.mean_E_DA = np.array([fret[i].mean_E_DA for i in range(n) if hasattr(fret[i], 'mean_E_DA')]).round(2)
        self.mean_R_DA = np.array([fret[i].mean_R_DA for i in range(n) if hasattr(fret[i], 'mean_R_DA')]).round(1)
        self.mean_R_DA_E = np.array([fret[i].mean_R_DA_E for i in range(n) if hasattr(fret[i], 'mean_R_DA_E')]).round(1)
        self.R_attach = np.array([fret[i].R_attach for i in range(n) if hasattr(fret[i], 'R_attach')]).round(1)
        self.R_mp = np.array([fret[i].R_mp for i in range(n) if hasattr(fret[i], 'R_mp')]).round(1)
        self.timestep = timestep

    @property
    def dataframe(self):
        """
        Dataframe view of the Trajectory object

        Returns
        -------
        df : pandas dataframe
        """
        df = pd.DataFrame((self.mean_R_DA, self.mean_E_DA, self.mean_R_DA_E, self.R_attach, self.R_mp),
                          index=['<R_DA> (A)', '<E_DA>', '<R_DA_E> (A)', 'R_attach (A)', 'R_mp (A)']).T
        if self.timestep:
            df = pd.concat((df, pd.Series(range(df.shape[0]), name='time (ps)')*self.timestep), axis=1)
        return df

    def save_traj(self, filename):
        """
        Save the trajectory as a .csv file

        Parameters
        ----------
        filename : str
        """
        with open(filename, 'w') as f:
            f.write(self.dataframe.to_csv(index=False))


class Volume:
    """
    Class object holding the accessible contact volume of a specific labeling position

    Parameters
    ----------
    structure : mdtraj.Trajectory
                trajectory of atom coordinates loaded from a pdb, xtc or other file
    site : str
           reference identifier for the labeling position
    labels : dict
             dye, linker and setup parameters for the accessible volume calculation

    Attributes
    ----------
    labeling_site : str
                    reference identifier for the labeling position
    structure : mdtraj.Trajectory
                trajectory of atom coordinates loaded from a pdb, xtc or other file
    state : int
            state in the pdb file (1 based indexing)
    frame_mdtraj : int
                   frame number of the mdtraj.Trajectory object (0 based indexing)
    attach_id : int
                serial atom id of the attachment point in the pdb file (1 based indexing)
    attach_id_mdtraj : int
                serial atom id of the attachment point in the mdtraj.Trajectory object (0 based indexing)
    resi_atom : str
                residue name, residue number and atom name of the attachment point
    mol_selection : str
                    atom selection expression compatible with `MDTraj <http://mdtraj.org/latest/atom_selection.html>`_
    linker_length : float
                    length of the dye linker in Angstrom
    linker_width : float
                   diameter of the dye linker in Angstrom
    simulation_type : {'AV1', 'AV3'}
                      type of accessible volume calculation
                      'AV1' parameterizes the dye as sphere with radius `dye_radii[0]`
                      'AV3' parameterizes the dye as an ellipsoid with three radii `dye_radii`
    dye_radii : ndarray([3,1])
                array of dye radii in Angstrom with shape [3,1]
    grid_spacing : float
    cv_thickness : float
                   width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
                   **Note:** good first approx.: 2*min(dye_radii)
    cv_fraction : float
                  fraction of dyes that are within the contact volume
                  (e.g. as determined by fluorescence anisotropy)
    frame : int
            frame number of the mdtraj.Trajectory object (0 based indexing)
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
        """
        Alternative constructor for the ft.cloud.Volume class by reading in one or multiple
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
        """
        Alternative constructor for the ft.cloud.Volume class by reading in one or multiple
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
        """
        Convert a 1D-grid to a 3D-Grid

        Parameters
        ----------
        grid : array_like
               one-dimensional grid array/list of length n_gridpts
        shape : list

        Returns
        -------
        grid_3d : ndarray
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
        """
        Convert 3D-grid with density values to xyz coordinates with a weight (q)

        Parameters
        ----------
        grid_3d : numpy.ndarray([nx,ny,nz])
                  3-dimensional array of grid points with a shape given by n_xyz
        xyz_min : list
                  origin coordinates of the grid
        d_xyz : list
                grid spacing in x-,y- and z-direction
        Returns
        -------
        cloud_xyzqt : ndarray
                     array of x-,y-,z-coordinates with corresponding weights and tags
                     with a shape [n_gridpts(+), 5]

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
        """
        Get coordinates and vdW radii of all selected atoms in the structure

        Returns
        -------
        xyzr : numpy.ndarray
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
        """
        Get coordinates of the dye attachment point

        Returns
        -------
        xyz : ndarray
              one-dimensional array of x-,y-,z-coordinates of the attachment point

        Examples
        --------

        >>> avobj.attach_xyz

        """
        xyz = self.structure.xyz[self.frame_mdtraj][self.attach_id_mdtraj] * 10
        return xyz.astype(np.float64)

    def save_acv(self, filename, format='xyz', **kwargs):
        """
        Write accessible contact volume to file

        Parameters
        ----------
        filename : str
        format : str
                 default: xyz
        **kwargs
            - write_weights : bool
            - encode_element : bool

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
            try:
                include_mdp = kwargs['include_mdp']
            except KeyError:
                include_mdp = False
            file_str = export.xyz(self.acv.cloud_xyzqt, self.acv.mp, self.acv.mdp, write_weights, encode_element,
                                  include_mdp)
        elif format == 'open_dx':
            d_xyz = [self.grid_spacing] * 3
            xyz_min = self.acv.originXYZ
            file_str = export.open_dx(self.acv.grid_3d, xyz_min, d_xyz)
        else:
            try:
                include_mdp = kwargs['include_mdp']
            except KeyError:
                include_mdp = False
            file_str = export.pdb(self.acv.cloud_xyzqt, self.acv.mp, self.acv.mdp, include_mdp=include_mdp)
        with open(filename, 'w') as fname:
            fname.write(file_str)

    def calc_av(self, use_LabelLib):
        """
        Calculate the dye accessible volume [#]_ [#]_

        Returns
        -------
        av : LabelLib.Grid3D
             the attributes of the av LabelLib.Grid3D object are:
                - discStep : float  (the grid spacing)
                - originXYZ : list  (x-/y-/z-coordinate of the grid origin)
                - shape : list      (number of grid points in x-/y-/z-direction)
                - grid : list       (flattened list of grid point values)

        See Also
        --------
        calc_acv : Calculate accessible contact volume

        References
        ----------
        .. [#] Kalinin, S. et al. "A toolkit and benchmark study for FRET-restrained high-precision \
        structural modeling", *Nat. Methods* **9**, 1218–1225 (2012).
        .. [#] Sindbert, S. et al. "Accurate distance determination of nucleic acids via Förster \
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
        """
        Calculate dye accessible surface by padding the vdW radius with the thickness of the contact volume

        Returns
        ------
-       das_xyzrm : numpy.ndarray([5,n_atoms])
                    array of marked coordinates and padded vdW radii (n_atoms = number of atoms in mdtraj.Trajectory)
        """
        cv_label = 2.0
        das_marker = np.full(self.mol_xyzr.shape[0], cv_label)
        das_xyzrm = np.vstack([self.mol_xyzr.T, das_marker])
        das_xyzrm[3] += self.cv_thickness
        return das_xyzrm

    def calc_acv(self, use_LabelLib):
        """
        Partition the accessible volume into a free and a contact volume [#]_ [#]_

        Returns
        -------
        acv : fretraj.cloud.ACV
              the attributes from the LabelLib.Grid3D object are:
                - discStep : float  (the grid spacing)
                - originXYZ : list  (x-/y-/z-coordinate of the grid origin)
                - shape : list      (number of grid points in x-/y-/z-direction)
                - grid : list       (flattened list of grid point values)

              additional attributes are:
                - grid_3d : numpy.ndarray([shape])
                - n_gridpts : int
                - tag : numpy.ndarray([1,n_gridpts])

        See Also
        --------
        calc_av : Calculate accessible volume

        References
        ----------
        .. [#] Steffen, F. D., Sigel, R. K. O. & Börner, R. "An atomistic view on carbocyanine \
        photophysics in the realm of RNA", *Phys. Chem. Chem. Phys.* **18**, 29045–29055 (2016).
        .. [#] Dimura, M. et al. Quantitative FRET studies and integrative modeling unravel the structure \
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
        """
        Calculate mean dye position in accessible contact volume

        Parameters
        ----------
        cloud_xyzq : ndarray
                     array of x-,y-,z-coordinates and corresponding weights
                     with a shape [n_gridpts(+), 4]

        Returns
        -------
        mp : ndarray
             mean position
        """
        x = np.dot(cloud_xyzqt[:, 0], cloud_xyzqt[:, 3])
        y = np.dot(cloud_xyzqt[:, 1], cloud_xyzqt[:, 3])
        z = np.dot(cloud_xyzqt[:, 2], cloud_xyzqt[:, 3])
        mp = np.array((x, y, z)) / cloud_xyzqt[:, 3].sum()
        return mp

    @staticmethod
    def median_pos(cloud_xyzqt):
        """
        Calculate median dye position in accessible contact volume

        Parameters
        ----------
        cloud_xyzq : ndarray
                     array of x-,y-,z-coordinates and corresponding weights
                     with a shape [n_gridpts(+), 4]

        Returns
        -------
        mdp : ndarray
              median position
        """
        quantile = 0.5
        mdp = np.apply_along_axis(Volume._weighted_quantile_1D, 0, cloud_xyzqt[:, 0:3], cloud_xyzqt[:, 3], quantile)
        return mdp

    @staticmethod
    def _weighted_quantile_1D(arr_1D, weights, quantile):
        """
        Compute the weighted quantile of a 1D-array

        Parameters
        ----------
        data : ndarray
        weights : ndarray
        quantile : float

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
        """
        Write mean dye position to file

        Parameters
        ----------
        filename : str
        format : ndarray
                 mean position
        units : str
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
