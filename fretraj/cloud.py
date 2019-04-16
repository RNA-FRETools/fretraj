#!/usr/bin/env python3

# Compute and analyze accessible volume clouds on an MD trajectory or a individual PDB structures

# F.Steffen, University of Zurich

import numpy as np
import os
import sys
import argparse
import json
import mdtraj as md
import LabelLib as ll
import numba as nb
import warnings
import copy

from fretraj import export


DISTANCE_SAMPLES = 200000
VERSION = 1.0

package_directory = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(package_directory, 'periodic_table.json')) as f:
    _periodic_table = json.load(f)

VDW_RADIUS = dict((key, _periodic_table[key]["vdw_radius"]) for key in _periodic_table.keys())

_default_struct = {'Position': {'pd_key': {'attach_id': ((int, float), None), 'linker_length': ((int, float), None),
                                           'linker_width': ((int, float), None), 'dye_radius1': ((int, float), None),
                                           'dye_radius2': ((int, float), None), 'dye_radius3': ((int, float), None),
                                           'cv_thickness': ((int, float), 0), 'grid_spacing': ((int, float), 0.5),
                                           'mol_selection': (str, 'all'), 'simulation_type': (str, 'AV3'),
                                           'cv_fraction': (float, 0)}},
                   'Distance': {'pd_key': {'R0': ((int, float), None)}}}

_default_params = {key: val[1] for key, val in _default_struct['Position']['pd_key'].items() if val[1] is not None}


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
                        version='%(prog)s ' + str(VERSION))
    parser.add_argument(
        '-i', '--input', help='Input PDB structure (.pdb)', required=True)
    parser.add_argument('-p', '--parameters',
                        help='Parameter file (.json)', required=True)
    args = parser.parse_args()
    pdb_file = args.input
    param_file = args.parameters
    return pdb_file, param_file


def labeling_params(param_file):
    """
    Parse the parameter file to get the configuration settings for acvCloud

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
        labels = check_labels(labels_json)
    except KeyError as e:
        key = e.args[0]
        pos = e.args[1]
        field = e.args[2]
        print('Missing Key: \'{}\' in {} {}. Exiting...'.format(key, field, pos))
    except TypeError as e:
        print('Wrong data type: {}'.format(e))
    else:
        return labels


def check_labels(labels):
    """
    Check integrity of parameter dictionary

    Parameters
    ----------
    labels : dict

    Returns
    -------
    labels_checked : dict
    """
    for field in _default_struct.keys():
        if field in labels.keys():
            for pos in labels[field].keys():
                if field == 'Position':
                    if 'simulation_type' in labels[field][pos]:
                        if labels[field][pos]['simulation_type'] == 'AV1':
                            labels[field][pos]['dye_radius2'] = 0
                            labels[field][pos]['dye_radius3'] = 0
                    else:
                        raise KeyError('simulation_type', pos, field)
                elif field == 'Distance':
                    pass
                for key, (t, d) in _default_struct[field]['pd_key'].items():
                    if key not in labels[field][pos]:
                        if key in _default_params.keys():
                            labels[field][pos][key] = _default_params[key]
                            print('Missing Key: \'{}\' in {} {}. Falling back to {}'.format(key, field, pos, _default_params[key]))
                        else:
                            raise KeyError(key, pos, field)
                    else:
                        if not isinstance(labels[field][pos][key], t):
                            raise TypeError('\'{}\' in {} {} must be of one of the following types: {}'.format(key, field, pos, t))
        else:
            labels[field] = None
    return labels


class ACV:
    """
    Reweighted accessible contact volume of a covalently linked fluorophore

    Parameters
    ----------
    ll_acv : LabelLib.Grid3D
             accessible volume with added weight labels for free and contact volume 
             (density label FV: 1.0, density label CV: 2.0); the two volumes are subsequently reweighted
    cv_thickness : float
                   width of the contact volume in Angstrom (default: 0, i.e. no CV is calculated)
                   **Note:** good first approx.: 2*min(dye_radii)
    cv_fraction : float
                  fraction of dyes that are within the contact volume
                  (e.g. as determined by fluorescence anisotropy)

    Attributes
    ----------
    grid : list
           flattened list of grid point values of the LabelLib.Grid3D
    shape : list
            shape of the LabelLib.Grid3D
    originXYZ : list
                origin of the LabelLib.Grid3D
    discStep : float
               spacing of the LabelLib.Grid3D

    n_gridpts : int
                total number of grid points

    grid_1d : ndarray
              one-dimensional array of grid points of length n_gridpts
    tag : ndarray
          one-dimensional array of length n_gridpts
    grid_3d : ndarray
              3-dimensional array of grid points with a shape given by n_xyz
    cloud_xyzq : ndarray
                 array of x-,y-,z-coordinates and corresponding weights
                 with a shape [n_gridpts(+), 4]
    """

    def __init__(self, ll_acv, cv_thickness, cv_fraction):
        self.grid = ll_acv.grid
        self.shape = ll_acv.shape
        self.originXYZ = ll_acv.originXYZ
        self.discStep = ll_acv.discStep

        self.n_gridpts = np.prod(self.shape)
        self.grid_1d = self._reweight_cv(cv_thickness, cv_fraction)
        self.tag = self._tag_volume()
        self.grid_3d = ACVolume.reshape_grid(self.grid_1d, self.shape)
        self.cloud_xyzq = ACVolume.grid2pts(self.grid_3d, self.originXYZ, [self.discStep] * 3)

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
        if cv_thickness != 0:
            weight_cv = self._weight_factor(grid_1d, cv_fraction)
            grid_1d[grid_1d > 1.0] = weight_cv
        grid_1d = np.clip(grid_1d, 0, None)
        return grid_1d

    @staticmethod
    def _weight_factor(grid_1d, cv_fraction):
        """
        Calculate the weight of contact volume grid points.

        Extended summary
        ----------------
        The weight factor corresponds to the factor by which the points of the contact volume (CV, dye trapped on surface)
        are favored over those belonging the free volume (FV, free dye diffusion). This accounts for the (experimentally)
        determined fraction of dyes populating the CV.

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
        n_FV = np.count_nonzero(grid_1d == 1.0)
        weight_cv = n_FV / n_CV * cv_fraction / (1.0 - cv_fraction)
        return weight_cv

    def _tag_volume(self):
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
        tag = np.full(self.n_gridpts, 1)
        tag[self.grid_1d > 1.0] = 2
        return tag


class ACVolume:
    """
    Class object holding the accessible contact volume of a specific labeling position

    Parameters
    ----------
    structure : mdtraj.Trajectory
                trajectory of atom coordinates loaded from a pdb, xtc or other file
    frame : int
            frame number of the mdtraj.Trajectory object
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
    frame : int
            frame number of the mdtraj.Trajectory object
    attach_id : int
                serial atom id of the attachment point in the mdtraj.Trajectory object
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
    """

    def __init__(self, structure, frame, site, labels):
        self.labeling_site = site
        self.structure = structure
        self.frame = frame
        self.attach_id = labels['Position'][site]['attach_id']
        self.mol_selection = labels['Position'][site]['mol_selection']
        self.linker_length = labels['Position'][site]['linker_length']
        self.linker_width = labels['Position'][site]['linker_width']
        self.simulation_type = labels['Position'][site]['simulation_type']
        self.dye_radii = np.array([labels['Position'][site]['dye_radius1'], labels['Position']
                                   [site]['dye_radius2'], labels['Position'][site]['dye_radius3']])
        self.grid_spacing = labels['Position'][site]['grid_spacing']
        self.cv_thickness = labels['Position'][site]["cv_thickness"]
        self.cv_fraction = labels['Position'][site]["cv_fraction"]
        self.av = self.calc_av()
        self.acv = self.calc_acv()

    @classmethod
    def from_attachID(cls, structure, frame, labels, attachID):
        """
        Alternative constructor for ACVolume

        Parameters
        ----------
        structure : mdtraj.Trajectory
        frame : int
        site : str
        labels : dict

        Returns
        -------
        ACVolume : class instance
                   Returns an instance of the ft.cloud.ACVolume class
        """
        site = 'global'
        try:
            global_params = labels['Position'][site]
        except KeyError as e:
            print('Site key missing: {}. Cannot initalize ACVolume.'.format(e))
        else:
            labels['Position'].clear()
            labels['Position'][attachID] = global_params
            labels['Position'][attachID]['attach_id'] = attachID
            return cls(structure, frame, attachID, labels)

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
        >>> ft.ACVolume.reshape_grid(grid, (3,3,3))
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
    @nb.jit
    def grid2pts(grid_3d, xyz_min, d_xyz):
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
        cloud_xyzq : ndarray
                     array of x-,y-,z-coordinates and corresponding weights
                     with a shape [n_gridpts(+), 4]

        Examples
        --------

        >>> grid_3d = numpy.zeros((3,3,3))
        >>> grid_3d[(1,1,1)] = 1
        >>> grid_3d[(1,2,1)] = 10
        >>> ft.ACVolume.grid2pts(grid_3d, xyz_min=[0,0,0], d_xyz=[1,1,1])
        array([[ 0.5,  0.5,  0.5,  1. ],
               [ 0.5,  1. ,  0.5, 10. ]])
        """
        xmin, ymin, zmin = xyz_min
        nx, ny, nz = grid_3d.shape
        dx, dy, dz = d_xyz

        cloud_xyzq = np.empty((nx * ny * nz, 4), dtype=np.float64, order='C')

        gdx = np.arange(0, nx, dtype=np.float64) * dx
        gdy = np.arange(0, ny, dtype=np.float64) * dy
        gdz = np.arange(0, nz, dtype=np.float64) * dz

        n = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    d = grid_3d[ix, iy, iz]
                    if d > 0:
                        cloud_xyzq[n, 0] = gdx[ix] + xmin
                        cloud_xyzq[n, 1] = gdy[iy] + ymin
                        cloud_xyzq[n, 2] = gdz[iz] + zmin
                        cloud_xyzq[n, 3] = d
                        n += 1
        return cloud_xyzq[:n]

    @property
    def mol_xyzr(self):
        """
        Get coordinates and vdW radii of all selected atoms in the structure

        Returns
        -------
        xyzr : ndarray
               array of x-,y-,z-coordinates and VdW radii with a shape [n_atoms, 4]

        Examples
        --------

        >>> avobj.mol_xyzr

        """
        frame = self.frame
        struct = self.structure
        # sele = struct.top.select(self.mol_selection)
        # cutoff = 3
        # sele = md.compute_neighbors(struct, cutoff, [self.attach_id])
        # xyz = struct.xyz[frame][sele[frame]]
        xyz = struct.xyz[frame] * 10

        radii = np.array([VDW_RADIUS[atom.element.symbol] /
                          100 for atom in struct.top.atoms], ndmin=2).T
        radii[self.attach_id] = 0.0
        xyzr = np.hstack((xyz, radii))

        # x0 = np.hstack([self.attach_xyz, 0.0])
        # asr = 4.5
        # sd = (xyzr - x0)**2
        # d2 = sd[:, 0] + sd[:, 1] + sd[:, 2]
        # xyzr[np.where(d2 < 0 ** 2)[0]] *= np.array([1.0, 1.0, 1.0, 0.0])
        return xyzr

    @property
    def attach_xyz(self):
        """
        Get coordinates of the dye attachment point

        Returns
        -------
        xyz : ndarray
              one-dimensional array of x-,y-,z-coordinates of length 3

        Examples
        --------

        >>> avobj.attach_xyz

        """
        xyz = self.structure.xyz[self.frame][self.attach_id] * 10
        return xyz

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

        Examples
        --------

        >>> obj.save_acv('cloud.xyz', format='xyz', write_weights=False)

        """
        if format == 'xyz':
            try:
                write_weights = kwargs['write_weights']
            except KeyError:
                write_weights = True
            file_str = export.xyz(self.acv.cloud_xyzq, write_weights)
        elif format == 'open_dx':
            d_xyz = [self.grid_spacing] * 3
            xyz_min = self.acv.originXYZ
            file_str = export.open_dx(self.acv.grid_3d, xyz_min, d_xyz)
        else:
            file_str = export.pdb(self.acv.cloud_xyzq, self.acv.tag)
        with open(filename, 'w') as fname:
            fname.write(file_str)

    def calc_av(self):
        """
        Calculate the dye accessible volume [1]_ [2]_

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
        .. [1] Kalinin, S. et al. "A toolkit and benchmark study for FRET-restrained high-precision \
        structural modeling", *Nat. Methods* **9**, 1218–1225 (2012).
        .. [2] Sindbert, S. et al. "Accurate distance determination of nucleic acids via Förster \
        resonance energy transfer: implications of dye linker length and rigidity", \
        *J. Am. Chem. Soc.* **133**, 2463–2480 (2011).

        Examples
        --------

        >>> avobj.calc_av()

        """
        mol_xyzr = self.mol_xyzr
        attach_xyz = self.attach_xyz
        if self.simulation_type == 'AV1':
            av = ll.dyeDensityAV1(mol_xyzr.T, attach_xyz, self.linker_length,
                                  self.linker_width, self.dye_radii[0], self.grid_spacing)
        else:
            av = ll.dyeDensityAV3(mol_xyzr.T, attach_xyz, self.linker_length,
                                  self.linker_width, self.dye_radii, self.grid_spacing)
        return av

    def _dye_acc_surf(self):
        """
        Calculate dye accessible surface by padding the vdW radius with the thickness of the contact volume

        Returns
        -------
        das_xyz : numpy.ndarray([5,n_atoms])
                   array of marked coordinates and padded vdW radii (n_atoms = number of atoms in mdtraj.Trajectory)
        """
        cv_label = 2.0
        das_marker = np.full(self.mol_xyzr.shape[0], cv_label)
        das_xyzm = np.vstack([self.mol_xyzr.T, das_marker])
        das_xyzm[3] += self.cv_thickness
        return das_xyzm

    def calc_acv(self):
        """
        Partition the accessible volume into a free and a contact volume [3]_ [4]_

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
        .. [3] Steffen, F. D., Sigel, R. K. O. & Börner, R. "An atomistic view on carbocyanine \
        photophysics in the realm of RNA", *Phys. Chem. Chem. Phys.* **18**, 29045–29055 (2016). 
        .. [4] Dimura, M. et al. Quantitative FRET studies and integrative modeling unravel the structure \
        and dynamics of biomolecular systems", *Curr. Opin. Struct. Biol.* **40**, 163–185 (2016).

        Examples
        --------

        >>> avobj.calc_acv()

        """
        if self.cv_thickness != 0:
            das_xyzm = self._dye_acc_surf()
            ll_acv = ll.addWeights(self.av, das_xyzm)
        else:
            ll_acv = copy.deepcopy(self.av)

        acv = ACV(ll_acv, self.cv_thickness, self.cv_fraction)
        return acv

        # print(acv)
        # acv.n_gridpts = np.prod(acv.shape)
        # acv.tag = self._tag_volume(acv.n_gridpts, grid_1d)
        # acv.grid = list(np.clip(grid_1d, 0, None))
        # acv.grid_3d = self._reshape_grid(acv)


if __name__ == "__main__":
    pdb_file, param_file = parseCmd()
    labels = labeling_params(param_file)
    struct = md.load_pdb(pdb_file)
    av1 = ACVolume(struct, 0, '344', labels)
    points = density2points(av1.av.discStep, density_av1, av1.av.originXYZ)
    writeXYZ('PDB', points, av1.attach_xyz)

    sys.exit()
    cloud = runACV(par, biomol)
    writeXYZ('PDB', cloud.bio.coords, cloud.coords_AttachPos)
    writeXYZ('AV', cloud.av.coords, cloud.coords_AttachPos, cloud.av.weights)
    writeXYZ('CV', cloud.cv.coords, cloud.coords_AttachPos)
    writeXYZ('FV', cloud.fv.coords, cloud.coords_AttachPos)
