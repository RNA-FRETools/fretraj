#!/usr/bin/env python3

NAME = 'D'
ELEMENT = 'D'
NAME_MP = 'Mn'
ELEMENT_MP = 'Mn'
OCCUPANCY_MP = -1

_pdb_format = '{:6}{:5d}{:>5}{:1}{:>3} {:1}{:4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}\n'


def open_dx(grid_3d, xyz_min, d_xyz):
    """Return an OpenDX formatted string

    Parameters
    ----------
    density : ndarray
        3-dimensional array of grid points with a shape given by n_xyz
    xyz_min : list
        origin coordinates of the grid
    d_xyz : list
        grid spacing in x-,y- and z-direction

    Returns
    -------
    str
        OpenDX formatted string

    """
    xmin, ymin, zmin = xyz_min
    nx, ny, nz = grid_3d.shape
    dx, dy, dz = d_xyz

    s = ''
    s += 'object 1 class gridpositions counts {:d} {:d} {:d}\n'.format(nx, ny, nz)
    s += 'origin {:0.3f} {:0.3f} {:0.3f}\n'.format(xmin, ymin, zmin)
    s += 'delta {:0.2f} 0 0\n'.format(dx)
    s += 'delta 0 {:0.2f} 0\n'.format(dy)
    s += 'delta 0 0 {:0.2f}\n'.format(dz)
    s += 'object 2 class gridconnections counts {:d} {:d} {:d}\n'.format(nx, ny, nz)
    s += 'object 3 class array type double rank 0 items {:d} data follows\n'.format(nx * ny * nz)
    k = 0
    for ix in range(0, nx):
        for iy in range(0, ny):
            for iz in range(0, nz):
                s += '{}'.format(grid_3d[ix, iy, iz])
                k += 1
                if k % 3 == 0:
                    s += '\n'
                else:
                    s += ' '
    s += 'attribute "dep" string "positions"\n'
    s += 'object "density" class field\n'
    s += 'component "positions" value 1\n'
    s += 'component "connections" value 2\n'
    s += 'component "data" value 3'

    return s


def xyz(cloud_xyzqt, mp, write_weights=True, encode_element=False):
    """Return an XYZ formatted string

    Parameters
    ----------
    cloud_xyzq : ndarray
        array of x-,y-,z-coordinates and corresponding weights with a shape [n_gridpts(+), 4]
    mp : ndarray
         mean position
    write_weights : bool, optional=True
        include weights in XYZ file (5th column)
    encode_element : bool, optional=False

    Returns
    -------
    str
        XYZ formatted string
    """
    n_points = cloud_xyzqt.shape[0]
    s = ''
    s += '{:d}\n'.format(n_points)
    s += '# Accessible contact volume\n'
    if write_weights:
        for k in range(n_points):
            if encode_element:
                if cloud_xyzqt[k, 4] == 2:
                    element = 'C'  # contact volume
                else:
                    element = 'F'  # free volume
            else:
                element = ELEMENT
            s += f'{element:1}\t{cloud_xyzqt[k, 0]:.3f}\t{cloud_xyzqt[k, 1]:.3f}\t{cloud_xyzqt[k, 2]:.3f}\t{cloud_xyzqt[k, 3]:.3f}\n'
        s += f'{ELEMENT_MP:1}\t{mp[0]:.3f}\t{mp[1]:.3f}\t{mp[2]:.3f}\t{-1:.3f}\n'

    else:
        for k in range(n_points):
            if encode_element:
                if cloud_xyzqt[k, 4] == 2:
                    element = 'C'
                else:
                    element = 'F'
            else:
                element = ELEMENT
            s += f'{element:1}\t{cloud_xyzqt[k, 0]:.3f}\t{cloud_xyzqt[k, 1]:.3f}\t{cloud_xyzqt[k, 2]:.3f}\n'
        s += f'{ELEMENT_MP:1}\t{mp[0]:.3f}\t{mp[1]:.3f}\t{mp[2]:.3f}\n'
    return s


def pdb(cloud_xyzqt, mp):
    """Returns a PDB formatted string

    Parameters
    ----------
    cloud_xyzq : ndarray
        array of x-,y-,z-coordinates and corresponding weights with a shape [n_gridpts(+), 4]
    tag : ndarray
        one-dimensional array of length n_gridpts

    Returns
    -------
    str
        PDB formatted string

    Notes
    -----
    the PDB properties are as follows:

    +--------------------------+------------+-----------+
    | entry                    | PyMOL      | pos.      |
    +==========================+============+===========+
    | ATOM / HETATM            |            | [1-6]     |
    +--------------------------+------------+-----------+
    | atom serial number       | (id)       | [7-11]    |
    +--------------------------+------------+-----------+
    | atom name                | (name)     | [12-16]   |
    +--------------------------+------------+-----------+
    | alternate loc. indicator | (alt)      | [17]      |
    +--------------------------+------------+-----------+
    | residue name             | (resn)     | [18-20]   |
    +--------------------------+------------+-----------+
    |                          |            | [21]      |
    +--------------------------+------------+-----------+
    | chain identifier         | (chain)    | [22]      |
    +--------------------------+------------+-----------+
    | residue seq. number      | (resi)     | [23-26]   |
    +--------------------------+------------+-----------+
    | residue insertion code   |            | [27]      |
    +--------------------------+------------+-----------+
    |                          |            | [28-30]   |
    +--------------------------+------------+-----------+
    | x coordinate (in A)      |            | [31-38]   |
    +--------------------------+------------+-----------+
    | y coordinate (in A)      |            | [39-46]   |
    +--------------------------+------------+-----------+
    | z coordinate (in A)      |            | [47-54]   |
    +--------------------------+------------+-----------+
    | occupancy                | (q)        | [55-60]   |
    +--------------------------+------------+-----------+
    | temperature factor       | (b)        | [61-66]   |
    +--------------------------+------------+-----------+
    |                          |            | [67-66]   |
    +--------------------------+------------+-----------+
    | element                  | (element)  | [77-78]   |
    +--------------------------+------------+-----------+
    | charge                   | (fc)       | [79-80]   |
    +--------------------------+------------+-----------+

    """
    n_points = cloud_xyzqt.shape[0]
    bfactor = 99
    s = ''
    for k in range(n_points):
        if cloud_xyzqt[k, 4] == 2:
            resn = 'CV'
        else:
            resn = 'FV'
        s += _pdb_format.format('ATOM', k+1, NAME, ' ', resn, ' ', int(cloud_xyzqt[k, 4]), ' ', cloud_xyzqt[k, 0],
                                cloud_xyzqt[k, 1], cloud_xyzqt[k, 2], bfactor, cloud_xyzqt[k, 3], ELEMENT, ' ')
    s += _pdb_format.format('ATOM', k+2, NAME_MP, ' ', 'MP', ' ', 0, ' ', mp[0], mp[1], mp[2], bfactor, -1,
                            ELEMENT_MP, ' ')
    return s
