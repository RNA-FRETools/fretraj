#!/usr/bin/env python3

from pymol import cmd


@cmd.extend
def smooth_map_from_xyz(name, selection, contour_level, grid_spacing, bfactor=100, gaussRes=3, grid_buffer=2):
    """Creates a map object from a selection with xyz coordinates (e.g. a PDB or XYZ object)
    and draws a smooth isosurface at the specified contour level.

    Parameters
    ----------

    name : str
        name of the map
    selection : xyz object
        loaded PDB or XYZ file
    contour_level : float
        contour level (in sigma units)
    grid_spacing : float
        spacing between grid points (in A)
    bfactor : int
        temperature factor; higher numbers generates smoother surfaces
    gaussRes : int
        Gaussian resolution; higher numbers generate smoother surfaces

    Notes
    ----
    If a map for each ACV of a multi-model PDB file should be generated,
    the PDB file must be loaded with the flag discrete=1 (load as discrete objects
    to allow each ACV to have different numbers of grid points.
    (see also https://www.pymolwiki.org/index.php/Discrete_objects)

    """
    name_surf = name + '_isosurf'
    name_map = name + '_map'
    bfactor_str = 'b={:d}'.format(int(bfactor))
    cmd.alter(selection, bfactor_str)
    cmd.alter(selection, bfactor_str)
    gaussRes_default = cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution', gaussRes)
    # Note: choose state=-3 if all ACVs should be combined into a single map
    cmd.map_new(name_map, 'gaussian', grid_spacing, selection, state=0)
    cmd.isosurface(name_surf, name_map, contour_level, selection, buffer=grid_buffer)
    cmd.set('gaussian_resolution', gaussRes_default)
    cmd.disable(selection)


@cmd.extend
def draw_map(name, isomap, contour_level):
    """
    Draws an isosurface of an open-dx/ccp4 map at the specified contour level.

    Parameters
    ----------
    name : str
        name of the map
    acv_map : map object
        an ccp4 or open-dx map
    contour_level : float
        contour level (in sigma units)
    """
    name_surf = name + '_isosurf'
    cmd.isosurface(name_surf, isomap, contour_level)


@cmd.extend
def set_acv_style(donor_name, acceptor_name, donor_site, acceptor_site, labels, volume_type='AV', transparency=False):
    """
    Set a default style for the ACV clouds

    Parameters
    ----------
    donor_name, acceptor_name : str
        names of the donor and acceptor ACV
    donor_site : str
        reference identifier for the donor labeling position
    acceptor_site : str
        reference identifier for the acceptor labeling position
    labels : dict
        dye, linker and setup parameters for the accessible volume calculation
    volume_type : {'AV', 'CV'}
        entire accessible volume or contact volume
    transparency : bool
        make volume transparent
    """
    cmd.set_color('donor_green', [108, 179, 129])
    cmd.set_color('acceptor_red', [194, 84, 73])
    if volume_type == 'CV':
        contour_level = 'contour_level_CV'
        donor_sele = '{} and resn CV'.format(donor_name)
        acceptor_sele = '{} and resn CV'.format(acceptor_name)
        donor_name = donor_name+'_CV'
        acceptor_name = acceptor_name+'_CV'
    else:
        contour_level = 'contour_level_AV'
        donor_sele = donor_name
        acceptor_sele = acceptor_name
    smooth_map_from_xyz(donor_name, donor_sele, labels['Position'][donor_site][contour_level],
                        labels['Position'][donor_site]['grid_spacing'], labels['Position'][donor_site]['state'])
    smooth_map_from_xyz(acceptor_name, acceptor_sele, labels['Position'][acceptor_site][contour_level],
                        labels['Position'][acceptor_site]['grid_spacing'], labels['Position'][acceptor_site]['state'])
    cmd.color('donor_green', donor_name+'_isosurf')
    cmd.color('acceptor_red', acceptor_name+'_isosurf')
    if transparency:
        cmd.set('transparency', 0.4, donor_name+'_isosurf')
        cmd.set('transparency', 0.4, acceptor_name+'_isosurf')
