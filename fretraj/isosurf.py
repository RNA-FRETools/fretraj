#!/usr/bin/env python3

from pymol import cmd


def smooth_map_from_xyz(name, selection, isolevel):
    """
    Creates a map object from a selection with xyz coordinates (e.g. a PDB or XYZ object) 
    and draws a smooth isosurface at the specified isolevel.

    Parameters
    ----------

    name : str
    selection : xyz object
              e.g. a loaded PDB or XYZ file
    level : float
    """
    grid_spacing = 0.9
    gaussRes = 3.0
    name_surf = name + '_isosurf'
    name_map = name + '_map'
    cmd.alter(selection, 'b=100')
    gaussRes_default = cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution', gaussRes)
    cmd.map_new(name_map, 'gaussian', grid_spacing, selection)
    cmd.isosurface(name_surf, name_map, isolevel, selection)
    cmd.set('gaussian_resolution', gaussRes_default)
    cmd.disable(selection)


def draw_map(name, isomap, isolevel):
    """
    Draws an isosurface of an open-dx/ccp4 map at the specified isolevel.

    Parameters
    ----------

    name : str
    acv_map : map object
              e.g. an ccp4 or open-dx map
    level : float
    """
    name_surf = name + '_isosurf'
    cmd.isosurface(name_surf, isomap, isolevel)


cmd.extend("show_acv", smooth_map_from_xyz)
cmd.extend("show_acv", draw_map)
