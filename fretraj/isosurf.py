#!/usr/bin/env python3

from pymol import cmd


def smooth_map_from_xyz(name, selection, contour_level, grid_spacing, bfactor=100, gaussRes=3, grid_buffer=2):
    """
    Creates a map object from a selection with xyz coordinates (e.g. a PDB or XYZ object) 
    and draws a smooth isosurface at the specified contour level.

    Parameters
    ----------

    name : str
    selection : xyz object
              e.g. a loaded PDB or XYZ file
    contour_level : float
    grid_spacing : float
    bfactor : int
              temperature factor; higher numbers generates smoother surfaces 
              (increasing the b-factor is more computationally efficient than increasing the gaussian resolution)
    gaussRes : int
               Gaussian resolution; higher numbers generate smoother surfaces
    """
    name_surf = name + '_isosurf'
    name_map = name + '_map'
    bfactor_str = 'b={:d}'.format(int(bfactor))
    cmd.alter(selection, bfactor_str)
    cmd.alter(selection, bfactor_str)
    gaussRes_default = cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution', gaussRes)
    cmd.map_new(name_map, 'gaussian', grid_spacing, selection)
    cmd.isosurface(name_surf, name_map, contour_level, selection, buffer=grid_buffer)
    cmd.set('gaussian_resolution', gaussRes_default)
    cmd.disable(selection)


def draw_map(name, isomap, contour_level):
    """
    Draws an isosurface of an open-dx/ccp4 map at the specified contour level.

    Parameters
    ----------

    name : str
    acv_map : map object
              e.g. an ccp4 or open-dx map
    contour_level : float
    """
    name_surf = name + '_isosurf'
    cmd.isosurface(name_surf, isomap, contour_level)


cmd.extend("smooth_map_from_xyz", smooth_map_from_xyz)
cmd.extend("draw_map", draw_map)
