#!/usr/bin/env python3


def show_acv2(name, selection, level):
    grid_spacing = 0.9
    gaussRes_acv = 3.0
    acv_surf = name + '_surf'
    acv_map = name + '_map'
    cmd.alter(selection, 'b=100')
    gaussRes_default = cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution', gaussRes_acv)
    cmd.map_new(acv_map, 'gaussian', grid_spacing, selection)
    cmd.isosurface(acv_surf, acv_map, level, selection)
    cmd.set('gaussian_resolution', gaussRes_default)
    cmd.disable(selection)


def show_acv(name, acv_map, level):
    grid_spacing = 0.9
    gaussRes_acv = 3.0
    acv_surf = name + '_surf'
    gaussRes_default = cmd.get('gaussian_resolution')
    cmd.set('gaussian_resolution', gaussRes_acv)
    cmd.isosurface(acv_surf, acv_map, level)
    cmd.set('gaussian_resolution', gaussRes_default)


cmd.extend("show_acv", show_acv)
