#!/usr/bin/env python3

# start a PyMOL server session from a terminal:
#   pymol -R

# On Windows you may create a shortcut that executes the following command:
#   C:\path\to\PyMOLWin.exe -R

import os
import re
import nglview
import ipywidgets


def connect2pymol():
    """
    Establish an RPC connection to run PyMOL commands from the command line or from a Jupyter notebook
    """
    import xmlrpc.client as xmlrpclib
    cmd = xmlrpclib.ServerProxy('http://localhost:9123')
    curr_wd = os.getcwd()
    try:
        cmd.cd(curr_wd)
    except:
        cmd.cd(re.sub(r'/mnt/([a-z])', r'\1:', curr_wd))
    return cmd


def nglview_trajectory(traj_biomol):
    view = nglview.NGLWidget()
    view.add_trajectory(traj_biomol)
    view.clear_representations(component=0)
    view.add_simplified_base(component=0, selection='/0', disablePicking=True, colorScheme='atomindex')
    view.add_cartoon(component=0, selection='/0', aspectRatio=4, disablePicking=True, colorScheme='atomindex')
    view.stage.set_parameters(mouse_preset='pymol')
    return view


def nglview_trajectory_AV(traj_biomol, traj_volume1, traj_volume2, surface_representation=False):
    """Create a nglview trajectory scene with donor and acceptor accessible volumes (AV)

    Parameters
    ----------
    traj_biomol : mdtraj.Trajectory
        trajectory of the biomolecule
    traj_volume1, traj_volume2 : mdtraj.Trajectory
        trajectory of the donor and acceptor accessible volume
    surface_representation : bool (default: False)
        show a surface representation (instead of spacefill, which is faster)

    Returns
    -------
    view : nglview.NGLWidget
    """
    view = nglview.NGLWidget()
    view.add_trajectory(traj_biomol)
    view.add_trajectory(traj_volume1)
    view.add_trajectory(traj_volume2)

    view.clear_representations(component=0)
    view.clear_representations(component=1)
    view.clear_representations(component=2)

    view.add_simplified_base(component=0, selection='/0', disablePicking=True, colorScheme='atomindex')
    view.add_cartoon(component=0, selection='/0', aspectRatio=4, disablePicking=True, colorScheme='atomindex')

    if surface_representation:
        view.add_surface(color='#6cb381', wireframe=True, opacity=0.4, isolevel=0, component=1,
                         disablePicking=True, selection='all')
        view.add_surface(color='#c25449', wireframe=True, opacity=0.4, isolevel=0, component=2,
                         disablePicking=True, selection='all')
    else:
        view.add_spacefill(color='#6cb381', component=1, disablePicking=True, selection='all')
        view.add_spacefill(color='#c25449', component=2, disablePicking=True, selection='all')
    view.stage.set_parameters(mouse_preset='pymol')
    return view


def nglview_trajectory_ACV(traj_biomol, traj_volume1_FV, traj_volume2_FV, traj_volume1_CV, traj_volume2_CV):
    """Create a nglview trajectory scene with donor and acceptor accessible-contact volumes (ACV)

    Parameters
    ----------
    traj_biomol : mdtraj.Trajectory
        trajectory of the biomolecule
    traj_volume1_FV, traj_volume2_FV : mdtraj.Trajectory
        trajectory of the donor and acceptor free volume
    traj_volume1_CV, traj_volume2_CV : mdtraj.Trajectory
        trajectory of the donor and acceptor contact volume

    Returns
    -------
    view : nglview.NGLWidget
    """
    view = nglview.NGLWidget()
    view.add_trajectory(traj_biomol)
    view.add_trajectory(traj_volume1_CV)
    view.add_trajectory(traj_volume2_CV)
    view.add_trajectory(traj_volume1_FV)
    view.add_trajectory(traj_volume2_FV)

    view.clear_representations(component=0)
    view.clear_representations(component=1)
    view.clear_representations(component=2)
    view.clear_representations(component=3)
    view.clear_representations(component=4)

    view.add_simplified_base(component=0, selection='/0', disablePicking=True, colorScheme='atomindex')
    view.add_cartoon(component=0, selection='/0', aspectRatio=4, disablePicking=True, colorScheme='atomindex')

    view.add_surface(color='#6cb381', wireframe=False, opacity=0.4, isolevel=0, component=1,
                     disablePicking=True, selection='all')
    view.add_surface(color='#c25449', wireframe=False, opacity=0.4, isolevel=0, component=2,
                     disablePicking=True, selection='all')
    view.add_surface(color='#6cb381', wireframe=True, opacity=0.4, isolevel=0, component=3,
                     disablePicking=True, selection='all')
    view.add_surface(color='#c25449', wireframe=True, opacity=0.4, isolevel=0, component=4,
                     disablePicking=True, selection='all')
    view.stage.set_parameters(mouse_preset='pymol')
    return view


def nglview_multimodel_ACV(biomol_filename, volume1_filename, volume2_filename):
    """Create a nglview multi-model structure scene with donor and acceptor accessible-contact volumes (ACV)

    Parameters
    ----------
    structure_filename : str
                         multi-model PDB of the biomolecule
    volume1_filename, volume2_filename : str
        multi-model PDB of the donor and acceptor accessible-contact volume

    Returns
    -------
    view : nglview.NGLWidget
    n_model : int
    """
    global view
    view = nglview.NGLWidget()
    struct = nglview.FileStructure(biomol_filename)
    acv_D = nglview.FileStructure(volume1_filename)
    acv_A = nglview.FileStructure(volume2_filename)

    struct_str = struct.get_structure_string()
    n_models = struct_str.count('MODEL')

    view.add_component(struct)
    view.add_component(acv_D)
    view.add_component(acv_A)

    view.clear_representations(component=0)
    view.clear_representations(component=1)
    view.clear_representations(component=2)

    view.add_simplified_base(component=0, selection='/0', disablePicking=True, colorScheme='atomindex')
    view.add_cartoon(component=0, selection='/0', aspectRatio=4, disablePicking=True, colorScheme='atomindex')

    view.add_surface(color='#6cb381', wireframe=False, opacity=0.4, isolevel=0, component=1,
                     disablePicking=True, selection='CV')
    view.add_surface(color='#c25449', wireframe=False, opacity=0.4, isolevel=0, component=2,
                     disablePicking=True, selection='CV')
    view.add_surface(color='#6cb381', wireframe=True, opacity=0.4, isolevel=0, component=1,
                     disablePicking=True, selection='all')
    view.add_surface(color='#c25449', wireframe=True, opacity=0.4, isolevel=0, component=2,
                     disablePicking=True, selection='all')
    view.stage.set_parameters(mouse_preset='pymol')
    return view, n_models


def _change_model(model):
    """Change the current model

    Parameters
    ----------
    model : int
    """
    for i in range(view.n_components):
        view._remote_call('setSelection', target='compList', args=[f'/{model}'], kwargs=dict(component_index=i))
        n_representations = 5  # more than 5 representations are unlikely
        for i in range(n_representations):
            view._remote_call('setSelection', target='Representation', args=[f'/{model}'],
                              kwargs=dict(component_index=0, repr_index=i))


def model_slider(n_models):
    """A slider widget for multi-model structures

    Parameters
    ----------
    n_models : int
               number of models in the structure

    Returns
    -------
    int_slider: ipywidgets.IntSlider

    Notes
    -----
    See also: https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20List.html
    """
    int_slider = ipywidgets.IntSlider(
        value=0,
        min=0,
        max=n_models-1,
        step=1,
        description='Frame:',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='d'
        )
    int_slider.observe(_on_slider_value_change, names='value')
    return int_slider


def _on_slider_value_change(change):
    """Updates the view upon changing the slider value

    Parameters
    ----------
    change : dict

    Notes
    -----
    See also: https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20Events.html
    """
    output = ipywidgets.Output()
    with output:
        _change_model(change['new'])


def render_view(view, gui=False):
    """Display an nglview scene with or without a GUI

    Parameters
    ----------
    view : nglview.NGLWidget
    gui : bool
          display the scene within a GUI window
          (useful for interacting with components and representations)
    """
    print('Building scene, please wait...')
    return view.display(gui)
