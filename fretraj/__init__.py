"""
FRETraj PyMOL Plugin 

Calculate accessible contact volumes and predict FRET efficiencies

(c) Fabio Steffen, University of Zurich, 2020
"""

from .__about__ import __version__, __author__ 

from . import cloud
from . import export
from . import fret
from . import grid
from . import restraints

# Import modules that are not strictly required to run FRETraj from the command line 
# (without the PyMOL GUI) and are therefore not listed as installation dependencies in setup.py
# They provide visualization functionalities for Jupyter notebooks (nglview) and PyMOL 
try:
    import pymol.cmd
except ModuleNotFoundError:
    print('Pymol is not installed. Submodules fretraj.gui and fretraj.isosurf will not be imported.')
else:
    from . import gui
    from . import isosurf

try:
    import nglview
except ModuleNotFoundError:
    print('nglview is not installed. Submodule fretraj.jupyter will not be imported.')
else:
    try:
        import ipywidgets
    except ModuleNotFoundError:
        print('ipywidgets is not installed. Submodule fretraj.jupyter will not be imported.')
    else:
        from . import jupyter

dialog = None


def __init_plugin__(app=None):
    """
    Add FRETraj plugin to the Plugins Menu
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('FRETraj', run_plugin_gui)


def run_plugin_gui():
    """
    Create the GUI Window
    """
    global dialog
    if dialog is None:
        dialog = gui.App(_pymol_running=True) # noqa: F821
    dialog.show()
