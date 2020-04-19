"""
FRETraj PyMOL Plugin 

Calculate accessible contact volumes and predict FRET efficiencies

(c) Fabio Steffen, University of Zurich, 2020
"""

from .__about__ import __version__, __author__ 

from . import cloud
from . import export
from . import fret
from . import isosurf
from . import gui
from . import grid
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
        dialog = gui.App(_pymol_running=True)
    dialog.show()
