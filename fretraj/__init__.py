from fretraj.cloud import *
from fretraj.export import *
from fretraj.fret import *
from fretraj.fit import *
from fretraj.isosurf import *
from fretraj.gui import *

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
        dialog = App(_pymol_running=True)
    dialog.show()
