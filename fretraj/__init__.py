try:
    import importlib.metadata as ilm # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as ilm # Python 3.7

metadata = ilm.metadata('fretraj')

__version__ = metadata['Version']
__author__ = metadata['Author']
__keywords__ = metadata['Keywords']
__license__ = metadata['License']

from . import cloud
from . import export
from . import fret
from . import grid
from . import restraints
from . import gui

# Import modules that are not strictly required to run FRETraj from the command line 
# (without the PyMOL GUI) and are therefore listed as optional installation dependencies pyproject.toml 
# They provide visualization functionalities for Jupyter notebooks (nglview) and PyMOL 
try:
    import pymol.cmd
except ModuleNotFoundError:
    print('Pymol is not installed. Submodules fretraj.isosurf will not be imported.')
else:
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
