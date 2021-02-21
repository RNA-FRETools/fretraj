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
from . import jupyter
