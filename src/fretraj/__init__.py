try:
    import importlib.metadata as meta  # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as meta  # Python 3.7

try:
    metadata = meta.metadata('fretraj')
except meta.PackageNotFoundError:
    print('FRETraj is not installed yet.')

__version__ = metadata['Version']
__author__ = metadata['Author']
__keywords__ = metadata['Keywords']
__license__ = metadata['License']

from . import cloud
from . import export
from . import fret
from . import grid
from . import jupyter

import warnings
import numba as nb

warnings.simplefilter('ignore', category=nb.errors.NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=nb.errors.NumbaPendingDeprecationWarning)
