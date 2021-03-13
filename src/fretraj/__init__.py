try:
    import importlib.metadata as meta  # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as meta  # Python 3.7

PACKAGE = 'fretraj'

try:
    metadata = meta.metadata(PACKAGE)
except meta.PackageNotFoundError:
    print(f'{PACKAGE} is not installed yet.')
            
__version__ = metadata['Version']
__author__ = metadata['Author']
__keywords__ = metadata['Keywords']
__license__ = metadata['License']

__urls__ = {}
for url in metadata.get_all('Project-URL'):
    __urls__.update(dict([url.replace(' ', '').split(',')]))

from . import cloud
from . import export
from . import fret
from . import grid
from . import jupyter
from . import burst

import warnings
import numba

warnings.simplefilter('ignore', category=numba.errors.NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=numba.errors.NumbaPendingDeprecationWarning)

