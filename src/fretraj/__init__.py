import sys

try:
    import importlib.metadata as meta  # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as meta  # Python 3.7

PACKAGE = 'fretraj'

try:
    metadata = meta.metadata(PACKAGE)
except meta.PackageNotFoundError:
    print(f'{PACKAGE} is not installed yet.')


def _get_urls():
    __urls__ = {}
    try:
        for _url in metadata.get_all('Project-URL'):
            __urls__.update(dict([_url.replace(' ', '').split(',')]))
    except TypeError:
        __urls__['Documentation'] = None
        __urls__['Repository'] = None
    return __urls__


__version__ = metadata['Version']
__author__ = metadata['Author']
__keywords__ = metadata['Keywords']
__license__ = metadata['License']
__urls__ = _get_urls()

from . import cloud
from . import export
from . import fret
from . import grid
from . import jupyter

if sys.platform.startswith('linux'):
    try:
        from . import burst
    except ModuleNotFoundError:
        print('The burst module could not be imported\n')
else:
    print('The burst module is currently only available on linux\n')

import warnings
import numba

warnings.simplefilter('ignore', category=numba.errors.NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=numba.errors.NumbaPendingDeprecationWarning)
