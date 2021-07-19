import sys

try:
    import importlib.metadata as meta  # Python >=3.8
except ModuleNotFoundError:
    import importlib_metadata as meta  # Python 3.7

PACKAGE = "fretraj"

try:
    metadata = meta.metadata(PACKAGE)
except meta.PackageNotFoundError:
    print(f"{PACKAGE} is not installed yet.")


def _get_urls():
    __urls__ = {}
    try:
        for _url in metadata.get_all("Project-URL"):
            __urls__.update(dict([_url.replace(" ", "").split(",")]))
    except TypeError:
        __urls__["Documentation"] = None
        __urls__["Repository"] = None
    return __urls__


__version__ = metadata["Version"]
__author__ = metadata["Author"]
__keywords__ = metadata["Keywords"]
__license__ = metadata["License"]
__urls__ = _get_urls()


try:
    import LabelLib as ll
except ModuleNotFoundError:
    _LabelLib_found = False
else:
    _LabelLib_found = True

try:
    import nglview
except ModuleNotFoundError:
    _nglview_found = False
    _jupyter_module_available = False
else:
    _nglview_found = True
    _jupyter_module_available = True
    from . import jupyter

try:
    from . import burst
except (ImportError, ModuleNotFoundError):
    _burst_module_available = False
    print("The burst submodule could not be imported\n")
else:
    _burst_module_available = True

from . import cloud
from . import export
from . import fret
from . import grid
