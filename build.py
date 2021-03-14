import os
import shutil


from distutils.command.build_ext import build_ext
from distutils.core import Distribution
from distutils.core import Extension
from distutils.errors import CCompilerError
from distutils.errors import DistutilsExecError
from distutils.errors import DistutilsPlatformError


class get_pybind_include(object):
    def __str__(self):
        import pybind11
        return pybind11.get_include()


extensions = [
    Extension('relaxation', sources=['src/fretraj/relaxation.cpp'], include_dirs=[get_pybind_include()], language='c++')
    ]


class ExtBuilder(build_ext):

    built_extensions = []

    def run(self):
        try:
            build_ext.run(self)
        except (DistutilsPlatformError, FileNotFoundError):
            print('Unable to build C++ extensions of FRETraj')

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError, ValueError):
            print(f'Unable to build the {ext.name} C++ extension, the burst module will not be available')


def build(setup_kwargs):

    distribution = Distribution({"name": "FRETraj", "ext_modules": extensions, "package_dir": "fretraj"})

    cmd = ExtBuilder(distribution)
    cmd.ensure_finalized()
    cmd.run()

    for output in cmd.get_outputs():
        if not os.path.exists(output):
            continue
        relative_extension = os.path.relpath(output, cmd.build_lib)
        destination = os.path.join('src', distribution.package_dir, relative_extension)
        shutil.copyfile(output, destination)
        mode = os.stat(destination).st_mode
        mode |= (mode & 0o444) >> 2
        os.chmod(destination, mode)

    return setup_kwargs


if __name__ == "__main__":
    build({})
