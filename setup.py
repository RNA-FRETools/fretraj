from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

# https://pybind11.readthedocs.io/en/stable/compiling.html
# https://github.com/python-poetry/poetry/issues/2740

ext_modules = [Pybind11Extension("relaxation", ["src/fretraj/relaxation.cpp"], extra_compile_args=["-std=c++11"])]


setup(ext_modules=ext_modules)
