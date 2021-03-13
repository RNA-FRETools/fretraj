import setuptools

class get_pybind_include(object):
    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [setuptools.Extension('relaxation', sources=['src/fretraj/relaxation.cpp'], include_dirs=[get_pybind_include()], language='c++')]

def build(setup_kwargs):
    setup_kwargs.update({
        'ext_modules': ext_modules,
    })
