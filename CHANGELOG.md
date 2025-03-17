## Changelog

### [0.2.11] 2025-03-16
- Switch to setuptools build
- Remove legacy Github workflows
- Add detection efficiency to burst simulation
- Use PyMOL's atom `id` / mdtraj's `serial` for `attach_id` and improve documentation

### [0.2.10] 2023-10-22
- Adjust numba, pandas dependencies for conda

### [0.2.9] 2023-10-22
- Loosen numba dependencies for conda

### [0.2.8] 2023-10-22
- Add Python 3.11 support
- Refactor build script
- wheels for Python 3.10, 3.11 for (Linux, Win, and macOS Intel)

### [0.2.7] 2023-10-16
- release yanked

### [0.2.6] 2021-07-19
- Remove nglview as dependency to make package installation faster
- Fix includes and remove tests from deployed wheel and sdist (test suite is available with the source code on Github)
- Fix pybind11 dependency
- Relocate the entrypoint scripts to console.py 

### [0.2.5] 2021-07-16
- release yanked

### [0.2.4] 2021-07-16
- Fix version dependencies of Numpy, Pandas and Matplotlib for Python 3.6

### [0.2.3] 2021-07-16
- Update versions of dependencies to ensure compatibility with Python 3.6
- Include .so and .pyd files in wheels only and .cpp in sdist only

### [0.2.2] 2021-07-16
- Add pybind11 to build-system and create wheels for linux-64/win-64 and Python versions 3.7, 3.8 and 3.9
- Use tox for testing
- Add docker image for PyMOL+FRETraj
- Add black code style
- Update docs for installation

### [0.2.1] 2021-03-19
- Create matrix deployments for different Python versions

### [0.2.0] 2021-03-19
- Add submodule to simulate photon burst to account for shot-noise and better comparison with single-molecule confocal experiments
- Add notebook documenting the burst simulation
- Add pybind11 and build.py to pyproject.toml to build burst.relaxation C++ extension
- setup.py is omitted when building locally (to avoid breaking the metadata) but generated when building via github actions since conda still depends on it. 
- Updated docstrings and docs

### [0.1.5] 2021-03-09
- Improve continuous integration with Github actions
- Add PyPI and conda deployment
- Add first unit and integration tests of ACV calculations
- NGLview for displaying ACVs in jupyter notebooks
- Updated docs with jupyter-books
- Add interactive notebooks with Binder

### [0.1.1] 2021-03-03
- Initial test release of FRETraj, a Python module to calculate 
accessible-contact volumes (ACV) and predict FRET efficiencies.