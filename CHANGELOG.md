## Changelog

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