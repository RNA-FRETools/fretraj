<img src="https://raw.githubusercontent.com/fdsteffen/fretraj/master/docs/images/fretraj_logo_readme.png">

[![Build Status](https://github.com/fdsteffen/fretraj/workflows/FRETraj%20test/badge.svg)](https://github.com/fdsteffen/fretraj/actions)
[![Docs Status](https://github.com/fdsteffen/fretraj/workflows/FRETraj%20docs/badge.svg)](https://github.com/fdsteffen/fretraj/actions)
[![PyPI](https://img.shields.io/pypi/v/fretraj)](https://pypi.org/project/fretraj/)
[![codecov](https://codecov.io/gh/fdsteffen/fretraj/branch/master/graph/badge.svg?token=A2E70FbycQ)](https://codecov.io/gh/fdsteffen/fretraj)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

*FRETraj* is a Python module for **predicting FRET efficiencies** by calculating multiple accessible-contact volumes (multi-ACV) to estimate donor and acceptor dye dynamics. The package features a user-friendly **PyMOL plugin**<sup>[1](#pymol)</sup> for for FRET-assisted, integrative structural modeling. It interfaces with the [*LabelLib*](https://github.com/Fluorescence-Tools/LabelLib) library for fast computation of ACVs. 
Specifically, *FRETraj* is designed to:
- plan (single-molecule) FRET experiments by optimizing **labeling positions**
- interpret FRET-based **distance measurements** on a biomolecule
- integrate FRET experiments with **molecular dynamics** simulations

<img src="https://raw.githubusercontent.com/fdsteffen/fretraj/master/docs/images/graphical_abstract.png">

## Installation and Documentation
Follow the instructions for your platform [here](https://rna-fretools.github.io/fretraj/getting_started/installation)

## References
If you use **FRETraj** in your work please refer to the following paper:

- F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. [![](https://img.shields.io/badge/DOI-10.1039/C6CP04277E-blue.svg)](https://doi.org/10.1039/C6CP04277E)

### Additional readings
- S. Kalinin, T. Peulen, C.A.M. Seidel et al. *Nat. Methods*, **2012**, *9*, 1218-1225.
- T. Eilert, M. Beckers, F. Drechsler, J. Michaelis, *Comput. Phys. Commun.*, **2017**, *219*, 377–389.
- M. Dimura, T. Peulen, C.A.M. Seidel et al. *Curr. Opin. Struct. Biol.* **2016**, *40*, 163-185.
- M. Dimura, T. Peulen, C.A.M Seidel et al. *Nat. Commun.* **2020**, *11*, 5394.

---

<sup><a name="pymol">1</a></sup> PyMOL is a trademark of Schrödinger, LLC.
