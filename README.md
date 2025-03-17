<img src="https://raw.githubusercontent.com/fdsteffen/fretraj/master/docs/images/fretraj_logo_readme.png">

[![Build and deploy](https://github.com/RNA-FRETools/fretraj/actions/workflows/build.yml/badge.svg)](https://github.com/RNA-FRETools/fretraj/actions/workflows/build.yml)
[![Documentation](https://github.com/RNA-FRETools/fretraj/actions/workflows/docs.yml/badge.svg)](https://github.com/RNA-FRETools/fretraj/actions/workflows/docs.yml)
[![PyPI](https://img.shields.io/pypi/v/fretraj)](https://pypi.org/project/fretraj/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/fretraj.svg)](https://anaconda.org/conda-forge/fretraj)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10898653.svg)](https://doi.org/10.5281/zenodo.10898653)

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

- F.D. Steffen, R.K.O. Sigel, R. Börner, *Bioinformatics* **2021**. [![](https://img.shields.io/badge/DOI-10.1093/bioinformatics/btab615-blue.svg)](https://doi.org/10.1093/bioinformatics/btab615)

### Additional readings
- F.D. Steffen, R.A. Cunha, R.K.O. Sigel, R. Börner, *Nucleic Acids Res.* **2024**, *52*, e59.
- F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055.
- S. Kalinin, T. Peulen, C.A.M. Seidel et al. *Nat. Methods*, **2012**, *9*, 1218-1225.
- T. Eilert, M. Beckers, F. Drechsler, J. Michaelis, *Comput. Phys. Commun.*, **2017**, *219*, 377–389.
- M. Dimura, T. Peulen, C.A.M. Seidel et al. *Curr. Opin. Struct. Biol.* **2016**, *40*, 163-185.
- M. Dimura, T. Peulen, C.A.M Seidel et al. *Nat. Commun.* **2020**, *11*, 5394.

---

<sup><a name="pymol">1</a></sup> PyMOL was developed by WarrenDeLano and is maintained by Schrödinger, LLC.
