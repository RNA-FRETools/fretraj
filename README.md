<img src="https://raw.githubusercontent.com/fdsteffen/fretraj/master/docs/images/fretraj_logo_readme.png">

[![Build Status](https://github.com/fdsteffen/fretraj/workflows/FRETraj%20build/badge.svg)](https://github.com/fdsteffen/fretraj/actions)
[![Docs Status](https://github.com/fdsteffen/fretraj/workflows/FRETraj%20docs/badge.svg)](https://github.com/fdsteffen/fretraj/actions)
[![PyPI](https://img.shields.io/pypi/v/fretraj)](https://pypi.org/project/fretraj/)
[![Anaconda-Server Badge](https://img.shields.io/badge/Install%20with-conda-green.svg)](https://anaconda.org/rna-fretools/fretraj)
[![codecov](https://codecov.io/gh/fdsteffen/fretraj/branch/master/graph/badge.svg?token=A2E70FbycQ)](https://codecov.io/gh/fdsteffen/fretraj)

*FRETraj* is a Python module for calculating **multiple accessible-contact volumes** (multi-ACV) and predicting **FRET efficiencies**. It provides an interface to the [*LabelLib*](https://github.com/Fluorescence-Tools/LabelLib) library to simulate fluorophore dynamics. The package features a user-friendly **PyMOL plugin**<sup>[1](#pymol)</sup> which can be used to explore different labeling positions when designing FRET experiments. In an AV simulation the fluorophore distribution is estimated by a shortest path search (Djikstra algorithm) using a coarse-grained dye probe. *FRETraj* further implements a **Python-only** version of the geometrical clash search used in *LabelLib*. This should facilitate prototyping of new features for the ACV algorithm.

<img src="https://raw.githubusercontent.com/fdsteffen/fretraj/master/docs/images/graphical_abstract.png">
     
A recent addition to the original AV model (Kalinin et al. *Nat. Methods*, 2012) is the so-called **contact volume** (Steffen et. al. *PCCP* 2016). Here, the accessible volume is split into a free volume (FV, transparent) where the dye is freely diffusing and a contact volume (CV, opaque) where the dye stacks to the biomolecular surface. Time-resolved anisotropy experiments suggest that certain fluorophores, among those the commonly used cyanines Cy3 and Cy5, are particularly prone to interact with both proteins and nucleic acids. The contact volume accounts for this effect by reweighting the point-cloud. By choosing different experimental weights for the free and contact component the AV dye model is refined, making *in silico* FRET predictions more reliable.

## Installation
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
