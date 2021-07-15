# Welcome to FRETraj
*FRETraj* is a Python module for **predicting FRET efficiencies** by calculating multiple accessible-contact volumes (multi-ACV) to estimate donor and acceptor dye dynamics. The package features a user-friendly **PyMOL plugin** [^PyMOL] for FRET-assisted, integrative structural modeling. It interfaces with the [*LabelLib*](https://github.com/Fluorescence-Tools/LabelLib) library for fast computation of ACVs. Specifically, *FRETraj* is designed to:
- plan (single-molecule) FRET experiments by optimizing **labeling positions**
- interpret FRET-based **distance measurements** on a biomolecule
- integrate FRET experiments with **molecular dynamics** simulations

```{figure} images/graphical_abstract.png
---
width: 100%
name: graphical_abstract
align: left
---
Schematic of *FRETraj*: A fluorophore (here Cy3) is parametrized by five distances ($R_1-R_3$, $L_\text{linker}$ and $W_\text{linker}$). The extent of surface stacking is assessed by fluorescence anisotropy decays. Accessible-contact volume calculations and FRET predictions are implemented in PyMOL (with a GUI) and Jupyter notebooks (for large-scale analyses).
```


[^PyMOL]: PyMOL is a trademark of Schrödinger, LLC.
