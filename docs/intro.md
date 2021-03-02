# Welcome to FRETraj

*FRETraj* is a Python module for calculating **multiple accessible-contact volumes** (multi-ACV) and predicting **FRET efficiencies**. It provides an interface to the [*LabelLib*](https://github.com/Fluorescence-Tools/LabelLib) library to simulate fluorophore dynamics. The package features a user-friendly **PyMOL plugin** [^PyMOL] which can be used to explore different labeling positions when designing FRET experiments. In an AV simulation the fluorophore distribution is estimated by a shortest path search (Djikstra algorithm) using a coarse-grained dye probe. *FRETraj* further implements a **Python-only** version of the geometrical clash search used in *LabelLib*. This should facilitate prototyping of new features for the ACV algorithm.

```{figure} images/graphical_abstract.png
---
width: 100%
name: graphical_abstract
align: left
---
Schematic of *FRETraj*: A fluorophore (here Cy3) is parametrized by five distances ($R_1-R_3$, $L_\text{linker}$ and $W_\text{linker}$). The extent of surface stacking is assessed by fluorescence anisotropy decays. Accessible-contact volume calculations and FRET predictions are implemented in PyMOL (with a GUI) and Jupyter notebooks (for large-scale analyses).
```


[^PyMOL]: PyMOL is a trademark of Schrödinger, LLC.
