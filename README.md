# <img src="docs/source/_static/fretraj_logo.png" width="200">

FRETraj is a high-level Python API to the **LabelLib** library (https://github.com/Fluorescence-Tools/LabelLib) to simulate fluorophores which are coupled to a biomolecule of interest. The package features a user-friendly **PyMOL plugin**<sup>[1](#pymol)</sup> which can be used to explore different labeling positions while designing new FRET experiments. In an AV simulation the fluorophore distribution is estimated by a shortest path search (Djikstra algorithm) using a coarse-grained dye probe. FRETraj further implements a **Python-only** version of the geometrical clash search used in LabelLib. This should facilitate prototyping of new features for the ACV algorithm.  

<img src="docs/source/_static/graphical_abstract.png" height="350">
     
A recent addition to the original AV model (Kalinin et al. *Nat. Methods*, 2012) is the so-called **contact volume** (Steffen et. al. *PCCP* 2016). Here, the accessible volume is split into a free volume (FV, transparent) where the dye is freely diffusing and a contact volume (CV, opaque) where the dye stacks to the biomolecular surface. Time-resolved anisotropy experiments suggest that certain fluorophores, among those the commonly used cyanine fluorophores Cy3 and Cy5, are particularly prone to interact with both proteins and nucleic acids. The contact volume accounts for this effect by reweighting the point-cloud. By choosing different experimental weights for the free and contact component the AV dye model is refined, making *in silico* FRET predictions more reliable.


### Plugin installation

You can get the latest version of PyMOL from [Schrödinger](https://pymol.org/). Start the **Anaconda prompt** which comes bundled with PyMOL 2.x and install the necessary dependencies:
```
conda install numpy "numba<=0.44" mdtraj packaging -c conda-forge
```
For a faster calculation of the AVs you may additionally install LabelLib, but this is not strictly required as FRETraj also runs a Python-only implementation of the AV algorithm.
```
conda install -c tpeulen labellib
```

To use the **FRETraj PyMOL plugin** simply download the .zip archive from Github and install it via PyMOL's Plugin manager: `Plugin` &rarr; `Plugin manager` &rarr; `Install New Plugin` &rarr; `Choose file...` and select the .zip archive. Upon first startup FRETraj will prompt you to select a root directory where to store the calculated ACVs and parameter files. You can load a **demo project** by going to ``Help`` :raw-html:`&rarr;` ``Load Example``. You may also want to have a look at this [step-by-step tutorial](https://fdsteffen.github.io/fretraj/pymol_plugin).

<img src="docs/source/_static/PyMOL_interface.PNG" height="250"> <img src="docs/source/_static/PyMOL_Plugin.PNG" height="250">


## References
If you use **FRETraj** in your work please refer to the following paper:
- F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. [![](https://img.shields.io/badge/DOI-10.1039/C6CP04277E-blue.svg)](https://doi.org/10.1039/C6CP04277E)

Here, we introduce the contact volume (CV) as an extension of the accessible volume (AV).

### Related projects
The accessible volume is described in:
- S. Kalinin, T. Peulen, C.A.M. Seidel et al. *Nat. Methods*, **2012**, *9*, 1218-1225. [![](https://img.shields.io/badge/DOI-10.1038/nmeth.2222-blue.svg)](https://doi.org/10.1038/nmeth.2222)

Fast-NPS (nano-positioning system) uses a Bayesian model to locate a dye with respect to a biomolecule: 
- T. Eilert, M. Beckers, F. Drechsler, J. Michaelis, *Comput. Phys. Commun.*, **2017**, *219*, 377–389. [![](https://img.shields.io/badge/DOI-10.1016/j.cpc.2017.05.027-blue.svg)](https://doi.org/10.1016/j.cpc.2017.05.027)

Various dye models have been reviewed in:
- M. Dimura, T. Peulen, C.A.M. Seidel et al. *Curr. Opin. Struct. Biol.* **2016**, *40*, 163-185. [![](https://img.shields.io/badge/DOI-10.1016/j.sbi.2016.11.012-blue.svg)](https://doi.org/10.1016/j.sbi.2016.11.012)
- T. Peulen, O. Opanasyuk, C.A.M Seidel, *J. Phys. Chem . B.*, **2017**, *121*, 8211-8241.[![](https://img.shields.io/badge/DOI-10.1021/acs.jpcb.7b03441-blue.svg)](https://doi.org/10.1021/acs.jpcb.7b03441)

<sup><a name="pymol">1</a></sup> PyMOL is a trademark of Schrödinger, LLC.