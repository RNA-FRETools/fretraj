---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.8.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# *FRETraj* and Jupyter

This notebook demonstrates how to interface *FRETraj* with Jupyter. We will
compute **ACVs** for the same DNA hairpin as in the previous PyMOL example and
predict a **mean FRET efficiency** between donor and acceptor.

```{code-cell} ipython3
import mdtraj as md
import fretraj as ft
import os
example_dir = '../../src/fretraj/examples/'
```


Let's first load a PDB file of a DNA.

```{code-cell} ipython3
struct = md.load(os.path.join(example_dir+'DNA_hairpin.pdb'))
```

Next, we load a [parameter file](../background/parameter_file) specifying the positions of the donor and acceptor labels.

```{code-cell} ipython3
:tags: [hide-output]
labels = ft.cloud.labeling_params(os.path.join(example_dir+'DNA_hairpin_labels.json'), verbose=False)
labels
```

We now calculate an ACV for the donor (`Cy3-20-C5`) and acceptor (`Cy5-44-P1`) position.

```{code-cell} ipython3
acv_D = ft.cloud.Volume(struct, 'Cy3-20-C5', labels)
acv_A = ft.cloud.Volume(struct, 'Cy5-44-P1', labels)
```

The transfer efficiency is calculated by randomly sampling distances between the two ACVs.

```{code-cell} ipython3
fret = ft.cloud.FRET(acv_D, acv_A, 'Cy3-Cy5', labels)
fret.values
```

Visualize the ACV together with the DNA directly in the notebook.

```{code-cell} ipython3
:tags: ['remove-output']
acv_D_traj = ft.cloud.create_acv_traj([acv_D])
acv_A_traj = ft.cloud.create_acv_traj([acv_A])
ft.jupyter.nglview_trajectory_ACV(struct, acv_D_traj['FV'], acv_A_traj['FV'], 
                                          acv_D_traj['CV'], acv_A_traj['CV'])
```

```{hint}
To see this notebook in action launch a **Binder** instance by clicking on the 🚀 at the top of the page.
```

