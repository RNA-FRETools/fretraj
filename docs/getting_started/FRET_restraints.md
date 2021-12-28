---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Introducing FRET restraints

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import mdtraj as md
import fretraj as ft
```

```{code-cell} ipython3
import importlib
```

### Accessible contact volume

Calculate ACVs for donor and acceptor dye on an **solvated and energy minimized** biomolecule.
> Note: the generated .gro file must first be converted to a .pdb file using `gmx trjconv`

```{code-cell} ipython3
labels = ft.cloud.labeling_params('gmx/mn_riboswitch_labels.json', verbose=False)
```

```{code-cell} ipython3
struct = md.load_pdb('gmx/mn_riboswitch.pdb')
struct_nosolv = struct.remove_solvent()
```

```{code-cell} ipython3
don = ft.cloud.Volume(struct_nosolv, '1-B-U1-O5\'', labels)
acc = ft.cloud.Volume(struct_nosolv, '1-A-A1-O5\'', labels)
```

Compute a FRET trajectory from the ACVs

```{code-cell} ipython3
fret_traj = ft.cloud.FRET(don, acc, 'A1-U1', labels)
```

```{code-cell} ipython3
fret_traj.values
```

Create a Plumed object from the loaded structure and the ACVs. Look for neighboring phosphorus atoms within a range of 15 A.

```{code-cell} ipython3
plumed = ft.restraints.Plumed(struct_nosolv, [don, acc], selection='name P', cutoff=15)
```

Set your desired FRET value that you want the simulation converge to. The corresponding $\langle R_{DA,E}\rangle$ will be computed from the $\langle E\rangle$. The $\langle R_{mp}\rangle$ will be calculated from $\langle E\rangle$ using a third order polynomial (Kalinin et al., *Nat. Methods*, 2012)

```{code-cell} ipython3
targetFRET = 0.6
mean_R_DA_E = ft.fret.mean_dist_DA_fromFRET(None,None,targetFRET,54)
targetRmp = ft.fret.R_DAE_to_Rmp(mean_R_DA_E)
```

Compile a plumed input file for the MD simulation

```{code-cell} ipython3
plumed.write_plumed('gmx/out/plumed_A1-U1.dat', targetRmp, 100, 100)
```

To illustrate the restaints save a visualization skript for VMD or PyMOL which can be run with either `vmd_vis -c filename.gro -x trajectory.xtc -v vis.vmd` or `pymol_vis -c filename.gro -x trajectory.xtc -v vis.py`

```{code-cell} ipython3
plumed.write_vmd('gmx/out/vis.vmd')
```

```{code-cell} ipython3
plumed.write_pymol('gmx/out/vis.py')
```

Save force field and topology files files of the mean position pseudoatoms

```{code-cell} ipython3
plumed.write_pseudo('gmx/out/MP.pdb', 'gmx/out/MP.itp')
```

```{code-cell} ipython3
don.save_mp('gmx/out/mn_riboswitch_D.dat', units='nm', format='plain')
acc.save_mp('gmx/out/mn_riboswitch_A.dat', units='nm', format='plain')
```

First, backup the original .gro file. Then, insert the mean position coordinates given in **MP.pdb** at positions specified by **riboswitch_D.dat** and **riboswitch_A.dat** 

```sh
mv em/mn_riboswitch.gro em/mn_riboswitch_original.gro
gmx insert-molecules -ci MP.pdb -f mn_riboswitch_original.gro -o mn_riboswitch.gro -ip out/mn_riboswitch_D.dat -replace water |& tee tmp_insert.dat
gmx insert-molecules -ci MP.pdb -f mn_riboswitch.gro -o mn_riboswitch.gro -ip out/mn_riboswitch_A.dat -replace water |& tee -a tmp_insert.dat
```

Run the following command will tell you how many solvent molecules have been replaced by the MP pseudoatoms and how to update the topology. 

```
awk '$1 == "Replaced" {sum += $2}; END {print "\nIn your topology.top file under the section [ molecules ] do the following:\n(1) decrease the number of solvent molecules by " sum "\n(2) add the following line:\nMP        2"}' tmp_insert.dat && rm tmp_insert.dat
```

+++

not needed
```sh
mkdir md0
mdp_dir=mdp
structureName=mn_riboswitch_s4_cleaned
plumedFile=plumed_A1-U1.dat

gmx grompp -f "$mdp_dir"/md0.mdp -c npt/"$structureName".gro -p "$structureName".top -o md0/"$structureName".tpr -po md0/"$structureName".mdp
gmx mdrun -v -s md0/"$structureName".tpr -c md0/"$structureName".gro -x md0/"$structureName".xtc -cpo md0/"$structureName".cpt -e md0/"$structureName".edr -g md0/"$structureName".log -plumed "$plumedFile"
```

```{code-cell} ipython3
single_run
```
