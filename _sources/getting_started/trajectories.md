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

# Working with trajectories

Running *FRETraj* with Jupyter has the advantage that we can easily work with
trajectories and other **multi-model structures**. As we saw in the previous 
notebook, *FRETraj* uses `mdtraj` to load single and multi-frame objects. 
For demonstrative purposes, we work here with an extract of an MD trajectory 
from our DNA hairpin. A more detailed notebook with the full trajectory is 
available at [here](https://github.com/RNA-FRETools/FRETraj-demo).

```{code-cell} ipython3
from matplotlib import pyplot as plt
import mdtraj as md
import fretraj as ft
import os
example_dir = '../../src/fretraj/examples/'
```

We load multiple snapshots of the DNA hairpin from a 1$\,\mu$s MD trajectory along with the labeling parameters.

```{code-cell} ipython3
traj = md.load(os.path.join(example_dir+'DNA_hairpin.xtc'), 
               top=os.path.join(example_dir+'DNA_hairpin.pdb'),
               stride=10)
labels = ft.cloud.labeling_params(os.path.join(example_dir+'DNA_hairpin_labels.json'), verbose=False)
print(f'timestep: {traj.timestep/1000 :.0f} ns')
print(f'length: {traj.time[-1]/1000 :.0f} ns')
```

Next, we calculate ACVs along the trajectory (here every 100 ns, for simplicity).

```{code-cell} ipython3
:tags: [remove-output]

frames = range(int(traj.n_frames/2))
acv_D = ft.cloud.Volume.from_frames(traj, 'Cy3-20-C5', labels, frames)
acv_A = ft.cloud.Volume.from_frames(traj, 'Cy5-44-P1', labels, frames)
```

Likewise, we compute a mean FRET efficiency per snapshot and combine them into a FRET trajectory.

```{code-cell} ipython3
:tags: [remove-output]

fret = ft.cloud.FRET.from_volumes(acv_D, acv_A, 'Cy3-Cy5', labels)
fret_traj = ft.cloud.Trajectory(fret, timestep=traj.timestep)
```

```{code-cell} ipython3
fret_traj.dataframe
```

Launch **Binder** 🚀 to visualize the multi-ACV trajectory.

```{code-cell} ipython3
:tags: ['remove-output']

acv_D_traj = ft.cloud.create_acv_traj(acv_D)
acv_A_traj = ft.cloud.create_acv_traj(acv_A)
ft.jupyter.nglview_trajectory_ACV(traj, acv_D_traj['FV'], acv_A_traj['FV'], acv_D_traj['CV'], acv_A_traj['CV'])
```
