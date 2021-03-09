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

```{code-cell} ipython3
from matplotlib import pyplot as plt
import mdtraj as md
import fretraj as ft
import os
example_dir = '../../src/fretraj/examples/'
```

We load again a DNA double helix but this time from a 1$\,\mu$s MD trajectory along with the labeling parameters.

```{code-cell} ipython3
traj = md.load(os.path.join(example_dir+'DNA.xtc'), 
               top=os.path.join(example_dir+'DNA.pdb'),
               stride=10)
labels = ft.cloud.labeling_params(os.path.join(example_dir+'DNA_labels.json'), verbose=False)
print(f'timestep: {traj.timestep/1000 :.0f} ns')
print(f'length: {traj.time[-1]/1000 :.0f} ns')
```

Next, we calculate ACVs along the trajectory (here every 100 ns, for simplicity).

```{code-cell} ipython3
:tags: [remove-output]

frames = range(traj.n_frames)
acv_D = ft.cloud.Volume.from_frames(traj, 'D-DT23-C7', labels, frames)
acv_A = ft.cloud.Volume.from_frames(traj, 'A-DT31-C7', labels, frames)
```

Likewise, we compute a mean FRET efficiency per snapshot and combine them into a FRET trajectory.

```{code-cell} ipython3
fret = ft.cloud.FRET.from_volumes(acv_D, acv_A, 'Cy3-Cy5', labels)
fret_traj = ft.cloud.Trajectory(fret, timestep=traj.timestep)
fret_traj.dataframe
```

Launch **Binder** to visualize the multi-ACV trajectory.

```{code-cell} ipython3
acv_D_traj = ft.cloud.create_acv_traj(acv_D)
acv_A_traj = ft.cloud.create_acv_traj(acv_A)
ft.jupyter.nglview_trajectory_ACV(traj, acv_D_traj['FV'], acv_A_traj['FV'], acv_D_traj['CV'], acv_A_traj['CV'])
```
