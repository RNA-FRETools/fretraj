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

# Generating photon bursts

*FRETraj* predicts mean FRET efficiencies and distributions thereof for dynamic biomolecules as outlined in the
previous sections. FRET histograms of single-molecule experiments are often broadened due to shot-noise. 
For better comparison of *in vitro* and *in silico* FRET measurements, *FRETraj* can take the **photon noise** into account by simulating fluorescence emission events. The probabilities of donor and acceptor emission are dependent on the quantum yields and fluroescence lifetimes of the dyes as well as the transfer efficiency and thus the distance between their ACVs {cite}`Hoefling.2011, Hoefling.2013`.
This notebook show how to simulate **photon bursts** similar to a confocal single-molecule experiment.

```{code-cell} ipython3
import fretraj as ft
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
```

First, we load a parameter file for the burst simulation. The format of this file is described [here](../background/parameter_file.html#burst-simulation)

```{code-cell} ipython3
:tags: [hide-output]

parameters = ft.burst.readParameters('burst_data/burst_parameters.json')
parameters
```

Importantly, key `species.unix_pattern_rkappa` in the parameter file points to any file matching the given regular expression. Here, the file `R_kappa.dat` is created from a `ft.cloud.Trajectory` object (see [Working with Trajectories](trajectories.md)) and contains inter-dye distance $R_\text{DA}$(t) and $\kappa^2$ values.

```{code-cell} ipython3
:hide-input: null

pd.read_csv('burst_data/R_kappa.dat', sep='\t', names=['R_DA (nm)', 'kappasquare']).head()
```

An analytical burst size distribution $P(x)$ is specified as a power law with a coefficient $\lambda$

$$P(x) = x^\lambda$$ 

Here we set the expoenent to $\lambda=-2.3$. We can now start a burst experiment.

```{code-cell} ipython3
experiment = ft.burst.Experiment('burst_data/', parameters)
```

The resulting FRET histogram is broadened by shot-noise.

```{code-cell} ipython3
:tags: [remove-cell]

sns.set_style('white')
sns.set_context('notebook')

def set_ticksStyle(x_size=4, y_size=4, x_dir='in', y_dir='in'):
    sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir, 'ytick.direction': y_dir})
```

```{code-cell} ipython3
with sns.axes_style('ticks'):
    set_ticksStyle()
    f, ax=plt.subplots(nrows=1, ncols=1, figsize=(3, 2), sharex=True, sharey=True, squeeze=False)
    ax[0, 0].hist(experiment.FRETefficiencies, bins=25, range=[0, 1], color=[0.75, 0.51, 0.38])
    ax[0, 0].set_xlabel('FRET')
    ax[0, 0].set_ylabel('occurence')
```

Launch **Binder** 🚀 to interact with this notebook.
