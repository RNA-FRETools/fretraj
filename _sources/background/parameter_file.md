# The parameter files

*FRETraj* uses two different json-formatted parameter files: one for **ACV/FRET calculation** and another to simulate **photon bursts**. Their general structure is outlined in the following:


## ACV/FRET calculation

```
{"Position":
    {"<string defining the dye position>": 
        {"attach_id": int,
            "linker_length": float,
            "linker_width": float,
            "dye_radius1": float,
            "dye_radius2": float,
            "dye_radius3": float},
    },
    "Distance": {"<string defining the fret pair>": 
        {"R0": float}
    }
}
```

The parameters above are mandatory and need to be defined by the user. Other parameters have default values:
```
{"Position":
    {"<string defining the dye position>": 
        {"cv_thickness": 0,
            "cv_fraction": 0.5,
            "simulation_type": "AV3",
            "grid_spacing": 1.0,
            "mol_selection": "all",
            "state": 1,
            "frame_mdtraj": 0,
            "use_LabelLib": True}
    },
    "Distance": {"<string defining the fret pair>": 
        {"n_dist": 10E+6}
    }
}
```

The parameters are described in more detail in the sections [Calculating ACVs](../getting_started/acv_calculation) and [making FRET predictions](../getting_started/acv_calculation)

```{hint}
The distance unit of FRETraj is **Angstrom** ([MDTraj](http://mdtraj.org/) in contrast uses nanometers). The time unit is **picosecond**.
```

```{note}
- If `simulation_type: "AV1"` the dye is only defined by a single radius `dye_radius1`. 
- `mol_selection` takes a valid atom selection expression compatible with [MDTraj](http://mdtraj.org/)
```

Finally, some parameters are solely used for visualization in PyMOL:

```
{"Position":
    {"<string defining the dye position>": 
        {"contour_level_AV": 0,
            "contour_level_CV": 0.7,
            "b_factor": 100,
            "gaussian_resolution": 2,
            "grid_buffer": 2.0,
            "transparent_AV": True}
    }
}
```

Further details about the simulation settings are provided in the [cloud](../module/cloud) and [fret](../module/fret) submodule descriptions.

<a id='burst-simulation'></a>

## Burst simulation

The structure of the second parameter file is as follows:

```
{"dyes": 
    {"tauD": float,
     "tauA": float,
     "QD": float,
     "QA": float
    },
"fret": 
    {"R0": 54},
"species":
    {"name": [string],
     "unix_pattern_rkappa": [<regex-string>],
     "probability": [float]
    },
"bursts": 
    {"lower_limit": 15,
     "upper_limit": 150,
     "lambda": -2.3,
     "averaging": "all"
    }
}
```

The above key-values are mandatory. The main keys define photopyhsical parameters (`dyes`), the FRET pair (`fret`), identifiers and weights of individual subpopulations (`species`) and characteristics of the burst size distribution (`bursts`). The individual parameters are briefly described below:
- `dyes`
    - `tauD` / `tauA`: fluorescence lifetime of the donor and acceptor dye
    - `QD` / `QA`: quantum yields of the donor and acceptor dye
- `fret`
    - `R0`: Förster radius in Angstrom
- `species`
    - `name`: list of names of the different subpopulations
    - `unix_pattern_rkappa` list of regular expressions identifying the R_kappa-files of the subspopulations 
    - `probability`: list of probabilities of the subpopulations
- `bursts`
    - `lower_limit` / `upper_limit`: lower and upper treshold of the burst size distribution
    - `lambda`: coefficient of the analytical burst size distribution (power law)
    - `averaging`: method of averaging over the MD simulations, available options are: `ensemble` (the entire ensemble is present in a single burst), 
    `species` (each burst only contains molecules of the same species) or `trajectory` (only a one trajectory per burst).

    ```{hint}
    The `averaging` parameter becomes important if multiple trajectories are run with **different starting conformations**, which are not expected to exchange on the timescale of the simulation
    (e.g. cis/trans isomerization of prolines, see {cite}`Hoefling.2011, Hoefling.2013`). If on the other hand, only a single trajectory is available all three settings should give the same results (although `ensemble` may run a bit faster than the other two). 
    ```

 Additional, optional parameters and their default values include:

```
{"dyes": 
    {dipole_angle_abs_em: 0},
"sampling":
    {"nbursts": 2000,
     "skipframesatstart": 0,
     "skipframesatend": 1000,
     "multiprocessing": True
    },
"fret": 
    {"kappasquare": 0.6666
     "no_gamma": False,
     "quenching_radius" 1,
    }
"species":
    {"unix_pattern_don_coords": None,
     "unix_pattern_acc_coords": None,
     "n_trajectory_splits": None
    },
"bursts": 
    {"burst_size_file": None,
     "QY_correction": False,
    }
}
```

```{note}
The `fretraj.burst` submodule can also be used to analyze MD simulations with **all-atom dyes**. Including the dyes explicitely in the simulations has the advantage that a **transition dipole vector** 
can be defined. This allows to additionally simulate **fluorescence anisotropy decays**, which can be activated by setting `compute_anisotropy=True` in `fretraj.burst.experiment()`.
```

- `dyes`
    - `dipole_angle_abs_em`: angle between the excitation and emission dipole in the absence of rotational diffusion (i.e. at time point 0; reduces the fundamental anisotropy $r_0$; for all-atom dye simulations only)
- `sampling`:
    - `nbursts`: number of bursts to generate
    - `skipframesatstart` / `skipframesatend`: number of frames (*not* time) to skip at the beginning and end of each trajectory when searching for excitation time points 
    - `multiprocessing`: distribute workload on multiple cores
- `fret`
    - `kappasquare`: $\kappa^2$ value for the given `R0` (usually isotropic averaging, i.e. 2/3)
    - `no_gamma`: mimic an uncorrected FRET experiment which is affected by the different quantum yields of donor and acceptor (i.e. before $\gamma$-correction).
    The detection efficiency ratio is always set to 1. If the simulation is compared to a $\gamma$-corrected experiment this parameter should be set to `False`.
    - `quenching_radius`: radius below which radiationless deactivation occurs due to contact quenching of the dyes
- `species`
    - `unix_pattern_don_coords` / `unix_pattern_acc_coords`: regular expression identifying the xyz-coordinates (xvg-files generated by GROMACS `gmx traj`) of the donor and acceptor dye (required for all-atom dye simulations)
    - `n_trajectory_splits` split the trajectory into $n$ parts (together with `averaging="trajectory"` this simulates non-interconverting species which leads to broadening; usually this should be left at the default: `None`)
- `bursts`
    - `burst_size_file`: relative path to an experimental burst size distribution file containing two columns: burst sizes and their probabilities). 
    If this setting is not specified, then an analytical power law distribution defined by the coefficient `lambda` is used.
    - `QY_correction`: correct the number of donor and acceptor photons by their respective quantum yield when evaluating whether the burstsize is reached. Usually experimental bursts are identified based on uncorrected photon numbers, so this setting defaults to `False`
        ```{admonition} Example for illustration
        Let's assume we collected 20 donor photons (with `QD=0.5`) and 25 acceptor photons (with `QA=0.75`) and the size of this particular burst is 50.
        If `QY_correction=False` the total photon count is 20+25=45 photons, so the required burstsize is not reached yet.
        if `QY_correction=True` the total photon count is (20/0.5 + 25/0.75 = 73) and so the burst is complete. 
        ```
    
