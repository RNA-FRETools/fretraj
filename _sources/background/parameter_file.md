# The parameter file

*FRETraj* uses a json-formatted parameter file to define the settings for the ACV calculation. The general structure is as follows:

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

The parameters above need to be defined by the user. Other parameters have default values:
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

Further details about the simulation settings are provided in the [code documentation](../module/cloud).
