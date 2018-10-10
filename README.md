# acv
Compute accessible volume clouds on PDB structures

## Dependencies
acvCloud depends on the following python modules
- numpy
- Biopython
- Scipy
- itertools
- heapq
- argparse

The modules can be installed via the package manager `pip`, e.g. for argparse
```
pip install argparse
```

## Run acvCloud
acvCloud can be run from the command line by specifiying a pdb structure and an associated parameter file (see [Parameter file](#parameter-file))
```
acvcloud.py -i structure.pdb -p parameters.dat
```

or from within PyMOL
```
run acvcloud.py -i structure.pdb -p parameters.dat
```

## Parameter file
A typical parameter file for acvCloud has the following entries

```
[dye parameters]
serialID=257
CVthick=3
n=18

[grid]
spacing=1

[dijkstra]
blub=test
```

To accurately visualize the accessible, contact and free volume in PyMOL, set the VdW size of the volume objects (not the pdb structure) to half of the grid size. In the Python API the following commands will do the job (assuming a grid spacing of 1 Angstrom and the respective volume objects are named `RNA`, `AV`, `CV` and `FV`)
```
cmd.alter("RNA", "vdw=0.5")
cmd.alter("AV", "vdw=0.5")
cmd.alter("CV", "vdw=0.5")
cmd.alter("FV", "vdw=0.5")
```
