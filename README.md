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
