# Installation

## 1. Install PyMOL and FRETraj
````{tabbed} For Windows
- Get the latest binaries of PyMOL from [Schrödinger](https://pymol.org/)
- Search for **Anaconda prompt** in the Windows start menu and run the following command to install *FRETraj* from the `RNA-FRETools` channel from anaconda.org.

    ```
    conda install fretraj -c rna-fretools -c conda-forge
    ```
- Locate the installation directory by running

    ```
    fretraj --path
    ```
````

``````{tabbed} For Linux and macOS

- Install the PyMOL binaries from [Schrödinger](https://pymol.org/) or directly from the Schrödinger Anaconda channel. 

    ```
    conda install -c schrodinger pymol-bundle
    ```

    ```{admonition} Open-source PyMOL
    You may also build PyMOL yourself from [source](https://github.com/schrodinger/pymol-open-source).
    ```
- Install *FRETraj* either from **PyPI** (with `pip`), **Anaconda** (with `conda`) or **Github** (clone or [download](https://github.com/RNA-FRETools/fretraj.git) and install from source with `poetry` or `pip`).
  

    ````{tabbed} pip
    Install from **PyPI**
    ```
    pip install fretraj
    ```
    ````

    ````{tabbed} conda
    Install from **Anaconda**
    ```
    conda install fretraj -c rna-fretools -c conda-forge
    ```
    ```` 

    `````{tabbed} from source
    Install the latest development version from **Github** with [Poetry](https://python-poetry.org/))
    

    ```
    curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
    git clone https://github.com/RNA-FRETools/fretraj.git
    cd fretraj/
    poetry install
    ```
    `````
- Locate the installation directory by running

    ```
    fretraj --path
    ```
``````


## 2. Register the Plugin
- Start PyMOL and install the *FRETraj* GUI with PyMOL's Plugin manager: `Plugin` &rarr; `Plugin manager` &rarr; `Install New Plugin` &rarr; `Choose file...` and select the `fretraj_gui.py` file located in the directory that was issued by `fretraj --path`. In the popup window select where you would like to install the plugin (default: `<PyMOL_path>/site-packages/pmg_tk/startup/`). Confirm with OK.

```{figure} ../images/pymol_plugin_manager.png
---
width: 450px
name: pymol_plugin_manager
---
PyMOL's plugin manager on Windows.
```

`````{admonition} Accelerated ACV calculation with LabelLib
For faster computation of the accessible-contact volumes you may additionally consider installing the C++ library [LabelLib](https://github.com/Fluorescence-Tools/LabelLib) with either `pip` or `conda`. However, is not strictly required since *FRETraj* features a built-in Python implementation of the ACV algorithm.

````{tabbed} pip
```
pip install labellib
```
````

````{tabbed} conda
```
conda install labellib -c tpeulen
```
````
`````

