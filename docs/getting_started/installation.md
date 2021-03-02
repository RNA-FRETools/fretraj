# Installation

## 1. Install PyMOL and FRETraj
````{tabbed} For Windows
- Get the latest binaries of PyMOL from [Schrödinger](https://pymol.org/)
- Start the **Anaconda prompt** which comes bundled with PyMOL and run the following command to install *FRETraj* from the `fdsteffen` channel on Anaconda.

    ```
    conda install fretraj -c fdsteffen
    ```

    ```{note}
    Installing FRETraj with `conda` will also install the C++ library [LabelLib](https://github.com/Fluorescence-Tools/LabelLib). This allows for faster computation of the ACVs. However, *FRETraj* also features a built-in Python implementation of the ACV algorithm.
    ```
````

`````{tabbed} For Linux and macOS

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
    conda install fretraj -c fdsteffen
    ```
    ```` 

    ````{tabbed} from source
    Install the latest development version from **Github** (preferentially, install with [Poetry](https://python-poetry.org/))
    ```
    git clone https://github.com/RNA-FRETools/fretraj.git
    cd fretraj/
    poetry install   # install in editable mode

    # alternatively, install in non-editable mode with pip
    pip install .  
    ```
    ````    
`````

## 2. Register the Plugin
Install the *FRETraj* GUI with PyMOL's Plugin manager: `Plugin` &rarr; `Plugin manager` &rarr; `Install New Plugin` &rarr; `Choose file...` and select the `fretraj_gui.py` file located under `src/fretraj/`. Upon first startup FRETraj will prompt you to select a root directory where to store the calculated ACVs and parameter files.

```{tip}
Within FRETraj you can load a **demo project** by going to `Help` &rarr; `Load Example`. You may also want to have a look at this [step-by-step tutorial](acv_calculation).
```

```{figure} ../images/pymol_plugin_manager.png
---
width: 450px
name: pymol_plugin_manager
---
PyMOL's plugin manager on Windows.
```


