What is FRETraj?
================

.. toctree::
   :maxdepth: 2

.. role:: raw-html(raw)
   :format: html

FRETraj is a high-level Python API to the **LabelLib** library (https://github.com/Fluorescence-Tools/LabelLib) to simulate fluorophores which are coupled to a biomolecule of interest. The package features a user-friendly **PyMOL plugin** which can be used to explore different labeling positions while designing new FRET experiments. In an AV simulation the fluorophore distribution is estimated by a shortest path search (Djikstra algorithm) using a coarse-grained dye probe. FRETraj further implements a **Python-only** version of the geometrical clash search used in LabelLib. This is particularly useful for trying out new features to improve the description of the dye distribution.  

.. image:: _static/preQ1_ACV.png
   :width: 300
   :align: center

Installation
************

FRETraj can be installed platform independently as a Python module or as a Plugin to PyMOL. If you would like to predict FRET efficiencies from a single PDB file (including multiple states), then the PyMOL Plugin provides you with the easiest and fastest access to FRETraj. If you like to integrate FRETraj in your Python analysis workflow (e.g. to calculate multiple accessible volumes along an molecular dynamics trajectory, hence the name), then install FRETraj as a module.

PyMOL plugin
------------

You can get the latest version of PyMOL from [Schrödinger](https://pymol.org/). Start the **Anaconda prompt** which comes bundled with PyMOL 2.x and install the necessary dependencies. ::

    conda install numpy "numba<=0.44" mdtraj packaging -c conda-forge

For a faster calculation of the AVs you may additionally install **LabelLib**, but this is not required as FRETraj also runs its own implementation of the AV algorithm. ::

    conda install -c tpeulen labellib
    or 
    pip install git+https://github.com/fdsteffen/labellib.git@src-neighbour-nodes


To use the **FRETraj PyMOL plugin** simply download the .zip archive from Github and install it via PyMOL's Plugin manager: ``Plugin`` :raw-html:`&rarr;` ``Plugin manager`` :raw-html:`&rarr;` ``Install New Plugin`` :raw-html:`&rarr;` ``Choose file...`` and select the .zip archive. Upon first startup FRETraj will prompt you to select a root directory where to store the calculated ACVs and parameter files.

.. image:: _static/PyMOL_interface.PNG
   :height: 200
   :align: left

.. image:: _static/PyMOL_Plugin.PNG
   :height: 200
   :align: right


Python module
-------------

You may install the latest development version of FRETraj as a Python module from Github using pip ::

    pip install -e git+https://github.com/fdsteffen/fretraj.git



Getting started
***************

A good way to start is to visit one of the :doc:`tutorials <tutorial/single/accessible_contact_volume>`