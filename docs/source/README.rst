.. toctree::
   :maxdepth: 2

.. role:: raw-html(raw)
   :format: html

FRETraj - Predicting FRET *in silico* 
*************************************

FRETraj is a high-level Python API to the **LabelLib** library (https://github.com/Fluorescence-Tools/LabelLib) to simulate fluorophores which are coupled to a biomolecule of interest. The package features a user-friendly **PyMOL plugin** which can be used to explore different labeling positions while designing new FRET experiments. In an AV simulation the fluorophore distribution is estimated by a shortest path search (Djikstra algorithm) using a coarse-grained dye probe. FRETraj further implements a **Python-only** version of the geometrical clash search used in LabelLib. This is particularly useful for rapid prototyping of new features of the ACV algorithm.

.. image:: _static/graphical_abstract.png
   :width: 90%
   :align: center

.. _pymol_plugin_installation:

Plugin Installation
*******************

You can get the latest version of PyMOL from `Schrödinger <https://pymol.org/>`_. Start the **Anaconda prompt** which comes bundled with PyMOL 2.x and install the necessary dependencies. ::

    conda install numpy "numba<=0.44" mdtraj packaging -c conda-forge

For a faster calculation of the AVs you may additionally install `LabelLib<https://github.com/Fluorescence-Tools/LabelLib>`_, but this is not required as FRETraj also runs its own implementation of the AV algorithm. ::

    conda install -c tpeulen labellib

To use the **FRETraj PyMOL plugin** simply download the .zip archive from Github and install it via PyMOL's Plugin manager: ``Plugin`` :raw-html:`&rarr;` ``Plugin manager`` :raw-html:`&rarr;` ``Install New Plugin`` :raw-html:`&rarr;` ``Choose file...`` and select the .zip archive. Upon first startup FRETraj will prompt you to select a root directory where to store the calculated ACVs and parameter files.

.. image:: _static/PyMOL_interface.PNG
   :height: 275
   :align: left

.. image:: _static/PyMOL_Plugin.PNG
   :height: 275
   :align: right


Getting started
***************

A good way to start is to have a look at the tutorial :doc:`tutorials <pymol_plugin>`

.. [#] PyMOL is a trademark of Schrodinger, LLC.