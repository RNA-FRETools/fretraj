PyMOL Plugin
============

The PyMOL plugin allows the user to interact with the molecular viewer while calculating accessible volumes. 

.. image:: _static/PyMOL_Plugin.PNG
   :width: 90%
   :align: center

The interface is divided into two sections (1) **Accessible contact clouds** and (2) **FRET trajectory**. Under the subpanel *PDB structure and position* the labeling positions are defined. The panels below define the *dye*, *linker*, *simulation* and *contact volume* parameters. All these settings are saved into the :doc:`parameter file <parameter_file>`.

The following steps will guide you through the workflow of calculating ACVs and predicting FRET efficiencies: 



ACV calculation
---------------

1. **Set folder**: define the root folder where the ACVs and the parameter file are saved
2. **Load PDB**: import a single or multi-state PDB/CIF file
3. **Show Text**: Show the PDB text file in a integrated viewer to help with the atom ID selection

.. Note ::
    
    If you already have a json formatted :doc:`parameter file <parameter_file>` you can load this file with the **Load Param** button proceed directly to step .


4. **atom ID**: select the attachment atom on the biomolecule by its serial ID
5. **state**: select a different state of the PDB file
6. Add the new position to the label dropdown by pressing on `>`

.. Note ::

    If you wish to remove a label position from the list press on **Delete label** (note: this does not remove any associated ACV-PDB file from the root folder)

7. **mol. selection**: consider only a subregion of your biomolecule while computing ACVs. While this may speed up the calculation, usually you want to leave select *all* for simplicity.
8. **Dye parameters**: Set the radii of the dye
9. **Linker parameters**: Set the length and width of the carbon linker
10. **Simulation**: 

    - **grid spacing**: distance between two grid points
    - **simulation type**: select between a single-radius (AV1) or a three-radii (AV3) dye probe
    - **AV library**: if the library has been :ref:`installed <pymol_plugin_installation>` you may select *use LabelLib* to speed up the calculations. Otherwise a pure :doc:`Python version<module_grid>` of the ACV algorithm is used.

11. **Contact volume**: Choose the thickness and relative fraction of the :doc:`contact volume<contact_volume>` 
12. **Compute ACV**: start the calculation (note: the first run can take slightly longer than subsequent ones). The calculated ACV is automatically saved to the selected root folder. Once the calculation has finished the ACV cloud is displayed in the molecular viewer and the *mean position (MP)* is indicated in the table

Predict FRET
------------

1. **R0**: Förster radius of the dye pair
2. **# distances**: number of distances to sample for the FRET calculation. The algorithm chooses *n* pairs of randomly distributed points in the donor and acceptor cloud and calculates their distance
3. choose the **donor** and **acceptor** cloud from the drop down
4. **Calculate FRET**: start the FRET calculation. The mean FRET efficiency :math:`\langle E\rangle`, the mean inter-dye distance :math:`\langle R_{DA}\rangle` as well as the distances between the mean dye positions :math:`\langle R_{MP}\rangle` or the two attachment sites are displayed in the table.


.. Note ::
    
    The style of the ACV clouds can be tuned with the *ACV visualization* settings:

        - **contour levels** of the accessible and contact volume
        - **b-factor** and the **gaussian resolution** define the smoothness of the cloud
        - **grid buffer** around the ACV
        - **transparency** of the cloud


.. toctree::
   :hidden:
