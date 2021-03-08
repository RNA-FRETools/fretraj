# Calculating ACVs

Accessible contact volumes and FRET calculations can be done directly from PyMOL. The PyMOL plugin allows the user to visually interact with the molecular viewer while exploring optimal sites for labeling a biomolecule of interest. This tutorial will guide you step-by-step through the workflow.

Start by opening the FRETraj GUI from the `Plugin Menu`.

```{note}
Upon first startup, *FRETraj* asks the user to define a root directory where all output files (ACVs, label and FRET parameters) are stored 
```

```{figure} ../images/PyMOL_interface.png
---
height: 240
name: pymol_interface
---
PyMOL with the *FRETraj*'s graphical user interface (left) and molecular viewer (right).
```

The graphical user interface (GUI) is divided into two sections (1) **Accessible contact clouds** and (2) **FRET trajectory**. The labeling positions are defined in the panel *PDB structure and position* . The panels below specify the *dye*, *linker*, *simulation* and *contact volume* parameters. All these settings are saved into the [parameter file](../background/parameter_file).

In the following the most important settings and functionalities of FRETraj for simulating accessible-contact volumes and predicting FRET efficiencies are described. 

## Load PDB and define dye positions 

1. **Set folder**: define the root folder where the ACVs and the parameter file are saved
2. **Load PDB**: import a single or multi-state PDB/CIF file
3. **Show Text**: Show the PDB text file in a integrated viewer to help with the atom ID selection
    ```{hint}
    If you already have a JSON-formatted [parameter file](../background/parameter_file) you can load this file with the **Load Param** button proceed directly to [step 12](#compute-acv).
    ```
4. **atom ID**: select the attachment atom on the biomolecule by its serial ID
5. **state**: select a different state of the PDB file
6. Add the new position to the label dropdown by pressing on ▶️
    ```{note}
    If you wish to remove a label position from the list press on **Delete label** (this does not remove any associated ACV PDB files from the root folder)
    ```
7. **mol. selection**: consider only a subregion of your biomolecule while computing ACVs. While this may speed up the calculation, usually you want to leave this at *all* (default) for simplicity.

## Define dye and linker dimensions

8. **Dye parameters**: Set the radii of the dye (see section [Accessible volume](../background/accessible_volume))
9. **Linker parameters**: Set the length and width of the carbon linker
10. **Simulation**: 
    - **grid spacing**: distance between two grid points
    - **simulation type**: select between a single-radius (AV1) or a three-radii (AV3) dye probe
    - **AV library**: if [LabelLib](https://github.com/Fluorescence-Tools/LabelLib) has been installed, you may select *use LabelLib* to speed up the calculations. Otherwise a pure [Python version](../module/grid) of the ACV algorithm is used.

11. **Contact volume**: Choose the thickness and relative fraction of the [contact volume](../background/contact_volume)
    <a name="compute-acv"></a>
    ```{note}
    The parameter file is dynamically updated whenever a new label is added/removed or a setting is changed.
    ```
12. **Compute ACV**: Start the calculation. The ACV is automatically saved to the selected root folder. Once the calculation has finished the ACV cloud is displayed in the molecular viewer and the *mean position (MP)* is indicated in the table