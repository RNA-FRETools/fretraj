The accessible volume
=====================

.. toctree::
   :hidden:

.. |Steffen2016| image:: https://img.shields.io/badge/DOI-10.1039/C6CP04277E-blue.svg
  :target: https://doi.org/10.1039/C6CP04277E

.. |Hoefling2011| image:: https://img.shields.io/badge/DOI-10.1371/journal.pone.0019791-blue.svg
  :target: https://doi.org/10.1371/journal.pone.0019791

.. |Kalinin2012| image:: https://img.shields.io/badge/DOI-10.1038/nmeth.2222-blue.svg
  :target: https://doi.org/10.1038/nmeth.2222

Fluorophores are usually attached to a biomolecule of interest via a relatively flexible, aliphatic spacer arm. This carbon linker endows the fluorophore with motional freedom which is the basis for the common assumption of near isotropic dye rotation, expressed by an orientation factor :math:`\kappa=2/3`. At the same time, the flexibility of the dye linker introduces some uncertainty in pinpointing the average position of the dye with respect to the biomolecule. Molecular dynamics simulations with explicitly included fluorophores can sample the dye distribution around a biomolecule.[1,2] Yet, these simulations are fairly computationally expensive. The accessible volume (AV) approach mitigates this problem by estimating the dye distribution with a purely geometrical search algorithm. For this purpose the dye is approximated by either a single or alternatively three radii. This dye probe is set to explore the space around the biomolecule without clashing with the surface.[3] The maximal distance from the attachment position given by the linker dimensions. The possible dye positions are updated iteratively using a Dijkstra shortest path algorithm. The resulting point cloud describes the location of the dye relative to a static biomolecule.

.. Todo::
    
    Graphic illustrating the shortest path algorithm



[1] F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. |Steffen2016|

[2] M. Hoefling, N. Lima, D. Haenni, C.A.M. Seidel, B. Schuler, H. Grubmüller, *Plos One.* **2011**, *6*, e19791. |Hoefling2011|

[3] S. Kalinin, T. Peulen, C.A.M. Seidel et al. *Nat. Methods*, **2012**, *9*, 1218-1225. |Kalinin2012|
