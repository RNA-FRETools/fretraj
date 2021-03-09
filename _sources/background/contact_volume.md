# Dye interactions with the biomolecule

Due to their conjugated $\pi$-system fluorophore often interact with the biomolecule. Such surface stacking changes the mean position of the dye with respect to the biomolecule.

```{figure} ../images/Cy3_MD_stacking.png
---
width: 300px
name: Cy3_MD_stacking
---
A Cy3 fluorophore is coupled covalently to a DNA double helix.
```

The fluorophore no longer rotates isotropically but displays an anisotropic distribution around the biomolecule. This is reflected in the fluorescence anisotropy by an incomplete depolarization of the emitted light resulting in a decay that levels off at $r_\infty>0$.

```{figure} ../images/Cy3_anisotropy.png
---
width: 450px
name: Cy3_anisotropy
---
Fluorescence anisotropy decays with (a) weak surface interaction and (b) strong stacking propensity.
```

 To account for the fraction of dyes that are sticking to the biomolecule we introduced the **contact volume** (CV) as an subdivision area of the accessible volume. {cite}`Steffen.2016` The **width** and **weight** of the CV are two additional parameters that go into the ACV calculation. The width is usually set to match one of the dye dimensions. The weight attributed to each point in the CV can be estimated from the residual anisotropy at $r_\infty$.

```{figure} ../images/Cy3_ACV.png
---
width: 300px
name: Cy3_ACV
---
The accessible-contact volume (ACV) model divides the AV into a contact and a free volume. The contact volume is occupied when the dye sticks to the surface of the biomolecule.
```
