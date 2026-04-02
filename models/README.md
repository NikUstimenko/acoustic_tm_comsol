# Instructions for calculating acoustic T-matrices of scatterers with COMSOL Multiphysics

In calculating T-matrices, one repeatedly illuminates the scatterer with many
different sources. Meanwhile, the FEM equation system is typically solved by
first factoring the system matrix and then solving the resulting linear system.
The `fastsweep` models implement the calculation in such a way that COMSOL only has to
factor the system matrix once for each frequency, and can then solve the linear
system for many different sources. This speeds up the calculation, often by a
factor of 10 or more compared to the `background` models.

**Important note 1:** The models support materials with isotropic density,
speed of sound (pressure waves) and speed of shear waves.

**Important note 2:** While COMSOL uses the `exp(i omega t)` time convention,
the T-matrix coefficients coming out of the model follow the `exp(-i omega t)`
time convention.

## Modelling your own scatterer

* Construct the geometry. Change the domain size parameters as needed, but keep
  `r_decomp` smaller than `r_domain`. For 2D axisymmetric models, keep the domain
  size fairly large, about 2 wavelengths or more.

* Assign materials. The first material (`mat1`) is the domain material; you can 
  change its properties but do not change its name.

* Put the same frequency sweep in both study 1 and study 2.

* Change the `lmax` parameter to change the order up to which the multipole 
  coefficients are calculated.

* Check the mesh. Recall that elements should be smaller than 
  lambda0/5/n_material in dielectric-type materials, and much smaller in metals.
  The mesh on the integration (decomposition) surface can be especially important
  for accuracy.

* In the `fastsweep` models, run study 1 (the simulations) and study 2 (the decomposition)
  one after another. Finally, evaluate `ap` in Global Evaluation.

When building a new model, please consider testing for convergence with respect
to mesh and domain size.

Non-reciprocity parameter and anisotropy are not currently supported.

Sanity checks in case of trouble:

* Sanity check for perfectly matched layer (PML): PML should be only in the 
  outermost layer of the model. This is under Definitions - Artificial Domains.
 
* Sanity check for integration surface definition: Definitions - Integration 1
  (intop1) should select the spherical surface that encloses the scatterer but
  does not touch it or the PML boundary (`r_decomp` should be sligtly bigger
  than the scatterer). 

* Sanity check for the sources (only for the `fastsweep` models): ensure that
  the Monopole Domain Source and the Dipole Domain Source selects the scatterer
  and (preferably, for speed) nothing else.

## Standard T-matrix HDF5 output

To make a 'standard' HDF5 T-matrix file, you can use the `readmph_tmatrix.py`
python script. It requires the `mph` package
(https://mph.readthedocs.io/en/stable/index.html)
as well as `tmatrix_tools.py`. This automatically saves the model parameters,
materials etc. into the hdf5 file. To use the script, run the COMSOL model first 
and save it. Then run the command

* `python readmph_tmatrix.py modelfile.mph tmatrixfile.hdf5`

Currently, mesh export does not work. For 

## Implementation notes of the `fastsweep` models

To avoid using COMSOL's own scattered-field formulation, which is extremely slow 
with multipole fields, the scattered field formulation is re-implemented using 
domain sources. In this case, the governing equations include sources

*  `curl E = -dB/dt + Jmag`

*  `curl H = dD/dt + J`

Then, substituting `p = p_inc + p_sc` and `v = v_inc + v_sc`, we can derive the 
expressions for the scattering currents

* `Q = omega^2 (beta_obj - beta_domain) p_inc`

* `q = -i omega (rho_obj - rho_domain) v_inc`

The incident fields `p_inc` and `v_inc` are the regular multipole fields.

In the 2D axisymmetric model, COMSOL's PML does not work well for azimuthal 
phase indices m != 0. For this reason, the domain needs to be fairly large 
(2 wavelengths or more) to achieve good convergence. The rectangular domain 
shape also tends to work better than the spherical for mysterious reasons.
