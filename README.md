# Acoustic T-Matrix Computation with Comsol

Compute the acoustic T-Matrix of arbitrarily shaped objects. The
embedding medium must be reciprocal. If the geometry allows it, for example
rotationally symmetric objects, computations in a reduced number of dimensions are
possible.

## Publication

When using these models please cite:

## Warning

COMSOL uses the time evolution $`\exp(i \omega t)`$. This means that all values,
especially the material property parameters, must be complex conjugated with respect to
the opposite convention $`\exp(-i \omega t)`$ that is often used in physics.

Please note that this difference has been taken into account, so that the models
calculate the T-matrix entries with the physics convention.

## Comsol Version Requirement

These models are implemented in COMSOL 6.2. They can be analogously impemented 
in all versions starting from 5.5. However, in version 5.4 and older one has to
define the associated Legendre polynomials manually. 

## General Usage

Define your object, your materials (see warning), and the parameters for the
decomposition. The latter includes the modes that are used in the calculation and
`r_decomp` the radius of the decomposition. The sphere (or cylinder for cylindrical
T-matrices) defined by this radius must completely enclose the object of interest.
Additionally, define the radius of the domain `r_domain >= r_decomp` and the thickness
of the perfectly matched layers (PML) `d_pml` around the domain. These three parameters
define the first object in the geometry, which uses the first material for its
properties. The first object and material are later used for the embedding material in
the decomposition of the scattered fields.

As computation object, all of the examples use either a sphere or an infinitely extended
cylinder defined by `r_obj`. Change the geometry or the parameters to your liking.

There are two types of the models based on the incident field formulation.

### Background field

This is the general formulation that apply to any boundary conditions at scatterer-embedding
interface. It requires you only to run Study 1, which you can adjust to your needs. Namely,
you can define a frequency range and other parametric sweeps. The study will sweep over 
all possible incident fields for a given `lmax` or `mmax`, being a single multipolar wave, 
and solve the scattering problem. Then you can get the results via *Global Evaluation* 
in the *Results* section, returning the variables `ap`, which are the desired T-matrix entries. 

### Fast sweep models

If the scatterer is homogeneous and fluid (no shear waves), it becomes possible to speed up
the calculations by defining the monopole and dipole sources instead of the background field.
In this case, one needs to run two studies.

The first study is the same as for the background field formulation. 
The second study is only used as an auxiliary for the post-processing. It does not solve
anything, but evaluates the necessary integrals and computes the T-matrix entries. Make
sure to mirror the frequency sweep of the first study and also all parameter sweeps, you
added to the first study. Run this second study. Then you can get the results via
*Global Evaluation* of the second study in the *Results* section, returning the
variables `ap`.

# Acoustic T-Matrix

The acoustic T-Matrix of a 3D finite object uses scalar spherical wave (SSW) functions as a basis set. For
non-reciprocal media, the eigensolutions get much more complicated, so they cannot be
included in a simple way (for the embedding medium). The same is true for anisotropy, so
the embedding medium must be fluid (no shear waves), homogeneous, isotropic and reciprocal. An additional
constraint is, that the object must have a finite size.

There are two types of geometries for the T-matrix calculation, `3D` and
`axisym` for general three-dimensional objects and for axisymmetric
objects, respectively. If applicable, the latter is more accurate, since the geometry
that needs to be solved is only two-dimensional and the diagonality of the T-matrix with
respect to `m` is enforced.

## Math

We begin with defining the spherical harmonics of a degree $`l`$ and an order $`m`$

```math
Y_{lm}(\theta, \varphi)
= \sqrt{\frac{2l + 1}{4 \pi} \frac{(l - m)!}{(l + m)!}}P_{l}^m(\cos \theta)\mathrm{e}^{\mathrm{i} m \varphi}\,,
```

which produce the scalar spherical waves as

```math
\Psi^{(n)}_{lm}(kr, \theta, \varphi) = z^{(n)}_l(k r) Y_{lm}(\theta, \varphi)\,.
```
where $`k = \omega/c`$ is the wavenumber of acoustic waves in the embedding.
The incident wave uses spherical Bessel functions of the first kind ($`n = 1`$). 
The scattered wave is expressed with spherical Hankel functions of the first kind ($`n = 3`$).

We can decompose the scattered pressure field as

```math
p_{\text{sca}}(kr, \theta, \varphi)
= \sum_{l,m} a_{lm} \Psi_{lm}^{(3)}(kr, \theta, \varphi)\,.
```

By using the following equation:

```math
\int \mathrm d\Omega Y_{lm}^\ast(\theta, \varphi)
\Psi_{lm}^{(n)}(kr, \theta, \varphi)
= z^{(n)}_l(k r)\,,
```

the coefficients $`a_{lm}`$ can be obtained by projecting the scattered field
onto different modes [[Tsimokha et al. PRB 105, 165311 (2022)](https://doi.org/10.1103/PhysRevB.105.165311)]

```math
a_{lm} = \frac{1}{r_d^2 h^{(1)}_l(k r_d)}\int \mathrm{d} S Y_{lm}^\ast(\theta, \varphi)
p_{\text{sca}}(kr_d, \theta, \varphi))\,,
```
where the integration is carried out over the spherical surface of a radius $`r_d`$.

In the case of an axisymmetric problem, the integral above can be simplified 
to a contour one [[Ustimenko et al., APL 126, 142201 (2025)](https://doi.org/10.1063/5.0257760)]
```math
a_{lm} = \frac{2 \pi}{r_d h^{(1)}_l(k r_d)}\int \mathrm d l \sin \theta Y_{lm}^\ast(\theta, 0)
p_{\text{sca}}(kr_d, \theta, 0)\,.
```

# Acoustic cylindrical T-Matrix

Two-dimensional (or cylindrical) acoustic T-matrices use scalar cylindrical wave (SCW) functions
as a basis set. The difference to the previous case is that the SCWs allow the objects to be infinitely 
extended in the z-direction, either being uniform or periodic along this axis.

Equivalently to the T-matrix case, the object can either have a general shape or be axisymmetric. The
file is for infinitely extended uniform objects. This has, like the axisymmetric
case, only a two-dimensional computation domain. Again, axisymmetry enforces diagonality
with respect to `m`.

## Math

The cylindrical vector waves are defined as

```math
\Psi^{(n)}_{k_z,m}(k_{\rho}\rho,\varphi,k_z z) = Z^{(n)}_m(k_{\rho}\rho) \mathrm{e}^{\mathrm{i} m \varphi + \mathrm{i} k_z z}\,,
```
where $`k_\rho = \sqrt{k^2 - k_z^2}`$ and $`Z_m^{(n)}`$ are the Bessel or first-kind Hankel
functions. Similar to the spherical-wave case, we expand the scattered field using the
Hankel functions.
