# Acoustic T-Matrix Computation with Comsol

Compute the acoustic T-Matrix of arbitrarily shaped objects. The
embedding medium must be reciprocal. If the geometry allows it, for example
rotationally symmetric objects, computations in a reduced number of dimensions are
possible.

## Publication

When using these models please cite:

## Warning

Comsol uses the time evolution $`\exp(i \omega t)`$. This means that all values,
especially the material property parameters, must be complex conjugated with respect to
the opposite convention $`\exp(-i \omega t)`$ that is often used in physics.

This difference has been taken into account in `readmph_tmatrix.py` using complex conjugation.

## Comsol Version Requirement

There is a major difference between Comsol 5.4 and Comsol 5.5 regarding these models: 
the associated Legendre polynomials are not implemented in the former one and should be defined manually. 
For Comsol 5.5 and above no such restrictions apply.

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

In the first study the actual solutions of the scattering problems are calculated. A
frequency sweep is assumed per default. Adjust it to your needs. You can define
additionally other parametric sweeps. Run the first study manually, if needed.

The second study is only used as an auxiliary for the post-processing. It does not solve
anything, but evaluates the necessary integrals and computes the T-matrix entries. Make
sure to mirror the frequency sweep of the first study and also all parameter sweeps, you
added to the first study. Run this second study. Then you can get the results via
*Global Evaluation* of the second study in the *Results* section, returning the
variables `ap` and `am`.

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
onto different modes

```math
a_{lm} = \frac{1}{h^{(1)}_l(k r_d)}\int \mathrm d\Omega Y_{lm}^\ast(\theta, \varphi)
p_{\text{sca}}(kr_d, \theta, \varphi))\,,
```
where the integration is carried out over the spherical surface of a radius $`r_d`$.

# Cylindrical T-Matrix

Two-dimensional (or cylindrical) acoustic T-matrices use scalar cylindrical wave (SCW) functions
as a basis set. The difference to the previous case is that the SCWs allow the objects to be infinitely 
extended in the z-direction, either being uniform or periodic along this axis.

There are three java files for the calculation, `tmatrixc.java`, `tmatrixc_axisym.java`,
and `tmatrixc_uni.java`. The first two of them are for periodic objects. Equivalently to
the T-matrix case, the object can either have a general shape or be axisymmetric. The
last file is for infinitely extended uniform objects. This has, like the axisymmetric
case, only a two-dimensional computation domain. Again, axisymmetry enforces diagonality
with respect to `m`.

In the periodic cases only `kz` values that differ by a multiple of the reciprocal
lattice vector are considered in the model, the number of reciprocal lattice vectors
included in each direction is `n_kz`. For uniform cases, different `kz` values do not
couple.

## Math

The cylindrical vector waves are defined as

```math
\boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z)
= \mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\left[
\mathrm i m \frac{Z_m^{(n)}(k_\rho \rho)}{k_\rho \rho} \boldsymbol{\hat\rho}
- Z_m^{(n)\prime}(k_\rho \rho) \boldsymbol{\hat\varphi}
\right] \\
\boldsymbol N_{mk_z}^{(n)}(\rho, \varphi, z)
= \frac{\nabla}{k} \times \boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z ) \\
\boldsymbol A_{m k_z p}^{(n)}(\rho, \varphi, z)
= \frac{1}{\sqrt{2}} \left(\boldsymbol N_{mk_z}^{(n)}(\rho, \varphi, z )
+ p \boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z ) \right)\,.
```

where $`k_\rho = \sqrt{k^2 - k_z^2}`$ and $`Z_m^{(n)}`$ are the Bessel or Hankel
functions. The projection onto different modes in the case of chiral media, and
therefore different $`k`$ and $`k_\rho`$ values, is not as straightforward as in the
cylindrical case. Similar to the spherical case we expand the scattered wave using
Hankel functions. We use the integrals

```math
I_1
= \frac{1}{2\pi a_z}
\int_0^{2\pi} \mathrm d \varphi
\int_0^{a_z} \mathrm d_z
\mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\boldsymbol{\hat z} \boldsymbol E_{\text{sca}}(\boldsymbol r)
= \frac{1}{\sqrt 2}
\left(
\frac{k_{\rho+} H_m^{(1)}(k_{\rho+} \rho)}{k_{\rho+} \rho} a_{mk_z +}
+ \frac{k_{\rho-} H_m^{(1)}(k_{\rho-} \rho)}{k_{\rho-} \rho} a_{mk_z -}
\right) \\
I_2
= \frac{1}{2\pi a_z}
\int_0^{2\pi} \mathrm d \varphi
\int_0^{a_z} \mathrm d_z
\mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\left(
-\mathrm i m \frac{H_m^{(1)}(k_{\rho+} \rho)}{k_{\rho+} \rho} \boldsymbol{\hat\rho}
- H_m^{(1)\prime}(k_{\rho+} \rho) \boldsymbol{\hat\varphi}\right)
\boldsymbol E_{\text{sca}}(\boldsymbol r) \\
= \frac{1}{\sqrt 2}
\left[
H_{m+1}^{(1)}(k_{\rho+} \rho) H_{m-1}^{(1)}(k_{\rho+} \rho) a_{mk_z +}
- \frac{H_{m+1}^{(1)}(k_{\rho+} \rho) H_{m-1}^{(1)}(k_{\rho-} \rho)(k_- - k_z)
+ H_{m+1}^{(1)}(k_{\rho-} \rho) H_{m-1}^{(1)}(k_{\rho+} \rho)(k_- + k_z)}{2k_-}a_{mk_z -}
\right]
```

to get a system of linear equations to determine $`a_{mk_z\pm}`$. The choice of these
integrals is somewhat arbitrary, e.g., for the second integral one could use equally
well the values for negative polarization. In the case of $`k_z=0`$ these integrals are
similar to the spherical case integral in the sense, that they separate
$`\boldsymbol M`$ and $`\boldsymbol N`$ modes.

# Closing remarks

The provided files can also be used as a starting point for any bi-isotropic calculation
and are not necessarily restricted to the computation of T-matrix coefficients.
