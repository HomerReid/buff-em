<h1> Spatially-varying permittivity tensors in
     <span class="SC">buff-em</span>
</h1>


[[buff-em]] supports objects with arbitrary user-specified
spatially-varying frequency-dependent permittivity tensors.
These tensors are described by simple text files conventionally
given the file extension `.SVTensor.`

The `.SVTensor` file is thought of as describing a
3$\times$3 matrix-valued function
$\mathbf{Q}$ of frequency and space:

$$\mathbf{Q}(\omega, \mathbf x)=\left(\begin{array}{ccc}
 Q_{xx}(\omega, \mathbf x) & 
 Q_{xy}(\omega, \mathbf x) & 
 Q_{xz}(\omega, \mathbf x) \\
 Q_{yx}(\omega, \mathbf x) & 
 Q_{yy}(\omega, \mathbf x) & 
 Q_{yz}(\omega, \mathbf x) \\
 Q_{zx}(\omega, \mathbf x) & 
 Q_{zy}(\omega, \mathbf x) & 
 Q_{zz}(\omega, \mathbf x)
\end{array}\right)$$

[TOC]

# 1. Syntax of the `.SVTensor` file

The `.SVTensor` file contains lines of the form 

+ `Q =` *function of space and frequency*

to define a spatially-varying but isotropic permittivity
(proportional to the $3 \times 3$ identity matrix), or

+ `Qxx =` *function of space and frequency*
+ `Qxy =` *function of space and frequency*
+ ...
+ `Qzz =` *function of space and frequency*

to define the individual cartesian components of the
permittivity tensor.

Note: Because `.SVTensor` files are most commonly
used to define permittivities $\boldsymbol{\epsilon}$, 
you can alternatively use the syntax

+ `Eps =` *function of space and frequency*

or 

+ `EpsXX =` *function of space and frequency*
+ `EpsXY =` *function of space and frequency*
+ ... 
+ `EpsZZ =` *function of space and frequency*

If you leave any off-diagonal components unspecified,
they will be assumed to be zero. If you do specify an 
off-diagonal component function, you only need to specify
*either* the above-diagonal *or* the below-diagonal component,
e.g. `Qxy` but not also `Qyx`; the code will automatically
set $Q_{yx}=Q_{xy}$. If you do specify two separate functions
for the two off-diagonal components of $\mathbf{Q}$, the code
will symmetrize by setting both components equal to 
their average.

For diagonal components, if you specify `Qxx` but omit
specifications for `Qyy` and `Qzz` then the code
will set $Q_{yy}=Q_{zz}=Q_{xx}.$

## Functions and variables

The user-defined *function of space and frequency* in the
definition of permittivity components is a character string
that may refer to any of the following variables:

+ `w`           (the angular frequency in units of 3e14 rad/sec)
+ `x,y,z`       (cartesian coordinates of points in space)
+ `r,Theta,Phi` (spherical coordinates of points in space)

### Referring to [[scuff-em]] material designations

In many cases your permittivity functions will want to refer to
[<span class="SC">scuff-em</span> material designations][scuffEMMaterials]
describing frequency-dependent homogeneous isotropic materials.
You can do this by including the string `MP_MATNAME`
in the function definition, where `MATNAME` is the name
of the [[scuff-em]] material.
For example, here is a description of a material
tensor that varies continuously from 100% gold to 
100% silicon dioxide as the *z* coordinate runs
from 0 to 1:

````bash
 Eps = MP_SIO2 * z + MP_GOLD*(1-z)
````

Note that there is no need to refer to the frequency `w` 
here; the `GOLD` and `SIO2` permittivities will automatically
by evaluated at the correct frequency.

The definitions for the [[scuff-em]] materials 
(`SIO2` and `GOLD` in this case) may appear in
`MATERIAL...ENDMATERIAL` sections within the `.SVTensor`
file; alternatively, they may be defined in the
global database file `${HOME}/.matprop.dat`
or the local data file `matprop.dat` in the current
working directory.


# 2. Examples of `.SVTensor` files

## 1. An isotropic but spatially-varying permittivity

Here's an isotropic material whose dielectric
constant varies linearly as a function of the *z*-coordinate,
from a value of $\epsilon=5+\frac{4i}{\omega}$ at
$z=-5$ to a value of $\epsilon=10+\frac{8i}{\omega}$
at $z=+5$:

````bash
 Eps = (5 + 4*i/w)*(1 + (z+5)/10 )
````

## 2. A non-isotropic permittivity

This example describes the constant permittivity tensor
$\boldsymbol{\epsilon}=
 \left(\begin{array}{ccc} 2+3i & 0.01 & 0    \\
                          0.01 & 2+3i & 0.00 \\
                          0    &    0 & 4.5i \\
       \end{array}\right)$:

````bash
EpsXX=2+3i
EpsXY=0.1
EpsZZ=4+5i
````

<a name="CoatedSphere"></a>
## 3. A gold-coated SiO2 sphere

This example is used to describe a silicon dioxide sphere of 
radius 2 um coated with a layer of gold.
Note that `step` refers to the Heaviside step function.

````bash
 Eps = step(2-r)*MP_SIO2 + step(r-2)*MP_GOLD
````

## 4. A gold/SiO2 Janus particle

In this case the material is gold at points above the *xy* plane
and SiO2 at points below the *xy* plane.

````bash
 Eps = step(z)*MP_GOLD + step(-z)*MP_SIO2
````

[scuffEMMaterials]:	http://homerreid.github.io/scuff-em-documentation/reference/Materials/
