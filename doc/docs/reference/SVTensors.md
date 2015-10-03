<h1> Spatially-varying permittivity tensors in
     <span class="SC">buff-em</span>
</h1>


[[buff-em]] supports objects with arbitrary user-specified
spatially-varying frequency-dependent permittivity tensors.
These tensors are described by simple text files conventionally
given the file extension `.SVTensor.`

[TOC]

# 1. Syntax of the `.SVTensor` file

The `.SVTensor` file contains lines of the form 

+ `Eps =` *function of space and frequency*

to define a spatially-varying but isotropic permittivity, or

+ `EpsXX =` *function of space and frequency*
+ `EpsXY =` *function of space and frequency*
+ ...
+ `EpsZZ =` *function of space and frequency*

to define the individual cartesian components of the
permittivity tensor.

If you leave any off-diagonal components unspecified,
they will be assumed to be zero. If you do specify an 
off-diagonal component function, you only need to specify
*either* the above-diagonal *or* the below-diagonal component,
e.g. `EpsXY` but not also `EpsYX`; the code will automatically
set $\epsilon_{yx}=\epsilon_{xy}$. 

For diagonal components, if you specify `EpsXX` but omit
specifications for `EpsYY` and `EpsZZ` then the code
will set $\epsilon_{yy}=\epsilon_{zz}=\epsilon_{xx}.$

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
EPSXX=2+3i
EPSXY=0.1
EPSZZ=4+5i
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
