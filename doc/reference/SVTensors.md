# Spatially-varying tensors in [[buff-em]]

[[buff-em]] supports objects with arbitrary user-specified
spatially-varying frequency-dependent permittivity tensors.
These tensors are described by simple text files conventionally
given the file extension `.SVTensor.`

The `.SVTensor` file contains lines of the form 

+ `Eps =` *function of space and frequency*

(to define a spatially-varying but isotropic permittivity)

+ `EpsXX =` *function of space and frequency*
+ `EpsXY =` *function of space and frequency*
+ ...
+ `EpsZZ =` *function of space and frequency*

(to define the individual cartesian components of the
permittivity tensor).

If you leave any off-diagonal components unspecified,
they will be assumed to be zero. If you do specify an 
off-diagonal component function, you only need to specify
*either* the above-diagonal *or* the below-diagonal component,
e.g. `EpsXY` but not also `EpsYX`; the code will automatically
set $\epsilon_{yx}=\epsilon_{xy}$. 

For diagonal components, if you specify `EpsXX` but omit
specifications for `EpsYY` and `EpsZZ` then the code
will set $\epsilon_{yy}=\epsilon_{zz}=\epsilon_{xx}.$

For example, the constant permittivity tensor
$\left(\begin{array}{ccc} 2+3i & 0.01 & 0    \\
                          0.01 & 2+3i & 0.00 \\
                          0    &    0 & 2.3i \\
 \end{array}\right)$
may be described like this:

````bash
EPSXX=2+3i
EPSXY=0.1
````

## Functions and variables

The user-defined *function of space and frequency* in the
definition of permittivity components is a character string
that may refer to any of the following variables:

| `w`          | angular frequency in units of 3e14 rad/sec
| `x,y,z`      | cartesian coordinates
| `r,Theta,Phi`| spherical coordinates

### Referring to [[scuff-em]] material designations

In many cases your permittivity functions will want to
refer to
[<span class="SC">scuff-em</span> material designations][scuffEMMaterials]
describing frequency-dependent homogeneous isotropic materials.
You can do this by including the string `MP_MATNAME`
in the function definition, where `MATNAME` is the name
of the [[scuff-em]] material.
For example, here

Note that there is no need to refer to the frequency `w` 
here; the `GOLD` and `SIO2` permittivities will automatically
by evaluated at the correct frequency.

## Isotropic but spatially-varying permittivities

If your material is isotropic (permittivity tensor
proportional to the 3x3 identity matrix), the `.SVTensor`
file may consist of just a single line of the
form `Eps = `<i>my permittivity function</i>`.

````bash
 Eps = (5 + 4*i/w)*(1 + (x+5)/10)
````

For example, here's an isotropic material whose dielectric
constant varies linearly as a function of the *x*-coordinate,
from a value of $\epsilon=5+\frac{4i}{\omega}$ at
$x=-5$ to a value of $\epsilon=10+\frac{8i}{\omega}$
at $x=+5$:

````bash
 Eps = (5 + 4*i/w)*(1 + (x+5)/10)
````

## Non-isotropic permittivities

````bash
 Eps = 
````

## Referring to

In many cases your spatially-varying permittivity

[[scuff-em]]-style
The functions specified

may 
`.SVTensor` file may

## Material definitions in `.SVTensor` files

Prior to the 

# Material definitions in `.SVTensor` files

[scuffEMMaterials]:	http://homerreid.com/
