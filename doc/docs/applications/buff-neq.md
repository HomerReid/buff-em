<h1> Modeling non-equilibrium electromagnetic fluctuations with
     <span class="SC">buff-neq</span>
</h1>

[[buff-neq]] is an application code in the [[buff-em]] suite for 
studying non-equilibrium (NEQ) electromagnetic-fluctuation-induced 
phenomena--specifically, for computing *radiative heat-transfer rates* 
and *non-equilibrium Casimir forces and torques* for bodies of 
arbitrary shapes and arbitrary (linear, isotropic, piecewise 
homogeneous) frequency-dependent permittivity and permeability.
[[buff-cas3d]] implements the 
*fluctuating-volume current (FSC)* approach
to numerical modeling of non-equilibrium fluctuation phenomena.

Mechanically, working with [[buff-neq]] is similar in many ways to
working with the non-equilibrium Casimir code 
[scuff-neq][scuff-neq]{.SC} in the [[scuff-em]] code suite.
In particular,

+ As in [[scuff-neq]], you can request either **(a)** frequency-resolved
information on heat-transfer rates and NEQ Casimir forces (in which case 
you will specify a list of frequencies and will get back a list of 
frequency-specific values of energy and momentum fluxes) or 
**(b)** frequency-integrated information, in which case you will assign 
temperatures to each body in your geometry and [[buff-neq]] will 
numerically integrate the fluxes, weighted by appropriate Bose-Einstein 
factors, to obtain the total heat-transfer rate or NEQ Casimir force. 
(For more details, see 
[What <span class="SC">buff-neq</span> actually computes](#WhatItComputes).)

+ As in [[scuff-neq]], you can specify an optional list of 
[geometrical transformations](../../reference/Transformations.md) 
describing various displacements and rotations of the bodies 
in your geometry; in this case you will get back values of the 
frequency-resolved or frequency-integrated quantities for each 
transformation you specify.

For Casimir forces and torques, the quantities computed by 
[[buff-neq]] are only the *non-equilibrium* contributions to 
the total force and torque---that
is, the contributions arising from the temperature 
differences between individual bodies and the surrounding environment.
To get the total force, these must be added to the *equilibrium*
contributions, which are the Casimir forces and torques for the 
case in which all bodies are at the temperature of the environment.
(For heat-transfer rates there is of course no equilibrium 
contribution, as there is no net power transfer between 
bodies at thermal equilibrium.)

[TOC]

# 1. What <span class="SC">buff-neq</span> actually computes

[[buff-neq]] implements the
FVC approach to non-equilibrium fluctuation phenomena.
This is 
an algorithm for computing the thermal averages of power,
force, and torque (PFT) quantities in geometries consisting
of compact material bodies at various temperatures embedded
in an finite-temperature or zero-temperature environment.
[[buff-neq]] can compute both
spatially-*resolved* and spatially-*integrated* PFT quantities.
(Examples of spatially-resolved quantities include components
of the average Poynting flux or Maxwell stress tensor at individual
points in space. Examples of spatially-integrated quantities 
include the total power absorbed by, or the total force or torque 
on, a compact homogeneous body. Spatially-integrated quantities 
are generally obtained by integrating spatially-resolved quantities
over closed bounding surfaces, although this is not necessarily
the way they are computed by [[buff-neq]].)

<a name="CommandLineOptions"></a>
# 2. <span class="SC">buff-neq</span> command-line options

### Common options

[[buff-neq]] recognizes the following subset of the 
[list of commonly accepted options to <span class="SC">buff-em</span> command-line codes][CommonOptions].

 
  ````
--geometry
--TransFile
--Omega
--OmegaFile
--OmegaQuadrature
--OmegaMin
--AbsTol
--RelTol
--FileBase
  ````
{.toc}

### Options requesting output quantities

  ````
--PAbs  
--PRad  
--XForce  
--YForce  
--ZForce
--XTorque  
--YTorque
--ZTorque
  ````
{.toc}

Specifies the quantities in which you are interested:
absorbed power (`--PAbs`), radiated power (`--PRad`),
Cartesian force components, or Cartesian torque components.
You may specify none, all, or any subset of these options,
but each option you specify will generally increase
the computation time (you can scrutinize the
[`.log` file](#LogFile) to see how *much* additional time each
extra output quantity takes to compute).

### Options specifying object temperatures

  ````
--Temperature     UpperSphere 300
--Temperature     LowerSphere 100

--TEnvironment    100

--TemperatureFile MyTemperatureFile.dat
  ````
{.toc}

> The first two options here set the temperatures
> of the objects labeled `UpperSphere` and
> `LowerSphere` in the `.buffgeo` file. 
> **Temperature specifications are interpreted in 
> units of Kelvin**, so `300` corresponds to 
> room temperature.
>
> The third option here sets the temperature of
> the environment in which the objects are embedded.
>
> Note that the temperatures of all objects, and of
> the environment, are zero by default. This means that,
> if you request a full frequency-integrated calculation
> (which you do by omitting the `--omega` or `--omegaFile`
> option) and you do not specify any `--temperature` 
> options, the code will chug for a while (computing 
> temperature-independent fluxes at various frequencies)
> before reporting strictly zero values for all
> quantities! This is probably not what you want.

### Options controlling the computation of power, force, and torque

[[buff-neq]] implements several algorithms for computing
powers, forces, and torques (PFTs). You may request the computation
of PFTs via more than one method.

  ````
 --JDEPFT
  ````

Requests that PFTs be computed using the $\mathbf{J} \cdot \mathbf{E}$ method.

  ````
 --MomentPFT
  ````

Requests that PFTs be computed using the moment method.

  ````
--DSIPoints  302
--DSIMesh BoundingMesh.msh

--DSIPoints2 590

--DSIRadius  5.0

  ````
{.toc}

These options request PFT calculations via the displaced-surface-integral
(DSI) method.

To request a DSIPFT calculation using `N` cubature points
a over a bounding sphere of radius `R,` say
`--DSIPoints N --DSIRadius R`. To see the allowed values
of `N`, type `buff-neq --help.`

Alternatively, you can say `--DSIMesh BoundingMesh.msh`
to request a DSIPFT calculation using 
a one-point cubature over each triangular panel of 
`BoundingMesh.msh`.

To get a "second opinion," you can request a second 
DSIPFT calculation by saying `--DSIPoints2 N`.
 This will perform a DSIPFT calculation over
bounding sphere with `N` points (with the sphere radius
set by `--DSIRadius`). This calculation will be 
performed in addition to the first DSIPFT calculation
you requested by specifying `--DSIPoints` or 
`--DSIMesh.`

--------------------------------------------------

<a name="Examples"></a>
# 3. Examples of calculations using <span class="SC">buff-neq</span>

[buff-cas3D]:    ../buff-cas3D/buff-cas3D.md
[EPFile]:        ../../applications/GeneralReference.md#EvaluationPoints
[LogFiles]:      ../GeneralReference.md#LogFiles
[SiO2Sphere]:    ../../examples/SiO2Spheres/SiO2Spheres.md
[SiO2Spheres]:   ../../examples/SiO2Spheres/SiO2Spheres.md
[TipSubstrate]:  ../../examples/TipSubstrate/TipSubstrate.md
[CommonOptions]: ../GeneralReference.md#CommonOptions
