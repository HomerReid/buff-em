<h1> Solving electromagnetic scattering problems with
     <span class="SC">buff-scatter</span>
</h1>

[[buff-scatter]] is a tool within the [[buff-em]] code suite
for solving classical scattering problems involving
user-specified incident fields impinging on a material
geometry.

To run a scattering calculation in [[buff-scatter]], you will

+ Create a [<span class="SC">buff-em</sc> geometry file][buffEMGeometries] describing your geometry. This will involve generating tetrahedral mesh representations of each object in your geometry and possibly writing a `.SVTensor` file to describe any objects with anisotropic and/or inhomogeneous permittivity.

+ Specify the incident field that will impinge on your objects:
  a plane wave, a gaussian beam, a point dipole radiator,
  or some combination.

+ Run [[buff-scatter]] with command-line options specifying 
  the geometry, the incident fields, the frequencies at which you
  wish to calculate, and the types of output you want to get back.

The types of output you can request from [[buff-scatter]] include

+ The components of the scattered and total electric and magnetic
fields at an arbitrary list of evaluation points you provide.

+ The absorbed power, scattered power, force, and torque for 
each object in the geometry.

+ Data on the electric and magnetic dipole moments induced
  by the incident field on each object.

+ Visualization files plotting the distribution of induced
current within each object.

[TOC]

<a name="CommandLineOptions"></a>
# 1. <span class="SC">buff-neq</span> command-line options

## Options defining the scattering problem

````
--geometry MyGeometry.buffgeo
````
{.toc}

Specifies the geometry input file.

````
--Omega      3.14
--OmegaFile  MyOmegaFile
````
{.toc}

Specifies the angular frequencies at which to
run calculations. (Angular frequencies are interpreted
in units of $c/1\,\mu\text{m}=3\cdot 10^{14} rad/sec.$)
The `--Omega` option may be used more than once 
to specify multiple frequencies. Alternatively,
the `--OmegaFile` option may be used to specify the
name of a file containing a list of frequencies (one per
line) at which to run calculations.h

## Options defining the incident field

The options for specifying incident fields in
[[buff-scatter]] are identical to those in [[scuff-scatter]],
as described in detail on the page
[Incident fields in <span class="SC">scuff-em</span>][IncidentFields].

we refer you to the
[<span class="SC">scuff-scatter</span> documentation][scuffScatter]
for details, and here quote only the
available options without commentary.

````
--pwDirection    nx ny nz
--pwPolarization Ex Ey Ez

--psStrength Px Py Pz
--psLocation xx yy zz

--gbDirection nx ny nz
--gbPolarization Ex Ey Ez
--gbCenter Cx Cy Cz
--gbWaist W
````
{.toc}

(As in [[scuff-scatter]], these options may occur multiple times 
to define superpositions of multiple types of incident field.)

## Options requesting scattered and total fields

````
 --EPFile MyEPFile
````
{.toc}

Specifies a list of evaluation points at which to
compute and report components of the scattered and total
fields. This option may be specified more than once to 
define multiple sets of field evaluation points. 

The `--EPFile` in [[buff-scatter]] is identical to 
that in [[scuff-scatter]]; for details, see the 
[<span class="SC">scuff-scatter</span> documentation][scuffScatter]

## Options requesting power, force, and torque (PFT) data

````
 --EPFile MyEPFile
````

````bash


<a name="Examples"></a>
# 2. <span class="SC">buff-neq</span> examples

buffEMGeometries:	       ../reference/Geometries.md
SVTensorFiles:   	       ../reference/SVTensors.md
[scuffScatter]:    	       http://homerreid.github.io/scuff-em-documentation/applications/scuff-scatter/scuff-scatter.md
[IncidentFields]:              ../../reference/IncidentFields.md
