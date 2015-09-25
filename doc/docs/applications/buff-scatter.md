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

### Common options

[[buff-scatter]] recognizes the following subset of the
list of commonly accepted options to [[buff-em]] command-line codes;
these all have the same meaning as the
[corresponding options in <span class="SC">scuff-em</span>][scuffOptions].
 (In general, almost everything mentioned in the
[General reference for <span class="SC">scuff-em</span> command-line applications][scuffGeneralReference] 
pertains to [[buff-em]] command-line applications as well.)

````bash
--geometry 
--Omega
--OmegaFile
````

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

````bash
 --PFTFile   MyFile.PFT
````

Requests that data on the absorbed and scattered power,
force, and torque for all bodies in the geometry be
written to the file `MyFile.PFT`. See the file header
in the output file for details on how to interpret its
contents.


````bash
 --OPFTFile     MyFile.OPFT
 --JDEPFTFile   MyFile.JDEPFT
 --DSIPFTFile   MyFile.DSIPFT
````

These options are similar to `--PFTFile`, but they 
request PFT calculations by specific algorithms: 
the "overlap" (OPFT) method, the "J dot E" (JDEPFT) 
method, or the "displaced surface-integral" (DSIPFT)
method. Note that JDEPFT is the default, so 
the options `--PFTFile` and `--JDEPFTFile` are
actually synonymous (if you specify both, only
one output file will be produced).

````bash
 --DSIRadius 5.0
 --DSIPoints 302

 --DSIMesh   MyBoundingMesh.msh
````

These options control the behavior of DSIPFT calculations.
The surface integral may be evaluated in one of two
ways: **(a)** by Lebedev cubature over a bounding sphere
surrounding the object, or **(b)** using a user-supplied
bounding surface mesh. 

For case **(a)**, the bounding sphere radius is set
by `--DSIRadius`, while the number of cubature points
is set by `--DSIPoints`. Type `buff-scatter --help`
to see list of allowed values for `--DSIPoints.`

For case **(b)**, `MyBoundingMesh.msh` should be a
[[gmsh]] mesh file describing a closed bounding
surface discretized into triangles.

Lebedev cubature generally yields better accuracy for 
the same number of cubature points; however, if your 
geometry contains closely-spaced objects
which cannot be enclosed in a bounding sphere without
also encompassing other objects, you will need to 
use `--DSIMesh.`

<a name="Examples"></a>
# 2. <span class="SC">buff-scatter</span> examples

+ [Power, force, and torque on spheres and Janus particles irradiated by plane waves](../examples/JanusParticles/index.md)

[buffEMGeometries]:	             ../reference/Geometries.md
[SVTensorFiles]:   	             ../reference/SVTensors.md
[scuffScatter]:    	             http://homerreid.github.io/scuff-em-documentation/applications/scuff-scatter/scuff-scatter.md
[IncidentFields]:  	             http://homerreid.github.io/scuff-em-documentation/reference/IncidentFields.md
[CommonOptions]:                     http://homerreid.github.io/scuff-em-documentation/applications/GeneralReference#CommonOptions
[scuffGeneralReference]:             http://homerreid.github.io/scuff-em-documentation/applications/GeneralReference
[scuffTransformations]:              http://homerreid.github.io/scuff-em-documentation/reference/Transformations
[scuffSIO2Spheres]:                  http://homerreid.github.io/scuff-em-documentation/examples/SiO2Spheres/SiO2Spheres/
