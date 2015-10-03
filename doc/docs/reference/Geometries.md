<h1> Geometry descriptions in 
     <span class="SC">buff-em</span>
</h1>

Geometries in [[buff-em]] are described by simple text files 
that are conventionally given the file extension `.buffgeo`.

[TOC]

# 1. Syntax of the `.buffgeo` file

The `.buffgeo` file consists of one or more sections,
delineated by the keywords `OBJECT...ENDOBJECT.`
Each section defines a single compact object in a scattering
geometry and takes the form

````bash
OBJECT Label
	KEYWORD argument
	KEYWORD argument
	...
ENDOBJECT
````

where the various possible `KEYWORD argument` pairs are
detailed below.

+ `MESHFILE MyMeshFile.vmsh`

Specifies the 3D (tetrahedra) volume mesh defining the object.
The file should be present either in the current working 
directory or in the search path specified by the environment
variable `BUFF_MESH_PATH`.

3D volume meshes are generated from a [[gmsh]] geometry 
file by running `%gmsh -3 MyObject.geo`; this will produce a 
file named `MyObject.msh`, which I typically rename to 
`MyObject.vmsh` to indicate that it is a *volume* mesh
as opposed to the surface mesh files used by [[scuff-em]].
 
+ `MATERIAL GOLD`
+ `SVTENSOR MyMaterial.SVTensor`

Specifies the material properties of the object. 

The `MATERIAL` keyword is used for homogeneous isotropic
bodies. In this case, `GOLD` should be a
[<span class="SC">scuff-em</span> material designation][scuffEMMaterials].
This may be defined either in a database file (such as the
file `${HOME}/.matprop.dat` or `./matprop.dat`) or on the fly
in your `.buffgeo` file by including a `MATERIAL...ENDMATERIAL`
section before your object definition.

The alternative `SVTENSOR` keyword is used for inhomogeneous 
and/or anisotropic materials. (It stands for "spatially-varying
tensor.") In this case, `MyMaterial.SVTensor`
should be a 
[<span class="SC">buff-em</span> spatially-varying tensor file][SVTensors]
present either in the current working directory or in the search
path defined by the environment variable `BUFF_SVTENSOR_PATH`.


+ `DISPLACED 2.3 4.5 6.7`
+ `ROTATED 90 ABOUT 0 0 1`

These options specify geometrical transformations to be
performed on the object after it has been read in from the
`.vmsh` file.

# 2. Examples of <span class="SC">buff-em</span> geometries

## a. An SiO2 sphere coated with a layer of gold

````bash
OBJECT CoatedSphere
	MESHFILE CoatedSphere_4017.vmsh
	SVTENSOR CoatedSphere.SVTensor
ENDOBJECT
````

The file `CoatedSphere.SVTensor` is described [here](SVTensors.md#CoatedSphere).

## b. Two gold spheres, described by the same volume mesh but separated 3 length units in the *z* direction

````bash
MATERIAL GOLD
    wp = 1.37e16; 
    gamma = 5.32e13;
    Eps(w) = 1 - wp^2 / (w * (w + i*gamma));
ENDMATERIAL

OBJECT UpperSphere
	MESHFILE Sphere_677.vmsh
	MATERIAL Gold
ENDOBJECT

OBJECT LowerSphere
	MESHFILE Sphere_677.vmsh
	MATERIAL Gold
	DISPLACED 0 0 3
ENDOBJECT
````


[scuffEMMaterials]:	http://homerreid.github.io/scuff-em-documentation/reference/Materials/
[SVTensors]:       	SVTensors.md
