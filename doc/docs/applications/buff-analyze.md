<h1> Analyzing objects and geometries with <span class="SC">buff-analyze</span></h1>

The [[buff-em]] suite comes with a simple standalone utility named 
[[buff-analyze]] that you can use to gather some quick statistics on
meshed objects and scattering geometries described by mesh files and 
geometry files.

There are several situations in which this can be useful:

-   You want to know how much memory will be occupied by the VIE matrix 
    for a geometry described by a `.buffgeo` file.
-   Your `.buffgeo` file contains multiple `OBJECTs`, each 
    described by a separate volume mesh and possibly displaced a
    and/or rotated, and you want to visualize the full geometry to 
    make sure the file you wrote actually describes what you want.
-   You have created a `.trans` file describing a list of 
    [geometrical transformations][Transformations]
    to be applied to your geometry, and before running a full calculation 
    you want to do a quick sanity check by visualizing the geometry under
    each of your transformations to make sure they are what you intended.
-   Your geometry involves one or more 
    [spatially inhomogeneous permittivities](../reference/SVTensors.md),
    and you want a plot of permittivity vs. space to double-check that 
    the permittivity you got is the one you wanted.
-   You have created a new volume mesh and you want to do a 
    preliminary calculation to compute the 
    [<span class="SC">buff-em</sc> cache](.././reference/BUFFvsSCUFF.md#Caching)
    before you start any actual calculations.

[TOC]

-------------------------------------
# buff-analyze Command-Line Options

## *Options specifying the file to analyze*

````bash
    --geometry MyGeometry.buffgeo
````

Analyze a full geometry described by a
[buff-em geometry.](../reference/Geometries.md)

````bash
    --mesh     MyObject.msh
    --meshFile MyObject.mesh
````

Analyze a single object described by a volume mesh. (The
two options are synonomous.)

## *Option specifying a list of geometrical transformations*

````bash
    --TransFile MyTransFile.trans
````

Specify a list of 
[geometrical transformations][Transformations]
to be applied to a geometry. This is useful for **(a)** checking
that your transformation file can be properly parsed by 
[[buff-em]], and **(b)** producing a visualization output file to 
confirm that the transformations you got are the ones you wanted.

## *Options controlling the generation of visualization files*

````bash
    --WriteGMSHFiles 
````

Append visualization data to [[gmsh]] visualization files that 
provides information on how the geometry is represented internally 
within [[buff-em]]. 

````bash
    --PlotPermittivity
````

Produces a visualization file showing the spatial variation
of the dielectric permittivity in your geometry. This is useful
if you specified any `.SVTensor` files and want to doublecheck
that your specification was interpreted correctly by [[buff-em]].

-------------------------------------
# [[buff-analyze]] console output

### Running [[buff-analyze]] on a geometry file

Running [[buff-analyze]] on a typical [[buff-em]] geometry file 
yields console output that looks like this:

````bash
% buff-analyze --geometry ThinPinwheel_N7_2190.buffgeo

*
* Geometry file MESH__ThinPinwheel_N7_2190__MAT__Gold.buffgeo: 1 objects 
*

**************************************************
* Object 0: 
**************************************************
Volume ThinPinwheel_N7_2190 (file ThinPinwheel_N7_2190.vmsh): 
 878 vertices 
 2190 tetrahedra 
 3568 interior faces 
 1624 exterior faces 

 Volume:          8.682680e-01 
 Surface area:    1.479837e+01 
 
 Moment of inertia: {+8.91e-01,+8.91e-01,+1.78e+00}
 
 Min/Max/Avg Quality factor: 1.376158e-01 / 9.928907e-01 / 7.479818e-01 

Thank you for your support.
````

-------------------------------------
# Running [[buff-analyze]] to precompute [[buff-em]] cache

[Transformations]:                   http://homerreid.github.io/scuff-em-documentation/reference/Transformations
