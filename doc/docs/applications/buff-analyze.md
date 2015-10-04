<h1> Analyzing objects and geometries with <span class="SC">buff-analyze</span></h1>

The [[buff-em]] suite comes with a simple standalone utility named 
[[buff-analyze]] that you can use to gather some quick statistics on
meshed objects and scattering geometries described by mesh files and 
geometry files.

There are several situations in which this can be useful:

-   You want to know how much memory will be occupied by the VIE matrix 
    for a geometry described by a `.buffgeo` file.
-   Your `.buffgeo` file contains multiple `OBJECTs`, each 
    described by separate volume meshes and possibly displaced
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
    [<span class="SC">buff-em</span> cache](.././reference/BUFFvsSCUFF.md#Caching)
    before you start any actual calculations.

[TOC]

-------------------------------------
# [[buff-analyze]] Command-Line Options

## *Options specifying the file to analyze*

  ````
    --geometry MyGeometry.buffgeo
  ````
{.toc}

Analyze a full geometry described by a
[buff-em geometry.](../reference/Geometries.md)

  ````
--mesh     MyObject.vmsh 
  ````
{.toc}

  ````
--meshFile MyObject.vmsh
  ````
{.toc}

Analyze a single object described by a volume mesh. (The
two options are synonymous.)

## *Option specifying a list of geometrical transformations*

  ````
    --TransFile MyTransFile.trans
  ````
{.toc}

Specify a list of 
[geometrical transformations][Transformations]
to be applied to a geometry. This is useful for **(a)** checking
that your transformation file can be properly parsed by 
[[buff-em]], and **(b)** producing a visualization output file to 
confirm that the transformations you got are the ones you wanted.

## *Options controlling the generation of visualization files*

  ````
    --WriteGMSHFiles 
  ````
{.toc}

Append visualization data to [[gmsh]] visualization files that 
provides information on how the geometry is represented internally 
within [[buff-em]]. 

  ````
    --PlotPermittivity
  ````
{.toc}

Produces a visualization file showing the spatial variation
of the dielectric permittivity in your geometry. This is useful
if you specified any `.SVTensor` files and want to doublecheck
that your specification was interpreted correctly by [[buff-em]].

## *Options controlling the generation of cache files*

  ````
 --WriteCache
  ````
{.toc}

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
# Running [[buff-analyze]] to precompute the [[buff-em]] cache

As discussed [here](.././reference/BUFFvsSCUFF.md#Caching),
the first calculation performed by [[buff-em]] on a given
volume mesh is considerably slower than all subsequent 
calculations, because on this first calculation the code 
computes and certain frequency- and material-independent 
quantities needed to assemble the system matrix, which 
it then writes to a binary data file named `Mesh.cache`
(where `Mesh.vmsh` is the name of the mesh file as 
specified using the
[`MESHFILE` keyword in a `.buffgeo` file][buffGeometries]).
In some case it is convenient to generate the `.cache`
file before running any actual scattering or 
non-equilibrium Casimir calculations. (For example,
if you are planning to launch a bunch of parallel 
[[buff-em]] jobs that will refer to the same `.vmsh` 
files---even if the corresponding objects have 
different material properties!---you will want to 
precompute the `.cache` file so that the jobs don't
all repeat the same calculation when they start up).

You can use [[buff-analyze]] to precompute the cache
file for a given volume mesh like this:

````bash
 % buff-analyze --mesh MyObject.vmsh --WriteCache
````

This will produce a file named `MyObject.cache`.
If the environment variable `BUFF_CACHE_DIR` is 
set, then this file will be written to the
directory it specifies; otherwise, the file will
be written to the current working directory.

Then all subsequent calculations that
refer to `MyObject.vmsh` (including calculations
performed by [[buff-scatter]], [[buff-neq]], or
API codes) will automatically read in the
`MyObject.cache` file (assuming they can find
it either in their current working directory
or in the directory specified by `BUFF_CACHE_DIR`),
thus bypassing the costly cache-computation
step and significantly accelerating calculations.

[Transformations]:                   http://homerreid.github.io/scuff-em-documentation/reference/Transformations
[buffGeometries]:                    ../reference/Geometries.md
