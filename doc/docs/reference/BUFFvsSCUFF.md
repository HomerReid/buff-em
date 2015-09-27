# Key differences between [[buff-em]] and [[scuff-em]]

## Differences in capabilities

### What [[buff-em]] can do that [[scuff-em]] can't do

+ [[buff-em]] can handle objects with anisotropic and/or
continuously spatially varying dielectric permittivity.
The components of the 3x3 dielectric tensor may be 
any arbitrary user-specified functions of space and
frequency.

+ The volume-integral technique used by [[buff-em]] to
handle non-equilibrium Casimir forces and torques is
significantly faster than the surface-integral
technique used by [[scuff-em.]]

### What [[buff-em]] *can't* do that [[suuff-em]] can do

+ [[buff-em]] does not currently support extended objects 
or periodic boundary conditions; all bodies must be compact 
objects. This may change in a later version of the code.

+ [[buff-em]] does not currently support magnetic materials 
($\mu$\ne 1$). This may change in a later version of the
code.

+ [[buff-em]] does not support perfectly electrically
conducting (PEC) bodies or bodies with surface
conductivity. This **will not change** in later versions
of the code, because these idealizations are not
supported by the underlying formulation of 
electromagnetism (the volume-integral-equation
approach) implemented by the [[buff-em]] core library.

## Differences in input files

+ Geometry files in [[buff-em]] are conventionally
given the extension `.buffgeo` instead of `.scuffgeo`
as in [[scuff-em]].

+ Mesh files---which describe collections of tetrahedra,
not collections of triangles---are conventionally given the 
extension `.vmsh` instead of `.msh`. (Tetrahedral meshes
may be generated in `gmsh` using the command-line 
argument `-3` as opposed to `-2` for triangular meshes.
This will produce a `.msh` file, which will be automatically
renamed to a `.vmsh` file by this [[bash]] script:
(`RenameMesh3D`)[RenameMesh3D].

+ [[buff-em]] looks for mesh files in the current
  working directory and in the directory specified
  by the environment variables `BUFF_MESH_DIR.`

<a name="Caching">
## Differences in caching

+ The discretized integral-equation formalisms implemented
by [[buff-em]] and [[scuff-em]] are similar in one key
respect: the assembly of the system matrix required
to solve scattering problems involves the computation of 
large numbers of singular multidimensional integrals. To
accelerate this task, both [[buff-em]] and [[scuff-em]]
implement a *caching* scheme in which certain 
frequency-independent contributions to matrix elements
for a given structure
are stored in binary data files to allow them to be
reused in subsequent calculations on the same structure.

+ However, whereas caching is considered as a somewhat
optional acceleration feature in [[scuff-em]],
it is treated as *mandatory* in [[buff-em]], because
the speedup afforded by caching is much greater
in this case (essentially because the singular
integrals in question are now 6-dimensional instead 
of 4-dimensional and thus much more expensive to 
compute from scratch; the acceleration enabled
by caching thus makes a greater difference).

+ For this reason, the caching process in [[buff-em]]
is designed to be largely transparent to the user;
in particular, there are no cache-related command-line
arguments to application codes. Instead, for each meshed 
object (each `.vmsh` file) specified in a `.buffgeo` 
file, [[buff-em]] looks **automatically** for the cache 
file; if no file is found, [[buff-em]] **automatically** 
writes this file to disk after the first matrix assembly.

+ The cache file associated to a meshed object
described by a [[gmsh]] file named `Object.vmsh` is
always named `Object.cache.` Before the first
matrix assembly, [[buff-em]] looks for
this file in a couple of different places:

> + The current working directory
> + The directory specified by the environment variable
`BUFF_CACHE_DIR`

+ If the `.cache` file is not found in either location,
[[buff-em]] computes all integrals from scratch
when it first assembles the system matrix (that is,
when it handles the first user-specified frequency)
and then writes the file to disk as soon as that
assembly is complete. If the environment variable
`BUFF_CACHE_DIR` is set, [[buff-em]] writes the
`.cache` file to the directory it specifies; otherwise,
the file is written to the current working directory.

+ As in [[scuff-em]], the `.cache` file is 
*independent of frequency and material properties*,
so it only needs to be computed **once** for
a given `.vmsh` file, after which it may be 
reused for different computations on that object
at different frequencies and even different 
material properties (including 
[anisotropic or inhomogeneous materials][SVTensors]).

+ The code [<span class="SC">buff-analyze</span>][buffAnalyze]
offers the command-line option `--WriteGCache` to precompute
and write to disk the `.cache` files for a given `.vmsh` file. 
Thus, before you start any [[buff-em]] calculation,
you can say 

````bash
 % buff-analyze --MeshFile Object.vmsh --WriteGCache
````

This will create a file named `Object.cache` in the
current working directory (or in the directory
specified by `BUFF_CACHE_DIR` if it is set). Although
this will take some time to complete, the advantage 
is that your calculations will run quickly already 
on the first frequency.

## Differences in the underlying formalism

+ The computational paradigm employed by `scuff-em`
separates space into contiguous homogeneous regions
bounded by closed surfaces. The fields in any
region are determined solely from knowledge of the 
surface currents on the surfaces bounding that surface
(together with the fields of any bulk sources that 
exist within the region).

In `buff-em` there is no such separation. Instead,
all objects exist in a single ginormous all-encompassing
homogeneous region (the vacuum) and the fields at any 
point receive contributions from all objects and from
any incident-field sources that may be present.

[SVTensors]:                          SVTensors.md
[buffAnalyze]:                        ../../applications/buff-analyze.md
