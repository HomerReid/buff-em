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

## Differences in caching

+ [[buff-em]] uses a caching scheme similar to that
used by [[scuff-em]] for storing the most 
computationally expensive contributions to
VIE matrix elements. However, unlike in [[scuff-em]],
in [[buff-em]] this process is designed to be 
entirely transparent to the user; in particular, 
there are no cache-related command-line arguments 
to application codes. 

Instead, for each `.vmsh`
file, [[buff-em]] looks **automatically** in
the directory specified by the environment variable
`BUFF_CACHE_DIR` for a file with the same base
filename and extension `.cache.' If this file is
 not found (or exists but is somehow incorrect),
[[buff-em]] will *automatically* write this file
(again to the directory specified by 
`BUFF_CACHE_DIR`)
after accumulating the appropriate data.
If the environment variable `BUFF_CACHE_DIR` is not 
set, `.cache` files are written to the current
working directory.

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
