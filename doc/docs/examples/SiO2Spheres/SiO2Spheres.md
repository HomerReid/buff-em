# Thermal radiation, heat transfer, and non-equilibrium Casimir forces between silicon dioxide spheres

In this example, we use [[buff-neq]] to reproduce the results
of [this example from the <span class="SC">scuff-neq</span> documentation][scuffSIO2Spheres]:
we will compute **(1)** the power radiated by a single SiO2 sphere, 
and **(2)** the heat transfer and non-equilibrium Casimir force
between two SiO2 spheres. 

The files for this example may be found in the 
`share/buff-em/examples/SiO2Spheres` subdirectory
of your [[scuff-em]] installation.

--------------------------------------------------
## [[gmsh]] geometry file and volume mesh for a single sphere

The [[gmsh]] geometry file [`Sphere.geo`](Sphere.geo)
describes a sphere of radius 1 micron. (This is the
same file used in [this example](../JanusParticles/index.md); 
as noted there, it is not the same `.geo` file that
was used in the
[<span class="SC">scuff-neq</span> version of this calculation][scuffSIO2Spheres],
because for [[buff-em]] we need volume meshes instead of 
surface meshes. This file may be meshed to create
coarse and fine volume meshes as follows:

````bash
% gmsh -3 -clscale 1 Sphere.geo
% RenameMesh3D Sphere.msh
% gmsh -3 -clscale 0.75 Sphere.geo
% RenameMesh3D Sphere.msh
````

(Here [`RenameMesh3D`][RenameMesh3D] is a simple `bash` script
that uses 
[<span class="SC">buff-analyze</span>][buffAnalyze].
to count the number of interior faces in a volume mesh and rename 
the mesh file accordingly; note that it also changes
the file extension from `.msh` to `.vmsh`, which I find
convenient for distinguishing volume mesh files from
surface mesh files. 

This produces the files `Sphere_637.msh` and `Sphere_1727.msh,`
which you can visualize by opening in [[gmsh]]::

````bash
% gmsh Sphere_637.msh
````
![Sphere_637 mesh image][Sphere637Image]

````bash
% gmsh Sphere_1727.msh
````
![Sphere_1727 mesh image][Sphere1727Image]

Note: As you can see from the first image here,
by default [[gmsh]]

--------------------------------------------------
## [[buff-em]] geometry files

The [[buff-em]] geometry file
[`SiO2Sphere_637.scuffgeo`](SiO2Sphere_637.scuffgeo)
describes a single SiO2 sphere.

The [[buff-em]] geometry files
[`SiO2Spheres_637.scuffgeo`](SiO2Spheres_637.scuffgeo)
and
[`SiO2Spheres_1727.scuffgeo`](SiO2Spheres_1727.scuffgeo)
each describe the same configuration: two SiO2 spheres
separated by a center--center distance of 10 microns.
You can visualize this configuration by typing e.g.

````bash
% buff-analyze --geometry SiO2Spheres_1727.scuffgeo --WriteGMSHFiles
% gmsh SiO2Spheres_1727.pp
````

![SiO2Spheres_1727 mesh image](SiO2Spheres_1727.png)

--------------------------------------------------
## Spectral density of radiated power

As described in the 
[<span class="CodeName">scuff-neq</span> documentation][scuff-neq],
[[scuff-neq]] computes the total power radiated by
finite-temperature objects as an integral over angular frequencies
$\omega,$ in which the integrand involves a
temperature-dependent Bose-Einstein factor 
and a temperature-independent dimensionless flux $\Phi.$ 
To calculate this radiated-power flux at a given set
of frequencies, we say

````bash
 % scuff-neq --geometry SiO2Sphere_501.scuffgeo --OmegaFile --PRad
````

where [`OmegaFile`](OmegaFile) is a list of
angular frequencies. (Here `--PRad` says that we 
are interested in the radiated power).
This produces the file
``SiO2Sphere_501.SiFlux``, which looks something
like this:

````bash
# scuff-neq run on superhr2 (07/11/15::00:31:36)
# data file columns: 
# 1 transform tag
# 2 omega 
# 3 (sourceObject,destObject) 
# 4 PRad flux spectral density
DEFAULT 1.000000e-01 11 4.18911788e-06 
DEFAULT 1.300000e-01 11 1.38869207e-05 
DEFAULT 1.600000e-01 11 3.93335327e-05 
DEFAULT 1.900000e-01 11 1.05263974e-04 
````

As the file header says, the second column here
is the angular frequency 
in units of $\omega_0=3\cdot 10^{14}$ rad/sec
and the fourth column is the dimensionless power
flux. (The first column lists the 
[geometrical transformation][Transformations]; since 
we didn't specify the `--transfile` option to 
[[scuff-neq]], we have just a single geometric
configuration, labeled `DEFAULT`. The third 
column identifies the source and destination objects;
since this geometry only has a single object,
the source and destination object are both 
always object 1 and this column always reads
`11`.)

Here's a plot of the data:

![Power radiation from an SiO2Sphere](SiO2Sphere_PowerRadiation.png)

In this plot, the solid line is the prediction of 
the [Krueger formalism][KruegerPaper], as computed
by a [[julia]] code called [`KruegerFormulas.jl`](KruegerFormulas.jl).

The plot is produced by [[gnuplot]] using 
[this script](Plotter.gp).

--------------------------------------------------
## Spectral density of power transfer and non-equilibrium force

Here's a [bash script](RunScript) that runs [[scuff-neq]]
for both the coarsely-meshed and finely-meshed two-sphere
geometry to compute the fluxes of power transfer
and nonequilibrium force between the spheres. 
Running the script produces files `SiO2Spheres_501.SIFlux`
and `SiO2Spheres_1479.SIFlux.` Here are plots (produced
by the same [[gnuplot]] script referenced above)
of the heat-transfer flux from sphere 1 to sphere 2,
and the force fluxes from sphere 1 to sphere 2 and
from sphere 2 to sphere 2, compared to the Krueger
T-matrix results (again computed using the [[julia]]
code referenced above).

![Power transfer between SiO2 Spheres](SiO2Spheres_PowerTransfer.png)

![Force between SiO2 Spheres](SiO2Spheres_F12.png)

![Force between SiO2 Spheres](SiO2Spheres_F22.png)

--------------------------------------------------

[scuff-neq]:                          ../../applications/scuff-neq/scuff-neq.md
[Transformations]:                    ../../reference/Transformations
[KruegerPaper]:                       http://dx.doi.org/10.1103/PhysRevB.86.115423
[scuffSIO2Spheres]:                   http://homerreid.github.io/scuff-em-documentation/examples/SiO2Spheres/SiO2Spheres/
[RenameMesh3D]:                       ../JanusParticles/JanusParticles.md#RenameMesh3D
[buffAnalyze]:                        ../../applications/buff-analyze.md
[Sphere637Image]:                     ../JanusParticles/Sphere637.png
[Sphere1727Image]:                    ../JanusParticles/Sphere1727.png
