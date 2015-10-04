# Top-level overview of [[buff-em]]

[[buff-em]] is a free, open-source software
implementation of the frequency-domain 
volume-integral-equation (VIE) method of 
classical electromagnetic
scattering using [SWG basis functions][SWGPaper].

[[buff-em]] is similar in many ways to 
[<span class="SC">scuff-em</span>][scuffEM],
which solves similar problems using the
alternative surface-integral-equation (SIE)
formalism. Some key differences between these
codes are

+ [[buff-em]] can handle bodies with inhomogeneous
and/or anisotropic dielectric permittivity. In
contrast, [[scuff-em]] can only handle homogeneous
isotropic materials.

+ [[buff-em]] requires scatterers to be described
by volume (tetrahedral) meshes instead of 
the surface (triangle) meshes used by [[scuff-em]].

+ [[buff-em]] does not aspire to solve the same
breadth of problems that [[scuff-em]] does. 
More specifically, [[buff-em]] is primarily
intended for just two classes of problem: 
**(a)** classical scattering of known incident
fields from compact bodies (implemented by the
command-line code 
[<span class="SC">buff-scatter</span>][buffScatter], 
and 
**(b)** non-equilibrium fluctuational
electrodynamics in the fluctuating-volume-current 
approach---including radiative heat transfer, thermal 
self propulsion/rotation, and non-equilibrium Casimir 
forces---for one or more compact bodies (implemented
by the command-line code
[<span class="SC">buff-neq</span>][buffNEQ]).

(For more on similarities and differences between
the two codes, see the document
[Key differences between <span class="SC">buff-em</span>
and <span class="SC">scuff-em</span>][BUFFvsSCUFF]).

Like [[scuff-em]], [[buff-em]] consists of a 
[core library](../API/libbuff.md),
implementing the basic VIE functionality, plus
the two specialized [application modules](#AvailableApplications) 
mentioned above for scattering and non-equilibrium
fluctuations.

[[buff-em]] stands for **BU**lk **F**ield **F**ormulation of 
**E**lectro**M**agnetism. This is a reference to the underlying solution 
methodology used by [[buff-em]] and other VIE solvers, in which
the primary goal of the solver is to compute the volume electric 
current distribution throughout a compact body (which is
locally proportional to the local bulk electric field,
whereupon the name).

Like [[scuff-em]], the entire [[buff-em]] suite is free software 
distributed under the [GNU GPL][GNUGPL]. The source code for
[[buff-em]] may be downloaded from the 
[<span class="SC">buff-em</span> GitHub page][GitHub]. 
**The GitHub page is also the right place for questions, 
bug reports, feature requests, and other discussion of [[buff-em]].**

Note that [[buff-em]] requires [[scuff-em]], so you will need a
working [[scuff-em]] installation on your system before you can 
start using [[buff-em]].

## Interfaces to [[buff-em]]

As is true for [[scuff-em]], the core computational engine 
in [[buff-em]] may be accessed via multiple interfaces.

The *command-line interface* consists of specialized
[command-line applications](#AvailableApplications) for
running specific calculations in computational
physics. Using [[buff-em]] in this way requires only
that you learn some basic command-line options;
it should be possible to come quickly up to speed
by following these
[tutorial examples](../index.md#Examples).

The *application programming interface* consists of 
a [C++ API](../API/libbuff.md)
that allows access to internal [[buff-em]] data structures
and methods for maximal flexibility in implementing your
own custom-designed physics codes.

## Inputs to [[buff-em]] calculations

Typical inputs to [[scuff-em]] calculations include

+ A [geometry file](Geometries.md) describing the scattering geometry

+ For anisotropic or inhomogeneous bodies, an 
  [`.SVTensor` file](SVTensors.md) describing
  the spatially-varying components of the permittivity
  tensor. ("SVTensor" stands for "spatially-varying tensor.")

+ Specification of the frequencies at which you want to 
  perform calculations.

+ An optional 
  [list of geometric transformations][Transformations]
  to be applied to the geometry, with calculations generally repeated
  at each transformation. (The usage and syntax of transformations
  in [[buff-em]] is identical to that in [[scuff-em]].)

+ For scattering codes: a specification of the 
  [Incident fields][IncidentFields]. (Incident fields in 
  in [[buff-em]] are handled the same way as in [[scuff-em]].)

+ Specifications of the output quantities you wish to get back: 
  field components at individual points in space, power/force/torque
  information, Casimir quantities, heat-transfer rates, etc.

## Outputs from [[buff-em]] calculations

Typical outputs from [[buff-em]] calculations include

+ text-based data files reporting output quantities

+ Visualization files written in 
  [<span class="SC">GMSH</span>][GMSH] post-processing
  format.

<a name="AvailableApplications"></a>
## Command-line Applications

### Nanophotonics / electromagnetic scattering 

 + [<span class="SC">buff-scatter</span>][buffScatter]
> A general-purpose solver for problems involving the 
> scattering of known incident fields from one or more
> compact objects.
> Available outputs include: scattered and total fields
> at arbitrary points in space; absorbed and scattered 
> power; force and torque (radiation pressure); and induced 
> multipole moments.

### Fluctuation-induced interactions

 + [<span class="SC">buff-neq</span>][buffNEQ]
> An implementation of the fluctuating-volume-current
> approach to non-equilibrium fluctuation-induced
> interactions among compact objects.
> Available outputs include: frequency-resolved or 
> frequency-integrated rates of heat radiation or 
> radiative heat transfer; non-equilibrium Casimir 
> forces; self-propulsion and self-rotation of 
> isolated bodies.

##Citing [[buff-em]]

If you find [[buff-em]] useful for generating
results included in publications, please consider citing both 
**(a)** one of the papers discussing the implementation of
[[buff-em]], and
**(b)** the URL for the code. For example, if you are writing
in LaTeX, you might write something like this:

````tex
Numerical computations were performed using {\sc buff-em}, a free,
open-source software implementation of the 
volume-integral-equation method~\cite{BUFF1, BUFF2}.
````

Here the ``BUFF1`` and ``BUFF2``
references refer to the following ``.bibtex`` entries:

````tex
@ARTICLE{BUFF1,
author = {{Homer Reid}, M.~T. and {Johnson}, S.~G.},
title = "{Efficient Computation of Power, Force, and Torque in 
BEM Scattering Calculations}",
journal = {ArXiv e-prints},
archivePrefix = "arXiv",
eprint = {1307.2966},
primaryClass = "physics.comp-ph",
keywords = {Physics - Computational Physics, Physics - Classical Physics},
year = 2013,
month = jul,
}

@ARTICLE{BUFF2,
note="\texttt{https://github.com/HomerReid/buff-em}"
}
````

[GMSH]: http://www.geuz.org/gmsh
[GNUGPL]:                            http://en.wikipedia.org/wiki/GNU_General_Public_License
[GitHub]:                            https://github.com/HomerReid/buff-em/
[scuffEM]:                           http://homerreid.github.io/scuff-em-documentation
[SWGPaper]:                          http://dx.doi.org/10.1109/TAP.1984.1143193
[buffScatter]:                       ../applications/buff-scatter.md
[buffNEQ]:                           ../applications/buff-neq.md
[BUFFvsSCUFF]:                       ../reference/BUFFvsSCUFF.md
[Examples]:                          ../../examples/index.md
[Transformations]:                   http://homerreid.github.io/scuff-em-documentation/reference/Transformations
[IncidentFields]:                    http://homerreid.github.io/scuff-em-documentation/reference/IncidentFields
