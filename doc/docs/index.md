<p align="center"><img align="center" src="img/buffEMLogo.png"></p>

<p align="center"><h1 align="center">
 <span class="SmallCaps">buff-em</span> documentation: Table of contents
</h1> 
</p>

## General reference

* [Top-level overview](reference/TopLevel.md)
* [Installation](reference/Installing.md)
* [Geometry files](reference/Geometries.md)
* [SVTensor files](reference/SVTensors.md)
* [Key differences between <span class="SC">buff-em</span> and <span class="SC">scuff-em</span>](reference/BUFFvsSCUFF.md)

<a name="Examples"></a>
## Tutorial examples

+ [Mie scattering](examples/MieScattering/index.md)
+ [Power, force, and torque on Janus particles irradiated by plane waves](examples/JanusParticles/index.md)
+ [Thermal radiation, heat transfer, and non-equilibrium Casimir forces between silicon dioxide spheres](examples/SiO2Spheres/index.md)
+ [Thermal radiation and self-propulsion of photon torpedoes](examples/PhotonTorpedoes/index.md)
+ [Thermal self-rotation of QED pinwheels](examples/Pinwheels/index.md)

## Command-line application reference

#### *Nanophotonics code*
- [buff-scatter][buff-scatter]    - general-purpose electromagnetic scattering
      
#### *Non-equilibrium Casimir/ heat-transfer code*
- [buff-neq][buff-neq]            - radiative heat transfer and non-equilibrium Casimir forces/torques

#### *Utility code*
- [buff-analyze][buff-analyze]    - diagnostic tool to print info on [[buff-em]] geometries

## API reference

* [libbuff][libbuff] - Accessing [[buff-em]] from C++ programs

## Technical memo

* [<span class="SC">buff-em</span> technical memo][memo]

[buffEMLogo]:         img/buffEMLogo.png
[buff-scatter]:       applications/buff-scatter.md
[buff-neq]:           applications/buff-neq.md
[buff-analyze]:       applications/buff-analyze.md
[libbuff]:            API/libbuff.md
[memo]:               tex/buff-em-tex.md
