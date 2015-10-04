[![Build Status](https://travis-ci.org/HomerReid/buff-em.svg?branch=master)](https://travis-ci.org/HomerReid/scuff-em)

BUFF-EM
========

A free, open-source software implementation of the 
volume-integral-equation method of computational electromagnetism
using SWG basis functions.
Includes a core library with a C++ API and command-line application codes
for 
electromagnetic scattering and non-equilibrium fluctuation-induced 
phenomena (Casimir forces and radiative heat transfer).

BUFF-EM stands for

 **BU**lk **F**ield **F**ormulation of **E**lectro**M**agnetism

For documentation and further information on BUFF-EM visit the 
[BUFF-EM documentation homepage:](http://homerreid.github.io/buff-em-documentation)

http://homerreid.github.io/buff-em-documentation

## General reference

* [Top-level overview](TopLevel)
* [Installation](Installing)
* [Geometry files](Geometries)
* [SVTensor files](SVTensors)
* [Key differences between BUFF-EM and SCUFF-EM](BUFFvsSCUFF)

<a name="Examples"></a>
## Tutorial examples

+ [Mie scattering](MieScattering)
+ [Power, force, and torque on Janus particles irradiated by plane waves](JanusParticles)
+ [Thermal radiation, heat transfer, and non-equilibrium Casimir forces between silicon dioxide spheres](SiO2Spheres)
+ [Thermal radiation and self-propulsion of photon torpedoes](PhotonTorpedoes)
+ [Thermal self-rotation of QED pinwheels](QEDPinwheels)

## Command-line application reference

#### *Nanophotonics code*
- [BUFF-SCATTER][buff-scatter]    - general-purpose electromagnetic scattering
      
#### *Non-equilibrium Casimir/ heat-transfer code*
- [BUFF-NEQ][buff-neq]            - radiative heat transfer and non-equilibrium Casimir forces/torques

#### *Utility code*
- [BUFF-ANALYZE][buff-analyze]    - diagnostic tool to print info on [[buff-em]] geometries

## API reference

* [LIBBUFF][libbuff] - Accessing [[buff-em]] from C++ programs

## Technical memo

* [BUFF-EM technical memo][memo]

[buff-scatter]:       http://homerreid.github.io/buff-em-documentation/applications/buff-scatter
[buff-neq]:           http://homerreid.github.io/buff-em-documentation/applications/buff-neq
[buff-analyze]:       http://homerreid.github.io/buff-em-documentation/applications/buff-analyze
[libbuff]:            http://homerreid.github.io/buff-em-documentation/API/libbuff
[memo]:               http://homerreid.github.io/buff-em-documentation/tex/buff-em-tex
[TopLevel]:           http://homerreid.github.io/buff-em-documentation/reference/TopLevel
[Installing]:         http://homerreid.github.io/buff-em-documentation/reference/Installing
[Geometries]:         http://homerreid.github.io/buff-em-documentation/reference/Geometries
[SVTensors]:          http://homerreid.github.io/buff-em-documentation/reference/SVTensors
[BUFFvsSCUFF]:        http://homerreid.github.io/buff-em-documentation/reference/BUFFvsSCUFF
[MieScattering]:      http://homerreid.github.io/buff-em-documentation/examples/MieScattering/MieScattering
[JanusParticles]:     http://homerreid.github.io/buff-em-documentation/examples/JanusParticles/JanusParticles
[SiO2Spheres]:        http://homerreid.github.io/buff-em-documentation/examples/SiO2Spheres/SiO2Spheres
[PhotonTorpedoes]:    http://homerreid.github.io/buff-em-documentation/examples/PhotonTorpedoes/PhotonTorpedoes
[QEDPinwheels]:       http://homerreid.github.io/buff-em-documentation/examples/Pinwheels/QEDPinwheels
