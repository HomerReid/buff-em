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

BUFF-EM is complementary to
[SCUFF-EM](http://homerreid.github.io/scuff-em-documentation),
which implements the
*surface*-integral-equation method of computational electromagnetism.

Compared to SCUFF-EM, BUFF-EM has the advantage of allowing
a wider class of electromagnetic materials, including

+ anisotropic materials (non-diagonal permittivity or permeability tensors)

+ continuously spatially varying materials (permittivity/permeability are continuous functions of **x**)

Disadvantages of BUFF-EM compared to SCUFF-EM include

+ BUFF-EM is less computationally efficient for geometries that can be handled by both codes

+ BUFF-EM does not support Bloch-periodic boundary conditions

+ The BUFF-EM suite of application codes is not as full-featured as
  the SCUFF-EM suite (although the core library offers equivalent functionality).

For documentation and further information on BUFF-EM visit the 
[BUFF-EM documentation homepage:](http://homerreid.github.io/buff-em-documentation)
