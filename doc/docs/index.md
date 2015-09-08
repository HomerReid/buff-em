<p align="center"><img align="center" src="img/buffEMLogo.png"></p>

<p align="center"><h1 align="center">
 <span class="SmallCaps">buff-em</span> documentation: Table of contents
</h1> 
</p>

## General reference

* [Top-level overview](reference/TopLevel.md)
* [Installation](reference/Installing.md)
* [Geometry files](reference/Geometries.md)
* [SVTensor files](reference/Materials.md)
* [Key differences between <span class="SC">buff-em</sc> and <span class="SC">scuff-em</sc>](reference/BUFFvsSCUFF.md)

## Command-line application reference

* *Nanophotonics code*
    * [buff-scatter][buff-scatter]    - general-purpose electromagnetic scattering
      
* *Non-equilibrium Casimir/ heat-transfer code*
    - [buff-neq][buff-neq]            - radiative heat transfer and non-equilibrium Casimir forces/torques
  
* *Utility codes*
    - [buff-analyze][buff-analyze]    - diagnostic tool to print info on [[buff-em]] geometries

## API reference

* [libbuff][libbuff] - Accessing [[buff-em]] from C++ programs

## Technical memo

* [<span class="SC">buff-em</span>technical memo][memo]

[buffEMLogo]:         img/buffEMLogo.png
[buff-scatter]:       applications/buff-scatter/buff-scatter.md
[buff-neq]:           applications/buff-neq/buff-neq.md
[buff-analyze]:       applications/buff-analyze/buff-analyze.md
[libbuff]:            API/libbuff.md
[memo]:               tex/buff-em-tex.md
