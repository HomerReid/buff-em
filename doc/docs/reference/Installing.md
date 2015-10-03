## Installing [[buff-em]]

## 1. Install <span class="CodeName">scuff-em</span>

[[buff-em]] requires [[scuff-em]], so if you don't have 
[[scuff-em]] installed you will need to do that first,
following the instructions [here][scuffEMInstallation].

## 2. External packages

[[buff-em]] doesn't require any external packages beyond
those required for [[scuff-em]].

Although not required to install, compile, or use
[[buff-em]],
[<span class="SC">gmsh</span>](http://geuz.org/gmsh)
is an extremely valuable open-source meshing and visualization
tool that is used throughout the
[[buff-em]] documentation.

On Debian/Ubuntu Linux systems, you can install [[gmsh]] by
doing a 

````bash
% sudo apt-get install gmsh
````

> Note: In some cases it seems the ``gmsh`` package conflicts with 
> the ``libhdf5-serial-dev`` package. In this case, just 
> remove ``gmsh`` from the above ``apt-get`` statement; you can 
> install it by hand following the instructions
> on the [GMSH website](http://geuz.org/gmsh).
> (Note that [[gmsh]], though very useful, 
> is not necessary to compile or run [[buff-em]].)

## 2. Cloning the GitHub repository and building the code

[[buff-em]] is hosted on [GitHub][GitHub].
To fetch and install the latest version of the 
code, execute the following steps. (Replace the string
``/path/to/buff-em-installation-directory``
with your desired installation directory.)

````bash
% git clone https://homerreid@github.com/HomerReid/buff-em.git
% cd buff-em
% sh autogen.sh --prefix=/path/to/buff-em-installation-directory
% make install
````

If this succeeds, the executable versions of the application
programs (such as ``buff-scatter``, ``buff-neq``, etc.) will be 
installed in the directory ``PREFIX/bin/`` 
and the demonstration examples for the various application programs 
will be in ``PREFIX/share/buff-em/examples``
(where ``PREFIX`` is the directory you specified using the 
``--prefix`` option above).

If you have trouble installing [[buff-em]],
please file an issue on the 
[<span class="SC">buff-em</span> GitHub page](GitHub).

### Build options

You may specify options to the ``autogen.sh``
(or ``configure``) command to guide the compilation process. 
For a full list of available options,
type ``configure --help.`` 

In some cases you may need to tweak certain environment 
variables to achieve maximal [[openmp]] performance
when running (not building) [[buff-em]].
For example, on my workstation (which has 8 CPU cores),
in order to get [[openmp]] codes
to use all 8 cores I need to set the following environment
variable:

````bash
% export GOMP_CPU_AFFINITY=0-7
````

[scuffEMInstallation]:               http://homerreid.github.io/scuff-em-documentation/reference/Installing
[GitHub]:                            http://github.com/HomerReid/buff-em
