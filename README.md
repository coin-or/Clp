# Clp

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Clp?sort=semver)](https://github.com/coin-or/Clp/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme).
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation script._

Clp (*C*oin-or *l*inear *p*rogramming) is an open-source linear programming solver.
It is primarily meant to be used as a callable library, but a basic, stand-alone executable version is also available.
It is designed to find solutions of mathematical optimization problems of the form

minimize   c'x
such that  lhs &le; Ax &le; rhs
and        lb &le; x &le; ub


CLP includes primal and dual Simplex solvers.
Both dual and primal algorithms can use matrix storage methods provided by the user (0-1 and network matrices are already supported in addition to the default sparse matrix).
The dual algorithm has Dantzig and Steepest edge row pivot choices; new ones may be provided by the user.
The same is true for the column pivot choice of the primal algorithm.
The primal can also use a non linear cost which should work for piecewise linear convex functions.
CLP also includes a barrier method for solving LPs.


Clp is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org)

The Clp website is https://github.com/coin-or/Clp.

## CITE

[![DOI](https://zenodo.org/badge/173496299.svg)](https://zenodo.org/badge/latestdoi/173496299)

## CURRENT BUILD STATUS

[![Build Status](https://travis-ci.org/coin-or/Clp.svg?branch=master)](https://travis-ci.org/coin-or/Clp)

[![Build status](https://ci.appveyor.com/api/projects/status/h3daf7woiig6n176/branch/master?svg=true)](https://ci.appveyor.com/project/tkralphs/clp-m0kud/branch/master)

## DOWNLOAD

Binaries for most platforms are available as part of [Clp](https://bintray.com/coin-or/download/Clp). 

 * *Linux*: On Debian/Ubuntu, Clp is available in the package `coinor-clp` and can be installed with apt. On Fedora, Clp is available in the package `coin-or-Clp`.
 * *Windows*: The easiest way to get Clp on Windows is to download from *[Bintray](https://bintray.com/coin-or/download/Clp)*.
 * *Mac OS X*: The easiest way to get Cbc on Mac OS X is through [Homebrew](https://brew.sh).
   * `brew tap coin-or-tools/coinor`
   * `brew install coin-or-tools/coinor/clp`

Due to license incompatibilities, pre-compiled binaries lack some functionality.
If binaries are not available for your platform for the latest version and you would like to request them to be built and posted, feel free to let us know on the mailing list.

*Source code* can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Clp from the
 [releases](https://github.com/coin-or/Clp/releases) page.
 * Cloning this repository from [Github](https://github.com/coin-or/Clp) or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html).

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew`

To build CoinUtils from source, obtain the `coinbrew` script, do
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Clp@master
./coinbrew build Clp
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Clp
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

### With Microsoft Visual Studio

For Microsoft Visual C++ users, there are project files for version 10
available in the `MSVisualStudio` directory. First, obtain the source code
using either a Windows git client or download a snapshot. In MSVC++ Version
10, open the solution file (this should be converted to whatever version of
MSVC+ you are using) and build the Clp project. The code should build out of
the box with default settings.

It is also possible to build Clp with the Visual Studio compiler from the
command line using the procedure for Unix-like environments, using the Msys2
shell or CYGWIN. This is the recommended and best-supported way of building
Clp in Windows from source. To do so, make sure the `cl` compiler is in your
path and add `--enable-msvc` to build command of `coinbrew`.  

## Quick start

Running clp gives you some hints.  It can do a unit test (`clp -unitTest`) and solve netlib 
problems (`-netlib` or `-netlibp` using primal).  It can also solve problems and set tolerances
etc.  Just do 
```
clp 
``` 
and then try `?` or setting various stuff.

```
clp filename                #read file, do presolve and dual algorithm
clp filename -primalsimplex #use primal instead
```
On Linux, clp can do file completion and line editing if it can find the 
history, readline, and termcap packages when building.

If you want to stress the code, you can set various stuff, e.g., dantzig pricing
and then go into netlib testing.  It is not guaranteed that it will solve all 
netlib instances if you get too creative.  For instance using presolve makes 
netlib solve faster - but pilot87 prefers a large infeasibility weight.  So
```
clp -presolve on -dualbound 1.0e10 -netlib
```
works well.

There are examples in [examples](examples).  To create an executable, build 
with `coinbrew` as above and then do
```
cd build/Cbc/master/examples
make DRIVER=minimum #build the driver minimum.cpp
```
or whichever driver you want.  A list is in [Makefile](Makefile.in).
Three useful samples are:

 * `minimum.cpp` This is the simplest possible program to read an mps file.

 * `defaults.cpp`.  This does not do much more, but it does it in much more 
complicated way by specifically setting defaults so it does give more
useful information.  It also prints a solution in a format similar to that
of MPSX.

 * `presolve.cpp`  This is a good driver for larger problems.

Other ones can get complicated so start simple and work your way up.

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxygen-docs` 

in the build directory. If Clp was built via `coinbrew`, then the build
directory will be `./build/Clp/master` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Clp/Doxygen).

## Project Links

Help:
 * [mailing list](http://list.coin-or.org/mailman/listinfo/clp)
 * [Report a bug](https://github.com/coin-or/Clp/issues/new)
 
Documentation:
 * [Doxygen-generated html documentation](http://www.coin-or.org/Doxygen/Clp)
 * Source code [examples](examples/)
 * [User's Guide](https://coin-or.github.io/Clp) (from 2004)

Interfaces:
 * [Matlab Interface + Windows x86 & x64 Interface Binaries (OPTI Toolbox)](https://www.inverseproblem.co.nz/OPTI/)
 * [Julia interface](https://github.com/JuliaOpt/Clp.jl)
 * [R and CLP - a quick start](https://cran.r-project.org/web/packages/clpAPI/vignettes/clpAPI.pdf)
 * [Java and CLP - performs well](http://orinanobworld.blogspot.co.uk/2016/06/using-clp-with-java.html)

