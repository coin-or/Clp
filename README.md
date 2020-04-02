# CLP

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


Clp is written in C++ and is released as open source code under the [Eclipse Public License (EPL)](http://www.opensource.org/licenses/eclipse-1.0).
It is available from the [COIN-OR initiative](http://www.coin-or.org/).
The code is written primarily by John J. Forrest, now retired from IBM Research.
The project is currently managed by John Forrest, Lou Hafer, [Julian Hall](https://www.maths.ed.ac.uk/school-of-mathematics/people?person=47), and Matthew Saltzman.

The Clp website is https://github.com/coin-or/Clp.

## CITE

[![DOI](https://zenodo.org/badge/173496299.svg)](https://zenodo.org/badge/latestdoi/173496299)

## CURRENT BUILD STATUS

[![Build Status](https://travis-ci.org/coin-or/Clp.svg?branch=master)](https://travis-ci.org/coin-or/Clp)

[![Build status](https://ci.appveyor.com/api/projects/status/h3daf7woiig6n176?svg=true)](https://ci.appveyor.com/project/tkralphs/clp-m0kud)

## DOWNLOAD

Binaries for most platforms are available for download from [Bintray](https://bintray.com/coin-or/download/Clp)

[ ![Download](https://api.bintray.com/packages/coin-or/download/Cbc/images/download.svg?version=1.17) ](https://bintray.com/coin-or/download/Clp/1.17/link)

 * *Linux*: On Debian/Ubuntu, Cbc is available in the package `coinor-clp` and can be installed with apt. On Fedora, Cbc is available in the package `coin-or-Clp`.
 * *Windows*: The easiest way to get Cbc on Windows is to download from *[Bintray](https://bintray.com/coin-or/download/Clp)*, although an old interactive installer for the [COIN-OR Optimization Suite](http://www.coin-or.org/download/binary/CoinAll) is also still available.
 * *Mac OS X*: The easiest way to get Cbc on Mac OS X is through [Homebrew](https://brew.sh).
   * `brew tap coin-or-tools/coinor`
   * `brew install clp`
 * AMPL also provides stand-alone [Clp executables](http://ampl.com/products/solvers/open-source/#clp) that can be used with (or without) AMPL.
 * The [GAMS](http://www.gams.com) distribution includes Clp.

Due to license incompatibilities, pre-compiled binaries lack some functionality.
If binaries are not available for your platform for the latest version and you would like to request them to be built and posted, feel free to let us know on the mailing list.

*Source code* can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Cbc from the [Clp source code download page](http://www.coin-or.org/download/source/Clp), or
 * Checking out the code from [Github](https://github.com/coin-or/Clp) or using the `coinbrew` script (recommended). 

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html) (this is a work in
progress).

## BUILDING from source

### Using CoinBrew

To build Clp from source, obtain the `coinbrew` script from
https://coin-or.github.io/coinbrew/
and run


    /path/to/coinbrew fetch Clp
    /path/to/coinbrew build Clp --prefix=/dir/to/install --test
    /path/to/coinbrew install Clp

The `coinbrew` script will fetch [these](.coin-or/Dependencies) additional projects.


### Without CoinBrew (Expert users)

 0. Install [these Dependencies](.coin-or/Dependencies)
 1. Obtain the source code, e.g., from https://github.com/coin-or/Clp
 2. Run `./configure -C` to generate makefiles
 3. Run `make` to build the CoinUtils library
 4. Run `make test` to build and run the CoinUtils unit test program
 5. Run `make install` to install library and header files.

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
path and add `--enable-msvc to build command of `coinbrew`.  

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

 `make doxydoc` 

in the build directory.
If Clp was build via `coinbrew`, then the build directory is `./build/Clp`.
The doxygen documentation main file is found at `./doxydoc/html/index.html` in the build directory.

If `Doxygen` is not available, you can use also use [this link](http://www.coin-or.org/Doxygen/Clp).

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
