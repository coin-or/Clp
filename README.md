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

Clp is available in [Debian](http://packages.debian.org/search?searchon=sourcenames&keywords=clp) and [Ubuntu](https://launchpad.net/ubuntu/+source/clp).


## Getting Started using CoinBrew

To build Clp from source, obtain the `coinbrew` script from
https://coin-or.github.io/coinbrew/
and run


    /path/to/coinbrew fetch --main-proj=Clp
    /path/to/coinbrew build --main-proj=Clp --test
    /path/to/coinbrew install --main-proj=Clp


The `coinbrew` script will fetch [these](Dependencies) additional projects.


## Getting Started without CoinBrew (Expert users)

 0. Install [these Dependencies](Dependencies)
 1. Obtain the source code, e.g., from https://github.com/coin-or/Clp
 2. Run `./configure -C` to generate makefiles
 3. Run `make` to build the CoinUtils library
 4. Run `make test` to build and run the CoinUtils unit test program
 5. Run `make install` to install library and header files.


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
 * Source code [examples](Clp/examples)
 * [User's Guide](https://coin-or.github.io/Clp) (from 2004)

Interfaces:
 * [Matlab Interface + Windows x86 & x64 Interface Binaries (OPTI Toolbox)](https://www.inverseproblem.co.nz/OPTI/)
 * [Julia interface](https://github.com/JuliaOpt/Clp.jl)
 * [R and CLP - a quick start](https://cran.r-project.org/web/packages/clpAPI/vignettes/clpAPI.pdf)
 * [Java and CLP - performs well](http://orinanobworld.blogspot.co.uk/2016/06/using-clp-with-java.html)
