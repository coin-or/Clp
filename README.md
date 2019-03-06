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
The project is currently managed by John Forrest, Lou Hafer, [Julian Hall](https://www.maths.ed.ac.uk/hall/), and Matthew Saltzman.

The Clp website is https://github.com/coin-or/Clp.

Clp is available in [Debian](http://packages.debian.org/search?searchon=sourcenames&keywords=clp) and [Ubuntu](https://launchpad.net/ubuntu/+source/clp).


## Getting Started using CoinBrew

To build Clp from source, obtain the `coinbrew` script from
https://coin-or.github.io/coinbrew/
and run


    /path/to/coinbrew fetch --mainProj=Clp
    /path/to/coinbrew build --mainProj=Clp --test
    /path/to/coinbrew install --mainProj=Clp


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

Support:
 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [mailing list](http://list.coin-or.org/mailman/listinfo/clp)
 * [Report a bug](https://github.com/coin-or/Clp/issues/new)
 
Documentation:
 * [User's Guide](http://www.coin-or.org/Clp/userguide/index.html) ([single page format](http://www.coin-or.org/Clp/userguide/clpuserguide.html))
 * [Doxygen-generated html documentation](http://www.coin-or.org/Doxygen/Clp)
 * Source code [examples](Clp/examples)

Interfaces:
 * [Matlab Interface + Windows x86 & x64 Interface Binaries (OPTI Toolbox)](https://www.inverseproblem.co.nz/OPTI/)
 * [Julia interface](https://github.com/JuliaOpt/Clp.jl)
 * [R and CLP - a quick start](https://cran.r-project.org/web/packages/clpAPI/vignettes/clpAPI.pdf)
 * [Java and CLP - performs well](http://orinanobworld.blogspot.co.uk/2016/06/using-clp-with-java.html)


## FAQ (from 2004)

### The barrier method sounds interesting, what are some of the details?

The CLP barrier method solves convex QPs as well as LPs.
In general, a barrier method requires implementation of the algorithm, as well as a fast Cholesky factorization.
CLP provides the algorithm, and is expected to have a reasonable factorization implementation.
However, the sparse factorization requires a good ordering algorithm, which the user is expected to provide (perhaps a better factorization code as well).

### Which Cholesky factorizations codes are supported by CLP's barrier method?

The Cholesky interface is flexible enough so that a variety of Cholesky ordering and factorization codes can be used.
Interfaces are provided to each of the following:
 * Anshul Gupta's WSSMP parallel enabled ordering and factorization code
 * Sivan Toledo's TAUCS parallel enabled factorization code (the package includes third party ordering codes)
 * University of Florida's Approximate Minimum Degree (AMD) ordering code (the CLP native factorization code is used with this ordering code)
 * CLP native code: very weak ordering but competitive nonparallel factorization
 * Fast dense factorization 


### When will CLP have a good native ordering?
The best outcome would be to have an existing ordering code available as part of the COIN-OR distribution under the EPL.
However, if this is not possible, the native ordering will be made respectable.


### Is the barrier code as mature as the simplex code?
The simplex code has been exposed to user testing for a while and the principal author, John Forrest, knows more about simplex algorithms than interior point algorithms, so the answer is "no".
However, it performs well on test sets and seems to be more reliable than some commercially available codes (including OSL).


### Which algorithm should I use for quadratic programming and should I keep an eye open for any issues?
The interior point algorithm for quadratic programming is much more elegant and normally much faster than the quadratic simplex code.
Caution is suggested with the presolve as not all bugs have been found and squashed when a quadratic objective is used.
One may wish to switch off the crossover to a basic feasible solution as the simplex code can be slow.
The sequential linear code is useful as a "crash" to the simplex code; its convergence is poor but, say, 100 iterations could set up the problem well for the simplex code.
