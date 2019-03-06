# CLP

Clp (*C*oin-or *l*inear *p*rogramming) is an open-source linear programming solver written in C++.
It is primarily meant to be used as a callable library, but a basic, stand-alone executable version is also available.
It is designed to find solutions of mathematical optimization problems of the from

min c<sup>t</sup>x

such that:
  row<sub>lower_bound</sub> &le; Ax &le; row<sub>upper_bound</sub>

  column<sub>lower_bound</sub> &le; x &le; column<sub>upper_bound</sub>


Clp is written in C++ and is released as open source code under the [Eclipse Public License (EPL)](http://www.opensource.org/licenses/eclipse-1.0).
It is available from the [COIN-OR initiative](http://www.coin-or.org/).
The code is written primarily by John J. Forrest, now retired from IBM Research.
The project is currently managed by John Forrest, Lou Hafer, [Julian Hall](https://www.maths.ed.ac.uk/hall/), and Matthew Saltzman.

The Clp website is https://github.com/coin-or/Clp.

## Getting Started using CoinBrew

To build Clp from source, obtain the `coinbrew` script from
https://coin-or.github.io/coinbrew/
and run


    /path/to/coinbrew fetch --mainProj=Clp
    /path/to/coinbrew build --mainProj=Clp --test
    /path/to/coinbrew install --mainProj=Clp


The `coinbrew` script will fetch [these](Dependencies) additional projects.


## Getting Started without CoinBrew (Expert users)

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

 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [mailing list](http://list.coin-or.org/mailman/listinfo/clp)
 * [Report a bug](https://github.com/coin-or/Clp/issues/new)
 * [Doxygen-generated html documentation](http://www.coin-or.org/Doxygen/Clp)
