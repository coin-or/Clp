# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Clp],[trunk],[clp@list.coin-or.org],[coin-or-clp],
        [https://projects.coin-or.org/Clp])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Clp which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there.

AC_CONFIG_SRCDIR(src/ClpSimplex.cpp)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE(Clp)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.

AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.

AC_LANG_PUSH(C++)

# Sort out what we need to do to link to Fortran libraries. Check if there's a
# Fortran compiler and call the appropriate macro based on the result. If
# there's no Fortran compiler, it's necessary to specify a library and
# function to use when attempting a link. For Clp, libblas (-lblas) and
# DAXPY are appropriate choices.
# Yes, we should be able to leverage this and avoid the check for BLAS below,
# but I haven't quite figured it out yet.

AC_COIN_PROG_F77
if test "$F77" != "unavailable" ; then
  AC_COIN_F77_SETUP
else
  AC_COIN_F77_WRAPPERS(blas,DAXPY)
fi

# Initialize libtool

AC_COIN_PROG_LIBTOOL

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h

AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

# Initialise variables for ClpLib and OsiClpLib

AC_COIN_INIT_LIBVARS([ClpLib],[-L\${libdir} -lClp],[-I\${includedir}])
AC_COIN_INIT_LIBVARS([OsiClpLib],[-L\${libdir} -lOsiClp],[-I\${includedir}])

AC_COIN_CHK_PKG(CoinUtils,[ClpLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHK_PKG(Osi,[OsiClpLib])
AC_COIN_CHK_PKG(OsiTests,[OsiClpUnitTest],[osi-unittests])
AC_COIN_CHK_PKG(Sample,,[coindatasample],[],dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],[],dataonly)

#############################################################################
#                                    Aboca                                  #
#############################################################################

AC_CLP_ABOCA(1)

#############################################################################
#                        Third party linear solvers                         #
#############################################################################

# AMD from UFL
# [amd.h] [amd_defaults]

AC_COIN_CHK_LIB(AMD,[ClpLib],[-lamd],[-I/usr/include/suitesparse])

# CHOLMOD from UFL (requires AMD).
# [cholmod.h] [cholmod_start]

AC_COIN_CHK_LIB(CHOLMOD,[ClpLib],[-lcholmod],[-I/usr/include/suitesparse])

# Glpk also brings AMD

if test $coin_has_cholmod = no -a $coin_has_amd = no ; then
  AC_COIN_CHK_PKG(Glpk,[ClpLib],[coinglpk])
  if test $coin_has_glpk = yes ; then
    AC_MSG_NOTICE([using AMD from GLPK package])
    AC_DEFINE(COIN_HAS_AMD,[1],[Define to 1 if the AMD package is available])
  fi
else
  # for configure
  AM_CONDITIONAL(COIN_HAS_GLPK, [test 0 = 1])
fi

# MUMPS
AC_COIN_CHK_PKG(Mumps,[ClpLib],[coinmumps])

# WSMP (Watson Sparse Matrix Package from IBM T.J. Watson)

AC_COIN_CHK_LIB(WSMP,[ClpLib],[],[],[],no)

if test "$use_wsmp" = yes ; then
  # AC_COIN_CHK_PKG_BLAS(ClpLib)
  use_wsmp=
  AC_COIN_CHK_PKG_FLINK([use_wsmp],wssmp,[$wsmp_libs])
  if test -n "$use_wsmp" ; then
    wsmp_libs=$use_wsmp
    AC_DEFINE_UNQUOTED([COIN_HAS_WSMP],[1],
        [Define to 1 if the WSMP package is available])
    CLPLIB_LIBS="$use_wsmp $CLPLIB_LIBS"
  else
    AC_MSG_WARN([WSMP failed to link with $wsmp_libs])
  fi
else
  # to please configure
  # AM_CONDITIONAL(COIN_HAS_BLAS, [test 0 = 1])
  :
fi
AM_CONDITIONAL(COIN_HAS_WSMP, [test "$use_wsmp" = yes])

#############################################################################
#                             Other dependencies                            #
#############################################################################

# AMPL solver library

AC_COIN_CHK_PKG(ASL,[ClpLib],[coinasl])

AC_COIN_CHK_GNU_READLINE([ClpLib])

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_EXAMPLE_FILES([examples/hello.mps examples/input.130])
AC_COIN_EXAMPLE_FILES([examples/g.tiny examples/gparm.tiny])

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils Osi])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Temporary hacks to shut down these two. This needs to be fixed.

AM_CONDITIONAL(COIN_CXX_IS_CL,[test 0 = 1])

# Finalize libs and cflags.
AC_COIN_FINALIZE_FLAGS([ClpLib OsiClpLib OsiClpUnitTest])
AC_COIN_CHK_HERE([ClpLib],[OsiClpLib],[clp])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/OsiClp/Makefile
                 test/Makefile
                 clp.pc])

if test $coin_has_osi = yes ; then
  AC_CONFIG_FILES([osi-clp.pc:src/OsiClp/osi-clp.pc.in])
fi

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_clp.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE

