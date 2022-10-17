# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Clp],[devel],[https://github.com/coin-or/Clp/issues/new],[coin-or-clp],[https://github.com/coin-or/Clp])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Clp which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there.
AC_CONFIG_SRCDIR(src/ClpSimplex.cpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h
AC_COIN_CHECK_MATH_HDRS

# Check for intrinsics header for ClpPackedMatrix.cpp
AC_CHECK_HEADERS([immintrin.h intrin.h arm_neon.h fmaintrin.h],[break],[])

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_CHK_PKG(CoinUtils,[ClpLib OsiClpUnitTest])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHK_PKG(Osi,[OsiClpLib OsiClpUnitTest])
AC_COIN_CHK_PKG(OsiTests,[OsiClpUnitTest],[osi-unittests])
AC_COIN_CHK_PKG(Sample,,[coindatasample],[],dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],[],dataonly)
AC_COIN_CHK_PKG(Glpk,[ClpLib],[coinglpk])

#############################################################################
#                                    Aboca                                  #
#############################################################################

AC_CLP_ABOCA(1)

#############################################################################
#                        Third party linear solvers                         #
#############################################################################

# AMD and CHOLMOD from UFL. Glpk is an alternative source for amd. CHOLMOD
# doc'n claims it requires METIS, but apparently that's not necessary for the
# version distributed with Fedora.

# amd.h
AC_COIN_CHK_LIBHDR(AMD,[ClpLib],[-lamd],[-I/usr/include/suitesparse],[],
  [[double Control[AMD_CONTROL];
   amd_defaults(Control);]],
  [#include "amd.h"])

# CHOLMOD from UFL
# Requires AMD, but apparently not for a link check with dynamic libraries.
# cholmod.h
AC_COIN_CHK_LIBHDR(CHOLMOD,[ClpLib],[-lcholmod],[-I/usr/include/suitesparse],[],
  [cholmod_start((cholmod_common*)0);],
  [#include "cholmod.h"])

# MUMPS. Check first for ThirdParty/MUMPS, on the assumption that if the user
# bothered to build it, we should use it. If it's not present, try for a
# system installation. If we find it, define CLP_HAS_MUMPS for export to code
# using clp.
AC_COIN_CHK_PKG(MUMPS,[ClpLib],[coinmumps])
if test $coin_has_mumps = no ; then
  AC_COIN_CHK_LIBHDR(MUMPS,[ClpLib],[-ldmumps],[-I/usr/include/MUMPS],[],
    [dmumps_c((DMUMPS_STRUC_C*)0)],
    [#include "dmumps_c.h"])
fi

# WSMP (Watson Sparse Matrix Package from IBM T.J. Watson)
AC_ARG_WITH([wsmp],
            AS_HELP_STRING([--with-wsmp],[specify WSMP library]),
            [have_wsmp=yes; wsmp_lflags=$withval],
            [have_wsmp=no])

if test "$have_wsmp" = "yes"; then
  AC_COIN_TRY_LINK([wssmp],[$wsmp_lflags],[],
    [AC_COIN_DEFINENAMEMANGLING([WSMP],[$wssmp_namemangling])
     CLPLIB_LFLAGS="$wsmp_lflags $CLPLIB_LFLAGS"
     AC_DEFINE(CLP_HAS_WSMP,1,[Define to 1 if WSMP is available])
    ],
    [AC_MSG_ERROR([Symbol wssmp not found with WSMP flags $wsmp_lflags.])])
fi

AM_CONDITIONAL([COIN_HAS_WSMP],[test $have_wsmp = yes])

#############################################################################
#                             Other dependencies                            #
#############################################################################

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
AC_CONFIG_HEADERS([src/config.h src/config_clp.h])

if test $coin_has_osi = yes ; then
  AC_CONFIG_HEADERS([src/OsiClp/config_osiclp.h])
fi

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
