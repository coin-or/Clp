# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Clp],[trunk],[clp@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/ClpSimplex.cpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Clp)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_HAS_MODULE(CoinUtils, [coinutils = trunk])
if test $coin_has_coinutils != yes ; then
  #if we have to use the fallback routines, then obviously it will not find CoinUtils dependencies
  #thus, we check whether ThirdParty/Blas and ThirdParty/Lapack have been setup and
  #assume that if present then they are dependencies of CoinUtils
  AC_COIN_HAS_MODULE_FALLBACK(Blas,   [coinblas],   [../ThirdParty/Blas],   [required=0 dodefine=0 doconditional=0])
  AC_COIN_HAS_MODULE_FALLBACK(Lapack, [coinlapack], [../ThirdParty/Lapack], [required=0 dodefine=0 doconditional=0])
fi
AC_COIN_HAS_MODULE_FALLBACK(CoinUtils, [coinutils], [../CoinUtils])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
#so if we had to do find Blas and Lapack separately, then add their linker flags to the CoinUtils flags
#actually the libs may already be included in libCoinUtils, but we may still need the fortran library flags
COINUTILS_LIBS="$COINUTILS_LIBS $BLAS_LIBS $LAPACK_LIBS"

AC_COIN_HAS_MODULE(Osi,      [osi],            [required=0], [../Osi])
AC_COIN_HAS_MODULE(OsiTests, [osi-unittests],  [required=0], [../Osi])
AC_COIN_HAS_MODULE(Sample,   [coindatasample], [required=0], [../Data/Sample])
AC_COIN_HAS_MODULE(Netlib,   [coindatanetlib], [required=0], [../Data/Netlib])

#############################################################################
#                        Third party linear solvers                         #
#############################################################################

# WSMP
AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],
                           [specify WSMP library]),
            [use_wsmp=$withval], [use_wsmp=no])

if test "$use_wsmp" != "no"; then
  LIBS="$use_wsmp $LIBS"
  ADDLIBS="$use_wsmp $ADDLIBS"

  defval=
  AC_MSG_CHECKING([whether user-supplied WSMP library \"$use_wsmp\" works (symbol wssmp_)])
  AC_TRY_LINK([],[wssmp_()],
              [AC_MSG_RESULT(yes)
               defval=2],
              [AC_MSG_RESULT(no)])
  if test "$defval" = ""; then
    AC_MSG_CHECKING([whether user-supplied WSMP library \"$use_wsmp\" works (symbol wssmp)])
    AC_TRY_LINK([],[wssmp()],
                [AC_MSG_RESULT(yes)
                 defval=3],
                [AC_MSG_RESULT(no)])
  fi
  if test "$defval" = ""; then
    AC_MSG_ERROR([WSMP library $use_wsmp does not seem to work])
  fi

  AC_DEFINE_UNQUOTED([WSSMP_BARRIER], [$defval],
                     [Define to 2 if wssmp_ is available, and to 3 if wssmp is available])
fi

#############################################################################
#                           CLP version number                              #
#############################################################################

AC_DEFINE_UNQUOTED(CLP_VERSION, "$PACKAGE_VERSION", "Version of Clp")

#############################################################################
#                             Other dependencies                            #
#############################################################################

AC_COIN_CHECK_GNU_READLINE

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_VPATH_LINK(examples/hello.mps)
AC_COIN_VPATH_LINK(examples/input.130)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/OsiClp/Makefile
                 test/Makefile
		 clp.pc
		 clp-uninstalled.pc])
		 
if test $coin_has_osi = yes ; then
  AC_CONFIG_FILES([src/OsiClp/osi-clp.pc
                   osi-clp-uninstalled.pc:src/OsiClp/osi-clp-uninstalled.pc.in])
fi

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_clp.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
