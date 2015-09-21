# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Clp],[1.16],[clp@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

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
AC_COIN_PROJECTDIR_INIT(Clp)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Clp)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinUtils, [coinutils], [ClpLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHECK_PACKAGE(Osi,       [osi],       [OsiClpLib])
AC_COIN_CHECK_PACKAGE(OsiTests,  [osi-unittests])
AC_COIN_CHECK_PACKAGE(Sample,    [coindatasample])
AC_COIN_CHECK_PACKAGE(Netlib,    [coindatanetlib])

#############################################################################
#                                    Aboca                                  #
#############################################################################

#  1 - build Abc serial but no inherit code
#  2 - build Abc serial and inherit code
#  3 - build Abc cilk parallel but no inherit code
#  4 - build Abc cilk parallel and inherit code
AC_ARG_ENABLE([aboca],
              [AC_HELP_STRING([--enable-aboca],[enables build of Aboca solver (set to 1,2,3,4)])],
              [use_aboca=$enableval],
              [use_aboca=no])
if test "$use_aboca" = yes ; then
  use_aboca=1
fi
case "$use_aboca" in 
  1 | 2 | 3 | 4)
    AC_DEFINE_UNQUOTED([CLP_HAS_ABC], [$use_aboca], [Define to 1, 2, 3, or 4 if Aboca should be build.])
    ;;
  no) ;;
  *)
    AC_MSG_ERROR([invalid argument for --enable-aboca: $use_aboca])
    ;;
esac
AM_CONDITIONAL(COIN_HAS_ABC, test ! "$use_aboca" = no)

#############################################################################
#                        Third party linear solvers                         #
#############################################################################

# AMD from UFL
AC_COIN_CHECK_USER_LIBRARY(AMD, AMD, [amd.h], [amd_defaults],,[ClpLib])

# CHOLMOD from UFL
# CHOLMOD requires AMD, which can be given by AMD_LIBS
# CHOLMOD requires Blas and Lapack?
coin_save_LIBS="$LIBS"
LIBS="$LIBS $AMD_LIBS"
AC_COIN_CHECK_USER_LIBRARY(CHOLMOD, CHOLMOD, [cholmod.h], [cholmod_start],,[ClpLib])
LIBS="$coin_save_LIBS"

# Glpk also brings AMD
if test $coin_has_cholmod = false -a $coin_has_amd = false ; then
  AC_COIN_CHECK_PACKAGE(Glpk, [coinglpk], [ClpLib])
  if test $coin_has_glpk = yes ; then
    AC_MSG_NOTICE([using AMD from GLPK package])
    AC_DEFINE(COIN_HAS_AMD,[1],[Define to 1 if the AMD package is available])
  fi
else
  # for configure
  AM_CONDITIONAL(COIN_HAS_GLPK, [test 0 = 1])
fi

# MUMPS
AC_COIN_CHECK_PACKAGE(Mumps, [coinmumps], [ClpLib])

# WSMP
AC_ARG_WITH([wsmp],
            AC_HELP_STRING([--with-wsmp],
                           [specify WSMP library]),
            [use_wsmp=$withval], [use_wsmp=no])

if test "$use_wsmp" != "no"; then
  # Check how to link against Fortran libraries from C
  AC_COIN_PROG_F77
  AC_COIN_F77_WRAPPERS
  
  # WSMP requires Blas
  AC_COIN_CHECK_PACKAGE_BLAS(ClpLib)

  coin_save_LIBS="$LIBS"
  LIBS="$LIBS $use_wsmp $BLAS_LIBS"
  AC_LANG_PUSH([Fortran 77])
  AC_MSG_CHECKING([whether user-supplied WSMP library \"$use_wsmp\" works])
  AC_TRY_LINK([],[      call WSSMP()],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_RESULT(no)
               AC_MSG_ERROR([WSMP library $use_wsmp does not seem to work])])
  AC_LANG_POP([Fortran 77])
  LIBS="$coin_save_LIBS"

  AC_DEFINE_UNQUOTED([COIN_HAS_WSMP], [1], [Define to 1 if the WSMP package is available])
  
  CLPLIB_LIBS="$use_wsmp $CLPLIB_LIBS"
  CLPLIB_LIBS_INSTALLED="$use_wsmp $CLPLIB_LIBS_INSTALLED"
  CLPLIB_PCLIBS="$use_wsmp $CLPLIB_PCLIBS"
else
  # to please configure
  AM_CONDITIONAL(COIN_HAS_BLAS, [test 0 = 1])
fi
AM_CONDITIONAL(COIN_HAS_WSMP, [test "$use_wsmp" != no])

#############################################################################
#                             Other dependencies                            #
#############################################################################
# Ampl Solver library
AC_COIN_CHECK_PACKAGE(ASL, [coinasl], [ClpLib])

AC_COIN_CHECK_GNU_READLINE(ClpLib)

AC_COIN_CHECK_CXX_CHEADER(math)
AC_COIN_CHECK_CXX_CHEADER(float)
AC_COIN_CHECK_CXX_CHEADER(ieeefp)

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################

# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.

AC_COIN_VPATH_LINK(examples/hello.mps)
AC_COIN_VPATH_LINK(examples/input.130)
AC_COIN_VPATH_LINK(examples/g.tiny)
AC_COIN_VPATH_LINK(examples/gparm.tiny)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi)

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
  AC_CONFIG_FILES([osi-clp.pc:src/OsiClp/osi-clp.pc.in
                   osi-clp-uninstalled.pc:src/OsiClp/osi-clp-uninstalled.pc.in])
fi

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_clp.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
