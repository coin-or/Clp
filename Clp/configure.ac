# Copyright (C) 2006, 2009 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Clp],[1.10.0],[coin-clp@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/ClpCholeskyBase.cpp)

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
#                              COIN components                              #
#############################################################################

#AC_COIN_HAS_MODULE_BLAS
#AC_COIN_HAS_MODULE_LAPACK
AC_COIN_HAS_MODULE(CoinUtils,  [coinutils = 2.5BSPsplit])
AC_COIN_HAS_MODULE(Sample,     [sample], [required=0])
AC_COIN_HAS_MODULE(Netlib,     [netlib], [required=0])

if test "$coin_has_coinutils" != yes ; then
  AC_COIN_HAS_PROJECT(CoinUtils)
  if test $coin_has_coinutils != notGiven && test $coin_has_coinutils != skipping ; then
    COINUTILS_CFLAGS="-I`${CYGPATH_W} ${COINUTILSSRCDIR}/src` -I`${CYGPATH_W} ${COINUTILSOBJDIR}/inc`"
    COINUTILS_LIBS="${COINUTILSOBJDIR}/src/libCoinUtils.la `cat ${COINUTILSOBJDIR}/coinutils_addlibs.txt`"
  else
    AC_MSG_ERROR("Required package CoinUtils unavailable.")
  fi
else
  #for configure
  AM_CONDITIONAL([COIN_HAS_COINUTILS_PREINSTALLED], [test 0 = 1])
fi

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
                 test/Makefile
		 clp_addlibs.txt
		 clp.pc
		 clp-uninstalled.pc])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_clp.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
