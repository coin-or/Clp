# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 1152 2007-12-28 03:59:07Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                             libOsiClp                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiClp.la

# List all source files for this library, including headers
libOsiClp_la_SOURCES = OsiClpSolverInterface.cpp OsiClpSolverInterface.hpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libOsiClp_la_LIBADD = $(OSICLPLIB_LIBS) ../libClp.la
endif

# This is for libtool (on Windows)
libOsiClp_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(COINUTILS_CFLAGS) $(OSI_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = OsiClpSolverInterface.hpp
