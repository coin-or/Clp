# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libClp                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libClp.la

# List all source files for this library, including headers
libClp_la_SOURCES = \
	ClpCholeskyBase.cpp ClpCholeskyBase.hpp \
	ClpCholeskyDense.cpp ClpCholeskyDense.hpp \
	ClpCholeskyUfl.cpp ClpCholeskyUfl.hpp \
	Clp_C_Interface.cpp Clp_C_Interface.h \
	ClpDualRowDantzig.cpp ClpDualRowDantzig.hpp \
	ClpDualRowPivot.cpp ClpDualRowPivot.hpp \
	ClpDualRowSteepest.cpp ClpDualRowSteepest.hpp \
	ClpDummyMatrix.cpp ClpDummyMatrix.hpp \
	ClpDynamicExampleMatrix.cpp ClpDynamicExampleMatrix.hpp \
	ClpDynamicMatrix.cpp ClpDynamicMatrix.hpp \
	ClpEventHandler.cpp ClpEventHandler.hpp \
	ClpFactorization.cpp ClpFactorization.hpp \
	ClpGubDynamicMatrix.cpp ClpGubDynamicMatrix.hpp \
	ClpGubMatrix.cpp ClpGubMatrix.hpp \
	ClpHelperFunctions.cpp ClpHelperFunctions.hpp \
	ClpInterior.cpp ClpInterior.hpp \
	ClpLinearObjective.cpp ClpLinearObjective.hpp \
	ClpMatrixBase.cpp ClpMatrixBase.hpp \
	ClpMessage.cpp ClpMessage.hpp \
	ClpModel.cpp ClpModel.hpp \
	ClpNetworkBasis.cpp ClpNetworkBasis.hpp \
	ClpNetworkMatrix.cpp ClpNetworkMatrix.hpp \
	ClpNonLinearCost.cpp ClpNonLinearCost.hpp \
	ClpObjective.cpp ClpObjective.hpp \
	ClpPackedMatrix.cpp ClpPackedMatrix.hpp \
	ClpParameters.hpp \
	ClpPlusMinusOneMatrix.cpp ClpPlusMinusOneMatrix.hpp \
	ClpPredictorCorrector.cpp ClpPredictorCorrector.hpp \
	ClpPresolve.cpp ClpPresolve.hpp \
	ClpPrimalColumnDantzig.cpp ClpPrimalColumnDantzig.hpp \
	ClpPrimalColumnPivot.cpp ClpPrimalColumnPivot.hpp \
	ClpPrimalColumnSteepest.cpp ClpPrimalColumnSteepest.hpp \
	ClpQuadraticObjective.cpp ClpQuadraticObjective.hpp \
	ClpSimplex.cpp ClpSimplex.hpp \
	ClpSimplexDual.cpp ClpSimplexDual.hpp \
	ClpSimplexNonlinear.cpp ClpSimplexNonlinear.hpp \
	ClpSimplexOther.cpp ClpSimplexOther.hpp \
	ClpSimplexPrimal.cpp ClpSimplexPrimal.hpp \
	ClpSolve.cpp ClpSolve.hpp \
	Idiot.cpp Idiot.hpp \
	IdiSolve.cpp

# This is for libtool (on Windows)
libClp_la_LDFLAGS = -no-undefined

########################################################################
#                        clp_solver program                            #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = clp_solver

# List all source files for this executable, including headers
clp_solver_SOURCES = \
	ClpMain.cpp \
	CbcOrClpParam.cpp CbcOrClpParam.hpp \
	MyEventHandler.cpp MyEventHandler.hpp \
	MyMessageHandler.cpp MyMessageHandler.hpp \
	unitTest.cpp

# List all additionally required COIN libraries
clp_solver_LDADD = libClp.la \
	$(COINBASEOBJDIR)/src/libCoinBase.la

# And finally we add all missing libraries.  ADDLIBS is set by configure to
# include all libraries that are determined to be necessary (such as -lgz).
# We also need to list the additional libraries for the other COIN libraries
# used.  Finally, the -rpath flag is used by libtool to make sure that
# the shared library is found (in the lib install directory) when we are
# using dynamic libraries.
clp_solver_LDFLAGS = $(ADDLIBS) `cat $(COINBASEOBJDIR)/coinbase_addlibs.txt` \
	-rpath $(libdir)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(COINBASESRCDIR)/src`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	ClpDualRowDantzig.hpp \
	ClpDualRowPivot.hpp \
	ClpDualRowSteepest.hpp \
	ClpEventHandler.hpp \
	ClpFactorization.hpp \
	ClpInterior.hpp \
	ClpLinearObjective.hpp \
	ClpMatrixBase.hpp \
	ClpMessage.hpp \
	ClpModel.hpp \
	ClpNetworkMatrix.hpp \
	ClpNonLinearCost.hpp \
	ClpObjective.hpp \
	ClpPackedMatrix.hpp \
	ClpPlusMinusOneMatrix.hpp \
	ClpParameters.hpp \
	ClpPresolve.hpp \
	ClpPrimalColumnDantzig.hpp \
	ClpPrimalColumnPivot.hpp \
	ClpPrimalColumnSteepest.hpp \
	ClpSimplex.hpp \
	ClpSimplexNonlinear.hpp \
	ClpSimplexOther.hpp \
	ClpSimplexPrimal.hpp \
	ClpSolve.hpp \
	Idiot.hpp

########################################################################
#     Generating the configure header that is to be distributed        #
########################################################################

# You should only need to adapt the following three lines
CONFIG_DIST = ClpConfig.hpp
CONFIG_H = config_clp.h
$(libClp_la_OBJECTS): ../inc/$(CONFIG_DIST)

nodist_include_HEADERS = ../inc/$(CONFIG_DIST)

../inc/$(CONFIG_DIST): ../inc/$(CONFIG_H)
	echo "/** Required defines from $(CONFIG_H) */" >bla.tmp
	echo "#ifndef __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	echo "#define __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	$(EGREP) 'COIN_' $< >>bla.tmp
	echo "#endif" >>bla.tmp
	if test -r ../inc/$(CONFIG_DIST); then \
	  if diff ./inc/$(CONFIG_DIST) bla.tmp >/dev/null 2>1 ; then \
	    rm -f bla.tmp; \
	  else \
	    rm -f ../inc/$(CONFIG_DIST); \
	    mv bla.tmp ../inc/$(CONFIG_DIST); \
	  fi; \
	else \
	  mv bla.tmp ../inc/$(CONFIG_DIST); \
	fi

CLEANFILES = bla.tmp

DISTCLEANFILES = ../inc/$(CONFIG_DIST)
