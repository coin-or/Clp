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
	ClpConfig.h \
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
libClp_la_LDFLAGS = $(LT_LDFLAGS) `$(CYGPATH_W) ../../CoinUtils/src/libCoinUtils.la`

########################################################################
#                            clp program                               #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = clp

# List all source files for this executable, including headers
clp_SOURCES = \
	ClpMain.cpp \
	CbcOrClpParam.cpp CbcOrClpParam.hpp \
	MyEventHandler.cpp MyEventHandler.hpp \
	MyMessageHandler.cpp MyMessageHandler.hpp \
	unitTest.cpp

# List all additionally required COIN libraries
clp_LDADD = libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# Here we add additional libraries
LIBS += $(ADDLIBS) `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
clp_LDFLAGS = -rpath $(libdir)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	../inc/config_clp.h \
	Clp_C_Interface.h \
	ClpCholeskyBase.hpp \
	ClpCholeskyDense.hpp \
	ClpCholeskyWssmp.hpp \
	ClpDualRowDantzig.hpp \
	ClpDualRowPivot.hpp \
	ClpDualRowSteepest.hpp \
	ClpDynamicExampleMatrix.hpp \
	ClpDynamicMatrix.hpp \
	ClpEventHandler.hpp \
	ClpFactorization.hpp \
	ClpGubDynamicMatrix.hpp \
	ClpGubMatrix.hpp \
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
	ClpQuadraticObjective.hpp \
	ClpSimplex.hpp \
	ClpSimplexNonlinear.hpp \
	ClpSimplexOther.hpp \
	ClpSimplexPrimal.hpp \
	ClpSolve.hpp \
	Idiot.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = ClpConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
