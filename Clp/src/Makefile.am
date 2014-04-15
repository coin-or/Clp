# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libClp                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libClp.la libClpSolver.la

# List all source files for this library, including headers
libClp_la_SOURCES = \
	ClpConfig.h \
	ClpCholeskyBase.cpp ClpCholeskyBase.hpp \
	ClpCholeskyDense.cpp ClpCholeskyDense.hpp \
	ClpConstraint.cpp ClpConstraint.hpp \
	ClpConstraintLinear.cpp ClpConstraintLinear.hpp \
	ClpConstraintQuadratic.cpp ClpConstraintQuadratic.hpp \
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
	ClpNode.cpp ClpNode.hpp \
	ClpObjective.cpp ClpObjective.hpp \
	ClpPackedMatrix.cpp ClpPackedMatrix.hpp \
	ClpParameters.hpp \
	ClpPlusMinusOneMatrix.cpp ClpPlusMinusOneMatrix.hpp \
	ClpPredictorCorrector.cpp ClpPredictorCorrector.hpp \
	ClpPdco.cpp ClpPdco.hpp \
	ClpPdcoBase.cpp ClpPdcoBase.hpp \
	ClpLsqr.cpp ClpLsqr.hpp \
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

if COIN_HAS_ABC
libClp_la_SOURCES += \
	AbcCommon.hpp \
	AbcDualRowDantzig.cpp AbcDualRowDantzig.hpp \
	AbcDualRowPivot.cpp AbcDualRowPivot.hpp \
	AbcDualRowSteepest.cpp AbcDualRowSteepest.hpp \
	AbcMatrix.cpp AbcMatrix.hpp \
	AbcNonLinearCost.cpp AbcNonLinearCost.hpp \
	AbcPrimalColumnDantzig.cpp AbcPrimalColumnDantzig.hpp \
	AbcPrimalColumnPivot.cpp AbcPrimalColumnPivot.hpp \
	AbcPrimalColumnSteepest.cpp AbcPrimalColumnSteepest.hpp \
	AbcSimplex.cpp AbcSimplex.hpp \
	AbcSimplexDual.cpp AbcSimplexDual.hpp \
	AbcSimplexPrimal.cpp AbcSimplexPrimal.hpp \
	AbcSimplexParallel.cpp \
	AbcSimplexFactorization.cpp AbcSimplexFactorization.hpp \
	AbcWarmStart.cpp AbcWarmStart.hpp \
	CoinAbcBaseFactorization.hpp \
	CoinAbcBaseFactorization1.cpp \
	CoinAbcBaseFactorization2.cpp \
	CoinAbcBaseFactorization3.cpp \
	CoinAbcBaseFactorization4.cpp \
	CoinAbcBaseFactorization5.cpp \
	CoinAbcCommonFactorization.hpp \
	CoinAbcCommon.hpp \
	CoinAbcDenseFactorization.cpp \
	CoinAbcDenseFactorization.hpp \
	CoinAbcFactorization.hpp \
	CoinAbcFactorization1.cpp \
	CoinAbcFactorization2.cpp \
	CoinAbcFactorization3.cpp \
	CoinAbcFactorization4.cpp \
	CoinAbcFactorization5.cpp \
	CoinAbcHelperFunctions.cpp CoinAbcHelperFunctions.hpp \
	CoinAbcOrderedFactorization1.cpp \
	CoinAbcOrderedFactorization2.cpp \
	CoinAbcOrderedFactorization3.cpp \
	CoinAbcOrderedFactorization4.cpp \
	CoinAbcOrderedFactorization5.cpp \
	CoinAbcSmallFactorization1.cpp \
	CoinAbcSmallFactorization2.cpp \
	CoinAbcSmallFactorization3.cpp \
	CoinAbcSmallFactorization4.cpp \
	CoinAbcSmallFactorization5.cpp
endif

libClpSolver_la_SOURCES = \
	ClpSolver.cpp \
	CbcOrClpParam.cpp CbcOrClpParam.hpp \
	Clp_ampl.cpp Clp_ampl.h \
	MyEventHandler.cpp MyEventHandler.hpp \
	MyMessageHandler.cpp MyMessageHandler.hpp \
	unitTest.cpp

# if AMD, CHOLMOD, or GLPK is available, then compile ClpCholeskyUfl.cpp
if COIN_HAS_AMD
libClp_la_SOURCES += ClpCholeskyUfl.cpp ClpCholeskyUfl.hpp
else
if COIN_HAS_CHOLMOD
libClp_la_SOURCES += ClpCholeskyUfl.cpp ClpCholeskyUfl.hpp
else
if COIN_HAS_GLPK
libClp_la_SOURCES += ClpCholeskyUfl.cpp ClpCholeskyUfl.hpp
endif
endif
endif

if COIN_HAS_MUMPS
libClp_la_SOURCES += ClpCholeskyMumps.cpp ClpCholeskyMumps.hpp
endif

if COIN_HAS_WSMP
libClp_la_SOURCES += ClpCholeskyWssmp.cpp ClpCholeskyWssmp.hpp ClpCholeskyWssmpKKT.cpp ClpCholeskyWssmpKKT.hpp
endif

# List all additionally required libraries
if DEPENDENCY_LINKING
libClp_la_LIBADD = $(CLPLIB_LIBS)
libClpSolver_la_LIBADD = $(CLPLIB_LIBS) libClp.la
endif

# This is for libtool
libClp_la_LDFLAGS = $(LT_LDFLAGS)
libClpSolver_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            clp program                               #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = clp

# List all source files for this executable, including headers
clp_SOURCES = ClpMain.cpp 

# List all additionally required libraries
clp_LDADD = libClpSolver.la libClp.la $(CLPLIB_LIBS) $(ASL_LIBS)

# List all dependency libraries (similar to LDADD, but without -l, -L flags)
clp_DEPENDENCIES = libClpSolver.la libClp.la $(CLPLIB_DEPENDENCIES)

########################################################################
#                            Additional flags                          #
########################################################################

# CbcOrClpParam seem to require COIN_HAS_CLP so that it knows that it is build within Clp
AM_CPPFLAGS = $(CLPLIB_CFLAGS) -DCOIN_HAS_CLP

if COIN_HAS_CHOLMOD
AM_CPPFLAGS += -I`$(CYGPATH_W) $(CHOLMODINCDIR)`
endif

if COIN_HAS_AMD
AM_CPPFLAGS += -I`$(CYGPATH_W) $(AMDINCDIR)`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	Clp_C_Interface.h \
	ClpCholeskyBase.hpp \
	ClpCholeskyDense.hpp \
	ClpConstraint.hpp \
	ClpConstraintLinear.hpp \
	ClpConstraintQuadratic.hpp \
	ClpDualRowDantzig.hpp \
	ClpDualRowPivot.hpp \
	ClpDualRowSteepest.hpp \
	ClpDummyMatrix.hpp \
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
	ClpNode.hpp \
	ClpObjective.hpp \
	ClpPackedMatrix.hpp \
	ClpPdcoBase.hpp \
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
	ClpSimplexDual.hpp \
	ClpSimplexPrimal.hpp \
	ClpSolve.hpp \
	CbcOrClpParam.hpp \
	Idiot.hpp

if COIN_HAS_ABC
includecoin_HEADERS += AbcSimplex.hpp CoinAbcCommon.hpp AbcCommon.hpp AbcNonLinearCost.hpp
endif


# if AMD, CHOLMOD, or GLPK is available, then install ClpCholeskyUfl.hpp (for advanced users)
if COIN_HAS_AMD
includecoin_HEADERS += ClpCholeskyUfl.hpp
else
if COIN_HAS_CHOLMOD
includecoin_HEADERS += ClpCholeskyUfl.hpp
else
if COIN_HAS_GLPK
includecoin_HEADERS += ClpCholeskyUfl.hpp
endif
endif
endif

if COIN_HAS_MUMPS
includecoin_HEADERS += ClpCholeskyMumps.hpp
endif

if COIN_HAS_WSMP
includecoin_HEADERS += ClpCholeskyWssmp.hpp ClpCholeskyWssmpKKT.hpp
endif

# needed by Clp
includecoin_HEADERS += CbcOrClpParam.cpp

########################################################################
#                     Installing manifest (MSVC++)                     #
########################################################################

if COIN_CXX_IS_CL
install-exec-hook:
	test -s clp.exe.manifest && \
	mt -manifest clp.exe.manifest -outputresource:clp.exe;\
	cp clp.exe $(bindir)/clp.exe
endif

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) config_clp.h $(DESTDIR)$(includecoindir)/ClpConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/ClpConfig.h
