# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED
#LibType := STATIC

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g
OptLevel := -O1

LIBNAME := Clp
LIBSRC :=
LIBSRC += ClpDualRowDantzig.cpp
LIBSRC += ClpDualRowPivot.cpp
LIBSRC += ClpDualRowSteepest.cpp
LIBSRC += ClpEventHandler.cpp
LIBSRC += ClpFactorization.cpp
LIBSRC += ClpMatrixBase.cpp
LIBSRC += ClpMessage.cpp
LIBSRC += ClpModel.cpp
LIBSRC += ClpNetworkBasis.cpp
LIBSRC += ClpNetworkMatrix.cpp
LIBSRC += ClpNonLinearCost.cpp
LIBSRC += ClpObjective.cpp
LIBSRC += ClpLinearObjective.cpp
LIBSRC += ClpPackedMatrix.cpp
LIBSRC += ClpPlusMinusOneMatrix.cpp
LIBSRC += ClpDummyMatrix.cpp
LIBSRC += ClpGubMatrix.cpp
LIBSRC += ClpDynamicMatrix.cpp
LIBSRC += ClpDynamicExampleMatrix.cpp
LIBSRC += ClpPrimalColumnDantzig.cpp
LIBSRC += ClpPrimalColumnPivot.cpp
LIBSRC += ClpPrimalColumnSteepest.cpp
LIBSRC += ClpQuadraticObjective.cpp
LIBSRC += ClpSimplex.cpp
LIBSRC += ClpSimplexDual.cpp
LIBSRC += ClpSimplexOther.cpp
LIBSRC += ClpSimplexPrimal.cpp
LIBSRC += ClpSimplexNonlinear.cpp
LIBSRC += ClpSolve.cpp
LIBSRC += ClpInterior.cpp
LIBSRC += ClpCholeskyBase.cpp
LIBSRC += ClpCholeskyDense.cpp
LIBSRC += ClpCholeskyUfl.cpp
#LIBSRC += ClpCholeskyWssmp.cpp
LIBSRC += ClpPredictorCorrector.cpp
LIBSRC += ClpHelperFunctions.cpp
LIBSRC += Clp_C_Interface.cpp
# and Presolve stuff
LIBSRC += ClpPresolve.cpp

LIBSRC += Idiot.cpp
LIBSRC += IdiSolve.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

# Say Idiot code can use Clp interface
#CXXFLAGS += -DCLP_IDIOT 
#CXXFLAGS += -DUSE_PRESOLVE
#CXXFLAGS += -DCORRECT_COLUMN_COUNTS
ifeq ($(OptLevel),-g)
#     CXXFLAGS += -DCLP_DEBUG
#CXXFLAGS += -DPRESOLVE_DEBUG=1
endif
ifeq ($(OptLevel),-O2)
#     CXXFLAGS += -DNDEBUG
endif

export ExtraIncDir  := ${CoinIncDir} ${zlibIncDir} ${bzlibIncDir}
export ExtraLibDir  := ${CoinLibDir} ${zlibLibDir} ${bzlibLibDir}
export ExtraLibName := ${CoinLibName} ${zlibLibName} ${bzlibLibName}
export ExtraDefine  := ${CoinDefine} ${zlibDefine} ${bzlibDefine}

export LibType OptLevel LIBNAME LIBSRC

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install library

default: install

install library:
	$(MAKE) -f ${MakefileDir}/Makefile.lib $@
