# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -O2
# I seem to need this at present
OPTFLAG := -O2
OptLevel := -g
ifeq ($(OptLevel),-g)
    CXXFLAGS += -DCLP_DEBUG
endif

LIBNAME := Clp
LIBSRC :=
LIBSRC += ClpDualRowDantzig.cpp
LIBSRC += ClpDualRowPivot.cpp
LIBSRC += ClpDualRowSteepest.cpp
LIBSRC += ClpFactorization.cpp
#LIBSRC += ClpMalloc.cpp
LIBSRC += ClpMatrixBase.cpp
LIBSRC += ClpMessage.cpp
LIBSRC += ClpModel.cpp
LIBSRC += ClpNonLinearCost.cpp
LIBSRC += ClpPackedMatrix.cpp
LIBSRC += ClpPrimalColumnDantzig.cpp
LIBSRC += ClpPrimalColumnPivot.cpp
LIBSRC += ClpPrimalColumnSteepest.cpp
LIBSRC += ClpSimplex.cpp
LIBSRC += ClpSimplexDual.cpp
LIBSRC += ClpSimplexPrimal.cpp

export CoinDir = $(shell cd ..; pwd)
##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export MakefileDir := ../Common/make
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location
#This modification seems to be needed
export ExtraIncDir := ../Osi/include
export ExtraLibDir := 
export ExtraLibName :=
export ExtraDefine := 

export LibType OptLevel LIBNAME LIBSRC

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libClp library clean doc

unitTest : install
	(cd Test && ${MAKE} unitTest)

default: install

install clean doc library: % :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

libClp:
	$(MAKE) -f ${MakefileDir}/Makefile.lib library
