# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g
OptLevel := -O2


LIBNAME := Clp
LIBSRC :=
LIBSRC += ClpDualRowDantzig.cpp
LIBSRC += ClpDualRowPivot.cpp
LIBSRC += ClpDualRowSteepest.cpp
LIBSRC += ClpFactorization.cpp
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
# and Presolve stuff
LIBSRC += Presolve.cpp
LIBSRC += PresolveDoubleton.cpp
LIBSRC += PresolveEmpty.cpp
LIBSRC += PresolveFixed.cpp
LIBSRC += PresolveMatrix.cpp
LIBSRC += PresolvePsdebug.cpp
LIBSRC += PresolveSingleton.cpp
LIBSRC += PresolveZeros.cpp
LIBSRC += PresolveDual.cpp		    
LIBSRC += PresolveDupcol.cpp		    
LIBSRC += PresolveForcing.cpp		    
LIBSRC += PresolveImpliedFree.cpp	    
LIBSRC += PresolveIsolated.cpp	    
LIBSRC += PresolveSubst.cpp		    
LIBSRC += PresolveTighten.cpp		    
LIBSRC += PresolveUseless.cpp             

LIBSRC += Idiot.cpp
LIBSRC += IdiSolve.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

# Say Idiot code can use Clp interface
CXXFLAGS += -DCLP_IDIOT
CXXFLAGS += -DUSE_PRESOLVE
ifeq ($(OptLevel),-g)
#     CXXFLAGS += -DCLP_DEBUG
#    CXXFLAGS += -DPRESOLVE_SUMMARY=1 -DDEBUG_PRESOLVE -DCHECK_CONSISTENCY=1
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
