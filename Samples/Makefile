OptLevel := -g

DRIVER ?= driver

# Possible choices for driver:
#   decompose
#   defaults
#   driver
#   driver2
#   dualCuts
#   ekk
#   ekk_interface
#   minimum
#   network
#   piece
#   sprint
#   sprint2
#   testBarrier
#   testBasis
#   useVolume


# The location of the customized Makefiles
CoinDir := ../../
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

IncDir :=
IncDir += ${ClpIncDir}
IncDir += ${CoinIncDir}
IncDir += ${zlibIncDir}
IncDir += ${bzlibIncDir}
IncDir += ${readlineIncDir}
IncDir += ${lapackIncDir}
IncDir += ${taucsIncDir}

LibDir :=
LibDir += ${ClpLibDir}
LibDir += ${CoinLibDir}
LibDir += ${zlibLibDir}
LibDir += ${bzlibLibDir}
LibDir += ${readlineLibDir}
LibDir += ${lapackLibDir}
#LibDir += ${taucsLibDir}
LibDir += $(GlpkLibDir)

LibName :=
LibName += ${ClpLibName}
LibName += ${CoinLibName}
LibName += ${zlibLibName}
LibName += ${bzlibLibName}
LibName += ${readlineLibName}
LibName += ${lapackLibName}
#LibName += ${taucsLibName}
LibName += $(GlpkLibName)

Define :=
Define += ${ClpDefine}
Define += ${CoinDefine}
Define += ${zlibDefine}
Define += ${bzlibDefine}
Define += ${readlineDefine}
Define += ${lapackDefine}
Define += ${taucsDefine}

CXXFLAGS += $(OPTFLAG)

###############################################################################

space:= $(empty) $(empty)
OptVersion := $(subst $(space),_,$(OptLevel))

TARGETDIR := $(UNAME)$(OptVersion)
DEPDIR := dep

VPATH := . : include : Junk : ${TARGETDIR} : ${DEPDIR}

#########################################################################

CXXFLAGS += $(addprefix -I,${IncDir})
CXXFLAGS += $(addprefix -D,${Define})

LIBDIRS := ${LibDir}
LIBS    := ${LibName}

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(call ifprefix,$(SHLINKPREFIX),$(LIBDIRS))
LDFLAGS += $(patsubst lib%,-l%,$(basename $(LIBS)))

#LDFLAGS += -lefence
#LDFLAGS += -Wl,-static
#LDFLAGS += -lwsmpP4 -lpthread
#OSL_LIB_DIR = /home/forrest/osl/lib
#LDFLAGS += -L $(OSL_LIB_DIR) -Wl,-rpath,$(OSL_LIB_DIR) -losl-O-nolic-native-32
#LDFLAGS += -pg
###############################################################################

ifeq ($(DRIVER),driverC)
    SRCEXT := c
else
    SRCEXT := cpp
endif

TESTSRC := $(DRIVER).$(SRCEXT)
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.$(SRCEXT)=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.$(SRCEXT)=.d))

###############################################################################
# This one must be right before the rules

include ${MakefileDir}/Makefile.rules

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean doc unitTest libdepend

default: testit

###############################################################################


testit : libdepend $(TESTOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating testit
	@rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJ) $(LDFLAGS) $(SYSLD) -lm

###############################################################################

-include $(TESTDEP)
