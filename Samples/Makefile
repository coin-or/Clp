OptLevel := -g

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
IncDir += ${lapackIncDir}

LibDir :=
LibDir += ${ClpLibDir}
LibDir += ${CoinLibDir}
LibDir += ${zlibLibDir}
LibDir += ${bzlibLibDir}
LibDir += ${lapackLibDir}

LibName :=
LibName += ${ClpLibName}
LibName += ${CoinLibName}
LibName += ${zlibLibName}
LibName += ${bzlibLibName}
LibName += ${lapackLibName}

Define :=
Define += ${ClpDefine}
Define += ${CoinDefine}
Define += ${zlibDefine}
Define += ${bzlibDefine}
Define += ${lapackDefine}

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

ifeq ($(findstring Linux,${UNAME}),Linux)
#	CXXFLAGS +=	-DREADLINE
#	LDFLAGS += -lhistory -lreadline -ltermcap
endif
#if DENSE and using given libraries
#LDFLAGS += -llapack -lblas -lg2c
#LDFLAGS += -lefence
###############################################################################

TESTSRC := $(DRIVER).cpp
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

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
