# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export MakefileDir := $(CoinDir)/Common/make
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

IncDir :=
IncDir += ${OsiSolverIncDir}
IncDir += ${OsiIncDir}
IncDir += ${ClpIncDir}
IncDir += ${zlibIncDir}
IncDir += ${bzlibIncDir}

LibDir :=
LibDir += ${OsiLibDir}
LibDir += ${ClpLibDir}
LibDir += ${zlibLibDir}
LibDir += ${bzlibLibDir}

LibName :=
LibName += ${OsiLibName}
LibName += ${ClpLibName}
LibName += ${zlibLibName}
LibName += ${bzlibLibName}

Define :=
Define += ${ClpDefine}
Define += ${OsiDefine}
Define += ${zlibDefine}
Define += ${bzlibDefine}

##############################################################################
# You should not need to edit below this line.
##############################################################################

CXXFLAGS += $(OPTFLAG)
ifeq ($(OptLevel),-g)
    CXXFLAGS += -DOSI_DEBUG
endif

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
	CXXFLAGS +=	-DREADLINE
	LDFLAGS += -lhistory -lreadline -ltermcap
endif

###############################################################################

TESTSRC := ClpMain.cpp
TESTSRC += unitTest.cpp
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

###############################################################################
# This one must be right before the rules

include ${MakefileDir}/Makefile.rules

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean doc unitTest libdepend

default: unitTest

###############################################################################

unitTest: $(TARGETDIR)/clp

libdepend :
	(cd .. && ${MAKE} install)

$(TARGETDIR)/clp : libdepend $(TESTOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating clp
	@mkdir -p $(TARGETDIR)
	@rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJ) $(LDFLAGS) $(SYSLD) -lm
	${CP} $@ ..

###############################################################################

doc:
	doxygen $(MakefileDir)/doxygen.conf

clean :
	rm -rf Junk
	@rm -rf $(DEPDIR)
	@rm -rf $(TARGETDIR)

###############################################################################

%::
	@mkdir -p Junk
	touch Junk/$(notdir $@)

###############################################################################

-include $(TESTDEP)
