# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean unitTest libdepend

default: unitTest

libdepend :
	(cd .. && ${MAKE} install)

unitTest: libdepend
	${MAKE} -f Makefile.test $*

clean:
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
