# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest libdepend libClp doc

default: install

libdepend:
	(cd $(CoinDir)/Coin && $(MAKE) install)

install library: libdepend
	${MAKE} -f Makefile.Clp $@

libClp: libdepend
	${MAKE} -f Makefile.Clp library

unitTest: 
	(cd Test && ${MAKE} unitTest)

clean: 
	rm -rf Junk
	@rm -rf $(DEPDIR)
	@rm -rf $(TARGETDIR)
