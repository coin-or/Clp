
# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Clp

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest libdepend libClp doc

default: install
libClp: library

libdepend:
	(cd $(CoinDir)/Coin && $(MAKE) install)

install library: libdepend
	${MAKE} -f Makefile.Clp $@

doc:
	doxygen $(MakefileDir)/doxygen.conf

unitTest: 
	(cd Test && ${MAKE} unitTest)

clean: 
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
