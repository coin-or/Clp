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

libClp: library

install library: libdepend
	${MAKE} -f Makefile.Clp $@

libdepend:
	(cd $(CoinDir)/Coin && $(MAKE) install)

unitTest: 
	(cd Test && ${MAKE} unitTest)

clean: 
	@rm -rf Junk
	@rm -rf $(UNAME)*
	@rm -rf dep

doc:
	doxygen $(MakefileDir)/doxygen.conf
