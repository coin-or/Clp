# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Clp

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

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
ifeq ($(VolDefine),COIN_HAS_VOL)
	(cd $(CoinDir)/Vol && $(MAKE) install)
endif 
unitTest: 
	(cd Test && ${MAKE} unitTest)

clean: 
	@rm -rf Junk
	@rm -rf $(UNAME)*
	@rm -rf dep
	@rm -rf Test/Junk
	@rm -rf Test/$(UNAME)*
	@rm -rf Test/dep
	@rm -f clp

doc:
	doxygen $(MakefileDir)/doxygen.conf
