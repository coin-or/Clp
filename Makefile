# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Clp

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location
###############################################################################
BACKUP ?= NULL_BACKUP
BACKUP_NAME := Clp_$(BACKUP).tar

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest libdepend libClp doc backup dist restore

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


backup: dist
dist:
 
ifeq ($(BACKUP),NULL_BACKUP)
	@echo ""
	@echo "Syntax is make dist BACKUP=name or make backup BACKUP=name,"
	@echo "where Clp_name.tar.gz and Clp_name.tar do not exist"
	@echo Current version seems to be $(word 3 ,$(shell grep define Test/ClpMain.cpp|grep CLPVERSION))
	@echo ""
else
ifeq	($(words $(wildcard $(BACKUP_NAME)*)),0)
	@echo "Creating tar file " $(BACKUP_NAME)
	(cd .. && \
	tar -cvf Clp/$(BACKUP_NAME) Coin/include/Coin*.hpp Coin/Coin*.cpp\
	Clp/include/Clp_C_Interface.h Coin/include/Coin_C_defines.h\
	Clp/include/Clp*.hpp Clp/include/Idiot.hpp Clp/Clp*.cpp Clp/Idiot.cpp Clp/IdiSolve.cpp\
	Clp/Test/CbcOr*.?pp Clp/Test/MyEvent*.?pp Clp/Test/MyMessage*.?pp Clp/Test/ClpMain.cpp\
	Clp/Test/unitTest.cpp Clp/Test/Makefile Clp/Test/Makefile.test Clp/Samples/Makefile\
	Clp/Makefile.Clp Clp/Makefile Makefiles/Makefile.* Coin/Makefile)
	@echo "Tar file " $(BACKUP_NAME) "created"
else
	@echo "Unable to create as " $(wildcard $(BACKUP_NAME)*) "exist(s)"
endif
endif

restore:

ifeq ($(BACKUP),NULL_BACKUP)
ifneq	($(words $(wildcard Clp_*.tar Clp_*tar.gz)),0)
	@ls -ltr  Clp_*.tar*
	@echo ""
	@echo "Syntax is make restore BACKUP=xxxx where Clp_xxxx.tar.gz or"
	@echo "Clp_xxxx.tar exists (see list above)"
	@echo "WARNING - All COIN/Clp/Clp* and COIN/Coin/Coin* will be replaced"
else
	@echo ""
	@echo "Syntax is make restore BACKUP=xxxx where Clp_xxxx.tar.gz or"
	@echo "Clp_xxxx.tar exists"
	@echo "No candidate files exist"
	@echo "WARNING - All COIN/Clp/Clp* and COIN/Coin/Coin* will be replaced"
endif
else
ifeq	($(words $(wildcard $(BACKUP_NAME) $(BACKUP_NAME).gz)),1)
ifeq	($(CONFIRM),YES)
#        could delete stuff here?
ifeq	($(words $(wildcard $(BACKUP_NAME))),1)
	(cd .. && \
	tar --exclude Makefiles --exclude Coin/Makefile -xvf Clp/$(BACKUP_NAME))
	@echo "Restored from tar file " $(BACKUP_NAME) 
endif
ifeq	($(words $(wildcard $(BACKUP_NAME).gz)),1)
	(cd .. && \
	tar --exclude Makefiles --exclude Coin/Makefile -xzvf Clp/$(BACKUP_NAME).gz)
	@echo "Restored from tar file " $(BACKUP_NAME).gz 
endif
	@rm -rf ../Coin/$(UNAME)* 
	@rm -rf ./$(UNAME)* 
	@rm -rf ./Samples//$(UNAME)* 
	@rm -rf ./Test/$(UNAME)* 
else
	@echo "Could restore from tar file " $(wildcard $(BACKUP_NAME) $(BACKUP_NAME).gz)
	@echo "Re-enter make with CONFIRM=YES"
endif
else
ifeq	($(words $(wildcard $(BACKUP_NAME) $(BACKUP_NAME).gz)),0)
	@echo "Unable to restore as none of" $(BACKUP_NAME) or $(BACKUP_NAME).gz exist
ifneq	($(words $(wildcard Clp_*.tar Clp_*tar.gz)),0)
	@echo candidates are
	@ls -ltr  Clp_*.tar*
	@echo ""
	@echo "Syntax is make restore BACKUP=name where name is of form Clp_*.tar.gz or"
	@echo "Clp_*.tar and exists (in given list)"
	@echo "WARNING - All COIN/Clp/Clp* and COIN/Coin/Coin* will be replaced"
else
	@echo ""
	@echo "No candidate files exist"
	@echo "Syntax is make restore BACKUP=name where name is of form Clp_*.tar.gz or"
	@echo "Clp_*.tar and exists"
	@echo "WARNING - All COIN/Clp/Clp* and COIN/Coin/Coin* will be replaced"
endif
else
	@echo Ambiguous backups $(wildcard $(BACKUP_NAME) $(BACKUP_NAME).gz)
endif
endif
endif

