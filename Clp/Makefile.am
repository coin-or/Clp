# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# Don't do src/OsiClp recursively, since src/OsiClp/libOsiClp depends on src/libClp
if COIN_HAS_OSI
  SUBDIRS += src/OsiClp
endif

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	doc/authors.xml \
	doc/basicmodelclasses.xml \
	doc/clpexe.xml \
	doc/clpuserguide.xml \
	doc/doxygen.xml \
	doc/faqcontent.xml \
	doc/faq.xml \
	doc/intro.xml \
	doc/legal.xml \
	doc/messages.xml \
	doc/moresamples.xml \
	doc/notsobasic.xml \
	doc/revhist.xml \
	doc/Howto/docbook4clp.xml \
	examples/addBits.cpp \
	examples/addColumns.cpp \
	examples/addRows.cpp \
	examples/decompose.cpp \
	examples/defaults.cpp \
	examples/driver2.cpp \
	examples/driverC.c \
	examples/driver.cpp \
	examples/dualCuts.cpp \
	examples/ekk.cpp \
	examples/ekk_interface.cpp \
	examples/hello.cpp \
	examples/hello.mps \
	examples/input.130 \
	examples/INSTALL \
	examples/makeDual.cpp \
	examples/Makefile.in \
	examples/minimum.cpp \
	examples/network.cpp \
	examples/piece.cpp \
	examples/rowColumn.cpp \
	examples/sprint2.cpp \
	examples/sprint.cpp \
	examples/testBarrier.cpp \
	examples/testBasis.cpp \
	examples/testGub2.cpp \
	examples/testGub.cpp \
	examples/testQP.cpp \
	examples/useVolume.cpp

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

# Doxygen documentation

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )


clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                    Creation of the addlibs file                      #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = clp.pc

if COIN_HAS_OSI
pkgconfiglib_DATA += osi-clp.pc
endif

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Clp

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@:$(DESTDIR)$(pkgconfiglibdir) \
	$(PKG_CONFIG) --libs clp > $(addlibsdir)/clp_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libClpSolver.lib libClp.lib @CLPLIB_LIBS_INSTALLED@" > $(addlibsdir)/clp_addlibs.txt
else
	echo -L@abs_lib_dir@ -lClp @CLPLIB_LIBS_INSTALLED@ > $(addlibsdir)/clp_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/clp_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
