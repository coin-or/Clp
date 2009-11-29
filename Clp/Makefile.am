# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

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

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

addlibsfile = clp_addlibs.txt
addlibsdir = $(prefix)/share/doc/coin/$(PACKAGE_NAME)/

install-exec-local: install-doc
	$(install_sh_DATA) $(addlibsfile) $(DESTDIR)$(addlibsdir)/$(addlibsfile)

uninstall-local: uninstall-doc
	rm -f $(DESTDIR)$(addlibsdir)/$(addlibsfile)

pkgconfiglibdir = $(prefix)/lib/pkgconfig
pkgconfiglib_DATA = clp.pc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include ../BuildTools/Makemain.inc
