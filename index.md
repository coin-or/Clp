CLP User's Guide
================

Authors:
- John Forrest, jjforre at us dot ibm dot com, IBM Research
- David de la Nuez, dmd57 at columbia dot edu, Columbia University & IBM Research
- Robin Lougee-Heimer, robinlh at us dot ibm dot com, IBM Research

Copyright: 2004 IBM Corporation

License: CLP and this documentation are provided under the terms of the
[Eclipse Public License](https://opensource.org/licenses/EPL-1.0) (EPL).
Any use, reproduction or distribution of the programs constitutes the
recipient's acceptance of the license.
The EPL is approved by the [Open Source Initiative](http://opensource.org/).
The Eclipse Foundation, the steward of the EPL, has an EPL FAQ available
which is based on the Eclipse Foundation's understanding of the EPL.

Welcome to CLP!
===============

The COIN-OR Linear Program code or CLP is an open-source simplex solver
written in C++. It is primarily meant to be used as a callable library,
but a basic, stand-alone [executable](./clpexe) version is also
available.

There are a number of resources available to help new CLP users get
started. This document is designed to be used in conjunction with the
files in the examples subdirectory of the main CLP directory.
The Examples illustrate how to use CLP and may also
serve as useful starting points for user projects. In the rare event
that either this document or the available [Doxygen content](https://www.coin-or.org/Doxygen/Clp/)
conflicts with the observed behavior of the source code, the comments in
the Clp header files are the ultimate reference.

Prerequisites
=============

CLP is written in C++, so it is expected that users of CLP will be
writing C++ programs which use CLP as a library. Thus a working
knowledge of [C++](http://www.cplusplus.com/doc/tutorial/), including
basic object-oriented programming terminology is assumed in this
document. In addition, the user should be familiar with the fundamental
concepts of [Linear Programming](http://carbon.cudenver.edu/~hgreenbe/courseware/LPshort/intro.html).

Content
=======

- [Basic Model Classes](./basicmodelclasses)
- [Not-Quite-So-Basic Model Classes](./notsobasic)
- [More Samples](./moresamples)
- [The CLP Executable](./clpexe)
- [Messages](./messages)
- [FAQ](./faq)
- [Doxygen Documentation](Doxygen)
- [Revision History](./revhist)
