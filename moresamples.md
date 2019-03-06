More Samples {#moreexamples}
============

CLP\'s Samples Directory
========================

The CLP dsitribution includes a number of `.cpp` sample files. Users are
encouraged to use them as starting points for their own CLP projects.
The files can be found in the `` directory. For the latest information
on compiling and running these samples, please see the file
`CLPSAMPLESDIRINSTALL`. Below is a list of some of the most useful
sample files with a short description for each file.

  Source file          Description
  -------------------- ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  MINIMUMCPP           This is a CLP \"Hello, world\" program. It reads a problem from an MPS file, and solves the problem. \[[More\...](#minimumcppdesc)\]
  DEFAULTSCPP          This is one of the simpler driver programs available. It sets tolerances to defaults and is a good place to find straightforward uses of \"set\" and \"get\" methods. It also prints out full MPS-like solutions. \[[More\...](#defaultscppdesc)\]
  DRIVERCPP            This is designed to be a file that a user could modify to get a useful driver program for his or her project. In particular, it demonstrates the use of CLP\'s presolve functionality. \[[More\...](#drivercppdesc)\]
  NETWORKCPP           This shows the use of non-standard matrices and how to load a problem without the use of MPS files. \[[More\...](#networkcppdesc)\]
  TESTBARRIERCPP       This is a basic driver file for the barrier method of CLP, similar to MINIMUMCPP. The barrier method is not currently addressed in this guide. \[[More\...](#testbarriercppdesc)\]

  : Basic Samples

  Source file          Description
  -------------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  DRIVER2CPP           This sample, in addition to some tasks common to other samples, does some advanced message handling and presolve.
  DUALCUTSCPP          This sample implements a method of treating a problem as a collection of cuts.
  DECOMPOSECPP         This does full Dantzig-Wolfe decomposition. It illustrates the use of many models, adding columns, et cetera.
  SPRINTCPP            This solves a long, thin problem by solving smaller subsets. It is a simplified version of work done by one of the authors on aircrew scheduling problems. It shows the use of two models and their synchronization. A more general version can be found in `COIN/Clp/ClpSolve.cpp`
  SPRINT2CPP           This is similar to `sprint.cpp` but is designed for solving large problems with little choice. The idea is that if relatively few variables are fixed, presolve can greatly reduce the problem size so that a series of solves can get close to the optimal solution much faster than would a naïve solve of the full problem.

  : Advanced Samples

The remaining Samples listed here are considered unsupported in that
they are of a more esoteric nature and are sometimes contributed as a
result of an individual\'s request. The are to be found in
`CLPSAMPLESDIRContributed`.

  Source file          Description
  -------------------- -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  TESTBASISCPP         This sample takes a problem, changes any inequality constraints to equality constraints, solves the problem, and creates the optimal basis.
  TESTGUBCPP           This sample illustrates the use of the GUB (\"Generalized Upper Bound\") technique.
  EKKCPP               This sample can be used to compare CLP and OSL. It uses an additional file in the Samples directory, `ekk_interface.cpp`. These sample files are not likely to be interesting to new CLP users who do not have experience with OSL.
  HELLOCPP             This sample creates a text-based picture of a matrix on screen (limited to an 80x80 matrix). It\'s not terribly useful but it does illustrate one way to step through the elements of a matrix.
  PIECECPP             This sample takes a matrix read in by `CoinMpsIo` (can be used to read in MPS files without a solver), deletes every second column and solves the resulting problem.
  USEVOLUMECPP         The Volume Algorithm is another solver available as part of the COIN-OR distribution. This sample shows how to use the Volume Algorithm with CLP.

  : Unsupported Samples

minimum.cpp {#minimumcppdesc}
-----------

This sample is examined in more detail in [???](#firstexample).

defaults.cpp {#defaultscppdesc}
------------

This sample begins by reading an MPS file. The default MPS file is
`COIN/Mps/Sample/p0033.mps`; this can be over-riden by a command-line
specification of a (path and) file name). The sample then sets the pivot
algorithm to be exact devex. It \"gets\" the default infeasibility cost
and \"sets\" it to that value (and prints it to standard out). This sort
of getting and setting of various parameters constitutes a common theme
in this sample, with the purpose of illustrating usage of some of the
more common get and set methods available in CLP.

At this point the model is solved by the primal method. A sequence of
sets, gets and prints is then followed by a number of calls to methods
which give specific information about the status of the problem (for
example, the code checks that the current solution has been proven to be
optimal by `assert(model.isProvenOptimal())`).

Next, a copy of the original model is made. More sets and gets are
performed to demonstrate the use of additional options (including the
setting of the default message handling as well as changing of the \"log
level\" (amount of output)). The model is solved again a number of times
between changes of the optimization direction (i.e. changing from min to
max or vice versa). The remaining lines of this sample serve to display
solution and problem information in much the same way as is done in
driver.cpp.

driver.cpp {#drivercppdesc}
----------

This sample begins by reading an MPS file. The default MPS file is
`COIN/Mps/Sample/p0033.mps`; this can be over-riden by a command-line
specification of a (path and) file name). A second command-line argument
can specify that either the \"primal\" or \"dual\" method (or even the
\"barrier\", see below) should be used by CLP.

Once the problem has been read, there are two options for how to solve
it, one of which must be chosen at compile-time (`STYLE1` being defined
or not determines this choice). The second manner is more flexible and
involves more specific directions being given to CLP, including the
ability to specify that the barrier method should be used.

At this point in the sample, the problem is solved by CLP, and some
basic ouput is generated. If more output is desired, at compile-time, an
`exit(0)` statement must either be removed or commented. There are two
levels of additional output, the first of which is suppressed by a
`#if 0` directive which may be modified at compile-time if desired. This
first level of output only involves non-zero columns, whereas the second
provides additional information.

network.cpp {#networkcppdesc}
-----------

This handy sample reads a network problem generated by
[netgen](http://www.netlib.org/lp/generators/netgen), converts it to an
LP using CLP\'s network matrix type, and solves. This entirely avoids
the use of an MPS file, as the LP is built in memory from the network
data file created by netgen. Also, the factorization frequency is
changed, and the problem is solved more than once (demonstrating the
change of optimization sense as well as switching from dual to primal
methods).

testBarrier.cpp {#testbarriercppdesc}
---------------

This straightfoward sample begins by reading a problem from an MPS file.
It then chooses a Cholesky factorization and solves the problem using
the predictor corrector barrier method. It then copies the problem and
performs a crossover to a simplex solution in the new copy.

dualCuts.cpp
------------

This sample begins with only the equality constraints of a problem. The
inequalities are considered to be part of a pool of available cuts in
much the same way as is done in integer programming. However, in this
case, the cuts are not \"generated\", they are simply the inequalities
of the problem.

decompose.cpp
-------------

More on this sample coming soon!

driver2.cpp
-----------

More on this sample coming soon!

Common CLP Tasks in the Samples
===============================

Below is a listing of a number of common CLP tasks, such as loading a
problem from an MPS file, matched with a list of each Sample file which
illustrates the performance of a given task.

  CLP Task(s)                                 Method(s)                                                                                                                                                                                                                                                                                           Sample(s)
  ------------------------------------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- ------------------------------------
  Read problem from MPS file                  `int readMps(const char *filename)`                                                                                                                                                                                                                                                                 DEFAULTSCPP, DRIVERCPP, MINIMUMCPP
  Solve by primal method                      `int primal()`                                                                                                                                                                                                                                                                                      DRIVERCPP
  Choose pivot rule                           `void setPrimalColumnPivotAlgorithm(ClpPrimalColumnPivot &choice)` `void setDualRowPivotAlgorithm(ClpDualRowPivot &choice)`                                                                                                                                                                         DEFAULTSCPP
  Get/set infeasibility cost                  `void setInfeasibilityCost(double value)` `void setInfeasibilityCost(double value)`                                                                                                                                                                                                                 DEFAULTSCPP
  Get string/\"double\"/integer information   `bool getStrParam(ClpStrParam key, std::string &value) const` `bool getDblParam(ClpDblParam key, double &value) const` `bool  getIntParam (ClpIntParam key, int &value) const `                                                                                                                     DEFAULTSCPP
  Set maximum number of iterations            `void setMaximumIterations(int value)`                                                                                                                                                                                                                                                              DEFAULTSCPP
  Check solution status                       `int status() const` `bool isAbandoned() const` `bool isProvenOptimal() const` `bool isProvenPrimalInfeasible() const` `bool isProvenDualInfeasible() const` `bool isPrimalObjectiveLimitReached() const` `bool isDualObjectiveLimitReached() const` `bool isIterationLimitReached() const` `` ``   
                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                  

  : Contents of the Samples directory
