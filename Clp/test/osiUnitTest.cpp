// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#include "OsiConfig.h"

#include <cstdio>
#include <iostream>

#include "OsiUnitTests.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#include "CoinError.hpp"

#include "OsiClpSolverInterface.hpp"

using namespace OsiUnitTest;

//----------------------------------------------------------------
// to see parameter list, call unitTest -usage
//----------------------------------------------------------------

int main (int argc, const char *argv[])
{
  bool exception = false;
/*
  Start off with various bits of initialisation that don't really belong
  anywhere else.
*/
  outcomes.clear();
/*
  Synchronise C++ stream i/o with C stdio. This makes debugging
  output a bit more comprehensible. It still suffers from interleave of cout
  (stdout) and cerr (stderr), but -nobuf deals with that.
*/
  std::ios::sync_with_stdio() ;
/*
  Suppress an popup window that Windows shows in response to a crash. See
  note at head of file.
*/
  WindowsErrorPopupBlocker();

/*
  Process command line parameters.
*/
  std::map<std::string,std::string> parms;
  if (processParameters(argc,argv,parms) == false)
  { return (1) ; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;

  try {
    /*
      Test Osi{Row,Col}Cut routines.
    */
    {
      OsiClpSolverInterface clpSi;
      testingMessage( "Testing OsiRowCut with OsiClpSolverInterface\n" );
      OsiRowCutUnitTest(&clpSi,mpsDir);
    }
    {
      OsiClpSolverInterface clpSi;
      testingMessage( "Testing OsiColCut with OsiClpSolverInterface\n" );
      OsiColCutUnitTest(&clpSi,mpsDir);
    }
    {
      OsiClpSolverInterface clpSi;
      testingMessage( "Testing OsiRowCutDebugger with OsiClpSolverInterface\n" );
      OsiRowCutDebuggerUnitTest(&clpSi,mpsDir);
    }

    /*
      Run the OsiClp class test. This will also call OsiSolverInterfaceCommonUnitTest.
    */
    testingMessage( "Testing OsiClpSolverInterface\n" );
    OsiClpSolverInterfaceUnitTest(mpsDir,netlibDir);

    /*
      We have run the specialised unit test.
      Check now to see if we need to run through the Netlib problems.
    */
    if (parms.find("-testOsiSolverInterface") != parms.end())
    {
      // Create vector of solver interfaces
      std::vector<OsiSolverInterface*> vecSi(1, new OsiClpSolverInterface);

      testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
      OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

      delete vecSi[0];
    }
    else {
      testingMessage( "***Skipped Testing of OsiClpSolverInterface on Netlib problems***\n" );
      testingMessage( "***use -testOsiSolverInterface to run them.***\n" );
    }
  } catch (CoinError& error) {
    std::cout.flush();
    std::cerr << "Caught CoinError exception: ";
    error.print(true);
    exception = true;
  }

  /*
    We're done. Report on the results.
  */
  std::cout.flush();
  outcomes.print();

  int nerrors;
  int nerrors_expected;
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);

  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected << " unexpected errors." << std::endl ;
  else if( exception )
    std::cerr << "Tests completed with exception\n";
  else
    std::cerr << "All tests completed successfully\n";

  return (nerrors - nerrors_expected) + (exception ? 1 : 0);
}
