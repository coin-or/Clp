/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "CoinStructuredModel.hpp"
#include <iomanip>

int main (int argc, const char *argv[])
{
  /* Create a structured model by reading mps file and trying
     Dantzig-Wolfe decomposition (that's the 1 parameter)
  */
  // At present D-W rows are hard coded - will move stuff from OSL
  CoinStructuredModel model((argc<2) ? "../../Data/Netlib/czprob.mps"
			: argv[1],1);
  if (!model.numberRows())
    exit(10);
  // Get default solver - could change stuff
  ClpSimplex solver;
  /*
    This driver does a simple Dantzig Wolfe decomposition
  */
  solver.solve(&model);
  // Double check
  solver.primal(1);
  return 0;
}    
