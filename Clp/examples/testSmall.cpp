// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "ClpSmallMatrix.hpp"
#include "CoinMpsIO.hpp"
#include "CoinTime.hpp"

int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  // Keep names
  if (argc<2) {
    status=model.readMps("../../Data/Netlib/czprob.mps",true);
  } else {
    status=model.readMps(argv[1],true);
  }
  if (status)
    exit(10);
  /*
    This driver checks small matrix works
  */
  int numberRows = model.numberRows();
  assert (numberRows<65536);
  // get small matrix
  ClpSmallMatrix * smallMatrix = new ClpSmallMatrix(*model.matrix());
  // replace and delete original
  model.setSpecialOptions(256); // to say no row copy
  model.replaceMatrix(smallMatrix,true);
  double time1 = CoinCpuTime();
  //model.setLogLevel(63);
  model.dual();
  printf("Dual took %g seconds\n",CoinCpuTime()-time1);
  model.dual();
  model.allSlackBasis();
  model.primal();
  model.primal(1);
  memset(model.statusArray(),0,model.numberColumns());
  model.primal(1);
  return 0;
}    
