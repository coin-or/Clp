// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpInterior.hpp"
#include "ClpSimplex.hpp"
//#include "ClpCholeskyWssmp.hpp"
int main (int argc, const char *argv[])
{
  ClpInterior  model;
  int status;
  if (argc<2)
    status=model.readMps("../../Mps/Sample/p0033.mps");
  else
    status=model.readMps(argv[1]);
  if (status) {
    printf("errors on input\n");
    exit(77);
  }
  // ** note this does not have presolve
  //ClpCholeskyWssmp * cholesky = new ClpCholeskyWssmp();
  //model.setCholesky(cholesky);
  model.primalDual();
  // Do crossover
  ClpSimplex model2(model);
  // make sure no status left
  model2.createStatus();
  model2.primal(1);
  return 0;
}    
