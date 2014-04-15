// $Id: driver3.cpp 1898 2013-04-09 18:06:04Z stefan $
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip> 

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#ifdef ABC_INHERIT
#include "AbcSimplex.hpp"
#endif
#include "CoinTime.hpp"

//#############################################################################


/************************************************************************

This main program shows how to take advantage of the standalone clp in your program.
It should perform very nearly the same as clp  
First it reads in a model from an mps file
Then it calls ClpMain1 passing all parameters apart from first
Finally it prints solution

************************************************************************/

#ifndef ABC_INHERIT
void ClpMain0(ClpSimplex * models);
int ClpMain1(int argc, const char *argv[],ClpSimplex * model);
#else
void ClpMain0(AbcSimplex * models);
int ClpMain1(int argc, const char *argv[],AbcSimplex * model);
#endif
int main (int argc, const char *argv[])
{

#ifndef ABC_INHERIT
  ClpSimplex model;
#else
  AbcSimplex model;
#endif
  if (argc > 1) {
    printf("command line - ");
    for (int i = 0; i < argc; i++)
      printf("%s ", argv[i]);
    printf("\n");
  }
  // change defaults to match standalone solver
  ClpMain0(&model);
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
#endif
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = model.readMps(mpsFileName.c_str());
  if( numMpsReadErrors != 0 )
  {
     printf("%d errors reading MPS file\n", numMpsReadErrors);
     return numMpsReadErrors;
  }
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc>2) {
    ClpMain1(argc-1,argv+1,&model);
  } else {
    const char * argv2[]={"driver3","-solve","-quit"};
    ClpMain1(3,argv2,&model);
  }

  // Print solution
    
  const double * solution = model.primalColumnSolution();
  int numberColumns = model.numberColumns();
  //const double * reducedCosts = model.dualColumnSolution();
    
  std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
  std::cout<<"--------------------------------------"<<std::endl;
    
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    double value=solution[iColumn];
    if (fabs(value)>1.0e-7) 
      std::cout<<std::setw(6)<<iColumn<<" "<<std::setw(8)<<setiosflags(std::ios::left)<<model.getColumnName(iColumn)
	       <<resetiosflags(std::ios::adjustfield)<<std::setw(14)<<" "<<value<<std::endl;
  }
  std::cout<<"--------------------------------------"<<std::endl;
  
  std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  return 0;
}    
