// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpFactorization.hpp"
#include "Presolve.hpp"
#include <iomanip>

// Time

#include  <time.h>
#if !defined(_MSC_VER)
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#endif
static double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow;        /* clock_t is same as int */
  
  ticksnow = (unsigned int)clock();
  
  cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}

int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  // Keep names
  if (argc<2)
    status=model.readMps("../../Mps/Sample/p0033.mps",true);
  else
    status=model.readMps(argv[1],true);
  /*
    This driver show how to do presolve
  */

  if (status) {
    fprintf(stderr,"Bad readMps %s\n",argv[1]);
    fprintf(stdout,"Bad readMps %s\n",argv[1]);
    exit(1);
  }

  double time1 = cpuTime();
  ClpSimplex * model2;
  Presolve pinfo;
  int numberPasses=5; // can change this
  model2 = pinfo.presolvedModel(model,1.0e-8,false,numberPasses);
  if (!model2) {
    fprintf(stderr,"Presolve says %s is infeasible with tolerance of %g\n",
	    argv[1],1.0e-8);
    fprintf(stdout,"Presolve says %s is infeasible with tolerance of %g\n",
	    argv[1],1.0e-8);
    model2 = pinfo.presolvedModel(model,1.0e-7,false,numberPasses);
    if (!model2) {
      fprintf(stderr,"Presolve says %s is infeasible with tolerance of %g\n",
	      argv[1],1.0e-7);
      fprintf(stdout,"Presolve says %s is infeasible with tolerance of %g\n",
	      argv[1],1.0e-7);
      exit(2);
    }
  }
  model2->checkSolution();
  // change from 200
  model2->factorization()->maximumPivots(100+model2->numberRows()/50);
  if (argc<3 ||!strstr(argv[2],"primal")) {
    // faster if bounds tightened
    int numberInfeasibilities = model2->tightenPrimalBounds();
    if (numberInfeasibilities)
      std::cout<<"** Analysis indicates model infeasible"
	       <<std::endl;
    // up dual bound for safety
    model2->setDualBound(1.0e10);
    model2->dual();
  } else {
    // up infeasibility cost for safety
    model2->setInfeasibilityCost(1.0e10);
    model2->primal();
  }
  pinfo.postsolve(true);
		
  int numberIterations=model2->numberIterations();;
  delete model2;
  printf("Resolving from postsolved model\n");
  model.primal(1);

  numberIterations += model.numberIterations();;
  // for running timing tests
  std::cout<<argv[1]<<" Objective "<<model.objectiveValue()<<" took "<<
    numberIterations<<" iterations and "<<cpuTime()-time1<<" seconds"<<std::endl;

  std::string modelName;
  model.getStrParam(ClpProbName,modelName);
  std::cout<<"Model "<<modelName<<" has "<<model.numberRows()<<" rows and "<<
    model.numberColumns()<<" columns"<<std::endl;

  // remove this to print solution

  exit(0);

  /*
    Now to print out solution.  The methods used return modifiable
    arrays while the alternative names return const pointers -
    which is of course much more virtuous.

    This version just does non-zero columns
  
   */
#if 0
  int numberRows = model.numberRows();

  // Alternatively getRowActivity()
  double * rowPrimal = model.primalRowSolution();
  // Alternatively getRowPrice()
  double * rowDual = model.dualRowSolution();
  // Alternatively getRowLower()
  double * rowLower = model.rowLower();
  // Alternatively getRowUpper()
  double * rowUpper = model.rowUpper();
  // Alternatively getRowObjCoefficients()
  double * rowObjective = model.rowObjective();
    
  // If we have not kept names (parameter to readMps) this will be 0
  assert(model.lengthNames());

  // Row names
  const std::vector<std::string> * rowNames = model.rowNames();


  int iRow;

  std::cout<<"                       Primal          Dual         Lower         Upper        (Cost)"
	   <<std::endl;

  for (iRow=0;iRow<numberRows;iRow++) {
    double value;
    std::cout<<std::setw(6)<<iRow<<" "<<std::setw(8)<<(*rowNames)[iRow];
    value = rowPrimal[iRow];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = rowDual[iRow];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = rowLower[iRow];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = rowUpper[iRow];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    if (rowObjective) {
      value = rowObjective[iRow];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    }
    std::cout<<std::endl;
  }
#endif
  std::cout<<"--------------------------------------"<<std::endl;

  // Columns

  int numberColumns = model.numberColumns();

  // Alternatively getColSolution()
  double * columnPrimal = model.primalColumnSolution();
  // Alternatively getReducedCost()
  double * columnDual = model.dualColumnSolution();
  // Alternatively getColLower()
  double * columnLower = model.columnLower();
  // Alternatively getColUpper()
  double * columnUpper = model.columnUpper();
  // Alternatively getObjCoefficients()
  double * columnObjective = model.objective();
    
  // If we have not kept names (parameter to readMps) this will be 0
  assert(model.lengthNames());

  // Column names
  const std::vector<std::string> * columnNames = model.columnNames();


  int iColumn;
  
  std::cout<<"                       Primal          Dual         Lower         Upper          Cost"
	   <<std::endl;

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double value;
    value = columnPrimal[iColumn];
    if (fabs(value)>1.0e-8) {
      std::cout<<std::setw(6)<<iColumn<<" "<<std::setw(8)<<(*columnNames)[iColumn];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
      value = columnDual[iColumn];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
      value = columnLower[iColumn];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
      value = columnUpper[iColumn];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
      value = columnObjective[iColumn];
      if (fabs(value)<1.0e5)
	std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
      else
	std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
      
      std::cout<<std::endl;
    }
  }
  std::cout<<"--------------------------------------"<<std::endl;

  return 0;
}    
