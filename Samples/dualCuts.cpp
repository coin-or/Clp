// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include <iomanip>

int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  // Keep names
  if (argc<2) {
    status=model.readMps("small.mps",true);
  } else {
    status=model.readMps(argv[1],true);
  }
  if (status)
    exit(10);
  /*
    This driver implements a method of treating a problem as all cuts.
    So it adds in all E rows, solves and then adds in violated rows.
  */

  double time1 = CoinCpuTime();
  
  
  int numberColumns = model.numberColumns();
  int originalNumberRows = model.numberRows();
  
  // We will need arrays to choose rows to add
  double * weight = new double [originalNumberRows];
  int * sort = new int [originalNumberRows];
  int numberSort=0;
  
  const double * rowLower = model.rowLower();
  const double * rowUpper = model.rowUpper();
  int iRow,iColumn;
  // Set up initial list
  numberSort=0;
  for (iRow=0;iRow<originalNumberRows;iRow++) {
    weight[iRow]=1.123e50;
    if (rowLower[iRow]==rowUpper[iRow]) {
      sort[numberSort++] = iRow;
      weight[iRow]=0.0;
    }
  }
  // Just add this number of rows each time in small problem
  int smallNumberRows = 2*numberColumns;
  smallNumberRows=min(smallNumberRows,originalNumberRows/20);
  // and pad out with random rows
  double ratio = ((double)(smallNumberRows-numberSort))/((double) originalNumberRows);
  for (iRow=0;iRow<originalNumberRows;iRow++) {
    if (weight[iRow]==1.123e50&&CoinDrand48()<ratio)
      sort[numberSort++] = iRow;
  }
  /* This is optional.
     The best thing to do is to miss out random rows and do a set which makes dual feasible.
     If that is not possible then make sure variables have bounds.

     One way that normally works is to automatically tighten bounds.
  */
#if 0
  model.tightenPrimalBounds();
  // However for some we need to do anyway
  double * columnLower = model.columnLower();
  double * columnUpper = model.columnUpper();
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    columnLower[iColumn]=max(-1.0e6,columnLower[iColumn]);
    columnUpper[iColumn]=min(1.0e6,columnUpper[iColumn]);
  }
#endif
  printf("%d rows in initial problem\n",numberSort);
  double * fullSolution = model.primalRowSolution();
  
  // Just do this number of passes
  int maxPass=100;
  // And take out slack rows until this pass
  int takeOutPass=20;
  int iPass;
  
  // We will be using all columns
  int * whichColumns = new int [numberColumns];
  for (int iColumn=0;iColumn<numberColumns;iColumn++)
    whichColumns[iColumn]=iColumn;
  
  for (iPass=0;iPass<maxPass;iPass++) {
    printf("Start of pass %d\n",iPass);
    // Cleaner this way
    std::sort(sort,sort+numberSort);
    // Create small problem
    ClpSimplex small(&model,numberSort,sort,numberColumns,whichColumns);
    small.setFactorizationFrequency(100+numberSort/100);
    //small.setPerturbation(50);
    //small.setLogLevel(63);
    // A variation is to just do N iterations
    // Solve 
    small.dual();
    // move solution back
    memcpy(model.primalColumnSolution(),small.primalColumnSolution(),numberColumns*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      model.setColumnStatus(iColumn,small.getColumnStatus(iColumn));
    
    for (iRow=0;iRow<numberSort;iRow++) {
      int kRow = sort[iRow];
      model.setRowStatus(kRow,small.getRowStatus(iRow));
    }
    // compute full solution
    memset(fullSolution,0,originalNumberRows*sizeof(double));
    model.times(1.0,model.primalColumnSolution(),fullSolution);
    if (iPass!=maxPass-1) {
      // Mark row as not looked at
      for (iRow=0;iRow<originalNumberRows;iRow++)
	weight[iRow]=1.123e50;
      // Look at rows already in small problem
      int iSort;
      int numberDropped=0;
      int numberKept=0;
      for (iSort=0;iSort<numberSort;iSort++) {
	iRow=sort[iSort];
	//printf("%d %g %g\n",iRow,fullSolution[iRow],small.primalRowSolution()[iSort]);
	if (model.getRowStatus(iRow)==ClpSimplex::basic) {
	  // Basic - we can get rid of if early on
	  if (iPass<takeOutPass) {
	    weight[iRow]=1.0;
	    numberDropped++;
	  } else {
	    // keep
	    weight[iRow]=-1.0e40;
	    numberKept++;
	  }
	} else {
	  // keep
	  weight[iRow]=-1.0e50;
	  numberKept++;
	}
      }
      int numberInfeasibilities=0;
      double sumInfeasibilities=0.0;
      // Now rest
      for (iRow=0;iRow<originalNumberRows;iRow++) {
	sort[iRow]=iRow;
	if (weight[iRow]==1.123e50) {
	  // not looked at yet
	  double infeasibility = max(fullSolution[iRow]-rowUpper[iRow],
				     rowLower[iRow]-fullSolution[iRow]);
	  weight[iRow]=-infeasibility;
	  if (infeasibility>1.0e-8) {
	    numberInfeasibilities++;
	    sumInfeasibilities += infeasibility;
	  }
	}
      }
      // sort
      CoinSort_2(weight,weight+originalNumberRows,sort);
      numberSort = min(originalNumberRows,smallNumberRows+numberKept);
      printf("%d rows kept, %d rows dropped - new size %d rows\n",
	     numberKept,numberDropped,numberSort);
      printf("%d rows are infeasible - sum is %g\n",numberInfeasibilities,
	       sumInfeasibilities);
      if (!numberInfeasibilities) {
	printf("Exiting as looks optimal\n");
	break;
      }
      numberInfeasibilities=0;
      sumInfeasibilities=0.0;
      for (iSort=0;iSort<numberSort;iSort++) {
	if (weight[iSort]>-1.0e30&&weight[iSort]<-1.0e-8) {
	  numberInfeasibilities++;
	  sumInfeasibilities += -weight[iSort];
	}
      }
      printf("in small model %d rows are infeasible - sum is %g\n",numberInfeasibilities,
	       sumInfeasibilities);
    }
  }
  delete [] weight;
  delete [] sort;
  delete [] whichColumns;
  // If problem is big you may wish to skip this
  model.dual();
  printf("Solve took %g seconds\n",CoinCpuTime()-time1);
  return 0;
}    
