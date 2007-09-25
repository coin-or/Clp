// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "ClpPresolve.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpNetworkMatrix.hpp"
#include "CoinTime.hpp"
#include <iomanip>
int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  // Keep names when reading an mps file
  if (argc<2)
    status=model.readMps("../../Data/Sample/p0033.mps",true);
  else
    status=model.readMps(argv[1],true);
  
  if (status) {
    fprintf(stderr,"Bad readMps %s\n",argv[1]);
    fprintf(stdout,"Bad readMps %s\n",argv[1]);
    exit(1);
  }
  ClpSimplex * model2;
  ClpPresolve pinfo;
  int numberPasses=5; // can change this
  /* Use a tolerance of 1.0e-8 for feasibility, treat problem as 
     not being integer, do "numberpasses" passes and throw away names
     in presolved model */
#define PRESOLVE
#ifdef PRESOLVE
  model2 = pinfo.presolvedModel(model,1.0e-8,false,numberPasses,false);
#else
  model2 = &model;
#endif
  if (!model2) {
    fprintf(stderr,"ClpPresolve says %s is infeasible with tolerance of %g\n",
	    argv[1],1.0e-8);
    fprintf(stdout,"ClpPresolve says %s is infeasible with tolerance of %g\n",
	    argv[1],1.0e-8);
    // model was infeasible - maybe try again with looser tolerances
    model2 = pinfo.presolvedModel(model,1.0e-7,false,numberPasses,false);
    if (!model2) {
      fprintf(stderr,"ClpPresolve says %s is infeasible with tolerance of %g\n",
	      argv[1],1.0e-7);
      fprintf(stdout,"ClpPresolve says %s is infeasible with tolerance of %g\n",
	      argv[1],1.0e-7);
      exit(2);
    }
  }
  // change factorization frequency from 200
  model2->setFactorizationFrequency(100+model2->numberRows()/50);
  model2->setPerturbation(50);
  double time1 = CoinCpuTime();
  int numberRows = model2->numberRows();
  int numberColumns = model2->numberColumns();
  char * type = new char[numberRows];
  if (false) {
    // scale rows
    const int * row = model2->matrix()->getIndices();
    const CoinBigIndex * columnStart = model2->matrix()->getVectorStarts();
    const int * columnLength = model2->matrix()->getVectorLengths(); 
    double * elementByColumn = model2->matrix()->getMutableElements();
    double * largest = new double [numberRows];
    double * smallest = new double [numberRows];
    int iRow,iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow = row[j];
	double value = fabs(elementByColumn[j]);
	largest[iRow] = CoinMax(largest[iRow],value);
	smallest[iRow] = CoinMin(smallest[iRow],value);
      }
    }
    int nScale=0;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (largest[iRow]==smallest[iRow]&&largest[iRow]!=1.0) {
	nScale++;
	largest[iRow]= 1.0/largest[iRow];
      } else {
	largest[iRow]=0.0;
      }
    }
    if (nScale) {
      printf("Scaling %d rows\n",nScale);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  iRow = row[j];
	  double value = largest[iRow];
	  if (value)
	    elementByColumn[j] *= value;
	  assert (fabs(elementByColumn[j])==1.0);
	}
      }  
      printf("would have to do rhs as well\n");
      abort();
    }
    delete [] largest;
    delete [] smallest;
  }
  int nFind=model2->findNetwork(type,0.001);
  printf("%d network rows found in %g seconds\n",nFind,CoinCpuTime()-time1);
  if (CoinAbs(nFind)>0.2*numberRows) {
    const int * row = model2->matrix()->getIndices();
    const CoinBigIndex * columnStart = model2->matrix()->getVectorStarts();
    const int * columnLength = model2->matrix()->getVectorLengths(); 
    const double * elementByColumn = model2->matrix()->getElements();
    int * whichRow = new int [numberRows];
    int * whichColumn = new int [numberColumns];
    int iRow,iColumn;
    double * newElement = new double[numberColumns];
    int * newColumn = new int[numberColumns];
    int nEl=0;
    int nSmallRows=0;
    int nSmallColumns=0;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (type[iRow]>=0)
	whichRow[nSmallRows++]=iRow;
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double product =1.0;
      int n=0;
      for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow = row[j];
	if (type[iRow]>=0) {
	  n++;
	  assert (fabs(elementByColumn[j])==1.0);
	  product *= elementByColumn[j];
	  if (type[iRow])
	    product = - product;
	}
      }
      assert (n<=2);
      if (n==2) {
	assert (product==-1.0);
      } else if (n==1) {
	newColumn[nEl]=nSmallColumns;
	newElement[nEl++] = -product;
      }
      if (n)
	whichColumn[nSmallColumns++]=iColumn;
    }
    // Create small problem
    ClpSimplex small(model2,nSmallRows,whichRow,nSmallColumns,whichColumn);
    printf("Small model has %d columns\n",nSmallColumns);
    // change elements and rhs
    double * lower = small.rowLower();
    double * upper = small.rowUpper();
    for (iRow=0;iRow<nSmallRows;iRow++) {
      if (type[whichRow[iRow]]) {
	double temp = -lower[iRow];
	lower[iRow]=-upper[iRow];
	upper[iRow]=temp;
      }
    }
    row = small.matrix()->getIndices();
    columnStart = small.matrix()->getVectorStarts();
    columnLength = small.matrix()->getVectorLengths(); 
    double * element = small.matrix()->getMutableElements();
    for (iColumn=0;iColumn<nSmallColumns;iColumn++) {
      for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow = row[j];
	if (type[whichRow[iRow]]) 
	  element[j] = - element[j];
      }
    }
    if (nEl) {
      // add row (probably should fix network)
      small.addRow(nEl,newColumn,newElement);
    }
    delete [] newElement;
    delete [] newColumn;
    ClpNetworkMatrix * network = new ClpNetworkMatrix(*small.matrix());
    small.replaceMatrix(network,true);
    small.initialSolve();
    double * fullSolution = model2->primalColumnSolution();
    // move solution back
    const double * solution = small.primalColumnSolution();
    for (iColumn=0;iColumn<nSmallColumns;iColumn++) {
      int kColumn = whichColumn[iColumn];
      model2->setColumnStatus(kColumn,small.getColumnStatus(iColumn));
      fullSolution[kColumn]=solution[iColumn];
    }
    for (iRow=0;iRow<nSmallRows;iRow++) {
      ClpSimplex::Status status=small.getRowStatus(iRow);
      if (type[whichRow[iRow]]) {
	if (status==ClpSimplex::atLowerBound)
	  status=ClpSimplex::atUpperBound;
	else if (status==ClpSimplex::atUpperBound)
	  status=ClpSimplex::atLowerBound;
      }
      model2->setRowStatus(whichRow[iRow],status);
    }
    double * rowSolution = model2->primalRowSolution();
    memset (rowSolution,0,numberRows*sizeof(double));
    model2->times(1.0,fullSolution,rowSolution);
    delete [] whichRow;
    delete [] whichColumn;
    model2->dual();
  } else {
    model2->initialSolve();
  }
  delete [] type;
  int numberIterations=model2->numberIterations();;
#ifdef PRESOLVE
  pinfo.postsolve(true);
  delete model2;
#endif
  /* After this postsolve model should be optimal.
     We can use checkSolution and test feasibility */
  model.checkSolution();
  if (model.numberDualInfeasibilities()||
      model.numberPrimalInfeasibilities()) 
    printf("%g dual %g(%d) Primal %g(%d)\n",
	   model.objectiveValue(),
	   model.sumDualInfeasibilities(),
	   model.numberDualInfeasibilities(),
	   model.sumPrimalInfeasibilities(),
	   model.numberPrimalInfeasibilities());
  // But resolve for safety
  model.primal(1);

  numberIterations += model.numberIterations();;
  // for running timing tests
  std::cout<<argv[1]<<" Objective "<<model.objectiveValue()<<" took "<<
    numberIterations<<" iterations and "<<
    CoinCpuTime()-time1<<" seconds"<<std::endl;
  return 0;
}    
