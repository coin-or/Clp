// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"
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
    This driver implements what I called Sprint.  Cplex calls it 
    "sifting" which is just as silly.  When I thought of this trivial idea
    it reminded me of an LP code of the 60's called sprint which after
    every factorization took a subset of the matrix into memory (all
    64K words!) and then iterated very fast on that subset.  On the
    problems of those days it did not work very well, but it worked very
    well on aircrew scheduling problems where there were very large numbers
    of columns all with the same flavor.
  */
  // Data in large model
  int numberColumns = model.numberColumns();
  int numberRows = model.numberRows();
  const double * columnLower = model.columnLower();
  const double * columnUpper = model.columnUpper();
  double * fullSolution = model.primalColumnSolution();
  double * rowSolution = model.primalRowSolution();
  const double * objective = model.objective();
  double originalOffset;
  model.getDblParam(ClpObjOffset,originalOffset);

  // We will need arrays to choose variables.
  double * weight = new double [numberColumns];
  int * sort = new int [numberColumns];
  // Just take this number of columns in small problem
  int smallNumberColumns = min(3*numberRows,numberColumns);
  smallNumberColumns = max(smallNumberColumns,3000);
  // Set up initial list 0,1,.... (assuming that has feasible solution)
  int numberSort=numberRows;
  CoinIotaN(sort, numberRows, 0);

  // We will be using all rows
  int * whichRows = new int [numberRows];
  CoinIotaN(whichRows, numberRows, 0);

  // Just do this number of passes
  int maxPass=100;
  double lastObjective=1.0e31;

  for (int iPass=0;iPass<maxPass;iPass++) {
    /* Create small problem.
       This constructor not in OsiSolverInterface - might be useful */
    ClpSimplex small(&model,numberRows,whichRows,numberSort,sort);
    // now see what variables left out do to row solution
    CoinZeroN (rowSolution,numberRows);
    int iRow,iColumn;
    // zero out ones in small problem
    for (iColumn=0;iColumn<numberSort;iColumn++) {
      int kColumn = sort[iColumn];
      fullSolution[kColumn]=0.0;
    }
    // Get objective offset for ones not in small problem
    double offset=0.0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      offset += fullSolution[iColumn]*objective[iColumn];
    small.setDblParam(ClpObjOffset,originalOffset-offset);

    // Modify RHS for ones not in small problem
    model.times(1.0,fullSolution,rowSolution);
    double * lower = small.rowLower();
    double * upper = small.rowUpper();
    for (iRow=0;iRow<numberRows;iRow++) {
      if (lower[iRow]>-1.0e50) 
	lower[iRow] -= rowSolution[iRow];
      if (upper[iRow]<1.0e50)
	upper[iRow] -= rowSolution[iRow];
    }
    // Solve 
    small.primal();
    // move solution back
    const double * smallSolution = small.primalColumnSolution();
    const double * smallRowSolution = small.primalRowSolution();
    for (iColumn=0;iColumn<numberSort;iColumn++) {
      int kColumn = sort[iColumn];
      model.setColumnStatus(kColumn,small.getColumnStatus(iColumn));
      fullSolution[kColumn]=smallSolution[iColumn];
    }
    for (iRow=0;iRow<numberRows;iRow++) {
      model.setRowStatus(iRow,small.getRowStatus(iRow));
      rowSolution[iRow] = smallRowSolution[iRow];
    }
    // Have we finished
    if ((small.objectiveValue()>lastObjective-1.0e-7&&iPass>5)||
	!small.numberIterations()||iPass==maxPass-1) {
      break; // finished
    } else {
      lastObjective = small.objectiveValue();
      // get reduced costs for large problem
      // this assumes minimization
      CoinMemcpyN(model.objective(),numberColumns,weight);
      model.transposeTimes(-1.0,small.dualRowSolution(),weight);
      // now massage weight so all basic in plus good djs
      int negativeDjs = 0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double dj = weight[iColumn];
	double value = fullSolution[iColumn];
	if (model.getColumnStatus(iColumn)==ClpSimplex::basic) {
	  dj = -1.0e50;
	} else if (dj<0.0&&value<columnUpper[iColumn]) {
	  dj = dj;
	  negativeDjs++;
	} else if (dj>0.0&&value>columnLower[iColumn]) {
	  dj = -dj;
	  negativeDjs++;
	} else if (columnUpper[iColumn]>columnLower[iColumn]) {
	  dj = fabs(dj);
	} else {
	  dj = 1.0e50;
	}
	weight[iColumn] = dj;
	sort[iColumn] = iColumn;
      }
      std::cout<<"After pass "<<iPass+1<<" there are "<<
	negativeDjs<<" negative reduced costs"<<std::endl;
      // sort (could do faster as exact order not important)
      CoinSort_2(weight,weight+numberColumns,sort);
      numberSort = smallNumberColumns;
    }
  }
  delete [] weight;
  delete [] sort;
  delete [] whichRows;
  // Just to show optimal - might not do in practice
  model.primal(1);
  return 0;
}    
