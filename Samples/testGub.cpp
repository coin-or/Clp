// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "ClpGubMatrix.hpp"
#include "ClpPrimalColumnSteepest.hpp"
int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  int maxIts=200;
  if (argc<2)
    status=model.readMps("../../Mps/Sample/p0033.mps");
  else
    status=model.readMps(argv[1]);
  if (status) {
    printf("errors on input\n");
    exit(77);
  }
  if (argc>2) {
    maxIts = atoi(argv[2]);
    printf("max its %d\n",maxIts);
  }
  // For now scaling off
  model.scaling(0);
  // Do partial dantzig
  ClpPrimalColumnSteepest dantzig(5);
  model.setPrimalColumnPivotAlgorithm(dantzig);
  model.messageHandler()->setLogLevel(63);
  //model.setFactorizationFrequency(1);
  model.setMaximumIterations(maxIts);
  model.primal();
  // find gub
  int numberRows = model.numberRows();
  int * gubStart = new int[numberRows];
  int * gubEnd = new int[numberRows];
  int * which = new int[numberRows];
  int * whichGub = new int[numberRows];
  double * lower = new double[numberRows];
  double * upper = new double[numberRows];
  const double * rowLower = model.rowLower();
  const double * rowUpper = model.rowUpper();
  int numberColumns = model.numberColumns();
  int * mark = new int[numberColumns];
  int iRow,iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++)
    mark[iColumn]=-1;
  CoinPackedMatrix * matrix = model.matrix();
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int * column = rowCopy.getIndices();
  const int * rowLength = rowCopy.getVectorLengths();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const double * element = rowCopy.getElements();
  int putGub=numberRows;
  int putNonGub=numberRows;
  for (iRow=numberRows-1;iRow>=0;iRow--) {
    bool gubRow=true;
    int first=numberColumns+1;
    int last=-1;
    for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      if (element[j]!=1.0) {
	gubRow=false;
	break;
      } else {
	int iColumn = column[j];
	if (mark[iColumn]>=0) {
	  gubRow=false;
	  break;
	} else {
	  last=max(last,iColumn);
	  first = min(first,iColumn);
	}
      }
    }
    if (last-first+1!=rowLength[iRow]||!gubRow) {
      which[--putNonGub]=iRow;
    } else {
      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	mark[iColumn]=iRow;
      }
      putGub--;
      gubStart[putGub]=first;
      gubEnd[putGub]=last+1;
      lower[putGub]=rowLower[iRow];
      upper[putGub]=rowUpper[iRow];
      whichGub[putGub]=iRow;
    }
  }
  int numberNonGub=numberRows-putNonGub;
  int numberGub = numberRows-putGub;
  if (numberGub>0) {
    printf("** %d gub rows\n",numberGub);
    for (iColumn=0;iColumn<numberColumns;iColumn++)
      mark[iColumn]=iColumn;
    ClpSimplex model2(&model,numberNonGub,which+putNonGub,numberColumns,mark);
    ClpMatrixBase * saveMatrix = model2.clpMatrix();
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
    model2.replaceMatrix(new ClpGubMatrix(clpMatrix,numberGub,
					 gubStart+putGub,gubEnd+putGub,
					 lower+putGub,upper+putGub));
#if 1
    saveMatrix = model2.clpMatrix();
    ClpGubMatrix* gubMatrix =
      dynamic_cast< ClpGubMatrix*>(saveMatrix);
    assert (gubMatrix);
    // deal with basis for gub
    for (int iSet=0;iSet<numberGub;iSet++) {
      int iRow = whichGub[putGub+iSet];
      gubMatrix->setStatus(iSet,model.getRowStatus(iRow));
    }
#endif
    // For now scaling off
    model2.scaling(0);
    // Do partial dantzig
    ClpPrimalColumnSteepest dantzig(5);
    model2.setPrimalColumnPivotAlgorithm(dantzig);
    model2.messageHandler()->setLogLevel(63);
    model2.setFactorizationFrequency(1);
    model2.setMaximumIterations(200);
    model2.primal();
  } else {
    // dummy gub
    ClpMatrixBase * saveMatrix = model.clpMatrix();
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
    model.replaceMatrix(new ClpGubMatrix(clpMatrix,0,NULL,NULL,NULL,NULL)); 
    // For now scaling off
    model.scaling(0);
    ClpPrimalColumnSteepest dantzig(5);
    model.setPrimalColumnPivotAlgorithm(dantzig);
    model.primal();
  }
  delete [] mark;
  delete [] gubStart;
  delete [] gubEnd;
  delete [] which;
  delete [] whichGub;
  delete [] lower;
  delete [] upper;
  return 0;
}    
