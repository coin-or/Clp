// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* This simple example takes a matrix read in by CoinMpsIo,
   deletes every second column and solves the resulting problem */

#include "ClpSimplex.hpp"
#include "CoinMpsIO.hpp"
#include <iomanip>


int main (int argc, const char *argv[])
{
  int status;
  CoinMpsIO m;
  if (argc<2)
    status=m.readMps("../../Mps/Sample/p0033.mps","");
  else
    status=m.readMps(argv[1],"");

  if (status) {
    fprintf(stdout,"Bad readMps %s\n",argv[1]);
    exit(1);
  }
  // Get data arrays
  const CoinPackedMatrix * matrix1 = m.getMatrixByCol();
  const int * start1 = matrix1->getVectorStarts();
  const int * length1 = matrix1->getVectorLengths();
  const int * row1 = matrix1->getIndices();
  const double * element1 = matrix1->getElements();

  const double * columnLower1 = m.getColLower();
  const double * columnUpper1 = m.getColUpper();
  const double * rowLower1 = m.getRowLower();
  const double * rowUpper1 = m.getRowUpper();
  const double * objective1 = m.getObjCoefficients();

  int numberColumns = m.getNumCols();
  int numberRows = m.getNumRows();
  int numberElements = m.getNumElements();

  // Get new arrays
  int numberColumns2 = (numberColumns+1)>>1;
  int * start2 = new int[numberColumns2+1];
  int * row2 = new int[numberElements];
  double * element2 = new double[numberElements];

  double * objective2 = new double[numberColumns2];
  double * columnLower2 = new double[numberColumns2];
  double * columnUpper2 = new double[numberColumns2];
  double * rowLower2 = new double[numberRows];
  double * rowUpper2 = new double[numberRows];

  // could have used old arrays
  memcpy(rowLower2,rowLower1,numberRows*sizeof(double));
  memcpy(rowUpper2,rowUpper1,numberRows*sizeof(double));

  int iColumn;
  numberColumns2=0;
  numberElements=0;
  start2[0]=0;
  for (iColumn=0;iColumn<numberColumns;iColumn+=2) {
    columnLower2[numberColumns2] = columnLower1[iColumn];
    columnUpper2[numberColumns2] = columnUpper1[iColumn];
    objective2[numberColumns2] = objective1[iColumn];
    int j;
    for (j=start1[iColumn];j<start1[iColumn]+length1[iColumn];j++) {
      row2[numberElements]=row1[j];
      element2[numberElements++] = element1[j];
    }
    start2[++numberColumns2]=numberElements;
  }

  ClpSimplex  model;

  // load up
  model.loadProblem(numberColumns2,numberRows,
		    start2,row2,element2,
		    columnLower2,columnUpper2,
		    objective2,
		    rowLower2,rowUpper2);
  // delete
  delete [] start2;
  delete [] row2 ;
  delete [] element2;

  delete [] objective2;
  delete [] columnLower2;
  delete [] columnUpper2;
  delete [] rowLower2;
  delete [] rowUpper2;

  // solve
  model.primal();
  return 0;
}    
