// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/* When playing around, it simplifies thinking if all lower bounds
   are zero and where possible coefficients are integral.  This
   helps do that. */
#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinMpsIO.hpp"
#include "CoinRational.hpp"
static long computeGcd(long a, long b) {
  // This is the standard Euclidean algorithm for gcd
  long remainder = 1;
  // Make sure a<=b (will always remain so)
  if (a > b) {
    // Swap a and b
    long temp = a;
    a = b;
    b = temp;
  }
  // If zero then gcd is nonzero
  if (!a) {
    if (b) {
      return b;
    } 
    else {
      printf("### WARNING: CglGMI::computeGcd() given two zeroes!\n");
      exit(1);
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
} /* computeGcd */
static bool scaleRowIntegral(double* rowElem, int rowNz)
{
  long gcd, lcm;
  double maxdelta = 1.0e-13;
  double maxscale = 1000; 
  long maxdnom = 1000;
  //long numerator = 0, denominator = 0;
  // Initialize gcd and lcm
  CoinRational r = CoinRational(rowElem[0], maxdelta, maxdnom);
  if (r.getNumerator() != 0){
    gcd = labs(r.getNumerator());
    lcm = r.getDenominator();
  } else {
    return false;
  } 
  for (int i = 1; i < rowNz; ++i) {
    if (rowElem[i]) {
      r = CoinRational(rowElem[i], maxdelta, maxdnom);
      if (r.getNumerator() != 0){
	gcd = computeGcd(gcd, r.getNumerator());
	lcm *= r.getDenominator()/(computeGcd(lcm,r.getDenominator()));
      } else {
	return false;
      }
    }
  }
  double scale = static_cast<double>(lcm)/static_cast<double>(gcd);
  if (fabs(scale) > maxscale) {
    return false;
  }
  scale = fabs(scale);
  // Looks like we have a good scaling factor; scale and return;
  for (int i = 0; i < rowNz; ++i) {
    double value = rowElem[i]*scale;
    rowElem[i] = floor(value+0.5);
    assert (fabs(rowElem[i]-value)<1.0e-9);
  }
  return true;
} /* scaleRowIntegral */
int main(int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find MPS file.\n");
    exit(1);
  } else {
    status = model.readMps(argv[1]);
  }
  if (status) {
    printf("errors on input\n");
    exit(77);
  }
  int cleanMode = 23;
  if (argc>2) {
    cleanMode = atoi(argv[2]);
    if (cleanMode<=0||cleanMode>31) {
      printf("Bad value for cleanMode - should be 1-31\n");
      exit(77);
    }
  }
  if ((cleanMode&1)!=0)
    printf("(1)moving a bound to zero ");
  if ((cleanMode&2)!=0)
    printf("(2)swap negative uppers ");
  if ((cleanMode&4)!=0)
    printf("(4)try scale to integer ");
  if ((cleanMode&8)!=0)
    printf("(8)majority coefficents positive");
  if ((cleanMode&16)!=0)
    printf("(16)G->L");
  printf("\n");
  model.dual();
  int nMoved = 0;
  int nSwappedColumn = 0;
  int nScaled = 0;
  int nSwappedRow = 0;
  /*
    1 - make one bound zero
    2 - swap all with negative upper bounds
    4 - try and make rows have integer coefficients
    8 - majority of coefficients positive
    16 - make G rows into L ones (overrides 8 - so 8 is for ranged or E)
   */
  int numberRows = model.numberRows();
  int numberColumns = model.numberColumns();
  double * columnLower = model.columnLower();
  double * columnUpper = model.columnUpper();
  double * columnActivity = model.primalColumnSolution();
  double offset = model.objectiveOffset();
  double * objective = model.objective();
  double * rowLower = model.rowLower();
  double * rowUpper = model.rowUpper();
  double * rowActivity = model.primalRowSolution();
  CoinPackedMatrix * matrixByColumn = model.matrix();
  const int *row = matrixByColumn->getIndices();
  const CoinBigIndex *columnStart = matrixByColumn->getVectorStarts();
  const int *columnLength = matrixByColumn->getVectorLengths();
  double *columnElements = matrixByColumn->getMutableElements();
  if ((cleanMode&2)!=0) {
    // swap columns
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (fabs(columnUpper[iColumn])<fabs(columnLower[iColumn])) {
	nSwappedColumn++;
	double value = -columnUpper[iColumn];
	columnUpper[iColumn] = -columnLower[iColumn];
	columnLower[iColumn] = value;
	objective[iColumn] = -objective[iColumn];
	for (CoinBigIndex j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++)
	  columnElements[j] = - columnElements[j];
      }
    }
  }
  if ((cleanMode&28)!=0) {
    CoinPackedMatrix matrixByRow(*matrixByColumn);
    matrixByRow.reverseOrdering();
    double *elementByRow = matrixByRow.getMutableElements();
    const int *column = matrixByRow.getIndices();
    const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
    const int *rowLength = matrixByRow.getVectorLengths();
    double * temp = new double [2*numberColumns+4+numberRows];
    double * tempSave = temp+numberColumns+2;
    double * rowScale = tempSave+numberColumns+2; 
    int nChanged = 0;
    for (int iRow=0;iRow<numberRows;iRow++) {
      int n = 0;
      // make majority positive?
      CoinBigIndex start = rowStart[iRow];
      CoinBigIndex end = start+rowLength[iRow];
      double multiplier = 1.0;
      for (CoinBigIndex j=start;j<end;j++) {
	if (elementByRow[j]<0)
	  n++;
      }
      if ((cleanMode&8)!=0 && ((cleanMode&16)==0 || rowUpper[iRow] < 1.0e30)) {
	if (2*n>end-start)
	  multiplier = -1.0;
      } else if ((cleanMode&16)!=0 && rowUpper[iRow] > 1.0e30) {
	multiplier = -1.0;
      }
      if ((cleanMode&4)!=0) {
	int nInRow = end-start;
	n = nInRow;
	for (int i=0;i<n;i++) 
	  temp[i] = multiplier*elementByRow[i+start];
	double lower = rowLower[iRow];
	double upper = rowUpper[iRow];
	if (lower>-1.0e20) 
	  temp[n++] = multiplier*lower;
	if (upper<1.0e20) 
	  temp[n++] = multiplier*upper;
	memcpy(tempSave,temp,n*sizeof(double));
	if (scaleRowIntegral(temp, n)) {
	  // double check
	  double largestError = 0.0;
	  double mult = temp[0]/elementByRow[start];
	  if (fabs(mult-floor(mult+0.1))<1.0e-12)
	    mult = floor(mult+0.1);
	  for (int i=0;i<n;i++) {
	    double value = mult*tempSave[i];
	    if (value) {
	      double vint = floor(value+0.01);
	      largestError = CoinMax(largestError,fabs(value-vint));
	      assert (fabs(vint)>0.9);
	    }
	  }
	  if (largestError<1.0e-9) {
	    multiplier = mult;
	  }
	}
	rowScale[iRow] = multiplier;
	if (multiplier!=1.0) {
	  nChanged++;
	  if (fabs(multiplier)!=1.0)
	    nScaled++;
	  //memcpy(elementByRow+start,temp,(end-start)*sizeof(double));
	  if (multiplier<0.0) {
	    double tempV = lower;
	    lower = -upper;
	    upper = -tempV;
	    nSwappedRow++;
	  }
	  if (lower>-1.0e20)
	    lower *= fabs(multiplier);
	  if (upper<1.0e20) 
	    upper *= fabs(multiplier);
	  rowLower[iRow] = lower;
	  rowUpper[iRow] = upper;
	}
      }
    }
    if (nChanged) {
      double largestDelta = 0.0;
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	for (CoinBigIndex j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  double value = columnElements[j];
	  if (fabs(rowScale[iRow])!=1.0) {
	    value *= rowScale[iRow];
	    double vint = floor(value+0.01);
	    largestDelta = CoinMax(largestDelta,fabs(value-vint));
	    assert (largestDelta<1.0e-9);
	    columnElements[j] = vint;
	    assert (fabs(vint)>0.9);
	  } else if (rowScale[iRow]==-1.0) {
	    columnElements[j] = -value;
	  }
	}
      }
      if (largestDelta)
	printf("largest error while scaling rows - %g\n", largestDelta);
    }
    delete [] temp;
  }
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    columnActivity[iColumn] = 0.0;
    if (columnLower[iColumn]<-1.0e20 &&
	columnUpper[iColumn]>1.0e20)
      continue;
    double move = 0.0;
    bool moveLower = fabs(columnLower[iColumn]) <= fabs(columnUpper[iColumn]);
    if (moveLower) 
      move = -columnLower[iColumn];
    else 
      move = -columnUpper[iColumn];
    if (move) {
      nMoved++;
      offset += move*objective[iColumn];
      if (columnLower[iColumn]>-1.0e20)
	columnLower[iColumn] += move;
      if (columnUpper[iColumn]<1.0e20)
	columnUpper[iColumn] += move;
      columnActivity[iColumn] = move;
    }
  }
  if (nMoved) {
    memset(rowActivity,0,numberRows*sizeof(double));
    matrixByColumn->times(columnActivity, rowActivity);
    for (int iRow=0;iRow<numberRows;iRow++) {
      double move = rowActivity[iRow];
      if (move) {
	if (rowLower[iRow]>-1.0e20)
	  rowLower[iRow] += move;
	if (rowUpper[iRow]<1.0e20)
	  rowUpper[iRow] += move;
      }
    }
  }
  if (nMoved || nSwappedColumn || nScaled || nSwappedRow) {
    model.setObjectiveOffset(offset);
    model.dual();
    model.writeMps("cleaned.mps");
    printf("%d columns moved, %d swapped - %d rows scaled, %d swapped\n",
	   nMoved, nSwappedColumn, nScaled, nSwappedRow);
  } else {
    printf("No changes were made\n");
  }
  return 0;
}
