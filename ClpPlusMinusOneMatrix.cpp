// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
// at end to get min/max!
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix () 
  : ClpMatrixBase()
{
  setType(12);
  elements_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  numberRows_=0;
  numberColumns_=0;
  columnOrdered_=true;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix (const ClpPlusMinusOneMatrix & rhs) 
: ClpMatrixBase(rhs)
{  
  elements_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  columnOrdered_=rhs.columnOrdered_;
  if (numberColumns_) {
    int numberElements = rhs.startPositive_[numberColumns_];
    indices_ = new int [ numberElements];
    memcpy(indices_,rhs.indices_,numberElements*sizeof(int));
    startPositive_ = new int [ numberColumns_+1];
    memcpy(startPositive_,rhs.startPositive_,(numberColumns_+1)*sizeof(int));
    startNegative_ = new int [ numberColumns_];
    memcpy(startNegative_,rhs.startNegative_,numberColumns_*sizeof(int));
  }
}

ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix (const CoinPackedMatrix & rhs) 
  : ClpMatrixBase()
{  
  setType(12);
  elements_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  int iColumn;
  assert (rhs.isColOrdered());
  // get matrix data pointers
  const int * row = rhs.getIndices();
  const CoinBigIndex * columnStart = rhs.getVectorStarts();
  const int * columnLength = rhs.getVectorLengths(); 
  const double * elementByColumn = rhs.getElements();
  numberColumns_ = rhs.getNumCols();
  bool goodPlusMinusOne=true;
  numberRows_=-1;
  indices_ = new int[rhs.getNumElements()];
  startPositive_ = new int [numberColumns_+1];
  startNegative_ = new int [numberColumns_];
  int * temp = new int [rhs.getNumRows()];
  CoinBigIndex j=0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    CoinBigIndex k;
    int iNeg=0;
    startPositive_[iColumn]=j;
    for (k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];
	 k++) {
      int iRow;
      if (fabs(elementByColumn[k]-1.0)<1.0e-10) {
	iRow = row[k];
	numberRows_ = max(numberRows_,iRow);
	indices_[j++]=iRow;
      } else if (fabs(elementByColumn[k]+1.0)<1.0e-10) {
	iRow = row[k];
	numberRows_ = max(numberRows_,iRow);
	temp[iNeg++]=iRow;
      } else {
	goodPlusMinusOne = false; // not a network
      }
    }
    if (goodPlusMinusOne) {
      // move negative
      startNegative_[iColumn]=j;
      for (k=0;k<iNeg;k++) {
	indices_[j++] = temp[k];
      }
    } else {
      break;
    }
  }
  startPositive_[numberColumns_]=j;
  delete [] temp;
  if (!goodPlusMinusOne) {
    delete [] indices_;
    // put in message
    indices_=NULL;
    numberRows_=0;
    numberColumns_=0;
    delete [] startPositive_;
    delete [] startNegative_;
    startPositive_ = NULL;
    startNegative_ = NULL;
  } else {
    numberRows_ ++; //  correct
    columnOrdered_ = true;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::~ClpPlusMinusOneMatrix ()
{
  delete [] elements_;
  delete [] startPositive_;
  delete [] startNegative_;
  delete [] lengths_;
  delete [] indices_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix &
ClpPlusMinusOneMatrix::operator=(const ClpPlusMinusOneMatrix& rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete [] elements_;
    delete [] startPositive_;
    delete [] startNegative_;
    delete [] lengths_;
    delete [] indices_;
    elements_ = NULL;
    startPositive_ = NULL;
    lengths_=NULL;
    indices_=NULL;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    columnOrdered_=rhs.columnOrdered_;
    if (numberColumns_) {
      int numberElements = rhs.startPositive_[numberColumns_];
      indices_ = new int [ numberElements];
      memcpy(indices_,rhs.indices_,numberElements*sizeof(int));
      startPositive_ = new int [ numberColumns_+1];
      memcpy(startPositive_,rhs.startPositive_,(numberColumns_+1)*sizeof(int));
      startNegative_ = new int [ numberColumns_];
      memcpy(startNegative_,rhs.startNegative_,numberColumns_*sizeof(int));
    }
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase * ClpPlusMinusOneMatrix::clone() const
{
  return new ClpPlusMinusOneMatrix(*this);
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase * 
ClpPlusMinusOneMatrix::subsetClone (int numberRows, const int * whichRows,
			      int numberColumns, 
			      const int * whichColumns) const 
{
  return new ClpPlusMinusOneMatrix(*this, numberRows, whichRows,
				   numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix (
		       const ClpPlusMinusOneMatrix & rhs,
		       int numberRows, const int * whichRow,
		       int numberColumns, const int * whichColumn)
: ClpMatrixBase(rhs)
{  
  elements_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  numberRows_=0;
  numberColumns_=0;
  columnOrdered_=rhs.columnOrdered_;
  if (numberRows<=0||numberColumns<=0) {
    startPositive_ = new int[1];
    startPositive_[0] = 0;
  } else {
    numberColumns_ = numberColumns;
    numberRows_ = numberRows;
    const int * index1 = rhs.indices_;
    int * startPositive1 = rhs.startPositive_;

    int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    int numberMinor1 = (!columnOrdered_) ? rhs.numberColumns_ : rhs.numberRows_;
    int numberMajor1 = (columnOrdered_) ? rhs.numberColumns_ : rhs.numberRows_;
    // Also swap incoming if not column ordered
    if (!columnOrdered_) {
      int temp1 = numberRows;
      numberRows = numberColumns;
      numberColumns = temp1;
      const int * temp2;
      temp2 = whichRow;
      whichRow = whichColumn;
      whichColumn = temp2;
    }
    // Throw exception if rhs empty
    if (numberMajor1 <= 0 || numberMinor1 <= 0)
      throw CoinError("empty rhs", "subset constructor", "ClpPlusMinusOneMatrix");
    // Array to say if an old row is in new copy
    int * newRow = new int [numberMinor1];
    int iRow;
    for (iRow=0;iRow<numberMinor1;iRow++) 
      newRow[iRow] = -1;
    // and array for duplicating rows
    int * duplicateRow = new int [numberMinor];
    int numberBad=0;
    for (iRow=0;iRow<numberMinor;iRow++) {
      duplicateRow[iRow] = -1;
      int kRow = whichRow[iRow];
      if (kRow>=0  && kRow < numberMinor1) {
	if (newRow[kRow]<0) {
	  // first time
	  newRow[kRow]=iRow;
	} else {
	  // duplicate
	  int lastRow = newRow[kRow];
	  newRow[kRow]=iRow;
	  duplicateRow[iRow] = lastRow;
	}
      } else {
	// bad row
	numberBad++;
      }
    }

    if (numberBad)
      throw CoinError("bad minor entries", 
		      "subset constructor", "ClpPlusMinusOneMatrix");
    // now get size and check columns
    int size = 0;
    int iColumn;
    numberBad=0;
    for (iColumn=0;iColumn<numberMajor;iColumn++) {
      int kColumn = whichColumn[iColumn];
      if (kColumn>=0  && kColumn <numberMajor1) {
	int i;
	for (i=startPositive1[kColumn];i<startPositive1[kColumn+1];i++) {
	  int kRow = index1[i];
	  kRow = newRow[kRow];
	  while (kRow>=0) {
	    size++;
	    kRow = duplicateRow[kRow];
	  }
	}
      } else {
	// bad column
	numberBad++;
	printf("%d %d %d %d\n",iColumn,numberMajor,numberMajor1,kColumn);
      }
    }
    if (numberBad)
      throw CoinError("bad major entries", 
		      "subset constructor", "ClpPlusMinusOneMatrix");
    // now create arrays
    startPositive_ = new int [numberMajor+1];
    startNegative_ = new int [numberMajor];
    indices_ = new int[size];
    // and fill them
    size = 0;
    startPositive_[0]=0;
    int * startNegative1 = rhs.startNegative_;
    for (iColumn=0;iColumn<numberMajor;iColumn++) {
      int kColumn = whichColumn[iColumn];
      int i;
      for (i=startPositive1[kColumn];i<startNegative1[kColumn];i++) {
	int kRow = index1[i];
	kRow = newRow[kRow];
	while (kRow>=0) {
	  indices_[size++] = kRow;
	  kRow = duplicateRow[kRow];
	}
      }
      startNegative_[iColumn] = size;
      for (;i<startPositive1[kColumn+1];i++) {
	int kRow = index1[i];
	kRow = newRow[kRow];
	while (kRow>=0) {
	  indices_[size++] = kRow;
	  kRow = duplicateRow[kRow];
	}
      }
      startPositive_[iColumn+1] = size;
    }
    delete [] newRow;
    delete [] duplicateRow;
  }
}


/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase * 
ClpPlusMinusOneMatrix::reverseOrderedCopy() const
{
  int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  // count number in each row/column
  int * tempP = new int [numberMinor];
  int * tempN = new int [numberMinor];
  memset(tempP,0,numberMinor*sizeof(int));
  memset(tempN,0,numberMinor*sizeof(int));
  CoinBigIndex j=0;
  int i;
  for (i=0;i<numberMajor;i++) {
    for (;j<startNegative_[i];j++) {
      int iRow = indices_[j];
      tempP[iRow]++;
    }
    for (;j<startPositive_[i+1];j++) {
      int iRow = indices_[j];
      tempN[iRow]++;
    }
  }
  int * newIndices = new int [startPositive_[numberMajor]];
  int * newP = new int [numberMinor+1];
  int * newN = new int[numberMinor];
  int iRow;
  j=0;
  // do starts
  for (iRow=0;iRow<numberMinor;iRow++) {
    newP[iRow]=j;
    j += tempP[iRow];
    tempP[iRow]=newP[iRow];
    newN[iRow] = j;
    j += tempN[iRow];
    tempN[iRow]=newN[iRow];
  }
  newP[numberMinor]=j;
  j=0;
  for (i=0;i<numberMajor;i++) {
    for (;j<startNegative_[i];j++) {
      int iRow = indices_[j];
      int put = tempP[iRow];
      newIndices[put++] = i;
      tempP[iRow] = put;
    }
    for (;j<startPositive_[i+1];j++) {
      int iRow = indices_[j];
      int put = tempN[iRow];
      newIndices[put++] = i;
      tempN[iRow] = put;
    }
  }
  delete [] tempP;
  delete [] tempN;
  ClpPlusMinusOneMatrix * newCopy = new ClpPlusMinusOneMatrix();
  newCopy->passInCopy(numberMinor, numberMajor,
		      !columnOrdered_,  newIndices, newP, newN);
  return newCopy;
}
//unscaled versions
void 
ClpPlusMinusOneMatrix::times(double scalar,
		   const double * x, double * y) const
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  int i;
  CoinBigIndex j;
  assert (columnOrdered_);
  for (i=0;i<numberMajor;i++) {
    double value = scalar*x[i];
    if (value) {
      for (j=startPositive_[i];j<startNegative_[i];j++) {
	int iRow = indices_[j];
	y[iRow] += value;
      }
      for (;j<startPositive_[i+1];j++) {
	int iRow = indices_[j];
	y[iRow] -= value;
      }
    }
  }
}
void 
ClpPlusMinusOneMatrix::transposeTimes(double scalar,
				const double * x, double * y) const
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  int i;
  CoinBigIndex j=0;
  assert (columnOrdered_);
  for (i=0;i<numberMajor;i++) {
    double value = 0.0;
    for (;j<startNegative_[i];j++) {
      int iRow = indices_[j];
      value += x[iRow];
    }
    for (;j<startPositive_[i+1];j++) {
      int iRow = indices_[j];
      value -= x[iRow];
    }
    y[i] += scalar*value;
  }
}
void 
ClpPlusMinusOneMatrix::times(double scalar,
		       const double * x, double * y,
		       const double * rowScale, 
		       const double * columnScale) const
{
  // we know it is not scaled 
  times(scalar, x, y);
}
void 
ClpPlusMinusOneMatrix::transposeTimes( double scalar,
				 const double * x, double * y,
				 const double * rowScale, 
				 const double * columnScale) const
{
  // we know it is not scaled 
  transposeTimes(scalar, x, y);
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpPlusMinusOneMatrix::transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  // we know it is not scaled 
  columnArray->clear();
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  int numberRows = model->numberRows();
  bool packed = rowArray->packedMode();
  ClpPlusMinusOneMatrix* rowCopy =
    dynamic_cast< ClpPlusMinusOneMatrix*>(model->rowCopy());
  double factor = 0.3;
  // We may not want to do by row if there may be cache problems
  int numberColumns = model->numberColumns();
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberColumns*sizeof(double)>1000000) {
    if (numberRows*10<numberColumns)
      factor=0.1;
    else if (numberRows*4<numberColumns)
      factor=0.15;
    else if (numberRows*2<numberColumns)
      factor=0.2;
  }
  if (numberInRowArray>factor*numberRows||!rowCopy) {
    assert (!y->getNumElements());
    // do by column
    // Need to expand if packed mode
    int iColumn;
    CoinBigIndex j=0;
    assert (columnOrdered_);
    if (packed) {
      // need to expand pi into y
      assert(y->capacity()>=numberRows);
      double * piOld = pi;
      pi = y->denseVector();
      const int * whichRow = rowArray->getIndices();
      int i;
      // modify pi so can collapse to one loop
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=scalar*piOld[i];
      }
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	double value = 0.0;
	for (;j<startNegative_[iColumn];j++) {
	  int iRow = indices_[j];
	  value += pi[iRow];
	}
	for (;j<startPositive_[iColumn+1];j++) {
	  int iRow = indices_[j];
	  value -= pi[iRow];
	}
	if (fabs(value)>zeroTolerance) {
	  array[numberNonZero]=value;
	  index[numberNonZero++]=iColumn;
	}
      }
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=0.0;
      }
    } else {
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	double value = 0.0;
	for (;j<startNegative_[iColumn];j++) {
	  int iRow = indices_[j];
	  value += pi[iRow];
	}
	for (;j<startPositive_[iColumn+1];j++) {
	  int iRow = indices_[j];
	  value -= pi[iRow];
	}
	value *= scalar;
	if (fabs(value)>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	  array[iColumn]=value;
	}
      }
    }
    columnArray->setNumElements(numberNonZero);
  } else {
    // do by row
    rowCopy->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  }
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpPlusMinusOneMatrix::transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  const int * column = indices_;
  const CoinBigIndex * startPositive = startPositive_;
  const CoinBigIndex * startNegative = startNegative_;
  const int * whichRow = rowArray->getIndices();
  bool packed = rowArray->packedMode();
  if (numberInRowArray>2||y->getNumElements()) {
    // do by rows
    int iRow;
    double * markVector = y->denseVector(); // probably empty .. but
    int * mark = y->getIndices();
    int numberOriginal=y->getNumElements();
    int i;
    if (packed) {
      assert(!numberOriginal);
      numberNonZero=0;
      // and set up mark as char array
      char * marked = (char *) (index+columnArray->capacity());
      double * array2 = y->denseVector();
#ifdef CLP_DEBUG
      int numberColumns = model->numberColumns();
      for (i=0;i<numberColumns;i++) {
	assert(!marked[i]);
	assert(!array2[i]);
      }
#endif
      for (i=0;i<numberInRowArray;i++) {
	iRow = whichRow[i]; 
	double value = pi[i]*scalar;
	CoinBigIndex j;
	for (j=startPositive[iRow];j<startNegative[iRow];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array2[iColumn] += value;
	}
	for (j=startNegative[iRow];j<startPositive[iRow+1];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array2[iColumn] -= value;
	}
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	if (marked[iColumn]) {
	  double value = array2[iColumn];
	  array2[iColumn]=0.0;
	  marked[iColumn]=0;
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
    } else {      
      for (i=0;i<numberOriginal;i++) {
	int iColumn = mark[i];
	index[i]=iColumn;
	array[iColumn]=markVector[iColumn];
	markVector[iColumn]=0.0;
      }
      numberNonZero=numberOriginal;
      // and set up mark as char array
      char * marked = (char *) markVector;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	marked[iColumn]=0;
      }
      for (i=0;i<numberInRowArray;i++) {
	iRow = whichRow[i]; 
	double value = pi[iRow]*scalar;
	CoinBigIndex j;
	for (j=startPositive[iRow];j<startNegative[iRow];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array[iColumn] += value;
	}
	for (j=startNegative[iRow];j<startPositive[iRow+1];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array[iColumn] -= value;
	}
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	marked[iColumn]=0;
	if (fabs(array[iColumn])>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	} else {
	  array[iColumn]=0.0;
	}
      }
    }
  } else if (numberInRowArray==2) {
    /* do by rows when two rows (do longer first when not packed
       and shorter first if packed */
    int iRow0 = whichRow[0];
    int iRow1 = whichRow[1];
    int j;
    if (packed) {
      double pi0 = pi[0];
      double pi1 = pi[1];
      if (startPositive[iRow0+1]-startPositive[iRow0]>
	  startPositive[iRow1+1]-startPositive[iRow1]) {
	int temp = iRow0;
	iRow0 = iRow1;
	iRow1 = temp;
	pi0=pi1;
	pi1=pi[0];
      }
      // and set up mark as char array
      char * marked = (char *) (index+columnArray->capacity());
      int * lookup = y->getIndices();
      double value = pi0*scalar;
      for (j=startPositive[iRow0];j<startNegative[iRow0];j++) {
	int iColumn = column[j];
	array[numberNonZero] = value;
	marked[iColumn]=1;
	lookup[iColumn]=numberNonZero;
	index[numberNonZero++]=iColumn;
      }
      for (j=startNegative[iRow0];j<startPositive[iRow0+1];j++) {
	int iColumn = column[j];
	array[numberNonZero] = -value;
	marked[iColumn]=1;
	lookup[iColumn]=numberNonZero;
	index[numberNonZero++]=iColumn;
      }
      int numberOriginal = numberNonZero;
      value = pi1*scalar;
      for (j=startPositive[iRow1];j<startNegative[iRow1];j++) {
	int iColumn = column[j];
	if (marked[iColumn]) {
	  int iLookup = lookup[iColumn];
	  array[iLookup] += value;
	} else {
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero] = value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
      for (j=startNegative[iRow1];j<startPositive[iRow1+1];j++) {
	int iColumn = column[j];
	if (marked[iColumn]) {
	  int iLookup = lookup[iColumn];
	  array[iLookup] -= value;
	} else {
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero] = -value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
      // get rid of tiny values and zero out marked
      int nDelete=0;
      for (j=0;j<numberOriginal;j++) {
	int iColumn = index[j];
	marked[iColumn]=0;
	if (fabs(array[j])<=zeroTolerance) 
	  nDelete++;
      }
      if (nDelete) {
	numberOriginal=numberNonZero;
	numberNonZero=0;
	for (j=0;j<numberOriginal;j++) {
	  int iColumn = index[j];
	  double value = array[j];
	  array[j]=0.0;
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
    } else {
      if (startPositive[iRow0+1]-startPositive[iRow0]<
	  startPositive[iRow1+1]-startPositive[iRow1]) {
	int temp = iRow0;
	iRow0 = iRow1;
	iRow1 = temp;
      }
      int numberOriginal;
      int i;
      numberNonZero=0;
      double value;
      value = pi[iRow0]*scalar;
      CoinBigIndex j;
      for (j=startPositive[iRow0];j<startNegative[iRow0];j++) {
	int iColumn = column[j];
	index[numberNonZero++]=iColumn;
	array[iColumn] = value;
      }
      for (j=startNegative[iRow0];j<startPositive[iRow0+1];j++) {
	int iColumn = column[j];
	index[numberNonZero++]=iColumn;
	array[iColumn] = -value;
      }
      value = pi[iRow1]*scalar;
      for (j=startPositive[iRow1];j<startNegative[iRow1];j++) {
	int iColumn = column[j];
	double value2= array[iColumn];
	if (value2) {
	  value2 += value;
	} else {
	  value2 = value;
	  index[numberNonZero++]=iColumn;
	}
	array[iColumn] = value2;
      }
      for (j=startNegative[iRow1];j<startPositive[iRow1+1];j++) {
	int iColumn = column[j];
	double value2= array[iColumn];
	if (value2) {
	  value2 -= value;
	} else {
	  value2 = -value;
	  index[numberNonZero++]=iColumn;
	}
	array[iColumn] = value2;
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	if (fabs(array[iColumn])>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	} else {
	  array[iColumn]=0.0;
	}
      }
    }
  } else if (numberInRowArray==1) {
    // Just one row
    int iRow=rowArray->getIndices()[0];
    numberNonZero=0;
    double value;
    iRow = whichRow[0]; 
    CoinBigIndex j;
    if (packed) {
      value = pi[0]*scalar;
      if (fabs(value)>zeroTolerance) {
	for (j=startPositive[iRow];j<startNegative[iRow];j++) {
	  int iColumn = column[j];
	  array[numberNonZero] = value;
	  index[numberNonZero++]=iColumn;
	}
	for (j=startNegative[iRow];j<startPositive[iRow+1];j++) {
	  int iColumn = column[j];
	  array[numberNonZero] = -value;
	  index[numberNonZero++]=iColumn;
	}
      }
    } else {
      value = pi[iRow]*scalar;
      if (fabs(value)>zeroTolerance) {
	for (j=startPositive[iRow];j<startNegative[iRow];j++) {
	  int iColumn = column[j];
	  array[iColumn] = value;
	  index[numberNonZero++]=iColumn;
	}
	for (j=startNegative[iRow];j<startPositive[iRow+1];j++) {
	  int iColumn = column[j];
	  array[iColumn] = -value;
	  index[numberNonZero++]=iColumn;
	}
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y.
   Squashes small elements and knows about ClpSimplex */
void 
ClpPlusMinusOneMatrix::subsetTransposeTimes(const ClpSimplex * model,
			      const CoinIndexedVector * rowArray,
			      const CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  double * array = columnArray->denseVector();
  int jColumn;
  int numberToDo = y->getNumElements();
  const int * which = y->getIndices();
  assert (!rowArray->packedMode());
  columnArray->setPacked();
  for (jColumn=0;jColumn<numberToDo;jColumn++) {
    int iColumn = which[jColumn];
    double value = 0.0;
    CoinBigIndex j=startPositive_[iColumn];
    for (;j<startNegative_[iColumn];j++) {
      int iRow = indices_[j];
      value += pi[iRow];
    }
    for (;j<startPositive_[iColumn+1];j++) {
      int iRow = indices_[j];
      value -= pi[iRow];
    }
    array[jColumn]=value;
  }
}
/* Returns number of elements in basis
   column is basic if entry >=0 */
CoinBigIndex 
ClpPlusMinusOneMatrix::numberInBasis(const int * columnIsBasic) const 
{
  int i;
  CoinBigIndex numberElements=0;
  assert (columnOrdered_);
  for (i=0;i<numberColumns_;i++) {
    if (columnIsBasic[i]>=0) 
      numberElements += startPositive_[i+1]-startPositive_[i];
  }
  return numberElements;
}
// Fills in basis (Returns number of elements and updates numberBasic)
CoinBigIndex 
ClpPlusMinusOneMatrix::fillBasis(const ClpSimplex * model,
				const int * columnIsBasic, int & numberBasic,
				int * indexRowU, int * indexColumnU,
				double * elementU) const 
{
#ifdef CLPDEBUG
  const double * rowScale = model->rowScale();
  assert (!rowScale);
#endif
  int i;
  CoinBigIndex numberElements=0;
  assert (columnOrdered_);
  for (i=0;i<numberColumns_;i++) {
    if (columnIsBasic[i]>=0) {
      CoinBigIndex j=startPositive_[i];
      for (;j<startNegative_[i];j++) {
	int iRow = indices_[j];
	indexRowU[numberElements]=iRow;
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements++]=1.0;
      }
      for (;j<startPositive_[i+1];j++) {
	int iRow = indices_[j];
	indexRowU[numberElements]=iRow;
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements++]=-1.0;
      }
      numberBasic++;
    }
  }
  return numberElements;
}
/* If element NULL returns number of elements in column part of basis,
   If not NULL fills in as well */
CoinBigIndex 
ClpPlusMinusOneMatrix::fillBasis(const ClpSimplex * model,
				 const int * whichColumn, 
				 int numberBasic,
				 int numberColumnBasic,
				 int * indexRowU, int * indexColumnU,
				 double * elementU) const 
{
  int i;
  CoinBigIndex numberElements=0;
  if (elementU!=NULL) {
    assert (columnOrdered_);
    for (i=0;i<numberColumnBasic;i++) {
      int iColumn = whichColumn[i];
      CoinBigIndex j=startPositive_[iColumn];
      for (;j<startNegative_[iColumn];j++) {
	int iRow = indices_[j];
	indexRowU[numberElements]=iRow;
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements++]=1.0;
      }
      for (;j<startPositive_[iColumn+1];j++) {
	int iRow = indices_[j];
	indexRowU[numberElements]=iRow;
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements++]=-1.0;
      }
      numberBasic++;
    }
  } else {
    for (i=0;i<numberColumnBasic;i++) {
      int iColumn = whichColumn[i];
      numberElements += startPositive_[iColumn+1]-startPositive_[iColumn];
    }
  }
  return numberElements;
}
/* Unpacks a column into an CoinIndexedvector
 */
void 
ClpPlusMinusOneMatrix::unpack(const ClpSimplex * model,
			      CoinIndexedVector * rowArray,
			      int iColumn) const 
{
  CoinBigIndex j=startPositive_[iColumn];
  for (;j<startNegative_[iColumn];j++) {
    int iRow = indices_[j];
    rowArray->add(iRow,1.0);
  }
  for (;j<startPositive_[iColumn+1];j++) {
    int iRow = indices_[j];
    rowArray->add(iRow,-1.0);
  }
}
/* Unpacks a column into an CoinIndexedvector
** in packed foramt
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void 
ClpPlusMinusOneMatrix::unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int iColumn) const
{
  int * index = rowArray->getIndices();
  double * array = rowArray->denseVector();
  int number = 0;
  CoinBigIndex j=startPositive_[iColumn];
  for (;j<startNegative_[iColumn];j++) {
    int iRow = indices_[j];
    array[number]=1.0;
    index[number++]=iRow;
  }
  for (;j<startPositive_[iColumn+1];j++) {
    int iRow = indices_[j];
    array[number]=-1.0;
    index[number++]=iRow;
  }
  rowArray->setNumElements(number);
  rowArray->setPackedMode(true);
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void 
ClpPlusMinusOneMatrix::add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn, double multiplier) const 
{
  CoinBigIndex j=startPositive_[iColumn];
  for (;j<startNegative_[iColumn];j++) {
    int iRow = indices_[j];
    rowArray->quickAdd(iRow,multiplier);
  }
  for (;j<startPositive_[iColumn+1];j++) {
    int iRow = indices_[j];
    rowArray->quickAdd(iRow,-multiplier);
  }
}

// Return a complete CoinPackedMatrix
CoinPackedMatrix * 
ClpPlusMinusOneMatrix::getPackedMatrix() const 
{
  int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  return new CoinPackedMatrix(columnOrdered_,numberMinor,numberMajor,
			      getNumElements(),
			      getElements(),indices_,
			      startPositive_,getVectorLengths());

}
/* A vector containing the elements in the packed matrix. Note that there
   might be gaps in this list, entries that do not belong to any
   major-dimension vector. To get the actual elements one should look at
   this vector together with vectorStarts and vectorLengths. */
const double * 
ClpPlusMinusOneMatrix::getElements() const 
{
  if (!elements_) {
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    int numberElements = startPositive_[numberMajor];
    elements_ = new double [numberElements];
    CoinBigIndex j=0;
    int i;
    for (i=0;i<numberMajor;i++) {
      for (;j<startNegative_[i];j++) {
	elements_[j]=1.0;
      }
      for (;j<startPositive_[i+1];j++) {
	elements_[j]=-1.0;
      }
    }
  }
  return elements_;
}

const CoinBigIndex * 
ClpPlusMinusOneMatrix::getVectorStarts() const 
{
  return startPositive_;
}
/* The lengths of the major-dimension vectors. */
const int * 
ClpPlusMinusOneMatrix::getVectorLengths() const
{
  if (!lengths_) {
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    lengths_ = new int [numberMajor];
    int i;
    for (i=0;i<numberMajor;i++) {
      lengths_[i]=startPositive_[i+1]-startPositive_[i];
    }
  }
  return lengths_;
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void 
ClpPlusMinusOneMatrix::deleteCols(const int numDel, const int * indDel) 
{
  int iColumn;
  int newSize=startPositive_[numberColumns_];;
  int numberBad=0;
  // Use array to make sure we can have duplicates
  int * which = new int[numberColumns_];
  memset(which,0,numberColumns_*sizeof(int));
  int nDuplicate=0;
  for (iColumn=0;iColumn<numDel;iColumn++) {
    int jColumn = indDel[iColumn];
    if (jColumn<0||jColumn>=numberColumns_) {
      numberBad++;
    } else {
      newSize -= startPositive_[jColumn+1]-startPositive_[jColumn];
      if (which[jColumn])
	nDuplicate++;
      else
	which[jColumn]=1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteCols", "ClpPlusMinusOneMatrix");
  int newNumber = numberColumns_-numDel+nDuplicate;
  // Get rid of temporary arrays
  delete [] lengths_;
  lengths_=NULL;
  delete [] elements_;
  elements_= NULL;
  int * newPositive = new int [newNumber+1];
  int * newNegative = new int [newNumber];
  int * newIndices = new int [newSize];
  newNumber=0;
  newSize=0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (!which[iColumn]) {
      int start,end;
      int i;
      start = startPositive_[iColumn];
      end=startNegative_[iColumn];
      newPositive[newNumber]=newSize;
      for (i=start;i<end;i++) 
	newIndices[newSize++]=indices_[i];
      start = startNegative_[iColumn];
      end=startPositive_[iColumn+1];
      newNegative[newNumber++]=newSize;
      for (i=start;i<end;i++) 
	newIndices[newSize++]=indices_[i];
    }
  }
  newPositive[newNumber]=newSize;
  delete [] which;
  delete [] startPositive_;
  startPositive_= newPositive;
  delete [] startNegative_;
  startNegative_= newNegative;
  delete [] indices_;
  indices_= newIndices;
  numberColumns_ = newNumber;
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void 
ClpPlusMinusOneMatrix::deleteRows(const int numDel, const int * indDel) 
{
  int iRow;
  int numberBad=0;
  // Use array to make sure we can have duplicates
  int * which = new int[numberRows_];
  memset(which,0,numberRows_*sizeof(int));
  int nDuplicate=0;
  for (iRow=0;iRow<numDel;iRow++) {
    int jRow = indDel[iRow];
    if (jRow<0||jRow>=numberRows_) {
      numberBad++;
    } else {
      if (which[jRow])
	nDuplicate++;
      else
	which[jRow]=1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteCols", "ClpPlusMinusOneMatrix");
  int iElement;
  int numberElements=startPositive_[numberColumns_];
  int newSize=0;
  for (iElement=0;iElement<numberElements;iElement++) {
    iRow = indices_[iElement];
    if (!which[iRow])
      newSize++;
  }
  int newNumber = numberRows_-numDel+nDuplicate;
  // Get rid of temporary arrays
  delete [] lengths_;
  lengths_=NULL;
  delete [] elements_;
  elements_= NULL;
  int * newIndices = new int [newSize];
  newSize=0;
  int iColumn;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    int start,end;
    int i;
    start = startPositive_[iColumn];
    end=startNegative_[iColumn];
    startPositive_[newNumber]=newSize;
    for (i=start;i<end;i++) {
	iRow = indices_[i];
	if (!which[iRow])
	  newIndices[newSize++]=iRow;
    }
    start = startNegative_[iColumn];
    end=startPositive_[iColumn+1];
    startNegative_[newNumber]=newSize;
    for (i=start;i<end;i++) {
	iRow = indices_[i];
	if (!which[iRow])
	  newIndices[newSize++]=iRow;
    }
  }
  startPositive_[numberColumns_]=newSize;
  delete [] which;
  delete [] indices_;
  indices_= newIndices;
  numberRows_ = newNumber;
}
bool 
ClpPlusMinusOneMatrix::isColOrdered() const 
{ 
  return columnOrdered_;
}
/* Number of entries in the packed matrix. */
CoinBigIndex 
ClpPlusMinusOneMatrix::getNumElements() const 
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  if (startPositive_) 
    return startPositive_[numberMajor];
  else
    return 0;
}
// pass in copy (object takes ownership)
void 
ClpPlusMinusOneMatrix::passInCopy(int numberRows, int numberColumns,
		  bool columnOrdered, int * indices,
		  int * startPositive, int * startNegative)
{
  columnOrdered_=columnOrdered;
  startPositive_ = startPositive;
  startNegative_ = startNegative;
  indices_ = indices;
  numberRows_=numberRows;
  numberColumns_=numberColumns;
}
/* Given positive integer weights for each row fills in sum of weights
   for each column (and slack).
   Returns weights vector
*/
CoinBigIndex * 
ClpPlusMinusOneMatrix::dubiousWeights(const ClpSimplex * model,int * inputWeights) const
{
  int numberRows = model->numberRows();
  int numberColumns =model->numberColumns();
  int number = numberRows+numberColumns;
  CoinBigIndex * weights = new CoinBigIndex[number];
  int i;
  for (i=0;i<numberColumns;i++) {
    CoinBigIndex j;
    CoinBigIndex count=0;
    for (j=startPositive_[i];j<startPositive_[i+1];j++) {
      int iRow=indices_[j];
      count += inputWeights[iRow];
    }
    weights[i]=count;
  }
  for (i=0;i<numberRows;i++) {
    weights[i+numberColumns]=inputWeights[i];
  }
  return weights;
}
// Append Columns
void 
ClpPlusMinusOneMatrix::appendCols(int number, const CoinPackedVectorBase * const * columns)
{
  int iColumn;
  int size=0;
  int numberBad=0;
  for (iColumn=0;iColumn<number;iColumn++) {
    int n=columns[iColumn]->getNumElements();
    const double * element = columns[iColumn]->getElements();
    size += n;
    int i;
    for (i=0;i<n;i++) {
      if (fabs(element[i])!=1.0)
	numberBad++;
    }
  }
  if (numberBad)
    throw CoinError("Not +- 1", "appendCols", "ClpPlusMinusOneMatrix");
  // Get rid of temporary arrays
  delete [] lengths_;
  lengths_=NULL;
  delete [] elements_;
  elements_= NULL;
  int numberNow = startPositive_[numberColumns_];
  int * temp;
  temp = new int [numberColumns_+1+number];
  memcpy(temp,startPositive_,(numberColumns_+1)*sizeof(int));
  delete [] startPositive_;
  startPositive_= temp;
  temp = new int [numberColumns_+number];
  memcpy(temp,startNegative_,numberColumns_*sizeof(int));
  delete [] startNegative_;
  startNegative_= temp;
  temp = new int [numberNow+size];
  memcpy(temp,indices_,numberNow*sizeof(int));
  delete [] indices_;
  indices_= temp;
  // now add
  size=numberNow;
  for (iColumn=0;iColumn<number;iColumn++) {
    int n=columns[iColumn]->getNumElements();
    const int * row = columns[iColumn]->getIndices();
    const double * element = columns[iColumn]->getElements();
    int i;
    for (i=0;i<n;i++) {
      if (element[i]==1.0) 
	indices_[size++]=row[i];
    }
    startNegative_[iColumn+numberColumns_]=size;
    for (i=0;i<n;i++) {
      if (element[i]==-1.0) 
	indices_[size++]=row[i];
    }
    startPositive_[iColumn+numberColumns_+1]=size;
  }
  
  numberColumns_ += number;
}
// Append Rows
void 
ClpPlusMinusOneMatrix::appendRows(int number, const CoinPackedVectorBase * const * rows)
{
  // Allocate arrays to use for counting
  int * countPositive = new int [numberColumns_+1];
  memset(countPositive,0,numberColumns_*sizeof(int));
  int * countNegative = new int [numberColumns_];
  memset(countNegative,0,numberColumns_*sizeof(int));
  int iRow;
  int size=0;
  int numberBad=0;
  for (iRow=0;iRow<number;iRow++) {
    int n=rows[iRow]->getNumElements();
    const int * column = rows[iRow]->getIndices();
    const double * element = rows[iRow]->getElements();
    size += n;
    int i;
    for (i=0;i<n;i++) {
      int iColumn = column[i];
      if (element[i]==1.0)
	countPositive[iColumn]++;
      else if (element[i]==-1.0)
	countNegative[iColumn]++;
      else
	numberBad++;
    }
  }
  if (numberBad)
    throw CoinError("Not +- 1", "appendRows", "ClpPlusMinusOneMatrix");
  // Get rid of temporary arrays
  delete [] lengths_;
  lengths_=NULL;
  delete [] elements_;
  elements_= NULL;
  int numberNow = startPositive_[numberColumns_];
  int * newIndices = new int [numberNow+size];
  // Update starts and turn counts into positions
  // also move current indices
  int iColumn;
  int numberAdded=0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    int n,now,move;
    now = startPositive_[iColumn];
    move = startNegative_[iColumn]-now;
    n = countPositive[iColumn];
    startPositive_[iColumn] += numberAdded;
    memcpy(indices_+now,newIndices+startPositive_[iColumn],move*sizeof(int));
    countPositive[iColumn]= startNegative_[iColumn]+numberAdded;
    numberAdded += n;
    now = startNegative_[iColumn];
    move = startPositive_[iColumn+1]-now;
    n = countNegative[iColumn];
    startNegative_[iColumn] += numberAdded;
    memcpy(indices_+now,newIndices+startNegative_[iColumn],move*sizeof(int));
    countNegative[iColumn]= startPositive_[iColumn+1]+numberAdded;
    numberAdded += n;
  }
  delete [] indices_;
  indices_ = newIndices;
  startPositive_[numberColumns_] += numberAdded;
  // Now put in
  for (iRow=0;iRow<number;iRow++) {
    int newRow = numberRows_+iRow;
    int n=rows[iRow]->getNumElements();
    const int * column = rows[iRow]->getIndices();
    const double * element = rows[iRow]->getElements();
    int i;
    for (i=0;i<n;i++) {
      int iColumn = column[i];
      int put;
      if (element[i]==1.0) {
	put = countPositive[iColumn];
	countPositive[iColumn] = put+1;
      } else {
	put = countNegative[iColumn];
	countNegative[iColumn] = put+1;
      }
      indices_[put]=newRow;
    }
  }
  delete [] countPositive;
  delete [] countNegative;
  numberRows_ += number;
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void 
ClpPlusMinusOneMatrix::rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive)
{
  int iColumn;
  bool plusOne=false;
  bool minusOne=false;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (startNegative_[iColumn]>startPositive_[iColumn])
      plusOne=true;
    if (startPositive_[iColumn+1]>startNegative_[iColumn])
      minusOne=true;
  }
  if (minusOne) {
    smallestNegative=-1.0;
    largestNegative=-1.0;
  } else {
    smallestNegative=0.0;
    largestNegative=0.0;
  }
  if (plusOne) {
    smallestPositive=1.0;
    largestPositive=1.0;
  } else {
    smallestPositive=0.0;
    largestPositive=0.0;
  }
}
// Says whether it can do partial pricing
bool 
ClpPlusMinusOneMatrix::canDoPartialPricing() const
{
  return true;
}
// Partial pricing 
void 
ClpPlusMinusOneMatrix::partialPricing(ClpSimplex * model, int start, int end,
			      int & bestSequence, int & numberWanted)
{
  int j;
  double tolerance=model->currentDualTolerance();
  double * reducedCost = model->djRegion();
  const double * duals = model->dualRowSolution();
  const double * cost = model->costRegion();
  double bestDj;
  if (bestSequence>=0)
    bestDj = fabs(reducedCost[bestSequence]);
  else
    bestDj=tolerance;
  int sequenceOut = model->sequenceOut();
  int saveSequence = bestSequence;
  int iSequence;
  for (iSequence=start;iSequence<end;iSequence++) {
    if (iSequence!=sequenceOut) {
      double value;
      ClpSimplex::Status status = model->getStatus(iSequence);
      
      switch(status) {
	
      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
	break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
	value=cost[iSequence];
	j=startPositive_[iSequence];
	for (;j<startNegative_[iSequence];j++) {
	  int iRow = indices_[j];
	  value -= duals[iRow];
	}
	for (;j<startPositive_[iSequence+1];j++) {
	  int iRow = indices_[j];
	  value += duals[iRow];
	}
	value = fabs(value);
	if (value>FREE_ACCEPT*tolerance) {
	  numberWanted--;
	  // we are going to bias towards free (but only if reasonable)
	  value *= FREE_BIAS;
	  if (value>bestDj) {
	    // check flagged variable and correct dj
	    if (!model->flagged(iSequence)) {
	      bestDj=value;
	      bestSequence = iSequence;
	    } else {
	      // just to make sure we don't exit before got something
	      numberWanted++;
	    }
	  }
	}
	break;
      case ClpSimplex::atUpperBound:
	value=cost[iSequence];
	j=startPositive_[iSequence];
	for (;j<startNegative_[iSequence];j++) {
	  int iRow = indices_[j];
	  value -= duals[iRow];
	}
	for (;j<startPositive_[iSequence+1];j++) {
	  int iRow = indices_[j];
	  value += duals[iRow];
	}
	if (value>tolerance) {
	  numberWanted--;
	  if (value>bestDj) {
	    // check flagged variable and correct dj
	    if (!model->flagged(iSequence)) {
	      bestDj=value;
	      bestSequence = iSequence;
	    } else {
	      // just to make sure we don't exit before got something
	      numberWanted++;
	    }
	  }
	}
	break;
      case ClpSimplex::atLowerBound:
	value=cost[iSequence];
	j=startPositive_[iSequence];
	for (;j<startNegative_[iSequence];j++) {
	  int iRow = indices_[j];
	  value -= duals[iRow];
	}
	for (;j<startPositive_[iSequence+1];j++) {
	  int iRow = indices_[j];
	  value += duals[iRow];
	}
	value = -value;
	if (value>tolerance) {
	  numberWanted--;
	  if (value>bestDj) {
	    // check flagged variable and correct dj
	    if (!model->flagged(iSequence)) {
	      bestDj=value;
	      bestSequence = iSequence;
	    } else {
	      // just to make sure we don't exit before got something
	      numberWanted++;
	    }
	  }
	}
	break;
      }
    }
    if (!numberWanted)
      break;
  }
  if (bestSequence!=saveSequence) {
    // recompute dj
    double value=cost[bestSequence];
    j=startPositive_[bestSequence];
    for (;j<startNegative_[bestSequence];j++) {
      int iRow = indices_[j];
      value -= duals[iRow];
    }
    for (;j<startPositive_[bestSequence+1];j++) {
      int iRow = indices_[j];
      value += duals[iRow];
    }
    reducedCost[bestSequence] = value;
  }
}
