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
    indices_ = new int [ 2*numberColumns_];
    memcpy(indices_,rhs.indices_,2*numberColumns_*sizeof(int));
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
    printf("Not all +-1 - can test if indices_ null\n");
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
      indices_ = new int [ 2*numberColumns_];
      memcpy(indices_,rhs.indices_,2*numberColumns_*sizeof(int));
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
  ClpPlusMinusOneMatrix* rowCopy =
    dynamic_cast< ClpPlusMinusOneMatrix*>(model->rowCopy());
  if (numberInRowArray>0.3*numberRows||!rowCopy) {
    // do by column
    int iColumn;
    double * markVector = y->denseVector(); // probably empty
    CoinBigIndex j=0;
    assert (columnOrdered_);
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      double value2 = 0.0;
      for (;j<startNegative_[iColumn];j++) {
	int iRow = indices_[j];
	value2 += pi[iRow];
      }
      for (;j<startPositive_[iColumn+1];j++) {
	int iRow = indices_[j];
	value2 -= pi[iRow];
      }
      double value = markVector[iColumn];
      markVector[iColumn]=0.0;
      value += scalar*value2;
      if (fabs(value)>zeroTolerance) {
	index[numberNonZero++]=iColumn;
	array[iColumn]=value;
      }
    }
    columnArray->setNumElements(numberNonZero);
    y->setNumElements(0);
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
  if (numberInRowArray>2||y->getNumElements()) {
    // do by rows
    int iRow;
    double * markVector = y->denseVector(); // probably empty .. but
    int * mark = y->getIndices();
    int numberOriginal=y->getNumElements();
    int i;
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
  } else if (numberInRowArray==2) {
    // do by rows when two rows (do longer first)
    int iRow0 = whichRow[0];
    int iRow1 = whichRow[1];
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
  } else if (numberInRowArray==1) {
    // Just one row
    int iRow=rowArray->getIndices()[0];
    numberNonZero=0;
    double value;
    iRow = whichRow[0]; 
    value = pi[iRow]*scalar;
    if (fabs(value)>zeroTolerance) {
      CoinBigIndex j;
      for (j=startPositive[iRow];j<startNegative[iRow];j++) {
	int iColumn = column[j];
	index[numberNonZero++]=iColumn;
	array[iColumn] = value;
      }
      for (j=startNegative[iRow];j<startPositive[iRow+1];j++) {
	int iColumn = column[j];
	index[numberNonZero++]=iColumn;
	array[iColumn] = -value;
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
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  int jColumn;
  int numberToDo = y->getNumElements();
  const int * which = y->getIndices();
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
    if (fabs(value)>zeroTolerance) {
      index[numberNonZero++]=iColumn;
      array[iColumn]=value;
    }
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
/* Unpacks a column into an CoinIndexedvector
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation */
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
  abort();
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void 
ClpPlusMinusOneMatrix::deleteRows(const int numDel, const int * indDel) 
{
  abort();
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
