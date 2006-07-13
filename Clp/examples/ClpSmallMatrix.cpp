// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"

#include "ClpSimplex.hpp"
#include "ClpSmallMatrix.hpp"
#include "ClpFactorization.hpp"
#include "ClpMessage.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpSmallMatrix::ClpSmallMatrix () 
  : ClpMatrixBase()
{
  setType(15);
  numberRows_=0;
  numberColumns_=0;
  columnStart_ = new int [1];
  columnStart_[0]=0;
  row_ = NULL;
  element_ = NULL;
}

/* Constructor from data */
ClpSmallMatrix::ClpSmallMatrix(int numberColumns, int numberRows,
			       int * starts, Int * row, Double * element)
  : ClpMatrixBase()
{
  setType(15);
  numberRows_=numberRows;
  numberColumns_=numberColumns;
  columnStart_ = starts;
  row_ = row;
  element_ = element;
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpSmallMatrix::ClpSmallMatrix (const ClpSmallMatrix & rhs) 
: ClpMatrixBase(rhs)
{  
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  columnStart_ = CoinCopyOfArray(rhs.columnStart_,numberColumns_+1);
  int numberElements = columnStart_[numberColumns_];
  row_ = CoinCopyOfArray(rhs.row_,numberElements);
  element_ = CoinCopyOfArray(rhs.element_,numberElements);
}

ClpSmallMatrix::ClpSmallMatrix (const CoinPackedMatrix & rhs) 
  : ClpMatrixBase()
{ 
  // This constructor is not strictly necessary and could abort
  numberRows_ = rhs.getNumRows();
  if (sizeof(Int)==2)
    assert (numberRows_<65536);
  numberColumns_ = rhs.getNumCols();
  int numberElements = rhs.getNumElements();
  assert (rhs.isColOrdered());
  // get matrix data pointers
  const int * row = rhs.getIndices();
  const CoinBigIndex * columnStart = rhs.getVectorStarts();
  const int * columnLength = rhs.getVectorLengths(); 
  const double * elementByColumn = rhs.getElements();
  columnStart_ = new int [numberColumns_+1];
  row_ = new Int [numberElements];
  element_ = new Double [numberElements];
  columnStart_[0]=0;
  numberElements=0;
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    for (CoinBigIndex j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];
	 j++) {
      row_[numberElements] = row[j];
      element_[numberElements++] = elementByColumn[j];
    }
    columnStart_[iColumn+1] = numberElements;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpSmallMatrix::~ClpSmallMatrix ()
{
  delete [] columnStart_;
  delete [] row_;
  delete [] element_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpSmallMatrix &
ClpSmallMatrix::operator=(const ClpSmallMatrix& rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    delete [] columnStart_;
    delete [] row_;
    delete [] element_;
    columnStart_ = CoinCopyOfArray(rhs.columnStart_,numberColumns_+1);
    int numberElements = columnStart_[numberColumns_];
    row_ = CoinCopyOfArray(rhs.row_,numberElements);
    element_ = CoinCopyOfArray(rhs.element_,numberElements);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase * ClpSmallMatrix::clone() const
{
  return new ClpSmallMatrix(*this);
}

//unscaled versions
void 
ClpSmallMatrix::times(double scalar,
		   const double * x, double * y) const
{
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    double value = scalar*x[iColumn];
    if (value) {
      for (CoinBigIndex j=columnStart_[iColumn];
	   j<columnStart_[iColumn+1];j++) {
	int iRow=row_[j];
	y[iRow] += value*element_[j];
      }
    }
  }
}
void 
ClpSmallMatrix::transposeTimes(double scalar,
				const double * x, double * y) const
{
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    double value = 0.0;
    for (CoinBigIndex j=columnStart_[iColumn];
	 j<columnStart_[iColumn+1];j++) {
      int iRow=row_[j];
      value += x[iRow] * element_[j];
    }
    y[iColumn] += value*scalar;
  }
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpSmallMatrix::transposeTimes(const ClpSimplex * model, double scalar,
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
  bool packed = rowArray->packedMode();
  const int * whichRow = rowArray->getIndices();
  if (packed) {
    // need to expand pi into y
    double * piOld = pi;
    pi = y->denseVector();
    // modify pi so can collapse to one loop
    for (int i=0;i<numberInRowArray;i++) {
      int iRow = whichRow[i];
      pi[iRow]=scalar*piOld[i];
    }
  }
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    double value = 0.0;
    for (CoinBigIndex j=columnStart_[iColumn];
	 j<columnStart_[iColumn+1];j++) {
      int iRow=row_[j];
      value += pi[iRow]*element_[j];
    }
    if (fabs(value)>zeroTolerance) {
      array[numberNonZero]=value;
      index[numberNonZero++]=iColumn;
    }
  }
  if (packed) {
    // zero out
    for (int i=0;i<numberInRowArray;i++) {
      int iRow = whichRow[i];
      pi[iRow]=0.0;
    }
    columnArray->setPackedMode(true);
  }
  columnArray->setNumElements(numberNonZero);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y */
void 
ClpSmallMatrix::subsetTransposeTimes(const ClpSimplex * model,
			      const CoinIndexedVector * rowArray,
			      const CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  double * array = columnArray->denseVector();
  int numberToDo = y->getNumElements();
  const int * which = y->getIndices();
  assert (!rowArray->packedMode());
  columnArray->setPacked();
  for (int jColumn=0;jColumn<numberToDo;jColumn++) {
    int iColumn = which[jColumn];
    double value = 0.0;
    for (CoinBigIndex j=columnStart_[iColumn];
	 j<columnStart_[iColumn+1];j++) {
      int iRow=row_[j];
      value += pi[iRow] * element_[j];
    }
    array[jColumn] = value;
  }
}
/// returns number of elements in column part of basis,
CoinBigIndex 
ClpSmallMatrix::countBasis(ClpSimplex * model,
				 const int * whichColumn, 
				 int numberBasic,
			   int & numberColumnBasic)
{
  CoinBigIndex numberElements=0;
  // just count - can be over so ignore zero problem
  for (int i=0;i<numberColumnBasic;i++) {
    int iColumn = whichColumn[i];
    numberElements += columnStart_[iColumn+1]-columnStart_[iColumn];
  }
  return numberElements;
}
void
ClpSmallMatrix::fillBasis(ClpSimplex * model,
			 const int * whichColumn, 
			 int & numberColumnBasic,
			 int * indexRowU, int * start,
			 int * rowCount, int * columnCount,
			 double * elementU)
{
  CoinBigIndex numberElements=start[0];
  assert(!model->rowScale()); // scaling not allowed
  // fill
  for (int i=0;i<numberColumnBasic;i++) {
    int iColumn = whichColumn[i];
    for (CoinBigIndex j=columnStart_[iColumn];
	 j<columnStart_[iColumn+1];j++) {
      int iRow=row_[j];
      indexRowU[numberElements]=iRow;
      rowCount[iRow]++;
      elementU[numberElements++]=element_[j];
    }
    start[i+1]=numberElements;
    columnCount[i]=columnStart_[iColumn+1]-columnStart_[iColumn];
  }
}
/* Unpacks a column into an CoinIndexedvector
 */
void 
ClpSmallMatrix::unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn) const 
{
  for (CoinBigIndex i=columnStart_[iColumn];
       i<columnStart_[iColumn+1];i++) {
    rowArray->add(row_[i],element_[i]);
  }
}
/* Unpacks a column into an CoinIndexedvector
** in packed format
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void 
ClpSmallMatrix::unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int iColumn) const
{
  int * index = rowArray->getIndices();
  double * array = rowArray->denseVector();
  int number = 0;
  for (CoinBigIndex i=columnStart_[iColumn];
       i<columnStart_[iColumn+1];i++) {
    array[number]=element_[i];
    index[number++]=row_[i];
  }
  rowArray->setNumElements(number);
  rowArray->setPackedMode(true);
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void 
ClpSmallMatrix::add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn, double multiplier) const 
{
  for (CoinBigIndex i=columnStart_[iColumn];
       i<columnStart_[iColumn+1];i++) {
    rowArray->quickAdd(row_[i],multiplier*element_[i]);
  }
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void 
ClpSmallMatrix::rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive)
{
  smallestNegative=-COIN_DBL_MAX;
  largestNegative=0.0;
  smallestPositive=COIN_DBL_MAX;
  largestPositive=0.0;
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    for (CoinBigIndex j=columnStart_[iColumn];
	 j<columnStart_[iColumn+1];j++) {
      double value = element_[j];
      if (value>0.0) {
	smallestPositive = CoinMin(smallestPositive,value);
	largestPositive = CoinMax(largestPositive,value);
      } else if (value<0.0) {
	smallestNegative = CoinMax(smallestNegative,value);
	largestNegative = CoinMin(largestNegative,value);
      }
    }
  }
}
/* Adds multiple of a column into an array */
void 
ClpSmallMatrix::add(const ClpSimplex * model,double * array,
		    int column, double multiplier) const
{
  std::cerr<<"add not supported - ClpSmallMatrix"<<std::endl;
  abort();
}

// Return a complete CoinPackedMatrix
CoinPackedMatrix * 
ClpSmallMatrix::getPackedMatrix() const 
{
  std::cerr<<"getPackedMatrix not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}
/* A vector containing the elements in the packed matrix. Note that there
   might be gaps in this list, entries that do not belong to any
   major-dimension vector. To get the actual elements one should look at
   this vector together with vectorStarts and vectorLengths. */
const double * 
ClpSmallMatrix::getElements() const 
{
  std::cerr<<"getElements not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}

const CoinBigIndex * 
ClpSmallMatrix::getVectorStarts() const 
{
  std::cerr<<"getVectorStarts not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}
/* The lengths of the major-dimension vectors. */
const int * 
ClpSmallMatrix::getVectorLengths() const
{
  std::cerr<<"get VectorLengths not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void ClpSmallMatrix::deleteCols(const int numDel, const int * indDel) 
{
  std::cerr<<"deleteCols not supported - ClpSmallMatrix"<<std::endl;
  abort();
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void ClpSmallMatrix::deleteRows(const int numDel, const int * indDel) 
{
  std::cerr<<"deleteRows not supported - ClpSmallMatrix"<<std::endl;
  abort();
}
const int *
ClpSmallMatrix::getIndices() const
{
  std::cerr<<"getIndices not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}
/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase * 
ClpSmallMatrix::reverseOrderedCopy() const
{
  std::cerr<<"reverseOrderedCopy not supported - ClpSmallMatrix"<<std::endl;
  abort();
  return NULL;
}
void 
ClpSmallMatrix::times(double scalar,
		       const double * x, double * y,
		       const double * rowScale, 
		       const double * columnScale) const
{
  std::cerr<<"timesnot supported - ClpSmallMatrix"<<std::endl;
  abort();
}
void 
ClpSmallMatrix::transposeTimes( double scalar,
				 const double * x, double * y,
				 const double * rowScale, 
				 const double * columnScale) const
{
  std::cerr<<"transposeTimesnot supported - ClpSmallMatrix"<<std::endl;
  abort();
}
