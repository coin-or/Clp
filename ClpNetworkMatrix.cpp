// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
// at end to get min/max!
#include "ClpNetworkMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNetworkMatrix::ClpNetworkMatrix () 
  : ClpMatrixBase()
{
  setType(11);
  elements_ = NULL;
  starts_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  numberRows_=0;
  numberColumns_=0;
  trueNetwork_=false;
}

/* Constructor from two arrays */
ClpNetworkMatrix::ClpNetworkMatrix(int numberColumns, const int * head,
				   const int * tail)
  : ClpMatrixBase()
{
  setType(11);
  elements_ = NULL;
  starts_ = NULL;
  lengths_=NULL;
  indices_=new int[2*numberColumns];;
  numberRows_=-1;
  numberColumns_=numberColumns;
  trueNetwork_=true;
  int iColumn;
  CoinBigIndex j=0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++, j+=2) {
    int iRow = head[iColumn];
    numberRows_ = max(numberRows_,iRow);
    indices_[j]=iRow;
    iRow = tail[iColumn];
    numberRows_ = max(numberRows_,iRow);
    indices_[j+1]=iRow;
  }
  numberRows_++;
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNetworkMatrix::ClpNetworkMatrix (const ClpNetworkMatrix & rhs) 
: ClpMatrixBase(rhs)
{  
  elements_ = NULL;
  starts_ = NULL;
  lengths_=NULL;
  indices_=NULL;
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  trueNetwork_=rhs.trueNetwork_;
  if (numberColumns_) {
    indices_ = new int [ 2*numberColumns_];
    memcpy(indices_,rhs.indices_,2*numberColumns_*sizeof(int));
  }
}

ClpNetworkMatrix::ClpNetworkMatrix (const CoinPackedMatrix & rhs) 
  : ClpMatrixBase()
{  
  setType(11);
  elements_ = NULL;
  starts_ = NULL;
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
  int goodNetwork=1;
  numberRows_=-1;
  indices_ = new int[2*numberColumns_];
  CoinBigIndex j=0;
  for (iColumn=0;iColumn<numberColumns_;iColumn++, j+=2) {
    CoinBigIndex k=columnStart[iColumn];
    int iRow;
    switch (columnLength[iColumn]) {
    case 0:
      goodNetwork=-1; // not classic network
      indices_[j]=-1;
      indices_[j+1]=-1;
      break;
      
    case 1:
      goodNetwork=-1; // not classic network
      if (fabs(elementByColumn[k]-1.0)<1.0e-10) {
	indices_[j] = -1;
	iRow = row[k];
	numberRows_ = max(numberRows_,iRow);
	indices_[j+1]=iRow;
      } else if (fabs(elementByColumn[k]+1.0)<1.0e-10) {
	indices_[j+1] = -1;
	iRow = row[k];
	numberRows_ = max(numberRows_,iRow);
	indices_[j]=iRow;
      } else {
	goodNetwork = 0; // not a network
      }
      break;
      
    case 2:
      if (fabs(elementByColumn[k]-1.0)<1.0e-10) {
	if (fabs(elementByColumn[k+1]+1.0)<1.0e-10) {
	  iRow = row[k];
	  numberRows_ = max(numberRows_,iRow);
	  indices_[j+1]=iRow;
	  iRow = row[k+1];
	  numberRows_ = max(numberRows_,iRow);
	  indices_[j] = iRow;
	} else {
	  goodNetwork = 0; // not a network
	}
      } else if (fabs(elementByColumn[k]+1.0)<1.0e-10) {
	if (fabs(elementByColumn[k+1]-1.0)<1.0e-10) {
	  iRow = row[k];
	  numberRows_ = max(numberRows_,iRow);
	  indices_[j]=iRow;
	  iRow = row[k+1];
	  numberRows_ = max(numberRows_,iRow);
	  indices_[j+1] = iRow;
	} else {
	  goodNetwork = 0; // not a network
	}
      } else {
	goodNetwork = 0; // not a network
      }
      break;

    default:
      goodNetwork = 0; // not a network
      break;
    }
    if (!goodNetwork)
      break;
  }
  if (!goodNetwork) {
    delete [] indices_;
    // put in message
    printf("Not a network - can test if indices_ null\n");
    indices_=NULL;
    numberRows_=0;
    numberColumns_=0;
  } else {
    numberRows_ ++; //  correct
    trueNetwork_ = goodNetwork>0;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNetworkMatrix::~ClpNetworkMatrix ()
{
  delete [] elements_;
  delete [] starts_;
  delete [] lengths_;
  delete [] indices_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNetworkMatrix &
ClpNetworkMatrix::operator=(const ClpNetworkMatrix& rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete [] elements_;
    delete [] starts_;
    delete [] lengths_;
    delete [] indices_;
    elements_ = NULL;
    starts_ = NULL;
    lengths_=NULL;
    indices_=NULL;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    trueNetwork_=rhs.trueNetwork_;
    if (numberColumns_) {
      indices_ = new int [ 2*numberColumns_];
      memcpy(indices_,rhs.indices_,2*numberColumns_*sizeof(int));
    }
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase * ClpNetworkMatrix::clone() const
{
  return new ClpNetworkMatrix(*this);
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase * 
ClpNetworkMatrix::reverseOrderedCopy() const
{
  // count number in each row
  int * tempP = new int [numberRows_];
  int * tempN = new int [numberRows_];
  memset(tempP,0,numberRows_*sizeof(int));
  memset(tempN,0,numberRows_*sizeof(int));
  CoinBigIndex j=0;
  int i;
  for (i=0;i<numberColumns_;i++,j+=2) {
    int iRow = indices_[j];
    tempN[iRow]++;
    iRow = indices_[j+1];
    tempP[iRow]++;
  }
  int * newIndices = new int [2*numberColumns_];
  int * newP = new int [numberRows_+1];
  int * newN = new int[numberRows_];
  int iRow;
  j=0;
  // do starts
  for (iRow=0;iRow<numberRows_;iRow++) {
    newP[iRow]=j;
    j += tempP[iRow];
    tempP[iRow]=newP[iRow];
    newN[iRow] = j;
    j += tempN[iRow];
    tempN[iRow]=newN[iRow];
  }
  newP[numberRows_]=j;
  j=0;
  for (i=0;i<numberColumns_;i++,j+=2) {
    int iRow = indices_[j];
    int put = tempN[iRow];
    newIndices[put++] = i;
    tempN[iRow] = put;
    iRow = indices_[j+1];
    put = tempP[iRow];
    newIndices[put++] = i;
    tempP[iRow] = put;
  }
  delete [] tempP;
  delete [] tempN;
  ClpPlusMinusOneMatrix * newCopy = new ClpPlusMinusOneMatrix();
  newCopy->passInCopy(numberRows_, numberColumns_,
		      false,  newIndices, newP, newN);
  return newCopy;
}
//unscaled versions
void 
ClpNetworkMatrix::times(double scalar,
		   const double * x, double * y) const
{
  int iColumn;
  CoinBigIndex j=0;
  if (trueNetwork_) {
    for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
      double value = scalar*x[iColumn];
      if (value) {
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	y[iRowM] -= value;
	y[iRowP] += value;
      }
    }
  } else {
    // skip negative rows
    for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
      double value = scalar*x[iColumn];
      if (value) {
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	if (iRowM>=0)
	  y[iRowM] -= value;
	if (iRowP>=0)
	  y[iRowP] += value;
      }
    }
  }
}
void 
ClpNetworkMatrix::transposeTimes(double scalar,
				const double * x, double * y) const
{
  int iColumn;
  CoinBigIndex j=0;
  if (trueNetwork_) {
    for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
      double value = y[iColumn];
      int iRowM = indices_[j];
      int iRowP = indices_[j+1];
      value -= scalar*x[iRowM];
      value += scalar*x[iRowP];
      y[iColumn] = value;
    }
  } else {
    // skip negative rows
    for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
      double value = y[iColumn];
      int iRowM = indices_[j];
      int iRowP = indices_[j+1];
      if (iRowM>=0)
	value -= scalar*x[iRowM];
      if (iRowP>=0)
	value += scalar*x[iRowP];
      y[iColumn] = value;
    }
  }
}
void 
ClpNetworkMatrix::times(double scalar,
		       const double * x, double * y,
		       const double * rowScale, 
		       const double * columnScale) const
{
  // we know it is not scaled 
  times(scalar, x, y);
}
void 
ClpNetworkMatrix::transposeTimes( double scalar,
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
ClpNetworkMatrix::transposeTimes(const ClpSimplex * model, double scalar,
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
    if (trueNetwork_) {
      for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
	double value = markVector[iColumn];
	markVector[iColumn]=0.0;
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	value -= scalar*pi[iRowM];
	value += scalar*pi[iRowP];
	if (fabs(value)>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	  array[iColumn]=value;
	}
      }
    } else {
      // skip negative rows
      for (iColumn=0;iColumn<numberColumns_;iColumn++,j+=2) {
	double value = markVector[iColumn];
	markVector[iColumn]=0.0;
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	if (iRowM>=0)
	  value -= scalar*pi[iRowM];
	if (iRowP>=0)
	  value += scalar*pi[iRowP];
	if (fabs(value)>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	  array[iColumn]=value;
	}
      }
    }
    columnArray->setNumElements(numberNonZero);
    y->setNumElements(0);
  } else {
    // do by row
    rowCopy->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  }
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y.
   Squashes small elements and knows about ClpSimplex */
void 
ClpNetworkMatrix::subsetTransposeTimes(const ClpSimplex * model,
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
  if (trueNetwork_) {
    for (jColumn=0;jColumn<numberToDo;jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex j=iColumn<<1;
      int iRowM = indices_[j];
      int iRowP = indices_[j+1];
      value -= pi[iRowM];
      value += pi[iRowP];
      if (fabs(value)>zeroTolerance) {
	index[numberNonZero++]=iColumn;
	array[iColumn]=value;
      }
    }
  } else {
    // skip negative rows
    for (jColumn=0;jColumn<numberToDo;jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex j=iColumn<<1;
      int iRowM = indices_[j];
      int iRowP = indices_[j+1];
      if (iRowM>=0)
	value -= pi[iRowM];
      if (iRowP>=0)
	value += pi[iRowP];
      if (fabs(value)>zeroTolerance) {
	index[numberNonZero++]=iColumn;
	array[iColumn]=value;
      }
    }
  }
}
/* Returns number of elements in basis
   column is basic if entry >=0 */
CoinBigIndex 
ClpNetworkMatrix::numberInBasis(const int * columnIsBasic) const 
{
  int i;
  CoinBigIndex numberElements=0;
  if (trueNetwork_) {
    for (i=0;i<numberColumns_;i++) {
      if (columnIsBasic[i]>=0) 
	numberElements ++;
    }
    numberElements *= 2;
  } else {
    for (i=0;i<numberColumns_;i++) {
      if (columnIsBasic[i]>=0) {
	CoinBigIndex j=i<<1;
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	if (iRowM>=0)
	  numberElements ++;
	if (iRowP>=0)
	  numberElements ++;
      }
    }
  }
  return numberElements;
}
// Fills in basis (Returns number of elements and updates numberBasic)
CoinBigIndex 
ClpNetworkMatrix::fillBasis(const ClpSimplex * model,
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
  if (trueNetwork_) {
    for (i=0;i<numberColumns_;i++) {
      if (columnIsBasic[i]>=0) {
	CoinBigIndex j=i<<1;
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	indexRowU[numberElements]=iRowM;
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements]=-1.0;
	indexRowU[numberElements+1]=iRowP;
	indexColumnU[numberElements+1]=numberBasic;
	elementU[numberElements+1]=1.0;
	numberElements+=2;
	numberBasic++;
      }
    }
  } else {
    for (i=0;i<numberColumns_;i++) {
      if (columnIsBasic[i]>=0) {
	CoinBigIndex j=i<<1;
	int iRowM = indices_[j];
	int iRowP = indices_[j+1];
	if (iRowM>=0) {
	  indexRowU[numberElements]=iRowM;
	  indexColumnU[numberElements]=numberBasic;
	  elementU[numberElements++]=-1.0;
	}
	if (iRowP>=0) {
	  indexRowU[numberElements]=iRowP;
	  indexColumnU[numberElements]=numberBasic;
	  elementU[numberElements++]=1.0;
	}
	numberBasic++;
      }
    }
  }
  return numberElements;
}
/* Unpacks a column into an CoinIndexedvector
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation */
void 
ClpNetworkMatrix::unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn) const 
{
  CoinBigIndex j=iColumn<<1;
  int iRowM = indices_[j];
  int iRowP = indices_[j+1];
  if (iRowM>=0) 
    rowArray->add(iRowM,-1.0);
  if (iRowP>=0) 
    rowArray->add(iRowP,1.0);
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void 
ClpNetworkMatrix::add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn, double multiplier) const 
{
  CoinBigIndex j=iColumn<<1;
  int iRowM = indices_[j];
  int iRowP = indices_[j+1];
  if (iRowM>=0) 
    rowArray->quickAdd(iRowM,-multiplier);
  if (iRowP>=0) 
    rowArray->quickAdd(iRowP,multiplier);
}

// Return a complete CoinPackedMatrix
CoinPackedMatrix * 
ClpNetworkMatrix::getPackedMatrix() const 
{
  return new CoinPackedMatrix(true,numberRows_,numberColumns_,
			      2*numberColumns_,
			      getElements(),indices_,
			      getVectorStarts(),getVectorLengths());

}
/* A vector containing the elements in the packed matrix. Note that there
   might be gaps in this list, entries that do not belong to any
   major-dimension vector. To get the actual elements one should look at
   this vector together with vectorStarts and vectorLengths. */
const double * 
ClpNetworkMatrix::getElements() const 
{
  assert (trueNetwork_); // fix later
  if (!elements_) {
    elements_ = new double [2*numberColumns_];
    int i;
    for (i=0;i<2*numberColumns_;i+=2) {
      elements_[i]=-1.0;
      elements_[i+1]=1.0;
    }
  }
  return elements_;
}

const CoinBigIndex * 
ClpNetworkMatrix::getVectorStarts() const 
{
  assert (trueNetwork_); // fix later
  if (!starts_) {
    starts_ = new int [numberColumns_+1];
    int i;
    for (i=0;i<numberColumns_+1;i++) {
      starts_[i]=i;
    }
  }
  return starts_;
}
/* The lengths of the major-dimension vectors. */
const int * 
ClpNetworkMatrix::getVectorLengths() const
{
  assert (trueNetwork_); // fix later
  if (!lengths_) {
    lengths_ = new int [numberColumns_];
    int i;
    for (i=0;i<numberColumns_;i++) {
      lengths_[i]=2;
    }
  }
  return lengths_;
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void ClpNetworkMatrix::deleteCols(const int numDel, const int * indDel) 
{
  abort();
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void ClpNetworkMatrix::deleteRows(const int numDel, const int * indDel) 
{
  abort();
}
