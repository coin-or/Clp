// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.



#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpInterior.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpCholeskyDense::ClpCholeskyDense () 
  : ClpCholeskyBase(),
    work_(NULL),
    rowCopy_(NULL)
{
  type_=11;;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpCholeskyDense::ClpCholeskyDense (const ClpCholeskyDense & rhs) 
: ClpCholeskyBase(rhs)
{
  type_=rhs.type_;
  work_ = ClpCopyOfArray(rhs.work_,numberRows_*numberRows_);
  rowCopy_ = rhs.rowCopy_->clone();
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpCholeskyDense::~ClpCholeskyDense ()
{
  delete [] work_;
  delete rowCopy_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpCholeskyDense &
ClpCholeskyDense::operator=(const ClpCholeskyDense& rhs)
{
  if (this != &rhs) {
    ClpCholeskyBase::operator=(rhs);
    delete [] work_;
    work_ = ClpCopyOfArray(rhs.work_,numberRows_*numberRows_);
    delete rowCopy_;
    rowCopy_ = rhs.rowCopy_->clone();
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase * ClpCholeskyDense::clone() const
{
  return new ClpCholeskyDense(*this);
}
/* Orders rows and saves pointer to matrix.and model */
int 
ClpCholeskyDense::order(ClpInterior * model) 
{
  numberRows_ = model->numberRows();
  work_ = new double [numberRows_*numberRows_];
  rowsDropped_ = new char [numberRows_];
  memset(rowsDropped_,0,numberRows_);
  numberRowsDropped_=0;
  model_=model;
  rowCopy_ = model->clpMatrix()->reverseOrderedCopy();
  return 0;
}
//#define CLP_DEBUG
/* Factorize - filling in rowsDropped and returning number dropped */
int 
ClpCholeskyDense::factorize(const double * diagonal, int * rowsDropped) 
{
  int iColumn;
  const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
  const int * columnLength = model_->clpMatrix()->getVectorLengths();
  const int * row = model_->clpMatrix()->getIndices();
  const double * element = model_->clpMatrix()->getElements();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  const double * elementByRow = rowCopy_->getElements();
  int numberColumns=model_->clpMatrix()->getNumCols();
  CoinZeroN(work_,numberRows_*numberRows_);
  int iRow;
  double * work = work_;
  const double * diagonalSlack = diagonal + numberColumns;
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      work[iRow] = diagonalSlack[iRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	CoinBigIndex start=columnStart[iColumn];
	CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	double multiplier = diagonal[iColumn]*elementByRow[k];
	for (CoinBigIndex j=start;j<end;j++) {
	  int jRow=row[j];
	  if (jRow>=iRow&&!rowsDropped_[jRow]) {
	    double value=element[j]*multiplier;
	    work[jRow] += value;
	  }
	}
      } 
    }
    work += numberRows_;
  }
  work = work_;
  // Get rid of any odd copies
  // model_->clpMatrix()->releasePackedMatrix();
#ifdef CLP_DEBUG
  double * save = NULL;
  if (numberRows_<20) {
    save = new double [ numberRows_*numberRows_];
    memcpy(save,work,numberRows_*numberRows_*sizeof(double));
  }
#endif
  int numberDropped=0;
  for (iColumn=0;iColumn<numberRows_;iColumn++) {
    int iRow;
    if (!rowsDropped_[iColumn]) {
      double diagonalValue = work[iColumn];
      if (diagonalValue>pivotTolerance_) {
	diagonalValue = sqrt(diagonalValue);
	work[iColumn]=diagonalValue;
	for (iRow=iColumn+1;iRow<numberRows_;iRow++)
	  work[iRow] /= diagonalValue;
	double * work2 = work;
	for (int jColumn=iColumn+1;jColumn<numberRows_;jColumn++) {
	  work2 += numberRows_;
	  double value = work[jColumn];
	  for (iRow=jColumn;iRow<numberRows_;iRow++)
	  work2[iRow] -= value*work[iRow];
	}
      } else {
	// drop column
	rowsDropped_[iColumn]=-1;
	rowsDropped[numberDropped++]=iColumn;
	numberRowsDropped_++;
	// clean up as this is a debug version
	for (iRow=0;iRow<iColumn;iRow++)
	  work_[iColumn+iRow*numberRows_]=0.0;
	for (iRow=iColumn;iRow<numberRows_;iRow++)
	  work[iRow]=0.0;
      }
    } else {
      // clean up as this is a debug version
      for (iRow=0;iRow<iColumn;iRow++)
	work_[iColumn+iRow*numberRows_]=0.0;
      for (iRow=iColumn;iRow<numberRows_;iRow++)
	work[iRow]=0.0;
    }
    work += numberRows_; // move on pointer
  }
#ifdef CLP_DEBUG
  if (save) {
    double * array = new double [numberRows_];
    int iColumn;
    for (iColumn=0;iColumn<numberRows_;iColumn++) {
      double * s = save;
      int iRow;
      for (iRow=0;iRow<iColumn;iRow++) {
	array[iRow]=s[iColumn];
	s += numberRows_;
      }
      for (iRow=iColumn;iRow<numberRows_;iRow++) {
	array[iRow]=s[iRow];
      }
      solve(array);
      for (iRow=0;iRow<numberRows_;iRow++) {
	if (iRow!=iColumn)
	  assert (fabs(array[iRow])<1.0e-7);
	else
	  assert (fabs(array[iRow]-1.0)<1.0e-7);
      }
    }
    delete [] array;
    delete [] save;
  }
#endif
  return numberDropped;
}
/* Uses factorization to solve. */
void 
ClpCholeskyDense::solve (double * region) 
{
  int iColumn;
  for (iColumn=0;iColumn<numberRows_;iColumn++) {
    if (!rowsDropped_[iColumn]) {
      double value = region[iColumn];
	  int iRow;
      for (iRow=0;iRow<iColumn;iRow++)
	value -= region[iRow]*work_[iColumn+iRow*numberRows_];
      for (iRow=0;iRow<iColumn;iRow++)
	if (rowsDropped_[iRow])
	  assert(!work_[iColumn+iRow*numberRows_]||!region[iRow]);
      region[iColumn]=value/work_[iColumn+iColumn*numberRows_];
    } else {
      region[iColumn]=0.0;
    }
  }
  double * work = work_ + numberRows_*numberRows_;
  for (iColumn=numberRows_-1;iColumn>=0;iColumn--) {
    work -= numberRows_;
    if (!rowsDropped_[iColumn]) {
      double value = region[iColumn];
      for (int iRow=iColumn+1;iRow<numberRows_;iRow++)
	value -= region[iRow]*work[iRow];
      for (int iRow=iColumn+1;iRow<numberRows_;iRow++)
	if (rowsDropped_[iRow])
	  assert(!work[iRow]||!region[iRow]);
      region[iColumn]=value/work[iColumn];
    } else {
      region[iColumn]=0.0;
    }
  }
}
