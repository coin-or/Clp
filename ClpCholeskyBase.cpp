// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "ClpCholeskyBase.hpp"
#include "ClpInterior.hpp"
#include "ClpHelperFunctions.hpp"
#include "CoinHelperFunctions.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpCholeskyBase::ClpCholeskyBase () :
  type_(-1),
  pivotTolerance_(1.0e-14),
  zeroTolerance_(1.0e-17),
  choleskyCondition_(0.0),
  model_(NULL),
  numberTrials_(),
  numberRows_(0),
  status_(0),
  rowsDropped_(NULL),
  permuteIn_(NULL),
  permuteOut_(NULL),
  numberRowsDropped_(0)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpCholeskyBase::ClpCholeskyBase (const ClpCholeskyBase & rhs) :
  type_(rhs.type_),
  pivotTolerance_(rhs.pivotTolerance_),
  zeroTolerance_(rhs.zeroTolerance_),
  choleskyCondition_(rhs.choleskyCondition_),
  model_(rhs.model_),
  numberTrials_(rhs.numberTrials_),
  numberRows_(rhs.numberRows_),
  status_(rhs.status_),
  numberRowsDropped_(rhs.numberRowsDropped_)
{  
  rowsDropped_ = ClpCopyOfArray(rhs.rowsDropped_,numberRows_);
  permuteIn_ = ClpCopyOfArray(rhs.permuteIn_,numberRows_);
  permuteOut_ = ClpCopyOfArray(rhs.permuteOut_,numberRows_);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpCholeskyBase::~ClpCholeskyBase ()
{
  delete [] rowsDropped_;
  delete [] permuteIn_;
  delete [] permuteOut_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpCholeskyBase &
ClpCholeskyBase::operator=(const ClpCholeskyBase& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    pivotTolerance_ = rhs.pivotTolerance_;
    zeroTolerance_ = rhs.zeroTolerance_;
    choleskyCondition_ = rhs.choleskyCondition_;
    model_ = rhs.model_;
    numberTrials_ = rhs.numberTrials_;
    numberRows_ = rhs.numberRows_;
    status_ = rhs.status_;
    numberRowsDropped_ = rhs.numberRowsDropped_;
    delete [] rowsDropped_;
    delete [] permuteIn_;
    delete [] permuteOut_;
    rowsDropped_ = ClpCopyOfArray(rhs.rowsDropped_,numberRows_);
    permuteIn_ = ClpCopyOfArray(rhs.permuteIn_,numberRows_);
    permuteOut_ = ClpCopyOfArray(rhs.permuteOut_,numberRows_);
  }
  return *this;
}
// reset numberRowsDropped and rowsDropped.
void 
ClpCholeskyBase::resetRowsDropped()
{
  numberRowsDropped_=0;
  memset(rowsDropped_,0,numberRows_);
}
/* Uses factorization to solve. - given as if KKT.
   region1 is rows+columns, region2 is rows */
void 
ClpCholeskyBase::solveKKT (double * region1, double * region2, const double * diagonal,
			   double diagonalScaleFactor)
{
  int iColumn;
  int numberColumns = model_->numberColumns();
  int numberTotal = numberRows_+numberColumns;
  double * region1Save = new double[numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    region1[iColumn] *= diagonal[iColumn];
    region1Save[iColumn]=region1[iColumn];
  }
  multiplyAdd(region1+numberColumns,numberRows_,-1.0,region2,1.0);
  model_->clpMatrix()->times(1.0,region1,region2);
  double maximumRHS = maximumAbsElement(region2,numberRows_);
  double scale=1.0;
  double unscale=1.0;
  if (maximumRHS>1.0e-30) {
    if (maximumRHS<=0.5) {
      double factor=2.0;
      while (maximumRHS<=0.5) {
	maximumRHS*=factor;
	scale*=factor;
      } /* endwhile */
    } else if (maximumRHS>=2.0&&maximumRHS<=COIN_DBL_MAX) {
      double factor=0.5;
      while (maximumRHS>=2.0) {
	maximumRHS*=factor;
	scale*=factor;
      } /* endwhile */
    } 
    unscale=diagonalScaleFactor/scale;
  } else {
    //effectively zero
    scale=0.0;
    unscale=0.0;
  } 
  multiplyAdd(NULL,numberRows_,0.0,region2,scale);
  solve(region2);
  multiplyAdd(NULL,numberRows_,0.0,region2,unscale);
  multiplyAdd(region2,numberRows_,-1.0,region1+numberColumns,0.0);
  CoinZeroN(region1,numberColumns);
  model_->clpMatrix()->transposeTimes(1.0,region2,region1);
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    region1[iColumn] = region1[iColumn]*diagonal[iColumn]-region1Save[iColumn];
  delete [] region1Save;
}
