// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "ClpCholeskyBase.hpp"
#include "ClpInterior.hpp"

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
