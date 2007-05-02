// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#include "ClpConstraint.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpConstraint::ClpConstraint () :
  lastGradient_(NULL),
  functionValue_(0.0),
  type_(-1),
  rowNumber_(-1)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpConstraint::ClpConstraint (const ClpConstraint & source) :
  lastGradient_(NULL),
  functionValue_(source.functionValue_),
  type_(source.type_),
  rowNumber_(source.rowNumber_)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpConstraint::~ClpConstraint ()
{
  delete [] lastGradient_;

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpConstraint &
ClpConstraint::operator=(const ClpConstraint& rhs)
{
  if (this != &rhs) {
    functionValue_ = rhs.functionValue_;
    type_ = rhs.type_;
    rowNumber_= rhs.rowNumber_;
    delete [] lastGradient_;
    lastGradient_ = NULL;
  }
  return *this;
}

