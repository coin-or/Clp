// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "ClpMatrixBase.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpMatrixBase::ClpMatrixBase () :
  type_(-1)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpMatrixBase::ClpMatrixBase (const ClpMatrixBase & source) :
  type_(source.type_)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpMatrixBase::~ClpMatrixBase ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpMatrixBase &
ClpMatrixBase::operator=(const ClpMatrixBase& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
  }
  return *this;
}
// And for scaling - default aborts for when scaling not supported
void 
ClpMatrixBase::times(double scalar,
		     const double * x, double * y,
		     const double * rowScale, 
		     const double * columnScale) const
{
  std::cerr<<"Scaling not supported - ClpMatrixBase"<<std::endl;
  abort();
}
// And for scaling - default aborts for when scaling not supported
void 
ClpMatrixBase::transposeTimes(double scalar,
				const double * x, double * y,
				const double * rowScale, 
				const double * columnScale) const
{
  std::cerr<<"Scaling not supported - ClpMatrixBase"<<std::endl;
  abort();
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given.
   Derived classes need not provide this as it may not always make
   sense */
ClpMatrixBase * 
ClpMatrixBase::subsetClone (
			    int numberRows, const int * whichRows,
			    int numberColumns, const int * whichColumns) const
 

{
  std::cerr<<"subsetClone not supported - ClpMatrixBase"<<std::endl;
  abort();
}
