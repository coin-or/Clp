// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "ClpMatrixBase.hpp"
#include "ClpSimplex.hpp"

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
  if (rowScale) {
    std::cerr<<"Scaling not supported - ClpMatrixBase"<<std::endl;
    abort();
  } else {
    times(scalar,x,y);
  }
}
// And for scaling - default aborts for when scaling not supported
void 
ClpMatrixBase::transposeTimes(double scalar,
				const double * x, double * y,
				const double * rowScale, 
				const double * columnScale) const
{
  if (rowScale) {
    std::cerr<<"Scaling not supported - ClpMatrixBase"<<std::endl;
    abort();
  } else {
    transposeTimes(scalar,x,y);
  }
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
  return NULL;
}
/* Given positive integer weights for each row fills in sum of weights
   for each column (and slack).
   Returns weights vector
   Default returns vector of ones
*/
CoinBigIndex * 
ClpMatrixBase::dubiousWeights(const ClpSimplex * model,int * inputWeights) const
{
  int number = model->numberRows()+model->numberColumns();
  CoinBigIndex * weights = new CoinBigIndex[number];
  int i;
  for (i=0;i<number;i++)
    weights[i]=1;
  return weights;
}
// Append Columns
void 
ClpMatrixBase::appendCols(int number, const CoinPackedVectorBase * const * columns)
{
  std::cerr<<"appendCols not supported - ClpMatrixBase"<<std::endl;
  abort();
}
// Append Rows
void 
ClpMatrixBase::appendRows(int number, const CoinPackedVectorBase * const * rows)
{
  std::cerr<<"appendRows not supported - ClpMatrixBase"<<std::endl;
  abort();
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void 
ClpMatrixBase::rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive)
{
  smallestNegative=0.0;
  largestNegative=0.0;
  smallestPositive=0.0;
  largestPositive=0.0;
}
// Says whether it can do partial pricing
bool 
ClpMatrixBase::canDoPartialPricing() const
{
  return false; // default is no
}
// Partial pricing 
void 
ClpMatrixBase::partialPricing(ClpSimplex * model, int start, int end,
			      int & bestSequence, int & numberWanted)
{
  std::cerr<<"partialPricing not supported - ClpMatrixBase"<<std::endl;
  abort();
}
