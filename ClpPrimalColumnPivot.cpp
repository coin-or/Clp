// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "ClpSimplex.hpp"
#include "ClpPrimalColumnPivot.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPrimalColumnPivot::ClpPrimalColumnPivot () :
  model_(NULL), 
  type_(-1)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPrimalColumnPivot::ClpPrimalColumnPivot (const ClpPrimalColumnPivot & source) :
  model_(source.model_),
  type_(source.type_)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPrimalColumnPivot::~ClpPrimalColumnPivot ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPrimalColumnPivot &
ClpPrimalColumnPivot::operator=(const ClpPrimalColumnPivot& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    model_ = rhs.model_;
  }
  return *this;
}
void 
ClpPrimalColumnPivot::saveWeights(ClpSimplex * model,int mode)
{
  model_=model;
}
// checks accuracy and may re-initialize (may be empty)

void 
ClpPrimalColumnPivot::updateWeights(CoinIndexedVector * input)
{
}

// Gets rid of all arrays
void 
ClpPrimalColumnPivot::clearArrays()
{
}
