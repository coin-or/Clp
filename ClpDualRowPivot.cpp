// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "ClpSimplex.hpp"
#include "ClpDualRowPivot.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpDualRowPivot::ClpDualRowPivot () :
  model_(NULL), 
  type_(-1)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpDualRowPivot::ClpDualRowPivot (const ClpDualRowPivot & source) :
  model_(source.model_),
  type_(source.type_)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpDualRowPivot::~ClpDualRowPivot ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpDualRowPivot &
ClpDualRowPivot::operator=(const ClpDualRowPivot& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    model_ = rhs.model_;
  }
  return *this;
}
void 
ClpDualRowPivot::saveWeights(ClpSimplex * model,int mode)
{
  model_=model;
}
// checks accuracy and may re-initialize (may be empty)
void 
ClpDualRowPivot::checkAccuracy()
{
}

void 
ClpDualRowPivot::updateWeights(CoinIndexedVector * input,
			       CoinIndexedVector * spare,
			       CoinIndexedVector * updatedColumn)
{
}
void 
ClpDualRowPivot::unrollWeights()
{
}
