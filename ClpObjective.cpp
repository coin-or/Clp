// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#include "ClpObjective.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpObjective::ClpObjective () :
  type_(-1)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpObjective::ClpObjective (const ClpObjective & source) :
  type_(source.type_)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpObjective::~ClpObjective ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpObjective &
ClpObjective::operator=(const ClpObjective& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
  }
  return *this;
}
/* Subset clone.  Duplicates are allowed
   and order is as given.
*/
ClpObjective * 
ClpObjective::subsetClone (int numberColumns, 
			   const int * whichColumns) const
{
  std::cerr<<"subsetClone not supported - ClpObjective"<<std::endl;
  abort();
  return (0) ;
}

