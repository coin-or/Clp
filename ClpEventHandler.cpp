// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include "ClpEventHandler.hpp"
#include "ClpSimplex.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpEventHandler::ClpEventHandler (ClpSimplex * model) :
  model_(model)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpEventHandler::ClpEventHandler (const ClpEventHandler & rhs) 
  : model_(rhs.model_)
{  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpEventHandler::~ClpEventHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpEventHandler &
ClpEventHandler::operator=(const ClpEventHandler& rhs)
{
  if (this != &rhs) {
    model_ = rhs.model_;
  }
  return *this;
}
// Clone
ClpEventHandler * 
ClpEventHandler::clone() const
{
  return new ClpEventHandler(*this);
}
// Event
int
ClpEventHandler::event(Event whichEvent)
{
  return -1; // do nothing
}
/* set model. */
void 
ClpEventHandler::setSimplex(ClpSimplex * model)
{
  model_= model;
}


