// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved. 

#if defined(_MSC_VER) 
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "MyEventHandler.hpp"


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler () 
  : ClpEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler (const MyEventHandler & rhs) 
: ClpEventHandler(rhs)
{  
}

// Constructor with pointer to model
MyEventHandler::MyEventHandler(ClpSimplex * model)
  : ClpEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyEventHandler::~MyEventHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyEventHandler &
MyEventHandler::operator=(const MyEventHandler& rhs)
{
  if (this != &rhs) {
    ClpEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpEventHandler * MyEventHandler::clone() const
{
  return new MyEventHandler(*this);
}

int 
MyEventHandler::event(Event whichEvent)
{
  if (whichEvent==endOfValuesPass)
    return 0; // say optimal
  else
    return -1; // carry on
}
