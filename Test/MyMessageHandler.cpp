// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved. 

#if defined(_MSC_VER) 
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <cstdio>

#include "ClpSimplex.hpp"
#include "MyMessageHandler.hpp"
#include "ClpMessage.hpp"


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyMessageHandler::MyMessageHandler () 
  : CoinMessageHandler(),
    model_(NULL),
    feasibleExtremePoints_(),
    iterationNumber_(-1)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyMessageHandler::MyMessageHandler (const MyMessageHandler & rhs) 
: CoinMessageHandler(rhs),
    model_(rhs.model_),
    feasibleExtremePoints_(rhs.feasibleExtremePoints_),
    iterationNumber_(rhs.iterationNumber_)
{  
}

MyMessageHandler::MyMessageHandler (const CoinMessageHandler & rhs) 
  : CoinMessageHandler(),
    model_(NULL),
    feasibleExtremePoints_(),
    iterationNumber_(-1)
{  
}

// Constructor with pointer to model
MyMessageHandler::MyMessageHandler(ClpSimplex * model,
               FILE * userPointer)
  : CoinMessageHandler(),
    model_(model),
    feasibleExtremePoints_(),
    iterationNumber_(-1)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyMessageHandler::~MyMessageHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyMessageHandler &
MyMessageHandler::operator=(const MyMessageHandler& rhs)
{
  if (this != &rhs) {
    CoinMessageHandler::operator=(rhs);
    model_ = rhs.model_;
    feasibleExtremePoints_ = rhs.feasibleExtremePoints_;
    iterationNumber_ = rhs.iterationNumber_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CoinMessageHandler * MyMessageHandler::clone() const
{
  return new MyMessageHandler(*this);
}

int 
MyMessageHandler::print()
{
  if (currentSource()=="Clp") {
    if (currentMessage().externalNumber()==5) {
      // Are we feasible
      // We can pick up two ways - check same
      assert (model_->numberPrimalInfeasibilities()==
        intValue(1));
      if (!intValue(1)&&intValue(0)!=iterationNumber_) {
        iterationNumber_ = intValue(0);
        // Column solution
        int numberColumns = model_->numberColumns();
        const double * solution = model_->solutionRegion(1);

        // Create vector to contain solution
        StdVectorDouble feasibleExtremePoint;

        if (!model_->columnScale()) {
          // No scaling
          for (int i=0;i<numberColumns;i++)
            feasibleExtremePoint.push_back(solution[i]);
        } else {
          // scaled
          const double * columnScale = model_->columnScale();
          for (int i=0;i<numberColumns;i++)
            feasibleExtremePoint.push_back(solution[i]*columnScale[i]);
        }
        // Save solution
        feasibleExtremePoints_.push_front(feasibleExtremePoint);

        // Want maximum of 10 solutions, so if more then 10 get rid of oldest
        int numExtremePointsSaved = feasibleExtremePoints_.size();
        if ( numExtremePointsSaved>=10 ) {
          feasibleExtremePoints_.pop_back();
          assert( feasibleExtremePoints_.size() == 
		  (unsigned int) numExtremePointsSaved-1 );
        };

      }
    }
  }
  return CoinMessageHandler::print();
}
const ClpSimplex *
MyMessageHandler::model() const
{
  return model_;
}
void 
MyMessageHandler::setModel(ClpSimplex * model)
{
  model_ = model;
}

const std::deque<StdVectorDouble> & MyMessageHandler::getFeasibleExtremePoints() const
{ 
  return feasibleExtremePoints_; 
}
void MyMessageHandler::clearFeasibleExtremePoints()
{
  feasibleExtremePoints_.clear();
}
