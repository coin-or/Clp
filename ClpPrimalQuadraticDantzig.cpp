// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <cstdio>

#include "CoinIndexedVector.hpp"

#include "ClpSimplex.hpp"
#include "ClpPrimalQuadraticDantzig.hpp"
#include "ClpSimplexPrimalQuadratic.hpp"
#include "ClpFactorization.hpp"
#include "ClpPackedMatrix.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPrimalQuadraticDantzig::ClpPrimalQuadraticDantzig () 
: ClpPrimalColumnPivot()
{
  type_=1;
  originalNumberRows_ = 0;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPrimalQuadraticDantzig::ClpPrimalQuadraticDantzig (const ClpPrimalQuadraticDantzig & source) 
: ClpPrimalColumnPivot(source)
{  
  originalNumberRows_ = source.originalNumberRows_;

}
// Constructor from model
ClpPrimalQuadraticDantzig::ClpPrimalQuadraticDantzig(
		  ClpSimplexPrimalQuadratic * model, int originalNumberRows)
  :ClpPrimalColumnPivot()
{
  type_=1;
  // Safe as has no extra data
  model_=(ClpSimplex *) model;
  originalNumberRows_ = originalNumberRows;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPrimalQuadraticDantzig::~ClpPrimalQuadraticDantzig ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPrimalQuadraticDantzig &
ClpPrimalQuadraticDantzig::operator=(const ClpPrimalQuadraticDantzig& rhs)
{
  if (this != &rhs) {
    ClpPrimalColumnPivot::operator=(rhs);
    originalNumberRows_ = rhs.originalNumberRows_;
  }
  return *this;
}

// Returns pivot column, -1 if none
int 
ClpPrimalQuadraticDantzig::pivotColumn(CoinIndexedVector * updates,
				    CoinIndexedVector * spareRow1,
				    CoinIndexedVector * spareRow2,
				    CoinIndexedVector * spareColumn1,
				    CoinIndexedVector * spareColumn2)
{
  assert(model_);
  // Find out where stuff is
  int originalNumberColumns = model_->numberRows()-originalNumberRows_;
  int start = model_->numberRows();
  double * solution = model_->solutionRegion();
  //double * lower = model_->lowerRegion();
  //double * upper = model_->upperRegion();
  //double * trueLower = model_->columnLower();
  //double * trueUpper = model_->columnUpper();
  
  //double tolerance=model_->currentPrimalTolerance();

  double bestDj = model_->dualTolerance();
  double bestInfeasibleDj = 10.0*model_->dualTolerance();
  int bestSequence=-1;
  double dualIn=0.0;

  int iSequence;
  assert (!model_->scalingFlag());

  for (iSequence=0;iSequence<originalNumberColumns;iSequence++) {
    int jSequence = iSequence + start;
    double value = solution[jSequence];
    ClpSimplex::Status status = model_->getStatus(iSequence);
    
    switch(status) {
      
    case ClpSimplex::basic:
      if (fabs(value)>bestInfeasibleDj) {
	// infeasible - high priority
	bestDj=COIN_DBL_MAX;
	bestInfeasibleDj = fabs(value);
	bestSequence = jSequence;
	dualIn = value;
      }
      break;
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      abort();
      break;
    case ClpSimplex::atUpperBound:
      if (value>bestDj) {
	bestDj = value;
	bestSequence = iSequence;
	dualIn = value;
      }
      break;
    case ClpSimplex::atLowerBound:
      if (value<-bestDj) {
	bestDj = -value;
	bestSequence = iSequence;
	dualIn = value;
      }
    }
  }
  // And slacks
  // Value of sj is - value of pi
  int firstSlack = model_->numberColumns();
  double * lower = model_->lowerRegion()+model_->numberColumns();
  double * upper = model_->upperRegion()+model_->numberColumns();
  int jSequence;
  for (jSequence=0;jSequence<originalNumberRows_;jSequence++) {
    int iSequence  = jSequence + firstSlack;
    int iPi = jSequence+originalNumberColumns;
    // for slacks either pi zero or row at bound
    // for L rows pi negative okay so choose if positive
    double value = solution[iPi];
    ClpSimplex::Status status = model_->getStatus(iSequence);
    switch(status) {
      
    case ClpSimplex::basic:
      if (fabs(value)>bestDj&&upper[jSequence]>lower[jSequence]) {
	bestDj = fabs(value);
	bestSequence = iPi; //?
	abort();
      }
      break;
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      if (fabs(value)>bestDj) {
	bestDj = fabs(value);
	bestSequence = iSequence;
	dualIn = value;
      }
      break;
    case ClpSimplex::atUpperBound:
      if (value>bestDj) {
	bestDj = value;
	bestSequence = iSequence;
	dualIn = value;
      }
      break;
    case ClpSimplex::atLowerBound:
      if (-value>bestDj) {
	bestDj = -value;
	bestSequence = iSequence;
	dualIn = value;
      }
    }
  }
  if (bestSequence>=0) {
    model_->djRegion()[bestSequence] = dualIn;
#if 0
    // free up depending which bound
    if (solution[iSequence]<trueLower[iSequence]+tolerance) {
      model_->setColumnStatus(iSequence,ClpSimplex::atLowerBound);
      upper[iSequence]=trueUpper[iSequence];
    } else if (solution[iSequence]>trueUpper[iSequence]-tolerance) {
      model_->setColumnStatus(iSequence,ClpSimplex::atUpperBound);
      lower[iSequence]=trueLower[iSequence];
    } else {
      abort();
    }
#endif
  }
  return bestSequence;
}
  
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpPrimalColumnPivot * ClpPrimalQuadraticDantzig::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPrimalQuadraticDantzig(*this);
  } else {
    return new ClpPrimalQuadraticDantzig();
  }
}

