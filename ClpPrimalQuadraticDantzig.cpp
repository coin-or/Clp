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
  quadraticInfo_=NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPrimalQuadraticDantzig::ClpPrimalQuadraticDantzig (const ClpPrimalQuadraticDantzig & source) 
: ClpPrimalColumnPivot(source)
{  
  quadraticInfo_=source.quadraticInfo_;

}
// Constructor from model
ClpPrimalQuadraticDantzig::ClpPrimalQuadraticDantzig(
		  ClpSimplexPrimalQuadratic * model, 
		  ClpQuadraticInfo * info)
  :ClpPrimalColumnPivot()
{
  type_=1;
  // Safe as has no extra data
  model_=(ClpSimplex *) model;
  quadraticInfo_ = info;
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
    quadraticInfo_ = rhs.quadraticInfo_;
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
  int originalNumberColumns = quadraticInfo_->numberXColumns();
  int start = model_->numberRows();
  double * solution = model_->solutionRegion();
  //double * lower = model_->lowerRegion();
  //double * upper = model_->upperRegion();
  //double * trueLower = model_->columnLower();
  //double * trueUpper = model_->columnUpper();
  
  //double tolerance=model_->currentPrimalTolerance();

  double dualTolerance = model_->dualTolerance()*1000.0;
  double bestDj = dualTolerance;
  double bestInfeasibleDj = 10.0*dualTolerance;
  int bestSequence=-1;
  double dualIn=0.0;

  int iSequence;
  assert (!model_->scalingFlag());
  const double * pi = solution+quadraticInfo_->numberXColumns();
  // Matrix for linear stuff
  CoinPackedMatrix * matrix = model_->matrix();
  const int * row = matrix->getIndices();
  const int * columnStart = matrix->getVectorStarts();
  const double * element =  matrix->getElements();
  const int * columnLength = matrix->getVectorLengths();
  const int * which = quadraticInfo_->quadraticSequence();
  const double * objective = quadraticInfo_->linearObjective();
  const double * djWeight = quadraticInfo_->djWeight();
  for (iSequence=0;iSequence<originalNumberColumns;iSequence++) {
    if (model_->flagged(iSequence))
      continue;
    int jSequence = which[iSequence];
    double value;
    if (jSequence>=0) {
      jSequence += start;
      value = solution[jSequence];
    } else {
      value=objective[iSequence];
      int j;
      for (j=columnStart[iSequence];j<columnStart[iSequence]+columnLength[iSequence]; j++) {
	int iRow = row[j];
	value -= element[j]*pi[iRow];
      }
    }
    double value2 = value*djWeight[iSequence];
    if (fabs(value2)>dualTolerance)
      value=value2;
    else if (value<-dualTolerance)
      value = -1.001*dualTolerance;
    else if (value>dualTolerance)
      value = 1.001*dualTolerance;
    if (djWeight[iSequence]<1.0e-6)
      value=value2;
    ClpSimplex::Status status = model_->getStatus(iSequence);
    
    switch(status) {
      
    case ClpSimplex::basic:
      if (fabs(value)>bestInfeasibleDj) {
	// infeasible - high priority
	bestDj=COIN_DBL_MAX;
	bestInfeasibleDj = fabs(value);
	bestSequence = jSequence;
	dualIn = value;
        abort();
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
  int originalNumberRows = quadraticInfo_->numberXRows();
  for (jSequence=0;jSequence<originalNumberRows;jSequence++) {
    int iSequence  = jSequence + firstSlack;
    if (model_->flagged(iSequence))
      continue;
    int iPi = jSequence+originalNumberColumns;
    // for slacks either pi zero or row at bound
    // for L rows pi negative okay so choose if positive
    double value = solution[iPi];
    double value2 = value*djWeight[iSequence];
    if (fabs(value2)>dualTolerance)
      value=value2;
    else if (value<-dualTolerance)
      value = -1.001*dualTolerance;
    else if (value>dualTolerance)
      value = 1.001*dualTolerance;
    if (djWeight[iSequence]<1.0e-6)
      value=value2;
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

