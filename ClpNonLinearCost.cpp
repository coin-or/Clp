// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>

#include "ClpNonLinearCost.hpp"
#include "ClpSimplex.hpp"
#include "OsiIndexedVector.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNonLinearCost::ClpNonLinearCost () :
  numberRows_(0),
  numberColumns_(0),
  start_(NULL),
  whichRange_(NULL),
  offset_(NULL),
  lower_(NULL),
  cost_(NULL),
  model_(NULL),
  numberInfeasibilities_(-1),
  changeCost_(0.0),
  largestInfeasibility_(0.0),
  sumInfeasibilities_(0.0),
  convex_(true)
{

}
/* Constructor from simplex.
   This will just set up wasteful arrays for linear, but
   later may do dual analysis and even finding duplicate columns 
*/
ClpNonLinearCost::ClpNonLinearCost ( ClpSimplex * model)
{
  model_ = model;
  numberRows_ = model_->numberRows();
  numberColumns_ = model_->numberColumns();
  int numberTotal = numberRows_+numberColumns_;
  convex_ = true;
  start_ = new int [numberTotal+1];
  whichRange_ = new int [numberTotal];
  offset_ = new int [numberTotal];
  memset(offset_,0,numberTotal*sizeof(int));

  numberInfeasibilities_=0;
  changeCost_=0.0;
  double infeasibilityCost = model_->infeasibilityCost();
  sumInfeasibilities_=0.0;
  largestInfeasibility_=0.0;

  // First see how much space we need
  int put=0;

  int iSequence;
  double * upper = model_->upperRegion();
  double * lower = model_->lowerRegion();
  double * cost = model_->costRegion();

  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    if (upper[iSequence]<1.0e20)
      put++;
    put += 3;
  }

  lower_ = new double [put];
  cost_ = new double [put];

  put=0;

  start_[0]=0;

  for (iSequence=0;iSequence<numberTotal;iSequence++) {

    lower_[put] = -DBL_MAX;
    cost_[put++] = cost[iSequence]-infeasibilityCost;
    whichRange_[iSequence]=put;
    lower_[put] = lower[iSequence];
    cost_[put++] = cost[iSequence];
    lower_[put] = upper[iSequence];
    cost_[put++] = cost[iSequence]+infeasibilityCost;
    if (upper[iSequence]<1.0e20) {
      lower_[put] = DBL_MAX;
      cost_[put++] = 1.0e50;
    }
    start_[iSequence+1]=put;
  }

}

ClpNonLinearCost::ClpNonLinearCost(ClpSimplex * model,const int * starts,
		   const double * lowerNon, const double * costNon)
{

  // what about scaling? - only try without it initially
  assert(!model->scalingFlag());
  model_ = model;
  numberRows_ = model_->numberRows();
  numberColumns_ = model_->numberColumns();
  int numberTotal = numberRows_+numberColumns_;
  convex_ = true;
  start_ = new int [numberTotal+1];
  whichRange_ = new int [numberTotal];
  offset_ = new int [numberTotal];
  memset(offset_,0,numberTotal*sizeof(int));

  numberInfeasibilities_=0;
  changeCost_=0.0;
  double infeasibilityCost = model_->infeasibilityCost();
  largestInfeasibility_=0.0;
  sumInfeasibilities_=0.0;

  int iSequence;
  double * upper = model_->upperRegion();
  double * lower = model_->lowerRegion();
  double * cost = model_->costRegion();
  // First see how much space we need - we know column part
  // but not infeasibilities - see how much extra
  // may be over-estimate
  int put=starts[numberColumns_]-numberColumns_;

  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    if (lower[iSequence]>-1.0e20)
      if (upper[iSequence]<1.0e20)
	put++;
    put += 3;
  }

  lower_ = new double [put];
  cost_ = new double [put];

  // now fill in 
  put=0;

  start_[0]=0;

  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    lower_[put] = -DBL_MAX;
    cost_[put++] = cost[iSequence]-infeasibilityCost;
    whichRange_[iSequence]=put;
    double thisCost;
    if (iSequence>=numberColumns_) {
      // rows
      lower_[put] = lower[iSequence];
      cost_[put++] = cost[iSequence];
      thisCost = cost[iSequence];
    } else {
      // columns - move costs and see if convex
      int iIndex = starts[iSequence];
      int end = starts[iSequence+1];
      assert (fabs(lower[iSequence]-lowerNon[iIndex])<1.0e-8);
      thisCost = -DBL_MAX;
      for (;iIndex<end;iIndex++) {
	if (lowerNon[iIndex]<upper[iSequence]-1.0e-8) {
	  lower_[put] = lowerNon[iIndex];
	  cost_[put++] = costNon[iIndex];
	  // check convexity
	  if (costNon[iIndex]<thisCost-1.0e-12)
	    convex_ = false;
	  thisCost = costNon[iIndex];
	} else {
	  break;
	}
      }
    }
    lower_[put] = upper[iSequence];
    cost_[put++] = thisCost+infeasibilityCost;
    if (upper[iSequence]<1.0e20) {
      lower_[put] = DBL_MAX;
      cost_[put++] = 1.0e50;
    }
    start_[iSequence+1]=put;
  }
  // can't handle non-convex at present
  assert(convex_);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNonLinearCost::ClpNonLinearCost (const ClpNonLinearCost & rhs) :
  numberRows_(rhs.numberRows_),
  numberColumns_(rhs.numberColumns_),
  start_(NULL),
  whichRange_(NULL),
  offset_(NULL),
  lower_(NULL),
  cost_(NULL),
  model_(NULL),
  numberInfeasibilities_(-1),
  changeCost_(0.0),
  largestInfeasibility_(0.0),
  sumInfeasibilities_(0.0),
  convex_(true)
{  
  if (numberRows_) {
    int numberTotal = numberRows_+numberColumns_;
    start_ = new int [numberTotal+1];
    memcpy(start_,rhs.start_,(numberTotal+1)*sizeof(int));
    whichRange_ = new int [numberTotal];
    memcpy(whichRange_,rhs.whichRange_,numberTotal*sizeof(int));
    offset_ = new int [numberTotal];
    memcpy(offset_,rhs.offset_,numberTotal*sizeof(int));
    int numberEntries = start_[numberTotal];
    lower_ = new double [numberEntries];
    memcpy(lower_,rhs.lower_,numberEntries*sizeof(double));
    cost_ = new double [numberEntries];
    memcpy(cost_,rhs.cost_,numberEntries*sizeof(double));
    model_ = rhs.model_;
    numberInfeasibilities_=rhs.numberInfeasibilities_;
    changeCost_ = rhs.changeCost_;
    largestInfeasibility_ = rhs.largestInfeasibility_;
    sumInfeasibilities_ = rhs.sumInfeasibilities_;
    convex_ = rhs.convex_;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNonLinearCost::~ClpNonLinearCost ()
{
  delete [] start_;
  delete [] whichRange_;
  delete [] offset_;
  delete [] lower_;
  delete [] cost_;

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNonLinearCost &
ClpNonLinearCost::operator=(const ClpNonLinearCost& rhs)
{
  if (this != &rhs) {
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    delete [] start_;
    delete [] whichRange_;
    delete [] offset_;
    delete [] lower_;
    delete []cost_;
    start_ = NULL;
    whichRange_ = NULL;
    lower_ = NULL;
    cost_= NULL;
    if (numberRows_) {
      int numberTotal = numberRows_+numberColumns_;
      start_ = new int [numberTotal+1];
      memcpy(start_,rhs.start_,(numberTotal+1)*sizeof(int));
      whichRange_ = new int [numberTotal];
      memcpy(whichRange_,rhs.whichRange_,numberTotal*sizeof(int));
      offset_ = new int [numberTotal];
      memcpy(offset_,rhs.offset_,numberTotal*sizeof(int));
      int numberEntries = start_[numberTotal];
      lower_ = new double [numberEntries];
      memcpy(lower_,rhs.lower_,numberEntries*sizeof(double));
      cost_ = new double [numberEntries];
      memcpy(cost_,rhs.cost_,numberEntries*sizeof(double));
    }
    model_ = rhs.model_;
    numberInfeasibilities_=rhs.numberInfeasibilities_;
    changeCost_ = rhs.changeCost_;
    largestInfeasibility_ = rhs.largestInfeasibility_;
    sumInfeasibilities_ = rhs.sumInfeasibilities_;
    convex_ = rhs.convex_;
  }
  return *this;
}

// Changes infeasible costs and computes number and cost of infeas
// We will need to re-think objective offsets later
// We will also need a 2 bit per variable array for some
// purpose which will come to me later
void 
ClpNonLinearCost::checkInfeasibilities(bool toNearest)
{
  numberInfeasibilities_=0;
  double infeasibilityCost = model_->infeasibilityCost();
  changeCost_=0.0;
  largestInfeasibility_ = 0.0;
  sumInfeasibilities_ = 0.0;
  double primalTolerance = model_->currentPrimalTolerance();
  
  int iSequence;
  double * solution = model_->solutionRegion();
  double * upper = model_->upperRegion();
  double * lower = model_->lowerRegion();
  double * cost = model_->costRegion();
    
  // nonbasic should be at a valid bound
  for (iSequence=0;iSequence<numberColumns_+numberRows_;iSequence++) {
    double lowerValue=lower[iSequence];
    double upperValue=upper[iSequence];
    double value=solution[iSequence];
    int iRange;
    // correct costs
    
    if (lowerValue>-1.0e20) {
      iRange=start_[iSequence];
      cost_[iRange] = cost[iSequence]-infeasibilityCost;
    }
    if (upperValue<1.0e20) {
      iRange=start_[iSequence+1]-2;
      cost_[iRange] = cost[iSequence]+infeasibilityCost;
    }
    // get correct place
    int start = start_[iSequence];
    int end = start_[iSequence+1]-1;
    for (iRange=start; iRange<end;iRange++) {
      if (value<lower_[iRange+1]+primalTolerance) {
	// put in better range if infeasible
	if (value>=lower_[iRange+1]-primalTolerance&&iRange==start) 
	  iRange++;
	whichRange_[iSequence]=iRange;
	break;
      }
    }
    assert(iRange<end);
    switch(model_->getStatus(iSequence)) {
      
    case ClpSimplex::basic:
    case ClpSimplex::superBasic:
      // iRange is in correct place
      // slot in here
      if (value<lower[iSequence]-primalTolerance) {
	value = lower[iSequence]-value;
	sumInfeasibilities_ += value;
	largestInfeasibility_ = max(largestInfeasibility_,value);
	changeCost_ -= lower[iSequence]*
	  (cost_[iRange]-cost[iSequence]);
	numberInfeasibilities_++;
      } else if (value>upper[iSequence]+primalTolerance) {
	value = value-upper[iSequence];
	sumInfeasibilities_ += value;
	largestInfeasibility_ = max(largestInfeasibility_,value);
	changeCost_ -= upper[iSequence]*
	  (cost_[iRange]-cost[iSequence]);
	numberInfeasibilities_++;
      }
      lower[iSequence] = lower_[iRange];
      upper[iSequence] = lower_[iRange+1];
      cost[iSequence] = cost_[iRange];
      break;
    case ClpSimplex::isFree:
      if (toNearest)
	solution[iSequence] = 0.0;
      break;
    case ClpSimplex::atUpperBound:
      if (!toNearest) {
	assert(fabs(value-upperValue)<=primalTolerance*1.0001) ;
      } else {
	solution[iSequence] = upperValue;
      }
      break;
    case ClpSimplex::atLowerBound:
      if (!toNearest) {
	assert(fabs(value-lowerValue)<=primalTolerance*1.0001); 
      } else {
	solution[iSequence] = lowerValue;
      }
      break;
    }
  }
}
/* Goes through one bound for each variable.
   If array[i]*multiplier>0 goes down, otherwise up.
   The indices are row indices and need converting to sequences
*/
void 
ClpNonLinearCost::goThru(int numberInArray, double multiplier,
	      const int * index, const double * array,
			 double * rhs)
{
  assert (model_!=NULL);
  const int * pivotVariable = model_->pivotVariable();
  int i;
  for (i=0;i<numberInArray;i++) {
    int iRow = index[i];
    int iPivot = pivotVariable[iRow];
    double alpha = multiplier*array[iRow];
    // get where in bound sequence
    int iRange = whichRange_[iPivot];
    iRange += offset_[iPivot]; //add temporary bias
    double value = model_->solution(iPivot);
    if (alpha>0.0) {
      // down one
      iRange--;
      assert(iRange>=start_[iPivot]);
      rhs[iRow] = value - lower_[iRange];
    } else {
      // up one
      iRange++;
      assert(iRange<start_[iPivot+1]-1);
      rhs[iRow] = lower_[iRange+1] - value;
    }
    offset_[iPivot] = iRange - whichRange_[iPivot];
  }
}
/* Takes off last iteration (i.e. offsets closer to 0)
 */
void 
ClpNonLinearCost::goBack(int numberInArray, const int * index, 
	      double * rhs)
{
  assert (model_!=NULL);
  const int * pivotVariable = model_->pivotVariable();
  int i;
  for (i=0;i<numberInArray;i++) {
    int iRow = index[i];
    int iPivot = pivotVariable[iRow];
    // get where in bound sequence
    int iRange = whichRange_[iPivot];
    bool down;
    // get closer to original
    if (offset_[iPivot]>0) {
      offset_[iPivot]--;
      assert (offset_[iPivot]>=0);
      down = false;
    } else {
      offset_[iPivot]++;
      assert (offset_[iPivot]<=0);
      down = true;
    }
    iRange += offset_[iPivot]; //add temporary bias
    double value = model_->solution(iPivot);
    if (down) {
      // down one
      assert(iRange>=start_[iPivot]);
      rhs[iRow] = value - lower_[iRange];
    } else {
      // up one
      assert(iRange<start_[iPivot+1]-1);
      rhs[iRow] = lower_[iRange+1] - value;
    }
  }
}
void 
ClpNonLinearCost::goBackAll(const OsiIndexedVector * update)
{
  assert (model_!=NULL);
  const int * pivotVariable = model_->pivotVariable();
  int i;
  int number = update->getNumElements();
  const int * index = update->getIndices();
  for (i=0;i<number;i++) {
    int iRow = index[i];
    int iPivot = pivotVariable[iRow];
    offset_[iPivot]=0;
  }
#ifdef CLP_DEBUG
  for (i=0;i<numberRows_+numberColumns_;i++) 
    assert(!offset_[i]);
#endif
}
void 
ClpNonLinearCost::checkInfeasibilities(int numberInArray, const int * index)
{
  assert (model_!=NULL);
  double primalTolerance = model_->currentPrimalTolerance();
  const int * pivotVariable = model_->pivotVariable();
  int i;
  for (i=0;i<numberInArray;i++) {
    int iRow = index[i];
    int iPivot = pivotVariable[iRow];
    // get where in bound sequence
    int iRange;
    double value = model_->solution(iPivot);
    int start = start_[iPivot];
    int end = start_[iPivot+1]-1;
    for (iRange=start; iRange<end;iRange++) {
      if (value<lower_[iRange+1]+primalTolerance) {
	// put in better range
	if (value>=lower_[iRange+1]-primalTolerance&&iRange==start) 
	  iRange++;
	break;
      }
    }
    assert(iRange<end);
    assert(model_->getStatus(iPivot)==ClpSimplex::basic);
    double & lower = model_->lowerAddress(iPivot);
    double & upper = model_->upperAddress(iPivot);
    double & cost = model_->costAddress(iPivot);
    whichRange_[iPivot]=iRange;
    lower = lower_[iRange];
    upper = lower_[iRange+1];
    cost = cost_[iRange];
  }
}
/* Puts back correct infeasible costs for each variable
   The input indices are row indices and need converting to sequences
   for costs.
   On input array is empty (but indices exist).  On exit just
   changed costs will be stored as normal OsiIndexedVector
*/
void 
ClpNonLinearCost::checkChanged(int numberInArray, OsiIndexedVector * update)
{
  assert (model_!=NULL);
  double primalTolerance = model_->currentPrimalTolerance();
  const int * pivotVariable = model_->pivotVariable();
  int number=0;
  int * index = update->getIndices();
  double * work = update->denseVector();
  int i;
  for (i=0;i<numberInArray;i++) {
    int iRow = index[i];
    int iPivot = pivotVariable[iRow];
    // get where in bound sequence
    int iRange;
    double value = model_->solution(iPivot);
    int start = start_[iPivot];
    int end = start_[iPivot+1]-1;
    for (iRange=start; iRange<end;iRange++) {
      if (value<lower_[iRange+1]+primalTolerance) {
	// put in better range
	if (value>=lower_[iRange+1]-primalTolerance&&iRange==start) 
	  iRange++;
	break;
      }
    }
    assert(iRange<end);
    assert(model_->getStatus(iPivot)==ClpSimplex::basic);
    int jRange = whichRange_[iPivot];
    if (iRange!=jRange) {
      // changed
      work[iRow] = cost_[jRange]-cost_[iRange];
      index[number++]=iRow;
      double & lower = model_->lowerAddress(iPivot);
      double & upper = model_->upperAddress(iPivot);
      double & cost = model_->costAddress(iPivot);
      whichRange_[iPivot]=iRange;
      lower = lower_[iRange];
      upper = lower_[iRange+1];
      cost = cost_[iRange];
    }
  }
  update->setNumElements(number);
}
/* Sets bounds and cost for one variable - returns change in cost*/
double 
ClpNonLinearCost::setOne(int iPivot, double value)
{
  assert (model_!=NULL);
  double primalTolerance = model_->currentPrimalTolerance();
  // get where in bound sequence
  int iRange;
  int start = start_[iPivot];
  int end = start_[iPivot+1]-1;
  for (iRange=start; iRange<end;iRange++) {
    if (value<lower_[iRange+1]+primalTolerance) {
      // put in better range
      if (value>=lower_[iRange+1]-primalTolerance&&iRange==start) 
	iRange++;
      break;
    }
  }
  assert(iRange<end);
  whichRange_[iPivot]=iRange;
  double & lower = model_->lowerAddress(iPivot);
  double & upper = model_->upperAddress(iPivot);
  double & cost = model_->costAddress(iPivot);
  lower = lower_[iRange];
  upper = lower_[iRange+1];
  switch(model_->getStatus(iPivot)) {
      
  case ClpSimplex::basic:
  case ClpSimplex::superBasic:
  case ClpSimplex::isFree:
    break;
  case ClpSimplex::atUpperBound:
  case ClpSimplex::atLowerBound:
    // set correctly
    if (fabs(value-lower)<=primalTolerance*1.001) 
      model_->setStatus(iPivot,ClpSimplex::atLowerBound);
    else if (fabs(value-upper)<=primalTolerance*1.001) 
      model_->setStatus(iPivot,ClpSimplex::atUpperBound);
    else
      assert(fabs(value-lower)<=primalTolerance*1.0001||
	     fabs(value-upper)<=primalTolerance*1.0001);
    break;
  }
  if (upper-lower<1.0e-8)
    model_->setFixed(iPivot);
  else
    model_->clearFixed(iPivot);
  double difference = cost-cost_[iRange]; 
  cost = cost_[iRange];
  changeCost_ += value*difference;
  return difference;
}
// Returns nearest bound
double 
ClpNonLinearCost::nearest(int sequence, double solutionValue)
{
  assert (model_!=NULL);
  // get where in bound sequence
  int iRange;
  int start = start_[sequence];
  int end = start_[sequence+1];
  int jRange=-1;
  double nearest=DBL_MAX;
  for (iRange=start; iRange<end;iRange++) {
    if (fabs(solutionValue-lower_[iRange])<nearest) {
      jRange=iRange;
      nearest=fabs(solutionValue-lower_[iRange]);
    }
  }
  assert(jRange<end);
  return lower_[jRange];
}

