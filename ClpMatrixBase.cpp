// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpSimplex.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpMatrixBase::ClpMatrixBase () :
  rhsOffset_(NULL),
  startFraction_(0.0),
  endFraction_(1.0),
  savedBestDj_(0.0),
  originalWanted_(0),
  currentWanted_(0),
  savedBestSequence_(-1),
  type_(-1),
  lastRefresh_(-1),
  refreshFrequency_(0),
  minimumObjectsScan_(-1),
  minimumGoodReducedCosts_(-1),
  trueSequenceIn_(-1),
  trueSequenceOut_(-1),
  skipDualCheck_(false)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpMatrixBase::ClpMatrixBase (const ClpMatrixBase & rhs) :
  type_(rhs.type_),
  skipDualCheck_(rhs.skipDualCheck_)
{  
  startFraction_ = rhs.startFraction_;
  endFraction_ = rhs.endFraction_;
  savedBestDj_ = rhs.savedBestDj_;
  originalWanted_ = rhs.originalWanted_;
  currentWanted_ = rhs.currentWanted_;
  savedBestSequence_ = rhs.savedBestSequence_;
  lastRefresh_ = rhs.lastRefresh_;
  refreshFrequency_ = rhs.refreshFrequency_;
  minimumObjectsScan_ = rhs.minimumObjectsScan_;
  minimumGoodReducedCosts_ = rhs.minimumGoodReducedCosts_;
  trueSequenceIn_ = rhs.trueSequenceIn_;
  trueSequenceOut_ = rhs.trueSequenceOut_;
  skipDualCheck_ = rhs.skipDualCheck_;
  int numberRows = rhs.getNumRows();
  if (rhs.rhsOffset_&&numberRows) {
    rhsOffset_ = ClpCopyOfArray(rhs.rhsOffset_,numberRows);
  } else {
    rhsOffset_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpMatrixBase::~ClpMatrixBase ()
{
  delete [] rhsOffset_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpMatrixBase &
ClpMatrixBase::operator=(const ClpMatrixBase& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    delete [] rhsOffset_;
    int numberRows = rhs.getNumRows();
    if (rhs.rhsOffset_&&numberRows) {
      rhsOffset_ = ClpCopyOfArray(rhs.rhsOffset_,numberRows);
    } else {
      rhsOffset_=NULL;
    }
    startFraction_ = rhs.startFraction_;
    endFraction_ = rhs.endFraction_;
    savedBestDj_ = rhs.savedBestDj_;
    originalWanted_ = rhs.originalWanted_;
    currentWanted_ = rhs.currentWanted_;
    savedBestSequence_ = rhs.savedBestSequence_;
    lastRefresh_ = rhs.lastRefresh_;
    refreshFrequency_ = rhs.refreshFrequency_;
    minimumObjectsScan_ = rhs.minimumObjectsScan_;
    minimumGoodReducedCosts_ = rhs.minimumGoodReducedCosts_;
    trueSequenceIn_ = rhs.trueSequenceIn_;
    trueSequenceOut_ = rhs.trueSequenceOut_;
    skipDualCheck_ = rhs.skipDualCheck_;
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
ClpMatrixBase::partialPricing(ClpSimplex * model, double start, double end,
			      int & bestSequence, int & numberWanted)
{
  std::cerr<<"partialPricing not supported - ClpMatrixBase"<<std::endl;
  abort();
}
/* expands an updated column to allow for extra rows which the main
   solver does not know about and returns number added.  If the arrays are NULL 
   then returns number of extra entries needed.
   
   This will normally be a no-op - it is in for GUB!
*/
int 
ClpMatrixBase::extendUpdated(ClpSimplex * model,CoinIndexedVector * update,int mode)
{
  return 0;
}
/*
     utility primal function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
void 
ClpMatrixBase::primalExpanded(ClpSimplex * model,int mode)
{
}
/*
     utility dual function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
void 
ClpMatrixBase::dualExpanded(ClpSimplex * model,
			    CoinIndexedVector * array,
			    double * other,int mode)
{
}
/*
     general utility function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
int
ClpMatrixBase::generalExpanded(ClpSimplex * model,int mode, int &number)
{
  int returnCode=0;
  switch (mode) {
    // Fill in pivotVariable but not for key variables
  case 0:
    {
      int i;
      int numberBasic=number;
      int numberColumns = model->numberColumns();
      // Use different array so can build from true pivotVariable_
      //int * pivotVariable = model->pivotVariable();
      int * pivotVariable = model->rowArray(0)->getIndices();
      for (i=0;i<numberColumns;i++) {
	if (model->getColumnStatus(i) == ClpSimplex::basic) 
	  pivotVariable[numberBasic++]=i;
      }
      number = numberBasic;
    }
    break;
    // Do initial extra rows + maximum basic
  case 2:
    {
      number = model->numberRows();
    }
    break;
    // To see if can dual or primal
  case 4:
    {
      returnCode= 3;
    }
    break;
  default:
    break;
  }
  return returnCode;
}
// Sets up an effective RHS
void 
ClpMatrixBase::useEffectiveRhs(ClpSimplex * model)
{
  std::cerr<<"useEffectiveRhs not supported - ClpMatrixBase"<<std::endl;
  abort();
}
/* Returns effective RHS if it is being used.  This is used for long problems
   or big gub or anywhere where going through full columns is
   expensive.  This may re-compute */
double * 
ClpMatrixBase::rhsOffset(ClpSimplex * model,bool forceRefresh,bool check)
{
  if (rhsOffset_) {
#ifdef CLP_DEBUG
    if (check) {
      // no need - but check anyway
      // zero out basic
      int numberRows = model->numberRows();
      int numberColumns = model->numberColumns();
      double * solution = new double [numberColumns];
      double * rhs = new double[numberRows];
      const double * solutionSlack = model->solutionRegion(0);
      CoinMemcpyN(model->solutionRegion(),numberColumns,solution);
      int iRow;
      for (iRow=0;iRow<numberRows;iRow++) {
	if (model->getRowStatus(iRow)!=ClpSimplex::basic)
	  rhs[iRow]=solutionSlack[iRow];
	else
	  rhs[iRow]=0.0;
      }
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (model->getColumnStatus(iColumn)==ClpSimplex::basic)
	  solution[iColumn]=0.0;
      }
      times(-1.0,solution,rhs);
      delete [] solution;
      for (iRow=0;iRow<numberRows;iRow++) {
	if (fabs(rhs[iRow]-rhsOffset_[iRow])>1.0e-3)
	  printf("** bad effective %d - true %g old %g\n",iRow,rhs[iRow],rhsOffset_[iRow]);
      }
    }
#endif
    if (forceRefresh||(refreshFrequency_&&model->numberIterations()>=
		       lastRefresh_+refreshFrequency_)) {
      // zero out basic
      int numberRows = model->numberRows();
      int numberColumns = model->numberColumns();
      double * solution = new double [numberColumns];
      const double * solutionSlack = model->solutionRegion(0);
      CoinMemcpyN(model->solutionRegion(),numberColumns,solution);
      for (int iRow=0;iRow<numberRows;iRow++) {
	if (model->getRowStatus(iRow)!=ClpSimplex::basic)
	  rhsOffset_[iRow]=solutionSlack[iRow];
	else
	  rhsOffset_[iRow]=0.0;
      }
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (model->getColumnStatus(iColumn)==ClpSimplex::basic)
	  solution[iColumn]=0.0;
      }
      times(-1.0,solution,rhsOffset_);
      delete [] solution;
      lastRefresh_ = model->numberIterations();
    }
  }
  return rhsOffset_;
}
/* 
   update information for a pivot (and effective rhs)
*/
int 
ClpMatrixBase::updatePivot(ClpSimplex * model,double oldInValue, double oldOutValue)
{
  if (rhsOffset_) {
    // update effective rhs
    int sequenceIn = model->sequenceIn();
    int sequenceOut = model->sequenceOut();
    double * solution = model->solutionRegion();
    int numberColumns = model->numberColumns();
    if (sequenceIn==sequenceOut) {
      if (sequenceIn<numberColumns)
	add(model,rhsOffset_,sequenceIn,oldInValue-solution[sequenceIn]);
    } else {
      if (sequenceIn<numberColumns)
	add(model,rhsOffset_,sequenceIn,oldInValue);
      if (sequenceOut<numberColumns)
	add(model,rhsOffset_,sequenceOut,-solution[sequenceOut]);
    }
  }
  return 0;
}
int 
ClpMatrixBase::hiddenRows() const
{ 
  return 0;
}
/* Creates a variable.  This is called after partial pricing and may modify matrix.
   May update bestSequence.
*/
void 
ClpMatrixBase::createVariable(ClpSimplex * model, int & bestSequence)
{
}
// Returns reduced cost of a variable
double 
ClpMatrixBase::reducedCost(ClpSimplex * model,int sequence) const
{
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  if (sequence<numberRows+numberColumns)
    return model->djRegion()[sequence];
  else
    return savedBestDj_;
}
/* Just for debug if odd type matrix.
   Returns number of primal infeasibilities.
*/
int 
ClpMatrixBase::checkFeasible() const 
{
  return 0;
}
// Correct sequence in and out to give true value
void 
ClpMatrixBase::correctSequence(int & sequenceIn, int & sequenceOut) const
{
}
// Really scale matrix
void 
ClpMatrixBase::reallyScale(const double * rowScale, const double * columnScale)
{
  std::cerr<<"reallyScale not supported - ClpMatrixBase"<<std::endl;
  abort();
}


