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
  effectiveRhs_(NULL),
  type_(-1),
  lastRefresh_(-1),
  refreshFrequency_(0),
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
  lastRefresh_ = rhs.lastRefresh_;
  refreshFrequency_ = rhs.refreshFrequency_;
  skipDualCheck_ = rhs.skipDualCheck_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpMatrixBase::~ClpMatrixBase ()
{
  delete [] effectiveRhs_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpMatrixBase &
ClpMatrixBase::operator=(const ClpMatrixBase& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    delete [] effectiveRhs_;
    int numberRows = getNumRows();
    if (rhs.effectiveRhs_&&numberRows) {
      effectiveRhs_ = ClpCopyOfArray(rhs.effectiveRhs_,numberRows);
    } else {
      effectiveRhs_=NULL;
    }
    lastRefresh_ = rhs.lastRefresh_;
    refreshFrequency_ = rhs.refreshFrequency_;
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
ClpMatrixBase::partialPricing(ClpSimplex * model, int start, int end,
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
  int numberColumns = model->numberColumns();
  switch (mode) {
    // Fill in pivotVariable but not for key variables
  case 0:
    {
      int i;
      int numberBasic=number;
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
    // Make all key variables basic
  case 1:
    {
    }
    break;
    // Do initial extra rows + maximum basic
  case 2:
    {
      returnCode= 0;
      number = model->numberRows();
    }
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
ClpMatrixBase::effectiveRhs(ClpSimplex * model,bool forceRefresh,bool check)
{
  if (effectiveRhs_) {
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
	if (fabs(rhs[iRow]-effectiveRhs_[iRow])>1.0e-3)
	  printf("** bad effective %d - true %g old %g\n",iRow,rhs[iRow],effectiveRhs_[iRow]);
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
	  effectiveRhs_[iRow]=solutionSlack[iRow];
	else
	  effectiveRhs_[iRow]=0.0;
      }
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (model->getColumnStatus(iColumn)==ClpSimplex::basic)
	  solution[iColumn]=0.0;
      }
      times(-1.0,solution,effectiveRhs_);
      delete [] solution;
      lastRefresh_ = model->numberIterations();
    }
  }
  return effectiveRhs_;
}
/* 
   update information for a pivot (and effective rhs)
*/
int 
ClpMatrixBase::updatePivot(ClpSimplex * model,double oldInValue, double oldOutValue)
{
  if (effectiveRhs_) {
    // update effective rhs
    int sequenceIn = model->sequenceIn();
    int sequenceOut = model->sequenceOut();
    double * solution = model->solutionRegion();
    int numberColumns = model->numberColumns();
    if (sequenceIn==sequenceOut) {
      if (sequenceIn<numberColumns)
	add(model,effectiveRhs_,sequenceIn,oldInValue-solution[sequenceIn]);
      else
	effectiveRhs_[sequenceIn-numberColumns] -= oldInValue-solution[sequenceIn];
    } else {
      if (sequenceIn<numberColumns)
	add(model,effectiveRhs_,sequenceIn,oldInValue);
      else
	effectiveRhs_[sequenceIn-numberColumns] -= oldInValue;
      if (sequenceOut<numberColumns)
	add(model,effectiveRhs_,sequenceOut,-solution[sequenceOut]);
      else
	effectiveRhs_[sequenceOut-numberColumns] -= -solution[sequenceOut];
    }
  }
  return 0;
}
int 
ClpMatrixBase::hiddenRows() const
{ 
  return 0;
}

