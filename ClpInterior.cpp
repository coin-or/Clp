// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.




#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpInterior.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpMessage.hpp"
#include "ClpLinearObjective.hpp"
#include <cfloat>

#include <string>
#include <stdio.h>
#include <iostream>
//#############################################################################

ClpInterior::ClpInterior () :

  ClpModel(),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL)
{
  solveType_=2; // say interior based life form
}

// Subproblem constructor
ClpInterior::ClpInterior ( const ClpModel * rhs,
		     int numberRows, const int * whichRow,
		     int numberColumns, const int * whichColumn,
		     bool dropNames, bool dropIntegers)
  : ClpModel(rhs, numberRows, whichRow,
	     numberColumns,whichColumn,dropNames,dropIntegers),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL)
{
  solveType_=2; // say interior based life form
  
}

//-----------------------------------------------------------------------------

ClpInterior::~ClpInterior ()
{
  gutsOfDelete();
}
//#############################################################################
/* 
   This does housekeeping
*/
int 
ClpInterior::housekeeping()
{
  numberIterations_++;
  return 0;
}
// Copy constructor. 
ClpInterior::ClpInterior(const ClpInterior &rhs) :
  ClpModel(rhs),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL)
{
  gutsOfDelete();
  gutsOfCopy(rhs);
  solveType_=2; // say interior based life form
}
// Copy constructor from model
ClpInterior::ClpInterior(const ClpModel &rhs) :
  ClpModel(rhs),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL)
{
  solveType_=2; // say interior based life form
}
// Assignment operator. This copies the data
ClpInterior & 
ClpInterior::operator=(const ClpInterior & rhs)
{
  if (this != &rhs) {
    gutsOfDelete();
    ClpModel::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
ClpInterior::gutsOfCopy(const ClpInterior & rhs)
{
  lower_ = ClpCopyOfArray(rhs.lower_,numberColumns_+numberRows_);
  rowLowerWork_ = lower_+numberColumns_;
  columnLowerWork_ = lower_;
  upper_ = ClpCopyOfArray(rhs.upper_,numberColumns_+numberRows_);
  rowUpperWork_ = upper_+numberColumns_;
  columnUpperWork_ = upper_;
  //cost_ = ClpCopyOfArray(rhs.cost_,2*(numberColumns_+numberRows_));
  cost_ = ClpCopyOfArray(rhs.cost_,numberColumns_);
  largestPrimalError_ = rhs.largestPrimalError_;
  largestDualError_ = rhs.largestDualError_;
  sumDualInfeasibilities_ = rhs.sumDualInfeasibilities_;
  sumPrimalInfeasibilities_ = rhs.sumPrimalInfeasibilities_;
  solveType_=rhs.solveType_;
}
// type == 0 do everything, most + pivot data, 2 factorization data as well
void 
ClpInterior::gutsOfDelete()
{
  delete [] lower_;
  lower_=NULL;
  rowLowerWork_=NULL;
  columnLowerWork_=NULL;
  delete [] upper_;
  upper_=NULL;
  rowUpperWork_=NULL;
  columnUpperWork_=NULL;
  delete [] cost_;
  cost_=NULL;
}
bool
ClpInterior::createWorkingData()
{
  bool goodMatrix=true;
  //check matrix
  if (!matrix_->allElementsInRange(this,1.0e-12,1.0e20)) {
    problemStatus_=4;
    goodMatrix= false;
  }
  // check rim of problem okay
  if (!sanityCheck())
    goodMatrix=false;
  return goodMatrix;
}
void
ClpInterior::deleteWorkingData()
{
}
// Sanity check on input data - returns true if okay
bool 
ClpInterior::sanityCheck()
{
  // bad if empty
  if (!numberRows_||!numberColumns_||!matrix_->getNumElements()) {
    handler_->message(CLP_EMPTY_PROBLEM,messages_)
      <<numberRows_
      <<numberColumns_
      <<matrix_->getNumElements()
      <<CoinMessageEol;
    problemStatus_=4;
    return false;
  }
  int numberBad ;
  double largestBound, smallestBound, minimumGap;
  double smallestObj, largestObj;
  int firstBad;
  int modifiedBounds=0;
  int i;
  numberBad=0;
  firstBad=-1;
  minimumGap=1.0e100;
  smallestBound=1.0e100;
  largestBound=0.0;
  smallestObj=1.0e100;
  largestObj=0.0;
  // If bounds are too close - fix
  double fixTolerance = 1.1*primalTolerance();
  for (i=numberColumns_;i<numberColumns_+numberRows_;i++) {
    double value;
    value = fabs(cost_[i]);
    if (value>1.0e50) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value) {
      if (value>largestObj)
	largestObj=value;
      if (value<smallestObj)
	smallestObj=value;
    }
    value=upper_[i]-lower_[i];
    if (value<-primalTolerance()) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value<=fixTolerance) {
      if (value) {
	// modify
	upper_[i] = lower_[i];
	modifiedBounds++;
      }
    } else {
      if (value<minimumGap)
	minimumGap=value;
    }
    if (lower_[i]>-1.0e100&&lower_[i]) {
      value = fabs(lower_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
    if (upper_[i]<1.0e100&&upper_[i]) {
      value = fabs(upper_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
  }
  if (largestBound)
    handler_->message(CLP_RIMSTATISTICS3,messages_)
      <<smallestBound
      <<largestBound
      <<minimumGap
      <<CoinMessageEol;
  minimumGap=1.0e100;
  smallestBound=1.0e100;
  largestBound=0.0;
  for (i=0;i<numberColumns_;i++) {
    double value;
    value = fabs(cost_[i]);
    if (value>1.0e50) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value) {
      if (value>largestObj)
	largestObj=value;
      if (value<smallestObj)
	smallestObj=value;
    }
    value=upper_[i]-lower_[i];
    if (value<-primalTolerance()) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value<=fixTolerance) {
      if (value) {
	// modify
	upper_[i] = lower_[i];
	modifiedBounds++;
      }
    } else {
      if (value<minimumGap)
	minimumGap=value;
    }
    if (lower_[i]>-1.0e100&&lower_[i]) {
      value = fabs(lower_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
    if (upper_[i]<1.0e100&&upper_[i]) {
      value = fabs(upper_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
  }
  char rowcol[]={'R','C'};
  if (numberBad) {
    handler_->message(CLP_BAD_BOUNDS,messages_)
      <<numberBad
      <<rowcol[isColumn(firstBad)]<<sequenceWithin(firstBad)
      <<CoinMessageEol;
    problemStatus_=4;
    return false;
  }
  if (modifiedBounds)
    handler_->message(CLP_MODIFIEDBOUNDS,messages_)
      <<modifiedBounds
      <<CoinMessageEol;
  handler_->message(CLP_RIMSTATISTICS1,messages_)
    <<smallestObj
    <<largestObj
    <<CoinMessageEol;  if (largestBound)
    handler_->message(CLP_RIMSTATISTICS2,messages_)
      <<smallestBound
      <<largestBound
      <<minimumGap
      <<CoinMessageEol;
  return true;
}
/* Loads a problem (the constraints on the
   rows are given by lower and upper bounds). If a pointer is 0 then the
   following values are the default:
   <ul>
   <li> <code>colub</code>: all columns have upper bound infinity
   <li> <code>collb</code>: all columns have lower bound 0 
   <li> <code>rowub</code>: all rows have upper bound infinity
   <li> <code>rowlb</code>: all rows have lower bound -infinity
   <li> <code>obj</code>: all variables have 0 objective coefficient
   </ul>
*/
void 
ClpInterior::loadProblem (  const ClpMatrixBase& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
}
void 
ClpInterior::loadProblem (  const CoinPackedMatrix& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
}

/* Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). */
void 
ClpInterior::loadProblem (  const int numcols, const int numrows,
		    const CoinBigIndex* start, const int* index,
		    const double* value,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(numcols, numrows, start, index, value,
			  collb, colub, obj, rowlb, rowub,
			  rowObjective);
}
void 
ClpInterior::loadProblem (  const int numcols, const int numrows,
			   const CoinBigIndex* start, const int* index,
			   const double* value,const int * length,
			   const double* collb, const double* colub,   
			   const double* obj,
			   const double* rowlb, const double* rowub,
			   const double * rowObjective)
{
  ClpModel::loadProblem(numcols, numrows, start, index, value, length,
			  collb, colub, obj, rowlb, rowub,
			  rowObjective);
}
// Read an mps file from the given filename
int 
ClpInterior::readMps(const char *filename,
	    bool keepNames,
	    bool ignoreErrors)
{
  int status = ClpModel::readMps(filename,keepNames,ignoreErrors);
  return status;
}
#include "ClpPdco.hpp"
/* Pdco algorithm - see ClpPdco.hpp for method */
int 
ClpInterior::pdco()
{
  return ((ClpPdco *) this)->pdco();
}
// ** Temporary version
int  
ClpInterior::pdco( Lsqr *lsqr, Options &options, Info &info, Outfo &outfo)
{
  return ((ClpPdco *) this)->pdco(lsqr,options,info,outfo);
}
