// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpFactorization.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpMessage.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream>
/* Dual ranging.
   This computes increase/decrease in cost for each given variable and corresponding
   sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
   and numberColumns.. for artificials/slacks.
   For non-basic variables the sequence number will be that of the non-basic variables.
   
   Up to user to provide correct length arrays.
   
*/
void ClpSimplexOther::dualRanging(int numberCheck,const int * which,
			    double * costIncreased, int * sequenceIncreased,
			    double * costDecreased, int * sequenceDecreased)
{
  rowArray_[0]->clear();
  rowArray_[1]->clear();
  rowArray_[3]->clear();
  columnArray_[0]->clear();
  for ( int i=0;i<numberCheck;i++) {
    int iSequence = which[i];
    double costIncrease=COIN_DBL_MAX;
    double costDecrease=COIN_DBL_MAX;
    int sequenceIncrease=-1;
    int sequenceDecrease=-1;
    
    switch(getStatus(iSequence)) {
      
    case basic:
      {
	// non-trvial
	// Find pivot row (could be faster)
	int iRow=-1;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (iSequence == pivotVariable_[iRow]) 
	    break;
	}
	assert (iRow>=0);
	double plusOne=1.0;
        rowArray_[0]->createPacked(1,&iRow,&plusOne);
	factorization_->updateColumnTranspose(rowArray_[1],rowArray_[0]);
	// put row of tableau in rowArray[0] and columnArray[0]
	matrix_->transposeTimes(this,-1.0,
				rowArray_[0],rowArray_[3],columnArray_[0]);
	// do ratio test up and down
	checkDualRatios(rowArray_[0],columnArray_[0],costIncrease,sequenceIncrease,
		    costDecrease,sequenceDecrease);
      }
      break;
    case isFixed:
      break;
    case isFree:
    case superBasic:
      costIncrease=0.0;
      costDecrease=0.0;
      sequenceIncrease=iSequence;
      sequenceDecrease=iSequence;
      break;
    case atUpperBound:
      costIncrease = max(0.0,-dj_[iSequence]);
      sequenceIncrease = iSequence;
      break;
    case atLowerBound:
      costDecrease = max(0.0,dj_[iSequence]);
      sequenceDecrease = iSequence;
      break;
    }
    double scaleFactor;
    if (rowScale_) {
      if (iSequence<numberColumns_) 
	scaleFactor = 1.0/(objectiveScale_*columnScale_[iSequence]);
      else
	scaleFactor = rowScale_[iSequence-numberColumns_]/objectiveScale_;
    } else {
      scaleFactor = 1.0/objectiveScale_;
    }
    if (costIncrease<1.0e30)
      costIncrease *= scaleFactor;
    if (costDecrease<1.0e30)
      costDecrease *= scaleFactor;
    if (optimizationDirection_==1.0) {
      costIncreased[i] = costIncrease;
      sequenceIncreased[i] = sequenceIncrease;
      costDecreased[i] = costDecrease;
      sequenceDecreased[i] = sequenceDecrease;
    } else if (optimizationDirection_==-1.0) {
      costIncreased[i] = costDecrease;
      sequenceIncreased[i] = sequenceDecrease;
      costDecreased[i] = costIncrease;
      sequenceDecreased[i] = sequenceIncrease;
    } else if (optimizationDirection_==0.0) {
      // !!!!!! ???
      costIncreased[i] = COIN_DBL_MAX;
      sequenceIncreased[i] = -1;
      costDecreased[i] = COIN_DBL_MAX;
      sequenceDecreased[i] = -1;
    } else {
      abort();
    }
  }
  if (!optimizationDirection_)
    printf("*** ????? Ranging with zero optimization costs\n");
}
/* 
   Row array has row part of pivot row
   Column array has column part.
   This is used in dual ranging
*/
void
ClpSimplexOther::checkDualRatios(CoinIndexedVector * rowArray,
			     CoinIndexedVector * columnArray,
			     double & costIncrease, int & sequenceIncrease,
			     double & costDecrease, int & sequenceDecrease)
{
  double acceptablePivot = 1.0e-7;
  double * work;
  int number;
  int * which;
  int iSection;

  double thetaDown = 1.0e31;
  double thetaUp = 1.0e31;
  int sequenceDown =-1;
  int sequenceUp = -1;

  int addSequence;

  for (iSection=0;iSection<2;iSection++) {

    int i;
    if (!iSection) {
      work = rowArray->denseVector();
      number = rowArray->getNumElements();
      which = rowArray->getIndices();
      addSequence = numberColumns_;
    } else {
      work = columnArray->denseVector();
      number = columnArray->getNumElements();
      which = columnArray->getIndices();
      addSequence = 0;
    }
    
    for (i=0;i<number;i++) {
      int iSequence = which[i];
      int iSequence2 = iSequence + addSequence;
      double alpha=work[i];
      if (fabs(alpha)<acceptablePivot)
	continue;
      double oldValue=dj_[iSequence2];

      switch(getStatus(iSequence2)) {
	  
      case basic: 
	break;
      case ClpSimplex::isFixed:
	break;
      case isFree:
      case superBasic:
	// treat dj as if zero
	thetaDown=0.0;
	thetaUp=0.0;
	sequenceDown=iSequence2;
	sequenceUp=iSequence2;
	break;
      case atUpperBound:
	if (alpha>0.0) {
	  // test up
	  if (oldValue + thetaUp*alpha > dualTolerance_) {
	    thetaUp = (dualTolerance_-oldValue)/alpha;
	    sequenceUp = iSequence2;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha > dualTolerance_) {
	    thetaDown = -(dualTolerance_-oldValue)/alpha;
	    sequenceDown = iSequence2;
	  }
	}
	break;
      case atLowerBound:
	if (alpha<0.0) {
	  // test up
	  if (oldValue + thetaUp*alpha <- dualTolerance_) {
	    thetaUp = -(dualTolerance_+oldValue)/alpha;
	    sequenceUp = iSequence2;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha < -dualTolerance_) {
	    thetaDown = (dualTolerance_+oldValue)/alpha;
	    sequenceDown = iSequence2;
	  }
	}
	break;
      }
    }
  }
  if (sequenceUp>=0) {
    costIncrease = thetaUp;
    sequenceIncrease = sequenceUp;
  }
  if (sequenceDown>=0) {
    costDecrease = thetaDown;
    sequenceDecrease = sequenceDown;
  }
}
/** Primal ranging.
    This computes increase/decrease in value for each given variable and corresponding
    sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
    and numberColumns.. for artificials/slacks.
    For basic variables the sequence number will be that of the basic variables.
    
    Up to user to provide correct length arrays.
    
    When here - guaranteed optimal
*/
void 
ClpSimplexOther::primalRanging(int numberCheck,const int * which,
		  double * valueIncreased, int * sequenceIncreased,
		  double * valueDecreased, int * sequenceDecreased)
{
  rowArray_[0]->clear();
  rowArray_[1]->clear();
  lowerIn_=-COIN_DBL_MAX;
  upperIn_=COIN_DBL_MAX;
  valueIn_ = 0.0;
  for ( int i=0;i<numberCheck;i++) {
    int iSequence = which[i];
    double valueIncrease=COIN_DBL_MAX;
    double valueDecrease=COIN_DBL_MAX;
    int sequenceIncrease=-1;
    int sequenceDecrease=-1;
    
    switch(getStatus(iSequence)) {
      
    case basic:
    case isFree:
    case superBasic:
      // Easy
      valueIncrease=max(0.0,upper_[iSequence]-solution_[iSequence]);
      valueDecrease=max(0.0,solution_[iSequence]-lower_[iSequence]);
      sequenceIncrease=iSequence;
      sequenceDecrease=iSequence;
      break;
    case isFixed:
    case atUpperBound:
    case atLowerBound:
      {
	// Non trivial
	// Other bound is ignored
	unpackPacked(rowArray_[1],iSequence);
	factorization_->updateColumn(rowArray_[2],rowArray_[1]);
	// Get extra rows
	matrix_->extendUpdated(this,rowArray_[1],0);
	// do ratio test
	checkPrimalRatios(rowArray_[1],-1);
	if (pivotRow_>=0) {
	  valueIncrease = theta_;
	  sequenceIncrease=pivotVariable_[pivotRow_];
	}
	directionIn_=-1; // down
	checkPrimalRatios(rowArray_[1],1);
	if (pivotRow_>=0) {
	  valueDecrease = theta_;
	  sequenceDecrease=pivotVariable_[pivotRow_];
	}
	rowArray_[1]->clear();
      }
      break;
    }
    double scaleFactor;
    if (rowScale_) {
      if (iSequence<numberColumns_) 
	scaleFactor = columnScale_[iSequence]/rhsScale_;
      else
	scaleFactor = 1.0/(rowScale_[iSequence-numberColumns_]*rhsScale_);
    } else {
      scaleFactor = 1.0/rhsScale_;
    }
    if (valueIncrease<1.0e30)
      valueIncrease *= scaleFactor;
    else
      valueIncrease = COIN_DBL_MAX;
    if (valueDecrease<1.0e30)
      valueDecrease *= scaleFactor;
    else
      valueDecrease = COIN_DBL_MAX;
    valueIncreased[i] = valueIncrease;
    sequenceIncreased[i] = sequenceIncrease;
    valueDecreased[i] = valueDecrease;
    sequenceDecreased[i] = sequenceDecrease;
  }
}
/* 
   Row array has pivot column
   This is used in primal ranging
*/
void 
ClpSimplexOther::checkPrimalRatios(CoinIndexedVector * rowArray,
			 int direction)
{
  // sequence stays as row number until end
  pivotRow_=-1;
  double acceptablePivot=1.0e-7;
  double * work=rowArray->denseVector();
  int number=rowArray->getNumElements();
  int * which=rowArray->getIndices();

  // we need to swap sign if going down
  double way = direction;
  theta_ = 1.0e30;
  for (int iIndex=0;iIndex<number;iIndex++) {
    
    int iRow = which[iIndex];
    double alpha = work[iIndex]*way;
    int iPivot=pivotVariable_[iRow];
    double oldValue = solution_[iPivot];
    if (fabs(alpha)>acceptablePivot) {
      if (alpha>0.0) {
	// basic variable going towards lower bound
	double bound = lower_[iPivot];
	oldValue -= bound;
	if (oldValue-theta_*alpha<0.0) {
	  pivotRow_ = iRow;
	  theta_ = max(0.0,oldValue/alpha);
	}
      } else {
	// basic variable going towards upper bound
	double bound = upper_[iPivot];
	oldValue = oldValue-bound;
	if (oldValue-theta_*alpha>0.0) {
	  pivotRow_ = iRow;
	  theta_ = max(0.0,oldValue/alpha);
	}
      }
    }
  }
}
