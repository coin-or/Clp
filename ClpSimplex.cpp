// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.



#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpMessage.hpp"
#include <cfloat>

#include <string>
#include <stdio.h>
#include <iostream>
// This returns a non const array filled with input from scalar
// or actual array
template <class T> inline T*
copyOfArray( const T * array, const int size, T value)
{
  T * arrayNew = new T[size];
  if (array) 
    CoinDisjointCopyN(array,size,arrayNew);
  else
    CoinFillN ( arrayNew, size,value);
  return arrayNew;
}

// This returns a non const array filled with actual array (or NULL)
template <class T> inline T*
copyOfArray( const T * array, const int size)
{
  if (array) {
    T * arrayNew = new T[size];
    CoinDisjointCopyN(array,size,arrayNew);
    return arrayNew;
  } else {
    return NULL;
  }
}

//#############################################################################

ClpSimplex::ClpSimplex () :

  ClpModel(),
  columnPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalInfeasibility_(0.0),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualInfeasibility_(0.0),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e7),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(0.0),
  dualIn_(0.0),
  sequenceIn_(-1),
  directionIn_(-1),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  status_(NULL),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  sumPrimalInfeasibilities_(0.0),
  numberPrimalInfeasibilities_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  numberRefinements_(0),
  rowScale_(NULL),
  savedSolution_(NULL),
  columnScale_(NULL),
  scalingFlag_(0),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  infeasibilityCost_(1.0e7),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  numberFake_(0)

{
  int i;
  for (i=0;i<6;i++) {
    rowArray_[i]=NULL;
    columnArray_[i]=NULL;
  }
  saveStatus_=NULL;
  // get an empty factorization so we can set tolerances etc
  factorization_ = new ClpFactorization();
  // Say sparse
  factorization_->sparseThreshold(1);
  // say Steepest pricing
  dualRowPivot_ = new ClpDualRowSteepest();
  // say Steepest pricing
  primalColumnPivot_ = new ClpPrimalColumnSteepest();
  
}


//-----------------------------------------------------------------------------

ClpSimplex::~ClpSimplex ()
{
  gutsOfDelete(0);
}
//#############################################################################
void ClpSimplex::setLargeValue( double value) 
{
  if (value>0.0&&value<DBL_MAX)
    largeValue_=value;
}
int
ClpSimplex::gutsOfSolution ( const double * rowActivities,
			     const double * columnActivities,
			     bool valuesPass)
{


  // if values pass, save values of basic variables
  double * save = NULL;
  double oldValue=0.0;
  if (valuesPass) {
    assert(algorithm_>0); // only primal at present
    assert(nonLinearCost_);
    int iRow;
    checkPrimalSolution( rowActivities, columnActivities);
    // get correct bounds on all variables
    nonLinearCost_->checkInfeasibilities();
    oldValue = nonLinearCost_->largestInfeasibility();
    save = new double[numberRows_];
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      save[iRow] = solution_[iPivot];
    }
  }
  // do work
  computePrimals(rowActivities, columnActivities);
  double objectiveModification = 0.0;
  if (algorithm_>0&&nonLinearCost_!=NULL) {
    // primal algorithm
    // get correct bounds on all variables
    nonLinearCost_->checkInfeasibilities();
    objectiveModification += nonLinearCost_->changeInCost();
    if (nonLinearCost_->numberInfeasibilities())
      handler_->message(CLP_SIMPLEX_NONLINEAR,messages_)
	<<nonLinearCost_->changeInCost()
	<<nonLinearCost_->numberInfeasibilities()
	<<CoinMessageEol;
  }
  if (valuesPass) {
#ifdef CLP_DEBUG
    std::cout<<"Largest given infeasibility "<<oldValue
	     <<" now "<<nonLinearCost_->largestInfeasibility()<<std::endl;
#endif
    int numberOut=0;
    if (oldValue<1.0
	&&nonLinearCost_->largestInfeasibility()>1.1*oldValue+1.0e-4||
	largestPrimalError_>1.0e-3) {
      // throw out up to 1000 structurals
      int iRow;
      int * sort = new int[numberRows_];
      // first put back solution and store difference
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];
	double difference = fabs(solution_[iPivot]=save[iRow]);
	solution_[iPivot]=save[iRow];
	save[iRow]=difference;
      }
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];

	if (iPivot<numberColumns_) {
	  // column
	  double difference= save[iRow];
	  if (difference>1.0e-4) {
	    sort[numberOut]=iPivot;
	    save[numberOut++]=difference;
	  }
	}
      }
      CoinSort_2(save, save + numberOut, sort,
		 CoinFirstGreater_2<double, int>());
      numberOut = min(1000,numberOut);
      for (iRow=0;iRow<numberOut;iRow++) {
	int iColumn=sort[iRow];
	setColumnStatus(iColumn,ClpSimplex::superBasic);

      }
      delete [] sort;
    }
    delete [] save;
    if (numberOut)
      return numberOut;
  }

  computeDuals();

  // now check solutions
  checkPrimalSolution( rowActivities, columnActivities);
  objectiveValue_ += objectiveModification;
  checkDualSolution();
  if (handler_->logLevel()>3||(largestPrimalError_>1.0e-2||
			       largestDualError_>1.0e-2)) 
    handler_->message(CLP_SIMPLEX_ACCURACY,messages_)
      <<largestPrimalError_
      <<largestDualError_
      <<CoinMessageEol;
  return 0;
}
void
ClpSimplex::computePrimals ( const double * rowActivities,
				     const double * columnActivities)
{

  //work space
  CoinIndexedVector  * workSpace = rowArray_[0];

  double * array = new double [numberRows_];
  double * save = new double [numberRows_];
  double * previous = new double [numberRows_];

  // accumulate non basic stuff 
  CoinFillN(array,numberRows_,0.0);

  int iRow;
  // order is this way for scaling
  // Use whole matrix every time to make it easier for ClpMatrixBase
  // So zero out basic
  if (columnActivities!=columnActivityWork_)
    CoinDisjointCopyN(columnActivities,numberColumns_,columnActivityWork_);
  if (rowActivities!=rowActivityWork_)
    CoinDisjointCopyN(rowActivities,numberRows_,rowActivityWork_);
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot=pivotVariable_[iRow];
    solution_[iPivot] = 0.0;
  }
  times(-1.0,columnActivityWork_,array);

  for (iRow=0;iRow<numberRows_;iRow++) {
    array[iRow] += rowActivityWork_[iRow];
  }

  // Ftran adjusted RHS and iterate to improve accuracy
  double lastError=DBL_MAX;
  int iRefine;
  double * work = workSpace->denseVector();
  factorization_->updateColumn(workSpace,array);

  for (iRefine=0;iRefine<numberRefinements_+1;iRefine++) {
    // save in case we want to restore
    CoinDisjointCopyN ( array, numberRows_ , save);

    // put solution in correct place
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      solution_[iPivot] = array[iRow];
    }
    // check Ax == b  (for all)
    times(-1.0,columnActivityWork_,work);
    for (iRow=0;iRow<numberRows_;iRow++) {
      work[iRow] += rowActivityWork_[iRow];
    }

    largestPrimalError_=0.0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (fabs(work[iRow])>largestPrimalError_) {
	largestPrimalError_=fabs(work[iRow]);
      }
      //work[iRow] -= save[iRow];
    }
    if (largestPrimalError_>=lastError) {
      // restore
      double * temp = array;
      array = previous;
      previous=temp;
      break;
    }
    if (iRefine<numberRefinements_&&largestPrimalError_>1.0e-12) {
      // try and make better
      // save this
      double * temp = array;
      array = previous;
      previous=temp;
      double multiplier = 131072.0;
      for (iRow=0;iRow<numberRows_;iRow++) {
	array[iRow] = multiplier*work[iRow];
	work[iRow]=0.0;
      }
      lastError=largestPrimalError_;
      factorization_->updateColumn(workSpace,array);
      multiplier = 1.0/multiplier;
      for (iRow=0;iRow<numberRows_;iRow++) {
	array[iRow] = previous[iRow] + multiplier*array[iRow];
	work[iRow]=0.0;
      }
    }
  }

  // solution as accurate as we are going to get
  CoinFillN(work,numberRows_,0.0);
  // put solution in correct place
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot=pivotVariable_[iRow];
    solution_[iPivot] = array[iRow];
  }
  delete [] previous;
  delete [] array;
  delete [] save;
}
// now dual side
void
ClpSimplex::computeDuals()
{
  double slackValue = factorization_->slackValue();
  //work space
  CoinIndexedVector  * workSpace = rowArray_[0];

  double * array = new double [numberRows_];
  double * save = new double [numberRows_];
  double * previous = new double [numberRows_];

  int iRow;
#ifdef CLP_DEBUG
  workSpace->checkClear();
#endif
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot=pivotVariable_[iRow];
    if (iPivot>=numberColumns_) {
      // slack
      array[iRow] = rowObjectiveWork_[iPivot-numberColumns_];
    } else {
      // column
      array[iRow]=objectiveWork_[iPivot];
    }
  }
  CoinDisjointCopyN ( array, numberRows_ , save);

  // Btran basic costs and get as accurate as possible
  double lastError=DBL_MAX;
  int iRefine;
  double * work = workSpace->denseVector();
  factorization_->updateColumnTranspose(workSpace,array);
  for (iRefine=0;iRefine<numberRefinements_+1;iRefine++) {
    // check basic reduced costs zero
    largestDualError_=0.0;
    // would be faster to do just for basic but this reduces code
    CoinDisjointCopyN(objectiveWork_,numberColumns_,reducedCostWork_);
    transposeTimes(-1.0,array,reducedCostWork_);
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      if (iPivot>=numberColumns_) {
	// slack
	//work[iRow] += slackValue*array[iPivot-numberColumns_];
	//work[iRow] += rowObjectiveWork_[iPivot-numberColumns_];
	work[iRow] = rowObjectiveWork_[iPivot-numberColumns_]
	- slackValue*array[iPivot-numberColumns_];
      } else {
	// column
	work[iRow] = reducedCostWork_[iPivot];
      }
      if (fabs(work[iRow])>largestDualError_) {
	largestDualError_=fabs(work[iRow]);
      }
    }
    if (largestDualError_>=lastError) {
      // restore
      double * temp = array;
      array = previous;
      previous=temp;
      break;
    }
    if (iRefine<numberRefinements_&&largestDualError_>1.0e-20) {
      // try and make better
      // save this
      double * temp = array;
      array = previous;
      previous=temp;
      double multiplier = 131072.0;
      for (iRow=0;iRow<numberRows_;iRow++) {
	array[iRow] = multiplier*work[iRow];
	work[iRow]=0.0;
      }
      lastError=largestDualError_;
      factorization_->updateColumnTranspose(workSpace,array);
      multiplier = 1.0/multiplier;
      for (iRow=0;iRow<numberRows_;iRow++) {
	array[iRow] = previous[iRow] + multiplier*array[iRow];
	work[iRow]=0.0;
      }
    }
  }
  CoinFillN(work,numberRows_,0.0);
  // now look at dual solution
  for (iRow=0;iRow<numberRows_;iRow++) {
    // slack
    double value = array[iRow];
    dual_[iRow]=value;
    value = rowObjectiveWork_[iRow] - value*slackValue;
    rowReducedCost_[iRow]=value;
  }
  CoinDisjointCopyN(objectiveWork_,numberColumns_,reducedCostWork_);
  transposeTimes(-1.0,dual_,reducedCostWork_);
  delete [] previous;
  delete [] array;
  delete [] save;
}
/* Given an existing factorization computes and checks 
   primal and dual solutions.  Uses input arrays for variables at
   bounds.  Returns feasibility states */
int ClpSimplex::getSolution ( const double * rowActivities,
			       const double * columnActivities)
{
  if (!factorization_->status()) {
    // put in standard form
    createRim(7+8+16);
    // do work
    gutsOfSolution ( rowActivities, columnActivities);
    // release extra memory
    deleteRim();
  }
  return factorization_->status();
}
/* Given an existing factorization computes and checks 
   primal and dual solutions.  Uses current problem arrays for
   bounds.  Returns feasibility states */
int ClpSimplex::getSolution ( )
{
  double * rowActivities = new double[numberRows_];
  double * columnActivities = new double[numberColumns_];
  CoinDisjointCopyN ( rowActivityWork_, numberRows_ , rowActivities);
  CoinDisjointCopyN ( columnActivityWork_, numberColumns_ , columnActivities);
  int status = getSolution( rowActivities, columnActivities);
  delete [] rowActivities;
  delete [] columnActivities;
  return status;
}
// Factorizes using current basis.  This is for external use
// Return codes are as from ClpFactorization
int ClpSimplex::factorize ()
{
  // put in standard form
  createRim(7+8+16);
  // do work
  int status = internalFactorize(-1);
  // release extra memory
  deleteRim();

  return status;
}

/* Factorizes using current basis.  
      solveType - 1 iterating, 0 initial, -1 external 
      If 10 added then in primal values pass
*/
/* Return codes are as from ClpFactorization unless initial factorization
   when total number of singularities is returned */
int ClpSimplex::internalFactorize ( int solveType)
{

  int iRow,iColumn;
  int totalSlacks=numberRows_;

  bool valuesPass=false;
  if (solveType>=10) {
    valuesPass=true;
    solveType -= 10;
  }
  if (!solveType) {
    if (!valuesPass) {
      // not values pass so set to bounds
      bool allSlack=true;
      if (status_) {
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (getRowStatus(iRow)!=ClpSimplex::basic) {
	    allSlack=false;
	    break;
	  }
	}
      }
      if (!allSlack) {
	// set values from warm start (if sensible)
	int numberBasic=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  switch(getRowStatus(iRow)) {
	    
	  case ClpSimplex::basic:
	    numberBasic++;
	    break;
	  case ClpSimplex::isFree:
	    assert(rowLowerWork_[iRow]<-largeValue_);
	    assert(rowUpperWork_[iRow]>largeValue_);
	    rowActivityWork_[iRow]=0.0;
	    break;
	  case ClpSimplex::atUpperBound:
	    rowActivityWork_[iRow]=rowUpperWork_[iRow];
	    if (rowActivityWork_[iRow]>largeValue_) {
	      assert(rowLowerWork_[iRow]>-largeValue_);
	      rowActivityWork_[iRow]=rowLowerWork_[iRow];
	      setRowStatus(iRow,ClpSimplex::atLowerBound);
	    }
	    break;
	  case ClpSimplex::atLowerBound:
	    rowActivityWork_[iRow]=rowLowerWork_[iRow];
	    if (rowActivityWork_[iRow]<-largeValue_) {
	      assert(rowUpperWork_[iRow]<largeValue_);
	      rowActivityWork_[iRow]=rowUpperWork_[iRow];
	      setRowStatus(iRow,ClpSimplex::atUpperBound);
	    }
	    break;
	  case ClpSimplex::superBasic:
	    abort();
	    break;
	  }
	}
	totalSlacks=numberBasic;
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  switch(getColumnStatus(iColumn)) {
	    
	  case ClpSimplex::basic:
	    if (numberBasic==numberRows_) {
	      // take out of basis
	      if (columnLowerWork_[iColumn]>-largeValue_) {
		if (columnActivityWork_[iColumn]-columnLowerWork_[iColumn]<
		    columnUpperWork_[iColumn]-columnActivityWork_[iColumn]) {
		  columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
		  setColumnStatus(iColumn,ClpSimplex::atLowerBound);
		} else {
		  columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		  setColumnStatus(iColumn,ClpSimplex::atUpperBound);
		}
	      } else if (columnUpperWork_[iColumn]<largeValue_) {
		columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		setColumnStatus(iColumn,ClpSimplex::atUpperBound);
	      } else {
		columnActivityWork_[iColumn]=0.0;
		setColumnStatus(iColumn,ClpSimplex::isFree);
	      }
	    } else {
	      numberBasic++;
	    }
	    break;
	  case ClpSimplex::isFree:
	    columnActivityWork_[iColumn]=0.0;
	    break;
	  case ClpSimplex::atUpperBound:
	    columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	    if (columnActivityWork_[iColumn]>largeValue_) {
	      assert(columnLowerWork_[iColumn]>-largeValue_);
	      columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	      setColumnStatus(iColumn,ClpSimplex::atLowerBound);
	    }
	    break;
	  case ClpSimplex::atLowerBound:
	    columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	    if (columnActivityWork_[iColumn]<-largeValue_) {
	      assert(columnUpperWork_[iColumn]<largeValue_);
	      columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	      setColumnStatus(iColumn,ClpSimplex::atUpperBound);
	    }
	    break;
	  case ClpSimplex::superBasic:
	    abort();
	    break;
	  }
	}
      } else {
	//#define TESTFREE
	// all slack basis
	int numberBasic=0;
	// changed to put free variables in basis
	if (!status_) {
	  status_ = new unsigned char [numberColumns_+numberRows_];
	  memset(status_,0,(numberColumns_+numberRows_)*sizeof(char));
	}
	for (iRow=0;iRow<numberRows_;iRow++) {
	  double lower=rowLowerWork_[iRow];
	  double upper=rowUpperWork_[iRow];
	  if (lower>-largeValue_||upper<largeValue_) {
	    if (fabs(lower)<=fabs(upper)) {
	      setRowStatus(iRow,ClpSimplex::atLowerBound);
	      rowActivityWork_[iRow]=lower;
	    } else {
	      setRowStatus(iRow,ClpSimplex::atUpperBound);
	      rowActivityWork_[iRow]=upper;
	    }
	  } else {
	    setRowStatus(iRow,ClpSimplex::isFree);
	    rowActivityWork_[iRow]=0.0;
	  }
#ifdef TESTFREE
	  setRowStatus(iRow,ClpSimplex::basic);
	  numberBasic++;
#endif
	}
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  double lower=columnLowerWork_[iColumn];
	  double upper=columnUpperWork_[iColumn];
	  double big_bound = largeValue_;
	  if (lower>-big_bound||upper<big_bound) {
	    if (fabs(lower)<=fabs(upper)) {
	      setColumnStatus(iColumn,ClpSimplex::atLowerBound);
	      columnActivityWork_[iColumn]=lower;
	    } else {
	      setColumnStatus(iColumn,ClpSimplex::atUpperBound);
	      columnActivityWork_[iColumn]=upper;
	    }
	  } else {
#ifndef TESTFREE
	    numberBasic++;
	    setColumnStatus(iColumn,ClpSimplex::basic);
#else
	    setColumnStatus(iColumn,ClpSimplex::isFree);
#endif
	    columnActivityWork_[iColumn]=0.0;
	  }
	}
	assert(numberBasic<=numberRows_); // problems if too many free
	if (!numberBasic) {
	  // might as well do all slack basis
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    setRowStatus(iRow,ClpSimplex::basic);
	  }
	}
      }
    } else {
      // values pass has less coding
      // make row activities correct
      memset(rowActivityWork_,0,numberRows_*sizeof(double));
      times(1.0,columnActivityWork_,rowActivityWork_);
      if (status_) {
	// set values from warm start (if sensible)
	int numberBasic=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (getRowStatus(iRow)==ClpSimplex::basic) 
	    numberBasic++;
	  else
	    setRowStatus(iRow,ClpSimplex::superBasic);
	}
	totalSlacks=numberBasic;
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  if (getColumnStatus(iColumn)==ClpSimplex::basic) {
	    if (numberBasic==numberRows_) 
	      // take out of basis
	      setColumnStatus(iColumn,ClpSimplex::superBasic);
	    else 
	      numberBasic++;
	  } else {
	    setColumnStatus(iColumn,ClpSimplex::superBasic);
	  }
	}
      } else {
	// all slack basis
	int numberBasic=0;
	if (!status_) {
	  status_ = new unsigned char [numberColumns_+numberRows_];
	  memset(status_,0,(numberColumns_+numberRows_)*sizeof(char));
	}
	for (iRow=0;iRow<numberRows_;iRow++) {
	  setRowStatus(iRow,ClpSimplex::basic);
	  numberBasic++;
	}
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  setColumnStatus(iColumn,ClpSimplex::superBasic);
	}
      }
    }
    numberRefinements_=1;
  }
  int status=-99;
  int * rowIsBasic = new int[numberRows_];
  int * columnIsBasic = new int[numberColumns_];
  //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  while (status<-98) {
  
    int i;
    int numberBasic=0;
    for (i=0;i<numberRows_;i++) {
      if (getRowStatus(i) == ClpSimplex::basic) {
	rowIsBasic[i]=1;
	numberBasic++;
      } else {
	rowIsBasic[i]=-1;
      }
    }
    for (i=0;i<numberColumns_;i++) {
      if (getColumnStatus(i) == ClpSimplex::basic) {
	columnIsBasic[i]=1;
	numberBasic++;
      } else {
	columnIsBasic[i]=-1;
      }
    }
    assert (numberBasic<=numberRows_);
    while (status==-99) {
      status =  factorization_->factorize(this,matrix_,
					  numberRows_,numberColumns_,
					  rowIsBasic, columnIsBasic,
					  0.0);
      if (status==-99) {
	// get more memory
	factorization_->areaFactor(2.0*factorization_->areaFactor());
      }
    }
    if (!status) {
      // do pivot information
      for (i=0;i<numberRows_;i++) {
	if (getRowStatus(i) == ClpSimplex::basic) {
	  pivotVariable_[rowIsBasic[i]]=i+numberColumns_;
	}
      }
      for (i=0;i<numberColumns_;i++) {
	if (getColumnStatus(i) == ClpSimplex::basic) {
	  pivotVariable_[columnIsBasic[i]]=i;
	}
      }
    } else {
      // leave pivotVariable_ in useful form for cleaning basis
      for (i=0;i<numberRows_;i++) {
	pivotVariable_[i]=-1;
      }
      for (i=0;i<numberRows_;i++) {
	if (getRowStatus(i) == ClpSimplex::basic) {
	  int iPivot = rowIsBasic[i];
	  if (iPivot>=0) 
	    pivotVariable_[iPivot]=i+numberColumns_;
	}
      }
      for (i=0;i<numberColumns_;i++) {
	if (getColumnStatus(i) == ClpSimplex::basic) {
	  int iPivot = columnIsBasic[i];
	  if (iPivot>=0) 
	    pivotVariable_[iPivot]=i;
	}
      }
    }
    if (status==-1) {
      if (!solveType) {
	//redo basis - first take ALL columns out
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  if (getColumnStatus(iColumn)==
	      ClpSimplex::basic) {
	    // take out
	    if (!valuesPass) {
	      double lower=columnLowerWork_[iColumn];
	      double upper=columnUpperWork_[iColumn];
	      double value=columnActivityWork_[iColumn];
	      if (lower>-largeValue_||upper<largeValue_) {
		if (fabs(value-lower)<fabs(value-upper)) {
		  setColumnStatus(iColumn,ClpSimplex::atLowerBound);
		  columnActivityWork_[iColumn]=lower;
		} else {
		  setColumnStatus(iColumn,ClpSimplex::atUpperBound);
		  columnActivityWork_[iColumn]=upper;
		}
	      } else {
		setColumnStatus(iColumn,ClpSimplex::isFree);
	      }
	    } else {
	      setColumnStatus(iColumn,ClpSimplex::superBasic);
	    }
	  }
	}
	for (iRow=0;iRow<numberRows_;iRow++) {
	  int iSequence=pivotVariable_[iRow];
	  if (iSequence>=0) {
	    // basic
	    if (iSequence>=numberColumns_) {
	      // slack in - leave
	      assert (iSequence-numberColumns_==iRow);
	    } else {
	      // put back structural
	      setColumnStatus(iSequence,ClpSimplex::basic);
	    }
	  } else {
	    // put in slack
	    setRowStatus(iRow,ClpSimplex::basic);
	  }
	}
	// signal repeat
	status=-99;
      }
    }
  } 
  delete [] rowIsBasic;
  delete [] columnIsBasic;
  if (status) {
    handler_->message(CLP_SIMPLEX_BADFACTOR,messages_)
      <<status
      <<CoinMessageEol;
    return -1;
  } else if (!solveType) {
    // Initial basis - return number of singularities
    int numberSlacks=0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (getRowStatus(iRow) == ClpSimplex::basic)
	numberSlacks++;
    }
    status= max(numberSlacks-totalSlacks,0);
  }

  // sparse methods
  if (factorization_->sparseThreshold()) {
    // get default value
    factorization_->sparseThreshold(0);
    factorization_->goSparse();
  }
  
  return status;
}
/* 
   This does basis housekeeping and does values for in/out variables.
   Can also decide to re-factorize
*/
int 
ClpSimplex::housekeeping(double objectiveChange)
{
  numberIterations_++;
  changeMade_++; // something has happened
  // incoming variable

  handler_->message(CLP_SIMPLEX_HOUSE1,messages_)
    <<directionOut_
    <<directionIn_<<theta_
    <<dualOut_<<dualIn_<<alpha_
    <<CoinMessageEol;
  if (getStatus(sequenceIn_)==ClpSimplex::isFree) {
    handler_->message(CLP_SIMPLEX_FREEIN,messages_)
      <<sequenceIn_
      <<CoinMessageEol;
  }
  // change of incoming
  if (algorithm_<0) {
    dualOut_ /= alpha_;
    dualOut_ *= -directionOut_;
  }
  char rowcol[]={'R','C'};
  double cost = cost_[sequenceIn_];
  double value=valueIn_;
  if (pivotRow_>=0)
    pivotVariable_[pivotRow_]=sequenceIn();
  if (directionIn_==-1) {
    // as if from upper bound
    if (sequenceIn_!=sequenceOut_) {
      // variable becoming basic
      setStatus(sequenceIn_,ClpSimplex::basic);
      if (algorithm_<0) {
	value = upperIn_+dualOut_;
	dj_[sequenceIn_]=0.0;
      } else {
	value = valueIn_-fabs(theta_);
      }
    } else {
      value=lowerIn_;
      setStatus(sequenceIn_, ClpSimplex::atLowerBound);
    }
  } else {
    // as if from lower bound
    if (sequenceIn_!=sequenceOut_) {
      // variable becoming basic
      setStatus(sequenceIn_,ClpSimplex::basic);
      if (algorithm_<0) {
	value = lowerIn_+dualOut_;
	dj_[sequenceIn_]=0.0;
      } else {
	value = valueIn_+fabs(theta_);
      }
    } else {
      value=upperIn_;
      setStatus(sequenceIn_, ClpSimplex::atUpperBound);
    }
  }
  if (algorithm_<0)
    objectiveChange += cost*(value-valueIn_);
  else
    objectiveChange += dualIn_*(value-valueIn_);
  solution_[sequenceIn_]=value;

  // outgoing
  if (sequenceIn_!=sequenceOut_) {
    assert( getStatus(sequenceOut_)== ClpSimplex::basic);
    if (algorithm_<0) {
      if (directionOut_>0) {
	value = lowerOut_;
	setStatus(sequenceOut_,ClpSimplex::atLowerBound);
	dj_[sequenceOut_]=theta_;
      } else {
	value = upperOut_;
	setStatus(sequenceOut_,ClpSimplex::atUpperBound);
	dj_[sequenceOut_]=-theta_;
      }
      solution_[sequenceOut_]=value;
    } else {
      if (directionOut_>0) {
	value = lowerOut_;
      } else {
	value = upperOut_;
      }
      nonLinearCost_->setOne(sequenceOut_,value);
      double lowerValue = lower_[sequenceOut_];
      double upperValue = upper_[sequenceOut_];
      assert(value>=lowerValue-primalTolerance_&&
	     value<=upperValue+primalTolerance_);
      // may not be exactly at bound and bounds may have changed
      if (value<=lowerValue+primalTolerance_) {
	setStatus(sequenceOut_,ClpSimplex::atLowerBound);
      } else if (value>=upperValue-primalTolerance_) {
	setStatus(sequenceOut_,ClpSimplex::atUpperBound);
      } else {
	printf("*** variable wandered off bound %g %g %g!\n",
	       lowerValue,value,upperValue);
	if (value-lowerValue<=upperValue-value) {
	  setStatus(sequenceOut_,ClpSimplex::atLowerBound);
	  value=lowerValue;
	} else {
	  setStatus(sequenceOut_,ClpSimplex::atUpperBound);
	  value=upperValue;
	}
      }
      solution_[sequenceOut_]=value;
    }
  }
  // change cost and bounds on incoming if primal
  if (algorithm_>0)
    nonLinearCost_->setOne(sequenceIn_,solution_[sequenceIn_]); 
  objectiveValue_ += objectiveChange;
  handler_->message(CLP_SIMPLEX_HOUSE2,messages_)
    <<numberIterations_<<objectiveValue_
    <<rowcol[isColumn(sequenceIn_)]<<sequenceWithin(sequenceIn_)
    <<rowcol[isColumn(sequenceOut_)]<<sequenceWithin(sequenceOut_);
  handler_->printing(algorithm_<0)<<theta_<<dualOut_;
  handler_->printing(algorithm_>0)<<dualIn_<<theta_;
  handler_->message()<<CoinMessageEol;

  if (numberIterations_>=maximumIterations_)
    return 2;
  // only time to re-factorize if one before real time
  // this is so user won't be surprised that maximumPivots has exact meaning
  if (factorization_->pivots()==factorization_->maximumPivots()) {
    return 1;
  } else {
    if (forceFactorization_>0&&
	factorization_->pivots()==forceFactorization_) {
      // relax
      forceFactorization_ *= 2;
      if (forceFactorization_>factorization_->maximumPivots())
	forceFactorization_ = -1; //off
      return 1;
    } else {
      // carry on iterating
      return 0;
    }
  }
}
// Copy constructor. 
ClpSimplex::ClpSimplex(const ClpSimplex &rhs) :
  ClpModel(rhs),
  columnPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalInfeasibility_(0.0),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualInfeasibility_(0.0),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e7),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(0.0),
  dualIn_(0.0),
  sequenceIn_(-1),
  directionIn_(-1),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  status_(NULL),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  sumPrimalInfeasibilities_(0.0),
  numberPrimalInfeasibilities_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  numberRefinements_(0),
  rowScale_(NULL),
  savedSolution_(NULL),
  columnScale_(NULL),
  scalingFlag_(0),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  infeasibilityCost_(1.0e7),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  numberFake_(0)
{
  int i;
  for (i=0;i<6;i++) {
    rowArray_[i]=NULL;
    columnArray_[i]=NULL;
  }
  saveStatus_=NULL;
  factorization_ = NULL;
  dualRowPivot_ = NULL;
  primalColumnPivot_ = NULL;
  gutsOfDelete(0);
  gutsOfCopy(rhs);
}
// Copy constructor from model
ClpSimplex::ClpSimplex(const ClpModel &rhs) :
  ClpModel(rhs),
  columnPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalInfeasibility_(0.0),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualInfeasibility_(0.0),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e7),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(0.0),
  dualIn_(0.0),
  sequenceIn_(-1),
  directionIn_(-1),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  status_(NULL),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  sumPrimalInfeasibilities_(0.0),
  numberPrimalInfeasibilities_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  numberRefinements_(0),
  rowScale_(NULL),
  savedSolution_(NULL),
  columnScale_(NULL),
  scalingFlag_(0),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  infeasibilityCost_(1.0e7),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  numberFake_(0)
{
  int i;
  for (i=0;i<6;i++) {
    rowArray_[i]=NULL;
    columnArray_[i]=NULL;
  }
  saveStatus_=NULL;
  // get an empty factorization so we can set tolerances etc
  factorization_ = new ClpFactorization();
  // say Steepest pricing
  dualRowPivot_ = new ClpDualRowSteepest();
  // say Dantzig pricing
  primalColumnPivot_ = new ClpPrimalColumnDantzig();
  
}
// Assignment operator. This copies the data
ClpSimplex & 
ClpSimplex::operator=(const ClpSimplex & rhs)
{
  if (this != &rhs) {
    gutsOfDelete(0);
    ClpModel::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
ClpSimplex::gutsOfCopy(const ClpSimplex & rhs)
{
  lower_ = copyOfArray(rhs.lower_,numberColumns_+numberRows_);
  rowLowerWork_ = lower_+numberColumns_;
  columnLowerWork_ = lower_;
  upper_ = copyOfArray(rhs.upper_,numberColumns_+numberRows_);
  rowUpperWork_ = upper_+numberColumns_;
  columnUpperWork_ = upper_;
  cost_ = copyOfArray(rhs.cost_,(numberColumns_+numberRows_));
  objectiveWork_ = cost_;
  rowObjectiveWork_ = cost_+numberColumns_;
  dj_ = copyOfArray(rhs.dj_,numberRows_+numberColumns_);
  if (dj_) {
    reducedCostWork_ = dj_;
    rowReducedCost_ = dj_+numberColumns_;
  }
  solution_ = copyOfArray(rhs.solution_,numberRows_+numberColumns_);
  if (solution_) {
    columnActivityWork_ = solution_;
    rowActivityWork_ = solution_+numberColumns_;
  }
  if (rhs.pivotVariable_) {
    pivotVariable_ = new int[numberRows_];
    CoinDisjointCopyN ( rhs.pivotVariable_, numberRows_ , pivotVariable_);
  } else {
    pivotVariable_=NULL;
  }
  if (rhs.factorization_) {
    factorization_ = new ClpFactorization(*rhs.factorization_);
  } else {
    factorization_=NULL;
  }
  rowScale_ = copyOfArray(rhs.rowScale_,numberRows_);
  savedSolution_ = copyOfArray(rhs.savedSolution_,numberColumns_+numberRows_);
  columnScale_ = copyOfArray(rhs.columnScale_,numberColumns_);
  int i;
  for (i=0;i<6;i++) {
    rowArray_[i]=NULL;
    if (rhs.rowArray_[i]) 
      rowArray_[i] = new CoinIndexedVector(*rhs.rowArray_[i]);
    columnArray_[i]=NULL;
    if (rhs.columnArray_[i]) 
      columnArray_[i] = new CoinIndexedVector(*rhs.columnArray_[i]);
  }
  if (rhs.status_) {
    status_ = copyOfArray( rhs.status_,numberColumns_+numberRows_);
  }
  if (rhs.saveStatus_) {
    saveStatus_ = copyOfArray( rhs.saveStatus_,numberColumns_+numberRows_);
  }
  columnPrimalInfeasibility_ = rhs.columnPrimalInfeasibility_;
  columnPrimalSequence_ = rhs.columnPrimalSequence_;
  rowPrimalInfeasibility_ = rhs.rowPrimalInfeasibility_;
  rowPrimalSequence_ = rhs.rowPrimalSequence_;
  columnDualInfeasibility_ = rhs.columnDualInfeasibility_;
  columnDualSequence_ = rhs.columnDualSequence_;
  rowDualInfeasibility_ = rhs.rowDualInfeasibility_;
  rowDualSequence_ = rhs.rowDualSequence_;
  primalToleranceToGetOptimal_ = rhs.primalToleranceToGetOptimal_;
  remainingDualInfeasibility_ = rhs.remainingDualInfeasibility_;
  largeValue_ = rhs.largeValue_;
  largestPrimalError_ = rhs.largestPrimalError_;
  largestDualError_ = rhs.largestDualError_;
  largestSolutionError_ = rhs.largestSolutionError_;
  dualBound_ = rhs.dualBound_;
  alpha_ = rhs.alpha_;
  theta_ = rhs.theta_;
  lowerIn_ = rhs.lowerIn_;
  valueIn_ = rhs.valueIn_;
  upperIn_ = rhs.upperIn_;
  dualIn_ = rhs.dualIn_;
  sequenceIn_ = rhs.sequenceIn_;
  directionIn_ = rhs.directionIn_;
  lowerOut_ = rhs.lowerOut_;
  valueOut_ = rhs.valueOut_;
  upperOut_ = rhs.upperOut_;
  dualOut_ = rhs.dualOut_;
  sequenceOut_ = rhs.sequenceOut_;
  directionOut_ = rhs.directionOut_;
  pivotRow_ = rhs.pivotRow_;
  numberRefinements_ = rhs.numberRefinements_;
  dualTolerance_ = rhs.dualTolerance_;
  primalTolerance_ = rhs.primalTolerance_;
  sumDualInfeasibilities_ = rhs.sumDualInfeasibilities_;
  numberDualInfeasibilities_ = rhs.numberDualInfeasibilities_;
  numberDualInfeasibilitiesWithoutFree_ = 
    rhs.numberDualInfeasibilitiesWithoutFree_;
  sumPrimalInfeasibilities_ = rhs.sumPrimalInfeasibilities_;
  numberPrimalInfeasibilities_ = rhs.numberPrimalInfeasibilities_;
  dualRowPivot_ = rhs.dualRowPivot_->clone(true);
  primalColumnPivot_ = rhs.primalColumnPivot_->clone(true);
  scalingFlag_ = rhs.scalingFlag_;
  numberTimesOptimal_ = rhs.numberTimesOptimal_;
  changeMade_ = rhs.changeMade_;
  algorithm_ = rhs.algorithm_;
  forceFactorization_ = rhs.forceFactorization_;
  perturbation_ = rhs.perturbation_;
  infeasibilityCost_ = rhs.infeasibilityCost_;
  specialOptions_ = rhs.specialOptions_;
  lastBadIteration_ = rhs.lastBadIteration_;
  numberFake_ = rhs.numberFake_;
  if (rhs.nonLinearCost_!=NULL)
    nonLinearCost_ = new ClpNonLinearCost(*rhs.nonLinearCost_);
}
// type == 0 do everything
void 
ClpSimplex::gutsOfDelete(int type)
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
  objectiveWork_=NULL;
  rowObjectiveWork_=NULL;
  delete [] dj_;
  dj_=NULL;
  reducedCostWork_=NULL;
  rowReducedCost_=NULL;
  delete [] solution_;
  solution_=NULL;
  rowActivityWork_=NULL;
  columnActivityWork_=NULL;
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] savedSolution_;
  savedSolution_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
  delete nonLinearCost_;
  nonLinearCost_ = NULL;
  int i;
  for (i=0;i<6;i++) {
    delete rowArray_[i];
    rowArray_[i]=NULL;
    delete columnArray_[i];
    columnArray_[i]=NULL;
  }
  delete rowCopy_;
  rowCopy_=NULL;
  if (!type) {
    // delete everything
    delete [] pivotVariable_;
    pivotVariable_=NULL;
    delete factorization_;
    factorization_ = NULL;
    delete dualRowPivot_;
    dualRowPivot_ = NULL;
    delete primalColumnPivot_;
    primalColumnPivot_ = NULL;
    delete [] saveStatus_;
    saveStatus_=NULL;
    delete status_;
    status_=NULL;
  }
}
// This sets largest infeasibility and most infeasible
void 
ClpSimplex::checkPrimalSolution(const double * rowActivities,
					const double * columnActivities)
{
  double * solution;
  int iRow,iColumn;

  objectiveValue_ = 0.0;
  // now look at primal solution
  columnPrimalInfeasibility_=0.0;
  columnPrimalSequence_=-1;
  rowPrimalInfeasibility_=0.0;
  rowPrimalSequence_=-1;
  largestSolutionError_=0.0;
  solution = rowActivityWork_;
  sumPrimalInfeasibilities_=0.0;
  numberPrimalInfeasibilities_=0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double infeasibility=0.0;
    objectiveValue_ += solution[iRow]*rowObjectiveWork_[iRow];
    if (solution[iRow]>rowUpperWork_[iRow]) {
      infeasibility=solution[iRow]-rowUpperWork_[iRow];
    } else if (solution[iRow]<rowLowerWork_[iRow]) {
      infeasibility=rowLowerWork_[iRow]-solution[iRow];
    }
    if (infeasibility>primalTolerance_) {
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
      numberPrimalInfeasibilities_ ++;
    }
    if (infeasibility>rowPrimalInfeasibility_) {
      rowPrimalInfeasibility_=infeasibility;
      rowPrimalSequence_=iRow;
    }
    infeasibility = fabs(rowActivities[iRow]-solution[iRow]);
    if (infeasibility>largestSolutionError_)
      largestSolutionError_=infeasibility;
    if (rowLowerWork_[iRow]!=rowUpperWork_[iRow])
      clearFixed(iRow+numberColumns_);
    else
      setFixed(iRow+numberColumns_);
  }
  solution = columnActivityWork_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double infeasibility=0.0;
    objectiveValue_ += objectiveWork_[iColumn]*solution[iColumn];
    if (solution[iColumn]>columnUpperWork_[iColumn]) {
      infeasibility=solution[iColumn]-columnUpperWork_[iColumn];
    } else if (solution[iColumn]<columnLowerWork_[iColumn]) {
      infeasibility=columnLowerWork_[iColumn]-solution[iColumn];
    }
    if (infeasibility>columnPrimalInfeasibility_) {
      columnPrimalInfeasibility_=infeasibility;
      columnPrimalSequence_=iColumn;
    }
    if (infeasibility>primalTolerance_) {
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
      numberPrimalInfeasibilities_ ++;
    }
    infeasibility = fabs(columnActivities[iColumn]-solution[iColumn]);
    if (infeasibility>largestSolutionError_)
      largestSolutionError_=infeasibility;
    if (columnUpperWork_[iColumn]-columnLowerWork_[iColumn]>primalTolerance_)
      clearFixed(iColumn);
    else
      setFixed(iColumn);
  }
}
void 
ClpSimplex::checkDualSolution()
{

  double * solution;
  int iRow,iColumn;
  sumDualInfeasibilities_=0.0;
  numberDualInfeasibilities_=0;
  numberDualInfeasibilitiesWithoutFree_=0;
  columnDualInfeasibility_=0.0;
  columnDualSequence_=-1;
  rowDualInfeasibility_=0.0;
  rowDualSequence_=-1;
  primalToleranceToGetOptimal_=max(rowPrimalInfeasibility_,
				   columnPrimalInfeasibility_);
  remainingDualInfeasibility_=0.0;
  solution = rowActivityWork_;
  double dualTolerance=dualTolerance_;
  if (algorithm_>0) {
    // primal
    // we can't really trust infeasibilities if there is dual error
    if (largestDualError_>1.0e-6)
      dualTolerance *= largestDualError_/1.0e-6;
  }
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (getRowStatus(iRow) != ClpSimplex::basic) {
      // not basic
      double value = rowReducedCost_[iRow];
      double distance;
      distance = rowUpperWork_[iRow]-solution[iRow];
      if (distance>primalTolerance_) {
	// should not be negative
	if (value<0.0) {
	  value = - value;
	  if (value>rowDualInfeasibility_) {
	    rowDualInfeasibility_=value;
	    rowDualSequence_=iRow;
	  }
	  if (value>dualTolerance) {
	    sumDualInfeasibilities_ += value-dualTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getRowStatus(iRow) != ClpSimplex::isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distance<largeValue_) {
	      if (distance>primalToleranceToGetOptimal_)
		primalToleranceToGetOptimal_=distance;
	    } else {
	      //gap too big for any tolerance
	      remainingDualInfeasibility_=
		max(remainingDualInfeasibility_,value);
	    }
	  }
	}
      }
      distance = solution[iRow] -rowLowerWork_[iRow];
      if (distance>primalTolerance_) {
	// should not be positive
	if (value>0.0) {
	  if (value>rowDualInfeasibility_) {
	    rowDualInfeasibility_=value;
	    rowDualSequence_=iRow;
	  }
	  if (value>dualTolerance) {
	    sumDualInfeasibilities_ += value-dualTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getRowStatus(iRow) != ClpSimplex::isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distance<largeValue_&&
		distance>primalToleranceToGetOptimal_)
	      primalToleranceToGetOptimal_=distance;
	  }
	}
      }
    }
  }
  solution = columnActivityWork_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (getColumnStatus(iColumn) != ClpSimplex::basic) {
      // not basic
      double value = reducedCostWork_[iColumn];
      double distance;
      distance = columnUpperWork_[iColumn]-solution[iColumn];
      if (distance>primalTolerance_) {
	// should not be negative
	if (value<0.0) {
	  value = - value;
	  if (value>columnDualInfeasibility_) {
	    columnDualInfeasibility_=value;
	    columnDualSequence_=iColumn;
	  }
	  if (value>dualTolerance) {
	    sumDualInfeasibilities_ += value-dualTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getColumnStatus(iColumn) != ClpSimplex::isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distance<largeValue_) {
	      if (distance>primalToleranceToGetOptimal_)
		primalToleranceToGetOptimal_=distance;
	    } else {
	      //gap too big for any tolerance
	      remainingDualInfeasibility_=
		max(remainingDualInfeasibility_,value);
	    }
	  }
	}
      }
      distance = solution[iColumn] -columnLowerWork_[iColumn];
      if (distance>primalTolerance_) {
	// should not be positive
	if (value>0.0) {
	  if (value>columnDualInfeasibility_) {
	    columnDualInfeasibility_=value;
	    columnDualSequence_=iColumn;
	  }
	  if (value>dualTolerance) {
	    sumDualInfeasibilities_ += value-dualTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getColumnStatus(iColumn) != ClpSimplex::isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distance<largeValue_&&
		distance>primalToleranceToGetOptimal_)
	      primalToleranceToGetOptimal_=distance;
	  }
	}
      }
    }
  }
}
/*
  Unpacks one column of the matrix into indexed array 
*/
void 
ClpSimplex::unpack(CoinIndexedVector * rowArray)
{
  rowArray->clear();
  if (sequenceIn_>=numberColumns_) {
    //slack
    rowArray->insert(sequenceIn_-numberColumns_,-1.0);
  } else {
    // column
    matrix_->unpack(this,rowArray,sequenceIn_);
  }
}
void 
ClpSimplex::unpack(CoinIndexedVector * rowArray,int sequence)
{
  rowArray->clear();
  if (sequence>=numberColumns_) {
    //slack
    rowArray->insert(sequence-numberColumns_,-1.0);
  } else {
    // column
    matrix_->unpack(this,rowArray,sequence);
  }
}
void
ClpSimplex::createRim(int what,bool makeRowCopy)
{
  if ((what&16)!=0) {
    // move information to work arrays
    // row reduced costs
    if (!dj_) {
      dj_ = new double[numberRows_+numberColumns_];
      reducedCostWork_ = dj_;
      rowReducedCost_ = dj_+numberColumns_;
      memcpy(reducedCostWork_,reducedCost_,numberColumns_*sizeof(double));
    }
    if (!solution_) {
      solution_ = new double[numberRows_+numberColumns_];
      columnActivityWork_ = solution_;
      rowActivityWork_ = solution_+numberColumns_;
      memcpy(columnActivityWork_,columnActivity_,
	     numberColumns_*sizeof(double));
      memcpy(rowActivityWork_,rowActivity_,
	     numberRows_*sizeof(double));
    }
    
    if (makeRowCopy) {
      delete rowCopy_;
      // may return NULL if can't give row copy
      rowCopy_ = matrix_->reverseOrderedCopy();
    }
  }
  int i;
  if ((what&4)!=0) {
    delete [] cost_;
    cost_ = new double[numberColumns_+numberRows_];
    objectiveWork_ = cost_;
    rowObjectiveWork_ = cost_+numberColumns_;
    memcpy(objectiveWork_,objective_,numberColumns_*sizeof(double));
    if (rowObjective_)
      memcpy(rowObjectiveWork_,rowObjective_,numberRows_*sizeof(double));
    else
      memset(rowObjectiveWork_,0,numberRows_*sizeof(double));
  }
  if ((what&1)!=0) {
    delete [] lower_;
    delete [] upper_;
    lower_ = new double[numberColumns_+numberRows_];
    upper_ = new double[numberColumns_+numberRows_];
    rowLowerWork_ = lower_+numberColumns_;
    columnLowerWork_ = lower_;
    rowUpperWork_ = upper_+numberColumns_;
    columnUpperWork_ = upper_;
    memcpy(rowLowerWork_,rowLower_,numberRows_*sizeof(double));
    memcpy(rowUpperWork_,rowUpper_,numberRows_*sizeof(double));
    memcpy(columnLowerWork_,columnLower_,numberColumns_*sizeof(double));
    memcpy(columnUpperWork_,columnUpper_,numberColumns_*sizeof(double));
    // clean up any mismatches on infinity
    for (i=0;i<numberColumns_;i++) {
      if (columnLowerWork_[i]<-CLP_INFINITY)
	columnLowerWork_[i] = -DBL_MAX;
      if (columnUpperWork_[i]>CLP_INFINITY)
	columnUpperWork_[i] = DBL_MAX;
    }
    // clean up any mismatches on infinity
    for (i=0;i<numberRows_;i++) {
      if (rowLowerWork_[i]<-1.0e30)
	rowLowerWork_[i] = -DBL_MAX;
      if (rowUpperWork_[i]>1.0e30)
	rowUpperWork_[i] = DBL_MAX;
    }
  }
  // do scaling if needed
  if (scalingFlag_>0&&!rowScale_) {
    if (matrix_->scale(this))
      scalingFlag_=-scalingFlag_; // not scaled after all
  }
  if ((what&4)!=0) {
    double direction = optimizationDirection_;
    // direction is actually scale out not scale in
    if (direction)
      direction = 1.0/direction;
    // but also scale by scale factors
    // not really sure about row scaling
    if (!rowScale_) {
      if (direction!=1.0) {
	for (i=0;i<numberRows_;i++)
	  rowObjectiveWork_[i] *= direction;
	for (i=0;i<numberColumns_;i++)
	  objectiveWork_[i] *= direction;
      }
    } else {
      for (i=0;i<numberRows_;i++)
	rowObjectiveWork_[i] *= direction/rowScale_[i];
      for (i=0;i<numberColumns_;i++)
	objectiveWork_[i] *= direction*columnScale_[i];
    }
  }
  if ((what&1)!=0&&rowScale_) {
    for (i=0;i<numberColumns_;i++) {
      double multiplier = 1.0/columnScale_[i];
      if (columnLowerWork_[i]>-1.0e50)
	columnLowerWork_[i] *= multiplier;
      if (columnUpperWork_[i]<1.0e50)
	columnUpperWork_[i] *= multiplier;
    }
    for (i=0;i<numberRows_;i++) {
      if (rowLowerWork_[i]>-1.0e50)
	rowLowerWork_[i] *= rowScale_[i];
      if (rowUpperWork_[i]<1.0e50)
	rowUpperWork_[i] *= rowScale_[i];
    }
  }
  if ((what&8)!=0&&rowScale_) {
    // on entry
    for (i=0;i<numberColumns_;i++) {
      columnActivityWork_[i] /= columnScale_[i];
      reducedCostWork_[i] *= columnScale_[i];
    }
    for (i=0;i<numberRows_;i++) {
      rowActivityWork_[i] *= rowScale_[i];
      dual_[i] /= rowScale_[i];
    }
  }
 
  // we need to treat matrix as if each element by rowScaleIn and columnScaleout??
  // maybe we need to move scales to SimplexModel for factorization?
  if ((what&8)!=0&&!pivotVariable_) {
    pivotVariable_=new int[numberRows_];
  }

  if ((what&16)!=0&&!rowArray_[2]) {
    // get some arrays
    int iRow,iColumn;
    // these are "indexed" arrays so we always know where nonzeros are
    /**********************************************************
      rowArray_[3] is long enough for columns as well
    *********************************************************/
    for (iRow=0;iRow<4;iRow++) {
      if (!rowArray_[iRow]) {
	rowArray_[iRow]=new CoinIndexedVector();
	int length =numberRows_+factorization_->maximumPivots();
	if (iRow==3)
	  length = max(length,numberColumns_);
	rowArray_[iRow]->reserve(length);
      }
    }
    
    for (iColumn=0;iColumn<2;iColumn++) {
      if (!columnArray_[iColumn]) {
	columnArray_[iColumn]=new CoinIndexedVector();
	columnArray_[iColumn]->reserve(numberColumns_);
      }
    }    
  }

}
void
ClpSimplex::deleteRim()
{
  int i;
  if (problemStatus_!=1&&problemStatus_!=2) {
    delete [] ray_;
    ray_=NULL;
  }
  if (rowScale_) {
    for (i=0;i<numberColumns_;i++) {
      columnActivity_[i] = columnActivityWork_[i]*columnScale_[i];
      reducedCost_[i] = reducedCostWork_[i]/columnScale_[i];
    }
    for (i=0;i<numberRows_;i++) {
      rowActivity_[i] = rowActivityWork_[i]/rowScale_[i];
      dual_[i] *= rowScale_[i];
    }
    if (problemStatus_==2) {
      for (i=0;i<numberColumns_;i++) {
	ray_[i] *= columnScale_[i];
      }
    } else if (problemStatus_==1) {
      for (i=0;i<numberRows_;i++) {
	ray_[i] *= rowScale_[i];
      }
    }
  } else {
    for (i=0;i<numberColumns_;i++) {
      columnActivity_[i] = columnActivityWork_[i];
      reducedCost_[i] = reducedCostWork_[i];
    }
    for (i=0;i<numberRows_;i++) {
      rowActivity_[i] = rowActivityWork_[i];
    }
  }
  // direction may have been modified by scaling - clean up
  if (optimizationDirection_>0.0)
    optimizationDirection_ = 1.0;
  else if (optimizationDirection_<0.0)
    optimizationDirection_ = -1.0;
  // scaling may have been turned off
  scalingFlag_ = abs(scalingFlag_);
  gutsOfDelete(1);
}
void 
ClpSimplex::setDualBound(double value)
{
  if (value>0.0)
    dualBound_=value;
}
void 
ClpSimplex::setInfeasibilityCost(double value)
{
  if (value>0.0)
    infeasibilityCost_=value;
}
void ClpSimplex::setnumberRefinements( int value) 
{
  if (value>=0&&value<10)
    numberRefinements_=value;
}
// Sets row pivot choice algorithm in dual
void 
ClpSimplex::setDualRowPivotAlgorithm(ClpDualRowPivot & choice)
{
  delete dualRowPivot_;
  dualRowPivot_ = choice.clone(true);
}
// Sets row pivot choice algorithm in dual
void 
ClpSimplex::setPrimalColumnPivotAlgorithm(ClpPrimalColumnPivot & choice)
{
  delete primalColumnPivot_;
  primalColumnPivot_ = choice.clone(true);
}
// Sets or unsets scaling, 0 -off, 1 on, 2 dynamic(later)
void 
ClpSimplex::scaling(int mode)
{
  if (mode>0&&mode<3) {
    scalingFlag_=mode;
  } else if (!mode) {
    scalingFlag_=0;
    delete [] rowScale_;
    rowScale_ = NULL;
    delete [] columnScale_;
    columnScale_ = NULL;
  }
}
// Sets up basis
void 
ClpSimplex::setBasis ( const CoinWarmStartBasis & basis)
{
  // transform basis to status arrays
  int iRow,iColumn;
  if (!status_) {
    /*
      get status arrays
      CoinWarmStartBasis would seem to have overheads and we will need
      extra bits anyway.
    */
    status_ = new unsigned char [numberColumns_+numberRows_];
    memset(status_,0,(numberColumns_+numberRows_)*sizeof(char));
  }
  CoinWarmStartBasis basis2 = basis;
  // resize if necessary
  basis2.resize(numberRows_,numberColumns_);
  // move status
  for (iRow=0;iRow<numberRows_;iRow++) {
    setRowStatus(iRow,
		 (Status) basis2.getArtifStatus(iRow));
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    setColumnStatus(iColumn,
		    (Status) basis2.getStructStatus(iColumn));
  }
}
// Passes in factorization
void 
ClpSimplex::setFactorization( ClpFactorization & factorization)
{
  delete factorization_;
  factorization_= new ClpFactorization(factorization);
}
// Warm start
CoinWarmStartBasis  
ClpSimplex::getBasis() const
{
  int iRow,iColumn;
  CoinWarmStartBasis basis;
  basis.setSize(numberColumns_,numberRows_);

  if(status_) {
    for (iRow=0;iRow<numberRows_;iRow++) {
      basis.setArtifStatus(iRow,
			   (CoinWarmStartBasis::Status) getRowStatus(iRow));
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      basis.setStructStatus(iColumn,
		       (CoinWarmStartBasis::Status) getColumnStatus(iColumn));
    }
  }
  return basis;
}
void 
ClpSimplex::times(double scalar,
		  const double * x, double * y) const
{
  if (rowScale_)
    matrix_->times(scalar,x,y,rowScale_,columnScale_);
  else
    matrix_->times(scalar,x,y);
}
void 
ClpSimplex::transposeTimes(double scalar,
			   const double * x, double * y) const 
{
  if (rowScale_)
    matrix_->transposeTimes(scalar,x,y,rowScale_,columnScale_);
  else
    matrix_->transposeTimes(scalar,x,y);
}
/* Perturbation:
   -50 to +50 - perturb by this power of ten (-6 sounds good)
   100 - auto perturb if takes too long (1.0e-6 largest nonzero)
   101 - we are perturbed
   102 - don't try perturbing again
   default is 100
*/
void 
ClpSimplex::setPerturbation(int value)
{
  if(value<=100&&value >=-50) {
    perturbation_=value;
  } 
}
// Sparsity on or off
bool 
ClpSimplex::sparseFactorization() const
{
  return factorization_->sparseThreshold();
}
void 
ClpSimplex::setSparseFactorization(bool value)
{
  if (value) {
    if (!factorization_->sparseThreshold())
      factorization_->goSparse();
  } else {
    factorization_->sparseThreshold(0);
  }
}
/* Tightens primal bounds to make dual faster.  Unless
   fixed, bounds are slightly looser than they could be.
   This is to make dual go faster and is probably not needed
   with a presolve.  Returns non-zero if problem infeasible
*/
int 
ClpSimplex::tightenPrimalBounds()
{
  
  // Get a row copy in standard format
  CoinPackedMatrix copy;
  copy.reverseOrderedCopyOf(*matrix());
  // get matrix data pointers
  const int * column = copy.getIndices();
  const int * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  const double * element = copy.getElements();
  int numberChanged=1,iPass=0;
  double large = largeValue(); // treat bounds > this as infinite
  int numberInfeasible=0;
  int totalTightened = 0;

  double tolerance = primalTolerance();


  // Save column bounds
  double * saveLower = new double [numberColumns_];
  memcpy(saveLower,columnLower_,numberColumns_*sizeof(double));
  double * saveUpper = new double [numberColumns_];
  memcpy(saveUpper,columnUpper_,numberColumns_*sizeof(double));
#define MAXPASS 10  

  // Loop round seeing if we can tighten bounds
  // Would be faster to have a stack of possible rows
  // and we put altered rows back on stack

  int iRow, iColumn;

  while(numberChanged) {

    numberChanged = 0; // Bounds tightened this pass
    
    if (iPass==MAXPASS) break;
    iPass++;
    
    for (iRow = 0; iRow < numberRows_; iRow++) {

      if (rowLower_[iRow]>-large||rowUpper_[iRow]<large) {

	// possible row
	int infiniteUpper = 0;
	int infiniteLower = 0;
	double maximumUp = 0.0;
	double maximumDown = 0.0;
	double newBound;
	int rStart = rowStart[iRow];
	int rEnd = rowStart[iRow]+rowLength[iRow];
	int j;

	// Compute possible lower and upper ranges

	for (j = rStart; j < rEnd; ++j) {
	  double value=element[j];
	  iColumn = column[j];
	  if (value > 0.0) {
	    if (columnUpper_[iColumn] >= large) {
	      maximumUp = 1e31;
	      ++infiniteUpper;
	    } else {
	      maximumUp += columnUpper_[iColumn] * value;
	    }
	    if (columnLower_[iColumn] <= -large) {
	      maximumDown = -1e31;
	      ++infiniteLower;
	    } else {
	      maximumDown += columnLower_[iColumn] * value;
	    }
	  } else if (value<0.0) {
	    if (columnUpper_[iColumn] >= large) {
	      maximumDown = -1e31;
	      ++infiniteLower;
	    } else {
	      maximumDown += columnUpper_[iColumn] * value;
	    }
	    if (columnLower_[iColumn] <= -large) {
	      maximumUp = 1e31;
	      ++infiniteUpper;
	    } else {
	      maximumUp += columnLower_[iColumn] * value;
	    }
	  }
	}
	if (maximumUp <= rowUpper_[iRow] + tolerance && 
	    maximumDown >= rowLower_[iRow] - tolerance) {

	  // Row is redundant - make totally free
	  rowLower_[iRow]=-DBL_MAX;
	  rowUpper_[iRow]=DBL_MAX;

	} else {
	  if (maximumUp < rowLower_[iRow] -tolerance ||
	      maximumDown > rowUpper_[iRow]+tolerance) {
	    // problem is infeasible - exit at once
	    numberInfeasible++;
	    break;
	  }

	  if (infiniteUpper == 0 && rowLower_[iRow] > -large) {
	    for (j = rStart; j < rEnd; ++j) {
	      double value=element[j];
	      iColumn = column[j];
	      if (value > 0.0) {
		if (columnUpper_[iColumn] < large) {
		  newBound = columnUpper_[iColumn] + 
		    (rowLower_[iRow] - maximumUp) / value;
		  if (newBound > columnLower_[iColumn] + 1.0e-12) {
		    // Tighten the lower bound 

		    columnLower_[iColumn] = newBound;
		    ++numberChanged;

		    // check infeasible (relaxed)
		    if (columnUpper_[iColumn] - columnLower_[iColumn] < 
			-100.0*tolerance) 
		      numberInfeasible++;
		    infiniteLower=1; // skip looking at other bound
		  }
		}
	      } else {
		if (columnLower_[iColumn] > -large) {
		  newBound = columnLower_[iColumn] + 
		    (rowLower_[iRow] - maximumUp) / value;
		  if (newBound < columnUpper_[iColumn] - 1.0e-12) {
		    // Tighten the upper bound 

		    columnUpper_[iColumn] = newBound;
		    ++numberChanged;

		    // check infeasible (relaxed)
		    if (columnUpper_[iColumn] - columnLower_[iColumn] < 
			-100.0*tolerance) 
		      numberInfeasible++;
		    infiniteLower=1; // skip looking at other bound
		  }
		}
	      }
	    }
	  }
	  
	  // Try other way
	  if (infiniteLower == 0 && rowUpper_[iRow] < large) {
	    for (j = rStart; j < rEnd; ++j) {
	      double value=element[j];
	      iColumn = column[j];
	      if (value < 0.0) {
		if (columnUpper_[iColumn] < large) {
		  newBound = columnUpper_[iColumn] + 
		    (rowUpper_[iRow] - maximumDown) / value;
		  if (newBound > columnLower_[iColumn] + 1.0e-12) {
		    // Tighten the lower bound 

		    columnLower_[iColumn] = newBound;
		    ++numberChanged;

		    // check infeasible (relaxed)
		    if (columnUpper_[iColumn] - columnLower_[iColumn] < 
			-100.0*tolerance) 
		      numberInfeasible++;
		  }
		} 
	      } else {
		if (columnLower_[iColumn] > -large) {
		  newBound = columnLower_[iColumn] + 
		    (rowUpper_[iRow] - maximumDown) / value;
		  if (newBound < columnUpper_[iColumn] - 1.0e-12) {
		    // Tighten the upper bound 

		    columnUpper_[iColumn] = newBound;
		    ++numberChanged;

		    // check infeasible (relaxed)
		    if (columnUpper_[iColumn] - columnLower_[iColumn] < 
			-100.0*tolerance) 
		      numberInfeasible++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    totalTightened += numberChanged;
    if (numberInfeasible) break;
  }
  if (!numberInfeasible) {
    handler_->message(CLP_SIMPLEX_BOUNDTIGHTEN,messages_)
      <<totalTightened
      <<CoinMessageEol;
    // Set bounds slightly loose
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (columnUpper_[iColumn]-columnLower_[iColumn]<tolerance) {
	// fix
	if (fabs(columnUpper_[iColumn])<fabs(columnLower_[iColumn])) 
	  columnLower_[iColumn]=columnUpper_[iColumn];
	else
	  columnUpper_[iColumn]=columnLower_[iColumn];
      } else {
	if (columnUpper_[iColumn]<saveUpper[iColumn]) {
	  // relax a bit
	  columnUpper_[iColumn] = min(columnUpper_[iColumn]+100.0*tolerance,
				      saveUpper[iColumn]);
	}
	if (columnLower_[iColumn]>saveLower[iColumn]) {
	  // relax a bit
	  columnLower_[iColumn] = max(columnLower_[iColumn]-100.0*tolerance,
				      saveLower[iColumn]);
	}
      }
    }
  } else {
    handler_->message(CLP_SIMPLEX_INFEASIBILITIES,messages_)
      <<numberInfeasible
      <<CoinMessageEol;
    // restore column bounds
    memcpy(columnLower_,saveLower,numberColumns_*sizeof(double));
    memcpy(columnUpper_,saveUpper,numberColumns_*sizeof(double));
  }
  delete [] saveLower;
  delete [] saveUpper;
  return (numberInfeasible);
}
// dual 
#include "ClpSimplexDual.hpp"
int ClpSimplex::dual ( )
{
  return ((ClpSimplexDual *) this)->dual();
}
// primal 
#include "ClpSimplexPrimal.hpp"
int ClpSimplex::primal (int ifValuesPass )
{
  return ((ClpSimplexPrimal *) this)->primal(ifValuesPass);
}
/* For strong branching.  On input lower and upper are new bounds
   while on output they are objective function values (>1.0e50 infeasible).
   Return code is 0 if nothing interesting, -1 if infeasible both
   ways and +1 if infeasible one way (check values to see which one(s))
*/
int ClpSimplex::strongBranching(int numberVariables,const int * variables,
				double * newLower, double * newUpper,
				bool stopOnFirstInfeasible,
				bool alwaysFinish)
{
  return ((ClpSimplexDual *) this)->strongBranching(numberVariables,variables,
						    newLower,  newUpper,
						    stopOnFirstInfeasible,
						    alwaysFinish);
}
/* Borrow model.  This is so we dont have to copy large amounts
   of data around.  It assumes a derived class wants to overwrite
   an empty model with a real one - while it does an algorithm.
   This is same as ClpModel one, but sets scaling on etc. */
void 
ClpSimplex::borrowModel(ClpModel & otherModel) 
{
  ClpModel::borrowModel(otherModel);
  scaling();
  ClpDualRowSteepest steep1;
  setDualRowPivotAlgorithm(steep1);
  ClpPrimalColumnSteepest steep2;
  setPrimalColumnPivotAlgorithm(steep2);
}
typedef struct {
  double optimizationDirection;
  double dblParam[ClpLastDblParam];
  double objectiveValue;
  double dualBound;
  double dualTolerance;
  double primalTolerance;
  double sumDualInfeasibilities;
  double sumPrimalInfeasibilities;
  double infeasibilityCost;
  int numberRows;
  int numberColumns;
  int intParam[ClpLastIntParam];
  int numberIterations;
  int problemStatus;
  int maximumIterations;
  int lengthNames;
  int numberDualInfeasibilities;
  int numberDualInfeasibilitiesWithoutFree;
  int numberPrimalInfeasibilities;
  int numberRefinements;
  int scalingFlag;
  int algorithm;
  unsigned int specialOptions;
  int dualPivotChoice;
  int primalPivotChoice;
  int matrixStorageChoice;
} Clp_scalars;

int outDoubleArray(double * array, int length, FILE * fp)
{
  int numberWritten;
  if (array&&length) {
    numberWritten = fwrite(&length,sizeof(int),1,fp);
    if (numberWritten!=1)
      return 1;
    numberWritten = fwrite(array,sizeof(double),length,fp);
    if (numberWritten!=length)
      return 1;
  } else {
    length = 0;
    numberWritten = fwrite(&length,sizeof(int),1,fp);
    if (numberWritten!=1)
      return 1;
  }
  return 0;
}
// Save model to file, returns 0 if success
int
ClpSimplex::saveModel(const char * fileName)
{
  FILE * fp = fopen(fileName,"wb");
  if (fp) {
    Clp_scalars scalars;
    int i, numberWritten;
    // Fill in scalars
    scalars.optimizationDirection = optimizationDirection_;
    memcpy(scalars.dblParam, dblParam_,ClpLastDblParam * sizeof(double));
    scalars.objectiveValue = objectiveValue_;
    scalars.dualBound = dualBound_;
    scalars.dualTolerance = dualTolerance_;
    scalars.primalTolerance = primalTolerance_;
    scalars.sumDualInfeasibilities = sumDualInfeasibilities_;
    scalars.sumPrimalInfeasibilities = sumPrimalInfeasibilities_;
    scalars.infeasibilityCost = infeasibilityCost_;
    scalars.numberRows = numberRows_;
    scalars.numberColumns = numberColumns_;
    memcpy(scalars.intParam, intParam_,ClpLastIntParam * sizeof(double));
    scalars.numberIterations = numberIterations_;
    scalars.problemStatus = problemStatus_;
    scalars.maximumIterations = maximumIterations_;
    scalars.lengthNames = lengthNames_;
    scalars.numberDualInfeasibilities = numberDualInfeasibilities_;
    scalars.numberDualInfeasibilitiesWithoutFree 
      = numberDualInfeasibilitiesWithoutFree_;
    scalars.numberPrimalInfeasibilities = numberPrimalInfeasibilities_;
    scalars.numberRefinements = numberRefinements_;
    scalars.scalingFlag = scalingFlag_;
    scalars.algorithm = algorithm_;
    scalars.specialOptions = specialOptions_;
    scalars.dualPivotChoice = dualRowPivot_->type();
    scalars.primalPivotChoice = primalColumnPivot_->type();
    scalars.matrixStorageChoice = matrix_->type();

    // put out scalars
    numberWritten = fwrite(&scalars,sizeof(Clp_scalars),1,fp);
    if (numberWritten!=1)
      return 1;
    // strings
    int length;
    for (i=0;i<ClpLastStrParam;i++) {
      length = strParam_[i].size();
      numberWritten = fwrite(&length,sizeof(int),1,fp);
      if (numberWritten!=1)
	return 1;
      numberWritten = fwrite(strParam_[i].c_str(),length,1,fp);
      if (numberWritten!=1)
	return 1;
    }
    // arrays - in no particular order
    if (outDoubleArray(rowActivity_,numberRows_,fp))
	return 1;
    if (outDoubleArray(columnActivity_,numberColumns_,fp))
	return 1;
    if (outDoubleArray(dual_,numberRows_,fp))
	return 1;
    if (outDoubleArray(reducedCost_,numberColumns_,fp))
	return 1;
    if (outDoubleArray(rowLower_,numberRows_,fp))
	return 1;
    if (outDoubleArray(rowUpper_,numberRows_,fp))
	return 1;
    if (outDoubleArray(objective_,numberColumns_,fp))
	return 1;
    if (outDoubleArray(rowObjective_,numberRows_,fp))
	return 1;
    if (outDoubleArray(columnLower_,numberColumns_,fp))
	return 1;
    if (outDoubleArray(columnUpper_,numberColumns_,fp))
	return 1;
    if (ray_) {
      if (problemStatus_==1)
	if (outDoubleArray(ray_,numberRows_,fp))
	  return 1;
      else if (problemStatus_==2)
	if (outDoubleArray(ray_,numberColumns_,fp))
	  return 1;
      else
	if (outDoubleArray(NULL,0,fp))
	  return 1;
    } else {
      if (outDoubleArray(NULL,0,fp))
	return 1;
    }
    if (status_&&(numberRows_+numberColumns_)>0) {
      length = numberRows_+numberColumns_;
      numberWritten = fwrite(&length,sizeof(int),1,fp);
      if (numberWritten!=1)
	return 1;
      numberWritten = fwrite(status_,sizeof(char),length, fp);
      if (numberWritten!=length)
	return 1;
    } else {
      length = 0;
      numberWritten = fwrite(&length,sizeof(int),1,fp);
      if (numberWritten!=1)
	return 1;
    }
    if (lengthNames_) {
      assert (numberRows_ == (int) rowNames_.size());
      for (i=0;i<numberRows_;i++) {
	length = rowNames_[i].size();
	numberWritten = fwrite(&length,sizeof(int),1,fp);
	if (numberWritten!=1)
	  return 1;
	numberWritten = fwrite(rowNames_[i].c_str(),length,1,fp);
	if (numberWritten!=1)
	  return 1;
      }
      assert (numberColumns_ == (int) columnNames_.size());
      for (i=0;i<numberColumns_;i++) {
	length = columnNames_[i].size();
	numberWritten = fwrite(&length,sizeof(int),1,fp);
	if (numberWritten!=1)
	  return 1;
	numberWritten = fwrite(columnNames_[i].c_str(),length,1,fp);
	if (numberWritten!=1)
	  return 1;
      }
    }
    // just standard type at present
    assert (matrix_->type()==1);
    assert (matrix_->getNumCols() == numberColumns_);
    assert (matrix_->getNumRows() == numberRows_);
    // we are going to save with gaps
    length = matrix_->getVectorStarts()[numberColumns_-1]
      + matrix_->getVectorLengths()[numberColumns_-1];
    numberWritten = fwrite(&length,sizeof(int),1,fp);
    if (numberWritten!=1)
      return 1;
    numberWritten = fwrite(matrix_->getElements(),
			   sizeof(double),length,fp);
    if (numberWritten!=length)
      return 1;
    numberWritten = fwrite(matrix_->getIndices(),
			   sizeof(int),length,fp);
    if (numberWritten!=length)
      return 1;
    numberWritten = fwrite(matrix_->getVectorStarts(),
			   sizeof(int),numberColumns_,fp);
    if (numberWritten!=numberColumns_)
      return 1;
    numberWritten = fwrite(matrix_->getVectorLengths(),
			   sizeof(int),numberColumns_,fp);
    if (numberWritten!=numberColumns_)
      return 1;
    // finished
    fclose(fp);
    return 0;
  } else {
    return -1;
  }
}

int inDoubleArray(double * &array, int length, FILE * fp)
{
  int numberRead;
  int length2;
  numberRead = fread(&length2,sizeof(int),1,fp);
  if (numberRead!=1)
    return 1;
  if (length2) {
    // lengths must match
    if (length!=length2)
      return 2;
    array = new double[length];
    numberRead = fread(array,sizeof(double),length,fp);
    if (numberRead!=length)
      return 1;
  } 
  return 0;
}
/* Restore model from file, returns 0 if success,
   deletes current model */
int 
ClpSimplex::restoreModel(const char * fileName)
{
  FILE * fp = fopen(fileName,"rb");
  if (fp) {
    // Get rid of current model
    ClpModel::gutsOfDelete();
    gutsOfDelete(0);
    int i;
    for (i=0;i<6;i++) {
      rowArray_[i]=NULL;
      columnArray_[i]=NULL;
    }
    // get an empty factorization so we can set tolerances etc
    factorization_ = new ClpFactorization();
    // Say sparse
    factorization_->sparseThreshold(1);
    Clp_scalars scalars;
    int numberRead;

    // get scalars
    numberRead = fread(&scalars,sizeof(Clp_scalars),1,fp);
    if (numberRead!=1)
      return 1;
    // Fill in scalars
    optimizationDirection_ = scalars.optimizationDirection;
    memcpy(dblParam_, scalars.dblParam, ClpLastDblParam * sizeof(double));
    objectiveValue_ = scalars.objectiveValue;
    dualBound_ = scalars.dualBound;
    dualTolerance_ = scalars.dualTolerance;
    primalTolerance_ = scalars.primalTolerance;
    sumDualInfeasibilities_ = scalars.sumDualInfeasibilities;
    sumPrimalInfeasibilities_ = scalars.sumPrimalInfeasibilities;
    infeasibilityCost_ = scalars.infeasibilityCost;
    numberRows_ = scalars.numberRows;
    numberColumns_ = scalars.numberColumns;
    memcpy(intParam_, scalars.intParam,ClpLastIntParam * sizeof(double));
    numberIterations_ = scalars.numberIterations;
    problemStatus_ = scalars.problemStatus;
    maximumIterations_ = scalars.maximumIterations;
    lengthNames_ = scalars.lengthNames;
    numberDualInfeasibilities_ = scalars.numberDualInfeasibilities;
    numberDualInfeasibilitiesWithoutFree_ 
      = scalars.numberDualInfeasibilitiesWithoutFree;
    numberPrimalInfeasibilities_ = scalars.numberPrimalInfeasibilities;
    numberRefinements_ = scalars.numberRefinements;
    scalingFlag_ = scalars.scalingFlag;
    algorithm_ = scalars.algorithm;
    specialOptions_ = scalars.specialOptions;
    // strings
    int length;
    for (i=0;i<ClpLastStrParam;i++) {
      numberRead = fread(&length,sizeof(int),1,fp);
      if (numberRead!=1)
	return 1;
      char * array = new char[length+1];
      numberRead = fread(array,length,1,fp);
      if (numberRead!=1)
	return 1;
      array[length]='\0';
      strParam_[i]=array;
      delete [] array;
    }
    // arrays - in no particular order
    if (inDoubleArray(rowActivity_,numberRows_,fp))
	return 1;
    if (inDoubleArray(columnActivity_,numberColumns_,fp))
	return 1;
    if (inDoubleArray(dual_,numberRows_,fp))
	return 1;
    if (inDoubleArray(reducedCost_,numberColumns_,fp))
	return 1;
    if (inDoubleArray(rowLower_,numberRows_,fp))
	return 1;
    if (inDoubleArray(rowUpper_,numberRows_,fp))
	return 1;
    if (inDoubleArray(objective_,numberColumns_,fp))
	return 1;
    if (inDoubleArray(rowObjective_,numberRows_,fp))
	return 1;
    if (inDoubleArray(columnLower_,numberColumns_,fp))
	return 1;
    if (inDoubleArray(columnUpper_,numberColumns_,fp))
	return 1;
    if (problemStatus_==1) {
      if (inDoubleArray(ray_,numberRows_,fp))
	return 1;
    } else if (problemStatus_==2) {
      if (inDoubleArray(ray_,numberColumns_,fp))
	return 1;
    } else {
      // ray should be null
      numberRead = fread(&length,sizeof(int),1,fp);
      if (numberRead!=1)
	return 1;
      if (length)
	return 2;
    }
    delete [] status_;
    status_=NULL;
    // status region
    numberRead = fread(&length,sizeof(int),1,fp);
    if (numberRead!=1)
	return 1;
    if (length) {
      if (length!=numberRows_+numberColumns_)
	return 1;
      status_ = new char unsigned[length];
      numberRead = fread(status_,sizeof(char),length, fp);
      if (numberRead!=length)
	return 1;
    }
    if (lengthNames_) {
      char * array = new char[lengthNames_+1];
      rowNames_.resize(0);
      for (i=0;i<numberRows_;i++) {
	numberRead = fread(&length,sizeof(int),1,fp);
	if (numberRead!=1)
	  return 1;
	numberRead = fread(array,length,1,fp);
	if (numberRead!=1)
	  return 1;
	rowNames_[i]=array;
      }
      columnNames_.resize(0);
      for (i=0;i<numberColumns_;i++) {
	numberRead = fread(&length,sizeof(int),1,fp);
	if (numberRead!=1)
	  return 1;
	numberRead = fread(array,length,1,fp);
	if (numberRead!=1)
	  return 1;
	columnNames_[i]=array;
      }
      delete [] array;
    }
    // Pivot choices
    assert(scalars.dualPivotChoice>0&&(scalars.dualPivotChoice&63)<3);
    delete dualRowPivot_;
    switch ((scalars.dualPivotChoice&63)) {
    case 1:
      // Dantzig
      dualRowPivot_ = new ClpDualRowDantzig();
      break;
    case 2:
      // Steepest - use mode
      dualRowPivot_ = new ClpDualRowSteepest(scalars.dualPivotChoice>>6);
      break;
    default:
      abort();
    }
    assert(scalars.primalPivotChoice>0&&(scalars.primalPivotChoice&63)<3);
    delete primalColumnPivot_;
    switch ((scalars.primalPivotChoice&63)) {
    case 1:
      // Dantzig
      primalColumnPivot_ = new ClpPrimalColumnDantzig();
      break;
    case 2:
      // Steepest - use mode
      primalColumnPivot_ 
	= new ClpPrimalColumnSteepest(scalars.primalPivotChoice>>6);
      break;
    default:
      abort();
    }
    assert(scalars.matrixStorageChoice==1);
    delete matrix_;
    // get arrays
    numberRead = fread(&length,sizeof(int),1,fp);
    if (numberRead!=1)
      return 1;
    double * elements = new double[length];
    int * indices = new int[length];
    int * starts = new int[numberColumns_];
    int * lengths = new int[numberColumns_];
    numberRead = fread(elements, sizeof(double),length,fp);
    if (numberRead!=length)
      return 1;
    numberRead = fread(indices, sizeof(int),length,fp);
    if (numberRead!=length)
      return 1;
    numberRead = fread(starts, sizeof(int),numberColumns_,fp);
    if (numberRead!=numberColumns_)
      return 1;
    numberRead = fread(lengths, sizeof(int),numberColumns_,fp);
    if (numberRead!=numberColumns_)
      return 1;
    // assign matrix
    CoinPackedMatrix * matrix = new CoinPackedMatrix();
    matrix->assignMatrix(true, numberRows_, numberColumns_,
			 length, elements, indices, starts, lengths);
    // and transfer to Clp
    matrix_ = new ClpPackedMatrix(matrix);
    // finished
    fclose(fp);
    return 0;
  } else {
    return -1;
  }
  return 0;
}
// value of incoming variable (in Dual)
double 
ClpSimplex::valueIncomingDual() const
{
  // Need value of incoming for list of infeasibilities as may be infeasible
  double valueIncoming = (dualOut_/alpha_)*directionOut_;
  if (directionIn_==-1)
    valueIncoming = upperIn_-valueIncoming;
  else
    valueIncoming = lowerIn_-valueIncoming;
  return valueIncoming;
}
