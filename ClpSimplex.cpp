// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.


#undef NDEBUG

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpMessage.hpp"
#include "ClpEventHandler.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpHelperFunctions.hpp"
#include "CoinModel.hpp"
#include <cfloat>

#include <string>
#include <stdio.h>
#include <iostream>
//#############################################################################

ClpSimplex::ClpSimplex () :

  ClpModel(),
  columnPrimalInfeasibility_(0.0),
  rowPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  rowDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e10),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(-COIN_DBL_MAX),
  dualIn_(0.0),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  infeasibilityCost_(1.0e10),
  sumOfRelaxedDualInfeasibilities_(0.0),
  sumOfRelaxedPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  sequenceIn_(-1),
  directionIn_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  lastGoodIteration_(-100),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  numberPrimalInfeasibilities_(100),
  numberRefinements_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  savedSolution_(NULL),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  lastFlaggedIteration_(-999999),
  numberFake_(0),
  numberChanged_(0),
  progressFlag_(0),
  firstFree_(-1),
  numberExtraRows_(0),
  maximumBasic_(0),
  incomingInfeasibility_(1.0),
  allowedInfeasibility_(10.0),
  automaticScale_(0),
  progress_(NULL)
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
  solveType_=1; // say simplex based life form
  
}

// Subproblem constructor
ClpSimplex::ClpSimplex ( const ClpModel * rhs,
		     int numberRows, const int * whichRow,
		     int numberColumns, const int * whichColumn,
		     bool dropNames, bool dropIntegers,bool fixOthers)
  : ClpModel(rhs, numberRows, whichRow,
	     numberColumns,whichColumn,dropNames,dropIntegers),
  columnPrimalInfeasibility_(0.0),
  rowPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  rowDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e10),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(-COIN_DBL_MAX),
  dualIn_(0.0),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  infeasibilityCost_(1.0e10),
  sumOfRelaxedDualInfeasibilities_(0.0),
  sumOfRelaxedPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  sequenceIn_(-1),
  directionIn_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  lastGoodIteration_(-100),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  numberPrimalInfeasibilities_(100),
  numberRefinements_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  savedSolution_(NULL),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  lastFlaggedIteration_(-999999),
  numberFake_(0),
  numberChanged_(0),
  progressFlag_(0),
  firstFree_(-1),
  numberExtraRows_(0),
  maximumBasic_(0),
  incomingInfeasibility_(1.0),
  allowedInfeasibility_(10.0),
  automaticScale_(0),
  progress_(NULL)
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
  // say Steepest pricing
  primalColumnPivot_ = new ClpPrimalColumnSteepest();
  solveType_=1; // say simplex based life form
  if (fixOthers) {
    int numberOtherColumns=rhs->numberColumns();
    int numberOtherRows=rhs->numberRows();
    double * solution = new double [numberOtherColumns];
    CoinZeroN(solution,numberOtherColumns);
    int i;
    for (i=0;i<numberColumns;i++) {
      int iColumn = whichColumn[i];
      if (solution[iColumn])
        fixOthers=false; // duplicates
      solution[iColumn]=1.0;
    }
    if (fixOthers) {
      const double * otherSolution = rhs->primalColumnSolution();
      const double * objective = rhs->objective();
      double offset=0.0;
      for (i=0;i<numberOtherColumns;i++) {
        if (solution[i]) {
          solution[i]=0.0; // in
        } else {
          solution[i] = otherSolution[i];
          offset += objective[i]*otherSolution[i];
        }
      }
      double * rhsModification = new double [numberOtherRows];
      CoinZeroN(rhsModification,numberOtherRows);
      rhs->matrix()->times(solution,rhsModification) ;
      for ( i=0;i<numberRows;i++) {
        int iRow = whichRow[i];
        if (rowLower_[i]>-1.0e20)
          rowLower_[i] -= rhsModification[iRow];
        if (rowUpper_[i]<1.0e20)
          rowUpper_[i] -= rhsModification[iRow];
      }
      delete [] rhsModification;
      setObjectiveOffset(rhs->objectiveOffset()-offset);
      // And set objective value to match
      setObjectiveValue(rhs->objectiveValue());
    }
    delete [] solution;
  }  
}
// Puts solution back small model
void 
ClpSimplex::getbackSolution(const ClpSimplex & smallModel,const int * whichRow, const int * whichColumn)
{
  setSumDualInfeasibilities(smallModel.sumDualInfeasibilities());
  setNumberDualInfeasibilities(smallModel.numberDualInfeasibilities());
  setSumPrimalInfeasibilities(smallModel.sumPrimalInfeasibilities());
  setNumberPrimalInfeasibilities(smallModel.numberPrimalInfeasibilities());
  setNumberIterations(smallModel.numberIterations());
  setProblemStatus(smallModel.status());
  setObjectiveValue(smallModel.objectiveValue());
  const double * solution2 = smallModel.primalColumnSolution();
  int i;
  int numberRows2 = smallModel.numberRows();
  int numberColumns2 = smallModel.numberColumns();
  const double * dj2 = smallModel.dualColumnSolution();
  for ( i=0;i<numberColumns2;i++) {
    int iColumn = whichColumn[i];
    columnActivity_[iColumn]=solution2[i];
    reducedCost_[iColumn]=dj2[i];
    setStatus(iColumn,smallModel.getStatus(i));
  }
  const double * dual2 = smallModel.dualRowSolution();
  memset(dual_,0,numberRows_*sizeof(double));
  for (i=0;i<numberRows2;i++) {
    int iRow=whichRow[i];
    setRowStatus(iRow,smallModel.getRowStatus(i));
    dual_[iRow]=dual2[i];
  }
  CoinZeroN(rowActivity_,numberRows_);
  matrix()->times(columnActivity_,rowActivity_) ;
}

//-----------------------------------------------------------------------------

ClpSimplex::~ClpSimplex ()
{
  gutsOfDelete(0);
  delete nonLinearCost_;
}
//#############################################################################
void ClpSimplex::setLargeValue( double value) 
{
  if (value>0.0&&value<COIN_DBL_MAX)
    largeValue_=value;
}
int
ClpSimplex::gutsOfSolution ( double * givenDuals,
			     const double * givenPrimals,
			     bool valuesPass)
{


  // if values pass, save values of basic variables
  double * save = NULL;
  double oldValue=0.0;
  if (valuesPass) {
    assert(algorithm_>0); // only primal at present
    assert(nonLinearCost_);
    int iRow;
    checkPrimalSolution( rowActivityWork_, columnActivityWork_);
    // get correct bounds on all variables
    nonLinearCost_->checkInfeasibilities(primalTolerance_);
    oldValue = nonLinearCost_->largestInfeasibility();
    save = new double[numberRows_];
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      save[iRow] = solution_[iPivot];
    }
  }
  // do work
  computePrimals(rowActivityWork_, columnActivityWork_);
  // If necessary - override results
  if (givenPrimals) {
    memcpy(columnActivityWork_,givenPrimals,numberColumns_*sizeof(double));
    memset(rowActivityWork_,0,numberRows_*sizeof(double));
    times(-1.0,columnActivityWork_,rowActivityWork_);
  }
  double objectiveModification = 0.0;
  if (algorithm_>0&&nonLinearCost_!=NULL) {
    // primal algorithm
    // get correct bounds on all variables
    // If  4 bit set - Force outgoing variables to exact bound (primal)
    if ((specialOptions_&4)==0)
      nonLinearCost_->checkInfeasibilities(primalTolerance_);
    else
      nonLinearCost_->checkInfeasibilities(0.0);
    objectiveModification += nonLinearCost_->changeInCost();
    if (nonLinearCost_->numberInfeasibilities())
      if (handler_->detail(CLP_SIMPLEX_NONLINEAR,messages_)<100) {
        handler_->message(CLP_SIMPLEX_NONLINEAR,messages_)
          <<nonLinearCost_->changeInCost()
          <<nonLinearCost_->numberInfeasibilities()
          <<CoinMessageEol;
      }
  }
  if (valuesPass) {
#ifdef CLP_DEBUG
    std::cout<<"Largest given infeasibility "<<oldValue
	     <<" now "<<nonLinearCost_->largestInfeasibility()<<std::endl;
#endif
    int numberOut=0;
    // But may be very large rhs etc
    double useError = CoinMin(largestPrimalError_,
                              1.0e5/maximumAbsElement(solution_,numberRows_+numberColumns_));  
    if (oldValue<incomingInfeasibility_
	&&nonLinearCost_->largestInfeasibility()>
	CoinMax(incomingInfeasibility_,allowedInfeasibility_)||
	useError>1.0e-3) {
      //printf("Original largest infeas %g, now %g, primalError %g\n",
      //     oldValue,nonLinearCost_->largestInfeasibility(),
      //     largestPrimalError_);
      // throw out up to 1000 structurals
      int iRow;
      int * sort = new int[numberRows_];
      // first put back solution and store difference
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];
	double difference = fabs(solution_[iPivot]-save[iRow]);
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
      numberOut = CoinMin(1000,numberOut);
      for (iRow=0;iRow<numberOut;iRow++) {
	int iColumn=sort[iRow];
	setColumnStatus(iColumn,superBasic);
	if (fabs(solution_[iColumn])>1.0e10) {
	  if (upper_[iColumn]<0.0) {
	    solution_[iColumn]=upper_[iColumn];
	  } else if (lower_[iColumn]>0.0) {
	    solution_[iColumn]=lower_[iColumn];
	  } else {
	    solution_[iColumn]=0.0;
	  }
	}
      }
      delete [] sort;
    }
    delete [] save;
    if (numberOut)
      return numberOut;
  }

  computeDuals(givenDuals);

  // now check solutions
  //checkPrimalSolution( rowActivityWork_, columnActivityWork_);
  //checkDualSolution();
  checkBothSolutions();
  objectiveValue_ += objectiveModification/(objectiveScale_*rhsScale_);
  if (handler_->logLevel()>3||(largestPrimalError_>1.0e-2||
			       largestDualError_>1.0e-2)) 
    handler_->message(CLP_SIMPLEX_ACCURACY,messages_)
      <<largestPrimalError_
      <<largestDualError_
      <<CoinMessageEol;
  // Switch off false values pass indicator
  if (!valuesPass&&algorithm_>0)
    firstFree_ = -1;
  return 0;
}
void
ClpSimplex::computePrimals ( const double * rowActivities,
				     const double * columnActivities)
{

  //work space
  CoinIndexedVector  * workSpace = rowArray_[0];

  CoinIndexedVector * arrayVector = rowArray_[1];
  arrayVector->clear();
  CoinIndexedVector * previousVector = rowArray_[2];
  previousVector->clear();
  // accumulate non basic stuff 

  int iRow;
  // order is this way for scaling
  if (columnActivities!=columnActivityWork_)
    ClpDisjointCopyN(columnActivities,numberColumns_,columnActivityWork_);
  if (rowActivities!=rowActivityWork_)
    ClpDisjointCopyN(rowActivities,numberRows_,rowActivityWork_);
  double * array = arrayVector->denseVector();
  int * index = arrayVector->getIndices();
  int number=0;
  const double * rhsOffset = matrix_->rhsOffset(this,false,true);
  if (!rhsOffset) {
    // Use whole matrix every time to make it easier for ClpMatrixBase
    // So zero out basic
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      assert (iPivot>=0);
      solution_[iPivot] = 0.0;
    }
    // Extended solution before "update"
    matrix_->primalExpanded(this,0);
    times(-1.0,columnActivityWork_,array);
    for (iRow=0;iRow<numberRows_;iRow++) {
      double value = array[iRow] + rowActivityWork_[iRow];
      if (value) {
	array[iRow]=value;
	index[number++]=iRow;
      } else {
	array[iRow]=0.0;
      }
    }
  } else {
    // we have an effective rhs lying around
    // zero out basic (really just for slacks)
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      solution_[iPivot] = 0.0;
    }
    for (iRow=0;iRow<numberRows_;iRow++) {
      double value = rhsOffset[iRow] + rowActivityWork_[iRow];
      if (value) {
	array[iRow]=value;
	index[number++]=iRow;
      } else {
	array[iRow]=0.0;
      }
    }
  }
  arrayVector->setNumElements(number);
#ifdef CLP_DEBUG
  if (numberIterations_==-3840) {
    int i;
    for (i=0;i<numberRows_+numberColumns_;i++)
      printf("%d status %d\n",i,status_[i]);
    printf("xxxxx1\n");
    for (i=0;i<numberRows_;i++)
      if (array[i])
	printf("%d rhs %g\n",i,array[i]);
    printf("xxxxx2\n");
    for (i=0;i<numberRows_+numberColumns_;i++)
      if (getStatus(i)!=basic)
	printf("%d non basic %g %g %g\n",i,lower_[i],solution_[i],upper_[i]);
    printf("xxxxx3\n");
  }
#endif
  // Ftran adjusted RHS and iterate to improve accuracy
  double lastError=COIN_DBL_MAX;
  int iRefine;
  double * work = workSpace->denseVector();
  CoinIndexedVector * thisVector = arrayVector;
  CoinIndexedVector * lastVector = previousVector;
  factorization_->updateColumn(workSpace,thisVector);
#ifdef CLP_DEBUG
  if (numberIterations_==-3840) {
    int i;
    for (i=0;i<numberRows_;i++)
      if (array[i])
	printf("%d after rhs %g\n",i,array[i]);
    printf("xxxxx4\n");
  }
#endif
  bool goodSolution=true;
  for (iRefine=0;iRefine<numberRefinements_+1;iRefine++) {

    int numberIn = thisVector->getNumElements();
    int * indexIn = thisVector->getIndices();
    double * arrayIn = thisVector->denseVector();
    // put solution in correct place
    if (!rhsOffset) {
      int j;
      for (j=0;j<numberIn;j++) {
	iRow = indexIn[j];
	int iPivot=pivotVariable_[iRow];
	solution_[iPivot] = arrayIn[iRow];
        //assert (fabs(solution_[iPivot])<1.0e100);
      }
    } else {
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];
	solution_[iPivot] = arrayIn[iRow];
        //assert (fabs(solution_[iPivot])<1.0e100);
      }
    }
    // Extended solution after "update"
    matrix_->primalExpanded(this,1);
    // check Ax == b  (for all)
    // signal column generated matrix to just do basic (and gub)
    unsigned int saveOptions = specialOptions();
    setSpecialOptions(16);
    times(-1.0,columnActivityWork_,work);
    setSpecialOptions(saveOptions);
    largestPrimalError_=0.0;
    double multiplier = 131072.0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      double value = work[iRow] + rowActivityWork_[iRow];
      work[iRow] = value*multiplier;
      if (fabs(value)>largestPrimalError_) {
	largestPrimalError_=fabs(value);
      }
    }
    if (largestPrimalError_>=lastError) {
      // restore
      CoinIndexedVector * temp = thisVector;
      thisVector = lastVector;
      lastVector=temp;
      goodSolution=false;
      break;
    }
    if (iRefine<numberRefinements_&&largestPrimalError_>1.0e-10) {
      // try and make better
      // save this
      CoinIndexedVector * temp = thisVector;
      thisVector = lastVector;
      lastVector=temp;
      int * indexOut = thisVector->getIndices();
      int number=0;
      array = thisVector->denseVector();
      thisVector->clear();
      for (iRow=0;iRow<numberRows_;iRow++) {
	double value = work[iRow];
	if (value) {
	  array[iRow]=value;
	  indexOut[number++]=iRow;
	  work[iRow]=0.0;
	}
      }
      thisVector->setNumElements(number);
      lastError=largestPrimalError_;
      factorization_->updateColumn(workSpace,thisVector);
      multiplier = 1.0/multiplier;
      double * previous = lastVector->denseVector();
      number=0;
      for (iRow=0;iRow<numberRows_;iRow++) {
	double value = previous[iRow] + multiplier*array[iRow];
	if (value) {
	  array[iRow]=value;
	  indexOut[number++]=iRow;
	} else {
	  array[iRow]=0.0;
	}
      }
      thisVector->setNumElements(number);
    } else {
      break;
    }
  }

  // solution as accurate as we are going to get
  ClpFillN(work,numberRows_,0.0);
  if (!goodSolution) {
    array = thisVector->denseVector();
    // put solution in correct place
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      solution_[iPivot] = array[iRow];
      //assert (fabs(solution_[iPivot])<1.0e100);
    }
  }
  arrayVector->clear();
  previousVector->clear();
#ifdef CLP_DEBUG
  if (numberIterations_==-3840) {
    exit(77);
  }
#endif
}
// now dual side
void
ClpSimplex::computeDuals(double * givenDjs)
{
  if (objective_->type()==1||!objective_->activated()) {
    // Linear
    //work space
    CoinIndexedVector  * workSpace = rowArray_[0];
    
    CoinIndexedVector * arrayVector = rowArray_[1];
    arrayVector->clear();
    CoinIndexedVector * previousVector = rowArray_[2];
    previousVector->clear();
    int iRow;
#ifdef CLP_DEBUG
    workSpace->checkClear();
#endif
    double * array = arrayVector->denseVector();
    int * index = arrayVector->getIndices();
    int number=0;
    if (!givenDjs) {
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];
	double value = cost_[iPivot];
	if (value) {
	  array[iRow]=value;
	  index[number++]=iRow;
	}
      }
    } else {
      // dual values pass - djs may not be zero
      for (iRow=0;iRow<numberRows_;iRow++) {
	int iPivot=pivotVariable_[iRow];
	// make sure zero if done
	if (!pivoted(iPivot))
	  givenDjs[iPivot]=0.0;
	double value =cost_[iPivot]-givenDjs[iPivot];
	if (value) {
	  array[iRow]=value;
	  index[number++]=iRow;
	}
      }
    }
    arrayVector->setNumElements(number);
    // Extended duals before "updateTranspose"
    matrix_->dualExpanded(this,arrayVector,givenDjs,0);
    
    // Btran basic costs and get as accurate as possible
    double lastError=COIN_DBL_MAX;
    int iRefine;
    double * work = workSpace->denseVector();
    CoinIndexedVector * thisVector = arrayVector;
    CoinIndexedVector * lastVector = previousVector;
    factorization_->updateColumnTranspose(workSpace,thisVector);
    
    for (iRefine=0;iRefine<numberRefinements_+1;iRefine++) {
      // check basic reduced costs zero
      largestDualError_=0.0;
      if (!numberExtraRows_) {
	// Just basic
	int * index2 = workSpace->getIndices();
	// use reduced costs for slacks as work array
	double * work2 = reducedCostWork_+numberColumns_;
	int numberStructurals=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  int iPivot=pivotVariable_[iRow];
	  if (iPivot<numberColumns_) 
	    index2[numberStructurals++]=iPivot;
	}
	matrix_->listTransposeTimes(this,array,index2,numberStructurals,work2);
	numberStructurals=0;
	if (!givenDjs) {
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    int iPivot=pivotVariable_[iRow];
	    double value;
	    if (iPivot>=numberColumns_) {
	      // slack
	      value = rowObjectiveWork_[iPivot-numberColumns_]
		+ array[iPivot-numberColumns_];
	    } else {
	      // column
	      value = objectiveWork_[iPivot]-work2[numberStructurals++];
	    }
	    work[iRow]=value;
	    if (fabs(value)>largestDualError_) {
	      largestDualError_=fabs(value);
	    }
	  }
	} else {
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    int iPivot=pivotVariable_[iRow];
	    if (iPivot>=numberColumns_) {
	      // slack
	      work[iRow] = rowObjectiveWork_[iPivot-numberColumns_]
		+ array[iPivot-numberColumns_]-givenDjs[iPivot];
	    } else {
	      // column
	      work[iRow] = objectiveWork_[iPivot]-work2[numberStructurals++]
		- givenDjs[iPivot];
	    }
	    if (fabs(work[iRow])>largestDualError_) {
	      largestDualError_=fabs(work[iRow]);
	      //assert (largestDualError_<1.0e-7);
	      //if (largestDualError_>1.0e-7)
	      //printf("large dual error %g\n",largestDualError_);
	    }
	  }
	}
      } else {
	// extra rows - be more careful
#if 1
	// would be faster to do just for basic but this reduces code
	ClpDisjointCopyN(objectiveWork_,numberColumns_,reducedCostWork_);
	transposeTimes(-1.0,array,reducedCostWork_);
#else
	// Just basic
	int * index2 = workSpace->getIndices();
	int numberStructurals=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  int iPivot=pivotVariable_[iRow];
	  if (iPivot<numberColumns_) 
	    index2[numberStructurals++]=iPivot;
	}
	matrix_->listTransposeTimes(this,array,index2,numberStructurals,work);
	for (iRow=0;iRow<numberStructurals;iRow++) {
	  int iPivot=index2[iRow];
	  reducedCostWork_[iPivot]=objectiveWork_[iPivot]-work[iRow];
	}
#endif
	// update by duals on sets
	matrix_->dualExpanded(this,NULL,NULL,1);
	if (!givenDjs) {
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    int iPivot=pivotVariable_[iRow];
	    double value;
	    if (iPivot>=numberColumns_) {
	      // slack
	      value = rowObjectiveWork_[iPivot-numberColumns_]
		+ array[iPivot-numberColumns_];
	    } else {
	      // column
	      value = reducedCostWork_[iPivot];
	    }
	    work[iRow]=value;
	    if (fabs(value)>largestDualError_) {
	      largestDualError_=fabs(value);
	    }
	  }
	} else {
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    int iPivot=pivotVariable_[iRow];
	    if (iPivot>=numberColumns_) {
	      // slack
	      work[iRow] = rowObjectiveWork_[iPivot-numberColumns_]
		+ array[iPivot-numberColumns_]-givenDjs[iPivot];
	    } else {
	      // column
	      work[iRow] = reducedCostWork_[iPivot]- givenDjs[iPivot];
	    }
	    if (fabs(work[iRow])>largestDualError_) {
	      largestDualError_=fabs(work[iRow]);
	      //assert (largestDualError_<1.0e-7);
	      //if (largestDualError_>1.0e-7)
	      //printf("large dual error %g\n",largestDualError_);
	    }
	  }
	}
      }
      if (largestDualError_>=lastError) {
	// restore
	CoinIndexedVector * temp = thisVector;
	thisVector = lastVector;
	lastVector=temp;
	break;
      }
      if (iRefine<numberRefinements_&&largestDualError_>1.0e-10
	  &&!givenDjs) {
	// try and make better
	// save this
	CoinIndexedVector * temp = thisVector;
	thisVector = lastVector;
	lastVector=temp;
	int * indexOut = thisVector->getIndices();
	int number=0;
	array = thisVector->denseVector();
	thisVector->clear();
	double multiplier = 131072.0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  double value = multiplier*work[iRow];
	  if (value) {
	    array[iRow]=value;
	    indexOut[number++]=iRow;
	    work[iRow]=0.0;
	  }
	  work[iRow]=0.0;
	}
	thisVector->setNumElements(number);
	lastError=largestDualError_;
	factorization_->updateColumnTranspose(workSpace,thisVector);
	multiplier = 1.0/multiplier;
	double * previous = lastVector->denseVector();
	number=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  double value = previous[iRow] + multiplier*array[iRow];
	  if (value) {
	    array[iRow]=value;
	    indexOut[number++]=iRow;
	  } else {
	    array[iRow]=0.0;
	  }
	}
	thisVector->setNumElements(number);
      } else {
	break;
      }
    }
    // now look at dual solution
    array = thisVector->denseVector();
    for (iRow=0;iRow<numberRows_;iRow++) {
      // slack
      double value = array[iRow];
      dual_[iRow]=value;
      value += rowObjectiveWork_[iRow];
      rowReducedCost_[iRow]=value;
    }
    ClpDisjointCopyN(objectiveWork_,numberColumns_,reducedCostWork_);
    // can use work if problem scaled (for better cache)
    if (numberRows_>10000)
      matrix_->transposeTimes(-1.0,dual_,reducedCostWork_,
			      rowScale_,columnScale_,work);
    else
      matrix_->transposeTimes(-1.0,dual_,reducedCostWork_,
			      rowScale_,columnScale_,NULL);
    ClpFillN(work,numberRows_,0.0);
    // Extended duals and check dual infeasibility
    if (!matrix_->skipDualCheck()||algorithm_<0||problemStatus_!=-2) 
      matrix_->dualExpanded(this,NULL,NULL,2);
    // If necessary - override results
    if (givenDjs) {
      // restore accurate duals
      memcpy(givenDjs,dj_,(numberRows_+numberColumns_)*sizeof(double));
    }
    arrayVector->clear();
    previousVector->clear();
  } else {
    // Nonlinear
    objective_->reducedGradient(this,dj_,false);
    // get dual_ by moving from reduced costs for slacks
    memcpy(dual_,dj_+numberColumns_,numberRows_*sizeof(double));
  }
}
/* Given an existing factorization computes and checks 
   primal and dual solutions.  Uses input arrays for variables at
   bounds.  Returns feasibility states */
int ClpSimplex::getSolution ( const double * rowActivities,
			       const double * columnActivities)
{
  if (!factorization_->status()) {
    // put in standard form
    createRim(7+8+16+32,false,-1);
    if (pivotVariable_[0]<0)
      internalFactorize(0);
    // do work
    gutsOfSolution ( NULL,NULL);
    // release extra memory
    deleteRim(0);
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
  ClpDisjointCopyN ( rowActivityWork_, numberRows_ , rowActivities);
  ClpDisjointCopyN ( columnActivityWork_, numberColumns_ , columnActivities);
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
  createRim(7+8+16+32,false);
  // do work
  int status = internalFactorize(-1);
  // release extra memory
  deleteRim(0);

  return status;
}
// Clean up status
void 
ClpSimplex::cleanStatus()
{
  int iRow,iColumn;
  int numberBasic=0;
  // make row activities correct
  memset(rowActivityWork_,0,numberRows_*sizeof(double));
  times(1.0,columnActivityWork_,rowActivityWork_);
  if (!status_)
    createStatus();
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (getRowStatus(iRow)==basic) 
      numberBasic++;
    else {
      setRowStatus(iRow,superBasic);
      // but put to bound if close
      if (fabs(rowActivityWork_[iRow]-rowLowerWork_[iRow])
	  <=primalTolerance_) {
	rowActivityWork_[iRow]=rowLowerWork_[iRow];
	setRowStatus(iRow,atLowerBound);
      } else if (fabs(rowActivityWork_[iRow]-rowUpperWork_[iRow])
		 <=primalTolerance_) {
	rowActivityWork_[iRow]=rowUpperWork_[iRow];
	setRowStatus(iRow,atUpperBound);
      }
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (getColumnStatus(iColumn)==basic) {
      if (numberBasic==numberRows_) {
	// take out of basis
	setColumnStatus(iColumn,superBasic);
	// but put to bound if close
	if (fabs(columnActivityWork_[iColumn]-columnLowerWork_[iColumn])
	    <=primalTolerance_) {
	  columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	  setColumnStatus(iColumn,atLowerBound);
	} else if (fabs(columnActivityWork_[iColumn]
			-columnUpperWork_[iColumn])
		   <=primalTolerance_) {
	  columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	  setColumnStatus(iColumn,atUpperBound);
	}
      } else 
	numberBasic++;
    } else {
      setColumnStatus(iColumn,superBasic);
      // but put to bound if close
      if (fabs(columnActivityWork_[iColumn]-columnLowerWork_[iColumn])
	  <=primalTolerance_) {
	columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	setColumnStatus(iColumn,atLowerBound);
      } else if (fabs(columnActivityWork_[iColumn]
		      -columnUpperWork_[iColumn])
		 <=primalTolerance_) {
	columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	setColumnStatus(iColumn,atUpperBound);
      }
    }
  }
}

/* Factorizes using current basis.  
   solveType - 1 iterating, 0 initial, -1 external 
   - 2 then iterating but can throw out of basis
   If 10 added then in primal values pass
   Return codes are as from ClpFactorization unless initial factorization
   when total number of singularities is returned.
   Special case is numberRows_+1 -> all slack basis.
*/
int ClpSimplex::internalFactorize ( int solveType)
{
  int iRow,iColumn;
  int totalSlacks=numberRows_;
  if (!status_)
    createStatus();

  bool valuesPass=false;
  if (solveType>=10) {
    valuesPass=true;
    solveType -= 10;
  }
#ifdef CLP_DEBUG
  if (solveType>0) {
    int numberFreeIn=0,numberFreeOut=0;
    double biggestDj=0.0;
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      switch(getColumnStatus(iColumn)) {

      case basic:
	if (columnLower_[iColumn]<-largeValue_
	    &&columnUpper_[iColumn]>largeValue_) 
	  numberFreeIn++;
	break;
      default:
	if (columnLower_[iColumn]<-largeValue_
	    &&columnUpper_[iColumn]>largeValue_) {
	  numberFreeOut++;
	  biggestDj = CoinMax(fabs(dj_[iColumn]),biggestDj);
	}
	break;
      }
    }
    if (numberFreeIn+numberFreeOut)
      printf("%d in basis, %d out - largest dj %g\n",
	     numberFreeIn,numberFreeOut,biggestDj);
  }
#endif
  if (solveType<=0) {
    // Make sure everything is clean
    for (iRow=0;iRow<numberRows_;iRow++) {
      if(getRowStatus(iRow)==isFixed) {
	// double check fixed
	if (rowUpperWork_[iRow]>rowLowerWork_[iRow])
	  setRowStatus(iRow,atLowerBound);
      } else if (getRowStatus(iRow)==isFree) {
	// may not be free after all
	if (rowLowerWork_[iRow]>-largeValue_||rowUpperWork_[iRow]<largeValue_)
	  setRowStatus(iRow,superBasic);
      }
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      if(getColumnStatus(iColumn)==isFixed) {
	// double check fixed
	if (columnUpperWork_[iColumn]>columnLowerWork_[iColumn])
	  setColumnStatus(iColumn,atLowerBound);
      } else if (getColumnStatus(iColumn)==isFree) {
	// may not be free after all
	if (columnLowerWork_[iColumn]>-largeValue_||columnUpperWork_[iColumn]<largeValue_)
	  setColumnStatus(iColumn,superBasic);
      }
    }
    if (!valuesPass) {
      // not values pass so set to bounds
      bool allSlack=true;
      if (status_) {
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (getRowStatus(iRow)!=basic) {
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
	    
	  case basic:
	    numberBasic++;
	    break;
	  case atUpperBound:
	    rowActivityWork_[iRow]=rowUpperWork_[iRow];
	    if (rowActivityWork_[iRow]>largeValue_) {
	      if (rowLowerWork_[iRow]>-largeValue_) {
		rowActivityWork_[iRow]=rowLowerWork_[iRow];
		setRowStatus(iRow,atLowerBound);
	      } else {
		// say free
		setRowStatus(iRow,isFree);
		rowActivityWork_[iRow]=0.0;
	      }
	    }
	    break;
	  case ClpSimplex::isFixed:
	  case atLowerBound:
	    rowActivityWork_[iRow]=rowLowerWork_[iRow];
	    if (rowActivityWork_[iRow]<-largeValue_) {
	      if (rowUpperWork_[iRow]<largeValue_) {
		rowActivityWork_[iRow]=rowUpperWork_[iRow];
		setRowStatus(iRow,atUpperBound);
	      } else {
		// say free
		setRowStatus(iRow,isFree);
		rowActivityWork_[iRow]=0.0;
	      }
	    }
	    break;
	  case isFree:
	      break;
	    // not really free - fall through to superbasic
	  case superBasic:
	    if (rowUpperWork_[iRow]>largeValue_) {
	      if (rowLowerWork_[iRow]>-largeValue_) {
		rowActivityWork_[iRow]=rowLowerWork_[iRow];
		setRowStatus(iRow,atLowerBound);
	      } else {
		// say free
		setRowStatus(iRow,isFree);
		rowActivityWork_[iRow]=0.0;
	      }
	    } else {
	      if (rowLowerWork_[iRow]>-largeValue_) {
		// set to nearest
		if (fabs(rowActivityWork_[iRow]-rowLowerWork_[iRow])
		    <fabs(rowActivityWork_[iRow]-rowLowerWork_[iRow])) {
		  rowActivityWork_[iRow]=rowLowerWork_[iRow];
		  setRowStatus(iRow,atLowerBound);
		} else {
		  rowActivityWork_[iRow]=rowUpperWork_[iRow];
		  setRowStatus(iRow,atUpperBound);
		}
	      } else {
		rowActivityWork_[iRow]=rowUpperWork_[iRow];
		setRowStatus(iRow,atUpperBound);
	      }
	    }
	    break;
	  }
	}
	totalSlacks=numberBasic;

	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  switch(getColumnStatus(iColumn)) {
	    
	  case basic:
	    if (numberBasic==maximumBasic_) {
	      // take out of basis
	      if (columnLowerWork_[iColumn]>-largeValue_) {
		if (columnActivityWork_[iColumn]-columnLowerWork_[iColumn]<
		    columnUpperWork_[iColumn]-columnActivityWork_[iColumn]) {
		  columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
		  setColumnStatus(iColumn,atLowerBound);
		} else {
		  columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		  setColumnStatus(iColumn,atUpperBound);
		}
	      } else if (columnUpperWork_[iColumn]<largeValue_) {
		columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		setColumnStatus(iColumn,atUpperBound);
	      } else {
		columnActivityWork_[iColumn]=0.0;
		setColumnStatus(iColumn,isFree);
	      }
	    } else {
	      numberBasic++;
	    }
	    break;
	  case atUpperBound:
	    columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	    if (columnActivityWork_[iColumn]>largeValue_) {
	      if (columnLowerWork_[iColumn]<-largeValue_) {
		columnActivityWork_[iColumn]=0.0;
		setColumnStatus(iColumn,isFree);
	      } else {
		columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
		setColumnStatus(iColumn,atLowerBound);
	      }
	    }
	    break;
	  case isFixed:
	  case atLowerBound:
	    columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	    if (columnActivityWork_[iColumn]<-largeValue_) {
	      if (columnUpperWork_[iColumn]>largeValue_) {
		columnActivityWork_[iColumn]=0.0;
		setColumnStatus(iColumn,isFree);
	      } else {
		columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		setColumnStatus(iColumn,atUpperBound);
	      }
	    }
	    break;
	  case isFree:
	      break;
	    // not really free - fall through to superbasic
	  case superBasic:
	    if (columnUpperWork_[iColumn]>largeValue_) {
	      if (columnLowerWork_[iColumn]>-largeValue_) {
		columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
		setColumnStatus(iColumn,atLowerBound);
	      } else {
		// say free
		setColumnStatus(iColumn,isFree);
		columnActivityWork_[iColumn]=0.0;
	      }
	    } else {
	      if (columnLowerWork_[iColumn]>-largeValue_) {
		// set to nearest
		if (fabs(columnActivityWork_[iColumn]-columnLowerWork_[iColumn])
		    <fabs(columnActivityWork_[iColumn]-columnLowerWork_[iColumn])) {
		  columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
		  setColumnStatus(iColumn,atLowerBound);
		} else {
		  columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		  setColumnStatus(iColumn,atUpperBound);
		}
	      } else {
		columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
		setColumnStatus(iColumn,atUpperBound);
	      }
	    }
	    break;
	  }
	}
      } else {
	// all slack basis
	int numberBasic=0;
	if (!status_) {
	  createStatus();
	}
	for (iRow=0;iRow<numberRows_;iRow++) {
	  double lower=rowLowerWork_[iRow];
	  double upper=rowUpperWork_[iRow];
	  if (lower>-largeValue_||upper<largeValue_) {
	    if (fabs(lower)<=fabs(upper)) {
	      rowActivityWork_[iRow]=lower;
	    } else {
	      rowActivityWork_[iRow]=upper;
	    }
	  } else {
	    rowActivityWork_[iRow]=0.0;
	  }
	  setRowStatus(iRow,basic);
	  numberBasic++;
	}
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  double lower=columnLowerWork_[iColumn];
	  double upper=columnUpperWork_[iColumn];
	  double big_bound = largeValue_;
	  if (lower>-big_bound||upper<big_bound) {
	    if ((getColumnStatus(iColumn)==atLowerBound&&
		 columnActivityWork_[iColumn]==lower)||
		(getColumnStatus(iColumn)==atUpperBound&&
		 columnActivityWork_[iColumn]==upper)) {
	      // status looks plausible
	    } else {
	      // set to sensible
	      if (fabs(lower)<=fabs(upper)) {
		setColumnStatus(iColumn,atLowerBound);
		columnActivityWork_[iColumn]=lower;
	      } else {
		setColumnStatus(iColumn,atUpperBound);
		columnActivityWork_[iColumn]=upper;
	      }
	    }
	  } else {
	    setColumnStatus(iColumn,isFree);
	    columnActivityWork_[iColumn]=0.0;
	  }
	}
      }
    } else {
      // values pass has less coding
      // make row activities correct and clean basis a bit
      cleanStatus();
      if (status_) {
	int numberBasic=0;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (getRowStatus(iRow)==basic) 
	    numberBasic++;
	}
	totalSlacks=numberBasic;
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  if (getColumnStatus(iColumn)==basic) 
	    numberBasic++;
	}
      } else {
	// all slack basis
	int numberBasic=0;
	if (!status_) {
	  createStatus();
	}
	for (iRow=0;iRow<numberRows_;iRow++) {
	  setRowStatus(iRow,basic);
	  numberBasic++;
	}
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  setColumnStatus(iColumn,superBasic);
	  // but put to bound if close
	  if (fabs(columnActivityWork_[iColumn]-columnLowerWork_[iColumn])
	      <=primalTolerance_) {
	    columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	    setColumnStatus(iColumn,atLowerBound);
	  } else if (fabs(columnActivityWork_[iColumn]
			-columnUpperWork_[iColumn])
		     <=primalTolerance_) {
	    columnActivityWork_[iColumn]=columnUpperWork_[iColumn];
	    setColumnStatus(iColumn,atUpperBound);
	  }
	}
      }
    }
    numberRefinements_=1;
    // set fixed if they are
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (getRowStatus(iRow)!=basic ) {
	if (rowLowerWork_[iRow]==rowUpperWork_[iRow]) {
	  rowActivityWork_[iRow]=rowLowerWork_[iRow];
	  setRowStatus(iRow,isFixed);
	}
      }
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (getColumnStatus(iColumn)!=basic ) {
	if (columnLowerWork_[iColumn]==columnUpperWork_[iColumn]) {
	  columnActivityWork_[iColumn]=columnLowerWork_[iColumn];
	  setColumnStatus(iColumn,isFixed);
	}
      }
    }
  }
  //for (iRow=0;iRow<numberRows_+numberColumns_;iRow++) {
  //if (fabs(solution_[iRow])>1.0e10) {
  //  printf("large %g at %d - status %d\n",
  //         solution_[iRow],iRow,status_[iRow]);
  //}
  //}
  if (0)  {
    static int k=0;
    printf("start basis\n");
    int i;
    for (i=0;i<numberRows_;i++)
      printf ("xx %d %d\n",i,pivotVariable_[i]);
    for (i=0;i<numberRows_+numberColumns_;i++)
      if (getColumnStatus(i)==basic)
	printf ("yy %d basic\n",i);
    if (k>20)
      exit(0);
    k++;
  }
  int status = factorization_->factorize(this, solveType,valuesPass);
  if (status) {
    handler_->message(CLP_SIMPLEX_BADFACTOR,messages_)
      <<status
      <<CoinMessageEol;
    return -1;
  } else if (!solveType) {
    // Initial basis - return number of singularities
    int numberSlacks=0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (getRowStatus(iRow) == basic)
	numberSlacks++;
    }
    status= CoinMax(numberSlacks-totalSlacks,0);
    // special case if all slack
    if (numberSlacks==numberRows_) {
      status=numberRows_+1;
    }
  }

  // sparse methods
  //if (factorization_->sparseThreshold()) {
    // get default value
    factorization_->sparseThreshold(0);
    factorization_->goSparse();
    //}
  
  return status;
}
/* 
   This does basis housekeeping and does values for in/out variables.
   Can also decide to re-factorize
*/
int 
ClpSimplex::housekeeping(double objectiveChange)
{
  // save value of incoming and outgoing
  double oldIn = solution_[sequenceIn_];
  double oldOut = solution_[sequenceOut_];
  numberIterations_++;
  changeMade_++; // something has happened
  // incoming variable
  if (handler_->detail(CLP_SIMPLEX_HOUSE1,messages_)<100) {
    handler_->message(CLP_SIMPLEX_HOUSE1,messages_)
      <<directionOut_
      <<directionIn_<<theta_
      <<dualOut_<<dualIn_<<alpha_
      <<CoinMessageEol;
    if (getStatus(sequenceIn_)==isFree) {
      handler_->message(CLP_SIMPLEX_FREEIN,messages_)
        <<sequenceIn_
        <<CoinMessageEol;
    }
  }
  // change of incoming
  char rowcol[]={'R','C'};
  if (pivotRow_>=0)
    pivotVariable_[pivotRow_]=sequenceIn();
  if (upper_[sequenceIn_]>1.0e20&&lower_[sequenceIn_]<-1.0e20)
    progressFlag_ |= 2; // making real progress
  solution_[sequenceIn_]=valueIn_;
  if (upper_[sequenceOut_]-lower_[sequenceOut_]<1.0e-12)
    progressFlag_ |= 1; // making real progress
  if (sequenceIn_!=sequenceOut_) {
    //assert( getStatus(sequenceOut_)== basic);
    setStatus(sequenceIn_,basic);
    if (upper_[sequenceOut_]-lower_[sequenceOut_]>0) {
      // As Nonlinear costs may have moved bounds (to more feasible)
      // Redo using value
      if (fabs(valueOut_-lower_[sequenceOut_])<fabs(valueOut_-upper_[sequenceOut_])) {
	// going to lower
	setStatus(sequenceOut_,atLowerBound);
	oldOut = lower_[sequenceOut_];
      } else {
	// going to upper
	setStatus(sequenceOut_,atUpperBound);
	oldOut = upper_[sequenceOut_];
      }
    } else {
      // fixed
      setStatus(sequenceOut_,isFixed);
    }
    solution_[sequenceOut_]=valueOut_;
  } else {
    //if (objective_->type()<2)
    assert (fabs(theta_)>1.0e-13);
    // flip from bound to bound
    // As Nonlinear costs may have moved bounds (to more feasible)
    // Redo using value
    if (fabs(valueIn_-lower_[sequenceIn_])<fabs(valueIn_-upper_[sequenceIn_])) {
      // as if from upper bound
      setStatus(sequenceIn_, atLowerBound);
    } else {
      // as if from lower bound
      setStatus(sequenceIn_, atUpperBound);
    }
  }
  
  // Update hidden stuff e.g. effective RHS and gub
  matrix_->updatePivot(this,oldIn,oldOut);
  objectiveValue_ += objectiveChange/(objectiveScale_*rhsScale_);
  if (handler_->detail(CLP_SIMPLEX_HOUSE2,messages_)<100) {
    handler_->message(CLP_SIMPLEX_HOUSE2,messages_)
      <<numberIterations_<<objectiveValue()
      <<rowcol[isColumn(sequenceIn_)]<<sequenceWithin(sequenceIn_)
      <<rowcol[isColumn(sequenceOut_)]<<sequenceWithin(sequenceOut_);
    handler_->printing(algorithm_<0)<<dualOut_<<theta_;
    handler_->printing(algorithm_>0)<<dualIn_<<theta_;
    handler_->message()<<CoinMessageEol;
  }
  if (hitMaximumIterations())
    return 2;
#if 1
  //if (numberIterations_>14000)
  //handler_->setLogLevel(63);
  //if (numberIterations_>24000)
  //exit(77);
  // check for small cycles
  int in = sequenceIn_;
  int out = sequenceOut_;
  matrix_->correctSequence(in,out);
  int cycle=progress_->cycle(in,out,
			    directionIn_,directionOut_);
  if (cycle>0&&objective_->type()<2) {
    //if (cycle>0) {
    if (handler_->logLevel()>=63)
      printf("Cycle of %d\n",cycle);
    // reset
    progress_->startCheck();
    double random = CoinDrand48();
    int extra = (int) (9.999*random);
    int off[]={1,1,1,1,2,2,2,3,3,4};
    if (factorization_->pivots()>cycle) {
      forceFactorization_=CoinMax(1,cycle-off[extra]);
    } else {
      // need to reject something
      int iSequence;
      if (algorithm_<0)
	iSequence = sequenceIn_;
      else
	iSequence = sequenceOut_;
      char x = isColumn(iSequence) ? 'C' :'R';
      if (handler_->logLevel()>=63)
	handler_->message(CLP_SIMPLEX_FLAG,messages_)
	  <<x<<sequenceWithin(iSequence)
	  <<CoinMessageEol;
      setFlagged(iSequence);
      //printf("flagging %d\n",iSequence);
    }
    return 1;
  }
#endif
  // only time to re-factorize if one before real time
  // this is so user won't be surprised that maximumPivots has exact meaning
  CoinBigIndex lengthL = factorization_->numberElementsL();
  CoinBigIndex lengthU = factorization_->numberElementsU();
  CoinBigIndex lengthR = factorization_->numberElementsR();
  int numberPivots=factorization_->pivots();
  int maximumPivots = factorization_->maximumPivots();
  int numberDense = factorization_->numberDense();
  if (numberPivots==maximumPivots||
      maximumPivots<2) {
    // If dense then increase
    if (maximumPivots>100&&numberDense>1.5*maximumPivots) {
      factorization_->maximumPivots(numberDense);
      dualRowPivot_->maximumPivotsChanged();
      primalColumnPivot_->maximumPivotsChanged();
      // and redo arrays
      for (int iRow=0;iRow<4;iRow++) {
	int length =rowArray_[iRow]->capacity()+numberDense-maximumPivots;
	rowArray_[iRow]->reserve(length);
      }
    }
    return 1;
  } else if (numberPivots*3>maximumPivots*2&&
             lengthR*3>(lengthL+lengthU)*2+1000&&!numberDense) {
    //printf("ret after %d pivots\n",factorization_->pivots());
    return 1;
  } else {
    if (forceFactorization_>0&&
	factorization_->pivots()==forceFactorization_) {
      // relax
      forceFactorization_ = (3+5*forceFactorization_)/4;
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
ClpSimplex::ClpSimplex(const ClpSimplex &rhs,int scalingMode) :
  ClpModel(rhs,scalingMode),
  columnPrimalInfeasibility_(0.0),
  rowPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  rowDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e10),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(-COIN_DBL_MAX),
  dualIn_(0.0),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  infeasibilityCost_(1.0e10),
  sumOfRelaxedDualInfeasibilities_(0.0),
  sumOfRelaxedPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  sequenceIn_(-1),
  directionIn_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  lastGoodIteration_(-100),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  numberPrimalInfeasibilities_(100),
  numberRefinements_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  savedSolution_(NULL),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  lastFlaggedIteration_(-999999),
  numberFake_(0),
  numberChanged_(0),
  progressFlag_(0),
  firstFree_(-1),
  numberExtraRows_(0),
  maximumBasic_(0),
  incomingInfeasibility_(1.0),
  allowedInfeasibility_(10.0),
  automaticScale_(0),
  progress_(NULL)
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
  specialOptions_ =0;
  delete nonLinearCost_;
  nonLinearCost_ = NULL;
  gutsOfCopy(rhs);
  solveType_=1; // say simplex based life form
}
// Copy constructor from model
ClpSimplex::ClpSimplex(const ClpModel &rhs, int scalingMode) :
  ClpModel(rhs,scalingMode),
  columnPrimalInfeasibility_(0.0),
  rowPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  rowDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  remainingDualInfeasibility_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  dualBound_(1.0e10),
  alpha_(0.0),
  theta_(0.0),
  lowerIn_(0.0),
  valueIn_(0.0),
  upperIn_(-COIN_DBL_MAX),
  dualIn_(0.0),
  lowerOut_(-1),
  valueOut_(-1),
  upperOut_(-1),
  dualOut_(-1),
  dualTolerance_(0.0),
  primalTolerance_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  infeasibilityCost_(1.0e10),
  sumOfRelaxedDualInfeasibilities_(0.0),
  sumOfRelaxedPrimalInfeasibilities_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rowObjectiveWork_(NULL),
  objectiveWork_(NULL),
  sequenceIn_(-1),
  directionIn_(-1),
  sequenceOut_(-1),
  directionOut_(-1),
  pivotRow_(-1),
  lastGoodIteration_(-100),
  dj_(NULL),
  rowReducedCost_(NULL),
  reducedCostWork_(NULL),
  solution_(NULL),
  rowActivityWork_(NULL),
  columnActivityWork_(NULL),
  numberDualInfeasibilities_(0),
  numberDualInfeasibilitiesWithoutFree_(0),
  numberPrimalInfeasibilities_(100),
  numberRefinements_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  savedSolution_(NULL),
  numberTimesOptimal_(0),
  changeMade_(1),
  algorithm_(0),
  forceFactorization_(-1),
  perturbation_(100),
  nonLinearCost_(NULL),
  specialOptions_(0),
  lastBadIteration_(-999999),
  lastFlaggedIteration_(-999999),
  numberFake_(0),
  numberChanged_(0),
  progressFlag_(0),
  firstFree_(-1),
  numberExtraRows_(0),
  maximumBasic_(0),
  incomingInfeasibility_(1.0),
  allowedInfeasibility_(10.0),
  automaticScale_(0),
  progress_(NULL)
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
  // say Steepest pricing
  primalColumnPivot_ = new ClpPrimalColumnSteepest();
  solveType_=1; // say simplex based life form
  
}
// Assignment operator. This copies the data
ClpSimplex & 
ClpSimplex::operator=(const ClpSimplex & rhs)
{
  if (this != &rhs) {
    gutsOfDelete(0);
    specialOptions_=0;
    delete nonLinearCost_;
    nonLinearCost_ = NULL;
    ClpModel::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
ClpSimplex::gutsOfCopy(const ClpSimplex & rhs)
{
  assert (numberRows_==rhs.numberRows_);
  assert (numberColumns_==rhs.numberColumns_);
  numberExtraRows_ = rhs.numberExtraRows_;
  maximumBasic_ = rhs.maximumBasic_;
  int numberRows2 = numberRows_+numberExtraRows_;
  if ((whatsChanged_&1)!=0) {
    lower_ = ClpCopyOfArray(rhs.lower_,numberColumns_+numberRows2);
    rowLowerWork_ = lower_+numberColumns_;
    columnLowerWork_ = lower_;
    upper_ = ClpCopyOfArray(rhs.upper_,numberColumns_+numberRows2);
    rowUpperWork_ = upper_+numberColumns_;
    columnUpperWork_ = upper_;
    //cost_ = ClpCopyOfArray(rhs.cost_,2*(numberColumns_+numberRows_));
    cost_ = ClpCopyOfArray(rhs.cost_,(numberColumns_+numberRows2));
    objectiveWork_ = cost_;
    rowObjectiveWork_ = cost_+numberColumns_;
    dj_ = ClpCopyOfArray(rhs.dj_,numberRows2+numberColumns_);
    if (dj_) {
      reducedCostWork_ = dj_;
      rowReducedCost_ = dj_+numberColumns_;
    }
    solution_ = ClpCopyOfArray(rhs.solution_,numberRows2+numberColumns_);
    if (solution_) {
      columnActivityWork_ = solution_;
      rowActivityWork_ = solution_+numberColumns_;
    }
    if (rhs.pivotVariable_) {
      pivotVariable_ = new int[numberRows2];
      CoinMemcpyN ( rhs.pivotVariable_, numberRows2 , pivotVariable_);
    } else {
      pivotVariable_=NULL;
    }
    savedSolution_ = ClpCopyOfArray(rhs.savedSolution_,numberColumns_+numberRows_);
    int i;
    for (i=0;i<6;i++) {
      rowArray_[i]=NULL;
      if (rhs.rowArray_[i]) 
        rowArray_[i] = new CoinIndexedVector(*rhs.rowArray_[i]);
      columnArray_[i]=NULL;
      if (rhs.columnArray_[i]) 
        columnArray_[i] = new CoinIndexedVector(*rhs.columnArray_[i]);
    }
    if (rhs.saveStatus_) {
      saveStatus_ = ClpCopyOfArray( rhs.saveStatus_,numberColumns_+numberRows_);
    }
  } else {
    lower_ = NULL;
    rowLowerWork_ = NULL;
    columnLowerWork_ = NULL; 
    upper_ = NULL; 
    rowUpperWork_ = NULL; 
    columnUpperWork_ = NULL;
    cost_ = NULL; 
    objectiveWork_ = NULL; 
    rowObjectiveWork_ = NULL; 
    dj_ = NULL; 
    reducedCostWork_ = NULL;
    rowReducedCost_ = NULL;
    solution_ = NULL; 
    columnActivityWork_ = NULL; 
    rowActivityWork_ = NULL; 
    pivotVariable_=NULL;
    savedSolution_ = NULL;
    int i;
    for (i=0;i<6;i++) {
      rowArray_[i]=NULL;
      columnArray_[i]=NULL;
    }
    saveStatus_ = NULL;
  }
  if (rhs.factorization_) {
    factorization_ = new ClpFactorization(*rhs.factorization_);
  } else {
    factorization_=NULL;
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
  lastGoodIteration_ = rhs.lastGoodIteration_;
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
  numberTimesOptimal_ = rhs.numberTimesOptimal_;
  changeMade_ = rhs.changeMade_;
  algorithm_ = rhs.algorithm_;
  forceFactorization_ = rhs.forceFactorization_;
  perturbation_ = rhs.perturbation_;
  infeasibilityCost_ = rhs.infeasibilityCost_;
  specialOptions_ = rhs.specialOptions_;
  lastBadIteration_ = rhs.lastBadIteration_;
  lastFlaggedIteration_ = rhs.lastFlaggedIteration_;
  numberFake_ = rhs.numberFake_;
  numberChanged_ = rhs.numberChanged_;
  progressFlag_ = rhs.progressFlag_;
  firstFree_ = rhs.firstFree_;
  incomingInfeasibility_ = rhs.incomingInfeasibility_;
  allowedInfeasibility_ = rhs.allowedInfeasibility_;
  automaticScale_ = rhs.automaticScale_;
  if (rhs.progress_)
    progress_ = new ClpSimplexProgress(*rhs.progress_);
  else
    progress_=NULL;
  sumOfRelaxedDualInfeasibilities_ = rhs.sumOfRelaxedDualInfeasibilities_;
  sumOfRelaxedPrimalInfeasibilities_ = rhs.sumOfRelaxedPrimalInfeasibilities_;
  if (rhs.nonLinearCost_!=NULL)
    nonLinearCost_ = new ClpNonLinearCost(*rhs.nonLinearCost_);
  else
    nonLinearCost_=NULL;
  solveType_=rhs.solveType_;
}
// type == 0 do everything, most + pivot data, 2 factorization data as well
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
  delete [] savedSolution_;
  savedSolution_ = NULL;
  if ((specialOptions_&2)==0) {
    delete nonLinearCost_;
    nonLinearCost_ = NULL;
  }
  int i;
  for (i=0;i<6;i++) {
    delete rowArray_[i];
    rowArray_[i]=NULL;
    delete columnArray_[i];
    columnArray_[i]=NULL;
  }
  delete rowCopy_;
  rowCopy_=NULL;
  delete [] saveStatus_;
  saveStatus_=NULL;
  if (!type) {
    // delete everything
    delete factorization_;
    factorization_ = NULL;
    delete [] pivotVariable_;
    pivotVariable_=NULL;
    delete dualRowPivot_;
    dualRowPivot_ = NULL;
    delete primalColumnPivot_;
    primalColumnPivot_ = NULL;
    delete progress_;
    progress_=NULL;
  } else {
    // delete any size information in methods
    if (type>1) {
      factorization_->clearArrays();
      delete [] pivotVariable_;
      pivotVariable_=NULL;
    }
    dualRowPivot_->clearArrays();
    primalColumnPivot_->clearArrays();
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
  double primalTolerance = primalTolerance_;
  double relaxedTolerance=primalTolerance_;
  // we can't really trust infeasibilities if there is primal error
  double error = CoinMin(1.0e-3,largestPrimalError_);
  // allow tolerance at least slightly bigger than standard
  relaxedTolerance = relaxedTolerance +  error;
  sumOfRelaxedPrimalInfeasibilities_ = 0.0;

  for (iRow=0;iRow<numberRows_;iRow++) {
    //assert (fabs(solution[iRow])<1.0e15||getRowStatus(iRow) == basic);
    double infeasibility=0.0;
    objectiveValue_ += solution[iRow]*rowObjectiveWork_[iRow];
    if (solution[iRow]>rowUpperWork_[iRow]) {
      infeasibility=solution[iRow]-rowUpperWork_[iRow];
    } else if (solution[iRow]<rowLowerWork_[iRow]) {
      infeasibility=rowLowerWork_[iRow]-solution[iRow];
    }
    if (infeasibility>primalTolerance) {
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
      if (infeasibility>relaxedTolerance) 
	sumOfRelaxedPrimalInfeasibilities_ += infeasibility-relaxedTolerance;
      numberPrimalInfeasibilities_ ++;
    }
    if (infeasibility>rowPrimalInfeasibility_) {
      rowPrimalInfeasibility_=infeasibility;
      rowPrimalSequence_=iRow;
    }
    infeasibility = fabs(rowActivities[iRow]-solution[iRow]);
    if (infeasibility>largestSolutionError_)
      largestSolutionError_=infeasibility;
  }
  // Check any infeasibilities from dynamic rows
  matrix_->primalExpanded(this,2);
  solution = columnActivityWork_;
  if (!matrix_->rhsOffset(this)) {
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      //assert (fabs(solution[iColumn])<1.0e15||getColumnStatus(iColumn) == basic);
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
      if (infeasibility>primalTolerance) {
	sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
	if (infeasibility>relaxedTolerance) 
	  sumOfRelaxedPrimalInfeasibilities_ += infeasibility-relaxedTolerance;
	numberPrimalInfeasibilities_ ++;
      }
      infeasibility = fabs(columnActivities[iColumn]-solution[iColumn]);
      if (infeasibility>largestSolutionError_)
	largestSolutionError_=infeasibility;
    }
  } else {
    // as we are using effective rhs we only check basics
    // But we do need to get objective
    objectiveValue_ += innerProduct(objectiveWork_,numberColumns_,solution);
    for (int j=0;j<numberRows_;j++) {
      int iColumn = pivotVariable_[j];
      //assert (fabs(solution[iColumn])<1.0e15||getColumnStatus(iColumn) == basic);
      double infeasibility=0.0;
      if (solution[iColumn]>columnUpperWork_[iColumn]) {
	infeasibility=solution[iColumn]-columnUpperWork_[iColumn];
      } else if (solution[iColumn]<columnLowerWork_[iColumn]) {
	infeasibility=columnLowerWork_[iColumn]-solution[iColumn];
      }
      if (infeasibility>columnPrimalInfeasibility_) {
	columnPrimalInfeasibility_=infeasibility;
	columnPrimalSequence_=iColumn;
      }
      if (infeasibility>primalTolerance) {
	sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
	if (infeasibility>relaxedTolerance) 
	  sumOfRelaxedPrimalInfeasibilities_ += infeasibility-relaxedTolerance;
	numberPrimalInfeasibilities_ ++;
      }
      infeasibility = fabs(columnActivities[iColumn]-solution[iColumn]);
      if (infeasibility>largestSolutionError_)
	largestSolutionError_=infeasibility;
    }
  }
  objectiveValue_ += objective_->nonlinearOffset();
  objectiveValue_ /= (objectiveScale_*rhsScale_);
}
void 
ClpSimplex::checkDualSolution()
{

  int iRow,iColumn;
  sumDualInfeasibilities_=0.0;
  numberDualInfeasibilities_=0;
  numberDualInfeasibilitiesWithoutFree_=0;
  columnDualInfeasibility_=0.0;
  columnDualSequence_=-1;
  if (matrix_->skipDualCheck()&&algorithm_>0&&problemStatus_==-2) {
    // pretend we found dual infeasibilities
    sumOfRelaxedDualInfeasibilities_ = 1.0;
    sumDualInfeasibilities_=1.0;
    numberDualInfeasibilities_=1;
    return;
  }
  rowDualInfeasibility_=0.0;
  rowDualSequence_=-1;
  int firstFreePrimal = -1;
  int firstFreeDual = -1;
  int numberSuperBasicWithDj=0;
  primalToleranceToGetOptimal_=CoinMax(rowPrimalInfeasibility_,
				   columnPrimalInfeasibility_);
  remainingDualInfeasibility_=0.0;
  double relaxedTolerance=dualTolerance_;
  // we can't really trust infeasibilities if there is dual error
  double error = CoinMin(1.0e-3,largestDualError_);
  // allow tolerance at least slightly bigger than standard
  relaxedTolerance = relaxedTolerance +  error;
  sumOfRelaxedDualInfeasibilities_ = 0.0;

  // Check any djs from dynamic rows
  matrix_->dualExpanded(this,NULL,NULL,3);
  numberDualInfeasibilitiesWithoutFree_= numberDualInfeasibilities_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (getColumnStatus(iColumn) != basic&&!flagged(iColumn)) {
      // not basic
      double distanceUp = columnUpperWork_[iColumn]-
	columnActivityWork_[iColumn];
      double distanceDown = columnActivityWork_[iColumn] -
	columnLowerWork_[iColumn];
      if (distanceUp>primalTolerance_) {
	double value = reducedCostWork_[iColumn];
	// Check if "free"
	if (distanceDown>primalTolerance_) {
	  if (fabs(value)>1.0e2*relaxedTolerance) {
	    numberSuperBasicWithDj++;
	    if (firstFreeDual<0)
	      firstFreeDual = iColumn;
	  }
	  if (firstFreePrimal<0)
	    firstFreePrimal = iColumn;
	}
	// should not be negative
	if (value<0.0) {
	  value = - value;
	  if (value>columnDualInfeasibility_) {
	    columnDualInfeasibility_=value;
	    columnDualSequence_=iColumn;
	  }
	  if (value>dualTolerance_) {
	    if (getColumnStatus(iColumn) != isFree) {
	      numberDualInfeasibilitiesWithoutFree_ ++;
	      sumDualInfeasibilities_ += value-dualTolerance_;
	      if (value>relaxedTolerance) 
		sumOfRelaxedDualInfeasibilities_ += value-relaxedTolerance;
	      numberDualInfeasibilities_ ++;
	    } else {
	      // free so relax a lot
	      value *= 0.01;
	      if (value>dualTolerance_) {
		sumDualInfeasibilities_ += value-dualTolerance_;
		if (value>relaxedTolerance) 
		  sumOfRelaxedDualInfeasibilities_ += value-relaxedTolerance;
		numberDualInfeasibilities_ ++;
	      }
	    }
	    // maybe we can make feasible by increasing tolerance
	    if (distanceUp<largeValue_) {
	      if (distanceUp>primalToleranceToGetOptimal_)
		primalToleranceToGetOptimal_=distanceUp;
	    } else {
	      //gap too big for any tolerance
	      remainingDualInfeasibility_=
		CoinMax(remainingDualInfeasibility_,value);
	    }
	  }
	}
      }
      if (distanceDown>primalTolerance_) {
	double value = reducedCostWork_[iColumn];
	// should not be positive
	if (value>0.0) {
	  if (value>columnDualInfeasibility_) {
	    columnDualInfeasibility_=value;
	    columnDualSequence_=iColumn;
	  }
	  if (value>dualTolerance_) {
	    sumDualInfeasibilities_ += value-dualTolerance_;
	    if (value>relaxedTolerance) 
	      sumOfRelaxedDualInfeasibilities_ += value-relaxedTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getColumnStatus(iColumn) != isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distanceDown<largeValue_&&
		distanceDown>primalToleranceToGetOptimal_)
	      primalToleranceToGetOptimal_=distanceDown;
	  }
	}
      }
    }
  }
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (getRowStatus(iRow) != basic&&!flagged(iRow+numberColumns_)) {
      // not basic
      double distanceUp = rowUpperWork_[iRow]-rowActivityWork_[iRow];
      double distanceDown = rowActivityWork_[iRow] -rowLowerWork_[iRow];
      if (distanceUp>primalTolerance_) {
	double value = rowReducedCost_[iRow];
	// Check if "free"
	if (distanceDown>primalTolerance_) {
	  if (fabs(value)>1.0e2*relaxedTolerance) {
	    numberSuperBasicWithDj++;
	    if (firstFreeDual<0)
	      firstFreeDual = iRow+numberColumns_;
	  }
	  if (firstFreePrimal<0)
	    firstFreePrimal = iRow+numberColumns_;
	}
	// should not be negative
	if (value<0.0) {
	  value = - value;
	  if (value>rowDualInfeasibility_) {
	    rowDualInfeasibility_=value;
	    rowDualSequence_=iRow;
	  }
	  if (value>dualTolerance_) {
	    sumDualInfeasibilities_ += value-dualTolerance_;
	    if (value>relaxedTolerance) 
	      sumOfRelaxedDualInfeasibilities_ += value-relaxedTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getRowStatus(iRow) != isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distanceUp<largeValue_) {
	      if (distanceUp>primalToleranceToGetOptimal_)
		primalToleranceToGetOptimal_=distanceUp;
	    } else {
	      //gap too big for any tolerance
	      remainingDualInfeasibility_=
		CoinMax(remainingDualInfeasibility_,value);
	    }
	  }
	}
      }
      if (distanceDown>primalTolerance_) {
	double value = rowReducedCost_[iRow];
	// should not be positive
	if (value>0.0) {
	  if (value>rowDualInfeasibility_) {
	    rowDualInfeasibility_=value;
	    rowDualSequence_=iRow;
	  }
	  if (value>dualTolerance_) {
	    sumDualInfeasibilities_ += value-dualTolerance_;
	    if (value>relaxedTolerance) 
	      sumOfRelaxedDualInfeasibilities_ += value-relaxedTolerance;
	    numberDualInfeasibilities_ ++;
	    if (getRowStatus(iRow) != isFree) 
	      numberDualInfeasibilitiesWithoutFree_ ++;
	    // maybe we can make feasible by increasing tolerance
	    if (distanceDown<largeValue_&&
		distanceDown>primalToleranceToGetOptimal_)
	      primalToleranceToGetOptimal_=distanceDown;
	  }
	}
      }
    }
  }
  if (algorithm_<0&&firstFreeDual>=0) {
    // dual
    firstFree_ = firstFreeDual;
  } else if (numberSuperBasicWithDj||
	     (progress_&&progress_->lastIterationNumber(0)<=0)) {
    firstFree_=firstFreePrimal;
  }
}
/* This sets sum and number of infeasibilities (Dual and Primal) */
void 
ClpSimplex::checkBothSolutions()
{
  if ((matrix_->skipDualCheck()&&algorithm_>0&&problemStatus_==-2)||
      matrix_->rhsOffset(this)) {
    // old way
    checkPrimalSolution(rowActivityWork_,columnActivityWork_);
    checkDualSolution();
    return;
  }
  int iSequence;

  objectiveValue_ = 0.0;
  sumPrimalInfeasibilities_=0.0;
  numberPrimalInfeasibilities_=0;
  double primalTolerance = primalTolerance_;
  double relaxedToleranceP=primalTolerance_;
  // we can't really trust infeasibilities if there is primal error
  double error = CoinMin(1.0e-3,largestPrimalError_);
  // allow tolerance at least slightly bigger than standard
  relaxedToleranceP = relaxedToleranceP +  error;
  sumOfRelaxedPrimalInfeasibilities_ = 0.0;
  sumDualInfeasibilities_=0.0;
  numberDualInfeasibilities_=0;
  double dualTolerance=dualTolerance_;
  double relaxedToleranceD=dualTolerance;
  // we can't really trust infeasibilities if there is dual error
  error = CoinMin(1.0e-3,largestDualError_);
  // allow tolerance at least slightly bigger than standard
  relaxedToleranceD = relaxedToleranceD +  error;
  sumOfRelaxedDualInfeasibilities_ = 0.0;

  // Check any infeasibilities from dynamic rows
  matrix_->primalExpanded(this,2);
  // Check any djs from dynamic rows
  matrix_->dualExpanded(this,NULL,NULL,3);
  int numberDualInfeasibilitiesFree= 0;
  int firstFreePrimal = -1;
  int firstFreeDual = -1;
  int numberSuperBasicWithDj=0;

  int numberTotal = numberRows_ + numberColumns_;
  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    double value = solution_[iSequence];
#ifdef COIN_DEBUG
    if (fabs(value)>1.0e20)
      printf("%d values %g %g %g - status %d\n",iSequence,lower_[iSequence],
             solution_[iSequence],upper_[iSequence],status_[iSequence]);
#endif
    objectiveValue_ += value*cost_[iSequence];
    double distanceUp =upper_[iSequence]-value;
    double distanceDown =value-lower_[iSequence];
    if (distanceUp<-primalTolerance) {
      double infeasibility=-distanceUp;
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
      if (infeasibility>relaxedToleranceP) 
	sumOfRelaxedPrimalInfeasibilities_ += infeasibility-relaxedToleranceP;
      numberPrimalInfeasibilities_ ++;
    } else if (distanceDown<-primalTolerance) {
      double infeasibility=-distanceDown;
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance_;
      if (infeasibility>relaxedToleranceP) 
	sumOfRelaxedPrimalInfeasibilities_ += infeasibility-relaxedToleranceP;
      numberPrimalInfeasibilities_ ++;
    } else {
      // feasible (so could be free)
      if (getStatus(iSequence) != basic&&!flagged(iSequence)) {
        // not basic
        double djValue = dj_[iSequence];
        if (distanceDown<primalTolerance) {
          if (distanceUp>primalTolerance&&djValue<-dualTolerance) {
	    sumDualInfeasibilities_ -= djValue+dualTolerance;
	    if (djValue<-relaxedToleranceD) 
	      sumOfRelaxedDualInfeasibilities_ -= djValue+relaxedToleranceD;
	    numberDualInfeasibilities_ ++;
          } 
        } else if (distanceUp<primalTolerance) {
          if (djValue>dualTolerance) {
	    sumDualInfeasibilities_ += djValue-dualTolerance;
	    if (djValue>relaxedToleranceD) 
	      sumOfRelaxedDualInfeasibilities_ += djValue-relaxedToleranceD;
	    numberDualInfeasibilities_ ++;
          }
        } else {
          // may be free
          djValue *= 0.01;
	  if (fabs(djValue)>dualTolerance) {
	    if (getStatus(iSequence) == isFree) 
              numberDualInfeasibilitiesFree++;
            sumDualInfeasibilities_ += fabs(djValue)-dualTolerance;
            numberDualInfeasibilities_ ++;
            if (fabs(djValue)>relaxedToleranceD) {
              sumOfRelaxedDualInfeasibilities_ += value-relaxedToleranceD;
              numberSuperBasicWithDj++;
              if (firstFreeDual<0)
                firstFreeDual = iSequence;
            }
          }
	  if (firstFreePrimal<0)
	    firstFreePrimal = iSequence;
	}
      }
    }
  }
  objectiveValue_ += objective_->nonlinearOffset();
  objectiveValue_ /= (objectiveScale_*rhsScale_);
  numberDualInfeasibilitiesWithoutFree_= numberDualInfeasibilities_-
    numberDualInfeasibilitiesFree;
  if (algorithm_<0&&firstFreeDual>=0) {
    // dual
    firstFree_ = firstFreeDual;
  } else if (numberSuperBasicWithDj||
	     (progress_&&progress_->lastIterationNumber(0)<=0)) {
    firstFree_=firstFreePrimal;
  }
}
/* Adds multiple of a column into an array */
void 
ClpSimplex::add(double * array,
		int sequence, double multiplier) const
{
  if (sequence>=numberColumns_&&sequence<numberColumns_+numberRows_) {
    //slack
    array [sequence-numberColumns_] -= multiplier;
  } else {
    // column
    matrix_->add(this,array,sequence,multiplier);
  }
}
/*
  Unpacks one column of the matrix into indexed array 
*/
void 
ClpSimplex::unpack(CoinIndexedVector * rowArray) const
{
  rowArray->clear();
  if (sequenceIn_>=numberColumns_&&sequenceIn_<numberColumns_+numberRows_) {
    //slack
    rowArray->insert(sequenceIn_-numberColumns_,-1.0);
  } else {
    // column
    matrix_->unpack(this,rowArray,sequenceIn_);
  }
}
void 
ClpSimplex::unpack(CoinIndexedVector * rowArray,int sequence) const
{
  rowArray->clear();
  if (sequence>=numberColumns_&&sequence<numberColumns_+numberRows_) {
    //slack
    rowArray->insert(sequence-numberColumns_,-1.0);
  } else {
    // column
    matrix_->unpack(this,rowArray,sequence);
  }
}
/*
  Unpacks one column of the matrix into indexed array 
*/
void 
ClpSimplex::unpackPacked(CoinIndexedVector * rowArray) 
{
  rowArray->clear();
  if (sequenceIn_>=numberColumns_&&sequenceIn_<numberColumns_+numberRows_) {
    //slack
    int * index = rowArray->getIndices();
    double * array = rowArray->denseVector();
    array[0]=-1.0;
    index[0]=sequenceIn_-numberColumns_;
    rowArray->setNumElements(1);
    rowArray->setPackedMode(true);
  } else {
    // column
    matrix_->unpackPacked(this,rowArray,sequenceIn_);
  }
}
void 
ClpSimplex::unpackPacked(CoinIndexedVector * rowArray,int sequence)
{
  rowArray->clear();
  if (sequence>=numberColumns_&&sequence<numberColumns_+numberRows_) {
    //slack
    int * index = rowArray->getIndices();
    double * array = rowArray->denseVector();
    array[0]=-1.0;
    index[0]=sequence-numberColumns_;
    rowArray->setNumElements(1);
    rowArray->setPackedMode(true);
  } else {
    // column
    matrix_->unpackPacked(this,rowArray,sequence);
  }
}
bool
ClpSimplex::createRim(int what,bool makeRowCopy, int startFinishOptions)
{
  bool goodMatrix=true;
  int saveLevel=handler_->logLevel();
  // Arrays will be there and correct size unless what is 63
  bool newArrays = (what==63);
  // We may be restarting with same size
  bool keepPivots=false;
  if (startFinishOptions==-1) {
    startFinishOptions=0;
    keepPivots=true;
  }
  bool oldMatrix = ((startFinishOptions&4)!=0&&(whatsChanged_&1)!=0);
  if (what==63) {
    if (!status_)
      createStatus();
    if (oldMatrix)
      newArrays=false;
    if (problemStatus_==10) {
      handler_->setLogLevel(0); // switch off messages
      if (rowArray_[0]) {
        // stuff is still there
        oldMatrix=true;
        newArrays=false;
        keepPivots=true;
      }
    } else if (factorization_) {
      // match up factorization messages
      if (handler_->logLevel()<3)
        factorization_->messageLevel(0);
      else
        factorization_->messageLevel(CoinMax(3,factorization_->messageLevel()));
      /* Faster to keep pivots rather than re-scan matrix.  Matrix may have changed
         i.e. oldMatrix false but okay as long as same number rows and status array exists
      */
      if ((startFinishOptions&2)!=0&&factorization_->numberRows()==numberRows_&&status_)
        keepPivots=true;
    }
    numberExtraRows_ = matrix_->generalExpanded(this,2,maximumBasic_);
    if (numberExtraRows_&&newArrays) {
      // make sure status array large enough
      assert (status_);
      int numberOld = numberRows_+numberColumns_;
      int numberNew = numberRows_+numberColumns_+numberExtraRows_;
      unsigned char * newStatus = new unsigned char [numberNew];
      memset(newStatus+numberOld,0,numberExtraRows_);
      memcpy(newStatus,status_,numberOld);
      delete [] status_;
      status_=newStatus;
    }
  }
  int numberRows2 = numberRows_+numberExtraRows_;
  //int numberTotal = numberRows_+numberColumns_;
  int i;
  bool doSanityCheck=true;
  if (what==63) {
    // We may want to switch stuff off for speed
    if ((specialOptions_&256)!=0)
      makeRowCopy=false; // no row copy
    if ((specialOptions_&128)!=0)
      doSanityCheck=false; // no sanity check
    //check matrix
    if (!matrix_)
      matrix_=new ClpPackedMatrix();
    int checkType=(doSanityCheck) ? 15 : 14;
    if (oldMatrix)
      checkType = 14;
    if (!matrix_->allElementsInRange(this,smallElement_,1.0e20,checkType)) {
      problemStatus_=4;
      goodMatrix= false;
    }
    if (makeRowCopy&&!oldMatrix) {
      delete rowCopy_;
      // may return NULL if can't give row copy
      rowCopy_ = matrix_->reverseOrderedCopy();
    }
    // do scaling if needed
    if (!oldMatrix) {
      assert (scalingFlag_>0||!rowScale_);
      delete [] rowScale_;
      delete [] columnScale_;
      rowScale_=NULL;
      columnScale_=NULL;
    }
    if (scalingFlag_>0&&!rowScale_) {
      if (matrix_->scale(this))
	scalingFlag_=-scalingFlag_; // not scaled after all
      if (rowScale_&&automaticScale_) {
	// try automatic scaling
	double smallestObj=1.0e100;
	double largestObj=0.0;
	double largestRhs=0.0;
        const double * obj = objective();
	for (i=0;i<numberColumns_;i++) {
	  double value = fabs(obj[i]);
          value *= columnScale_[i];
	  if (value&&columnLower_[i]!=columnUpper_[i]) {
	    smallestObj = CoinMin(smallestObj,value);
	    largestObj = CoinMax(largestObj,value);
	  }
	  if (columnLower_[i]>0.0||columnUpper_[i]<0.0) {
	    double scale = 1.0/columnScale_[i];
	    //printf("%d %g %g %g %g\n",i,scale,lower_[i],upper_[i],largestRhs);
	    if (columnLower_[i]>0)
	      largestRhs=CoinMax(largestRhs,columnLower_[i]*scale);
	    if (columnUpper_[i]<0.0)
	      largestRhs=CoinMax(largestRhs,-columnUpper_[i]*scale);
	  }
	}
	for (i=0;i<numberRows_;i++) {
	  if (rowLower_[i]>0.0||rowUpper_[i]<0.0) {
	    double scale = rowScale_[i];
	    //printf("%d %g %g %g %g\n",i,scale,lower_[i],upper_[i],largestRhs);
	    if (rowLower_[i]>0)
	      largestRhs=CoinMax(largestRhs,rowLower_[i]*scale);
	    if (rowUpper_[i]<0.0)
	      largestRhs=CoinMax(largestRhs,-rowUpper_[i]*scale);
	  }
	}
	printf("small obj %g, large %g - rhs %g\n",smallestObj,largestObj,largestRhs);
	bool scalingDone=false;
	// look at element range
	double smallestNegative;
	double largestNegative;
	double smallestPositive;
	double largestPositive;
	matrix_->rangeOfElements(smallestNegative, largestNegative,
				 smallestPositive, largestPositive);
	smallestPositive = CoinMin(fabs(smallestNegative),smallestPositive);
	largestPositive = CoinMax(fabs(largestNegative),largestPositive);
	if (largestObj) {
	  double ratio = largestObj/smallestObj;
	  double scale=1.0;
	  if (ratio<1.0e8) {
	    // reasonable
	    if (smallestObj<1.0e-4) {
	      // may as well scale up
	      scalingDone=true;
	      scale=1.0e-3/smallestObj;
	    } else if (largestObj<1.0e6||(algorithm_>0&&largestObj<1.0e-4*infeasibilityCost_)) {
	      //done=true;
	    } else {
	      scalingDone=true;
	      if (algorithm_<0) {
		scale = 1.0e6/largestObj;
	      } else {
		scale = CoinMax(1.0e6,1.0e-4*infeasibilityCost_)/largestObj;
	      }
	    }
	  } else if (ratio<1.0e12) {
	    // not so good
	    if (smallestObj<1.0e-7) {
	      // may as well scale up
	      scalingDone=true;
	      scale=1.0e-6/smallestObj;
	    } else if (largestObj<1.0e7||(algorithm_>0&&largestObj<1.0e-3*infeasibilityCost_)) {
	      //done=true;
	    } else {
	      scalingDone=true;
	      if (algorithm_<0) {
		scale = 1.0e7/largestObj;
	      } else {
		scale = CoinMax(1.0e7,1.0e-3*infeasibilityCost_)/largestObj;
	      }
	    }
	  } else {
	    // Really nasty problem
	    if (smallestObj<1.0e-8) {
	      // may as well scale up
	      scalingDone=true;
	      scale=1.0e-7/smallestObj;
	      largestObj *= scale;
	    } 
	    if (largestObj<1.0e7||(algorithm_>0&&largestObj<1.0e-3*infeasibilityCost_)) {
	      //done=true;
	    } else {
	      scalingDone=true;
	      if (algorithm_<0) {
		scale = 1.0e7/largestObj;
	      } else {
		scale = CoinMax(1.0e7,1.0e-3*infeasibilityCost_)/largestObj;
	      }
	    }
	  }
	  objectiveScale_=scale;
	}
	if (largestRhs>1.0e12) {
	  scalingDone=true;
	  rhsScale_=1.0e9/largestRhs;
	} else if (largestPositive>1.0e-14*smallestPositive&&largestRhs>1.0e6) {
	  scalingDone=true;
	  rhsScale_=1.0e6/largestRhs;
	} else {
	  rhsScale_=1.0;
	}
	if (scalingDone) {
	  handler_->message(CLP_RIM_SCALE,messages_)
	    <<objectiveScale_<<rhsScale_
	    <<CoinMessageEol;
	}
      }
    } else if (makeRowCopy&&scalingFlag_>0&&!oldMatrix) {
      matrix_->scaleRowCopy(this);
    }
  }
  if (what==63) {
    if (newArrays) {
      delete [] cost_;
      // extra copy with original costs
      int nTotal = numberRows2+numberColumns_;
      //cost_ = new double[2*nTotal];
      cost_ = new double[nTotal];
      delete [] lower_;
      delete [] upper_;
      lower_ = new double[numberColumns_+numberRows2];
      upper_ = new double[numberColumns_+numberRows2];
      delete [] dj_;
      dj_ = new double[numberRows2+numberColumns_];
      reducedCostWork_ = dj_;
      rowReducedCost_ = dj_+numberColumns_;
      delete [] solution_;
      solution_ = new double[numberRows2+numberColumns_];
      columnActivityWork_ = solution_;
      rowActivityWork_ = solution_+numberColumns_;
    }
    objectiveWork_ = cost_;
    rowObjectiveWork_ = cost_+numberColumns_;
    rowLowerWork_ = lower_+numberColumns_;
    columnLowerWork_ = lower_;
    rowUpperWork_ = upper_+numberColumns_;
    columnUpperWork_ = upper_;
  }
  if ((what&4)!=0) {
    double direction = optimizationDirection_*objectiveScale_;
    const double * obj = objective();
    // and also scale by scale factors
    if (rowScale_) {
      if (rowObjective_) {
        for (i=0;i<numberRows_;i++)
          rowObjectiveWork_[i] = rowObjective_[i]*direction/rowScale_[i];
      } else {
        memset(rowObjectiveWork_,0,numberRows_*sizeof(double));
      }
      // If scaled then do all columns later in one loop 
      if (what!=63||!rowScale_) {
        for (i=0;i<numberColumns_;i++) {
          CoinAssert(fabs(obj[i])<1.0e25);
          objectiveWork_[i] = obj[i]*direction*columnScale_[i];
        }
      }
    } else {
      if (rowObjective_) {
        for (i=0;i<numberRows_;i++)
          rowObjectiveWork_[i] = rowObjective_[i]*direction;
      } else {
        memset(rowObjectiveWork_,0,numberRows_*sizeof(double));
      }
      for (i=0;i<numberColumns_;i++) {
        CoinAssert(fabs(obj[i])<1.0e25);
        objectiveWork_[i] = obj[i]*direction;
      }
    }
  }
  if ((what&1)!=0) {
    // clean up any mismatches on infinity
    // and fix any variables with tiny gaps
    double primalTolerance=dblParam_[ClpPrimalTolerance];
    if(rowScale_) {
      // If scaled then do all columns later in one loop 
      if (what!=63||!rowScale_) {
        for (i=0;i<numberColumns_;i++) {
          double multiplier = rhsScale_/columnScale_[i];
          double lowerValue = columnLower_[i];
          double upperValue = columnUpper_[i];
          if (lowerValue>-1.0e20) {
            columnLowerWork_[i]=lowerValue*multiplier;
            if (upperValue>=1.0e20) {
              columnUpperWork_[i]=COIN_DBL_MAX;
            } else {
              columnUpperWork_[i]=upperValue*multiplier;
              if (columnUpperWork_[i]-columnLowerWork_[i]<=primalTolerance) {
                if (columnLowerWork_[i]>=0.0) {
                  columnUpperWork_[i] = columnLowerWork_[i];
                } else if (columnUpperWork_[i]<=0.0) {
                  columnLowerWork_[i] = columnUpperWork_[i];
                } else {
                  columnUpperWork_[i] = 0.0;
                  columnLowerWork_[i] = 0.0;
                }
              }
            }
          } else if (upperValue<1.0e20) {
            columnLowerWork_[i]=-COIN_DBL_MAX;
            columnUpperWork_[i]=upperValue*multiplier;
          } else {
            // free
            columnLowerWork_[i]=-COIN_DBL_MAX;
            columnUpperWork_[i]=COIN_DBL_MAX;
          }
        }
      }
      for (i=0;i<numberRows_;i++) {
	double multiplier = rhsScale_*rowScale_[i];
        double lowerValue = rowLower_[i];
        double upperValue = rowUpper_[i];
        if (lowerValue>-1.0e20) {
          rowLowerWork_[i]=lowerValue*multiplier;
          if (upperValue>=1.0e20) {
            rowUpperWork_[i]=COIN_DBL_MAX;
          } else {
            rowUpperWork_[i]=upperValue*multiplier;
            if (rowUpperWork_[i]-rowLowerWork_[i]<=primalTolerance) {
              if (rowLowerWork_[i]>=0.0) {
                rowUpperWork_[i] = rowLowerWork_[i];
              } else if (rowUpperWork_[i]<=0.0) {
                rowLowerWork_[i] = rowUpperWork_[i];
              } else {
                rowUpperWork_[i] = 0.0;
                rowLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=upperValue*multiplier;
        } else {
          // free
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=COIN_DBL_MAX;
        }
      }
    } else if (rhsScale_!=1.0) {
      for (i=0;i<numberColumns_;i++) {
        double lowerValue = columnLower_[i];
        double upperValue = columnUpper_[i];
        if (lowerValue>-1.0e20) {
          columnLowerWork_[i]=lowerValue*rhsScale_;
          if (upperValue>=1.0e20) {
            columnUpperWork_[i]=COIN_DBL_MAX;
          } else {
            columnUpperWork_[i]=upperValue*rhsScale_;
            if (columnUpperWork_[i]-columnLowerWork_[i]<=primalTolerance) {
              if (columnLowerWork_[i]>=0.0) {
                columnUpperWork_[i] = columnLowerWork_[i];
              } else if (columnUpperWork_[i]<=0.0) {
                columnLowerWork_[i] = columnUpperWork_[i];
              } else {
                columnUpperWork_[i] = 0.0;
                columnLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=upperValue*rhsScale_;
        } else {
          // free
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=COIN_DBL_MAX;
        }
      }
      for (i=0;i<numberRows_;i++) {
        double lowerValue = rowLower_[i];
        double upperValue = rowUpper_[i];
        if (lowerValue>-1.0e20) {
          rowLowerWork_[i]=lowerValue*rhsScale_;
          if (upperValue>=1.0e20) {
            rowUpperWork_[i]=COIN_DBL_MAX;
          } else {
            rowUpperWork_[i]=upperValue*rhsScale_;
            if (rowUpperWork_[i]-rowLowerWork_[i]<=primalTolerance) {
              if (rowLowerWork_[i]>=0.0) {
                rowUpperWork_[i] = rowLowerWork_[i];
              } else if (rowUpperWork_[i]<=0.0) {
                rowLowerWork_[i] = rowUpperWork_[i];
              } else {
                rowUpperWork_[i] = 0.0;
                rowLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=upperValue*rhsScale_;
        } else {
          // free
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=COIN_DBL_MAX;
        }
      }
    } else {
      for (i=0;i<numberColumns_;i++) {
        double lowerValue = columnLower_[i];
        double upperValue = columnUpper_[i];
        if (lowerValue>-1.0e20) {
          columnLowerWork_[i]=lowerValue;
          if (upperValue>=1.0e20) {
            columnUpperWork_[i]=COIN_DBL_MAX;
          } else {
            columnUpperWork_[i]=upperValue;
            if (columnUpperWork_[i]-columnLowerWork_[i]<=primalTolerance) {
              if (columnLowerWork_[i]>=0.0) {
                columnUpperWork_[i] = columnLowerWork_[i];
              } else if (columnUpperWork_[i]<=0.0) {
                columnLowerWork_[i] = columnUpperWork_[i];
              } else {
                columnUpperWork_[i] = 0.0;
                columnLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=upperValue;
        } else {
          // free
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=COIN_DBL_MAX;
        }
      }
      for (i=0;i<numberRows_;i++) {
        double lowerValue = rowLower_[i];
        double upperValue = rowUpper_[i];
        if (lowerValue>-1.0e20) {
          rowLowerWork_[i]=lowerValue;
          if (upperValue>=1.0e20) {
            rowUpperWork_[i]=COIN_DBL_MAX;
          } else {
            rowUpperWork_[i]=upperValue;
            if (rowUpperWork_[i]-rowLowerWork_[i]<=primalTolerance) {
              if (rowLowerWork_[i]>=0.0) {
                rowUpperWork_[i] = rowLowerWork_[i];
              } else if (rowUpperWork_[i]<=0.0) {
                rowLowerWork_[i] = rowUpperWork_[i];
              } else {
                rowUpperWork_[i] = 0.0;
                rowLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=upperValue;
        } else {
          // free
          rowLowerWork_[i]=-COIN_DBL_MAX;
          rowUpperWork_[i]=COIN_DBL_MAX;
        }
      }
    }
  }
  if (what==63) {
    // move information to work arrays
    double direction = optimizationDirection_;
    // direction is actually scale out not scale in
    if (direction)
      direction = 1.0/direction;
    if (direction!=1.0) {
      // reverse all dual signs
      for (i=0;i<numberColumns_;i++) 
	reducedCost_[i] *= direction;
      for (i=0;i<numberRows_;i++) 
	dual_[i] *= direction;
    }
    for (i=0;i<numberRows_+numberColumns_;i++) {
      setFakeBound(i,noFake);
    }
    if (rowScale_) {
      const double * obj = objective();
      double direction = optimizationDirection_*objectiveScale_;
      // clean up any mismatches on infinity
      // and fix any variables with tiny gaps
      double primalTolerance=dblParam_[ClpPrimalTolerance];
      // on entry
      for (i=0;i<numberColumns_;i++) {
        CoinAssert(fabs(obj[i])<1.0e25);
        double scaleFactor = columnScale_[i];
	double multiplier = rhsScale_/scaleFactor;
        scaleFactor *= direction;
	objectiveWork_[i] = obj[i]*scaleFactor;
	reducedCostWork_[i] = reducedCost_[i]*scaleFactor;
        double lowerValue = columnLower_[i];
        double upperValue = columnUpper_[i];
        if (lowerValue>-1.0e20) {
          columnLowerWork_[i]=lowerValue*multiplier;
          if (upperValue>=1.0e20) {
            columnUpperWork_[i]=COIN_DBL_MAX;
          } else {
            columnUpperWork_[i]=upperValue*multiplier;
            if (columnUpperWork_[i]-columnLowerWork_[i]<=primalTolerance) {
              if (columnLowerWork_[i]>=0.0) {
                columnUpperWork_[i] = columnLowerWork_[i];
              } else if (columnUpperWork_[i]<=0.0) {
                columnLowerWork_[i] = columnUpperWork_[i];
              } else {
                columnUpperWork_[i] = 0.0;
                columnLowerWork_[i] = 0.0;
              }
            }
          }
        } else if (upperValue<1.0e20) {
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=upperValue*multiplier;
        } else {
          // free
          columnLowerWork_[i]=-COIN_DBL_MAX;
          columnUpperWork_[i]=COIN_DBL_MAX;
        }
        double value = columnActivity_[i] * multiplier;
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for column %d\n",value,i);
          setColumnStatus(i,superBasic);
          if (columnUpperWork_[i]<0.0) {
            value=columnUpperWork_[i];
          } else if (columnLowerWork_[i]>0.0) {
            value=columnLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	columnActivityWork_[i] = value;
      }
      for (i=0;i<numberRows_;i++) {
	dual_[i] /= rowScale_[i];
	dual_[i] *= objectiveScale_;
	rowReducedCost_[i] = dual_[i];
	double multiplier = rhsScale_*rowScale_[i];
        double value = rowActivity_[i] * multiplier;
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for row %d\n",value,i);
          setRowStatus(i,superBasic);
          if (rowUpperWork_[i]<0.0) {
            value=rowUpperWork_[i];
          } else if (rowLowerWork_[i]>0.0) {
            value=rowLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	rowActivityWork_[i] = value;
      }
    } else if (objectiveScale_!=1.0||rhsScale_!=1.0) {
      // on entry
      for (i=0;i<numberColumns_;i++) {
        double value = columnActivity_[i];
	value *= rhsScale_;
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for column %d\n",value,i);
          setColumnStatus(i,superBasic);
          if (columnUpperWork_[i]<0.0) {
            value=columnUpperWork_[i];
          } else if (columnLowerWork_[i]>0.0) {
            value=columnLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	columnActivityWork_[i] = value;
	reducedCostWork_[i] = reducedCost_[i]*objectiveScale_;
      }
      for (i=0;i<numberRows_;i++) {
        double value = rowActivity_[i];
	value *= rhsScale_;
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for row %d\n",value,i);
          setRowStatus(i,superBasic);
          if (rowUpperWork_[i]<0.0) {
            value=rowUpperWork_[i];
          } else if (rowLowerWork_[i]>0.0) {
            value=rowLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	rowActivityWork_[i] = value;
	dual_[i] *= objectiveScale_;
	rowReducedCost_[i] = dual_[i];
      }
    } else {
      // on entry
      for (i=0;i<numberColumns_;i++) {
        double value = columnActivity_[i];
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for column %d\n",value,i);
          setColumnStatus(i,superBasic);
          if (columnUpperWork_[i]<0.0) {
            value=columnUpperWork_[i];
          } else if (columnLowerWork_[i]>0.0) {
            value=columnLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	columnActivityWork_[i] = value;
	reducedCostWork_[i] = reducedCost_[i];
      }
      for (i=0;i<numberRows_;i++) {
        double value = rowActivity_[i];
        if (fabs(value)>1.0e20) {
          //printf("bad value of %g for row %d\n",value,i);
          setRowStatus(i,superBasic);
          if (rowUpperWork_[i]<0.0) {
            value=rowUpperWork_[i];
          } else if (rowLowerWork_[i]>0.0) {
            value=rowLowerWork_[i];
          } else {
            value=0.0;
          }
        }
	rowActivityWork_[i] = value;
	rowReducedCost_[i] = dual_[i];
      }
    }
  }
 
  if (what==63&&doSanityCheck) {
    // check rim of problem okay
    if (!sanityCheck())
      goodMatrix=false;
  } 
  // we need to treat matrix as if each element by rowScaleIn and columnScaleout??
  // maybe we need to move scales to SimplexModel for factorization?
  if ((what==63&&!pivotVariable_)||(newArrays&&!keepPivots)) {
    delete [] pivotVariable_;
    pivotVariable_=new int[numberRows2];
    for (int i=0;i<numberRows2;i++)
      pivotVariable_[i]=-1;
  } else if (what==63&&!keepPivots) {
    // just reset
    for (int i=0;i<numberRows2;i++)
      pivotVariable_[i]=-1;
  } else if (what==63) {
    // check pivots
    for (int i=0;i<numberRows2;i++) {
      int iSequence = pivotVariable_[i];
      if (iSequence<numberRows_+numberColumns_&&
          getStatus(iSequence)!=basic) {
	keepPivots =false;
	break;
      }
    }
    if (!keepPivots) {
      // reset
      for (int i=0;i<numberRows2;i++)
	pivotVariable_[i]=-1;
    } else {
      // clean
      for (int i=0;i<numberColumns_+numberRows_;i++) {
	Status status =getStatus(i);
	if (status!=basic) {
	  if (upper_[i]==lower_[i]) {
	    setStatus(i,isFixed);
	    solution_[i]=lower_[i];
	  } else if (status==atLowerBound) {
	    if (lower_[i]>-1.0e20) {
              solution_[i]=lower_[i];
            } else {
              //printf("seq %d at lower of %g\n",i,lower_[i]);
              if (upper_[i]<1.0e20) {
                solution_[i]=upper_[i];
                setStatus(i,atUpperBound);
              } else {
                setStatus(i,isFree);
              }
            }
	  } else if (status==atUpperBound) {
	    if (upper_[i]<1.0e20) {
              solution_[i]=upper_[i];
            } else {
              //printf("seq %d at upper of %g\n",i,upper_[i]);
              if (lower_[i]>-1.0e20) {
                solution_[i]=lower_[i];
                setStatus(i,atLowerBound);
              } else {
                setStatus(i,isFree);
              }
            }
	  }
	}
      }
    }
  }

  if (what==63) {
    if (newArrays) {
      // get some arrays
      int iRow,iColumn;
      // these are "indexed" arrays so we always know where nonzeros are
      /**********************************************************
      rowArray_[3] is long enough for rows+columns
      *********************************************************/
      for (iRow=0;iRow<4;iRow++) {
	delete rowArray_[iRow];
	rowArray_[iRow]=new CoinIndexedVector();
	int length =numberRows2+factorization_->maximumPivots();
	if (iRow==3||objective_->type()>1)
	  length += numberColumns_;
	rowArray_[iRow]->reserve(length);
      }
      
      for (iColumn=0;iColumn<2;iColumn++) {
	delete columnArray_[iColumn];
	columnArray_[iColumn]=new CoinIndexedVector();
	if (!iColumn)
	  columnArray_[iColumn]->reserve(numberColumns_);
	else
	  columnArray_[iColumn]->reserve(CoinMax(numberRows2,numberColumns_));
      }
    }    
  }
  if (problemStatus_==10) {
    problemStatus_=-1;
    handler_->setLogLevel(saveLevel); // switch back messages
  }
  if ((what&5)!=0) 
    matrix_->generalExpanded(this,9,what); // update costs and bounds if necessary
  return goodMatrix;
}
void
ClpSimplex::deleteRim(int getRidOfFactorizationData)
{
  int i;
  if (problemStatus_!=1&&problemStatus_!=2) {
    delete [] ray_;
    ray_=NULL;
  }
  // set upperOut_ to furthest away from bound so can use in dual for dualBound_
  upperOut_=1.0;
#if 0
  {
    int nBad=0;
    for (i=0;i<numberColumns_;i++) {
      if (lower_[i]==upper_[i]&&getColumnStatus(i)==basic)
        nBad++;
    }
    if (nBad)
      printf("yy %d basic fixed\n",nBad);
  }
#endif
  // ray may be null if in branch and bound
  if (rowScale_) {
    // Collect infeasibilities
    int numberPrimalScaled=0;
    int numberPrimalUnscaled=0;
    int numberDualScaled=0;
    int numberDualUnscaled=0;
    double scaleC = 1.0/objectiveScale_;
    double scaleR = 1.0/rhsScale_;
    for (i=0;i<numberColumns_;i++) {
      double scaleFactor = columnScale_[i];
      double valueScaled = columnActivityWork_[i];
      double lowerScaled = columnLowerWork_[i];
      double upperScaled = columnUpperWork_[i];
      if (lowerScaled>-1.0e20||upperScaled<1.0e20) {
	if (valueScaled<lowerScaled-primalTolerance_||
	    valueScaled>upperScaled+primalTolerance_)
	  numberPrimalScaled++;
	else
	  upperOut_ = CoinMax(upperOut_,CoinMin(valueScaled-lowerScaled,upperScaled-valueScaled));
      }
      columnActivity_[i] = valueScaled*scaleFactor*scaleR;
      double value = columnActivity_[i];
      if (value<columnLower_[i]-primalTolerance_)
	numberPrimalUnscaled++;
      else if (value>columnUpper_[i]+primalTolerance_)
	numberPrimalUnscaled++;
      double valueScaledDual = reducedCostWork_[i];
      if (valueScaled>columnLowerWork_[i]+primalTolerance_&&valueScaledDual>dualTolerance_)
	numberDualScaled++;
      if (valueScaled<columnUpperWork_[i]-primalTolerance_&&valueScaledDual<-dualTolerance_)
	numberDualScaled++;
      reducedCost_[i] = (valueScaledDual*scaleC)/scaleFactor;
      double valueDual = reducedCost_[i];
      if (value>columnLower_[i]+primalTolerance_&&valueDual>dualTolerance_)
	numberDualUnscaled++;
      if (value<columnUpper_[i]-primalTolerance_&&valueDual<-dualTolerance_)
	numberDualUnscaled++;
    }
    for (i=0;i<numberRows_;i++) {
      double scaleFactor = rowScale_[i];
      double valueScaled = rowActivityWork_[i];
      double lowerScaled = rowLowerWork_[i];
      double upperScaled = rowUpperWork_[i];
      if (lowerScaled>-1.0e20||upperScaled<1.0e20) {
	if (valueScaled<lowerScaled-primalTolerance_||
	    valueScaled>upperScaled+primalTolerance_)
	  numberPrimalScaled++;
	else
	  upperOut_ = CoinMax(upperOut_,CoinMin(valueScaled-lowerScaled,upperScaled-valueScaled));
      }
      rowActivity_[i] = (valueScaled*scaleR)/scaleFactor;
      double value = rowActivity_[i];
      if (value<rowLower_[i]-primalTolerance_)
	numberPrimalUnscaled++;
      else if (value>rowUpper_[i]+primalTolerance_)
	numberPrimalUnscaled++;
      double valueScaledDual = dual_[i]+rowObjectiveWork_[i];;
      if (valueScaled>rowLowerWork_[i]+primalTolerance_&&valueScaledDual>dualTolerance_)
	numberDualScaled++;
      if (valueScaled<rowUpperWork_[i]-primalTolerance_&&valueScaledDual<-dualTolerance_)
	numberDualScaled++;
      dual_[i] *= scaleFactor*scaleC;
      double valueDual = dual_[i]; 
      if (rowObjective_)
	valueDual += rowObjective_[i];
      if (value>rowLower_[i]+primalTolerance_&&valueDual>dualTolerance_)
	numberDualUnscaled++;
      if (value<rowUpper_[i]-primalTolerance_&&valueDual<-dualTolerance_)
	numberDualUnscaled++;
    }
    if (!problemStatus_&&!secondaryStatus_) {
      // See if we need to set secondary status
      if (numberPrimalUnscaled) {
	if (numberDualUnscaled) 
	  secondaryStatus_=4;
	else
	  secondaryStatus_=2;
      } else {
	if (numberDualUnscaled) 
	  secondaryStatus_=3;
      }
    }
    if (problemStatus_==2) {
      for (i=0;i<numberColumns_;i++) {
	ray_[i] *= columnScale_[i];
      }
    } else if (problemStatus_==1&&ray_) {
      for (i=0;i<numberRows_;i++) {
	ray_[i] *= rowScale_[i];
      }
    }
  } else if (rhsScale_!=1.0||objectiveScale_!=1.0) {
    // Collect infeasibilities
    int numberPrimalScaled=0;
    int numberPrimalUnscaled=0;
    int numberDualScaled=0;
    int numberDualUnscaled=0;
    double scaleC = 1.0/objectiveScale_;
    double scaleR = 1.0/rhsScale_;
    for (i=0;i<numberColumns_;i++) {
      double valueScaled = columnActivityWork_[i];
      double lowerScaled = columnLowerWork_[i];
      double upperScaled = columnUpperWork_[i];
      if (lowerScaled>-1.0e20||upperScaled<1.0e20) {
	if (valueScaled<lowerScaled-primalTolerance_||
	    valueScaled>upperScaled+primalTolerance_)
	  numberPrimalScaled++;
	else
	  upperOut_ = CoinMax(upperOut_,CoinMin(valueScaled-lowerScaled,upperScaled-valueScaled));
      }
      columnActivity_[i] = valueScaled*scaleR;
      double value = columnActivity_[i];
      if (value<columnLower_[i]-primalTolerance_)
	numberPrimalUnscaled++;
      else if (value>columnUpper_[i]+primalTolerance_)
	numberPrimalUnscaled++;
      double valueScaledDual = reducedCostWork_[i];
      if (valueScaled>columnLowerWork_[i]+primalTolerance_&&valueScaledDual>dualTolerance_)
	numberDualScaled++;
      if (valueScaled<columnUpperWork_[i]-primalTolerance_&&valueScaledDual<-dualTolerance_)
	numberDualScaled++;
      reducedCost_[i] = valueScaledDual*scaleC;
      double valueDual = reducedCost_[i];
      if (value>columnLower_[i]+primalTolerance_&&valueDual>dualTolerance_)
	numberDualUnscaled++;
      if (value<columnUpper_[i]-primalTolerance_&&valueDual<-dualTolerance_)
	numberDualUnscaled++;
    }
    for (i=0;i<numberRows_;i++) {
      double valueScaled = rowActivityWork_[i];
      double lowerScaled = rowLowerWork_[i];
      double upperScaled = rowUpperWork_[i];
      if (lowerScaled>-1.0e20||upperScaled<1.0e20) {
	if (valueScaled<lowerScaled-primalTolerance_||
	    valueScaled>upperScaled+primalTolerance_)
	  numberPrimalScaled++;
	else
	  upperOut_ = CoinMax(upperOut_,CoinMin(valueScaled-lowerScaled,upperScaled-valueScaled));
      }
      rowActivity_[i] = valueScaled*scaleR;
      double value = rowActivity_[i];
      if (value<rowLower_[i]-primalTolerance_)
	numberPrimalUnscaled++;
      else if (value>rowUpper_[i]+primalTolerance_)
	numberPrimalUnscaled++;
      double valueScaledDual = dual_[i]+rowObjectiveWork_[i];;
      if (valueScaled>rowLowerWork_[i]+primalTolerance_&&valueScaledDual>dualTolerance_)
	numberDualScaled++;
      if (valueScaled<rowUpperWork_[i]-primalTolerance_&&valueScaledDual<-dualTolerance_)
	numberDualScaled++;
      dual_[i] *= scaleC;
      double valueDual = dual_[i]; 
      if (rowObjective_)
	valueDual += rowObjective_[i];
      if (value>rowLower_[i]+primalTolerance_&&valueDual>dualTolerance_)
	numberDualUnscaled++;
      if (value<rowUpper_[i]-primalTolerance_&&valueDual<-dualTolerance_)
	numberDualUnscaled++;
    }
    if (!problemStatus_&&!secondaryStatus_) {
      // See if we need to set secondary status
      if (numberPrimalUnscaled) {
	if (numberDualUnscaled) 
	  secondaryStatus_=4;
	else
	  secondaryStatus_=2;
      } else {
	if (numberDualUnscaled) 
	  secondaryStatus_=3;
      }
    }
  } else {
    if (columnActivityWork_) {
      for (i=0;i<numberColumns_;i++) {
	double value = columnActivityWork_[i];
	double lower = columnLowerWork_[i];
	double upper = columnUpperWork_[i];
	if (lower>-1.0e20||upper<1.0e20) {
	  if (value>lower&&value<upper)
	    upperOut_ = CoinMax(upperOut_,CoinMin(value-lower,upper-value));
	}
	columnActivity_[i] = columnActivityWork_[i];
	reducedCost_[i] = reducedCostWork_[i];
      }
      for (i=0;i<numberRows_;i++) {
	double value = rowActivityWork_[i];
	double lower = rowLowerWork_[i];
	double upper = rowUpperWork_[i];
	if (lower>-1.0e20||upper<1.0e20) {
	  if (value>lower&&value<upper)
	    upperOut_ = CoinMax(upperOut_,CoinMin(value-lower,upper-value));
	}
	rowActivity_[i] = rowActivityWork_[i];
      }
    }
  }
  // switch off scalefactor if auto
  if (automaticScale_) {
    rhsScale_=1.0;
    objectiveScale_=1.0;
  }
  if (optimizationDirection_!=1.0) {
    // and modify all dual signs
    for (i=0;i<numberColumns_;i++) 
      reducedCost_[i] *= optimizationDirection_;
    for (i=0;i<numberRows_;i++) 
      dual_[i] *= optimizationDirection_;
  }
  // scaling may have been turned off
  scalingFlag_ = abs(scalingFlag_);
  if(getRidOfFactorizationData>0) {
    gutsOfDelete(getRidOfFactorizationData+1);
  } else {
    // at least get rid of nonLinearCost_
    delete nonLinearCost_;
    nonLinearCost_=NULL;
  }
  // get rid of data
  matrix_->generalExpanded(this,13,scalingFlag_);
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
void ClpSimplex::setNumberRefinements( int value) 
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
// Passes in factorization
void 
ClpSimplex::setFactorization( ClpFactorization & factorization)
{
  delete factorization_;
  factorization_= new ClpFactorization(factorization);
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
  if(value<=100&&value >=-1000) {
    perturbation_=value;
  } 
}
// Sparsity on or off
bool 
ClpSimplex::sparseFactorization() const
{
  return factorization_->sparseThreshold()!=0;
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
void checkCorrect(ClpSimplex * model,int iRow,
		  const double * element,const int * rowStart,const int * rowLength,
		  const int * column,
		  const double * columnLower_, const double * columnUpper_,
		  int infiniteUpperC,
		  int infiniteLowerC,
		  double &maximumUpC,
		  double &maximumDownC)
{
  int infiniteUpper = 0;
  int infiniteLower = 0;
  double maximumUp = 0.0;
  double maximumDown = 0.0;
  CoinBigIndex rStart = rowStart[iRow];
  CoinBigIndex rEnd = rowStart[iRow]+rowLength[iRow];
  CoinBigIndex j;
  double large=1.0e15;
  int iColumn;
  // Compute possible lower and upper ranges
  
  for (j = rStart; j < rEnd; ++j) {
    double value=element[j];
    iColumn = column[j];
    if (value > 0.0) {
      if (columnUpper_[iColumn] >= large) {
	++infiniteUpper;
      } else {
	maximumUp += columnUpper_[iColumn] * value;
      }
      if (columnLower_[iColumn] <= -large) {
	++infiniteLower;
      } else {
	maximumDown += columnLower_[iColumn] * value;
      }
    } else if (value<0.0) {
      if (columnUpper_[iColumn] >= large) {
	++infiniteLower;
      } else {
	maximumDown += columnUpper_[iColumn] * value;
      }
      if (columnLower_[iColumn] <= -large) {
	++infiniteUpper;
      } else {
	maximumUp += columnLower_[iColumn] * value;
      }
    }
  }
  assert (infiniteLowerC==infiniteLower);
  assert (infiniteUpperC==infiniteUpper);
  if (fabs(maximumUp-maximumUpC)>1.0e-12*CoinMax(fabs(maximumUp),fabs(maximumUpC)))
    printf("row %d comp up %g, true up %g\n",iRow,
	   maximumUpC,maximumUp);
  if (fabs(maximumDown-maximumDownC)>1.0e-12*CoinMax(fabs(maximumDown),fabs(maximumDownC)))
    printf("row %d comp down %g, true down %g\n",iRow,
	   maximumDownC,maximumDown);
  maximumUpC=maximumUp;
  maximumDownC=maximumDown;
}

/* Tightens primal bounds to make dual faster.  Unless
   fixed, bounds are slightly looser than they could be.
   This is to make dual go faster and is probably not needed
   with a presolve.  Returns non-zero if problem infeasible

   Fudge for branch and bound - put bounds on columns of factor *
   largest value (at continuous) - should improve stability
   in branch and bound on infeasible branches (0.0 is off)
*/
int 
ClpSimplex::tightenPrimalBounds(double factor,int doTight)
{
  
  // Get a row copy in standard format
  CoinPackedMatrix copy;
  copy.reverseOrderedCopyOf(*matrix());
  // get matrix data pointers
  const int * column = copy.getIndices();
  const CoinBigIndex * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  const double * element = copy.getElements();
  int numberChanged=1,iPass=0;
  double large = largeValue(); // treat bounds > this as infinite
#ifndef NDEBUG
  double large2= 1.0e10*large;
#endif
  int numberInfeasible=0;
  int totalTightened = 0;

  double tolerance = primalTolerance();


  // Save column bounds
  double * saveLower = new double [numberColumns_];
  memcpy(saveLower,columnLower_,numberColumns_*sizeof(double));
  double * saveUpper = new double [numberColumns_];
  memcpy(saveUpper,columnUpper_,numberColumns_*sizeof(double));

  int iRow, iColumn;

  // If wanted - tighten column bounds using solution
  if (factor) {
    double largest=0.0;
    if (factor>0.0) {
      assert (factor>1.0);
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
        if (columnUpper_[iColumn]-columnLower_[iColumn]>tolerance) {
          largest = CoinMax(largest,fabs(columnActivity_[iColumn]));
        }
      }
      largest *= factor;
    } else {
      // absolute
       largest = - factor;
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (columnUpper_[iColumn]-columnLower_[iColumn]>tolerance) {
        columnUpper_[iColumn] = CoinMin(columnUpper_[iColumn],largest);
        columnLower_[iColumn] = CoinMax(columnLower_[iColumn],-largest);
      }
    }
  }
#define MAXPASS 10

  // Loop round seeing if we can tighten bounds
  // Would be faster to have a stack of possible rows
  // and we put altered rows back on stack
  int numberCheck=-1;
  while(numberChanged>numberCheck) {

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
	CoinBigIndex rStart = rowStart[iRow];
	CoinBigIndex rEnd = rowStart[iRow]+rowLength[iRow];
	CoinBigIndex j;
	// Compute possible lower and upper ranges
      
	for (j = rStart; j < rEnd; ++j) {
	  double value=element[j];
	  iColumn = column[j];
	  if (value > 0.0) {
	    if (columnUpper_[iColumn] >= large) {
	      ++infiniteUpper;
	    } else {
	      maximumUp += columnUpper_[iColumn] * value;
	    }
	    if (columnLower_[iColumn] <= -large) {
	      ++infiniteLower;
	    } else {
	      maximumDown += columnLower_[iColumn] * value;
	    }
	  } else if (value<0.0) {
	    if (columnUpper_[iColumn] >= large) {
	      ++infiniteLower;
	    } else {
	      maximumDown += columnUpper_[iColumn] * value;
	    }
	    if (columnLower_[iColumn] <= -large) {
	      ++infiniteUpper;
	    } else {
	      maximumUp += columnLower_[iColumn] * value;
	    }
	  }
	}
	// Build in a margin of error
	maximumUp += 1.0e-8*fabs(maximumUp);
	maximumDown -= 1.0e-8*fabs(maximumDown);
	double maxUp = maximumUp+infiniteUpper*1.0e31;
	double maxDown = maximumDown-infiniteLower*1.0e31;
	if (maxUp <= rowUpper_[iRow] + tolerance && 
	    maxDown >= rowLower_[iRow] - tolerance) {
	  
	  // Row is redundant - make totally free
	  // NO - can't do this for postsolve
	  // rowLower_[iRow]=-COIN_DBL_MAX;
	  // rowUpper_[iRow]=COIN_DBL_MAX;
	  //printf("Redundant row in presolveX %d\n",iRow);

	} else {
	  if (maxUp < rowLower_[iRow] -100.0*tolerance ||
	      maxDown > rowUpper_[iRow]+100.0*tolerance) {
	    // problem is infeasible - exit at once
	    numberInfeasible++;
	    break;
	  }
	  double lower = rowLower_[iRow];
	  double upper = rowUpper_[iRow];
	  for (j = rStart; j < rEnd; ++j) {
	    double value=element[j];
	    iColumn = column[j];
	    double nowLower = columnLower_[iColumn];
	    double nowUpper = columnUpper_[iColumn];
	    if (value > 0.0) {
	      // positive value
	      if (lower>-large) {
		if (!infiniteUpper) {
		  assert(nowUpper < large2);
		  newBound = nowUpper + 
		    (lower - maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumUp);
		} else if (infiniteUpper==1&&nowUpper>large) {
		  newBound = (lower -maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumUp);
		} else {
		  newBound = -COIN_DBL_MAX;
		}
		if (newBound > nowLower + 1.0e-12&&newBound>-large) {
		  // Tighten the lower bound 
		  columnLower_[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (nowUpper - newBound < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowLower<-large) {
		    now=0.0;
		    infiniteLower--;
		  } else {
		    now = nowLower;
		  }
		  maximumDown += (newBound-now) * value;
		  nowLower = newBound;
#ifdef DEBUG
		  checkCorrect(this,iRow,
			       element, rowStart, rowLength,
			       column,
			       columnLower_,  columnUpper_,
			       infiniteUpper,
			       infiniteLower,
			       maximumUp,
			       maximumDown);
#endif
		}
	      } 
	      if (upper <large) {
		if (!infiniteLower) {
		  assert(nowLower >- large2);
		  newBound = nowLower + 
		    (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumDown);
		} else if (infiniteLower==1&&nowLower<-large) {
		  newBound =   (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumDown);
		} else {
		  newBound = COIN_DBL_MAX;
		}
		if (newBound < nowUpper - 1.0e-12&&newBound<large) {
		  // Tighten the upper bound 
		  columnUpper_[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (newBound - nowLower < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust 
		  double now;
		  if (nowUpper>large) {
		    now=0.0;
		    infiniteUpper--;
		  } else {
		    now = nowUpper;
		  }
		  maximumUp += (newBound-now) * value;
		  nowUpper = newBound;
#ifdef DEBUG
		  checkCorrect(this,iRow,
			       element, rowStart, rowLength,
			       column,
			       columnLower_,  columnUpper_,
			       infiniteUpper,
			       infiniteLower,
			       maximumUp,
			       maximumDown);
#endif
		}
	      }
	    } else {
	      // negative value
	      if (lower>-large) {
		if (!infiniteUpper) {
		  assert(nowLower < large2);
		  newBound = nowLower + 
		    (lower - maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumUp);
		} else if (infiniteUpper==1&&nowLower<-large) {
		  newBound = (lower -maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumUp);
		} else {
		  newBound = COIN_DBL_MAX;
		}
		if (newBound < nowUpper - 1.0e-12&&newBound<large) {
		  // Tighten the upper bound 
		  columnUpper_[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (newBound - nowLower < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowUpper>large) {
		    now=0.0;
		    infiniteLower--;
		  } else {
		    now = nowUpper;
		  }
		  maximumDown += (newBound-now) * value;
		  nowUpper = newBound;
#ifdef DEBUG
		  checkCorrect(this,iRow,
			       element, rowStart, rowLength,
			       column,
			       columnLower_,  columnUpper_,
			       infiniteUpper,
			       infiniteLower,
			       maximumUp,
			       maximumDown);
#endif
		}
	      }
	      if (upper <large) {
		if (!infiniteLower) {
		  assert(nowUpper < large2);
		  newBound = nowUpper + 
		    (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumDown);
		} else if (infiniteLower==1&&nowUpper>large) {
		  newBound =   (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumDown);
		} else {
		  newBound = -COIN_DBL_MAX;
		}
		if (newBound > nowLower + 1.0e-12&&newBound>-large) {
		  // Tighten the lower bound 
		  columnLower_[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (nowUpper - newBound < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowLower<-large) {
		    now=0.0;
		    infiniteUpper--;
		  } else {
		    now = nowLower;
		  }
		  maximumUp += (newBound-now) * value;
		  nowLower = newBound;
#ifdef DEBUG
		  checkCorrect(this,iRow,
			       element, rowStart, rowLength,
			       column,
			       columnLower_,  columnUpper_,
			       infiniteUpper,
			       infiniteLower,
			       maximumUp,
			       maximumDown);
#endif
		}
	      }
	    }
	  }
	}
      }
    }
    totalTightened += numberChanged;
    if (iPass==1)
      numberCheck=numberChanged>>4;
    if (numberInfeasible) break;
  }
  if (!numberInfeasible) {
    handler_->message(CLP_SIMPLEX_BOUNDTIGHTEN,messages_)
      <<totalTightened
      <<CoinMessageEol;
    // Set bounds slightly loose
    double useTolerance = 1.0e-3;
    if (doTight>0) {
      if (doTight>10) { 
        useTolerance=0.0;
      } else {
        while (doTight) {
          useTolerance *= 0.1;
          doTight--;
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (saveUpper[iColumn]>saveLower[iColumn]+useTolerance) {
        // Make large bounds stay infinite
        if (saveUpper[iColumn]>1.0e30&&columnUpper_[iColumn]>1.0e10) {
          columnUpper_[iColumn]=COIN_DBL_MAX;
        }
        if (saveLower[iColumn]<-1.0e30&&columnLower_[iColumn]<-1.0e10) {
          columnLower_[iColumn]=-COIN_DBL_MAX;
        }
	if (columnUpper_[iColumn]-columnLower_[iColumn]<useTolerance+1.0e-8) {
	  // relax enough so will have correct dj
#if 1
	  columnLower_[iColumn]=CoinMax(saveLower[iColumn],
				    columnLower_[iColumn]-100.0*useTolerance);
	  columnUpper_[iColumn]=CoinMin(saveUpper[iColumn],
				    columnUpper_[iColumn]+100.0*useTolerance);
#else
	  if (fabs(columnUpper_[iColumn])<fabs(columnLower_[iColumn])) {
	    if (columnUpper_[iColumn]- 100.0*useTolerance>saveLower[iColumn]) {
	      columnLower_[iColumn]=columnUpper_[iColumn]-100.0*useTolerance;
	    } else {
	      columnLower_[iColumn]=saveLower[iColumn];
	      columnUpper_[iColumn]=CoinMin(saveUpper[iColumn],
					saveLower[iColumn]+100.0*useTolerance);
	    }
	  } else {
	    if (columnLower_[iColumn]+100.0*useTolerance<saveUpper[iColumn]) {
	      columnUpper_[iColumn]=columnLower_[iColumn]+100.0*useTolerance;
	    } else {
	      columnUpper_[iColumn]=saveUpper[iColumn];
	      columnLower_[iColumn]=CoinMax(saveLower[iColumn],
					saveUpper[iColumn]-100.0*useTolerance);
	    }
	  }
#endif
	} else {
	  if (columnUpper_[iColumn]<saveUpper[iColumn]) {
	    // relax a bit
	    columnUpper_[iColumn] = CoinMin(columnUpper_[iColumn]+100.0*useTolerance,
					saveUpper[iColumn]);
	  }
	  if (columnLower_[iColumn]>saveLower[iColumn]) {
	    // relax a bit
	    columnLower_[iColumn] = CoinMax(columnLower_[iColumn]-100.0*useTolerance,
					saveLower[iColumn]);
	  }
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
#include "ClpSimplexPrimal.hpp"
int ClpSimplex::dual (int ifValuesPass , int startFinishOptions)
{
  //double savedPivotTolerance = factorization_->pivotTolerance();
  int saveQuadraticActivated = objective_->activated();
  objective_->setActivated(0);
  CoinAssert (ifValuesPass>=0&&ifValuesPass<3);
  /*  Note use of "down casting".  The only class the user sees is ClpSimplex.
      Classes ClpSimplexDual, ClpSimplexPrimal, (ClpSimplexNonlinear) 
      and ClpSimplexOther all exist and inherit from ClpSimplex but have no
      additional data and have no destructor or (non-default) constructor.

      This is to stop classes becoming too unwieldy and so I (JJF) can use e.g. "perturb"
      in primal and dual.

      As far as I can see this is perfectly safe.
  */
  int returnCode = ((ClpSimplexDual *) this)->dual(ifValuesPass, startFinishOptions);
  if ((specialOptions_&2048)!=0&&problemStatus_==10&&!numberPrimalInfeasibilities_
      &&sumDualInfeasibilities_<1000.0*dualTolerance_&&perturbation_>=100)
    problemStatus_=0; // ignore
  if (problemStatus_==10) {
    //printf("Cleaning up with primal\n");
    int savePerturbation = perturbation_;
    int saveLog = handler_->logLevel();
    //handler_->setLogLevel(63);
    perturbation_=100;
    bool denseFactorization = initialDenseFactorization();
    // It will be safe to allow dense
    setInitialDenseFactorization(true);
    // Allow for catastrophe
    int saveMax = intParam_[ClpMaxNumIteration];
    if (intParam_[ClpMaxNumIteration]>100000+numberIterations_)
      intParam_[ClpMaxNumIteration] = numberIterations_ + 1000 + 2*numberRows_+numberColumns_;
    // check which algorithms allowed
    int dummy;
    if ((matrix_->generalExpanded(this,4,dummy)&1)!=0)
      returnCode = ((ClpSimplexPrimal *) this)->primal(1,startFinishOptions);
    else
      returnCode = ((ClpSimplexDual *) this)->dual(0,startFinishOptions);
    if (problemStatus_==3&&numberIterations_<saveMax) {
      if (handler_->logLevel()==63)
	printf("looks like trouble - too many iterations in clean up - trying again\n");
      // flatten solution and try again
      int iRow,iColumn;
      for (iRow=0;iRow<numberRows_;iRow++) {
	if (getRowStatus(iRow)!=basic) {
	  setRowStatus(iRow,superBasic);
	  // but put to bound if close
	  if (fabs(rowActivity_[iRow]-rowLower_[iRow])
	      <=primalTolerance_) {
	    rowActivity_[iRow]=rowLower_[iRow];
	    setRowStatus(iRow,atLowerBound);
	  } else if (fabs(rowActivity_[iRow]-rowUpper_[iRow])
		     <=primalTolerance_) {
	    rowActivity_[iRow]=rowUpper_[iRow];
	    setRowStatus(iRow,atUpperBound);
	  }
	}
      }
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	if (getColumnStatus(iColumn)!=basic) {
	  setColumnStatus(iColumn,superBasic);
	  // but put to bound if close
	  if (fabs(columnActivity_[iColumn]-columnLower_[iColumn])
	      <=primalTolerance_) {
	    columnActivity_[iColumn]=columnLower_[iColumn];
	    setColumnStatus(iColumn,atLowerBound);
	  } else if (fabs(columnActivity_[iColumn]
			  -columnUpper_[iColumn])
		     <=primalTolerance_) {
	    columnActivity_[iColumn]=columnUpper_[iColumn];
	    setColumnStatus(iColumn,atUpperBound);
	  }
	}
      }
      problemStatus_=-1;
      intParam_[ClpMaxNumIteration] = CoinMin(numberIterations_ + 1000 + 
					  2*numberRows_+numberColumns_,saveMax);
      perturbation_=savePerturbation;
      returnCode = ((ClpSimplexPrimal *) this)->primal(0);
      if (problemStatus_==3&&numberIterations_<saveMax&& 
	  handler_->logLevel()>0)
	printf("looks like real trouble - too many iterations in second clean up - giving up\n");
    }
    intParam_[ClpMaxNumIteration] = saveMax;

    setInitialDenseFactorization(denseFactorization);
    perturbation_=savePerturbation;
    if (problemStatus_==10) 
      problemStatus_=0;
    handler_->setLogLevel(saveLog);
  }
  objective_->setActivated(saveQuadraticActivated);
  //factorization_->pivotTolerance(savedPivotTolerance);
  return returnCode;
}
#include "ClpQuadraticObjective.hpp"
// primal 
int ClpSimplex::primal (int ifValuesPass , int startFinishOptions)
{
  //double savedPivotTolerance = factorization_->pivotTolerance();
  // See if nonlinear
  if (objective_->type()>1&&objective_->activated()) 
    return reducedGradient();
  
  CoinAssert (ifValuesPass>=0&&ifValuesPass<3);
  /*  Note use of "down casting".  The only class the user sees is ClpSimplex.
      Classes ClpSimplexDual, ClpSimplexPrimal, (ClpSimplexNonlinear) 
      and ClpSimplexOther all exist and inherit from ClpSimplex but have no
      additional data and have no destructor or (non-default) constructor.

      This is to stop classes becoming too unwieldy and so I (JJF) can use e.g. "perturb"
      in primal and dual.

      As far as I can see this is perfectly safe.
  */
  int returnCode = ((ClpSimplexPrimal *) this)->primal(ifValuesPass,startFinishOptions);
  if (problemStatus_==10) {
    //printf("Cleaning up with dual\n");
    int savePerturbation = perturbation_;
    perturbation_=100;
    bool denseFactorization = initialDenseFactorization();
    // It will be safe to allow dense
    setInitialDenseFactorization(true);
    // check which algorithms allowed
    int dummy;
    if ((matrix_->generalExpanded(this,4,dummy)&2)!=0&&(specialOptions_&8192)==0) {
      double saveBound = dualBound_;
      // upperOut_ has largest away from bound
      dualBound_=CoinMin(2.0*upperOut_,dualBound_);
      returnCode = ((ClpSimplexDual *) this)->dual(0,startFinishOptions);
      dualBound_=saveBound;
    } else {
      returnCode = ((ClpSimplexPrimal *) this)->primal(0,startFinishOptions);
    }
    setInitialDenseFactorization(denseFactorization);
    perturbation_=savePerturbation;
    if (problemStatus_==10) 
      problemStatus_=0;
  }
  //factorization_->pivotTolerance(savedPivotTolerance);
  return returnCode;
}
/* Dual ranging.
   This computes increase/decrease in cost for each given variable and corresponding
   sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
   and numberColumns.. for artificials/slacks.
   For non-basic variables the sequence number will be that of the non-basic variables.
   
   Up to user to provide correct length arrays.
   
   Returns non-zero if infeasible unbounded etc
*/
#include "ClpSimplexOther.hpp"
int ClpSimplex::dualRanging(int numberCheck,const int * which,
			    double * costIncrease, int * sequenceIncrease,
			    double * costDecrease, int * sequenceDecrease)
{
  int savePerturbation = perturbation_;
  perturbation_=100;
  int returnCode = ((ClpSimplexPrimal *) this)->primal(0,1);
  if (problemStatus_==10) {
    //printf("Cleaning up with dual\n");
    bool denseFactorization = initialDenseFactorization();
    // It will be safe to allow dense
    setInitialDenseFactorization(true);
    // check which algorithms allowed
    int dummy;
    if ((matrix_->generalExpanded(this,4,dummy)&2)!=0) {
      // upperOut_ has largest away from bound
      double saveBound=dualBound_;
      if (upperOut_>0.0)
	dualBound_=2.0*upperOut_;
      returnCode = ((ClpSimplexDual *) this)->dual(0,1);
      dualBound_=saveBound;
    } else {
      returnCode = ((ClpSimplexPrimal *) this)->primal(0,1);
    }
    setInitialDenseFactorization(denseFactorization);
    if (problemStatus_==10) 
      problemStatus_=0;
  }
  perturbation_=savePerturbation;
  if (problemStatus_||secondaryStatus_==6) {
    finish(); // get rid of arrays
    return 1; // odd status
  }
  ((ClpSimplexOther *) this)->dualRanging(numberCheck,which,
					  costIncrease,sequenceIncrease,
					  costDecrease,sequenceDecrease);
  finish(); // get rid of arrays
  return 0;
}
/* Primal ranging.
   This computes increase/decrease in value for each given variable and corresponding
   sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
   and numberColumns.. for artificials/slacks.
   For basic variables the sequence number will be that of the basic variables.
   
   Up to user to providen correct length arrays.
   
   Returns non-zero if infeasible unbounded etc
*/
int ClpSimplex::primalRanging(int numberCheck,const int * which,
		  double * valueIncrease, int * sequenceIncrease,
		  double * valueDecrease, int * sequenceDecrease)
{
  int savePerturbation = perturbation_;
  perturbation_=100;
  int returnCode = ((ClpSimplexPrimal *) this)->primal(0,1);
  if (problemStatus_==10) {
    //printf("Cleaning up with dual\n");
    bool denseFactorization = initialDenseFactorization();
    // It will be safe to allow dense
    setInitialDenseFactorization(true);
    // check which algorithms allowed
    int dummy;
    if ((matrix_->generalExpanded(this,4,dummy)&2)!=0) {
      // upperOut_ has largest away from bound
      double saveBound=dualBound_;
      if (upperOut_>0.0)
	dualBound_=2.0*upperOut_;
      returnCode = ((ClpSimplexDual *) this)->dual(0,1);
      dualBound_=saveBound;
    } else {
      returnCode = ((ClpSimplexPrimal *) this)->primal(0,1);
    }
    setInitialDenseFactorization(denseFactorization);
    if (problemStatus_==10) 
      problemStatus_=0;
  }
  perturbation_=savePerturbation;
  if (problemStatus_||secondaryStatus_==6) {
    finish(); // get rid of arrays
    return 1; // odd status
  }
  ((ClpSimplexOther *) this)->primalRanging(numberCheck,which,
					  valueIncrease,sequenceIncrease,
					  valueDecrease,sequenceDecrease);
  finish(); // get rid of arrays
  return 0;
}
/* Write the basis in MPS format to the specified file.
   If writeValues true writes values of structurals
   (and adds VALUES to end of NAME card)
   
   Row and column names may be null.
   formatType is
   <ul>
   <li> 0 - normal
   <li> 1 - extra accuracy 
   <li> 2 - IEEE hex (later)
   </ul>
   
   Returns non-zero on I/O error
*/
int 
ClpSimplex::writeBasis(const char *filename,
			    bool writeValues,
			    int formatType) const
{
  return ((const ClpSimplexOther *) this)->writeBasis(filename,writeValues,
					 formatType);
}
// Read a basis from the given filename
int 
ClpSimplex::readBasis(const char *filename)
{
  return ((ClpSimplexOther *) this)->readBasis(filename);
}
#include "ClpSimplexNonlinear.hpp"
/* Solves nonlinear problem using SLP - may be used as crash
   for other algorithms when number of iterations small
*/
int 
ClpSimplex::nonlinearSLP(int numberPasses, double deltaTolerance)
{
  return ((ClpSimplexNonlinear *) this)->primalSLP(numberPasses,deltaTolerance);
}
// Solves non-linear using reduced gradient
int ClpSimplex::reducedGradient(int phase)
{
  if (objective_->type()<2||!objective_->activated()) {
    // no quadratic part
    return primal(0);
  }
  // get feasible
  if ((this->status()<0||numberPrimalInfeasibilities())&&phase==0) {
    objective_->setActivated(0);
    double saveDirection = optimizationDirection();
    setOptimizationDirection(0.0);
    primal(1);
    setOptimizationDirection(saveDirection);
    objective_->setActivated(1);
    // still infeasible
    if (numberPrimalInfeasibilities())
      return 0;
  }
  // Now enter method
  int returnCode = ((ClpSimplexNonlinear *) this)->primal();
  return returnCode;
}
#include "ClpPredictorCorrector.hpp"
#include "ClpCholeskyBase.hpp"
// Preference is WSSMP, TAUCS, UFL (just ordering) then base
#if WSSMP_BARRIER
#include "ClpCholeskyWssmp.hpp"
#include "ClpCholeskyWssmpKKT.hpp"
#elif UFL_BARRIER
#include "ClpCholeskyUfl.hpp"
#elif TAUCS_BARRIER
#include "ClpCholeskyTaucs.hpp"
#endif
#include "ClpPresolve.hpp"
/* Solves using barrier (assumes you have good cholesky factor code).
   Does crossover to simplex if asked*/
int 
ClpSimplex::barrier(bool crossover)
{
  ClpSimplex * model2 = this;
  int savePerturbation=perturbation_;
  ClpInterior barrier;
  barrier.borrowModel(*model2);
// Preference is WSSMP, UFL, TAUCS then base
#ifdef WSSMP_BARRIER
  ClpCholeskyWssmp * cholesky = new ClpCholeskyWssmp(CoinMax(100,model2->numberRows()/10));
  //ClpCholeskyWssmp * cholesky = new ClpCholeskyWssmp();
  //ClpCholeskyWssmpKKT * cholesky = new ClpCholeskyWssmpKKT(CoinMax(100,model2->numberRows()/10));
#elif UFL_BARRIER
  ClpCholeskyUfl * cholesky = new ClpCholeskyUfl();
#elif TAUCS_BARRIER
  ClpCholeskyTaucs * cholesky = new ClpCholeskyTaucs();
#else
  ClpCholeskyBase * cholesky = new ClpCholeskyBase();
#endif
  barrier.setCholesky(cholesky);
  int numberRows = model2->numberRows();
  int numberColumns = model2->numberColumns();
  int saveMaxIts = model2->maximumIterations();
  if (saveMaxIts<1000) {
    barrier.setMaximumBarrierIterations(saveMaxIts);
    model2->setMaximumIterations(1000000);
  }
  barrier.primalDual();
  int barrierStatus=barrier.status();
  double gap = barrier.complementarityGap();
  // get which variables are fixed
  double * saveLower=NULL;
  double * saveUpper=NULL;
  ClpPresolve pinfo2;
  ClpSimplex * saveModel2=NULL;
  int numberFixed = barrier.numberFixed();
  if (numberFixed*20>barrier.numberRows()&&numberFixed>5000&&crossover&&0) {
    // may as well do presolve
    int numberRows = barrier.numberRows();
    int numberColumns = barrier.numberColumns();
    int numberTotal = numberRows+numberColumns;
    saveLower = new double [numberTotal];
    saveUpper = new double [numberTotal];
    memcpy(saveLower,barrier.columnLower(),numberColumns*sizeof(double));
    memcpy(saveLower+numberColumns,barrier.rowLower(),numberRows*sizeof(double));
    memcpy(saveUpper,barrier.columnUpper(),numberColumns*sizeof(double));
    memcpy(saveUpper+numberColumns,barrier.rowUpper(),numberRows*sizeof(double));
    barrier.fixFixed();
    saveModel2=model2;
  }
  barrier.returnModel(*model2);
  double * rowPrimal = new double [numberRows];
  double * columnPrimal = new double [numberColumns];
  double * rowDual = new double [numberRows];
  double * columnDual = new double [numberColumns];
  // move solutions other way
  CoinMemcpyN(model2->primalRowSolution(),
	      numberRows,rowPrimal);
  CoinMemcpyN(model2->dualRowSolution(),
	      numberRows,rowDual);
  CoinMemcpyN(model2->primalColumnSolution(),
	      numberColumns,columnPrimal);
  CoinMemcpyN(model2->dualColumnSolution(),
	      numberColumns,columnDual);
  if (saveModel2) {
    // do presolve
    model2 = pinfo2.presolvedModel(*model2,1.0e-8,
				   false,5,true);
  }
  if (barrierStatus<4&&crossover) {
    // make sure no status left
    model2->createStatus();
    // solve
    model2->setPerturbation(100);
    // throw some into basis 
    {
      int numberRows = model2->numberRows();
      int numberColumns = model2->numberColumns();
      double * dsort = new double[numberColumns];
      int * sort = new int[numberColumns];
      int n=0;
      const double * columnLower = model2->columnLower();
      const double * columnUpper = model2->columnUpper();
      const double * primalSolution = model2->primalColumnSolution();
      double tolerance = 10.0*primalTolerance_;
      int i;
      for ( i=0;i<numberRows;i++) 
	model2->setRowStatus(i,superBasic);
      for ( i=0;i<numberColumns;i++) {
	double distance = CoinMin(columnUpper[i]-primalSolution[i],
			      primalSolution[i]-columnLower[i]);
	if (distance>tolerance) {
	  dsort[n]=-distance;
	  sort[n++]=i;
	  model2->setStatus(i,superBasic);
	} else if (distance>primalTolerance_) {
	  model2->setStatus(i,superBasic);
	} else if (primalSolution[i]<=columnLower[i]+primalTolerance_) {
	  model2->setStatus(i,atLowerBound);
	} else {
	  model2->setStatus(i,atUpperBound);
	}
      }
      CoinSort_2(dsort,dsort+n,sort);
      n = CoinMin(numberRows,n);
      for ( i=0;i<n;i++) {
	int iColumn = sort[i];
	model2->setStatus(iColumn,basic);
      }
      delete [] sort;
      delete [] dsort;
    }
    if (gap<1.0e-3*((double) (numberRows+numberColumns))) {
      int numberRows = model2->numberRows();
      int numberColumns = model2->numberColumns();
      // just primal values pass
      double saveScale = model2->objectiveScale();
      model2->setObjectiveScale(1.0e-3);
      model2->primal(2);
      model2->setObjectiveScale(saveScale);
      // save primal solution and copy back dual
      CoinMemcpyN(model2->primalRowSolution(),
		  numberRows,rowPrimal);
      CoinMemcpyN(rowDual,
		  numberRows,model2->dualRowSolution());
      CoinMemcpyN(model2->primalColumnSolution(),
		  numberColumns,columnPrimal);
      CoinMemcpyN(columnDual,
		  numberColumns,model2->dualColumnSolution());
      //model2->primal(1);
      // clean up reduced costs and flag variables
      {
	double * dj = model2->dualColumnSolution();
	double * cost = model2->objective();
	double * saveCost = new double[numberColumns];
	memcpy(saveCost,cost,numberColumns*sizeof(double));
	double * saveLower = new double[numberColumns];
	double * lower = model2->columnLower();
	memcpy(saveLower,lower,numberColumns*sizeof(double));
	double * saveUpper = new double[numberColumns];
	double * upper = model2->columnUpper();
	memcpy(saveUpper,upper,numberColumns*sizeof(double));
	int i;
	double tolerance = 10.0*dualTolerance_;
	for ( i=0;i<numberColumns;i++) {
	  if (model2->getStatus(i)==basic) {
	    dj[i]=0.0;
	  } else if (model2->getStatus(i)==atLowerBound) {
	    if (optimizationDirection_*dj[i]<tolerance) {
	      if (optimizationDirection_*dj[i]<0.0) {
		//if (dj[i]<-1.0e-3)
		//printf("bad dj at lb %d %g\n",i,dj[i]);
		cost[i] -= dj[i];
		dj[i]=0.0;
	      }
	    } else {
	      upper[i]=lower[i];
	    }
	  } else if (model2->getStatus(i)==atUpperBound) {
	    if (optimizationDirection_*dj[i]>tolerance) {
	      if (optimizationDirection_*dj[i]>0.0) {
		//if (dj[i]>1.0e-3)
		//printf("bad dj at ub %d %g\n",i,dj[i]);
		cost[i] -= dj[i];
		dj[i]=0.0;
	      }
	    } else {
	      lower[i]=upper[i];
	    }
	  }
	}
	// just dual values pass
	//model2->setLogLevel(63);
	//model2->setFactorizationFrequency(1);
	model2->dual(2);
	memcpy(cost,saveCost,numberColumns*sizeof(double));
	delete [] saveCost;
	memcpy(lower,saveLower,numberColumns*sizeof(double));
	delete [] saveLower;
	memcpy(upper,saveUpper,numberColumns*sizeof(double));
	delete [] saveUpper;
      }
      // and finish
      // move solutions
      CoinMemcpyN(rowPrimal,
		  numberRows,model2->primalRowSolution());
      CoinMemcpyN(columnPrimal,
		  numberColumns,model2->primalColumnSolution());
    }
//     double saveScale = model2->objectiveScale();
//     model2->setObjectiveScale(1.0e-3);
//     model2->primal(2);
//    model2->setObjectiveScale(saveScale);
    model2->primal(1);
  } else if (barrierStatus==4&&crossover) {
    // memory problems
    model2->setPerturbation(savePerturbation);
    model2->createStatus();
    model2->dual();
  }
  model2->setMaximumIterations(saveMaxIts);
  delete [] rowPrimal;
  delete [] columnPrimal;
  delete [] rowDual;
  delete [] columnDual;
  if (saveLower) {
    pinfo2.postsolve(true);
    delete model2;
    model2=saveModel2;
    int numberRows = model2->numberRows();
    int numberColumns = model2->numberColumns();
    memcpy(model2->columnLower(),saveLower,numberColumns*sizeof(double));
    memcpy(model2->rowLower(),saveLower+numberColumns,numberRows*sizeof(double));
    delete [] saveLower;
    memcpy(model2->columnUpper(),saveUpper,numberColumns*sizeof(double));
    memcpy(model2->rowUpper(),saveUpper+numberColumns,numberRows*sizeof(double));
    delete [] saveUpper;
    model2->primal(1);
  }
  model2->setPerturbation(savePerturbation);
  return model2->status();
}
/* For strong branching.  On input lower and upper are new bounds
   while on output they are objective function values (>1.0e50 infeasible).
   Return code is 0 if nothing interesting, -1 if infeasible both
   ways and +1 if infeasible one way (check values to see which one(s))
*/
int ClpSimplex::strongBranching(int numberVariables,const int * variables,
				double * newLower, double * newUpper,
				double ** outputSolution,
				int * outputStatus, int * outputIterations,
				bool stopOnFirstInfeasible,
				bool alwaysFinish,
				int startFinishOptions)
{
  return ((ClpSimplexDual *) this)->strongBranching(numberVariables,variables,
						    newLower,  newUpper,outputSolution,
						    outputStatus, outputIterations,
						    stopOnFirstInfeasible,
						    alwaysFinish,startFinishOptions);
}
/* Borrow model.  This is so we dont have to copy large amounts
   of data around.  It assumes a derived class wants to overwrite
   an empty model with a real one - while it does an algorithm.
   This is same as ClpModel one, but sets scaling on etc. */
void 
ClpSimplex::borrowModel(ClpModel & otherModel) 
{
  ClpModel::borrowModel(otherModel);
  createStatus();
  //ClpDualRowSteepest steep1;
  //setDualRowPivotAlgorithm(steep1);
  //ClpPrimalColumnSteepest steep2;
  //setPrimalColumnPivotAlgorithm(steep2);
}
void 
ClpSimplex::borrowModel(ClpSimplex & otherModel) 
{
  ClpModel::borrowModel(otherModel);
  createStatus();
  dualBound_ = otherModel.dualBound_;
  dualTolerance_ = otherModel.dualTolerance_;
  primalTolerance_ = otherModel.primalTolerance_;
  delete dualRowPivot_;
  dualRowPivot_ = otherModel.dualRowPivot_->clone(false);
  delete primalColumnPivot_;
  primalColumnPivot_ = otherModel.primalColumnPivot_->clone(false);
  perturbation_ = otherModel.perturbation_;
  specialOptions_ = otherModel.specialOptions_;
  automaticScale_ = otherModel.automaticScale_;
}
/// Saves scalars for ClpSimplex
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
    int i;
    CoinBigIndex numberWritten;
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
    scalars.maximumIterations = maximumIterations();
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
    CoinBigIndex length;
    for (i=0;i<ClpLastStrParam;i++) {
      length = strParam_[i].size();
      numberWritten = fwrite(&length,sizeof(int),1,fp);
      if (numberWritten!=1)
	return 1;
      if (length) {
	numberWritten = fwrite(strParam_[i].c_str(),length,1,fp);
	if (numberWritten!=1)
	  return 1;
      }
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
    if (outDoubleArray(objective(),numberColumns_,fp))
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
      char * array = 
	new char[CoinMax(numberRows_,numberColumns_)*(lengthNames_+1)];
      char * put = array;
      CoinAssert (numberRows_ == (int) rowNames_.size());
      for (i=0;i<numberRows_;i++) {
	assert((int)rowNames_[i].size()<=lengthNames_);
	strcpy(put,rowNames_[i].c_str());
	put += lengthNames_+1;
      }
      numberWritten = fwrite(array,lengthNames_+1,numberRows_,fp);
      if (numberWritten!=numberRows_)
	return 1;
      put=array;
      CoinAssert (numberColumns_ == (int) columnNames_.size());
      for (i=0;i<numberColumns_;i++) {
	assert((int) columnNames_[i].size()<=lengthNames_);
	strcpy(put,columnNames_[i].c_str());
	put += lengthNames_+1;
      }
      numberWritten = fwrite(array,lengthNames_+1,numberColumns_,fp);
      if (numberWritten!=numberColumns_)
	return 1;
      delete [] array;
    }
    // just standard type at present
    assert (matrix_->type()==1);
    CoinAssert (matrix_->getNumCols() == numberColumns_);
    CoinAssert (matrix_->getNumRows() == numberRows_);
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
			   sizeof(int),numberColumns_+1,fp);
    if (numberWritten!=numberColumns_+1)
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
    // save event handler in case already set
    ClpEventHandler * handler = eventHandler_->clone();
    ClpModel::gutsOfDelete();
    eventHandler_ = handler;
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
    CoinBigIndex numberRead;

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
    setMaximumIterations(scalars.maximumIterations);
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
    CoinBigIndex length;
    for (i=0;i<ClpLastStrParam;i++) {
      numberRead = fread(&length,sizeof(int),1,fp);
      if (numberRead!=1)
	return 1;
      if (length) {
	char * array = new char[length+1];
	numberRead = fread(array,length,1,fp);
	if (numberRead!=1)
	  return 1;
	array[length]='\0';
	strParam_[i]=array;
	delete [] array;
      }
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
    double * objective;
    if (inDoubleArray(objective,numberColumns_,fp))
	return 1;
    delete objective_;
    objective_ = new ClpLinearObjective(objective,numberColumns_);
    delete [] objective;
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
      char * array = 
	new char[CoinMax(numberRows_,numberColumns_)*(lengthNames_+1)];
      char * get = array;
      numberRead = fread(array,lengthNames_+1,numberRows_,fp);
      if (numberRead!=numberRows_)
	return 1;
      rowNames_ = std::vector<std::string> ();
      rowNames_.resize(numberRows_);
      for (i=0;i<numberRows_;i++) {
	rowNames_.push_back(get);
	get += lengthNames_+1;
      }
      get = array;
      numberRead = fread(array,lengthNames_+1,numberColumns_,fp);
      if (numberRead!=numberColumns_)
	return 1;
      columnNames_ = std::vector<std::string> ();
      columnNames_.resize(numberColumns_);
      for (i=0;i<numberColumns_;i++) {
	columnNames_.push_back(get);
	get += lengthNames_+1;
      }
      delete [] array;
    }
    // Pivot choices
    assert(scalars.dualPivotChoice>0&&(scalars.dualPivotChoice&63)<3);
    delete dualRowPivot_;
    switch ((scalars.dualPivotChoice&63)) {
    default:
      printf("Need another dualPivot case %d\n",scalars.dualPivotChoice&63);
    case 1:
      // Dantzig
      dualRowPivot_ = new ClpDualRowDantzig();
      break;
    case 2:
      // Steepest - use mode
      dualRowPivot_ = new ClpDualRowSteepest(scalars.dualPivotChoice>>6);
      break;
    }
    assert(scalars.primalPivotChoice>0&&(scalars.primalPivotChoice&63)<3);
    delete primalColumnPivot_;
    switch ((scalars.primalPivotChoice&63)) {
    default:
      printf("Need another primalPivot case %d\n",
	     scalars.primalPivotChoice&63);
    case 1:
      // Dantzig
      primalColumnPivot_ = new ClpPrimalColumnDantzig();
      break;
    case 2:
      // Steepest - use mode
      primalColumnPivot_ 
	= new ClpPrimalColumnSteepest(scalars.primalPivotChoice>>6);
      break;
    }
    assert(scalars.matrixStorageChoice==1);
    delete matrix_;
    // get arrays
    numberRead = fread(&length,sizeof(int),1,fp);
    if (numberRead!=1)
      return 1;
    double * elements = new double[length];
    int * indices = new int[length];
    CoinBigIndex * starts = new CoinBigIndex[numberColumns_+1];
    int * lengths = new int[numberColumns_];
    numberRead = fread(elements, sizeof(double),length,fp);
    if (numberRead!=length)
      return 1;
    numberRead = fread(indices, sizeof(int),length,fp);
    if (numberRead!=length)
      return 1;
    numberRead = fread(starts, sizeof(int),numberColumns_+1,fp);
    if (numberRead!=numberColumns_+1)
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
// Sanity check on input data - returns true if okay
bool 
ClpSimplex::sanityCheck()
{
  // bad if empty
  if (!numberColumns_||((!numberRows_||!matrix_->getNumElements())&&objective_->type()<2)) {
    int infeasNumber[2];
    double infeasSum[2];
    problemStatus_=emptyProblem(infeasNumber,infeasSum,false);
    numberDualInfeasibilities_=infeasNumber[0];
    sumDualInfeasibilities_=infeasSum[0];
    numberPrimalInfeasibilities_=infeasNumber[1];
    sumPrimalInfeasibilities_=infeasSum[1];
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
  double fixTolerance = 10.0*primalTolerance_;
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
    if (value<-primalTolerance_) {
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
    if (value<-primalTolerance_) {
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
    <<CoinMessageEol; 
  if (largestBound)
    handler_->message(CLP_RIMSTATISTICS2,messages_)
      <<smallestBound
      <<largestBound
      <<minimumGap
      <<CoinMessageEol;
  return true;
}
// Set up status array (for OsiClp)
void 
ClpSimplex::createStatus() 
{
  if(!status_)
    status_ = new unsigned char [numberColumns_+numberRows_];
  memset(status_,0,(numberColumns_+numberRows_)*sizeof(char));
  int i;
  // set column status to one nearest zero
  for (i=0;i<numberColumns_;i++) {
#if 0
    if (columnLower_[i]>=0.0) {
      setColumnStatus(i,atLowerBound);
    } else if (columnUpper_[i]<=0.0) {
      setColumnStatus(i,atUpperBound);
    } else if (columnLower_[i]<-1.0e20&&columnUpper_[i]>1.0e20) {
      // free
      setColumnStatus(i,isFree);
    } else if (fabs(columnLower_[i])<fabs(columnUpper_[i])) {
      setColumnStatus(i,atLowerBound);
    } else {
      setColumnStatus(i,atUpperBound);
    }
#else
    setColumnStatus(i,atLowerBound);
#endif
  }
  for (i=0;i<numberRows_;i++) {
    setRowStatus(i,basic);
  }
}
/* Sets up all slack basis and resets solution to 
   as it was after initial load or readMps */
void ClpSimplex::allSlackBasis(bool resetSolution)
{
  createStatus();
  if (resetSolution) {
    // put back to as it was originally
    int i;
    // set column status to one nearest zero
    // But set value to zero if lb <0.0 and ub>0.0
    for (i=0;i<numberColumns_;i++) {
      if (columnLower_[i]>=0.0) {
	columnActivity_[i]=columnLower_[i];
	setColumnStatus(i,atLowerBound);
      } else if (columnUpper_[i]<=0.0) {
	columnActivity_[i]=columnUpper_[i];
	setColumnStatus(i,atUpperBound);
      } else if (columnLower_[i]<-1.0e20&&columnUpper_[i]>1.0e20) {
	// free
	columnActivity_[i]=0.0;
	setColumnStatus(i,isFree);
      } else if (fabs(columnLower_[i])<fabs(columnUpper_[i])) {
	columnActivity_[i]=0.0;
	setColumnStatus(i,atLowerBound);
      } else {
	columnActivity_[i]=0.0;
	setColumnStatus(i,atUpperBound);
      }
    }
  }
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
ClpSimplex::loadProblem (  const ClpMatrixBase& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
  createStatus();
}
void 
ClpSimplex::loadProblem (  const CoinPackedMatrix& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
  createStatus();
}

/* Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). */
void 
ClpSimplex::loadProblem (  const int numcols, const int numrows,
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
  createStatus();
}
// This loads a model from a coinModel object - returns number of errors
int 
ClpSimplex::loadProblem (  CoinModel & modelObject, bool keepSolution)
{
  unsigned char * status = NULL;
  double * psol = NULL;
  double * dsol = NULL;
  if (status_&&numberRows_&&numberRows_==modelObject.numberRows()&&
      numberColumns_==modelObject.numberColumns()) {
    status = new unsigned char [numberRows_+numberColumns_];
    memcpy(status,status_,numberRows_+numberColumns_);
    psol = new double [numberRows_+numberColumns_];
    memcpy(psol,columnActivity_,numberColumns_*sizeof(double));
    memcpy(psol+numberColumns_,rowActivity_,numberRows_*sizeof(double));
    dsol = new double [numberRows_+numberColumns_];
    memcpy(dsol,reducedCost_,numberColumns_*sizeof(double));
    memcpy(dsol+numberColumns_,dual_,numberRows_*sizeof(double));
  }
  int returnCode = ClpModel::loadProblem(modelObject);
  createStatus();
  if (status) {
    // copy back
    memcpy(status_,status,numberRows_+numberColumns_);
    memcpy(columnActivity_,psol,numberColumns_*sizeof(double));
    memcpy(rowActivity_,psol+numberColumns_,numberRows_*sizeof(double));
    memcpy(reducedCost_,dsol,numberColumns_*sizeof(double));
    memcpy(dual_,dsol+numberColumns_,numberRows_*sizeof(double));
    delete [] status;
    delete [] psol;
    delete [] dsol;
  }
  return returnCode;
}
void 
ClpSimplex::loadProblem (  const int numcols, const int numrows,
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
  createStatus();
}
// Read an mps file from the given filename
int 
ClpSimplex::readMps(const char *filename,
	    bool keepNames,
	    bool ignoreErrors)
{
  int status = ClpModel::readMps(filename,keepNames,ignoreErrors);
  createStatus();
  return status;
}
// Read GMPL files from the given filenames
int 
ClpSimplex::readGMPL(const char *filename,const char * dataName,
                     bool keepNames)
{
  int status = ClpModel::readGMPL(filename,dataName,keepNames);
  createStatus();
  return status;
}
// Just check solution (for external use)
void 
ClpSimplex::checkSolution()
{
  // put in standard form
  createRim(7+8+16+32);
  dualTolerance_=dblParam_[ClpDualTolerance];
  primalTolerance_=dblParam_[ClpPrimalTolerance];
  checkPrimalSolution( rowActivityWork_, columnActivityWork_);
  checkDualSolution();
  if (!numberDualInfeasibilities_&&
      !numberPrimalInfeasibilities_)
    problemStatus_=0;
  else
    problemStatus_=-1;
#ifdef CLP_DEBUG
  int i;
  double value=0.0;
  for (i=0;i<numberRows_+numberColumns_;i++)
    value += dj_[i]*solution_[i];
  printf("dual value %g, primal %g\n",value,objectiveValue());
#endif
  // release extra memory
  deleteRim(0);
}
/* Crash - at present just aimed at dual, returns
   -2 if dual preferred and crash basis created
   -1 if dual preferred and all slack basis preferred
   0 if basis going in was not all slack
   1 if primal preferred and all slack basis preferred
   2 if primal preferred and crash basis created.
   
   if gap between bounds <="gap" variables can be flipped
   ( If pivot -1 then can be made super basic!)
   
   If "pivot" is
   -1 No pivoting - always primal
   0 No pivoting (so will just be choice of algorithm)
   1 Simple pivoting e.g. gub
   2 Mini iterations
   3 Just throw all free variables in basis
*/
int 
ClpSimplex::crash(double gap,int pivot)
{
  //CoinAssert(!rowObjective_); // not coded
  int iColumn;
  int numberBad=0;
  int numberBasic=0;
  double dualTolerance=dblParam_[ClpDualTolerance];
  //double primalTolerance=dblParam_[ClpPrimalTolerance];
  int returnCode=0;
  // If no basis then make all slack one
  if (!status_)
    createStatus();

  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (getColumnStatus(iColumn)==basic)
      numberBasic++;
  }
  if (!numberBasic||pivot==3) {
    if (pivot==3) {
      // Just throw all free variables in basis
      for (int iRow=0;iRow<numberRows_;iRow++) {
        if (fabs(rowLower_[iRow])<fabs(rowUpper_[iRow]))
          setRowStatus(iRow,atLowerBound);
        else
          setRowStatus(iRow,atUpperBound);
      }
      for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
        if (columnLower_[iColumn]<-1.0e20&&columnUpper_[iColumn]>1.0e20)
          setColumnStatus(iColumn,basic);
      }
      return 0;
    }
    // all slack
    double * dj = new double [numberColumns_];
    double * solution = columnActivity_;
    const double * linearObjective = objective();
    //double objectiveValue=0.0;
    int iColumn;
    double direction = optimizationDirection_;
    // direction is actually scale out not scale in
    if (direction)
      direction = 1.0/direction;
    for (iColumn=0;iColumn<numberColumns_;iColumn++)
      dj[iColumn] = direction*linearObjective[iColumn];
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      // assume natural place is closest to zero
      double lowerBound = columnLower_[iColumn];
      double upperBound = columnUpper_[iColumn];
      if (lowerBound>-1.0e20||upperBound<1.0e20) {
	bool atLower;
	if (fabs(upperBound)<fabs(lowerBound)) {
	  atLower=false;
	  setColumnStatus(iColumn,atUpperBound);
	  solution[iColumn]=upperBound;
	} else {
	  atLower=true;
	  setColumnStatus(iColumn,atLowerBound);
	  solution[iColumn]=lowerBound;
	}
	if (dj[iColumn]<-dualTolerance_) {
	  // should be at upper bound
	  if (atLower) {
	    // can we flip
	    if (upperBound-lowerBound<=gap) {
	      columnActivity_[iColumn]=upperBound;
	      setColumnStatus(iColumn,atUpperBound);
	    } else if (pivot<0) {
	      // set superbasic
	      columnActivity_[iColumn]=lowerBound+gap;
	      setColumnStatus(iColumn,superBasic);
	    } else if (dj[iColumn]<-dualTolerance) {
	      numberBad++;
	    }
	  }
	} else if (dj[iColumn]>dualTolerance_) {
	  // should be at lower bound
	  if (!atLower) {
	    // can we flip
	    if (upperBound-lowerBound<=gap) {
	      columnActivity_[iColumn]=lowerBound;
	      setColumnStatus(iColumn,atLowerBound);
	    } else if (pivot<0) {
	      // set superbasic
	      columnActivity_[iColumn]=upperBound-gap;
	      setColumnStatus(iColumn,superBasic);
	    } else if (dj[iColumn]>dualTolerance) {
	      numberBad++;
	    }
	  }
	}
      } else {
	// free
	setColumnStatus(iColumn,isFree);
	if (fabs(dj[iColumn])>dualTolerance) 
	  numberBad++;
      }
    }
    if (numberBad||pivot) {
      if (pivot<=0) {
	delete [] dj;
	returnCode = 1;
      } else {
	// see if can be made dual feasible with gubs etc
	double * pi = new double[numberRows_];
	memset (pi,0,numberRows_*sizeof(double));
	int * way = new int[numberColumns_];
	int numberIn = 0;

	// Get column copy
	CoinPackedMatrix * columnCopy = matrix();
	// Get a row copy in standard format
	CoinPackedMatrix copy;
	copy.reverseOrderedCopyOf(*columnCopy);
	// get matrix data pointers
	const int * column = copy.getIndices();
	const CoinBigIndex * rowStart = copy.getVectorStarts();
	const int * rowLength = copy.getVectorLengths(); 
	const double * elementByRow = copy.getElements();
	//const int * row = columnCopy->getIndices();
	//const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
	//const int * columnLength = columnCopy->getVectorLengths(); 
	//const double * element = columnCopy->getElements();


	// if equality row and bounds mean artificial in basis bad
	// then do anyway

	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  // - if we want to reduce dj, + if we want to increase
	  int thisWay = 100;
	  double lowerBound = columnLower_[iColumn];
	  double upperBound = columnUpper_[iColumn];
	  if (upperBound>lowerBound) {
	    switch(getColumnStatus(iColumn)) {
	      
	    case basic:
	      thisWay=0;
	    case ClpSimplex::isFixed:
	      break;
	    case isFree:
	    case superBasic:
	      if (dj[iColumn]<-dualTolerance) 
		thisWay = 1;
	      else if (dj[iColumn]>dualTolerance) 
		thisWay = -1;
	      else
		thisWay =0;
	      break;
	    case atUpperBound:
	      if (dj[iColumn]>dualTolerance) 
		thisWay = -1;
	      else if (dj[iColumn]<-dualTolerance) 
		thisWay = -3;
	      else
		thisWay = -2;
	      break;
	    case atLowerBound:
	      if (dj[iColumn]<-dualTolerance) 
		thisWay = 1;
	      else if (dj[iColumn]>dualTolerance) 
		thisWay = 3;
	      else
		thisWay = 2;
	      break;
	    }
	  }
	  way[iColumn] = thisWay;
	}
	/*if (!numberBad)
	  printf("Was dual feasible before passes - rows %d\n",
	  numberRows_);*/
	int lastNumberIn = -100000;
	int numberPasses=5;
	while (numberIn>lastNumberIn+numberRows_/100) {
	  lastNumberIn = numberIn;
	  // we need to maximize chance of doing good
	  int iRow;
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    double lowerBound = rowLower_[iRow];
	    double upperBound = rowUpper_[iRow];
	    if (getRowStatus(iRow)==basic) {
	      // see if we can find a column to pivot on
	      int j;
	      // down is amount pi can go down
	      double maximumDown = COIN_DBL_MAX;
	      double maximumUp = COIN_DBL_MAX;
	      double minimumDown =0.0;
	      double minimumUp =0.0;
	      int iUp=-1;
	      int iDown=-1;
	      int iUpB=-1;
	      int iDownB=-1;
	      if (lowerBound<-1.0e20)
		maximumUp = -1.0;
	      if (upperBound>1.0e20)
		maximumDown = -1.0;
	      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		int iColumn = column[j];
		double value = elementByRow[j];
		double djValue = dj[iColumn];
		/* way -
		   -3 - okay at upper bound with negative dj
		   -2 - marginal at upper bound with zero dj - can only decrease
		   -1 - bad at upper bound
		   0 - we can never pivot on this row
		   1 - bad at lower bound
		   2 - marginal at lower bound with zero dj - can only increase
		   3 - okay at lower bound with positive dj
		   100 - fine we can just ignore
		*/
		if (way[iColumn]!=100) {
		  switch(way[iColumn]) {
		    
		  case -3:
		    if (value>0.0) {
		      if (maximumDown*value>-djValue) {
			maximumDown = -djValue/value;
			iDown = iColumn;
		      }
		    } else {
		      if (-maximumUp*value>-djValue) {
			maximumUp = djValue/value;
			iUp = iColumn;
		      }
		    }
		    break;
		  case -2:
		    if (value>0.0) {
		      maximumDown = 0.0;
		    } else {
		      maximumUp = 0.0;
		    }
		    break;
		  case -1:
		    // see if could be satisfied
		    // dj value > 0
		    if (value>0.0) {
		      maximumDown=0.0;
		      if (maximumUp*value<djValue-dualTolerance) {
			maximumUp = 0.0; // would improve but not enough
		      } else {
			if (minimumUp*value<djValue) {
			  minimumUp = djValue/value;
			  iUpB = iColumn;
			}
		      }
		    } else {
		      maximumUp=0.0;
		      if (-maximumDown*value<djValue-dualTolerance) {
			maximumDown = 0.0; // would improve but not enough
		      } else {
			if (-minimumDown*value<djValue) {
			  minimumDown = -djValue/value;
			  iDownB = iColumn;
			}
		      }
		    }
		    
		    break;
		  case 0:
		    maximumDown = -1.0;
		    maximumUp=-1.0;
		    break;
		  case 1:
		    // see if could be satisfied
		    // dj value < 0
		    if (value>0.0) {
		      maximumUp=0.0;
		      if (maximumDown*value<-djValue-dualTolerance) {
			maximumDown = 0.0; // would improve but not enough
		      } else {
			if (minimumDown*value<-djValue) {
			  minimumDown = -djValue/value;
			  iDownB = iColumn;
			}
		      }
		    } else {
		      maximumDown=0.0;
		      if (-maximumUp*value<-djValue-dualTolerance) {
			maximumUp = 0.0; // would improve but not enough
		      } else {
			if (-minimumUp*value<-djValue) {
			  minimumUp = djValue/value;
			  iUpB = iColumn;
			}
		      }
		    }
		    
		    break;
		  case 2:
		    if (value>0.0) {
		      maximumUp = 0.0;
		    } else {
		      maximumDown = 0.0;
		    }
		    
		    break;
		  case 3:
		    if (value>0.0) {
		      if (maximumUp*value>djValue) {
			maximumUp = djValue/value;
			iUp = iColumn;
		      }
		    } else {
		      if (-maximumDown*value>djValue) {
			maximumDown = -djValue/value;
			iDown = iColumn;
		      }
		    }
		    
		    break;
		  default:
		    break;
		  }
		}
	      }
	      if (iUpB>=0)
		iUp=iUpB;
	      if (maximumUp<=dualTolerance||maximumUp<minimumUp)
		iUp=-1;
	      if (iDownB>=0)
		iDown=iDownB;
	      if (maximumDown<=dualTolerance||maximumDown<minimumDown)
		iDown=-1;
	      if (iUp>=0||iDown>=0) {
		// do something
		if (iUp>=0&&iDown>=0) {
		  if (maximumDown>maximumUp)
		    iUp=-1;
		}
		double change;
		int kColumn;
		if (iUp>=0) {
		  kColumn=iUp;
		  change=maximumUp;
		  // just do minimum if was dual infeasible
		  // ? only if maximum large?
		  if (minimumUp>0.0)
		    change=minimumUp;
		  setRowStatus(iRow,atUpperBound);
		} else {
		  kColumn=iDown;
		  change=-maximumDown;
		  // just do minimum if was dual infeasible
		  // ? only if maximum large?
		  if (minimumDown>0.0)
		    change=-minimumDown;
		  setRowStatus(iRow,atLowerBound);
		}
		assert (fabs(change)<1.0e20);
		setColumnStatus(kColumn,basic);
		numberIn++;
		pi[iRow]=change;
		for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		  int iColumn = column[j];
		  double value = elementByRow[j];
		  double djValue = dj[iColumn]-change*value;
		  dj[iColumn]=djValue;
		  if (abs(way[iColumn])==1) {
		    numberBad--;
		    /*if (!numberBad)
		      printf("Became dual feasible at row %d out of %d\n",
		      iRow, numberRows_);*/
		    lastNumberIn=-1000000;
		  }
		  int thisWay = 100;
		  double lowerBound = columnLower_[iColumn];
		  double upperBound = columnUpper_[iColumn];
		  if (upperBound>lowerBound) {
		    switch(getColumnStatus(iColumn)) {
		      
		    case basic:
		      thisWay=0;
		    case isFixed:
		      break;
		    case isFree:
		    case superBasic:
		      if (djValue<-dualTolerance) 
			thisWay = 1;
		      else if (djValue>dualTolerance) 
			thisWay = -1;
		      else
			{ thisWay =0; abort();}
		      break;
		    case atUpperBound:
		      if (djValue>dualTolerance) 
			{ thisWay =-1; abort();}
		      else if (djValue<-dualTolerance) 
			thisWay = -3;
		      else
			thisWay = -2;
		      break;
		    case atLowerBound:
		      if (djValue<-dualTolerance) 
			{ thisWay =1; abort();}
		      else if (djValue>dualTolerance) 
			thisWay = 3;
		      else
			thisWay = 2;
		      break;
		    }
		  }
		  way[iColumn] = thisWay;
		}
	      }
	    }
	  }
	  if (numberIn==lastNumberIn||numberBad||pivot<2)
	    break;
	  if (!(--numberPasses))
	    break;
	  //printf("%d put in so far\n",numberIn);
	}
	// last attempt to flip
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  double lowerBound = columnLower_[iColumn];
	  double upperBound = columnUpper_[iColumn];
	  if (upperBound-lowerBound<=gap&&upperBound>lowerBound) {
	    double djValue=dj[iColumn];
	    switch(getColumnStatus(iColumn)) {
	      
	    case basic:
	    case ClpSimplex::isFixed:
	      break;
	    case isFree:
	    case superBasic:
	      break;
	    case atUpperBound:
	      if (djValue>dualTolerance) {
		setColumnStatus(iColumn,atUpperBound);
		solution[iColumn]=upperBound;
	      } 
	      break;
	    case atLowerBound:
	      if (djValue<-dualTolerance) {
		setColumnStatus(iColumn,atUpperBound);
		solution[iColumn]=upperBound;
	      }
	      break;
	    }
	  }
	}
	delete [] pi;
	delete [] dj;
	delete [] way;
	handler_->message(CLP_CRASH,messages_)
	  <<numberIn
	  <<numberBad
	  <<CoinMessageEol;
	returnCode =  -1;
      }
    } else {
      delete [] dj;
      returnCode =  -1;
    }
    //cleanStatus();
  }
  return returnCode;
}
/* Pivot in a variable and out a variable.  Returns 0 if okay,
   1 if inaccuracy forced re-factorization, -1 if would be singular.
   Also updates primal/dual infeasibilities.
   Assumes sequenceIn_ and pivotRow_ set and also directionIn and Out.
*/
int ClpSimplex::pivot()
{
  // scaling not allowed
  assert (!scalingFlag_);
  // assume In_ and Out_ are correct and directionOut_ set
  // (or In_ if flip
  lowerIn_ = lower_[sequenceIn_];
  valueIn_ = solution_[sequenceIn_];
  upperIn_ = upper_[sequenceIn_];
  dualIn_ = dj_[sequenceIn_];
  lowerOut_ = lower_[sequenceOut_];
  valueOut_ = solution_[sequenceOut_];
  upperOut_ = upper_[sequenceOut_];
  // for now assume primal is feasible (or in dual)
  dualOut_ = dj_[sequenceOut_];
  assert(fabs(dualOut_)<1.0e-6);
  bool roundAgain = true;
  int returnCode=0;
  while (roundAgain) {
    roundAgain=false;
    unpack(rowArray_[1]);
    factorization_->updateColumnFT(rowArray_[2],rowArray_[1]);
    alpha_=0.0;
    int i;
    int * index = rowArray_[1]->getIndices();
    int number = rowArray_[1]->getNumElements();
    double * element = rowArray_[1]->denseVector();
    for (i=0;i<number;i++) {
      int ii = index[i];
      if ( pivotVariable_[ii]==sequenceOut_) {
	pivotRow_=ii;
	alpha_=element[pivotRow_];
	break;
      }
    }
    assert (fabs(alpha_)>1.0e-12);
    // we are going to subtract movement from current basic
    double movement;
    // see where incoming will go to
    if (sequenceOut_<0||sequenceIn_==sequenceOut_) {
      // flip so go to bound
      movement  = ((directionIn_>0) ? upperIn_ : lowerIn_) - valueIn_;
    } else {
      // get where outgoing needs to get to
      double outValue = (directionOut_<0) ? upperOut_ : lowerOut_;
      // solutionOut_ - movement*alpha_ == outValue
      movement = (valueOut_-outValue)/alpha_;
      // set directionIn_ correctly
      directionIn_ = (movement>0) ? 1 :-1;
    }
    // update primal solution
    for (i=0;i<number;i++) {
      int ii = index[i];
      // get column
      int ij = pivotVariable_[ii];
      solution_[ij] -= movement*element[ii];
      element[ii]=0.0;
    }
    rowArray_[1]->setNumElements(0);
    // see where something went to
    if (sequenceOut_<0) {
      if (directionIn_<0) {
	assert (fabs(solution_[sequenceIn_]-upperIn_)<1.0e-7);
	solution_[sequenceIn_]=upperIn_;
      } else {
	assert (fabs(solution_[sequenceIn_]-lowerIn_)<1.0e-7);
	solution_[sequenceIn_]=lowerIn_;
      }
    } else {
      if (directionOut_<0) {
	assert (fabs(solution_[sequenceOut_]-upperOut_)<1.0e-7);
	solution_[sequenceOut_]=upperOut_;
      } else {
	assert (fabs(solution_[sequenceOut_]-lowerOut_)<1.0e-7);
	solution_[sequenceOut_]=lowerOut_;
      }
      solution_[sequenceIn_]=valueIn_+movement;
    }
    double objectiveChange = dualIn_*movement;
    // update duals
    if (pivotRow_>=0) {
      assert (fabs(alpha_)>1.0e-8);
      double multiplier = dualIn_/alpha_;
      rowArray_[0]->insert(pivotRow_,multiplier);
      factorization_->updateColumnTranspose(rowArray_[2],rowArray_[0]);
      // put row of tableau in rowArray[0] and columnArray[0]
      matrix_->transposeTimes(this,-1.0,
			      rowArray_[0],columnArray_[1],columnArray_[0]);
      // update column djs
      int i;
      int * index = columnArray_[0]->getIndices();
      int number = columnArray_[0]->getNumElements();
      double * element = columnArray_[0]->denseVector();
      for (i=0;i<number;i++) {
	int ii = index[i];
	dj_[ii] += element[ii];
	element[ii]=0.0;
      }
      columnArray_[0]->setNumElements(0);
      // and row djs
      index = rowArray_[0]->getIndices();
      number = rowArray_[0]->getNumElements();
      element = rowArray_[0]->denseVector();
      for (i=0;i<number;i++) {
	int ii = index[i];
	dj_[ii+numberColumns_] += element[ii];
	dual_[ii] = dj_[ii+numberColumns_];
	element[ii]=0.0;
      }
      rowArray_[0]->setNumElements(0);
      // check incoming
      assert (fabs(dj_[sequenceIn_])<1.0e-6);
    }
    
    // if stable replace in basis
    int updateStatus = factorization_->replaceColumn(this,
						     rowArray_[2],
						     rowArray_[1],
						     pivotRow_,
						     alpha_);
    bool takePivot=true;
    // if no pivots, bad update but reasonable alpha - take and invert
    if (updateStatus==2&&
	lastGoodIteration_==numberIterations_&&fabs(alpha_)>1.0e-5)
      updateStatus=4;
    if (updateStatus==1||updateStatus==4) {
      // slight error
      if (factorization_->pivots()>5||updateStatus==4) {
	returnCode=-1;
      }
    } else if (updateStatus==2) {
      // major error
      rowArray_[1]->clear();
      takePivot=false;
      if (factorization_->pivots()) {
	// refactorize here
	statusOfProblem();
	roundAgain=true;
      } else {
	returnCode=1;
      }
    } else if (updateStatus==3) {
      // out of memory
      // increase space if not many iterations
      if (factorization_->pivots()<
	  0.5*factorization_->maximumPivots()&&
	  factorization_->pivots()<200)
	factorization_->areaFactor(
				   factorization_->areaFactor() * 1.1);
      returnCode =-1; // factorize now
    }
    if (takePivot) {
      int save = algorithm_;
      // make simple so always primal
      algorithm_=1;
      housekeeping(objectiveChange);
      algorithm_=save;
    }
  }
  if (returnCode == -1) {
    // refactorize here
    statusOfProblem();
  } else {
    // just for now - recompute anyway
    gutsOfSolution(NULL,NULL);
  }
  return returnCode;
}

/* Pivot in a variable and choose an outgoing one.  Assumes primal
   feasible - will not go through a bound.  Returns step length in theta
   Returns ray in ray_ (or NULL if no pivot)
   Return codes as before but -1 means no acceptable pivot
*/
int ClpSimplex::primalPivotResult()
{
  assert (sequenceIn_>=0);
  valueIn_=solution_[sequenceIn_];
  lowerIn_=lower_[sequenceIn_];
  upperIn_=upper_[sequenceIn_];
  dualIn_=dj_[sequenceIn_];

  int returnCode = ((ClpSimplexPrimal *) this)->pivotResult();
  if (returnCode<0&&returnCode>-4) {
    return 0;
  } else {
    printf("Return code of %d from ClpSimplexPrimal::pivotResult\n",
	   returnCode);
    return -1;
  }
}
  
/* Pivot out a variable and choose an incoing one.  Assumes dual
   feasible - will not go through a reduced cost.  
   Returns step length in theta
   Returns ray in ray_ (or NULL if no pivot)
   Return codes as before but -1 means no acceptable pivot
*/
int 
ClpSimplex::dualPivotResult()
{
  return ((ClpSimplexDual *) this)->pivotResult();
}
// Factorization frequency
int 
ClpSimplex::factorizationFrequency() const
{
  if (factorization_)
    return factorization_->maximumPivots();
  else 
    return -1;
}
void 
ClpSimplex::setFactorizationFrequency(int value)
{
  if (factorization_)
    factorization_->maximumPivots(value);
}
// Common bits of coding for dual and primal
int 
ClpSimplex::startup(int ifValuesPass, int startFinishOptions)
{
  // Get rid of some arrays and empty factorization
  int useFactorization=false;
  if ((startFinishOptions&2)!=0&&(whatsChanged_&(2+512))==2+512)
    useFactorization=true; // Keep factorization if possible
  // sanity check
  // bad if empty (trap here to avoid using bad matrix_)
  if (!matrix_||(!matrix_->getNumElements()&&objective_->type()<2)) {
    int infeasNumber[2];
    double infeasSum[2];
    bool printIt = (specialOptions_&32768)==0 ? true : false; // no message if from Osi
    problemStatus_=emptyProblem(infeasNumber,infeasSum,printIt);
    numberDualInfeasibilities_=infeasNumber[0];
    sumDualInfeasibilities_=infeasSum[0];
    numberPrimalInfeasibilities_=infeasNumber[1];
    sumPrimalInfeasibilities_=infeasSum[1];
    return 2;
  }
  pivotRow_=-1;
  sequenceIn_=-1;
  sequenceOut_=-1;
  secondaryStatus_=0;

  primalTolerance_=dblParam_[ClpPrimalTolerance];
  dualTolerance_=dblParam_[ClpDualTolerance];
  if (problemStatus_!=10)
    numberIterations_=0;

  // put in standard form (and make row copy)
  // create modifiable copies of model rim and do optional scaling
  bool goodMatrix=createRim(7+8+16+32,true,startFinishOptions);

  if (goodMatrix) {
    // switch off factorization if bad
    if (pivotVariable_[0]<0)
      useFactorization=false;
    // Model looks okay
    // Do initial factorization
    // and set certain stuff
    // We can either set increasing rows so ...IsBasic gives pivot row
    // or we can just increment iBasic one by one
    // for now let ...iBasic give pivot row
    int saveThreshold = factorization_->denseThreshold();
    if (!useFactorization||factorization_->numberRows()!=numberRows_) {
      useFactorization = false;
      factorization_->increasingRows(2);
      // row activities have negative sign
      factorization_->slackValue(-1.0);
      factorization_->zeroTolerance(1.0e-13);
      // Switch off dense (unless special option set)
      if ((specialOptions_&8)==0)
        factorization_->setDenseThreshold(0);
    }
    // If values pass then perturb (otherwise may be optimal so leave a bit)
    if (ifValuesPass) {
      // do perturbation if asked for
      
      if (perturbation_<100) {
	if (algorithm_>0&&(objective_->type()<2||!objective_->activated())) {
	  ((ClpSimplexPrimal *) this)->perturb(0);
	} else if (algorithm_<0) {
	((ClpSimplexDual *) this)->perturb();
	}
      }
    }
    // for primal we will change bounds using infeasibilityCost_
    if (nonLinearCost_==NULL&&algorithm_>0) {
      // get a valid nonlinear cost function
      nonLinearCost_= new ClpNonLinearCost(this);
    }
    
    // loop round to clean up solution if values pass
    int numberThrownOut = -1;
    int totalNumberThrownOut=0;
    // see if we are re-using factorization
    if (!useFactorization) {
      while(numberThrownOut) {
	int status = internalFactorize(ifValuesPass ? 10 : 0);
	if (status<0)
	  return 1; // some error
	else
	  numberThrownOut = status;
	
	// for this we need clean basis so it is after factorize
	if (!numberThrownOut||numberThrownOut==numberRows_+1) {
	  // solution will be done again - skip if absolutely sure
	  if ((specialOptions_&512)==0||numberThrownOut==numberRows_+1) {
	    numberThrownOut = gutsOfSolution(  NULL,NULL,
					       ifValuesPass!=0);
	  } else {
	    // make sure not optimal at once
	    numberPrimalInfeasibilities_=1;
	    numberThrownOut=0;
	  }
	} else {
	  matrix_->rhsOffset(this,true); // redo rhs offset
	}
	totalNumberThrownOut+= numberThrownOut;
	
      }
    } else {
      // using previous factorization - we assume fine
      // but we need to say not optimal
      numberPrimalInfeasibilities_=1;
      numberDualInfeasibilities_=1;
      matrix_->rhsOffset(this,true); // redo rhs offset
    }
    if (totalNumberThrownOut)
      handler_->message(CLP_SINGULARITIES,messages_)
	<<totalNumberThrownOut
	<<CoinMessageEol;
    // Switch back dense
    factorization_->setDenseThreshold(saveThreshold);
    
    problemStatus_ = -1;
    
    // number of times we have declared optimality
    numberTimesOptimal_=0;

    return 0;
  } else {
    // bad matrix
    return 2;
  }
    
}


void 
ClpSimplex::finish(int startFinishOptions)
{
  // Get rid of some arrays and empty factorization
  int getRidOfData=1;
  if ((startFinishOptions&1)!=0||problemStatus_==10) {
    getRidOfData=0; // Keep stuff
    // mark all as current
    whatsChanged_ = 0xffff;
  } else {
    whatsChanged_=0;
  }
  deleteRim(getRidOfData);
  // Skip message if changing algorithms
  if (problemStatus_!=10) {
    if (problemStatus_==-1)
      problemStatus_=4;
    assert(problemStatus_>=0&&problemStatus_<6);
    if (handler_->detail(CLP_SIMPLEX_FINISHED,messages_)<100) {
      handler_->message(CLP_SIMPLEX_FINISHED+problemStatus_,messages_)
        <<objectiveValue()
        <<CoinMessageEol;
    }
  }
  factorization_->relaxAccuracyCheck(1.0);
  // get rid of any network stuff - could do more
  factorization_->cleanUp();
}
// Save data
ClpDataSave 
ClpSimplex::saveData() 
{
  ClpDataSave saved;
  saved.dualBound_ = dualBound_;
  saved.infeasibilityCost_ = infeasibilityCost_;
  saved.sparseThreshold_ = factorization_->sparseThreshold();
  saved.pivotTolerance_ = factorization_->pivotTolerance();
  saved.perturbation_ = perturbation_;
  saved.forceFactorization_=forceFactorization_;
  // Progress indicator
  delete progress_;
  progress_ = new ClpSimplexProgress (this);
  return saved;
}
// Restore data
void 
ClpSimplex::restoreData(ClpDataSave saved)
{
  factorization_->sparseThreshold(saved.sparseThreshold_);
  factorization_->pivotTolerance(saved.pivotTolerance_);
  perturbation_ = saved.perturbation_;
  infeasibilityCost_ = saved.infeasibilityCost_;
  dualBound_ = saved.dualBound_;
  forceFactorization_=saved.forceFactorization_;
  delete progress_;
  progress_=NULL;
}
// To flag a variable (not inline to allow for column generation)
void 
ClpSimplex::setFlagged( int sequence)
{
  status_[sequence] |= 64;
  matrix_->generalExpanded(this,7,sequence);
  lastFlaggedIteration_=numberIterations_;
}
/* Factorizes and returns true if optimal.  Used by user */
bool
ClpSimplex::statusOfProblem(bool initial)
{
  createRim(7+8+16+32);
  // is factorization okay?
  if (initial) {
    // First time - allow singularities
    int numberThrownOut = -1;
    int totalNumberThrownOut=0;
    while(numberThrownOut) {
      int status = internalFactorize(0);
      if (status==numberRows_+1)
	status=0; // all slack
      if (status<0) {
        deleteRim(-1);
	return false; // some error
      } else {
	numberThrownOut = status;
      }
      
      // for this we need clean basis so it is after factorize
      //if (!numberThrownOut)
      //numberThrownOut = gutsOfSolution(  NULL,NULL,
      //				   false);
      //else
      //matrix_->rhsOffset(this,true); // redo rhs offset
      totalNumberThrownOut+= numberThrownOut;
      
    }
    
    if (totalNumberThrownOut)
      handler_->message(CLP_SINGULARITIES,messages_)
	<<totalNumberThrownOut
	<<CoinMessageEol;
  } else {
#ifndef NDEBUG    
    int returnCode=internalFactorize(1);
    assert (!returnCode);
#else
    internalFactorize(1);
#endif
  }
  //memcpy(rowActivityWork_,rowActivity_,numberRows_*sizeof(double));
  //memcpy(columnActivityWork_,columnActivity_,numberColumns_*sizeof(double));
  gutsOfSolution(NULL,NULL);
  //memcpy(rowActivity_,rowActivityWork_,numberRows_*sizeof(double));
  //memcpy(columnActivity_,columnActivityWork_,numberColumns_*sizeof(double));
  //memcpy(reducedCost_,dj_,numberColumns_*sizeof(double));
  deleteRim(-1);
  return (primalFeasible()&&dualFeasible());
}
/* Return model - updates any scalars */
void 
ClpSimplex::returnModel(ClpSimplex & otherModel)
{
  ClpModel::returnModel(otherModel);
  otherModel.columnPrimalInfeasibility_ = columnPrimalInfeasibility_;
  otherModel.columnPrimalSequence_ = columnPrimalSequence_;
  otherModel.rowPrimalInfeasibility_ = rowPrimalInfeasibility_;
  otherModel.rowPrimalSequence_ = rowPrimalSequence_;
  otherModel.columnDualInfeasibility_ = columnDualInfeasibility_;
  otherModel.columnDualSequence_ = columnDualSequence_;
  otherModel.rowDualInfeasibility_ = rowDualInfeasibility_;
  otherModel.rowDualSequence_ = rowDualSequence_;
  otherModel.primalToleranceToGetOptimal_ = primalToleranceToGetOptimal_;
  otherModel.remainingDualInfeasibility_ = remainingDualInfeasibility_;
  otherModel.largestPrimalError_ = largestPrimalError_;
  otherModel.largestDualError_ = largestDualError_;
  otherModel.largestSolutionError_ = largestSolutionError_;
  otherModel.alpha_ = alpha_;
  otherModel.theta_ = theta_;
  otherModel.lowerIn_ = lowerIn_;
  otherModel.valueIn_ = valueIn_;
  otherModel.upperIn_ = upperIn_;
  otherModel.dualIn_ = dualIn_;
  otherModel.sequenceIn_ = sequenceIn_;
  otherModel.directionIn_ = directionIn_;
  otherModel.lowerOut_ = lowerOut_;
  otherModel.valueOut_ = valueOut_;
  otherModel.upperOut_ = upperOut_;
  otherModel.dualOut_ = dualOut_;
  otherModel.sequenceOut_ = sequenceOut_;
  otherModel.directionOut_ = directionOut_;
  otherModel.pivotRow_ = pivotRow_;
  otherModel.sumDualInfeasibilities_ = sumDualInfeasibilities_;
  otherModel.numberDualInfeasibilities_ = numberDualInfeasibilities_;
  otherModel.numberDualInfeasibilitiesWithoutFree_ = 
    numberDualInfeasibilitiesWithoutFree_;
  otherModel.sumPrimalInfeasibilities_ = sumPrimalInfeasibilities_;
  otherModel.numberPrimalInfeasibilities_ = numberPrimalInfeasibilities_;
  otherModel.numberTimesOptimal_ = numberTimesOptimal_;
  otherModel.sumOfRelaxedDualInfeasibilities_ = sumOfRelaxedDualInfeasibilities_;
  otherModel.sumOfRelaxedPrimalInfeasibilities_ = sumOfRelaxedPrimalInfeasibilities_;
}
/* Constructs a non linear cost from list of non-linearities (columns only)
   First lower of each column is taken as real lower
   Last lower is taken as real upper and cost ignored
   
   Returns nonzero if bad data e.g. lowers not monotonic
*/
int 
ClpSimplex::createPiecewiseLinearCosts(const int * starts,
				       const double * lower, const double * gradient)
{
  delete nonLinearCost_;
  // Set up feasible bounds and check monotonicity
  int iColumn;
  int returnCode=0;

  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    int iIndex = starts[iColumn];
    int end = starts[iColumn+1]-1;
    columnLower_[iColumn] = lower[iIndex];
    columnUpper_[iColumn] = lower[end];
    double value = columnLower_[iColumn];
    iIndex++;
    for (;iIndex<end;iIndex++) {
      if (lower[iIndex]<value)
	returnCode++; // not monotonic
      value=lower[iIndex];
    }
  }
  nonLinearCost_ = new ClpNonLinearCost(this,starts,lower,gradient);
  specialOptions_ |= 2; // say keep
  return returnCode;
}
/* For advanced use.  When doing iterative solves things can get
   nasty so on values pass if incoming solution has largest
   infeasibility < incomingInfeasibility throw out variables
   from basis until largest infeasibility < allowedInfeasibility
   or incoming largest infeasibility.
   If allowedInfeasibility>= incomingInfeasibility this is
   always possible altough you may end up with an all slack basis.
   
   Defaults are 1.0,10.0
*/
void 
ClpSimplex::setValuesPassAction(float incomingInfeasibility,
				float allowedInfeasibility)
{
  incomingInfeasibility_=incomingInfeasibility;
  allowedInfeasibility_=allowedInfeasibility;
  CoinAssert(incomingInfeasibility_>=0.0);
  CoinAssert(allowedInfeasibility_>=incomingInfeasibility_);
}
//#############################################################################

ClpSimplexProgress::ClpSimplexProgress () 
{
  int i;
  for (i=0;i<CLP_PROGRESS;i++) {
    objective_[i] = COIN_DBL_MAX;
    infeasibility_[i] = -1.0; // set to an impossible value
    realInfeasibility_[i] = COIN_DBL_MAX;
    numberInfeasibilities_[i]=-1; 
    iterationNumber_[i]=-1;
  }
  for (i=0;i<CLP_CYCLE;i++) {
    //obj_[i]=COIN_DBL_MAX;
    in_[i]=-1;
    out_[i]=-1;
    way_[i]=0;
  }
  numberTimes_ = 0;
  numberBadTimes_ = 0;
  model_ = NULL;
  oddState_=0;
}


//-----------------------------------------------------------------------------

ClpSimplexProgress::~ClpSimplexProgress ()
{
}
// Copy constructor. 
ClpSimplexProgress::ClpSimplexProgress(const ClpSimplexProgress &rhs) 
{
  int i;
  for (i=0;i<CLP_PROGRESS;i++) {
    objective_[i] = rhs.objective_[i];
    infeasibility_[i] = rhs.infeasibility_[i];
    realInfeasibility_[i] = rhs.realInfeasibility_[i];
    numberInfeasibilities_[i]=rhs.numberInfeasibilities_[i]; 
    iterationNumber_[i]=rhs.iterationNumber_[i];
  }
  for (i=0;i<CLP_CYCLE;i++) {
    //obj_[i]=rhs.obj_[i];
    in_[i]=rhs.in_[i];
    out_[i]=rhs.out_[i];
    way_[i]=rhs.way_[i];
  }
  numberTimes_ = rhs.numberTimes_;
  numberBadTimes_ = rhs.numberBadTimes_;
  model_ = rhs.model_;
  oddState_=rhs.oddState_;
}
// Copy constructor.from model
ClpSimplexProgress::ClpSimplexProgress(ClpSimplex * model) 
{
  model_ = model;
  int i;
  for (i=0;i<CLP_PROGRESS;i++) {
    if (model_->algorithm()>=0)
      objective_[i] = COIN_DBL_MAX;
    else
      objective_[i] = -COIN_DBL_MAX;
    infeasibility_[i] = -1.0; // set to an impossible value
    realInfeasibility_[i] = COIN_DBL_MAX;
    numberInfeasibilities_[i]=-1; 
    iterationNumber_[i]=-1;
  }
  for (i=0;i<CLP_CYCLE;i++) {
    //obj_[i]=COIN_DBL_MAX;
    in_[i]=-1;
    out_[i]=-1;
    way_[i]=0;
  }
  numberTimes_ = 0;
  numberBadTimes_ = 0;
  oddState_=0;
}
// Assignment operator. This copies the data
ClpSimplexProgress & 
ClpSimplexProgress::operator=(const ClpSimplexProgress & rhs)
{
  if (this != &rhs) {
    int i;
    for (i=0;i<CLP_PROGRESS;i++) {
      objective_[i] = rhs.objective_[i];
      infeasibility_[i] = rhs.infeasibility_[i];
      realInfeasibility_[i] = rhs.realInfeasibility_[i];
      numberInfeasibilities_[i]=rhs.numberInfeasibilities_[i]; 
      iterationNumber_[i]=rhs.iterationNumber_[i];
    }
    for (i=0;i<CLP_CYCLE;i++) {
      //obj_[i]=rhs.obj_[i];
      in_[i]=rhs.in_[i];
      out_[i]=rhs.out_[i];
      way_[i]=rhs.way_[i];
    }
    numberTimes_ = rhs.numberTimes_;
    numberBadTimes_ = rhs.numberBadTimes_;
    model_ = rhs.model_;
    oddState_=rhs.oddState_;
  }
  return *this;
}
// Seems to be something odd about exact comparison of doubles on linux
static bool equalDouble(double value1, double value2)
{
  unsigned int *i1 = (unsigned int *) &value1;
  unsigned int *i2 = (unsigned int *) &value2;
  if (sizeof(unsigned int)*2==sizeof(double)) 
    return (i1[0]==i2[0]&&i1[1]==i2[1]);
  else
    return (i1[0]==i2[0]);
}
int
ClpSimplexProgress::looping()
{
  if (!model_)
    return -1;
  double objective = model_->rawObjectiveValue();
  double infeasibility;
  double realInfeasibility=0.0;
  int numberInfeasibilities;
  int iterationNumber = model_->numberIterations();
  if (model_->algorithm()<0) {
    // dual
    infeasibility = model_->sumPrimalInfeasibilities();
    numberInfeasibilities = model_->numberPrimalInfeasibilities();
  } else {
    //primal
    infeasibility = model_->sumDualInfeasibilities();
    realInfeasibility = model_->nonLinearCost()->sumInfeasibilities();
    numberInfeasibilities = model_->numberDualInfeasibilities();
  }
  int i;
  int numberMatched=0;
  int matched=0;
  int nsame=0;
  for (i=0;i<CLP_PROGRESS;i++) {
    bool matchedOnObjective = equalDouble(objective,objective_[i]);
    bool matchedOnInfeasibility = equalDouble(infeasibility,infeasibility_[i]);
    bool matchedOnInfeasibilities = 
      (numberInfeasibilities==numberInfeasibilities_[i]);
    
    if (matchedOnObjective&&matchedOnInfeasibility&&matchedOnInfeasibilities) {
      matched |= (1<<i);
      // Check not same iteration
      if (iterationNumber!=iterationNumber_[i]) {
	numberMatched++;
	// here mainly to get over compiler bug?
	if (model_->messageHandler()->logLevel()>10)
	  printf("%d %d %d %d %d loop check\n",i,numberMatched,
		 matchedOnObjective, matchedOnInfeasibility, 
		 matchedOnInfeasibilities);
      } else {
	// stuck but code should notice
	nsame++;
      }
    }
    if (i) {
      objective_[i-1] = objective_[i];
      infeasibility_[i-1] = infeasibility_[i];
      realInfeasibility_[i-1] = realInfeasibility_[i];
      numberInfeasibilities_[i-1]=numberInfeasibilities_[i]; 
      iterationNumber_[i-1]=iterationNumber_[i];
    }
  }
  objective_[CLP_PROGRESS-1] = objective;
  infeasibility_[CLP_PROGRESS-1] = infeasibility;
  realInfeasibility_[CLP_PROGRESS-1] = realInfeasibility;
  numberInfeasibilities_[CLP_PROGRESS-1]=numberInfeasibilities;
  iterationNumber_[CLP_PROGRESS-1]=iterationNumber;
  if (nsame==CLP_PROGRESS)
    numberMatched=CLP_PROGRESS; // really stuck
  if (model_->progressFlag())
    numberMatched=0;
  numberTimes_++;
  if (numberTimes_<10)
    numberMatched=0;
  // skip if just last time as may be checking something
  if (matched==(1<<(CLP_PROGRESS-1)))
    numberMatched=0;
  if (numberMatched) {
    model_->messageHandler()->message(CLP_POSSIBLELOOP,model_->messages())
      <<numberMatched
      <<matched
      <<numberTimes_
      <<CoinMessageEol;
    numberBadTimes_++;
    if (numberBadTimes_<10) {
      // make factorize every iteration
      model_->forceFactorization(1);
      if (numberBadTimes_<2) {
	startCheck(); // clear other loop check
	if (model_->algorithm()<0) {
	  // dual - change tolerance
	  model_->setCurrentDualTolerance(model_->currentDualTolerance()*1.05);
	  // if infeasible increase dual bound
	  if (model_->dualBound()<1.0e17) {
	    model_->setDualBound(model_->dualBound()*1.1);
	  }
	} else {
	  // primal - change tolerance	
	  if (numberBadTimes_>3)
	    model_->setCurrentPrimalTolerance(model_->currentPrimalTolerance()*1.05);
	  // if infeasible increase infeasibility cost
	  if (model_->nonLinearCost()->numberInfeasibilities()&&
	      model_->infeasibilityCost()<1.0e17) {
	    model_->setInfeasibilityCost(model_->infeasibilityCost()*1.1);
	  }
	}
      } else {
	// flag
	int iSequence;
	if (model_->algorithm()<0) {
	  // dual
	  if (model_->dualBound()>1.0e14) 
	    model_->setDualBound(1.0e14);
	  iSequence=in_[CLP_CYCLE-1];
	} else {
	  // primal 
	  if (model_->infeasibilityCost()>1.0e14) 
	    model_->setInfeasibilityCost(1.0e14);
	  iSequence=out_[CLP_CYCLE-1];
	}
	if (iSequence>=0) {
	  char x = model_->isColumn(iSequence) ? 'C' :'R';
	  if (model_->messageHandler()->logLevel()>=63)
	    model_->messageHandler()->message(CLP_SIMPLEX_FLAG,model_->messages())
	      <<x<<model_->sequenceWithin(iSequence)
	      <<CoinMessageEol;
	  model_->setFlagged(iSequence);
	  //printf("flagging %d from loop\n",iSequence);
	  startCheck();
	} else {
	  // Give up
	  if (model_->messageHandler()->logLevel()>=63)
	    printf("***** All flagged?\n");
	  return 4;
	}
	// reset
	numberBadTimes_=2;
      }
      return -2;
    } else {
      // look at solution and maybe declare victory
      if (infeasibility<1.0e-4) {
	return 0;
      } else {
	model_->messageHandler()->message(CLP_LOOP,model_->messages())
	  <<CoinMessageEol;
#ifndef NDEBUG
        printf("debug loop ClpSimplex A\n");
	abort();
#endif
	return 3;
      }
    }
  }
  return -1;
}
// Returns previous objective (if -1) - current if (0)
double 
ClpSimplexProgress::lastObjective(int back) const
{
  return objective_[CLP_PROGRESS-1-back];
}
// Returns previous infeasibility (if -1) - current if (0)
double 
ClpSimplexProgress::lastInfeasibility(int back) const
{
  return realInfeasibility_[CLP_PROGRESS-1-back];
}
// Sets real primal infeasibility
void
ClpSimplexProgress::setInfeasibility(double value)
{
  for (int i=1;i<CLP_PROGRESS;i++) 
    realInfeasibility_[i-1] = realInfeasibility_[i];
  realInfeasibility_[CLP_PROGRESS-1]=value;
}
// Modify objective e.g. if dual infeasible in dual
void 
ClpSimplexProgress::modifyObjective(double value)
{
  objective_[CLP_PROGRESS-1]=value;
}
// Returns previous iteration number (if -1) - current if (0)
int 
ClpSimplexProgress::lastIterationNumber(int back) const
{
  return iterationNumber_[CLP_PROGRESS-1-back];
}
// clears iteration numbers (to switch off panic)
void 
ClpSimplexProgress::clearIterationNumbers()
{
  for (int i=0;i<CLP_PROGRESS;i++) 
    iterationNumber_[i]=-1;
}
// Start check at beginning of whileIterating
void 
ClpSimplexProgress::startCheck()
{
  int i;
  for (i=0;i<CLP_CYCLE;i++) {
    //obj_[i]=COIN_DBL_MAX;
    in_[i]=-1;
    out_[i]=-1;
    way_[i]=0;
  }
}
// Returns cycle length in whileIterating
int 
ClpSimplexProgress::cycle(int in, int out,int wayIn,int wayOut)
{
  int i;
#if 0
  if (model_->numberIterations()>206571) {
    printf("in %d out %d\n",in,out);
    for (i=0;i<CLP_CYCLE;i++) 
      printf("cy %d in %d out %d\n",i,in_[i],out_[i]);
  }
#endif
  int matched=0;
  // first see if in matches any out
  for (i=1;i<CLP_CYCLE;i++) {
    if (in==out_[i]) {
      // even if flip then suspicious
      matched=-1;
      break;
    }
  }
#if 0
  if (!matched||in_[0]<0) {
    // can't be cycle
    for (i=0;i<CLP_CYCLE-1;i++) {
      //obj_[i]=obj_[i+1];
      in_[i]=in_[i+1];
      out_[i]=out_[i+1];
      way_[i]=way_[i+1];
    }
  } else {
    // possible cycle
    matched=0;
    for (i=0;i<CLP_CYCLE-1;i++) {
      int k;
      char wayThis = way_[i];
      int inThis = in_[i];
      int outThis = out_[i];
      //double objThis = obj_[i];
      for(k=i+1;k<CLP_CYCLE;k++) {
	if (inThis==in_[k]&&outThis==out_[k]&&wayThis==way_[k]) {
	  int distance = k-i;
	  if (k+distance<CLP_CYCLE) {
	    // See if repeats
	    int j=k+distance;
	    if (inThis==in_[j]&&outThis==out_[j]&&wayThis==way_[j]) {
	      matched=distance;
	      break;
	    }
	  } else {
	    matched=distance;
	    break;
	  }
	}
      }
      //obj_[i]=obj_[i+1];
      in_[i]=in_[i+1];
      out_[i]=out_[i+1];
      way_[i]=way_[i+1];
    }
  }
#else
  if (matched&&in_[0]>=0) {
    // possible cycle - only check [0] against all
    matched=0;
    int nMatched=0;
    char way0 = way_[0];
    int in0 = in_[0];
    int out0 = out_[0];
    //double obj0 = obj_[i];
    for(int k=1;k<CLP_CYCLE-4;k++) {
      if (in0==in_[k]&&out0==out_[k]&&way0==way_[k]) {
	nMatched++;
	// See if repeats
	int end = CLP_CYCLE-k;
	int j;
	for ( j=1;j<end;j++) {
	  if (in_[j+k]!=in_[j]||out_[j+k]!=out_[j]||way_[j+k]!=way_[j]) 
	    break;
	}
	if (j==end) {
	  matched=k;
	  break;
	}
      }
    }
    // If three times then that is too much even if not regular
    if (matched<=0&&nMatched>1)
      matched=100;
  }
  for (i=0;i<CLP_CYCLE-1;i++) {
    //obj_[i]=obj_[i+1];
    in_[i]=in_[i+1];
    out_[i]=out_[i+1];
    way_[i]=way_[i+1];
  }
#endif
  char way = 1-wayIn+4*(1-wayOut);
  //obj_[i]=model_->objectiveValue();
  in_[CLP_CYCLE-1]=in;
  out_[CLP_CYCLE-1]=out;
  way_[CLP_CYCLE-1]=way;
  return matched;
}
// Allow initial dense factorization
void 
ClpSimplex::setInitialDenseFactorization(bool onOff)
{
  if (onOff)
    specialOptions_ |= 8;
  else
    specialOptions_ &= ~8;
}
bool 
ClpSimplex::initialDenseFactorization() const
{
  return (specialOptions_&8)!=0;
}
/* This constructor modifies original ClpSimplex and stores
   original stuff in created ClpSimplex.  It is only to be used in
   conjunction with originalModel */
ClpSimplex::ClpSimplex (ClpSimplex * wholeModel,
			int numberColumns, const int * whichColumns)
{

  // Set up dummy row selection list
  numberRows_ = wholeModel->numberRows_;
  int * whichRow = new int [numberRows_];
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++)
    whichRow[iRow]=iRow;
  // ClpModel stuff (apart from numberColumns_)
  matrix_ = wholeModel->matrix_;
  rowCopy_ = wholeModel->rowCopy_;
  if (wholeModel->rowCopy_) {
    // note reversal of order
    wholeModel->rowCopy_ = wholeModel->rowCopy_->subsetClone(numberRows_,whichRow,
							     numberColumns,whichColumns);
  } else {
    wholeModel->rowCopy_=NULL;
  }
  whatsChanged_=0;
  CoinAssert (wholeModel->matrix_);
  wholeModel->matrix_ = wholeModel->matrix_->subsetClone(numberRows_,whichRow,
					numberColumns,whichColumns);
  delete [] whichRow;
  numberColumns_ = wholeModel->numberColumns_;
  // Now ClpSimplex stuff and status_
  delete  wholeModel->primalColumnPivot_;
  wholeModel->primalColumnPivot_ = new ClpPrimalColumnSteepest(0);
  nonLinearCost_ = wholeModel->nonLinearCost_;

  // Now main arrays
  int iColumn;
  int numberTotal = numberRows_+numberColumns;
  printf("%d %d %d\n",numberTotal,numberRows_,numberColumns);
  // mapping 
  int * mapping = new int[numberRows_+numberColumns_];
  for (iColumn=0;iColumn<numberColumns_;iColumn++) 
    mapping[iColumn]=-1;
  for (iRow=0;iRow<numberRows_;iRow++) 
    mapping[iRow+numberColumns_] = iRow+numberColumns;
  // Redo costs and bounds of whole model
  wholeModel->createRim(1+4,false);
  lower_ = wholeModel->lower_;
  wholeModel->lower_ = new double [numberTotal];
  memcpy(wholeModel->lower_+numberColumns,lower_+numberColumns_,numberRows_*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int jColumn = whichColumns[iColumn];
    wholeModel->lower_[iColumn]=lower_[jColumn];
    // and pointer back 
    mapping[jColumn]=iColumn;
  }
#ifdef CLP_DEBUG
  for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) 
    printf("mapx %d %d\n",iColumn,mapping[iColumn]);
#endif
  // Re-define pivotVariable_
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot = wholeModel->pivotVariable_[iRow];
    wholeModel->pivotVariable_[iRow]=mapping[iPivot];
#ifdef CLP_DEBUG
    printf("p row %d, pivot %d -> %d\n",iRow,iPivot,mapping[iPivot]);
#endif
    assert (wholeModel->pivotVariable_[iRow]>=0);
  }
  // Reverse mapping (so extended version of whichColumns)
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    mapping[iColumn]=whichColumns[iColumn];
  for (;iColumn<numberRows_+numberColumns;iColumn++) 
    mapping[iColumn] = iColumn + (numberColumns_-numberColumns);
#ifdef CLP_DEBUG
  for (iColumn=0;iColumn<numberRows_+numberColumns;iColumn++) 
    printf("map %d %d\n",iColumn,mapping[iColumn]);
#endif
  // Save mapping somewhere - doesn't matter
  rowUpper_ = (double *) mapping;
  upper_ = wholeModel->upper_;
  wholeModel->upper_ = new double [numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    wholeModel->upper_[iColumn]=upper_[jColumn];
  }
  cost_ = wholeModel->cost_;
  wholeModel->cost_ = new double [numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    wholeModel->cost_[iColumn]=cost_[jColumn];
  }
  dj_ = wholeModel->dj_;
  wholeModel->dj_ = new double [numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    wholeModel->dj_[iColumn]=dj_[jColumn];
  }
  solution_ = wholeModel->solution_;
  wholeModel->solution_ = new double [numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    wholeModel->solution_[iColumn]=solution_[jColumn];
  }
  // now see what variables left out do to row solution
  double * rowSolution = wholeModel->solution_+numberColumns;
  double * fullSolution = solution_;
  double * sumFixed = new double[numberRows_];
  memset (sumFixed,0,numberRows_*sizeof(double));
  // zero out ones in small problem
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int jColumn = mapping[iColumn];
    fullSolution[jColumn]=0.0;
  }
  // Get objective offset
  double originalOffset;
  wholeModel->getDblParam(ClpObjOffset,originalOffset);
  double offset=0.0;
  const double * cost = cost_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) 
    offset += fullSolution[iColumn]*cost[iColumn];
  wholeModel->setDblParam(ClpObjOffset,originalOffset-offset);
  setDblParam(ClpObjOffset,originalOffset);
  matrix_->times(1.0,fullSolution,sumFixed,wholeModel->rowScale_,wholeModel->columnScale_);
      
  double * lower = lower_+numberColumns;
  double * upper = upper_+numberColumns;
  double fixed=0.0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    fixed += fabs(sumFixed[iRow]);
    if (lower[iRow]>-1.0e50) 
      lower[iRow] -= sumFixed[iRow];
    if (upper[iRow]<1.0e50)
      upper[iRow] -= sumFixed[iRow];
    rowSolution[iRow] -= sumFixed[iRow];
  }
  printf("offset %g sumfixed %g\n",offset,fixed);
  delete [] sumFixed;
  columnScale_ = wholeModel->columnScale_;
  if (columnScale_) {
    wholeModel->columnScale_ = new double [numberTotal];
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int jColumn = mapping[iColumn];
      wholeModel->columnScale_[iColumn]=columnScale_[jColumn];
    }
  }
  status_ = wholeModel->status_;
  wholeModel->status_ = new unsigned char [numberTotal];
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    wholeModel->status_[iColumn]=status_[jColumn];
  }
  savedSolution_ = wholeModel->savedSolution_;
  if (savedSolution_) {
    wholeModel->savedSolution_ = new double [numberTotal];
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      int jColumn = mapping[iColumn];
      wholeModel->savedSolution_[iColumn]=savedSolution_[jColumn];
    }
  }
  saveStatus_ = wholeModel->saveStatus_;
  if (saveStatus_) {
    wholeModel->saveStatus_ = new unsigned char [numberTotal];
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      int jColumn = mapping[iColumn];
      wholeModel->saveStatus_[iColumn]=saveStatus_[jColumn];
    }
  }
  
  wholeModel->numberColumns_ = numberColumns;
  // Initialize weights
  wholeModel->primalColumnPivot_->saveWeights(wholeModel,2);
  // Costs
  wholeModel->nonLinearCost_ = new ClpNonLinearCost(wholeModel);
  wholeModel->nonLinearCost_->checkInfeasibilities();
  printf("after contraction %d infeasibilities summing to %g\n",
	 nonLinearCost_->numberInfeasibilities(),nonLinearCost_->sumInfeasibilities());
  // Redo some stuff
  wholeModel->reducedCostWork_ = wholeModel->dj_;
  wholeModel->rowReducedCost_ = wholeModel->dj_+wholeModel->numberColumns_;
  wholeModel->columnActivityWork_ = wholeModel->solution_;
  wholeModel->rowActivityWork_ = wholeModel->solution_+wholeModel->numberColumns_;
  wholeModel->objectiveWork_ = wholeModel->cost_;
  wholeModel->rowObjectiveWork_ = wholeModel->cost_+wholeModel->numberColumns_;
  wholeModel->rowLowerWork_ = wholeModel->lower_+wholeModel->numberColumns_;
  wholeModel->columnLowerWork_ = wholeModel->lower_;
  wholeModel->rowUpperWork_ = wholeModel->upper_+wholeModel->numberColumns_;
  wholeModel->columnUpperWork_ = wholeModel->upper_;
#ifndef NDEBUG
  // Check status
  ClpSimplex * xxxx = wholeModel;
  int nBasic=0;
  for (iColumn=0;iColumn<xxxx->numberRows_+xxxx->numberColumns_;iColumn++)
    if (xxxx->getStatus(iColumn)==basic)
      nBasic++;
  assert (nBasic==xxxx->numberRows_);
  for (iRow=0;iRow<xxxx->numberRows_;iRow++) {
    int iPivot=xxxx->pivotVariable_[iRow];
    assert (xxxx->getStatus(iPivot)==basic);
  }
#endif
}
/* This copies back stuff from miniModel and then deletes miniModel.
   Only to be used with mini constructor */
void 
ClpSimplex::originalModel(ClpSimplex * miniModel)
{
  int numberSmall = numberColumns_;
  numberColumns_ = miniModel->numberColumns_;
  int numberTotal = numberSmall+numberRows_;
  // copy back
  int iColumn;
  int * mapping = (int *) miniModel->rowUpper_;
#ifdef CLP_DEBUG
  for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) 
    printf("mapb %d %d\n",iColumn,mapping[iColumn]);
#endif
  // miniModel actually has full arrays
  // now see what variables left out do to row solution
  double * fullSolution = miniModel->solution_;
  double * sumFixed = new double[numberRows_];
  memset (sumFixed,0,numberRows_*sizeof(double));
  miniModel->matrix_->times(1.0,fullSolution,sumFixed,rowScale_,miniModel->columnScale_);
      
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    int jColumn = mapping[iColumn];
    miniModel->lower_[jColumn]=lower_[iColumn];
    miniModel->upper_[jColumn]=upper_[iColumn];
    miniModel->cost_[jColumn]=cost_[iColumn];
    miniModel->dj_[jColumn]=dj_[iColumn];
    miniModel->solution_[jColumn]=solution_[iColumn];
    miniModel->status_[jColumn]=status_[iColumn];
#ifdef CLP_DEBUG
    printf("%d in small -> %d in original\n",iColumn,jColumn);
#endif
  }
  delete [] lower_;
  lower_ =  miniModel->lower_;
  delete [] upper_;
  upper_ =  miniModel->upper_;
  delete [] cost_;
  cost_ =  miniModel->cost_;
  delete [] dj_;
  dj_ =  miniModel->dj_;
  delete [] solution_;
  solution_ =  miniModel->solution_;
  delete [] status_;
  status_ =  miniModel->status_;
  if (columnScale_) {
    for (iColumn=0;iColumn<numberSmall;iColumn++) {
      int jColumn = mapping[iColumn];
      miniModel->columnScale_[jColumn]=columnScale_[iColumn];
    }
    delete [] columnScale_;
    columnScale_ =  miniModel->columnScale_;
  }
  if (savedSolution_) {
    if (!miniModel->savedSolution_) {
      miniModel->savedSolution_ = ClpCopyOfArray(solution_,numberColumns_+numberRows_);
    } else {
      for (iColumn=0;iColumn<numberTotal;iColumn++) {
	int jColumn = mapping[iColumn];
	miniModel->savedSolution_[jColumn]=savedSolution_[iColumn];
      }
    }
    delete [] savedSolution_;
    savedSolution_ =  miniModel->savedSolution_;
  }
  if (saveStatus_) {
    if (!miniModel->saveStatus_) {
      miniModel->saveStatus_ = ClpCopyOfArray(status_,numberColumns_+numberRows_);
    } else {
      for (iColumn=0;iColumn<numberTotal;iColumn++) {
	int jColumn = mapping[iColumn];
	miniModel->saveStatus_[jColumn]=saveStatus_[iColumn];
      }
    }
    delete [] saveStatus_;
    saveStatus_ =  miniModel->saveStatus_;
  }
  // Re-define pivotVariable_
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot = pivotVariable_[iRow];
#ifdef CLP_DEBUG
    printf("pb row %d, pivot %d -> %d\n",iRow,iPivot,mapping[iPivot]);
#endif
    pivotVariable_[iRow]=mapping[iPivot];
    assert (pivotVariable_[iRow]>=0);
  }
  // delete stuff and move back
  delete matrix_;
  delete rowCopy_;
  delete primalColumnPivot_;
  delete nonLinearCost_;
  matrix_ = miniModel->matrix_;
  rowCopy_ = miniModel->rowCopy_;
  nonLinearCost_ = miniModel->nonLinearCost_;
  double originalOffset;
  miniModel->getDblParam(ClpObjOffset,originalOffset);
  setDblParam(ClpObjOffset,originalOffset);
  // Redo some stuff
  reducedCostWork_ = dj_;
  rowReducedCost_ = dj_+numberColumns_;
  columnActivityWork_ = solution_;
  rowActivityWork_ = solution_+numberColumns_;
  objectiveWork_ = cost_;
  rowObjectiveWork_ = cost_+numberColumns_;
  rowLowerWork_ = lower_+numberColumns_;
  columnLowerWork_ = lower_;
  rowUpperWork_ = upper_+numberColumns_;
  columnUpperWork_ = upper_;
  // Cleanup
  for (iRow=0;iRow<numberRows_;iRow++) {
    double value = rowActivityWork_[iRow] + sumFixed[iRow];
    rowActivityWork_[iRow] = value;
    switch(getRowStatus(iRow)) {
      
    case basic:
      break;
    case atUpperBound:
      //rowActivityWork_[iRow]=rowUpperWork_[iRow];
      break;
    case ClpSimplex::isFixed:
    case atLowerBound:
      //rowActivityWork_[iRow]=rowLowerWork_[iRow];
      break;
    case isFree:
      break;
      // superbasic
    case superBasic:
      break;
    }
  }
  delete [] sumFixed;
  nonLinearCost_->checkInfeasibilities();
  printf("in original %d infeasibilities summing to %g\n",
	 nonLinearCost_->numberInfeasibilities(),nonLinearCost_->sumInfeasibilities());
  // Initialize weights
  primalColumnPivot_ = new ClpPrimalColumnSteepest(10);
  primalColumnPivot_->saveWeights(this,2);
#ifndef NDEBUG
  // Check status
  ClpSimplex * xxxx = this;
  int nBasic=0;
  for (iColumn=0;iColumn<xxxx->numberRows_+xxxx->numberColumns_;iColumn++)
    if (xxxx->getStatus(iColumn)==basic)
      nBasic++;
  assert (nBasic==xxxx->numberRows_);
  for (iRow=0;iRow<xxxx->numberRows_;iRow++) {
    int iPivot=xxxx->pivotVariable_[iRow];
    assert (xxxx->getStatus(iPivot)==basic);
  }
#endif
}
// Pass in Event handler (cloned and deleted at end)
void 
ClpSimplex::passInEventHandler(const ClpEventHandler * eventHandler)
{
  delete eventHandler_;
  eventHandler_ = eventHandler->clone();
  eventHandler_->setSimplex(this);
}
#ifndef NDEBUG
// For errors to make sure print to screen
// only called in debug mode
static void indexError(int index,
			std::string methodName)
{
  std::cerr<<"Illegal index "<<index<<" in ClpSimplex::"<<methodName<<std::endl;
  throw CoinError("Illegal index",methodName,"ClpSimplex");
}
#endif
// These are only to be used using startFinishOptions (ClpSimplexDual, ClpSimplexPrimal)
//Get a row of the tableau (slack part in slack if not NULL)
void 
ClpSimplex::getBInvARow(int row, double* z, double * slack)
{
#ifndef NDEBUG
  int n = numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvARow");
  }
#endif
  if (!rowArray_[0]) {
    printf("ClpSimplexPrimal or ClpSimplexDual must have been called with correct startFinishOption\n");
    abort();
  }
  CoinIndexedVector * rowArray0 = rowArray(0);
  CoinIndexedVector * rowArray1 = rowArray(1);
  CoinIndexedVector * columnArray0 = columnArray(0);
  CoinIndexedVector * columnArray1 = columnArray(1);
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
  // put +1 in row 
  // But swap if pivot variable was slack as clp stores slack as -1.0
  int pivot = pivotVariable_[row];
  double value;
  // And if scaled then adjust
  if (!rowScale_) {
    if (pivot<numberColumns_)
      value = 1.0;
    else
      value = -1.0;
  } else {
    if (pivot<numberColumns_)
      value = columnScale_[pivot];
    else
      value = -1.0/rowScale_[pivot-numberColumns_];
  }
  rowArray1->insert(row,value);
  factorization_->updateColumnTranspose(rowArray0,rowArray1);
  // put row of tableau in rowArray1 and columnArray0
  clpMatrix()->transposeTimes(this,1.0,
			    rowArray1,columnArray1,columnArray0);
  if (!rowScale_) {
    memcpy(z,columnArray0->denseVector(),
	   numberColumns_*sizeof(double));
  } else {
    double * array = columnArray0->denseVector();
    for (int i=0;i<numberColumns_;i++)
      z[i] = array[i]/columnScale_[i];
  }
  if (slack) {
    if (!rowScale_) {
      memcpy(slack,rowArray1->denseVector(),
	     numberRows_*sizeof(double));
    } else {
      double * array = rowArray1->denseVector();
      for (int i=0;i<numberRows_;i++)
	slack[i] = array[i]*rowScale_[i];
    }
  }
  // don't need to clear everything always, but doesn't cost
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
}

//Get a row of the basis inverse
void 
ClpSimplex::getBInvRow(int row, double* z)

{
#ifndef NDEBUG
  int n = numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvRow");
  }
#endif
  if (!rowArray_[0]) {
    printf("ClpSimplexPrimal or ClpSimplexDual must have been called with correct startFinishOption\n");
    abort();
  }
  ClpFactorization * factorization = factorization_;
  CoinIndexedVector * rowArray0 = rowArray(0);
  CoinIndexedVector * rowArray1 = rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // put +1 in row
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double value = (pivotVariable_[row]<numberColumns_) ? 1.0 : -1.0;
  // What about scaling ?
  rowArray1->insert(row,value);
  factorization->updateColumnTranspose(rowArray0,rowArray1);
  memcpy(z,rowArray1->denseVector(),numberRows()*sizeof(double));
  rowArray1->clear();
}

//Get a column of the tableau
void 
ClpSimplex::getBInvACol(int col, double* vec)
{
  if (!rowArray_[0]) {
    printf("ClpSimplexPrimal or ClpSimplexDual should have been called with correct startFinishOption\n");
    abort();
  }
  CoinIndexedVector * rowArray0 = rowArray(0);
  CoinIndexedVector * rowArray1 = rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // get column of matrix
#ifndef NDEBUG
  int n = numberColumns_+numberRows_;
  if (col<0||col>=n) {
    indexError(col,"getBInvACol");
  }
#endif
  if (!rowScale_) {
    if (col<numberColumns_) {
      unpack(rowArray1,col);
    } else {
      rowArray1->insert(col-numberColumns_,1.0);
    }
  } else {
    if (col<numberColumns_) {
      unpack(rowArray1,col);
      double multiplier = 1.0/columnScale_[col];
      int number = rowArray1->getNumElements();
      int * index = rowArray1->getIndices();
      double * array = rowArray1->denseVector();
      for (int i=0;i<number;i++) {
	int iRow = index[i];
	// make sure not packed
	assert (array[iRow]);
	array[iRow] *= multiplier;
      }
    } else {
      rowArray1->insert(col-numberColumns_,rowScale_[col-numberColumns_]);
    }
  }
  factorization_->updateColumn(rowArray0,rowArray1,false);
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double * array = rowArray1->denseVector();
  if (!rowScale_) {
    for (int i=0;i<numberRows_;i++) {
      double multiplier = (pivotVariable_[i]<numberColumns_) ? 1.0 : -1.0;
      vec[i] = multiplier * array[i];
    }
  } else {
    for (int i=0;i<numberRows_;i++) {
      int pivot = pivotVariable_[i];
      if (pivot<numberColumns_)
	vec[i] = array[i] * columnScale_[pivot];
      else
	vec[i] = - array[i] / rowScale_[pivot-numberColumns_];
    }
  }
  rowArray1->clear();
}

//Get a column of the basis inverse
void 
ClpSimplex::getBInvCol(int col, double* vec)
{
  if (!rowArray_[0]) {
    printf("ClpSimplexPrimal or ClpSimplexDual must have been called with correct startFinishOption\n");
    abort();
  }
  ClpFactorization * factorization = factorization_;
  CoinIndexedVector * rowArray0 = rowArray(0);
  CoinIndexedVector * rowArray1 = rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
#ifndef NDEBUG
  int n = numberRows();
  if (col<0||col>=n) {
    indexError(col,"getBInvCol");
  }
#endif
  // put +1 in row
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double value = (pivotVariable_[col]<numberColumns_) ? 1.0 : -1.0;
  // What about scaling ?
  rowArray1->insert(col,value);
  factorization->updateColumn(rowArray0,rowArray1,false);
  memcpy(vec,rowArray1->denseVector(),numberRows()*sizeof(double));
  rowArray1->clear();
}

/* Get basic indices (order of indices corresponds to the
   order of elements in a vector retured by getBInvACol() and
   getBInvCol()).
*/
void 
ClpSimplex::getBasics(int* index)
{
  if (!rowArray_[0]) {
    printf("ClpSimplexPrimal or ClpSimplexDual must have been called with correct startFinishOption\n");
    abort();
  }
  CoinAssert (index);
  memcpy(index,pivotVariable(),
	 numberRows()*sizeof(int));
}
/* Set an objective function coefficient */
void 
ClpSimplex::setObjectiveCoefficient( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  if (elementIndex<0||elementIndex>=numberColumns_) {
    indexError(elementIndex,"setObjectiveCoefficient");
  }
#endif
  if (objective()[elementIndex] != elementValue) {
    objective()[elementIndex] = elementValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~64;
      double direction = optimizationDirection_*objectiveScale_;
      if (!rowScale_) {
	objectiveWork_[elementIndex] = direction*elementValue;
      } else {
	objectiveWork_[elementIndex] = direction*elementValue
	  *columnScale_[elementIndex];
      }
    }
  }
}
/* Set a single row lower bound<br>
   Use -DBL_MAX for -infinity. */
void 
ClpSimplex::setRowLower( int elementIndex, double elementValue ) {
#ifndef NDEBUG
  int n = numberRows_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowLower");
  }
#endif
  if (elementValue<-1.0e27)
    elementValue=-COIN_DBL_MAX;
  if (rowLower_[elementIndex] != elementValue) {
    rowLower_[elementIndex] = elementValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~16;
      if (rowLower_[elementIndex]==-COIN_DBL_MAX) {
	rowLowerWork_[elementIndex] = -COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowLowerWork_[elementIndex] = elementValue * rhsScale_;
      } else {
	rowLowerWork_[elementIndex] = elementValue * rhsScale_ 
	  * rowScale_[elementIndex];
      }
    }
  }
}
      
/* Set a single row upper bound<br>
   Use DBL_MAX for infinity. */
void 
ClpSimplex::setRowUpper( int elementIndex, double elementValue ) {
#ifndef NDEBUG
  int n = numberRows_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowUpper");
  }
#endif
  if (elementValue>1.0e27)
    elementValue=COIN_DBL_MAX;
  if (rowUpper_[elementIndex] != elementValue) {
    rowUpper_[elementIndex] = elementValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~32;
      if (rowUpper_[elementIndex]==COIN_DBL_MAX) {
	rowUpperWork_[elementIndex] = COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowUpperWork_[elementIndex] = elementValue * rhsScale_;
      } else {
	rowUpperWork_[elementIndex] = elementValue * rhsScale_ 
	  * rowScale_[elementIndex];
      }
    }
  }
}
    
/* Set a single row lower and upper bound */
void 
ClpSimplex::setRowBounds( int elementIndex,
	      double lowerValue, double upperValue ) {
#ifndef NDEBUG
  int n = numberRows_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowBounds");
  }
#endif
  if (lowerValue<-1.0e27)
    lowerValue=-COIN_DBL_MAX;
  if (upperValue>1.0e27)
    upperValue=COIN_DBL_MAX;
  CoinAssert (upperValue>=lowerValue);
  if (rowLower_[elementIndex] != lowerValue) {
    rowLower_[elementIndex] = lowerValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~16;
      if (rowLower_[elementIndex]==-COIN_DBL_MAX) {
	rowLowerWork_[elementIndex] = -COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowLowerWork_[elementIndex] = lowerValue * rhsScale_;
      } else {
	rowLowerWork_[elementIndex] = lowerValue * rhsScale_ 
	  * rowScale_[elementIndex];
      }
    }
  }
  if (rowUpper_[elementIndex] != upperValue) {
    rowUpper_[elementIndex] = upperValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~32;
      if (rowUpper_[elementIndex]==COIN_DBL_MAX) {
	rowUpperWork_[elementIndex] = COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowUpperWork_[elementIndex] = upperValue * rhsScale_;
      } else {
	rowUpperWork_[elementIndex] = upperValue * rhsScale_ 
	  * rowScale_[elementIndex];
      }
    }
  }
}
void ClpSimplex::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
#ifndef NDEBUG
  int n = numberRows_;
#endif
  int numberChanged=0;
  const int * saveFirst=indexFirst;
  while (indexFirst != indexLast) {
    const int iRow=*indexFirst++;
#ifndef NDEBUG
    if (iRow<0||iRow>=n) {
      indexError(iRow,"setRowSetBounds");
    }
#endif
    double lowerValue = *boundList++;
    double upperValue = *boundList++;
    if (lowerValue<-1.0e27)
      lowerValue=-COIN_DBL_MAX;
    if (upperValue>1.0e27)
      upperValue=COIN_DBL_MAX;
    CoinAssert (upperValue>=lowerValue);
    if (rowLower_[iRow] != lowerValue) {
      rowLower_[iRow] = lowerValue;
      whatsChanged_ &= ~16;
      numberChanged++;
    }
    if (rowUpper_[iRow] != upperValue) {
      rowUpper_[iRow] = upperValue;
      whatsChanged_ &= ~32;
      numberChanged++;
    }
  }
  if (numberChanged&&(whatsChanged_&1)!=0) {
    indexFirst = saveFirst;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      if (rowLower_[iRow]==-COIN_DBL_MAX) {
	rowLowerWork_[iRow] = -COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowLowerWork_[iRow] = rowLower_[iRow] * rhsScale_;
      } else {
	rowLowerWork_[iRow] = rowLower_[iRow] * rhsScale_ 
	  * rowScale_[iRow];
      }
      if (rowUpper_[iRow]==COIN_DBL_MAX) {
	rowUpperWork_[iRow] = COIN_DBL_MAX;
      } else if (!rowScale_) {
	rowUpperWork_[iRow] = rowUpper_[iRow] * rhsScale_;
      } else {
	rowUpperWork_[iRow] = rowUpper_[iRow] * rhsScale_ 
	  * rowScale_[iRow];
      }
    }
  }
}
//-----------------------------------------------------------------------------
/* Set a single column lower bound<br>
   Use -DBL_MAX for -infinity. */
void 
ClpSimplex::setColumnLower( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  int n = numberColumns_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColumnLower");
  }
#endif
  if (elementValue<-1.0e27)
    elementValue=-COIN_DBL_MAX;
  if (columnLower_[elementIndex] != elementValue) {
    columnLower_[elementIndex] = elementValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~128;
      if (columnLower_[elementIndex]==-COIN_DBL_MAX) {
	columnLowerWork_[elementIndex] = -COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnLowerWork_[elementIndex] = elementValue * rhsScale_;
      } else {
	columnLowerWork_[elementIndex] = elementValue * rhsScale_ 
	  / columnScale_[elementIndex];
      }
    }
  }
}
      
/* Set a single column upper bound<br>
   Use DBL_MAX for infinity. */
void 
ClpSimplex::setColumnUpper( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  int n = numberColumns_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColumnUpper");
  }
#endif
  if (elementValue>1.0e27)
    elementValue=COIN_DBL_MAX;
  if (columnUpper_[elementIndex] != elementValue) {
    columnUpper_[elementIndex] = elementValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~256;
      if (columnUpper_[elementIndex]==COIN_DBL_MAX) {
	columnUpperWork_[elementIndex] = COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnUpperWork_[elementIndex] = elementValue * rhsScale_;
      } else {
	columnUpperWork_[elementIndex] = elementValue * rhsScale_ 
	  / columnScale_[elementIndex];
      }
    }
  }
}

/* Set a single column lower and upper bound */
void 
ClpSimplex::setColumnBounds( int elementIndex,
				     double lowerValue, double upperValue )
{
#ifndef NDEBUG
  int n = numberColumns_;
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColumnBounds");
  }
#endif
  if (lowerValue<-1.0e27)
    lowerValue=-COIN_DBL_MAX;
  if (columnLower_[elementIndex] != lowerValue) {
    columnLower_[elementIndex] = lowerValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~128;
      if (columnLower_[elementIndex]==-COIN_DBL_MAX) {
	columnLowerWork_[elementIndex] = -COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnLowerWork_[elementIndex] = lowerValue * rhsScale_;
      } else {
	columnLowerWork_[elementIndex] = lowerValue * rhsScale_ 
	  / columnScale_[elementIndex];
      }
    }
  }
  if (upperValue>1.0e27)
    upperValue=COIN_DBL_MAX;
  CoinAssert (upperValue>=lowerValue);
  if (columnUpper_[elementIndex] != upperValue) {
    columnUpper_[elementIndex] = upperValue;
    if ((whatsChanged_&1)!=0) {
      // work arrays exist - update as well
      whatsChanged_ &= ~256;
      if (columnUpper_[elementIndex]==COIN_DBL_MAX) {
	columnUpperWork_[elementIndex] = COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnUpperWork_[elementIndex] = upperValue * rhsScale_;
      } else {
	columnUpperWork_[elementIndex] = upperValue * rhsScale_ 
	  / columnScale_[elementIndex];
      }
    }
  }
}
void ClpSimplex::setColumnSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
#ifndef NDEBUG
  int n = numberColumns_;
#endif
  int numberChanged=0;
  const int * saveFirst=indexFirst;
  while (indexFirst != indexLast) {
    const int iColumn=*indexFirst++;
#ifndef NDEBUG
    if (iColumn<0||iColumn>=n) {
      indexError(iColumn,"setColumnSetBounds");
    }
#endif
    double lowerValue = *boundList++;
    double upperValue = *boundList++;
    if (lowerValue<-1.0e27)
      lowerValue=-COIN_DBL_MAX;
    if (upperValue>1.0e27)
      upperValue=COIN_DBL_MAX;
    CoinAssert (upperValue>=lowerValue);
    if (columnLower_[iColumn] != lowerValue) {
      columnLower_[iColumn] = lowerValue;
      whatsChanged_ &= ~16;
      numberChanged++;
    }
    if (columnUpper_[iColumn] != upperValue) {
      columnUpper_[iColumn] = upperValue;
      whatsChanged_ &= ~32;
      numberChanged++;
    }
  }
  if (numberChanged&&(whatsChanged_&1)!=0) {
    indexFirst = saveFirst;
    while (indexFirst != indexLast) {
      const int iColumn=*indexFirst++;
      if (columnLower_[iColumn]==-COIN_DBL_MAX) {
	columnLowerWork_[iColumn] = -COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnLowerWork_[iColumn] = columnLower_[iColumn] * rhsScale_;
      } else {
	columnLowerWork_[iColumn] = columnLower_[iColumn] * rhsScale_ 
	  / columnScale_[iColumn];
      }
      if (columnUpper_[iColumn]==COIN_DBL_MAX) {
	columnUpperWork_[iColumn] = COIN_DBL_MAX;
      } else if (!columnScale_) {
	columnUpperWork_[iColumn] = columnUpper_[iColumn] * rhsScale_;
      } else {
	columnUpperWork_[iColumn] = columnUpper_[iColumn] * rhsScale_ 
	  / columnScale_[iColumn];
      }
    }
  }
}
/* Just check solution (for internal use) - sets sum of
   infeasibilities etc. */
void 
ClpSimplex::checkSolutionInternal()
{
  double dualTolerance=dblParam_[ClpDualTolerance];
  double primalTolerance=dblParam_[ClpPrimalTolerance];
  double nonLinearOffset=0.0;
  const double * objective = objective_->gradient(this,columnActivity_,
                                                  nonLinearOffset,true);
  int iRow,iColumn;
  assert (!rowObjective_);

  objectiveValue_ = 0.0;
  // now look at solution
  sumPrimalInfeasibilities_=0.0;
  numberPrimalInfeasibilities_=0;

  sumDualInfeasibilities_=0.0;
  numberDualInfeasibilities_=0;
  double maxmin = optimizationDirection_;

  for (iRow=0;iRow<numberRows_;iRow++) {
    double dualValue = dual_[iRow]*maxmin;
    double primalValue = rowActivity_[iRow];
    double lower = rowLower_[iRow];
    double upper = rowUpper_[iRow];
    if (primalValue>upper+primalTolerance) {
      sumPrimalInfeasibilities_ += primalValue-upper-primalTolerance;
      numberPrimalInfeasibilities_ ++;
    } else if (primalValue<lower-primalTolerance) {
      sumPrimalInfeasibilities_ += lower-primalValue-primalTolerance;
      numberPrimalInfeasibilities_ ++;
    } else if (getRowStatus(iRow) != basic) {
      // not basic
      if (primalValue<upper-primalTolerance) {
	// dual should not be negative
	if (dualValue<-dualTolerance) {
          sumDualInfeasibilities_ -= dualValue+dualTolerance_;
          numberDualInfeasibilities_ ++;
        }
      }
      if (primalValue>lower+primalTolerance) {
	// dual should not be positive
	if (dualValue>dualTolerance) {
          sumDualInfeasibilities_ += dualValue-dualTolerance_;
          numberDualInfeasibilities_ ++;
        }
      }
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double dualValue = reducedCost_[iColumn]*maxmin;
    double primalValue = columnActivity_[iColumn];
    objectiveValue_ += objective[iColumn]*primalValue;
    double lower = columnLower_[iColumn];
    double upper = columnUpper_[iColumn];
    if (primalValue>upper+primalTolerance) {
      sumPrimalInfeasibilities_ += primalValue-upper-primalTolerance;
      numberPrimalInfeasibilities_ ++;
    } else if (primalValue<lower-primalTolerance) {
      sumPrimalInfeasibilities_ += lower-primalValue-primalTolerance;
      numberPrimalInfeasibilities_ ++;
    } else if (getColumnStatus(iColumn) != basic) {
      // not basic
      if (primalValue<upper-primalTolerance) {
	// dual should not be negative
	if (dualValue<-dualTolerance) {
          sumDualInfeasibilities_ -= dualValue+dualTolerance_;
          numberDualInfeasibilities_ ++;
        }
      }
      if (primalValue>lower+primalTolerance) {
	// dual should not be positive
	if (dualValue>dualTolerance) {
          sumDualInfeasibilities_ += dualValue-dualTolerance_;
          numberDualInfeasibilities_ ++;
        }
      }
    }
  }
  objectiveValue_ += objective_->nonlinearOffset();

  if (!numberDualInfeasibilities_&&
      !numberPrimalInfeasibilities_)
    problemStatus_=0;
  else
    problemStatus_=-1;
}
/*
  When scaling is on it is possible that the scaled problem
  is feasible but the unscaled is not.  Clp returns a secondary
  status code to that effect.  This option allows for a cleanup.
  If you use it I would suggest 1.
  This only affects actions when scaled optimal
  0 - no action
  1 - clean up using dual if primal infeasibility
  2 - clean up using dual if dual infeasibility
  3 - clean up using dual if primal or dual infeasibility
  11,12,13 - as 1,2,3 but use primal
*/
#ifdef COUNT_CLEANUPS
static int n1=0;
static int n2=0;
static int n3=0;
#endif
int
ClpSimplex::cleanup(int cleanupScaling)
{
#ifdef COUNT_CLEANUPS
  n1++;
#endif
  int returnCode=0;
  if (!problemStatus_&&cleanupScaling) {
    int check = cleanupScaling%10;
    bool primal = (secondaryStatus_==2||secondaryStatus_==4);
    bool dual = (secondaryStatus_==3||secondaryStatus_==4);
    if (((check&1)!=0&&primal)||((check&2)!=0)&&dual) {
      // need cleanup
      int saveScalingFlag = scalingFlag_;
      // say matrix changed
      whatsChanged_ |= 1;
      scaling(0);
      if (cleanupScaling<10) {
        // dual
        returnCode = this->dual();
      } else {
        // primal
        returnCode = this->primal();
      }
#ifdef COUNT_CLEANUPS
      n2++;
      n3+= numberIterations_;
      //printf("**cleanup took %d iterations\n",numberIterations_);
#endif
      scaling(saveScalingFlag);
    }
  }
  return returnCode;
}
#ifdef COUNT_CLEANUPS
void printHowMany()
{
  printf("There were %d cleanups out of %d solves and %d iterations\n",
         n2,n1,n3);
}
#endif
