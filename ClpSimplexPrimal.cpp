// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

/* Notes on implementation of primal simplex algorithm.

   When primal feasible(A):

   If dual feasible, we are optimal.  Otherwise choose an infeasible
   basic variable to enter basis from a bound (B).  We now need to find an 
   outgoing variable which will leave problem primal feasible so we get
   the column of the tableau corresponding to the incoming variable
   (with the correct sign depending if variable will go up or down).

   We now perform a ratio test to determine which outgoing variable will
   preserve primal feasibility (C).  If no variable found then problem
   is unbounded (in primal sense).  If there is a variable, we then 
   perform pivot and repeat.  Trivial?

   -------------------------------------------

   A) How do we get primal feasible?  All variables have fake costs
   outside their feasible region so it is trivial to declare problem
   feasible.  OSL did not have a phase 1/phase 2 approach but
   instead effectively put an extra cost on infeasible basic variables
   I am taking the same approach here, although it is generalized
   to allow for non-linear costs and dual information.

   In OSL, this weight was changed heuristically, here at present
   it is only increased if problem looks finished.  If problem is
   feasible I check for unboundedness.  If not unbounded we 
   could play with going into dual.  As long as weights increase
   any algorithm would be finite.
   
   B) Which incoming variable to choose is a virtual base class.
   For difficult problems steepest edge is preferred while for
   very easy (large) problems we will need partial scan.

   C) Sounds easy, but this is hardest part of algorithm.
      1) Instead of stopping at first choice, we may be able
      to allow that variable to go through bound and if objective
      still improving choose again.  These mini iterations can
      increase speed by orders of magnitude but we may need to
      go to more of a bucket choice of variable rather than looking
      at them one by one (for speed).
      2) Accuracy.  Basic infeasibilities may be less than
      tolerance.  Pivoting on these makes objective go backwards.
      OSL modified cost so a zero move was made, Gill et al
      modified so a strictly positive move was made.
      The two problems are that re-factorizations can
      change rinfeasibilities above and below tolerances and that when
      finished we need to reset costs and try again.
      3) Degeneracy.  Gill et al helps but may not be enough.  We
      may need more.  Also it can improve speed a lot if we perturb
      the rhs and bounds significantly.  

  References:
     Forrest and Goldfarb, Steepest-edge simplex algorithms for
       linear programming - Mathematical Programming 1992
     Forrest and Tomlin, Implementing the simplex method for
       the Optimization Subroutine Library - IBM Systems Journal 1992
     Gill, Murray, Saunders, Wright A Practical Anti-Cycling
       Procedure for Linear and Nonlinear Programming SOL report 1988


  TODO:
 
  a) Better recovery procedures.  At present I never check on forward
     progress.  There is checkpoint/restart with reducing 
     re-factorization frequency, but this is only on singular 
     factorizations.
  b) Fast methods for large easy problems (and also the option for
     the code to automatically choose which method).
  c) We need to be able to stop in various ways for OSI - this
     is fairly easy.

 */


#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexPrimal.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpPrimalColumnPivot.hpp"
#include "ClpMessage.hpp"
#include "ClpEventHandler.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream>
// primal 
int ClpSimplexPrimal::primal (int ifValuesPass , int startFinishOptions)
{

  /*
      Method

     It tries to be a single phase approach with a weight of 1.0 being
     given to getting optimal and a weight of infeasibilityCost_ being
     given to getting primal feasible.  In this version I have tried to
     be clever in a stupid way.  The idea of fake bounds in dual
     seems to work so the primal analogue would be that of getting
     bounds on reduced costs (by a presolve approach) and using
     these for being above or below feasible region.  I decided to waste
     memory and keep these explicitly.  This allows for non-linear
     costs!

     The code is designed to take advantage of sparsity so arrays are
     seldom zeroed out from scratch or gone over in their entirety.
     The only exception is a full scan to find incoming variable for 
     Dantzig row choice.  For steepest edge we keep an updated list 
     of dual infeasibilities (actually squares).  
     On easy problems we don't need full scan - just
     pick first reasonable.

     One problem is how to tackle degeneracy and accuracy.  At present
     I am using the modification of costs which I put in OSL and which was
     extended by Gill et al.  I am still not sure of the exact details.

     The flow of primal is three while loops as follows:

     while (not finished) {

       while (not clean solution) {

          Factorize and/or clean up solution by changing bounds so
	  primal feasible.  If looks finished check fake primal bounds.
	  Repeat until status is iterating (-1) or finished (0,1,2)

       }

       while (status==-1) {

         Iterate until no pivot in or out or time to re-factorize.

         Flow is:

         choose pivot column (incoming variable).  if none then
	 we are primal feasible so looks as if done but we need to
	 break and check bounds etc.

	 Get pivot column in tableau

         Choose outgoing row.  If we don't find one then we look
	 primal unbounded so break and check bounds etc.  (Also the
	 pivot tolerance is larger after any iterations so that may be
	 reason)

         If we do find outgoing row, we may have to adjust costs to
	 keep going forwards (anti-degeneracy).  Check pivot will be stable
	 and if unstable throw away iteration and break to re-factorize.
	 If minor error re-factorize after iteration.

	 Update everything (this may involve changing bounds on 
	 variables to stay primal feasible.

       }

     }

     At present we never check we are going forwards.  I overdid that in
     OSL so will try and make a last resort.

     Needs partial scan pivot in option.

     May need other anti-degeneracy measures, especially if we try and use
     loose tolerances as a way to solve in fewer iterations.

     I like idea of dynamic scaling.  This gives opportunity to decouple
     different implications of scaling for accuracy, iteration count and
     feasibility tolerance.

  */

  algorithm_ = +1;
  //specialOptions_ |= 4;

  // save data
  ClpDataSave data = saveData();
  matrix_->refresh(this); // make sure matrix okay
  
  // Save so can see if doing after dual
  int initialStatus=problemStatus_;
  // initialize - maybe values pass and algorithm_ is +1
  if (!startup(ifValuesPass)) {
    
    // Set average theta
    nonLinearCost_->setAverageTheta(1.0e3);
    int lastCleaned=0; // last time objective or bounds cleaned up
    
    // Say no pivot has occurred (for steepest edge and updates)
    pivotRow_=-2;
    
    // This says whether to restore things etc
    int factorType=0;
    if (problemStatus_<0&&perturbation_<100) {
      perturb(0);
      // Can't get here if values pass
      assert (!ifValuesPass);
      gutsOfSolution(NULL,NULL);
      if (handler_->logLevel()>2) {
	handler_->message(CLP_SIMPLEX_STATUS,messages_)
	  <<numberIterations_<<objectiveValue();
	handler_->printing(sumPrimalInfeasibilities_>0.0)
	  <<sumPrimalInfeasibilities_<<numberPrimalInfeasibilities_;
	handler_->printing(sumDualInfeasibilities_>0.0)
	  <<sumDualInfeasibilities_<<numberDualInfeasibilities_;
	handler_->printing(numberDualInfeasibilitiesWithoutFree_
			   <numberDualInfeasibilities_)
			     <<numberDualInfeasibilitiesWithoutFree_;
	handler_->message()<<CoinMessageEol;
      }
    }
    ClpSimplex * saveModel=NULL;
    int stopSprint=-1;
    int sprintPass=0;
    int reasonableSprintIteration=0;
    int lastSprintIteration=0;
    double lastObjectiveValue=COIN_DBL_MAX;
    // Start check for cycles
    progress_->startCheck();
    /*
      Status of problem:
      0 - optimal
      1 - infeasible
      2 - unbounded
      -1 - iterating
      -2 - factorization wanted
      -3 - redo checking without factorization
      -4 - looks infeasible
      -5 - looks unbounded
    */
    while (problemStatus_<0) {
      int iRow,iColumn;
      // clear
      for (iRow=0;iRow<4;iRow++) {
	rowArray_[iRow]->clear();
      }    
      
      for (iColumn=0;iColumn<2;iColumn++) {
	columnArray_[iColumn]->clear();
      }    
      
      // give matrix (and model costs and bounds a chance to be
      // refreshed (normally null)
      matrix_->refresh(this);
      // If getting nowhere - why not give it a kick
#if 1
      if (perturbation_<101&&numberIterations_>2*(numberRows_+numberColumns_)&&(specialOptions_&4)==0
	  &&initialStatus!=10) {
	perturb(1);
	matrix_->rhsOffset(this,true,false);
      }
#endif
      // If we have done no iterations - special
      if (lastGoodIteration_==numberIterations_&&factorType)
	factorType=3;
      if (saveModel) {
	// Doing sprint
	if (sequenceIn_<0||numberIterations_>=stopSprint) {
	  problemStatus_=-1;
	  originalModel(saveModel);
	  saveModel=NULL;
	  if (sequenceIn_<0&&numberIterations_<reasonableSprintIteration&&
	      sprintPass>100)
	    primalColumnPivot_->switchOffSprint();
	  //lastSprintIteration=numberIterations_;
	  printf("End small model\n");
	}
      }
	  
      // may factorize, checks if problem finished
      statusOfProblemInPrimal(lastCleaned,factorType,progress_,true,ifValuesPass,saveModel);
      // See if sprint says redo beacuse of problems
      if (numberDualInfeasibilities_==-776) {
	// Need new set of variables
	problemStatus_=-1;
	originalModel(saveModel);
	saveModel=NULL;
	//lastSprintIteration=numberIterations_;
	printf("End small model after\n");
	statusOfProblemInPrimal(lastCleaned,factorType,progress_,true,ifValuesPass,saveModel);
      } 
      int numberSprintIterations=0;
      int numberSprintColumns = primalColumnPivot_->numberSprintColumns(numberSprintIterations);
      if (problemStatus_==777) {
	// problems so do one pass with normal
	problemStatus_=-1;
	originalModel(saveModel);
	saveModel=NULL;
	// Skip factorization
	//statusOfProblemInPrimal(lastCleaned,factorType,progress_,false,saveModel);
	statusOfProblemInPrimal(lastCleaned,factorType,progress_,true,ifValuesPass,saveModel);
      } else if (problemStatus_<0&&!saveModel&&numberSprintColumns&&firstFree_<0) {
	int numberSort=0;
	int numberFixed=0;
	int numberBasic=0;
	reasonableSprintIteration = numberIterations_ + 100;
	int * whichColumns = new int[numberColumns_];
	double * weight = new double[numberColumns_];
	int numberNegative=0;
	double sumNegative = 0.0;
	// now massage weight so all basic in plus good djs
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  double dj = dj_[iColumn];
	  switch(getColumnStatus(iColumn)) {
	    
	  case basic:
	    dj = -1.0e50;
	    numberBasic++;
	    break;
	  case atUpperBound:
	    dj = -dj;
	    break;
	  case isFixed:
	    dj=1.0e50;
	    numberFixed++;
	    break;
	  case atLowerBound:
	    dj = dj;
	    break;
	  case isFree:
	    dj = -100.0*fabs(dj);
	      break;
	  case superBasic:
	    dj = -100.0*fabs(dj);
	    break;
	  }
	  if (dj<-dualTolerance_&&dj>-1.0e50) {
	    numberNegative++;
	    sumNegative -= dj;
	  }
	  weight[iColumn]=dj;
	  whichColumns[iColumn] = iColumn;
	}
	handler_->message(CLP_SPRINT,messages_)
	  <<sprintPass<<numberIterations_-lastSprintIteration<<objectiveValue()<<sumNegative
	  <<numberNegative
	  <<CoinMessageEol;
	sprintPass++;
	lastSprintIteration=numberIterations_;
	if (objectiveValue()*optimizationDirection_>lastObjectiveValue-1.0e-7&&sprintPass>5) {
	  // switch off
	  printf("Switching off sprint\n");
	  primalColumnPivot_->switchOffSprint();
	} else {
	  lastObjectiveValue = objectiveValue()*optimizationDirection_;
	  // sort
	  CoinSort_2(weight,weight+numberColumns_,whichColumns);
	  numberSort = CoinMin(numberColumns_-numberFixed,numberBasic+numberSprintColumns);
	  // Sort to make consistent ?
	  std::sort(whichColumns,whichColumns+numberSort);
	  saveModel = new ClpSimplex(this,numberSort,whichColumns);
	  delete [] whichColumns;
	  delete [] weight;
	  // Skip factorization
	  //statusOfProblemInPrimal(lastCleaned,factorType,progress_,false,saveModel);
	  //statusOfProblemInPrimal(lastCleaned,factorType,progress_,true,saveModel);
	  stopSprint = numberIterations_+numberSprintIterations;
	  printf("Sprint with %d columns for %d iterations\n",
		 numberSprintColumns,numberSprintIterations);
	}
      }
      
      // Say good factorization
      factorType=1;
      
      // Say no pivot has occurred (for steepest edge and updates)
      pivotRow_=-2;

      // exit if victory declared
      if (problemStatus_>=0)
	break;
      
      // test for maximum iterations
      if (hitMaximumIterations()||(ifValuesPass==2&&firstFree_<0)) {
	problemStatus_=3;
	break;
      }

      if (firstFree_<0) {
	if (ifValuesPass) {
	  // end of values pass
	  ifValuesPass=0;
	  int status = eventHandler_->event(ClpEventHandler::endOfValuesPass);
	  if (status>=0) {
	    problemStatus_=5;
	    secondaryStatus_=ClpEventHandler::endOfValuesPass;
	    break;
	  }
	}
      }
      // Check event
      {
	int status = eventHandler_->event(ClpEventHandler::endOfFactorization);
	if (status>=0) {
	  problemStatus_=5;
	  secondaryStatus_=ClpEventHandler::endOfFactorization;
	  break;
	}
      }
      // Iterate
      whileIterating(ifValuesPass ? 1 : 0);
    }
  }
  // if infeasible get real values
  if (problemStatus_==1) {
    infeasibilityCost_=0.0;
    createRim(7);
    nonLinearCost_->checkInfeasibilities(0.0);
    sumPrimalInfeasibilities_=nonLinearCost_->sumInfeasibilities();
    numberPrimalInfeasibilities_= nonLinearCost_->numberInfeasibilities();
    // and get good feasible duals
    computeDuals(NULL);
  }
  // clean up
  unflag();
  if (!startFinishOptions)
    finish();
  restoreData(data);
  return problemStatus_;
}
/*
  Reasons to come out:
  -1 iterations etc
  -2 inaccuracy 
  -3 slight inaccuracy (and done iterations)
  -4 end of values pass and done iterations
  +0 looks optimal (might be infeasible - but we will investigate)
  +2 looks unbounded
  +3 max iterations 
*/
int
ClpSimplexPrimal::whileIterating(int valuesOption)
{
  // Say if values pass
  int ifValuesPass=(firstFree_>=0) ? 1 : 0;
  int returnCode=-1;
  int superBasicType=1;
  if (valuesOption>1)
    superBasicType=3;
  // status stays at -1 while iterating, >=0 finished, -2 to invert
  // status -3 to go to top without an invert
  while (problemStatus_==-1) {
    //#define CLP_DEBUG 1
#ifdef CLP_DEBUG
    {
      int i;
      // not [1] as has information
      for (i=0;i<4;i++) {
	if (i!=1)
	  rowArray_[i]->checkClear();
      }    
      for (i=0;i<2;i++) {
	columnArray_[i]->checkClear();
      }    
    }      
#endif
#if 0
    {
      int iPivot;
      double * array = rowArray_[3]->denseVector();
      int * index = rowArray_[3]->getIndices();
      int i;
      for (iPivot=0;iPivot<numberRows_;iPivot++) {
	int iSequence = pivotVariable_[iPivot];
	unpackPacked(rowArray_[3],iSequence);
	factorization_->updateColumn(rowArray_[2],rowArray_[3]);
	int number = rowArray_[3]->getNumElements();
	for (i=0;i<number;i++) {
	  int iRow = index[i];
	  if (iRow==iPivot)
	    assert (fabs(array[i]-1.0)<1.0e-4);
	  else
	    assert (fabs(array[i])<1.0e-4);
	}
	rowArray_[3]->clear();
      }
    }
#endif
#if CLP_DEBUG>2
    // very expensive
    if (numberIterations_>0&&numberIterations_<-2534) {
      handler_->setLogLevel(63);
      double saveValue = objectiveValue_;
      double * saveRow1 = new double[numberRows_];
      double * saveRow2 = new double[numberRows_];
      memcpy(saveRow1,rowReducedCost_,numberRows_*sizeof(double));
      memcpy(saveRow2,rowActivityWork_,numberRows_*sizeof(double));
      double * saveColumn1 = new double[numberColumns_];
      double * saveColumn2 = new double[numberColumns_];
      memcpy(saveColumn1,reducedCostWork_,numberColumns_*sizeof(double));
      memcpy(saveColumn2,columnActivityWork_,numberColumns_*sizeof(double));
      createRim(7);
      gutsOfSolution(NULL,NULL);
      printf("xxx %d old obj %g, recomputed %g, sum primal inf %g\n",
	     numberIterations_,
	     saveValue,objectiveValue_,sumPrimalInfeasibilities_);
      memcpy(rowReducedCost_,saveRow1,numberRows_*sizeof(double));
      memcpy(rowActivityWork_,saveRow2,numberRows_*sizeof(double));
      memcpy(reducedCostWork_,saveColumn1,numberColumns_*sizeof(double));
      memcpy(columnActivityWork_,saveColumn2,numberColumns_*sizeof(double));
      delete [] saveRow1;
      delete [] saveRow2;
      delete [] saveColumn1;
      delete [] saveColumn2;
      objectiveValue_=saveValue;
    }
#endif
    if (!ifValuesPass) {
      // choose column to come in
      // can use pivotRow_ to update weights
      // pass in list of cost changes so can do row updates (rowArray_[1])
      // NOTE rowArray_[0] is used by computeDuals which is a 
      // slow way of getting duals but might be used 
      primalColumn(rowArray_[1],rowArray_[2],rowArray_[3],
		   columnArray_[0],columnArray_[1]);
    } else {
      // in values pass
      int sequenceIn=nextSuperBasic(superBasicType,columnArray_[0]);
      if (valuesOption>1)
	superBasicType=2;
      if (sequenceIn<0) {
	// end of values pass - initialize weights etc
	handler_->message(CLP_END_VALUES_PASS,messages_)
	  <<numberIterations_;
	primalColumnPivot_->saveWeights(this,5);
	problemStatus_=-2; // factorize now
	pivotRow_=-1; // say no weights update
	returnCode=-4;
	// Clean up
	int i;
	for (i=0;i<numberRows_+numberColumns_;i++) {
	  if (getColumnStatus(i)==atLowerBound||getColumnStatus(i)==isFixed)
	    solution_[i]=lower_[i];
	  else if (getColumnStatus(i)==atUpperBound)
	    solution_[i]=upper_[i];
	}
	break;
      } else {
	// normal
	sequenceIn_ = sequenceIn;
	valueIn_=solution_[sequenceIn_];
	lowerIn_=lower_[sequenceIn_];
	upperIn_=upper_[sequenceIn_];
	dualIn_=dj_[sequenceIn_];
      }
    }
    pivotRow_=-1;
    sequenceOut_=-1;
    rowArray_[1]->clear();
    if (sequenceIn_>=0) {
      // we found a pivot column
      assert (!flagged(sequenceIn_));
#ifdef CLP_DEBUG
      if ((handler_->logLevel()&32)) {
	char x = isColumn(sequenceIn_) ? 'C' :'R';
	std::cout<<"pivot column "<<
	  x<<sequenceWithin(sequenceIn_)<<std::endl;
      }
#endif
      // do second half of iteration
      returnCode = pivotResult(ifValuesPass);
      if (returnCode<-1&&returnCode>-5) {
	problemStatus_=-2; // 
      } else if (returnCode==-5) {
	// something flagged - continue;
      } else if (returnCode==2) {
	problemStatus_=-5; // looks unbounded
      } else if (returnCode==4) {
	problemStatus_=-2; // looks unbounded but has iterated
      } else if (returnCode!=-1) {
	assert(returnCode==3);
	problemStatus_=3;
      }
    } else {
      // no pivot column
#ifdef CLP_DEBUG
      if (handler_->logLevel()&32)
	printf("** no column pivot\n");
#endif
      if (nonLinearCost_->numberInfeasibilities())
	problemStatus_=-4; // might be infeasible 
      // Force to re-factorize early next time
      int numberPivots = factorization_->pivots();
      forceFactorization_=CoinMin(forceFactorization_,(numberPivots+1)>>1);
      returnCode=0;
      break;
    }
  }
  if (valuesOption>1) 
    columnArray_[0]->setNumElements(0);
  return returnCode;
}
/* Checks if finished.  Updates status */
void 
ClpSimplexPrimal::statusOfProblemInPrimal(int & lastCleaned,int type,
					  ClpSimplexProgress * progress,
					  bool doFactorization,
					  int ifValuesPass,
					  ClpSimplex * originalModel)
{
  int dummy; // for use in generalExpanded
  if (type==2) {
    // trouble - restore solution
    memcpy(status_ ,saveStatus_,(numberColumns_+numberRows_)*sizeof(char));
    memcpy(rowActivityWork_,savedSolution_+numberColumns_ ,
	   numberRows_*sizeof(double));
    memcpy(columnActivityWork_,savedSolution_ ,
	   numberColumns_*sizeof(double));
    // restore extra stuff
    matrix_->generalExpanded(this,6,dummy);
    forceFactorization_=1; // a bit drastic but ..
    pivotRow_=-1; // say no weights update
    changeMade_++; // say change made
  }
  int saveThreshold = factorization_->sparseThreshold();
  int tentativeStatus = problemStatus_;
  int numberThrownOut=1; // to loop round on bad factorization in values pass
  double lastSumInfeasibility=COIN_DBL_MAX;
  if (numberIterations_)
    lastSumInfeasibility=nonLinearCost_->sumInfeasibilities();
  while (numberThrownOut) {
    if (problemStatus_>-3||problemStatus_==-4) {
      // factorize
      // later on we will need to recover from singularities
      // also we could skip if first time
      // do weights
      // This may save pivotRow_ for use 
      if (doFactorization)
	primalColumnPivot_->saveWeights(this,1);
      
      if (type&&doFactorization) {
	// is factorization okay?
	int factorStatus = internalFactorize(1);
	if (factorStatus) {
	  if (solveType_==2+8) {
	    // say odd
	    problemStatus_=5;
	    return;
	  }
	  if (type!=1||largestPrimalError_>1.0e3
	      ||largestDualError_>1.0e3) {
	    // switch off dense
	    int saveDense = factorization_->denseThreshold();
	    factorization_->setDenseThreshold(0);
	    // Go to safe
	    factorization_->pivotTolerance(0.99);
	    // make sure will do safe factorization
	    pivotVariable_[0]=-1;
	    internalFactorize(2);
	    factorization_->setDenseThreshold(saveDense);
	    // restore extra stuff
	    matrix_->generalExpanded(this,6,dummy);
	  } else {
	    // no - restore previous basis
	    memcpy(status_ ,saveStatus_,(numberColumns_+numberRows_)*sizeof(char));
	    memcpy(rowActivityWork_,savedSolution_+numberColumns_ ,
		   numberRows_*sizeof(double));
	    memcpy(columnActivityWork_,savedSolution_ ,
		   numberColumns_*sizeof(double));
	    // restore extra stuff
	    matrix_->generalExpanded(this,6,dummy);
	    matrix_->generalExpanded(this,5,dummy);
	    forceFactorization_=1; // a bit drastic but ..
	    type = 2;
	    // Go to safe 
	    factorization_->pivotTolerance(0.99);
	    if (internalFactorize(1)!=0)
	       largestPrimalError_=1.0e4; // force other type
	  }
	  changeMade_++; // say change made
	}
      }
      if (problemStatus_!=-4)
	problemStatus_=-3;
    }
    // at this stage status is -3 or -5 if looks unbounded
    // get primal and dual solutions
    // put back original costs and then check
    createRim(4);
    // May need to do more if column generation
    dummy=4;
    matrix_->generalExpanded(this,9,dummy);
    numberThrownOut=gutsOfSolution(NULL,NULL,(firstFree_>=0));
    double sumInfeasibility =  nonLinearCost_->sumInfeasibilities();
    if (numberThrownOut||
	(sumInfeasibility>1.0e7&&sumInfeasibility>100.0*lastSumInfeasibility
	 &&factorization_->pivotTolerance()<0.11)) {
      problemStatus_=tentativeStatus;
      doFactorization=true;
    }
  }
  // Double check reduced costs if no action
  if (progress->lastIterationNumber(0)==numberIterations_) {
    if (primalColumnPivot_->looksOptimal()) {
      numberDualInfeasibilities_ = 0;
      sumDualInfeasibilities_ = 0.0;
    }
  }
  // Check if looping
  int loop;
  if (type!=2&&!ifValuesPass) 
    loop = progress->looping();
  else
    loop=-1;
  if (loop>=0) {
    if (!problemStatus_) {
      // declaring victory
      numberPrimalInfeasibilities_ = 0;
      sumPrimalInfeasibilities_ = 0.0;
    } else {
      problemStatus_ = loop; //exit if in loop 
      problemStatus_ = 10; // instead - try other algorithm
    }
    problemStatus_ = 10; // instead - try other algorithm
    return ;
  } else if (loop<-1) {
    // Is it time for drastic measures
    if (nonLinearCost_->numberInfeasibilities()&&progress->badTimes()>5&&
	progress->oddState()<10&&progress->oddState()>=0) {
      progress->newOddState();
      nonLinearCost_->zapCosts();
    }
    // something may have changed
    gutsOfSolution(NULL,NULL,ifValuesPass!=0);
  }
  // If progress then reset costs
  if (loop==-1&&!nonLinearCost_->numberInfeasibilities()&&progress->oddState()<0) {
    createRim(4,false); // costs back
    delete nonLinearCost_;
    nonLinearCost_ = new ClpNonLinearCost(this);
    progress->endOddState();
    gutsOfSolution(NULL,NULL,ifValuesPass!=0);
  }
  // Flag to say whether to go to dual to clean up
  bool goToDual=false;
  // really for free variables in
  //if((progressFlag_&2)!=0)
  //problemStatus_=-1;;
  progressFlag_ = 0; //reset progress flag

  handler_->message(CLP_SIMPLEX_STATUS,messages_)
    <<numberIterations_<<nonLinearCost_->feasibleReportCost();
  handler_->printing(nonLinearCost_->numberInfeasibilities()>0)
    <<nonLinearCost_->sumInfeasibilities()<<nonLinearCost_->numberInfeasibilities();
  handler_->printing(sumDualInfeasibilities_>0.0)
    <<sumDualInfeasibilities_<<numberDualInfeasibilities_;
  handler_->printing(numberDualInfeasibilitiesWithoutFree_
		     <numberDualInfeasibilities_)
		       <<numberDualInfeasibilitiesWithoutFree_;
  handler_->message()<<CoinMessageEol;
  if (!primalFeasible()) {
    nonLinearCost_->checkInfeasibilities(primalTolerance_);
    gutsOfSolution(NULL,NULL,ifValuesPass!=0);
    nonLinearCost_->checkInfeasibilities(primalTolerance_);
  }
  double trueInfeasibility =nonLinearCost_->sumInfeasibilities();
  if (trueInfeasibility>1.0) {
    // If infeasibility going up may change weights
    double testValue = trueInfeasibility-1.0e-4*(10.0+trueInfeasibility);
    if(progress->lastInfeasibility()<testValue) {
      if (infeasibilityCost_<1.0e14) {
	infeasibilityCost_ *= 1.5;
	printf("increasing weight to %g\n",infeasibilityCost_);
	gutsOfSolution(NULL,NULL,ifValuesPass!=0);
      }
    }
  }
  // we may wish to say it is optimal even if infeasible
  bool alwaysOptimal = (specialOptions_&1)!=0;
  // give code benefit of doubt
  if (sumOfRelaxedDualInfeasibilities_ == 0.0&&
      sumOfRelaxedPrimalInfeasibilities_ == 0.0) {
    // say optimal (with these bounds etc)
    numberDualInfeasibilities_ = 0;
    sumDualInfeasibilities_ = 0.0;
    numberPrimalInfeasibilities_ = 0;
    sumPrimalInfeasibilities_ = 0.0;
    // But check if in sprint
    if (originalModel) {
      // Carry on and re-do
      numberDualInfeasibilities_ = -776;
    }
  }
  // had ||(type==3&&problemStatus_!=-5) -- ??? why ????
  if ((dualFeasible()||problemStatus_==-4)&&!ifValuesPass) {
    // see if extra helps
    if (nonLinearCost_->numberInfeasibilities()&&
	 (nonLinearCost_->sumInfeasibilities()>1.0e-3||sumOfRelaxedPrimalInfeasibilities_)
	&&!alwaysOptimal) {
      //may need infeasiblity cost changed
      // we can see if we can construct a ray
      // make up a new objective
      double saveWeight = infeasibilityCost_;
      // save nonlinear cost as we are going to switch off costs
      ClpNonLinearCost * nonLinear = nonLinearCost_;
      // do twice to make sure Primal solution has settled
      // put non-basics to bounds in case tolerance moved
      // put back original costs
      createRim(4);
      nonLinearCost_->checkInfeasibilities(0.0);
      gutsOfSolution(NULL,NULL,ifValuesPass!=0);

      infeasibilityCost_=1.0e100;
      // put back original costs
      createRim(4);
      nonLinearCost_->checkInfeasibilities(primalTolerance_);
      // may have fixed infeasibilities - double check
      if (nonLinearCost_->numberInfeasibilities()==0) {
	// carry on
	problemStatus_ = -1;
	infeasibilityCost_=saveWeight;
	nonLinearCost_->checkInfeasibilities(primalTolerance_);
      } else {
	nonLinearCost_=NULL;
	// scale
	int i;
	for (i=0;i<numberRows_+numberColumns_;i++) 
	  cost_[i] *= 1.0e-95;
	gutsOfSolution(NULL,NULL,ifValuesPass!=0);
	nonLinearCost_=nonLinear;
	infeasibilityCost_=saveWeight;
	if ((infeasibilityCost_>=1.0e18||
	     numberDualInfeasibilities_==0)&&perturbation_==101) {
	  goToDual=unPerturb(); // stop any further perturbation
	  if (nonLinearCost_->sumInfeasibilities()>1.0e-1)
	    goToDual=false;
	  nonLinearCost_->checkInfeasibilities(primalTolerance_);
	  numberDualInfeasibilities_=1; // carry on
	  problemStatus_=-1;
	}
	if (infeasibilityCost_>=1.0e20||
	    numberDualInfeasibilities_==0) {
	  // we are infeasible - use as ray
	  delete [] ray_;
	  ray_ = new double [numberRows_];
	  memcpy(ray_,dual_,numberRows_*sizeof(double));
	  // and get feasible duals
	  infeasibilityCost_=0.0;
	  createRim(4);
	  nonLinearCost_->checkInfeasibilities(primalTolerance_);
	  gutsOfSolution(NULL,NULL,ifValuesPass!=0);
	  // so will exit
	  infeasibilityCost_=1.0e30;
	  // reset infeasibilities
	  sumPrimalInfeasibilities_=nonLinearCost_->sumInfeasibilities();;
	  numberPrimalInfeasibilities_=
	    nonLinearCost_->numberInfeasibilities();
	}
	if (infeasibilityCost_<1.0e20) {
	  infeasibilityCost_ *= 5.0;
	  changeMade_++; // say change made
	  handler_->message(CLP_PRIMAL_WEIGHT,messages_)
	    <<infeasibilityCost_
	    <<CoinMessageEol;
	  // put back original costs and then check
	  createRim(4);
	  nonLinearCost_->checkInfeasibilities(0.0);
	  gutsOfSolution(NULL,NULL,ifValuesPass!=0);
	  problemStatus_=-1; //continue
	  goToDual=false;
	} else {
	  // say infeasible
	  problemStatus_ = 1;
	}
      }
    } else {
      // may be optimal
      if (perturbation_==101) {
	goToDual=unPerturb(); // stop any further perturbation
	lastCleaned=-1; // carry on
      }
      bool unflagged = unflag();
      if ( lastCleaned!=numberIterations_||unflagged) {
	handler_->message(CLP_PRIMAL_OPTIMAL,messages_)
	  <<primalTolerance_
	  <<CoinMessageEol;
	if (numberTimesOptimal_<4) {
	  numberTimesOptimal_++;
	  changeMade_++; // say change made
	  if (numberTimesOptimal_==1) {
	    // better to have small tolerance even if slower
	    factorization_->zeroTolerance(1.0e-15);
	  }
	  lastCleaned=numberIterations_;
	  if (primalTolerance_!=dblParam_[ClpPrimalTolerance])
	    handler_->message(CLP_PRIMAL_ORIGINAL,messages_)
	      <<CoinMessageEol;
	  double oldTolerance = primalTolerance_;
	  primalTolerance_=dblParam_[ClpPrimalTolerance];
#if 0
	  double * xcost = new double[numberRows_+numberColumns_];
	  double * xlower = new double[numberRows_+numberColumns_];
	  double * xupper = new double[numberRows_+numberColumns_];
	  double * xdj = new double[numberRows_+numberColumns_];
	  double * xsolution = new double[numberRows_+numberColumns_];
	  memcpy(xcost,cost_,(numberRows_+numberColumns_)*sizeof(double));
	  memcpy(xlower,lower_,(numberRows_+numberColumns_)*sizeof(double));
	  memcpy(xupper,upper_,(numberRows_+numberColumns_)*sizeof(double));
	  memcpy(xdj,dj_,(numberRows_+numberColumns_)*sizeof(double));
	  memcpy(xsolution,solution_,(numberRows_+numberColumns_)*sizeof(double));
#endif
	  // put back original costs and then check
	  createRim(4);
	  nonLinearCost_->checkInfeasibilities(oldTolerance);
#if 0
	  int i;
	  for (i=0;i<numberRows_+numberColumns_;i++) {
	    if (cost_[i]!=xcost[i])
	      printf("** %d old cost %g new %g sol %g\n",
		     i,xcost[i],cost_[i],solution_[i]);
	    if (lower_[i]!=xlower[i])
	      printf("** %d old lower %g new %g sol %g\n",
		     i,xlower[i],lower_[i],solution_[i]);
	    if (upper_[i]!=xupper[i])
	      printf("** %d old upper %g new %g sol %g\n",
		     i,xupper[i],upper_[i],solution_[i]);
	    if (dj_[i]!=xdj[i])
	      printf("** %d old dj %g new %g sol %g\n",
		     i,xdj[i],dj_[i],solution_[i]);
	    if (solution_[i]!=xsolution[i])
	      printf("** %d old solution %g new %g sol %g\n",
		     i,xsolution[i],solution_[i],solution_[i]);
	  }
	  delete [] xcost;
	  delete [] xupper;
	  delete [] xlower;
	  delete [] xdj;
	  delete [] xsolution;
#endif
	  gutsOfSolution(NULL,NULL,ifValuesPass!=0);
	  if (sumOfRelaxedDualInfeasibilities_ == 0.0&&
	      sumOfRelaxedPrimalInfeasibilities_ == 0.0) {
	    // say optimal (with these bounds etc)
	    numberDualInfeasibilities_ = 0;
	    sumDualInfeasibilities_ = 0.0;
	    numberPrimalInfeasibilities_ = 0;
	    sumPrimalInfeasibilities_ = 0.0;
	  }
	  if (dualFeasible()&&!nonLinearCost_->numberInfeasibilities()&&lastCleaned>=0)
	    problemStatus_=0;
	  else
	    problemStatus_ = -1;
	} else {
	  problemStatus_=0; // optimal
	  if (lastCleaned<numberIterations_) {
	    handler_->message(CLP_SIMPLEX_GIVINGUP,messages_)
	      <<CoinMessageEol;
	  }
	}
      } else {
	problemStatus_=0; // optimal
      }
    }
  } else {
    // see if looks unbounded
    if (problemStatus_==-5) {
      if (nonLinearCost_->numberInfeasibilities()) {
	if (infeasibilityCost_>1.0e18&&perturbation_==101) {
	  // back off weight
	  infeasibilityCost_ = 1.0e13;
	  unPerturb(); // stop any further perturbation
	}
	//we need infeasiblity cost changed
	if (infeasibilityCost_<1.0e20) {
	  infeasibilityCost_ *= 5.0;
	  changeMade_++; // say change made
	  handler_->message(CLP_PRIMAL_WEIGHT,messages_)
	    <<infeasibilityCost_
	    <<CoinMessageEol;
	  // put back original costs and then check
	  createRim(4);
	  gutsOfSolution(NULL,NULL,ifValuesPass!=0);
	  problemStatus_=-1; //continue
	} else {
	  // say unbounded
	  problemStatus_ = 2;
	}
      } else {
	// say unbounded
	problemStatus_ = 2;
      } 
    } else {
      if(type==3&&problemStatus_!=-5)
	unflag(); // odd
      // carry on
      problemStatus_ = -1;
    }
  }
  // save extra stuff
  matrix_->generalExpanded(this,5,dummy);
  if (type==0||type==1) {
    if (type!=1||!saveStatus_) {
      // create save arrays
      delete [] saveStatus_;
      delete [] savedSolution_;
      saveStatus_ = new unsigned char [numberRows_+numberColumns_];
      savedSolution_ = new double [numberRows_+numberColumns_];
    }
    // save arrays
    memcpy(saveStatus_,status_,(numberColumns_+numberRows_)*sizeof(char));
    memcpy(savedSolution_+numberColumns_ ,rowActivityWork_,
	   numberRows_*sizeof(double));
    memcpy(savedSolution_ ,columnActivityWork_,numberColumns_*sizeof(double));
  }
  if (doFactorization) {
    // restore weights (if saved) - also recompute infeasibility list
    if (tentativeStatus>-3) 
      primalColumnPivot_->saveWeights(this,(type <2) ? 2 : 4);
    else
      primalColumnPivot_->saveWeights(this,3);
    if (saveThreshold) {
      // use default at present
      factorization_->sparseThreshold(0);
      factorization_->goSparse();
    }
  }
  if (problemStatus_<0&&!changeMade_) {
    problemStatus_=4; // unknown
  }
  lastGoodIteration_ = numberIterations_;
  if (goToDual) 
    problemStatus_=10; // try dual
#if 0
  double thisObj = progress->lastObjective(0);
  double lastObj = progress->lastObjective(1);
  if (lastObj<thisObj-1.0e-7*max(fabs(thisObj),fabs(lastObj))-1.0e-8
      &&firstFree_<0) {
    int maxFactor = factorization_->maximumPivots();
    if (maxFactor>10) {
      if (forceFactorization_<0)
	forceFactorization_= maxFactor;
      forceFactorization_ = CoinMax(1,(forceFactorization_>>1));
      printf("Reducing factorization frequency\n");
    } 
  }
#endif
}
/* 
   Row array has pivot column
   This chooses pivot row.
   For speed, we may need to go to a bucket approach when many
   variables go through bounds
   On exit rhsArray will have changes in costs of basic variables
*/
void 
ClpSimplexPrimal::primalRow(CoinIndexedVector * rowArray,
			    CoinIndexedVector * rhsArray,
			    CoinIndexedVector * spareArray,
			    CoinIndexedVector * spareArray2,
			    int valuesPass)
{
  double saveDj = dualIn_;
  if (valuesPass&&objective_->type()<2) {
    dualIn_ = cost_[sequenceIn_];

    double * work=rowArray->denseVector();
    int number=rowArray->getNumElements();
    int * which=rowArray->getIndices();

    int iIndex;
    for (iIndex=0;iIndex<number;iIndex++) {
      
      int iRow = which[iIndex];
      double alpha = work[iIndex];
      int iPivot=pivotVariable_[iRow];
      dualIn_ -= alpha*cost_[iPivot];
    }
    // determine direction here
    if (dualIn_<-dualTolerance_) {
      directionIn_=1;
    } else if (dualIn_>dualTolerance_) {
      directionIn_=-1;
    } else {
      // towards nearest bound
      if (valueIn_-lowerIn_<upperIn_-valueIn_) {
	directionIn_=-1;
	dualIn_=dualTolerance_;
      } else {
	directionIn_=1;
	dualIn_=-dualTolerance_;
      }
    }
  }

  // sequence stays as row number until end
  pivotRow_=-1;
  int numberRemaining=0;

  double totalThru=0.0; // for when variables flip
  double acceptablePivot=1.0e-7;
  if (factorization_->pivots())
    acceptablePivot=1.0e-5; // if we have iterated be more strict
  double bestEverPivot=acceptablePivot;
  int lastPivotRow = -1;
  double lastPivot=0.0;
  double lastTheta=1.0e50;

  // use spareArrays to put ones looked at in
  // First one is list of candidates
  // We could compress if we really know we won't need any more
  // Second array has current set of pivot candidates
  // with a backup list saved in double * part of indexed vector

  // pivot elements
  double * spare;
  // indices
  int * index;
  spareArray->clear();
  spare = spareArray->denseVector();
  index = spareArray->getIndices();

  // we also need somewhere for effective rhs
  double * rhs=rhsArray->denseVector();
  // and we can use indices to point to alpha
  // that way we can store fabs(alpha)
  int * indexPoint = rhsArray->getIndices();
  //int numberFlip=0; // Those which may change if flips

  /*
    First we get a list of possible pivots.  We can also see if the
    problem looks unbounded.

    At first we increase theta and see what happens.  We start
    theta at a reasonable guess.  If in right area then we do bit by bit.
    We save possible pivot candidates

   */

  // do first pass to get possibles 
  // We can also see if unbounded

  double * work=rowArray->denseVector();
  int number=rowArray->getNumElements();
  int * which=rowArray->getIndices();

  // we need to swap sign if coming in from ub
  double way = directionIn_;
  double maximumMovement;
  if (way>0.0) 
    maximumMovement = CoinMin(1.0e30,upperIn_-valueIn_);
  else
    maximumMovement = CoinMin(1.0e30,valueIn_-lowerIn_);

  double averageTheta = nonLinearCost_->averageTheta();
  double tentativeTheta = CoinMin(10.0*averageTheta,maximumMovement);
  double upperTheta = maximumMovement;
  if (tentativeTheta>0.5*maximumMovement)
    tentativeTheta=maximumMovement;

  double dualCheck = fabs(dualIn_);
  // but make a bit more pessimistic
  dualCheck=max(dualCheck-100.0*dualTolerance_,0.99*dualCheck);

  int iIndex;
  bool gotList=false;
  int pivotOne=-1;
  //#define CLP_DEBUG
#ifdef CLP_DEBUG
  if (numberIterations_==-3839||numberIterations_==-3840) {
    double dj=cost_[sequenceIn_];
    printf("cost in on %d is %g, dual in %g\n",sequenceIn_,dj,dualIn_);
    for (iIndex=0;iIndex<number;iIndex++) {

      int iRow = which[iIndex];
      double alpha = work[iIndex];
      int iPivot=pivotVariable_[iRow];
      dj -= alpha*cost_[iPivot];
      printf("row %d var %d current %g %g %g, alpha %g so sol => %g (cost %g, dj %g)\n",
	     iRow,iPivot,lower_[iPivot],solution_[iPivot],upper_[iPivot],
	     alpha, solution_[iPivot]-1.0e9*alpha,cost_[iPivot],dj);
    }
  }
#endif
  while (!gotList) {
    pivotOne=-1;
    totalThru=0.0;
    // We also re-compute reduced cost
    numberRemaining=0;
    dualIn_ = cost_[sequenceIn_];
    double tolerance = primalTolerance_*1.002;
    for (iIndex=0;iIndex<number;iIndex++) {

      int iRow = which[iIndex];
      double alpha = work[iIndex];
      int iPivot=pivotVariable_[iRow];
      dualIn_ -= alpha*cost_[iPivot];
      alpha *= way;
      double oldValue = solution_[iPivot];
      // get where in bound sequence
      // note that after this alpha is actually fabs(alpha)
      if (alpha>0.0) {
	// basic variable going towards lower bound
	double bound = lower_[iPivot];
	oldValue -= bound;
      } else {
	// basic variable going towards upper bound
	double bound = upper_[iPivot];
	oldValue = bound-oldValue;
	alpha = - alpha;
      }
      
      double value = oldValue-tentativeTheta*alpha;
      assert (oldValue>=-tolerance);
      if (value<=tolerance) {
	value=oldValue-upperTheta*alpha;
	if (value<-primalTolerance_&&alpha>=acceptablePivot) {
	  upperTheta = (oldValue+primalTolerance_)/alpha;
	  pivotOne=numberRemaining;
	}
	// add to list
	spare[numberRemaining]=alpha;
	rhs[numberRemaining]=oldValue;
	indexPoint[numberRemaining]=iIndex;
	index[numberRemaining++]=iRow;
	totalThru += alpha;
	setActive(iRow);
	//} else if (value<primalTolerance_*1.002) {
	// May change if is a flip
	//indexRhs[numberFlip++]=iRow;
      }
    }
    if (upperTheta<maximumMovement&&totalThru*infeasibilityCost_>=1.0001*dualCheck) {
      // Can pivot here
      gotList=true;
    } else if (tentativeTheta<maximumMovement) {
      //printf("Going round with average theta of %g\n",averageTheta);
      tentativeTheta=maximumMovement;
    } else {
      gotList=true;
    }
  }
  totalThru=0.0;

  theta_=maximumMovement;

  bool goBackOne = false;
  if (objective_->type()>1) 
    dualIn_=saveDj;

  //printf("%d remain out of %d\n",numberRemaining,number);
  int iTry=0;
#define MAXTRY 1000
  if (numberRemaining&&upperTheta<maximumMovement) {
    // First check if previously chosen one will work
    if (pivotOne>=0&&0) {
      double thruCost = infeasibilityCost_*spare[pivotOne];
      if (thruCost>=0.99*fabs(dualIn_))
	printf("Could pivot on %d as change %g dj %g\n",
	       index[pivotOne],thruCost,dualIn_);
      double alpha = spare[pivotOne];
      double oldValue = rhs[pivotOne];
      theta_ = oldValue/alpha;
      pivotRow_=pivotOne;
      // Stop loop
      iTry=MAXTRY;
    }

    // first get ratio with tolerance
    for ( ;iTry<MAXTRY;iTry++) {
      
      upperTheta=maximumMovement;
      int iBest=-1;
      for (iIndex=0;iIndex<numberRemaining;iIndex++) {
	
	double alpha = spare[iIndex];
	double oldValue = rhs[iIndex];
	double value = oldValue-upperTheta*alpha;
	
	if (value<-primalTolerance_ && alpha>=acceptablePivot) {
	  upperTheta = (oldValue+primalTolerance_)/alpha;
	  iBest=iIndex; // just in case weird numbers
	}
      }
      
      // now look at best in this lot
      double bestPivot=acceptablePivot;
      pivotRow_=-1;
      for (iIndex=0;iIndex<numberRemaining;iIndex++) {
	
	int iRow = index[iIndex];
	double alpha = spare[iIndex];
	double oldValue = rhs[iIndex];
	double value = oldValue-upperTheta*alpha;
	
	if (value<=0||iBest==iIndex) {
	  // how much would it cost to go thru and modify bound
	  double trueAlpha=way*work[indexPoint[iIndex]];
	  totalThru += nonLinearCost_->changeInCost(pivotVariable_[iRow],trueAlpha,rhs[iIndex]);
	  setActive(iRow);
	  if (alpha>bestPivot) {
	    bestPivot=alpha;
	    theta_ = oldValue/bestPivot;
	    pivotRow_=iIndex;
	  }
	}
      }
      if (bestPivot<0.1*bestEverPivot&&
	  bestEverPivot>1.0e-6&& bestPivot<1.0e-3) {
	// back to previous one
	goBackOne = true;
	break;
      } else if (pivotRow_==-1&&upperTheta>largeValue_) {
	if (lastPivot>acceptablePivot) {
	  // back to previous one
	  goBackOne = true;
	} else {
	  // can only get here if all pivots so far too small
	}
	break;
      } else if (totalThru>=dualCheck) {
	break; // no point trying another loop
      } else {
	lastPivotRow=pivotRow_;
	lastTheta = theta_;
	if (bestPivot>bestEverPivot)
	  bestEverPivot=bestPivot;
      }    }
    // can get here without pivotRow_ set but with lastPivotRow
    if (goBackOne||(pivotRow_<0&&lastPivotRow>=0)) {
      // back to previous one
      pivotRow_=lastPivotRow;
      theta_ = lastTheta;
    }
  } else {
    // mark ones which may move
    //for (int i=0;i<numberFlip;i++) {
    //int iRow= indexRhs[i];
    //setActive(iRow);
    //}
  }
  //if (iTry>50)
  //printf("** %d tries\n",iTry);
  if (pivotRow_>=0) {
    int position=pivotRow_; // position in list
    pivotRow_=index[position];
    alpha_=work[indexPoint[position]];
    // translate to sequence
    sequenceOut_ = pivotVariable_[pivotRow_];
    valueOut_ = solution(sequenceOut_);
    lowerOut_=lower_[sequenceOut_];
    upperOut_=upper_[sequenceOut_];
#define MINIMUMTHETA 1.0e-12
    // Movement should be minimum for anti-degeneracy - unless
    // fixed variable out
    double minimumTheta;
    if (upperOut_>lowerOut_)
      minimumTheta=MINIMUMTHETA;
    else
      minimumTheta=0.0;
    // will we need to increase tolerance
    //#define CLP_DEBUG
    double largestInfeasibility = primalTolerance_;
    if (theta_<minimumTheta&&(specialOptions_&4)==0&&!valuesPass) {
      theta_=minimumTheta;
      for (iIndex=0;iIndex<numberRemaining-numberRemaining;iIndex++) {
	largestInfeasibility = CoinMax(largestInfeasibility,
				    -(rhs[iIndex]-spare[iIndex]*theta_));
      }
//#define CLP_DEBUG
#ifdef CLP_DEBUG
      if (largestInfeasibility>primalTolerance_&&(handler_->logLevel()&32)>-1)
	printf("Primal tolerance increased from %g to %g\n",
	       primalTolerance_,largestInfeasibility);
#endif
//#undef CLP_DEBUG
      primalTolerance_ = max(primalTolerance_,largestInfeasibility);
    }
    // Need to look at all in some cases
    if (theta_>tentativeTheta) {
      for (iIndex=0;iIndex<number;iIndex++) 
	setActive(which[iIndex]);
    }
    if (way<0.0) 
      theta_ = - theta_;
    double newValue = valueOut_ - theta_*alpha_;
    // If  4 bit set - Force outgoing variables to exact bound (primal)
    if (alpha_*way<0.0) {
      directionOut_=-1;      // to upper bound
      if (fabs(theta_)>1.0e-6||(specialOptions_&4)!=0) {
	upperOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
      } else {
	  upperOut_ = newValue;
      }
    } else {
      directionOut_=1;      // to lower bound
      if (fabs(theta_)>1.0e-6||(specialOptions_&4)!=0) {
	lowerOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
      } else {
	lowerOut_ = newValue;
      }
    }
    dualOut_ = reducedCost(sequenceOut_);
  } else if (maximumMovement<1.0e20) {
    // flip
    pivotRow_ = -2; // so we can tell its a flip
    sequenceOut_ = sequenceIn_;
    valueOut_ = valueIn_;
    dualOut_ = dualIn_;
    lowerOut_ = lowerIn_;
    upperOut_ = upperIn_;
    alpha_ = 0.0;
    if (way<0.0) {
      directionOut_=1;      // to lower bound
      theta_ = lowerOut_ - valueOut_;
    } else {
      directionOut_=-1;      // to upper bound
      theta_ = upperOut_ - valueOut_;
    }
  }

  double theta1 = max(theta_,1.0e-12);
  double theta2 = numberIterations_*nonLinearCost_->averageTheta();
  // Set average theta
  nonLinearCost_->setAverageTheta((theta1+theta2)/((double) (numberIterations_+1)));
  //if (numberIterations_%1000==0)
  //printf("average theta is %g\n",nonLinearCost_->averageTheta());

  // clear arrays

  memset(spare,0,numberRemaining*sizeof(double));

  // put back original bounds etc
  memcpy(rhsArray->getIndices(),index,numberRemaining*sizeof(int));
  rhsArray->setNumElements(numberRemaining);
  rhsArray->setPacked();
  nonLinearCost_->goBackAll(rhsArray);
  rhsArray->clear();

}
/* 
   Chooses primal pivot column
   updateArray has cost updates (also use pivotRow_ from last iteration)
   Would be faster with separate region to scan
   and will have this (with square of infeasibility) when steepest
   For easy problems we can just choose one of the first columns we look at
*/
void 
ClpSimplexPrimal::primalColumn(CoinIndexedVector * updates,
			       CoinIndexedVector * spareRow1,
			       CoinIndexedVector * spareRow2,
			       CoinIndexedVector * spareColumn1,
			       CoinIndexedVector * spareColumn2)
{
  sequenceIn_ = primalColumnPivot_->pivotColumn(updates,spareRow1,
					       spareRow2,spareColumn1,
					       spareColumn2);
  if (sequenceIn_>=0) {
    valueIn_=solution_[sequenceIn_];
    dualIn_=dj_[sequenceIn_];
    if (nonLinearCost_->lookBothWays()) {
      // double check 
      ClpSimplex::Status status = getStatus(sequenceIn_);
      
      switch(status) {
      case ClpSimplex::atUpperBound:
	if (dualIn_<0.0) {
	  // move to other side
	  printf("For %d U (%g, %g, %g) dj changed from %g",
		 sequenceIn_,lower_[sequenceIn_],solution_[sequenceIn_],
		 upper_[sequenceIn_],dualIn_);
	  dualIn_ -= nonLinearCost_->changeUpInCost(sequenceIn_);
	  printf(" to %g\n",dualIn_);
	  nonLinearCost_->setOne(sequenceIn_,upper_[sequenceIn_]+2.0*currentPrimalTolerance());
	  setStatus(sequenceIn_,ClpSimplex::atLowerBound);
	}
	break;
      case ClpSimplex::atLowerBound:
	if (dualIn_>0.0) {
	  // move to other side
	  printf("For %d L (%g, %g, %g) dj changed from %g",
		 sequenceIn_,lower_[sequenceIn_],solution_[sequenceIn_],
		 upper_[sequenceIn_],dualIn_);
	  dualIn_ -= nonLinearCost_->changeDownInCost(sequenceIn_);
	  printf(" to %g\n",dualIn_);
	  nonLinearCost_->setOne(sequenceIn_,lower_[sequenceIn_]-2.0*currentPrimalTolerance());
	  setStatus(sequenceIn_,ClpSimplex::atUpperBound);
	}
	break;
      default:
	break;
      }
    }
    lowerIn_=lower_[sequenceIn_];
    upperIn_=upper_[sequenceIn_];
    if (dualIn_>0.0)
      directionIn_ = -1;
    else 
      directionIn_ = 1;
  } else {
    sequenceIn_ = -1;
  }
}
/* The primals are updated by the given array.
   Returns number of infeasibilities.
   After rowArray will have list of cost changes
*/
int 
ClpSimplexPrimal::updatePrimalsInPrimal(CoinIndexedVector * rowArray,
					double theta,
					double & objectiveChange,
					int valuesPass)
{
  // Cost on pivot row may change - may need to change dualIn
  double oldCost=0.0;
  if (pivotRow_>=0)
    oldCost = cost_[sequenceOut_];
  //rowArray->scanAndPack();
  double * work=rowArray->denseVector();
  int number=rowArray->getNumElements();
  int * which=rowArray->getIndices();

  int newNumber = 0;
  int pivotPosition = -1;
  nonLinearCost_->setChangeInCost(0.0);
  int iIndex;
  if (!valuesPass) {
    for (iIndex=0;iIndex<number;iIndex++) {
      
      int iRow = which[iIndex];
      double alpha = work[iIndex];
      work[iIndex]=0.0;
      int iPivot=pivotVariable_[iRow];
      double change = theta*alpha;
      double value = solution_[iPivot] - change;
      solution_[iPivot]=value;
#ifndef NDEBUG
      // check if not active then okay
      if (!active(iRow)&&(specialOptions_&4)==0) {
	// But make sure one going out is feasible
	if (change>0.0) {
	  // going down
	  if (value<lower_[iPivot]+primalTolerance_) {
	    if (iPivot==sequenceOut_&&value>lower_[iPivot]-1.001*primalTolerance_)
	      value=lower_[iPivot];
	    double difference = nonLinearCost_->setOne(iPivot,value);
	    assert (!difference);
	  }
	} else {
	  // going up
	  if (value>upper_[iPivot]-primalTolerance_) {
	    if (iPivot==sequenceOut_&&value<upper_[iPivot]+1.001*primalTolerance_)
	      value=upper_[iPivot];
	    double difference = nonLinearCost_->setOne(iPivot,value);
	    assert (!difference);
	  }
	}
      }
#endif    
      if (active(iRow)) {
	clearActive(iRow);
	// But make sure one going out is feasible
	if (change>0.0) {
	  // going down
	  if (value<lower_[iPivot]+primalTolerance_) {
	    if (iPivot==sequenceOut_&&value>lower_[iPivot]-1.001*primalTolerance_)
	      value=lower_[iPivot];
	    double difference = nonLinearCost_->setOne(iPivot,value);
	    if (difference) {
	      if (iRow==pivotRow_)
		pivotPosition=newNumber;
	      work[newNumber] = difference;
	      //change reduced cost on this
	      dj_[iPivot] = -difference;
	      which[newNumber++]=iRow;
	    }
	  }
	} else {
	  // going up
	  if (value>upper_[iPivot]-primalTolerance_) {
	    if (iPivot==sequenceOut_&&value<upper_[iPivot]+1.001*primalTolerance_)
	      value=upper_[iPivot];
	    double difference = nonLinearCost_->setOne(iPivot,value);
	    if (difference) {
	      if (iRow==pivotRow_)
		pivotPosition=newNumber;
	      work[newNumber] = difference;
	      //change reduced cost on this
	      dj_[iPivot] = -difference;
	      which[newNumber++]=iRow;
	    }
	  }
	}
      }
    }
  } else {
    // values pass so look at all
    for (iIndex=0;iIndex<number;iIndex++) {
      
      int iRow = which[iIndex];
      double alpha = work[iIndex];
      work[iIndex]=0.0;
      int iPivot=pivotVariable_[iRow];
      double change = theta*alpha;
      double value = solution_[iPivot] - change;
      solution_[iPivot]=value;
      clearActive(iRow);
      // But make sure one going out is feasible
      if (change>0.0) {
	// going down
	if (value<lower_[iPivot]+primalTolerance_) {
	  if (iPivot==sequenceOut_&&value>lower_[iPivot]-1.001*primalTolerance_)
	    value=lower_[iPivot];
	  double difference = nonLinearCost_->setOne(iPivot,value);
	  if (difference) {
	    if (iRow==pivotRow_)
	      pivotPosition=newNumber;
	    work[newNumber] = difference;
	    //change reduced cost on this
	    dj_[iPivot] = -difference;
	    which[newNumber++]=iRow;
	  }
	}
      } else {
	// going up
	if (value>upper_[iPivot]-primalTolerance_) {
	  if (iPivot==sequenceOut_&&value<upper_[iPivot]+1.001*primalTolerance_)
	    value=upper_[iPivot];
	  double difference = nonLinearCost_->setOne(iPivot,value);
	  if (difference) {
	    if (iRow==pivotRow_)
	      pivotPosition=newNumber;
	    work[newNumber] = difference;
	    //change reduced cost on this
	    dj_[iPivot] = -difference;
	    which[newNumber++]=iRow;
	  }
	}
      }
    }
  }
  objectiveChange += nonLinearCost_->changeInCost();
  rowArray->setPacked();
#if 0
  rowArray->setNumElements(newNumber);
  rowArray->expand();
  if (pivotRow_>=0) {
    dualIn_ += (oldCost-cost_[sequenceOut_]);
    // update change vector to include pivot
    rowArray->add(pivotRow_,-dualIn_);
    // and convert to packed
    rowArray->scanAndPack();
  } else {
    // and convert to packed
    rowArray->scanAndPack();
  }
#else
  if (pivotRow_>=0) {
    double dualIn = dualIn_+(oldCost-cost_[sequenceOut_]);
    // update change vector to include pivot
    if (pivotPosition>=0) {
      work[pivotPosition] -= dualIn;
    } else {
      work[newNumber]=-dualIn;
      which[newNumber++]=pivotRow_;
    }
  }
  rowArray->setNumElements(newNumber);
#endif
  return 0;
}
// Perturbs problem
void 
ClpSimplexPrimal::perturb(int type)
{
  if (perturbation_>100)
    return; //perturbed already
  if (perturbation_==100)
    perturbation_=50; // treat as normal
  int savePerturbation = perturbation_;
  int i;
  if (!numberIterations_)
    cleanStatus(); // make sure status okay
  // look at element range
  double smallestNegative;
  double largestNegative;
  double smallestPositive;
  double largestPositive;
  matrix_->rangeOfElements(smallestNegative, largestNegative,
			   smallestPositive, largestPositive);
  smallestPositive = CoinMin(fabs(smallestNegative),smallestPositive);
  largestPositive = max(fabs(largestNegative),largestPositive);
  double elementRatio = largestPositive/smallestPositive;
  if (!numberIterations_&&perturbation_==50) {
    // See if we need to perturb
    double * sort = new double[numberRows_];
    for (i=0;i<numberRows_;i++) {
      double lo = fabs(lower_[i]);
      double up = fabs(upper_[i]);
      double value=0.0;
      if (lo&&lo<1.0e20) {
	if (up&&up<1.0e20)
	  value = 0.5*(lo+up);
	else
	  value=lo;
      } else {
	if (up&&up<1.0e20)
	  value = up;
      }
      sort[i]=value;
    }
    std::sort(sort,sort+numberRows_);
    int number=1;
    double last = sort[0];
    for (i=1;i<numberRows_;i++) {
      if (last!=sort[i])
	number++;
      last=sort[i];
    }
    delete [] sort;
    //printf("ratio number diff rhs %g, element ratio %g\n",((double)number)/((double) numberRows_),
    //								      elementRatio);
    if (number*4>numberRows_||elementRatio>1.0e12) {
      perturbation_=100;
      return; // good enough
    }
  }
  // primal perturbation
  double perturbation=1.0e-20;
  int numberNonZero=0;
  // maximum fraction of rhs/bounds to perturb
  double maximumFraction = 1.0e-5;
  if (perturbation_>=50) {
    perturbation = 1.0e-4;
    for (i=0;i<numberColumns_+numberRows_;i++) {
      if (upper_[i]>lower_[i]+primalTolerance_) {
	double lowerValue, upperValue;
	if (lower_[i]>-1.0e20)
	  lowerValue = fabs(lower_[i]);
	else
	  lowerValue=0.0;
	if (upper_[i]<1.0e20)
	  upperValue = fabs(upper_[i]);
	else
	  upperValue=0.0;
	double value = max(fabs(lowerValue),fabs(upperValue));
	value = CoinMin(value,upper_[i]-lower_[i]);
#if 1
	if (value) {
	  perturbation += value;
	  numberNonZero++;
	}
#else
	perturbation = max(perturbation,value);
#endif
      }
    }
    if (numberNonZero) 
      perturbation /= (double) numberNonZero;
    else
      perturbation = 1.0e-1;
  } else if (perturbation_<100) {
    perturbation = pow(10.0,perturbation_);
    // user is in charge
    maximumFraction = 1.0;
  }
  double largestZero=0.0;
  double largest=0.0;
  double largestPerCent=0.0;
  bool printOut=(handler_->logLevel()==63);
  printOut=false; //off
  // Check if all slack
  int number=0;
  int iSequence;
  for (iSequence=0;iSequence<numberRows_;iSequence++) {
    if (getRowStatus(iSequence)==basic) 
      number++;
  }
  if (rhsScale_>100.0) {
    // tone down perturbation
    maximumFraction *= 0.1;
  }
  if (number!=numberRows_)
    type=1;
  // modify bounds
  // Change so at least 1.0e-5 and no more than 0.1
  // For now just no more than 0.1
  // printf("Pert type %d perturbation %g, maxF %g\n",type,perturbation,maximumFraction);
  if (type==1) {
    //double multiplier = perturbation*maximumFraction;
    for (iSequence=0;iSequence<numberRows_+numberColumns_;iSequence++) {
      if (getStatus(iSequence)==basic) {
	double solutionValue = solution_[iSequence];
	double lowerValue = lower_[iSequence];
	double upperValue = upper_[iSequence];
	double difference = upperValue-lowerValue;
	difference = CoinMin(difference,perturbation);
	difference = CoinMin(difference,fabs(solutionValue)+1.0);
	double value = maximumFraction*(difference+1.0);
	value = CoinMin(value,0.1);
	value *= CoinDrand48();
	if (solutionValue-lowerValue<=primalTolerance_) {
	  lower_[iSequence] -= value;
	} else if (upperValue-solutionValue<=primalTolerance_) {
	  upper_[iSequence] += value;
	} else {
#if 0
	  if (iSequence>=numberColumns_) {
	    // may not be at bound - but still perturb (unless free)
	    if (upperValue>1.0e30&&lowerValue<-1.0e30)
	      value=0.0;
	    else
	      value = - value; // as -1.0 in matrix
	  } else {
	    value = 0.0;
	  }
#else
	  value=0.0;
#endif
	}
	if (value) {
	  if (printOut)
	    printf("col %d lower from %g to %g, upper from %g to %g\n",
		   iSequence,lower_[iSequence],lowerValue,upper_[iSequence],upperValue);
	  if (solutionValue) {
	    largest = max(largest,value);
	    if (value>(fabs(solutionValue)+1.0)*largestPerCent)
	      largestPerCent=value/(fabs(solutionValue)+1.0);
	  } else {
	    largestZero = max(largestZero,value);
	  } 
	}
      }
    }
  } else {
    double tolerance = 100.0*primalTolerance_;
    for (i=0;i<numberColumns_;i++) {
      double lowerValue=lower_[i], upperValue=upper_[i];
      if (upperValue>lowerValue+primalTolerance_) {
	double value = perturbation*maximumFraction;
	value = CoinMin(value,0.1);
	value *= CoinDrand48();
	if (savePerturbation!=50) {
	  if (fabs(value)<=primalTolerance_)
	    value=0.0;
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += value * (max(1.0e-2,1.0e-5*fabs(upperValue))); 
	} else if (value) {
	  double valueL =value *(max(1.0e-2,1.0e-5*fabs(lowerValue)));
	  // get in range 
	  if (valueL<=tolerance) {
	    valueL *= 10.0;
	    while (valueL<=tolerance) 
	      valueL *= 10.0;
	  } else if (valueL>1.0) {
	    valueL *= 0.1;
	    while (valueL>1.0) 
	      valueL *= 0.1;
	  }
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= valueL;
	  double valueU =value *(max(1.0e-2,1.0e-5*fabs(upperValue)));
	  // get in range 
	  if (valueU<=tolerance) {
	    valueU *= 10.0;
	    while (valueU<=tolerance) 
	      valueU *= 10.0;
	  } else if (valueU>1.0) {
	    valueU *= 0.1;
	    while (valueU>1.0) 
	      valueU *= 0.1;
	  }
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += valueU;
	}
	if (lowerValue!=lower_[i]) {
	  double difference = fabs(lowerValue-lower_[i]);
	  largest = max(largest,difference);
	  if (difference>fabs(lower_[i])*largestPerCent)
	    largestPerCent=fabs(difference/lower_[i]);
	} 
	if (upperValue!=upper_[i]) {
	  double difference = fabs(upperValue-upper_[i]);
	  largest = max(largest,difference);
	  if (difference>fabs(upper_[i])*largestPerCent)
	    largestPerCent=fabs(difference/upper_[i]);
	} 
	if (printOut)
	  printf("col %d lower from %g to %g, upper from %g to %g\n",
		 i,lower_[i],lowerValue,upper_[i],upperValue);
      }
      lower_[i]=lowerValue;
      upper_[i]=upperValue;
    }
    for (;i<numberColumns_+numberRows_;i++) {
      double lowerValue=lower_[i], upperValue=upper_[i];
      double value = perturbation*maximumFraction;
      value = CoinMin(value,0.1);
      value *= CoinDrand48();
      if (upperValue>lowerValue+tolerance) {
	if (savePerturbation!=50) {
	  if (fabs(value)<=primalTolerance_)
	    value=0.0;
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += value * (max(1.0e-2,1.0e-5*fabs(upperValue))); 
	} else if (value) {
	  double valueL =value *(max(1.0e-2,1.0e-5*fabs(lowerValue)));
	  // get in range 
	  if (valueL<=tolerance) {
	    valueL *= 10.0;
	    while (valueL<=tolerance) 
	      valueL *= 10.0;
	  } else if (valueL>1.0) {
	    valueL *= 0.1;
	    while (valueL>1.0) 
	      valueL *= 0.1;
	  }
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= valueL;
	  double valueU =value *(max(1.0e-2,1.0e-5*fabs(upperValue)));
	  // get in range 
	  if (valueU<=tolerance) {
	    valueU *= 10.0;
	    while (valueU<=tolerance) 
	      valueU *= 10.0;
	  } else if (valueU>1.0) {
	    valueU *= 0.1;
	    while (valueU>1.0) 
	      valueU *= 0.1;
	  }
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += valueU;
	}
      } else if (upperValue>0.0) {
	upperValue -= value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
	lowerValue -= value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
      } else if (upperValue<0.0) {
	upperValue += value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
	lowerValue += value * (max(1.0e-2,1.0e-5*fabs(lowerValue))); 
      } else {
      }
      if (lowerValue!=lower_[i]) {
	double difference = fabs(lowerValue-lower_[i]);
	largest = max(largest,difference);
	if (difference>fabs(lower_[i])*largestPerCent)
	  largestPerCent=fabs(difference/lower_[i]);
      } 
      if (upperValue!=upper_[i]) {
	double difference = fabs(upperValue-upper_[i]);
	largest = max(largest,difference);
	if (difference>fabs(upper_[i])*largestPerCent)
	  largestPerCent=fabs(difference/upper_[i]);
      } 
      if (printOut)
	printf("row %d lower from %g to %g, upper from %g to %g\n",
	       i-numberColumns_,lower_[i],lowerValue,upper_[i],upperValue);
      lower_[i]=lowerValue;
      upper_[i]=upperValue;
    }
  }
  // Clean up
  for (i=0;i<numberColumns_+numberRows_;i++) {
    switch(getStatus(i)) {
      
    case basic:
      break;
    case atUpperBound:
      solution_[i]=upper_[i];
      break;
    case isFixed:
    case atLowerBound:
      solution_[i]=lower_[i];
      break;
    case isFree:
      break;
    case superBasic:
      break;
    }
  }
  handler_->message(CLP_SIMPLEX_PERTURB,messages_)
    <<100.0*maximumFraction<<perturbation<<largest<<100.0*largestPerCent<<largestZero
    <<CoinMessageEol;
  // redo nonlinear costs
  // say perturbed
  perturbation_=101;
}
// un perturb
bool
ClpSimplexPrimal::unPerturb()
{
  if (perturbation_!=101)
    return false;
  // put back original bounds and costs
  createRim(7);
  sanityCheck();
  // unflag
  unflag();
  // get a valid nonlinear cost function
  delete nonLinearCost_;
  nonLinearCost_= new ClpNonLinearCost(this);
  perturbation_ = 102; // stop any further perturbation
  // move non basic variables to new bounds
  nonLinearCost_->checkInfeasibilities(0.0);
#if 1
  // Try using dual
  return true;
#else
  gutsOfSolution(NULL,NULL,ifValuesPass!=0);
  return false;
#endif
  
}
// Unflag all variables and return number unflagged
int 
ClpSimplexPrimal::unflag()
{
  int i;
  int number = numberRows_+numberColumns_;
  int numberFlagged=0;
  for (i=0;i<number;i++) {
    if (flagged(i)) {
      clearFlagged(i);
      numberFlagged++;
    }
  }
  numberFlagged += matrix_->generalExpanded(this,8,i);
  if (handler_->logLevel()>2&&numberFlagged&&objective_->type()>1)
    printf("%d unflagged\n",numberFlagged);
  return numberFlagged;
}
// Do not change infeasibility cost and always say optimal
void 
ClpSimplexPrimal::alwaysOptimal(bool onOff)
{
  if (onOff)
    specialOptions_ |= 1;
  else
    specialOptions_ &= ~1;
}
bool 
ClpSimplexPrimal::alwaysOptimal() const
{
  return (specialOptions_&1)!=0;
}
// Flatten outgoing variables i.e. - always to exact bound
void 
ClpSimplexPrimal::exactOutgoing(bool onOff)
{
  if (onOff)
    specialOptions_ |= 4;
  else
    specialOptions_ &= ~4;
}
bool 
ClpSimplexPrimal::exactOutgoing() const
{
  return (specialOptions_&4)!=0;
}
/*
  Reasons to come out (normal mode/user mode):
  -1 normal
  -2 factorize now - good iteration/ NA
  -3 slight inaccuracy - refactorize - iteration done/ same but factor done
  -4 inaccuracy - refactorize - no iteration/ NA
  -5 something flagged - go round again/ pivot not possible
  +2 looks unbounded
  +3 max iterations (iteration done)
*/
int
ClpSimplexPrimal::pivotResult(int ifValuesPass)
{

  bool roundAgain=true;
  int returnCode=-1;

  // loop round if user setting and doing refactorization
  while (roundAgain) {
    roundAgain=false;
    returnCode=-1;
    pivotRow_=-1;
    sequenceOut_=-1;
    rowArray_[1]->clear();
#if 0
    {
      int seq[]={612,643};
      int k;
      for (k=0;k<sizeof(seq)/sizeof(int);k++) {
	int iSeq=seq[k];
	if (getColumnStatus(iSeq)!=basic) {
	  double djval;
	  double * work;
	  int number;
	  int * which;
	  
	  int iIndex;
	  unpack(rowArray_[1],iSeq);
	  factorization_->updateColumn(rowArray_[2],rowArray_[1]);
	  djval = cost_[iSeq];
	  work=rowArray_[1]->denseVector();
	  number=rowArray_[1]->getNumElements();
	  which=rowArray_[1]->getIndices();
	  
	  for (iIndex=0;iIndex<number;iIndex++) {
	    
	    int iRow = which[iIndex];
	    double alpha = work[iRow];
	    int iPivot=pivotVariable_[iRow];
	    djval -= alpha*cost_[iPivot];
	  }
	  double comp = 1.0e-8 + 1.0e-7*(max(fabs(dj_[iSeq]),fabs(djval)));
	  if (fabs(djval-dj_[iSeq])>comp)
	    printf("Bad dj %g for %d - true is %g\n",
		   dj_[iSeq],iSeq,djval);
	  assert (fabs(djval)<1.0e-3||djval*dj_[iSeq]>0.0);
	  rowArray_[1]->clear();
	}
      }
    }
#endif
	
    // we found a pivot column
    // update the incoming column
    unpackPacked(rowArray_[1]);
    // save reduced cost
    double saveDj = dualIn_;
    factorization_->updateColumnFT(rowArray_[2],rowArray_[1]);
    // Get extra rows
    matrix_->extendUpdated(this,rowArray_[1],0);
    // do ratio test and re-compute dj
    primalRow(rowArray_[1],rowArray_[3],rowArray_[2],rowArray_[0],
	      ifValuesPass);
    if (ifValuesPass) {
      saveDj=dualIn_;
      //assert (fabs(alpha_)>=1.0e-5||(objective_->type()<2||!objective_->activated())||pivotRow_==-2);
      if (pivotRow_==-1||(pivotRow_>=0&&fabs(alpha_)<1.0e-5)) {
	if(fabs(dualIn_)<1.0e2*dualTolerance_&&objective_->type()<2) {
	  // try other way
	  directionIn_=-directionIn_;
	  primalRow(rowArray_[1],rowArray_[3],rowArray_[2],rowArray_[0],
		    0);
	}
	if (pivotRow_==-1||(pivotRow_>=0&&fabs(alpha_)<1.0e-5)) {
	  if (solveType_==1) {
	    // reject it
	    char x = isColumn(sequenceIn_) ? 'C' :'R';
	    handler_->message(CLP_SIMPLEX_FLAG,messages_)
	      <<x<<sequenceWithin(sequenceIn_)
	      <<CoinMessageEol;
	    setFlagged(sequenceIn_);
	    progress_->clearBadTimes();
	    lastBadIteration_ = numberIterations_; // say be more cautious
	    clearAll();
	    pivotRow_=-1;
	  }
	  returnCode=-5;
	  break;
	}
      }
    }
    // need to clear toIndex_ in gub
    // ? when can I clear stuff
    // Clean up any gub stuff
    matrix_->extendUpdated(this,rowArray_[1],1);
    double checkValue=1.0e-2;
    if (largestDualError_>1.0e-5)
      checkValue=1.0e-1;
    if (solveType_==1&&((saveDj*dualIn_<1.0e-20&&!ifValuesPass)||
	fabs(saveDj-dualIn_)>checkValue*(1.0+fabs(saveDj)))) {
      char x = isColumn(sequenceIn_) ? 'C' :'R';
      handler_->message(CLP_PRIMAL_DJ,messages_)
	<<x<<sequenceIn_<<saveDj<<dualIn_
	<<CoinMessageEol;
      if(lastGoodIteration_ != numberIterations_) {
	clearAll();
	pivotRow_=-1; // say no weights update
	returnCode=-4;
	if(lastGoodIteration_+1 == numberIterations_) {
	  // not looking wonderful - try cleaning bounds
	  // put non-basics to bounds in case tolerance moved
	  nonLinearCost_->checkInfeasibilities(0.0);
	}
	sequenceOut_=-1;
	break;
      } else {
	// take on more relaxed criterion
	if (saveDj*dualIn_<1.0e-20||
	    fabs(saveDj-dualIn_)>2.0e-1*(1.0+fabs(dualIn_))) {
	  // need to reject something
	  char x = isColumn(sequenceIn_) ? 'C' :'R';
	  handler_->message(CLP_SIMPLEX_FLAG,messages_)
	    <<x<<sequenceWithin(sequenceIn_)
	    <<CoinMessageEol;
	  setFlagged(sequenceIn_);
	  progress_->clearBadTimes();
	  lastBadIteration_ = numberIterations_; // say be more cautious
	  clearAll();
	  pivotRow_=-1;
	  returnCode=-5;
	  sequenceOut_=-1;
	  break;
	}
      }
    }
    if (pivotRow_>=0) {
      if (solveType_==2) {
	// **** Coding for user interface
	// do ray
	primalRay(rowArray_[1]);
	// update duals
	// as packed need to find pivot row
	//assert (rowArray_[1]->packedMode());
	//int i;
	
	//alpha_ = rowArray_[1]->denseVector()[pivotRow_];
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
	assert (fabs(dj_[sequenceIn_])<1.0e-1);
      }
      // if stable replace in basis
      // If gub or odd then alpha and pivotRow may change
      int updateType=0;
      int updateStatus = matrix_->generalExpanded(this,3,updateType);
      if (updateType>=0)
	updateStatus = factorization_->replaceColumn(this,
						     rowArray_[2],
						     rowArray_[1],
						     pivotRow_,
						     alpha_);

      // if no pivots, bad update but reasonable alpha - take and invert
      if (updateStatus==2&&
	  lastGoodIteration_==numberIterations_&&fabs(alpha_)>1.0e-5)
	updateStatus=4;
      if (updateStatus==1||updateStatus==4) {
	// slight error
	if (factorization_->pivots()>5||updateStatus==4) {
	  returnCode=-3;
	}
      } else if (updateStatus==2) {
	// major error
	// better to have small tolerance even if slower
	factorization_->zeroTolerance(1.0e-15);
	int maxFactor = factorization_->maximumPivots();
	if (maxFactor>10) {
	  if (forceFactorization_<0)
	    forceFactorization_= maxFactor;
	  forceFactorization_ = CoinMax(1,(forceFactorization_>>1));
	} 
	// later we may need to unwind more e.g. fake bounds
	if(lastGoodIteration_ != numberIterations_) {
	  clearAll();
	  pivotRow_=-1;
	  if (solveType_==1) {
	    returnCode=-4;
	    break;
	  } else {
	    // user in charge - re-factorize
	    int lastCleaned;
	    ClpSimplexProgress dummyProgress;
	    if (saveStatus_)
	      statusOfProblemInPrimal(lastCleaned,1,&dummyProgress,true,ifValuesPass);
	    else
	      statusOfProblemInPrimal(lastCleaned,0,&dummyProgress,true,ifValuesPass);
	    roundAgain=true;
	    continue;
	  }
	} else {
	  // need to reject something
	  if (solveType_==1) {
	    char x = isColumn(sequenceIn_) ? 'C' :'R';
	    handler_->message(CLP_SIMPLEX_FLAG,messages_)
	      <<x<<sequenceWithin(sequenceIn_)
	      <<CoinMessageEol;
	    setFlagged(sequenceIn_);
	    progress_->clearBadTimes();
	  }
	  lastBadIteration_ = numberIterations_; // say be more cautious
	  clearAll();
	  pivotRow_=-1;
	  sequenceOut_=-1;
	  returnCode = -5;
	  break;

	}
      } else if (updateStatus==3) {
	// out of memory
	// increase space if not many iterations
	if (factorization_->pivots()<
	    0.5*factorization_->maximumPivots()&&
	    factorization_->pivots()<200)
	  factorization_->areaFactor(
				     factorization_->areaFactor() * 1.1);
	returnCode =-2; // factorize now
      } else if (updateStatus==5) {
	problemStatus_=-2; // factorize now
      }
      // here do part of steepest - ready for next iteration
      if (!ifValuesPass)
	primalColumnPivot_->updateWeights(rowArray_[1]);
    } else {
      if (pivotRow_==-1) {
	// no outgoing row is valid
	rowArray_[0]->clear();
	if (!factorization_->pivots()) {
	  returnCode = 2; //say looks unbounded
	  // do ray
	  primalRay(rowArray_[1]);
	} else if (solveType_==2) {
	  // refactorize
	  int lastCleaned;
	  ClpSimplexProgress dummyProgress;
	  if (saveStatus_)
	    statusOfProblemInPrimal(lastCleaned,1,&dummyProgress,true,ifValuesPass);
	  else
	    statusOfProblemInPrimal(lastCleaned,0,&dummyProgress,true,ifValuesPass);
	  roundAgain=true;
	  continue;
	} else {
	  returnCode = 4; //say looks unbounded but has iterated
	}
	break;
      } else {
	// flipping from bound to bound
      }
    }

    
    // update primal solution
    
    double objectiveChange=0.0;
    // after this rowArray_[1] is not empty - used to update djs
    // If pivot row >= numberRows then may be gub
    int savePivot = pivotRow_;
    if (pivotRow_>=numberRows_)
      pivotRow_=-1;
    updatePrimalsInPrimal(rowArray_[1],theta_, objectiveChange,ifValuesPass);
    pivotRow_=savePivot;
    
    double oldValue = valueIn_;
    if (directionIn_==-1) {
      // as if from upper bound
      if (sequenceIn_!=sequenceOut_) {
	// variable becoming basic
	valueIn_ -= fabs(theta_);
      } else {
	valueIn_=lowerIn_;
      }
    } else {
      // as if from lower bound
      if (sequenceIn_!=sequenceOut_) {
	// variable becoming basic
	valueIn_ += fabs(theta_);
      } else {
	valueIn_=upperIn_;
      }
    }
    objectiveChange += dualIn_*(valueIn_-oldValue);
    // outgoing
    if (sequenceIn_!=sequenceOut_) {
      if (directionOut_>0) {
	valueOut_ = lowerOut_;
      } else {
	valueOut_ = upperOut_;
      }
      if(valueOut_<lower_[sequenceOut_]-primalTolerance_)
	valueOut_=lower_[sequenceOut_]-0.9*primalTolerance_;
      else if (valueOut_>upper_[sequenceOut_]+primalTolerance_)
	valueOut_=upper_[sequenceOut_]+0.9*primalTolerance_;
      // may not be exactly at bound and bounds may have changed
      // Make sure outgoing looks feasible
      directionOut_=nonLinearCost_->setOneOutgoing(sequenceOut_,valueOut_);
      solution_[sequenceOut_]=valueOut_;
    }
    // change cost and bounds on incoming if primal
    nonLinearCost_->setOne(sequenceIn_,valueIn_); 
    int whatNext=housekeeping(objectiveChange);
#if CLP_DEBUG >1
    {
      int ninf= matrix_->checkFeasible(this);
      if (ninf)
	printf("infeas %d\n",ninf);
    }
#endif
    if (whatNext==1) {
	returnCode =-2; // refactorize
    } else if (whatNext==2) {
      // maximum iterations or equivalent
      returnCode=3;
    } else if(numberIterations_ == lastGoodIteration_
	      + 2 * factorization_->maximumPivots()) {
      // done a lot of flips - be safe
      returnCode =-2; // refactorize
    }
    // Check event
    {
      int status = eventHandler_->event(ClpEventHandler::endOfIteration);
      if (status>=0) {
	problemStatus_=5;
	secondaryStatus_=ClpEventHandler::endOfIteration;
	returnCode=4;
      }
    }
  }
  if (solveType_==2&&(returnCode == -2||returnCode==-3)) {
    // refactorize here
    int lastCleaned;
    ClpSimplexProgress dummyProgress;
    if (saveStatus_)
      statusOfProblemInPrimal(lastCleaned,1,&dummyProgress,true,ifValuesPass);
    else
      statusOfProblemInPrimal(lastCleaned,0,&dummyProgress,true,ifValuesPass);
    if (problemStatus_==5) {
      printf("Singular basis\n");
      problemStatus_=-1;
      returnCode=5;
    }
  }
#ifdef CLP_DEBUG
  {
    int i;
    // not [1] as may have information
    for (i=0;i<4;i++) {
      if (i!=1)
	rowArray_[i]->checkClear();
    }    
    for (i=0;i<2;i++) {
      columnArray_[i]->checkClear();
    }    
  }      
#endif
  return returnCode;
}
// Create primal ray
void 
ClpSimplexPrimal::primalRay(CoinIndexedVector * rowArray)
{
  delete [] ray_;
  ray_ = new double [numberColumns_];
  ClpFillN(ray_,numberColumns_,0.0);
  int number=rowArray->getNumElements();
  int * index = rowArray->getIndices();
  double * array = rowArray->denseVector();
  double way=-directionIn_;
  int i;
  double zeroTolerance=1.0e-12;
  if (sequenceIn_<numberColumns_)
    ray_[sequenceIn_]=directionIn_;
  if (!rowArray->packedMode()) {
    for (i=0;i<number;i++) {
      int iRow=index[i];
      int iPivot=pivotVariable_[iRow];
      double arrayValue = array[iRow];
      if (iPivot<numberColumns_&&fabs(arrayValue)>=zeroTolerance)
	ray_[iPivot] = way* arrayValue;
    }
  } else {
    for (i=0;i<number;i++) {
      int iRow=index[i];
      int iPivot=pivotVariable_[iRow];
      double arrayValue = array[i];
      if (iPivot<numberColumns_&&fabs(arrayValue)>=zeroTolerance)
	ray_[iPivot] = way* arrayValue;
    }
  }
}
/* Get next superbasic -1 if none,
   Normal type is 1
   If type is 3 then initializes sorted list
   if 2 uses list.
*/
int 
ClpSimplexPrimal::nextSuperBasic(int superBasicType,CoinIndexedVector * columnArray)
{
  if (firstFree_>=0&&superBasicType) {
    int returnValue=-1;
    bool finished=false;
    while (!finished) {
      returnValue=firstFree_;
      int iColumn=firstFree_+1;
      if (superBasicType>1) {
	if (superBasicType>2) {
	  // Initialize list
	  // Wild guess that lower bound more natural than upper
	  int number=0;
	  double * work=columnArray->denseVector();
	  int * which=columnArray->getIndices();
	  for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
	    if (!flagged(iColumn)) {
	      if (getStatus(iColumn)==superBasic) {
		if (fabs(solution_[iColumn]-lower_[iColumn])<=primalTolerance_) {
		  solution_[iColumn]=lower_[iColumn];
		  setStatus(iColumn,atLowerBound);
		} else if (fabs(solution_[iColumn]-upper_[iColumn])
			   <=primalTolerance_) {
		  solution_[iColumn]=upper_[iColumn];
		  setStatus(iColumn,atUpperBound);
		} else if (lower_[iColumn]<-1.0e20&&upper_[iColumn]>1.0e20) {
		  setStatus(iColumn,isFree);
		  break;
		} else if (!flagged(iColumn)) {
		  // put ones near bounds at end after sorting
		  work[number]= - CoinMin(0.1*(solution_[iColumn]-lower_[iColumn]),
				      upper_[iColumn]-solution_[iColumn]);
		  which[number++] = iColumn;
		}
	      }
	    }
	  }
	  CoinSort_2(work,work+number,which);
	  columnArray->setNumElements(number);
	  memset(work,0,number*sizeof(double));
	}
	int * which=columnArray->getIndices();
	int number = columnArray->getNumElements();
	if (!number) {
	  // finished
	  iColumn = numberRows_+numberColumns_;
	  returnValue=-1;
	} else {
	  number--;
	  returnValue=which[number];
	  iColumn=returnValue;
	  columnArray->setNumElements(number);
	}      
      } else {
	for (;iColumn<numberRows_+numberColumns_;iColumn++) {
	  if (!flagged(iColumn)) {
	    if (getStatus(iColumn)==superBasic) {
	      if (fabs(solution_[iColumn]-lower_[iColumn])<=primalTolerance_) {
		solution_[iColumn]=lower_[iColumn];
		setStatus(iColumn,atLowerBound);
	      } else if (fabs(solution_[iColumn]-upper_[iColumn])
			 <=primalTolerance_) {
		solution_[iColumn]=upper_[iColumn];
		setStatus(iColumn,atUpperBound);
	      } else if (lower_[iColumn]<-1.0e20&&upper_[iColumn]>1.0e20) {
		setStatus(iColumn,isFree);
		break;
	      } else {
		break;
	      }
	    }
	  }
	}
      }
      firstFree_ = iColumn;
      finished=true;
      if (firstFree_==numberRows_+numberColumns_)
	firstFree_=-1;
      if (returnValue>=0&&getStatus(returnValue)!=superBasic)
	finished=false; // somehow picked up odd one
    }
    return returnValue;
  } else {
    return -1;
  }
}
void
ClpSimplexPrimal::clearAll()
{
  // Clean up any gub stuff
  matrix_->extendUpdated(this,rowArray_[1],1);
  int number=rowArray_[1]->getNumElements();
  int * which=rowArray_[1]->getIndices();
  
  int iIndex;
  for (iIndex=0;iIndex<number;iIndex++) {
    
    int iRow = which[iIndex];
    clearActive(iRow);
  }
  rowArray_[1]->clear();
  // make sure any gub sets are clean
  matrix_->generalExpanded(this,11,sequenceIn_);
  
}

