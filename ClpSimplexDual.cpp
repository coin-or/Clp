// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.


/* Notes on implementation of dual simplex algorithm.

   When dual feasible:

   If primal feasible, we are optimal.  Otherwise choose an infeasible
   basic variable to leave basis (normally going to nearest bound) (B).  We
   now need to find an incoming variable which will leave problem
   dual feasible so we get the row of the tableau corresponding to
   the basic variable (with the correct sign depending if basic variable
   above or below feasibility region - as that affects whether reduced
   cost on outgoing variable has to be positive or negative).

   We now perform a ratio test to determine which incoming variable will
   preserve dual feasibility (C).  If no variable found then problem
   is infeasible (in primal sense).  If there is a variable, we then 
   perform pivot and repeat.  Trivial?

   -------------------------------------------

   A) How do we get dual feasible?  If all variables have bounds then
   it is trivial to get feasible by putting non-basic variables to
   correct bounds.  OSL did not have a phase 1/phase 2 approach but
   instead effectively put fake bounds on variables and this is the
   approach here, although I had hoped to make it cleaner.

   If there is a weight of X on getting dual feasible:
     Non-basic variables with negative reduced costs are put to
     lesser of their upper bound and their lower bound + X.
     Similarly, mutatis mutandis, for positive reduced costs.

   Free variables should normally be in basis, otherwise I have
   coding which may be able to come out (and may not be correct).

   In OSL, this weight was changed heuristically, here at present
   it is only increased if problem looks finished.  If problem is
   feasible I check for unboundedness.  If not unbounded we 
   could play with going into primal.  As long as weights increase
   any algorithm would be finite.
   
   B) Which outgoing variable to choose is a virtual base class.
   For difficult problems steepest edge is preferred while for
   very easy (large) problems we will need partial scan.

   C) Sounds easy, but this is hardest part of algorithm.
      1) Instead of stopping at first choice, we may be able
      to flip that variable to other bound and if objective
      still improving choose again.  These mini iterations can
      increase speed by orders of magnitude but we may need to
      go to more of a bucket choice of variable rather than looking
      at them one by one (for speed).
      2) Accuracy.  Reduced costs may be of wrong sign but less than
      tolerance.  Pivoting on these makes objective go backwards.
      OSL modified cost so a zero move was made, Gill et al
      (in primal analogue) modified so a strictly positive move was 
      made.  It is not quite as neat in dual but that is what we
      try and do.  The two problems are that re-factorizations can
      change reduced costs above and below tolerances and that when
      finished we need to reset costs and try again.
      3) Degeneracy.  Gill et al helps but may not be enough.  We
      may need more.  Also it can improve speed a lot if we perturb
      the costs significantly.  

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
#include "ClpSimplexDual.hpp"
#include "ClpFactorization.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpMessage.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream>
#define CHECK_DJ
// dual 
int ClpSimplexDual::dual ( )
{

  /* *** Method
     This is a vanilla version of dual simplex.

     It tries to be a single phase approach with a weight of 1.0 being
     given to getting optimal and a weight of dualBound_ being
     given to getting dual feasible.  In this version I have used the
     idea that this weight can be thought of as a fake bound.  If the
     distance between the lower and upper bounds on a variable is less
     than the feasibility weight then we are always better off flipping
     to other bound to make dual feasible.  If the distance is greater
     then we make up a fake bound dualBound_ away from one bound.
     If we end up optimal or primal infeasible, we check to see if
     bounds okay.  If so we have finished, if not we increase dualBound_
     and continue (after checking if unbounded). I am undecided about
     free variables - there is coding but I am not sure about it.  At
     present I put them in basis anyway.

     The code is designed to take advantage of sparsity so arrays are
     seldom zeroed out from scratch or gone over in their entirety.
     The only exception is a full scan to find outgoing variable.  This
     will be changed to keep an updated list of infeasibilities (or squares
     if steepest edge).  Also on easy problems we don't need full scan - just
     pick first reasonable.

     One problem is how to tackle degeneracy and accuracy.  At present
     I am using the modification of costs which I put in OSL and which was
     extended by Gill et al.  I am still not sure of the exact details.

     The flow of dual is three while loops as follows:

     while (not finished) {

       while (not clean solution) {

          Factorize and/or clean up solution by flipping variables so
	  dual feasible.  If looks finished check fake dual bounds.
	  Repeat until status is iterating (-1) or finished (0,1,2)

       }

       while (status==-1) {

         Iterate until no pivot in or out or time to re-factorize.

         Flow is:

         choose pivot row (outgoing variable).  if none then
	 we are primal feasible so looks as if done but we need to
	 break and check bounds etc.

	 Get pivot row in tableau

         Choose incoming column.  If we don't find one then we look
	 primal infeasible so break and check bounds etc.  (Also the
	 pivot tolerance is larger after any iterations so that may be
	 reason)

         If we do find incoming column, we may have to adjust costs to
	 keep going forwards (anti-degeneracy).  Check pivot will be stable
	 and if unstable throw away iteration (we will need to implement
	 flagging of basic variables sometime) and break to re-factorize.
	 If minor error re-factorize after iteration.

	 Update everything (this may involve flipping variables to stay
	 dual feasible.

       }

     }

     At present we never check we are going forwards.  I overdid that in
     OSL so will try and make a last resort.

     Needs partial scan pivot out option.
     Needs dantzig, uninitialized and full steepest edge options (can still
     use partial scan)

     May need other anti-degeneracy measures, especially if we try and use
     loose tolerances as a way to solve in fewer iterations.

     I like idea of dynamic scaling.  This gives opportunity to decouple
     different implications of scaling for accuracy, iteration count and
     feasibility tolerance.

  */


  // sanity check
  assert (numberRows_==matrix_->getNumRows());
  assert (numberColumns_==matrix_->getNumCols());
  // for moment all arrays must exist
  assert(columnLower_);
  assert(columnUpper_);
  assert(rowLower_);
  assert(rowUpper_);


  algorithm_ = -1;
  dualTolerance_=dblParam_[ClpDualTolerance];
  primalTolerance_=dblParam_[ClpPrimalTolerance];

  // put in standard form (and make row copy)
  // create modifiable copies of model rim and do optional scaling
  bool goodMatrix = createRim(7+8+16,true);

  // save dual bound
  double saveDualBound_ = dualBound_;
  // save perturbation
  int savePerturbation = perturbation_;
  // save if sparse factorization wanted
  int saveSparse = factorization_->sparseThreshold();

  if (goodMatrix) {
    // Problem looks okay
    int iRow,iColumn;
    // Do initial factorization
    // and set certain stuff
    // We can either set increasing rows so ...IsBasic gives pivot row
    // or we can just increment iBasic one by one
    // for now let ...iBasic give pivot row
    factorization_->increasingRows(2);
    // row activities have negative sign
    factorization_->slackValue(-1.0);
    factorization_->zeroTolerance(1.0e-13);
    
    int factorizationStatus = internalFactorize(0);
    if (factorizationStatus<0)
      return 1; // some error
    else if (factorizationStatus)
      handler_->message(CLP_SINGULARITIES,messages_)
	<<factorizationStatus
	<<CoinMessageEol;
    
    // If user asked for perturbation - do it
    
    if (perturbation_<100) 
      perturb();
    
    double objectiveChange;
    // for dual we will change bounds using dualBound_
    // for this we need clean basis so it is after factorize
    gutsOfSolution(rowActivityWork_,columnActivityWork_);
    
    numberFake_ =0; // Number of variables at fake bounds
    changeBounds(true,NULL,objectiveChange);
    
    problemStatus_ = -1;
    numberIterations_=0;
    
    int lastCleaned=0; // last time objective or bounds cleaned up
    
    // number of times we have declared optimality
    numberTimesOptimal_=0;
    
    // Progress indicator
    ClpSimplexProgress progress(this);
    
    // This says whether to restore things etc
    int factorType=0;
    /*
      Status of problem:
      0 - optimal
      1 - infeasible
      2 - unbounded
      -1 - iterating
      -2 - factorization wanted
      -3 - redo checking without factorization
      -4 - looks infeasible
    */
    while (problemStatus_<0) {
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
#if 0
      // does not seem to work too well - do some more work
      if (perturbation_<101&&numberIterations_>2*(numberRows_+numberColumns_)) 
	perturb();
#endif
      // may factorize, checks if problem finished
      statusOfProblemInDual(lastCleaned,factorType,progress);
      
      // Say good factorization
      factorType=1;
      if (saveSparse) {
	// use default at present
	factorization_->sparseThreshold(0);
	factorization_->goSparse();
      }
      
      // Do iterations
      whileIterating();
    }
  }

  //assert(!numberFake_||problemStatus_); // all bounds should be okay
  factorization_->sparseThreshold(saveSparse);
  // Get rid of some arrays and empty factorization
  deleteRim();
  handler_->message(CLP_SIMPLEX_FINISHED+problemStatus_,messages_)
    <<objectiveValue()
    <<CoinMessageEol;
  // Restore any saved stuff
  perturbation_ = savePerturbation;
  dualBound_ = saveDualBound_;
  return problemStatus_;
}
/* Reasons to come out:
   -1 iterations etc
   -2 inaccuracy 
   -3 slight inaccuracy (and done iterations)
   +0 looks optimal (might be unbounded - but we will investigate)
   +1 looks infeasible
   +3 max iterations 
 */
int
ClpSimplexDual::whileIterating()
{
#ifdef CLP_DEBUG
  int debugIteration=-1;
#endif
  // status stays at -1 while iterating, >=0 finished, -2 to invert
  // status -3 to go to top without an invert
  int returnCode = -1;
  while (problemStatus_==-1) {
#ifdef CLP_DEBUG
    {
      int i;
      for (i=0;i<4;i++) {
	rowArray_[i]->checkClear();
      }    
      for (i=0;i<2;i++) {
	columnArray_[i]->checkClear();
      }    
    }      
#endif
#if CLP_DEBUG>2
    // very expensive
    if (numberIterations_>0&&numberIterations_<-801) {
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
      gutsOfSolution(rowActivityWork_,columnActivityWork_);
      printf("xxx %d old obj %g, recomputed %g, sum dual inf %g\n",
	     numberIterations_,
	     saveValue,objectiveValue_,sumDualInfeasibilities_);
      if (saveValue>objectiveValue_+1.0e-2)
	printf("**bad**\n");
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
#ifdef CLP_DEBUG
    {
      int iSequence, number=numberRows_+numberColumns_;
      for (iSequence=0;iSequence<number;iSequence++) {
	double lowerValue=lower_[iSequence];
	double upperValue=upper_[iSequence];
	double value=solution_[iSequence];
	if(getStatus(iSequence)!=basic) {
	  assert(lowerValue>-1.0e20);
	  assert(upperValue<1.0e20);
	}
	switch(getStatus(iSequence)) {
	  
	case basic:
	  break;
	case isFree:
	case superBasic:
	  break;
	case atUpperBound:
	  assert (fabs(value-upperValue)<=primalTolerance_) ;
	  break;
	case atLowerBound:
	  assert (fabs(value-lowerValue)<=primalTolerance_) ;
	  break;
	}
      }
    }
    if(numberIterations_==debugIteration) {
      printf("dodgy iteration coming up\n");
    }
#endif
    // choose row to go out
    // dualRow will go to virtual row pivot choice algorithm
    dualRow();
    if (pivotRow_>=0) {
      // we found a pivot row
      handler_->message(CLP_SIMPLEX_PIVOTROW,messages_)
	<<pivotRow_
	<<CoinMessageEol;
      // check accuracy of weights
      dualRowPivot_->checkAccuracy();
      // get sign for finding row of tableau
      rowArray_[0]->insert(pivotRow_,directionOut_);
      factorization_->updateColumnTranspose(rowArray_[1],rowArray_[0]);
      // put row of tableau in rowArray[0] and columnArray[0]
      matrix_->transposeTimes(this,-1.0,
			      rowArray_[0],columnArray_[1],columnArray_[0]);
      // rowArray has pi equivalent
      // do ratio test
      dualColumn(rowArray_[0],columnArray_[0],columnArray_[1],
		 rowArray_[3]);
      if (sequenceIn_>=0) {
	// normal iteration
	// update the incoming column
	unpack(rowArray_[1]);
	factorization_->updateColumn(rowArray_[2],rowArray_[1],true);
	// and update dual weights (can do in parallel - with extra array)
	dualRowPivot_->updateWeights(rowArray_[0],rowArray_[2],
				     rowArray_[1]);
	// see if update stable
	double btranAlpha = -alpha_*directionOut_; // for check
	alpha_=(*rowArray_[1])[pivotRow_];
#ifdef CLP_DEBUG
	if ((handler_->logLevel()&32))
	  printf("btran alpha %g, ftran alpha %g\n",btranAlpha,alpha_);
#endif
	if (fabs(btranAlpha)<1.0e-12||fabs(alpha_)<1.0e-12||
	    fabs(btranAlpha-alpha_)>1.0e-7*(1.0+fabs(alpha_))) {
	  handler_->message(CLP_DUAL_CHECK,messages_)
	    <<btranAlpha
	    <<alpha_
	    <<CoinMessageEol;
	  dualRowPivot_->unrollWeights();
	  if (factorization_->pivots()) {
	    problemStatus_=-2; // factorize now
	    rowArray_[0]->clear();
	    rowArray_[1]->clear();
	    columnArray_[0]->clear();
	    returnCode=-2;
	    break;
	  } else {
	    // take on more relaxed criterion
	    if (fabs(btranAlpha)<1.0e-12||fabs(alpha_)<1.0e-12||
		fabs(btranAlpha-alpha_)>1.0e-4*(1.0+fabs(alpha_))) {
	      // need to reject something
	      char x = isColumn(sequenceOut_) ? 'C' :'R';
	      handler_->message(CLP_SIMPLEX_FLAG,messages_)
		<<x<<sequenceWithin(sequenceOut_)
		<<CoinMessageEol;
	      setFlagged(sequenceOut_);
	      lastBadIteration_ = numberIterations_; // say be more cautious
	      rowArray_[0]->clear();
	      rowArray_[1]->clear();
	      columnArray_[0]->clear();
	      continue;
	    }
	  }
	}
	// update duals BEFORE replaceColumn so can do updateColumn
	double objectiveChange=0.0;
	// do duals first as variables may flip bounds
	// rowArray_[0] and columnArray_[0] may have flips
	// so use rowArray_[3] for work array from here on
	int nswapped = 
	  updateDualsInDual(rowArray_[0],columnArray_[0],rowArray_[2],theta_,
			    objectiveChange);
	// which will change basic solution
	if (nswapped) {
#ifdef CLP_DEBUG
	  if ((handler_->logLevel()&16))
	    printf("old dualOut_ %g, v %g, l %g, u %g - new ",
		   dualOut_,valueOut_,lowerOut_,upperOut_);
	  double oldOut=dualOut_;
#endif
	  factorization_->updateColumn(rowArray_[3],rowArray_[2],false);
	  dualRowPivot_->updatePrimalSolution(rowArray_[2],
					      1.0,objectiveChange);
	  
	  // recompute dualOut_
	  valueOut_ = solution_[sequenceOut_];
	  if (directionOut_<0) {
	    dualOut_ = valueOut_ - upperOut_;
	  } else {
	    dualOut_ = lowerOut_ - valueOut_;
	  }
#ifdef CLP_DEBUG
	  if ((handler_->logLevel()&16))
	    printf("%g\n",dualOut_);
	  assert(dualOut_<=oldOut);
#endif
	  if(dualOut_<0.0&&factorization_->pivots()) {
	    // going backwards - factorize
	    dualRowPivot_->unrollWeights();
	    problemStatus_=-2; // factorize now
	    returnCode=-2;
	    break;
	  }
	}
	// amount primal will move
	double movement = -dualOut_*directionOut_/alpha_;
	// so objective should increase by fabs(dj)*movement
	// but we already have objective change - so check will be good
	if (objectiveChange+fabs(movement*dualIn_)<-1.0e-5) {
#ifdef CLP_DEBUG
	  if (handler_->logLevel()&32)
	    printf("movement %g, swap change %g, rest %g  * %g\n",
		   objectiveChange+fabs(movement*dualIn_),
		   objectiveChange,movement,dualIn_);
#endif
	  if(factorization_->pivots()>5) {
	    // going backwards - factorize
	    dualRowPivot_->unrollWeights();
	    problemStatus_=-2; // factorize now
	    returnCode=-2;
	    break;
	  }
	}
	// if stable replace in basis
	int updateStatus = factorization_->replaceColumn(rowArray_[2],
							 pivotRow_,
							 alpha_);
	if (updateStatus==1) {
	  // slight error
	  if (factorization_->pivots()>5) {
	    problemStatus_=-2; // factorize now
	    returnCode=-3;
	  }
	} else if (updateStatus==2) {
	  // major error
	  dualRowPivot_->unrollWeights();
	  // later we may need to unwind more e.g. fake bounds
	  if (factorization_->pivots()) {
	    problemStatus_=-2; // factorize now
	    returnCode=-2;
	    break;
	  } else {
	    // need to reject something
	    char x = isColumn(sequenceOut_) ? 'C' :'R';
	    handler_->message(CLP_SIMPLEX_FLAG,messages_)
	      <<x<<sequenceWithin(sequenceOut_)
	      <<CoinMessageEol;
	    setFlagged(sequenceOut_);
	    lastBadIteration_ = numberIterations_; // say be more cautious
	    rowArray_[0]->clear();
	    rowArray_[1]->clear();
	    columnArray_[0]->clear();
	    // make sure dual feasible
	    // look at all rows and columns
	    CoinIotaN(rowArray_[0]->getIndices(),numberRows_,0);
	    rowArray_[0]->setNumElements(numberRows_);
	    CoinIotaN(columnArray_[0]->getIndices(),numberColumns_,0);
	    columnArray_[0]->setNumElements(numberColumns_);
	    double objectiveChange=0.0;
	    updateDualsInDual(rowArray_[0],columnArray_[0],rowArray_[1],
			      0.0,objectiveChange);
	    continue;
	  }
	} else if (updateStatus==3) {
	  // out of memory
	  // increase space if not many iterations
	  if (factorization_->pivots()<
	      0.5*factorization_->maximumPivots()&&
	      factorization_->pivots()<200)
	    factorization_->areaFactor(
				       factorization_->areaFactor() * 1.1);
	  problemStatus_=-2; // factorize now
	} 
	// update primal solution
	if (theta_<0.0) {
#ifdef CLP_DEBUG
	  if (handler_->logLevel()&32)
	    printf("negative theta %g\n",theta_);
#endif
	  theta_=0.0;
	}
	// do actual flips
	flipBounds(rowArray_[0],columnArray_[0],theta_);
	dualRowPivot_->updatePrimalSolution(rowArray_[1],
					    movement,
					    objectiveChange);
#ifdef CLP_DEBUG
	double oldobj=objectiveValue_;
#endif
	int whatNext=housekeeping(objectiveChange);
	// and set bounds correctly
	originalBound(sequenceIn_); 
	changeBound(sequenceOut_);
#ifdef CLP_DEBUG
	if (objectiveValue_<oldobj-1.0e-5&&(handler_->logLevel()&16))
	  printf("obj backwards %g %g\n",objectiveValue_,oldobj);
#endif
	if (whatNext==1) {
	  problemStatus_ =-2; // refactorize
	} else if (whatNext==2) {
	  // maximum iterations or equivalent
	  problemStatus_= 3;
	  returnCode=3;
	  break;
	}
      } else {
	// no incoming column is valid
#ifdef CLP_DEBUG
	if (handler_->logLevel()&32)
	  printf("** no column pivot\n");
#endif
	if (factorization_->pivots()<5) {
	  problemStatus_=-4; //say looks infeasible
	  // create ray anyway
	  delete [] ray_;
	  ray_ = new double [ numberRows_];
	  ClpDisjointCopyN(rowArray_[0]->denseVector(),numberRows_,ray_);
	}
	rowArray_[0]->clear();
	columnArray_[0]->clear();
	returnCode=1;
	break;
      }
    } else {
      // no pivot row
#ifdef CLP_DEBUG
      if (handler_->logLevel()&32)
	printf("** no row pivot\n");
#endif
      if (!factorization_->pivots()) {
	// may have crept through - so may be optimal
	//problemStatus_=-5; //say looks unbounded
	problemStatus_=0;
	// check any flagged variables
	int iRow;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  int iPivot=pivotVariable_[iRow];
	  if (flagged(iPivot))
	    break;
	}
	if (iRow<numberRows_) {
#ifdef CLP_DEBUG
	  std::cerr<<"Flagged variables at end - infeasible?"<<std::endl;
#endif
	  problemStatus_=-4; //say looks infeasible
	  // create ray anyway
	  delete [] ray_;
	  ray_ = new double [ numberRows_];
	  ClpDisjointCopyN(rowArray_[0]->denseVector(),numberRows_,ray_);
	}
      }
      returnCode=0;
      break;
    }
  }
  return returnCode;
}
/* The duals are updated by the given arrays.
   Returns number of infeasibilities.
   rowArray and columnarray will have flipped
   The output vector has movement (row length array) */
int
ClpSimplexDual::updateDualsInDual(CoinIndexedVector * rowArray,
				  CoinIndexedVector * columnArray,
				  CoinIndexedVector * outputArray,
				  double theta,
				  double & objectiveChange)
{
  
  outputArray->clear();
  
  double * work;
  int number;
  int * which;
  
  int numberInfeasibilities=0;
  int numberRowInfeasibilities=0;
  
  // see whether we will be doing full recompute
  bool fullRecompute= (rowArray->getNumElements()==numberRows_&&
		       columnArray->getNumElements()==numberColumns_);
  int numberAtFake=0;
  
  // use a tighter tolerance except for all being okay
  double tolerance = dualTolerance_;
  
  double changeObj=0.0;
  
  int iSection;
  
  for (iSection=0;iSection<2;iSection++) {
    int i;
    double * solution = solutionRegion(iSection);
    double * reducedCost = djRegion(iSection);
    double * lower = lowerRegion(iSection);
    double * upper = upperRegion(iSection);
    double * cost = costRegion(iSection);
    int addSequence;
    if (!iSection) {
      addSequence = numberColumns_;
      work = rowArray->denseVector();
      number = rowArray->getNumElements();
      which = rowArray->getIndices();
    } else {
      // set number of infeasibilities in row array
      addSequence=0;
      numberRowInfeasibilities=numberInfeasibilities;
      rowArray->setNumElements(numberInfeasibilities);
      numberInfeasibilities=0;
      work = columnArray->denseVector();
      number = columnArray->getNumElements();
      which = columnArray->getIndices();
    }
    
    for (i=0;i<number;i++) {
      int iSequence = which[i];
      double alphaI = work[iSequence];
      double value = reducedCost[iSequence]-theta*alphaI;
      work[iSequence]=0.0;
      reducedCost[iSequence]=value;
      
      if (!fixed(iSequence+addSequence)) {
	double movement=0.0;
	FakeBound bound = getFakeBound(iSequence+addSequence);
	Status status = getStatus(iSequence+addSequence);

	switch(status) {
	  
	case basic:
	case superBasic:
	  break;
	case isFree:
	  if (fabs(value)>tolerance) { 
#ifdef CLP_DEBUG
	    if (handler_->logLevel()&32)
	      printf("%d %d, free has dj of %g, alpha %g\n",
		     iSection,iSequence,value,alphaI);
#endif
	  }
	  break;
	case atUpperBound:
	  if (value>tolerance) {
	    // to lower bound (if swap)
	    // put back alpha
	    which[numberInfeasibilities++]=iSequence;
	    work[iSequence]=alphaI;
	    movement = lower[iSequence]-upper[iSequence];
#ifdef CLP_DEBUG
	    if ((handler_->logLevel()&32))
	      printf("%d %d, new dj %g, alpha %g, movement %g\n",
		   iSection,iSequence,value,alphaI,movement);
#endif
	    changeObj += movement*cost[iSequence];
	    if (bound==ClpSimplexDual::bothFake||
		bound==ClpSimplexDual::lowerFake) 
	      numberAtFake++;
	  } else if (fullRecompute) {
	    // at correct bound
	    if (bound==ClpSimplexDual::bothFake||
		bound==ClpSimplexDual::upperFake) {
	      // but flip if dj would allow
	      if (bound==ClpSimplexDual::upperFake&&
		  value>=-tolerance) {
		movement = lower[iSequence]-upper[iSequence];
		setStatus(iSequence+addSequence,atLowerBound);
		solution[iSequence] = lower[iSequence];
		changeObj += movement*cost[iSequence];
	      } else {
		numberAtFake++;
	      }
	    }
	  }
	  break;
	case atLowerBound:
	  if (value<-tolerance) {
	    // to upper bound 
	    // put back alpha
	    which[numberInfeasibilities++]=iSequence;
	    work[iSequence]=alphaI;
	    movement = upper[iSequence] - lower[iSequence];
#ifdef CLP_DEBUG
	    if ((handler_->logLevel()&32))
	      printf("%d %d, new dj %g, alpha %g, movement %g\n",
		   iSection,iSequence,value,alphaI,movement);
#endif
	    changeObj += movement*cost[iSequence];
	    if (bound==ClpSimplexDual::bothFake||
		bound==ClpSimplexDual::upperFake) 
	      numberAtFake++;
	  } else if (fullRecompute) {
	    // at correct bound
	    if (bound==ClpSimplexDual::bothFake||
		bound==ClpSimplexDual::lowerFake) {
	      // but flip if dj would allow
	      if (bound==ClpSimplexDual::lowerFake&&
		  value<=tolerance) {
		movement = upper[iSequence] - lower[iSequence];
		setStatus(iSequence+addSequence,atUpperBound);
		solution[iSequence] = upper[iSequence];
		changeObj += movement*cost[iSequence];
	      } else {
		numberAtFake++;
	      }
	    }
	  }
	  break;
	}
	if (!fullRecompute) {
	  if (movement) {
	    if (!iSection) {
	      // row (sign ?)
	      outputArray->quickAdd(iSequence,-movement);
	    } else {
	      matrix_->add(this,outputArray,iSequence,movement);
	    }
	  }
	}
      }
    }
  }
#ifdef CLP_DEBUG
  if (fullRecompute&&numberAtFake&&(handler_->logLevel()&16)!=0) 
    printf("%d fake after full update\n",numberAtFake);
#endif
  outputArray->stopQuickAdd();
  // set number of infeasibilities
  columnArray->setNumElements(numberInfeasibilities);
  numberInfeasibilities += numberRowInfeasibilities;
  if (fullRecompute) {
    // do actual flips
    flipBounds(rowArray,columnArray,theta);
    numberFake_ = numberAtFake;
  }
  objectiveChange += changeObj;
  return numberInfeasibilities;
}
/* 
   Chooses dual pivot row
   Would be faster with separate region to scan
   and will have this (with square of infeasibility) when steepest
   For easy problems we can just choose one of the first rows we look at
*/
void
ClpSimplexDual::dualRow()
{
  // get pivot row using whichever method it is
  pivotRow_=dualRowPivot_->pivotRow();
  if (pivotRow_>=0) {
    sequenceOut_ = pivotVariable_[pivotRow_];
    valueOut_ = solution_[sequenceOut_];
    lowerOut_ = lower_[sequenceOut_];
    upperOut_ = upper_[sequenceOut_];
    // if we have problems we could try other way and hope we get a 
    // zero pivot?
    if (valueOut_>upperOut_) {
      directionOut_ = -1;
      dualOut_ = valueOut_ - upperOut_;
    } else {
      directionOut_ = 1;
      dualOut_ = lowerOut_ - valueOut_;
    }
#ifdef CLP_DEBUG
    assert(dualOut_>=0.0);
#endif
  }
  return ;
}
// Checks if any fake bounds active - if so returns number and modifies
// dualBound_ and everything.
// Free variables will be left as free
// Returns number of bounds changed if >=0
// Returns -1 if not initialize and no effect
// Fills in changeVector which can be used to see if unbounded
// and cost of change vector
int
ClpSimplexDual::changeBounds(bool initialize,
				 CoinIndexedVector * outputArray,
				 double & changeCost)
{
  if (!initialize) {
    int numberInfeasibilities;
    double newBound;
    newBound = 5.0*dualBound_;
    numberInfeasibilities=0;
    changeCost=0.0;
    // put back original bounds and then check
    createRim(3);
    int iSequence;
    // bounds will get bigger - just look at ones at bounds
    for (iSequence=0;iSequence<numberRows_+numberColumns_;iSequence++) {
      double lowerValue=lower_[iSequence];
      double upperValue=upper_[iSequence];
      double value=solution_[iSequence];
      setFakeBound(iSequence,ClpSimplexDual::noFake);
      switch(getStatus(iSequence)) {

      case basic:
	break;
      case isFree:
      case superBasic:
	break;
      case atUpperBound:
	if (fabs(value-upperValue)>primalTolerance_) 
	  numberInfeasibilities++;
	break;
      case atLowerBound:
	if (fabs(value-lowerValue)>primalTolerance_) 
	  numberInfeasibilities++;
	break;
      }
    }
    if (numberInfeasibilities) {
      int iSequence;
      for (iSequence=0;iSequence<numberRows_+numberColumns_;iSequence++) {
	double lowerValue=lower_[iSequence];
	double upperValue=upper_[iSequence];
	double newLowerValue;
	double newUpperValue;
	Status status = getStatus(iSequence);
	if (status==atUpperBound||
	    status==atLowerBound) {
	  double value = solution_[iSequence];
	  if (value-lowerValue<=upperValue-value) {
	    newLowerValue = max(lowerValue,value-0.666667*newBound);
	    newUpperValue = min(upperValue,newLowerValue+newBound);
	  } else {
	    newUpperValue = min(upperValue,value+0.666667*newBound);
	    newLowerValue = max(lowerValue,newUpperValue-newBound);
	  }
	  lower_[iSequence]=newLowerValue;
	  upper_[iSequence]=newUpperValue;
	  if (newLowerValue > lowerValue) {
	    if (newUpperValue < upperValue) 
	      setFakeBound(iSequence,ClpSimplexDual::bothFake);
	    else
	      setFakeBound(iSequence,ClpSimplexDual::lowerFake);
	  } else {
	    if (newUpperValue < upperValue) 
	      setFakeBound(iSequence,ClpSimplexDual::upperFake);
	  }
	  if (status==atUpperBound)
	    solution_[iSequence] = newUpperValue;
	  else
	    solution_[iSequence] = newLowerValue;
	  double movement = solution_[iSequence] - value;
	  if (movement&&outputArray) {
	    if (iSequence>=numberColumns_) {
	      outputArray->quickAdd(iSequence,-movement);
	      changeCost += movement*cost_[iSequence];
	    } else {
	      matrix_->add(this,outputArray,iSequence,movement);
	      changeCost += movement*cost_[iSequence];
	    }
	  }
	}
      }
      dualBound_ = newBound;
      if (outputArray)
	outputArray->stopQuickAdd();
    } else {
      numberInfeasibilities=-1;
    }
    return numberInfeasibilities;
  } else {
    int iSequence;
      
    for (iSequence=0;iSequence<numberRows_+numberColumns_;iSequence++) {
      Status status = getStatus(iSequence);
      if (status==atUpperBound||
	  status==atLowerBound) {
	double lowerValue=lower_[iSequence];
	double upperValue=upper_[iSequence];
	double value = solution_[iSequence];
	if (lowerValue>-largeValue_||upperValue<largeValue_) {
	  if (lowerValue-value>-0.5*dualBound_||
	      upperValue-value<0.5*dualBound_) {
	    if (fabs(lowerValue-value)<=fabs(upperValue-value)) {
	      if (upperValue > lowerValue + dualBound_) {
		upper_[iSequence]=lowerValue+dualBound_;
		setFakeBound(iSequence,ClpSimplexDual::upperFake);
	      }
	    } else {
	      if (lowerValue < upperValue - dualBound_) {
		lower_[iSequence]=upperValue-dualBound_;
		setFakeBound(iSequence,ClpSimplexDual::lowerFake);
	      }
	    }
	  } else {
	    lower_[iSequence]=-0.5*dualBound_;
	    upper_[iSequence]= 0.5*dualBound_;
	    setFakeBound(iSequence,ClpSimplexDual::bothFake);
	  }
	}
      }
    }
    return 1;
  }
}
/* 
   Row array has row part of pivot row (as duals so sign may be switched)
   Column array has column part.
   This chooses pivot column.
   Spare array will be needed when we start getting clever.
   We will check for basic so spare array will never overflow.
   If necessary will modify costs
*/
void
ClpSimplexDual::dualColumn(CoinIndexedVector * rowArray,
		       CoinIndexedVector * columnArray,
		       CoinIndexedVector * spareArray,
		       CoinIndexedVector * spareArray2)
{
  double * work;
  int number;
  int * which;
  double * reducedCost;

  int iSection;

  sequenceIn_=-1;
  int numberPossiblySwapped=0;
  int numberRemaining=0;
  
  double totalThru=0.0; // for when variables flip
  double acceptablePivot=1.0e-7;
  if (factorization_->pivots())
    acceptablePivot=1.0e-5; // if we have iterated be more strict
  double bestEverPivot=acceptablePivot;
  int lastSequence = -1;
  double lastPivot=0.0;
  double upperTheta;
  double newTolerance = dualTolerance_;
  // will we need to increase tolerance
  bool thisIncrease=false;
  // If we think we need to modify costs (not if something from broad sweep)
  bool modifyCosts=false;
  // Increase in objective due to swapping bounds (may be negative)
  double increaseInObjective=0.0;

  // use spareArrays to put ones looked at in
  // we are going to flip flop between
  int iFlip = 0;
  // Possible list of pivots
  int interesting[2];
  // where possible swapped ones are
  int swapped[2];
  // for zeroing out arrays after
  int marker[2][2];
  // pivot elements
  double * array[2], * spare, * spare2;
  // indices
  int * indices[2], * index, * index2;
  spareArray->clear();
  spareArray2->clear();
  array[0] = spareArray->denseVector();
  indices[0] = spareArray->getIndices();
  spare = array[0];
  index = indices[0];
  array[1] = spareArray2->denseVector();
  indices[1] = spareArray2->getIndices();
  int i;
  double * lower;
  double * upper;

  // initialize lists
  for (i=0;i<2;i++) {
    interesting[i]=0;
    swapped[i]=numberColumns_;
    marker[i][0]=0;
    marker[i][1]=numberColumns_;
  }

  /*
    First we get a list of possible pivots.  We can also see if the
    problem looks infeasible or whether we want to pivot in free variable.
    This may make objective go backwards but can only happen a finite
    number of times and I do want free variables basic.

    Then we flip back and forth.  At the start of each iteration
    interesting[iFlip] should have possible candidates and swapped[iFlip]
    will have pivots if we decide to take a previous pivot.
    At end of each iteration interesting[1-iFlip] should have
    candidates if we go through this theta and swapped[1-iFlip]
    pivots if we don't go through.

    At first we increase theta and see what happens.  We start
    theta at a reasonable guess.  If in right area then we do bit by bit.

   */

  // do first pass to get possibles 
  // We can also see if infeasible or pivoting on free
  double tentativeTheta = 1.0e22;
  upperTheta = 1.0e31;
  double freePivot = acceptablePivot;
  for (iSection=0;iSection<2;iSection++) {

    int addSequence;

    if (!iSection) {
      lower = rowLowerWork_;
      upper = rowUpperWork_;
      work = rowArray->denseVector();
      number = rowArray->getNumElements();
      which = rowArray->getIndices();
      reducedCost = rowReducedCost_;
      addSequence = numberColumns_;
    } else {
      lower = columnLowerWork_;
      upper = columnUpperWork_;
      work = columnArray->denseVector();
      number = columnArray->getNumElements();
      which = columnArray->getIndices();
      reducedCost = reducedCostWork_;
      addSequence = 0;
    }
    
    for (i=0;i<number;i++) {
      int iSequence = which[i];
      double alpha = work[iSequence];
      if (fixed(iSequence+addSequence)||!alpha) 
	continue; // skip fixed ones or (zeroed out)
      double oldValue = reducedCost[iSequence];
      double value = oldValue-tentativeTheta*alpha;
      int keep = 0;

      switch(getStatus(iSequence+addSequence)) {
	  
      case basic:
	break;
      case isFree:
      case superBasic:
	if (oldValue>dualTolerance_) {
	  if (value<-newTolerance) 
	    keep = 2;
	} else if (oldValue<-dualTolerance_) {
	  if (value>newTolerance) 
	    keep = 2;
	} else {
	  if (alpha>=acceptablePivot) 
	    keep = 2;
	  else if (-alpha>=acceptablePivot) 
	    keep = 2;
	}
	break;
      case atUpperBound:
	assert (oldValue<=dualTolerance_*1.0001);
	if (value>newTolerance) {
	  keep = 1;
	  value = oldValue-upperTheta*alpha;
	  if (value>newTolerance && -alpha>=acceptablePivot) 
	    upperTheta = (oldValue-newTolerance)/alpha;
	}
	break;
      case atLowerBound:
	assert (oldValue>=-dualTolerance_*1.0001);
	if (value<-newTolerance) {
	  keep = 1;
	  value = oldValue-upperTheta*alpha;
	  if (value<-newTolerance && alpha>=acceptablePivot) 
	    upperTheta = (oldValue+newTolerance)/alpha;
	}
	break;
      }
      if (keep) {
	if (keep==2) {
	  // free - choose largest
	  if (fabs(alpha)>freePivot) {
	    freePivot=fabs(alpha);
	    sequenceIn_ = iSequence + addSequence;
	    theta_=oldValue/alpha;
	  }
	} else {
	  // add to list
	  spare[numberRemaining]=alpha;
	  index[numberRemaining++]=iSequence+addSequence;
	}
      }
    }
  }
  interesting[0]=numberRemaining;
  marker[0][0] = numberRemaining;

  if (!numberRemaining)
    return; // Looks infeasible

  if (sequenceIn_>=0) {
    // free variable - always choose
  } else {

    theta_=1.0e50;
    // now flip flop between spare arrays until reasonable theta
    tentativeTheta = max(10.0*upperTheta,1.0e-7);
    
    // loops increasing tentative theta until can't go through
    
    while (tentativeTheta < 1.0e22) {
      double thruThis = 0.0;
      
      double bestPivot=acceptablePivot;
      int bestSequence=-1;
      
      numberPossiblySwapped = numberColumns_;
      numberRemaining = 0;

      upperTheta = 1.0e50;

      spare = array[iFlip];
      index = indices[iFlip];
      spare2 = array[1-iFlip];
      index2 = indices[1-iFlip];

      // try 3 different ways
      // 1 bias increase by ones with slightly wrong djs
      // 2 bias by all
      // 3 bias by all - tolerance (doesn't seem very good)
#define TRYBIAS 1


      double increaseInThis=0.0; //objective increase in this loop
      
      for (i=0;i<interesting[iFlip];i++) {
	int iSequence = index[i];
	double alpha = spare[i];
	double oldValue = dj_[iSequence];
	double value = oldValue-tentativeTheta*alpha;

	if (alpha < 0.0) {
	  //at upper bound
	  if (value>newTolerance) {
	    double range = upper_[iSequence] - lower_[iSequence];
	    thruThis -= range*alpha;
#if TRYBIAS==1
	    if (oldValue>0.0)
	      increaseInThis -= oldValue*range;
#elif TRYBIAS==2
	    increaseInThis -= oldValue*range;
#else
	    increaseInThis -= (oldValue+dualTolerance_)*range;
#endif
	    // goes on swapped list (also means candidates if too many)
	    spare2[--numberPossiblySwapped]=alpha;
	    index2[numberPossiblySwapped]=iSequence;
	    if (fabs(alpha)>bestPivot) {
	      bestPivot=fabs(alpha);
	      bestSequence=numberPossiblySwapped;
	    }
	  } else {
	    value = oldValue-upperTheta*alpha;
	    if (value>newTolerance && -alpha>=acceptablePivot) 
	      upperTheta = (oldValue-newTolerance)/alpha;
	    spare2[numberRemaining]=alpha;
	    index2[numberRemaining++]=iSequence;
	  }
	} else {
	  // at lower bound
	  if (value<-newTolerance) {
	    double range = upper_[iSequence] - lower_[iSequence];
	    thruThis += range*alpha;
	    //?? is this correct - and should we look at good ones
#if TRYBIAS==1
	    if (oldValue<0.0)
	      increaseInThis += oldValue*range;
#elif TRYBIAS==2
	    increaseInThis += oldValue*range;
#else
	    increaseInThis += (oldValue-dualTolerance_)*range;
#endif
	    // goes on swapped list (also means candidates if too many)
	    spare2[--numberPossiblySwapped]=alpha;
	    index2[numberPossiblySwapped]=iSequence;
	    if (fabs(alpha)>bestPivot) {
	      bestPivot=fabs(alpha);
	      bestSequence=numberPossiblySwapped;
	    }
	  } else {
	    value = oldValue-upperTheta*alpha;
	    if (value<-newTolerance && alpha>=acceptablePivot) 
	      upperTheta = (oldValue+newTolerance)/alpha;
	    spare2[numberRemaining]=alpha;
	    index2[numberRemaining++]=iSequence;
	  }
	}
      }
      swapped[1-iFlip]=numberPossiblySwapped;
      interesting[1-iFlip]=numberRemaining;
      marker[1-iFlip][0]= max(marker[1-iFlip][0],numberRemaining);
      marker[1-iFlip][1]= min(marker[1-iFlip][1],numberPossiblySwapped);
      
      if (totalThru+thruThis>=fabs(dualOut_)||
	  increaseInObjective+increaseInThis<0.0) {
	// We should be pivoting in this batch
	// so compress down to this lot
	numberRemaining=0;
	for (i=numberColumns_-1;i>=swapped[1-iFlip];i--) {
	  spare[numberRemaining]=spare2[i];
	  index[numberRemaining++]=index2[i];
	}
	interesting[iFlip]=numberRemaining;
	int iTry;
#define MAXTRY 100
	// first get ratio with tolerance
	for (iTry=0;iTry<MAXTRY;iTry++) {
	  
	  upperTheta=1.0e50;
	  numberPossiblySwapped = numberColumns_;
	  numberRemaining = 0;

	  increaseInThis=0.0; //objective increase in this loop

	  thruThis=0.0;

	  spare = array[iFlip];
	  index = indices[iFlip];
	  spare2 = array[1-iFlip];
	  index2 = indices[1-iFlip];
      
	  for (i=0;i<interesting[iFlip];i++) {
	    int iSequence=index[i];
	    double alpha=spare[i];
	    double oldValue = dj_[iSequence];
	    double value = oldValue-upperTheta*alpha;
	    
	    if (alpha < 0.0) {
	      //at upper bound
	      if (value>newTolerance) {
		if (-alpha>=acceptablePivot) {
		  upperTheta = (oldValue-newTolerance)/alpha;
		  // recompute value and make sure works
		  value = oldValue-upperTheta*alpha;
		  if (value<0)
		    upperTheta *= 1.0 +1.0e-11; // must be large
		}
	      }
	    } else {
	      // at lower bound
	      if (value<-newTolerance) {
		if (alpha>=acceptablePivot) {
		  upperTheta = (oldValue+newTolerance)/alpha;
		  // recompute value and make sure works
		  value = oldValue-upperTheta*alpha;
		  if (value>0)
		    upperTheta *= 1.0 +1.0e-11; // must be large
		}
	      }
	    }
	  }
	  bestPivot=acceptablePivot;
	  sequenceIn_=-1;
	  // now choose largest and sum all ones which will go through
#define MINIMUMTHETA 1.0e-12
	  for (i=0;i<interesting[iFlip];i++) {
	    int iSequence=index[i];
	    double alpha=spare[i];
	    double value = dj_[iSequence]-upperTheta*alpha;
	    double badDj=0.0;
	    
	    bool addToSwapped=false;
	    
	    if (alpha < 0.0) {
	      //at upper bound
	      if (value>=0.0) { 
		addToSwapped=true;
#if TRYBIAS==1
		badDj = -max(dj_[iSequence],0.0);
#elif TRYBIAS==2
		badDj = -dj_[iSequence];
#else
		badDj = -dj_[iSequence]-dualTolerance_;
#endif
	      }
	    } else {
	      // at lower bound
	      if (value<=0.0) {
		addToSwapped=true;
#if TRYBIAS==1
		badDj = min(dj_[iSequence],0.0);
#elif TRYBIAS==2
		badDj = dj_[iSequence];
#else
		badDj = dj_[iSequence]-dualTolerance_;
#endif
	      }
	    }
	    if (!addToSwapped) {
	      // add to list of remaining
	      spare2[numberRemaining]=alpha;
	      index2[numberRemaining++]=iSequence;
	    } else {
	      // add to list of swapped
	      spare2[--numberPossiblySwapped]=alpha;
	      index2[numberPossiblySwapped]=iSequence;
	      // select if largest pivot
	      if (fabs(alpha)>bestPivot) {
		sequenceIn_ = numberPossiblySwapped;
		bestPivot =  fabs(alpha);
		theta_ = dj_[iSequence]/alpha;
	      }
	      double range = upper[iSequence] - lower[iSequence];
	      thruThis += range*fabs(alpha);
	      increaseInThis += badDj*range;
	    }
	  }
	  swapped[1-iFlip]=numberPossiblySwapped;
	  interesting[1-iFlip]=numberRemaining;
	  marker[1-iFlip][0]= max(marker[1-iFlip][0],numberRemaining);
	  marker[1-iFlip][1]= min(marker[1-iFlip][1],numberPossiblySwapped);
	  // If we stop now this will be increase in objective (I think)
	  double increase = (fabs(dualOut_)-totalThru)*theta_;
	  increase += increaseInObjective;
	  if (theta_<0.0)
	    thruThis += fabs(dualOut_); // force using this one
	  if (increaseInObjective<0.0&&increase<0.0&&lastSequence>=0) {
	    // back
	    // We may need to be more careful - we could do by
	    // switch so we always do fine grained?
	    bestPivot=0.0;
	  } else {
	    // add in
	    totalThru += thruThis;
	    increaseInObjective += increaseInThis;
	  }
	  if (bestPivot<0.1*bestEverPivot&&
	      bestEverPivot>1.0e-6&&bestPivot<1.0e-3) {
	    // back to previous one
	    sequenceIn_=lastSequence;
	    // swap regions
	    iFlip = 1-iFlip;
	    break;
	  } else if (sequenceIn_==-1&&upperTheta>largeValue_) {
	    if (lastPivot>acceptablePivot) {
	      // back to previous one
	      sequenceIn_=lastSequence;
	      // swap regions
	      iFlip = 1-iFlip;
	    } else {
	      // can only get here if all pivots too small
	    }
	    break;
	  } else if (totalThru>=fabs(dualOut_)) {
	    modifyCosts=true; // fine grain - we can modify costs
	    break; // no point trying another loop
	  } else {
	    lastSequence=sequenceIn_;
	    if (bestPivot>bestEverPivot)
	      bestEverPivot=bestPivot;
	    iFlip = 1 -iFlip;
	    modifyCosts=true; // fine grain - we can modify costs
	}
	}
	if (iTry==MAXTRY)
	  iFlip = 1-iFlip; // flip back
	break;
      } else {
	// skip this lot
	if (bestPivot>1.0e-3||bestPivot>bestEverPivot) {
	  bestEverPivot=bestPivot;
	  lastSequence=bestSequence;
	} else {
	  // keep old swapped
	  memcpy(array[1-iFlip]+swapped[iFlip],
		 array[iFlip]+swapped[iFlip],
		 (numberColumns_-swapped[iFlip])*sizeof(double));
	  memcpy(indices[1-iFlip]+swapped[iFlip],
		 indices[iFlip]+swapped[iFlip],
		 (numberColumns_-swapped[iFlip])*sizeof(int));
	  marker[1-iFlip][1] = min(marker[1-iFlip][1],swapped[iFlip]);
	  swapped[1-iFlip]=swapped[iFlip];
	}
	increaseInObjective += increaseInThis;
	iFlip = 1 - iFlip; // swap regions
	tentativeTheta = 2.0*upperTheta;
	totalThru += thruThis;
      }
    }
    
    // can get here without sequenceIn_ set but with lastSequence
    if (sequenceIn_<0&&lastSequence>=0) {
      // back to previous one
      sequenceIn_=lastSequence;
      // swap regions
      iFlip = 1-iFlip;
    }
    
    if (sequenceIn_>=0) {
      // at this stage sequenceIn_ is just pointer into index array
      // flip just so we can use iFlip
      iFlip = 1 -iFlip;
      spare = array[iFlip];
      index = indices[iFlip];
      double oldValue;
      double alpha = spare[sequenceIn_];
      sequenceIn_ = indices[iFlip][sequenceIn_];
      oldValue = dj_[sequenceIn_];
      theta_ = oldValue/alpha;
#if 0
      if (numberIterations_>2000)
	exit(99);
      if (numberIterations_>2000-20) 
	handler_->setLogLevel(63);
      if (numberIterations_>2000-20) 
	printf("theta %g oldValue %g tol %g %g\n",theta_,oldValue,dualTolerance_,
	       newTolerance);
#endif
      if (theta_<MINIMUMTHETA) {
	// can't pivot to zero
	if (oldValue-MINIMUMTHETA*alpha>=-dualTolerance_) {
	  theta_=MINIMUMTHETA;
	} else if (oldValue-MINIMUMTHETA*alpha>=-newTolerance) {
	  theta_=MINIMUMTHETA;
	  thisIncrease=true;
	} else {
	  theta_=(oldValue+newTolerance)/alpha;
	  assert(theta_>=0.0);
	  thisIncrease=true;
	} 
      }
      // may need to adjust costs so all dual feasible AND pivoted is exactly 0
      if (modifyCosts) {
	int i;
	double * workRow = rowArray->denseVector();
	double * workColumn = columnArray->denseVector();
	for (i=numberColumns_-1;i>=swapped[iFlip];i--) {
	  int iSequence=index[i];
	  double alpha;
	  if (iSequence>=numberColumns_)
	    alpha=workRow[iSequence-numberColumns_];
	  else
	    alpha=workColumn[iSequence];
	  double value = dj_[iSequence]-theta_*alpha;
#if 0
	  if (numberIterations_>2000-20) 
	    printf("%d alpha %g value %g\n",iSequence,alpha,value);
#endif
	    
	  // can't be free here
	  
	  if (alpha < 0.0) {
	    //at upper bound
	    if (value>dualTolerance_) {
	      thisIncrease=true;
#define MODIFYCOST 2
#if MODIFYCOST
	      // modify cost to hit new tolerance
	      double modification = alpha*theta_-dj_[iSequence]
		+newTolerance;
	      //modification = min(modification,dualTolerance_);
	      //assert (fabs(modification)<1.0e-7);
	      dj_[iSequence] += modification;
	      cost_[iSequence] +=  modification;
#if 0
	      if (numberIterations_>2000-20) 
		printf("%d acost %g mod %g\n",iSequence,cost_[iSequence],
		       modification);
#endif
#endif
	    }
	  } else {
	    // at lower bound
	    if (-value>dualTolerance_) {
	      thisIncrease=true;
#if MODIFYCOST
	      // modify cost to hit new tolerance
	      double modification = alpha*theta_-dj_[iSequence]
		-newTolerance;
	      //modification = max(modification,-dualTolerance_);
	      //assert (fabs(modification)<1.0e-7);
	      dj_[iSequence] += modification;
	      cost_[iSequence] +=  modification;
#if 0
	      if (numberIterations_>2000-20) 
		printf("%d bcost %g mod %g\n",iSequence,cost_[iSequence],
		       modification);
#endif
#endif
	    }
	  }
	}
      }
    }
  }

  if (sequenceIn_>=0) {
    if (sequenceIn_>=numberColumns_) {
      //slack
      alpha_ = rowArray->denseVector()[sequenceIn_-numberColumns_];
    } else {
      // column
      alpha_ = columnArray->denseVector()[sequenceIn_];
    }
    lowerIn_ = lower_[sequenceIn_];
    upperIn_ = upper_[sequenceIn_];
    valueIn_ = solution_[sequenceIn_];
    dualIn_ = dj_[sequenceIn_];

    if (numberTimesOptimal_) {
      // can we adjust cost back closer to original
      //*** add coding
    }
#if MODIFYCOST>1    
    // modify cost to hit zero exactly
    // so (dualIn_+modification)==theta_*alpha_
    double modification = theta_*alpha_-dualIn_;
    dualIn_ += modification;
    dj_[sequenceIn_]=dualIn_;
    cost_[sequenceIn_] += modification;
    //assert (fabs(modification)<1.0e-6);
#ifdef CLP_DEBUG
    if ((handler_->logLevel()&32)&&fabs(modification)>1.0e-15)
      printf("exact %d new cost %g, change %g\n",sequenceIn_,
	     cost_[sequenceIn_],modification);
#endif
#endif
    
    if (alpha_<0.0) {
      // as if from upper bound
      directionIn_=-1;
      upperIn_=valueIn_;
    } else {
      // as if from lower bound
      directionIn_=1;
      lowerIn_=valueIn_;
    }
  }
  if (thisIncrease) {
    newTolerance = dualTolerance_+1.0e-4*dblParam_[ClpDualTolerance];
  }

  // clear arrays

  for (i=0;i<2;i++) {
    memset(array[i],0,marker[i][0]*sizeof(double));
    memset(array[i]+marker[i][1],0,
	   (numberColumns_-marker[i][1])*sizeof(double));
  }
}
/* Checks if tentative optimal actually means unbounded
   Returns -3 if not, 2 if is unbounded */
int 
ClpSimplexDual::checkUnbounded(CoinIndexedVector * ray,
				   CoinIndexedVector * spare,
				   double changeCost)
{
  int status=2; // say unbounded
  factorization_->updateColumn(spare,ray);
  // get reduced cost
  int i;
  int number=ray->getNumElements();
  int * index = ray->getIndices();
  double * array = ray->denseVector();
  for (i=0;i<number;i++) {
    int iRow=index[i];
    int iPivot=pivotVariable_[iRow];
    changeCost -= cost(iPivot)*array[iRow];
  }
  double way;
  if (changeCost>0.0) {
    //try going down
    way=1.0;
  } else if (changeCost<0.0) {
    //try going up
    way=-1.0;
  } else {
#ifdef CLP_DEBUG
    printf("can't decide on up or down\n");
#endif
    way=0.0;
    status=-3;
  }
  double movement=1.0e10*way; // some largish number
  double zeroTolerance = 1.0e-14*dualBound_;
  for (i=0;i<number;i++) {
    int iRow=index[i];
    int iPivot=pivotVariable_[iRow];
    double arrayValue = array[iRow];
    if (fabs(arrayValue)<zeroTolerance)
      arrayValue=0.0;
    double newValue=solution(iPivot)+movement*arrayValue;
    if (newValue>upper(iPivot)+primalTolerance_||
	newValue<lower(iPivot)-primalTolerance_)
      status=-3; // not unbounded
  }
  if (status==2) {
    // create ray
    delete [] ray_;
    ray_ = new double [numberColumns_];
    ClpFillN(ray_,numberColumns_,0.0);
    for (i=0;i<number;i++) {
      int iRow=index[i];
      int iPivot=pivotVariable_[iRow];
      double arrayValue = array[iRow];
      if (iPivot<numberColumns_&&fabs(arrayValue)>=zeroTolerance)
	ray_[iPivot] = way* array[iRow];
    }
  }
  ray->clear();
  return status;
}
/* Checks if finished.  Updates status */
void 
ClpSimplexDual::statusOfProblemInDual(int & lastCleaned,int type,
			       ClpSimplexProgress &progress)
{
  if (type==2) {
    // trouble - restore solution
    memcpy(status_ ,saveStatus_,(numberColumns_+numberRows_)*sizeof(char));
    memcpy(rowActivityWork_,savedSolution_+numberColumns_ ,
	   numberRows_*sizeof(double));
    memcpy(columnActivityWork_,savedSolution_ ,
	   numberColumns_*sizeof(double));
    forceFactorization_=1; // a bit drastic but ..
    changeMade_++; // say something changed
  }
  int tentativeStatus = problemStatus_;
  double changeCost;

  if (problemStatus_>-3) {
    // factorize
    // later on we will need to recover from singularities
    // also we could skip if first time
    // save dual weights
    dualRowPivot_->saveWeights(this,1);
    if (type) {
      // is factorization okay?
      if (internalFactorize(1)) {
	// no - restore previous basis
	assert (type==1);
	changeMade_++; // say something changed
	memcpy(status_ ,saveStatus_,(numberColumns_+numberRows_)*sizeof(char));
	memcpy(rowActivityWork_,savedSolution_+numberColumns_ ,
	       numberRows_*sizeof(double));
	memcpy(columnActivityWork_,savedSolution_ ,
	       numberColumns_*sizeof(double));
	// get correct bounds on all variables
	double dummyChangeCost=0.0;
	changeBounds(true,rowArray_[2],dummyChangeCost);
	// throw away change
	rowArray_[2]->clear();
	forceFactorization_=1; // a bit drastic but ..
	type = 2;
	assert (internalFactorize(1)==0);
      }
    }
    problemStatus_=-3;
  }
  // at this stage status is -3 or -4 if looks infeasible
  // get primal and dual solutions
  gutsOfSolution(rowActivityWork_,columnActivityWork_);
  // Check if looping
  int loop = progress.looping();
  bool situationChanged=false;
  if (loop>=0) {
    problemStatus_ = loop; //exit if in loop
    return;
  } else if (loop<-1) {
    // something may have changed
    gutsOfSolution(rowActivityWork_,columnActivityWork_);
    situationChanged=true;
  }
  progressFlag_ = 0; //reset progress flag
#ifdef CLP_DEBUG
  if (!rowScale_&&(handler_->logLevel()&32)) {
    double * objectiveSimplex 
      = ClpCopyOfArray(objective_,numberColumns_,0.0);
    double * rowObjectiveSimplex 
      = ClpCopyOfArray(rowObjective_,numberRows_,0.0);
    int i;
    double largest;
    largest=0.0;
    for (i=0;i<numberRows_;i++) {
      rowObjectiveSimplex[i] *= optimizationDirection_;
      double difference = fabs(rowObjectiveWork_[i]-rowObjectiveSimplex[i]);
      if (difference>largest)
	largest=difference;
    }
    for (i=0;i<numberColumns_;i++) {
      objectiveSimplex[i] *= optimizationDirection_;
      double difference = fabs(objectiveWork_[i]-objectiveSimplex[i]);
      if (difference>largest)
	largest=difference;
    }
    if ((handler_->logLevel()&16))
      printf("difference in obj %g\n",largest);
    delete [] objectiveSimplex;
    delete [] rowObjectiveSimplex;
  }
#endif
  handler_->message(CLP_SIMPLEX_STATUS,messages_)
    <<numberIterations_<<objectiveValue();
  handler_->printing(sumPrimalInfeasibilities_>0.0)
    <<sumPrimalInfeasibilities_<<numberPrimalInfeasibilities_;
  handler_->printing(sumDualInfeasibilities_>0.0)
    <<sumDualInfeasibilities_<<numberDualInfeasibilities_;
  handler_->printing(numberDualInfeasibilitiesWithoutFree_
		     <numberDualInfeasibilities_)
		       <<numberDualInfeasibilities_-
    numberDualInfeasibilitiesWithoutFree_;
  handler_->message()<<CoinMessageEol;

  while (problemStatus_<=-3) {
    bool cleanDuals=situationChanged;
    int numberChangedBounds=0;
    int doOriginalTolerance=0;
    if ( lastCleaned==numberIterations_)
      doOriginalTolerance=1;
    // check optimal
    // give code benefit of doubt
    if (sumOfRelaxedDualInfeasibilities_ == 0.0&&
	sumOfRelaxedPrimalInfeasibilities_ == 0.0) {
      // say optimal (with these bounds etc)
      numberDualInfeasibilities_ = 0;
      sumDualInfeasibilities_ = 0.0;
      numberPrimalInfeasibilities_ = 0;
      sumPrimalInfeasibilities_ = 0.0;
    }
    if (dualFeasible()||problemStatus_==-4) {
      if (primalFeasible()) {
	// may be optimal - or may be bounds are wrong
	handler_->message(CLP_DUAL_BOUNDS,messages_)
	  <<dualBound_
	  <<CoinMessageEol;
	// save solution in case unbounded
	ClpDisjointCopyN(columnActivityWork_,numberColumns_,
			  columnArray_[0]->denseVector());
	ClpDisjointCopyN(rowActivityWork_,numberRows_,
			  rowArray_[2]->denseVector());
	numberChangedBounds=changeBounds(false,rowArray_[0],changeCost);
	if (numberChangedBounds<=0) {
	  //looks optimal - do we need to reset tolerance
	  if (perturbation_==101) {
	    perturbation_=102; // stop any perturbations
	    cleanDuals=true;
	    createRim(4);
	  }
	  if (lastCleaned<numberIterations_&&numberTimesOptimal_<4) {
	    doOriginalTolerance=2;
	    numberTimesOptimal_++;
	    changeMade_++; // say something changed
	    if (numberTimesOptimal_==1) {
	      dualTolerance_ = min(dualTolerance_,1.0e-8);
	      // better to have small tolerance even if slower
	      factorization_->zeroTolerance(1.0e-15);
	    }
	    cleanDuals=true;
	  } else {
	    problemStatus_=0; // optimal
	    if (lastCleaned<numberIterations_) {
	      handler_->message(CLP_SIMPLEX_GIVINGUP,messages_)
		<<CoinMessageEol;
	    }
	  }
	} else {
	  cleanDuals=true;
	  if (doOriginalTolerance==1) {
	    // check unbounded
	    problemStatus_ = checkUnbounded(rowArray_[0],rowArray_[1],
					    changeCost);
	    if (problemStatus_==2&&perturbation_==101) {
	      perturbation_=102; // stop any perturbations
	      cleanDuals=true;
	      createRim(4);
	      problemStatus_=-1;
	    }
	    if (problemStatus_==2) {
	      // it is unbounded - restore solution
	      // but first add in changes to non-basic
	      int iColumn;
	      double * original = columnArray_[0]->denseVector();
	      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
		if(getColumnStatus(iColumn)!= basic)
		  ray_[iColumn] += 
		    columnActivityWork_[iColumn]-original[iColumn];
		columnActivityWork_[iColumn] = original[iColumn];
	      }
	      ClpDisjointCopyN(rowArray_[2]->denseVector(),numberRows_,
				rowActivityWork_);
	    }
	  } else {
	    doOriginalTolerance=2;
	    rowArray_[0]->clear();
	  }
	}
	ClpFillN(columnArray_[0]->denseVector(),numberColumns_,0.0);
	ClpFillN(rowArray_[2]->denseVector(),numberRows_,0.0);
      } 
      if (problemStatus_==-4) {
	// may be infeasible - or may be bounds are wrong
	handler_->message(CLP_DUAL_CHECKB,messages_)
	  <<dualBound_
	  <<CoinMessageEol;
	numberChangedBounds=changeBounds(false,NULL,changeCost);
	if (perturbation_==101) {
	  perturbation_=102; // stop any perturbations
	  cleanDuals=true;
	  numberChangedBounds=1;
	  createRim(4);
	}
	if (numberChangedBounds<=0||dualBound_>1.0e20||
	    (largestPrimalError_>1.0&&dualBound_>1.0e17)) {
	  problemStatus_=1; // infeasible
	} else {
	  problemStatus_=-1; //iterate
	  cleanDuals=true;
	  doOriginalTolerance=2;
	  // and delete ray which has been created
	  delete [] ray_;
	  ray_ = NULL;
	}

      }
    } else {
      cleanDuals=true;
    }
    if (problemStatus_<0) {
      if (doOriginalTolerance==2) {
	// put back original tolerance
	lastCleaned=numberIterations_;
	handler_->message(CLP_DUAL_ORIGINAL,messages_)
	  <<CoinMessageEol;

	perturbation_=102; // stop any perturbations
	createRim(4);
	// make sure duals are current
	computeDuals();
	// put back bounds as they were if was optimal
	if (doOriginalTolerance==2) {
	  changeMade_++; // say something changed
	  changeBounds(true,NULL,changeCost);
	  cleanDuals=true;
	}
      }
      if (cleanDuals) {
	// make sure dual feasible
	// look at all rows and columns
	rowArray_[0]->clear();
	CoinIotaN(rowArray_[0]->getIndices(),numberRows_,0);
	rowArray_[0]->setNumElements(numberRows_);
	columnArray_[0]->clear();
	CoinIotaN(columnArray_[0]->getIndices(),numberColumns_,0);
	columnArray_[0]->setNumElements(numberColumns_);
	double objectiveChange=0.0;
	updateDualsInDual(rowArray_[0],columnArray_[0],rowArray_[1],
			  0.0,objectiveChange);
	// for now - recompute all
	gutsOfSolution(rowActivityWork_,columnActivityWork_);
	assert(numberDualInfeasibilitiesWithoutFree_==0);
	if (numberDualInfeasibilities_) {
	  // bad free variables
	  if (primalFeasible()) {
	    std::cerr<<"Free variable problem?"<<std::endl;
	    abort(); // what now
	  }
	  problemStatus_=-1; // carry on as normal
	}
      } else {
	// iterate
	problemStatus_=-1;
      }
    }
  }
  if (type==0||type==1) {
    if (!type) {
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

  // restore weights (if saved) - also recompute infeasibility list
  if (tentativeStatus>-3) 
    dualRowPivot_->saveWeights(this,(type <2) ? 2 : 4);
  else
    dualRowPivot_->saveWeights(this,3);
  // unflag all variables (we may want to wait a bit?)
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot=pivotVariable_[iRow];
    clearFlagged(iPivot);
  }
  if (problemStatus_<0&&!changeMade_) {
    problemStatus_=4; // unknown
  }

}
/* While updateDualsInDual sees what effect is of flip
   this does actual flipping.
   If change >0.0 then value in array >0.0 => from lower to upper
*/
void 
ClpSimplexDual::flipBounds(CoinIndexedVector * rowArray,
		  CoinIndexedVector * columnArray,
		  double change)
{
  double * work;
  int number;
  int * which;
  
  int iSection;

  for (iSection=0;iSection<2;iSection++) {
    int i;
    double * solution = solutionRegion(iSection);
    double * lower = lowerRegion(iSection);
    double * upper = upperRegion(iSection);
    int addSequence;
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
#ifndef NDEBUG
      double value = work[iSequence]*change;
#endif
      work[iSequence]=0.0;
      Status status = getStatus(iSequence+addSequence);

      switch(status) {

      case basic:
      case isFree:
      case superBasic:
	break;
      case atUpperBound:
	assert (value<=0.0);
	// to lower bound
	setStatus(iSequence+addSequence,atLowerBound);
	solution[iSequence] = lower[iSequence];
	break;
      case atLowerBound:
	assert (value>=0.0);
	// to upper bound 
	setStatus(iSequence+addSequence,atUpperBound);
	solution[iSequence] = upper[iSequence];
	break;
      }
    }
  }
  rowArray->setNumElements(0);
  columnArray->setNumElements(0);
}
// Restores bound to original bound
void 
ClpSimplexDual::originalBound( int iSequence)
{
  if (iSequence>=numberColumns_) {
    // rows
    int iRow = iSequence-numberColumns_;
    rowLowerWork_[iRow]=rowLower_[iRow];
    rowUpperWork_[iRow]=rowUpper_[iRow];
    if (rowScale_) {
      if (rowLowerWork_[iRow]>-1.0e50)
	rowLowerWork_[iRow] *= rowScale_[iRow];
      if (rowUpperWork_[iRow]<1.0e50)
	rowUpperWork_[iRow] *= rowScale_[iRow];
    }
  } else {
    // columns
    columnLowerWork_[iSequence]=columnLower_[iSequence];
    columnUpperWork_[iSequence]=columnUpper_[iSequence];
    if (rowScale_) {
      double multiplier = 1.0/columnScale_[iSequence];
      if (columnLowerWork_[iSequence]>-1.0e50)
	columnLowerWork_[iSequence] *= multiplier;
      if (columnUpperWork_[iSequence]<1.0e50)
	columnUpperWork_[iSequence] *= multiplier;
    }
  }
  setFakeBound(iSequence,noFake);
}
/* As changeBounds but just changes new bounds for a single variable.
   Returns true if change */
bool 
ClpSimplexDual::changeBound( int iSequence)
{
  // old values 
  double oldLower=lower_[iSequence];
  double oldUpper=upper_[iSequence];
  double value=solution_[iSequence];
  bool modified=false;
  originalBound(iSequence);
  // original values
  double lowerValue=lower_[iSequence];
  double upperValue=upper_[iSequence];
  // back to altered values
  lower_[iSequence] = oldLower;
  upper_[iSequence] = oldUpper;
  if (value==oldLower) {
    if (upperValue > oldLower + dualBound_) {
      upper_[iSequence]=oldLower+dualBound_;
      setFakeBound(iSequence,upperFake);
      modified=true;
    }
  } else if (value==oldUpper) {
    if (lowerValue < oldUpper - dualBound_) {
      lower_[iSequence]=oldUpper-dualBound_;
      setFakeBound(iSequence,lowerFake);
      modified=true;
    }
  } else {
    assert(value==oldLower||value==oldUpper);
  }
  return modified;
}
// Perturbs problem
void 
ClpSimplexDual::perturb()
{
  if (perturbation_>100)
    return; //perturbed already
  int iRow,iColumn;
  // dual perturbation
  double perturbation=1.0e-20;
  // maximum fraction of cost to perturb
  double maximumFraction = 1.0e-4;
  if (perturbation_>=50) {
    perturbation = 1.0e-4;
    for (iRow=0;iRow<numberRows_;iRow++) {
      double value = fabs(rowObjectiveWork_[iRow]);
      perturbation = max(perturbation,value);
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) { 
      double value = 
	fabs(objectiveWork_[iColumn]);
      perturbation = max(perturbation,value);
    }
    perturbation *= 1.0e-8;
  } else if (perturbation_<100) {
    perturbation = pow(10.0,perturbation_);
    // user is in charge
    maximumFraction = 1.0e100;
  }
  // modify costs
  handler_->message(CLP_SIMPLEX_PERTURB,messages_)
    <<perturbation
    <<CoinMessageEol;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double value = perturbation;
    double currentValue = rowObjectiveWork_[iRow];
    value = min(value,maximumFraction*fabs(currentValue)+1.0e-6);
    if (rowLowerWork_[iRow]>-largeValue_) {
      if (fabs(rowLowerWork_[iRow])<fabs(rowUpperWork_[iRow])) 
	value *= CoinDrand48();
      else
	value *= -CoinDrand48();
    } else if (rowUpperWork_[iRow]<largeValue_) {
      value *= -CoinDrand48();
    } else {
      value=0.0;
    }
    rowObjectiveWork_[iRow] += value;
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double value = perturbation;
    double currentValue = objectiveWork_[iColumn];
    value = min(value,maximumFraction*fabs(currentValue)+1.0e-6);
    if (columnLowerWork_[iColumn]>-largeValue_) {
      if (fabs(columnLowerWork_[iColumn])<
	  fabs(columnUpperWork_[iColumn])) 
	value *= CoinDrand48();
      else
	value *= -CoinDrand48();
    } else if (columnUpperWork_[iColumn]<largeValue_) {
      value *= -CoinDrand48();
    } else {
      value=0.0;
    }
    objectiveWork_[iColumn] += value;
  }
  // say perturbed
  perturbation_=101;

}
/* For strong branching.  On input lower and upper are new bounds
   while on output they are change in objective function values 
   (>1.0e50 infeasible).
   Return code is 0 if nothing interesting, -1 if infeasible both
   ways and +1 if infeasible one way (check values to see which one(s))
*/
int ClpSimplexDual::strongBranching(int numberVariables,const int * variables,
				    double * newLower, double * newUpper,
				    bool stopOnFirstInfeasible,
				    bool alwaysFinish)
{
  int i;
  int returnCode=0;
  double saveObjectiveValue = objectiveValue_;
#if 1
  algorithm_ = -1;

  //scaling(false);

  // put in standard form (and make row copy)
  // create modifiable copies of model rim and do optional scaling
  createRim(7+8+16,true);

  // change newLower and newUpper if scaled

  // Do initial factorization
  // and set certain stuff
  // We can either set increasing rows so ...IsBasic gives pivot row
  // or we can just increment iBasic one by one
  // for now let ...iBasic give pivot row
  factorization_->increasingRows(2);
  // row activities have negative sign
  factorization_->slackValue(-1.0);
  factorization_->zeroTolerance(1.0e-13);
  // save if sparse factorization wanted
  int saveSparse = factorization_->sparseThreshold();

  int factorizationStatus = internalFactorize(0);
  if (factorizationStatus<0)
    return 1; // some error
  else if (factorizationStatus)
    handler_->message(CLP_SINGULARITIES,messages_)
      <<factorizationStatus
      <<CoinMessageEol;
  if (saveSparse) {
    // use default at present
    factorization_->sparseThreshold(0);
    factorization_->goSparse();
  }
  
  // save stuff
  ClpFactorization saveFactorization(*factorization_);
  // save basis and solution 
  double * saveSolution = new double[numberRows_+numberColumns_];
  memcpy(saveSolution,solution_,
	 (numberRows_+numberColumns_)*sizeof(double));
  unsigned char * saveStatus = 
    new unsigned char [numberRows_+numberColumns_];
  memcpy(saveStatus,status_,(numberColumns_+numberRows_)*sizeof(char));
  // save bounds as createRim makes clean copies
  double * saveLower = new double[numberRows_+numberColumns_];
  memcpy(saveLower,lower_,
	 (numberRows_+numberColumns_)*sizeof(double));
  double * saveUpper = new double[numberRows_+numberColumns_];
  memcpy(saveUpper,upper_,
	 (numberRows_+numberColumns_)*sizeof(double));
  double * saveObjective = new double[numberRows_+numberColumns_];
  memcpy(saveObjective,cost_,
	 (numberRows_+numberColumns_)*sizeof(double));
  int * savePivot = new int [numberRows_];
  memcpy(savePivot, pivotVariable_, numberRows_*sizeof(int));
  // need to save/restore weights.

  for (i=0;i<numberVariables;i++) {
    int iColumn = variables[i];
    double objectiveChange;
    double saveBound;
    
    // try down
    
    saveBound = columnUpper_[iColumn];
    // external view - in case really getting optimal 
    columnUpper_[iColumn] = newUpper[i];
    if (scalingFlag_<=0) 
      upper_[iColumn] = newUpper[i];
    else 
      upper_[iColumn] = newUpper[i]/columnScale_[iColumn]; // scale
    // Start of fast iterations
    int status = fastDual(alwaysFinish);

    // restore
    memcpy(solution_,saveSolution,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(status_,saveStatus,(numberColumns_+numberRows_)*sizeof(char));
    memcpy(lower_,saveLower,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(upper_,saveUpper,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(cost_,saveObjective,
	 (numberRows_+numberColumns_)*sizeof(double));
    columnUpper_[iColumn] = saveBound;
    memcpy(pivotVariable_, savePivot, numberRows_*sizeof(int));
    delete factorization_;
    factorization_ = new ClpFactorization(saveFactorization);

    if (status||(problemStatus_==0&&!isDualObjectiveLimitReached())) {
      objectiveChange = objectiveValue_-saveObjectiveValue;
    } else {
      objectiveChange = 1.0e100;
    }
    newUpper[i]=objectiveChange;
#ifdef CLP_DEBUG
    printf("down on %d costs %g\n",iColumn,objectiveChange);
#endif
	  
    // try up
    
    saveBound = columnLower_[iColumn];
    // external view - in case really getting optimal 
    columnLower_[iColumn] = newLower[i];
    if (scalingFlag_<=0) 
      lower_[iColumn] = newLower[i];
    else 
      lower_[iColumn] = newLower[i]/columnScale_[iColumn]; // scale
    // Start of fast iterations
    status = fastDual(alwaysFinish);

    // restore
    memcpy(solution_,saveSolution,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(status_,saveStatus,(numberColumns_+numberRows_)*sizeof(char));
    memcpy(lower_,saveLower,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(upper_,saveUpper,
	   (numberRows_+numberColumns_)*sizeof(double));
    memcpy(cost_,saveObjective,
	 (numberRows_+numberColumns_)*sizeof(double));
    columnLower_[iColumn] = saveBound;
    memcpy(pivotVariable_, savePivot, numberRows_*sizeof(int));
    delete factorization_;
    factorization_ = new ClpFactorization(saveFactorization);

    if (status||(problemStatus_==0&&!isDualObjectiveLimitReached())) {
      objectiveChange = objectiveValue_-saveObjectiveValue;
    } else {
      objectiveChange = 1.0e100;
    }
    newLower[i]=objectiveChange;
#ifdef CLP_DEBUG
    printf("up on %d costs %g\n",iColumn,objectiveChange);
#endif
	  
    /* Possibilities are:
       Both sides feasible - store
       Neither side feasible - set objective high and exit
       One side feasible - change bounds and resolve
    */
    if (newUpper[i]<1.0e100) {
      if(newLower[i]<1.0e100) {
	// feasible - no action
      } else {
	// up feasible, down infeasible
	returnCode=1;
	if (stopOnFirstInfeasible)
	  break;
      }
    } else {
      if(newLower[i]<1.0e100) {
	// down feasible, up infeasible
	returnCode=1;
	if (stopOnFirstInfeasible)
	  break;
      } else {
	// neither side feasible
	returnCode=-1;
	break;
      }
    }
  }
  delete [] saveSolution;
  delete [] saveLower;
  delete [] saveUpper;
  delete [] saveObjective;
  delete [] saveStatus;
  delete [] savePivot;

  factorization_->sparseThreshold(saveSparse);
  // Get rid of some arrays and empty factorization
  deleteRim();
#else
  // save basis and solution 
  double * rowSolution = new double[numberRows_];
  memcpy(rowSolution,rowActivity_,numberRows_*sizeof(double));
  double * columnSolution = new double[numberColumns_];
  memcpy(columnSolution,columnActivity_,numberColumns_*sizeof(double));
  unsigned char * saveStatus = 
    new unsigned char [numberRows_+numberColumns_];
  if (!status_) {
    // odd but anyway setup all slack basis
    status_ = new unsigned char [numberColumns_+numberRows_];
    memset(status_,0,(numberColumns_+numberRows_)*sizeof(char));
  }
  memcpy(saveStatus,status_,(numberColumns_+numberRows_)*sizeof(char));
  for (i=0;i<numberVariables;i++) {
    int iColumn = variables[i];
    double objectiveChange;
    
    // try down
    
    double saveUpper = columnUpper_[iColumn];
    columnUpper_[iColumn] = newUpper[i];
    dual();

    // restore
    columnUpper_[iColumn] = saveUpper;
    memcpy(rowActivity_,rowSolution,numberRows_*sizeof(double));
    memcpy(columnActivity_,columnSolution,numberColumns_*sizeof(double));
    memcpy(status_,saveStatus,(numberColumns_+numberRows_)*sizeof(char));

    if (problemStatus_==0&&!isDualObjectiveLimitReached()) {
      objectiveChange = objectiveValue_-saveObjectiveValue;
    } else {
      objectiveChange = 1.0e100;
    }
    newUpper[i]=objectiveChange;
#ifdef CLP_DEBUG
    printf("down on %d costs %g\n",iColumn,objectiveChange);
#endif
	  
    // try up
    
    double saveLower = columnLower_[iColumn];
    columnLower_[iColumn] = newLower[i];
    dual();

    // restore
    columnLower_[iColumn] = saveLower;
    memcpy(rowActivity_,rowSolution,numberRows_*sizeof(double));
    memcpy(columnActivity_,columnSolution,numberColumns_*sizeof(double));
    memcpy(status_,saveStatus,(numberColumns_+numberRows_)*sizeof(char));

    if (problemStatus_==0&&!isDualObjectiveLimitReached()) {
      objectiveChange = objectiveValue_-saveObjectiveValue;
    } else {
      objectiveChange = 1.0e100;
    }
    newLower[i]=objectiveChange;
#ifdef CLP_DEBUG
    printf("up on %d costs %g\n",iColumn,objectiveChange);
#endif
	  
    /* Possibilities are:
       Both sides feasible - store
       Neither side feasible - set objective high and exit
       One side feasible - change bounds and resolve
    */
    if (newUpper[i]<1.0e100) {
      if(newLower[i]<1.0e100) {
	// feasible - no action
      } else {
	// up feasible, down infeasible
	returnCode=1;
	if (stopOnFirstInfeasible)
	  break;
      }
    } else {
      if(newLower[i]<1.0e100) {
	// down feasible, up infeasible
	returnCode=1;
	if (stopOnFirstInfeasible)
	  break;
      } else {
	// neither side feasible
	returnCode=-1;
	break;
      }
    }
  }
  delete [] rowSolution;
  delete [] columnSolution;
  delete [] saveStatus;
#endif
  objectiveValue_ = saveObjectiveValue;
  return returnCode;
}
// treat no pivot as finished (unless interesting)
int ClpSimplexDual::fastDual(bool alwaysFinish)
{
  algorithm_ = -1;
  dualTolerance_=dblParam_[ClpDualTolerance];
  primalTolerance_=dblParam_[ClpPrimalTolerance];

  // save dual bound
  double saveDualBound = dualBound_;

  double objectiveChange;
  // for dual we will change bounds using dualBound_
  // for this we need clean basis so it is after factorize
  gutsOfSolution(rowActivityWork_,columnActivityWork_);
  numberFake_ =0; // Number of variables at fake bounds
  changeBounds(true,NULL,objectiveChange);

  problemStatus_ = -1;
  numberIterations_=0;

  int lastCleaned=0; // last time objective or bounds cleaned up

  // number of times we have declared optimality
  numberTimesOptimal_=0;

  // Progress indicator
  ClpSimplexProgress progress(this);

  // This says whether to restore things etc
  int factorType=0;
  /*
    Status of problem:
    0 - optimal
    1 - infeasible
    2 - unbounded
    -1 - iterating
    -2 - factorization wanted
    -3 - redo checking without factorization
    -4 - looks infeasible

    BUT also from whileIterating return code is:

   -1 iterations etc
   -2 inaccuracy 
   -3 slight inaccuracy (and done iterations)
   +0 looks optimal (might be unbounded - but we will investigate)
   +1 looks infeasible
   +3 max iterations 

  */

  int returnCode = 0;

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
    // may factorize, checks if problem finished
    // should be able to speed this up on first time
    statusOfProblemInDual(lastCleaned,factorType,progress);

    // Say good factorization
    factorType=1;

    // Do iterations
    if (problemStatus_<0) {
      returnCode = whileIterating();
      if (!alwaysFinish&&returnCode<1) {
	double limit = 0.0;
	getDblParam(ClpDualObjectiveLimit, limit);
	if(fabs(limit)>1.0e30||objectiveValue()*optimizationDirection_<
	   optimizationDirection_*limit|| 
	   numberAtFakeBound()) {
	  // can't say anything interesting - might as well return
#ifdef CLP_DEBUG
	  printf("returning from fastDual after %d iterations with code %d\n",
		 numberIterations_,returnCode);
#endif
	  returnCode=1;
	  break;
	}
      }
      returnCode=0;
    }
  }

  assert(!numberFake_||returnCode||problemStatus_); // all bounds should be okay
  dualBound_ = saveDualBound;
  return returnCode;
}
/* Checks number of variables at fake bounds.  This is used by fastDual
   so can exit gracefully before end */
int 
ClpSimplexDual::numberAtFakeBound()
{
  int iSequence;
  int numberFake=0;
  
  for (iSequence=0;iSequence<numberRows_+numberColumns_;iSequence++) {
    FakeBound bound = getFakeBound(iSequence);
    switch(getStatus(iSequence)) {

    case basic:
      break;
    case isFree:
    case superBasic:
      break;
    case atUpperBound:
      if (bound==upperFake||bound==bothFake)
	numberFake++;
      break;
    case atLowerBound:
      if (bound==lowerFake||bound==bothFake)
	numberFake++;
      break;
    }
  }
  numberFake_ = numberFake;
  return numberFake;
}
