// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Implements crude primal dual predictor corrector algorithm

 */
//#define SOME_DEBUG

#include "CoinPragma.hpp"
#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpPredictorCorrector.hpp"
#include "CoinPackedMatrix.hpp"
#include "ClpMessage.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpHelperFunctions.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <cstdio>
#include <iostream>

static double eScale=1.0e57;
static double eBaseCaution=1.0e-12;
static double eBase=1.0e-12;
static double eRatio=1.0e40;
static double eRatioCaution=1.0e25;
static double eDiagonal=1.0e25;
static double eDiagonalCaution=1.0e18;
static double eExtra=1.0e-12;

// main function

int ClpPredictorCorrector::solve ( )
{
  problemStatus_=-1;
  algorithm_=1;
  //create all regions
  if (!createWorkingData()) {
    problemStatus_=4;
    return 2;
  }
  ClpMatrixBase * saveMatrix = NULL;
  // If scaled then really scale matrix
  if (scalingFlag_>0&&rowScale_) {
    saveMatrix = matrix_;
    matrix_ = matrix_->scaledColumnCopy(this);
  }
  //initializeFeasible(); - this just set fixed flag
  smallestInfeasibility_=COIN_DBL_MAX;
  int i;
  for (i=0;i<LENGTH_HISTORY;i++) 
    historyInfeasibility_[i]=COIN_DBL_MAX;

  //bool firstTime=true;
  //firstFactorization(true);
  if (cholesky_->order(this)) {
    printf("Not enough memory\n");
    problemStatus_=4;
    //delete all temporary regions
    deleteWorkingData();
    if (saveMatrix) {
      // restore normal copy
      delete matrix_;
      matrix_ = saveMatrix;
    }
    return -1;
  }
  mu_=1.0e10;
  diagonalScaleFactor_=1.0;
  //set iterations
  numberIterations_=-1;
  int numberTotal = numberRows_+numberColumns_;
  //initialize solution here
  if(createSolution()<0) {
    printf("Not enough memory\n");
    problemStatus_=4;
    //delete all temporary regions
    deleteWorkingData();
    if (saveMatrix) {
      // restore normal copy
      delete matrix_;
      matrix_ = saveMatrix;
    }
    return -1;
  }
  // Could try centering steps without any original step i.e. just center
  //firstFactorization(false);
  CoinZeroN(dual_,numberRows_);
  multiplyAdd(solution_+numberColumns_,numberRows_,-1.0,errorRegion_,0.0);
  matrix_->times(1.0,solution_,errorRegion_);
  maximumRHSError_=maximumAbsElement(errorRegion_,numberRows_);
  maximumBoundInfeasibility_=maximumRHSError_;
  //double maximumDualError_=COIN_DBL_MAX;
  //initialize
  actualDualStep_=0.0;
  actualPrimalStep_=0.0;
  gonePrimalFeasible_=false;
  goneDualFeasible_=false;
  //bool hadGoodSolution=false;
  diagonalNorm_=solutionNorm_;
  mu_=solutionNorm_;
  int numberFixed=updateSolution(-COIN_DBL_MAX);
  int numberFixedTotal=numberFixed;
  //int numberRows_DroppedBefore=0;
  //double extra=eExtra;
  //double maximumPerturbation=COIN_DBL_MAX;
  //constants for infeas interior point
  const double beta2 = 0.99995;
  const double tau   = 0.00002;
  double lastComplementarityGap=COIN_DBL_MAX;
  // use to see if to take affine
  double checkGap = COIN_DBL_MAX;
  int lastGoodIteration=0;
  double bestObjectiveGap=COIN_DBL_MAX;
  int saveIteration=-1;
  bool sloppyOptimal=false;
  double * savePi=NULL;
  double * savePrimal=NULL;
  // Extra regions for centering
  double * saveX = new double[numberTotal];
  double * saveY = new double[numberRows_];
  double * saveZ = new double[numberTotal];
  double * saveT = new double[numberTotal];
  // Save smallest mu used in primal dual moves
  double smallestPrimalDualMu=COIN_DBL_MAX;
  while (problemStatus_<0) {
    //#define FULL_DEBUG
#ifdef FULL_DEBUG
    {
      int i;
      printf("row    pi          artvec       rhsfx\n");
      for (i=0;i<numberRows_;i++) {
	printf("%d %g %g %g\n",i,dual_[i],errorRegion_[i],rhsFixRegion_[i]);
      }
      printf(" col  dsol  ddiag  dwvec  dzvec dbdslu dbdsll\n");
      for (i=0;i<numberColumns_+numberRows_;i++) {
	printf(" %d %g %g %g %g %g %g\n",i,solution_[i],diagonal_[i],tVec_[i],
	       zVec_[i],upperSlack_[i],lowerSlack_[i]);
      }
    }
#endif
    complementarityGap_=complementarityGap(numberComplementarityPairs_,
					   numberComplementarityItems_,0);
      handler_->message(CLP_BARRIER_ITERATION,messages_)
    <<numberIterations_
    <<primalObjective_*optimizationDirection_- dblParam_[ClpObjOffset]
    << dualObjective_*optimizationDirection_- dblParam_[ClpObjOffset]
    <<complementarityGap_
    <<numberFixedTotal
    <<cholesky_->rank()
    <<CoinMessageEol;
    double goodGapChange;
    if (!sloppyOptimal) {
      goodGapChange=0.93;
    } else {
      goodGapChange=0.7;
    } 
    double gapO;
    double lastGood=bestObjectiveGap;
    if (gonePrimalFeasible_&&goneDualFeasible_) {
      double largestObjective;
      if (fabs(primalObjective_)>fabs(dualObjective_)) {
        largestObjective = fabs(primalObjective_);
      } else {
        largestObjective = fabs(dualObjective_);
      } 
      if (largestObjective<1.0) {
        largestObjective=1.0;
      } 
      gapO=fabs(primalObjective_-dualObjective_)/largestObjective;
      handler_->message(CLP_BARRIER_OBJECTIVE_GAP,messages_)
	<<gapO
	<<CoinMessageEol;
      //start saving best
      if (gapO<bestObjectiveGap) {
        saveIteration=numberIterations_;
        bestObjectiveGap=gapO;
        if (!savePi) {
          savePi=new double[numberRows_];
          savePrimal = new double [numberTotal];
        } 
        CoinMemcpyN(dual_,numberRows_,savePi);
	CoinMemcpyN(solution_,numberTotal,savePrimal);
      } else if(gapO>2.0*bestObjectiveGap) {
        //maybe be more sophisticated e.g. re-initialize having
        //fixed variables and dropped rows
        //std::cout <<" gap increasing "<<std::endl;
      } 
      //std::cout <<"could stop"<<std::endl;
      //gapO=0.0;
      if (fabs(primalObjective_-dualObjective_)<dualTolerance()) {
        gapO=0.0;
      } 
    } else {
      gapO=COIN_DBL_MAX;
      if (saveIteration>=0) {
	handler_->message(CLP_BARRIER_GONE_INFEASIBLE,messages_)
	  <<CoinMessageEol;
	if (sloppyOptimal) {
	  // vaguely optimal
	  double scaledRHSError=maximumRHSError_/solutionNorm_;
	  if (maximumBoundInfeasibility_>1.0e-2||
	      scaledRHSError>1.0e-2||
	      maximumDualError_>objectiveNorm_*1.0e-2) {
	    handler_->message(CLP_BARRIER_EXIT2,messages_)
	      <<saveIteration
	      <<CoinMessageEol;
	    break;
	  }
	} else {
	  // not close to optimal but check if getting bad
	  double scaledRHSError=maximumRHSError_/solutionNorm_;
	  if (maximumBoundInfeasibility_>1.0e-1||
	      scaledRHSError>1.0e-1||
	      maximumDualError_>objectiveNorm_*1.0e-1) {
	    handler_->message(CLP_BARRIER_EXIT2,messages_)
	      <<saveIteration
	      <<CoinMessageEol;
	    break;
	  }
	}
      } 
    } 
    if ((gapO<1.0e-6||(gapO<1.0e-4&&complementarityGap_<0.1))&&!sloppyOptimal) {
      sloppyOptimal=true;
      handler_->message(CLP_BARRIER_CLOSE_TO_OPTIMAL,messages_)
	<<numberIterations_<<complementarityGap_
	<<CoinMessageEol;
    } 
    //if (complementarityGap_>=0.98*lastComplementarityGap) {
    //tryJustPredictor=true;
    //printf("trying just predictor\n");
    //}
    if (complementarityGap_>=1.05*lastComplementarityGap) {
      handler_->message(CLP_BARRIER_COMPLEMENTARITY,messages_)
	<<complementarityGap_<<"increasing"
	<<CoinMessageEol;
      if (saveIteration>=0&&sloppyOptimal) {
	handler_->message(CLP_BARRIER_EXIT2,messages_)
	  <<saveIteration
	  <<CoinMessageEol;
        break;
      } else if (numberIterations_-lastGoodIteration>=5&&
		 complementarityGap_<1.0e-6) {
	break; // not doing very well - give up
      } 
    } else if (complementarityGap_<goodGapChange*lastComplementarityGap) {
      lastGoodIteration=numberIterations_;
      lastComplementarityGap=complementarityGap_;
    } else if (numberIterations_-lastGoodIteration>=5&&
	       complementarityGap_<1.0e-3) {
      handler_->message(CLP_BARRIER_COMPLEMENTARITY,messages_)
	<<complementarityGap_<<"not decreasing"
	<<CoinMessageEol;
      if (gapO>0.75*lastGood) {
        break;
      } 
    } else if (numberIterations_-lastGoodIteration>=2&&
	       complementarityGap_<1.0e-6) {
      handler_->message(CLP_BARRIER_COMPLEMENTARITY,messages_)
	<<complementarityGap_<<"not decreasing"
	<<CoinMessageEol;
      break;
    } 
    if (numberIterations_>maximumBarrierIterations_) {
      handler_->message(CLP_BARRIER_STOPPING,messages_)
	<<CoinMessageEol;
      break;
    } 
    if (gapO<targetGap_) {
      problemStatus_=0;
      handler_->message(CLP_BARRIER_EXIT,messages_)
	<<" "
	<<CoinMessageEol;
        break;//finished
    } 
    if (complementarityGap_<1.0e-12) {
      problemStatus_=0;
      handler_->message(CLP_BARRIER_EXIT,messages_)
        <<"- small complementarity gap"
	<<CoinMessageEol;
        break;//finished
    } 
    if (complementarityGap_<1.0e-10&&gapO<1.0e-10) {
      problemStatus_=0;
      handler_->message(CLP_BARRIER_EXIT,messages_)
        <<"- objective gap and complementarity gap both small"
	<<CoinMessageEol;
        break;//finished
    } 
    if (gapO<1.0e-9) {
      double value=gapO*complementarityGap_;
      value*=actualPrimalStep_;
      value*=actualDualStep_;
      //std::cout<<value<<std::endl;
      if (value<1.0e-17&&numberIterations_>lastGoodIteration) {
	problemStatus_=0;
	handler_->message(CLP_BARRIER_EXIT,messages_)
	  <<"- objective gap and complementarity gap both smallish and small steps"
	  <<CoinMessageEol;
        break;//finished
      } 
    } 
    double nextGap=COIN_DBL_MAX;
    int nextNumber=0;
    int nextNumberItems=0;
    worstDirectionAccuracy_=0.0;
    int newDropped=0;
    //Predictor step
    //prepare for cholesky.  Set up scaled diagonal in deltaX
    //  ** for efficiency may be better if scale factor known before
    double norm2=0.0;
    double maximumValue;
    getNorms(diagonal_,numberTotal,maximumValue,norm2);
    diagonalNorm_ = sqrt(norm2/numberComplementarityPairs_);
    diagonalScaleFactor_=1.0;
    double maximumAllowable=eScale;
    //scale so largest is less than allowable ? could do better
    double factor=0.5;
    while (maximumValue>maximumAllowable) {
      diagonalScaleFactor_*=factor;
      maximumValue*=factor;
    } /* endwhile */
    if (diagonalScaleFactor_!=1.0) {
      handler_->message(CLP_BARRIER_SCALING,messages_)
	<<"diagonal"<<diagonalScaleFactor_
	<<CoinMessageEol;
      diagonalNorm_*=diagonalScaleFactor_;
    } 
    multiplyAdd(NULL,numberTotal,0.0,diagonal_,
		diagonalScaleFactor_);
    int * rowsDroppedThisTime = new int [numberRows_];
    newDropped=cholesky_->factorize(diagonal_,rowsDroppedThisTime);
    if (newDropped) {
      if (newDropped==-1) {
	printf("Out of memory\n");
	problemStatus_=4;
	//delete all temporary regions
	deleteWorkingData();
	if (saveMatrix) {
	  // restore normal copy
	  delete matrix_;
	  matrix_ = saveMatrix;
	}
	return -1;
      } else {
	//int newDropped2=cholesky_->factorize(diagonal_,rowsDroppedThisTime);
	//assert(!newDropped2);
	if (newDropped<0&&0) {
	  //replace dropped
	  newDropped=-newDropped;
	  //off 1 to allow for reset all
	  newDropped--;
	  //set all bits false
	  cholesky_->resetRowsDropped();
	} 
      }
    } 
    delete [] rowsDroppedThisTime;
    if (cholesky_->status()) {
      std::cout << "bad cholesky?" <<std::endl;
      abort();
    }
    int phase=0; // predictor, corrector , primal dual
    double directionAccuracy=0.0;
    bool doCorrector=true;
    bool goodMove=true;
    //set up for affine direction
    setupForSolve(phase);
    directionAccuracy=findDirectionVector(phase);
    if (directionAccuracy>worstDirectionAccuracy_) {
      worstDirectionAccuracy_=directionAccuracy;
    } 
    findStepLength(phase);
    nextGap=complementarityGap(nextNumber,nextNumberItems,1);
    double affineGap=nextGap;
    int bestPhase=0;
    double bestNextGap=nextGap;
    // ?
    bestNextGap=max(nextGap,0.8*complementarityGap_);
    bestNextGap=max(nextGap,0.99*complementarityGap_);
    if (complementarityGap_>1.0e-4*numberComplementarityPairs_) {
      //std::cout <<"predicted duality gap "<<nextGap<<std::endl;
      double part1=nextGap/numberComplementarityPairs_;
      part1=nextGap/numberComplementarityItems_;
      double part2=nextGap/complementarityGap_;
      mu_=part1*part2*part2;
#if 0
      double papermu =complementarityGap_/numberComplementarityPairs_;
      double affmu = nextGap/nextNumber;
      double sigma = pow(affmu/papermu,3);
      printf("mu %g, papermu %g, affmu %g, sigma %g sigmamu %g\n",
	     mu_,papermu,affmu,sigma,sigma*papermu);
#endif	
      //printf("paper mu %g\n",(nextGap*nextGap*nextGap)/(complementarityGap_*complementarityGap_*
      //					    (double) numberComplementarityPairs_));
    } else {
      double phi;
      if (numberComplementarityPairs_<=500) {
	phi=pow((double) numberComplementarityPairs_,2.0);
      } else {
	phi=pow((double) numberComplementarityPairs_,1.5);
	if (phi<500.0*500.0) {
	  phi=500.0*500.0;
	} 
      }
      mu_=complementarityGap_/phi;
    } 
    //save information
    double product=affineProduct();
    //#define ALWAYS
#ifndef ALWAYS
#if 0
    //can we do corrector step?
    double xx= complementarityGap_*(beta2-tau) +product;
    if (xx>0.0) {
      double saveMu = mu_;
      double mu2=numberComplementarityPairs_;
      mu2=xx/mu2;
      if (mu2>mu_) {
	//std::cout<<" could increase to "<<mu2<<std::endl;
	//was mu2=mu2*0.25;
	mu2=mu2*0.99;
	if (mu2<mu_) {
	  mu_=mu2;
	  //std::cout<<" changing mu to "<<mu_<<std::endl;
	} else {
	  //std::cout<<std::endl;
	} 
      } else {
	//std::cout<<" should decrease to "<<mu2<<std::endl;
	mu_=0.5*mu2;
	//std::cout<<" changing mu to "<<mu_<<std::endl;
      } 
      handler_->message(CLP_BARRIER_MU,messages_)
	<<saveMu<<mu_
	<<CoinMessageEol;
    } else {
      //std::cout<<" bad by any standards"<<std::endl;
    }
#endif 
    if (complementarityGap_*(beta2-tau)+product-mu_*numberComplementarityPairs_<0.0) {
#ifdef SOME_DEBUG
      printf("failed 1 product %.18g mu %.18g - %.18g < 0.0, nextGap %.18g\n",product,mu_,
	     complementarityGap_*(beta2-tau)+product-mu_*numberComplementarityPairs_,
	     nextGap);
#endif
      doCorrector=false;
      if (nextGap>0.9*complementarityGap_||1) {
	goodMove=false;
	bestNextGap=COIN_DBL_MAX;
      }
      //double floatNumber = 2.0*numberComplementarityPairs_;
      //floatNumber = 1.0*numberComplementarityItems_;
      //mu_=nextGap/floatNumber;
      handler_->message(CLP_BARRIER_INFO,messages_)
	<<"no corrector step"
	<<CoinMessageEol;
    } else {
      phase=1;
    }
#else
    phase=1;
#endif
    if (goodMove&&doCorrector) {
      //set up for next step
      setupForSolve(phase);
      double directionAccuracy2=findDirectionVector(phase);
      if (directionAccuracy2>worstDirectionAccuracy_) {
	worstDirectionAccuracy_=directionAccuracy2;
      } 
      double testValue=1.0e2*directionAccuracy;
      if (1.0e2*projectionTolerance_>testValue) {
	testValue=1.0e2*projectionTolerance_;
      } 
      if (primalTolerance()>testValue) {
	testValue=primalTolerance();
      } 
      if (maximumRHSError_>testValue) {
	testValue=maximumRHSError_;
      } 
      if (directionAccuracy2>testValue&&numberIterations_>=-77) {
	goodMove=false;
#ifdef SOME_DEBUG
	printf("accuracy %g phase 1 failed, test value %g\n",
	       directionAccuracy2,testValue);
#endif
	bestNextGap=COIN_DBL_MAX;
      } 
      if (goodMove) {
	findStepLength(phase);
	nextGap = complementarityGap(nextNumber,nextNumberItems,1);
#ifndef ALWAYS
	if (numberIterations_>=-77) {
	  goodMove=checkGoodMove(true,bestNextGap);
	  if (!goodMove) {
#ifdef SOME_DEBUG
	    printf("checkGoodMove failed\n");
#endif
	    if ((affineGap<0.5*complementarityGap_&&complementarityGap_>0.9*checkGap)
	      ||(affineGap<0.95*complementarityGap_&&
		 complementarityGap_>0.99*checkGap)) {
	      // Back to affine
	      phase=0;
  	      // Try primal dual step instead - but with small mu
	      phase=2;
	      double floatNumber;
	      floatNumber = 2.0*numberComplementarityPairs_;
	      mu_=complementarityGap_/floatNumber;
	      double mu1=mu_;
	      double phi;
	      if (numberComplementarityPairs_<=500) {
		phi=pow((double) numberComplementarityPairs_,2.0);
	      } else {
		phi=pow((double) numberComplementarityPairs_,1.5);
		if (phi<500.0*500.0) {
		  phi=500.0*500.0;
		} 
	      }
	      mu_=complementarityGap_/phi;
	      //printf("pd mu %g, alternate %g, smallest %g\n",
	      //     mu_,mu1,smallestPrimalDualMu);
	      mu_ = sqrt(mu_*mu1);
	      mu_=mu1*0.8;
              if (numberIterations_>100) 
                mu_ *=0.1; // so close to affine
	      setupForSolve(phase);
	      directionAccuracy=findDirectionVector(phase);
	      findStepLength(phase);
	      nextGap=complementarityGap(nextNumber,nextNumberItems,1);
	      goodMove=true;
	    }
	  }
	} else {
	  goodMove=true;
	} 
#endif
      }
    }
    //bestPhase=-1;
    //goodMove=false;
    if (!goodMove) {
      // Just primal dual step
      double floatNumber;
      floatNumber = 2.5*numberComplementarityPairs_;
      //floatNumber = 1.5*numberComplementarityItems_;
      double saveMu=mu_; // use one from predictor corrector
      mu_=complementarityGap_/floatNumber;
      double mu1=mu_;
      double phi;
      if (numberComplementarityPairs_<=500) {
	phi=pow((double) numberComplementarityPairs_,2.0);
      } else {
	phi=pow((double) numberComplementarityPairs_,1.5);
	if (phi<500.0*500.0) {
	  phi=500.0*500.0;
	} 
      }
      mu_=complementarityGap_/phi;
      //printf("pd mu %g, alternate %g, smallest %g\n",
      //     mu_,mu1,smallestPrimalDualMu);
      mu_ = sqrt(mu_*mu1);
      mu_=mu1;
      if ((numberIterations_&1)==0||numberIterations_<10)
	mu_=saveMu;
      //mu_=min(smallestPrimalDualMu*0.95,mu_);
      smallestPrimalDualMu = mu_;
      //set up for next step
      setupForSolve(2);
      findDirectionVector(2);
      findStepLength(2);
      // just for debug
      nextGap=complementarityGap(nextNumber,nextNumberItems,2);
      if (nextGap>0.99*complementarityGap_&&bestPhase==0&&affineGap<nextGap) {
	// Back to affine
	phase=0;
	// no
	phase=2;
	mu_ *= 0.5;
	setupForSolve(phase);
	directionAccuracy=findDirectionVector(phase);
	findStepLength(phase);
	nextGap=complementarityGap(nextNumber,nextNumberItems,1);
      }
    }
    if (numberIterations_==0)
      smallestPrimalDualMu=mu_;
    if (!goodMove)
      mu_=nextGap / ((double) 1.1*nextNumber);
    goodMove=true;
    //goodMove=false; //TEMP
    // Do centering steps
    int numberTries=0;
    double nextCenterGap=0.0;
    int numberGoodTries=0;
    double originalDualStep=actualDualStep_;
    double originalPrimalStep=actualPrimalStep_;
    if (actualDualStep_>0.9&&actualPrimalStep_>0.9)
      goodMove=false; // don't bother
    while (goodMove&&numberTries<5) {
      goodMove=false;
      numberTries++;
      memcpy(saveX,deltaX_,numberTotal*sizeof(double));
      memcpy(saveY,deltaY_,numberRows_*sizeof(double));
      memcpy(saveZ,deltaZ_,numberTotal*sizeof(double));
      memcpy(saveT,deltaT_,numberTotal*sizeof(double));
      double savePrimalStep = actualPrimalStep_;
      double saveDualStep = actualDualStep_;
      double saveMu = mu_;
      setupForSolve(3);
      findDirectionVector(3);
      findStepLength(3);
      double xGap = complementarityGap(nextNumber,nextNumberItems,3);
      // If one small then that's the one that counts
      double checkDual=saveDualStep;
      double checkPrimal=savePrimalStep;
      if (checkDual>5.0*checkPrimal) {
	checkDual=2.0*checkPrimal;
      } else if (checkPrimal>5.0*checkDual) {
	checkPrimal=2.0*checkDual;
      }
      if (actualPrimalStep_<=checkPrimal||
	  actualDualStep_<=checkDual||
	  (xGap>nextGap&&xGap>0.9*complementarityGap_)) {
	//if (actualPrimalStep_<=checkPrimal||
	//actualDualStep_<=checkDual) {
#ifdef SOME_DEBUG
	printf("PP rejected gap %.18g, steps %.18g %.18g, 2 gap %.18g, steps %.18g %.18g\n",xGap,
	       actualPrimalStep_,actualDualStep_,nextGap,savePrimalStep,saveDualStep);
#endif
	mu_=saveMu;
	actualPrimalStep_ = savePrimalStep;
	actualDualStep_ = saveDualStep;
	memcpy(deltaX_,saveX,numberTotal*sizeof(double));
	memcpy(deltaY_,saveY,numberRows_*sizeof(double));
	memcpy(deltaZ_,saveZ,numberTotal*sizeof(double));
	memcpy(deltaT_,saveT,numberTotal*sizeof(double));
      } else {
#ifdef SOME_DEBUG
	printf("PPphase 3 gap %.18g, steps %.18g %.18g, 2 gap %.18g, steps %.18g %.18g\n",xGap,
	       actualPrimalStep_,actualDualStep_,nextGap,savePrimalStep,saveDualStep);
#endif
	numberGoodTries++;
	nextCenterGap=xGap;
	// See if big enough change
	if (actualPrimalStep_<1.01*checkPrimal||
	    actualDualStep_<1.01*checkDual) {
	  // stop now
	} else {
	  // carry on
	  goodMove=true;
	}
      }
    }
    if (numberGoodTries&&handler_->logLevel()>1) {
      printf("%d centering steps moved from (gap %.18g, dual %.18g, primal %.18g) to (gap %.18g, dual %.18g, primal %.18g)\n",
	     numberGoodTries,nextGap,originalDualStep,originalPrimalStep,
	     nextCenterGap, actualDualStep_,actualPrimalStep_);
    }
    // save last gap
    checkGap = complementarityGap_;
    numberFixed=updateSolution(nextGap);
    numberFixedTotal+=numberFixed;
  } /* endwhile */
  delete [] saveX;
  delete [] saveY;
  delete [] saveZ;
  delete [] saveT;
  if (savePi) {
    //std::cout<<"Restoring from iteration "<<saveIteration<<std::endl;
    CoinMemcpyN(savePi,numberRows_,dual_);
    CoinMemcpyN(savePrimal,numberTotal,solution_);
    delete [] savePi;
    delete [] savePrimal;
  } 
  //recompute slacks
  // Split out solution
  CoinZeroN(rowActivity_,numberRows_);
  CoinMemcpyN(solution_,numberColumns_,columnActivity_);
  matrix_->times(1.0,columnActivity_,rowActivity_);
  //unscale objective
  multiplyAdd(NULL,numberTotal,0.0,cost_,scaleFactor_);
  multiplyAdd(NULL,numberRows_,0,dual_,scaleFactor_);
  CoinMemcpyN(cost_,numberColumns_,reducedCost_);
  matrix_->transposeTimes(-1.0,dual_,reducedCost_);
  CoinMemcpyN(reducedCost_,numberColumns_,dj_);
  checkSolution();
  handler_->message(CLP_BARRIER_END,messages_)
    <<sumPrimalInfeasibilities_
    <<sumDualInfeasibilities_
    <<complementarityGap_
    <<objectiveValue()
    <<CoinMessageEol;
  //#ifdef SOME_DEBUG
  if (handler_->logLevel()>1) 
    printf("ENDRUN status %d after %d iterations\n",problemStatus_,numberIterations_);
  //#endif
  //std::cout<<"Absolute primal infeasibility at end "<<sumPrimalInfeasibilities_<<std::endl;
  //std::cout<<"Absolute dual infeasibility at end "<<sumDualInfeasibilities_<<std::endl;
  //std::cout<<"Absolute complementarity at end "<<complementarityGap_<<std::endl;
  //std::cout<<"Primal objective "<<objectiveValue()<<std::endl;
  //std::cout<<"maximum complementarity "<<worstComplementarity_<<std::endl;
  //delete all temporary regions
  deleteWorkingData();
  if (saveMatrix) {
    // restore normal copy
    delete matrix_;
    matrix_ = saveMatrix;
  }
  return problemStatus_;
}
// findStepLength.
//phase  - 0 predictor
//         1 corrector
//         2 primal dual
double ClpPredictorCorrector::findStepLength( int phase)
{
  double directionNorm=0.0;
  double maximumPrimalStep=COIN_DBL_MAX;
  double maximumDualStep=COIN_DBL_MAX;
  int numberTotal = numberRows_+numberColumns_;
  double tolerance = 1.0e-12;
  int chosenPrimalSequence=-1;
  int chosenDualSequence=-1;
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  //direction vector in deltaX
  double * deltaX = deltaX_;
  int iColumn;
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double directionElement=deltaX[iColumn];
      if (directionNorm<fabs(directionElement)) {
	directionNorm=fabs(directionElement);
      } 
      if (lowerBound(iColumn)) {
	double delta = - deltaSL_[iColumn];
	double z1 = deltaZ_[iColumn];
	if (lowerSlack[iColumn]<maximumPrimalStep*delta) {
	  maximumPrimalStep=lowerSlack[iColumn]/delta;
	  chosenPrimalSequence=iColumn;
	} 
	if (zVec[iColumn]>tolerance) {
	  if (zVec[iColumn]<-z1*maximumDualStep) {
	    maximumDualStep=-zVec[iColumn]/z1;
	    chosenDualSequence=iColumn;
	  } 
	} 
      }
      if (upperBound(iColumn)) {
	double delta = - deltaSU_[iColumn];;
	double t1 = deltaT_[iColumn];
	if (upperSlack[iColumn]<maximumPrimalStep*delta) {
	  maximumPrimalStep=upperSlack[iColumn]/delta;
	  chosenPrimalSequence=iColumn;
	} 
	if (tVec[iColumn]>tolerance) {
	  if (tVec[iColumn]<-t1*maximumDualStep) {
	    maximumDualStep=-tVec[iColumn]/t1;
	    chosenDualSequence=iColumn;
	  } 
	} 
      } 
    } 
  }
#ifdef SOME_DEBUG
  printf("new step - phase %d, norm %.18g, dual step %.18g, primal step %.18g\n",
	 phase,directionNorm,maximumDualStep,maximumPrimalStep);
#endif
  actualPrimalStep_=stepLength_*maximumPrimalStep;
  if (phase>=0&&actualPrimalStep_>1.0) {
    actualPrimalStep_=1.0;
  } 
  actualDualStep_=stepLength_*maximumDualStep;
  if (phase>=0&&actualDualStep_>1.0) {
    actualDualStep_=1.0;
  } 
#ifdef FULL_DEBUG
  if (phase==3){
    double minBeta = 0.1*mu_;
    double maxBeta = 10.0*mu_;
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change = rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (delta2Z_[iColumn]<minBeta||delta2Z_[iColumn]>maxBeta)
	    printf("3lower %d primal %g, dual %g, gap %g, old gap %g\n",
		   iColumn,primalValue,dualValue,gapProduct,delta2Z_[iColumn]);
	}  
	if (upperBound(iColumn)) {
	  double change = rhsU_[iColumn]-deltaX_[iColumn];
	  double dualValue=tVec[iColumn]+actualDualStep_*deltaT_[iColumn];
	  double primalValue=upperSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (delta2T_[iColumn]<minBeta||delta2T_[iColumn]>maxBeta)
	    printf("3upper %d primal %g, dual %g, gap %g, old gap %g\n",
		 iColumn,primalValue,dualValue,gapProduct,delta2T_[iColumn]);
	} 
      } 
    }
  }
#endif
#ifdef SOME_DEBUG
  {
    double largestL=0.0;
    double smallestL=COIN_DBL_MAX;
    double largestU=0.0;
    double smallestU=COIN_DBL_MAX;
    double sumL=0.0;
    double sumU=0.0;
    int nL=0;
    int nU=0;
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change = rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  largestL = max(largestL,gapProduct);
	  smallestL = min(smallestL,gapProduct);
	  nL++;
	  sumL += gapProduct;
	}  
	if (upperBound(iColumn)) {
	  double change = rhsU_[iColumn]-deltaX_[iColumn];
	  double dualValue=tVec[iColumn]+actualDualStep_*deltaT_[iColumn];
	  double primalValue=upperSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  largestU = max(largestU,gapProduct);
	  smallestU = min(smallestU,gapProduct);
	  nU++;
	  sumU += gapProduct;
	} 
      } 
    }
    double mu = (sumL+sumU)/((double) (nL+nU));

    double minBeta = 0.1*mu;
    double maxBeta = 10.0*mu;
    int nBL=0;
    int nAL=0;
    int nBU=0;
    int nAU=0;
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change = rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (gapProduct<minBeta)
	    nBL++;
	  else if (gapProduct>maxBeta)
	    nAL++;
	  //if (gapProduct<0.1*minBeta)
	  //printf("Lsmall one %d dual %g primal %g\n",iColumn,
	  //   dualValue,primalValue);
	}  
	if (upperBound(iColumn)) {
	  double change = rhsU_[iColumn]-deltaX_[iColumn];
	  double dualValue=tVec[iColumn]+actualDualStep_*deltaT_[iColumn];
	  double primalValue=upperSlack[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (gapProduct<minBeta)
	    nBU++;
	  else if (gapProduct>maxBeta)
	    nAU++;
	  //if (gapProduct<0.1*minBeta)
	  //printf("Usmall one %d dual %g primal %g\n",iColumn,
	  //   dualValue,primalValue);
	} 
      } 
    }
    printf("phase %d new mu %.18g new gap %.18g\n",phase,mu,sumL+sumU);
    printf("          %d lower, smallest %.18g, %d below - largest %.18g, %d above\n",
	   nL,smallestL,nBL,largestL,nAL);
    printf("          %d upper, smallest %.18g, %d below - largest %.18g, %d above\n",
	   nU,smallestU,nBU,largestU,nAU);
  }
#endif
  return directionNorm;
}
//#define KKT 2
/* Does solve. region1 is for deltaX (columns+rows), region2 for deltaPi (rows) */
void 
ClpPredictorCorrector::solveSystem(double * region1, double * region2,
				   const double * region1In, const double * region2In,
				   const double * saveRegion1, const double * saveRegion2,
				   bool gentleRefine)
{
  int iRow;
  int numberTotal = numberRows_+numberColumns_;
  if (region2In) {
    // normal
    for (iRow=0;iRow<numberRows_;iRow++)
      region2[iRow] = region2In[iRow];
  } else {
    // initial solution - (diagonal is 1 or 0)
    CoinZeroN(region2,numberRows_);
  }
  int iColumn;
#if KKT<2
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    region1[iColumn] = region1In[iColumn]*diagonal_[iColumn];
  multiplyAdd(region1+numberColumns_,numberRows_,-1.0,region2,1.0);
  matrix_->times(1.0,region1,region2);
  double maximumRHS = maximumAbsElement(region2,numberRows_);
  double scale=1.0;
  double unscale=1.0;
  if (maximumRHS>1.0e-30) {
    if (maximumRHS<=0.5) {
      double factor=2.0;
      while (maximumRHS<=0.5) {
	maximumRHS*=factor;
	scale*=factor;
      } /* endwhile */
    } else if (maximumRHS>=2.0&&maximumRHS<=COIN_DBL_MAX) {
      double factor=0.5;
      while (maximumRHS>=2.0) {
	maximumRHS*=factor;
	scale*=factor;
      } /* endwhile */
    } 
    unscale=diagonalScaleFactor_/scale;
  } else {
    //effectively zero
    scale=0.0;
    unscale=0.0;
  } 
  multiplyAdd(NULL,numberRows_,0.0,region2,scale);
  cholesky_->solve(region2);
  multiplyAdd(NULL,numberRows_,0.0,region2,unscale);
  multiplyAdd(region2,numberRows_,-1.0,region1+numberColumns_,0.0);
  CoinZeroN(region1,numberColumns_);
  matrix_->transposeTimes(1.0,region2,region1);
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    region1[iColumn] = (region1[iColumn]-region1In[iColumn])*diagonal_[iColumn];
#else
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    region1[iColumn] = region1In[iColumn];
  cholesky_->solveKKT(region1,region2,diagonal_,diagonalScaleFactor_);
#endif
  if (saveRegion2) {
    //refine?
    double scaleX=1.0;
    if (gentleRefine) 
      scaleX=0.8;
    multiplyAdd(saveRegion2,numberRows_,1.0,region2,scaleX);
    assert (saveRegion1);
    multiplyAdd(saveRegion1,numberTotal,1.0,region1,scaleX);
  } 
}
// findDirectionVector.
double ClpPredictorCorrector::findDirectionVector(const int phase)
{
  double projectionTolerance=projectionTolerance_;
  //temporary
  //projectionTolerance=1.0e-15;
  double errorCheck=0.9*maximumRHSError_/solutionNorm_;
  if (errorCheck>primalTolerance()) {
    if (errorCheck<projectionTolerance) {
      projectionTolerance=errorCheck;
    } 
  } else {
    if (primalTolerance()<projectionTolerance) {
      projectionTolerance=primalTolerance();
    } 
  } 
  double * newError = new double [numberRows_];
  double * workArray = workArray_;
  int numberTotal = numberRows_+numberColumns_;
  //if flagged then entries zero so can do
  // For KKT separate out
#ifndef KKT
  int iColumn;
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    deltaX_[iColumn] = workArray[iColumn] - solution_[iColumn];
  multiplyAdd(deltaX_+numberColumns_,numberRows_,-1.0,deltaY_,0.0);
  matrix_->times(1.0,deltaX_,deltaY_);
#else
  // regions in will be workArray and newError
  // regions out deltaX_ and deltaY_
  multiplyAdd(solution_+numberColumns_,numberRows_,1.0,newError,0.0);
  matrix_->times(-1.0,solution_,newError);
  double * region1Save=NULL;//for refinement
#endif
  bool goodSolve=false;
  double * regionSave=NULL;//for refinement
  int numberTries=0;
  double relativeError=COIN_DBL_MAX;
  double tryError=1.0e31;
  while (!goodSolve&&numberTries<30) {
    double lastError=relativeError;
    goodSolve=true;
    double maximumRHS = maximumAbsElement(deltaY_,numberRows_);
    double saveMaximum = maximumRHS;
#ifndef KKT
    double scale=1.0;
    double unscale=1.0;
    if (maximumRHS>1.0e-30) {
      if (maximumRHS<=0.5) {
        double factor=2.0;
        while (maximumRHS<=0.5) {
          maximumRHS*=factor;
          scale*=factor;
        } /* endwhile */
      } else if (maximumRHS>=2.0&&maximumRHS<=COIN_DBL_MAX) {
        double factor=0.5;
        while (maximumRHS>=2.0) {
          maximumRHS*=factor;
          scale*=factor;
        } /* endwhile */
      } 
      unscale=diagonalScaleFactor_/scale;
    } else {
      //effectively zero
      scale=0.0;
      unscale=0.0;
    } 
    //printf("--putting scales to 1.0\n");
    //scale=1.0;
    //unscale=1.0;
    multiplyAdd(NULL,numberRows_,0.0,deltaY_,scale);
    cholesky_->solve(deltaY_);
    multiplyAdd(NULL,numberRows_,0.0,deltaY_,unscale);
#if 0
 {
   printf("deltay\n");
   for (int i=0;i<numberRows_;i++) 
     printf("%d %.18g\n",i,deltaY_[i]);
   }
    exit(66);
#endif
    if (numberTries) {
      //refine?
      double scaleX=1.0;
      if (lastError>1.0e-5) 
        scaleX=0.8;
      multiplyAdd(regionSave,numberRows_,1.0,deltaY_,scaleX);
    } 
    //CoinZeroN(newError,numberRows_);
    multiplyAdd(deltaY_,numberRows_,-1.0,deltaX_+numberColumns_,0.0);
    CoinZeroN(deltaX_,numberColumns_);
    matrix_->transposeTimes(1.0,deltaY_,deltaX_);
    //if flagged then entries zero so can do
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      deltaX_[iColumn] = deltaX_[iColumn]*diagonal_[iColumn]
	-workArray[iColumn];
#else
    solveSystem(deltaX_, deltaY_,
		workArray,newError,region1Save,regionSave,lastError>1.0e-5);
#endif
    multiplyAdd(deltaX_+numberColumns_,numberRows_,-1.0,newError,0.0);
    matrix_->times(1.0,deltaX_,newError);
    numberTries++;
    
    //now add in old Ax - doing extra checking
    double maximumRHSError=0.0;
    double maximumRHSChange=0.0;
    int iRow;
    char * dropped = cholesky_->rowsDropped();
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (!dropped[iRow]) {
	double newValue=newError[iRow];
	double oldValue=errorRegion_[iRow];
	//severity of errors depend on signs
	//**later                                                             */
	if (fabs(newValue)>maximumRHSChange) {
	  maximumRHSChange=fabs(newValue);
	} 
	double result=newValue+oldValue;
	if (fabs(result)>maximumRHSError) {
	  maximumRHSError=fabs(result);
	} 
	newError[iRow]=result;
      } else {
	double newValue=newError[iRow];
	double oldValue=errorRegion_[iRow];
	if (fabs(newValue)>maximumRHSChange) {
	  maximumRHSChange=fabs(newValue);
	} 
	double result=newValue+oldValue;
	newError[iRow]=result;
	//newError[iRow]=0.0;
	//assert(deltaY_[iRow]==0.0);
	deltaY_[iRow]=0.0;
      } 
    } 
    relativeError = maximumRHSError/solutionNorm_;
    relativeError = maximumRHSError/saveMaximum;
    if (relativeError>tryError) 
      relativeError=tryError;
    if (relativeError<lastError) {
      maximumRHSChange_= maximumRHSChange;
      if (relativeError>1.0e-9
	  ||numberTries>1) {
	handler_->message(CLP_BARRIER_ACCURACY,messages_)
	  <<phase<<numberTries<<relativeError
	  <<CoinMessageEol;
      } 
      if (relativeError>projectionTolerance&&numberTries<=3) {
        //try and refine
        goodSolve=false;
      } 
      //*** extra test here
      if (!goodSolve) {
        if (!regionSave) {
          regionSave = new double [numberRows_];
#ifdef KKT
	  region1Save = new double [numberTotal];
#endif
        } 
        CoinMemcpyN(deltaY_,numberRows_,regionSave);
#ifndef KKT
	multiplyAdd(newError,numberRows_,-1.0,deltaY_,0.0);
#else
        CoinMemcpyN(deltaX_,numberTotal,region1Save);
	// and back to input region
        CoinMemcpyN(deltaY_,numberRows_,newError);
#endif
      } 
    } else {
      //std::cout <<" worse residual = "<<relativeError;
      //bring back previous
      relativeError=lastError;
      CoinMemcpyN(regionSave,numberRows_,deltaY_);
#ifndef KKT
      multiplyAdd(deltaY_,numberRows_,-1.0,deltaX_+numberColumns_,0.0);
      CoinZeroN(deltaX_,numberColumns_);
      matrix_->transposeTimes(1.0,deltaY_,deltaX_);
      //if flagged then entries zero so can do
      for (iColumn=0;iColumn<numberTotal;iColumn++)
	deltaX_[iColumn] = deltaX_[iColumn]*diagonal_[iColumn]
	  -workArray[iColumn];
#else
        CoinMemcpyN(region1Save,numberTotal,deltaX_);
#endif
    } 
  } /* endwhile */
  delete [] regionSave;
#ifdef KKT
  delete [] region1Save;
  int iColumn;
#endif
  delete [] newError;
  // now rest
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double extra=eExtra;

  for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    deltaSU_[iColumn]=0.0;
    deltaSL_[iColumn]=0.0;
    deltaZ_[iColumn]=0.0;
    deltaT_[iColumn]=0.0;
    if (!flagged(iColumn)) {
      double deltaX = deltaX_[iColumn];
      if (lowerBound(iColumn)) {
	double deltaSL = rhsL_[iColumn]+deltaX;
	double slack = lowerSlack[iColumn]+extra;
	deltaSL_[iColumn]=deltaSL;
	deltaZ_[iColumn]=(rhsZ_[iColumn]-zVec[iColumn]*deltaSL)/slack;
      } 
      if (upperBound(iColumn)) {
	double deltaSU = rhsU_[iColumn]-deltaX;
	double slack = upperSlack[iColumn]+extra;
	deltaSU_[iColumn]=deltaSU;
	deltaT_[iColumn]=(rhsT_[iColumn]-tVec[iColumn]*deltaSU)/slack;
      }
    }
  } 
  return relativeError;
}
// createSolution.  Creates solution from scratch
int ClpPredictorCorrector::createSolution()
{
  int numberTotal = numberRows_+numberColumns_;
  int iColumn;
  double tolerance = primalTolerance();
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    if (upper_[iColumn]-lower_[iColumn]>tolerance) 
      clearFixed(iColumn);
    else
      setFixed(iColumn);
  }
  double initialValue =1.0e-12;

  double maximumObjective=1.0e-12;
  double objectiveNorm2=0.0;
  getNorms(cost_,numberTotal,maximumObjective,objectiveNorm2);
  if (maximumObjective<1.0e-12) {
    maximumObjective=1.0e-12;
  } 
  objectiveNorm2=sqrt(objectiveNorm2)/(double) numberTotal;
  objectiveNorm_=maximumObjective;
  scaleFactor_=1.0;
  if (maximumObjective>0.0) {
    if (maximumObjective<1.0) {
      scaleFactor_=maximumObjective;
    } else if (maximumObjective>1.0e4) {
      scaleFactor_=maximumObjective/1.0e4;
    } 
  } 
  if (scaleFactor_!=1.0) {
    objectiveNorm2*=scaleFactor_;
    multiplyAdd(NULL,numberTotal,0.0,cost_,1.0/scaleFactor_);
    objectiveNorm_=maximumObjective/scaleFactor_;
  } 
  baseObjectiveNorm_=objectiveNorm_;
  //accumulate fixed in dj region (as spare)
  //accumulate primal solution in primal region
  //DZ in lowerDual
  //DW in upperDual
  double infiniteCheck=1.0e40;
  //double     fakeCheck=1.0e10;
  //use deltaX region for work region
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    double primalValue = solution_[iColumn];
    clearFlagged(iColumn);
    clearFixedOrFree(iColumn);
    clearLowerBound(iColumn);
    clearUpperBound(iColumn);
    clearFakeLower(iColumn);
    clearFakeUpper(iColumn);
    if (!fixed(iColumn)) {
      dj_[iColumn]=0.0;
      diagonal_[iColumn]=1.0;
      deltaX_[iColumn]=1.0;
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      if (lowerValue>-infiniteCheck) {
        if (upperValue<infiniteCheck) {
          //upper and lower bounds
          setLowerBound(iColumn);
          setUpperBound(iColumn);
          if (lowerValue>=0.0) {
            solution_[iColumn]=lowerValue;
          } else if (upperValue<=0.0) {
            solution_[iColumn]=upperValue;
          } else {
            solution_[iColumn]=0.0;
          } 
        } else {
          //just lower bound
          setLowerBound(iColumn);
          if (lowerValue>=0.0) {
            solution_[iColumn]=lowerValue;
          } else {
            solution_[iColumn]=0.0;
          } 
        } 
      } else {
        if (upperValue<infiniteCheck) {
          //just upper bound
          setUpperBound(iColumn);
          if (upperValue<=0.0) {
            solution_[iColumn]=upperValue;
          } else {
            solution_[iColumn]=0.0;
          } 
        } else {
          //free
          setFixedOrFree(iColumn);
          solution_[iColumn]=0.0;
          //std::cout<<" free "<<i<<std::endl;
        } 
      } 
    } else {
      setFlagged(iColumn);
      setFixedOrFree(iColumn);
      setLowerBound(iColumn);
      setUpperBound(iColumn);
      dj_[iColumn]=primalValue;;
      solution_[iColumn]=lower_[iColumn];
      diagonal_[iColumn]=0.0;
      deltaX_[iColumn]=0.0;
    } 
  } 
  //   modify fixed RHS
  multiplyAdd(dj_+numberColumns_,numberRows_,-1.0,rhsFixRegion_,0.0);
  //   create plausible RHS?
  matrix_->times(-1.0,dj_,rhsFixRegion_);
  multiplyAdd(solution_+numberColumns_,numberRows_,1.0,errorRegion_,0.0);
  matrix_->times(-1.0,solution_,errorRegion_);
  rhsNorm_=maximumAbsElement(errorRegion_,numberRows_);
  if (rhsNorm_<1.0) {
    rhsNorm_=1.0;
  } 
  int * rowsDropped = new int [numberRows_];
  int returnCode=cholesky_->factorize(diagonal_,rowsDropped);
  if (returnCode==-1) {
    printf("Out of memory\n");
    problemStatus_=4;
    return -1;
  }
  if (cholesky_->status()) {
    std::cout << "singular on initial cholesky?" <<std::endl;
    cholesky_->resetRowsDropped();
    //cholesky_->factorize(rowDropped_);
    //if (cholesky_->status()) {
      //std::cout << "bad cholesky??? (after retry)" <<std::endl;
      //abort();
    //} 
  } 
  delete [] rowsDropped;
#ifndef KKT
  cholesky_->solve(errorRegion_);
  //create information for solution
  multiplyAdd(errorRegion_,numberRows_,-1.0,deltaX_+numberColumns_,0.0);
  CoinZeroN(deltaX_,numberColumns_);
  matrix_->transposeTimes(1.0,errorRegion_,deltaX_);
#else
  // reverse sign on solution
  multiplyAdd(NULL,numberRows_+numberColumns_,0.0,solution_,-1.0);
  solveSystem(deltaX_,errorRegion_,solution_,NULL,NULL,NULL,false);
#endif
  //do reduced costs
  CoinZeroN(dj_+numberColumns_,numberRows_);
  for ( iColumn=0;iColumn<numberColumns_;iColumn++) {
    double value=cost_[iColumn];
    if (lowerBound(iColumn)) {
      value+=linearPerturbation_;
    } 
    if (upperBound(iColumn)) {
      value-=linearPerturbation_;
    } 
    dj_[iColumn]=value;
  } 
  initialValue=1.0e2;
  if (rhsNorm_*1.0e-2>initialValue) {
    initialValue=rhsNorm_*1.0e-2;
  } 
  double smallestBoundDifference=COIN_DBL_MAX;
  double * fakeSolution = deltaX_;
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      if (lower_[iColumn]-fakeSolution[iColumn]>initialValue) {
        initialValue=lower_[iColumn]-fakeSolution[iColumn];
      } 
      if (fakeSolution[iColumn]-upper_[iColumn]>initialValue) {
        initialValue=fakeSolution[iColumn]-upper_[iColumn];
      } 
      if (upper_[iColumn]-lower_[iColumn]<smallestBoundDifference) {
        smallestBoundDifference=upper_[iColumn]-lower_[iColumn];
      } 
    } 
  } 
  solutionNorm_=1.0e-12;
  handler_->message(CLP_BARRIER_SAFE,messages_)
    <<initialValue<<objectiveNorm_
    <<CoinMessageEol;
  int strategy=0;
  double extra=1.0e-10;
  double largeGap=1.0e15;
  double safeObjectiveValue=2.0*objectiveNorm_;
  double safeFree=1.0e-1*initialValue;
  //printf("normal safe dual value of %g, primal value of %g\n",
  // safeObjectiveValue,initialValue);
  //safeObjectiveValue=max(2.0,1.0e-1*safeObjectiveValue);
  //initialValue=max(100.0,1.0e-1*initialValue);;
  //printf("temp safe dual value of %g, primal value of %g\n",
  // safeObjectiveValue,initialValue);
  double zwLarge=1.0e2*initialValue;
  //zwLarge=1.0e40;
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      double objectiveValue = cost_[iColumn];
      double newValue;
      double low;
      double high;
      double randomZ =0.1*CoinDrand48()+0.9;
      double randomW =0.1*CoinDrand48()+0.9;
      double setPrimal=initialValue;
      if (strategy==0) {
        randomZ=1.0;
        randomW=1.0;
      } 
      if (lowerBound(iColumn)) {
        if (upperBound(iColumn)) {
          //upper and lower bounds
          if (upperValue-lowerValue>2.0*setPrimal) {
            double fakeValue=fakeSolution[iColumn];
            if (fakeValue<lowerValue+setPrimal) {
              fakeValue=lowerValue+setPrimal;
            } 
            if (fakeValue>upperValue-setPrimal) {
              fakeValue=upperValue-setPrimal;
            } 
            newValue=fakeValue;
            low=fakeValue-lowerValue;
            high=upperValue-fakeValue;
          } else {
            newValue=0.5*(upperValue+lowerValue);
            low=setPrimal;
            high=setPrimal;
          } 
          low *= randomZ;
          high *= randomW;
          double s = low+extra;
          double ratioZ;
          if (s<zwLarge) {
            ratioZ=1.0;
          } else {
            ratioZ=sqrt(zwLarge/s);
          } 
          double t = high+extra;
          double ratioW;
          if (t<zwLarge) {
            ratioW=1.0;
          } else {
            ratioW=sqrt(zwLarge/t);
          } 
          //modify s and t
          if (s>largeGap) {
            s=largeGap;
          } 
          if (t>largeGap) {
            t=largeGap;
          } 
          //modify if long long way away from bound
          if (objectiveValue>=0.0) {
            zVec_[iColumn]=objectiveValue + safeObjectiveValue*ratioZ;
            tVec_[iColumn]=safeObjectiveValue*ratioW;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            tVec_[iColumn]=-objectiveValue + safeObjectiveValue*ratioW;
          } 
          diagonal_[iColumn] = (t*s)/(s*tVec_[iColumn]+t*zVec_[iColumn]);
        } else {
          //just lower bound
          double fakeValue=fakeSolution[iColumn];
          if (fakeValue<lowerValue+setPrimal) {
            fakeValue=lowerValue+setPrimal;
          } 
          newValue=fakeValue;
          low=fakeValue-lowerValue;
          low *= randomZ;
          high=0.0;
          double s = low+extra;
          double ratioZ;
          if (s<zwLarge) {
            ratioZ=1.0;
          } else {
            ratioZ=sqrt(zwLarge/s);
          } 
          //modify s
          if (s>largeGap) {
            s=largeGap;
          } 
          if (objectiveValue>=0.0) {
            zVec_[iColumn]=objectiveValue + safeObjectiveValue*ratioZ;
            tVec_[iColumn]=0.0;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            tVec_[iColumn]=0.0;
          } 
          diagonal_[iColumn] = s/zVec_[iColumn];
        } 
      } else {
        if (upperBound(iColumn)) {
          //just upper bound
          double fakeValue=fakeSolution[iColumn];
          if (fakeValue>upperValue-setPrimal) {
            fakeValue=upperValue-setPrimal;
          } 
          newValue=fakeValue;
          low=0.0;
          high=upperValue-fakeValue;
          high*=randomW;
          double t = high+extra;
          double ratioW;
          if (t<zwLarge) {
            ratioW=1.0;
          } else {
            ratioW=sqrt(zwLarge/t);
          } 
          //modify t
          if (t>largeGap) {
            t=largeGap;
          } 
          if (objectiveValue>=0.0) {
            zVec_[iColumn]=0.0;
            tVec_[iColumn]=safeObjectiveValue*ratioW;
          } else {
            zVec_[iColumn]=0.0;
            tVec_[iColumn]=-objectiveValue + safeObjectiveValue*ratioW;
          } 
          diagonal_[iColumn] =  t/tVec_[iColumn];
        } else {
          //free
          zVec_[iColumn]=safeObjectiveValue;
          tVec_[iColumn]=safeObjectiveValue;
          newValue=fakeSolution[iColumn];
          if (newValue>=0.0) {
            if (newValue<safeFree) {
              newValue=safeFree;
            } 
          } else {
            if (newValue>-safeFree) {
              newValue=-safeFree;
            } 
          } 
          if (fabs(newValue)>1.0) {
            diagonal_[iColumn]=fabs(newValue)/(zVec_[iColumn]+tVec_[iColumn]);
          } else {
            diagonal_[iColumn]=1.0/(zVec_[iColumn]+tVec_[iColumn]);
          } 
          low=0.0;
          high=0.0;
        } 
      } 
      lowerSlack_[iColumn]=low;
      upperSlack_[iColumn]=high;
      solution_[iColumn]=newValue;
    } else {
      lowerSlack_[iColumn]=0.0;
      upperSlack_[iColumn]=0.0;
      solution_[iColumn]=lower_[iColumn];
      zVec_[iColumn]=0.0;
      tVec_[iColumn]=0.0;
      diagonal_[iColumn]=0.0;
    } 
  }
  solutionNorm_ =  maximumAbsElement(solution_,numberTotal);
  // Set bounds
  largeGap = max(1.0e7,1.02*solutionNorm_);
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double newValue=solution_[iColumn];
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      if (newValue>lowerValue+largeGap&&newValue<upperValue-largeGap) {
	clearFixedOrFree(iColumn);
	setLowerBound(iColumn);
	setUpperBound(iColumn);
	double objectiveValue = cost_[iColumn];
	lowerValue=max(lowerValue,newValue-largeGap);
	upperValue=min(upperValue,newValue+largeGap);
	lower_[iColumn]=lowerValue;
	upper_[iColumn]=upperValue;
	double low;
	double high;
	double setPrimal=initialValue;
	//upper and lower bounds
	if (upperValue-lowerValue>2.0*setPrimal) {
	  low=newValue-lowerValue;
	  high=upperValue-newValue;
	} else {
	  newValue=0.5*(upperValue+lowerValue);
	  low=setPrimal;
	  high=setPrimal;
	} 
	double s = low+extra;
	double ratioZ;
	if (s<zwLarge) {
	  ratioZ=1.0;
	} else {
	  ratioZ=sqrt(zwLarge/s);
	} 
	double t = high+extra;
	double ratioW;
	if (t<zwLarge) {
	  ratioW=1.0;
	} else {
	  ratioW=sqrt(zwLarge/t);
	} 
	//modify s and t
	if (s>largeGap) {
	  s=largeGap;
	} 
	if (t>largeGap) {
	  t=largeGap;
	} 
	//modify if long long way away from bound
	if (objectiveValue>=0.0) {
	  zVec_[iColumn]=objectiveValue + safeObjectiveValue*ratioZ;
	  tVec_[iColumn]=safeObjectiveValue*ratioW;
	} else {
	  zVec_[iColumn]=safeObjectiveValue*ratioZ;
	  tVec_[iColumn]=-objectiveValue + safeObjectiveValue*ratioW;
	} 
	diagonal_[iColumn] = (t*s)/(s*tVec_[iColumn]+t*zVec_[iColumn]);
	lowerSlack_[iColumn]=low;
	upperSlack_[iColumn]=high;
      } 
    }
  }
#if 0
  if (solution_[0]>0.0) {
    for (int i=0;i<numberTotal;i++)
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	     diagonal_[i],fabs(dj_[i]),
	     lowerSlack_[i],zVec_[i],
	     upperSlack_[i],tVec_[i]);
  } else {
    for (int i=0;i<numberTotal;i++)
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	     diagonal_[i],fabs(dj_[i]),
	     upperSlack_[i],tVec_[i],
	     lowerSlack_[i],zVec_[i] );
  }
  exit(66);
#endif
  return 0;
}
// complementarityGap.  Computes gap
//phase 0=as is , 1 = after predictor , 2 after corrector
double ClpPredictorCorrector::complementarityGap(int & numberComplementarityPairs,
						 int & numberComplementarityItems,
						 const int phase)
{
  double gap=0.0;
  //seems to be same coding for phase = 1 or 2
  numberComplementarityPairs=0;
  numberComplementarityItems=0;
  double toleranceGap=0.0;
  double largestGap=0.0;
  double smallestGap=COIN_DBL_MAX;
  //seems to be same coding for phase = 1 or 2
  int numberNegativeGaps=0;
  double sumNegativeGap=0.0;
  double largeGap=1.0e2*solutionNorm_;
  if (largeGap<1.0e10) {
    largeGap=1.0e10;
  } 
  largeGap=1.0e30;
  double dualTolerance =  dblParam_[ClpDualTolerance];
  double primalTolerance =  dblParam_[ClpPrimalTolerance];
  dualTolerance=dualTolerance/scaleFactor_;
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * primal = solution_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * deltaZ = deltaZ_;
  double * deltaT = deltaT_;
  double * deltaX = deltaX_;
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!fixedOrFree(iColumn)) {
      numberComplementarityPairs++;
      //can collapse as if no lower bound both zVec and deltaZ 0.0
      if (lowerBound(iColumn)) {
	numberComplementarityItems++;
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=zVec[iColumn];
          primalValue=lowerSlack[iColumn];
        } else {
          double change;
	  change =primal[iColumn]+deltaX[iColumn]-lowerSlack[iColumn]-lower[iColumn];
          dualValue=zVec[iColumn]+actualDualStep_*deltaZ[iColumn];
          primalValue=lowerSlack[iColumn]+actualPrimalStep_*change;
        } 
        //reduce primalValue
        if (primalValue>largeGap) {
          primalValue=largeGap;
        } 
        double gapProduct=dualValue*primalValue;
        if (gapProduct<0.0) {
          //cout<<"negative gap component "<<iColumn<<" "<<dualValue<<" "<<
              //primalValue<<endl;
          numberNegativeGaps++;
          sumNegativeGap-=gapProduct;
          gapProduct=0.0;
        } 
        gap+=gapProduct;
        if (gapProduct>largestGap) {
          largestGap=gapProduct;
        }
	smallestGap = min(smallestGap,gapProduct);
        if (dualValue>dualTolerance&&primalValue>primalTolerance) {
          toleranceGap+=dualValue*primalValue;
        } 
      } 
      if (upperBound(iColumn)) {
	numberComplementarityItems++;
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=tVec[iColumn];
          primalValue=upperSlack[iColumn];
        } else {
          double change;
	  change =upper[iColumn]-primal[iColumn]-deltaX[iColumn]-upperSlack[iColumn];
          dualValue=tVec[iColumn]+actualDualStep_*deltaT[iColumn];
          primalValue=upperSlack[iColumn]+actualPrimalStep_*change;
        } 
        //reduce primalValue
        if (primalValue>largeGap) {
          primalValue=largeGap;
        } 
        double gapProduct=dualValue*primalValue;
        if (gapProduct<0.0) {
          //cout<<"negative gap component "<<iColumn<<" "<<dualValue<<" "<<
              //primalValue<<endl;
          numberNegativeGaps++;
          sumNegativeGap-=gapProduct;
          gapProduct=0.0;
        } 
        gap+=gapProduct;
        if (gapProduct>largestGap) {
          largestGap=gapProduct;
        } 
        if (dualValue>dualTolerance&&primalValue>primalTolerance) {
          toleranceGap+=dualValue*primalValue;
        } 
      } 
    } 
  } 
  if (!phase&&numberNegativeGaps) {
      handler_->message(CLP_BARRIER_NEGATIVE_GAPS,messages_)
    <<numberNegativeGaps<<sumNegativeGap
    <<CoinMessageEol;
  } 
  
  //in case all free!
  if (!numberComplementarityPairs) {
    numberComplementarityPairs=1;
  } 
  //printf("gap %g - smallest %g, largest %g, pairs %d\n",
  // gap,smallestGap,largestGap,numberComplementarityPairs);
  return gap;
}
// setupForSolve.
//phase 0=affine , 1 = corrector , 2 = primal-dual
void ClpPredictorCorrector::setupForSolve(const int phase)
{
  double extra =eExtra;
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * primal = solution_;
  double * dj = dj_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * workArray = workArray_;
  int numberTotal = numberRows_ + numberColumns_;
  int iColumn;
#ifdef SOME_DEBUG
  printf("phase %d in setupForSolve, mu %.18g\n",phase,mu_);
#endif
  switch (phase) {
  case 0:
    CoinMemcpyN(errorRegion_,numberRows_,rhsB_);
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsC_[iColumn]=0.0;
      rhsU_[iColumn]=0.0;
      rhsL_[iColumn]=0.0;
      rhsZ_[iColumn]=0.0;
      rhsT_[iColumn]=0.0;
      if (!flagged(iColumn)) {
        rhsC_[iColumn] = dj[iColumn]-zVec[iColumn]+tVec[iColumn];
        if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = -zVec[iColumn]*(lowerSlack[iColumn]+extra);
	  rhsL_[iColumn] = min(0.0,primal[iColumn]-lowerSlack[iColumn]-lower[iColumn]);
        } 
        if (upperBound(iColumn)) {
	  rhsT_[iColumn] = -tVec[iColumn]*(upperSlack[iColumn]+extra);
	  rhsU_[iColumn] = min(0.0,-(primal[iColumn] + upperSlack[iColumn] - upper[iColumn])); 
        }
      }
    } 
    if (0) {
      int i=1324;
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,solution_[i],
	     diagonal_[i],dj_[i],
	     lowerSlack_[i],zVec_[i],
	     upperSlack_[i],tVec_[i]);
      printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	     rhsZ_[i],rhsL_[i],
	     rhsT_[i],rhsU_[i]);
    }
#if 0
    for (int i=0;i<3;i++) {
      if (!fabs(rhsZ_[i]))
	rhsZ_[i]=0.0;
      if (!fabs(rhsT_[i]))
	rhsT_[i]=0.0;
      if (!fabs(rhsU_[i]))
	rhsU_[i]=0.0;
      if (!fabs(rhsL_[i]))
	rhsL_[i]=0.0;
    }
    if (solution_[0]>0.0) {
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,solution_[i],
	       diagonal_[i],dj_[i],
	       lowerSlack_[i],zVec_[i],
	       upperSlack_[i],tVec_[i]);
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	       rhsZ_[i],rhsL_[i],
	       rhsT_[i],rhsU_[i]);
    } else {
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,solution_[i],
	       diagonal_[i],dj_[i],
	       lowerSlack_[i],zVec_[i],
	       upperSlack_[i],tVec_[i]);
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	       rhsZ_[i],rhsL_[i],
	       rhsT_[i],rhsU_[i]);
    }
#endif
    break;
  case 1:
    // could be stored in delta2?
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsC_[iColumn]=0.0;
      //rhsU_[iColumn]=0.0;
      //rhsL_[iColumn]=0.0;
      rhsZ_[iColumn]=0.0;
      rhsT_[iColumn]=0.0;
      if (!flagged(iColumn)) {
        rhsC_[iColumn] = dj[iColumn]-zVec[iColumn]+tVec[iColumn];
        if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = mu_ -zVec[iColumn]*(lowerSlack[iColumn]+extra)
	    - deltaZ_[iColumn]*deltaX_[iColumn];
	  // To bring in line with OSL
	  rhsZ_[iColumn] -= deltaZ_[iColumn]*rhsL_[iColumn];
	  //rhsZ_[iColumn] = mu_ -deltaZ_[iColumn]*rhsL_[iColumn]
	  //- deltaZ_[iColumn]*deltaX_[iColumn];
        } 
        if (upperBound(iColumn)) {
	  rhsT_[iColumn] = mu_ -tVec[iColumn]*(upperSlack[iColumn]+extra)
	    +deltaT_[iColumn]*deltaX_[iColumn];
	  // To bring in line with OSL
	  rhsT_[iColumn] -= deltaT_[iColumn]*rhsU_[iColumn];
	  //rhsT_[iColumn] = mu_ -deltaT_[iColumn]*rhsU_[iColumn]
	  //+deltaT_[iColumn]*deltaX_[iColumn];
        } 
      } 
    } 
#if 0
    for (int i=0;i<numberTotal;i++) {
      if (!fabs(rhsZ_[i]))
	rhsZ_[i]=0.0;
      if (!fabs(rhsT_[i]))
	rhsT_[i]=0.0;
      if (!fabs(rhsU_[i]))
	rhsU_[i]=0.0;
      if (!fabs(rhsL_[i]))
	rhsL_[i]=0.0;
    }
    if (solution_[0]>0.0) {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	       diagonal_[i],fabs(dj_[i]),
	       lowerSlack_[i],zVec_[i],
	       upperSlack_[i],tVec_[i]);
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,fabs(rhsC_[i]),
	       rhsZ_[i],rhsL_[i],
	       rhsT_[i],rhsU_[i]);
    } else {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	       diagonal_[i],fabs(dj_[i]),
	       upperSlack_[i],tVec_[i],
	       lowerSlack_[i],zVec_[i] );
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,fabs(rhsC_[i]),
	       rhsT_[i],rhsU_[i],
	       rhsZ_[i],rhsL_[i]);
    }
    exit(66);
#endif
    break;
  case 2:
    CoinMemcpyN(errorRegion_,numberRows_,rhsB_);
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsC_[iColumn]=0.0;
      rhsZ_[iColumn]=0.0;
      rhsT_[iColumn]=0.0;
      if (!flagged(iColumn)) {
        rhsC_[iColumn] = dj[iColumn]-zVec[iColumn]+tVec[iColumn];
        if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = mu_ - zVec[iColumn]*(lowerSlack[iColumn]+extra);
        } 
        if (upperBound(iColumn)) {
	  rhsT_[iColumn] = mu_ - tVec[iColumn]*(upperSlack[iColumn]+extra);
        }
      }
    } 
    break;
  case 3:
    {
      double minBeta = 0.1*mu_;
      double maxBeta = 10.0*mu_;
      double dualStep = min(1.0,actualDualStep_+0.1);
      double primalStep = min(1.0,actualPrimalStep_+0.1);
#ifdef SOME_DEBUG
      printf("good complementarity range %g to %g\n",minBeta,maxBeta);
#endif
      //minBeta=0.0;
      //maxBeta=COIN_DBL_MAX;
      for (iColumn=0;iColumn<numberTotal;iColumn++) {
	if (!flagged(iColumn)) {
	  if (lowerBound(iColumn)) {
	    double change = rhsL_[iColumn] + deltaX_[iColumn];
	    double dualValue=zVec[iColumn]+dualStep*deltaZ_[iColumn];
	    double primalValue=lowerSlack[iColumn]+primalStep*change;
	    double gapProduct=dualValue*primalValue;
	    if (gapProduct>0.0&&dualValue<0.0)
	      gapProduct = - gapProduct;
#ifdef FULL_DEBUG
	    delta2Z_[iColumn]=gapProduct;
	    if (delta2Z_[iColumn]<minBeta||delta2Z_[iColumn]>maxBeta)
	      printf("lower %d primal %g, dual %g, gap %g\n",
		     iColumn,primalValue,dualValue,gapProduct);
#endif
	    double value= 0.0;
	    if (gapProduct<minBeta) {
	      value= 2.0*(minBeta-gapProduct);
	      value= (mu_-gapProduct);
	      value= (minBeta-gapProduct);
	      assert (value>0.0);
	    } else if (gapProduct>maxBeta) {
	      value= max(maxBeta-gapProduct,-maxBeta);
	      assert (value<0.0);
	    }
	    //#define AGAIN
#ifdef AGAIN
	    //rhsZ_[iColumn] = mu_ -zVec[iColumn]*(lowerSlack[iColumn]+extra)
	    //- deltaZ_[iColumn]*deltaX_[iColumn];
	    // To bring in line with OSL
	    rhsZ_[iColumn] += deltaZ_[iColumn]*rhsL_[iColumn];
#endif
	    rhsZ_[iColumn] += value;
	  }  
	  if (upperBound(iColumn)) {
	    double change = rhsU_[iColumn]-deltaX_[iColumn];
	    double dualValue=tVec[iColumn]+dualStep*deltaT_[iColumn];
	    double primalValue=upperSlack[iColumn]+primalStep*change;
	    double gapProduct=dualValue*primalValue;
	    if (gapProduct>0.0&&dualValue<0.0)
	      gapProduct = - gapProduct;
#ifdef FULL_DEBUG
	    delta2T_[iColumn]=gapProduct;
	    if (delta2T_[iColumn]<minBeta||delta2T_[iColumn]>maxBeta)
	      printf("upper %d primal %g, dual %g, gap %g\n",
		     iColumn,primalValue,dualValue,gapProduct);
#endif
	    double value= 0.0;
	    if (gapProduct<minBeta) {
	      value= (minBeta-gapProduct);
	      assert (value>0.0);
	    } else if (gapProduct>maxBeta) {
	      value= max(maxBeta-gapProduct,-maxBeta);
	      assert (value<0.0);
	    }
#ifdef AGAIN
	    //rhsT_[iColumn] = mu_ -tVec[iColumn]*(upperSlack[iColumn]+extra)
	    //+deltaT_[iColumn]*deltaX_[iColumn];
	    // To bring in line with OSL
	    rhsT_[iColumn] -= deltaT_[iColumn]*rhsU_[iColumn];
#endif
	    rhsT_[iColumn] += value;
	  } 
	} 
      }
    }
    break;
  } /* endswitch */
#ifndef KKT
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    double value = rhsC_[iColumn];
    double zValue = rhsZ_[iColumn];
    double tValue = rhsT_[iColumn];
#if 1
    if (phase==0) {
      // more accurate
      value = dj[iColumn];
      zValue=0.0;
      tValue=0.0;
    } else if (phase==2) {
      // more accurate
      value = dj[iColumn];
      zValue=mu_;
      tValue=mu_;
    }
#endif
    assert (rhsL_[iColumn]<=0.0);
    assert (rhsU_[iColumn]<=0.0);
    if (lowerBound(iColumn)) {
      value += (zVec[iColumn]*rhsL_[iColumn]-zValue)/
	(lowerSlack[iColumn]+extra);
    }
    if (upperBound(iColumn)) {
      value += (tValue-tVec[iColumn]*rhsU_[iColumn])/
	(upperSlack[iColumn]+extra);
    }
    workArray[iColumn]=diagonal_[iColumn]*value;
  } 
#if 0
    if (solution_[0]>0.0) {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g\n",i,workArray[i]);
    } else {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g\n",i,workArray[i]);
    }
    exit(66);
#endif
#else
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    double value = rhsC_[iColumn];
    if (lowerBound(iColumn)) {
      value -= rhsZ_[iColumn]/(lowerSlack[iColumn]+extra);
      value += zVec[iColumn]*rhsL_[iColumn]/(lowerSlack[iColumn]+extra);
    }
    if (upperBound(iColumn)) {
      value += rhsT_[iColumn]/(upperSlack[iColumn]+extra);
      value -= tVec[iColumn]*rhsU_[iColumn]/(upperSlack[iColumn]+extra);
    }
      workArray[iColumn]=value;
  }
#endif 
}
//method: sees if looks plausible change in complementarity
bool ClpPredictorCorrector::checkGoodMove(const bool doCorrector,double & bestNextGap)
{
  const double beta3 = 0.99997;
  bool goodMove=false;
  int nextNumber;
  int nextNumberItems;
  int numberTotal = numberRows_+numberColumns_;
  double returnGap=bestNextGap;
  double nextGap=complementarityGap(nextNumber,nextNumberItems,2);
  if (nextGap>bestNextGap&&nextGap>-0.9*complementarityGap_) {
#ifdef SOME_DEBUG
    printf("checkGood phase 1 next gap %.18g, phase 0 %.18g, old gap %.18g\n",
	   nextGap,bestNextGap,complementarityGap_);
#endif
    return false;
  } else {
    returnGap=nextGap;
  }
  double step;
  if (actualDualStep_>actualPrimalStep_) {
    step=actualDualStep_;
  } else {
    step=actualPrimalStep_;
  } 
  double testValue=1.0-step*(1.0-beta3);
  //testValue=0.0;
  testValue*=complementarityGap_;
  if (nextGap<testValue) {
    //std::cout <<"predicted duality gap "<<nextGap<<std::endl;
    goodMove=true;
  } else if(doCorrector) {
    //if (actualDualStep_<actualPrimalStep_) {
      //step=actualDualStep_;
    //} else {
      //step=actualPrimalStep_;
    //} 
    double gap = bestNextGap;
    goodMove=checkGoodMove2(step,gap);
    if (goodMove)
      returnGap=gap;
  } else {
    goodMove=true;
  } 
  if (!goodMove) {
    //try smaller of two
    if (actualDualStep_<actualPrimalStep_) {
      step=actualDualStep_;
    } else {
      step=actualPrimalStep_;
    } 
    if (step>1.0) {
      step=1.0;
    } 
    actualPrimalStep_=step;
    actualDualStep_=step;
    while (!goodMove) {
      double gap = bestNextGap;
      goodMove=checkGoodMove2(step,gap);
      if (goodMove)
	returnGap=gap;
      if (step<1.0e-10) {
        break;
      } 
      step*=0.5;
      actualPrimalStep_=step;
      actualDualStep_=step;
    } /* endwhile */
    if (doCorrector) {
      //say bad move if both small
      if (numberIterations_&1) {
        if (actualPrimalStep_<1.0e-2&&actualDualStep_<1.0e-2) {
          goodMove=false;
        } 
      } else {
        if (actualPrimalStep_<1.0e-5&&actualDualStep_<1.0e-5) {
          goodMove=false;
        } 
        if (actualPrimalStep_*actualDualStep_<1.0e-20) {
          goodMove=false;
        } 
      } 
    } 
  } 
  if (goodMove) {
    //compute delta in objectives
    double deltaObjectivePrimal=0.0;
    double deltaObjectiveDual=
      innerProduct(deltaY_,numberRows_,
		   rhsFixRegion_);
    double error=0.0;
    double * workArray = workArray_;
    CoinZeroN(workArray,numberColumns_);
    CoinMemcpyN(deltaY_,numberRows_,workArray+numberColumns_);
    matrix_->transposeTimes(-1.0,deltaY_,workArray);
    double * deltaZ = deltaZ_;
    double * deltaT = deltaT_;
    double * lower = lower_;
    double * upper = upper_;
    //direction vector in deltaX
    double * deltaX = deltaX_;
    double * cost = cost_;
    //double sumPerturbCost=0.0;
    for (int iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  //sumPerturbCost+=deltaX[iColumn];
	  deltaObjectiveDual+=deltaZ[iColumn]*lower[iColumn];
	} 
	if (upperBound(iColumn)) {
	  //sumPerturbCost-=deltaX[iColumn];
	  deltaObjectiveDual-=deltaT[iColumn]*upper[iColumn];
	} 
	double change = fabs(workArray[iColumn]-deltaZ[iColumn]+deltaT[iColumn]);
	error = max (change,error);
      } 
      deltaObjectivePrimal += cost[iColumn] * deltaX[iColumn];
    } 
    //deltaObjectivePrimal+=sumPerturbCost*linearPerturbation_;
    double testValue;
    if (error>0.0) {
      testValue=1.0e1*maximumDualError_/error;
    } else {
      testValue=1.0e1;
    } 
    if (testValue<actualDualStep_) {
      handler_->message(CLP_BARRIER_REDUCING,messages_)
      <<"dual"<<actualDualStep_
      << testValue
      <<CoinMessageEol;
      actualDualStep_=testValue;
    } 
  } 
  if (maximumRHSError_<1.0e1*solutionNorm_*primalTolerance()
                            &&maximumRHSChange_>1.0e-16*solutionNorm_) {
    //check change in AX not too much
    //??? could be dropped row going infeasible
    double ratio = 1.0e1*maximumRHSError_/maximumRHSChange_;
    if (ratio<actualPrimalStep_) {
      handler_->message(CLP_BARRIER_REDUCING,messages_)
      <<"primal"<<actualPrimalStep_
      <<ratio
      <<CoinMessageEol;
      if (ratio>1.0e-6) {
        actualPrimalStep_=ratio;
      } else {
        actualPrimalStep_=ratio;
        //std::cout <<"sign we should be stopping"<<std::endl;
      } 
    } 
  }
  if (goodMove)
    bestNextGap=returnGap;
  return goodMove;
}
//:  checks for one step size
bool ClpPredictorCorrector::checkGoodMove2(const double move,double & bestNextGap)
{
  double complementarityMultiplier =1.0/numberComplementarityPairs_;
  const double gamma = 1.0e-8;
  const double gammap = 1.0e-8;
  const double gammad = 1.0e-8;
  int nextNumber;
  int nextNumberItems;
  double nextGap=complementarityGap(nextNumber,nextNumberItems,2);
  if (nextGap>bestNextGap)
    return false;
  double lowerBoundGap = gamma*nextGap*complementarityMultiplier;
  bool goodMove=true;
  double * deltaZ = deltaZ_;
  double * deltaT = deltaT_;
  double * deltaSL = deltaSL_;
  double * deltaSU = deltaSU_;
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!flagged(iColumn)) {
      if (lowerBound(iColumn)) {
	double part1=lowerSlack[iColumn]+actualPrimalStep_*deltaSL[iColumn];
	double part2=zVec[iColumn]+actualDualStep_*deltaZ[iColumn];
	if (part1*part2<lowerBoundGap) {
	  goodMove=false;
	  break;
	} 
      } 
      if (upperBound(iColumn)) {
	double part1=upperSlack[iColumn]+actualPrimalStep_*deltaT[iColumn];
	double part2=tVec[iColumn]+actualDualStep_*deltaSU[iColumn];
	if (part1*part2<lowerBoundGap) {
	  goodMove=false;
	  break;
	} 
      } 
    } 
  } 
  //      Satisfy g_p(alpha)?
  if (rhsNorm_>solutionNorm_) {
    solutionNorm_=rhsNorm_;
  } 
  double errorCheck=maximumRHSError_/solutionNorm_;
  if (errorCheck<maximumBoundInfeasibility_) {
    errorCheck=maximumBoundInfeasibility_;
  } 
  //scale
  if ((1.0-move)*errorCheck>primalTolerance()) {
    if (nextGap<gammap*(1.0-move)*errorCheck) {
      goodMove=false;
    } 
  } 
  //      Satisfy g_d(alpha)?
  errorCheck=maximumDualError_/objectiveNorm_;
  if ((1.0-move)*errorCheck>dualTolerance()) {
    if (nextGap<gammad*(1.0-move)*errorCheck) {
      goodMove=false;
    } 
  }
  if (goodMove)
    bestNextGap=nextGap;
  return goodMove;
}
// updateSolution.  Updates solution at end of iteration
//returns number fixed
int ClpPredictorCorrector::updateSolution(double nextGap)
{
  //update pi
  multiplyAdd(deltaY_,numberRows_,actualDualStep_,dual_,1.0);
  CoinZeroN(errorRegion_,numberRows_);
  CoinZeroN(rhsFixRegion_,numberRows_);
  double maximumRhsInfeasibility=0.0;
  double maximumBoundInfeasibility=0.0;
  double maximumDualError=1.0e-12;
  double primalObjectiveValue=0.0;
  double dualObjectiveValue=0.0;
  double solutionNorm=1.0e-12;
  int numberKilled=0;
  double freeMultiplier=1.0e6;
  double trueNorm =diagonalNorm_/diagonalScaleFactor_;
  if (freeMultiplier<trueNorm) {
    freeMultiplier=trueNorm;
  } 
  if (freeMultiplier>1.0e12) {
    freeMultiplier=1.0e12;
  } 
  freeMultiplier=0.5/freeMultiplier;
  double condition = fabs(cholesky_->choleskyCondition());
  bool caution;
  if ((condition<1.0e10&&trueNorm<1.0e12)||numberIterations_<20) {
    caution=false;
  } else {
    caution=true;
  } 
  // do reduced costs
  CoinMemcpyN(dual_,numberRows_,dj_+numberColumns_);
  CoinMemcpyN(cost_,numberColumns_,dj_);
  matrix_->transposeTimes(-1.0,dual_,dj_);
  double extra=eExtra;
  const double largeFactor=1.0e2;
  double largeGap=largeFactor*solutionNorm_;
  if (largeGap<largeFactor) {
    largeGap=largeFactor;
  } 
  double dualFake=0.0;
  double dualTolerance =  dblParam_[ClpDualTolerance];
  dualTolerance=dualTolerance/scaleFactor_;
  if (dualTolerance<1.0e-12) {
    dualTolerance=1.0e-12;
  } 
  double offsetObjective=0.0;
  const double killTolerance=primalTolerance();
  double qDiagonal;
  if (mu_<1.0) {
    qDiagonal=1.0e-8;
  } else {
    qDiagonal=1.0e-8*mu_;
  } 
  //qDiagonal *= 1.0e2;
  //largest allowable ratio of lowerSlack/zVec (etc)
  double largestRatio;
  double epsilonBase;
  double diagonalLimit;
  if (!caution) {
    epsilonBase=eBase;
    largestRatio=eRatio;
    diagonalLimit=eDiagonal;
  } else {
    epsilonBase=eBaseCaution;
    largestRatio=eRatioCaution;
    diagonalLimit=eDiagonalCaution;
  } 
  double smallGap=1.0e2;
  double maximumDJInfeasibility=0.0;
  int numberIncreased=0;
  int numberDecreased=0;
  double largestDiagonal=0.0;
  double smallestDiagonal=1.0e50;
  double * zVec = zVec_;
  double * tVec = tVec_;
  double * primal = solution_;
  double * dual = dj_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * deltaZ = deltaZ_;
  double * deltaT = deltaT_;
  double * diagonal = diagonal_;
  //direction vector in deltaX
  double * deltaX = deltaX_;
  double * cost = cost_;
  double largeGap2 = max(1.0e7,1.0e2*solutionNorm_);
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!flagged(iColumn)) {
      double reducedCost=dual[iColumn];
      bool thisKilled=false;
      double zValue = zVec[iColumn] + actualDualStep_*deltaZ[iColumn];
      double tValue = tVec[iColumn] + actualDualStep_*deltaT[iColumn];
      zVec[iColumn]=zValue;
      tVec[iColumn]=tValue;
      double thisWeight=deltaX[iColumn];
      double oldPrimal = primal[iColumn];
      double newPrimal=primal[iColumn]+actualPrimalStep_*thisWeight;
      double primalObjectiveThis=newPrimal*cost[iColumn];
      double dualObjectiveThis=0.0;
      double t=extra;
      double s=extra;
      double kill;
      if (fabs(newPrimal)>1.0e4) {
        kill=killTolerance*1.0e-4*newPrimal;
        //kill=killTolerance;
        //kill*=1.0e-3;//???????
      } else {
        kill=killTolerance;
        //kill*=1.0e-3;//???????
      } 
      kill*=1.0e-3;//be conservative
      double smallerSlack=COIN_DBL_MAX;
      double largerzw;
      if (zValue>tValue) {
        largerzw=zValue;
      } else {
        largerzw=tValue;
      } 
      bool fakeOldBounds=false;
      bool fakeNewBounds=false;
      double trueLower;
      double trueUpper;
      if (iColumn<numberColumns_) {
	trueLower = columnLower_[iColumn];
	trueUpper = columnUpper_[iColumn];
      } else {
	trueLower = rowLower_[iColumn-numberColumns_];
	trueUpper = rowUpper_[iColumn-numberColumns_];
      }
      if (oldPrimal>trueLower+largeGap2&&
	  oldPrimal<trueUpper-largeGap2)
	fakeOldBounds=true;
      if (newPrimal>trueLower+largeGap2&&
	  newPrimal<trueUpper-largeGap2)
	fakeNewBounds=true;
      if (fakeOldBounds) {
	if (fakeNewBounds) {
	  lower_[iColumn]=newPrimal-largeGap2;
	  lowerSlack[iColumn] = largeGap2;
	  upper_[iColumn]=newPrimal+largeGap2;
	  upperSlack[iColumn] = largeGap2;
	} else {
	  lower_[iColumn]=trueLower;
	  setLowerBound(iColumn);
	  lowerSlack[iColumn] = max(newPrimal-trueLower,1.0);
	  upper_[iColumn]=trueUpper;
	  setUpperBound(iColumn);
	  upperSlack[iColumn] = max(trueUpper-newPrimal,1.0);
	}
      } else if (fakeNewBounds) {
	lower_[iColumn]=newPrimal-largeGap2;
	lowerSlack[iColumn] = largeGap2;
	upper_[iColumn]=newPrimal+largeGap2;
	upperSlack[iColumn] = largeGap2;
	// so we can just have one test
	fakeOldBounds=true;
      }
      if (lowerBound(iColumn)) {
        double oldSlack = lowerSlack[iColumn];
        double newSlack;
	newSlack=
	  lowerSlack[iColumn]+actualPrimalStep_*(oldPrimal-oldSlack
						 + thisWeight-lower[iColumn]);
	if (fakeOldBounds)
	  newSlack = lowerSlack[iColumn];
        double epsilon = fabs(newSlack)*epsilonBase;
        if (epsilon>1.0e-5) {
          //cout<<"bad"<<endl;
          epsilon=1.0e-5;
        } 
        //for changing slack
        double zValue2 = zValue;
        if (zValue2<epsilonBase) {
          zValue2=epsilonBase;
        } 
        //make sure reasonable
        if (zValue<epsilon) {
          zValue=epsilon;
        } 
        //store modified zVec
        //no zVec[iColumn]=zValue;
        double feasibleSlack=newPrimal-lower[iColumn];
        if (feasibleSlack>0.0&&newSlack>0.0) {
          double smallGap2=smallGap;
          if (fabs(0.1*newPrimal)>smallGap) {
            smallGap2=0.1*fabs(newPrimal);
          } 
	  double larger;
	  if (newSlack>feasibleSlack) {
	    larger=newSlack;
	  } else {
	    larger=feasibleSlack;
	  } 
	  if (fabs(feasibleSlack-newSlack)<1.0e-6*larger) {
	    newSlack=feasibleSlack;
	  } 
        } 
        if (zVec[iColumn]>dualTolerance) {
          dualObjectiveThis+=lower[iColumn]*zVec[iColumn];
        } 
        lowerSlack[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
	double infeasibility = fabs(newPrimal-lowerSlack[iColumn]-lower[iColumn]);
	if (infeasibility>maximumBoundInfeasibility) {
	  maximumBoundInfeasibility=infeasibility;
	} 
        if (lowerSlack[iColumn]<=1.0e5*kill&&fabs(newPrimal-lower[iColumn])<=1.0e5*kill) {
	  double step = min(actualPrimalStep_*1.1,1.0);
	  double newPrimal2=primal[iColumn]+step*thisWeight;
	  if (newPrimal2<newPrimal&&dj_[iColumn]>1.0e-5&&numberIterations_>50-40) {
	    newPrimal=lower[iColumn];
	    lowerSlack[iColumn]=0.0;
	    printf("fixing %d to lower\n",iColumn);
	  }
	}
        if (lowerSlack[iColumn]<=kill&&fabs(newPrimal-lower[iColumn])<=kill) {
          //may be better to leave at value?
          newPrimal=lower[iColumn];
          lowerSlack[iColumn]=0.0;
          thisKilled=true;
          //cout<<j<<" l "<<reducedCost<<" "<<zVec[iColumn]<<endl;
        } else {
          s+=lowerSlack[iColumn];
        } 
      } 
      if (upperBound(iColumn)) {
        //reducedCost-=perturbation;
        //primalObjectiveThis-=perturbation*newPrimal;
        double oldSlack = upperSlack[iColumn];
        double newSlack;
	newSlack=
	  upperSlack[iColumn]+actualPrimalStep_*(-oldPrimal-oldSlack
						 - thisWeight+upper[iColumn]);
	if (fakeOldBounds)
	  newSlack = upperSlack[iColumn];
        double epsilon = fabs(newSlack)*epsilonBase;
        if (epsilon>1.0e-5) {
          //cout<<"bad"<<endl;
          epsilon=1.0e-5;
        } 
        //for changing slack
        double tValue2 = tValue;
        if (tValue2<epsilonBase) {
          tValue2=epsilonBase;
        } 
        //make sure reasonable
        if (tValue<epsilon) {
          tValue=epsilon;
        } 
        //store modified tVec
        //no tVec[iColumn]=tValue;
        double feasibleSlack=upper[iColumn]-newPrimal;
        if (feasibleSlack>0.0&&newSlack>0.0) {
          double smallGap2=smallGap;
          if (fabs(0.1*newPrimal)>smallGap) {
            smallGap2=0.1*fabs(newPrimal);
          } 
	  double larger;
	  if (newSlack>feasibleSlack) {
	    larger=newSlack;
	  } else {
	    larger=feasibleSlack;
	  } 
	  if (fabs(feasibleSlack-newSlack)<1.0e-6*larger) {
	    newSlack=feasibleSlack;
	  } 
        } 
        if (tVec[iColumn]>dualTolerance) {
          dualObjectiveThis-=upper[iColumn]*tVec[iColumn];
        } 
        upperSlack[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
	double infeasibility = fabs(newPrimal+upperSlack[iColumn]-upper[iColumn]);
	if (infeasibility>maximumBoundInfeasibility) {
	  maximumBoundInfeasibility=infeasibility;
        } 
        if (upperSlack[iColumn]<=1.0e5*kill&&fabs(newPrimal-upper[iColumn])<=1.0e5*kill) {
	  double step = min(actualPrimalStep_*1.1,1.0);
	  double newPrimal2=primal[iColumn]+step*thisWeight;
	  if (newPrimal2>newPrimal&&dj_[iColumn]<-1.0e-5&&numberIterations_>50-40) {
	    newPrimal=upper[iColumn];
	    upperSlack[iColumn]=0.0;
	    printf("fixing %d to upper\n",iColumn);
	  }
	}
        if (upperSlack[iColumn]<=kill&&fabs(newPrimal-upper[iColumn])<=kill) {
          //may be better to leave at value?
          newPrimal=upper[iColumn];
          upperSlack[iColumn]=0.0;
          thisKilled=true;
        } else {
          t+=upperSlack[iColumn];
        } 
      } 
      primal[iColumn]=newPrimal;
      if (fabs(newPrimal)>solutionNorm) {
        solutionNorm=fabs(newPrimal);
      } 
      if (!thisKilled) {
        primalObjectiveValue+=primalObjectiveThis;
        dualObjectiveValue+=dualObjectiveThis;
        if (s>largeGap) {
          s=largeGap;
        } 
        if (t>largeGap) {
          t=largeGap;
        } 
        double divisor = s*tValue+t*zValue;
        double diagonalValue=(t*s)/divisor;
        diagonal[iColumn]=diagonalValue;
        //FUDGE
        if (diagonalValue>diagonalLimit) {
          //cout<<"large diagonal "<<diagonalValue<<endl;
          diagonal[iColumn]=diagonalLimit;
        } 
        if (diagonalValue<1.0e-10) {
          //cout<<"small diagonal "<<diagonalValue<<endl;
        } 
        if (diagonalValue>largestDiagonal) {
          largestDiagonal=diagonalValue;
        } 
        if (diagonalValue<smallestDiagonal) {
          smallestDiagonal=diagonalValue;
        } 
        double dualInfeasibility=reducedCost-zVec[iColumn]+tVec[iColumn];
        if (fabs(dualInfeasibility)>dualTolerance) {
	  dualFake+=newPrimal*dualInfeasibility;
        } 
        dualInfeasibility=fabs(dualInfeasibility);
        if (dualInfeasibility>maximumDualError) {
          maximumDualError=dualInfeasibility;
        } 
        deltaZ[iColumn]=0.0;
      } else {
        numberKilled++;
        diagonal[iColumn]=0.0;
        zVec[iColumn]=0.0;
        tVec[iColumn]=0.0;
        setFlagged(iColumn);
        setFixedOrFree(iColumn);
        deltaZ[iColumn]=newPrimal;
        offsetObjective+=newPrimal*cost[iColumn];
      } 
    } else {
      deltaZ[iColumn]=primal[iColumn];
      diagonal[iColumn]=0.0;
      offsetObjective+=primal[iColumn]*cost[iColumn];
      if (upper[iColumn]-lower[iColumn]>1.0e-5) {
        if (primal[iColumn]<lower[iColumn]+1.0e-8&&dual[iColumn]<-1.0e-8) {
          if (-dual[iColumn]>maximumDJInfeasibility) {
            maximumDJInfeasibility=-dual[iColumn];
          } 
        } 
        if (primal[iColumn]>upper[iColumn]-1.0e-8&&dual[iColumn]>1.0e-8) {
          if (dual[iColumn]>maximumDJInfeasibility) {
            maximumDJInfeasibility=dual[iColumn];
          } 
        } 
      } 
    } 
  } 
  handler_->message(CLP_BARRIER_DIAGONAL,messages_)
    <<largestDiagonal<<smallestDiagonal
    <<CoinMessageEol;
#if 0
  // If diagonal wild - kill some
  if (largestDiagonal>1.0e17*smallestDiagonal) {
    double killValue =largestDiagonal*1.0e-17;
    for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (fabs(diagonal_[iColumn])<killValue)
	diagonal_[iolumn]=0.0;
    }
  }
#endif
  // update rhs region
  multiplyAdd(deltaZ+numberColumns_,numberRows_,-1.0,rhsFixRegion_,1.0);
  matrix_->times(1.0,deltaZ,rhsFixRegion_);
  primalObjectiveValue=primalObjectiveValue+offsetObjective;
  dualObjectiveValue+=offsetObjective+dualFake;
  if (numberIncreased||numberDecreased) {
    handler_->message(CLP_BARRIER_SLACKS,messages_)
      <<numberIncreased<<numberDecreased
      <<CoinMessageEol;
  } 
  if (maximumDJInfeasibility) {
    handler_->message(CLP_BARRIER_DUALINF,messages_)
      <<maximumDJInfeasibility
      <<CoinMessageEol;
  } 
  // Need to rethink (but it is only for printing)
  sumPrimalInfeasibilities_=maximumRhsInfeasibility;
  sumDualInfeasibilities_=maximumDualError;
  maximumBoundInfeasibility_ = maximumBoundInfeasibility;
  //compute error and fixed RHS
  multiplyAdd(solution_+numberColumns_,numberRows_,-1.0,errorRegion_,0.0);
  matrix_->times(1.0,solution_,errorRegion_);
  maximumDualError_=maximumDualError;
  maximumBoundInfeasibility_=maximumBoundInfeasibility;
  solutionNorm_=solutionNorm;
  //finish off objective computation
  primalObjective_=primalObjectiveValue*scaleFactor_;
  double dualValue2=innerProduct(dual_,numberRows_,
                rhsFixRegion_);
  dualObjectiveValue-=dualValue2;
  dualObjective_=dualObjectiveValue*scaleFactor_;
  if (numberKilled) {
      handler_->message(CLP_BARRIER_KILLED,messages_)
	<<numberKilled
	<<CoinMessageEol;
  } 
  double maximumRHSError1=0.0;
  double maximumRHSError2=0.0;
  double primalOffset=0.0;
  char * dropped = cholesky_->rowsDropped();
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double value=errorRegion_[iRow];
    if (!dropped[iRow]) {
      if (fabs(value)>maximumRHSError1) {
	maximumRHSError1=fabs(value);
      } 
    } else {
      if (fabs(value)>maximumRHSError2) {
	maximumRHSError2=fabs(value);
      } 
      primalOffset+=value*dual_[iRow];
    } 
  } 
  primalObjective_-=primalOffset*scaleFactor_;
  if (maximumRHSError1>maximumRHSError2) {
    maximumRHSError_=maximumRHSError1;
  } else {
    maximumRHSError_=maximumRHSError1; //note change
    if (maximumRHSError2>primalTolerance()) {
      handler_->message(CLP_BARRIER_ABS_DROPPED,messages_)
	<<maximumRHSError2
	<<CoinMessageEol;
    } 
    //if (maximumRHSError2>1.0e-1) {
      //set all bits false
      //cholesky_->resetRowsDropped();
    //}
  } 
  objectiveNorm_=maximumAbsElement(dual_,numberRows_);
  if (objectiveNorm_<1.0e-12) {
    objectiveNorm_=1.0e-12;
  } 
  if (objectiveNorm_<baseObjectiveNorm_) {
    //std::cout<<" base "<<baseObjectiveNorm_<<" "<<objectiveNorm_<<std::endl;
    if (objectiveNorm_<baseObjectiveNorm_*1.0e-4) {
      objectiveNorm_=baseObjectiveNorm_*1.0e-4;
    } 
  } 
  bool primalFeasible=true;
  if (maximumRHSError_>primalTolerance()||
      maximumDualError_>dualTolerance/scaleFactor_) {
    handler_->message(CLP_BARRIER_ABS_ERROR,messages_)
      <<maximumRHSError_<<maximumDualError_
      <<CoinMessageEol;
  } 
  if (rhsNorm_>solutionNorm_) {
    solutionNorm_=rhsNorm_;
  } 
  double scaledRHSError=maximumRHSError_/solutionNorm_;
  bool dualFeasible=true;
  if (maximumBoundInfeasibility_>primalTolerance()||
      scaledRHSError>primalTolerance())
    primalFeasible=false;
  if (maximumDualError_>objectiveNorm_*dualTolerance) 
    dualFeasible=false;
  if (!primalFeasible||!dualFeasible) {
    handler_->message(CLP_BARRIER_FEASIBLE,messages_)
      <<maximumBoundInfeasibility_<<scaledRHSError
      <<maximumDualError_/objectiveNorm_
      <<CoinMessageEol;
  } 
  if (!gonePrimalFeasible_) {
    gonePrimalFeasible_=primalFeasible;
  } else if (!primalFeasible) {
    gonePrimalFeasible_=primalFeasible;
    if (!numberKilled) {
      handler_->message(CLP_BARRIER_GONE_INFEASIBLE,messages_)
      <<CoinMessageEol;
    } 
  } 
  if (!goneDualFeasible_) {
    goneDualFeasible_=dualFeasible;
  } else if (!dualFeasible) {
    handler_->message(CLP_BARRIER_GONE_INFEASIBLE,messages_)
      <<CoinMessageEol;
    goneDualFeasible_=dualFeasible;
  } 
  //objectiveValue();
  if (solutionNorm_>1.0e40) {
    std::cout <<"primal off to infinity"<<std::endl;
    abort();
  } 
  if (objectiveNorm_>1.0e40) {
    std::cout <<"dual off to infinity"<<std::endl;
    abort();
  } 
  handler_->message(CLP_BARRIER_STEP,messages_)
    <<actualPrimalStep_
    <<actualDualStep_
    <<mu_
    <<CoinMessageEol;
  numberIterations_++;
  return numberKilled;
}
//  Save info on products of affine deltaT*deltaW and deltaS*deltaZ
double 
ClpPredictorCorrector::affineProduct()
{
  double * deltaZ = deltaZ_;
  double * deltaT = deltaT_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * primal = solution_;
  //direction vector in deltaX
  double * deltaX = deltaX_;
  double product = 0.0;
  //IF zVec starts as 0 then deltaZ always zero
  //(remember if free then zVec not 0)
  //I think free can be done with careful use of boundSlacks to zero
  //out all we want
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    double w3=deltaZ[iColumn]*deltaX[iColumn];
    double w4=-deltaT[iColumn]*deltaX[iColumn];
    if (lowerBound(iColumn)) {
      w3+=deltaZ[iColumn]*(primal[iColumn]-lowerSlack[iColumn]-lower[iColumn]);
      product+=w3;
    } 
    if (upperBound(iColumn)) {
      w4+=deltaT[iColumn]*(-primal[iColumn]-upperSlack[iColumn]+upper[iColumn]);
      product+=w4;
    } 
  } 
  return product;
}
