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
#include "ClpQuadraticObjective.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <cstdio>
#include <iostream>
#if 0
static int yyyyyy=0;
void ClpPredictorCorrector::saveSolution(std::string fileName)
{
  FILE * fp=fopen(fileName.c_str(),"wb");
  if (fp) {
    int numberRows=numberRows_;
    int numberColumns=numberColumns_;
    fwrite(&numberRows,sizeof(int),1,fp);
    fwrite(&numberColumns,sizeof(int),1,fp);
    double dsave[20];
    memset(dsave,0,sizeof(dsave));
    fwrite(dsave,sizeof(double),20,fp);
    int msave[20];
    memset(msave,0,sizeof(msave));
    msave[0]=numberIterations_;
    fwrite(msave,sizeof(int),20,fp);
    fwrite(dual_,sizeof(double),numberRows,fp);
    fwrite(errorRegion_,sizeof(double),numberRows,fp);
    fwrite(rhsFixRegion_,sizeof(double),numberRows,fp);
    fwrite(solution_,sizeof(double),numberColumns,fp);
    fwrite(solution_+numberColumns,sizeof(double),numberRows,fp);
    fwrite(diagonal_,sizeof(double),numberColumns,fp);
    fwrite(diagonal_+numberColumns,sizeof(double),numberRows,fp);
    fwrite(wVec_,sizeof(double),numberColumns,fp);
    fwrite(wVec_+numberColumns,sizeof(double),numberRows,fp);
    fwrite(zVec_,sizeof(double),numberColumns,fp);
    fwrite(zVec_+numberColumns,sizeof(double),numberRows,fp);
    fwrite(upperSlack_,sizeof(double),numberColumns,fp);
    fwrite(upperSlack_+numberColumns,sizeof(double),numberRows,fp);
    fwrite(lowerSlack_,sizeof(double),numberColumns,fp);
    fwrite(lowerSlack_+numberColumns,sizeof(double),numberRows,fp);
    fclose(fp);
  } else {
    std::cout<<"Unable to open file "<<fileName<<std::endl;
  }
}
#endif
static double eScale=1.0e27;
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
  //diagonalPerturbation_=1.0e-25;
  ClpMatrixBase * saveMatrix = NULL;
  // If quadratic then make copy so we can actually scale or normalize
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  /* If modeSwitch is :
     0 - normal
     1 - bit switch off centering
     2 - bit always do type 2
     4 - accept corrector nearly always
  */
  int modeSwitch=0;
  //if (quadraticObj)
    //modeSwitch |= 1; // switch off centring for now
  //if (quadraticObj)
  //modeSwitch |=4;
  ClpObjective * saveObjective = NULL;
  if (quadraticObj) {
    // check KKT is on
    if (!cholesky_->kkt()) {
      //No!
      handler_->message(CLP_BARRIER_KKT,messages_)
	<<CoinMessageEol;
      return -1;
    }
    saveObjective=objective_;
    // We are going to make matrix full rather half
    objective_ = new ClpQuadraticObjective(*quadraticObj,1);
  }
  bool allowIncreasingGap = (modeSwitch&4)!=0;
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
  int returnCode = cholesky_->order(this);
  if (returnCode||cholesky_->symbolic()) {
    printf("Error return from symbolic - probably not enough memory\n");
    problemStatus_=4;
    //delete all temporary regions
    deleteWorkingData();
    if (saveMatrix) {
      // restore normal copy
      delete matrix_;
      matrix_ = saveMatrix;
    }
    // Restore quadratic objective if necessary
    if (saveObjective) {
      delete objective_;
      objective_=saveObjective;
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
  double lastComplementarityGap=COIN_DBL_MAX*1.0e20;
  double lastStep=1.0;
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
  double * saveW = new double[numberTotal];
  double * saveSL = new double[numberTotal];
  double * saveSU = new double[numberTotal];
  // Save smallest mu used in primal dual moves
  double smallestPrimalDualMu=COIN_DBL_MAX;
  double objScale = optimizationDirection_/
    (rhsScale_*objectiveScale_);
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
	printf(" %d %g %g %g %g %g %g\n",i,solution_[i],diagonal_[i],wVec_[i],
	       zVec_[i],upperSlack_[i],lowerSlack_[i]);
      }
    }
#endif
    complementarityGap_=complementarityGap(numberComplementarityPairs_,
					   numberComplementarityItems_,0);
    handler_->message(CLP_BARRIER_ITERATION,messages_)
      <<numberIterations_
      <<primalObjective_*objScale- dblParam_[ClpObjOffset]
      << dualObjective_*objScale- dblParam_[ClpObjOffset]
      <<complementarityGap_
      <<numberFixedTotal
      <<cholesky_->rank()
      <<CoinMessageEol;
#if 0
    if (numberIterations_==-1) {
      saveSolution("xxx.sav");
      if (yyyyyy)
        exit(99);
    }
#endif
    // move up history
    for (i=1;i<LENGTH_HISTORY;i++) 
      historyInfeasibility_[i-1]=historyInfeasibility_[i];
    historyInfeasibility_[LENGTH_HISTORY-1]=complementarityGap_;
    // switch off saved if changes
    if (saveIteration+10<numberIterations_&&
	complementarityGap_*2.0<historyInfeasibility_[0])
      saveIteration=-1;
    lastStep = CoinMin(actualPrimalStep_,actualDualStep_);
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
	    problemStatus_=0; // benefit of doubt
	    break;
	  }
	} else {
	  // not close to optimal but check if getting bad
	  double scaledRHSError=maximumRHSError_/solutionNorm_;
	  if ((maximumBoundInfeasibility_>1.0e-1||
	      scaledRHSError>1.0e-1||
	       maximumDualError_>objectiveNorm_*1.0e-1)
	      &&(numberIterations_>50
		 &&complementarityGap_>0.9*historyInfeasibility_[0])) {
	    handler_->message(CLP_BARRIER_EXIT2,messages_)
	      <<saveIteration
	      <<CoinMessageEol;
	    break;
	  }
	  if (complementarityGap_>0.95*checkGap&&bestObjectiveGap<1.0e-3&&
	      (numberIterations_>saveIteration+5||numberIterations_>100)) {
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
    int numberBack = quadraticObj ? 10 : 5;
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
      } else if (numberIterations_-lastGoodIteration>=numberBack&&
		 complementarityGap_<1.0e-6) {
	break; // not doing very well - give up
      } 
    } else if (complementarityGap_<goodGapChange*lastComplementarityGap) {
      lastGoodIteration=numberIterations_;
      lastComplementarityGap=complementarityGap_;
    } else if (numberIterations_-lastGoodIteration>=numberBack&&
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
    if (numberIterations_>maximumBarrierIterations_||hitMaximumIterations()) {
      handler_->message(CLP_BARRIER_STOPPING,messages_)
	<<CoinMessageEol;
      problemStatus_=3;
      onStopped(); // set secondary status
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
#ifndef NDEBUG
	//int newDropped2=cholesky_->factorize(diagonal_,rowsDroppedThisTime);
	//assert(!newDropped2);
#endif
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
    if ((modeSwitch&2)==0) {
      directionAccuracy=findDirectionVector(phase);
      if (directionAccuracy>worstDirectionAccuracy_) {
	worstDirectionAccuracy_=directionAccuracy;
      }
      if (saveIteration>0&&directionAccuracy>1.0) {
	handler_->message(CLP_BARRIER_EXIT2,messages_)
	  <<saveIteration
	  <<CoinMessageEol;
	break;
      }
      findStepLength(phase);
      nextGap=complementarityGap(nextNumber,nextNumberItems,1);
      debugMove(0,actualPrimalStep_,actualDualStep_);
      debugMove(0,1.0e-2,1.0e-2);
    }
    double affineGap=nextGap;
    int bestPhase=0;
    double bestNextGap=nextGap;
    // ?
    bestNextGap=CoinMax(nextGap,0.8*complementarityGap_);
    if (quadraticObj)
      bestNextGap=CoinMax(nextGap,0.99*complementarityGap_);
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
      if (numberComplementarityPairs_<=5000) {
	phi=pow(static_cast<double> (numberComplementarityPairs_),2.0);
      } else {
	phi=pow(static_cast<double> (numberComplementarityPairs_),1.5);
	if (phi<500.0*500.0) {
	  phi=500.0*500.0;
	} 
      }
      mu_=complementarityGap_/phi;
    } 
    //save information
    double product=affineProduct();
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
    if (complementarityGap_*(beta2-tau)+product-mu_*numberComplementarityPairs_<0.0&&0) {
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
    // If bad gap - try standard primal dual
    if (nextGap>complementarityGap_*1.001)
      goodMove=false;
    if ((modeSwitch&2)!=0)
      goodMove=false;
    if (goodMove&&doCorrector) {
      CoinMemcpyN(deltaX_,numberTotal,saveX);
      CoinMemcpyN(deltaY_,numberRows_,saveY);
      CoinMemcpyN(deltaZ_,numberTotal,saveZ);
      CoinMemcpyN(deltaW_,numberTotal,saveW);
      CoinMemcpyN(deltaSL_,numberTotal,saveSL);
      CoinMemcpyN(deltaSU_,numberTotal,saveSU);
#ifdef HALVE
      double savePrimalStep = actualPrimalStep_;
      double saveDualStep = actualDualStep_;
      double saveMu = mu_;
#endif
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
      } 
      if (goodMove) {
	phase=1;
	double norm = findStepLength(phase);
	nextGap = complementarityGap(nextNumber,nextNumberItems,1);
	debugMove(1,actualPrimalStep_,actualDualStep_);
	//debugMove(1,1.0e-7,1.0e-7);
	goodMove=checkGoodMove(true,bestNextGap,allowIncreasingGap);
	if (norm<0)
	  goodMove=false;
	if (!goodMove) {
#ifdef SOME_DEBUG
	  printf("checkGoodMove failed\n");
#endif
	}
      }
#ifdef HALVE
      int nHalve=0;
      // relax test
      bestNextGap=CoinMax(bestNextGap,0.9*complementarityGap_);
      while (!goodMove) {
	mu_=saveMu;
	actualPrimalStep_ = savePrimalStep;
	actualDualStep_ = saveDualStep;
	int i;
	//printf("halve %d\n",nHalve);
	nHalve++;
	const double lambda=0.5;
	for (i=0;i<numberRows_;i++)
	  deltaY_[i] = lambda*deltaY_[i]+(1.0-lambda)*saveY[i];
	for (i=0;i<numberTotal;i++) {
	  deltaX_[i] = lambda*deltaX_[i]+(1.0-lambda)*saveX[i];
	  deltaZ_[i] = lambda*deltaZ_[i]+(1.0-lambda)*saveZ[i];
	  deltaW_[i] = lambda*deltaW_[i]+(1.0-lambda)*saveW[i];
	  deltaSL_[i] = lambda*deltaSL_[i]+(1.0-lambda)*saveSL[i];
	  deltaSU_[i] = lambda*deltaSU_[i]+(1.0-lambda)*saveSU[i];
	}
   CoinMemcpyN(saveX,numberTotal,deltaX_);
   CoinMemcpyN(saveY,numberRows_,deltaY_);
   CoinMemcpyN(saveZ,numberTotal,deltaZ_);
   CoinMemcpyN(saveW,numberTotal,deltaW_);
   CoinMemcpyN(saveSL,numberTotal,deltaSL_);
   CoinMemcpyN(saveSU,numberTotal,deltaSU_);
	findStepLength(1);
	nextGap = complementarityGap(nextNumber,nextNumberItems,1);
	goodMove=checkGoodMove(true,bestNextGap,allowIncreasingGap);
	if (nHalve>10)
	  break;
	//assert (goodMove);
      }
      if (nHalve&&handler_->logLevel()>1) 
	printf("halved %d times\n",nHalve);
#endif
    }
    //bestPhase=-1;
    //goodMove=false;
    if (!goodMove) {
      // Just primal dual step
      double floatNumber;
      floatNumber = 2.0*numberComplementarityPairs_;
      //floatNumber = numberComplementarityItems_;
      double saveMu=mu_; // use one from predictor corrector
      mu_=complementarityGap_/floatNumber;
      // If going well try small mu
      mu_ *= sqrt((1.0-lastStep)/(1.0+10.0*lastStep));
      double mu1=mu_;
      double phi;
      if (numberComplementarityPairs_<=500) {
	phi=pow(static_cast<double> (numberComplementarityPairs_),2.0);
      } else {
	phi=pow(static_cast<double> (numberComplementarityPairs_),1.5);
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
      //mu_=CoinMin(smallestPrimalDualMu*0.95,mu_);
      smallestPrimalDualMu = mu_;
      // Try simpler
      floatNumber = numberComplementarityItems_;
      mu_=0.5*complementarityGap_/floatNumber;
      //if ((modeSwitch&2)==0) {
      //if ((numberIterations_&1)==0)
      //  mu_ *= 0.5;
      //} else {
      //mu_ *= 0.8;
      //}
      //set up for next step
      setupForSolve(2);
      findDirectionVector(2);
      double norm=findStepLength(2);
      // just for debug
      bestNextGap = complementarityGap_*1.0005;
      //bestNextGap=COIN_DBL_MAX;
      nextGap=complementarityGap(nextNumber,nextNumberItems,2);
      debugMove(2,actualPrimalStep_,actualDualStep_);
      //debugMove(2,1.0e-7,1.0e-7);
      checkGoodMove(false, bestNextGap,allowIncreasingGap);
      if ((nextGap>0.9*complementarityGap_&&bestPhase==0&&affineGap<nextGap
	  &&(numberIterations_>80||(numberIterations_>20&&quadraticObj)))||norm<0.0) {
	// Back to affine
	phase=0;
	setupForSolve(phase);
	directionAccuracy=findDirectionVector(phase);
	findStepLength(phase);
	nextGap=complementarityGap(nextNumber,nextNumberItems,1);
	bestNextGap = complementarityGap_;
	//checkGoodMove(false, bestNextGap,allowIncreasingGap);
      }
    }
    if (numberIterations_==0)
      smallestPrimalDualMu=mu_;
    if (!goodMove)
      mu_=nextGap / (static_cast<double> (nextNumber)*1.1);
    //if (quadraticObj)
    //goodMove=true; 
    //goodMove=false; //TEMP
    // Do centering steps
    int numberTries=0;
    double nextCenterGap=0.0;
    int numberGoodTries=0;
    double originalDualStep=actualDualStep_;
    double originalPrimalStep=actualPrimalStep_;
    if (actualDualStep_>0.9&&actualPrimalStep_>0.9)
      goodMove=false; // don't bother
    if ((modeSwitch&1)!=0)
      goodMove=false;
    while (goodMove&&numberTries<5) {
      goodMove=false;
      numberTries++;
      CoinMemcpyN(deltaX_,numberTotal,saveX);
      CoinMemcpyN(deltaY_,numberRows_,saveY);
      CoinMemcpyN(deltaZ_,numberTotal,saveZ);
      CoinMemcpyN(deltaW_,numberTotal,saveW);
      double savePrimalStep = actualPrimalStep_;
      double saveDualStep = actualDualStep_;
      double saveMu = mu_;
      setupForSolve(3);
      findDirectionVector(3);
      findStepLength(3);
      debugMove(3,actualPrimalStep_,actualDualStep_);
      //debugMove(3,1.0e-7,1.0e-7);
      double xGap = complementarityGap(nextNumber,nextNumberItems,3);
      // If one small then that's the one that counts
      double checkDual=saveDualStep;
      double checkPrimal=savePrimalStep;
      if (checkDual>5.0*checkPrimal) {
	checkDual=2.0*checkPrimal;
      } else if (checkPrimal>5.0*checkDual) {
	checkPrimal=2.0*checkDual;
      }
      if (actualPrimalStep_<checkPrimal||
	  actualDualStep_<checkDual||
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
 CoinMemcpyN(saveX,numberTotal,deltaX_);
 CoinMemcpyN(saveY,numberRows_,deltaY_);
 CoinMemcpyN(saveZ,numberTotal,deltaZ_);
 CoinMemcpyN(saveW,numberTotal,deltaW_);
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
  delete [] saveW;
  delete [] saveSL;
  delete [] saveSU;
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
  checkSolution();
  CoinMemcpyN(reducedCost_,numberColumns_,dj_);
  // If quadratic use last solution 
  // Restore quadratic objective if necessary
  if (saveObjective) {
    delete objective_;
    objective_=saveObjective;
    objectiveValue_=0.5*(primalObjective_+dualObjective_);
  }
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
  double maximumPrimalStep=COIN_DBL_MAX*1.0e-20;
  double maximumDualStep=COIN_DBL_MAX;
  int numberTotal = numberRows_+numberColumns_;
  double tolerance = 1.0e-12;
  int chosenPrimalSequence=-1;
  int chosenDualSequence=-1;
  bool lowPrimal=false;
  bool lowDual=false;
  // If done many iterations then allow to hit boundary
  double hitTolerance;
  //printf("objective norm %g\n",objectiveNorm_);
  if (numberIterations_<80||!gonePrimalFeasible_)
    hitTolerance = COIN_DBL_MAX;
  else
    hitTolerance = CoinMax(1.0e3,1.0e-3*objectiveNorm_);
  int iColumn;
  //printf("dual value %g\n",dual_[0]);
  //printf("     X     dX      lS     dlS     uS     dUs    dj    Z dZ     t   dT\n"); 
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double directionElement=deltaX_[iColumn];
      if (directionNorm<fabs(directionElement)) {
	directionNorm=fabs(directionElement);
      }
      if (0) 
      printf("%d %g %g %g %g %g %g %g %g %g %g %g\n",
	     iColumn,solution_[iColumn],deltaX_[iColumn],
	     lowerSlack_[iColumn],deltaSL_[iColumn],
	     upperSlack_[iColumn],deltaSU_[iColumn],
	     dj_[iColumn],
	     zVec_[iColumn],deltaZ_[iColumn],
	     wVec_[iColumn],deltaW_[iColumn]);
      if (lowerBound(iColumn)) {
	double delta = - deltaSL_[iColumn];
	double z1 = deltaZ_[iColumn];
	double newZ = zVec_[iColumn]+z1;
	if (zVec_[iColumn]>tolerance) {
	  if (zVec_[iColumn]<-z1*maximumDualStep) {
	    maximumDualStep=-zVec_[iColumn]/z1;
	    chosenDualSequence=iColumn;
	    lowDual=true;
	  } 
	} 
	if (lowerSlack_[iColumn]<maximumPrimalStep*delta) {
	  double newStep=lowerSlack_[iColumn]/delta;
	  if (newStep>0.2||newZ<hitTolerance||delta>1.0e3||delta<=1.0e-6||dj_[iColumn]<hitTolerance) {
	    maximumPrimalStep = newStep;
	    chosenPrimalSequence=iColumn;
	    lowPrimal=true;
	  } else {
	    //printf("small %d delta %g newZ %g step %g\n",iColumn,delta,newZ,newStep); 
	  }
	} 
      }
      if (upperBound(iColumn)) {
	double delta = - deltaSU_[iColumn];;
	double w1 = deltaW_[iColumn];
	double newT = wVec_[iColumn]+w1;
	if (wVec_[iColumn]>tolerance) {
	  if (wVec_[iColumn]<-w1*maximumDualStep) {
	    maximumDualStep=-wVec_[iColumn]/w1;
	    chosenDualSequence=iColumn;
	    lowDual=false;
	  } 
	} 
	if (upperSlack_[iColumn]<maximumPrimalStep*delta) {
	  double newStep=upperSlack_[iColumn]/delta;
	  if (newStep>0.2||newT<hitTolerance||delta>1.0e3||delta<=1.0e-6||dj_[iColumn]>-hitTolerance) {
	    maximumPrimalStep = newStep;
	    chosenPrimalSequence=iColumn;
	    lowPrimal=false;
	  } else {
	    //printf("small %d delta %g newT %g step %g\n",iColumn,delta,newT,newStep); 
	  }
	} 
      } 
    } 
  }
#ifdef SOME_DEBUG
  printf("new step - phase %d, norm %.18g, dual step %.18g, primal step %.18g\n",
	 phase,directionNorm,maximumDualStep,maximumPrimalStep);
  if (lowDual) 
    printf("ld %d %g %g => %g (dj %g,sol %g) ",
	   chosenDualSequence,zVec_[chosenDualSequence],
	   deltaZ_[chosenDualSequence],zVec_[chosenDualSequence]+
	   maximumDualStep*deltaZ_[chosenDualSequence],dj_[chosenDualSequence],
	   solution_[chosenDualSequence]);
  else
    printf("ud %d %g %g => %g (dj %g,sol %g) ",
	   chosenDualSequence,wVec_[chosenDualSequence],
	   deltaW_[chosenDualSequence],wVec_[chosenDualSequence]+
	   maximumDualStep*deltaW_[chosenDualSequence],dj_[chosenDualSequence],
	   solution_[chosenDualSequence]);
  if (lowPrimal) 
    printf("lp %d %g %g => %g (dj %g,sol %g)\n",
	   chosenPrimalSequence,lowerSlack_[chosenPrimalSequence],
	   deltaSL_[chosenPrimalSequence],lowerSlack_[chosenPrimalSequence]+
	   maximumPrimalStep*deltaSL_[chosenPrimalSequence],
	   dj_[chosenPrimalSequence],solution_[chosenPrimalSequence]);
  else
    printf("up %d %g %g => %g (dj %g,sol %g)\n",
	   chosenPrimalSequence,upperSlack_[chosenPrimalSequence],
	   deltaSU_[chosenPrimalSequence],upperSlack_[chosenPrimalSequence]+
	   maximumPrimalStep*deltaSU_[chosenPrimalSequence],
	   dj_[chosenPrimalSequence],solution_[chosenPrimalSequence]);
#endif
  actualPrimalStep_=stepLength_*maximumPrimalStep;
  if (phase>=0&&actualPrimalStep_>1.0) {
    actualPrimalStep_=1.0;
  } 
  actualDualStep_=stepLength_*maximumDualStep;
  if (phase>=0&&actualDualStep_>1.0) {
    actualDualStep_=1.0;
  } 
  // See if quadratic objective
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (quadraticObj) {
    // Use smaller unless very small
    double smallerStep=CoinMin(actualDualStep_,actualPrimalStep_);
    if (smallerStep>0.0001) {
      actualDualStep_=smallerStep;
      actualPrimalStep_=smallerStep;
    }
  }
#define OFFQ
#ifndef OFFQ
  if (quadraticObj) {
    // Don't bother if phase 0 or 3 or large gap
    //if ((phase==1||phase==2||phase==0)&&maximumDualError_>0.1*complementarityGap_
    //&&smallerStep>0.001) {
    if ((phase==1||phase==2||phase==0||phase==3)) {
      // minimize complementarity + norm*dual inf ? primal inf
      // at first - just check better - if not 
      // Complementarity gap will be a*change*change + b*change +c
      double a=0.0;
      double b=0.0;
      double c=0.0;
      /* SQUARE of dual infeasibility will be:
	 square of dj - ......
      */
      double aq=0.0;
      double bq=0.0;
      double cq=0.0;
      double gamma2 = gamma_*gamma_; // gamma*gamma will be added to diagonal
      double * linearDjChange = new double[numberTotal];
      CoinZeroN(linearDjChange,numberColumns_);
      multiplyAdd(deltaY_,numberRows_,1.0,linearDjChange+numberColumns_,0.0);
      matrix_->transposeTimes(-1.0,deltaY_,linearDjChange);
      CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
      const int * columnQuadratic = quadratic->getIndices();
      const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
      const int * columnQuadraticLength = quadratic->getVectorLengths();
      double * quadraticElement = quadratic->getMutableElements();
      for (iColumn=0;iColumn<numberTotal;iColumn++) {
	double oldPrimal = solution_[iColumn];
	if (!flagged(iColumn)) {
	  if (lowerBound(iColumn)) {
	    double change =oldPrimal+deltaX_[iColumn]-lowerSlack_[iColumn]-lower_[iColumn];
	    c += lowerSlack_[iColumn]*zVec_[iColumn];
	    b += lowerSlack_[iColumn]*deltaZ_[iColumn]+zVec_[iColumn]*change;
	    a += deltaZ_[iColumn]*change;
	  }
	  if (upperBound(iColumn)) {
	    double change =upper_[iColumn]-oldPrimal-deltaX_[iColumn]-upperSlack_[iColumn];
	    c += upperSlack_[iColumn]*wVec_[iColumn];
	    b += upperSlack_[iColumn]*deltaW_[iColumn]+wVec_[iColumn]*change;
	    a += deltaW_[iColumn]*change;
	  } 
	  // new djs are dj_ + change*value
	  double djChange = linearDjChange[iColumn];
	  if (iColumn<numberColumns_) {
	    for (CoinBigIndex j=columnQuadraticStart[iColumn];
		 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	      int jColumn = columnQuadratic[j];
	      double changeJ = deltaX_[jColumn];
	      double elementValue = quadraticElement[j];
	      djChange += changeJ*elementValue;
	    }
	  }
	  double gammaTerm = gamma2;
	  if (primalR_) {
	    gammaTerm += primalR_[iColumn];
	  }
	  djChange += gammaTerm;
	  // and dual infeasibility
	  double oldInf = dj_[iColumn]-zVec_[iColumn]+wVec_[iColumn]+
	    gammaTerm*solution_[iColumn];
	  double changeInf = djChange-deltaZ_[iColumn]+deltaW_[iColumn];
	  cq += oldInf*oldInf;
	  bq += 2.0*oldInf*changeInf;
	  aq += changeInf*changeInf;
	} else {
	  // fixed
	  if (lowerBound(iColumn)) {
	    c += lowerSlack_[iColumn]*zVec_[iColumn];
	  }
	  if (upperBound(iColumn)) {
	    c += upperSlack_[iColumn]*wVec_[iColumn];
	  } 
	  // new djs are dj_ + change*value
	  double djChange = linearDjChange[iColumn];
	  if (iColumn<numberColumns_) {
	    for (CoinBigIndex j=columnQuadraticStart[iColumn];
		 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	      int jColumn = columnQuadratic[j];
	      double changeJ = deltaX_[jColumn];
	      double elementValue = quadraticElement[j];
	      djChange += changeJ*elementValue;
	    }
	  }
	  double gammaTerm = gamma2;
	  if (primalR_) {
	    gammaTerm += primalR_[iColumn];
	  }
	  djChange += gammaTerm;
	  // and dual infeasibility
	  double oldInf = dj_[iColumn]-zVec_[iColumn]+wVec_[iColumn]+
	    gammaTerm*solution_[iColumn];
	  double changeInf = djChange-deltaZ_[iColumn]+deltaW_[iColumn];
	  cq += oldInf*oldInf;
	  bq += 2.0*oldInf*changeInf;
	  aq += changeInf*changeInf;
	}
      }
      delete [] linearDjChange;
      // ? We want to minimize complementarityGap + solutionNorm_*square of inf ??
      // maybe use inf and do line search
      // To check see if matches at current step
      double step=actualPrimalStep_;
      //Current gap + solutionNorm_ * sqrt (sum square inf)
      double multiplier = solutionNorm_;
      multiplier *= 0.01;
      multiplier=1.0;
      double currentInf =  multiplier*sqrt(cq);
      double nextInf = 	multiplier*sqrt(CoinMax(cq+step*bq+step*step*aq,0.0));
      double allowedIncrease=1.4;
#ifdef SOME_DEBUG
      printf("lin %g %g %g -> %g\n",a,b,c,
	     c+b*step+a*step*step);
      printf("quad %g %g %g -> %g\n",aq,bq,cq,
	     cq+bq*step+aq*step*step);
      debugMove(7,step,step);
      printf ("current dualInf %g, with step of %g is %g\n",
	      currentInf,step,nextInf);
#endif
      if (b>-1.0e-6) {
	if (phase!=0)
	  directionNorm=-1.0;
      }
      if ((phase==1||phase==2||phase==0||phase==3)&&nextInf>0.1*complementarityGap_&&
	  nextInf>currentInf*allowedIncrease) {
	//cq = CoinMax(cq,10.0);
	// convert to (x+q)*(x+q) = w
	double q = bq/(1.0*aq);
	double w = CoinMax(q*q + (cq/aq)*(allowedIncrease-1.0),0.0);
	w = sqrt(w);
	double stepX = w-q;
	step=stepX;
	nextInf = 
	  multiplier*sqrt(CoinMax(cq+step*bq+step*step*aq,0.0));
#ifdef SOME_DEBUG
	printf ("with step of %g dualInf is %g\n",
		step,nextInf);
#endif
	actualDualStep_=CoinMin(step,actualDualStep_);
	actualPrimalStep_=CoinMin(step,actualPrimalStep_);
      }
    }
  } else {
    // probably pointless as linear
    // minimize complementarity 
    // Complementarity gap will be a*change*change + b*change +c
    double a=0.0;
    double b=0.0;
    double c=0.0;
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      double oldPrimal = solution_[iColumn];
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change =oldPrimal+deltaX_[iColumn]-lowerSlack_[iColumn]-lower_[iColumn];
	  c += lowerSlack_[iColumn]*zVec_[iColumn];
	  b += lowerSlack_[iColumn]*deltaZ_[iColumn]+zVec_[iColumn]*change;
	  a += deltaZ_[iColumn]*change;
	}
	if (upperBound(iColumn)) {
	  double change =upper_[iColumn]-oldPrimal-deltaX_[iColumn]-upperSlack_[iColumn];
	  c += upperSlack_[iColumn]*wVec_[iColumn];
	  b += upperSlack_[iColumn]*deltaW_[iColumn]+wVec_[iColumn]*change;
	  a += deltaW_[iColumn]*change;
	} 
      } else {
	// fixed
	if (lowerBound(iColumn)) {
	  c += lowerSlack_[iColumn]*zVec_[iColumn];
	}
	if (upperBound(iColumn)) {
	  c += upperSlack_[iColumn]*wVec_[iColumn];
	} 
      }
    }
    // ? We want to minimize complementarityGap;
    // maybe use inf and do line search
    // To check see if matches at current step
    double step=CoinMin(actualPrimalStep_,actualDualStep_);
    double next = c+b*step+a*step*step;
#ifdef SOME_DEBUG
    printf("lin %g %g %g -> %g\n",a,b,c,
	   c+b*step+a*step*step);
    debugMove(7,step,step);
#endif
    if (b>-1.0e-6) {
      if (phase==0) {
#ifdef SOME_DEBUG
	printf("*** odd phase 0 direction\n");
#endif
      } else {
	directionNorm=-1.0;
      }
    }
    // and with ratio
    a=0.0;
    b=0.0;
    double ratio = actualDualStep_/actualPrimalStep_;
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      double oldPrimal = solution_[iColumn];
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change =oldPrimal+deltaX_[iColumn]-lowerSlack_[iColumn]-lower_[iColumn];
	  b += lowerSlack_[iColumn]*deltaZ_[iColumn]*ratio+zVec_[iColumn]*change;
	  a += deltaZ_[iColumn]*change*ratio;
	}
	if (upperBound(iColumn)) {
	  double change =upper_[iColumn]-oldPrimal-deltaX_[iColumn]-upperSlack_[iColumn];
	  b += upperSlack_[iColumn]*deltaW_[iColumn]*ratio+wVec_[iColumn]*change;
	  a += deltaW_[iColumn]*change*ratio;
	} 
      }
    }
    // ? We want to minimize complementarityGap;
    // maybe use inf and do line search
    // To check see if matches at current step
    step=actualPrimalStep_;
    double next2 = c+b*step+a*step*step;
    if (next2>next) {
      actualPrimalStep_=CoinMin(actualPrimalStep_,actualDualStep_);
      actualDualStep_=actualPrimalStep_;
    }
#ifdef SOME_DEBUG
    printf("linb %g %g %g -> %g\n",a,b,c,
	   c+b*step+a*step*step);
    debugMove(7,actualPrimalStep_,actualDualStep_);
#endif
    if (b>-1.0e-6) {
      if (phase==0) {
#ifdef SOME_DEBUG
	printf("*** odd phase 0 direction\n");
#endif
      } else {
	directionNorm=-1.0;
      }
    }
  }
#else
  //actualPrimalStep_ =0.5*actualDualStep_;
#endif
#ifdef FULL_DEBUG
  if (phase==3){
    double minBeta = 0.1*mu_;
    double maxBeta = 10.0*mu_;
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change = -rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack_[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (delta2Z_[iColumn]<minBeta||delta2Z_[iColumn]>maxBeta)
	    printf("3lower %d primal %g, dual %g, gap %g, old gap %g\n",
		   iColumn,primalValue,dualValue,gapProduct,delta2Z_[iColumn]);
	}  
	if (upperBound(iColumn)) {
	  double change = rhsU_[iColumn]-deltaX_[iColumn];
	  double dualValue=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	  double primalValue=upperSlack_[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  if (delta2W_[iColumn]<minBeta||delta2W_[iColumn]>maxBeta)
	    printf("3upper %d primal %g, dual %g, gap %g, old gap %g\n",
		 iColumn,primalValue,dualValue,gapProduct,delta2W_[iColumn]);
	} 
      } 
    }
  }
#endif
#ifdef SOME_DEBUG_not
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
	  double change = -rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack_[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  largestL = CoinMax(largestL,gapProduct);
	  smallestL = CoinMin(smallestL,gapProduct);
	  nL++;
	  sumL += gapProduct;
	}  
	if (upperBound(iColumn)) {
	  double change = rhsU_[iColumn]-deltaX_[iColumn];
	  double dualValue=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	  double primalValue=upperSlack_[iColumn]+actualPrimalStep_*change;
	  double gapProduct=dualValue*primalValue;
	  largestU = CoinMax(largestU,gapProduct);
	  smallestU = CoinMin(smallestU,gapProduct);
	  nU++;
	  sumU += gapProduct;
	} 
      } 
    }
    double mu = (sumL+sumU)/(static_cast<double> (nL+nU));

    double minBeta = 0.1*mu;
    double maxBeta = 10.0*mu;
    int nBL=0;
    int nAL=0;
    int nBU=0;
    int nAU=0;
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  double change = -rhsL_[iColumn] + deltaX_[iColumn];
	  double dualValue=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  double primalValue=lowerSlack_[iColumn]+actualPrimalStep_*change;
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
	  double dualValue=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	  double primalValue=upperSlack_[iColumn]+actualPrimalStep_*change;
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
  if (cholesky_->type()<20) {
    // not KKT
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
  } else {
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      region1[iColumn] = region1In[iColumn];
    cholesky_->solveKKT(region1,region2,diagonal_,diagonalScaleFactor_);
  }
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
  int numberTotal = numberRows_+numberColumns_;
  //if flagged then entries zero so can do
  // For KKT separate out
  double * region1Save=NULL;//for refinement
  int iColumn;
  if (cholesky_->type()<20) {
    int iColumn;
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      deltaX_[iColumn] = workArray_[iColumn] - solution_[iColumn];
    multiplyAdd(deltaX_+numberColumns_,numberRows_,-1.0,deltaY_,0.0);
    matrix_->times(1.0,deltaX_,deltaY_);
  } else {
    // regions in will be workArray and newError
    // regions out deltaX_ and deltaY_
    multiplyAdd(solution_+numberColumns_,numberRows_,1.0,newError,0.0);
    matrix_->times(-1.0,solution_,newError);
    // This is inefficient but just for now get values which will be in deltay
    int iColumn;
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      deltaX_[iColumn] = workArray_[iColumn] - solution_[iColumn];
    multiplyAdd(deltaX_+numberColumns_,numberRows_,-1.0,deltaY_,0.0);
    matrix_->times(1.0,deltaX_,deltaY_);
  }
  bool goodSolve=false;
  double * regionSave=NULL;//for refinement
  int numberTries=0;
  double relativeError=COIN_DBL_MAX;
  double tryError=1.0e31;
  while (!goodSolve&&numberTries<30) {
    double lastError=relativeError;
    goodSolve=true;
    double maximumRHS;
    double saveMaximum;
    maximumRHS = CoinMax(maximumAbsElement(deltaY_,numberRows_),1.0e-12);
    saveMaximum = maximumRHS;
    if (cholesky_->type()<20) {
      // no kkt
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
	  -workArray_[iColumn];
    } else {
      // KKT
      solveSystem(deltaX_, deltaY_,
		  workArray_,newError,region1Save,regionSave,lastError>1.0e-5);
    }
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
	  if (cholesky_->type()>=20) 
	    region1Save = new double [numberTotal];
        } 
        CoinMemcpyN(deltaY_,numberRows_,regionSave);
	if (cholesky_->type()<20) {
	  // not KKT
	  multiplyAdd(newError,numberRows_,-1.0,deltaY_,0.0);
	} else {
	  // KKT
	  CoinMemcpyN(deltaX_,numberTotal,region1Save);
	  // and back to input region
	  CoinMemcpyN(deltaY_,numberRows_,newError);
	}
      } 
    } else {
      //std::cout <<" worse residual = "<<relativeError;
      //bring back previous
      relativeError=lastError;
      if (regionSave) {
	CoinMemcpyN(regionSave,numberRows_,deltaY_);
	if (cholesky_->type()<20) {
	  // not KKT
	  multiplyAdd(deltaY_,numberRows_,-1.0,deltaX_+numberColumns_,0.0);
	  CoinZeroN(deltaX_,numberColumns_);
	  matrix_->transposeTimes(1.0,deltaY_,deltaX_);
	  //if flagged then entries zero so can do
	  for (iColumn=0;iColumn<numberTotal;iColumn++)
	    deltaX_[iColumn] = deltaX_[iColumn]*diagonal_[iColumn]
	      -workArray_[iColumn];
	} else {
	  // KKT
	  CoinMemcpyN(region1Save,numberTotal,deltaX_);
	}
      } else {
	// disaster
	CoinFillN(deltaX_,numberTotal,1.0);
	CoinFillN(deltaY_,numberRows_,1.0);
	printf("bad cholesky\n");
      }
    }
  } /* endwhile */
  delete [] regionSave;
  delete [] region1Save;
  delete [] newError;
  // now rest
  double extra=eExtra;
  //multiplyAdd(deltaY_,numberRows_,1.0,deltaW_+numberColumns_,0.0);
  //CoinZeroN(deltaW_,numberColumns_);
  //matrix_->transposeTimes(-1.0,deltaY_,deltaW_);
  
  for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    deltaSU_[iColumn]=0.0;
    deltaSL_[iColumn]=0.0;
    deltaZ_[iColumn]=0.0;
    double dd=deltaW_[iColumn];
    deltaW_[iColumn]=0.0;
    if (!flagged(iColumn)) {
      double deltaX = deltaX_[iColumn];
      if (lowerBound(iColumn)) {
	double zValue = rhsZ_[iColumn];
	double gHat = zValue + zVec_[iColumn]*rhsL_[iColumn];
	double slack = lowerSlack_[iColumn]+extra;
	deltaSL_[iColumn] = -rhsL_[iColumn]+deltaX;
	deltaZ_[iColumn]=(gHat-zVec_[iColumn]*deltaX)/slack;
      } 
      if (upperBound(iColumn)) {
	double wValue = rhsW_[iColumn];
	double hHat = wValue - wVec_[iColumn]*rhsU_[iColumn];
	double slack = upperSlack_[iColumn]+extra;
	deltaSU_[iColumn] = rhsU_[iColumn]-deltaX;
	deltaW_[iColumn]=(hHat+wVec_[iColumn]*deltaX)/slack;
      }
      if (0) {
	// different way of calculating
	double gamma2 = gamma_*gamma_;
	double dZ=0.0;
	double dW=0.0;
	double zValue = rhsZ_[iColumn];
	double gHat = zValue + zVec_[iColumn]*rhsL_[iColumn];
	double slackL = lowerSlack_[iColumn]+extra;
	double wValue = rhsW_[iColumn];
	double hHat = wValue - wVec_[iColumn]*rhsU_[iColumn];
	double slackU = upperSlack_[iColumn]+extra;
	double q = rhsC_[iColumn]+gamma2 * deltaX +dd;
	if (primalR_)
	  q += deltaX*primalR_[iColumn];
	dW = (gHat+hHat -slackL*q + (wValue-zValue)*deltaX)/(slackL+slackU);
	dZ = dW + q;
	//printf("B %d old %g %g new %g %g\n",iColumn,deltaZ_[iColumn],
	//deltaW_[iColumn],dZ,dW);
	if (lowerBound(iColumn)) {
	  if (upperBound(iColumn)) {
	    //printf("B %d old %g %g new %g %g\n",iColumn,deltaZ_[iColumn],
	    //deltaW_[iColumn],dZ,dW);
	    deltaW_[iColumn]=dW;
	    deltaZ_[iColumn]=dZ;
	  } else {
	    // just lower
	    //printf("L %d old %g new %g\n",iColumn,deltaZ_[iColumn],
	    //dZ);
	  }
	} else {
	  assert (upperBound(iColumn));
	  //printf("U %d old %g new %g\n",iColumn,deltaW_[iColumn],
	  //dW);
	}
      }
    }
  }
#if 0
  double * check = new double[numberTotal];  
  // Check out rhsC_
  multiplyAdd(deltaY_,numberRows_,-1.0,check+numberColumns_,0.0);
  CoinZeroN(check,numberColumns_);
  matrix_->transposeTimes(1.0,deltaY_,check);
  quadraticDjs(check,deltaX_,-1.0);
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    check[iColumn] += deltaZ_[iColumn]-deltaW_[iColumn];
    if (fabs(check[iColumn]-rhsC_[iColumn])>1.0e-3)
      printf("rhsC %d %g %g\n",iColumn,check[iColumn],rhsC_[iColumn]);
  }
  // Check out rhsZ_
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    check[iColumn] += lowerSlack_[iColumn]*deltaZ_[iColumn]+
      zVec_[iColumn]*deltaSL_[iColumn];
    if (fabs(check[iColumn]-rhsZ_[iColumn])>1.0e-3)
      printf("rhsZ %d %g %g\n",iColumn,check[iColumn],rhsZ_[iColumn]);
  }
  // Check out rhsW_
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    check[iColumn] += upperSlack_[iColumn]*deltaW_[iColumn]+
      wVec_[iColumn]*deltaSU_[iColumn];
    if (fabs(check[iColumn]-rhsW_[iColumn])>1.0e-3)
      printf("rhsW %d %g %g\n",iColumn,check[iColumn],rhsW_[iColumn]);
  }
  delete [] check;
#endif
  return relativeError;
}
// createSolution.  Creates solution from scratch
int ClpPredictorCorrector::createSolution()
{
  int numberTotal = numberRows_+numberColumns_;
  int iColumn;
  double tolerance = primalTolerance();
  // See if quadratic objective
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (!quadraticObj) {
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      if (upper_[iColumn]-lower_[iColumn]>tolerance) 
	clearFixed(iColumn);
      else
	setFixed(iColumn);
    }
  } else {
    // try leaving fixed
    for (iColumn=0;iColumn<numberTotal;iColumn++) 
      clearFixed(iColumn);
  }
  double initialValue =1.0e-12;

  double maximumObjective=0.0;
  double objectiveNorm2=0.0;
  getNorms(cost_,numberTotal,maximumObjective,objectiveNorm2);
  if (!maximumObjective) {
    maximumObjective=1.0; // objective all zero
  } 
  objectiveNorm2=sqrt(objectiveNorm2)/static_cast<double> (numberTotal);
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
  // See if quadratic objective
  if (quadraticObj) {
    // If scaled then really scale matrix
    double scaleFactor = 
      scaleFactor_*optimizationDirection_*objectiveScale_*
      rhsScale_;
    if ((scalingFlag_>0&&rowScale_)||scaleFactor != 1.0) {
      CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
      const int * columnQuadratic = quadratic->getIndices();
      const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
      const int * columnQuadraticLength = quadratic->getVectorLengths();
      double * quadraticElement = quadratic->getMutableElements();
      int numberColumns = quadratic->getNumCols();
      double scale = 1.0/scaleFactor;
      if (scalingFlag_>0&&rowScale_) {
	for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	  double scaleI = columnScale_[iColumn]*scale;
	  for (CoinBigIndex j=columnQuadraticStart[iColumn];
	       j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	    int jColumn = columnQuadratic[j];
	    double scaleJ = columnScale_[jColumn];
	    quadraticElement[j] *= scaleI*scaleJ;
	    objectiveNorm_ = CoinMax(objectiveNorm_,fabs(quadraticElement[j]));
	  }
	}
      } else {
	// not scaled
	for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	  for (CoinBigIndex j=columnQuadraticStart[iColumn];
	       j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	    quadraticElement[j] *= scale;
	    objectiveNorm_ = CoinMax(objectiveNorm_,fabs(quadraticElement[j]));
	  }
	}
      }
    }
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
  if (cholesky_->type()<20) {
    // not KKT
    cholesky_->solve(errorRegion_);
    //create information for solution
    multiplyAdd(errorRegion_,numberRows_,-1.0,deltaX_+numberColumns_,0.0);
    CoinZeroN(deltaX_,numberColumns_);
    matrix_->transposeTimes(1.0,errorRegion_,deltaX_);
  } else {
    // KKT
    // reverse sign on solution
    multiplyAdd(NULL,numberRows_+numberColumns_,0.0,solution_,-1.0);
    solveSystem(deltaX_,errorRegion_,solution_,NULL,NULL,NULL,false);
  }
  initialValue=1.0e2;
  if (rhsNorm_*1.0e-2>initialValue) {
    initialValue=rhsNorm_*1.0e-2;
  } 
  //initialValue = CoinMax(1.0,rhsNorm_);
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
  double extra=1.0e-10;
  double largeGap=1.0e15;
  //double safeObjectiveValue=2.0*objectiveNorm_;
  double safeObjectiveValue=objectiveNorm_+1.0;
  double safeFree=1.0e-1*initialValue;
  //printf("normal safe dual value of %g, primal value of %g\n",
  // safeObjectiveValue,initialValue);
  //safeObjectiveValue=CoinMax(2.0,1.0e-1*safeObjectiveValue);
  //initialValue=CoinMax(100.0,1.0e-1*initialValue);;
  //printf("temp safe dual value of %g, primal value of %g\n",
  // safeObjectiveValue,initialValue);
  double zwLarge=1.0e2*initialValue;
  //zwLarge=1.0e40;
  if (cholesky_->choleskyCondition()<0.0&&cholesky_->type()<20) {
    // looks bad - play safe
    initialValue *=10.0;
    safeObjectiveValue *= 10.0;
    safeFree *= 10.0;
  }
  double gamma2 = gamma_*gamma_; // gamma*gamma will be added to diagonal
#if 0
  fakeSolution[0 ] =   0.072310129 ;
  fakeSolution[1 ] =   0.053083871; 
  fakeSolution[2 ] =      0.178127; 
  fakeSolution[3 ] =    0.13215151; 
  fakeSolution[4 ] =   0.072715642; 
  fakeSolution[5 ] =    0.15680727; 
  fakeSolution[6 ] =    0.16841689; 
  fakeSolution[7 ] =   0.093612798 ;
  fakeSolution[8 ] =   0.072774891 ;
  fakeSolution[9]=1.0;
  initialValue=1.0e-5;
  safeObjectiveValue=1.0e-5;
#endif
  // First do primal side
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      double newValue;
      double setPrimal=initialValue;
      if (quadraticObj) {
	// perturb primal solution a bit
	//fakeSolution[iColumn]  *= 0.002*CoinDrand48()+0.999;
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
          } else {
            newValue=0.5*(upperValue+lowerValue);
          }
        } else {
          //just lower bound
          double fakeValue=fakeSolution[iColumn];
          if (fakeValue<lowerValue+setPrimal) {
            fakeValue=lowerValue+setPrimal;
          } 
          newValue=fakeValue;
        } 
      } else {
        if (upperBound(iColumn)) {
          //just upper bound
          double fakeValue=fakeSolution[iColumn];
          if (fakeValue>upperValue-setPrimal) {
            fakeValue=upperValue-setPrimal;
          } 
          newValue=fakeValue;
        } else {
          //free
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
        } 
      } 
      solution_[iColumn]=newValue;
    } else {
      // fixed
      lowerSlack_[iColumn]=0.0;
      upperSlack_[iColumn]=0.0;
      solution_[iColumn]=lower_[iColumn];
      zVec_[iColumn]=0.0;
      wVec_[iColumn]=0.0;
      diagonal_[iColumn]=0.0;
    } 
  }
  solutionNorm_ =  maximumAbsElement(solution_,numberTotal);
  // Set bounds and do dj including quadratic
  largeGap = CoinMax(1.0e7,1.02*solutionNorm_);
  CoinPackedMatrix * quadratic = NULL;
  const int * columnQuadratic = NULL;
  const CoinBigIndex * columnQuadraticStart = NULL;
  const int * columnQuadraticLength = NULL;
  const double * quadraticElement = NULL;
  if (quadraticObj) {
    quadratic = quadraticObj->quadraticObjective();
    columnQuadratic = quadratic->getIndices();
    columnQuadraticStart = quadratic->getVectorStarts();
    columnQuadraticLength = quadratic->getVectorLengths();
    quadraticElement = quadratic->getElements();
  }
  double quadraticNorm=0.0;
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double primalValue=solution_[iColumn];
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      // Do dj
      double reducedCost=cost_[iColumn];
      if (lowerBound(iColumn)) {
	reducedCost+=linearPerturbation_;
      } 
      if (upperBound(iColumn)) {
	reducedCost-=linearPerturbation_;
      }
      if (quadraticObj&&iColumn<numberColumns_) {
	for (CoinBigIndex j=columnQuadraticStart[iColumn];
	     j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	  int jColumn = columnQuadratic[j];
	  double valueJ = solution_[jColumn];
	  double elementValue = quadraticElement[j];
	  reducedCost += valueJ*elementValue;
	}
	quadraticNorm = CoinMax(quadraticNorm,fabs(reducedCost));
      }
      dj_[iColumn]=reducedCost;
      if (primalValue>lowerValue+largeGap&&primalValue<upperValue-largeGap) {
	clearFixedOrFree(iColumn);
	setLowerBound(iColumn);
	setUpperBound(iColumn);
	lowerValue=CoinMax(lowerValue,primalValue-largeGap);
	upperValue=CoinMin(upperValue,primalValue+largeGap);
	lower_[iColumn]=lowerValue;
	upper_[iColumn]=upperValue;
      }
    }
  }
  safeObjectiveValue=CoinMax(safeObjectiveValue,quadraticNorm);
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double primalValue=solution_[iColumn];
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
      double reducedCost=dj_[iColumn];
      double low=0.0;
      double high=0.0;
      if (lowerBound(iColumn)) {
        if (upperBound(iColumn)) {
          //upper and lower bounds
          if (upperValue-lowerValue>2.0*initialValue) {
            low=primalValue-lowerValue;
            high=upperValue-primalValue;
          } else {
            low=initialValue;
            high=initialValue;
          }
          double s = low+extra;
          double ratioZ;
          if (s<zwLarge) {
            ratioZ=1.0;
          } else {
            ratioZ=sqrt(zwLarge/s);
          } 
          double t = high+extra;
          double ratioT;
          if (t<zwLarge) {
            ratioT=1.0;
          } else {
            ratioT=sqrt(zwLarge/t);
          } 
          //modify s and t
          if (s>largeGap) {
            s=largeGap;
          } 
          if (t>largeGap) {
            t=largeGap;
          } 
          //modify if long long way away from bound
          if (reducedCost>=0.0) {
            zVec_[iColumn]=reducedCost + safeObjectiveValue*ratioZ;
            zVec_[iColumn]=CoinMax(reducedCost, safeObjectiveValue*ratioZ);
            wVec_[iColumn]=safeObjectiveValue*ratioT;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            wVec_[iColumn]=-reducedCost + safeObjectiveValue*ratioT;
            wVec_[iColumn]=CoinMax(-reducedCost , safeObjectiveValue*ratioT);
          }
	  double gammaTerm = gamma2;
	  if (primalR_)
	    gammaTerm += primalR_[iColumn];
          diagonal_[iColumn] = (t*s)/
	    (s*wVec_[iColumn]+t*zVec_[iColumn]+gammaTerm*t*s);
        } else {
          //just lower bound
          low=primalValue-lowerValue;
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
          if (reducedCost>=0.0) {
            zVec_[iColumn]=reducedCost + safeObjectiveValue*ratioZ;
            zVec_[iColumn]=CoinMax(reducedCost , safeObjectiveValue*ratioZ);
            wVec_[iColumn]=0.0;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            wVec_[iColumn]=0.0;
          } 
	  double gammaTerm = gamma2;
	  if (primalR_)
	    gammaTerm += primalR_[iColumn];
          diagonal_[iColumn] = s/(zVec_[iColumn]+s*gammaTerm);
        } 
      } else {
        if (upperBound(iColumn)) {
          //just upper bound
          low=0.0;
          high=upperValue-primalValue;
          double t = high+extra;
          double ratioT;
          if (t<zwLarge) {
            ratioT=1.0;
          } else {
            ratioT=sqrt(zwLarge/t);
          } 
          //modify t
          if (t>largeGap) {
            t=largeGap;
          } 
          if (reducedCost>=0.0) {
            zVec_[iColumn]=0.0;
            wVec_[iColumn]=safeObjectiveValue*ratioT;
          } else {
            zVec_[iColumn]=0.0;
            wVec_[iColumn]=-reducedCost + safeObjectiveValue*ratioT;
            wVec_[iColumn]=CoinMax(-reducedCost , safeObjectiveValue*ratioT);
          } 
	  double gammaTerm = gamma2;
	  if (primalR_)
	    gammaTerm += primalR_[iColumn];
          diagonal_[iColumn] =  t/(wVec_[iColumn]+t*gammaTerm);
        } 
      } 
      lowerSlack_[iColumn]=low;
      upperSlack_[iColumn]=high;
    }
  }
#if 0
  if (solution_[0]>0.0) {
    for (int i=0;i<numberTotal;i++)
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	     diagonal_[i],fabs(dj_[i]),
	     lowerSlack_[i],zVec_[i],
	     upperSlack_[i],wVec_[i]);
  } else {
    for (int i=0;i<numberTotal;i++)
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	     diagonal_[i],fabs(dj_[i]),
	     upperSlack_[i],wVec_[i],
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
  int numberTotal = numberRows_+numberColumns_;
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
  for (int iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!fixedOrFree(iColumn)) {
      numberComplementarityPairs++;
      //can collapse as if no lower bound both zVec and deltaZ 0.0
      double newZ=0.0;
      double newW=0.0;
      if (lowerBound(iColumn)) {
	numberComplementarityItems++;
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=zVec_[iColumn];
          primalValue=lowerSlack_[iColumn];
        } else {
          double change;
	  change =solution_[iColumn]+deltaX_[iColumn]-lowerSlack_[iColumn]-lower_[iColumn];
          dualValue=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
	  newZ=dualValue;
          primalValue=lowerSlack_[iColumn]+actualPrimalStep_*change;
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
	//printf("l %d prim %g dual %g totalGap %g\n",
	//   iColumn,primalValue,dualValue,gap);
        if (gapProduct>largestGap) {
          largestGap=gapProduct;
        }
	smallestGap = CoinMin(smallestGap,gapProduct);
        if (dualValue>dualTolerance&&primalValue>primalTolerance) {
          toleranceGap+=dualValue*primalValue;
        } 
      } 
      if (upperBound(iColumn)) {
	numberComplementarityItems++;
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=wVec_[iColumn];
          primalValue=upperSlack_[iColumn];
        } else {
          double change;
	  change =upper_[iColumn]-solution_[iColumn]-deltaX_[iColumn]-upperSlack_[iColumn];
          dualValue=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	  newW=dualValue;
          primalValue=upperSlack_[iColumn]+actualPrimalStep_*change;
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
	//printf("u %d prim %g dual %g totalGap %g\n",
	//   iColumn,primalValue,dualValue,gap);
        if (gapProduct>largestGap) {
          largestGap=gapProduct;
        } 
        if (dualValue>dualTolerance&&primalValue>primalTolerance) {
          toleranceGap+=dualValue*primalValue;
        } 
      } 
    } 
  }
  //if (numberIterations_>4)
  //exit(9);
  if (!phase&&numberNegativeGaps) {
      handler_->message(CLP_BARRIER_NEGATIVE_GAPS,messages_)
    <<numberNegativeGaps<<sumNegativeGap
    <<CoinMessageEol;
  } 
  
  //in case all free!
  if (!numberComplementarityPairs) {
    numberComplementarityPairs=1;
  } 
#ifdef SOME_DEBUG
  printf("with d,p steps %g,%g gap %g - smallest %g, largest %g, pairs %d\n",
	 actualDualStep_,actualPrimalStep_,
	 gap,smallestGap,largestGap,numberComplementarityPairs);
#endif
  return gap;
}
// setupForSolve.
//phase 0=affine , 1 = corrector , 2 = primal-dual
void ClpPredictorCorrector::setupForSolve(const int phase)
{
  double extra =eExtra;
  int numberTotal = numberRows_ + numberColumns_;
  int iColumn;
#ifdef SOME_DEBUG
  printf("phase %d in setupForSolve, mu %.18g\n",phase,mu_);
#endif
  double gamma2 = gamma_*gamma_; // gamma*gamma will be added to diagonal
  switch (phase) {
  case 0:
    CoinMemcpyN(errorRegion_,numberRows_,rhsB_);
    if (delta_||dualR_) {
      // add in regularization
      double delta2 = delta_*delta_;
      for (int iRow=0;iRow<numberRows_;iRow++) {
	rhsB_[iRow] -= delta2*dual_[iRow];
	if (dualR_)
	  rhsB_[iRow] -= dualR_[iRow]*dual_[iRow];
      }
    }
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsC_[iColumn]=0.0;
      rhsU_[iColumn]=0.0;
      rhsL_[iColumn]=0.0;
      rhsZ_[iColumn]=0.0;
      rhsW_[iColumn]=0.0;
      if (!flagged(iColumn)) {
	rhsC_[iColumn] = dj_[iColumn]-zVec_[iColumn]+wVec_[iColumn];
	rhsC_[iColumn] += gamma2*solution_[iColumn];
	if (primalR_)
	  rhsC_[iColumn] += primalR_[iColumn]*solution_[iColumn];
	if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = -zVec_[iColumn]*(lowerSlack_[iColumn]+extra);
	  rhsL_[iColumn] = CoinMax(0.0,(lower_[iColumn]+lowerSlack_[iColumn])-solution_[iColumn]);
	} 
	if (upperBound(iColumn)) {
	  rhsW_[iColumn] = -wVec_[iColumn]*(upperSlack_[iColumn]+extra);
	  rhsU_[iColumn] = CoinMin(0.0,(upper_[iColumn]-upperSlack_[iColumn])-solution_[iColumn]);
	}
      }
    } 
    if (0) {
      int i=1324;
      printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,solution_[i],
	     diagonal_[i],dj_[i],
	     lowerSlack_[i],zVec_[i],
	     upperSlack_[i],wVec_[i]);
      printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	     rhsZ_[i],rhsL_[i],
	     rhsW_[i],rhsU_[i]);
    }
#if 0
    for (int i=0;i<3;i++) {
      if (!fabs(rhsZ_[i]))
	rhsZ_[i]=0.0;
      if (!fabs(rhsW_[i]))
	rhsW_[i]=0.0;
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
	       upperSlack_[i],wVec_[i]);
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	       rhsZ_[i],rhsL_[i],
	       rhsW_[i],rhsU_[i]);
    } else {
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,solution_[i],
	       diagonal_[i],dj_[i],
	       lowerSlack_[i],zVec_[i],
	       upperSlack_[i],wVec_[i]);
      for (int i=0;i<3;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,rhsC_[i],
	       rhsZ_[i],rhsL_[i],
	       rhsW_[i],rhsU_[i]);
    }
#endif
    break;
  case 1:
    // could be stored in delta2?
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsZ_[iColumn]=0.0;
      rhsW_[iColumn]=0.0;
      if (!flagged(iColumn)) {
        if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = mu_ -zVec_[iColumn]*(lowerSlack_[iColumn]+extra)
	    - deltaZ_[iColumn]*deltaX_[iColumn];
	  // To bring in line with OSL
	  rhsZ_[iColumn] += deltaZ_[iColumn]*rhsL_[iColumn];
        } 
        if (upperBound(iColumn)) {
	  rhsW_[iColumn] = mu_ -wVec_[iColumn]*(upperSlack_[iColumn]+extra)
	    +deltaW_[iColumn]*deltaX_[iColumn];
	  // To bring in line with OSL
	  rhsW_[iColumn] -= deltaW_[iColumn]*rhsU_[iColumn];
        } 
      } 
    } 
#if 0
    for (int i=0;i<numberTotal;i++) {
      if (!fabs(rhsZ_[i]))
	rhsZ_[i]=0.0;
      if (!fabs(rhsW_[i]))
	rhsW_[i]=0.0;
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
	       upperSlack_[i],wVec_[i]);
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,fabs(rhsC_[i]),
	       rhsZ_[i],rhsL_[i],
	       rhsW_[i],rhsU_[i]);
    } else {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g %.18g %.18g\n",i,fabs(solution_[i]),
	       diagonal_[i],fabs(dj_[i]),
	       upperSlack_[i],wVec_[i],
	       lowerSlack_[i],zVec_[i] );
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g %.18g %.18g %.18g %.18g\n",i,fabs(rhsC_[i]),
	       rhsW_[i],rhsU_[i],
	       rhsZ_[i],rhsL_[i]);
    }
    exit(66);
#endif
    break;
  case 2:
    CoinMemcpyN(errorRegion_,numberRows_,rhsB_);
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      rhsZ_[iColumn]=0.0;
      rhsW_[iColumn]=0.0;
      if (!flagged(iColumn)) {
        if (lowerBound(iColumn)) {
	  rhsZ_[iColumn] = mu_ - zVec_[iColumn]*(lowerSlack_[iColumn]+extra);
        } 
        if (upperBound(iColumn)) {
	  rhsW_[iColumn] = mu_ - wVec_[iColumn]*(upperSlack_[iColumn]+extra);
        }
      }
    } 
    break;
  case 3:
    {
      double minBeta = 0.1*mu_;
      double maxBeta = 10.0*mu_;
      double dualStep = CoinMin(1.0,actualDualStep_+0.1);
      double primalStep = CoinMin(1.0,actualPrimalStep_+0.1);
#ifdef SOME_DEBUG
      printf("good complementarity range %g to %g\n",minBeta,maxBeta);
#endif
      //minBeta=0.0;
      //maxBeta=COIN_DBL_MAX;
      for (iColumn=0;iColumn<numberTotal;iColumn++) {
	if (!flagged(iColumn)) {
	  if (lowerBound(iColumn)) {
	    double change = -rhsL_[iColumn] + deltaX_[iColumn];
	    double dualValue=zVec_[iColumn]+dualStep*deltaZ_[iColumn];
	    double primalValue=lowerSlack_[iColumn]+primalStep*change;
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
	      value= CoinMax(maxBeta-gapProduct,-maxBeta);
	      assert (value<0.0);
	    }
	    rhsZ_[iColumn] += value;
	  }  
	  if (upperBound(iColumn)) {
	    double change = rhsU_[iColumn]-deltaX_[iColumn];
	    double dualValue=wVec_[iColumn]+dualStep*deltaW_[iColumn];
	    double primalValue=upperSlack_[iColumn]+primalStep*change;
	    double gapProduct=dualValue*primalValue;
	    if (gapProduct>0.0&&dualValue<0.0)
	      gapProduct = - gapProduct;
#ifdef FULL_DEBUG
	    delta2W_[iColumn]=gapProduct;
	    if (delta2W_[iColumn]<minBeta||delta2W_[iColumn]>maxBeta)
	      printf("upper %d primal %g, dual %g, gap %g\n",
		     iColumn,primalValue,dualValue,gapProduct);
#endif
	    double value= 0.0;
	    if (gapProduct<minBeta) {
	      value= (minBeta-gapProduct);
	      assert (value>0.0);
	    } else if (gapProduct>maxBeta) {
	      value= CoinMax(maxBeta-gapProduct,-maxBeta);
	      assert (value<0.0);
	    }
	    rhsW_[iColumn] += value;
	  } 
	} 
      }
    }
    break;
  } /* endswitch */
  if (cholesky_->type()<20) {
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      double value = rhsC_[iColumn];
      double zValue = rhsZ_[iColumn];
      double wValue = rhsW_[iColumn];
#if 0
#if 1
      if (phase==0) {
	// more accurate
	value = dj[iColumn];
	zValue=0.0;
	wValue=0.0;
      } else if (phase==2) {
	// more accurate
	value = dj[iColumn];
	zValue=mu_;
	wValue=mu_;
      }
#endif
      assert (rhsL_[iColumn]>=0.0);
      assert (rhsU_[iColumn]<=0.0);
      if (lowerBound(iColumn)) {
	value += (-zVec_[iColumn]*rhsL_[iColumn]-zValue)/
	  (lowerSlack_[iColumn]+extra);
      }
      if (upperBound(iColumn)) {
	value += (wValue-wVec_[iColumn]*rhsU_[iColumn])/
	  (upperSlack_[iColumn]+extra);
      }
#else
      if (lowerBound(iColumn)) {
	double gHat = zValue + zVec_[iColumn]*rhsL_[iColumn];
	value -= gHat/(lowerSlack_[iColumn]+extra);
      }
      if (upperBound(iColumn)) {
	double hHat = wValue - wVec_[iColumn]*rhsU_[iColumn];
	value += hHat/(upperSlack_[iColumn]+extra);
      }
#endif
      workArray_[iColumn]=diagonal_[iColumn]*value;
    } 
#if 0
    if (solution_[0]>0.0) {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g\n",i,workArray_[i]);
    } else {
      for (int i=0;i<numberTotal;i++)
	printf("%d %.18g\n",i,workArray_[i]);
    }
    exit(66);
#endif
  } else {
    // KKT
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      double value = rhsC_[iColumn];
      double zValue = rhsZ_[iColumn];
      double wValue = rhsW_[iColumn];
      if (lowerBound(iColumn)) {
	double gHat = zValue + zVec_[iColumn]*rhsL_[iColumn];
	value -= gHat/(lowerSlack_[iColumn]+extra);
      }
      if (upperBound(iColumn)) {
	double hHat = wValue - wVec_[iColumn]*rhsU_[iColumn];
	value += hHat/(upperSlack_[iColumn]+extra);
      }
      workArray_[iColumn]=value;
    }
  }
}
//method: sees if looks plausible change in complementarity
bool ClpPredictorCorrector::checkGoodMove(const bool doCorrector,
					  double & bestNextGap,
					  bool allowIncreasingGap)
{
  const double beta3 = 0.99997;
  bool goodMove=false;
  int nextNumber;
  int nextNumberItems;
  int numberTotal = numberRows_+numberColumns_;
  double returnGap=bestNextGap;
  double nextGap=complementarityGap(nextNumber,nextNumberItems,2);
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (nextGap>bestNextGap&&nextGap>0.9*complementarityGap_&&doCorrector
      &&!quadraticObj&&!allowIncreasingGap) {
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
    goodMove=checkGoodMove2(step,gap,allowIncreasingGap);
    if (goodMove)
      returnGap=gap;
  } else {
    goodMove=true;
  } 
  if (goodMove)
    goodMove=checkGoodMove2(step,bestNextGap,allowIncreasingGap);
  // Say good if small
  //if (quadraticObj) {
  if (CoinMax(actualDualStep_,actualPrimalStep_)<1.0e-6)
    goodMove=true;
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
    //if (quadraticObj)
    //actualPrimalStep_ *=0.5;
    actualDualStep_=step;
    goodMove=checkGoodMove2(step,bestNextGap,allowIncreasingGap);
    int pass=0;
    while (!goodMove) {
      pass++;
      double gap = bestNextGap;
      goodMove=checkGoodMove2(step,gap,allowIncreasingGap);
      if (goodMove||pass>3) {
	returnGap=gap;
	break;
      }
      if (step<1.0e-4) {
        break;
      } 
      step*=0.5;
      actualPrimalStep_=step;
      //if (quadraticObj)
      //actualPrimalStep_ *=0.5;
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
    //double sumPerturbCost=0.0;
    for (int iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  //sumPerturbCost+=deltaX_[iColumn];
	  deltaObjectiveDual+=deltaZ_[iColumn]*lower_[iColumn];
	} 
	if (upperBound(iColumn)) {
	  //sumPerturbCost-=deltaX_[iColumn];
	  deltaObjectiveDual-=deltaW_[iColumn]*upper_[iColumn];
	} 
	double change = fabs(workArray_[iColumn]-deltaZ_[iColumn]+deltaW_[iColumn]);
	error = CoinMax(change,error);
      } 
      deltaObjectivePrimal += cost_[iColumn] * deltaX_[iColumn];
    } 
    //deltaObjectivePrimal+=sumPerturbCost*linearPerturbation_;
    double testValue;
    if (error>0.0) {
      testValue=1.0e1*CoinMax(maximumDualError_,1.0e-12)/error;
    } else {
      testValue=1.0e1;
    } 
    // If quadratic then primal step may compensate
    if (testValue<actualDualStep_&&!quadraticObj) {
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
    double ratio = 1.0e1*CoinMax(maximumRHSError_,1.0e-12)/maximumRHSChange_;
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
bool ClpPredictorCorrector::checkGoodMove2(double move,
					   double & bestNextGap,
					   bool allowIncreasingGap)
{
  double complementarityMultiplier =1.0/numberComplementarityPairs_;
  const double gamma = 1.0e-8;
  const double gammap = 1.0e-8;
  double gammad = 1.0e-8;
  int nextNumber;
  int nextNumberItems;
  double nextGap=complementarityGap(nextNumber,nextNumberItems,2);
  if (nextGap>bestNextGap&&!allowIncreasingGap)
    return false;
  double lowerBoundGap = gamma*nextGap*complementarityMultiplier;
  bool goodMove=true;
  int iColumn;
  for ( iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!flagged(iColumn)) {
      if (lowerBound(iColumn)) {
	double part1=lowerSlack_[iColumn]+actualPrimalStep_*deltaSL_[iColumn];
	double part2=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
	if (part1*part2<lowerBoundGap) {
	  goodMove=false;
	  break;
	} 
      } 
      if (upperBound(iColumn)) {
	double part1=upperSlack_[iColumn]+actualPrimalStep_*deltaSU_[iColumn];
	double part2=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	if (part1*part2<lowerBoundGap) {
	  goodMove=false;
	  break;
	} 
      } 
    } 
  } 
   double * nextDj=NULL;
   double maximumDualError = maximumDualError_;
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (quadraticObj) {
    // change gammad
    gammad=1.0e-4;
    double gamma2 = gamma_*gamma_;
    nextDj = new double [numberColumns_];
    double * nextSolution = new double [numberColumns_];
    // put next primal into nextSolution
    for ( iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
	nextSolution[iColumn]=solution_[iColumn]+
	  actualPrimalStep_*deltaX_[iColumn];
      } else {
	nextSolution[iColumn]=solution_[iColumn];
      }
    }
    // do reduced costs
    CoinMemcpyN(cost_,numberColumns_,nextDj);
    matrix_->transposeTimes(-1.0,dual_,nextDj);
    matrix_->transposeTimes(-actualDualStep_,deltaY_,nextDj);
    quadraticDjs(nextDj,nextSolution,1.0);
    delete [] nextSolution;
    CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
    const int * columnQuadraticLength = quadratic->getVectorLengths();
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (!fixedOrFree(iColumn)) {
	double newZ=0.0;
	double newW=0.0;
	if (lowerBound(iColumn)) {
          newZ=zVec_[iColumn]+actualDualStep_*deltaZ_[iColumn];
        } 
	if (upperBound(iColumn)) {
          newW=wVec_[iColumn]+actualDualStep_*deltaW_[iColumn];
	} 
	if (columnQuadraticLength[iColumn]) {
	  double gammaTerm = gamma2;
	  if (primalR_)
	    gammaTerm += primalR_[iColumn];
	  //double dualInfeasibility=
	  //dj_[iColumn]-zVec_[iColumn]+wVec_[iColumn]
	  //+gammaTerm*solution_[iColumn];
	  double newInfeasibility=
	    nextDj[iColumn]-newZ+newW
	    +gammaTerm*(solution_[iColumn]+actualPrimalStep_*deltaX_[iColumn]);
	  maximumDualError = CoinMax(maximumDualError,newInfeasibility);
	  //if (fabs(newInfeasibility)>CoinMax(2000.0*maximumDualError_,1.0e-2)) {
	  //if (dualInfeasibility*newInfeasibility<0.0) {
	  //  printf("%d current %g next %g\n",iColumn,dualInfeasibility,
	  //       newInfeasibility);
	  //  goodMove=false;
	  //}
	  //}
	}
      } 
    }
    delete [] nextDj;
  }
 //      Satisfy g_p(alpha)?
  if (rhsNorm_>solutionNorm_) {
    solutionNorm_=rhsNorm_;
  } 
  double errorCheck=maximumRHSError_/solutionNorm_;
  if (errorCheck<maximumBoundInfeasibility_) {
    errorCheck=maximumBoundInfeasibility_;
  } 
  // scale back move
  move = CoinMin(move,0.95);
  //scale
  if ((1.0-move)*errorCheck>primalTolerance()) {
    if (nextGap<gammap*(1.0-move)*errorCheck) {
      goodMove=false;
    } 
  } 
  //      Satisfy g_d(alpha)?
  errorCheck=maximumDualError/objectiveNorm_;
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
  int numberTotal = numberRows_+numberColumns_;
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
  //double nextMu = nextGap/(static_cast<double>(2*numberComplementarityPairs_));
  //printf("using gap of %g\n",nextMu);
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
  double largeGap2 = CoinMax(1.0e7,1.0e2*solutionNorm_);
  //largeGap2 = 1.0e9;
  // When to start looking at killing (factor0
  double killFactor;
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (!quadraticObj||1) {
    if (numberIterations_<50) {
      killFactor = 1.0;
    } else if (numberIterations_<100) {
      killFactor = 10.0;
      stepLength_=CoinMax(stepLength_,0.9995);
    }else if (numberIterations_<150) {
      killFactor = 100.0;
      stepLength_=CoinMax(stepLength_,0.99995);
    } else {
      killFactor = 1.0e5;
      stepLength_=CoinMax(stepLength_,0.999995);
    }
  } else {
    killFactor=1.0;
  }
  // put next primal into deltaSL_
  int iColumn;
  int iRow;
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    double thisWeight=deltaX_[iColumn];
    double newPrimal=solution_[iColumn]+1.0*actualPrimalStep_*thisWeight;
    deltaSL_[iColumn]=newPrimal;
  }
#if 0
  // nice idea but doesn't work
  multiplyAdd(solution_+numberColumns_,numberRows_,-1.0,errorRegion_,0.0);
  matrix_->times(1.0,solution_,errorRegion_);
  multiplyAdd(deltaSL_+numberColumns_,numberRows_,-1.0,rhsFixRegion_,0.0);
  matrix_->times(1.0,deltaSL_,rhsFixRegion_);
  double newNorm =  maximumAbsElement(deltaSL_,numberTotal);
  double tol = newNorm*primalTolerance();
  bool goneInf=false;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double value=errorRegion_[iRow];
    double valueNew=rhsFixRegion_[iRow];
    if (fabs(value)<tol&&fabs(valueNew)>tol) {
      printf("row %d old %g new %g\n",iRow,value,valueNew);
      goneInf=true;
    }
  }
  if (goneInf) {
    actualPrimalStep_ *= 0.5;
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      double thisWeight=deltaX_[iColumn];
      double newPrimal=solution_[iColumn]+1.0*actualPrimalStep_*thisWeight;
      deltaSL_[iColumn]=newPrimal;
    }
  }
  CoinZeroN(errorRegion_,numberRows_);
  CoinZeroN(rhsFixRegion_,numberRows_);
#endif
  // do reduced costs
  CoinMemcpyN(dual_,numberRows_,dj_+numberColumns_);
  CoinMemcpyN(cost_,numberColumns_,dj_);
  double quadraticOffset=quadraticDjs(dj_,deltaSL_,1.0);
  // Save modified costs for fixed variables
  CoinMemcpyN(dj_,numberColumns_,deltaSU_);
  matrix_->transposeTimes(-1.0,dual_,dj_);
  double gamma2 = gamma_*gamma_; // gamma*gamma will be added to diagonal
  double gammaOffset=0.0;
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double reducedCost=dj_[iColumn];
      bool thisKilled=false;
      double zValue = zVec_[iColumn] + actualDualStep_*deltaZ_[iColumn];
      double wValue = wVec_[iColumn] + actualDualStep_*deltaW_[iColumn];
      zVec_[iColumn]=zValue;
      wVec_[iColumn]=wValue;
      double thisWeight=deltaX_[iColumn];
      double oldPrimal = solution_[iColumn];
      double newPrimal=solution_[iColumn]+actualPrimalStep_*thisWeight;
      double dualObjectiveThis=0.0;
      double sUpper=extra;
      double sLower=extra;
      double kill;
      if (fabs(newPrimal)>1.0e4) {
        kill=killTolerance*1.0e-4*newPrimal;
      } else {
        kill=killTolerance;
      } 
      kill*=1.0e-3;//be conservative
      double smallerSlack=COIN_DBL_MAX;
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
	  lowerSlack_[iColumn] = largeGap2;
	  upper_[iColumn]=newPrimal+largeGap2;
	  upperSlack_[iColumn] = largeGap2;
	} else {
	  lower_[iColumn]=trueLower;
	  setLowerBound(iColumn);
	  lowerSlack_[iColumn] = CoinMax(newPrimal-trueLower,1.0);
	  upper_[iColumn]=trueUpper;
	  setUpperBound(iColumn);
	  upperSlack_[iColumn] = CoinMax(trueUpper-newPrimal,1.0);
	}
      } else if (fakeNewBounds) {
	lower_[iColumn]=newPrimal-largeGap2;
	lowerSlack_[iColumn] = largeGap2;
	upper_[iColumn]=newPrimal+largeGap2;
	upperSlack_[iColumn] = largeGap2;
	// so we can just have one test
	fakeOldBounds=true;
      }
      double lowerBoundInfeasibility=0.0;
      double upperBoundInfeasibility=0.0;
      if (lowerBound(iColumn)) {
        double oldSlack = lowerSlack_[iColumn];
        double newSlack;
	newSlack=
	  lowerSlack_[iColumn]+actualPrimalStep_*(oldPrimal-oldSlack
						 + thisWeight-lower_[iColumn]);
	if (fakeOldBounds)
	  newSlack = lowerSlack_[iColumn];
        double epsilon = fabs(newSlack)*epsilonBase;
        if (epsilon>1.0e-5) {
          //cout<<"bad"<<endl;
          epsilon=1.0e-5;
        } 
	//epsilon=1.0e-14;
        //make sure reasonable
        if (zValue<epsilon) {
          zValue=epsilon;
        } 
        double feasibleSlack=newPrimal-lower_[iColumn];
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
        if (zVec_[iColumn]>dualTolerance) {
          dualObjectiveThis+=lower_[iColumn]*zVec_[iColumn];
        } 
        lowerSlack_[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
	lowerBoundInfeasibility = fabs(newPrimal-lowerSlack_[iColumn]-lower_[iColumn]);
        if (lowerSlack_[iColumn]<=kill*killFactor&&fabs(newPrimal-lower_[iColumn])<=kill*killFactor) {
	  double step = CoinMin(actualPrimalStep_*1.1,1.0);
	  double newPrimal2=solution_[iColumn]+step*thisWeight;
	  if (newPrimal2<newPrimal&&dj_[iColumn]>1.0e-5&&numberIterations_>50-40) {
	    newPrimal=lower_[iColumn];
	    lowerSlack_[iColumn]=0.0;
	    //printf("fixing %d to lower\n",iColumn);
	  }
	}
        if (lowerSlack_[iColumn]<=kill&&fabs(newPrimal-lower_[iColumn])<=kill) {
          //may be better to leave at value?
          newPrimal=lower_[iColumn];
          lowerSlack_[iColumn]=0.0;
          thisKilled=true;
          //cout<<j<<" l "<<reducedCost<<" "<<zVec_[iColumn]<<endl;
        } else {
          sLower+=lowerSlack_[iColumn];
        } 
      } 
      if (upperBound(iColumn)) {
        double oldSlack = upperSlack_[iColumn];
        double newSlack;
	newSlack=
	  upperSlack_[iColumn]+actualPrimalStep_*(-oldPrimal-oldSlack
						 - thisWeight+upper_[iColumn]);
	if (fakeOldBounds)
	  newSlack = upperSlack_[iColumn];
        double epsilon = fabs(newSlack)*epsilonBase;
        if (epsilon>1.0e-5) {
          //cout<<"bad"<<endl;
          epsilon=1.0e-5;
        } 
        //make sure reasonable
	//epsilon=1.0e-14;
        if (wValue<epsilon) {
          wValue=epsilon;
        } 
        double feasibleSlack=upper_[iColumn]-newPrimal;
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
        if (wVec_[iColumn]>dualTolerance) {
          dualObjectiveThis-=upper_[iColumn]*wVec_[iColumn];
        } 
        upperSlack_[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
	upperBoundInfeasibility = fabs(newPrimal+upperSlack_[iColumn]-upper_[iColumn]);
        if (upperSlack_[iColumn]<=kill*killFactor&&fabs(newPrimal-upper_[iColumn])<=kill*killFactor) {
	  double step = CoinMin(actualPrimalStep_*1.1,1.0);
	  double newPrimal2=solution_[iColumn]+step*thisWeight;
	  if (newPrimal2>newPrimal&&dj_[iColumn]<-1.0e-5&&numberIterations_>50-40) {
	    newPrimal=upper_[iColumn];
	    upperSlack_[iColumn]=0.0;
	    //printf("fixing %d to upper\n",iColumn);
	  }
	}
        if (upperSlack_[iColumn]<=kill&&fabs(newPrimal-upper_[iColumn])<=kill) {
          //may be better to leave at value?
          newPrimal=upper_[iColumn];
          upperSlack_[iColumn]=0.0;
          thisKilled=true;
        } else {
          sUpper+=upperSlack_[iColumn];
        } 
      } 
      solution_[iColumn]=newPrimal;
      if (fabs(newPrimal)>solutionNorm) {
        solutionNorm=fabs(newPrimal);
      } 
      if (!thisKilled) {
	double gammaTerm = gamma2;
	if (primalR_) {
	  gammaTerm += primalR_[iColumn];
	  quadraticOffset += newPrimal*newPrimal*primalR_[iColumn];
	}
        double dualInfeasibility=
	  reducedCost-zVec_[iColumn]+wVec_[iColumn]+gammaTerm*newPrimal;
        if (fabs(dualInfeasibility)>dualTolerance) {
#if 0
	  if (dualInfeasibility>0.0) {
	    // To improve we could reduce t and/or increase z
	    if (lowerBound(iColumn)) {
	      double complementarity =zVec_[iColumn]*lowerSlack_[iColumn];
	      if (complementarity<nextMu) {
		double change=
		  CoinMin(dualInfeasibility,
		      (nextMu-complementarity)/lowerSlack_[iColumn]);
		dualInfeasibility -= change;
		printf("%d lb locomp %g - dual inf from %g to %g\n",
		       iColumn,complementarity,dualInfeasibility+change,
		       dualInfeasibility);
		zVec_[iColumn] += change;
		zValue = CoinMax(zVec_[iColumn],1.0e-12);
	      }
	    }
	    if (upperBound(iColumn)) {
	      double complementarity =wVec_[iColumn]*upperSlack_[iColumn];
	      if (complementarity>nextMu) {
		double change=
		  CoinMin(dualInfeasibility,
		      (complementarity-nextMu)/upperSlack_[iColumn]);
		dualInfeasibility -= change;
		printf("%d ub hicomp %g - dual inf from %g to %g\n",
		       iColumn,complementarity,dualInfeasibility+change,
		       dualInfeasibility);
		wVec_[iColumn] -= change;
		wValue = CoinMax(wVec_[iColumn],1.0e-12);
	      }
	    }
	  } else {
	    // To improve we could reduce z and/or increase t
	    if (lowerBound(iColumn)) {
	      double complementarity =zVec_[iColumn]*lowerSlack_[iColumn];
	      if (complementarity>nextMu) {
		double change=
		  CoinMax(dualInfeasibility,
		      (nextMu-complementarity)/lowerSlack_[iColumn]);
		dualInfeasibility -= change;
		printf("%d lb hicomp %g - dual inf from %g to %g\n",
		       iColumn,complementarity,dualInfeasibility+change,
		       dualInfeasibility);
		zVec_[iColumn] += change;
		zValue = CoinMax(zVec_[iColumn],1.0e-12);
	      }
	    }
	    if (upperBound(iColumn)) {
	      double complementarity =wVec_[iColumn]*upperSlack_[iColumn];
	      if (complementarity<nextMu) {
		double change=
		  CoinMax(dualInfeasibility,
		      (complementarity-nextMu)/upperSlack_[iColumn]);
		dualInfeasibility -= change;
		printf("%d ub locomp %g - dual inf from %g to %g\n",
		       iColumn,complementarity,dualInfeasibility+change,
		       dualInfeasibility);
		wVec_[iColumn] -= change;
		wValue = CoinMax(wVec_[iColumn],1.0e-12);
	      }
	    }
	  }
#endif
	  dualFake+=newPrimal*dualInfeasibility;
        } 
	if (lowerBoundInfeasibility>maximumBoundInfeasibility) {
	  maximumBoundInfeasibility=lowerBoundInfeasibility;
	} 
	if (upperBoundInfeasibility>maximumBoundInfeasibility) {
	  maximumBoundInfeasibility=upperBoundInfeasibility;
	} 
        dualInfeasibility=fabs(dualInfeasibility);
        if (dualInfeasibility>maximumDualError) {
	  //printf("bad dual %d %g\n",iColumn,
	  // reducedCost-zVec_[iColumn]+wVec_[iColumn]+gammaTerm*newPrimal);
          maximumDualError=dualInfeasibility;
        } 
        dualObjectiveValue+=dualObjectiveThis;
	gammaOffset += newPrimal*newPrimal;
        if (sLower>largeGap) {
          sLower=largeGap;
        } 
        if (sUpper>largeGap) {
          sUpper=largeGap;
        } 
#if 1
        double divisor = sLower*wValue+sUpper*zValue+gammaTerm*sLower*sUpper;
        double diagonalValue=(sUpper*sLower)/divisor;
#else
        double divisor = sLower*wValue+sUpper*zValue+gammaTerm*sLower*sUpper;
        double diagonalValue2=(sUpper*sLower)/divisor;
	double diagonalValue;
	if (!lowerBound(iColumn)) {
	  diagonalValue = wValue/sUpper + gammaTerm;
	} else if (!upperBound(iColumn)) {
	  diagonalValue = zValue/sLower + gammaTerm;
	} else {
	  diagonalValue = zValue/sLower + wValue/sUpper + gammaTerm;
	}
	diagonalValue = 1.0/diagonalValue;
#endif
        diagonal_[iColumn]=diagonalValue;
        //FUDGE
        if (diagonalValue>diagonalLimit) {
#ifdef COIN_DEVELOP
          std::cout<<"large diagonal "<<diagonalValue<<std::endl;
#endif
          diagonal_[iColumn]=diagonalLimit;
        } 
        if (diagonalValue<1.0e-10) {
          //std::cout<<"small diagonal "<<diagonalValue<<std::endl;
        } 
        if (diagonalValue>largestDiagonal) {
          largestDiagonal=diagonalValue;
        } 
        if (diagonalValue<smallestDiagonal) {
          smallestDiagonal=diagonalValue;
        } 
        deltaX_[iColumn]=0.0;
      } else {
        numberKilled++;
        diagonal_[iColumn]=0.0;
        zVec_[iColumn]=0.0;
        wVec_[iColumn]=0.0;
        setFlagged(iColumn);
        setFixedOrFree(iColumn);
        deltaX_[iColumn]=newPrimal;
	offsetObjective+=newPrimal*deltaSU_[iColumn];
      } 
    } else {
      deltaX_[iColumn]=solution_[iColumn];
      diagonal_[iColumn]=0.0;
      offsetObjective+=solution_[iColumn]*deltaSU_[iColumn];
      if (upper_[iColumn]-lower_[iColumn]>1.0e-5) {
        if (solution_[iColumn]<lower_[iColumn]+1.0e-8&&dj_[iColumn]<-1.0e-8) {
          if (-dj_[iColumn]>maximumDJInfeasibility) {
            maximumDJInfeasibility=-dj_[iColumn];
          } 
        } 
        if (solution_[iColumn]>upper_[iColumn]-1.0e-8&&dj_[iColumn]>1.0e-8) {
          if (dj_[iColumn]>maximumDJInfeasibility) {
            maximumDJInfeasibility=dj_[iColumn];
          } 
        } 
      } 
    } 
    primalObjectiveValue+=solution_[iColumn]*cost_[iColumn];
  }
  handler_->message(CLP_BARRIER_DIAGONAL,messages_)
    <<largestDiagonal<<smallestDiagonal
    <<CoinMessageEol;
#if 0
  // If diagonal wild - kill some
  if (largestDiagonal>1.0e17*smallestDiagonal) {
    double killValue =largestDiagonal*1.0e-17;
    for (int iColumn=0;iColumn<numberTotal;iColumn++) {
      if (fabs(diagonal_[iColumn])<killValue)
	diagonal_[iolumn]=0.0;
    }
  }
#endif
  // update rhs region
  multiplyAdd(deltaX_+numberColumns_,numberRows_,-1.0,rhsFixRegion_,1.0);
  matrix_->times(1.0,deltaX_,rhsFixRegion_);
  primalObjectiveValue += 0.5*gamma2*gammaOffset+0.5*quadraticOffset; 
  if (quadraticOffset) {
    //  printf("gamma offset %g %g, quadoffset %g\n",gammaOffset,gamma2*gammaOffset,quadraticOffset);
  }
  
  dualObjectiveValue+=offsetObjective+dualFake;
  dualObjectiveValue -= 0.5*gamma2*gammaOffset+0.5*quadraticOffset; 
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
  // relax dual test if obj big and gap smallish
  double gap=fabs(primalObjective_-dualObjective_);
  double sizeObj = CoinMin(fabs(primalObjective_),fabs(dualObjective_))+1.0e-50;
  //printf("gap %g sizeObj %g ratio %g comp %g\n",
  //     gap,sizeObj,gap/sizeObj,complementarityGap_);
  if (numberIterations_>100&&gap/sizeObj<1.0e-9&&complementarityGap_<1.0e-7*sizeObj)
    dualTolerance *= 1.0e2;
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
//  Save info on products of affine deltaSU*deltaW and deltaSL*deltaZ
double 
ClpPredictorCorrector::affineProduct()
{
  double product = 0.0;
  //IF zVec starts as 0 then deltaZ always zero
  //(remember if free then zVec not 0)
  //I think free can be done with careful use of boundSlacks to zero
  //out all we want
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    double w3=deltaZ_[iColumn]*deltaX_[iColumn];
    double w4=-deltaW_[iColumn]*deltaX_[iColumn];
    if (lowerBound(iColumn)) {
      w3+=deltaZ_[iColumn]*(solution_[iColumn]-lowerSlack_[iColumn]-lower_[iColumn]);
      product+=w3;
    } 
    if (upperBound(iColumn)) {
      w4+=deltaW_[iColumn]*(-solution_[iColumn]-upperSlack_[iColumn]+upper_[iColumn]);
      product+=w4;
    } 
  } 
  return product;
}
//See exactly what would happen given current deltas
void 
ClpPredictorCorrector::debugMove(int phase,double primalStep, double dualStep)
{
#ifndef SOME_DEBUG
  return;
#endif
  int numberTotal = numberRows_+numberColumns_;
  double * dualNew = ClpCopyOfArray(dual_,numberRows_);
  double * errorRegionNew = new double [numberRows_];
  double * rhsFixRegionNew = new double [numberRows_];
  double * primalNew = ClpCopyOfArray(solution_,numberTotal);
  double * djNew = new double[numberTotal];
  //update pi
  multiplyAdd(deltaY_,numberRows_,dualStep,dualNew,1.0);
  // do reduced costs
  CoinMemcpyN(dualNew,numberRows_,djNew+numberColumns_);
  CoinMemcpyN(cost_,numberColumns_,djNew);
  matrix_->transposeTimes(-1.0,dualNew,djNew);
  // update x
  int iColumn;
  for (iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) 
      primalNew[iColumn] +=primalStep*deltaX_[iColumn];
  }
  double quadraticOffset=quadraticDjs(djNew,primalNew,1.0);
  CoinZeroN(errorRegionNew,numberRows_);
  CoinZeroN(rhsFixRegionNew,numberRows_);
  double maximumBoundInfeasibility=0.0;
  double maximumDualError=1.0e-12;
  double primalObjectiveValue=0.0;
  double dualObjectiveValue=0.0;
  double solutionNorm=1.0e-12;
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
  double newGap=0.0;
  double offsetObjective=0.0;
  double gamma2 = gamma_*gamma_; // gamma*gamma will be added to diagonal
  double gammaOffset=0.0;
  double maximumDjInfeasibility=0.0;
  for ( iColumn=0;iColumn<numberTotal;iColumn++) {
    if (!flagged(iColumn)) {
      double reducedCost=djNew[iColumn];
      double zValue = zVec_[iColumn] + dualStep*deltaZ_[iColumn];
      double wValue = wVec_[iColumn] + dualStep*deltaW_[iColumn];
      double thisWeight=deltaX_[iColumn];
      double oldPrimal = solution_[iColumn];
      double newPrimal=primalNew[iColumn];
      double lowerBoundInfeasibility=0.0;
      double upperBoundInfeasibility=0.0;
      if (lowerBound(iColumn)) {
        double oldSlack = lowerSlack_[iColumn];
        double newSlack=
	  lowerSlack_[iColumn]+primalStep*(oldPrimal-oldSlack
						 + thisWeight-lower_[iColumn]);
        if (zValue>dualTolerance) {
          dualObjectiveValue+=lower_[iColumn]*zVec_[iColumn];
        } 
	lowerBoundInfeasibility = fabs(newPrimal-newSlack-lower_[iColumn]);
	newGap += newSlack*zValue;
      } 
      if (upperBound(iColumn)) {
        double oldSlack = upperSlack_[iColumn];
        double newSlack=
	  upperSlack_[iColumn]+primalStep*(-oldPrimal-oldSlack
						 - thisWeight+upper_[iColumn]);
        if (wValue>dualTolerance) {
          dualObjectiveValue-=upper_[iColumn]*wVec_[iColumn];
        } 
	upperBoundInfeasibility = fabs(newPrimal+newSlack-upper_[iColumn]);
	newGap += newSlack*wValue;
      } 
      if (fabs(newPrimal)>solutionNorm) {
        solutionNorm=fabs(newPrimal);
      } 
      double gammaTerm = gamma2;
      if (primalR_) {
	gammaTerm += primalR_[iColumn];
	quadraticOffset += newPrimal*newPrimal*primalR_[iColumn];
      }
      double dualInfeasibility=
	reducedCost-zValue+wValue+gammaTerm*newPrimal;
      if (fabs(dualInfeasibility)>dualTolerance) {
	dualFake+=newPrimal*dualInfeasibility;
      } 
      if (lowerBoundInfeasibility>maximumBoundInfeasibility) {
	maximumBoundInfeasibility=lowerBoundInfeasibility;
      } 
      if (upperBoundInfeasibility>maximumBoundInfeasibility) {
	maximumBoundInfeasibility=upperBoundInfeasibility;
      } 
      dualInfeasibility=fabs(dualInfeasibility);
      if (dualInfeasibility>maximumDualError) {
	//printf("bad dual %d %g\n",iColumn,
	// reducedCost-zVec_[iColumn]+wVec_[iColumn]+gammaTerm*newPrimal);
	maximumDualError=dualInfeasibility;
      } 
      gammaOffset += newPrimal*newPrimal;
      djNew[iColumn]=0.0;
    } else {
      offsetObjective+=primalNew[iColumn]*cost_[iColumn];
      if (upper_[iColumn]-lower_[iColumn]>1.0e-5) {
        if (primalNew[iColumn]<lower_[iColumn]+1.0e-8&&djNew[iColumn]<-1.0e-8) {
          if (-djNew[iColumn]>maximumDjInfeasibility) {
            maximumDjInfeasibility=-djNew[iColumn];
          } 
        } 
        if (primalNew[iColumn]>upper_[iColumn]-1.0e-8&&djNew[iColumn]>1.0e-8) {
          if (djNew[iColumn]>maximumDjInfeasibility) {
            maximumDjInfeasibility=djNew[iColumn];
          } 
        } 
      } 
      djNew[iColumn]=primalNew[iColumn];
    } 
    primalObjectiveValue+=solution_[iColumn]*cost_[iColumn];
  }
  // update rhs region
  multiplyAdd(djNew+numberColumns_,numberRows_,-1.0,rhsFixRegionNew,1.0);
  matrix_->times(1.0,djNew,rhsFixRegionNew);
  primalObjectiveValue += 0.5*gamma2*gammaOffset+0.5*quadraticOffset; 
  dualObjectiveValue+=offsetObjective+dualFake;
  dualObjectiveValue -= 0.5*gamma2*gammaOffset+0.5*quadraticOffset; 
  // Need to rethink (but it is only for printing)
  //compute error and fixed RHS
  multiplyAdd(primalNew+numberColumns_,numberRows_,-1.0,errorRegionNew,0.0);
  matrix_->times(1.0,primalNew,errorRegionNew);
  //finish off objective computation
  double primalObjectiveNew=primalObjectiveValue*scaleFactor_;
  double dualValue2=innerProduct(dualNew,numberRows_,
                rhsFixRegionNew);
  dualObjectiveValue-=dualValue2;
  double dualObjectiveNew=dualObjectiveValue*scaleFactor_;
  double maximumRHSError1=0.0;
  double maximumRHSError2=0.0;
  double primalOffset=0.0;
  char * dropped = cholesky_->rowsDropped();
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double value=errorRegionNew[iRow];
    if (!dropped[iRow]) {
      if (fabs(value)>maximumRHSError1) {
	maximumRHSError1=fabs(value);
      } 
    } else {
      if (fabs(value)>maximumRHSError2) {
	maximumRHSError2=fabs(value);
      } 
      primalOffset+=value*dualNew[iRow];
    } 
  } 
  primalObjectiveNew-=primalOffset*scaleFactor_;
  double maximumRHSError;
  if (maximumRHSError1>maximumRHSError2) {
    maximumRHSError=maximumRHSError1;
  } else {
    maximumRHSError=maximumRHSError1; //note change
    if (maximumRHSError2>primalTolerance()) {
      handler_->message(CLP_BARRIER_ABS_DROPPED,messages_)
	<<maximumRHSError2
	<<CoinMessageEol;
    } 
  }
  printf("PH %d %g, %g new comp %g, b %g, p %g, d %g\n",phase,
	 primalStep,dualStep,newGap,maximumBoundInfeasibility,
	 maximumRHSError,maximumDualError);
  if (handler_->logLevel()>1) 
    printf("       objs %g %g\n",
	   primalObjectiveNew,dualObjectiveNew);
  if (maximumDjInfeasibility) {
    printf(" max dj error on fixed %g\n",
	   maximumDjInfeasibility);
  } 
  delete [] dualNew;
  delete [] errorRegionNew;
  delete [] rhsFixRegionNew;
  delete [] primalNew;
  delete [] djNew;
}
