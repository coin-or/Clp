// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Implements crude primal dual predictor corrector algorithm

 */


#include "CoinPragma.hpp"
static int xxxxxx=-1;
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
static double eFree =1.0e3;

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
    return -1;
  }
  mu_=1.0e10;
  diagonalScaleFactor_=1.0;
  //set iterations
  numberIterations_=-1;
  //initialize solution here
  if(createSolution()<0) {
    printf("Not enough memory\n");
    problemStatus_=4;
    return -1;
  }
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
  int numberFixed=updateSolution();
  int numberFixedTotal=numberFixed;
  //int numberRows_DroppedBefore=0;
  //double extra=eExtra;
  //double maximumPerturbation=COIN_DBL_MAX;
  //constants for infeas interior point
  const double beta2 = 0.99995;
  const double tau   = 0.00002;
  double lastComplementarityGap=COIN_DBL_MAX;
  int lastGoodIteration=0;
  double bestObjectiveGap=COIN_DBL_MAX;
  int saveIteration=-1;
  bool sloppyOptimal=false;
  double * savePi=NULL;
  double * savePrimal=NULL;
  int numberTotal = numberRows_+numberColumns_;
  // Extra regions for centering
  double * saveWeights = new double[numberTotal];
  double * saveZ = new double[numberTotal];
  double * saveW = new double[numberTotal];
  while (problemStatus_<0) {
    if (numberIterations_==xxxxxx) {
      FILE *fp = fopen("yy.yy","wb");
      if (fp) {
	int nrow=numberRows_;
	int ncol=numberColumns_;
	// and flip
	int i;
	for (i=0;i<nrow;i++) {
	  solution_[ncol+i]= -solution_[ncol+i];
	  double value;
	  value=lowerSlack_[ncol+i];
	  lowerSlack_[ncol+i]=upperSlack_[ncol+i];
	  upperSlack_[ncol+i]=value;
	  value=zVec_[ncol+i];
	  zVec_[ncol+i]=wVec_[ncol+i];
	  wVec_[ncol+i]=value;
	}
	fwrite(dual_,sizeof(double),nrow,fp);
	fwrite(errorRegion_,sizeof(double),nrow,fp);
	fwrite(rhsFixRegion_,sizeof(double),nrow,fp);
	fwrite(solution_+ncol,sizeof(double),nrow,fp);
	fwrite(diagonal_+ncol,sizeof(double),nrow,fp);
	fwrite(wVec_+ncol,sizeof(double),nrow,fp);
	fwrite(zVec_+ncol,sizeof(double),nrow,fp);
	fwrite(upperSlack_+ncol,sizeof(double),nrow,fp);
	fwrite(lowerSlack_+ncol,sizeof(double),nrow,fp);
	fwrite(solution_,sizeof(double),ncol,fp);
	fwrite(diagonal_,sizeof(double),ncol,fp);
	fwrite(wVec_,sizeof(double),ncol,fp);
	fwrite(zVec_,sizeof(double),ncol,fp);
	fwrite(upperSlack_,sizeof(double),ncol,fp);
	fwrite(lowerSlack_,sizeof(double),ncol,fp);
	fclose(fp);
	actualDualStep_=0.0;
	actualPrimalStep_=0.0;
      } else {
	printf("no file\n");
      }
    }
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
    complementarityGap_=complementarityGap(numberComplementarityPairs_,0);
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
      } else {
        //lastComplementarityGap=complementarityGap_;
      } 
    } else if (complementarityGap_<goodGapChange*lastComplementarityGap) {
      lastGoodIteration=numberIterations_;
      lastComplementarityGap=complementarityGap_;
    } else if (numberIterations_-lastGoodIteration>=5&&complementarityGap_<1.0e-3) {
      handler_->message(CLP_BARRIER_COMPLEMENTARITY,messages_)
	<<complementarityGap_<<"not decreasing"
	<<CoinMessageEol;
      if (gapO>0.75*lastGood) {
        break;
      } 
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
    bool goodMove=false;
    double bestNextGap=COIN_DBL_MAX;
    double nextGap=bestNextGap;
    int nextNumber=0;
    worstDirectionAccuracy_=0.0;
    while (!goodMove) {
      goodMove=true;
      int newDropped=0;
      //Predictor step
      //prepare for cholesky.  Set up scaled diagonal in weights
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
      if (newDropped>0) {
	//int newDropped2=cholesky_->factorize(diagonal_,rowsDroppedThisTime);
	//assert(!newDropped2);
	if (newDropped<0) {
	  //replace dropped
	  newDropped=-newDropped;
	  //off 1 to allow for reset all
	  newDropped--;
	  //set all bits false
	  cholesky_->resetRowsDropped();
	} 
      } else if (newDropped==-1) {
	printf("Out of memory\n");
	problemStatus_=4;
	return -1;
      } 
      delete [] rowsDroppedThisTime;
      if (cholesky_->status()) {
	std::cout << "bad cholesky?" <<std::endl;
	abort();
      }
      int phase=0; // predictor, corrector , primal dual
      double directionAccuracy=0.0;
      bool doCorrector=true;
      //goodMove=false;
      if (goodMove) {
	//set up for affine direction
	setupForSolve(phase);
	directionAccuracy=findDirectionVector(phase);
	if (directionAccuracy>worstDirectionAccuracy_) {
	  worstDirectionAccuracy_=directionAccuracy;
	} 
	findStepLength(phase,NULL);
	nextGap=complementarityGap(nextNumber,1);
	
	bestNextGap=max(nextGap,0.8*complementarityGap_);
	if (complementarityGap_>1.0e-5*numberComplementarityPairs_) {
	  //std::cout <<"predicted duality gap "<<nextGap<<std::endl;
	  double part1=nextGap/numberComplementarityPairs_;
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
	  //? could use gapO?
	  //if small then should be stopping
	  //if (mu_<1.0e-4/(numberComplementarityPairs_*qqqq)) {
	  //mu_=1.0e-4/(numberComplementarityPairs_*qqqq);
	  //? better to skip corrector?
	  //} 
	} 
	//save information
	double product=affineProduct();
	//can we do corrector step?
	double xx= complementarityGap_*(beta2-tau) +product;
	//std::cout<<"gap part "<<
	//complementarityGap_*(beta2-tau)<<" after adding product = "<<xx<<std::endl;
	//xx=-1.0;
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
	//printf("product %g mu %g\n",product,mu_);
	if (complementarityGap_*(beta2-tau)+product-mu_*numberComplementarityPairs_<0.0) {
	  doCorrector=false;
	  goodMove=false;
	  bestNextGap=COIN_DBL_MAX;
	  handler_->message(CLP_BARRIER_INFO,messages_)
	    <<"no corrector step"
	    <<CoinMessageEol;
	} else {
	  phase=1;
	} 
      }
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
	  bestNextGap=COIN_DBL_MAX;
	} 
	if (goodMove) {
	  findStepLength(phase,NULL);
	  // just for debug
	  nextGap = complementarityGap(nextNumber,1);
	  if (numberIterations_>=-77) {
	    goodMove=checkGoodMove(true,bestNextGap);
	  } else {
	    goodMove=true;
	  } 
	}
      } 
      if (!goodMove) {
	// Just primal dual step
	goodMove=true;
	double floatNumber;
	floatNumber = 2.0*numberComplementarityPairs_;
	//floatNumber = 1.1*numberComplementarityPairs_;
	mu_=complementarityGap_/floatNumber;
	//set up for next step
	setupForSolve(2);
	double directionAccuracy=findDirectionVector(2);
	if (directionAccuracy>worstDirectionAccuracy_) {
	  worstDirectionAccuracy_=directionAccuracy;
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
	findStepLength(2,NULL);
	// just for debug
	nextGap=complementarityGap(nextNumber,2);
      }
    } /* endwhile */
    //numberFixed=updateSolution();
    //numberFixedTotal+=numberFixed;
    if (numberIterations_==1) {
      FILE *fp = fopen("xx.xx","rb");
      if (fp) {
	int nrow=numberRows_;
	int ncol=numberColumns_;
	fread(dual_,sizeof(double),nrow,fp);
	fread(errorRegion_,sizeof(double),nrow,fp);
	fread(rhsFixRegion_,sizeof(double),nrow,fp);
	fread(solution_+ncol,sizeof(double),nrow,fp);
	fread(diagonal_+ncol,sizeof(double),nrow,fp);
	fread(wVec_+ncol,sizeof(double),nrow,fp);
	fread(zVec_+ncol,sizeof(double),nrow,fp);
	fread(upperSlack_+ncol,sizeof(double),nrow,fp);
	fread(lowerSlack_+ncol,sizeof(double),nrow,fp);
	fread(solution_,sizeof(double),ncol,fp);
	fread(diagonal_,sizeof(double),ncol,fp);
	fread(wVec_,sizeof(double),ncol,fp);
	fread(zVec_,sizeof(double),ncol,fp);
	fread(upperSlack_,sizeof(double),ncol,fp);
	fread(lowerSlack_,sizeof(double),ncol,fp);
	fclose(fp);
	// and flip
	int i;
	for (i=0;i<nrow;i++) {
	  solution_[ncol+i]= -solution_[ncol+i];
	  double value;
	  value=lowerSlack_[ncol+i];
	  lowerSlack_[ncol+i]=upperSlack_[ncol+i];
	  upperSlack_[ncol+i]=value;
	  value=zVec_[ncol+i];
	  zVec_[ncol+i]=wVec_[ncol+i];
	  wVec_[ncol+i]=value;
	}
	actualDualStep_=0.0;
	actualPrimalStep_=0.0;
      } else {
	printf("no file\n");
      }
    }
#if 0
    // See if we can do centering steps
    memcpy(saveWeights,weights_,numberTotal*sizeof(double));
    memcpy(saveZ,deltaZ_,numberTotal*sizeof(double));
    memcpy(saveW,deltaW_,numberTotal*sizeof(double));
    double savePrimalStep = actualPrimalStep_;
    double saveDualStep = actualDualStep_;
    double saveMu = mu_;
    mu_=nextGap / (double) nextNumber;
    setupForSolve(3);
    findDirectionVector(3);
    memcpy(deltaZ_,saveZ,numberTotal*sizeof(double));
    memcpy(deltaW_,saveW,numberTotal*sizeof(double));
    findStepLength(3,saveWeights);
    double xGap = complementarityGap(nextNumber,3);
#if 0
    goodMove=checkGoodMove(true,bestNextGap);
#endif
    if (xGap>nextGap) {
      mu_=saveMu;
      actualPrimalStep_ = savePrimalStep;
      actualDualStep_ = saveDualStep;
      memcpy(weights_,saveWeights,numberTotal*sizeof(double));
      memcpy(deltaZ_,saveZ,numberTotal*sizeof(double));
      memcpy(deltaW_,saveW,numberTotal*sizeof(double));
    }
#endif
    numberFixed=updateSolution();
    numberFixedTotal+=numberFixed;
  } /* endwhile */
  delete [] saveWeights;
  delete [] saveZ;
  delete [] saveW;
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
  //std::cout<<"Absolute primal infeasibility at end "<<sumPrimalInfeasibilities_<<std::endl;
  //std::cout<<"Absolute dual infeasibility at end "<<sumDualInfeasibilities_<<std::endl;
  //std::cout<<"Absolute complementarity at end "<<complementarityGap_<<std::endl;
  //std::cout<<"Primal objective "<<objectiveValue()<<std::endl;
  //std::cout<<"maximum complementarity "<<worstComplementarity_<<std::endl;
  //delete all temporary regions
  deleteWorkingData();
  return problemStatus_;
}
// findStepLength.
//phase  - 0 predictor
//         1 corrector
//         2 primal dual
double ClpPredictorCorrector::findStepLength(const int phase,const double * oldWeight)
{
  double directionNorm=0.0;
  double maximumPrimalStep=COIN_DBL_MAX;
  double maximumDualStep=COIN_DBL_MAX;
  int numberTotal = numberRows_+numberColumns_;
  double tolerance = 1.0e-12;
  int chosenPrimalSequence=-1;
  int chosenDualSequence=-1;
  double extra=eExtra;
  double * zVec = zVec_;
  double * wVec = wVec_;
  double * primal = solution_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * work1 = deltaZ_;
  double * work2 = deltaW_;
  double * work3 = deltaS_;
  double * work4 = deltaT_;
  //direction vector in weights
  double * weights = weights_;
  int iColumn;
  switch (phase) {
  case 0:
    //Now get affine deltas for Z(duals on LBds) and W (duals on UBds)
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
        double z1=-zVec[iColumn];
        double w1=-wVec[iColumn];
        double work3Value=0.0;
        double work4Value=0.0;
        double directionElement=weights[iColumn];
        double value=primal[iColumn];
        if (directionNorm<fabs(directionElement)) {
          directionNorm=fabs(directionElement);
        } 
        //below does not feel right - can it be simplified because
        //of zero values for zVec and wVec
        if (lowerBound(iColumn)||
                        upperBound(iColumn)) {
          if (lowerBound(iColumn)) {
            double gap=lowerSlack[iColumn]+extra;
            double delta;
            if (!fakeLower(iColumn)) {
              delta = lower[iColumn]+lowerSlack[iColumn]
             -  value - directionElement;
            } else {
              delta = - directionElement;
            } 
            z1+=(zVec[iColumn]*delta)/gap;
            work3Value=-delta;
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
            double gap=upperSlack[iColumn]+extra;
            double delta;
            if (!fakeUpper(iColumn)) {
              delta = -upper[iColumn]+upperSlack[iColumn]
                        + value + directionElement;
            } else {
              delta = directionElement;
            } 
            w1+=(wVec[iColumn]*delta)/gap;
            work4Value=-delta;
            if (upperSlack[iColumn]<maximumPrimalStep*delta) {
              maximumPrimalStep=upperSlack[iColumn]/delta;
              chosenPrimalSequence=iColumn;
            } 
            if (wVec[iColumn]>tolerance) {
              if (wVec[iColumn]<-w1*maximumDualStep) {
                maximumDualStep=-wVec[iColumn]/w1;
                chosenDualSequence=iColumn;
              } 
            } 
          } 
        } else {
          //free
          double gap=fabs(value);
          double multiplier=1.0/gap;
          if (gap<1.0) {
            multiplier=1,0;
          } 
          z1-=multiplier*directionElement*zVec[iColumn];
          w1+=multiplier*directionElement*wVec[iColumn];
        } 
        work1[iColumn]=z1;
        work2[iColumn]=w1;
        work3[iColumn]=work3Value;
        work4[iColumn]=work4Value;
      } else {
        work1[iColumn]=0.0;
        work2[iColumn]=0.0;
        work3[iColumn]=0.0;
        work4[iColumn]=0.0;
      } 
    }
    break;
  case 1:
    //corrector step
    for ( iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
        double z1=-zVec[iColumn];
        double w1=-wVec[iColumn];
        double work3Value=0.0;
        double work4Value=0.0;
        double directionElement=weights[iColumn];
        double value=primal[iColumn];
        if (directionNorm<fabs(directionElement)) {
          directionNorm=fabs(directionElement);
        } 
        //below does not feel right - can it be simplified because
        //of zero values for zVec and wVec
        if (lowerBound(iColumn)||
                        upperBound(iColumn)) {
          if (lowerBound(iColumn)) {
            double gap=lowerSlack[iColumn]+extra;
            double delta;
            if (!fakeLower(iColumn)) {
              delta = lower[iColumn]+lowerSlack[iColumn]
             -  value - directionElement;
            } else {
              delta = - directionElement;
            } 
            z1+=(mu_-work3[iColumn]+zVec[iColumn]*delta)/gap;
            work3Value=-delta;
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
            double gap=upperSlack[iColumn]+extra;
            double delta;
            if (!fakeUpper(iColumn)) {
              delta = -upper[iColumn]+upperSlack[iColumn]
                        + value + directionElement;
            } else {
              delta = directionElement;
            } 
            //double delta = -upper[iColumn]+upperSlack[iColumn]-
              //+ value + directionElement;
            w1+=(mu_-work4[iColumn]+wVec[iColumn]*delta)/gap;
            work4Value=-delta;
            if (upperSlack[iColumn]<maximumPrimalStep*delta) {
              maximumPrimalStep=upperSlack[iColumn]/delta;
              chosenPrimalSequence=iColumn;
            } 
            if (wVec[iColumn]>tolerance) {
              if (wVec[iColumn]<-w1*maximumDualStep) {
                maximumDualStep=-wVec[iColumn]/w1;
                chosenDualSequence=iColumn;
              } 
            } 
          } 
        } else {
          //free
          double gap=fabs(value);
          double multiplier=1.0/gap;
          if (gap<1.0) {
            multiplier=1,0;
          } 
          z1+=multiplier*(mu_-work3[iColumn]-directionElement*zVec[iColumn]);
          w1+=multiplier*(mu_-work4[iColumn]+directionElement*wVec[iColumn]);
        } 
        work1[iColumn]=z1;
        work2[iColumn]=w1;
        work3[iColumn]=work3Value;
        work4[iColumn]=work4Value;
      } else {
        work1[iColumn]=0.0;
        work2[iColumn]=0.0;
        work3[iColumn]=0.0;
        work4[iColumn]=0.0;
      } 
    }
    break;
  case 2:
    //just primal dual
    for ( iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
        double z1=-zVec[iColumn];
        double w1=-wVec[iColumn];
        double work3Value=0.0;
        double work4Value=0.0;
        double directionElement=weights[iColumn];
        double value=primal[iColumn];
        if (directionNorm<fabs(directionElement)) {
          directionNorm=fabs(directionElement);
        } 
        //below does not feel right - can it be simplified because
        //of zero values for zVec and wVec
        if (lowerBound(iColumn)||
                        upperBound(iColumn)) {
          if (lowerBound(iColumn)) {
            double gap=lowerSlack[iColumn]+extra;
            double delta;
            if (!fakeLower(iColumn)) {
              delta = lower[iColumn]+lowerSlack[iColumn]
             -  value - directionElement;
            } else {
              delta = - directionElement;
            } 
            z1+=(mu_+zVec[iColumn]*delta)/gap;
            work3Value=-delta;
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
            double gap=upperSlack[iColumn]+extra;
            double delta;
            if (!fakeUpper(iColumn)) {
              delta = -upper[iColumn]+upperSlack[iColumn]
                        + value + directionElement;
            } else {
              delta = directionElement;
            } 
            w1+=(mu_+wVec[iColumn]*delta)/gap;
            work4Value=-delta;
            if (upperSlack[iColumn]<maximumPrimalStep*delta) {
              maximumPrimalStep=upperSlack[iColumn]/delta;
              chosenPrimalSequence=iColumn;
            } 
            if (wVec[iColumn]>tolerance) {
              if (wVec[iColumn]<-w1*maximumDualStep) {
                maximumDualStep=-wVec[iColumn]/w1;
                chosenDualSequence=iColumn;
              } 
            } 
          } 
        } else {
          //free
          double gap=fabs(value);
          double multiplier=1.0/gap;
          if (gap<1.0) {
            multiplier=1,0;
          } 
          z1+=multiplier*(mu_-directionElement*zVec[iColumn]);
          w1+=multiplier*(mu_+directionElement*wVec[iColumn]);
        } 
        work1[iColumn]=z1;
        work2[iColumn]=w1;
        work3[iColumn]=work3Value;
        work4[iColumn]=work4Value;
      } else {
        work1[iColumn]=0.0;
        work2[iColumn]=0.0;
        work3[iColumn]=0.0;
        work4[iColumn]=0.0;
      } 
    }
    break;
  case 3:
    //centering step
    for ( iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
	// old deltas from previous step
        double z1=work1[iColumn];
        double w1=work2[iColumn];
        double work3Value=0.0;
        double work4Value=0.0;
        double directionElement=weights[iColumn];
	double oldElement = oldWeight[iColumn];
	double newElement = oldElement+directionElement;
	weights[iColumn]=newElement;
        double value=primal[iColumn];
        if (directionNorm<fabs(directionElement)) {
          directionNorm=fabs(directionElement);
        } 
        if (lowerBound(iColumn)||
                        upperBound(iColumn)) {
          if (lowerBound(iColumn)) {
            double gap=lowerSlack[iColumn]+extra;
            double delta;
	    double zDelta;
            if (!fakeLower(iColumn)) {
              zDelta = lower[iColumn]+lowerSlack[iColumn]
             -  value - directionElement;
            } else {
              zDelta = - directionElement;
            } 
	    delta = zDelta - oldElement;
            z1+=(work3[iColumn]-zVec[iColumn]*zDelta)/gap;
            work3Value=-delta;
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
            double gap=upperSlack[iColumn]+extra;
            double delta;
	    double zDelta;
            if (!fakeUpper(iColumn)) {
              zDelta = -upper[iColumn]+upperSlack[iColumn]
                        + value + directionElement;
            } else {
              zDelta = directionElement;
            } 
	    delta = zDelta + oldElement;
            //double delta = -upper[iColumn]+upperSlack[iColumn]-
              //+ value + directionElement;
            w1+=(work4[iColumn]+wVec[iColumn]*zDelta)/gap;
            work4Value=-delta;
            if (upperSlack[iColumn]<maximumPrimalStep*delta) {
              maximumPrimalStep=upperSlack[iColumn]/delta;
              chosenPrimalSequence=iColumn;
            } 
            if (wVec[iColumn]>tolerance) {
              if (wVec[iColumn]<-w1*maximumDualStep) {
                maximumDualStep=-wVec[iColumn]/w1;
                chosenDualSequence=iColumn;
              } 
            } 
          } 
        } else {
          //free
          double gap=fabs(value);
          double multiplier=1.0/gap;
          if (gap<1.0) {
            multiplier=1,0;
          } 
          z1+=multiplier*(work3[iColumn]-directionElement*zVec[iColumn]);
          w1+=multiplier*(work4[iColumn]+directionElement*wVec[iColumn]);
        } 
        work1[iColumn]=z1;
        work2[iColumn]=w1;
        work3[iColumn]=work3Value;
        work4[iColumn]=work4Value;
      } else {
        work1[iColumn]=0.0;
        work2[iColumn]=0.0;
        work3[iColumn]=0.0;
        work4[iColumn]=0.0;
      } 
    }
    break;
  } 
  //printf("step - phase %d, norm %g, dual step %g, primal step %g\n",
  // phase,directionNorm,maximumPrimalStep,maximumDualStep);
  actualPrimalStep_=stepLength_*maximumPrimalStep;
  if (phase>=0&&actualPrimalStep_>1.0) {
    actualPrimalStep_=1.0;
  } 
  actualDualStep_=stepLength_*maximumDualStep;
  if (phase>=0&&actualDualStep_>1.0) {
    actualDualStep_=1.0;
  } 
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
  if (phase==3) {
    // centering step
    double * work2 = deltaW_;
    int numberTotal = numberRows_+numberColumns_;
    double * tempRegion = new double [numberRows_];
    double * newError = new double [numberRows_];
    // For KKT separate out
#ifndef KKT
    int iColumn;
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      weights_[iColumn] = work2[iColumn];
    multiplyAdd(weights_+numberColumns_,numberRows_,-1.0,tempRegion,0.0);
    matrix_->times(1.0,weights_,tempRegion);
#else
    // regions in will be work2 and newError
    // regions out weights_ and tempRegion
    memset(newError,0,numberRows_*sizeof(double));
    double * region1Save=NULL;//for refinement
#endif
#ifndef KKT
    double maximumRHS = maximumAbsElement(tempRegion,numberRows_);
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
    multiplyAdd(NULL,numberRows_,0.0,tempRegion,scale);
    cholesky_->solve(tempRegion);
    multiplyAdd(NULL,numberRows_,0.0,tempRegion,unscale);
    //CoinZeroN(newError,numberRows_);
    multiplyAdd(tempRegion,numberRows_,-1.0,weights_+numberColumns_,0.0);
    CoinZeroN(weights_,numberColumns_);
    matrix_->transposeTimes(1.0,tempRegion,weights_);
    //if flagged then entries zero so can do
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      weights_[iColumn] = weights_[iColumn]*diagonal_[iColumn]
	-work2[iColumn];
#else
    solveSystem(weights_, tempRegion,
		work2,newError,region1Save,regionSave,lastError>1.0e-5);
#endif
    multiplyAdd(weights_+numberColumns_,numberRows_,-1.0,newError,0.0);
    matrix_->times(1.0,weights_,newError);
    delete [] newError;
    delete [] tempRegion;
    return 0.0;
  }
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
  double * work2 = deltaW_;
  int numberTotal = numberRows_+numberColumns_;
  //if flagged then entries zero so can do
  // For KKT separate out
#ifndef KKT
  int iColumn;
  for (iColumn=0;iColumn<numberTotal;iColumn++)
    weights_[iColumn] = work2[iColumn] - solution_[iColumn];
  multiplyAdd(weights_+numberColumns_,numberRows_,-1.0,updateRegion_,0.0);
  matrix_->times(1.0,weights_,updateRegion_);
#else
  // regions in will be work2 and newError
  // regions out weights_ and updateRegion_
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
#ifndef KKT
    double maximumRHS = maximumAbsElement(updateRegion_,numberRows_);
    double saveMaximum = maximumRHS;
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
    multiplyAdd(NULL,numberRows_,0.0,updateRegion_,scale);
    cholesky_->solve(updateRegion_);
    multiplyAdd(NULL,numberRows_,0.0,updateRegion_,unscale);
    if (numberTries) {
      //refine?
      double scaleX=1.0;
      if (lastError>1.0e-5) 
        scaleX=0.8;
      multiplyAdd(regionSave,numberRows_,1.0,updateRegion_,scaleX);
    } 
    //CoinZeroN(newError,numberRows_);
    multiplyAdd(updateRegion_,numberRows_,-1.0,weights_+numberColumns_,0.0);
    CoinZeroN(weights_,numberColumns_);
    matrix_->transposeTimes(1.0,updateRegion_,weights_);
    //if flagged then entries zero so can do
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      weights_[iColumn] = weights_[iColumn]*diagonal_[iColumn]
	-work2[iColumn];
#else
    solveSystem(weights_, updateRegion_,
		work2,newError,region1Save,regionSave,lastError>1.0e-5);
#endif
    multiplyAdd(weights_+numberColumns_,numberRows_,-1.0,newError,0.0);
    matrix_->times(1.0,weights_,newError);
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
	//assert(updateRegion_[iRow]==0.0);
	updateRegion_[iRow]=0.0;
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
        CoinMemcpyN(updateRegion_,numberRows_,regionSave);
#ifndef KKT
	multiplyAdd(newError,numberRows_,-1.0,updateRegion_,0.0);
#else
        CoinMemcpyN(weights_,numberTotal,region1Save);
	// and back to input region
        CoinMemcpyN(updateRegion_,numberRows_,newError);
#endif
      } 
    } else {
      //std::cout <<" worse residual = "<<relativeError;
      //bring back previous
      relativeError=lastError;
      CoinMemcpyN(regionSave,numberRows_,updateRegion_);
#ifndef KKT
      multiplyAdd(updateRegion_,numberRows_,-1.0,weights_+numberColumns_,0.0);
      CoinZeroN(weights_,numberColumns_);
      matrix_->transposeTimes(1.0,updateRegion_,weights_);
      //if flagged then entries zero so can do
      for (iColumn=0;iColumn<numberTotal;iColumn++)
	weights_[iColumn] = weights_[iColumn]*diagonal_[iColumn]
	  -work2[iColumn];
#else
        CoinMemcpyN(region1Save,numberTotal,weights_);
#endif
    } 
  } /* endwhile */
  delete [] regionSave;
#ifdef KKT
  delete [] region1Save;
#endif
  delete [] newError;
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
  //use weights region for work region
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
      weights_[iColumn]=1.0;
      double lowerValue=lower_[iColumn];
      double upperValue=upper_[iColumn];
#if 0
      //fake it
      if (lowerValue<-fakeCheck&&upperValue>fakeCheck) {
        lowerValue=-1.0e5;
        upperValue=1.0e5;
        lower_[iColumn]=lowerValue;
        upper_[iColumn]=upperValue;
        std::cout<<"faking free variable "<<iColumn<<std::endl;
      } 
#endif
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
      weights_[iColumn]=0.0;
    } 
  } 
  //   modify fixed RHS
  multiplyAdd(dj_+numberColumns_,numberRows_,1.0,rhsFixRegion_,0.0);
  //   create plausible RHS?
  matrix_->times(-1.0,dj_,rhsFixRegion_);
  multiplyAdd(solution_+numberColumns_,numberRows_,-1.0,errorRegion_,0.0);
  matrix_->times(1.0,solution_,errorRegion_);
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
  multiplyAdd(errorRegion_,numberRows_,-1.0,weights_+numberColumns_,0.0);
  CoinZeroN(weights_,numberColumns_);
  matrix_->transposeTimes(1.0,errorRegion_,weights_);
#else
  // reverse sign on solution
  multiplyAdd(NULL,numberRows_+numberColumns_,0.0,solution_,-1.0);
  solveSystem(weights_,errorRegion_,solution_,NULL,NULL,NULL,false);
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
  double * fakeSolution = weights_;
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
  //safeObjectiveValue=1.0+objectiveNorm_;
  //printf("temp safe obj value of %g\n",safeObjectiveValue);
  double safeFree=1.0e-1*initialValue;
  safeFree=1.0;
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
            wVec_[iColumn]=safeObjectiveValue*ratioW;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            wVec_[iColumn]=-objectiveValue + safeObjectiveValue*ratioW;
          } 
          diagonal_[iColumn] = (t*s)/(s*wVec_[iColumn]+t*zVec_[iColumn]);
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
            wVec_[iColumn]=0.0;
          } else {
            zVec_[iColumn]=safeObjectiveValue*ratioZ;
            wVec_[iColumn]=0.0;
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
            wVec_[iColumn]=safeObjectiveValue*ratioW;
          } else {
            zVec_[iColumn]=0.0;
            wVec_[iColumn]=-objectiveValue + safeObjectiveValue*ratioW;
          } 
          diagonal_[iColumn] =  t/wVec_[iColumn];
        } else {
          //free
          zVec_[iColumn]=safeObjectiveValue;
          wVec_[iColumn]=safeObjectiveValue;
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
            diagonal_[iColumn]=fabs(newValue)/(zVec_[iColumn]+wVec_[iColumn]);
          } else {
            diagonal_[iColumn]=1.0/(zVec_[iColumn]+wVec_[iColumn]);
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
      wVec_[iColumn]=0.0;
      diagonal_[iColumn]=0.0;
    } 
  } 
  solutionNorm_ =  maximumAbsElement(solution_,numberTotal);
  return 0;
}
// complementarityGap.  Computes gap
//phase 0=as is , 1 = after predictor , 2 after corrector
double ClpPredictorCorrector::complementarityGap(int & numberComplementarityPairs,
			  const int phase)
{
  double gap=0.0;
  //seems to be same coding for phase = 1 or 2
  numberComplementarityPairs=0;
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
  double * wVec = wVec_;
  double * primal = solution_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * work1 = deltaZ_;
  double * work2 = deltaW_;
  double * weights = weights_;
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!fixedOrFree(iColumn)) {
      numberComplementarityPairs++;
      //can collapse as if no lower bound both zVec and work1 0.0
      if (lowerBound(iColumn)) {
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=zVec[iColumn];
          primalValue=lowerSlack[iColumn];
        } else {
          double change;
          if (!fakeLower(iColumn)) {
            change =primal[iColumn]+weights[iColumn]-lowerSlack[iColumn]-lower[iColumn];
          } else {
            change =weights[iColumn];
          } 
          dualValue=zVec[iColumn]+actualDualStep_*work1[iColumn];
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
        double dualValue;
        double primalValue;
        if (!phase) {
          dualValue=wVec[iColumn];
          primalValue=upperSlack[iColumn];
        } else {
          double change;
          if (!fakeUpper(iColumn)) {
            change =upper[iColumn]-primal[iColumn]-weights[iColumn]-upperSlack[iColumn];
          } else {
            change =-weights[iColumn];
          } 
          dualValue=wVec[iColumn]+actualDualStep_*work2[iColumn];
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
  double * wVec = wVec_;
  double * primal = solution_;
  double * dual = dj_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * work2 = deltaW_;
  double * work3 = deltaS_;
  double * work4 = deltaT_;

  int iColumn;
  //printf("phase %d in setupForSolve, mu %g\n",phase,mu_);
  switch (phase) {
  case 0:
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
        double value= dual[iColumn];
        if (lowerBound(iColumn)) {
          if (!fakeLower(iColumn)) {
            value+=zVec[iColumn]*
              (primal[iColumn]-lowerSlack[iColumn]-lower[iColumn])/
                 (lowerSlack[iColumn]+extra);
          } 
        } 
        if (upperBound(iColumn)) {
          if (!fakeUpper(iColumn)) {
            value+=wVec[iColumn]*
              (primal[iColumn]+upperSlack[iColumn]-upper[iColumn])/
                 (upperSlack[iColumn]+extra);
          } 
        } 
#ifndef KKT
        work2[iColumn]=diagonal_[iColumn]*value;
#else
        work2[iColumn]=value;
#endif
      } else {
        work2[iColumn]=0.0;
      } 
    } 
    break;
  case 1:
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
        double value= 0.0;
        if (lowerBound(iColumn)) {
          if (!fakeLower(iColumn)) {
            value-=(mu_-work3[iColumn]-zVec[iColumn]*
              (primal[iColumn]-lowerSlack[iColumn]-lower[iColumn]))/
                 (lowerSlack[iColumn]+extra);
          } else {
            value-=(mu_-work3[iColumn])/(lowerSlack[iColumn]+extra);
          } 
        } 
        if (upperBound(iColumn)) {
          if (!fakeUpper(iColumn)) {
            value+=(mu_-work4[iColumn]+wVec[iColumn]*
              (primal[iColumn]+upperSlack[iColumn]-upper[iColumn]))/
                 (upperSlack[iColumn]+extra);
          } else {
            value+=(mu_-work4[iColumn])/(upperSlack[iColumn]+extra);
          } 
        } 
#ifndef KKT
        work2[iColumn]=diagonal_[iColumn]*(dual[iColumn]+value);
#else
        work2[iColumn]=dual[iColumn]+value;
#endif
      } else {
        work2[iColumn]=0.0;
      } 
    } 
    break;
  case 2:
    for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
      if (!flagged(iColumn)) {
        double value= 0.0;
        if (lowerBound(iColumn)) {
          if (!fakeLower(iColumn)) {
            value-=(mu_-zVec[iColumn]*
              (primal[iColumn]-lowerSlack[iColumn]-lower[iColumn]))/
                 (lowerSlack[iColumn]+extra);
          } else {
            value-=(mu_-zVec[iColumn])/ (lowerSlack[iColumn]+extra);
          } 
        } 
        if (upperBound(iColumn)) {
          if (!fakeUpper(iColumn)) {
            value+=(mu_+wVec[iColumn]*
              (primal[iColumn]+upperSlack[iColumn]-upper[iColumn]))/
                 (upperSlack[iColumn]+extra);
          } else {
            value+=(mu_+wVec[iColumn])/ (upperSlack[iColumn]+extra);
          } 
        } 
#ifndef KKT
        work2[iColumn]=diagonal_[iColumn]*(dual[iColumn]+value);
#else
        work2[iColumn]=dual[iColumn]+value;
#endif
      } else {
        work2[iColumn]=0.0;
      } 
    } 
    break;
  case 3:
    {
      double minBeta = 0.1*mu_;
      double maxBeta = 10.0*mu_;
      double * work1 = deltaZ_;
      double * weights = weights_;
      for (iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
	work3[iColumn]=0.0;
	work4[iColumn]=0.0;
	if (!flagged(iColumn)) {
	  double value= 0.0;
	  if (lowerBound(iColumn)) {
	    double change;
	    if (!fakeLower(iColumn)) {
	      change =primal[iColumn]+weights[iColumn]-lowerSlack[iColumn]-lower[iColumn];
	    } else {
	      change =weights[iColumn];
	    } 
	    double dualValue=zVec[iColumn]+actualDualStep_*work1[iColumn];
	    double primalValue=lowerSlack[iColumn]+actualPrimalStep_*change;
	    double gapProduct=dualValue*primalValue;
	    if (gapProduct<minBeta) {
	      value-= (minBeta-gapProduct)/(lowerSlack[iColumn]+extra);
	    } else if (gapProduct>maxBeta) {
	      value-= max(maxBeta-gapProduct,-maxBeta)/(lowerSlack[iColumn]+extra);
	    }
	    value=-value;
	    work3[iColumn]=value;
	  }  
	  if (upperBound(iColumn)) {
	    double change;
	    if (!fakeUpper(iColumn)) {
	      change =upper[iColumn]-primal[iColumn]-weights[iColumn]-upperSlack[iColumn];
	    } else {
	      change =-weights[iColumn];
	    } 
	    double dualValue=wVec[iColumn]+actualDualStep_*work2[iColumn];
	    double primalValue=upperSlack[iColumn]+actualPrimalStep_*change;
	    double gapProduct=dualValue*primalValue;
	    if (gapProduct<minBeta) {
	      value+= (minBeta-gapProduct)/(upperSlack[iColumn]+extra);
	    } else if (gapProduct>maxBeta) {
	      value+= max(maxBeta-gapProduct,-maxBeta)/(upperSlack[iColumn]+extra);
	    }
	    value=-value;
	    work4[iColumn]=value-work3[iColumn];
	  } 
#ifndef KKT
	  work2[iColumn]=diagonal_[iColumn]*value;
#else
	  work2[iColumn]=value;
#endif
	} else {
	  work2[iColumn]=0.0;
	} 
      }
    }
    break;
  } /* endswitch */
}
//method: sees if looks plausible change in complementarity
bool ClpPredictorCorrector::checkGoodMove(const bool doCorrector,double & bestNextGap)
{
  const double beta3 = 0.99997;
  bool goodMove=false;
  int nextNumber;
  int numberTotal = numberRows_+numberColumns_;
  double returnGap=bestNextGap;
  double nextGap=complementarityGap(nextNumber,2);
  if (nextGap>bestNextGap)
    return false;
  else
    returnGap=nextGap;
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
      innerProduct(updateRegion_,numberRows_,
		   rhsFixRegion_);
    double error=0.0;
    double * work3 = deltaS_;
    CoinZeroN(work3,numberColumns_);
    CoinMemcpyN(updateRegion_,numberRows_,work3+numberColumns_);
    matrix_->transposeTimes(-1.0,updateRegion_,work3);
    double * work1 = deltaZ_;
    double * work2 = deltaW_;
    double * lower = lower_;
    double * upper = upper_;
    //direction vector in weights
    double * weights = weights_;
    double * cost = cost_;
    //double sumPerturbCost=0.0;
    for (int iColumn=0;iColumn<numberTotal;iColumn++) {
      if (!flagged(iColumn)) {
	if (lowerBound(iColumn)) {
	  //sumPerturbCost+=weights[iColumn];
	  deltaObjectiveDual+=work1[iColumn]*lower[iColumn];
	} 
	if (upperBound(iColumn)) {
	  //sumPerturbCost-=weights[iColumn];
	  deltaObjectiveDual-=work2[iColumn]*upper[iColumn];
	} 
	double change = fabs(work3[iColumn]-work1[iColumn]+work2[iColumn]);
	error = max (change,error);
      } 
      deltaObjectivePrimal += cost[iColumn] * weights[iColumn];
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
  double nextGap=complementarityGap(nextNumber,2);
  if (nextGap>bestNextGap)
    return false;
  double lowerBoundGap = gamma*nextGap*complementarityMultiplier;
  bool goodMove=true;
  double * deltaZ = deltaZ_;
  double * deltaW = deltaW_;
  double * deltaS = deltaS_;
  double * deltaT = deltaT_;
  double * zVec = zVec_;
  double * wVec = wVec_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!flagged(iColumn)) {
      if (lowerBound(iColumn)) {
	double part1=lowerSlack[iColumn]+actualPrimalStep_*deltaS[iColumn];
	double part2=zVec[iColumn]+actualDualStep_*deltaZ[iColumn];
	if (part1*part2<lowerBoundGap) {
	  goodMove=false;
	  break;
	} 
      } 
      if (upperBound(iColumn)) {
	double part1=upperSlack[iColumn]+actualPrimalStep_*deltaT[iColumn];
	double part2=wVec[iColumn]+actualDualStep_*deltaW[iColumn];
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
int ClpPredictorCorrector::updateSolution()
{
  //update pi
  multiplyAdd(updateRegion_,numberRows_,actualDualStep_,dual_,1.0);
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
  double condition = cholesky_->choleskyCondition();
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
  double largeGap=1.0e2*solutionNorm_;
  if (largeGap<1.0e2) {
    largeGap=1.0e2;
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
  double norm=1.0e-12;
  double widenGap=1.0e1;
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
  double * wVec = wVec_;
  double * primal = solution_;
  double * dual = dj_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * work1 = deltaZ_;
  double * work2 = deltaW_;
  double * diagonal = diagonal_;
  //direction vector in weights
  double * weights = weights_;
  double * cost = cost_;
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    if (!flagged(iColumn)) {
      double reducedCost=dual[iColumn];
      bool thisKilled=false;
      double zValue = zVec[iColumn] + actualDualStep_*work1[iColumn];
      double wValue = wVec[iColumn] + actualDualStep_*work2[iColumn];
      zVec[iColumn]=zValue;
      wVec[iColumn]=wValue;
      double thisWeight=weights[iColumn];
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
      if (zValue>wValue) {
        largerzw=zValue;
      } else {
        largerzw=wValue;
      } 
      if (lowerBound(iColumn)) {
        double oldSlack = lowerSlack[iColumn];
        double newSlack;
        if (!fakeLower(iColumn)) {
          newSlack=
                   lowerSlack[iColumn]+actualPrimalStep_*(oldPrimal-oldSlack
                + thisWeight-lower[iColumn]);
        } else {
          newSlack= lowerSlack[iColumn]+actualPrimalStep_*thisWeight;
          if (newSlack<0.0) {
            abort();
          } 
        } 
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
          if (!fakeLower(iColumn)) {
            double larger;
            if (newSlack>feasibleSlack) {
              larger=newSlack;
            } else {
              larger=feasibleSlack;
            } 
            if (fabs(feasibleSlack-newSlack)<1.0e-6*larger) {
              newSlack=feasibleSlack;
            } 
            //set FAKE here
            if (newSlack>zValue2*largestRatio&&newSlack>smallGap2) {
              setFakeLower(iColumn);
              newSlack=zValue2*largestRatio;
              if (newSlack<smallGap2) {
                newSlack=smallGap2;
              } 
              numberDecreased++;
            } 
          } else {
            newSlack=zValue2*largestRatio;
            if (newSlack<smallGap2) {
              newSlack=smallGap2;
            } 
            if (newSlack>largeGap) {
              //increase up to smaller of z.. and largeGap
              newSlack=largeGap;
            } 
            if (newSlack>widenGap*oldSlack) {
              newSlack=widenGap*oldSlack;
              numberIncreased++;
              //cout<<"wider "<<j<<" "<<newSlack<<" "<<oldSlack<<" "<<feasibleSlack<<endl;
            } 
            if (newSlack>feasibleSlack) {
              newSlack=feasibleSlack;
              clearFakeLower(iColumn);
            } 
          } 
        } 
        if (zVec[iColumn]>dualTolerance) {
          dualObjectiveThis+=lower[iColumn]*zVec[iColumn];
        } 
        lowerSlack[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
        if (!fakeLower(iColumn)) {
          double infeasibility = fabs(newPrimal-lowerSlack[iColumn]-lower[iColumn]);
          if (infeasibility>maximumBoundInfeasibility) {
            maximumBoundInfeasibility=infeasibility;
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
        if (!fakeUpper(iColumn)) {
          newSlack=
                   upperSlack[iColumn]+actualPrimalStep_*(-oldPrimal-oldSlack
                - thisWeight+upper[iColumn]);
        } else {
          newSlack= upperSlack[iColumn]-actualPrimalStep_*thisWeight;
        } 
        double epsilon = fabs(newSlack)*epsilonBase;
        if (epsilon>1.0e-5) {
          //cout<<"bad"<<endl;
          epsilon=1.0e-5;
        } 
        //for changing slack
        double wValue2 = wValue;
        if (wValue2<epsilonBase) {
          wValue2=epsilonBase;
        } 
        //make sure reasonable
        if (wValue<epsilon) {
          wValue=epsilon;
        } 
        //store modified wVec
        //no wVec[iColumn]=wValue;
        double feasibleSlack=upper[iColumn]-newPrimal;
        if (feasibleSlack>0.0&&newSlack>0.0) {
          double smallGap2=smallGap;
          if (fabs(0.1*newPrimal)>smallGap) {
            smallGap2=0.1*fabs(newPrimal);
          } 
          if (!fakeUpper(iColumn)) {
            double larger;
            if (newSlack>feasibleSlack) {
              larger=newSlack;
            } else {
              larger=feasibleSlack;
            } 
            if (fabs(feasibleSlack-newSlack)<1.0e-6*larger) {
              newSlack=feasibleSlack;
            } 
            //set FAKE here
            if (newSlack>wValue2*largestRatio&&newSlack>smallGap2) {
              setFakeUpper(iColumn);
              newSlack=wValue2*largestRatio;
              if (newSlack<smallGap2) {
                newSlack=smallGap2;
              } 
              numberDecreased++;
            } 
          } else {
            newSlack=wValue2*largestRatio;
            if (newSlack<smallGap2) {
              newSlack=smallGap2;
            } 
            if (newSlack>largeGap) {
              //increase up to smaller of w.. and largeGap
              newSlack=largeGap;
            } 
            if (newSlack>widenGap*oldSlack) {
              numberIncreased++;
              //cout<<"wider "<<-j<<" "<<newSlack<<" "<<oldSlack<<" "<<feasibleSlack<<endl;
            } 
            if (newSlack>feasibleSlack) {
              newSlack=feasibleSlack;
              clearFakeUpper(iColumn);
            } 
          } 
        } 
        if (wVec[iColumn]>dualTolerance) {
          dualObjectiveThis-=upper[iColumn]*wVec[iColumn];
        } 
        upperSlack[iColumn]=newSlack;
        if (newSlack<smallerSlack) {
          smallerSlack=newSlack;
        } 
        if (!fakeUpper(iColumn)) {
          double infeasibility = fabs(newPrimal+upperSlack[iColumn]-upper[iColumn]);
          if (infeasibility>maximumBoundInfeasibility) {
            maximumBoundInfeasibility=infeasibility;
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
      if ((!lowerBound(iColumn))&&
               (!upperBound(iColumn))) {
        double gap;
        if (fabs(newPrimal)>eFree) {
          gap=fabs(newPrimal);
        } else {
          gap=eFree;
        } 
        zVec[iColumn]=gap*freeMultiplier;
        wVec[iColumn]=zVec[iColumn];
        //fake to give correct result
        s=1.0;
        t=1.0;
        wValue=0.0;
        zValue=2.0*freeMultiplier+qDiagonal;
      } 
      primal[iColumn]=newPrimal;
      if (fabs(newPrimal)>norm) {
        norm=fabs(newPrimal);
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
        double divisor = s*wValue+t*zValue;
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
        double dualInfeasibility=reducedCost-zVec[iColumn]+wVec[iColumn];
        if (fabs(dualInfeasibility)>dualTolerance) {
          if ((!lowerBound(iColumn))&&
                 (!upperBound(iColumn))) {
            dualFake+=newPrimal*dualInfeasibility;
          } else {
            dualFake+=newPrimal*dualInfeasibility;
          } 
        } 
        dualInfeasibility=fabs(dualInfeasibility);
        if (dualInfeasibility>maximumDualError) {
          maximumDualError=dualInfeasibility;
        } 
        work1[iColumn]=0.0;
      } else {
        numberKilled++;
        diagonal[iColumn]=0.0;
        zVec[iColumn]=0.0;
        wVec[iColumn]=0.0;
        setFlagged(iColumn);
        setFixedOrFree(iColumn);
        work1[iColumn]=newPrimal;
        offsetObjective+=newPrimal*cost[iColumn];
      } 
    } else {
      work1[iColumn]=primal[iColumn];
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
  multiplyAdd(work1+numberColumns_,numberRows_,-1.0,rhsFixRegion_,1.0);
  matrix_->times(1.0,work1,rhsFixRegion_);
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
  solutionNorm_ = norm;
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
  double * work1 = deltaZ_;
  double * work2 = deltaW_;
  double * work3 = deltaS_;
  double * work4 = deltaT_;
  double * lower = lower_;
  double * upper = upper_;
  double * lowerSlack = lowerSlack_;
  double * upperSlack = upperSlack_;
  double * primal = solution_;
  //direction vector in weights
  double * weights = weights_;
  double product = 0.0;
  //IF zVec starts as 0 then work1 always zero
  //(remember if free then zVec not 0)
  //I think free can be done with careful use of boundSlacks to zero
  //out all we want
  for (int iColumn=0;iColumn<numberRows_+numberColumns_;iColumn++) {
    double w3=work1[iColumn]*weights[iColumn];
    double w4=-work2[iColumn]*weights[iColumn];
    if (lowerBound(iColumn)) {
      if (!fakeLower(iColumn)) {
	w3+=work1[iColumn]*(primal[iColumn]-lowerSlack[iColumn]-lower[iColumn]);
      } 
      product+=w3;
    } 
    if (upperBound(iColumn)) {
      if (!fakeUpper(iColumn)) {
	w4+=work2[iColumn]*(-primal[iColumn]-upperSlack[iColumn]+upper[iColumn]);
      } 
      product+=w4;
    } 
    work3[iColumn]=w3;
    work4[iColumn]=w4;
  } 
  return product;
}
