// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

// This file has higher level solve functions


#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpMessage.hpp"

#include "ClpPresolve.hpp"
#ifdef CLP_IDIOT
#include "Idiot.hpp"
#endif
//#############################################################################
// Allow for interrupts
// But is this threadsafe ? (so switched off by option
#include <signal.h>
static ClpSimplex * currentModel = NULL;
static void signal_handler(int whichSignal)
{
  if (currentModel!=NULL) 
    currentModel->setMaximumIterations(0); // stop at next iterations
  return;
}

/** General solve algorithm which can do presolve
    special options (bits)
    1 - do not perturb
    2 - do not scale
    4 - use crash (default allslack in dual, idiot in primal)
    8 - all slack basis in primal
    16 - switch off interrupt handling
    32 - do not try and make plus minus one matrix
 */
int 
ClpSimplex::initialSolve(SolveType method, PresolveType presolve,
			 int specialOption)
{
  int saveMaxIterations = maximumIterations();
  int finalStatus=-1;
  int numberIterations=0;
  double time1 = CoinCpuTime();
  double timeX = time1;
  double time2;
  ClpMatrixBase * saveMatrix=NULL;
  ClpSimplex * model2 = this;
  bool interrupt = ((specialOption&16)==0);
  sighandler_t saveSignal=SIG_DFL;
  if (interrupt) {
    currentModel = model2;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
  }
  ClpPresolve pinfo;
  double timePresolve=0.0;
  double timeIdiot=0.0;
  double timeCore=0.0;
  double timeSimplex=0.0;
  assert (method!=automatic); // later
  int savePerturbation=perturbation_;
  if ((specialOption&1)!=0||method==usePrimal)
    perturbation_=100;
  int saveScaling = scalingFlag_;
  if ((specialOption&2)!=0)
    scalingFlag_=0;

  if (presolve!=presolveOff) {
    int numberPasses=5;
    if (presolve==presolveMaximum)
      numberPasses=25;
    model2 = pinfo.presolvedModel(*this,1.0e-8,
				  false,numberPasses,true);
    time2 = CoinCpuTime();
    timePresolve = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Presolve"<<timePresolve<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    if (model2) {
      if (method==useDual) {
	int numberInfeasibilities = model2->tightenPrimalBounds();
	if (numberInfeasibilities) {
	  handler_->message(CLP_INFEASIBLE,messages_)
	    <<CoinMessageEol;
	  model2 = this;
	  presolve=presolveOff;
	}
      }
    } else {
      handler_->message(CLP_INFEASIBLE,messages_)
	<<CoinMessageEol;
      model2 = this;
      presolve=presolveOff;
    }
  }
  if (interrupt)
    currentModel = model2;
  // See if worth trying +- one matrix
  bool plusMinus=false;
  if ((specialOption&32)==0) {
    int numberElements=model2->getNumElements();
    if(numberElements>100000)
      plusMinus=true;
    if(numberElements>10000&&(specialOption&12)==0&&method==usePrimal) 
      plusMinus=true;
  }
  if (plusMinus) {
    saveMatrix = model2->clpMatrix();
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
    if (clpMatrix) {
      ClpPlusMinusOneMatrix * newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
      if (newMatrix->getIndices()) {
	model2->replaceMatrix(newMatrix);
      } else {
	handler_->message(CLP_MATRIX_CHANGE,messages_)
	  <<"+- 1"
	  <<CoinMessageEol;
	saveMatrix=NULL;
	plusMinus=false;
	delete newMatrix;
      }
    } else {
      saveMatrix=NULL;
      plusMinus=false;
    }
  }
  if (model2->factorizationFrequency()==200) {
    // User did not touch preset
    model2->setFactorizationFrequency(100+model2->numberRows()/100);
  }
  int numberColumns = model2->numberColumns();
  int numberRows = model2->numberRows();
  if (method==useDual) {
    if ((specialOption&4)!=0)
      model2->crash(1000,1);
    model2->dual();
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Dual"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  } else {
#ifdef CLP_IDIOT
    if ((specialOption&12)==0) {
      int nPasses=0;
      if (numberRows>2000&&numberColumns>2*numberRows) {
	if (plusMinus) {
	  nPasses = 10+numberColumns/1000;
	  nPasses = min(nPasses,100);
	} else {
	  nPasses = 10+numberColumns/100000;
	  if (numberColumns>4*numberRows) 
	    nPasses = min(nPasses,50);
	  else
	    nPasses=5;
	}
      }
      if (nPasses) {
	Idiot info(*model2);
	info.crash(nPasses,model2->messageHandler(),model2->messagesPointer());
	time2 = CoinCpuTime();
	timeIdiot = time2-timeX;
	handler_->message(CLP_INTERVAL_TIMING,messages_)
	  <<"Idiot Crash"<<timeIdiot<<time2-time1
	  <<CoinMessageEol;
	timeX=time2;
      }
    }
#endif
    if ((specialOption&4)!=0)
      model2->crash(1000,1);
    model2->primal(1);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    timeSimplex = timeCore;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Primal"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  }
  if (saveMatrix) {
    // delete and replace
    delete model2->clpMatrix();
    model2->replaceMatrix(saveMatrix);
  }
  numberIterations = model2->numberIterations();
  finalStatus=model2->status();
  if (presolve==presolveOn) {
    int saveLevel = logLevel();
    setLogLevel(1);
    pinfo.postsolve(true);
    time2 = CoinCpuTime();
    timePresolve += time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Postsolve"<<time2-timeX<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    delete model2;
    if (interrupt)
      currentModel = this;
    checkSolution();
    setLogLevel(saveLevel);
    if (finalStatus!=3&&(finalStatus||status()==-1)) {
      int savePerturbation = perturbation();
      setPerturbation(100);
      primal(1);
      setPerturbation(savePerturbation);
      numberIterations += this->numberIterations();
      finalStatus=status();
      time2 = CoinCpuTime();
      timeSimplex += time2-timeX;
      handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Cleanup"<<time2-timeX<<time2-time1
      <<CoinMessageEol;
      timeX=time2;
    }
  }
  setMaximumIterations(saveMaxIterations);
  std::string statusMessage[]={"Unknown","Optimal","PrimalInfeasible","DualInfeasible","Stopped"};
  assert (finalStatus>=-1&&finalStatus<=3);
  handler_->message(CLP_TIMING,messages_)
    <<statusMessage[finalStatus+1]<<objectiveValue()<<numberIterations<<time2-time1;
  handler_->printing(presolve==presolveOn)
    <<timePresolve;
  handler_->printing(timeIdiot)
    <<timeIdiot;
  handler_->message()<<CoinMessageEol;
  if (interrupt) 
    signal(SIGINT,saveSignal);
  perturbation_=savePerturbation;
  scalingFlag_=saveScaling;
  return finalStatus;
}
