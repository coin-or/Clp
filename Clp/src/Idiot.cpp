// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "ClpPresolve.hpp"
#include "Idiot.hpp"
#include "CoinTime.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinHelperFunctions.hpp"
// Redefine stuff for Clp
#ifndef OSI_IDIOT
#include "ClpMessage.hpp"
#define OsiObjOffset ClpObjOffset
#endif
/**** strategy 4 - drop, exitDrop and djTolerance all relative:
For first two major iterations these are small.  Then:

drop - exit a major iteration if drop over 5*checkFrequency < this is
used as info->drop*(10.0+fabs(last weighted objective))

exitDrop - exit idiot if feasible and drop < this is
used as info->exitDrop*(10.0+fabs(last objective))

djExit - exit a major iteration if largest dj (averaged over 5 checks)
drops below this - used as info->djTolerance*(10.0+fabs(last weighted objective)

djFlag - mostly skip variables with bad dj worse than this => 2*djExit

djTol - only look at variables with dj better than this => 0.01*djExit
****************/

#define IDIOT_FIX_TOLERANCE 1e-6
#define SMALL_IDIOT_FIX_TOLERANCE 1e-10
int 
Idiot::dropping(IdiotResult result,
		    double tolerance,
		    double small,
                    int *nbad)
{
  if (result.infeas<=small) {
    double value=CoinMax(fabs(result.objval),fabs(result.dropThis))+1.0;
    if (result.dropThis>tolerance*value) {
      *nbad=0;
      return 1;
    } else {
      (*nbad)++;
      if (*nbad>4) {
	return 0;
      } else {
	return 1;
      }
    }
  } else {
    *nbad=0;
    return 1;
  }
}

/* returns -1 or start of costed slacks */
 static int countCostedSlacks(OsiSolverInterface * model)
{
#ifdef OSI_IDIOT
  const CoinPackedMatrix * matrix = model->getMatrixByCol();
#else
  ClpMatrixBase * matrix = model->clpMatrix();
#endif
  const int * row = matrix->getIndices();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const int * columnLength = matrix->getVectorLengths(); 
  const double * element = matrix->getElements();
  const double * rowupper = model->getRowUpper();
  int nrows=model->getNumRows();
  int ncols=model->getNumCols();
  int slackStart = ncols-nrows;
  int nSlacks=nrows;
  int i;
  
  if (ncols<=nrows) return -1;
  while (1) {
    for (i=0;i<nrows;i++) {
      int j=i+slackStart;
      CoinBigIndex k=columnStart[j];
      if (columnLength[j]==1) {
	if (row[k]!=i||element[k]!=1.0) {
	  nSlacks=0;
	  break;
	}
      } else {
	nSlacks=0;
	break;
      }
      if (rowupper[i]<=0.0) {
	nSlacks=0;
	break;
      }
    }
    if (nSlacks||!slackStart) break;
    slackStart=0;
  }
  if (!nSlacks) slackStart=-1;
  return slackStart;
}
void
Idiot::crash(int numberPass, CoinMessageHandler * handler,const CoinMessages *messages)
{
  // lightweight options
  int numberColumns = model_->getNumCols();
  const double * objective = model_->getObjCoefficients();
  int nnzero=0;
  double sum=0.0;
  int i;
  for (i=0;i<numberColumns;i++) {
    if (objective[i]) {
      sum += fabs(objective[i]);
      nnzero++;
    }
  }
  sum /= (double) (nnzero+1);
  maxIts_=2;
  if (numberPass<=0)
    // Cast to double to avoid VACPP complaining
    majorIterations_=(int)(2+log10((double)(numberColumns+1)));
  else
    majorIterations_=numberPass;
  // If mu not changed then compute
  if (mu_==1e-4)
    mu_= CoinMax(1.0e-3,sum*1.0e-5);
  if (!lightWeight_) {
    maxIts2_=105;
  } else if (lightWeight_==1) {
    mu_ *= 1000.0;
    maxIts2_=23;
  } else if (lightWeight_==2) {
    maxIts2_=11;
  } else {
    maxIts2_=23;
  }
  //printf("setting mu to %g and doing %d passes\n",mu_,majorIterations_);
  solve2(handler,messages);
#ifndef OSI_IDIOT
  double averageInfeas = model_->sumPrimalInfeasibilities()/((double) model_->numberRows());
  if (averageInfeas<0.01&&(strategy_&512)!=0) 
    crossOver(16+1); 
  else
    crossOver(3);
#endif
}
void
Idiot::solve()
{
  CoinMessages dummy;
  solve2(NULL,&dummy);
}
void
Idiot::solve2(CoinMessageHandler * handler,const CoinMessages * messages)
{
  int strategy=0;
  double d2;
  int i,n;
  int allOnes=1;
  int iteration=0;
  int iterationTotal=0;
  int nTry=0; /* number of tries at same weight */
  double fixTolerance=IDIOT_FIX_TOLERANCE;
  int maxBigIts=maxBigIts_;
  int maxIts=maxIts_;
  int logLevel=logLevel_;
  int saveMajorIterations = majorIterations_;
  if (handler) {
    if (handler->logLevel()>0&&handler->logLevel()<3)
      logLevel=1;
    else if (!handler->logLevel())
      logLevel=0;
    else
      logLevel=7;
  }
  double djExit=djTolerance_;
  double djFlag=1.0+100.0*djExit;
  double djTol=0.00001;
  double mu =mu_;
  double drop=drop_;
  int maxIts2=maxIts2_;
  double factor=muFactor_;
  double smallInfeas=smallInfeas_;
  double reasonableInfeas=reasonableInfeas_;
  double stopMu=stopMu_;
  double maxmin,offset;
  double lastWeighted=1.0e50;
  double exitDrop=exitDrop_;
  double fakeSmall=smallInfeas;
  double firstInfeas;
  int badIts=0;
  int slackStart,slackEnd,ordStart,ordEnd;
  int checkIteration=0;
  int lambdaIteration=0;
  int belowReasonable=0; /* set if ever gone below reasonable infeas */
  double bestWeighted=1.0e60;
  double bestFeasible=1.0e60; /* best solution while feasible */
  IdiotResult result,lastResult;
  int saveStrategy=strategy_;
  const int strategies[]={0,2,128};
  double saveLambdaScale=0.0;
  if ((saveStrategy&128)!=0) {
    fixTolerance=SMALL_IDIOT_FIX_TOLERANCE;
  }
#ifdef OSI_IDIOT
  const CoinPackedMatrix * matrix = model_->getMatrixByCol();
#else
  ClpMatrixBase * matrix = model_->clpMatrix();
#endif
  const int * row = matrix->getIndices();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const int * columnLength = matrix->getVectorLengths(); 
  const double * element = matrix->getElements();
  int nrows=model_->getNumRows();
  int ncols=model_->getNumCols();
  double * rowsol, * colsol;
  double * pi, * dj;
#ifndef OSI_IDIOT
  double * cost = model_->objective();
  double * lower = model_->columnLower();
  double * upper = model_->columnUpper();
#else
  double * cost = new double [ncols];
  memcpy( cost, model_->getObjCoefficients(), ncols*sizeof(double));
  const double * lower = model_->getColLower();
  const double * upper = model_->getColUpper();
#endif
  const double *elemXX;
  double * saveSol;
  double * rowupper= new double[nrows]; // not const as modified
  memcpy(rowupper,model_->getRowUpper(),nrows*sizeof(double));
  double * rowlower= new double[nrows]; // not const as modified
  memcpy(rowlower,model_->getRowLower(),nrows*sizeof(double));
  int * whenUsed;
  double * lambda;
  saveSol=new double[ncols];
  lambda=new double [nrows];
  rowsol= new double[nrows];
  colsol= new double [ncols];
  memcpy(colsol,model_->getColSolution(),ncols*sizeof(double));
  pi= new double[nrows];
  dj=new double[ncols];
  delete [] whenUsed_;
  whenUsed=whenUsed_=new int[ncols];
  if (model_->getObjSense()==-1.0) {
    maxmin=-1.0;
  } else {
    maxmin=1.0;
  }
  model_->getDblParam(OsiObjOffset,offset);
  if (!maxIts2) maxIts2=maxIts;
  strategy=strategy_;
  strategy &= 3;
  memset(lambda,0,nrows*sizeof(double));
  slackStart=countCostedSlacks(model_);
  if (slackStart>=0) {
    printf("This model has costed slacks\n");
    slackEnd=slackStart+nrows;
    if (slackStart) {
      ordStart=0;
      ordEnd=slackStart;
    } else {
      ordStart=nrows;
      ordEnd=ncols;
    }
  } else {
    slackEnd=slackStart;
    ordStart=0;
    ordEnd=ncols;
  }
  if (offset&&logLevel>2) {
    printf("** Objective offset is %g\n",offset);
  }
  /* compute reasonable solution cost */
  for (i=0;i<nrows;i++) {
    rowsol[i]=1.0e31;
  }
  for (i=0;i<ncols;i++) {
    CoinBigIndex j;
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      if (element[j]!=1.0) {
	allOnes=0;
	break;
      }
    }
  }
  if (allOnes) {
    elemXX=NULL;
  } else {
    elemXX=element;
  }
  // Do scaling if wanted
  bool scaled=false;
#ifndef OSI_IDIOT
  if ((strategy_&32)!=0&&!allOnes) {
    if (model_->scalingFlag()>0)
      scaled = model_->clpMatrix()->scale(model_)==0;
    if (scaled) {
      const double * rowScale = model_->rowScale();
      const double * columnScale = model_->columnScale();
      double * oldLower = lower;
      double * oldUpper = upper;
      double * oldCost = cost;
      lower = new double[ncols];
      upper = new double[ncols];
      cost = new double[ncols];
      memcpy(lower,oldLower,ncols*sizeof(double));
      memcpy(upper,oldUpper,ncols*sizeof(double));
      memcpy(cost,oldCost,ncols*sizeof(double));
      int icol,irow;
      for (icol=0;icol<ncols;icol++) {
	double multiplier = 1.0/columnScale[icol];
	if (lower[icol]>-1.0e50)
	  lower[icol] *= multiplier;
	if (upper[icol]<1.0e50)
	  upper[icol] *= multiplier;
	colsol[icol] *= multiplier;
	cost[icol] *= columnScale[icol];
      }
      memcpy(rowlower,model_->rowLower(),nrows*sizeof(double));
      for (irow=0;irow<nrows;irow++) {
	double multiplier = rowScale[irow];
	if (rowlower[irow]>-1.0e50)
	  rowlower[irow] *= multiplier;
	if (rowupper[irow]<1.0e50)
	  rowupper[irow] *= multiplier;
	rowsol[irow] *= multiplier;
      }
      int length = columnStart[ncols-1]+columnLength[ncols-1];
      double * elemYY = new double[length];
      for (i=0;i<ncols;i++) {
	CoinBigIndex j;
	double scale = columnScale[i];
	for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	  int irow=row[j];
	  elemYY[j] = element[j]*scale*rowScale[irow];
	}
      }
      elemXX=elemYY;
    }
  }
#endif
  for (i=0;i<ncols;i++) {
    CoinBigIndex j;
    double dd=columnLength[i];
    dd=cost[i]/dd;
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      int irow=row[j];
      if (dd<rowsol[irow]) {
	rowsol[irow]=dd;
      }
    }
  }
  d2=0.0;
  for (i=0;i<nrows;i++) {
    d2+=rowsol[i];
  }
  d2*=2.0; /* for luck */
  
  d2=d2/((double) (4*nrows+8000));
  d2*=0.5; /* halve with more flexible method */
  if (d2<5.0) d2=5.0;
  if (djExit==0.0) {
    djExit=d2;
  }
  if ((saveStrategy&4)!=0) {
    /* go to relative tolerances - first small */
    djExit=1.0e-10;
    djFlag=1.0e-5;
    drop=1.0e-10;
  }
  memset(whenUsed,0,ncols*sizeof(int));
  strategy=strategies[strategy];
  if ((saveStrategy&8)!=0) strategy |= 64; /* don't allow large theta */
  memcpy(saveSol,colsol,ncols*sizeof(double));
  
  lastResult=IdiSolve(nrows,ncols,rowsol ,colsol,pi,
		       dj,cost,rowlower,rowupper,
		       lower,upper,elemXX,row,columnStart,columnLength,lambda,
		       0,mu,drop,
		       maxmin,offset,strategy,djTol,djExit,djFlag);
  n=0;
  for (i=ordStart;i<ordEnd;i++) {
    if (colsol[i]>lower[i]+fixTolerance) {
      if (colsol[i]<upper[i]-fixTolerance) {
	n++;
      } else {
	colsol[i]=upper[i];
      }
      whenUsed[i]=iteration;
    } else {
      colsol[i]=lower[i];
    }
  }
  if ((logLevel_&1)!=0) {
#ifndef OSI_IDIOT
    if (!handler) {
#endif
      printf("Iteration %d infeasibility %g, objective %g - mu %g, its %d, %d interior\n", 
	     iteration,lastResult.infeas,lastResult.objval,mu,lastResult.iteration,n);
#ifndef OSI_IDIOT
    } else {
      handler->message(CLP_IDIOT_ITERATION,*messages)
	<<iteration<<lastResult.infeas<<lastResult.objval<<mu<<lastResult.iteration<<n
	<<CoinMessageEol;
    }
#endif
  }
  int numberBaseTrys=0; // for first time
  int numberAway=-1;
  iterationTotal = lastResult.iteration;
  firstInfeas=lastResult.infeas;
  if ((strategy_&1024)!=0) reasonableInfeas=0.5*firstInfeas;
  if (lastResult.infeas<reasonableInfeas) lastResult.infeas=reasonableInfeas;
  double keepinfeas=1.0e31;
  double lastInfeas=1.0e31;
  double bestInfeas=1.0e31;
  while ((mu>stopMu&&lastResult.infeas>smallInfeas)||
         (lastResult.infeas<=smallInfeas&&
	 dropping(lastResult,exitDrop,smallInfeas,&badIts))||
	 checkIteration<2||lambdaIteration<lambdaIterations_) {
    if (lastResult.infeas<=exitFeasibility_)
      break; 
    iteration++;
    checkIteration++;
    if (lastResult.infeas<=smallInfeas&&lastResult.objval<bestFeasible) {
      bestFeasible=lastResult.objval;
    }
    if ((saveStrategy&4096)) strategy |=256;
    if ((saveStrategy&4)!=0&&iteration>2) {
      /* go to relative tolerances */
      double weighted=10.0+fabs(lastWeighted);
      djExit=djTolerance_*weighted;
      djFlag=2.0*djExit;
      drop=drop_*weighted;
      djTol=0.01*djExit;
    }
    result=IdiSolve(nrows,ncols,rowsol ,colsol,pi,dj,
		     cost,rowlower,rowupper,
		     lower,upper,elemXX,row,columnStart,columnLength,lambda,
		     maxIts,mu,drop,
		     maxmin,offset,strategy,djTol,djExit,djFlag);
    n=0;
    for (i=ordStart;i<ordEnd;i++) {
      if (colsol[i]>lower[i]+fixTolerance) {
	if (colsol[i]<upper[i]-fixTolerance) {
	  n++;
	} else {
	  colsol[i]=upper[i];
	}
	whenUsed[i]=iteration;
      } else {
	colsol[i]=lower[i];
      }
    }
    if ((logLevel_&1)!=0) {
#ifndef OSI_IDIOT
      if (!handler) {
#endif
	printf("Iteration %d infeasibility %g, objective %g - mu %g, its %d, %d interior\n", 
	       iteration,result.infeas,result.objval,mu,result.iteration,n);
#ifndef OSI_IDIOT
      } else {
	handler->message(CLP_IDIOT_ITERATION,*messages)
	  <<iteration<<result.infeas<<result.objval<<mu<<result.iteration<<n
	  <<CoinMessageEol;
      }
#endif
    }
    if (iteration>50&&n==numberAway&&result.infeas<1.0e-4)
      break; // not much happening
    if (lightWeight_==1&&iteration>10&&result.infeas>1.0&&maxIts!=7) {
      if (lastInfeas!=bestInfeas&&CoinMin(result.infeas,lastInfeas)>0.95*bestInfeas)
	majorIterations_ = CoinMin(majorIterations_,iteration); // not getting feasible
    }
    lastInfeas = result.infeas;
    numberAway=n;
    keepinfeas = result.infeas;
    lastWeighted=result.weighted;
    iterationTotal += result.iteration;
    if (iteration==1) {
      if ((strategy_&1024)!=0&&mu<1.0e-10) 
	result.infeas=firstInfeas*0.8;
      if (majorIterations_>=50)
        result.infeas *= 0.8;
      if (result.infeas>firstInfeas*0.9
	  &&result.infeas>reasonableInfeas) {
	iteration--;
	if (majorIterations_<50)
	  mu*=1.0e-1;
	else
	  mu*=0.7;
	bestFeasible=1.0e31;
        bestWeighted=1.0e60;
        numberBaseTrys++;
        if (mu<1.0e-30||(numberBaseTrys>10&&lightWeight_)) {
          // back to all slack basis
          lightWeight_=2;
          break;
        }
	memcpy(colsol,saveSol,ncols*sizeof(double));
      } else {
        // Save best solution
        memcpy(saveSol,colsol,ncols*sizeof(double));
	maxIts=maxIts2;
	checkIteration=0;
	if ((strategy_&1024)!=0) mu *= 1.0e-1;
      }
    } else if (result.infeas<bestInfeas) {
      // Save best solution
      memcpy(saveSol,colsol,ncols*sizeof(double));
    }
    bestInfeas=CoinMin(bestInfeas,result.infeas);
    if (iteration) {
      /* this code is in to force it to terminate sometime */
      double changeMu=factor;
      if ((saveStrategy&64)!=0) {
	keepinfeas=0.0; /* switch off ranga's increase */
	fakeSmall=smallInfeas;
      } else {
	fakeSmall=-1.0;
      }
      saveLambdaScale=0.0;
      if (result.infeas>reasonableInfeas||
	  (nTry+1==maxBigIts&&result.infeas>fakeSmall)) {
	if (result.infeas>lastResult.infeas*(1.0-dropEnoughFeasibility_)||
	    nTry+1==maxBigIts||
	    (result.infeas>lastResult.infeas*0.9
	     &&result.weighted>lastResult.weighted
	     -dropEnoughWeighted_*fabs(lastResult.weighted))) {
	  mu*=changeMu;
          if ((saveStrategy&32)!=0&&result.infeas<reasonableInfeas&&0) {
	    reasonableInfeas=CoinMax(smallInfeas,reasonableInfeas*sqrt(changeMu));
	    printf("reasonable infeas now %g\n",reasonableInfeas);
	  }
	  nTry=0;
	  bestFeasible=1.0e31;
	  bestWeighted=1.0e60;
	  checkIteration=0;
	  lambdaIteration=0;
#define LAMBDA
#ifdef LAMBDA
	  if ((saveStrategy&2048)==0) {
	    memset(lambda,0,nrows*sizeof(double));
	  }
#else
	  memset(lambda,0,nrows*sizeof(double));
#endif
	} else {
	  nTry++;
	}
      } else if (lambdaIterations_>=0) {
	/* update lambda  */
	double scale=1.0/mu;
	int i,nnz=0;
	saveLambdaScale=scale;
         lambdaIteration++;
         if ((saveStrategy&4)==0) drop = drop_/50.0;
         if (lambdaIteration>4 && 
            (((lambdaIteration%10)==0 && smallInfeas<keepinfeas) ||
             (lambdaIteration%5)==0 && 1.5*smallInfeas<keepinfeas)) {
           //printf(" Increasing smallInfeas from %f to %f\n",smallInfeas,1.5*smallInfeas);
           smallInfeas *= 1.5;
         }
         if ((saveStrategy&2048)==0) {
	   for (i=0;i<nrows;i++) {
	     if (lambda[i]) nnz++;
	     lambda[i]+= scale*rowsol[i];
	   }
	 } else {
	   nnz=1;
#ifdef LAMBDA
	   for (i=0;i<nrows;i++) {
	     lambda[i]+= scale*rowsol[i];
	   }
#else
	   for (i=0;i<nrows;i++) {
	     lambda[i] = scale*rowsol[i];
	   }
	   for (i=0;i<ncols;i++) {
	     CoinBigIndex j;
	     double value=cost[i]*maxmin;
	     for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	       int irow=row[j];
	       value+=element[j]*lambda[irow];
	     }
	     cost[i]=value*maxmin;
	   }
	   for (i=0;i<nrows;i++) {
	     offset+=lambda[i]*rowupper[i];
	     lambda[i]=0.0;
	   }
#ifdef DEBUG
           printf("offset %g\n",offset);
#endif
	   model_->setDblParam(OsiObjOffset,offset);
#endif
	 }
        nTry++;
	if (!nnz) {
	  bestFeasible=1.0e32;
	  bestWeighted=1.0e60;
	  checkIteration=0;
	  result.weighted=1.0e31;
	}
#ifdef DEBUG
        double trueCost=0.0;
	for (i=0;i<ncols;i++) {
	  int j;
	  trueCost+=cost[i]*colsol[i];
	}
	printf("True objective %g\n",trueCost-offset);
#endif
      } else {
        nTry++;
      }
      lastResult=result;
      if (result.infeas<reasonableInfeas&&!belowReasonable) {
	belowReasonable=1;
	bestFeasible=1.0e32;
        bestWeighted=1.0e60;
	checkIteration=0;
	result.weighted=1.0e31;
      }
    }
    if (iteration>=majorIterations_) {
      // If small and not feasible and crash then dive dive dive
      if (0&&result.infeas>1.0&&majorIterations_<30&&(maxIts2_==11||maxIts2_==23)) {
	maxIts=7;
	majorIterations_=100;
      } else {
	if (logLevel>2) 
	  printf("Exiting due to number of major iterations\n");
	break;
      }
    }
  }
  majorIterations_ = saveMajorIterations;
  // put back best solution
  memcpy(colsol,saveSol,ncols*sizeof(double));
#ifndef OSI_IDIOT
  if (scaled) {
    // Scale solution and free arrays
    const double * rowScale = model_->rowScale();
    const double * columnScale = model_->columnScale();
    int icol,irow;
    for (icol=0;icol<ncols;icol++) {
      colsol[icol] *= columnScale[icol];
      dj[icol] /= columnScale[icol];
    }
    for (irow=0;irow<nrows;irow++) {
      rowsol[irow] /= rowScale[irow];
      pi[irow] *= rowScale[irow];
    }
    // Don't know why getting Microsoft problems
#if defined (_MSC_VER)
    delete [] ( double *) elemXX;
#else
    delete [] elemXX;
#endif
    model_->setRowScale(NULL);
    model_->setColumnScale(NULL);
    delete [] lower;
    delete [] upper;
    delete [] cost;
    lower = model_->columnLower();
    upper = model_->columnUpper();
    cost = model_->objective();
    //rowlower = model_->rowLower();
  }
#endif
#define TRYTHIS
#ifdef TRYTHIS
  if ((saveStrategy&2048)!=0) {
    double offset;
    model_->getDblParam(OsiObjOffset,offset);
    for (i=0;i<ncols;i++) {
      CoinBigIndex j;
      double djval=cost[i]*maxmin;
      for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	int irow=row[j];
	djval -= element[j]*lambda[irow];
      }
      cost[i]=djval;
    }
    for (i=0;i<nrows;i++) {
      offset+=lambda[i]*rowupper[i];
    }
    model_->setDblParam(OsiObjOffset,offset);
  }
#endif
  if (saveLambdaScale) {
    /* back off last update */
    for (i=0;i<nrows;i++) {
      lambda[i]-= saveLambdaScale*rowsol[i];
    }
  }
  muAtExit_=mu;
  n=0;
  for (i=ordStart;i<ordEnd;i++) {
    if (colsol[i]>lower[i]+fixTolerance) {
      n++;
      whenUsed[i]=iteration;
    } else {
      colsol[i]=lower[i];
    }
  }
  if ((logLevel&1)==0) {
    printf(
	    "%d - mu %g, infeasibility %g, objective %g, %d interior\n",
	    iteration,mu,lastResult.infeas,lastResult.objval,n);
  }
#ifndef OSI_IDIOT
  model_->setSumPrimalInfeasibilities(lastResult.infeas);
#endif
  {
    double large=0.0;
    int i;
    memset(rowsol,0,nrows*sizeof(double));
    for (i=0;i<ncols;i++) {
      CoinBigIndex j;
      double value=colsol[i];
      for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	int irow=row[j];
	rowsol[irow] += element[j]*value;
      }
    }
    for (i=0;i<nrows;i++) {
      if (rowsol[i] > rowupper[i]) {
	double diff=rowsol[i]-rowupper[i];
	if (diff>large) 
	  large=diff;
      } else if (rowsol[i] < rowlower[i]) {
	double diff=rowlower[i]-rowsol[i];
	if (diff>large) 
	  large=diff;
      } 
    }
    if (logLevel>2)
      printf("largest infeasibility is %g\n",large);
  }
  /* subtract out lambda */
  for (i=0;i<nrows;i++) {
    pi[i]-=lambda[i];
  }
  for (i=0;i<ncols;i++) {
    CoinBigIndex j;
    double djval=cost[i]*maxmin;
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      int irow=row[j];
      djval -= element[j]*pi[irow];
    }
    dj[i]=djval;
  }
  if ((strategy_&1024)!=0) {
    double ratio = ((double) ncols)/((double) nrows);
    printf("col/row ratio %g infeas ratio %g\n",ratio,lastResult.infeas/firstInfeas);
    if (lastResult.infeas>0.01*firstInfeas*ratio) {
      strategy_ &= (~1024);
      printf(" - layer off\n");
    } else {
      printf(" - layer on\n");
    }
  }
  delete [] saveSol;
  delete [] lambda;
  // save solution 
  // duals not much use - but save anyway
#ifndef OSI_IDIOT
  memcpy(model_->primalRowSolution(),rowsol,nrows*sizeof(double));
  memcpy(model_->primalColumnSolution(),colsol,ncols*sizeof(double));
  memcpy(model_->dualRowSolution(),pi,nrows*sizeof(double));
  memcpy(model_->dualColumnSolution(),dj,ncols*sizeof(double));
#else
  model_->setColSolution(colsol);
  model_->setRowPrice(pi);
  delete [] cost;
#endif
  delete [] rowsol;
  delete [] colsol;
  delete [] pi;
  delete [] dj;
  delete [] rowlower;
  delete [] rowupper;
  return ;
}
#ifndef OSI_IDIOT
void
Idiot::crossOver(int mode)
{
  if (lightWeight_==2) {
    // total failure
    model_->allSlackBasis();
    return;
  }
  double fixTolerance=IDIOT_FIX_TOLERANCE;
  double startTime = CoinCpuTime();
  ClpSimplex * saveModel=NULL;
  ClpMatrixBase * matrix = model_->clpMatrix();
  const int * row = matrix->getIndices();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const int * columnLength = matrix->getVectorLengths(); 
  const double * element = matrix->getElements();
  const double * rowupper = model_->getRowUpper();
  int nrows=model_->getNumRows();
  int ncols=model_->getNumCols();
  double * rowsol, * colsol;
  // different for Osi
  double * lower = model_->columnLower();
  double * upper = model_->columnUpper();
  const double * rowlower= model_->getRowLower();
  int * whenUsed=whenUsed_;
  rowsol= model_->primalRowSolution();
  colsol= model_->primalColumnSolution();;
  double * cost=model_->objective();

  int slackEnd,ordStart,ordEnd;
  int slackStart = countCostedSlacks(model_);

  int addAll = mode&7;
  int presolve=0;

  double djTolerance = djTolerance_;
  if (djTolerance>0.0&&djTolerance<1.0)
    djTolerance=1.0;
  int iteration;
  int i, n=0;
  double ratio=1.0;
  double objValue=0.0;
  if ((strategy_&128)!=0) {
    fixTolerance=SMALL_IDIOT_FIX_TOLERANCE;
  }
  if ((mode&16)!=0&&addAll<3) presolve=1;
  double * saveUpper = NULL;
  double * saveLower = NULL;
  if (addAll<3) {
    saveUpper = new double [ncols];
    saveLower = new double [ncols];
    memcpy(saveUpper,upper,ncols*sizeof(double));
    memcpy(saveLower,lower,ncols*sizeof(double));
  }
  if (slackStart>=0) {
    slackEnd=slackStart+nrows;
    if (slackStart) {
      ordStart=0;
      ordEnd=slackStart;
    } else {
      ordStart=nrows;
      ordEnd=ncols;
    }
  } else {
    slackEnd=slackStart;
    ordStart=0;
    ordEnd=ncols;
  }
  /* get correct rowsol (without known slacks) */
  memset(rowsol,0,nrows*sizeof(double));
  for (i=ordStart;i<ordEnd;i++) {
    CoinBigIndex j;
    double value=colsol[i];
    if (value<lower[i]+fixTolerance) {
      value=lower[i];
      colsol[i]=value;
    }
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      int irow=row[j];
      rowsol[irow]+=value*element[j];
    }
  }
  if (slackStart>=0) {
    for (i=0;i<nrows;i++) {
      if (ratio*rowsol[i]>rowlower[i]&&rowsol[i]>1.0e-8) {
	ratio=rowlower[i]/rowsol[i];
      }
    }
    for (i=0;i<nrows;i++) {
      rowsol[i]*=ratio;
    }
    for (i=ordStart;i<ordEnd;i++) {
      double value=colsol[i]*ratio;
      colsol[i]=value;
      objValue+=value*cost[i];
    }
    for (i=0;i<nrows;i++) {
      double value=rowlower[i]-rowsol[i];
      colsol[i+slackStart]=value;
      objValue+=value*cost[i+slackStart];
    }
    printf("New objective after scaling %g\n",objValue);
  }
#if 0
   maybe put back - but just get feasible ?
  // If not many fixed then just exit
  int numberFixed=0;
  for (i=ordStart;i<ordEnd;i++) {
    if (colsol[i]<lower[i]+fixTolerance) 
      numberFixed++;
    else if (colsol[i]>upper[i]-fixTolerance) 
      numberFixed++;
  }
  if (numberFixed<ncols/2) {
    addAll=3;
    presolve=0;
  }
#endif
  model_->createStatus();
  /* addAll
     0 - chosen,all used, all 
     1 - chosen, all
     2 - all
     3 - do not do anything  - maybe basis
  */
  for (i=ordStart;i<ordEnd;i++) {
    if (addAll<2) {
      if (colsol[i]<lower[i]+fixTolerance) {
	upper[i]=lower[i];
	colsol[i]=lower[i];
      } else if (colsol[i]>upper[i]-fixTolerance) {
	lower[i]=upper[i];
	colsol[i]=upper[i];
      }
    }
    model_->setColumnStatus(i,ClpSimplex::superBasic);
  }
  double maxmin;
  if (model_->getObjSense()==-1.0) {
    maxmin=-1.0;
  } else {
    maxmin=1.0;
  }
  if (slackStart>=0) {
    for (i=0;i<nrows;i++) {
      model_->setRowStatus(i,ClpSimplex::superBasic);
    }
    for (i=slackStart;i<slackEnd;i++) {
      model_->setColumnStatus(i,ClpSimplex::basic);
    }
  } else {
    /* still try and put singletons rather than artificials in basis */
    int ninbas=0;
    for (i=0;i<nrows;i++) {
      model_->setRowStatus(i,ClpSimplex::basic);
    }
    for (i=0;i<ncols;i++) {
      if (columnLength[i]==1&&upper[i]>lower[i]+1.0e-5) {
	CoinBigIndex j =columnStart[i];
	double value=element[j];
	int irow=row[j];
	double rlo=rowlower[irow];
	double rup=rowupper[irow];
	double clo=lower[i];
	double cup=upper[i];
	double csol=colsol[i];
	/* adjust towards feasibility */
	double move=0.0;
	if (rowsol[irow]>rup) {
	  move=(rup-rowsol[irow])/value;
	  if (value>0.0) {
	    /* reduce */
	    if (csol+move<clo) move=CoinMin(0.0,clo-csol);
	  } else {
	    /* increase */
	    if (csol+move>cup) move=CoinMax(0.0,cup-csol);
	  }
	} else if (rowsol[irow]<rlo) {
	  move=(rlo-rowsol[irow])/value;
	  if (value>0.0) {
	    /* increase */
	    if (csol+move>cup) move=CoinMax(0.0,cup-csol);
	  } else {
	    /* reduce */
	    if (csol+move<clo) move=CoinMin(0.0,clo-csol);
	  }
	} else {
	  /* move to improve objective */
	  if (cost[i]*maxmin>0.0) {
	    if (value>0.0) {
	      move=(rlo-rowsol[irow])/value;
	      /* reduce */
	      if (csol+move<clo) move=CoinMin(0.0,clo-csol);
	    } else {
	      move=(rup-rowsol[irow])/value;
	      /* increase */
	      if (csol+move>cup) move=CoinMax(0.0,cup-csol);
	    }
	  } else if (cost[i]*maxmin<0.0) {
	    if (value>0.0) {
	      move=(rup-rowsol[irow])/value;
	      /* increase */
	      if (csol+move>cup) move=CoinMax(0.0,cup-csol);
	    } else {
	      move=(rlo-rowsol[irow])/value;
	      /* reduce */
	      if (csol+move<clo) move=CoinMin(0.0,clo-csol);
	    }
	  }
	}
	rowsol[irow] +=move*value;
	colsol[i]+=move;
	/* put in basis if row was artificial */
	if (rup-rlo<1.0e-7&&model_->getRowStatus(irow)==ClpSimplex::basic) {
	  model_->setRowStatus(irow,ClpSimplex::superBasic);
	  model_->setColumnStatus(i,ClpSimplex::basic);
	  ninbas++;
	}
      }
    }
    /*printf("%d in basis\n",ninbas);*/
  }
  bool wantVector=false;
  if (dynamic_cast< ClpPackedMatrix*>(model_->clpMatrix())) {
    // See if original wanted vector
    ClpPackedMatrix * clpMatrixO = dynamic_cast< ClpPackedMatrix*>(model_->clpMatrix());
    wantVector = clpMatrixO->wantsSpecialColumnCopy();
  }
  if (addAll<3) {
    ClpPresolve pinfo;
    if (presolve) {
      saveModel = model_;
      model_ = pinfo.presolvedModel(*model_,1.0e-8,false,5);
    }
    if (model_) {
      if (!wantVector) {
	model_->primal(1);
      } else {
	ClpMatrixBase * matrix = model_->clpMatrix();
	ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
	assert (clpMatrix);
	clpMatrix->makeSpecialColumnCopy();
	model_->primal(1);
	clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
	pinfo.postsolve(true);
	delete model_;
	model_ = saveModel;
	saveModel=NULL;
      }
    } else {
      // not feasible
      addAll=1;
      presolve=0;
      model_ = saveModel;
      saveModel=NULL;
    }
    if (addAll<2) {
      n=0;
      if (!addAll ) {
	/* could do scans to get a good number */
	iteration=1;
	for (i=ordStart;i<ordEnd;i++) {
	  if (whenUsed[i]>=iteration) {
	    if (upper[i]-lower[i]<1.0e-5&&saveUpper[i]-saveLower[i]>1.0e-5) {
	      n++;
	      upper[i]=saveUpper[i];
	      lower[i]=saveLower[i];
	    }
	  }
	}
      } else {
	for (i=ordStart;i<ordEnd;i++) {
	  if (upper[i]-lower[i]<1.0e-5&&saveUpper[i]-saveLower[i]>1.0e-5) {
	    n++;
	    upper[i]=saveUpper[i];
	    lower[i]=saveLower[i];
	  }
	}
      }
      printf("Time so far %g, %d now added from previous iterations\n",
	     CoinCpuTime()-startTime,n);
      if (addAll)
	presolve=0;
      if (presolve) {
	saveModel = model_;
	model_ = pinfo.presolvedModel(*model_,1.0e-8,false,5);
      } else {
	presolve=0;
      }
      if (!wantVector) {
	model_->primal(1);
      } else {
	ClpMatrixBase * matrix = model_->clpMatrix();
	ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
	assert (clpMatrix);
	clpMatrix->makeSpecialColumnCopy();
	model_->primal(1);
	clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
	pinfo.postsolve(true);
	delete model_;
	model_ = saveModel;
	saveModel=NULL;
      }
      if (!addAll) {
	n=0;
	for (i=ordStart;i<ordEnd;i++) {
	  if (upper[i]-lower[i]<1.0e-5&&saveUpper[i]-saveLower[i]>1.0e-5) {
	    n++;
	    upper[i]=saveUpper[i];
	    lower[i]=saveLower[i];
	  }
	}
	printf("Time so far %g, %d now added from previous iterations\n",
	       CoinCpuTime()-startTime,n);
      }
      if (presolve) {
	saveModel = model_;
	model_ = pinfo.presolvedModel(*model_,1.0e-8,false,5);
      } else {
	presolve=0;
      }
      if (!wantVector) {
	model_->primal(1);
      } else {
	ClpMatrixBase * matrix = model_->clpMatrix();
	ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
	assert (clpMatrix);
	clpMatrix->makeSpecialColumnCopy();
	model_->primal(1);
	clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
	pinfo.postsolve(true);
	delete model_;
	model_ = saveModel;
	saveModel=NULL;
      }
    }
    printf("Total time in crossover %g\n", CoinCpuTime()-startTime);
    delete [] saveUpper;
    delete [] saveLower;
  }
  return ;
}
#endif
/*****************************************************************************/

// Default contructor
Idiot::Idiot()
{
  model_ = NULL;
  maxBigIts_ = 3;
  maxIts_ = 5;
  logLevel_ = 1; 
  logFreq_ = 100;
  maxIts2_ = 100;
  djTolerance_ = 1e-1;
  mu_ = 1e-4;
  drop_ = 5.0;
  exitDrop_=-1.0e20;
  muFactor_ = 0.3333;
  stopMu_ = 1e-12;
  smallInfeas_ = 1e-1;
  reasonableInfeas_ = 1e2;
  muAtExit_ =1.0e31;
  strategy_ =8;
  lambdaIterations_ =0;
  checkFrequency_ =100;
  whenUsed_ = NULL;
  majorIterations_ =30;
  exitFeasibility_ =-1.0;
  dropEnoughFeasibility_ =0.02;
  dropEnoughWeighted_ =0.01;
  // adjust
  double nrows=10000.0;
  int baseIts =(int) sqrt((double)nrows);
  baseIts =baseIts/10;
  baseIts *= 10;
  maxIts2_ =200+baseIts+5;
  reasonableInfeas_ =((double) nrows)*0.05;
  lightWeight_=0;
}
// Constructor from model
Idiot::Idiot(OsiSolverInterface &model)
{
  model_ = & model;
  maxBigIts_ = 3;
  maxIts_ = 5;
  logLevel_ = 1; 
  logFreq_ = 100;
  maxIts2_ = 100;
  djTolerance_ = 1e-1;
  mu_ = 1e-4;
  drop_ = 5.0;
  exitDrop_=-1.0e20;
  muFactor_ = 0.3333;
  stopMu_ = 1e-12;
  smallInfeas_ = 1e-1;
  reasonableInfeas_ = 1e2;
  muAtExit_ =1.0e31;
  strategy_ =8;
  lambdaIterations_ =0;
  checkFrequency_ =100;
  whenUsed_ = NULL;
  majorIterations_ =30;
  exitFeasibility_ =-1.0;
  dropEnoughFeasibility_ =0.02;
  dropEnoughWeighted_ =0.01;
  // adjust
  double nrows;
  if (model_)
    nrows=model_->getNumRows();
  else
    nrows=10000.0;
  int baseIts =(int) sqrt((double)nrows);
  baseIts =baseIts/10;
  baseIts *= 10;
  maxIts2_ =200+baseIts+5;
  reasonableInfeas_ =((double) nrows)*0.05;
  lightWeight_=0;
}
// Copy constructor. 
Idiot::Idiot(const Idiot &rhs)
{
  model_ = rhs.model_;
  if (model_&&rhs.whenUsed_) {
    int numberColumns = model_->getNumCols();
    whenUsed_ = new int [numberColumns];
    memcpy(whenUsed_,rhs.whenUsed_,numberColumns*sizeof(int));
  } else {
    whenUsed_=NULL;
  }
  djTolerance_ = rhs.djTolerance_;
  mu_ = rhs.mu_;
  drop_ = rhs.drop_;
  muFactor_ = rhs.muFactor_;
  stopMu_ = rhs.stopMu_;
  smallInfeas_ = rhs.smallInfeas_;
  reasonableInfeas_ = rhs.reasonableInfeas_;
  exitDrop_ = rhs.exitDrop_;
  muAtExit_ = rhs.muAtExit_;
  exitFeasibility_ = rhs.exitFeasibility_;
  dropEnoughFeasibility_ = rhs.dropEnoughFeasibility_;
  dropEnoughWeighted_ = rhs.dropEnoughWeighted_;
  maxBigIts_ = rhs.maxBigIts_;
  maxIts_ = rhs.maxIts_;
  majorIterations_ = rhs.majorIterations_;
  logLevel_ = rhs.logLevel_;
  logFreq_ = rhs.logFreq_;
  checkFrequency_ = rhs.checkFrequency_;
  lambdaIterations_ = rhs.lambdaIterations_;
  maxIts2_ = rhs.maxIts2_;
  strategy_ = rhs.strategy_;
  lightWeight_=rhs.lightWeight_;
}
// Assignment operator. This copies the data
Idiot & 
Idiot::operator=(const Idiot & rhs)
{
  if (this != &rhs) {
    delete [] whenUsed_;
    model_ = rhs.model_;
    if (model_&&rhs.whenUsed_) {
      int numberColumns = model_->getNumCols();
      whenUsed_ = new int [numberColumns];
      memcpy(whenUsed_,rhs.whenUsed_,numberColumns*sizeof(int));
    } else {
      whenUsed_=NULL;
    }
    djTolerance_ = rhs.djTolerance_;
    mu_ = rhs.mu_;
    drop_ = rhs.drop_;
    muFactor_ = rhs.muFactor_;
    stopMu_ = rhs.stopMu_;
    smallInfeas_ = rhs.smallInfeas_;
    reasonableInfeas_ = rhs.reasonableInfeas_;
    exitDrop_ = rhs.exitDrop_;
    muAtExit_ = rhs.muAtExit_;
    exitFeasibility_ = rhs.exitFeasibility_;
    dropEnoughFeasibility_ = rhs.dropEnoughFeasibility_;
    dropEnoughWeighted_ = rhs.dropEnoughWeighted_;
    maxBigIts_ = rhs.maxBigIts_;
    maxIts_ = rhs.maxIts_;
    majorIterations_ = rhs.majorIterations_;
    logLevel_ = rhs.logLevel_;
    logFreq_ = rhs.logFreq_;
    checkFrequency_ = rhs.checkFrequency_;
    lambdaIterations_ = rhs.lambdaIterations_;
    maxIts2_ = rhs.maxIts2_;
    strategy_ = rhs.strategy_;
    lightWeight_=rhs.lightWeight_;
  }
  return *this;
}
Idiot::~Idiot()
{
  delete [] whenUsed_;
}
