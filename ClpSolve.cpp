// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

// This file has higher level solve functions


#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpSolve.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpMessage.hpp"
#include "CoinTime.hpp"

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
ClpSimplex::initialSolve(ClpSolve & options)
{
  ClpSolve::SolveType method=options.getSolveType();
  ClpSolve::PresolveType presolve = options.getPresolveType();
  int saveMaxIterations = maximumIterations();
  int finalStatus=-1;
  int numberIterations=0;
  double time1 = CoinCpuTime();
  double timeX = time1;
  double time2;
  ClpMatrixBase * saveMatrix=NULL;
  ClpSimplex * model2 = this;
  bool interrupt = (options.getSpecialOption(2)==0);
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
  int savePerturbation=perturbation_;
  int saveScaling = scalingFlag_;
  if (dynamic_cast< ClpNetworkMatrix*>(matrix_)) {
    // network - switch off stuff
    presolve = ClpSolve::presolveOff;
  }
  // For below >0 overrides
  // 0 means no, -1 means maybe
  int doIdiot=0;
  int doCrash=0;
  int doSprint=0;
  if (method!=ClpSolve::useDual) {
    switch (options.getSpecialOption(1)) {
    case 0:
      doIdiot=-1;
      doCrash=-1;
      doSprint=-1;
      break;
    case 1:
      doIdiot=0;
      doCrash=1;
      doSprint=0;
      break;
    case 2:
      doIdiot=1;
      doCrash=0;
      doSprint=0;
      break;
    case 3:
      doIdiot=0;
      doCrash=0;
      doSprint=1;
      break;
    case 4:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      break;
    case 5:
      doIdiot=0;
      doCrash=-1;
      doSprint=-1;
      break;
    case 6:
      doIdiot=-1;
      doCrash=-1;
      doSprint=0;
      break;
    case 7:
      doIdiot=-1;
      doCrash=0;
      doSprint=-1;
      break;
    case 8:
      doIdiot=-1;
      doCrash=0;
      doSprint=0;
      break;
    case 9:
      doIdiot=0;
      doCrash=0;
      doSprint=-1;
      break;
    default:
      abort();
    }
  } else {
    // Dual
    switch (options.getSpecialOption(0)) {
    case 0:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      break;
    case 1:
      doIdiot=0;
      doCrash=1;
      doSprint=0;
      break;
    case 2:
      doIdiot=-1;
      if (options.getExtraInfo(0)>0)
	doIdiot = options.getExtraInfo(0);
      doCrash=0;
      doSprint=0;
      break;
    default:
      abort();
    }
  }
  // Just do this number of passes in Sprint
  int maxSprintPass=100;

  if (presolve!=ClpSolve::presolveOff) {
    int numberPasses=5;
    if (presolve==ClpSolve::presolveNumber) {
      numberPasses=options.getPresolvePasses();
      presolve = ClpSolve::presolveOn;
    }
    model2 = pinfo.presolvedModel(*this,1.0e-8,
				  false,numberPasses,true);
    time2 = CoinCpuTime();
    timePresolve = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Presolve"<<timePresolve<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    if (model2) {
    } else {
      handler_->message(CLP_INFEASIBLE,messages_)
	<<CoinMessageEol;
      model2 = this;
      presolve=ClpSolve::presolveOff;
    }
    // We may be better off using original (but if dual leave because of bounds)
    if (presolve!=ClpSolve::presolveOff&&
	numberRows_<1.01*model2->numberRows_&&numberColumns_<1.01*model2->numberColumns_
	&&method!=ClpSolve::useDual) {
      delete model2;
      model2 = this;
      presolve=ClpSolve::presolveOff;
    }
  }
  if (interrupt)
    currentModel = model2;
  // See if worth trying +- one matrix
  bool plusMinus=false;
  int numberElements=model2->getNumElements();
  if (dynamic_cast< ClpNetworkMatrix*>(matrix_)) {
    // network - switch off stuff
    doIdiot=0;
    doSprint=0;
  }
  int numberColumns = model2->numberColumns();
  int numberRows = model2->numberRows();
  // If not all slack basis - switch off all
  int number=0;
  int iRow;
  for (iRow=0;iRow<numberRows;iRow++)
    if (model2->getRowStatus(iRow)==basic)
      number++;
  if (number<numberRows) {
    doIdiot=0;
    doCrash=0;
    doSprint=0;
  }
  if (options.getSpecialOption(3)==0) {
    if(numberElements>100000)
      plusMinus=true;
    if(numberElements>10000&&(doIdiot||doSprint)) 
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
    model2->setFactorizationFrequency(100+model2->numberRows()/200);
  }
  if (method==ClpSolve::usePrimalorSprint) {
    if (doSprint<0) { 
      if (numberElements<500000) {
	// Small problem
	if(numberRows*10>numberColumns||numberColumns<6000
	   ||(numberRows*20>numberColumns&&!plusMinus))
	  method=ClpSolve::usePrimal; // switch off sprint
      } else {
	// larger problem
	if(numberRows*8>numberColumns)
	  method=ClpSolve::usePrimal; // switch off sprint
	// but make lightweight
	if(numberRows*10>numberColumns||numberColumns<6000
	   ||(numberRows*20>numberColumns&&!plusMinus))
	  maxSprintPass=5;
      }
    } else if (doSprint==0) {
      method=ClpSolve::usePrimal; // switch off sprint
    }
  }
  if (method==ClpSolve::useDual) {
    // pick up number passes
    int nPasses=0;
    int numberNotE=0;
    if ((doIdiot<0&&plusMinus)||doIdiot>0) {
      // See if candidate for idiot
      nPasses=0;
      Idiot info(*model2);
      // Get average number of elements per column
      double ratio  = ((double) numberElements/(double) numberColumns);
      // look at rhs
      int iRow;
      double largest=0.0;
      double smallest = 1.0e30;
      double largestGap=0.0;
      for (iRow=0;iRow<numberRows;iRow++) {
	double value1 = model2->rowLower_[iRow];
	if (value1&&value1>-1.0e31) {
	  largest = max(largest,fabs(value1));
	  smallest=min(smallest,fabs(value1));
	}
	double value2 = model2->rowUpper_[iRow];
	if (value2&&value2<1.0e31) {
	  largest = max(largest,fabs(value2));
	  smallest=min(smallest,fabs(value2));
	}
	if (value2>value1) {
	  numberNotE++;
	  if (value2>1.0e31||value1<-1.0e31)
	    largestGap = COIN_DBL_MAX;
	  else
	    largestGap = value2-value1;
	}
      }
      if (doIdiot<0) {
	if (numberRows>200&&numberColumns>5000&&ratio>=3.0&&
	    largest/smallest<1.1&&!numberNotE) {
	  nPasses = 71;
	}
      } 
      if (doIdiot>0) {
	nPasses=max(nPasses,doIdiot);
	if (nPasses>70) 
	  info.setStartingWeight(1.0e3);
      }
      if (nPasses) {
	info.setReduceIterations(5);
	doCrash=0;
	info.crash(nPasses,model2->messageHandler(),model2->messagesPointer());
	time2 = CoinCpuTime();
	timeIdiot = time2-timeX;
	handler_->message(CLP_INTERVAL_TIMING,messages_)
	  <<"Idiot Crash"<<timeIdiot<<time2-time1
	  <<CoinMessageEol;
	timeX=time2;
      }
    }
    int numberInfeasibilities = model2->tightenPrimalBounds();
    if (numberInfeasibilities) {
      handler_->message(CLP_INFEASIBLE,messages_)
	<<CoinMessageEol;
      model2 = this;
      presolve=ClpSolve::presolveOff;
    }
    if (options.getSpecialOption(0)==1)
      model2->crash(1000,1);
    if (!nPasses) {
      model2->dual(0);
    } else if (!numberNotE&&0) {
      // E so we can do in another way
      double * pi = model2->dualRowSolution();
      int i;
      int numberColumns = model2->numberColumns();
      int numberRows = model2->numberRows();
      double * saveObj = new double[numberColumns];
      memcpy(saveObj,model2->objective(),numberColumns*sizeof(double));
      memcpy(model2->dualColumnSolution(),model2->objective(),
	     numberColumns*sizeof(double));
      model2->clpMatrix()->transposeTimes(-1.0,pi,model2->dualColumnSolution());
      memcpy(model2->objective(),model2->dualColumnSolution(),
	     numberColumns*sizeof(double));
      const double * rowsol = model2->primalRowSolution();
      double offset=0.0;
      for (i=0;i<numberRows;i++) {
	offset += pi[i]*rowsol[i];
      }
      double value2;
      model2->getDblParam(ClpObjOffset,value2);
      printf("Offset %g %g\n",offset,value2);
      model2->setRowObjective(pi);
      // zero out pi
      memset(pi,0,numberRows*sizeof(double));
      // Could put some in basis - only partially tested
      model2->allSlackBasis(); 
      model2->factorization()->maximumPivots(200);
      //model2->setLogLevel(63);
      // solve
      model2->dual(1);
      memcpy(model2->objective(),saveObj,numberColumns*sizeof(double));
      // zero out pi
      memset(pi,0,numberRows*sizeof(double));
      model2->setRowObjective(pi);
      delete [] saveObj;
      model2->primal();
    } else {
      // solve
      model2->dual(1);
    }
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Dual"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  } else if (method==ClpSolve::usePrimal) {
#ifdef CLP_IDIOT
    if (doIdiot) {
      int nPasses=0;
      Idiot info(*model2);
      // Get average number of elements per column
      double ratio  = ((double) numberElements/(double) numberColumns);
      // look at rhs
      int iRow;
      double largest=0.0;
      double smallest = 1.0e30;
      double largestGap=0.0;
      int numberNotE=0;
      for (iRow=0;iRow<numberRows;iRow++) {
	double value1 = model2->rowLower_[iRow];
	if (value1&&value1>-1.0e31) {
	  largest = max(largest,fabs(value1));
	  smallest=min(smallest,fabs(value1));
	}
	double value2 = model2->rowUpper_[iRow];
	if (value2&&value2<1.0e31) {
	  largest = max(largest,fabs(value2));
	  smallest=min(smallest,fabs(value2));
	}
	if (value2>value1) {
	  numberNotE++;
	  if (value2>1.0e31||value1<-1.0e31)
	    largestGap = COIN_DBL_MAX;
	  else
	    largestGap = value2-value1;
	}
      }
      if (numberRows>200&&numberColumns>5000&&numberColumns>2*numberRows) {
	if (plusMinus) {
	  if (largest/smallest>2.0) {
	    nPasses = 10+numberColumns/100000;
	    nPasses = min(nPasses,50);
	    nPasses = max(nPasses,15);
	    if (numberRows>25000&&nPasses>5) {
	      // Might as well go for it
	      nPasses = max(nPasses,71);
	    } else if (numberElements<3*numberColumns) {
	      nPasses=min(nPasses,10); // probably not worh it
	    }
	    if (doIdiot<0)
	      info.setLightweight(1); // say lightweight idiot
	  } else if (largest/smallest>1.01||numberElements<=3*numberColumns) {
	    nPasses = 10+numberColumns/1000;
	    nPasses = min(nPasses,100);
	    nPasses = max(nPasses,30);
	    if (numberRows>25000) {
	      // Might as well go for it
	      nPasses = max(nPasses,71);
	    }
	    if (!largestGap)
	      nPasses *= 2;
	  } else {
	    nPasses = 10+numberColumns/1000;
	    nPasses = min(nPasses,200);
	    nPasses = max(nPasses,100);
	    info.setStartingWeight(1.0e-1);
	    info.setReduceIterations(6);
	    if (!largestGap)
	      nPasses *= 2;
	    //info.setFeasibilityTolerance(1.0e-7);
	  }
	  // If few passes - don't bother
	  if (nPasses<=5)
	    nPasses=0;
	} else {
	  if (doIdiot<0)
	    info.setLightweight(1); // say lightweight idiot
	  if (largest/smallest>1.01||numberNotE) {
	    if (numberRows>25000||numberColumns>5*numberRows) {
	      nPasses = 50;
	    } else if (numberColumns>4*numberRows) {
	      nPasses = 20;
	    } else {
	      nPasses=5;
	    }
	  } else {
	    if (numberRows>25000||numberColumns>5*numberRows) {
	      nPasses = 50;
	      info.setLightweight(0); // say not lightweight idiot
	    } else if (numberColumns>4*numberRows) {
	      nPasses = 20;
	    } else {
	      nPasses=15;
	    }
	  }
	  if (numberElements<3*numberColumns) { 
	    nPasses=(int) ((2.0*(double) nPasses)/ratio); // probably not worh it
	  } else {
	    nPasses = max(nPasses,5);
	    nPasses = (int) (((double) nPasses)*5.0/ratio); // reduce if lots of elements per column
	  }
	  if (numberRows>25000&&nPasses>5) {
	    // Might as well go for it
	    nPasses = max(nPasses,71);
	  } else if (plusMinus) {
	    nPasses *= 2;
	    nPasses=min(nPasses,71);
	  }
	  if (nPasses<=5)
	    nPasses=0;
	  //info.setStartingWeight(1.0e-1);
	}
      }
      if (doIdiot>0) {
	// pick up number passes
	nPasses=options.getExtraInfo(1);
	if (nPasses>70) {
	  info.setStartingWeight(1.0e3);
	  info.setReduceIterations(6);
	} else if (nPasses>=50) {
	  info.setStartingWeight(1.0e3);
	  //info.setReduceIterations(6);
	} 
	// For experimenting
	if (nPasses<70&&(nPasses%10)>0&&(nPasses%10)<4) {
	  info.setStartingWeight(1.0e3);
	  info.setLightweight(nPasses%10); // special testing
	  //info.setReduceIterations(6);
	}
      }
      if (nPasses) {
	doCrash=0;
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
    // ?
    if (doCrash)
      model2->crash(1000,1);
    model2->primal(1);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    timeSimplex = timeCore;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Primal"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  } else if (method==ClpSolve::usePrimalorSprint) {
    // Sprint
    /*
      This driver implements what I called Sprint when I introduced the idea
      many years ago.  Cplex calls it "sifting" which I think is just as silly.
      When I thought of this trivial idea
      it reminded me of an LP code of the 60's called sprint which after
      every factorization took a subset of the matrix into memory (all
      64K words!) and then iterated very fast on that subset.  On the
      problems of those days it did not work very well, but it worked very
      well on aircrew scheduling problems where there were very large numbers
      of columns all with the same flavor.
    */
    
    /* The idea works best if you can get feasible easily.  To make it
       more general we can add in costed slacks */
    
    int originalNumberColumns = model2->numberColumns();
    int numberRows = model2->numberRows();
    
    // We will need arrays to choose variables.  These are too big but ..
    double * weight = new double [numberRows+originalNumberColumns];
    int * sort = new int [numberRows+originalNumberColumns];
    int numberSort=0;
    // We are going to add slacks to get feasible.
    // initial list will just be artificials
    // first we will set all variables as close to zero as possible
    int iColumn;
    const double * columnLower = model2->columnLower();
    const double * columnUpper = model2->columnUpper();
    double * columnSolution = model2->primalColumnSolution();
    
    for (iColumn=0;iColumn<originalNumberColumns;iColumn++) {
      double value =0.0;
      if (columnLower[iColumn]>0.0)
	value = columnLower[iColumn];
      else if (columnUpper[iColumn]<0.0)
	value = columnUpper[iColumn];
      columnSolution[iColumn]=value;
    }
    // now see what that does to row solution
    double * rowSolution = model2->primalRowSolution();
    memset (rowSolution,0,numberRows*sizeof(double));
    model2->times(1.0,columnSolution,rowSolution);
    
    int * addStarts = new int [numberRows+1];
    int * addRow = new int[numberRows];
    double * addElement = new double[numberRows];
    const double * lower = model2->rowLower();
    const double * upper = model2->rowUpper();
    addStarts[0]=0;
    int numberArtificials=0;
    double * addCost = new double [numberRows];
    const double penalty=1.0e8;
    int iRow;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (lower[iRow]>rowSolution[iRow]) {
	addRow[numberArtificials]=iRow;
	addElement[numberArtificials]=1.0;
	addCost[numberArtificials]=penalty;
	numberArtificials++;
	addStarts[numberArtificials]=numberArtificials;
      } else if (upper[iRow]<rowSolution[iRow]) {
	addRow[numberArtificials]=iRow;
	addElement[numberArtificials]=-1.0;
	addCost[numberArtificials]=penalty;
	numberArtificials++;
	addStarts[numberArtificials]=numberArtificials;
      }
    }
    model2->addColumns(numberArtificials,NULL,NULL,addCost,
		       addStarts,addRow,addElement);
    delete [] addStarts;
    delete [] addRow;
    delete [] addElement;
    delete [] addCost;
    // look at rhs to see if to perturb
    double largest=0.0;
    double smallest = 1.0e30;
    for (iRow=0;iRow<numberRows;iRow++) {
      double value;
      value = fabs(model2->rowLower_[iRow]);
      if (value&&value<1.0e30) {
	largest = max(largest,value);
	smallest=min(smallest,value);
      }
      value = fabs(model2->rowUpper_[iRow]);
      if (value&&value<1.0e30) {
	largest = max(largest,value);
	smallest=min(smallest,value);
      }
    }
    double * saveLower = NULL;
    double * saveUpper = NULL;
    if (largest<2.01*smallest) {
      // perturb - so switch off standard
      model2->setPerturbation(100);
      saveLower = new double[numberRows];
      memcpy(saveLower,model2->rowLower_,numberRows*sizeof(double));
      saveUpper = new double[numberRows];
      memcpy(saveUpper,model2->rowUpper_,numberRows*sizeof(double));
      double * lower = model2->rowLower();
      double * upper = model2->rowUpper();
      for (iRow=0;iRow<numberRows;iRow++) {
	double lowerValue=lower[iRow], upperValue=upper[iRow];
	double value = CoinDrand48();
	if (upperValue>lowerValue+primalTolerance_) {
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= value * 1.0e-4*fabs(lowerValue); 
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += value * 1.0e-4*fabs(upperValue); 
	} else if (upperValue>0.0) {
	  upperValue -= value * 1.0e-4*fabs(lowerValue); 
	  lowerValue -= value * 1.0e-4*fabs(lowerValue); 
	} else if (upperValue<0.0) {
	  upperValue += value * 1.0e-4*fabs(lowerValue); 
	  lowerValue += value * 1.0e-4*fabs(lowerValue); 
	} else {
	}
	lower[iRow]=lowerValue;
	upper[iRow]=upperValue;
      }
    }
    int i;
    // Set up initial list
    if (numberArtificials) {
      numberSort=numberArtificials;
      for (i=0;i<numberSort;i++)
	sort[i] = i+originalNumberColumns;
    } else {
      numberSort = min(numberRows_,numberColumns_);
      for (i=0;i<numberSort;i++)
	sort[i] = i;
    }
    
    // redo as will have changed
    columnLower = model2->columnLower();
    columnUpper = model2->columnUpper();
    int numberColumns = model2->numberColumns();
    double * fullSolution = model2->primalColumnSolution();
    
    // Just do this number of passes in Sprint
    if (doSprint>0)
      maxSprintPass=options.getExtraInfo(1);
    int iPass;
    double lastObjective=1.0e31;
    // It will be safe to allow dense
    model2->setInitialDenseFactorization(true);
    
    // Just take this number of columns in small problem
    int smallNumberColumns = min(3*numberRows,numberColumns);
    smallNumberColumns = max(smallNumberColumns,3000);
    //int smallNumberColumns = min(12*numberRows/10,numberColumns);
    //smallNumberColumns = max(smallNumberColumns,3000);
    //smallNumberColumns = max(smallNumberColumns,numberRows+1000);
    // We will be using all rows
    int * whichRows = new int [numberRows];
    for (int iRow=0;iRow<numberRows;iRow++)
      whichRows[iRow]=iRow;
    double originalOffset;
    model2->getDblParam(ClpObjOffset,originalOffset);
    int totalIterations=0;
    for (iPass=0;iPass<maxSprintPass;iPass++) {
      //printf("Bug until submodel new version\n");
      //CoinSort_2(sort,sort+numberSort,weight);
      // Create small problem
      ClpSimplex small(model2,numberRows,whichRows,numberSort,sort);
      small.setPerturbation(model2->perturbation());
      // now see what variables left out do to row solution
      double * rowSolution = model2->primalRowSolution();
      double * sumFixed = new double[numberRows];
      memset (sumFixed,0,numberRows*sizeof(double));
      int iRow,iColumn;
      // zero out ones in small problem
      for (iColumn=0;iColumn<numberSort;iColumn++) {
	int kColumn = sort[iColumn];
	fullSolution[kColumn]=0.0;
      }
      // Get objective offset
      double offset=0.0;
      const double * objective = model2->objective();
      for (iColumn=0;iColumn<numberColumns;iColumn++) 
	offset += fullSolution[iColumn]*objective[iColumn];
      small.setDblParam(ClpObjOffset,originalOffset-offset);
      model2->times(1.0,fullSolution,sumFixed);
      
      double * lower = small.rowLower();
      double * upper = small.rowUpper();
      for (iRow=0;iRow<numberRows;iRow++) {
	if (lower[iRow]>-1.0e50) 
	  lower[iRow] -= sumFixed[iRow];
	if (upper[iRow]<1.0e50)
	  upper[iRow] -= sumFixed[iRow];
	rowSolution[iRow] -= sumFixed[iRow];
      }
      delete [] sumFixed;
      // Solve 
      if (interrupt)
	currentModel = &small;
      small.primal();
      totalIterations += small.numberIterations();
      // move solution back
      const double * solution = small.primalColumnSolution();
      for (iColumn=0;iColumn<numberSort;iColumn++) {
	int kColumn = sort[iColumn];
	model2->setColumnStatus(kColumn,small.getColumnStatus(iColumn));
	fullSolution[kColumn]=solution[iColumn];
      }
      for (iRow=0;iRow<numberRows;iRow++) 
	model2->setRowStatus(iRow,small.getRowStatus(iRow));
      memcpy(model2->primalRowSolution(),small.primalRowSolution(),
	     numberRows*sizeof(double));
      // get reduced cost for large problem
      memcpy(weight,model2->objective(),numberColumns*sizeof(double));
      model2->transposeTimes(-1.0,small.dualRowSolution(),weight);
      int numberNegative=0;
      double sumNegative = 0.0;
      // now massage weight so all basic in plus good djs
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double dj = weight[iColumn]*optimizationDirection_;
	double value = fullSolution[iColumn];
	if (model2->getColumnStatus(iColumn)==ClpSimplex::basic) 
	  dj = -1.0e50;
	else if (dj<0.0&&value<columnUpper[iColumn])
	  dj = dj;
	else if (dj>0.0&&value>columnLower[iColumn])
	  dj = -dj;
	else if (columnUpper[iColumn]>columnLower[iColumn])
	  dj = fabs(dj);
	else
	  dj = 1.0e50;
	weight[iColumn] = dj;
	if (dj<-dualTolerance_&&dj>-1.0e50) {
	  numberNegative++;
	  sumNegative -= dj;
	}
	sort[iColumn] = iColumn;
      }
      handler_->message(CLP_SPRINT,messages_)
	<<iPass+1<<small.numberIterations()<<small.objectiveValue()<<sumNegative
	<<numberNegative
	<<CoinMessageEol;
      if ((small.objectiveValue()*optimizationDirection_>lastObjective-1.0e-7&&iPass>5)||
	  !small.numberIterations()||
	  iPass==maxSprintPass-1||small.status()==3) {
	
	break; // finished
      } else {
	lastObjective = small.objectiveValue()*optimizationDirection_;
	// sort
	CoinSort_2(weight,weight+numberColumns,sort);
	numberSort = smallNumberColumns;
      }
    }
    if (interrupt) 
      currentModel = model2;
    for (i=0;i<numberArtificials;i++)
      sort[i] = i + originalNumberColumns;
    model2->deleteColumns(numberArtificials,sort);
    delete [] weight;
    delete [] sort;
    delete [] whichRows;
    if (saveLower) {
      // unperturb and clean
      for (iRow=0;iRow<numberRows;iRow++) {
	double diffLower = saveLower[iRow]-model2->rowLower_[iRow];
	double diffUpper = saveUpper[iRow]-model2->rowUpper_[iRow];
	model2->rowLower_[iRow]=saveLower[iRow];
	model2->rowUpper_[iRow]=saveUpper[iRow];
	if (diffLower) 
	  assert (!diffUpper||fabs(diffLower-diffUpper)<1.0e-5);
	else
	  diffLower = diffUpper;
	model2->rowActivity_[iRow] += diffLower;
      }
      delete [] saveLower;
      delete [] saveUpper;
    }
    model2->primal(0);
    model2->setPerturbation(savePerturbation);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    timeSimplex = timeCore;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Sprint"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    model2->setNumberIterations(model2->numberIterations()+totalIterations);
  } else {
    assert (method!=ClpSolve::automatic); // later
    time2=0.0;
  }
  if (saveMatrix) {
    // delete and replace
    delete model2->clpMatrix();
    model2->replaceMatrix(saveMatrix);
  }
  numberIterations = model2->numberIterations();
  finalStatus=model2->status();
  if (presolve==ClpSolve::presolveOn) {
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
  handler_->printing(presolve==ClpSolve::presolveOn)
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
// General solve
int 
ClpSimplex::initialSolve()
{
  // Default so use dual
  ClpSolve options;
  return initialSolve(options);
}
// General dual solve
int 
ClpSimplex::initialDualSolve()
{
  ClpSolve options;
  // Use dual
  options.setSolveType(ClpSolve::useDual);
  return initialSolve(options);
}
// General dual solve
int 
ClpSimplex::initialPrimalSolve()
{
  ClpSolve options;
  // Use primal
  options.setSolveType(ClpSolve::usePrimal);
  return initialSolve(options);
}
// Default constructor
ClpSolve::ClpSolve (  )
{
  method_ = useDual;
  presolveType_=presolveOn;
  numberPasses_=5;
  for (int i=0;i<4;i++)
    options_[i]=0;
  for (int i=0;i<4;i++)
    extraInfo_[i]=-1;
}

// Copy constructor. 
ClpSolve::ClpSolve(const ClpSolve & rhs)
{
  method_ = rhs.method_;
  presolveType_=rhs.presolveType_;
  numberPasses_=rhs.numberPasses_;
  for (int i=0;i<4;i++)
    options_[i]=rhs.options_[i];
  for (int i=0;i<4;i++)
    extraInfo_[i]=rhs.extraInfo_[i];
}
// Assignment operator. This copies the data
ClpSolve & 
ClpSolve::operator=(const ClpSolve & rhs)
{
  if (this != &rhs) {
    method_ = rhs.method_;
    presolveType_=rhs.presolveType_;
    numberPasses_=rhs.numberPasses_;
    for (int i=0;i<4;i++)
      options_[i]=rhs.options_[i];
    for (int i=0;i<4;i++)
      extraInfo_[i]=rhs.extraInfo_[i];
  }
  return *this;

}
// Destructor
ClpSolve::~ClpSolve (  )
{
}
/*   which translation is:
     which:
     0 - startup in Dual  (nothing if basis exists).:
             0 - no basis, 1 crash
     1 - startup in Primal (nothing if basis exists):
        0 - use initiative
        1 - use crash
        2 - use idiot and look at further info
        3 - use sprint and look at further info
        4 - use all slack
        5 - use initiative but no idiot
        6 - use initiative but no sprint
        7 - use initiative but no crash
        8 - do allslack or idiot
        9 - do allslack or sprint
     2 - interrupt handling - 0 yes, 1 no (for threadsafe)
     3 - whether to make +- 1matrix - 0 yes, 1 no
*/
void 
ClpSolve::setSpecialOption(int which,int value,int extraInfo)
{
  options_[which]=value;
  extraInfo_[which]=extraInfo;
}
int 
ClpSolve::getSpecialOption(int which) const
{
  return options_[which];
}

// Solve types
void 
ClpSolve::setSolveType(SolveType method, int extraInfo)
{
  method_=method;
}

ClpSolve::SolveType 
ClpSolve::getSolveType()
{
  return method_;
}

// Presolve types
void 
ClpSolve::setPresolveType(PresolveType amount, int extraInfo)
{
  presolveType_=amount;
  numberPasses_=extraInfo;
}
ClpSolve::PresolveType 
ClpSolve::getPresolveType()
{
  return presolveType_;
}
// Extra info for idiot (or sprint)
int 
ClpSolve::getExtraInfo(int which) const
{
  return extraInfo_[which];
}
int 
ClpSolve::getPresolvePasses() const
{
  return numberPasses_;
}
