// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#include "ClpDualRowSteepest.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "CoinHelperFunctions.hpp"
#include <cstdio>
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpDualRowSteepest::ClpDualRowSteepest (int mode) 
  : ClpDualRowPivot(),
    state_(-1),
    mode_(mode),
    persistence_(normal),
    weights_(NULL),
    infeasible_(NULL),
    alternateWeights_(NULL),
    savedWeights_(NULL),
    dubiousWeights_(NULL)
{
  type_=2+64*mode;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpDualRowSteepest::ClpDualRowSteepest (const ClpDualRowSteepest & rhs) 
: ClpDualRowPivot(rhs)
{  
  state_=rhs.state_;
  mode_ = rhs.mode_;
  persistence_ = rhs.persistence_;
  model_ = rhs.model_;
  if (rhs.infeasible_) {
    infeasible_= new CoinIndexedVector(rhs.infeasible_);
  } else {
    infeasible_=NULL;
  }
  if (rhs.weights_) {
    assert(model_);
    int number = model_->numberRows();
    weights_= new double[number];
    ClpDisjointCopyN(rhs.weights_,number,weights_);
  } else {
    weights_=NULL;
  }
  if (rhs.alternateWeights_) {
    alternateWeights_= new CoinIndexedVector(rhs.alternateWeights_);
  } else {
    alternateWeights_=NULL;
  }
  if (rhs.savedWeights_) {
    savedWeights_= new CoinIndexedVector(rhs.savedWeights_);
  } else {
    savedWeights_=NULL;
  }
  if (rhs.dubiousWeights_) {
    assert(model_);
    int number = model_->numberRows();
    dubiousWeights_= new int[number];
    ClpDisjointCopyN(rhs.dubiousWeights_,number,dubiousWeights_);
  } else {
    dubiousWeights_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpDualRowSteepest::~ClpDualRowSteepest ()
{
  delete [] weights_;
  delete [] dubiousWeights_;
  delete infeasible_;
  delete alternateWeights_;
  delete savedWeights_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpDualRowSteepest &
ClpDualRowSteepest::operator=(const ClpDualRowSteepest& rhs)
{
  if (this != &rhs) {
    ClpDualRowPivot::operator=(rhs);
    state_=rhs.state_;
    mode_ = rhs.mode_;
    persistence_ = rhs.persistence_;
    model_ = rhs.model_;
    delete [] weights_;
    delete [] dubiousWeights_;
    delete infeasible_;
    delete alternateWeights_;
    delete savedWeights_;
    if (rhs.infeasible_!=NULL) {
      infeasible_= new CoinIndexedVector(rhs.infeasible_);
    } else {
      infeasible_=NULL;
    }
    if (rhs.weights_!=NULL) {
      assert(model_);
      int number = model_->numberRows();
      weights_= new double[number];
      ClpDisjointCopyN(rhs.weights_,number,weights_);
    } else {
      weights_=NULL;
    }
    if (rhs.alternateWeights_!=NULL) {
      alternateWeights_= new CoinIndexedVector(rhs.alternateWeights_);
    } else {
      alternateWeights_=NULL;
    }
    if (rhs.savedWeights_!=NULL) {
      savedWeights_= new CoinIndexedVector(rhs.savedWeights_);
    } else {
      savedWeights_=NULL;
    }
    if (rhs.dubiousWeights_) {
      assert(model_);
      int number = model_->numberRows();
      dubiousWeights_= new int[number];
      ClpDisjointCopyN(rhs.dubiousWeights_,number,dubiousWeights_);
    } else {
      dubiousWeights_=NULL;
    }
  }
  return *this;
}
// Returns pivot row, -1 if none
int 
ClpDualRowSteepest::pivotRow()
{
  assert(model_);
  int i,iRow;
  double * infeas = infeasible_->denseVector();
  double largest=0.0;
  int * index = infeasible_->getIndices();
  int number = infeasible_->getNumElements();
  const int * pivotVariable =model_->pivotVariable();
  int chosenRow=-1;
  int lastPivotRow = model_->pivotRow();
  double tolerance=model_->currentPrimalTolerance();
  // we can't really trust infeasibilities if there is primal error
  // this coding has to mimic coding in checkPrimalSolution
  double error = CoinMin(1.0e-3,model_->largestPrimalError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance +  error;
  // But cap
  tolerance = CoinMin(1000.0,tolerance);
  tolerance *= tolerance; // as we are using squares
  double * solution = model_->solutionRegion();
  double * lower = model_->lowerRegion();
  double * upper = model_->upperRegion();
  // do last pivot row one here
  //#define COLUMN_BIAS 4.0
  //#define FIXED_BIAS 10.0
  if (lastPivotRow>=0) {
#ifdef COLUMN_BIAS 
    int numberColumns = model_->numberColumns();
#endif
    int iPivot=pivotVariable[lastPivotRow];
    double value = solution[iPivot];
    double lower = model_->lower(iPivot);
    double upper = model_->upper(iPivot);
    if (value>upper+tolerance) {
      value -= upper;
      value *= value;
#ifdef COLUMN_BIAS 
      if (iPivot<numberColumns)
	value *= COLUMN_BIAS; // bias towards columns
k
#endif
      // store square in list
      if (infeas[lastPivotRow])
	infeas[lastPivotRow] = value; // already there
      else
	infeasible_->quickAdd(lastPivotRow,value);
    } else if (value<lower-tolerance) {
      value -= lower;
      value *= value;
#ifdef COLUMN_BIAS 
      if (iPivot<numberColumns)
	value *= COLUMN_BIAS; // bias towards columns
#endif
      // store square in list
      if (infeas[lastPivotRow])
	infeas[lastPivotRow] = value; // already there
      else
	infeasible_->add(lastPivotRow,value);
    } else {
      // feasible - was it infeasible - if so set tiny
      if (infeas[lastPivotRow])
	infeas[lastPivotRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
    }
    number = infeasible_->getNumElements();
  }
  if(model_->numberIterations()<model_->lastBadIteration()+200) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (!model_->factorization()->pivots())
      checkTolerance = 1.0e-6;
    if (model_->largestPrimalError()>checkTolerance)
      tolerance *= model_->largestPrimalError()/checkTolerance;
  }
  int numberWanted;
  if (mode_<2 ) {
    numberWanted = number+1;
  } else if (mode_==2) {
    numberWanted = CoinMax(2000,number/8);
  } else {
    int numberElements = model_->factorization()->numberElements();
    double ratio = (double) numberElements/(double) model_->numberRows();
    numberWanted = CoinMax(2000,number/8);
    if (ratio<1.0) {
      numberWanted = CoinMax(2000,number/20);
    } else if (ratio>10.0) {
      ratio = number * (ratio/80.0);
      if (ratio>number)
	numberWanted=number+1;
      else
	numberWanted = CoinMax(2000,(int) ratio);
    }
  }
  int iPass;
  // Setup two passes
  int start[4];
  start[1]=number;
  start[2]=0;
  double dstart = ((double) number) * CoinDrand48();
  start[0]=(int) dstart;
  start[3]=start[0];
  //double largestWeight=0.0;
  //double smallestWeight=1.0e100;
  for (iPass=0;iPass<2;iPass++) {
    int end = start[2*iPass+1];
    for (i=start[2*iPass];i<end;i++) {
      iRow = index[i];
      double value = infeas[iRow];
      if (value>tolerance) {
	//#define OUT_EQ
#ifdef OUT_EQ
	{
	  int iSequence = pivotVariable[iRow];
	  if (upper[iSequence]==lower[iSequence])
	    value *= 2.0;
	}
#endif
	double weight = CoinMin(weights_[iRow],1.0e50);
	//largestWeight = CoinMax(largestWeight,weight);
	//smallestWeight = CoinMin(smallestWeight,weight);
	//double dubious = dubiousWeights_[iRow];
	//weight *= dubious;
	//if (value>2.0*largest*weight||(value>0.5*largest*weight&&value*largestWeight>dubious*largest*weight)) {
	if (value>largest*weight) {
	  // make last pivot row last resort choice
	  if (iRow==lastPivotRow) {
	    if (value*1.0e-10<largest*weight) 
	      continue;
	    else 
	      value *= 1.0e-10;
	  }
	  int iSequence = pivotVariable[iRow];
	  if (!model_->flagged(iSequence)) {
	    //#define CLP_DEBUG 3
#ifdef CLP_DEBUG
	    double value2=0.0;
	    if (solution[iSequence]>upper[iSequence]+tolerance) 
	      value2=solution[iSequence]-upper[iSequence];
	    else if (solution[iSequence]<lower[iSequence]-tolerance) 
	      value2=solution[iSequence]-lower[iSequence];
	    assert(fabs(value2*value2-infeas[iRow])<1.0e-8*CoinMin(value2*value2,infeas[iRow]));
#endif
	    if (solution[iSequence]>upper[iSequence]+tolerance||
		solution[iSequence]<lower[iSequence]-tolerance) {
	      chosenRow=iRow;
	      largest=value/weight;
	      //largestWeight = dubious;
	    }
	  } else {
	    // just to make sure we don't exit before got something
	    numberWanted++;
	  }
	}
	numberWanted--;
	if (!numberWanted)
	  break;
      }
    }
    if (!numberWanted)
      break;
  }
  //printf("smallest %g largest %g\n",smallestWeight,largestWeight);
  return chosenRow;
}
// Updates weights and returns pivot alpha
double
ClpDualRowSteepest::updateWeights(CoinIndexedVector * input,
				  CoinIndexedVector * spare,
				  CoinIndexedVector * updatedColumn)
{
#if CLP_DEBUG>2
  // Very expensive debug
  {
    int numberRows = model_->numberRows();
    CoinIndexedVector * temp = new CoinIndexedVector();
    temp->reserve(numberRows+
		  model_->factorization()->maximumPivots());
    double * array = alternateWeights_->denseVector();
    int * which = alternateWeights_->getIndices();
    int i;
    for (i=0;i<numberRows;i++) {
      double value=0.0;
      array[i] = 1.0;
      which[0] = i;
      alternateWeights_->setNumElements(1);
      model_->factorization()->updateColumnTranspose(temp,
						     alternateWeights_);
      int number = alternateWeights_->getNumElements();
      int j;
      for (j=0;j<number;j++) {
	int iRow=which[j];
	value+=array[iRow]*array[iRow];
	array[iRow]=0.0;
      }
      alternateWeights_->setNumElements(0);
      double w = CoinMax(weights_[i],value)*.1;
      if (fabs(weights_[i]-value)>w) {
	printf("%d old %g, true %g\n",i,weights_[i],value);
	weights_[i]=value; // to reduce printout
      }
      //else 
      //printf("%d matches %g\n",i,value);
    }
    delete temp;
  }
#endif
  assert (input->packedMode());
  assert (updatedColumn->packedMode());
  double alpha=0.0;
  if (!model_->factorization()->networkBasis()) {
    // clear other region
    alternateWeights_->clear();
    double norm = 0.0;
    int i;
    double * work = input->denseVector();
    int numberNonZero = input->getNumElements();
    int * which = input->getIndices();
    double * work2 = spare->denseVector();
    int * which2 = spare->getIndices();
    // ftran
    //permute and move indices into index array
    //also compute norm
    //int *regionIndex = alternateWeights_->getIndices (  );
    const int *permute = model_->factorization()->permute();
    //double * region = alternateWeights_->denseVector();
    for ( i = 0; i < numberNonZero; i ++ ) {
      int iRow = which[i];
      double value = work[i];
      norm += value*value;
      iRow = permute[iRow];
      work2[iRow] = value;
      which2[i] = iRow;
    }
    spare->setNumElements ( numberNonZero );
    // Only one array active as already permuted
    model_->factorization()->updateColumn(spare,spare,true);
    // permute back
    numberNonZero = spare->getNumElements();
    // alternateWeights_ should still be empty
    int pivotRow = model_->pivotRow();
#ifdef CLP_DEBUG
    if ( model_->logLevel (  ) >4  && 
	 fabs(norm-weights_[pivotRow])>1.0e-3*(1.0+norm)) 
      printf("on row %d, true weight %g, old %g\n",
	     pivotRow,sqrt(norm),sqrt(weights_[pivotRow]));
#endif
    // could re-initialize here (could be expensive)
    norm /= model_->alpha() * model_->alpha();
    
    assert(norm);
    // pivot element
    alpha=0.0;
    double multiplier = 2.0 / model_->alpha();
    // look at updated column
    work = updatedColumn->denseVector();
    numberNonZero = updatedColumn->getNumElements();
    which = updatedColumn->getIndices();
    
    int nSave=0;
    double * work3 = alternateWeights_->denseVector();
    int * which3 = alternateWeights_->getIndices();
    const int * pivotColumn = model_->factorization()->pivotColumn();
    for (i =0; i < numberNonZero; i++) {
      int iRow = which[i];
      double theta = work[i];
      if (iRow==pivotRow)
	alpha = theta;
      double devex = weights_[iRow];
      work3[nSave]=devex; // save old
      which3[nSave++]=iRow;
      // transform to match spare
      int jRow = pivotColumn[iRow];
      double value = work2[jRow];
      devex +=  theta * (theta*norm+value * multiplier);
      if (devex < DEVEX_TRY_NORM) 
	devex = DEVEX_TRY_NORM;
      weights_[iRow]=devex;
    }
    alternateWeights_->setPackedMode(true);
    alternateWeights_->setNumElements(nSave);
    if (norm < DEVEX_TRY_NORM) 
      norm = DEVEX_TRY_NORM;
    // Try this to make less likely will happen again and stop cycling
    //norm *= 1.02;
    weights_[pivotRow] = norm;
    spare->clear();
#ifdef CLP_DEBUG
    spare->checkClear();
#endif
  } else {
    // clear other region
    alternateWeights_->clear();
    double norm = 0.0;
    int i;
    double * work = input->denseVector();
    int number = input->getNumElements();
    int * which = input->getIndices();
    double * work2 = spare->denseVector();
    int * which2 = spare->getIndices();
    for (i=0;i<number;i++) {
      int iRow = which[i];
      double value = work[i];
      norm += value*value;
      work2[iRow]=value;
      which2[i]=iRow;
    }
    spare->setNumElements(number);
    // ftran
    model_->factorization()->updateColumn(alternateWeights_,spare);
    // alternateWeights_ should still be empty
    int pivotRow = model_->pivotRow();
#ifdef CLP_DEBUG
    if ( model_->logLevel (  ) >4  && 
	 fabs(norm-weights_[pivotRow])>1.0e-3*(1.0+norm)) 
      printf("on row %d, true weight %g, old %g\n",
	     pivotRow,sqrt(norm),sqrt(weights_[pivotRow]));
#endif
    // could re-initialize here (could be expensive)
    norm /= model_->alpha() * model_->alpha();
    
    assert(norm);
    //if (norm < DEVEX_TRY_NORM) 
    //norm = DEVEX_TRY_NORM;
    // pivot element
    alpha=0.0;
    double multiplier = 2.0 / model_->alpha();
    // look at updated column
    work = updatedColumn->denseVector();
    number = updatedColumn->getNumElements();
    which = updatedColumn->getIndices();
    
    int nSave=0;
    double * work3 = alternateWeights_->denseVector();
    int * which3 = alternateWeights_->getIndices();
    for (i =0; i < number; i++) {
      int iRow = which[i];
      double theta = work[i];
      if (iRow==pivotRow)
	alpha = theta;
      double devex = weights_[iRow];
      work3[nSave]=devex; // save old
      which3[nSave++]=iRow;
      double value = work2[iRow];
      devex +=  theta * (theta*norm+value * multiplier);
      if (devex < DEVEX_TRY_NORM) 
	devex = DEVEX_TRY_NORM;
      weights_[iRow]=devex;
    }
    assert (alpha);
    alternateWeights_->setPackedMode(true);
    alternateWeights_->setNumElements(nSave);
    if (norm < DEVEX_TRY_NORM) 
      norm = DEVEX_TRY_NORM;
    weights_[pivotRow] = norm;
    spare->clear();
  }
#ifdef CLP_DEBUG
  spare->checkClear();
#endif
  return alpha;
}
  
/* Updates primal solution (and maybe list of candidates)
   Uses input vector which it deletes
   Computes change in objective function
*/
void 
ClpDualRowSteepest::updatePrimalSolution(
					CoinIndexedVector * primalUpdate,
					double primalRatio,
					double & objectiveChange)
{
  double * work = primalUpdate->denseVector();
  int number = primalUpdate->getNumElements();
  int * which = primalUpdate->getIndices();
  int i;
  double changeObj=0.0;
  double tolerance=model_->currentPrimalTolerance();
  const int * pivotVariable = model_->pivotVariable();
  double * infeas = infeasible_->denseVector();
  int pivotRow = model_->pivotRow();
  double * solution = model_->solutionRegion();
#ifdef COLUMN_BIAS 
  int numberColumns = model_->numberColumns();
#endif
  if (primalUpdate->packedMode()) {
    for (i=0;i<number;i++) {
      int iRow=which[i];
      int iPivot=pivotVariable[iRow];
      double value = solution[iPivot];
      double cost = model_->cost(iPivot);
      double change = primalRatio*work[i];
      work[i]=0.0;
      value -= change;
      changeObj -= change*cost;
      solution[iPivot] = value;
      double lower = model_->lower(iPivot);
      double upper = model_->upper(iPivot);
      // But if pivot row then use value of incoming
      // Although it is safer to recompute before next selection
      // in case something odd happens
      if (iRow==pivotRow) {
	iPivot = model_->sequenceIn();
	lower = model_->lower(iPivot);
	upper = model_->upper(iPivot);
	value = model_->valueIncomingDual();
      }
      if (value<lower-tolerance) {
	value -= lower;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	if (infeas[iRow])
	  infeas[iRow] = value; // already there
	else
	  infeasible_->quickAdd(iRow,value);
      } else if (value>upper+tolerance) {
	value -= upper;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	if (infeas[iRow])
	  infeas[iRow] = value; // already there
	else
	  infeasible_->quickAdd(iRow,value);
      } else {
	// feasible - was it infeasible - if so set tiny
	if (infeas[iRow])
	  infeas[iRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
      }
    }
  } else {
    for (i=0;i<number;i++) {
      int iRow=which[i];
      int iPivot=pivotVariable[iRow];
      double value = solution[iPivot];
      double cost = model_->cost(iPivot);
      double change = primalRatio*work[iRow];
      value -= change;
      changeObj -= change*cost;
      solution[iPivot] = value;
      double lower = model_->lower(iPivot);
      double upper = model_->upper(iPivot);
      // But if pivot row then use value of incoming
      // Although it is safer to recompute before next selection
      // in case something odd happens
      if (iRow==pivotRow) {
	iPivot = model_->sequenceIn();
	lower = model_->lower(iPivot);
	upper = model_->upper(iPivot);
	value = model_->valueIncomingDual();
      }
      if (value<lower-tolerance) {
	value -= lower;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	if (infeas[iRow])
	  infeas[iRow] = value; // already there
	else
	  infeasible_->quickAdd(iRow,value);
      } else if (value>upper+tolerance) {
	value -= upper;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	if (infeas[iRow])
	  infeas[iRow] = value; // already there
	else
	  infeasible_->quickAdd(iRow,value);
      } else {
	// feasible - was it infeasible - if so set tiny
	if (infeas[iRow])
	  infeas[iRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
      }
      work[iRow]=0.0;
    }
  }
  primalUpdate->setNumElements(0);
  objectiveChange += changeObj;
}
/* Saves any weights round factorization as pivot rows may change
   1) before factorization
   2) after factorization
   3) just redo infeasibilities
   4) restore weights
*/
void 
ClpDualRowSteepest::saveWeights(ClpSimplex * model,int mode)
{
  // alternateWeights_ is defined as indexed but is treated oddly
  model_ = model;
  int numberRows = model_->numberRows();
  int numberColumns = model_->numberColumns();
  const int * pivotVariable = model_->pivotVariable();
  int i;
  if (mode==1) {
    if(weights_) {
      // Check if size has changed
      if (infeasible_->capacity()==numberRows) {
	alternateWeights_->clear();
	// change from row numbers to sequence numbers
	int * which = alternateWeights_->getIndices();
	for (i=0;i<numberRows;i++) {
	  int iPivot=pivotVariable[i];
	  which[i]=iPivot;
	}
	state_=1;
      } else {
	// size has changed - clear everything
	delete [] weights_;
	weights_=NULL;
	delete [] dubiousWeights_;
	dubiousWeights_=NULL;
	delete infeasible_;
	infeasible_=NULL;
	delete alternateWeights_;
	alternateWeights_=NULL;
	delete savedWeights_;
	savedWeights_=NULL;
	state_=-1;
      }
    }
  } else if (mode==2||mode==4||mode==5) {
    // restore
    if (!weights_||state_==-1||mode==5) {
      // initialize weights
      delete [] weights_;
      delete alternateWeights_;
      weights_ = new double[numberRows];
      alternateWeights_ = new CoinIndexedVector();
      // enough space so can use it for factorization
      alternateWeights_->reserve(numberRows+
				 model_->factorization()->maximumPivots());
      if (mode_!=1||mode==5) {
	// initialize to 1.0 (can we do better?)
	for (i=0;i<numberRows;i++) {
	  weights_[i]=1.0;
	}
      } else {
	CoinIndexedVector * temp = new CoinIndexedVector();
	temp->reserve(numberRows+
		      model_->factorization()->maximumPivots());
	double * array = alternateWeights_->denseVector();
	int * which = alternateWeights_->getIndices();
	for (i=0;i<numberRows;i++) {
	  double value=0.0;
	  array[0] = 1.0;
	  which[0] = i;
	  alternateWeights_->setNumElements(1);
	  alternateWeights_->setPackedMode(true);
	  model_->factorization()->updateColumnTranspose(temp,
							 alternateWeights_);
	  int number = alternateWeights_->getNumElements();
	  int j;
	  for (j=0;j<number;j++) {
	    value+=array[j]*array[j];
	    array[j]=0.0;
	  }
	  alternateWeights_->setNumElements(0);
	  weights_[i] = value;
	}
	delete temp;
      }
      // create saved weights (not really indexedvector)
      savedWeights_ = new CoinIndexedVector();
      savedWeights_->reserve(numberRows);
      
      double * array = savedWeights_->denseVector();
      int * which = savedWeights_->getIndices();
      for (i=0;i<numberRows;i++) {
	array[i]=weights_[i];
	which[i]=pivotVariable[i];
      }
    } else {
      int * which = alternateWeights_->getIndices();
      CoinIndexedVector * rowArray3 = model_->rowArray(3);
      assert (!rowArray3->getNumElements());
      int * back = rowArray3->getIndices();
      // In case something went wrong
      for (i=0;i<numberRows+numberColumns;i++)
	back[i]=-1;
      if (mode!=4) {
	// save
	memcpy(savedWeights_->getIndices(),which,
	       numberRows*sizeof(int));
	memcpy(savedWeights_->denseVector(),weights_,
	       numberRows*sizeof(double));
      } else {
	// restore
	//memcpy(which,savedWeights_->getIndices(),
	//     numberRows*sizeof(int));
	//memcpy(weights_,savedWeights_->denseVector(),
	//     numberRows*sizeof(double));
	which = savedWeights_->getIndices();
      }
      // restore (a bit slow - but only every re-factorization)
      double * array = savedWeights_->denseVector();
      for (i=0;i<numberRows;i++) {
	int iSeq=which[i];
	back[iSeq]=i;
      }
      for (i=0;i<numberRows;i++) {
	int iPivot=pivotVariable[i];
	iPivot = back[iPivot];
	if (iPivot>=0) {
	  weights_[i]=array[iPivot];
	  if (weights_[i]<DEVEX_TRY_NORM)
	    weights_[i] = DEVEX_TRY_NORM; // may need to check more
	} else {
	  // odd
	  weights_[i]=1.0;
	}
      }
    }
    state_=0;
    // set up infeasibilities
    if (!infeasible_) {
      infeasible_ = new CoinIndexedVector();
      infeasible_->reserve(numberRows);
    }
  }
  if (mode>=2) {
    // Get dubious weights
    //if (!dubiousWeights_)
    //dubiousWeights_=new int[numberRows];
    //model_->factorization()->getWeights(dubiousWeights_);
    infeasible_->clear();
    int iRow;
    const int * pivotVariable = model_->pivotVariable();
    double tolerance=model_->currentPrimalTolerance();
    for (iRow=0;iRow<numberRows;iRow++) {
      int iPivot=pivotVariable[iRow];
      double value = model_->solution(iPivot);
      double lower = model_->lower(iPivot);
      double upper = model_->upper(iPivot);
      if (value<lower-tolerance) {
	value -= lower;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	infeasible_->quickAdd(iRow,value);
      } else if (value>upper+tolerance) {
	value -= upper;
	value *= value;
#ifdef COLUMN_BIAS 
	if (iPivot<numberColumns)
	  value *= COLUMN_BIAS; // bias towards columns
#endif
#ifdef FIXED_BIAS 
	if (lower==upper)
	  value *= FIXED_BIAS; // bias towards taking out fixed variables
#endif
	// store square in list
	infeasible_->quickAdd(iRow,value);
      }
    }
  }
}
// Gets rid of last update
void 
ClpDualRowSteepest::unrollWeights()
{
  double * saved = alternateWeights_->denseVector();
  int number = alternateWeights_->getNumElements();
  int * which = alternateWeights_->getIndices();
  int i;
  if (alternateWeights_->packedMode()) {
    for (i=0;i<number;i++) {
      int iRow = which[i];
      weights_[iRow]=saved[i];
      saved[i]=0.0;
    }
  } else {
    for (i=0;i<number;i++) {
      int iRow = which[i];
      weights_[iRow]=saved[iRow];
      saved[iRow]=0.0;
    }
  }
  alternateWeights_->setNumElements(0);
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpDualRowPivot * ClpDualRowSteepest::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpDualRowSteepest(*this);
  } else {
    return new ClpDualRowSteepest();
  }
}
// Gets rid of all arrays
void 
ClpDualRowSteepest::clearArrays()
{
  if (persistence_==normal) {
    delete [] weights_;
    weights_=NULL;
    delete [] dubiousWeights_;
    dubiousWeights_=NULL;
    delete infeasible_;
    infeasible_ = NULL;
    delete alternateWeights_;
    alternateWeights_ = NULL;
    delete savedWeights_;
    savedWeights_ = NULL;
  }
  state_ =-1;
}
// Returns true if would not find any row
bool 
ClpDualRowSteepest::looksOptimal() const
{
  int iRow;
  const int * pivotVariable = model_->pivotVariable();
  double tolerance=model_->currentPrimalTolerance();
  // we can't really trust infeasibilities if there is primal error
  // this coding has to mimic coding in checkPrimalSolution
  double error = CoinMin(1.0e-3,model_->largestPrimalError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance +  error;
  // But cap
  tolerance = CoinMin(1000.0,tolerance);
  int numberRows = model_->numberRows();
  int numberInfeasible=0;
  for (iRow=0;iRow<numberRows;iRow++) {
    int iPivot=pivotVariable[iRow];
    double value = model_->solution(iPivot);
    double lower = model_->lower(iPivot);
    double upper = model_->upper(iPivot);
    if (value<lower-tolerance) {
      numberInfeasible++;
    } else if (value>upper+tolerance) {
      numberInfeasible++;
    }
  }
  return (numberInfeasible==0);
}

