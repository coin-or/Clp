// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "ClpSimplex.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpMessage.hpp"
#include "CoinHelperFunctions.hpp"
#include <stdio.h>
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::ClpPrimalColumnSteepest (int mode) 
  : ClpPrimalColumnPivot(),
    state_(-1),
    mode_(mode),
    weights_(NULL),
    infeasible_(NULL),
    alternateWeights_(NULL),
    savedWeights_(NULL),
    pivotSequence_(-1),
    savedPivotSequence_(-1),
    savedSequenceOut_(-1),  
    reference_(NULL),
    devex_(0.0)
{
  type_=2+64*mode;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::ClpPrimalColumnSteepest (const ClpPrimalColumnSteepest & rhs) 
: ClpPrimalColumnPivot(rhs)
{  
  state_=rhs.state_;
  mode_ = rhs.mode_;
  model_ = rhs.model_;
  pivotSequence_ = rhs.pivotSequence_;
  savedPivotSequence_ = rhs.savedPivotSequence_;
  savedSequenceOut_ = rhs.savedSequenceOut_;
  devex_ = rhs.devex_;
  if (rhs.infeasible_) {
    infeasible_= new CoinIndexedVector(rhs.infeasible_);
  } else {
    infeasible_=NULL;
  }
  reference_=NULL;
  if (rhs.weights_) {
    assert(model_);
    int number = model_->numberRows()+model_->numberColumns();
    weights_= new double[number];
    CoinDisjointCopyN(rhs.weights_,number,weights_);
    savedWeights_= new double[number];
    CoinDisjointCopyN(rhs.savedWeights_,number,savedWeights_);
    if (!mode_) {
      reference_ = new unsigned int[(number+31)>>5];
      memcpy(reference_,rhs.reference_,((number+31)>>5)*sizeof(unsigned int));
    }
  } else {
    weights_=NULL;
    savedWeights_=NULL;
  }
  if (rhs.alternateWeights_) {
    alternateWeights_= new CoinIndexedVector(rhs.alternateWeights_);
  } else {
    alternateWeights_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPrimalColumnSteepest::~ClpPrimalColumnSteepest ()
{
  delete [] weights_;
  delete infeasible_;
  delete alternateWeights_;
  delete [] savedWeights_;
  delete [] reference_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPrimalColumnSteepest &
ClpPrimalColumnSteepest::operator=(const ClpPrimalColumnSteepest& rhs)
{
  if (this != &rhs) {
    ClpPrimalColumnPivot::operator=(rhs);
    state_=rhs.state_;
    mode_ = rhs.mode_;
    model_ = rhs.model_;
    pivotSequence_ = rhs.pivotSequence_;
    savedPivotSequence_ = rhs.savedPivotSequence_;
    savedSequenceOut_ = rhs.savedSequenceOut_;
    devex_ = rhs.devex_;
    delete [] weights_;
    delete [] reference_;
    reference_=NULL;
    delete infeasible_;
    delete alternateWeights_;
    delete [] savedWeights_;
    savedWeights_ = NULL;
    if (rhs.infeasible_!=NULL) {
      infeasible_= new CoinIndexedVector(rhs.infeasible_);
    } else {
      infeasible_=NULL;
    }
    if (rhs.weights_!=NULL) {
      assert(model_);
      int number = model_->numberRows()+model_->numberColumns();
      weights_= new double[number];
      CoinDisjointCopyN(rhs.weights_,number,weights_);
      savedWeights_= new double[number];
      CoinDisjointCopyN(rhs.savedWeights_,number,savedWeights_);
      if (!mode_) {
	reference_ = new unsigned int[(number+31)>>5];
	memcpy(reference_,rhs.reference_,
	       ((number+31)>>5)*sizeof(unsigned int));
      }
    } else {
      weights_=NULL;
    }
    if (rhs.alternateWeights_!=NULL) {
      alternateWeights_= new CoinIndexedVector(rhs.alternateWeights_);
    } else {
      alternateWeights_=NULL;
    }
  }
  return *this;
}

#define TRY_NORM 1.0e-4
#define ADD_ONE 1.0
// Returns pivot column, -1 if none
int 
ClpPrimalColumnSteepest::pivotColumn(CoinIndexedVector * updates,
				    CoinIndexedVector * spareRow1,
				    CoinIndexedVector * spareRow2,
				    CoinIndexedVector * spareColumn1,
				    CoinIndexedVector * spareColumn2)
{
  assert(model_);
  int iSection,j;
  int number;
  int * index;
  double * updateBy;
  double * reducedCost;
  // dj could be very small (or even zero - take care)
  double dj = model_->dualIn();
  double * infeas = infeasible_->denseVector();
  double tolerance=model_->currentDualTolerance();
  int pivotRow = model_->pivotRow();

  int anyUpdates;

  if (updates->getNumElements()) {
    // would have to have two goes for devex, three for steepest
    anyUpdates=2;
    // add in pivot contribution
    if (pivotRow>=0) 
      updates->add(pivotRow,-dj);
  } else if (pivotRow>=0) {
    if (fabs(dj)>1.0e-15) {
      // some dj
      updates->insert(pivotRow,-dj);
      if (fabs(dj)>1.0e-6) {
	// reasonable size
	anyUpdates=1;
      } else {
	// too small
	anyUpdates=2;
      }
    } else {
      // zero dj
      anyUpdates=-1;
    }
  } else if (pivotSequence_>=0){
    // just after re-factorization
    anyUpdates=-1;
  } else {
    // sub flip - nothing to do
    anyUpdates=0;
  }

  if (anyUpdates>0) {
    model_->factorization()->updateColumnTranspose(spareRow2,updates);
    // put row of tableau in rowArray and columnArray
    model_->clpMatrix()->transposeTimes(model_,-1.0,
					updates,spareColumn2,spareColumn1);
    for (iSection=0;iSection<2;iSection++) {
      
      reducedCost=model_->djRegion(iSection);
      int addSequence;
      
      if (!iSection) {
	number = updates->getNumElements();
	index = updates->getIndices();
	updateBy = updates->denseVector();
	addSequence = model_->numberColumns();
      } else {
	number = spareColumn1->getNumElements();
	index = spareColumn1->getIndices();
	updateBy = spareColumn1->denseVector();
	addSequence = 0;
      }

      for (j=0;j<number;j++) {
	int iSequence = index[j];
	double value = reducedCost[iSequence];
	value -= updateBy[iSequence];
	reducedCost[iSequence] = value;
	if (model_->fixed(iSequence+addSequence))
	  continue;
	ClpSimplex::Status status = model_->getStatus(iSequence+addSequence);

	switch(status) {
	  
	case ClpSimplex::basic:
	  infeasible_->zero(iSequence+addSequence);
	  break;
	case ClpSimplex::isFree:
	case ClpSimplex::superBasic:
	  if (fabs(value)>tolerance) {
	    // we are going to bias towards free
	    value *= 10.0;
	    // store square in list
	    if (infeas[iSequence+addSequence])
	      infeas[iSequence+addSequence] = value*value; // already there
	    else
	      infeasible_->quickAdd(iSequence+addSequence,value*value);
	  } else {
	    infeasible_->zero(iSequence+addSequence);
	  }
	  break;
	case ClpSimplex::atUpperBound:
	  if (value>tolerance) {
	    // store square in list
	    if (infeas[iSequence+addSequence])
	      infeas[iSequence+addSequence] = value*value; // already there
	    else
	      infeasible_->quickAdd(iSequence+addSequence,value*value);
	  } else {
	    infeasible_->zero(iSequence+addSequence);
	  }
	  break;
	case ClpSimplex::atLowerBound:
	  if (value<-tolerance) {
	    // store square in list
	    if (infeas[iSequence+addSequence])
	      infeas[iSequence+addSequence] = value*value; // already there
	    else
	      infeasible_->quickAdd(iSequence+addSequence,value*value);
	  } else {
	    infeasible_->zero(iSequence+addSequence);
	  }
	}
      }
    }
    if (anyUpdates==2) {
      // we can zero out as will have to get pivot row
      updates->clear();
      spareColumn1->clear();
    }
    if (pivotRow>=0) {
      // make sure infeasibility on incoming is 0.0
      int sequenceIn = model_->sequenceIn();
      infeasible_->zero(sequenceIn);
    }
  }
  // make sure outgoing from last iteration okay
  int sequenceOut = model_->sequenceOut();
  if (sequenceOut>=0) {
    if (!model_->fixed(sequenceOut)) {
      ClpSimplex::Status status = model_->getStatus(sequenceOut);
      double value = model_->reducedCost(sequenceOut);
      
      switch(status) {

      case ClpSimplex::basic:
	break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
	if (fabs(value)>tolerance) { 
	  // we are going to bias towards free
	  value *= 10.0;
	  // store square in list
	  if (infeas[sequenceOut])
	    infeas[sequenceOut] = value*value; // already there
	  else
	    infeasible_->quickAdd(sequenceOut,value*value);
	} else {
	  infeasible_->zero(sequenceOut);
	}
	break;
      case ClpSimplex::atUpperBound:
	if (value>tolerance) {
	  // store square in list
	  if (infeas[sequenceOut])
	    infeas[sequenceOut] = value*value; // already there
	  else
	    infeasible_->quickAdd(sequenceOut,value*value);
	} else {
	  infeasible_->zero(sequenceOut);
	}
	break;
      case ClpSimplex::atLowerBound:
	if (value<-tolerance) {
	  // store square in list
	  if (infeas[sequenceOut])
	    infeas[sequenceOut] = value*value; // already there
	  else
	    infeasible_->quickAdd(sequenceOut,value*value);
	} else {
	  infeasible_->zero(sequenceOut);
	}
      }
    }
  }

#ifdef RECOMPUTE_ALL_DJS
  for (iSection=0;iSection<2;iSection++) {
      
    reducedCost=model_->djRegion(iSection);
    int addSequence;
    int iSequence;
      
    if (!iSection) {
      number = model_->numberRows();
      addSequence = model_->numberColumns();
    } else {
      number = model_->numberColumns();
      addSequence = 0;
    }

    for (iSequence=0;iSequence<number;iSequence++) {
      double value = reducedCost[iSequence];
      if (model_->fixed(iSequence+addSequence))
	continue;
      ClpSimplex::Status status = model_->getStatus(iSequence+addSequence);
      
      switch(status) {

      case ClpSimplex::basic:
	infeasible_->zero(iSequence+addSequence);
	break;
      case ClpSimplex::isFree:
	if (fabs(value)>tolerance) {
	  // we are going to bias towards free
	  value *= 10.0;
	  // store square in list
	  if (infeas[iSequence+addSequence])
	    infeas[iSequence+addSequence] = value*value; // already there
	  else
	    infeasible_->quickAdd(iSequence+addSequence,value*value);
	  } else {
	    infeasible_->zero(iSequence+addSequence);
	  }
	break;
      case ClpSimplex::atUpperBound:
	if (value>tolerance) {
	    // store square in list
	  if (infeas[iSequence+addSequence])
	    infeas[iSequence+addSequence] = value*value; // already there
	  else
	    infeasible_->quickAdd(iSequence+addSequence,value*value);
	} else {
	  infeasible_->zero(iSequence+addSequence);
	}
	break;
      case ClpSimplex::atLowerBound:
	if (value<-tolerance) {
	  // store square in list
	  if (infeas[iSequence+addSequence])
	    infeas[iSequence+addSequence] = value*value; // already there
	  else
	    infeasible_->quickAdd(iSequence+addSequence,value*value);
	} else {
	  infeasible_->zero(iSequence+addSequence);
	}
      }
    }
  }
#endif
  // for weights update we use pivotSequence
  pivotRow = pivotSequence_;
  // unset in case sub flip
  pivotSequence_=-1;
  if (pivotRow>=0) {
    // make sure infeasibility on incoming is 0.0
    const int * pivotVariable = model_->pivotVariable();
    int sequenceIn = pivotVariable[pivotRow];
    infeasible_->zero(sequenceIn);
    // and we can see if reference
    double referenceIn=0.0;
    if (!mode_&&reference(sequenceIn))
      referenceIn=1.0;
    // save outgoing weight round update
    double outgoingWeight=weights_[sequenceOut];
    // update weights
    if (anyUpdates!=1) {
      updates->setNumElements(0);
      spareColumn1->setNumElements(0);
      // might as well set dj to 1
      dj = 1.0;
      updates->insert(pivotRow,-dj);
      model_->factorization()->updateColumnTranspose(spareRow2,updates);
      // put row of tableau in rowArray and columnArray
      model_->clpMatrix()->transposeTimes(model_,-1.0,
					  updates,spareColumn2,spareColumn1);
    }
    // now update weight update array
    model_->factorization()->updateColumnTranspose(spareRow2,
						   alternateWeights_);
    double * other = alternateWeights_->denseVector();
    double * weight;
    int numberColumns = model_->numberColumns();
    double scaleFactor = -1.0/dj; // as formula is with 1.0
    // rows
    number = updates->getNumElements();
    index = updates->getIndices();
    updateBy = updates->denseVector();
    weight = weights_+numberColumns;
    
    for (j=0;j<number;j++) {
      int iSequence = index[j];
      double thisWeight = weight[iSequence];
      // row has -1 
      double pivot = updateBy[iSequence]*scaleFactor;
      updateBy[iSequence]=0.0;
      double modification = other[iSequence];
      double pivotSquared = pivot * pivot;
      
      thisWeight += pivotSquared * devex_ + pivot * modification;
      if (thisWeight<TRY_NORM) {
	if (mode_) {
	  // steepest
	  thisWeight = max(TRY_NORM,ADD_ONE+pivotSquared);
	} else {
	  // exact
	  thisWeight = referenceIn*pivotSquared;
	  if (reference(iSequence+numberColumns))
	    thisWeight += 1.0;
	  thisWeight = max(thisWeight,TRY_NORM);
	}
      }
      weight[iSequence] = thisWeight;
    }
    
    // columns
    weight = weights_;
    
    scaleFactor = -scaleFactor;
    
    number = spareColumn1->getNumElements();
    index = spareColumn1->getIndices();
    updateBy = spareColumn1->denseVector();
    // get subset which have nonzero tableau elements
    model_->clpMatrix()->subsetTransposeTimes(model_,alternateWeights_,
					      spareColumn1,
					      spareColumn2);
    double * updateBy2 = spareColumn2->denseVector();
    for (j=0;j<number;j++) {
      int iSequence = index[j];
      double thisWeight = weight[iSequence];
      double pivot = updateBy[iSequence]*scaleFactor;
      updateBy[iSequence]=0.0;
      double modification = updateBy2[iSequence];
      updateBy2[iSequence]=0.0;
      double pivotSquared = pivot * pivot;
      
      thisWeight += pivotSquared * devex_ + pivot * modification;
      if (thisWeight<TRY_NORM) {
	if (mode_) {
	  // steepest
	  thisWeight = max(TRY_NORM,ADD_ONE+pivotSquared);
	} else {
	  // exact
	  thisWeight = referenceIn*pivotSquared;
	  if (reference(iSequence))
	    thisWeight += 1.0;
	  thisWeight = max(thisWeight,TRY_NORM);
	}
      }
      weight[iSequence] = thisWeight;
    }
    // restore outgoing weight
    weights_[sequenceOut]=outgoingWeight;
    alternateWeights_->clear();
    spareColumn2->setNumElements(0);
#ifdef SOME_DEBUG_1
#if 1
    // check for accuracy
    int iCheck=229;
    printf("weight for iCheck is %g\n",weights_[iCheck]);
    if (model_->getStatus(iCheck)!=ClpSimplex::basic&&
	!model_->fixed(iCheck))
      checkAccuracy(iCheck,1.0e-1,updates,spareRow2);
#else
    // check for accuracy
    number = updates->getNumElements();
    index = updates->getIndices();
    for (j=0;j<number;j++) {
      if (model_->getStatus(index[j])!=ClpSimplex::basic&&
	  !model_->fixed(index[j]))
	checkAccuracy(index[j],1.0e-1,updates,spareRow2);
    }
#endif
#endif
    updates->setNumElements(0);
    spareColumn1->setNumElements(0);
  }
  infeasible_->stopQuickAdd();


  // update of duals finished - now do pricing


  double bestDj = 1.0e-30;
  int bestSequence=-1;

  int i,iSequence;
  index = infeasible_->getIndices();
  number = infeasible_->getNumElements();
  if(model_->numberIterations()<model_->lastBadIteration()+200) {
    // we can't really trust infeasibilities if there is dual error
    double checkTolerance = 1.0e-8;
    if (!model_->factorization()->pivots())
      checkTolerance = 1.0e-6;
    if (model_->largestDualError()>checkTolerance)
      tolerance *= model_->largestDualError()/checkTolerance;
  }
  tolerance *= tolerance; // as we are using squares
  for (i=0;i<number;i++) {
    iSequence = index[i];
    double value = infeas[iSequence];
    double weight = weights_[iSequence];
    /*if (model_->numberIterations()%100==0)
      printf("%d inf %g wt %g\n",iSequence,value,weight);*/
    //weight=1.0;
    if (value>bestDj*weight&&value>tolerance) {
      // check flagged variable
      if (!model_->flagged(iSequence)) {
	/*if (model_->numberIterations()%100==0)
	  printf("chosen %g => %g\n",bestDj,value/weight);*/
	bestDj=value/weight;
	bestSequence = iSequence;
      }
    }
  }
  /*if (model_->numberIterations()%100==0)
    printf("%d best %g\n",bestSequence,bestDj);*/
  
#ifdef CLP_DEBUG
  if (bestSequence>=0) {
    if (model_->getStatus(bestSequence)==ClpSimplex::atLowerBound)
      assert(model_->reducedCost(bestSequence)<0.0);
    if (model_->getStatus(bestSequence)==ClpSimplex::atUpperBound)
      assert(model_->reducedCost(bestSequence)>0.0);
  }
#endif
  return bestSequence;
}
/* 
   1) before factorization
   2) after factorization
   3) just redo infeasibilities
   4) restore weights
*/
void 
ClpPrimalColumnSteepest::saveWeights(ClpSimplex * model,int mode)
{
  // alternateWeights_ is defined as indexed but is treated oddly
  // at times
  model_ = model;
  int numberRows = model_->numberRows();
  int numberColumns = model_->numberColumns();
  const int * pivotVariable = model_->pivotVariable();
  if (mode==1&&weights_) {
    if (pivotSequence_>=0) {
      // save pivot order
      memcpy(alternateWeights_->getIndices(),pivotVariable,
	     numberRows*sizeof(int));
      // change from pivot row number to sequence number
      pivotSequence_=pivotVariable[pivotSequence_];
    } else {
      alternateWeights_->clear();
    }
    state_=1;
  } else if (mode==2||mode==4||mode==5) {
    // restore
    if (!weights_||state_==-1||mode==5) {
      // initialize weights
      delete [] weights_;
      delete alternateWeights_;
      weights_ = new double[numberRows+numberColumns];
      alternateWeights_ = new CoinIndexedVector();
      // enough space so can use it for factorization
      alternateWeights_->reserve(numberRows+
				 model_->factorization()->maximumPivots());
      initializeWeights();
      // create saved weights 
      savedWeights_ = new double[numberRows+numberColumns];
      memcpy(savedWeights_,weights_,(numberRows+numberColumns)*
	     sizeof(double));
      savedPivotSequence_=-2;
      savedSequenceOut_ = -2;
      
    } else {
      if (mode!=4) {
	// save
	memcpy(savedWeights_,weights_,(numberRows+numberColumns)*
	       sizeof(double));
	savedPivotSequence_= pivotSequence_;
	savedSequenceOut_ = model_->sequenceOut();
      } else {
	// restore
	memcpy(weights_,savedWeights_,(numberRows+numberColumns)*
	       sizeof(double));
	pivotSequence_= savedPivotSequence_;
	model_->setSequenceOut(savedSequenceOut_); 
      }
    }
    state_=0;
    // set up infeasibilities
    if (!infeasible_) {
      infeasible_ = new CoinIndexedVector();
      infeasible_->reserve(numberColumns+numberRows);
    }
  }
  if (mode>=2&&mode!=5) {
    if (mode!=3&&pivotSequence_>=0) {
      // restore pivot row
      int iRow;
      // permute alternateWeights
      double * temp = new double[numberRows+numberColumns];
      double * work = alternateWeights_->denseVector();
      int * oldPivotOrder = alternateWeights_->getIndices();
      for (iRow=0;iRow<numberRows;iRow++) {
	int iPivot=oldPivotOrder[iRow];
	temp[iPivot]=work[iRow];
      }
      int number=0;
      int found=-1;
      int * which = oldPivotOrder;
      // find pivot row and re-create alternateWeights
      for (iRow=0;iRow<numberRows;iRow++) {
	int iPivot=pivotVariable[iRow];
	if (iPivot==pivotSequence_) 
	  found=iRow;
	work[iRow]=temp[iPivot];
	if (work[iRow])
	  which[number++]=iRow;
      }
      alternateWeights_->setNumElements(number);
#ifdef CLP_DEBUG
      // Can happen but I should clean up
      assert(found>=0);
#endif
      pivotSequence_ = found;
      delete [] temp;
    }
    infeasible_->clear();
    double tolerance=model_->currentDualTolerance();
    int number = model_->numberRows() + model_->numberColumns();
    int iSequence;

    double * reducedCost = model_->djRegion();
      
    for (iSequence=0;iSequence<number;iSequence++) {
      if (model_->fixed(iSequence))
	continue;
      double value = reducedCost[iSequence];
      ClpSimplex::Status status = model_->getStatus(iSequence);
      
      switch(status) {

      case ClpSimplex::basic:
	break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
	if (fabs(value)>tolerance) { 
	  // store square in list
	  infeasible_->quickAdd(iSequence,value*value);
	}
	break;
      case ClpSimplex::atUpperBound:
	if (value>tolerance) {
	  infeasible_->quickAdd(iSequence,value*value);
	}
	break;
      case ClpSimplex::atLowerBound:
	if (value<-tolerance) {
	  infeasible_->quickAdd(iSequence,value*value);
	}
      }
    }
    infeasible_->stopQuickAdd();
  }
}
// Gets rid of last update
void 
ClpPrimalColumnSteepest::unrollWeights()
{
  double * saved = alternateWeights_->denseVector();
  int number = alternateWeights_->getNumElements();
  int * which = alternateWeights_->getIndices();
  int i;
  for (i=0;i<number;i++) {
    int iRow = which[i];
    weights_[iRow]=saved[iRow];
    saved[iRow]=0.0;
  }
  alternateWeights_->setNumElements(0);
}
  
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpPrimalColumnPivot * ClpPrimalColumnSteepest::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpPrimalColumnSteepest(*this);
  } else {
    return new ClpPrimalColumnSteepest();
  }
}
void
ClpPrimalColumnSteepest::updateWeights(CoinIndexedVector * input)
{
  int number=input->getNumElements();
  int * which = input->getIndices();
  double * work = input->denseVector();
  int newNumber = 0;
  int * newWhich = alternateWeights_->getIndices();
  double * newWork = alternateWeights_->denseVector();
  int i;
  int sequenceIn = model_->sequenceIn();
  int sequenceOut = model_->sequenceOut();
  const int * pivotVariable = model_->pivotVariable();

  int pivotRow = model_->pivotRow();
  pivotSequence_ = pivotRow;

  devex_ =0.0;

  if (pivotRow>=0) {
    if (mode_) {
      for (i=0;i<number;i++) {
	int iRow = which[i];
	devex_ += work[iRow]*work[iRow];
	newWork[iRow] = -2.0*work[iRow];
      }
      newWork[pivotRow] = -2.0*max(devex_,0.0);
      devex_+=ADD_ONE;
      weights_[sequenceOut]=1.0+ADD_ONE;
      memcpy(newWhich,which,number*sizeof(int));
      alternateWeights_->setNumElements(number);
    } else {
      for (i=0;i<number;i++) {
	int iRow = which[i];
	int iPivot = pivotVariable[iRow];
	if (reference(iPivot)) {
	  devex_ += work[iRow]*work[iRow];
	  newWork[iRow] = -2.0*work[iRow];
	  newWhich[newNumber++]=iRow;
	}
      }
      if (!newWork[pivotRow])
	newWhich[newNumber++]=pivotRow; // add if not already in
      newWork[pivotRow] = -2.0*max(devex_,0.0);
      if (reference(sequenceIn)) {
	devex_+=1.0;
      } else {
      }
      if (reference(sequenceOut)) {
	weights_[sequenceOut]=1.0+1.0;
      } else {
	weights_[sequenceOut]=1.0;
      }
      alternateWeights_->setNumElements(newNumber);
    }
  } else {
    if (mode_) {
      for (i=0;i<number;i++) {
	int iRow = which[i];
	devex_ += work[iRow]*work[iRow];
      }
      devex_ += ADD_ONE;
    } else {
      for (i=0;i<number;i++) {
	int iRow = which[i];
	int iPivot = pivotVariable[iRow];
	if (reference(iPivot)) {
	  devex_ += work[iRow]*work[iRow];
	}
      }
      if (reference(sequenceIn)) 
	devex_+=1.0;
    }
  }
  double oldDevex = weights_[sequenceIn];
#ifdef CLP_DEBUG
  if ((model_->messageHandler()->logLevel()&32))
    printf("old weight %g, new %g\n",oldDevex,devex_);
#endif
  double check = max(devex_,oldDevex);;
  weights_[sequenceIn] = devex_;
  if ( fabs ( devex_ - oldDevex ) > 1.0e-1 * check ) {
#ifdef CLP_DEBUG
    if ((model_->messageHandler()->logLevel()&48)==16)
      printf("old weight %g, new %g\n",oldDevex,devex_);
#endif
    if ( fabs ( devex_ - oldDevex ) > 1.0e5 * check ) {
      // need to redo
      model_->messageHandler()->message(CLP_INITIALIZE_STEEP,
					*model_->messagesPointer())
					  <<oldDevex<<devex_
					  <<CoinMessageEol;
      initializeWeights();
    }
  }
  if (pivotRow>=0) {
    // set outgoing weight here
    weights_[model_->sequenceOut()]=devex_/(model_->alpha()*model_->alpha());
  }
}
// Checks accuracy - just for debug
void 
ClpPrimalColumnSteepest::checkAccuracy(int sequence,
				       double relativeTolerance,
				       CoinIndexedVector * rowArray1,
				       CoinIndexedVector * rowArray2)
{
  model_->unpack(rowArray1,sequence);
  model_->factorization()->updateColumn(rowArray2,rowArray1);
  int number=rowArray1->getNumElements();
  int * which = rowArray1->getIndices();
  double * work = rowArray1->denseVector();
  const int * pivotVariable = model_->pivotVariable();
  
  double devex =0.0;
  int i;

  if (mode_) {
    for (i=0;i<number;i++) {
      int iRow = which[i];
      devex += work[iRow]*work[iRow];
      work[iRow]=0.0;
    }
    devex += ADD_ONE;
  } else {
    for (i=0;i<number;i++) {
      int iRow = which[i];
      int iPivot = pivotVariable[iRow];
      if (reference(iPivot)) {
	devex += work[iRow]*work[iRow];
      }
      work[iRow]=0.0;
    }
    if (reference(sequence)) 
      devex+=1.0;
  }

  double oldDevex = weights_[sequence];
  double check = max(devex,oldDevex);;
  if ( fabs ( devex - oldDevex ) > relativeTolerance * check ) {
    printf("check %d old weight %g, new %g\n",sequence,oldDevex,devex);
    // update so won't print again
    weights_[sequence]=devex;
  }
  rowArray1->setNumElements(0);
}

// Initialize weights
void 
ClpPrimalColumnSteepest::initializeWeights()
{
  int numberRows = model_->numberRows();
  int numberColumns = model_->numberColumns();
  int number = numberRows + numberColumns;
  int iSequence;
  if (!mode_) {
    // initialize to 1.0 
    // and set reference framework
    if (!reference_) {
      int nWords = (number+31)>>5;
      reference_ = new unsigned int[nWords];
      assert (sizeof(unsigned int)==4);
      // tiny overhead to zero out (stops valgrind complaining)
      memset(reference_,0,nWords*sizeof(int));
    }
    
    for (iSequence=0;iSequence<number;iSequence++) {
      weights_[iSequence]=1.0;
      if (model_->getStatus(iSequence)==ClpSimplex::basic) {
	setReference(iSequence,false);
      } else {
	setReference(iSequence,true);
      }
    }
  } else {
    CoinIndexedVector * temp = new CoinIndexedVector();
    temp->reserve(numberRows+
		  model_->factorization()->maximumPivots());
    double * array = alternateWeights_->denseVector();
    int * which = alternateWeights_->getIndices();
      
    for (iSequence=0;iSequence<number;iSequence++) {
      weights_[iSequence]=1.0+ADD_ONE;
      if (model_->fixed(iSequence))
	continue;
      if (model_->getStatus(iSequence)!=ClpSimplex::basic) {
	model_->unpack(alternateWeights_,iSequence);
	double value=ADD_ONE;
	model_->factorization()->updateColumn(temp,alternateWeights_);
	int number = alternateWeights_->getNumElements();
	int j;
	for (j=0;j<number;j++) {
	  int iRow=which[j];
	  value+=array[iRow]*array[iRow];
	  array[iRow]=0.0;
	}
	alternateWeights_->setNumElements(0);
	weights_[iSequence] = value;
      }
    }
    delete temp;
  }
}
// Gets rid of all arrays
void 
ClpPrimalColumnSteepest::clearArrays()
{
  delete [] weights_;
  weights_=NULL;
  delete infeasible_;
  infeasible_ = NULL;
  delete alternateWeights_;
  alternateWeights_ = NULL;
  delete [] savedWeights_;
  savedWeights_ = NULL;
  delete [] reference_;
  reference_ = NULL;
}
