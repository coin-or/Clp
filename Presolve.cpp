// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

//#define	CHECK_CONSISTENCY	1

#include <stdio.h>

#include <assert.h>
#include <iostream>

#include "CoinHelperFunctions.hpp"

#include "CoinPackedMatrix.hpp"
#include "ClpSimplex.hpp"

#include "Presolve.hpp"
#include "PresolveMatrix.hpp"

#include "PresolveEmpty.hpp"
#include "PresolveFixed.hpp"
#include "PresolvePsdebug.hpp"
#include "PresolveSingleton.hpp"
#include "PresolveDoubleton.hpp"
#include "PresolveZeros.hpp"
#include "PresolveSubst.hpp"
#include "PresolveForcing.hpp"
#include "PresolveDual.hpp"
#include "PresolveTighten.hpp"
#include "PresolveUseless.hpp"
#include "PresolveDupcol.hpp"
#include "PresolveImpliedFree.hpp"
#include "PresolveIsolated.hpp"
#include "ClpMessage.hpp"



Presolve::Presolve() :
  originalModel_(NULL),
  presolvedModel_(NULL),
  nonLinearValue_(0.0),
  originalColumn_(NULL),
  paction_(0),
  ncols_(0),
  nrows_(0),
  nelems_(0),
  numberPasses_(5)
{
}

Presolve::~Presolve()
{
  gutsOfDestroy();
}
// Gets rid of presolve actions (e.g.when infeasible)
void 
Presolve::gutsOfDestroy()
{
 const PresolveAction *paction = paction_;
  while (paction) {
    const PresolveAction *next = paction->next;
    delete paction;
    paction = next;
  }
  delete [] originalColumn_;
  paction_=NULL;
  originalColumn_=NULL;
}

#if 0
void Presolve::presolve(ClpSimplex& si)
{
  ncols_ = si.getNumCols();
  nrows_ = si.getNumRows();
  nelems_ = si.getNumElements();

  double maxmin = si.getObjSense();

  delete presolvedModel_;
  presolvedModel_ = NULL;
  originalModel_ = NULL;

  PresolveMatrix prob(ncols_,
			 maxmin,
			 si,
			 nrows_, nelems_,false);
  prob.originalModel_ = originalModel_;

  paction_ = presolve(&prob);

  prob.update_model(si, nrows_, ncols_, nelems_);
  int ncolsNow = si.getNumCols();
  originalColumn_ = new int[ncolsNow];
  memcpy(originalColumn_,prob.originalColumn_,ncolsNow*sizeof(int));
}
  

// This is the toplevel postsolving routine.
void Presolve::postsolve(ClpSimplex& si,
			    unsigned char *colstat,
			    unsigned char *rowstat)
{
  if (!si.isProvenOptimal()) {
    throw CoinError("solution not optimal",
		    "postsolve",
		    "Presolve");
  }

  // this is the size of the original problem
  const int ncols0  = ncols_;
  const int nrows0  = nrows_;
  const CoinBigIndex nelems0 = nelems_;

  double *acts = new double[nrows0];
  double *sol  = new double[ncols0];

  {
    // put this in a block to make sure it gets deleted early
    PostsolveMatrix prob(si,
			    ncols0,
			    nrows0,
			    nelems0,
			    si.getObjSense(),
			    // end prepost

			    sol, acts,
			 colstat, rowstat);
    prob.originalModel_ = originalModel_;

    postsolve(prob);

    // Postsolve does not use a standard matrix format, since row
    // and col numbers change rapidly.
    // The col rep is just a linked list.
    int *hrow = new int[nelems0];
    int *mcstrt = new int[ncols0+1];
    double *colels = new double[nelems0];

    {
      const int *link = prob.link_;
      CoinBigIndex start = 0;
      
      for (int j=0; j<ncols0; j++) {
	CoinBigIndex k = prob.mcstrt_[j];
	int n = prob.hincol_[j];

	mcstrt[j] = start;
	for (int i=0; i<n; i++) {
	  hrow[start] = prob.hrow_[k];
	  colels[start] = prob.colels_[k];
	  k = link[k];
	  start++;
	}
      }
      mcstrt[ncols0] = start;
    }
    
    si.loadProblem(ncols0, nrows0,
		   mcstrt, hrow, colels,
		   
		   prob.clo_, prob.cup_,
		   prob.cost_,
		   prob.rlo_, prob.rup_);

    delete[]hrow;
    delete[]mcstrt;
    delete[]colels;
  }


  si.setColSolution(sol);
  delete[]acts;

  delete[]sol;
}
#endif
/* This version of presolve returns a pointer to a new presolved 
   model.  NULL if infeasible
*/
ClpSimplex * 
Presolve::presolvedModel(ClpSimplex & si,
			 double feasibilityTolerance,
			 bool keepIntegers,
			 int numberPasses)
{
  ncols_ = si.getNumCols();
  nrows_ = si.getNumRows();
  nelems_ = si.getNumElements();
  numberPasses_ = numberPasses;

  double maxmin = si.getObjSense();
  originalModel_ = &si;
  delete [] originalColumn_;
  originalColumn_ = new int[ncols_];

  // Check matrix
  if (!si.clpMatrix()->allElementsInRange(&si,1.0e-20,1.0e20))
    return NULL;

  // result is 0 - okay, 1 infeasible, -1 go round again
  int result = -1;
  
  // User may have deleted - its their responsibility
  presolvedModel_=NULL;

  while (result==-1) {

    // make new copy
    delete presolvedModel_;
    presolvedModel_ = new ClpSimplex(si);
    

    // drop integer information if wanted
    if (!keepIntegers)
      presolvedModel_->deleteIntegerInformation();

    
    PresolveMatrix prob(ncols_,
			maxmin,
			*presolvedModel_,
			nrows_, nelems_,true,nonLinearValue_);
    // make sure row solution correct
    {
      double *colels	= prob.colels_;
      int *hrow		= prob.hrow_;
      CoinBigIndex *mcstrt		= prob.mcstrt_;
      int *hincol		= prob.hincol_;
      int ncols		= prob.ncols_;
      
      
      double * csol = prob.sol_;
      double * acts = prob.acts_;
      int nrows = prob.nrows_;

      int colx;

      memset(acts,0,nrows*sizeof(double));
      
      for (colx = 0; colx < ncols; ++colx) {
	double solutionValue = csol[colx];
	for (int i=mcstrt[colx]; i<mcstrt[colx]+hincol[colx]; ++i) {
	  int row = hrow[i];
	  double coeff = colels[i];
	  acts[row] += solutionValue*coeff;
	}
      }
    }
    prob.handler_ = presolvedModel_->messageHandler();
    prob.messages_ = presolvedModel_->messages();

    // move across feasibility tolerance
    prob.feasibilityTolerance_ = feasibilityTolerance;

    // Do presolve

    paction_ = presolve(&prob);

    result =0; 

    if (prob.status_==0&&paction_) {
      // feasible
    
      prob.update_model(*presolvedModel_, nrows_, ncols_, nelems_);
      // copy status and solution
      memcpy(presolvedModel_->primalColumnSolution(),
	     prob.sol_,prob.ncols_*sizeof(double));
      memcpy(presolvedModel_->primalRowSolution(),
	     prob.acts_,prob.nrows_*sizeof(double));
      memcpy(presolvedModel_->statusArray(),
	     prob.colstat_,prob.ncols_*sizeof(unsigned char));
      memcpy(presolvedModel_->statusArray()+prob.ncols_,
	     prob.rowstat_,prob.nrows_*sizeof(unsigned char));
      delete [] prob.sol_;
      delete [] prob.acts_;
      delete [] prob.colstat_;
      prob.sol_=NULL;
      prob.acts_=NULL;
      prob.colstat_=NULL;
      
      int ncolsNow = presolvedModel_->getNumCols();
      memcpy(originalColumn_,prob.originalColumn_,ncolsNow*sizeof(int));
      // now clean up integer variables.  This can modify original
      int i;
      const char * information = presolvedModel_->integerInformation();
      if (information) {
	int numberChanges=0;
	double * lower0 = originalModel_->columnLower();
	double * upper0 = originalModel_->columnUpper();
	double * lower = presolvedModel_->columnLower();
	double * upper = presolvedModel_->columnUpper();
	for (i=0;i<ncolsNow;i++) {
	  if (!information[i])
	    continue;
	  int iOriginal = originalColumn_[i];
	  double lowerValue0 = lower0[iOriginal];
	  double upperValue0 = upper0[iOriginal];
	  double lowerValue = ceil(lower[i]-1.0e-5);
	  double upperValue = floor(upper[i]+1.0e-5);
	  lower[i]=lowerValue;
	  upper[i]=upperValue;
	  if (lowerValue>upperValue) {
	    numberChanges++;
	    presolvedModel_->messageHandler()->message(CLP_PRESOLVE_COLINFEAS,
						       presolvedModel_->messages())
							 <<iOriginal
							 <<lowerValue
							 <<upperValue
							 <<CoinMessageEol;
	    result=1;
	  } else {
	    if (lowerValue>lowerValue0+1.0e-8) {
	      lower0[iOriginal] = lowerValue;
	      numberChanges++;
	    }
	    if (upperValue<upperValue0-1.0e-8) {
	      upper0[iOriginal] = upperValue;
	      numberChanges++;
	    }
	  }	  
	}
	if (numberChanges) {
	  presolvedModel_->messageHandler()->message(CLP_PRESOLVE_INTEGERMODS,
						     presolvedModel_->messages())
						       <<numberChanges
						       <<CoinMessageEol;
	  if (!result) {
	    result = -1; // round again
	  }
	}
      }
    } else {
      // infeasible
      delete [] prob.sol_;
      delete [] prob.acts_;
      delete [] prob.colstat_;
      result=1;
    }
  }
  if (!result) {
    int nrowsAfter = presolvedModel_->getNumRows();
    int ncolsAfter = presolvedModel_->getNumCols();
    CoinBigIndex nelsAfter = presolvedModel_->getNumElements();
    presolvedModel_->messageHandler()->message(CLP_PRESOLVE_STATS,
					       presolvedModel_->messages())
						 <<nrowsAfter<< -(nrows_ - nrowsAfter)
						 << ncolsAfter<< -(ncols_ - ncolsAfter)
						 <<nelsAfter<< -(nelems_ - nelsAfter)
						 <<CoinMessageEol;
  } else {
    gutsOfDestroy();
    delete presolvedModel_;
    presolvedModel_=NULL;
  }
  return presolvedModel_;
}

// Return pointer to presolved model
ClpSimplex * 
Presolve::model() const
{
  return presolvedModel_;
}
// Return pointer to original model
ClpSimplex * 
Presolve::originalModel() const
{
  return originalModel_;
}
void 
Presolve::postsolve(bool updateStatus)
{
  if (!presolvedModel_->isProvenOptimal()) {
    presolvedModel_->messageHandler()->message(CLP_PRESOLVE_NONOPTIMAL,
					     presolvedModel_->messages())
					       <<CoinMessageEol;
  }

  // this is the size of the original problem
  const int ncols0  = ncols_;
  const int nrows0  = nrows_;
  const CoinBigIndex nelems0 = nelems_;

  // reality check
  assert(ncols0==originalModel_->getNumCols());
  assert(nrows0==originalModel_->getNumRows());

  // this is the reduced problem
  int ncols = presolvedModel_->getNumCols();
  int nrows = presolvedModel_->getNumRows();

  double *acts = originalModel_->primalRowSolution();
  double *sol  = originalModel_->primalColumnSolution();

  unsigned char * rowstat=NULL;
  unsigned char * colstat = NULL;
  if (updateStatus) {
    unsigned char *status = originalModel_->statusArray();
    rowstat = status + ncols0;
    colstat = status;
    memcpy(colstat, presolvedModel_->statusArray(), ncols);
    memcpy(rowstat, presolvedModel_->statusArray()+ncols, nrows);
  }


  PostsolveMatrix prob(*presolvedModel_,
		       ncols0,
		       nrows0,
		       nelems0,
		       presolvedModel_->getObjSense(),
		       // end prepost
		       
		       sol, acts,
		       colstat, rowstat);
    
  postsolve(prob);

}

// return pointer to original columns
const int * 
Presolve::originalColumns() const
{
  return originalColumn_;
}
// Set pointer to original model
void 
Presolve::setOriginalModel(ClpSimplex * model)
{
  originalModel_=model;
}
#if 0
// A lazy way to restrict which transformations are applied
// during debugging.
static int ATOI(const char *name)
{
 return true;
#if	DEBUG_PRESOLVE || PRESOLVE_SUMMARY
  if (getenv(name)) {
    int val = atoi(getenv(name));
    printf("%s = %d\n", name, val);
    return (val);
  } else {
    if (strcmp(name,"off"))
      return (true);
    else
      return (false);
  }
#else
  return (true);
#endif
}
#endif
//#define DEBUG_PRESOLVE 1
#if DEBUG_PRESOLVE
void check_sol(PresolveMatrix *prob,double tol)
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  int *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *hinrow		= prob->hinrow_;
  int ncols		= prob->ncols_;


  double * csol = prob->sol_;
  double * acts = prob->acts_;
  double * clo = prob->clo_;
  double * cup = prob->cup_;
  int nrows = prob->nrows_;
  double * rlo = prob->rlo_;
  double * rup = prob->rup_;

  int colx;

  double * rsol = new double[nrows];
  memset(rsol,0,nrows*sizeof(double));

  for (colx = 0; colx < ncols; ++colx) {
    if (1) {
      CoinBigIndex k = mcstrt[colx];
      int nx = hincol[colx];
      double solutionValue = csol[colx];
      for (int i=0; i<nx; ++i) {
	int row = hrow[k];
	double coeff = colels[k];
	k++;
	rsol[row] += solutionValue*coeff;
      }
      if (csol[colx]<clo[colx]-tol) {
	printf("low CSOL:  %d  - %g %g %g\n",
		   colx, clo[colx], csol[colx], cup[colx]);
      } else if (csol[colx]>cup[colx]+tol) {
	printf("high CSOL:  %d  - %g %g %g\n",
		   colx, clo[colx], csol[colx], cup[colx]);
      } 
    }
  }
  int rowx;
  for (rowx = 0; rowx < nrows; ++rowx) {
    if (hinrow[rowx]) {
      if (fabs(rsol[rowx]-acts[rowx])>tol)
	printf("inacc RSOL:  %d - %g %g (acts_ %g) %g\n",
		   rowx,  rlo[rowx], rsol[rowx], acts[rowx], rup[rowx]);
      if (rsol[rowx]<rlo[rowx]-tol) {
	printf("low RSOL:  %d - %g %g %g\n",
		   rowx,  rlo[rowx], rsol[rowx], rup[rowx]);
      } else if (rsol[rowx]>rup[rowx]+tol ) {
	printf("high RSOL:  %d - %g %g %g\n",
		   rowx,  rlo[rowx], rsol[rowx], rup[rowx]);
      } 
    }
  }
  delete [] rsol;
}
#endif
// This is the presolve loop.
// It is a separate virtual function so that it can be easily
// customized by subclassing Presolve.
const PresolveAction *Presolve::presolve(PresolveMatrix *prob)
{
  paction_ = 0;

  prob->status_=0; // say feasible

  paction_ = make_fixed(prob, paction_);
  // if integers then switch off dual stuff
  // later just do individually
  bool doDualStuff = (presolvedModel_->integerInformation()==NULL);

#if	CHECK_CONSISTENCY
  presolve_links_ok(prob->rlink_, prob->mrstrt_, prob->hinrow_, prob->nrows_);
#endif

  if (!prob->status_) {
#if 0
    const bool slackd = ATOI("SLACKD")!=0;
    //const bool forcing = ATOI("FORCING")!=0;
    const bool doubleton = ATOI("DOUBLETON")!=0;
    const bool forcing = ATOI("off")!=0;
    const bool ifree = ATOI("off")!=0;
    const bool zerocost = ATOI("off")!=0;
    const bool dupcol = ATOI("off")!=0;
    const bool duprow = ATOI("off")!=0;
    const bool dual = ATOI("off")!=0;
#else
    const bool slackd = true;
    const bool doubleton = true;
    const bool forcing = true;
    const bool ifree = true;
    const bool zerocost = true;
    const bool dupcol = true;
    const bool duprow = true;
    const bool dual = doDualStuff;
#endif
    
    // some things are expensive so just do once (normally)

    int i;
    // say look at all
    if (!prob->anyProhibited()) {
      for (i=0;i<nrows_;i++) 
	prob->rowsToDo_[i]=i;
      prob->numberRowsToDo_=nrows_;
      for (i=0;i<ncols_;i++) 
	prob->colsToDo_[i]=i;
      prob->numberColsToDo_=ncols_;
    } else {
      // some stuff must be left alone
      prob->numberRowsToDo_=0;
      for (i=0;i<nrows_;i++) 
	if (!prob->rowProhibited(i))
	    prob->rowsToDo_[prob->numberRowsToDo_++]=i;
      prob->numberColsToDo_=0;
      for (i=0;i<ncols_;i++) 
	if (!prob->colProhibited(i))
	    prob->colsToDo_[prob->numberColsToDo_++]=i;
    }
	

    int iLoop;
#if	DEBUG_PRESOLVE
    check_sol(prob,1.0e0);
#endif

    // Check number rows dropped
    int lastDropped=0;
    for (iLoop=0;iLoop<numberPasses_;iLoop++) {
#ifdef PRESOLVE_SUMMARY
      printf("Starting major pass %d\n",iLoop+1);
#endif
      const PresolveAction * const paction0 = paction_;
      // look for substitutions with no fill
      int fill_level=2;
      int whichPass=0;
      while (1) {
	whichPass++;
	const PresolveAction * const paction1 = paction_;
	
	if (slackd) {
	  bool notFinished = true;
	  while (notFinished) 
	    paction_ = slack_doubleton_action::presolve(prob, paction_,
							notFinished);
	  if (prob->status_)
	    break;
	}
	
	if (doubleton) {
	  paction_ = doubleton_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	}
	
	if (zerocost) {
	  paction_ = do_tighten_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	}
	
	if (forcing) {
	  paction_ = forcing_constraint_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	}
	
	if (ifree) {
	  paction_ = implied_free_action::presolve(prob, paction_,fill_level);
	  if (prob->status_)
	    break;
	}
	
#if	DEBUG_PRESOLVE
	check_sol(prob,1.0e0);
#endif

#if	CHECK_CONSISTENCY
	presolve_links_ok(prob->rlink_, prob->mrstrt_, prob->hinrow_, 
			  prob->nrows_);
#endif
	
#if	DEBUG_PRESOLVE
	presolve_no_zeros(prob->mcstrt_, prob->colels_, prob->hincol_, 
			  prob->ncols_);
#endif
#if	CHECK_CONSISTENCY
	prob->consistent();
#endif
	
	  
	// set up for next pass
	// later do faster if many changes i.e. memset and memcpy
	prob->numberRowsToDo_ = prob->numberNextRowsToDo_;
	int kcheck;
	bool found;
	kcheck=-1;
	found=false;
	for (i=0;i<prob->numberNextRowsToDo_;i++) {
	  int index = prob->nextRowsToDo_[i];
	  prob->unsetRowChanged(index);
	  prob->rowsToDo_[i] = index;
	  if (index==kcheck) {
	    printf("row %d on list after pass %d\n",kcheck,
		   whichPass);
	    found=true;
	  }
	}
	if (!found&&kcheck>=0)
	  prob->rowsToDo_[prob->numberRowsToDo_++]=kcheck;
	prob->numberNextRowsToDo_=0;
	prob->numberColsToDo_ = prob->numberNextColsToDo_;
	kcheck=-1;
	found=false;
	for (i=0;i<prob->numberNextColsToDo_;i++) {
	  int index = prob->nextColsToDo_[i];
	  prob->unsetColChanged(index);
	  prob->colsToDo_[i] = index;
	  if (index==kcheck) {
	    printf("col %d on list after pass %d\n",kcheck,
		   whichPass);
	    found=true;
	  }
	}
	if (!found&&kcheck>=0)
	  prob->colsToDo_[prob->numberColsToDo_++]=kcheck;
	prob->numberNextColsToDo_=0;
	if (paction_ == paction1&&fill_level>0)
	  break;
      }
      // say look at all
      int i;
      if (!prob->anyProhibited()) {
	for (i=0;i<nrows_;i++) 
	  prob->rowsToDo_[i]=i;
	prob->numberRowsToDo_=nrows_;
	for (i=0;i<ncols_;i++) 
	  prob->colsToDo_[i]=i;
	prob->numberColsToDo_=ncols_;
      } else {
	// some stuff must be left alone
	prob->numberRowsToDo_=0;
	for (i=0;i<nrows_;i++) 
	  if (!prob->rowProhibited(i))
	    prob->rowsToDo_[prob->numberRowsToDo_++]=i;
	prob->numberColsToDo_=0;
	for (i=0;i<ncols_;i++) 
	  if (!prob->colProhibited(i))
	    prob->colsToDo_[prob->numberColsToDo_++]=i;
      }
      // now expensive things
      // this caused world.mps to run into numerical difficulties
#ifdef PRESOLVE_SUMMARY
      printf("Starting expensive\n");
#endif

      if (dual) {
	int itry;
	for (itry=0;itry<5;itry++) {
	  const PresolveAction * const paction2 = paction_;
	  paction_ = remove_dual_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	  if (ifree) {
	    int fill_level=0; // switches off substitution
	    paction_ = implied_free_action::presolve(prob, paction_,fill_level);
	    if (prob->status_)
	      break;
	  }
	  if (paction_ == paction2)
	    break;
	}
      }
#if	DEBUG_PRESOLVE
      check_sol(prob,1.0e0);
#endif
      if (dupcol) {
	paction_ = dupcol_action::presolve(prob, paction_);
	if (prob->status_)
	  break;
      }
#if	DEBUG_PRESOLVE
	check_sol(prob,1.0e0);
#endif
      
      if (duprow) {
	paction_ = duprow_action::presolve(prob, paction_);
	if (prob->status_)
	  break;
      }
#if	DEBUG_PRESOLVE
      check_sol(prob,1.0e0);
#endif
      {
	int * hinrow = prob->hinrow_;
	int numberDropped=0;
	for (i=0;i<nrows_;i++) 
	  if (!hinrow[i])
	    numberDropped++;
	//printf("%d rows dropped after pass %d\n",numberDropped,
	//     iLoop+1);
	if (numberDropped==lastDropped)
	  break;
	else
	  lastDropped = numberDropped;
      }
      if (paction_ == paction0)
	break;
	  
    }
  }
  if (!prob->status_) {
    paction_ = drop_zero_coefficients(prob, paction_);
#if	DEBUG_PRESOLVE
	check_sol(prob,1.0e0);
#endif

    paction_ = drop_empty_cols_action::presolve(prob, paction_);
    paction_ = drop_empty_rows_action::presolve(prob, paction_);
#if	DEBUG_PRESOLVE
	check_sol(prob,1.0e0);
#endif
  }
  
  if (prob->status_) {
    if (prob->status_==1)
	  presolvedModel_->messageHandler()->message(CLP_PRESOLVE_INFEAS,
					     presolvedModel_->messages())
					       <<prob->feasibilityTolerance_
					       <<CoinMessageEol;
    else if (prob->status_==2)
	  presolvedModel_->messageHandler()->message(CLP_PRESOLVE_UNBOUND,
					     presolvedModel_->messages())
					       <<CoinMessageEol;
    else
	  presolvedModel_->messageHandler()->message(CLP_PRESOLVE_INFEASUNBOUND,
					     presolvedModel_->messages())
					       <<CoinMessageEol;
    // get rid of data
    gutsOfDestroy();
  }
  return (paction_);
}

void check_djs(PostsolveMatrix *prob);


// We could have implemented this by having each postsolve routine
// directly call the next one, but this may make it easier to add debugging checks.
void Presolve::postsolve(PostsolveMatrix &prob)
{
  const PresolveAction *paction = paction_;

  if (prob.colstat_)
    prob.check_nbasic();
  
#if	DEBUG_PRESOLVE
  check_djs(&prob);
#endif
  
  
  while (paction) {
#if	DEBUG_PRESOLVE
    printf("POSTSOLVING %s\n", paction->name());
#endif

    paction->postsolve(&prob);
    
#if	DEBUG_PRESOLVE
    if (prob.colstat_)
      prob.check_nbasic();
#endif
    paction = paction->next;
#if	DEBUG_PRESOLVE
    check_djs(&prob);
#endif
  }    
  
#if	0 && DEBUG_PRESOLVE
  for (i=0; i<ncols0; i++) {
    if (!cdone[i]) {
      printf("!cdone[%d]\n", i);
      abort();
    }	
  }
  
  for (i=0; i<nrows0; i++) {
    if (!rdone[i]) {
      printf("!rdone[%d]\n", i);
      abort();
    }
  }
  
  
  for (i=0; i<ncols0; i++) {
    if (sol[i] < -1e10 || sol[i] > 1e10)
      printf("!!!%d %g\n", i, sol[i]);
    
  }
  
  
#endif
  
#if	0 && DEBUG_PRESOLVE
  // debug check:  make sure we ended up with same original matrix
  {
    int identical = 1;
    
    for (int i=0; i<ncols0; i++) {
      PRESOLVEASSERT(hincol[i] == &prob->mcstrt0[i+1] - &prob->mcstrt0[i]);
      CoinBigIndex kcs0 = &prob->mcstrt0[i];
      CoinBigIndex kcs = mcstrt[i];
      int n = hincol[i];
      for (int k=0; k<n; k++) {
	CoinBigIndex k1 = presolve_find_row1(&prob->hrow0[kcs0+k], kcs, kcs+n, hrow);
	
	if (k1 == kcs+n) {
	  printf("ROW %d NOT IN COL %d\n", &prob->hrow0[kcs0+k], i);
	  abort();
	}
	
	if (colels[k1] != &prob->dels0[kcs0+k])
	  printf("BAD COLEL[%d %d %d]:  %g\n",
		 k1, i, &prob->hrow0[kcs0+k], colels[k1] - &prob->dels0[kcs0+k]);
	
	if (kcs0+k != k1)
	  identical=0;
      }
    }
    printf("identical? %d\n", identical);
  }
#endif
  // put back duals
  memcpy(originalModel_->dualRowSolution(),prob.rowduals_,
	 nrows_*sizeof(double));
  double maxmin = originalModel_->getObjSense();
  if (maxmin<0.0) {
    // swap signs
    int i;
    double * pi = originalModel_->dualRowSolution();
    for (i=0;i<nrows_;i++)
      pi[i] = -pi[i];
  }
  // Now check solution
  memcpy(originalModel_->dualColumnSolution(),
	 originalModel_->objective(),ncols_*sizeof(double));
  originalModel_->transposeTimes(-1.0,
				 originalModel_->dualRowSolution(),
				 originalModel_->dualColumnSolution());
  memset(originalModel_->primalRowSolution(),0,nrows_*sizeof(double));
  originalModel_->times(1.0,originalModel_->primalColumnSolution(),
			originalModel_->primalRowSolution());
  originalModel_->checkSolution();
  presolvedModel_->messageHandler()->message(CLP_PRESOLVE_POSTSOLVE,
					    presolvedModel_->messages())
					      <<originalModel_->objectiveValue()
					      <<originalModel_->sumDualInfeasibilities()
					      <<originalModel_->numberDualInfeasibilities()
					      <<originalModel_->sumPrimalInfeasibilities()
					      <<originalModel_->numberPrimalInfeasibilities()
					       <<CoinMessageEol;
  
  //originalModel_->objectiveValue_=objectiveValue_;
  originalModel_->setNumberIterations(presolvedModel_->numberIterations());
  if (!presolvedModel_->status()) {
    if (!originalModel_->numberDualInfeasibilities()&&
	!originalModel_->numberPrimalInfeasibilities()) {
      originalModel_->setProblemStatus( 0);
    } else {
      originalModel_->setProblemStatus( -1);
      presolvedModel_->messageHandler()->message(CLP_PRESOLVE_NEEDS_CLEANING,
					    presolvedModel_->messages())
					      <<CoinMessageEol;
    }
  } else {
    originalModel_->setProblemStatus( presolvedModel_->status());
  }
}








#if	DEBUG_PRESOLVE
void check_djs(PostsolveMatrix *prob)
{
  return;
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  int *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;
  int ncols		= prob->ncols_;

  double *dcost	= prob->cost_;

  double *rcosts	= prob->rcosts_;

  double *rowduals = prob->rowduals_;

  const double maxmin	= prob->maxmin_;

  char *cdone	= prob->cdone_;

  double * csol = prob->sol_;
  double * clo = prob->clo_;
  double * cup = prob->cup_;
  int nrows = prob->nrows_;
  double * rlo = prob->rlo_;
  double * rup = prob->rup_;
  char *rdone	= prob->rdone_;

  int colx;

  double * rsol = new double[nrows];
  memset(rsol,0,nrows*sizeof(double));

  for (colx = 0; colx < ncols; ++colx) {
    if (cdone[colx]) {
      CoinBigIndex k = mcstrt[colx];
      int nx = hincol[colx];
      double dj = maxmin * dcost[colx];
      double solutionValue = csol[colx];
      for (int i=0; i<nx; ++i) {
	int row = hrow[k];
	double coeff = colels[k];
	k = link[k];
	
	dj -= rowduals[row] * coeff;
	rsol[row] += solutionValue*coeff;
      }
      if (! (fabs(rcosts[colx] - dj) < 100*ZTOLDP))
	printf("BAD DJ:  %d %g %g\n",
	       colx, rcosts[colx], dj);
      if (cup[colx]-clo[colx]>1.0e-6) {
	if (csol[colx]<clo[colx]+1.0e-6) {
	  if (dj <-1.0e-6)
	    printf("neg DJ:  %d %g  - %g %g %g\n",
		   colx, dj, clo[colx], csol[colx], cup[colx]);
	} else if (csol[colx]>cup[colx]-1.0e-6) {
	  if (dj > 1.0e-6)
	    printf("pos DJ:  %d %g  - %g %g %g\n",
		   colx, dj, clo[colx], csol[colx], cup[colx]);
	} else {
	  if (fabs(dj) >1.0e-6)
	    printf("nonzero DJ:  %d %g  - %g %g %g\n",
		   colx, dj, clo[colx], csol[colx], cup[colx]);
	}
      } 
    }
  }
  int rowx;
  for (rowx = 0; rowx < nrows; ++rowx) {
    if (rdone[rowx]) {
      if (rup[rowx]-rlo[rowx]>1.0e-6) {
	double dj = rowduals[rowx];
	if (rsol[rowx]<rlo[rowx]+1.0e-6) {
	  if (dj <-1.0e-6)
	    printf("neg rDJ:  %d %g  - %g %g %g\n",
		   rowx, dj, rlo[rowx], rsol[rowx], rup[rowx]);
	} else if (rsol[rowx]>rup[rowx]-1.0e-6) {
	  if (dj > 1.0e-6)
	    printf("pos rDJ:  %d %g  - %g %g %g\n",
		   rowx, dj, rlo[rowx], rsol[rowx], rup[rowx]);
	} else {
	  if (fabs(dj) >1.0e-6)
	    printf("nonzero rDJ:  %d %g  - %g %g %g\n",
		   rowx, dj, rlo[rowx], rsol[rowx], rup[rowx]);
	}
      } 
    }
  }
  delete [] rsol;
}
#endif


