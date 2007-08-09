// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

//#define	PRESOLVE_CONSISTENCY	1
//#define	PRESOLVE_DEBUG	1

#include <stdio.h>

#include <assert.h>
#include <iostream>

#include "CoinHelperFunctions.hpp"

#include "CoinPackedMatrix.hpp"
#include "ClpSimplex.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#endif

#include "ClpPresolve.hpp"
#include "CoinPresolveMatrix.hpp"

#include "CoinPresolveEmpty.hpp"
#include "CoinPresolveFixed.hpp"
#include "CoinPresolvePsdebug.hpp"
#include "CoinPresolveSingleton.hpp"
#include "CoinPresolveDoubleton.hpp"
#include "CoinPresolveTripleton.hpp"
#include "CoinPresolveZeros.hpp"
#include "CoinPresolveSubst.hpp"
#include "CoinPresolveForcing.hpp"
#include "CoinPresolveDual.hpp"
#include "CoinPresolveTighten.hpp"
#include "CoinPresolveUseless.hpp"
#include "CoinPresolveDupcol.hpp"
#include "CoinPresolveImpliedFree.hpp"
#include "CoinPresolveIsolated.hpp"
#include "CoinMessage.hpp"



ClpPresolve::ClpPresolve() :
  originalModel_(NULL),
  presolvedModel_(NULL),
  nonLinearValue_(0.0),
  originalColumn_(NULL),
  originalRow_(NULL),
  rowObjective_(NULL),
  paction_(0),
  ncols_(0),
  nrows_(0),
  nelems_(0),
  numberPasses_(5),
  substitution_(3),
#ifndef CLP_NO_STD
  saveFile_(""),
#endif
  presolveActions_(0)
{
}

ClpPresolve::~ClpPresolve()
{
  destroyPresolve();
}
// Gets rid of presolve actions (e.g.when infeasible)
void 
ClpPresolve::destroyPresolve()
{
 const CoinPresolveAction *paction = paction_;
  while (paction) {
    const CoinPresolveAction *next = paction->next;
    delete paction;
    paction = next;
  }
  delete [] originalColumn_;
  delete [] originalRow_;
  paction_=NULL;
  originalColumn_=NULL;
  originalRow_=NULL;
  delete [] rowObjective_;
  rowObjective_=NULL;
}

/* This version of presolve returns a pointer to a new presolved 
   model.  NULL if infeasible
*/
ClpSimplex * 
ClpPresolve::presolvedModel(ClpSimplex & si,
                            double feasibilityTolerance,
                            bool keepIntegers,
                            int numberPasses,
                            bool dropNames,
                            bool doRowObjective)
{
  // Check matrix
  if (!si.clpMatrix()->allElementsInRange(&si,si.getSmallElementValue(),
					  1.0e20))
    return NULL;
  else
    return gutsOfPresolvedModel(&si,feasibilityTolerance,keepIntegers,numberPasses,dropNames,
                                doRowObjective);
}
#ifndef CLP_NO_STD
/* This version of presolve updates
   model and saves original data to file.  Returns non-zero if infeasible
*/
int
ClpPresolve::presolvedModelToFile(ClpSimplex &si,std::string fileName,
                                  double feasibilityTolerance,
                                  bool keepIntegers,
                                  int numberPasses,
                                  bool doRowObjective)
{
  // Check matrix
  if (!si.clpMatrix()->allElementsInRange(&si,si.getSmallElementValue(),
					  1.0e20))
    return 2;
  saveFile_=fileName;
  si.saveModel(saveFile_.c_str());
  ClpSimplex * model = gutsOfPresolvedModel(&si,feasibilityTolerance,keepIntegers,numberPasses,true,
                                            doRowObjective);
  if (model==&si) {
    return 0;
  } else {
    si.restoreModel(saveFile_.c_str());
    remove(saveFile_.c_str());
    return 1;
  }
}
#endif
// Return pointer to presolved model
ClpSimplex * 
ClpPresolve::model() const
{
  return presolvedModel_;
}
// Return pointer to original model
ClpSimplex * 
ClpPresolve::originalModel() const
{
  return originalModel_;
}
void 
ClpPresolve::postsolve(bool updateStatus)
{
  // Return at once if no presolved model
  if (!presolvedModel_)
    return;
  // Messages
  CoinMessages messages = originalModel_->coinMessages();
  if (!presolvedModel_->isProvenOptimal()) {
    presolvedModel_->messageHandler()->message(COIN_PRESOLVE_NONOPTIMAL,
					     messages)
					       <<CoinMessageEol;
  }

  // this is the size of the original problem
  const int ncols0  = ncols_;
  const int nrows0  = nrows_;
  const CoinBigIndex nelems0 = nelems_;
  
  // this is the reduced problem
  int ncols = presolvedModel_->getNumCols();
  int nrows = presolvedModel_->getNumRows();

  double * acts=NULL;
  double * sol =NULL;
  unsigned char * rowstat=NULL;
  unsigned char * colstat = NULL;
#ifndef CLP_NO_STD
  if (saveFile_=="") {
#endif
    // reality check
    assert(ncols0==originalModel_->getNumCols());
    assert(nrows0==originalModel_->getNumRows());
    acts = originalModel_->primalRowSolution();
    sol  = originalModel_->primalColumnSolution();
    if (updateStatus) {
      // postsolve does not know about fixed
      int i;
      for (i=0;i<nrows+ncols;i++) {
        if (presolvedModel_->getColumnStatus(i)==ClpSimplex::isFixed)
          presolvedModel_->setColumnStatus(i,ClpSimplex::atLowerBound);
      }
      unsigned char *status = originalModel_->statusArray();
      if (!status) {
        originalModel_->createStatus();
        status = originalModel_->statusArray();
      }
      rowstat = status + ncols0;
      colstat = status;
      memcpy(colstat, presolvedModel_->statusArray(), ncols);
      memcpy(rowstat, presolvedModel_->statusArray()+ncols, nrows);
    }
#ifndef CLP_NO_STD
  } else {
    // from file
    acts = new double[nrows0];
    sol  = new double[ncols0];
    CoinZeroN(acts,nrows0);
    CoinZeroN(sol,ncols0);
    if (updateStatus) {
      unsigned char *status = new unsigned char [nrows0+ncols0];
      rowstat = status + ncols0;
      colstat = status;
      memcpy(colstat, presolvedModel_->statusArray(), ncols);
      memcpy(rowstat, presolvedModel_->statusArray()+ncols, nrows);
    }
  }
#endif

  // CoinPostsolveMatrix object assumes ownership of sol, acts, colstat;
  // will be deleted by ~CoinPostsolveMatrix. delete[] operations below
  // cause duplicate free. In case where saveFile == "", as best I can see
  // arrays are owned by originalModel_. fix is to
  // clear fields in prob to avoid delete[] in ~CoinPostsolveMatrix.
  CoinPostsolveMatrix prob(presolvedModel_,
		       ncols0,
		       nrows0,
		       nelems0,
		       presolvedModel_->getObjSense(),
		       // end prepost
		       
		       sol, acts,
		       colstat, rowstat);
    
  postsolve(prob);

#ifndef CLP_NO_STD
  if (saveFile_!="") {
    // From file
    assert (originalModel_==presolvedModel_);
    originalModel_->restoreModel(saveFile_.c_str());
    remove(saveFile_.c_str());
    memcpy(originalModel_->primalRowSolution(),acts,nrows0*sizeof(double));
    // delete [] acts;
    memcpy(originalModel_->primalColumnSolution(),sol,ncols0*sizeof(double));
    // delete [] sol;
    if (updateStatus) {
      memcpy(originalModel_->statusArray(),colstat,nrows0+ncols0);
      // delete [] colstat;
    }
  } else {
#endif
    prob.sol_ = 0 ;
    prob.acts_ = 0 ;
    prob.colstat_ = 0 ;
#ifndef CLP_NO_STD
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
  double offset;
  memcpy(originalModel_->dualColumnSolution(),
	 originalModel_->objectiveAsObject()->gradient(originalModel_,
						       originalModel_->primalColumnSolution(),
						       offset,true),ncols_*sizeof(double));
  originalModel_->transposeTimes(-1.0,
				 originalModel_->dualRowSolution(),
				 originalModel_->dualColumnSolution());
  memset(originalModel_->primalRowSolution(),0,nrows_*sizeof(double));
  originalModel_->times(1.0,originalModel_->primalColumnSolution(),
			originalModel_->primalRowSolution());
  originalModel_->checkSolutionInternal();
  // Messages
  presolvedModel_->messageHandler()->message(COIN_PRESOLVE_POSTSOLVE,
					    messages)
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
      // Say not optimal after presolve
      originalModel_->setSecondaryStatus(7);
      presolvedModel_->messageHandler()->message(COIN_PRESOLVE_NEEDS_CLEANING,
					    messages)
					      <<CoinMessageEol;
    }
  } else {
    originalModel_->setProblemStatus( presolvedModel_->status());
  }
#ifndef CLP_NO_STD
  if (saveFile_!="") 
    presolvedModel_=NULL;
#endif
}

// return pointer to original columns
const int * 
ClpPresolve::originalColumns() const
{
  return originalColumn_;
}
// return pointer to original rows
const int * 
ClpPresolve::originalRows() const
{
  return originalRow_;
}
// Set pointer to original model
void 
ClpPresolve::setOriginalModel(ClpSimplex * model)
{
  originalModel_=model;
}
#if 0
// A lazy way to restrict which transformations are applied
// during debugging.
static int ATOI(const char *name)
{
 return true;
#if	PRESOLVE_DEBUG || PRESOLVE_SUMMARY
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
//#define PRESOLVE_DEBUG 1
#if PRESOLVE_DEBUG
void check_sol(CoinPresolveMatrix *prob,double tol)
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
// customized by subclassing CoinPresolve.
const CoinPresolveAction *ClpPresolve::presolve(CoinPresolveMatrix *prob)
{
  // Messages
  CoinMessages messages = CoinMessage(prob->messages().language());
  paction_ = 0;

  prob->status_=0; // say feasible
  paction_ = make_fixed(prob, paction_);
  // if integers then switch off dual stuff
  // later just do individually
  bool doDualStuff = (presolvedModel_->integerInformation()==NULL);
  // but allow in some cases
  if ((presolveActions_&512)!=0)
    doDualStuff=true;
  if (prob->anyProhibited()) 
    doDualStuff=false;
  if (!doDual())
    doDualStuff=false;
#if	PRESOLVE_CONSISTENCY
//  presolve_links_ok(prob->rlink_, prob->mrstrt_, prob->hinrow_, prob->nrows_);
    presolve_links_ok(prob,false,true) ;
#endif

  if (!prob->status_) {
    bool slackSingleton = doSingletonColumn();
    slackSingleton = false;
    const bool slackd = doSingleton();
    const bool doubleton = doDoubleton();
    const bool tripleton = doTripleton();
    const bool forcing = doForcing();
    const bool ifree = doImpliedFree();
    const bool zerocost = doTighten();
    const bool dupcol = doDupcol();
    const bool duprow = doDuprow();
    const bool dual = doDualStuff;
    
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
#if	PRESOLVE_DEBUG
    check_sol(prob,1.0e0);
#endif
    if (dupcol) {
      //paction_ = dupcol_action::presolve(prob, paction_);
    }

    // Check number rows dropped
    int lastDropped=0;
    prob->pass_=0;
    for (iLoop=0;iLoop<numberPasses_;iLoop++) {
    // See if we want statistics
    if ((presolveActions_&0x80000000)!=0)
      printf("Starting major pass %d after %g seconds\n",iLoop+1,CoinCpuTime()-prob->startTime_);
#ifdef PRESOLVE_SUMMARY
      printf("Starting major pass %d\n",iLoop+1);
#endif
      const CoinPresolveAction * const paction0 = paction_;
      // look for substitutions with no fill
      int fill_level=2;
      //fill_level=10;
      //printf("** fill_level == 10 !\n");
      int whichPass=0;
      while (1) {
	whichPass++;
	prob->pass_++;
	const CoinPresolveAction * const paction1 = paction_;

	if (slackd) {
	  bool notFinished = true;
	  while (notFinished) 
	    paction_ = slack_doubleton_action::presolve(prob, paction_,
							notFinished);
	  if (prob->status_)
	    break;
	}
	if (dual&&whichPass==1) {
	  // this can also make E rows so do one bit here
	  paction_ = remove_dual_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	}

	if (doubleton) {
	  paction_ = doubleton_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	}
	if (tripleton) {
	  paction_ = tripleton_action::presolve(prob, paction_);
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

	if (ifree&&(whichPass%5)==1) {
	paction_ = implied_free_action::presolve(prob, paction_,fill_level);
	if (prob->status_)
	  break;
	}

#if	PRESOLVE_DEBUG
	check_sol(prob,1.0e0);
#endif

#if	PRESOLVE_CONSISTENCY
//	presolve_links_ok(prob->rlink_, prob->mrstrt_, prob->hinrow_, 
//			  prob->nrows_);
	presolve_links_ok(prob,false,true) ;
#endif

//#if	PRESOLVE_DEBUG
//	presolve_no_zeros(prob->mcstrt_, prob->colels_, prob->hincol_, 
//			  prob->ncols_);
//#endif
//#if	PRESOLVE_CONSISTENCY
//	prob->consistent();
//#endif
#if	PRESOLVE_CONSISTENCY
	presolve_no_zeros(prob,true,false) ;
	presolve_consistent(prob,true) ;
#endif

	  
	// set up for next pass
	// later do faster if many changes i.e. memset and memcpy
	prob->numberRowsToDo_ = prob->numberNextRowsToDo_;
	int kcheck;
	bool found=false;
	kcheck=-1;
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
	  paction_ = remove_dual_action::presolve(prob, paction_);
	  if (prob->status_)
	    break;
	  const CoinPresolveAction * const paction2 = paction_;
	  if (ifree) {
	    //int fill_level=0; // switches off substitution
	    paction_ = implied_free_action::presolve(prob, paction_,fill_level);
	    if (prob->status_)
	      break;
	  }
	  if (paction_ == paction2)
	    break;
	}
      } else if (ifree) {
	// just free
	int fill_level=0; // switches off substitution
	paction_ = implied_free_action::presolve(prob, paction_,fill_level);
	if (prob->status_)
	  break;
      }
#if	PRESOLVE_DEBUG
      check_sol(prob,1.0e0);
#endif
      if (dupcol) {
        // maybe allow integer columns to be checked
        if ((presolveActions_&512)!=0)
          prob->setPresolveOptions(prob->presolveOptions()|1);
	paction_ = dupcol_action::presolve(prob, paction_);
	if (prob->status_)
	  break;
      }
#if	PRESOLVE_DEBUG
	check_sol(prob,1.0e0);
#endif
      
      if (duprow) {
	paction_ = duprow_action::presolve(prob, paction_);
	if (prob->status_)
	  break;
      }
#if	PRESOLVE_DEBUG
      check_sol(prob,1.0e0);
#endif
      bool stopLoop=false;
      {
	int * hinrow = prob->hinrow_;
	int numberDropped=0;
	for (i=0;i<nrows_;i++) 
	  if (!hinrow[i])
	    numberDropped++;
	
	prob->messageHandler()->message(COIN_PRESOLVE_PASS,
					messages)
					  <<numberDropped<<iLoop+1
					  <<CoinMessageEol;
	//printf("%d rows dropped after pass %d\n",numberDropped,
	//     iLoop+1);
	if (numberDropped==lastDropped)
	  stopLoop=true;
	else
	  lastDropped = numberDropped;
      }
      // Do this here as not very loopy
      if (slackSingleton) {
        // On most passes do not touch costed slacks
        if (paction_ != paction0&&!stopLoop) {
          paction_ = slack_singleton_action::presolve(prob, paction_,NULL);
        } else {
          // do costed if Clp (at end as ruins rest of presolve)
          paction_ = slack_singleton_action::presolve(prob, paction_,rowObjective_);
          stopLoop=true;
        }
      }
#if	PRESOLVE_DEBUG
      check_sol(prob,1.0e0);
#endif
      if (paction_ == paction0||stopLoop)
	break;
	  
    }
  }
  if (!prob->status_) {
    paction_ = drop_zero_coefficients(prob, paction_);
#if	PRESOLVE_DEBUG
	check_sol(prob,1.0e0);
#endif

    paction_ = drop_empty_cols_action::presolve(prob, paction_);
    paction_ = drop_empty_rows_action::presolve(prob, paction_);
#if	PRESOLVE_DEBUG
	check_sol(prob,1.0e0);
#endif
  }
  
  if (prob->status_) {
    if (prob->status_==1) 
	  prob->messageHandler()->message(COIN_PRESOLVE_INFEAS,
					     messages)
					       <<prob->feasibilityTolerance_
					       <<CoinMessageEol;
    else if (prob->status_==2)
	  prob->messageHandler()->message(COIN_PRESOLVE_UNBOUND,
					     messages)
					       <<CoinMessageEol;
    else
	  prob->messageHandler()->message(COIN_PRESOLVE_INFEASUNBOUND,
					     messages)
					       <<CoinMessageEol;
    // get rid of data
    destroyPresolve();
  }
  return (paction_);
}

void check_djs(CoinPostsolveMatrix *prob);


// We could have implemented this by having each postsolve routine
// directly call the next one, but this may make it easier to add debugging checks.
void ClpPresolve::postsolve(CoinPostsolveMatrix &prob)
{
  {
    // Check activities
    double *colels	= prob.colels_;
    int *hrow		= prob.hrow_;
    CoinBigIndex *mcstrt		= prob.mcstrt_;
    int *hincol		= prob.hincol_;
    int *link		= prob.link_;
    int ncols		= prob.ncols_;

    char *cdone	= prob.cdone_;

    double * csol = prob.sol_;
    int nrows = prob.nrows_;

    int colx;

    double * rsol = prob.acts_;
    memset(rsol,0,nrows*sizeof(double));

    for (colx = 0; colx < ncols; ++colx) {
      if (cdone[colx]) {
	CoinBigIndex k = mcstrt[colx];
	int nx = hincol[colx];
	double solutionValue = csol[colx];
	for (int i=0; i<nx; ++i) {
	  int row = hrow[k];
	  double coeff = colels[k];
	  k = link[k];
	  rsol[row] += solutionValue*coeff;
	}
      }
    }
  }
  const CoinPresolveAction *paction = paction_;
  //#define PRESOLVE_DEBUG 1
#if	PRESOLVE_DEBUG
  // Check only works after first one
  int checkit=-1;
#endif
  
  while (paction) {
#if	PRESOLVE_DEBUG
    printf("POSTSOLVING %s\n", paction->name());
#endif

    paction->postsolve(&prob);
    
#if	PRESOLVE_DEBUG
  {
    int nr=0;
    int i;
    for (i=0;i<prob.nrows_;i++) {
      if ((prob.rowstat_[i]&7)==1)
        nr++;
    }
    int nc=0;
    for (i=0;i<prob.ncols_;i++) {
      if ((prob.colstat_[i]&7)==1)
        nc++;
    }
    printf("%d rows (%d basic), %d cols (%d basic)\n",prob.nrows_,nr,prob.ncols_,nc);
  }
    checkit++;
    if (prob.colstat_&&checkit>0) {
      presolve_check_nbasic(&prob) ;
      presolve_check_sol(&prob,2,2,1) ;
    }
#endif
    paction = paction->next;
#if	PRESOLVE_DEBUG
//  check_djs(&prob);
    if (checkit>0)
      presolve_check_reduced_costs(&prob) ;
#endif
  }    
#if	PRESOLVE_DEBUG
  if (prob.colstat_) {
    presolve_check_nbasic(&prob) ;
    presolve_check_sol(&prob,2,2,1) ;
  }
#endif
#undef PRESOLVE_DEBUG
  
#if	0 && PRESOLVE_DEBUG
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
  
#if	0 && PRESOLVE_DEBUG
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
}








static inline double getTolerance(const ClpSimplex  *si, ClpDblParam key)
{
  double tol;
  if (! si->getDblParam(key, tol)) {
    CoinPresolveAction::throwCoinError("getDblParam failed",
				      "CoinPrePostsolveMatrix::CoinPrePostsolveMatrix");
  }
  return (tol);
}


// Assumptions:
// 1. nrows>=m.getNumRows()
// 2. ncols>=m.getNumCols()
//
// In presolve, these values are equal.
// In postsolve, they may be inequal, since the reduced problem
// may be smaller, but we need room for the large problem.
// ncols may be larger than si.getNumCols() in postsolve,
// this at that point si will be the reduced problem,
// but we need to reserve enough space for the original problem.
CoinPrePostsolveMatrix::CoinPrePostsolveMatrix(const ClpSimplex * si,
					     int ncols_in,
					     int nrows_in,
					     CoinBigIndex nelems_in,
                                               double bulkRatio) 
  : ncols_(si->getNumCols()),
    nrows_(si->getNumRows()),
    nelems_(si->getNumElements()),
    ncols0_(ncols_in),
    nrows0_(nrows_in),
    bulkRatio_(bulkRatio),
    mcstrt_(new CoinBigIndex[ncols_in+1]),
    hincol_(new int[ncols_in+1]),
    cost_(new double[ncols_in]),
    clo_(new double[ncols_in]),
    cup_(new double[ncols_in]),
    rlo_(new double[nrows_in]),
    rup_(new double[nrows_in]),
    originalColumn_(new int[ncols_in]),
    originalRow_(new int[nrows_in]),
    ztolzb_(getTolerance(si, ClpPrimalTolerance)),
    ztoldj_(getTolerance(si, ClpDualTolerance)),
    maxmin_(si->getObjSense()),
    sol_(NULL),
    rowduals_(NULL),
    acts_(NULL),
    rcosts_(NULL),
    colstat_(NULL),
    rowstat_(NULL),
    handler_(NULL),
    defaultHandler_(false),
    messages_()

{
  bulk0_ = (CoinBigIndex) (bulkRatio_*nelems_in);
  hrow_  = new int   [bulk0_];
  colels_ = new double[bulk0_];
  si->getDblParam(ClpObjOffset,originalOffset_);
  int ncols = si->getNumCols();
  int nrows = si->getNumRows();

  setMessageHandler(si->messageHandler()) ;

  ClpDisjointCopyN(si->getColLower(), ncols, clo_);
  ClpDisjointCopyN(si->getColUpper(), ncols, cup_);
  //ClpDisjointCopyN(si->getObjCoefficients(), ncols, cost_);
  double offset;
  ClpDisjointCopyN(si->objectiveAsObject()->gradient(si,si->getColSolution(),offset,true), ncols, cost_);
  ClpDisjointCopyN(si->getRowLower(), nrows,  rlo_);
  ClpDisjointCopyN(si->getRowUpper(), nrows,  rup_);
  int i;
  for (i=0;i<ncols_in;i++) 
    originalColumn_[i]=i;
  for (i=0;i<nrows_in;i++) 
    originalRow_[i]=i;
  sol_=NULL;
  rowduals_=NULL;
  acts_=NULL;

  rcosts_=NULL;
  colstat_=NULL;
  rowstat_=NULL;
}

// I am not familiar enough with CoinPackedMatrix to be confident
// that I will implement a row-ordered version of toColumnOrderedGapFree
// properly.
static bool isGapFree(const CoinPackedMatrix& matrix)
{
  const CoinBigIndex * start = matrix.getVectorStarts();
  const int * length = matrix.getVectorLengths();
  int i = matrix.getSizeVectorLengths() - 1;
  // Quick check
  if (matrix.getNumElements()==start[i]) {
    return true;
  } else {
    for (i = matrix.getSizeVectorLengths() - 1; i >= 0; --i) {
      if (start[i+1] - start[i] != length[i])
        break;
    }
    return (! (i >= 0));
  }
}
#if	PRESOLVE_DEBUG
static void matrix_bounds_ok(const double *lo, const double *up, int n)
{
  int i;
  for (i=0; i<n; i++) {
    PRESOLVEASSERT(lo[i] <= up[i]);
    PRESOLVEASSERT(lo[i] < PRESOLVE_INF);
    PRESOLVEASSERT(-PRESOLVE_INF < up[i]);
  }
}
#endif
CoinPresolveMatrix::CoinPresolveMatrix(int ncols0_in,
				     double maxmin,
				     // end prepost members

				     ClpSimplex * si,

				     // rowrep
				     int nrows_in,
				     CoinBigIndex nelems_in,
                                       bool doStatus,
                                       double nonLinearValue,
                                       double bulkRatio) :

  CoinPrePostsolveMatrix(si,
			ncols0_in, nrows_in, nelems_in,bulkRatio),
  clink_(new presolvehlink[ncols0_in+1]),
  rlink_(new presolvehlink[nrows_in+1]),

  dobias_(0.0),


  // temporary init
  integerType_(new unsigned char[ncols0_in]),
  tuning_(false),
  startTime_(0.0),
  feasibilityTolerance_(0.0),
  status_(-1),
  colsToDo_(new int [ncols0_in]),
  numberColsToDo_(0),
  nextColsToDo_(new int[ncols0_in]),
  numberNextColsToDo_(0),
  rowsToDo_(new int [nrows_in]),
  numberRowsToDo_(0),
  nextRowsToDo_(new int[nrows_in]),
  numberNextRowsToDo_(0),
  presolveOptions_(0)

{
  const int bufsize = bulk0_;

  nrows_ = si->getNumRows() ;

  // Set up change bits etc
  rowChanged_ = new unsigned char[nrows_];
  memset(rowChanged_,0,nrows_);
  colChanged_ = new unsigned char[ncols_];
  memset(colChanged_,0,ncols_);
  CoinPackedMatrix * m = si->matrix();

  // The coefficient matrix is a big hunk of stuff.
  // Do the copy here to try to avoid running out of memory.

  const CoinBigIndex * start = m->getVectorStarts();
  const int * length = m->getVectorLengths();
  const int * row = m->getIndices();
  const double * element = m->getElements();
  int icol,nel=0;
  mcstrt_[0]=0;
  for (icol=0;icol<ncols_;icol++) {
    int j;
    for (j=start[icol];j<start[icol]+length[icol];j++) {
      hrow_[nel]=row[j];
      colels_[nel++]=element[j];
    }
    mcstrt_[icol+1]=nel;
  }
  assert(mcstrt_[ncols_] == nelems_);
  ClpDisjointCopyN(m->getVectorLengths(),ncols_,  hincol_);

  // same thing for row rep
  CoinPackedMatrix * mRow = new CoinPackedMatrix();
  mRow->reverseOrderedCopyOf(*m);
  mRow->removeGaps();
  mRow->setExtraGap(0.0);

  // Now get rid of matrix
  si->createEmptyMatrix();

  double * el = mRow->getMutableElements();
  int * ind = mRow->getMutableIndices();
  CoinBigIndex * strt = mRow->getMutableVectorStarts();
  int * len = mRow->getMutableVectorLengths();
  // Do carefully to save memory
  rowels_ = new double[bulk0_];
  ClpDisjointCopyN(el,      nelems_, rowels_);
  mRow->nullElementArray();
  delete [] el;
  hcol_ = new int[bulk0_];
  ClpDisjointCopyN(ind,       nelems_, hcol_);
  mRow->nullIndexArray();
  delete [] ind;
  mrstrt_ = new CoinBigIndex[nrows_in+1];
  ClpDisjointCopyN(strt,  nrows_,  mrstrt_);
  mRow->nullStartArray();
  mrstrt_[nrows_] = nelems_;
  delete [] strt;
  hinrow_ = new int[nrows_in+1];
  ClpDisjointCopyN(len, nrows_,  hinrow_);

  delete mRow;
  if (si->integerInformation()) {
    memcpy(integerType_,si->integerInformation(),ncols_*sizeof(char));
  } else {
    ClpFillN<unsigned char>(integerType_, ncols_, (unsigned char) 0);
  }

#ifndef SLIM_CLP
#ifndef NO_RTTI
    ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(si->objectiveAsObject()));
#else
    ClpQuadraticObjective * quadraticObj = NULL;
    if (si->objectiveAsObject()->type()==2)
      quadraticObj = (static_cast< ClpQuadraticObjective*>(si->objectiveAsObject()));
#endif
#endif
  // Set up prohibited bits if needed
  if (nonLinearValue) {
    anyProhibited_ = true;
    for (icol=0;icol<ncols_;icol++) {
      int j;
      bool nonLinearColumn = false;
      if (cost_[icol]==nonLinearValue)
	nonLinearColumn=true;
      for (j=mcstrt_[icol];j<mcstrt_[icol+1];j++) {
	if (colels_[j]==nonLinearValue) {
	  nonLinearColumn=true;
	  setRowProhibited(hrow_[j]);
	}
      }
      if (nonLinearColumn)
	setColProhibited(icol);
    }
#ifndef SLIM_CLP
  } else if (quadraticObj) {
    CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
    //const int * columnQuadratic = quadratic->getIndices();
    //const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
    const int * columnQuadraticLength = quadratic->getVectorLengths();
    //double * quadraticElement = quadratic->getMutableElements();
    int numberColumns = quadratic->getNumCols();
    anyProhibited_ = true;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnQuadraticLength[iColumn]) {
	setColProhibited(iColumn);
	//printf("%d prohib\n",iColumn);
      }
    }
#endif
  } else {
    anyProhibited_ = false;
  }

  if (doStatus) {
    // allow for status and solution
    sol_ = new double[ncols_];
    memcpy(sol_,si->primalColumnSolution(),ncols_*sizeof(double));;
    acts_ = new double [nrows_];
    memcpy(acts_,si->primalRowSolution(),nrows_*sizeof(double));
    if (!si->statusArray())
      si->createStatus();
    colstat_ = new unsigned char [nrows_+ncols_];
    memcpy(colstat_,si->statusArray(),
	   (nrows_+ncols_)*sizeof(unsigned char));
    rowstat_ = colstat_+ncols_;
  }

  // the original model's fields are now unneeded - free them
  
  si->resize(0,0);

#if	PRESOLVE_DEBUG
  matrix_bounds_ok(rlo_, rup_, nrows_);
  matrix_bounds_ok(clo_, cup_, ncols_);
#endif

#if 0
  for (i=0; i<nrows; ++i)
    printf("NR: %6d\n", hinrow[i]);
  for (int i=0; i<ncols; ++i)
    printf("NC: %6d\n", hincol[i]);
#endif

  presolve_make_memlists(mcstrt_, hincol_, clink_, ncols_);
  presolve_make_memlists(mrstrt_, hinrow_, rlink_, nrows_);

  // this allows last col/row to expand up to bufsize-1 (22);
  // this must come after the calls to presolve_prefix
  mcstrt_[ncols_] = bufsize-1;
  mrstrt_[nrows_] = bufsize-1;

#if	PRESOLVE_CONSISTENCY
//consistent(false);
  presolve_consistent(this,false) ;
#endif
}

void CoinPresolveMatrix::update_model(ClpSimplex * si,
				     int nrows0,
				     int ncols0,
				     CoinBigIndex nelems0)
{
  si->loadProblem(ncols_, nrows_, mcstrt_, hrow_, colels_, hincol_,
		 clo_, cup_, cost_, rlo_, rup_);
  //delete [] si->integerInformation();
  int numberIntegers=0;
  for (int i=0; i<ncols_; i++) {
    if (integerType_[i])
      numberIntegers++;
  }
  if (numberIntegers) 
    si->copyInIntegerInformation((const char *) integerType_);
  else
    si->copyInIntegerInformation(NULL);

#if	PRESOLVE_SUMMARY
  printf("NEW NCOL/NROW/NELS:  %d(-%d) %d(-%d) %d(-%d)\n",
	 ncols_, ncols0-ncols_,
	 nrows_, nrows0-nrows_,
	 si->getNumElements(), nelems0-si->getNumElements());
#endif
  si->setDblParam(ClpObjOffset,originalOffset_-dobias_);

}











////////////////  POSTSOLVE

CoinPostsolveMatrix::CoinPostsolveMatrix(ClpSimplex*  si,
				       int ncols0_in,
				       int nrows0_in,
				       CoinBigIndex nelems0,
				   
				       double maxmin,
				       // end prepost members

				       double *sol_in,
				       double *acts_in,

				       unsigned char *colstat_in,
				       unsigned char *rowstat_in) :
  CoinPrePostsolveMatrix(si,
			ncols0_in, nrows0_in, nelems0,2.0),

  free_list_(0),
  // link, free_list, maxlink
  maxlink_(bulk0_),
  link_(new int[/*maxlink*/ bulk0_]),
      
  cdone_(new char[ncols0_]),
  rdone_(new char[nrows0_in])

{
  bulk0_ = maxlink_ ;
  nrows_ = si->getNumRows() ;
  ncols_ = si->getNumCols() ;

  sol_=sol_in;
  rowduals_=NULL;
  acts_=acts_in;

  rcosts_=NULL;
  colstat_=colstat_in;
  rowstat_=rowstat_in;

  // this is the *reduced* model, which is probably smaller
  int ncols1 = ncols_ ;
  int nrows1 = nrows_ ;

  const CoinPackedMatrix * m = si->matrix();

  const CoinBigIndex nelemsr = m->getNumElements();
  if (m->getNumElements()&&!isGapFree(*m)) {
    // Odd - gaps
    CoinPackedMatrix mm(*m);
    mm.removeGaps();
    mm.setExtraGap(0.0);
    
    ClpDisjointCopyN(mm.getVectorStarts(), ncols1, mcstrt_);
    CoinZeroN(mcstrt_+ncols1,ncols0_-ncols1);
    mcstrt_[ncols1] = nelems0;	// ??    (should point to end of bulk store   -- lh --)
    ClpDisjointCopyN(mm.getVectorLengths(),ncols1,  hincol_);
    ClpDisjointCopyN(mm.getIndices(),      nelemsr, hrow_);
    ClpDisjointCopyN(mm.getElements(),     nelemsr, colels_);
  } else {
    // No gaps
    
    ClpDisjointCopyN(m->getVectorStarts(), ncols1, mcstrt_);
    CoinZeroN(mcstrt_+ncols1,ncols0_-ncols1);
    mcstrt_[ncols1] = nelems0;	// ??    (should point to end of bulk store   -- lh --)
    ClpDisjointCopyN(m->getVectorLengths(),ncols1,  hincol_);
    ClpDisjointCopyN(m->getIndices(),      nelemsr, hrow_);
    ClpDisjointCopyN(m->getElements(),     nelemsr, colels_);
  }



#if	0 && PRESOLVE_DEBUG
  presolve_check_costs(model, &colcopy);
#endif

  // This determines the size of the data structure that contains
  // the matrix being postsolved.  Links are taken from the free_list
  // to recreate matrix entries that were presolved away,
  // and links are added to the free_list when entries created during
  // presolve are discarded.  There is never a need to gc this list.
  // Naturally, it should contain
  // exactly nelems0 entries "in use" when postsolving is done,
  // but I don't know whether the matrix could temporarily get
  // larger during postsolving.  Substitution into more than two
  // rows could do that, in principle.  I am being very conservative
  // here by reserving much more than the amount of space I probably need.
  // If this guess is wrong, check_free_list may be called.
  //  int bufsize = 2*nelems0;

  memset(cdone_, -1, ncols0_);
  memset(rdone_, -1, nrows0_);

  rowduals_ = new double[nrows0_];
  ClpDisjointCopyN(si->getRowPrice(), nrows1, rowduals_);

  rcosts_ = new double[ncols0_];
  ClpDisjointCopyN(si->getReducedCost(), ncols1, rcosts_);
  if (maxmin<0.0) {
    // change so will look as if minimize
    int i;
    for (i=0;i<nrows1;i++)
      rowduals_[i] = - rowduals_[i];
    for (i=0;i<ncols1;i++) {
      rcosts_[i] = - rcosts_[i];
    }
  }

  //ClpDisjointCopyN(si->getRowUpper(), nrows1, rup_);
  //ClpDisjointCopyN(si->getRowLower(), nrows1, rlo_);

  ClpDisjointCopyN(si->getColSolution(), ncols1, sol_);
  si->setDblParam(ClpObjOffset,originalOffset_);

  for (int j=0; j<ncols1; j++) {
    CoinBigIndex kcs = mcstrt_[j];
    CoinBigIndex kce = kcs + hincol_[j];
    for (CoinBigIndex k=kcs; k<kce; ++k) {
      link_[k] = k+1;
    }
    link_[kce-1] = NO_LINK ;
  }
  {
    int ml = maxlink_;
    for (CoinBigIndex k=nelemsr; k<ml; ++k)
      link_[k] = k+1;
    if (ml)
      link_[ml-1] = NO_LINK;
  }
  free_list_ = nelemsr;
# if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
/*
  These are used to track the action of postsolve transforms during debugging.
*/
  CoinFillN(cdone_,ncols1,PRESENT_IN_REDUCED) ;
  CoinZeroN(cdone_+ncols1,ncols0_in-ncols1) ;
  CoinFillN(rdone_,nrows1,PRESENT_IN_REDUCED) ;
  CoinZeroN(rdone_+nrows1,nrows0_in-nrows1) ;
# endif
}
/* This is main part of Presolve */
ClpSimplex * 
ClpPresolve::gutsOfPresolvedModel(ClpSimplex * originalModel,
				  double feasibilityTolerance,
				  bool keepIntegers,
				  int numberPasses,
				  bool dropNames,
                                  bool doRowObjective)
{
  ncols_ = originalModel->getNumCols();
  nrows_ = originalModel->getNumRows();
  nelems_ = originalModel->getNumElements();
  numberPasses_ = numberPasses;

  double maxmin = originalModel->getObjSense();
  originalModel_ = originalModel;
  delete [] originalColumn_;
  originalColumn_ = new int[ncols_];
  delete [] originalRow_;
  originalRow_ = new int[nrows_];
  // and fill in case returns early
  int i;
  for (i=0;i<ncols_;i++)
    originalColumn_[i]=i;
  for (i=0;i<nrows_;i++)
    originalRow_[i]=i;
  delete [] rowObjective_;
  if (doRowObjective) {
    rowObjective_ = new double [nrows_];
    memset(rowObjective_,0,nrows_*sizeof(double));
  } else {
    rowObjective_=NULL;
  }

  // result is 0 - okay, 1 infeasible, -1 go round again, 2 - original model
  int result = -1;
  
  // User may have deleted - its their responsibility
  presolvedModel_=NULL;
  // Messages
  CoinMessages messages = originalModel->coinMessages();
  while (result==-1) {

#ifndef CLP_NO_STD
    // make new copy
    if (saveFile_=="") {
#endif
      delete presolvedModel_;
#ifndef CLP_NO_STD
      // So won't get names
      int lengthNames = originalModel->lengthNames();
      originalModel->setLengthNames(0);
#endif
      presolvedModel_ = new ClpSimplex(*originalModel);
#ifndef CLP_NO_STD
      originalModel->setLengthNames(lengthNames);
    } else {
      presolvedModel_=originalModel;
    }
    presolvedModel_->dropNames();
#endif

    // drop integer information if wanted
    if (!keepIntegers)
      presolvedModel_->deleteIntegerInformation();

    double ratio=2.0;
    if (substitution_>3)
      ratio=substitution_;
    else if (substitution_==2)
      ratio=1.5;
    CoinPresolveMatrix prob(ncols_,
			maxmin,
			presolvedModel_,
			nrows_, nelems_,true,nonLinearValue_,ratio);
    prob.setMaximumSubstitutionLevel(substitution_);
    if (doRowObjective) 
      memset(rowObjective_,0,nrows_*sizeof(double));
    // See if we want statistics
    if ((presolveActions_&0x80000000)!=0)
      prob.statistics();
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

    // move across feasibility tolerance
    prob.feasibilityTolerance_ = feasibilityTolerance;

    // Do presolve
    paction_ = presolve(&prob);

    result =0; 

    if (prob.status_==0&&paction_) {
      // Looks feasible but double check to see if anything slipped through
      int n		= prob.ncols_;
      double * lo = prob.clo_;
      double * up = prob.cup_;
      int i;
      
      for (i=0;i<n;i++) {
	if (up[i]<lo[i]) {
	  if (up[i]<lo[i]-1.0e-8) {
	    // infeasible
	    prob.status_=1;
	  } else {
	    up[i]=lo[i];
	  }
	}
      }
      
      n = prob.nrows_;
      lo = prob.rlo_;
      up = prob.rup_;

      for (i=0;i<n;i++) {
	if (up[i]<lo[i]) {
	  if (up[i]<lo[i]-1.0e-8) {
	    // infeasible
	    prob.status_=1;
	  } else {
	    up[i]=lo[i];
	  }
	}
      }
    }
    if (prob.status_==0&&paction_) {
      // feasible
    
      prob.update_model(presolvedModel_, nrows_, ncols_, nelems_);
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
#ifndef SLIM_CLP
#ifndef NO_RTTI
      ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(originalModel->objectiveAsObject()));
#else
      ClpQuadraticObjective * quadraticObj = NULL;
      if (originalModel->objectiveAsObject()->type()==2)
	quadraticObj = (static_cast< ClpQuadraticObjective*>(originalModel->objectiveAsObject()));
#endif
      if (quadraticObj) {
	// set up for subset
	char * mark = new char [ncols_];
	memset(mark,0,ncols_);
	CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
	//const int * columnQuadratic = quadratic->getIndices();
	//const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	const int * columnQuadraticLength = quadratic->getVectorLengths();
	//double * quadraticElement = quadratic->getMutableElements();
	int numberColumns = quadratic->getNumCols();
	ClpQuadraticObjective * newObj = new ClpQuadraticObjective(*quadraticObj,
								   ncolsNow,
								   originalColumn_);
	// and modify linear and check
	double * linear = newObj->linearObjective();
	memcpy(linear,presolvedModel_->objective(),ncolsNow*sizeof(double));
	int iColumn;
	for ( iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (columnQuadraticLength[iColumn])
	    mark[iColumn]=1;
	}
	// and new
	quadratic = newObj->quadraticObjective();
	columnQuadraticLength = quadratic->getVectorLengths();
	int numberColumns2 = quadratic->getNumCols();
	for ( iColumn=0;iColumn<numberColumns2;iColumn++) {
	  if (columnQuadraticLength[iColumn])
	    mark[originalColumn_[iColumn]]=0;
	}
	presolvedModel_->setObjective(newObj);
	delete newObj;
	// final check
	for ( iColumn=0;iColumn<numberColumns;iColumn++) 
	  if (mark[iColumn])
	    printf("Quadratic column %d modified - may be okay\n",iColumn);
	delete [] mark;
      }
#endif
      delete [] prob.originalColumn_;
      prob.originalColumn_=NULL;
      int nrowsNow = presolvedModel_->getNumRows();
      memcpy(originalRow_,prob.originalRow_,nrowsNow*sizeof(int));
      delete [] prob.originalRow_;
      prob.originalRow_=NULL;
#ifndef CLP_NO_STD
      if (!dropNames&&originalModel->lengthNames()) {
	// Redo names
	int iRow;
	std::vector<std::string> rowNames;
	rowNames.reserve(nrowsNow);
	for (iRow=0;iRow<nrowsNow;iRow++) {
	  int kRow = originalRow_[iRow];
	  rowNames.push_back(originalModel->rowName(kRow));
	}
      
	int iColumn;
	std::vector<std::string> columnNames;
	columnNames.reserve(ncolsNow);
	for (iColumn=0;iColumn<ncolsNow;iColumn++) {
	  int kColumn = originalColumn_[iColumn];
	  columnNames.push_back(originalModel->columnName(kColumn));
	}
	presolvedModel_->copyNames(rowNames,columnNames);
      }
#endif
      if (rowObjective_) {
	int iRow;
        int k=-1;
        int nObj=0;
	for (iRow=0;iRow<nrowsNow;iRow++) {
	  int kRow = originalRow_[iRow];
          assert (kRow>k);
          k=kRow;
          rowObjective_[iRow]=rowObjective_[kRow];
          if (rowObjective_[iRow])
            nObj++;
	}
        if (nObj) {
          printf("%d costed slacks\n",nObj);
          presolvedModel_->setRowObjective(rowObjective_);
        }
      }
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
	    presolvedModel_->messageHandler()->message(COIN_PRESOLVE_COLINFEAS,
						       messages)
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
	  presolvedModel_->messageHandler()->message(COIN_PRESOLVE_INTEGERMODS,
						     messages)
						       <<numberChanges
						       <<CoinMessageEol;
	  if (!result) {
	    result = -1; // round again
	  }
	}
      }
    } else if (prob.status_) {
      // infeasible or unbounded
      result=1;
      originalModel->setProblemStatus(prob.status_);
    } else {
      // no changes - model needs restoring after Lou's changes
#ifndef CLP_NO_STD
      if (saveFile_=="") {
#endif
	delete presolvedModel_;
	presolvedModel_ = new ClpSimplex(*originalModel);
	// but we need to remove gaps
	ClpPackedMatrix* clpMatrix =
	  dynamic_cast< ClpPackedMatrix*>(presolvedModel_->clpMatrix());
	if (clpMatrix) {
	  clpMatrix->getPackedMatrix()->removeGaps();
	}
#ifndef CLP_NO_STD
      } else {
	presolvedModel_=originalModel;
      }
      presolvedModel_->dropNames();
#endif
      
      // drop integer information if wanted
      if (!keepIntegers)
	presolvedModel_->deleteIntegerInformation();
      result=2;
    }
  }
  if (result==0||result==2) {
    int nrowsAfter = presolvedModel_->getNumRows();
    int ncolsAfter = presolvedModel_->getNumCols();
    CoinBigIndex nelsAfter = presolvedModel_->getNumElements();
    presolvedModel_->messageHandler()->message(COIN_PRESOLVE_STATS,
					       messages)
						 <<nrowsAfter<< -(nrows_ - nrowsAfter)
						 << ncolsAfter<< -(ncols_ - ncolsAfter)
						 <<nelsAfter<< -(nelems_ - nelsAfter)
						 <<CoinMessageEol;
  } else {
    destroyPresolve();
    if (presolvedModel_!=originalModel_)
      delete presolvedModel_;
    presolvedModel_=NULL;
  }
  return presolvedModel_;
}


