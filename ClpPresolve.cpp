// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

//#define	CHECK_CONSISTENCY	1

#include <stdio.h>

#include <assert.h>
#include <iostream>

#include "CoinHelperFunctions.hpp"

#include "CoinPackedMatrix.hpp"
#include "ClpSimplex.hpp"

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
  paction_(0),
  ncols_(0),
  nrows_(0),
  nelems_(0),
  numberPasses_(5),
  saveFile_("")
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
}

/* This version of presolve returns a pointer to a new presolved 
   model.  NULL if infeasible
*/
ClpSimplex * 
ClpPresolve::presolvedModel(ClpSimplex & si,
			 double feasibilityTolerance,
			 bool keepIntegers,
			 int numberPasses,
			 bool dropNames)
{
  // Check matrix
  if (!si.clpMatrix()->allElementsInRange(&si,si.getSmallElementValue(),
					  1.0e20))
    return NULL;
  else
    return gutsOfPresolvedModel(&si,feasibilityTolerance,keepIntegers,numberPasses,dropNames);
}
/* This version of presolve updates
   model and saves original data to file.  Returns non-zero if infeasible
*/
int
ClpPresolve::presolvedModelToFile(ClpSimplex &si,std::string fileName,
			    double feasibilityTolerance,
			    bool keepIntegers,
			    int numberPasses)
{
  // Check matrix
  if (!si.clpMatrix()->allElementsInRange(&si,si.getSmallElementValue(),
					  1.0e20))
    return 2;
  saveFile_=fileName;
  si.saveModel(saveFile_.c_str());
  ClpSimplex * model = gutsOfPresolvedModel(&si,feasibilityTolerance,keepIntegers,numberPasses,true);
  if (model==&si) {
    return 0;
  } else {
    si.restoreModel(saveFile_.c_str());
    return 1;
  }
}

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
  // Messages
  CoinMessages messages = CoinMessage(presolvedModel_->messages().language());
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
  if (saveFile_=="") {
    // reality check
    assert(ncols0==originalModel_->getNumCols());
    assert(nrows0==originalModel_->getNumRows());
    acts = originalModel_->primalRowSolution();
    sol  = originalModel_->primalColumnSolution();
    if (updateStatus) {
      unsigned char *status = originalModel_->statusArray();
      rowstat = status + ncols0;
      colstat = status;
      memcpy(colstat, presolvedModel_->statusArray(), ncols);
      memcpy(rowstat, presolvedModel_->statusArray()+ncols, nrows);
    }
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


  CoinPostsolveMatrix prob(presolvedModel_,
		       ncols0,
		       nrows0,
		       nelems0,
		       presolvedModel_->getObjSense(),
		       // end prepost
		       
		       sol, acts,
		       colstat, rowstat);
    
  postsolve(prob);

  if (saveFile_!="") {
    // From file
    assert (originalModel_==presolvedModel_);
    originalModel_->restoreModel(saveFile_.c_str());
    memcpy(originalModel_->primalRowSolution(),acts,nrows0*sizeof(double));
    delete [] acts;
    memcpy(originalModel_->primalColumnSolution(),sol,ncols0*sizeof(double));
    delete [] sol;
    if (updateStatus) {
      memcpy(originalModel_->statusArray(),colstat,nrows0+ncols0);
      delete [] colstat;
    }
  }
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
      presolvedModel_->messageHandler()->message(COIN_PRESOLVE_NEEDS_CLEANING,
					    messages)
					      <<CoinMessageEol;
    }
  } else {
    originalModel_->setProblemStatus( presolvedModel_->status());
  }
  if (saveFile_!="") 
    presolvedModel_=NULL;
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
    // normal
#if 0
    const bool slackd = true;
    const bool doubleton = false;
    const bool tripleton = true;
    const bool forcing = false;
    const bool ifree = true;
    const bool zerocost = true;
    const bool dupcol = false;
    const bool duprow = false;
    const bool dual = doDualStuff;
#else
    const bool slackd = true;
    const bool doubleton = true;
    const bool tripleton = true;
    const bool forcing = true;
    const bool ifree = true;
    const bool zerocost = true;
    const bool dupcol = true;
    const bool duprow = true;
    const bool dual = doDualStuff;
#endif
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
    prob->pass_=0;
    for (iLoop=0;iLoop<numberPasses_;iLoop++) {
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
	  const CoinPresolveAction * const paction2 = paction_;
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
	
	prob->messageHandler()->message(COIN_PRESOLVE_PASS,
					messages)
					  <<numberDropped<<iLoop+1
					  <<CoinMessageEol;
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
					     CoinBigIndex nelems_in) :
  ncols_(si->getNumCols()),
  ncols0_(ncols_in),
  nelems_(si->getNumElements()),

  mcstrt_(new CoinBigIndex[ncols_in+1]),
  hincol_(new int[ncols_in+1]),
  hrow_  (new int   [2*nelems_in]),
  colels_(new double[2*nelems_in]),

  cost_(new double[ncols_in]),
  clo_(new double[ncols_in]),
  cup_(new double[ncols_in]),
  rlo_(new double[nrows_in]),
  rup_(new double[nrows_in]),
  originalColumn_(new int[ncols_in]),
  originalRow_(new int[nrows_in]),

  ztolzb_(getTolerance(si, ClpPrimalTolerance)),
  ztoldj_(getTolerance(si, ClpDualTolerance)),

  maxmin_(si->getObjSense())

{
  si->getDblParam(ClpObjOffset,originalOffset_);
  int ncols = si->getNumCols();
  int nrows = si->getNumRows();

  ClpDisjointCopyN(si->getColLower(), ncols, clo_);
  ClpDisjointCopyN(si->getColUpper(), ncols, cup_);
  ClpDisjointCopyN(si->getObjCoefficients(), ncols, cost_);
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
  int i;
  for (i = matrix.getSizeVectorLengths() - 1; i >= 0; --i) {
    if (start[i+1] - start[i] != length[i])
      break;
  }
  return (! (i >= 0));
}
#if	DEBUG_PRESOLVE
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
				     double maxmin_,
				     // end prepost members

				     ClpSimplex * si,

				     // rowrep
				     int nrows_in,
				     CoinBigIndex nelems_in,
			       bool doStatus,
			       double nonLinearValue) :

  CoinPrePostsolveMatrix(si,
			ncols0_in, nrows_in, nelems_in),
  clink_(new presolvehlink[ncols0_in+1]),
  rlink_(new presolvehlink[nrows_in+1]),

  dobias_(0.0),

  nrows_(si->getNumRows()),

  // temporary init
  mrstrt_(new CoinBigIndex[nrows_in+1]),
  hinrow_(new int[nrows_in+1]),
  rowels_(new double[2*nelems_in]),
  hcol_(new int[2*nelems_in]),
  integerType_(new char[ncols0_in]),
  feasibilityTolerance_(0.0),
  status_(-1),
  rowsToDo_(new int [nrows_in]),
  numberRowsToDo_(0),
  nextRowsToDo_(new int[nrows_in]),
  numberNextRowsToDo_(0),
  colsToDo_(new int [ncols0_in]),
  numberColsToDo_(0),
  nextColsToDo_(new int[ncols0_in]),
  numberNextColsToDo_(0)

{
  const int bufsize = 2*nelems_in;

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
  m = new CoinPackedMatrix();
  m->reverseOrderedCopyOf(*si->matrix());
  m->removeGaps();


  ClpDisjointCopyN(m->getVectorStarts(),  nrows_,  mrstrt_);
  mrstrt_[nrows_] = nelems_;
  ClpDisjointCopyN(m->getVectorLengths(), nrows_,  hinrow_);
  ClpDisjointCopyN(m->getIndices(),       nelems_, hcol_);
  ClpDisjointCopyN(m->getElements(),      nelems_, rowels_);

  delete m;
  if (si->integerInformation()) {
    memcpy(integerType_,si->integerInformation(),ncols_*sizeof(char));
  } else {
    ClpFillN<char>(integerType_, ncols_, 0);
  }

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

#if	DEBUG_PRESOLVE
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

#if	CHECK_CONSISTENCY
  consistent(false);
#endif
}

void CoinPresolveMatrix::update_model(ClpSimplex * si,
				     int nrows0,
				     int ncols0,
				     CoinBigIndex nelems0)
{
  si->loadProblem(ncols_, nrows_, mcstrt_, hrow_, colels_, hincol_,
		 clo_, cup_, cost_, rlo_, rup_);

  delete [] si->integerInformation();
  int numberIntegers=0;
  for (int i=0; i<ncols_; i++) {
    if (integerType_[i])
      numberIntegers++;
  }
  if (numberIntegers) 
    si->copyInIntegerInformation(integerType_);
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
				   
				       double maxmin_,
				       // end prepost members

				       double *sol_in,
				       double *acts_in,

				       unsigned char *colstat_in,
				       unsigned char *rowstat_in) :
  CoinPrePostsolveMatrix(si,
			ncols0_in, nrows0_in, nelems0),

  free_list_(0),
  // link, free_list, maxlink
  maxlink_(2*nelems0),
  link_(new int[/*maxlink*/ 2*nelems0]),
      
  cdone_(new char[ncols0_]),
  rdone_(new char[nrows0_in]),

  nrows_(si->getNumRows()),
  nrows0_(nrows0_in)
{

  sol_=sol_in;
  rowduals_=NULL;
  acts_=acts_in;

  rcosts_=NULL;
  colstat_=colstat_in;
  rowstat_=rowstat_in;

  // this is the *reduced* model, which is probably smaller
  int ncols1 = si->getNumCols();
  int nrows1 = si->getNumRows();

  const CoinPackedMatrix * m = si->matrix();

  if (! isGapFree(*m)) {
    CoinPresolveAction::throwCoinError("Matrix not gap free",
				      "CoinPostsolveMatrix");
  }

  const CoinBigIndex nelemsr = m->getNumElements();

  ClpDisjointCopyN(m->getVectorStarts(), ncols1, mcstrt_);
  CoinZeroN(mcstrt_+ncols1,ncols0_-ncols1);
  mcstrt_[ncols_] = nelems0;	// ??
  ClpDisjointCopyN(m->getVectorLengths(),ncols1,  hincol_);
  ClpDisjointCopyN(m->getIndices(),      nelemsr, hrow_);
  ClpDisjointCopyN(m->getElements(),     nelemsr, colels_);


#if	0 && DEBUG_PRESOLVE
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
  if (maxmin_<0.0) {
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
  }
  {
    int ml = maxlink_;
    for (CoinBigIndex k=nelemsr; k<ml; ++k)
      link_[k] = k+1;
    if (ml)
      link_[ml-1] = NO_LINK;
  }
  free_list_ = nelemsr;
}
/* This is main part of Presolve */
ClpSimplex * 
ClpPresolve::gutsOfPresolvedModel(ClpSimplex * originalModel,
				  double feasibilityTolerance,
				  bool keepIntegers,
				  int numberPasses,
				  bool dropNames)
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

  // result is 0 - okay, 1 infeasible, -1 go round again
  int result = -1;
  
  // User may have deleted - its their responsibility
  presolvedModel_=NULL;
  // Messages
  CoinMessages messages = CoinMessage(originalModel->messages().language());
  while (result==-1) {

    // make new copy
    if (saveFile_=="") {
      delete presolvedModel_;
      presolvedModel_ = new ClpSimplex(*originalModel);
    } else {
      presolvedModel_=originalModel;
    }
    presolvedModel_->dropNames();

    // drop integer information if wanted
    if (!keepIntegers)
      presolvedModel_->deleteIntegerInformation();

    
    CoinPresolveMatrix prob(ncols_,
			maxmin,
			presolvedModel_,
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
    prob.messages_ = CoinMessage(presolvedModel_->messages().language());

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
      delete [] prob.originalColumn_;
      prob.originalColumn_=NULL;
      int nrowsNow = presolvedModel_->getNumRows();
      memcpy(originalRow_,prob.originalRow_,nrowsNow*sizeof(int));
      delete [] prob.originalRow_;
      prob.originalRow_=NULL;
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
    presolvedModel_->messageHandler()->message(COIN_PRESOLVE_STATS,
					       messages)
						 <<nrowsAfter<< -(nrows_ - nrowsAfter)
						 << ncolsAfter<< -(ncols_ - ncolsAfter)
						 <<nelsAfter<< -(nelems_ - nelsAfter)
						 <<CoinMessageEol;
  } else {
    destroyPresolve();
    delete presolvedModel_;
    presolvedModel_=NULL;
  }
  return presolvedModel_;
}


