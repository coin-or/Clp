// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <stdio.h>
#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "PresolveMatrix.hpp"

#include "PresolveEmpty.hpp"	// for DROP_COL/DROP_ROW
#include "PresolveFixed.hpp"
#include "PresolveSingleton.hpp"
#include "ClpMessage.hpp"

/*
 * Transfers singleton row bound information to the corresponding column bounds.
 * What I refer to as a row singleton would be called a doubleton
 * in the paper, since my terminology doesn't refer to the slacks.
 * In terms of the paper, we transfer the bounds of the slack onto
 * the variable (vii) and then "substitute" the slack out of the problem 
 * (which is a noop).
 */
const PresolveAction *
slack_doubleton_action::presolve(PresolveMatrix *prob,
				 const PresolveAction *next,
				 bool & notFinished)
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int ncols		= prob->ncols_;

  double *clo		= prob->clo_;
  double *cup		= prob->cup_;

  double *rowels	= prob->rowels_;
  const int *hcol	= prob->hcol_;
  const CoinBigIndex *mrstrt	= prob->mrstrt_;
  int *hinrow		= prob->hinrow_;
  //  int nrows		= prob->nrows_;

  double *rlo		= prob->rlo_;
  double *rup		= prob->rup_;

  // If rowstat exists then all do
  unsigned char *rowstat	= prob->rowstat_;
  //  double *acts	= prob->acts_;
  double * sol = prob->sol_;
  //  unsigned char * colstat = prob->colstat_;

  //  const char *integerType = prob->integerType_;

  const double ztolzb	= prob->ztolzb_;

  int numberLook = prob->numberRowsToDo_;
  int iLook;
  int * look = prob->rowsToDo_;

  action actions[MAX_SLACK_DOUBLETONS];
  int nactions = 0;
  notFinished = false;

  int * fixed_cols = new int [ncols];
  int nfixed_cols	= 0;

  // this loop can apparently create new singleton rows;
  // I haven't bothered to detect this.
  // wasfor (int irow=0; irow<nrows; irow++)
  for (iLook=0;iLook<numberLook;iLook++) {
    int irow = look[iLook];
    if (hinrow[irow] == 1) {
      int jcol = hcol[mrstrt[irow]];
      double coeff = rowels[mrstrt[irow]];
      double lo = rlo[irow];
      double up = rup[irow];
      double acoeff = fabs(coeff);
      //      const bool singleton_col = (hincol[jcol] == 1);

      if (acoeff < ZTOLDP)
	continue;

      // don't bother with fixed cols
      if (fabs(cup[jcol] - clo[jcol]) < ztolzb)
	continue;

      {
	// put column on stack of things to do next time
	prob->addCol(jcol);
	action *s = &actions[nactions];
	nactions++;

	s->col = jcol;
	s->clo = clo[jcol];
	s->cup = cup[jcol];

	s->row = irow;
	s->rlo = rlo[irow];
	s->rup = rup[irow];

	s->coeff = coeff;
      }

      if (coeff < 0.0) {
	swap(lo, up);
	lo = -lo;
	up = -up;
      }
      
      if (lo <= -PRESOLVE_INF)
	lo = -PRESOLVE_INF;
      else {
	lo /= acoeff;
	if (lo <= -PRESOLVE_INF)
	  lo = -PRESOLVE_INF;
      }

      if (up > PRESOLVE_INF)
	up = PRESOLVE_INF;
      else {
	up /= acoeff;
	if (up > PRESOLVE_INF)
	  up = PRESOLVE_INF;
      }
      
      if (clo[jcol] < lo)
	clo[jcol] = lo;
	
      if (cup[jcol] > up) 
	cup[jcol] = up;

      if (fabs(cup[jcol] - clo[jcol]) < ZTOLDP) {
	fixed_cols[nfixed_cols++] = jcol;
      }

      if (lo > up) {
	if (lo <= up + prob->feasibilityTolerance_) {
	  // If close to integer then go there
	  double nearest = floor(lo+0.5);
	  if (fabs(nearest-lo)<2.0*prob->feasibilityTolerance_) {
	    lo = nearest;
	    up = nearest;
	  } else {
	    lo = up;
	  }
	  clo[jcol] = lo;
	  cup[jcol] = up;
	} else {
	  prob->status_ |= 1;
	  prob->messageHandler()->message(CLP_PRESOLVE_COLINFEAS,
					     prob->messages())
					       <<jcol
					       <<lo
					       <<up
					       <<CoinMessageEol;
	  break;
	}
      }

#if	DEBUG_PRESOLVE
      printf("SINGLETON R-%d C-%d\n", irow, jcol);
#endif

      // eliminate the row entirely from the row rep
      // rlinks don't seem to be used PRESOLVE_REMOVE_LINK(rlink, irow);
      hinrow[irow] = 0;

      // just to make things squeeky
      rlo[irow] = 0.0;
      rup[irow] = 0.0;

      if (rowstat) {
	// update solution and basis
	int basisChoice=0;
	int numberBasic=0;
	if (prob->columnIsBasic(jcol))
	  numberBasic++;
	if (prob->rowIsBasic(irow))
	  numberBasic++;
	if (sol[jcol]<=clo[jcol]+ztolzb) {
	  sol[jcol] =clo[jcol];
	  prob->setColumnStatus(jcol,PrePostsolveMatrix::atLowerBound);
	} else if (sol[jcol]>=cup[jcol]-ztolzb) {
	  sol[jcol] =cup[jcol];
	  prob->setColumnStatus(jcol,PrePostsolveMatrix::atUpperBound);
	} else {
	  basisChoice=1;
	}
	if (numberBasic>1||basisChoice)
	  prob->setColumnStatus(jcol,PrePostsolveMatrix::basic);
      }

      // remove the row from this col in the col rep
      presolve_delete_from_row(jcol, irow, mcstrt, hincol, hrow, colels);

      if (nactions >= MAX_SLACK_DOUBLETONS) {
	notFinished=true;
	break;
      }
    }
  }

  if (nactions) {
#if	PRESOLVE_SUMMARY
    printf("SINGLETON ROWS:  %d\n", nactions);
#endif
    action *save_actions = new action[nactions];
    ClpDisjointCopyN(actions, nactions, save_actions);
    next = new slack_doubleton_action(nactions, save_actions, next);

    if (nfixed_cols)
      next = make_fixed_action::presolve(prob, fixed_cols, nfixed_cols,
					 true, // arbitrary
					 next);
  }
  delete [] fixed_cols;
  return (next);
}

void slack_doubleton_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions = nactions_;

  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;
  //  int ncols		= prob->ncols_;

  double *clo		= prob->clo_;
  double *cup		= prob->cup_;

  double *rlo		= prob->rlo_;
  double *rup		= prob->rup_;

  double *sol		= prob->sol_;
  double *rcosts	= prob->rcosts_;

  double *acts		= prob->acts_;
  double *rowduals 	= prob->rowduals_;

  unsigned char *colstat		= prob->colstat_;
  //  unsigned char *rowstat		= prob->rowstat_;

  //  char *cdone		= prob->cdone_;
  char *rdone		= prob->rdone_;
  CoinBigIndex free_list		= prob->free_list_;

  const double ztolzb	= prob->ztolzb_;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {
    int irow = f->row;
    double lo0 = f->clo;
    double up0 = f->cup;
    double coeff = f->coeff;
    int jcol = f->col;

    rlo[irow] = f->rlo;
    rup[irow] = f->rup;

    clo[jcol] = lo0;
    cup[jcol] = up0;

    acts[irow] = coeff * sol[jcol];

    // copy col to end to make room for new element
    {
      CoinBigIndex k = free_list;
      free_list = link[free_list];

      check_free_list(free_list);

      hrow[k] = irow;
      colels[k] = coeff;
      link[k] = mcstrt[jcol];
      mcstrt[jcol] = k;
    }
    hincol[jcol]++;	// right?

    /*
     * Have to compute rowstat and rowduals, since we are adding the row.
     * that satisfy complentarity slackness.
     *
     * We may also have to modify colstat and rcosts if bounds
     * have been relaxed.
     */
    if (!colstat) {
      // ????
      rowduals[irow] = 0.0;
    } else {
      if (prob->columnIsBasic(jcol)) {
	/* variable is already basic, so slack in this row must be basic */
	prob->setRowStatus(irow,PrePostsolveMatrix::basic);
	
	rowduals[irow] = 0.0;
	/* hence no reduced costs change */
	/* since the column was basic, it doesn't matter if it is now
	   away from its bounds. */
	/* the slack is basic and its reduced cost is 0 */
      } else if ((fabs(sol[jcol] - lo0) <= ztolzb &&
		  rcosts[jcol] >= 0) ||
		 
		 (fabs(sol[jcol] - up0) <= ztolzb &&
		  rcosts[jcol] <= 0)) {
	/* up against its bound but the rcost is still ok */
	prob->setRowStatus(irow,PrePostsolveMatrix::basic);
	
	rowduals[irow] = 0.0;
	/* hence no reduced costs change */
      } else if (! (fabs(sol[jcol] - lo0) <= ztolzb) &&
		 ! (fabs(sol[jcol] - up0) <= ztolzb)) {
	/* variable not marked basic,
	 * so was at a bound in the reduced problem,
	 * but its bounds were tighter in the reduced problem,
	 * so now it has to be made basic.
	 */
	prob->setColumnStatus(jcol,PrePostsolveMatrix::basic);
	prob->setRowStatusUsingValue(irow);

	/* choose a value for rowduals[irow] that forces rcosts[jcol] to 0.0 */
	/* new reduced cost = 0 = old reduced cost - new dual * coeff */
	rowduals[irow] = rcosts[jcol] / coeff;
	rcosts[jcol] = 0.0;

	/* this value is provably of the right sign for the slack */
	/* SHOULD SHOW THIS */
      } else {
	/* the solution is at a bound, but the rcost is wrong */

	prob->setColumnStatus(jcol,PrePostsolveMatrix::basic);
	prob->setRowStatusUsingValue(irow);

	/* choose a value for rowduals[irow] that forcesd rcosts[jcol] to 0.0 */
	rowduals[irow] = rcosts[jcol] / coeff;
	rcosts[jcol] = 0.0;

	/* this value is provably of the right sign for the slack */
	/*rcosts[irow] = -rowduals[irow];*/
      }
    }

    rdone[irow] = SLACK_DOUBLETON;
  }
  prob->free_list_ = free_list;
}

