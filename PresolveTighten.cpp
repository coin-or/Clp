#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "PresolveMatrix.hpp"
#include "PresolveFixed.hpp"
#include "PresolveTighten.hpp"
#include "PresolveUseless.hpp"

const char *do_tighten_action::name() const
{
  return ("do_tighten_action");
}

// This is ekkredc2.
// This fairly simple transformation is not mentioned in the paper.
// Say there is a costless variable such all its constraints
// would be satisfied as it approaches plus or minus infinity,
// because all its constraints have only one bound, and increasing/decreasing
// the variable makes the row activity grow away from the bound
// (in the right direction).
//
// If the variable is unbounded in that direction,
// that means we can determine right now how large it needs
// to get in order to satisfy the constraints, so we can
// just drop the variable and those constraints from the problem.
//
// If the variable *is* bounded in that direction,
// there is no reason not to set it to that bound.
// This effectively weakens the constraints, and in fact 
// may be subsequently presolved away.
//
// Note that none of the constraints may be bounded both above and below,
// since then we don't know which way to move the variable in order
// to satisfy the constraint.
//
// To drop constraints, we just make them useless and let other
// transformations take care of the rest.
//
// Note that more than one such costless unbounded variable
// may be part of a given constraint.
// In that case, the first one processed will make the
// constraint useless, and the second will ignore it.
// In postsolve, the first will be responsible for satisfying
// the constraint.
//
// Note that if the constraints are dropped (as in the first case),
// then we just make them useless.  It is subsequently discovered
// the the variable does not appear in any constraints, and since it
// has no cost it is just set to some value (either zero or a bound)
// and removed (by remove_empty_cols).
//
// oddly, pilots and baxter do *worse* when this transform is applied.
const PresolveAction *do_tighten_action::presolve(PresolveMatrix *prob,
						   const PresolveAction *next)
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int ncols		= prob->ncols_;

  int nrows		= prob->nrows_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double *dcost	= prob->cost_;

  // NEED TO THINK MORE ABOUT INTEGER VARS
  const char *integerType = prob->integerType_;

  int *fixup_cols	= new int[ncols];
  int nfixup_cols	= 0;

  int *fixdown_cols	= new int[ncols];
  int nfixdown_cols	= 0;

  int *useless_rows	= new int[nrows];
  int nuseless_rows	= 0;
  
  action *actions	= new action [ncols];
  int nactions		= 0;

  int numberLook = prob->numberColsToDo_;
  int iLook;
  int * look = prob->colsToDo_;

  // singleton columns are especially likely to be caught here
  for (iLook=0;iLook<numberLook;iLook++) {
    int j = look[iLook];
    if (dcost[j]==0.0) {
      int iflag=0; /* 1 - up is towards feasibility, -1 down is towards */

      CoinBigIndex kcs = mcstrt[j];
      CoinBigIndex kce = kcs + hincol[j];

      // check constraints
      for (CoinBigIndex k=kcs; k<kce; ++k) {
	int i = hrow[k];
	double coeff = colels[k];
	double rlb = rlo[i];
	double rub = rup[i];

	if (-1.0e28 < rlb && rub < 1.0e28) {
	  // bounded - we lose
	  iflag=0;
	  break;
	}

	PRESOLVEASSERT(fabs(coeff) > ZTOLDP);

	// see what this particular row says
	// jflag == 1 ==> up is towards feasibility
	int jflag = (coeff > 0.0
		     ? (rub >  1.0e28 ? 1 : -1)
		     : (rlb < -1.0e28 ? 1 : -1));

	if (iflag) {
	  // check that it agrees with iflag.
	  if (iflag!=jflag) {
	    iflag=0;
	    break;
	  }
	} else {
	  // first row -- initialize iflag
	  iflag=jflag;
	}
      }
      // done checking constraints

      if (iflag) {
	if (iflag==1 && cup[j]<1.0e10) {
#if	DEBUG_PRESOLVE
	  printf("TIGHTEN UP:  %d\n", j);
#endif
	  fixup_cols[nfixup_cols] = j;
	  ++nfixup_cols;

	} else if (iflag==-1&&clo[j]>-1.0e10) {
	  // symmetric case
	  //mpre[j] = PRESOLVE_XUP;

#if	DEBUG_PRESOLVE
	  printf("TIGHTEN DOWN:  %d\n", j);
#endif

	  fixdown_cols[nfixdown_cols] = j;
	  ++nfixdown_cols;

	} else {
#if 0
	  static int limit;
	  static int which = atoi(getenv("WZ"));
	  if (which == -1)
	    ;
	  else if (limit != which) {
	    limit++;
	    continue;
	  } else
	    limit++;

	  printf("TIGHTEN STATS %d %g %g %d:  \n", j, clo[j], cup[j], integerType[j]); 
  double *rowels	= prob->rowels_;
  int *hcol		= prob->hcol_;
  int *mrstrt		= prob->mrstrt_;
  int *hinrow		= prob->hinrow_;
	  for (CoinBigIndex k=kcs; k<kce; ++k) {
	    int irow = hrow[k];
	    CoinBigIndex krs = mrstrt[irow];
	    CoinBigIndex kre = krs + hinrow[irow];
	    printf("%d  %g %g %g:  ",
		   irow, rlo[irow], rup[irow], colels[irow]);
	    for (CoinBigIndex kk=krs; kk<kre; ++kk)
	      printf("%d(%g) ", hcol[kk], rowels[kk]);
	    printf("\n");
	  }
#endif

	  {
	    action *s = &actions[nactions];	  
	    nactions++;
	    s->col = j;
	    s->direction = iflag;

	    s->rows =   new int[hincol[j]];
	    s->lbound = new double[hincol[j]];
	    s->ubound = new double[hincol[j]];
#ifdef DEBUG_PRESOLVE
	    printf("TIGHTEN FREE:  %d   ", j);
#endif
	    int nr = 0;
            prob->addCol(j);
	    for (CoinBigIndex k=kcs; k<kce; ++k) {
	      int irow = hrow[k];
	      // ignore this if we've already made it useless
	      if (! (rlo[irow] == -PRESOLVE_INF && rup[irow] == PRESOLVE_INF)) {
		prob->addRow(irow);
		s->rows  [nr] = irow;
		s->lbound[nr] = rlo[irow];
		s->ubound[nr] = rup[irow];
		nr++;

		useless_rows[nuseless_rows++] = irow;

		rlo[irow] = -PRESOLVE_INF;
		rup[irow] = PRESOLVE_INF;

#ifdef DEBUG_PRESOLVE
		printf("%d ", irow);
#endif
	      }
	    }
	    s->nrows = nr;

#ifdef DEBUG_PRESOLVE
	    printf("\n");
#endif
	  }
	}
      }
    }
  }


#if	PRESOLVE_SUMMARY
  if (nfixdown_cols || nfixup_cols || nuseless_rows) {
    printf("NTIGHTENED:  %d %d %d\n", nfixdown_cols, nfixup_cols, nuseless_rows);
  }
#endif

  if (nuseless_rows) {
    next = new do_tighten_action(nactions, copyOfArray(actions,nactions), next);

    next = useless_constraint_action::presolve(prob,
					       useless_rows, nuseless_rows,
					       next);
  }
  delete [] actions;
  delete[]useless_rows;

  if (nfixdown_cols) {
    next = make_fixed_action::presolve(prob, fixdown_cols, nfixdown_cols,
				       true,
				       next);
  }
  delete[]fixdown_cols;

  if (nfixup_cols) {
    next = make_fixed_action::presolve(prob, fixup_cols, nfixup_cols,
				       false,
				       next);
  }
  delete[]fixup_cols;

  return (next);
}

void do_tighten_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions	= nactions_;

  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;
  int ncols		= prob->ncols_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;
  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double *sol	= prob->sol_;
  double *dcost	= prob->cost_;
  double *rcosts	= prob->rcosts_;

  double *acts	= prob->acts_;
  double *rowduals = prob->rowduals_;


  const double ztolzb	= prob->ztolzb_;

  char *cdone	= prob->cdone_;
  char *rdone	= prob->rdone_;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {
    int jcol = f->col;
    int iflag = f->direction;
    int nr   = f->nrows;
    const int *rows = f->rows;
    const double *lbound = f->lbound;
    const double *ubound = f->ubound;

    PRESOLVEASSERT(prob->getColumnStatus(jcol)!=PrePostsolveMatrix::basic);
    for (int i=0;i<nr; ++i) {
      int irow = rows[i];

      rlo[irow] = lbound[i];
      rup[irow] = ubound[i];

     PRESOLVEASSERT(prob->getRowStatus(irow)==PrePostsolveMatrix::basic);
    }

    // We have just tightened the row bounds.
    // That means we'll have to compute a new value
    // for this variable that will satisfy everybody.
    // We are supposed to be in a position where this
    // is always possible.

    // Each constraint has exactly one bound.
    // The correction should only ever be forced to move in one direction.
    double orig_sol = sol[jcol];
    double correction = 0.0;
    
    int last_corrected = -1;
    CoinBigIndex k = mcstrt[jcol];
    int nk = hincol[jcol];
    for (int i=0; i<nk; ++i) {
      int irow = hrow[k];
      double coeff = colels[k];
      k = link[k];
      double newrlo = rlo[irow];
      double newrup = rup[irow];
      double activity = acts[irow];

      if (activity + correction * coeff < newrlo) {
	// only one of these two should fire
	PRESOLVEASSERT( ! (activity + correction * coeff > newrup) );

	last_corrected = irow;

	// adjust to just meet newrlo (solve for correction)
	double new_correction = (newrlo - activity) / coeff;
	PRESOLVEASSERT((iflag == 1) == (correction < new_correction));
	correction = new_correction;
      } else if (activity + correction * coeff > newrup) {
	last_corrected = irow;

	double new_correction = (newrup - activity) / coeff;
	PRESOLVEASSERT((iflag == 1) == (correction < new_correction));
	correction = new_correction;
      }
    }

    sol[jcol] += correction;
    PRESOLVEASSERT(clo[jcol] - ztolzb < sol[jcol] &&
	      sol[jcol] <= cup[jcol] + ztolzb);

    // by construction, the last row corrected (if there was one)
    // must be at its bound, so it can be non-basic.
    // All other rows may not be at a bound (but may if the difference
    // is very small, causing a new correction by a tiny amount).

    // now adjust the activities
    k = mcstrt[jcol];
    for (int i=0; i<nk; ++i) {
      int irow = hrow[k];
      double coeff = colels[k];
      k = link[k];
      double activity = acts[irow];

      acts[irow] += correction * coeff;

      // ? there isn't really a tolerance for tihs
      PRESOLVEASSERT(rlo[irow] - ztolzb <= acts[irow] &&
		acts[irow] <= rup[irow] + ztolzb);
    }

    // if the col happens to get pushed to its bound,
    // we may as well leave it non-basic.
    if (fabs(sol[jcol] - clo[jcol]) > ZTOLDP &&
	fabs(sol[jcol] - cup[jcol]) > ZTOLDP) {

      PRESOLVEASSERT(last_corrected != -1);
      prob->setRowStatus(last_corrected,PrePostsolveMatrix::atLowerBound);
      prob->setColumnStatus(jcol,PrePostsolveMatrix::basic);
    }
  }
}
