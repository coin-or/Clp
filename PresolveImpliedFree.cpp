#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "PresolveMatrix.hpp"
#include "PresolveSubst.hpp"
#include "PresolveIsolated.hpp"
#include "PresolveImpliedFree.hpp"
#include "ClpMessage.hpp"



// If there is a row with a singleton column such that no matter what
// the values of the other variables are, the constraint forces the singleton
// column to have a feasible value, then we can drop the column and row,
// since we just compute the value of the column from the row in postsolve.
// This seems vaguely similar to the case of a useless constraint, but it
// is different.  For example, if the singleton column is already free,
// then this operation will eliminate it, but the constraint is not useless
// (assuming the constraint is not trivial), since the variables do not imply an
// upper or lower bound.
//
// If the column is not singleton, we can still do something similar if the
// constraint is an equality constraint.  In that case, we substitute away
// the variable in the other constraints it appears in.  This introduces
// new coefficients, but the total number of coefficients never increases
// if the column has only two constraints, and may not increase much even
// if there are more.
//
// There is nothing to prevent us from substituting away a variable
// in an equality from the other constraints it appears in, but since
// that causes fill-in, it wouldn't make sense unless we could then
// drop the equality itself.  We can't do that if the bounds on the
// variable in equation aren't implied by the equality.
// Another way of thinking of this is that there is nothing special
// about an equality; just like one can't always drop an inequality constraint
// with a column singleton, one can't always drop an equality.
//
// It is possible for two singleton columns to be in the same row.
// In that case, the other one will become empty.  If it's bounds and
// costs aren't just right, this signals an unbounded problem.
// We don't need to check that specially here.
//
// invariant:  loosely packed
const PresolveAction *implied_free_action::presolve(PresolveMatrix *prob,
						     const PresolveAction *next,
						    int & fill_level)
{
  double *colels	= prob->colels_;
  int *hrow	= prob->hrow_;
  const CoinBigIndex *mcstrt	= prob->mcstrt_;
  int *hincol	= prob->hincol_;
  const int ncols	= prob->ncols_;

  const double *clo	= prob->clo_;
  const double *cup	= prob->cup_;

  const double *rowels	= prob->rowels_;
  const int *hcol	= prob->hcol_;
  const CoinBigIndex *mrstrt	= prob->mrstrt_;
  int *hinrow	= prob->hinrow_;
  const int nrows	= prob->nrows_;

  /*const*/ double *rlo	= prob->rlo_;
  /*const*/ double *rup	= prob->rup_;

  double *cost		= prob->cost_;

  presolvehlink *rlink = prob->rlink_;
  presolvehlink *clink = prob->clink_;

  const char *integerType = prob->integerType_;

  const double tol = prob->feasibilityTolerance_;

  int nbounds = 0;

  action *actions	= new action [ncols];
  int nactions = 0;

  char *implied_free = new char[ncols];
  bzero(implied_free, ncols*sizeof(char));

  double *ilbound = new double[ncols];
  double *iubound = new double[ncols];

  double *tclo = new double[ncols];
  double *tcup = new double[ncols];

#if 	PRESOLVE_TRY_TIGHTEN
  for (int j=0; j<ncols; j++) {
    ilbound[j] = -PRESOLVE_INF;
    iubound[j] =  PRESOLVE_INF;
    tclo[j] = clo[j];
    tcup[j] = cup[j];
  }

  {
    int ntightened;
    do {
      implied_bounds1(rowels, mrstrt, hrow, hinrow, clo, cup, hcol, ncols,
		      rlo, rup, integerType, nrows,
		      ilbound, iubound);

      ntightened = 0;
      for (int j=0; j<ncols; j++) {
	if (tclo[j] < ilbound[j]) {
	  tclo[j] = ilbound[j];
	  ntightened++;
	}
	if (tcup[j] > iubound[j]) {
	  tcup[j] = iubound[j];
	  ntightened++;
	}
      }
      printf("NTIGHT: %d\n", ntightened);
    } while (ntightened);
  }
#endif

  int numberLook = prob->numberColsToDo_;
  int iLook;
  int * look = prob->colsToDo_;
  int * look2 = NULL;
  // if gone from 2 to 3 look at all
  if (fill_level<0) {
    look2 = new int[ncols];
    look=look2;
    for (iLook=0;iLook<ncols;iLook++) 
      look[iLook]=iLook;
    numberLook=ncols;
 }

  for (iLook=0;iLook<numberLook;iLook++) {
    int j=look[iLook];
    if ((hincol[j] >= 1 && hincol[j] <= 3) &&
	!integerType[j]) {
      CoinBigIndex kcs = mcstrt[j];
      CoinBigIndex kce = kcs + hincol[j];

      for (CoinBigIndex k=kcs; k<kce; ++k) {
	int row = hrow[k];
	double coeffj = colels[k];

	// if its row is an equality constraint...
	if (hinrow[row] > 1 &&	// don't bother with singleton rows

	    // either this is a singleton col,
	    // or this particular row is an equality
	    (hincol[j] == 1 || fabs(rlo[row] - rup[row]) < tol) &&

	    fabs(coeffj) > ZTOLDP) {

	  CoinBigIndex krs = mrstrt[row];
	  CoinBigIndex kre = krs + hinrow[row];

	  double maxup, maxdown, ilow, iup;
	  implied_bounds(rowels, clo, cup, hcol,
			 krs, kre,
			 &maxup, &maxdown,
			 j, rlo[row], rup[row], &ilow, &iup);

#if	PRESOLVE_TRY_TIGHTEN
	  if ((clo[j] <= ilbound[j] && iubound[j] <= cup[j]) &&
	      ! (clo[j] <= ilow && iup <= cup[j]))
	    printf("TIGHTER:  %6d %6d\n", row, j);
#endif

	  if (maxup < PRESOLVE_INF && maxup + tol < rlo[row]) {
	    /* there is an upper bound and it can't be reached */
	    prob->status_|= 1;
	    prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
					     prob->messages())
					       <<row
					       <<rlo[row]
					       <<rup[row]
					       <<CoinMessageEol;
	    break;
	  } else if (-PRESOLVE_INF < maxdown && rup[row] < maxdown - tol) {
	    /* there is a lower bound and it can't be reached */
	    prob->status_|= 1;
	    prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
					     prob->messages())
					       <<row
					       <<rlo[row]
					       <<rup[row]
					       <<CoinMessageEol;
	    break;
	  } else if (clo[j] <= ilow && iup <= cup[j]) {

	    // both column bounds implied by the constraints of the problem
	    implied_free[j] = hincol[j];
	    break;
	  }
	}
      }
    }
  }
  // implied_free[j] == hincol[j] && hincol[j] > 0 ==> j is implied free

#if 0
  // DEBUG
  static int nfree = 0;
  static int maxfree = atoi(getenv("MAXFREE"));
#endif

  int isolated_row = -1;

  // first pick off the easy ones
  // note that this will only deal with columns that were originally
  // singleton; it will not deal with doubleton columns that become
  // singletons as a result of dropping rows.
  for (iLook=0;iLook<numberLook;iLook++) {
    int j=look[iLook];
    if (hincol[j] == 1 && implied_free[j] == 1) {
      CoinBigIndex kcs = mcstrt[j];
      int row = hrow[kcs];
      double coeffj = colels[kcs];

      CoinBigIndex krs = mrstrt[row];
      CoinBigIndex kre = krs + hinrow[row];

#if 0
      if (nfree >= maxfree)
	continue;
      nfree++;
#endif

      // isolated rows are weird
      {
	int n = 0;
	for (CoinBigIndex k=krs; k<kre; ++k)
	  n += hincol[hcol[k]];
	if (n==hinrow[row]) {
	  isolated_row = row;
	  break;
	}
      }

      const bool nonzero_cost = (cost[j] != 0.0&&fabs(rup[row]-rlo[row])<=tol);

      double *save_costs = nonzero_cost ? new double[hinrow[row]] : NULL;

      {
	action *s = &actions[nactions++];
	      
	s->row = row;
	s->col = j;

	s->clo = clo[j];
	s->cup = cup[j];
	s->rlo = rlo[row];
	s->rup = rup[row];

	s->ninrow = hinrow[row];
	s->rowels = presolve_duparray(&rowels[krs], hinrow[row]);
	s->rowcols = presolve_duparray(&hcol[krs], hinrow[row]);
	s->costs = save_costs;
      }

      if (nonzero_cost) {
	double rhs = rlo[row];
	double costj = cost[j];

#if	DEBUG_PRESOLVE
	printf("FREE COSTS:  %g  ", costj);
#endif
	for (CoinBigIndex k=krs; k<kre; k++) {
	  int jcol = hcol[k];
	  save_costs[k-krs] = cost[jcol];

	  if (jcol != j) {
	    double coeff = rowels[k];

#if	DEBUG_PRESOLVE
	    printf("%g %g   ", cost[jcol], coeff/coeffj);
#endif
	    /*
	     * Similar to eliminating doubleton:
	     *   cost1 x = cost1 (c - b y) / a = (c cost1)/a - (b cost1)/a
	     *   cost[icoly] += cost[icolx] * (-coeff2 / coeff1);
	     */
	    cost[jcol] += costj * (-coeff / coeffj);
	  }
	}
#if	DEBUG_PRESOLVE
	printf("\n");

	/* similar to doubleton */
	printf("BIAS??????? %g %g %g %g\n",
	       costj * rhs / coeffj,
	       costj, rhs, coeffj);
#endif
	prob->change_bias(costj * rhs / coeffj);
	// ??
	cost[j] = 0.0;
      }

      /* remove the row from the columns in the row */
      for (CoinBigIndex k=krs; k<kre; k++) {
	int jcol=hcol[k];
	prob->addCol(jcol);
	presolve_delete_from_row(jcol, row, mcstrt, hincol, hrow, colels);
      }
      PRESOLVE_REMOVE_LINK(rlink, row);
      hinrow[row] = 0;

      // just to make things squeeky
      rlo[row] = 0.0;
      rup[row] = 0.0;

      PRESOLVE_REMOVE_LINK(clink, j);
      hincol[j] = 0;

      implied_free[j] = 0;	// probably unnecessary
    }
  }

  delete [] look2;
  if (nactions) {
#if	PRESOLVE_SUMMARY
    printf("NIMPLIED FREE:  %d\n", nactions);
#endif
    next = new implied_free_action(nactions, copyOfArray(actions,nactions), next);
  }
  delete [] actions;

  delete[]ilbound;
  delete[]iubound;
  delete[]tclo;
  delete[]tcup;

  if (isolated_row != -1) {
    const PresolveAction *nextX = isolated_constraint_action::presolve(prob, 
						isolated_row, next);
    if (nextX)
      next = nextX; // may fail
  }

  // try more complex ones
  if (fill_level)
    next = subst_constraint_action::presolve(prob, implied_free, next,fill_level);
  delete[]implied_free;

  return (next);
}



const char *implied_free_action::name() const
{
  return ("implied_free_action");
}

void implied_free_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions = nactions_;

  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double *sol	= prob->sol_;

  double *rcosts	= prob->rcosts_;
  double *dcost		= prob->cost_;

  double *acts	= prob->acts_;
  double *rowduals = prob->rowduals_;

  const double ztoldj	= prob->ztoldj_;
  const double ztolzb	= prob->ztolzb_;

  const double maxmin	= prob->maxmin_;

  char *cdone	= prob->cdone_;
  char *rdone	= prob->rdone_;
  CoinBigIndex free_list = prob->free_list_;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {

    int irow = f->row;
    int icol = f->col;
	  
    int ninrow = f->ninrow;
    const double *rowels = f->rowels;
    const int *rowcols = f->rowcols;
    const double *save_costs = f->costs;

    // put back coefficients in the row
    // this includes recreating the singleton column
    {
      for (int k = 0; k<ninrow; k++) {
	int jcol = rowcols[k];
	double coeff = rowels[k];

	if (save_costs)
	  dcost[jcol] = save_costs[k];

	{
	  CoinBigIndex kk = free_list;
	  free_list = link[free_list];

	  check_free_list(free_list);

	  link[kk] = mcstrt[jcol];
	  mcstrt[jcol] = kk;
	  colels[kk] = coeff;
	  hrow[kk] = irow;
	}

	if (jcol == icol) {
	  // initialize the singleton column
	  hincol[jcol] = 1;
	  clo[icol] = f->clo;
	  cup[icol] = f->cup;

	  cdone[icol] = IMPLIED_FREE;
	} else {
	  hincol[jcol]++;
	}
      }
      rdone[irow] = IMPLIED_FREE;

      rlo[irow] = f->rlo;
      rup[irow] = f->rup;
    }
    delete [] save_costs;
    // coeff has now been initialized

    // compute solution
    {
      double act = 0.0;
      double coeff;

      for (int k = 0; k<ninrow; k++)
	if (rowcols[k] == icol)
	  coeff = rowels[k];
	else {
	  int jcol = rowcols[k];
	  PRESOLVE_STMT(CoinBigIndex kk = presolve_find_row2(irow, mcstrt[jcol], hincol[jcol], hrow, link));
	  PRESOLVEASSERT(colels[kk] == rowels[k]);
	  act += rowels[k] * sol[jcol];
	}
	    
      PRESOLVEASSERT(fabs(coeff) > ZTOLDP);
      // choose rowdual to make this col basic
      rowduals[irow] = maxmin*dcost[icol] / coeff;
      rcosts[icol] = 0.0;
      if ((rlo[irow] < rup[irow] && rowduals[irow] < 0.0)
	  || rlo[irow]< -1.0e20) {
	if (rlo[irow]<-1.0e20&&rowduals[irow]>=0.0)
	  printf("IMP %g %g %g\n",rlo[irow],rup[irow],rowduals[irow]);
	sol[icol] = (rup[irow] - act) / coeff;
	acts[irow] = rup[irow];
      } else {
	sol[icol] = (rlo[irow] - act) / coeff;
	acts[irow] = rlo[irow];
      }


      prob->setRowStatus(irow,PrePostsolveMatrix::atLowerBound);
      prob->setColumnStatus(icol,PrePostsolveMatrix::basic);
    }
  }
  prob->free_list_ = free_list;
}
