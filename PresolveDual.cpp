#include <stdio.h>
#include <math.h>

#include "PresolveMatrix.hpp"
#include "PresolveFixed.hpp"
#include "PresolveDual.hpp"
#include "ClpMessage.hpp"

// this looks for "dominated columns"
// following ekkredc
// rdmin == dwrow2
// rdmax == dwrow
// this transformation is applied in: k4.mps, lp22.mps

// make inferences using this constraint on reduced cost:
//	at optimality	dj>0 ==> var must be at lb
//			dj<0 ==> var must be at ub
//
// This implies:
//	there is no lb ==> dj<=0 at optimality
//	there is no ub ==> dj>=0 at optimality
const PresolveAction *remove_dual_action::presolve(PresolveMatrix *prob,
				   const PresolveAction *next)
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int ncols		= prob->ncols_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  //  double *rowels	= prob->rowels_;
  //  int *hcol		= prob->hcol_;
  //  CoinBigIndex *mrstrt		= prob->mrstrt_;
  //  int *hinrow		= prob->hinrow_;
  int nrows		= prob->nrows_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double *dcost	= prob->cost_;

  const double maxmin	= prob->maxmin_;

  const double ekkinf = 1e28;
  const double ekkinf2 = 1e20;

  double *rdmin	= new double[nrows];
  double *rdmax	= new double[nrows];

  // This combines the initialization of rdmin/rdmax to extreme values
  // (PRESOLVE_INF/-PRESOLVE_INF) with a version of the next loop specialized
  // for row slacks.
  // In this case, it is always the case that dprice==0.0 and coeff==1.0.
  for (int i = 0; i < nrows; i++) {
    double sup = -rlo[i];	// slack ub; corresponds to cup[j]
    double slo = -rup[i];	// slack lb; corresponds to clo[j]
    bool no_lb = (slo <= -ekkinf);
    bool no_ub = (sup >= ekkinf);

    // here, dj==-row dual
    // the row slack has no lower bound, but it does have an upper bound,
    // then the reduced cost must be <= 0, so the row dual must be >= 0
    rdmin[i] = (no_lb && !no_ub) ? 0.0 : -PRESOLVE_INF;

    rdmax[i] = (no_ub && !no_lb) ? 0.0 :  PRESOLVE_INF;
  }

  // Look for col singletons and update bounds on dual costs
  // Take the min of the maxes and max of the mins
  for (int j = 0; j<ncols; j++) {
    bool no_ub = (cup[j] >= ekkinf);
    bool no_lb = (clo[j] <= -ekkinf);
    
    if (hincol[j] == 1 &&

	// we need singleton cols that have exactly one bound
	(no_ub ^ no_lb)) {
      int row = hrow[mcstrt[j]];
      double coeff = colels[mcstrt[j]];

      PRESOLVEASSERT(fabs(coeff) > ZTOLDP);

      // I don't think the sign of dcost[j] matters

      // this row dual would make this col's reduced cost be 0
      double dprice = maxmin * dcost[j] / coeff;

      // no_ub == !no_lb
      // no_ub ==> !(dj<0)
      // no_lb ==> !(dj>0)
      // I don't think the case where !no_ub has actually been tested
      if ((coeff > 0.0) == no_ub) {
	// coeff>0 ==> making the row dual larger would make dj *negative*
	//		==> dprice is an upper bound on dj if no *ub*
	// coeff<0 ==> making the row dual larger would make dj *positive*
	//		==> dprice is an upper bound on dj if no *lb*
	if (rdmax[row] > dprice)	// reduced cost may be positive
	  rdmax[row] = dprice;
      } else {			// no lower bound
	if (rdmin[row] < dprice) 	// reduced cost may be negative
	  rdmin[row] = dprice;
      }
    }
  }

  int *fixup_cols	= new int[ncols];

  int *fixdown_cols	= new int[ncols];

#if	PRESOLVE_TIGHTEN_DUALS
  double *djmin	= new double[ncols];
  double *djmax	= new double[ncols];
#endif
  int nfixup_cols	= 0;
  int nfixdown_cols	= 0;

  while (1) {
    /* Perform duality tests */
    for (int j = 0; j<ncols; j++)
      if (hincol[j] > 0) {

	//if (cup[j] - clo[j] <= ztolzb && cup[j] - clo[j] > ZTOLDP) {
	//  cup[j] = clo[j];
	//}
	CoinBigIndex kcs = mcstrt[j];
	CoinBigIndex kce = kcs + hincol[j];
	int iflagu = 0;
	int iflagl = 0;
	double ddjlo = maxmin * dcost[j];
	double ddjhi = ddjlo;

	for (CoinBigIndex k = kcs; k < kce; k++) {
	  int i = hrow[k];
	  double coeff = colels[k];

	  if (coeff > 0.0) {
	    ddjhi -= coeff * rdmin[i];
	    if (rdmin[i] >= ekkinf2) {
	      iflagu = 1;
	      ddjhi = 0.0;
	    }
	    ddjlo -= coeff * rdmax[i];
	    if (rdmax[i] <= -ekkinf2) {
	      iflagl = 1;
	      ddjlo = 0.0;
	    }
	  } else {
	    ddjhi -= coeff * rdmax[i];
	    if (rdmax[i] <= -ekkinf2) {
	      iflagu = 1;
	      ddjhi = 0.0;
	    }
	    ddjlo -= coeff * rdmin[i];
	    if (rdmin[i] >= ekkinf2) {
	      iflagl = 1;
	      ddjlo = 0.0;
	    }
	  }
	}

#if	PRESOLVE_TIGHTEN_DUALS
	djmin[j] = (iflagl ?  PRESOLVE_INF : ddjlo);
	djmax[j] = (iflagu ? -PRESOLVE_INF : ddjhi);
#endif

	if (ddjlo > ZTOLDP && iflagl == 0) {
	  // dj>0 at optimality ==> must be at lower bound
	  if (clo[j] <= -ekkinf) {
	    prob->messageHandler()->message(CLP_PRESOLVE_COLUMNBOUNDB,
					     prob->messages())
					       <<j
					       <<CoinMessageEol;
	    prob->status_ |= 2;
	    break;
	  } else {
	    fixdown_cols[nfixdown_cols++] = j;
	  }
	} else if (ddjhi < -ZTOLDP && iflagu == 0) {
	  // dj<0 at optimality ==> must be at upper bound
	  if (cup[j] >= ekkinf) {
	    prob->messageHandler()->message(CLP_PRESOLVE_COLUMNBOUNDA,
					     prob->messages())
					       <<j
					       <<CoinMessageEol;
	    prob->status_ |= 2;
	    break;
	  } else {
	    fixup_cols[nfixup_cols++] = j;
	  }
	}
      }

    int tightened = 0;

    // I don't know why I stopped doing this.
#if	PRESOLVE_TIGHTEN_DUALS
    // tighten row dual bounds, as described on p. 229
    for (int i = 0; i < nrows; i++) {
      bool no_ub = (rup[i] >= ekkinf);
      bool no_lb = (rlo[i] <= -ekkinf);
      
      if ((no_ub ^ no_lb) == true) {
	CoinBigIndex krs = mrstrt[i];
	CoinBigIndex kre = krs + hinrow[i];
	double rmax  = rdmax[i];
	double rmin  = rdmin[i];

	// it may be that we could just use rmin/rmax, but I'm not sure,
	// since I don't update djmax0/djmin0
	double rmax0 = rmax;
	double rmin0 = rmin;

	// all row columns are non-empty
	for (CoinBigIndex k=krs; k<kre; k++) {
	  double coeff = rowels[k];
	  double djmax0 = djmax[hcol[k]];
	  double djmin0 = djmin[hcol[k]];

	  if (no_ub) {
	    if (coeff > ZTOLDP && djmax0 > -PRESOLVE_INF && rmin0 < PRESOLVE_INF) {
	      double bnd = djmax0 / coeff + rmin0;
	      if (rmax < bnd) {
#if	DEBUG_PRESOLVE
		printf("MAX TIGHT[%d,%d]:  %g --> %g\n", i,hrow[k], rdmax[i], bnd);
#endif
		rdmax[i] = rmax = bnd;
		tightened = 1;
	      }
	    } else if (coeff < -ZTOLDP && djmax0 > -PRESOLVE_INF && rmax0 > -PRESOLVE_INF) {
	      double bnd = djmax0 / coeff + rmax0;
	      if (rmin > bnd) {
#if	DEBUG_PRESOLVE
		printf("MIN TIGHT[%d,%d]:  %g --> %g\n", i, hrow[k], rdmin[i], bnd);
#endif
		rdmin[i] = rmin = bnd;
		tightened = 1;
	      }		
	    }
	  } else {	// no_lb
	    if (coeff > ZTOLDP && djmin0 < PRESOLVE_INF && rmax0 > -PRESOLVE_INF) {
	      double bnd = djmin0 / coeff + rmax0;
	      if (rmin > bnd) {
#if	DEBUG_PRESOLVE
		printf("MIN1 TIGHT[%d,%d]:  %g --> %g\n", i, hrow[k], rdmin[i], bnd);
#endif
		rdmin[i] = rmin = bnd;
		tightened = 1;
	      }
	    } else if (coeff < -ZTOLDP && djmin0 < PRESOLVE_INF && rmin0 < PRESOLVE_INF) {
	      double bnd = djmin0 / coeff + rmin0;
	      if (rmax < bnd) {
#if	DEBUG_PRESOLVE
		printf("MAX TIGHT1[%d,%d]:  %g --> %g\n", i,hrow[k], rdmax[i], bnd);
#endif
		rdmax[i] = rmax = bnd;
		tightened = 1;
	      }
	    }
	  }
	}
      }
    }
#endif

    if (!tightened)
      break;
#if	PRESOLVE_TIGHTEN_DUALS
    else
      printf("DUAL TIGHTENED!  %d\n", nactions);
#endif
  }

  if (nfixup_cols) {
#if	DEBUG_PRESOLVE
    printf("NDUAL:  %d\n", nfixup_cols);
#endif
    next = make_fixed_action::presolve(prob, fixup_cols, nfixup_cols, false, next);
  }

  if (nfixdown_cols) {
#if	DEBUG_PRESOLVE
    printf("NDUAL:  %d\n", nfixdown_cols);
#endif
    next = make_fixed_action::presolve(prob, fixdown_cols, nfixdown_cols, true, next);
  }

  delete[]rdmin;
  delete[]rdmax;

  delete[]fixup_cols;
  delete[]fixdown_cols;

#if	PRESOLVE_TIGHTEN_DUALS
  delete[]djmin;
  delete[]djmax;
#endif

  return (next);
}
