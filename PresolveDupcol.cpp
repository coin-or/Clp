#include <stdio.h>
#include <math.h>

#include "PresolveMatrix.hpp"
#include "PresolveFixed.hpp"
#include "PresolveDupcol.hpp"
#include "CoinSort.hpp"
#include "PresolveUseless.hpp"
#include "ClpMessage.hpp"

#define DSEED2 2147483647.0

const char *dupcol_action::name() const
{
  return ("dupcol_action");
}

const char *duprow_action::name() const
{
  return ("duprow_action");
}

void init_random_vec(double *work, int n)
{
  double deseed = 12345678.0;

  for (int i = 0; i < n; ++i) {
    deseed *= 16807.;
    int jseed = (int) (deseed /    DSEED2);
    deseed -= (double) jseed * DSEED2;
    double random = deseed /  DSEED2;

    work[i]=random;
  }
}

int compute_sums(const int *len, const CoinBigIndex *starts,
		 /*const*/ int *index, /*const*/ double *elems,	// index/elems are sorted
		 const double *work,
		 double *sums, int *sorted, int n)
{
  int nlook=0;
  for (int i = 0; i < n; ++i)
    if (len[i] > 0 /*1?*/) {
      CoinBigIndex kcs = starts[i];
      CoinBigIndex kce = kcs + len[i];

      double value=0.0;

      // sort all columns for easy comparison in next loop
      CoinSort_2(index+kcs,index+kcs+len[i],elems+kcs);
      //ekk_sort2(index+kcs, elems+kcs, len[i]);

      for (CoinBigIndex k=kcs;k<kce;k++) {
	int irow=index[k];
	value += work[irow]*elems[k];
      }
      sums[nlook]=value;
      sorted[nlook]=i;
      ++nlook;
    }
  return (nlook);
}

dupcol_action::dupcol_action(int nactions,
		const action *actions,
		const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions)
{
}


// This is just ekkredc5, adapted into the new framework.
// the datasets scorpion.mps and allgrade.mps have duplicate columns.
const PresolveAction *dupcol_action::presolve(PresolveMatrix *prob,
					       const PresolveAction *next)
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int ncols		= prob->ncols_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;
  double * sol = prob->sol_;

  double *rowels	= prob->rowels_;
  int *hcol	= prob->hcol_;
  const CoinBigIndex *mrstrt	= prob->mrstrt_;
  int *hinrow		= prob->hinrow_;
  int nrows		= prob->nrows_;

  //  double *rlo	= prob->rlo_;
  //  double *rup	= prob->rup_;

  double *dcost	= prob->cost_;

  const char *integerType = prob->integerType_;

  double maxmin	= prob->maxmin_;

  action *actions	= new action [ncols];
  int nactions = 0;

  int *fixed_down	= new int[ncols];
  int nfixed_down	= 0;
  int *fixed_up		= new int[ncols];
  int nfixed_up		= 0;

  double * workcol = new double[ncols];
  double * workrow = new double[nrows];
  int * sort = new int[ncols];

  // initialize workrow to have "random" values
  init_random_vec(workrow, nrows);

  // If we sum all the coefficients of a column,
  // then duplicate columns must have the same sum.
  // However, other columns may also.
  // To weed most of them out, we multiply the coefficients
  // by the random value associated with the row.
  int nlook = compute_sums(hincol, mcstrt, hrow, colels, workrow, workcol, sort, ncols);

#if 1
  // not tested yet
  if (maxmin < 0.0)
    return (next);
#endif

  CoinSort_2(workcol,workcol+nlook,sort);
  //ekk_sortonDouble(workcol,sort,nlook);

#if 0
  // It may be the case that several columns are duplicate.
  // If not all have the same cost, then we have to make sure
  // that we set the most expensive one to its minimum
  // now sort in each class by cost
  {
    double dval = workcol[0];
    int first = 0;
    for (int jj = 1; jj < nlook; jj++) {
      while (workcol[jj]==dval) 
	jj++;

      if (first + 1 < jj) {
	double buf[jj - first];
	for (int i=first; i<jj; ++i)
	  buf[i-first] = dcost[sort[i]]*maxmin;

	CoinSort_2(buf,buf+jj-first,sort+first);
	//ekk_sortonDouble(buf,&sort[first],jj-first);
      }
    }
  }
#endif

  // it appears to be the case that this loop is finished,
  // there may still be duplicate cols left.
  // I haven't done anything about that yet.
  for (int jj = 1; jj < nlook; jj++) {
    if (workcol[jj]==workcol[jj-1]) {
      int ithis=sort[jj];
      int ilast=sort[jj-1];
      CoinBigIndex kcs = mcstrt[ithis];
      CoinBigIndex kce = kcs + hincol[ithis];

      if (hincol[ithis] == hincol[ilast]&&!prob->colProhibited2(ithis)) {
	int ishift = mcstrt[ilast] - kcs;
	CoinBigIndex k;
	for (k=kcs;k<kce;k++) {
	  if (hrow[k] != hrow[k+ishift] ||
	      colels[k] != colels[k+ishift]) {
	    break;
	  }
	}
	if (k == kce) {
	  // these really are duplicate columns

	  /* now check bounds and costs to see what is what */
	  double clo1=clo[ilast];
	  double cup1=cup[ilast];
	  double clo2=clo[ithis];
	  double cup2=cup[ithis];
	  double dcost1=dcost[ilast]*maxmin;
	  double dcost2=dcost[ithis]*maxmin;
          double newSolution = sol[ilast]+sol[ithis];
	  //int itype=27;
	  if (clo1==cup1||clo2==cup2)
	    abort();

	  if (/*ddjs[ilast]||ddjs[ithis]*/
	      integerType[ilast] || integerType[ithis]) {
#ifdef PRINT_DEBUG
	    printf("at least one of %d or %d integer\n",ilast-nrowsmx,
		   ithis-nrowsmx);
#endif
	    continue;
	  } else if (dcost1==dcost2) {
	    //
	    // SAME COSTS
	    //
	    double l_j = clo[ithis];
	    double u_j = cup[ithis];
	    double l_k = clo[ilast];
	    double u_k = cup[ilast];

	    if (! (l_j + u_k <= l_k + u_j)) {
	      swap(ilast, ithis);
	      swap(clo1, clo2);
	      swap(cup1, cup2);
	      //swap(dcost1, dcost2);
	    }

	    //PRESOLVE_STMT(printf("DUPCOL (%d,%d)\n", ithis, ilast));

	    /* identical cost so can add ranges */
	    clo1 += clo2;
	    cup1 += cup2;
	    if (clo1<-1.0e20) {
	      clo1=-1.0e31;
	    }
	    if (cup1>1.0e20) {
	      cup1=1.0e31;
	    }
	    /*dfix=0.0;*/
	    //itype=26;

	    {
	      action *s = &actions[nactions];	  
	      nactions++;

	      s->thislo = clo[ithis];
	      s->thisup = cup[ithis];
	      s->lastlo = clo[ilast];
	      s->lastup = cup[ilast];
	      s->ithis  = ithis;
	      s->ilast  = ilast;

	      s->nincol	= hincol[ithis];
	      s->colels	= copyOfArray(&colels[mcstrt[ithis]], hincol[ithis]);
	      s->colrows= copyOfArray(&hrow[mcstrt[ithis]], hincol[ithis]);
	    }

	    // adjust ilast
	    clo[ilast]=clo1;
	    cup[ilast]=cup1;
	    sol[ilast] = newSolution;
	    if (prob->getColumnStatus(ilast)==PrePostsolveMatrix::basic||
		prob->getColumnStatus(ithis)==PrePostsolveMatrix::basic)
	      prob->setColumnStatus(ilast,PrePostsolveMatrix::basic);

	    // delete ithis
	    dcost[ithis] = 0.0;
	    sol[ithis]=clo2;
	    {
	      CoinBigIndex kcs = mcstrt[ithis];
	      CoinBigIndex kce = kcs + hincol[ithis];
	      for (CoinBigIndex k=kcs; k<kce; ++k)
		// delete ithis from its rows
		presolve_delete_from_row(hrow[k], ithis, mrstrt, hinrow, hcol, rowels);
	    }
	    hincol[ithis] = 0;

	    /* make sure fixed one out of way */
	    sort[jj] = ilast;
	    //if (ithis==sort[jj]) {
	    //sort[jj]=sort[jj-1];
	    //}
	  } else {
	    //
	    // (dcost1!=dcost2) - DIFFERENT COSTS
	    //
#define	PRINT_DEBUG
	    /* look to see whether we can drop one ? */
	    if (clo1<-1.0e20 && clo2<-1.0e20) {
	      // both unbounded below 
	      /* can't cope for now */
#ifdef PRINT_DEBUG
	      printf("** Odd columns %d and %d\n", ilast,ithis);
#endif
	      continue;
	    } else if (clo1<-1.0e20) {
	      printf("ILAST:  %d %d\n", ilast, ithis);
	      swap(ilast, ithis);
	      printf("ILAST1:  %d %d\n", ilast, ithis);
	      swap(clo1, clo2);
	      swap(cup1, cup2);
	      swap(dcost1, dcost2);
	    }
	    // first one (ilast) bounded below - PRESOLVEASSERT(! (clo1<-1.0e20) );

	    if (cup1>1.0e20) {
	      /* ilast can go to plus infinity */
	      if (clo2<-1.0e20) {
		if (dcost1<dcost2) {
#ifdef PRINT_DEBUG
		  printf("** Problem seems unbounded %d and %d\n", ilast,ithis);
#endif
		  break;
		} else {
		  continue;
		}
	      } else if (cup2>1.0e20) {
		//PRESOLVE_STMT(printf("DUPCOL1 (%d,%d) %g\n", ithis, ilast, workcol[jj]));
		// both have an lb, both have no ub
		// the more expensive one would end up at its lb
		if (dcost1>dcost2) {
		  swap(ilast, ithis);
		  swap(clo1, clo2);
		  swap(cup1, cup2);
		  swap(dcost1, dcost2);
		}
		sol[ilast] = newSolution;
		if (prob->getColumnStatus(ilast)==PrePostsolveMatrix::basic||
		    prob->getColumnStatus(ithis)==PrePostsolveMatrix::basic)
		  prob->setColumnStatus(ilast,PrePostsolveMatrix::basic);
		
		sol[ithis]=clo2;
		fixed_down[nfixed_down++] = ithis ;
		sort[jj]                  = ilast;
	      } else {
		/* ithis ranged - last not */
		if (dcost2>dcost1) {
		  //PRESOLVE_STMT(printf("DUPCOL2 (%d,%d)\n", ithis, ilast));
		  /* set ithis to lower bound */
		  fixed_down[nfixed_down++] = ithis;
		  sort[jj] = ilast;
		  sol[ilast] = newSolution;
		  if (prob->getColumnStatus(ilast)==PrePostsolveMatrix::basic||
		      prob->getColumnStatus(ithis)==PrePostsolveMatrix::basic)
		    prob->setColumnStatus(ilast,PrePostsolveMatrix::basic);
		  
		  sol[ithis]=clo2;
		} else {
		  /* no good */
		  continue;
		}
	      }
	    } else {
	      /* ilast ranged */
	      if (clo2<-1.0e20) {
		// ithis has no lb
		if (dcost2>dcost1) {
		  //PRESOLVE_STMT(printf("DUPCOL3 (%d,%d)\n", ithis, ilast));
		  /* set ilast to upper bound */
		  fixed_up[nfixed_up++] = ilast;
		  sort[jj] = ithis;
		  sol[ithis] = newSolution;
		  if (prob->getColumnStatus(ithis)==PrePostsolveMatrix::basic||
		      prob->getColumnStatus(ilast)==PrePostsolveMatrix::basic)
		    prob->setColumnStatus(ithis,PrePostsolveMatrix::basic);
		  
		  sol[ilast]=clo1;
		} else {
		  /* need both */
		  continue;
		}
	      } else if (cup2>1.0e20) {
		if (dcost2<dcost1) {
		  //PRESOLVE_STMT(printf("DUPCOL4 (%d,%d)\n", ithis, ilast));
		  /* set ilast to lower bound */
		  //SWAP(int, ilast, ithis);
		  fixed_down[nfixed_down++] = ilast;
		  sort[jj] = ithis;
		  sol[ithis] = newSolution;
		  if (prob->getColumnStatus(ithis)==PrePostsolveMatrix::basic||
		      prob->getColumnStatus(ilast)==PrePostsolveMatrix::basic)
		    prob->setColumnStatus(ithis,PrePostsolveMatrix::basic);
		  
		  sol[ilast]=clo1;
		} else {
		  /* need both */
		  continue;
		}
	      } else {
		/* both ranged - no hope */
		continue;
	      }
	    }
	  }
	}
      }
    }
  }

  delete[]workrow;
  delete[]workcol;
  delete[]sort;


  if (nactions) {
#if	PRESOLVE_SUMMARY
    printf("DUPLICATE COLS:  %d\n", nactions);
#endif
    next = new dupcol_action(nactions, copyOfArray(actions,nactions), next);
  }
  deleteAction(actions,action*);
  if (nfixed_down)
    next = make_fixed_action::presolve(prob, fixed_down, nfixed_down,
				       true,
				       next);

  if (nfixed_up)
    next = make_fixed_action::presolve(prob, fixed_up, nfixed_up,
				       false,
				       next);

  delete[]fixed_down;
  delete[]fixed_up;

  return (next);
}

void create_col(int col, int n, int *rows, double *els,
		CoinBigIndex *mcstrt, double *colels, int *hrow, int *link,
		CoinBigIndex *free_listp)
{
  CoinBigIndex free_list = *free_listp;
  int xstart = NO_LINK;
  for (int i=0; i<n; ++i) {
    CoinBigIndex k = free_list;
    free_list = link[free_list];

    check_free_list(free_list);

    hrow[k]   = rows[i];
    colels[k] = els[i];
    link[k] = xstart;
    xstart = k;
  }
  mcstrt[col] = xstart;
  *free_listp = free_list;
}


void dupcol_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions = nactions_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  double *sol	= prob->sol_;
  double *dcost	= prob->cost_;
  
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;

  double *rcosts	= prob->rcosts_;

  CoinBigIndex free_list		= prob->free_list_;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {
    int icol  = f->ithis;	// was fixed
    int icol2 = f->ilast;	// was kept

    dcost[icol] = dcost[icol2];
    clo[icol] = f->thislo;
    cup[icol] = f->thisup;
    clo[icol2] = f->lastlo;
    cup[icol2] = f->lastup;

    create_col(icol, f->nincol, f->colrows, f->colels,
	       mcstrt, colels, hrow, link, &free_list);
    hincol[icol] = hincol[icol2]; // right?

    int nfinite = ((fabs(f->thislo) < PRESOLVE_INF) +
		   (fabs(f->thisup) < PRESOLVE_INF) +
		   (fabs(f->lastlo) < PRESOLVE_INF) +
		   (fabs(f->lastup) < PRESOLVE_INF));

    if (nfinite > 1) {
       //      double l_j = f->thislo;
      double u_j = f->thisup;
      double l_k = f->lastlo;
      double u_k = f->lastup;
      double x_k_sol = sol[icol2];

      prob->setColumnStatus(icol,prob->getColumnStatus(icol2));
      if (x_k_sol <= l_k + u_j) {
	sol[icol2] = l_k;
	sol[icol] = x_k_sol - l_k;
	prob->setColumnStatus(icol2,PrePostsolveMatrix::atLowerBound);
      } else {
	sol[icol2] = u_k;
	sol[icol] = x_k_sol - u_k;
	prob->setColumnStatus(icol2,PrePostsolveMatrix::atLowerBound);
      }
    } else if (nfinite == 1) {
      double x_k_sol = sol[icol2];

      if (fabs(f->thislo) < PRESOLVE_INF) {
	prob->setColumnStatus(icol,PrePostsolveMatrix::atLowerBound);
	sol[icol] = f->thislo;
	sol[icol2] = x_k_sol - sol[icol];
      } else if (fabs(f->thisup) < PRESOLVE_INF) {
	prob->setColumnStatus(icol,PrePostsolveMatrix::atUpperBound);
	sol[icol] = f->thisup;
	sol[icol2] = x_k_sol - sol[icol];
      } else if (fabs(f->lastlo) < PRESOLVE_INF) {
	prob->setColumnStatus(icol,prob->getColumnStatus(icol2));
	prob->setColumnStatus(icol2,PrePostsolveMatrix::atLowerBound);
	sol[icol2] = f->lastlo;
	sol[icol] = x_k_sol - sol[icol2];
      } else {
	// (fabs(f->lastup) < PRESOLVE_INF)
	prob->setColumnStatus(icol,prob->getColumnStatus(icol2));
	prob->setColumnStatus(icol2,PrePostsolveMatrix::atUpperBound);
	sol[icol2] = f->lastup;
	sol[icol] = x_k_sol - sol[icol2];
      }
    } else {
      // both free!  superbasic time
      sol[icol] = 0.0;	// doesn't matter
      prob->setColumnStatus(icol2,PrePostsolveMatrix::isFree);
    }

    // row activity doesn't change
    // dj of both variables is the same
    rcosts[icol] = rcosts[icol2];

#ifdef DEBUG_PRESOLVE
    const double ztolzb = prob->ztolzb_;
    if (! (clo[icol] - ztolzb <= sol[icol] && sol[icol] <= cup[icol] + ztolzb))
	     printf("BAD DUPCOL BOUNDS:  %g %g %g\n", clo[icol], sol[icol], cup[icol]);
    if (! (clo[icol2] - ztolzb <= sol[icol2] && sol[icol2] <= cup[icol2] + ztolzb))
	     printf("BAD DUPCOL BOUNDS:  %g %g %g\n", clo[icol2], sol[icol2], cup[icol2]);
#endif
  }
  prob->free_list_ = free_list;
}





// This is just ekkredc4, adapted into the new framework.
const PresolveAction *duprow_action::presolve(PresolveMatrix *prob,
					       const PresolveAction *next)
{
  //  CoinBigIndex *mcstrt		= prob->mcstrt_;
  //  int *hincol		= prob->hincol_;
  int ncols		= prob->ncols_;

  double *rowels	= prob->rowels_;
  /*const*/ int *hcol	= prob->hcol_;
  const CoinBigIndex *mrstrt	= prob->mrstrt_;
  int *hinrow		= prob->hinrow_;
  int nrows		= prob->nrows_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  int nuseless_rows = 0;

  double * workcol = new double[ncols];
  double * workrow = new double[nrows];
  int * sort = new int[nrows];

  init_random_vec(workcol, ncols);

  int nlook = compute_sums(hinrow, mrstrt, hcol, rowels, workcol, workrow, sort, nrows);
  CoinSort_2(workrow,workrow+nlook,sort);

  double dval = workrow[0];
  for (int jj = 1; jj < nlook; jj++) {
    if (workrow[jj]==dval) {
      int ithis=sort[jj];
      int ilast=sort[jj-1];
      CoinBigIndex krs = mrstrt[ithis];
      CoinBigIndex kre = krs + hinrow[ithis];
      if (hinrow[ithis] == hinrow[ilast]) {
	int ishift = mrstrt[ilast] - krs;
	CoinBigIndex k;
	for (k=krs;k<kre;k++) {
	  if (hcol[k] != hcol[k+ishift] ||
	      rowels[k] != rowels[k+ishift]) {
	    break;
	  }
	}
	if (k == kre) {
	  /* now check rhs to see what is what */
	  double rlo1=rlo[ilast];
	  double rup1=rup[ilast];
	  double rlo2=rlo[ithis];
	  double rup2=rup[ithis];

	  int idelete=-1;
	  if (rlo1<=rlo2) {
	    if (rup2<=rup1) {
	      /* this is strictly tighter than last */
	      idelete=ilast;
	    } else if (fabs(rlo1-rlo2)<1.0e-12) {
	      /* last is strictly tighter than this */
	      idelete=ithis;
	      // swap so can carry on deleting
	      sort[jj-1]=ithis;
	      sort[jj]=ilast;
	    } else {
	      /* overlapping - could merge */
#ifdef PRINT_DEBUG
	      printf("overlapping duplicate row %g %g, %g %g\n",
		     rlo1,rup1,rlo2,rup2);
#endif
	      if (rup1<rlo2) {
		// infeasible
		prob->status_|= 1;
		// wrong message - correct if works
		prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
						prob->messages())
						  <<ithis
						  <<rlo[ithis]
						  <<rup[ithis]
						  <<CoinMessageEol;
		break;
	      }
	    }
	  } else {
	    // rlo2<=rlo1
	    if (rup1<=rup2) {
	      /* last is strictly tighter than this */
	      idelete=ithis;
	      // swap so can carry on deleting
	      sort[jj-1]=ithis;
	      sort[jj]=ilast;
	    } else {
	      /* overlapping - could merge */
#ifdef PRINT_DEBUG
	      printf("overlapping duplicate row %g %g, %g %g\n",
		     rlo1,rup1,rlo2,rup2);
#endif
	      if (rup2<rlo1) {
		// infeasible
		prob->status_|= 1;
		// wrong message - correct if works
		prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
						prob->messages())
						  <<ithis
						  <<rlo[ithis]
						  <<rup[ithis]
						  <<CoinMessageEol;
		break;
	      }
	    }
	  }
	  if (idelete>=0) 
	    sort[nuseless_rows++]=idelete;
	}
      }
    }
    dval=workrow[jj];
  }

  delete[]workrow;
  delete[]workcol;


  if (nuseless_rows) {
#if	PRESOLVE_SUMMARY
    printf("DUPLICATE ROWS:  %d\n", nuseless_rows);
#endif
    next = useless_constraint_action::presolve(prob,
					       sort, nuseless_rows,
					       next);
  }
  delete[]sort;

  return (next);
}

void duprow_action::postsolve(PostsolveMatrix *prob) const
{
  printf("STILL NO POSTSOLVE FOR DUPROW!\n");
  abort();
}
