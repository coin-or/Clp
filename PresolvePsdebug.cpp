// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <new>
#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "PresolveMatrix.hpp"


static inline const double *
lookup_by_col(const int *mcstrt, const int *hrow,
	      const int *hincol, const int *hinrow, const double *colels,
	      int row, int col)
{
  int kcs = mcstrt[col];
  int kce = kcs + hincol[col];
  int k;

  if (hinrow[row] <= 0)
    return (0);

  for (k=kcs; k<kce; k++) {
    int r = hrow[k];
    if (r == row)
      return (&colels[k]);
  }
  return (0);
}

static inline void
no_dups(const char *done,
	const int *mcstrt, const int *hrow, const int *hincol, int ncols)
{
#if	DEBUG_PRESOLVE
  for (int jcol=0; jcol<ncols; jcol++) 
    if ((!done || done[jcol]) && hincol[jcol] > 0) {
      int kcs = mcstrt[jcol];
      int kce = kcs + hincol[jcol];

      for (int k=kcs; k<kce; k++) {
	int row = hrow[k];

	PRESOLVEASSERT(presolve_find_row1(row, k+1, kce, hrow) == kce);
      }
    }
#endif
}

void presolve_no_zeros(const int *mcstrt, const double *colels, const int *hincol, int ncols)
{
#if	DEBUG_PRESOLVE
  for (int jcol=0; jcol<ncols; jcol++) 
    if (hincol[jcol] > 0) {
      int kcs = mcstrt[jcol];
      int kce = kcs + hincol[jcol];

      for (int k=kcs; k<kce; k++) {
	PRESOLVEASSERT(fabs(colels[k]) > ZTOLDP);
      }
    }
#endif
}


void presolve_hincol_ok(const int *mcstrt, const int *hincol,
	       const int *hinrow,
	       const int *hrow, int ncols)
{
#if	CHECK_CONSISTENCY
  int jcol;

  for (jcol=0; jcol<ncols; jcol++) 
    if (hincol[jcol] > 0) {
      int kcs = mcstrt[jcol];
      int kce = kcs + hincol[jcol];
      int n=0;
      
      int k;
      for (k=kcs; k<kce; k++) {
	int row = hrow[k];
	if (hinrow[row] > 0)
	  n++;
      }
      if (n != hincol[jcol])
	abort();
    }
#endif
}

/*
 * The linked list 
 */
void presolve_links_ok(presolvehlink *link, int *starts, int *lengths, int n)
{
#if	CHECK_CONSISTENCY
  int i;

  for (i=0; i<n; i++) {
    int pre = link[i].pre;
    int suc = link[i].suc;

    if (pre != NO_LINK) {
      PRESOLVEASSERT(0 <= pre && pre <= n);
      PRESOLVEASSERT(link[pre].suc == i);
    }
    if (suc != NO_LINK) {
      PRESOLVEASSERT(0 <= suc && suc <= n);
      PRESOLVEASSERT(link[suc].pre == i);
    }
  }

  for (i=0; i<n; i++) 
    if (link[i].pre == NO_LINK)
      break;
  PRESOLVEASSERT(i<n);

  while (i != NO_LINK) {
    if (link[i].suc != NO_LINK) 
      PRESOLVEASSERT(starts[i] + lengths[i] <= starts[link[i].suc]);
    i = link[i].suc;
  }
#endif
}

// I've forgotton what this is all about
void check_pivots(const int *mrstrt, const int *hinrow, const int *hcol, int nrows,
		  const unsigned char *colstat, const unsigned char *rowstat,
		  int ncols)
{
#if 0
  int i;
  int nbasic = 0;
  int gotone = 1;
  int stillmore;

  return;

  int *bcol = new int[nrows];
  memset(bcol, -1, nrows*sizeof(int));

  char *coldone = new char[ncols];
  memset(coldone, 0, ncols);

  while (gotone) {
    gotone = 0;
    stillmore = 0;
    for (i=0; i<nrows; i++)
      if (!prob->rowIsBasic(i)) {
	int krs = mrstrt[i];
	int kre = mrstrt[i] + hinrow[i];
	int nb = 0;
	int kk;
	for (int k=krs; k<kre; k++)
	  if (prob->columnIsBasic(hcol[k]) && !coldone[hcol[k]]) {
	    nb++;
	    kk = k;
	    if (nb > 1)
	      break;
	  }
	if (nb == 1) {
	  PRESOLVEASSERT(bcol[i] == -1);
	  bcol[i] = hcol[kk];
	  coldone[hcol[kk]] = 1;
	  nbasic++;
	  gotone = 1;
	}
	else
	  stillmore = 1;
      }
  }
  PRESOLVEASSERT(!stillmore);

  for (i=0; i<nrows; i++)
    if (prob->rowIsBasic(i)) {
      int krs = mrstrt[i];
      int kre = mrstrt[i] + hinrow[i];
      for (int k=krs; k<kre; k++)
	PRESOLVEASSERT(!prob->columnIsBasic(hcol[k]) || coldone[hcol[k]]);
      nbasic++;
    }
  PRESOLVEASSERT(nbasic == nrows);
#endif
}



static inline void a_ok(PostsolveMatrix *prob)
{
#if 0
  static int warned = 0;
  if (!warned) {
    warned = 1;
    printf("********** NO CHECKING!!!!!!! \n");
  }

  double *colels	= prob->colels;
  int *hrow		= prob->hrow;
  int *mcstrt		= prob->mcstrt;
  int *hincol		= prob->hincol;

  int *colstat	= prob->colstat;
  int *rowstat	= prob->rowstat;

  double *clo	= prob->clo;
  double *cup	= prob->cup;

  double *dcost	= prob->cost;

  double *sol	= prob->sol;
  double *rcosts	= prob->rcosts;

  char *cdone	= prob->cdone;
  char *rdone	= prob->rdone;

  const double ztoldj	= prob->ztoldj;
  const double ztolzb	= prob->ztolzb;

  double *rowduals = prob->rowduals;

  int ncols0		= prob->ncols0;

#if 0
  {
    int ncols		= prob->ncols;
    int nrows		= prob->nrows;
    int *mrstrt		= prob->mrstrt;
    int *hinrow		= prob->hinrow;
    int *hcol		= prob->hcol;

    no_dups(cdone, mcstrt, hrow, hincol, ncols);
    no_dups(rdone, mrstrt, hcol, hinrow, nrows);
  }
#endif

  /* DEBUG HACK */
  static double *warned_rcosts;
  if (!warned_rcosts) {
    warned_rcosts = new double[ncols0];

    for (int i=0; i<ncols0; i++) 
      warned_rcosts[i] = rcosts[i];
  }

  for (int j=0; j<ncols0; j++)
    if (cdone[j]) {
      //int i = &paction[npaction-1]-pa;
      int i = -666;

      if (prob->columnIsBasic(j)) {
	if (! (fabs(rcosts[j]) < ztoldj) &&
	    warned_rcosts[j] != rcosts[j])
	  printf("ITER[%d]: basic rcosts[%d]:  %g\n", i, j, rcosts[j]);
      } else if (fabs(sol[j] - cup[j]) < ztolzb &&
		 ! (fabs(sol[j] - clo[j]) < ztolzb)) {
	if (! (rcosts[j] <= ztoldj) &&
	    warned_rcosts[j] != rcosts[j])
	  printf("ITER[%d]: ub rcosts[%d]:  %g\n", i, j, rcosts[j]);
      } else if (fabs(sol[j] - clo[j]) < ztolzb &&
		 ! (fabs(sol[j] - cup[j]) < ztolzb)){
	if (! (rcosts[j] >= -ztoldj) &&
	    warned_rcosts[j] != rcosts[j])
	  printf("ITER[%d]: lb rcosts[%d]:  %g\n", i, j, rcosts[j]);
      } else if (! (fabs(sol[j] - clo[j]) < ztolzb) &&
		 ! (fabs(sol[j] - cup[j]) < ztolzb)) {
	printf("SUPERBASIC (cost=%g):  %d %g < %g < %g\n",
	       dcost[j], j, clo[j], sol[j], cup[j]);
      }

      {
	int kcs = mcstrt[j];
	int kce = mcstrt[j] + hincol[j];
	int k;
	double dj = dcost[j];
	int row0 = (kcs<kce ? hrow[kcs] : 0);

	for (k=kcs; k<kce; k++) {
	  int row = hrow[k];
	  PRESOLVEASSERT(rdone[row]);
	  PRESOLVEASSERT(row != row0 || k==kcs);
	  if (rdone[row])
	    dj -= rowduals[row] * colels[k];
	}
	if (! (fabs(rcosts[j] - dj) < ztoldj) &&
	    warned_rcosts[j] != rcosts[j])
	  printf("ITER[%d]: rcosts[%d]:  %g <-> %g\n", i, j, rcosts[j], dj);
      }
    }

  {
    int i;
    for (i=0; i<ncols0; i++) 
      warned_rcosts[i] = rcosts[i];
  }
#endif
}


