#include <stdio.h>
#include <math.h>

#include "PresolveMatrix.hpp"
#include "PresolveEmpty.hpp"	// for DROP_COL/DROP_ROW
#include "PresolveFixed.hpp"
#include "PresolveZeros.hpp"
#include "PresolveSubst.hpp"
#include "ClpMessage.hpp"
#include "CoinSort.hpp"

inline void prepend_elem(int jcol, double coeff, int irow,
		    CoinBigIndex *mcstrt,
		    double *colels,
		    int *hrow,
		    int *link, CoinBigIndex *free_listp)
{
  CoinBigIndex kk = *free_listp;
  *free_listp = link[*free_listp];
  check_free_list(*free_listp);

  link[kk] = mcstrt[jcol];
  mcstrt[jcol] = kk;
  colels[kk] = coeff;
  hrow[kk] = irow;
}

const char *subst_constraint_action::name() const
{
  return ("subst_constraint_action");
}

void compact_rep(double *elems, int *indices, CoinBigIndex *starts, const int *lengths, int n,
		 const presolvehlink *link);

inline int min(int x, int y)
{
  return (x < y) ? x : y;
}


// copy of expand_col; have to rename params
static void expand_row(CoinBigIndex *mcstrt, 
		    double *colels,
		       int *hrow, // int *hcol,
		       //int *hinrow,
		       int *hincol,
		    presolvehlink *clink, int ncols,

		       int icolx
		       ///, int icoly
		       )
{
  /////CoinBigIndex kcs = mcstrt[icoly];
  /////CoinBigIndex kce = kcs + hincol[icoly];
  CoinBigIndex kcsx = mcstrt[icolx];
  CoinBigIndex kcex = kcsx + hincol[icolx];

  const int maxk = mcstrt[ncols];	// (22)

  // update col rep - need to expand the column, though.
  int nextcol = clink[icolx].suc;

  // (22)
  if (kcex + 1 < mcstrt[nextcol] || nextcol == ncols) {
    if (! (kcex + 1 < mcstrt[nextcol])) {
      // nextcol==ncols and no space - must compact
      compact_rep(colels, hrow, mcstrt, hincol, ncols, clink);

      // update vars
      kcsx = mcstrt[icolx];
      kcex = kcsx + hincol[icolx];

      if (! (kcex + 1 < mcstrt[nextcol])) {
	abort();
      }
    }
  } else {
    // this is not the last col 
    // fetch last non-empty col (presolve_make_memlists-1)
    int lastcol = clink[ncols].pre;
    // (clink[icolx].suc != ncols) ==> (icolx != lastcol)

    // put it directly after the last column 
    int newkcsx = mcstrt[lastcol] + hincol[lastcol];

    // well, pad it a bit
    newkcsx += min(hincol[icolx], 5); // slack

    //printf("EXPAND_ROW:  %d %d %d\n", newkcsx, maxk, icolx);

    if (newkcsx + hincol[icolx] + 1 >= maxk) {
      compact_rep(colels, hrow, mcstrt, hincol, ncols, clink);

      // update vars
      kcsx = mcstrt[icolx];
      kcex = kcsx + hincol[icolx];

      newkcsx = mcstrt[lastcol] + hincol[lastcol];

      if (newkcsx + hincol[icolx] + 1 >= maxk) {
	abort();
      }
      // have to adjust various induction variables
      ////kcoly = mcstrt[icoly] + (kcoly - kcs);
      /////kcs = mcstrt[icoly];			// do this for ease of debugging
      /////kce = mcstrt[icoly] + hincol[icoly];
	    
      /////kcolx = mcstrt[icolx] + (kcolx - kcs);	// don't really need to do this
      kcsx = mcstrt[icolx];
      kcex = mcstrt[icolx] + hincol[icolx];
    }

    // move the column - 1:  copy the entries
    memcpy((void*)&hrow[newkcsx], (void*)&hrow[kcsx], hincol[icolx] * sizeof(int));
    memcpy((void*)&colels[newkcsx], (void*)&colels[kcsx], hincol[icolx] * sizeof(double));

    // move the column - 2:  update the memory-order linked list
    PRESOLVE_REMOVE_LINK(clink, icolx);
    PRESOLVE_INSERT_LINK(clink, icolx, lastcol);

    // move the column - 3:  update loop variables to maintain invariant
    mcstrt[icolx] = newkcsx;
    kcsx = newkcsx;
    kcex = newkcsx + hincol[icolx];

#if 0
    hincol[icolx]++;
    kcex = newkcsx + hincol[icolx];

    // move the column - 4:  add the new entry
    hrow[kcex-1] = row;
    colels[kcex-1] = colels[kcoly] * coeff_factor;
#endif
  }
}

// add coeff_factor * rowy to rowx
void add_row(CoinBigIndex *mrstrt, 
	     double *rlo, double * acts, double *rup,
	     double *rowels,
	     int *hcol,
	     int *hinrow,
	     presolvehlink *rlink, int nrows,
	     double coeff_factor,
	     int irowx, int irowy,
	     int *x_to_y)
{
  CoinBigIndex krs = mrstrt[irowy];
  CoinBigIndex kre = krs + hinrow[irowy];
  CoinBigIndex krsx = mrstrt[irowx];
  CoinBigIndex krex = krsx + hinrow[irowx];
  //  const int maxk = mrstrt[nrows];	// (22)

  // if irowx is very long, the searching gets very slow,
  // so we always sort.
  // whatever sorts rows should handle almost-sorted data efficiently
  // (quicksort may not)
  CoinSort_2(hcol+krsx,hcol+krsx+hinrow[irowx],rowels+krsx);
  CoinSort_2(hcol+krs,hcol+krs+hinrow[irowy],rowels+krs);
  //ekk_sort2(hcol+krsx, rowels+krsx, hinrow[irowx]);
  //ekk_sort2(hcol+krs,  rowels+krs,  hinrow[irowy]);

  //printf("%s x=%d y=%d cf=%g nx=%d ny=%d\n",
  // "ADD_ROW:",
  //  irowx, irowy, coeff_factor, hinrow[irowx], hinrow[irowy]);

#if	DEBUG_PRESOLVE
  printf("%s x=%d y=%d cf=%g nx=%d ycols=(",
	 "ADD_ROW:",
	  irowx, irowy, coeff_factor, hinrow[irowx]);
#endif

  // adjust row bounds of rowx;
  // analogous to adjusting bounds info of colx in doubleton,
  // or perhaps adjustment to rlo/rup in elim_doubleton
  //
  // I believe that since we choose a column that is implied free,
  // no other column bounds need to be updated.
  // This is what would happen in doubleton if y's bounds were implied free;
  // in that case,
  // lo1 would never improve clo[icolx] and
  // up1 would never improve cup[icolx].
  {
    double rhsy = rlo[irowy];

    // (1)
    if (-PRESOLVE_INF < rlo[irowx]) {
#if	DEBUG_PRESOLVE
      if (rhsy * coeff_factor)
	printf("ELIM_ROW RLO:  %g -> %g\n",
	       rlo[irowx],
	       rlo[irowx] + rhsy * coeff_factor);
#endif
      rlo[irowx] += rhsy * coeff_factor;
    }
    // (2)
    if (rup[irowx] < PRESOLVE_INF) {
#if	DEBUG_PRESOLVE
      if (rhsy * coeff_factor)
	printf("ELIM_ROW RUP:  %g -> %g\n",
	       rup[irowx],
	       rup[irowx] + rhsy * coeff_factor);
#endif
      rup[irowx] += rhsy * coeff_factor;
    }
    acts[irowx] += rhsy * coeff_factor;
  }

  CoinBigIndex kcolx = krsx;
  CoinBigIndex krex0 = krex;
  int x_to_y_i = 0;

  for (CoinBigIndex krowy=krs; krowy<kre; krowy++) {
    int jcol = hcol[krowy];

    // even though these values are updated, they remain consistent
    PRESOLVEASSERT(krex == krsx + hinrow[irowx]);

    // see if row appears in colx
    // do NOT look beyond the original elements of rowx
    //CoinBigIndex kcolx = presolve_find_row1(jcol, krsx, krex, hcol);
    while (kcolx < krex0 && hcol[kcolx] < jcol)
      kcolx++;

#if	DEBUG_PRESOLVE
    printf("%d%s ", jcol, (kcolx < krex0 && hcol[kcolx] == jcol) ? "+" : "");
#endif

    if (kcolx < krex0 && hcol[kcolx] == jcol) {
      // before:  both x and y are in the jcol
      // after:   only x is in the jcol
      // so: number of elems in col x unchanged, and num elems in jcol is one less

      // update row rep - just modify coefficent
      // column y is deleted as a whole at the end of the loop
#if	DEBUG_PRESOLVE
      printf("CHANGING %g + %g -> %g\n",
	     rowels[kcolx],
	     rowels[krowy],
	     rowels[kcolx] + rowels[krowy] * coeff_factor);
#endif
      rowels[kcolx] += rowels[krowy] * coeff_factor;

      // this is where this element in rowy ended up
      x_to_y[x_to_y_i++] = kcolx - krsx;
      kcolx++;
    } else {
      // before:  only y is in the jcol
      // after:   only x is in the jcol
      // so: number of elems in col x is one greater, but num elems in jcol remains same
      {
	expand_row(mrstrt, rowels, hcol, hinrow, rlink, nrows, irowx);
	// this may force a compaction
	// this will be called excessively if the rows are packed too tightly

	// have to adjust various induction variables
	krowy = mrstrt[irowy] + (krowy - krs);
	krs = mrstrt[irowy];			// do this for ease of debugging
	kre = mrstrt[irowy] + hinrow[irowy];
	    
	kcolx = mrstrt[irowx] + (kcolx - krsx);	// don't really need to do this
	krex0 = mrstrt[irowx] + (krex0 - krsx);
	krsx = mrstrt[irowx];
	krex = mrstrt[irowx] + hinrow[irowx];
      }

      // this is where this element in rowy ended up
      x_to_y[x_to_y_i++] = krex - krsx;

      // there is now an unused entry in the memory after the column - use it
      // mrstrt[nrows] == penultimate index of arrays hcol/rowels
      hcol[krex] = jcol;
      rowels[krex] = rowels[krowy] * coeff_factor;
      hinrow[irowx]++, krex++;	// expand the col

      // do NOT increment kcolx
    }
  }

#if	DEBUG_PRESOLVE
  printf(")\n");
#endif
}


// It is common in osl to copy from one representation to another
// (say from a col rep to a row rep).
// One such routine is ekkclcp.
// This is similar, except that it does not assume that the
// representation is packed, and it adds some slack space
// in the target rep.
// It assumes both hincol/hinrow are correct.
// Note that such routines automatically sort the target rep by index,
// because they sweep the rows in ascending order.
void copyrep(const int * mrstrt, const int *hcol, const double *rowels, 
	     const int *hinrow, int nrows,
	     int *mcstrt, int *hrow, double *colels, 
	     int *hincol, int ncols)
{
  int pos = 0;
  for (int j = 0; j < ncols; ++j) {
    mcstrt[j] = pos;
    pos += hincol[j];
    pos += min(hincol[j], 10); // slack
    hincol[j] = 0;
  }

  for (int i = 0; i < nrows; ++i) {
    CoinBigIndex krs = mrstrt[i];
    CoinBigIndex kre = krs + hinrow[i];
    for (CoinBigIndex kr = krs; kr < kre; ++kr) {
      int icol = hcol[kr];
      int iput = hincol[icol];
      hincol[icol] = iput + 1;
      iput += mcstrt[icol];
      
      hrow[iput] = i;
      colels[iput] = rowels[kr];
    }
  }
}

#if 0
// variant of function of same name in ekk_7.c
// exactly the same, with some code #if'ed out
int ekk_makeColumnOrdered(int * indexRow , int * indexColumn , double * element ,
			  int * rowCount , int * columnCount , CoinBigIndex * startColumn ,
			  int numberRows , int numberColumns,
			  CoinBigIndex numberElements, double tolerance)
{
  int iColumn,i,k;
#if 0
  for (i=0;i<numberRows;i++) {
    rowCount[i]=0;
  }
  for (i=0;i<numberColumns;i++) {
    columnCount[i]=0;
  }
  for (i=0;i<numberElements;i++) {
    int iRow=indexRow[i];
    int iColumn=indexColumn[i];
    rowCount[iRow]++;
    columnCount[iColumn]++;
  }
#endif

  i=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    /* position after end of Column */
    i+=columnCount[iColumn];
    startColumn[iColumn]=i;
  } /* endfor */
  startColumn[iColumn]=i;

  for (k=numberElements-1;k>=0;k--) {
    iColumn=indexColumn[k];
    if (iColumn>=0) {
      /* pick up the entry with your right hand */
      double value = element[k];
      int iRow=indexRow[k];
      int iColumnSave=0;
      indexColumn[k]=-2;	/* the hole */

      while (1) {
	/* pick this up with your left */
        int iLook=startColumn[iColumn]-1;
        double valueSave=element[iLook];
        int iColumnSave=indexColumn[iLook];
        int iRowSave=indexRow[iLook];

	/* put the right-hand entry where it wanted to go */
        startColumn[iColumn]=iLook;
        element[iLook]=value;
        indexRow[iLook]=iRow;
        indexColumn[iLook]=-1;	/* mark it as being where it wants to be */
	
	/* there was something there */
        if (iColumnSave>=0) {
          iColumn=iColumnSave;
          value=valueSave;
          iRow=iRowSave;
	} else if (iColumnSave = -2)	/* that was the hole */
          break;
	else
	  ekkmesg_no(158);	/* should never happen */
	/* endif */
      } /* endwhile */
    } /* endif */
  } /* endfor */

#if 0
  /* now pack the elements and combine entries with the same row and column */
  /* also, drop entries with "small" coefficients */
  numberElements=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex start=startColumn[iColumn];
    CoinBigIndex end =startColumn[iColumn+1];
    startColumn[iColumn]=numberElements;
    if (end>start) {
      int lastRow;
      double lastValue;
      CoinSort_2(indexRow+start,indexRow+end,element+start);
      //ekk_sort2(indexRow+start,element+start,end-start);
      lastRow=indexRow[start];
      lastValue=element[start];
      for (i=start+1;i<end;i++) {
        int iRow=indexRow[i];
        double value=element[i];
        if (iRow>lastRow) {
          if(fabs(lastValue)>tolerance) {
            indexRow[numberElements]=lastRow;
            element[numberElements]=lastValue;
            numberElements++;
          }
          lastRow=iRow;
          lastValue=value;
        } else {
          lastValue+=value;
        } /* endif */
      } /* endfor */
      if(fabs(lastValue)>tolerance) {
        indexRow[numberElements]=lastRow;
        element[numberElements]=lastValue;
        numberElements++;
      }
    }
  } /* endfor */
#endif

  startColumn[numberColumns]=numberElements;
  return numberElements;
}
#endif


// add -x/y times row y to row x, thus cancelling out one column of rowx;
// afterwards, that col will be singleton for rowy, so we drop the row.
//
// This no longer maintains the col rep as it goes along.
// Instead, it reconstructs it from scratch afterward.
//
// This implements the functionality of ekkrdc3.
const PresolveAction *subst_constraint_action::presolve(PresolveMatrix *prob,
					 char *implied_free,
					const PresolveAction *next,
							int &try_fill_level)
{
  double *colels	= prob->colels_;
  int *hrow	= prob->hrow_;
  CoinBigIndex *mcstrt	= prob->mcstrt_;
  int *hincol	= prob->hincol_;
  const int ncols	= prob->ncols_;

  const double *clo	= prob->clo_;
  const double *cup	= prob->cup_;

  double *rowels	= prob->rowels_;
  int *hcol	= prob->hcol_;
  CoinBigIndex *mrstrt	= prob->mrstrt_;
  int *hinrow	= prob->hinrow_;
  const int nrows	= prob->nrows_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;
  double *acts	= prob->acts_;

  double *dcost		= prob->cost_;

  presolvehlink *clink = prob->clink_;
  presolvehlink *rlink = prob->rlink_;

  const double tol = prob->feasibilityTolerance_;

  action *actions	= new action [ncols];
  int nactions = 0;

  int *zerocols	= new int[ncols];
  int nzerocols	= 0;

  int *x_to_y	= new int[ncols];

#if 0
  // follmer.mps presents a challenge, since it has some very
  // long rows.  I started experimenting with how to deal with it,
  // but haven't yet finished.
  // The idea was to space out the rows to add some padding between them.
  // Ideally, we shouldn't have to do this just here, but could try to
  // do it a little everywhere.

  // sort the row rep by reconstructing from col rep
  copyrep(mcstrt, hrow, colels, hincol, ncols,
	  mrstrt, hcol, rowels, hinrow, nrows);
  presolve_make_memlists(mrstrt, hinrow, rlink, nrows);
  // NEED SOME ASSERTION ABOUT NELEMS

  copyrep(mrstrt, hcol, rowels, hinrow, nrows,
	  mcstrt, hrow, colels, hincol, ncols);
  presolve_make_memlists(mcstrt, hincol, clink, ncols);
#endif

  // in the original presolve, I don't think the two representations were
  // kept in sync.  It may be useful not to do that here, either,
  // but rather just keep the columns with nfill_level rows accurate
  // and resync at the end of the function.

  // DEBUGGING
  //  int nt = 0;
  int ngood = 0;
  int nsubst = 0;
#ifdef	DEBUG_PRESOLVEx
  int maxsubst = atoi(getenv("MAXSUBST"));
#else
  const int maxsubst = 1000000;
#endif

  // This loop does very nearly the same thing as
  // the first loop in implied_free_action::presolve.
  // We have to do it again in case constraints change while we process them (???).
  int numberLook = prob->numberColsToDo_;
  int iLook;
  int * look = prob->colsToDo_;
  int fill_level = try_fill_level;
  int * look2 = NULL;
  // if gone from 2 to 3 look at all
  if (fill_level<0) {
    fill_level=-fill_level;
    try_fill_level=fill_level;
    look2 = new int[ncols];
    look=look2;
    for (iLook=0;iLook<ncols;iLook++) 
      look[iLook]=iLook;
    numberLook=ncols;
 }
 

  for (iLook=0;iLook<numberLook;iLook++) {
    int jcoly=look[iLook];
    if (hincol[jcoly] > 1 && hincol[jcoly] <= fill_level &&
	implied_free[jcoly] == hincol[jcoly]) {
      CoinBigIndex kcs = mcstrt[jcoly];
      CoinBigIndex kce = kcs + hincol[jcoly];

      int bestrowy_size = 0;
      int bestrowy_row=-1;
      int bestrowy_k=-1;
      double bestrowy_coeff=0.0;

      for (CoinBigIndex k=kcs; k<kce; ++k) {
	int row = hrow[k];
	double coeffj = colels[k];

	// we don't clean up zeros in the middle of the routine.
	// if there is one, skip this candidate.
	if (fabs(coeffj) <= ZTOLDP) {
	  bestrowy_size = 0;
	  break;
	}

	// if its row is an equality constraint...
	if (hinrow[row] > 1 &&	// don't bother with singleton rows

	    fabs(rlo[row] - rup[row]) < tol) {
	  CoinBigIndex krs = mrstrt[row];
	  CoinBigIndex kre = krs + hinrow[row];

	  double maxup, maxdown, ilow, iup;

	  implied_bounds(rowels, clo, cup, hcol,
			 krs, kre,
			 &maxup, &maxdown,
			 jcoly, rlo[row], rup[row], &ilow, &iup);

	  if (maxup < PRESOLVE_INF && maxup + tol < rlo[row]) {
	    prob->status_|= 1;
	    prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
					     prob->messages())
					       <<row
					       <<rlo[row]
					       <<rup[row]
					       <<CoinMessageEol;
	    break;
	  } else if (-PRESOLVE_INF < maxdown && rup[row] < maxdown - tol) {
	    prob->status_|= 1;
	    prob->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
					     prob->messages())
					       <<row
					       <<rlo[row]
					       <<rup[row]
					       <<CoinMessageEol;
	    break;
	  } else {
	    // the row has an implied upper or a lower bound 

	    if (clo[jcoly] <= ilow && iup <= cup[jcoly]) {
	      // both column bounds implied by the constraint bounds

	      // we want coeffy to be smaller than x, BACKWARDS from in doubleton
	      if (bestrowy_size == 0 ||
		  fabs(coeffj) > fabs(bestrowy_coeff) ||
		  (fabs(coeffj) == fabs(bestrowy_coeff) &&
		   hinrow[row] < bestrowy_size)) {
		bestrowy_size = hinrow[row];
		bestrowy_row = row;
		bestrowy_coeff = coeffj;
		bestrowy_k = k;
	      }
	    }
	  }
	}
      }

      if (bestrowy_size == 0)
	continue;

      bool all_ok = true;
      for (CoinBigIndex k=kcs; k<kce; ++k) {
	double coeff_factor = fabs(colels[k] / bestrowy_coeff);
	if (fabs(coeff_factor) > 1.3)
	  all_ok = false;
      }

      // check fill-in
      if (all_ok && hincol[jcoly] == 3) {
	// compute fill-in
	int row1 = -1;
	int row2=-1;
	for (CoinBigIndex kk=kcs; kk<kce; ++kk) 
	  if (kk != bestrowy_k) {
	    if (row1 == -1)
	      row1 = hrow[kk];
	    else
	      row2 = hrow[kk];
	  }


	CoinBigIndex krs = mrstrt[bestrowy_row];
	CoinBigIndex kre = krs + hinrow[bestrowy_row];
	CoinBigIndex krs1 = mrstrt[row1];
	CoinBigIndex kre1 = krs + hinrow[row1];
	CoinBigIndex krs2 = mrstrt[row2];
	CoinBigIndex kre2 = krs + hinrow[row2];

	CoinSort_2(hcol+krs,hcol+krs+hinrow[bestrowy_row],rowels+krs);
	CoinSort_2(hcol+krs1,hcol+krs1+hinrow[row1],rowels+krs1);
	CoinSort_2(hcol+krs2,hcol+krs2+hinrow[row2],rowels+krs2);
	//ekk_sort2(hcol+krs,  rowels+krs,  hinrow[bestrowy_row]);
	//ekk_sort2(hcol+krs1, rowels+krs1, hinrow[row1]);
	//ekk_sort2(hcol+krs2, rowels+krs2, hinrow[row2]);

	int nfill = -hinrow[bestrowy_row];
	CoinBigIndex kcol1 = krs1;
	CoinBigIndex kk;
	for (kk=krs; kk<kre; ++kk) {
	  int jcol = hcol[kk];

	  while (kcol1 < kre1 && hcol[kcol1] < jcol)
	    kcol1++;
	  if (! (kcol1 < kre1 && hcol[kcol1] == jcol))
	    nfill++;
	}
	CoinBigIndex kcol2 = krs2;
	for (kk=krs; kk<kre; ++kk) {
	  int jcol = hcol[kk];

	  while (kcol2 < kre2 && hcol[kcol2] < jcol)
	    kcol2++;
	  if (! (kcol2 < kre2 && hcol[kcol2] == jcol))
	    nfill++;
	}
#if	DEBUG_PRESOLVE
	printf("FILL:  %d\n", nfill);
#endif

#if 0
	static int maxfill = atoi(getenv("MAXFILL"));

	if (nfill > maxfill)
	  all_ok = false;
#endif

	// not too much
	if (nfill <= 0)
	  ngood++;

#if 0
	static int nts = 0;
	if (++nts > atoi(getenv("NTS")))
	  all_ok = false;
	else
	  nt++;
#endif
      }

      // probably never happens
      if (all_ok && nzerocols + hinrow[bestrowy_row] >= ncols)
	all_ok = false;

      if (nsubst >= maxsubst) {
	all_ok = false;
      } 

      if (all_ok) {
	nsubst++;
#if 0
	// debug
	if (numberLook<ncols&&iLook==numberLook-1) {
	  printf("found last one?? %d\n", jcoly);
	}
#endif

	CoinBigIndex kcs = mcstrt[jcoly];
	int rowy = bestrowy_row;
	double coeffy = bestrowy_coeff;

	PRESOLVEASSERT(fabs(colels[kcs]) > ZTOLDP);
	PRESOLVEASSERT(fabs(colels[kcs+1]) > ZTOLDP);

	PRESOLVEASSERT(hinrow[rowy] > 1);

	const bool nonzero_cost = (fabs(dcost[jcoly]) > tol);

	double *costsx = nonzero_cost ? new double[hinrow[rowy]] : 0;

	int ntotels = 0;
	for (CoinBigIndex k=kcs; k<kce; ++k) {
	  int irow = hrow[k];
	  ntotels += hinrow[irow];
	}

	{
	  action *ap = &actions[nactions++];
	  int nincol = hincol[jcoly];

	  ap->col = jcoly;
	  ap->rowy = rowy;

	  ap->nincol = nincol;
	  ap->rows = new int[nincol];
	  ap->rlos = new double[nincol];
	  ap->rups = new double[nincol];

	  // coefficients in deleted col
	  ap->coeffxs = new double[nincol];

	  ap->ninrowxs = new int[nincol];
	  ap->rowcolsxs = new int[ntotels];
	  ap->rowelsxs = new double[ntotels];

	  ap->costsx = costsx;

	  // copy all the rows for restoring later - wasteful
	  {
	    int nel = 0;
	    for (CoinBigIndex k=kcs; k<kce; ++k) {
	      int irow = hrow[k];
	      CoinBigIndex krs = mrstrt[irow];
	      //	      CoinBigIndex kre = krs + hinrow[irow];

	      prob->addRow(irow);
	      ap->rows[k-kcs] = irow;
	      ap->ninrowxs[k-kcs] = hinrow[irow];
	      ap->rlos[k-kcs] = rlo[irow];
	      ap->rups[k-kcs] = rup[irow];

	      ap->coeffxs[k-kcs] = colels[k];

	      memcpy( &ap->rowcolsxs[nel], &hcol[krs],hinrow[irow]*sizeof(int));
	      memcpy( &ap->rowelsxs[nel], &rowels[krs],hinrow[irow]*sizeof(double));
	      nel += hinrow[irow];
	    }
	  }
	}

	// rowy is supposed to be an equality row
	PRESOLVEASSERT(fabs(rup[rowy] - rlo[rowy]) < ZTOLDP);

	// now adjust for the implied free row - COPIED
	if (nonzero_cost) {
#ifdef	DEBUG_PRESOLVE
	  printf("NONZERO SUBST COST:  %d %g\n", jcoly, dcost[jcoly]);
#endif
	  double *cost = dcost;
	  double *save_costs = costsx;
	  double coeffj = coeffy;
	  CoinBigIndex krs = mrstrt[rowy];
	  CoinBigIndex kre = krs + hinrow[rowy];

	  double rhs = rlo[rowy];
	  double costj = cost[jcoly];

	  for (CoinBigIndex k=krs; k<kre; k++) {
	    int jcol = hcol[k];
	    prob->addCol(jcol);
	    save_costs[k-krs] = cost[jcol];

	    if (jcol != jcoly) {
	      double coeff = rowels[k];

	      /*
	       * Similar to eliminating doubleton:
	       *   cost1 x = cost1 (c - b y) / a = (c cost1)/a - (b cost1)/a
	       *   cost[icoly] += cost[icolx] * (-coeff2 / coeff1);
	       */
	      cost[jcol] += costj * (-coeff / coeffj);
	    }
	  }

	  // I'm not sure about this
	  prob->change_bias(costj * rhs / coeffj);

	  // ??
	  cost[jcoly] = 0.0;
	}

#if	DEBUG_PRESOLVE
	    if (hincol[jcoly] == 3) {
	      CoinBigIndex krs = mrstrt[rowy];
	      CoinBigIndex kre = krs + hinrow[rowy];
	      printf("HROW0 (%d):  ", rowy);
	      for (CoinBigIndex k=krs; k<kre; ++k) {
		int jcol = hcol[k];
		double coeff = rowels[k];
		printf("%d:%g (%d) ", jcol, coeff, hincol[jcol]);
	      }
	      printf("\n");
	    }
#endif

	    if (hincol[jcoly] != 2) {
	      CoinBigIndex krs = mrstrt[rowy];
	      //	      CoinBigIndex kre = krs + hinrow[rowy];
	      CoinSort_2(hcol+krs,hcol+krs+hinrow[rowy],rowels+krs);
	      //ekk_sort2(hcol+krs,  rowels+krs,  hinrow[rowy]);
	    }

	    // substitute away jcoly in the other rows
	    // Use ap as mcstrt etc may move if compacted
	    kce = hincol[jcoly];
	    CoinBigIndex k;
	    action *ap = &actions[nactions-1];
	    for (k=0; k<kce; ++k) {
	      int rowx = ap->rows[k];
	      //assert(rowx==hrow[k+kcs]);
	      //assert(ap->coeffxs[k]==colels[k+kcs]);
	      if (rowx != rowy) {
		double coeffx = ap->coeffxs[k];
		double coeff_factor = -coeffx / coeffy;	// backwards from doubleton

#if	DEBUG_PRESOLVE
		{
		  CoinBigIndex krs = mrstrt[rowx];
		  CoinBigIndex kre = krs + hinrow[rowx];
		  printf("HROW (%d %d %d):  ", rowx, hinrow[rowx], jcoly);
		  for (CoinBigIndex k=krs; k<kre; ++k) {
		    int jcol = hcol[k];
		    double coeff = rowels[k];
		    printf("%d ", jcol);
		  }
		  printf("\n");
#if 0
		  for (CoinBigIndex k=krs; k<kre; ++k) {
		    int jcol = hcol[k];
		    prob->addCol(jcol);
		    double coeff = rowels[k];
		    printf("%g ", coeff);
		  }
		  printf("\n");
#endif
		}
#endif
		{
		  CoinBigIndex krsx = mrstrt[rowx];
		  CoinBigIndex krex = krsx + hinrow[rowx];
		  int i;
		  for (i=krsx;i<krex;i++) 
		    prob->addCol(hcol[i]);
		  if (hincol[jcoly] != 2) 
		    CoinSort_2(hcol+krsx,hcol+krsx+hinrow[rowx],rowels+krsx);
		  //ekk_sort2(hcol+krsx, rowels+krsx, hinrow[rowx]);
		}
		
		// add (coeff_factor * <rowy>) to rowx
		// does not affect rowy
		// may introduce (or cancel) elements in rowx
		add_row(mrstrt,
			rlo, acts, rup,
			rowels, hcol,
			hinrow,
			rlink, nrows,
			coeff_factor,
			rowx, rowy,
			x_to_y);
		
		// update col rep of rowx from row rep:
		// for every col in rowy, copy the elem for that col in rowx
		// from the row rep to the col rep
		{
		  CoinBigIndex krs = mrstrt[rowy];
		  //		  CoinBigIndex kre = krs + hinrow[rowy];
		  int niny = hinrow[rowy];
		  
		  CoinBigIndex krsx = mrstrt[rowx];
		  //		  CoinBigIndex krex = krsx + hinrow[rowx];
		  for (CoinBigIndex ki=0; ki<niny; ++ki) {
		    CoinBigIndex k = krs + ki;
		    int jcol = hcol[k];
		    prob->addCol(jcol);
		    CoinBigIndex kcs = mcstrt[jcol];
		    CoinBigIndex kce = kcs + hincol[jcol];
		    
		    //double coeff = rowels[presolve_find_row(jcol, krsx, krex, hcol)];
		    if (hcol[krsx + x_to_y[ki]] != jcol)
		      abort();
		    double coeff = rowels[krsx + x_to_y[ki]];
		    
		    // see if rowx appeared in jcol in the col rep
		    CoinBigIndex k2 = presolve_find_row1(rowx, kcs, kce, hrow);
		    
		    //PRESOLVEASSERT(fabs(coeff) > ZTOLDP);
		    
		    if (k2 < kce) {
		      // yes - just update the entry
		      colels[k2] = coeff;
		    } else {
		      // no - make room, then append
		      expand_row(mcstrt, colels, hrow, hincol, clink, ncols, jcol);
		      kcs = mcstrt[jcol];
		      kce = kcs + hincol[jcol];
		      
		      hrow[kce] = rowx;
		      colels[kce] = coeff;
		      hincol[jcol]++;
		    }
		  }
		}
		// now colels[k] == 0.0

#if 1
		// now remove jcoly from rowx in the row rep
		// better if this were first
		presolve_delete_from_row(rowx, jcoly, mrstrt, hinrow, hcol, rowels);
#endif
#if	DEBUG_PRESOLVE
		{
		  CoinBigIndex krs = mrstrt[rowx];
		  CoinBigIndex kre = krs + hinrow[rowx];
		  printf("HROW (%d %d %d):  ", rowx, hinrow[rowx], jcoly);
		  for (CoinBigIndex k=krs; k<kre; ++k) {
		    int jcol = hcol[k];
		    double coeff = rowels[k];
		    printf("%d ", jcol);
		  }
		  printf("\n");
#if 0
		  for (CoinBigIndex k=krs; k<kre; ++k) {
		    int jcol = hcol[k];
		    double coeff = rowels[k];
		    printf("%g ", coeff);
		  }
		  printf("\n");
#endif
		}
#endif
		
		// don't have to update col rep, since entire col deleted later
	      }
	    }

#if	DEBUG_PRESOLVE
	    printf("\n");
#endif

	    // the addition of rows may have created zero coefficients
	    memcpy( &zerocols[nzerocols], &hcol[mrstrt[rowy]],hinrow[rowy]*sizeof(int));
	    nzerocols += hinrow[rowy];
	    
	    // delete rowy in col rep
	    {
	      CoinBigIndex krs = mrstrt[rowy];
	      CoinBigIndex kre = krs + hinrow[rowy];
	      for (CoinBigIndex k=krs; k<kre; ++k) {
		int jcol = hcol[k];
		
		// delete rowy from the jcol
		presolve_delete_from_row(jcol, rowy, mcstrt, hincol, hrow, colels);
	      }
	    }
	    // delete rowy in row rep
	    hinrow[rowy] = 0;
	    
	    // This last is entirely dual to doubleton, but for the cost adjustment
	    
	    // eliminate col entirely from the col rep
	    PRESOLVE_REMOVE_LINK(clink, jcoly);
	    hincol[jcoly] = 0;
	    
	    // eliminate rowy entirely from the row rep
	    PRESOLVE_REMOVE_LINK(rlink, rowy);
	    //cost[irowy] = 0.0;
	    
	    rlo[rowy] = 0.0;
	    rup[rowy] = 0.0;
	    
#if	0 && DEBUG_PRESOLVE
	    printf("ROWY COLS:  ");
	    for (CoinBigIndex k=0; k<save_ninrowy; ++k)
	      if (rowycols[k] != col) {
		printf("%d ", rowycols[k]);
		(void)presolve_find_row(rowycols[k], mrstrt[rowx], mrstrt[rowx]+hinrow[rowx],
					hcol);
	      }
	    printf("\n");
#endif
#if 0
	presolve_links_ok(clink, mcstrt, hincol, ncols);
	presolve_links_ok(rlink, mrstrt, hinrow, nrows);
	prob->consistent();
#endif
      }
    }
  }

  // general idea - only do doubletons until there are almost none left
  if (nactions < 30&&fill_level==2)
    try_fill_level = -3;

  if (nactions) {
#if	PRESOLVE_SUMMARY
    printf("NSUBSTS:  %d\n", nactions);
    printf("NT: %d  NGOOD:  %d FILL_LEVEL:  %d\n", nt, ngood, fill_level);
#endif
    next = new subst_constraint_action(nactions, copyOfArray(actions,nactions), next);

    next = drop_zero_coefficients_action::presolve(prob, zerocols, nzerocols, next);
  }
  delete [] look2;
  deleteAction(actions);

  delete[]x_to_y;
  delete[]zerocols;

  return (next);
}


void subst_constraint_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions	= nactions_;

  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *link		= prob->link_;
  //  int ncols		= prob->ncols_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double *dcost	= prob->cost_;

  double *sol	= prob->sol_;
  double *rcosts	= prob->rcosts_;

  double *acts	= prob->acts_;
  double *rowduals = prob->rowduals_;

  char *cdone	= prob->cdone_;
  char *rdone	= prob->rdone_;
  CoinBigIndex free_list = prob->free_list_;

  const double ztolzb	= prob->ztolzb_;
  //  const double ztoldj	= prob->ztoldj_;
  const double maxmin = prob->maxmin_;
  int k;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {
    int icol = f->col;

    int nincoly = f->nincol;
    double *rlos = f->rlos;
    double *rups = f->rups;
    int *rows = f->rows;

    double *coeffxs = f->coeffxs;

    int jrowy = f->rowy;

    int *ninrowxs = f->ninrowxs;
    const int *rowcolsxs = f->rowcolsxs;
    const double *rowelsxs = f->rowelsxs;
    
    /* the row was in the reduced problem */
    for (int i=0; i<nincoly; ++i) {
      if (rows[i] != jrowy)
	PRESOLVEASSERT(rdone[rows[i]]);
    }
    PRESOLVEASSERT(cdone[icol]==DROP_COL);
    PRESOLVEASSERT(rdone[jrowy]==DROP_ROW);

    // DEBUG CHECK
#if	0 && DEBUG_PRESOLVE
    {
      double actx = 0.0;
      for (int j=0; j<ncols; ++j)
	if (hincol[j] > 0 && cdone[j]) {
	  CoinBigIndex krow = presolve_find_row1(jrowx, mcstrt[j], mcstrt[j] + hincol[j], hrow);
	  if (krow < mcstrt[j] + hincol[j])
	    actx += colels[krow] * sol[j];
      }
      if (! (fabs(acts[jrowx] - actx) < 100*ztolzb))
	printf("BAD ACTSX:  acts[%d]==%g != %g\n",
	       jrowx, acts[jrowx], actx);
      if (! (rlo[jrowx] - 100*ztolzb <= actx && actx <= rup[jrowx] + 100*ztolzb))
	printf("ACTSX NOT IN RANGE:  %d %g %g %g\n",
	       jrowx, rlo[jrowx], actx, rup[jrowx]);
    }
#endif

    int ninrowy=-1;
    const int *rowcolsy=NULL;
    const double *rowelsy=NULL;
    double coeffy=0.0;

    double rloy=1.0e50;
    {
      int nel = 0;
      for (int i=0; i<nincoly; ++i) {
	int row = rows[i];
	rlo[row] = rlos[i];
	rup[row] = rups[i];
	if (row == jrowy) {
	  ninrowy = ninrowxs[i];
	  rowcolsy = &rowcolsxs[nel];
	  rowelsy  = &rowelsxs[nel];

	  coeffy = coeffxs[i];
	  rloy = rlo[row];

	}
	nel += ninrowxs[i];
      }
    }
    double rhsy = rloy;

    // restore costs
    {
      const double *costs = f->costsx;
      if (costs)
	for (int i = 0; i<ninrowy; ++i) {
	  dcost[rowcolsy[i]] = costs[i];
	}
    }

    // solve for the equality to find the solution for the eliminated col
    // this is why we want coeffx < coeffy (55)
    {
      double sol0 = rloy;
      sol[icol] = 0.0;	// to avoid condition in loop
      for (k = 0; k<ninrowy; ++k) {
	int jcolx = rowcolsy[k];
	double coeffx = rowelsy[k];
	sol0 -= coeffx * sol[jcolx];
      }
      sol[icol] = sol0 / coeffy;

      if (! (sol[icol] > clo[icol] - ztolzb &&
	     cup[icol] + ztolzb > sol[icol]))
	printf("NEW SOL OUT-OF-TOL:  %g %g %g\n", clo[icol], sol[icol], cup[icol]);
    }

    // since this row is fixed 
    acts[jrowy] = rloy;

    // acts[irow] always ok, since slack is fixed
    prob->setRowStatus(jrowy,PrePostsolveMatrix::atLowerBound);

    // remove old rowx from col rep
    // we don't explicitly store what the current rowx is;
    // however, after the presolve, rowx contains a col for every
    // col in either the original rowx or the original rowy.
    // If there were cancellations, those were handled in subsequent
    // presolves.
    {
      // erase those cols in the other rows that occur in rowy
      // (with the exception of icol, which was deleted);
      // the other rows *must* contain these cols
      for (k = 0; k<ninrowy; ++k) {
	int col = rowcolsy[k];

	// remove jrowx from col in the col rep
	// icol itself was deleted, so won't be there
	if (col != icol)
	  for (int i = 0; i<nincoly; ++i) {
	    if (rows[i] != jrowy)
	      presolve_delete_from_row2(col, rows[i], mcstrt, hincol, hrow, colels, link, &free_list);
	  }
      }

      // initialize this for loops below
      hincol[icol] = 0;

      // now restore the original rows (other than rowy).
      // those cols that were also in rowy were just removed;
      // otherwise, they *must* already be there.
      // This loop and the next automatically create the rep for the new col.
      {
	const int *rowcolsx = rowcolsxs;
	const double *rowelsx = rowelsxs;

	for (int i = 0; i<nincoly; ++i) {
	  int ninrowx = ninrowxs[i];
	  int jrowx = rows[i];

	  if (jrowx != jrowy)
	    for (k = 0; k<ninrowx; ++k) {
	      int col = rowcolsx[k];
	      CoinBigIndex kcolx = presolve_find_row3(jrowx, mcstrt[col], hincol[col], hrow, link);

	      if (kcolx != -1) {
		PRESOLVEASSERT(presolve_find_row1(col, 0, ninrowy, rowcolsy) == ninrowy);
		// overwrite the existing entry
		colels[kcolx] = rowelsx[k];
	      } else {
		PRESOLVEASSERT(presolve_find_row1(col, 0, ninrowy, rowcolsy) < ninrowy);

		{
		  CoinBigIndex kk = free_list;
		  free_list = link[free_list];
		  check_free_list(free_list);

		  link[kk] = mcstrt[col];
		  mcstrt[col] = kk;
		  colels[kk] = rowelsx[k];
		  hrow[kk] = jrowx;
		}
		++hincol[col];
	      }
	    }
	  rowcolsx += ninrowx;
	  rowelsx += ninrowx;
	}
      }

      // finally, add original rowy elements
      for (k = 0; k<ninrowy; ++k) {
	int col = rowcolsy[k];

	{
	  prepend_elem(col, rowelsy[k], jrowy, mcstrt, colels, hrow, link, &free_list);
	  ++hincol[col];
	}
      }
    }

    // my guess is that the CLAIM in doubleton generalizes to
    // equations with more than one x-style variable.
    // Since I can't see how to distinguish among them,
    // I assume that any of them will do.

    {
       //      CoinBigIndex k;
      double dj = maxmin*dcost[icol];
      double bounds_factor = rhsy/coeffy;
      for (int i=0; i<nincoly; ++i)
	if (rows[i] != jrowy) {
	  int row = rows[i];
	  double coeff = coeffxs[i];

	  // PROBABLY DOESN'T MAKE SENSE
	  acts[row] += coeff * bounds_factor;

	  dj -= rowduals[row] * coeff;
	}

      // DEBUG CHECK
      double acty = 0.0;
      for (k = 0; k<ninrowy; ++k) {
	int col = rowcolsy[k];
	acty += rowelsy[k] * sol[col];
      }
      PRESOLVEASSERT(fabs(acty - acts[jrowy]) < 100*ZTOLDP);

      // RECOMPUTING
      {
	const int *rowcolsx = rowcolsxs;
	const double *rowelsx = rowelsxs;

	for (int i=0; i<nincoly; ++i) {
	  int ninrowx = ninrowxs[i];

	  if (rows[i] != jrowy) {
	    int jrowx = rows[i];

	    double actx = 0.0;
	    for (k = 0; k<ninrowx; ++k) {
	      int col = rowcolsx[k];
	      actx += rowelsx[k] * sol[col];
	    }
#if	DEBUG_PRESOLVE
	    PRESOLVEASSERT(rlo[jrowx] - ztolzb <= actx && actx <= rup[jrowx] + ztolzb);
#endif
	    acts[jrowx] = actx;
	  }
	  rowcolsx += ninrowx;
	  rowelsx += ninrowx;
	}
      }

      // this is the coefficient we need to force col y's reduced cost to 0.0;
      // for example, this is obviously true if y is a singleton column
      rowduals[jrowy] = dj / coeffy;
      rcosts[icol] = 0.0;

      // furthermore, this should leave rcosts[colx] for all colx
      // in jrowx unchanged (????).
    }

    // Unlike doubleton, there should never be a problem with keeping
    // the reduced costs the way they were, because the other
    // variable's bounds are never changed, since col was implied free.
    //rowstat[jrowy] = 0;
    prob->setColumnStatus(icol,PrePostsolveMatrix::basic);

    cdone[icol] = SUBST_ROW;
    rdone[jrowy] = SUBST_ROW;
  }

  prob->free_list_ = free_list;
}



subst_constraint_action::~subst_constraint_action()
{
  const action *actions = actions_;

  for (int i=0; i<nactions_; ++i) {
    delete[]actions[i].rows;
    delete[]actions[i].rlos;
    delete[]actions[i].rups;
    delete[]actions[i].coeffxs;
    delete[]actions[i].ninrowxs;
    delete[]actions[i].rowcolsxs;
    delete[]actions[i].rowelsxs;
    delete[]actions[i].costsx;
  }

  deleteAction(actions_);
}
