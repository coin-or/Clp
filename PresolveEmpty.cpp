// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "PresolveMatrix.hpp"

#include "PresolveEmpty.hpp"
#include "ClpMessage.hpp"


const PresolveAction *drop_empty_cols_action::presolve(PresolveMatrix *prob,
							int *ecols,
							int necols,
							const PresolveAction *next)
{
  int ncols		= prob->ncols_;
  int *mcstrt		= prob->mcstrt_;
  int *hincol		= prob->hincol_;
  int *hcol		= prob->hcol_;

  // We know we are not going to need row copy again
  //presolvehlink *clink	= prob->clink_;

  double *clo		= prob->clo_;
  double *cup		= prob->cup_;
  double *dcost		= prob->cost_;

  //int *mrstrt		= prob->mrstrt_;
  //int *hinrow		= prob->hinrow_;
  //int *hrow		= prob->hrow_;

  char * integerType     = prob->integerType_;
  int * originalColumn  = prob->originalColumn_;

  const double maxmin	= prob->maxmin_;

  double * sol = prob->sol_;
  unsigned char * colstat = prob->colstat_;

  action *actions 	= new action[necols];
  int * colmapping = new int [ncols];

  memset(colmapping,0,ncols*sizeof(int));
  int i;
  for (i=necols-1; i>=0; i--) {
    int jcol = ecols[i];
    colmapping[jcol]=-1;
    action &e = actions[i];

    e.jcol	= jcol;
    e.clo	= clo[jcol];
    e.cup	= cup[jcol];
    e.cost	= dcost[jcol];

    // there are no more constraints on this variable, 
    // so we had better be able to compute the answer now

    if (dcost[jcol] * maxmin == 0.0) {
      // hopefully, we can make this non-basic
      // what does OSL currently do in this case???
      e.sol = (-PRESOLVE_INF < clo[jcol]
	       ? clo[jcol]
	       : cup[jcol] < PRESOLVE_INF
	       ? cup[jcol]
	       : 0.0);
    } else if (dcost[jcol] * maxmin > 0.0) {
      if (-PRESOLVE_INF < clo[jcol])
	e.sol = clo[jcol];
      else {
	  prob->originalModel_->messageHandler()->message(CLP_PRESOLVE_COLUMNBOUNDB,
					     prob->originalModel_->messages())
					       <<jcol
					       <<CoinMessageEol;
	prob->status_ |= 2;
	break;
      }
    } else {
      if (cup[jcol] < PRESOLVE_INF)
	e.sol = cup[jcol];
      else {
	  prob->originalModel_->messageHandler()->message(CLP_PRESOLVE_COLUMNBOUNDA,
					     prob->originalModel_->messages())
					       <<jcol
					       <<CoinMessageEol;
	prob->status_ |= 2;
	break;
      }
    }
	
#if	DEBUG_PRESOLVE
    if (e.sol * dcost[jcol]) {
      //printf("\a>>>NON-ZERO COST FOR EMPTY COL %d:  %g\n", jcol, dcost[jcol]);
    }
#endif
    prob->change_bias(e.sol * dcost[jcol]);


  }
  int ncols2=0;

  // now move remaining ones down
  for (i=0;i<ncols;i++) {
    if (!colmapping[i]) {
      mcstrt[ncols2] = mcstrt[i];
      hincol[ncols2] = hincol[i];
    
      clo[ncols2]   = clo[i];
      cup[ncols2]   = cup[i];

      dcost[ncols2] = dcost[i];
      if (sol) {
	sol[ncols2] = sol[i];
	colstat[ncols2] = colstat[i];
      }

      integerType[ncols2] = integerType[i];
      originalColumn[ncols2] = originalColumn[i];
      ncols2++;
    }
  }
  
  delete [] colmapping;
  prob->ncols_ = ncols2;

  return (new drop_empty_cols_action(necols, actions, next));
}


const PresolveAction *drop_empty_cols_action::presolve(PresolveMatrix *prob,
							  const PresolveAction *next)
{
  const int *hincol	= prob->hincol_;
  const double *clo	= prob->clo_;
  const double *cup	= prob->cup_;
  int ncols		= prob->ncols_;
  int i;
  int nempty		= 0;
  int * empty = new int [ncols];

  // count empty cols
  for (i=0; i<ncols; i++)
    if (hincol[i] == 0) {
#if	DEBUG_PRESOLVE
      if (nempty==0)
	printf("UNUSED COLS:  ");
      printf("%d ", i);
#endif
      empty[nempty++] = i;
    }

  if (nempty) {
#if	DEBUG_PRESOLVE
      printf("\ndropped %d cols\n", nempty);
#endif
    next = drop_empty_cols_action::presolve(prob, empty, nempty, next);
  }
  delete [] empty;
  return (next);
}


void drop_empty_cols_action::postsolve(PostsolveMatrix *prob) const
{
  const int nactions	= nactions_;
  const action *const actions = actions_;

  int ncols		= prob->ncols_;
  int nelems		= prob->nelems_;

  int *mcstrt	= prob->mcstrt_;
  int *hincol	= prob->hincol_;
  int *hrow	= prob->hrow_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  double *sol	= prob->sol_;
  double *cost	= prob->cost_;
  double *rcosts	= prob->rcosts_;
  unsigned char *colstat	= prob->colstat_;
  const double maxmin = prob->maxmin_;

  int ncols2 = ncols+nactions;
  int * colmapping = new int [ncols2];

  memset(colmapping,0,ncols2*sizeof(int));
  char *cdone	= prob->cdone_;

  for (int action_i = 0; action_i < nactions; action_i++) {
    const action *e = &actions[action_i];
    int jcol = e->jcol;
    colmapping[jcol]=-1;
  }

  int i;

  // now move remaining ones up
  for (i=ncols2-1;i>=0;i--) {
    if (!colmapping[i]) {
      ncols--;
      mcstrt[i] = mcstrt[ncols];
      hincol[i] = hincol[ncols];
    
      clo[i]   = clo[ncols];
      cup[i]   = cup[ncols];

      cost[i] = cost[ncols];

      sol[i] = sol[ncols];
      
      rcosts[i] = rcosts[ncols];
      
      if (colstat) 
	colstat[i] = colstat[ncols];
      cdone[i]   = cdone[ncols];
    }
  }
  assert (!ncols);
  
  delete [] colmapping;

  for (int action_i = 0; action_i < nactions; action_i++) {
    const action *e = &actions[action_i];
    int jcol = e->jcol;
    
    // now recreate jcol
    clo[jcol] = e->clo;
    cup[jcol] = e->cup;
    sol[jcol] = e->sol;
    cost[jcol] = e->cost;

    rcosts[jcol] = maxmin*cost[jcol];

    hincol[jcol] = 0;

    cdone[jcol] = DROP_COL;
    if (colstat) 
      prob->setColumnStatusUsingValue(jcol);
  }

  prob->ncols_ += nactions;
}



const PresolveAction *drop_empty_rows_action::presolve(PresolveMatrix *prob,
				       const PresolveAction *next)
{
  int ncols	= prob->ncols_;
  int *mcstrt	= prob->mcstrt_;
  int *hincol	= prob->hincol_;
  int *hrow	= prob->hrow_;

  int nrows	= prob->nrows_;
  // This is done after row copy needed
  //int *mrstrt	= prob->mrstrt_;
  int *hinrow	= prob->hinrow_;
  //int *hcol	= prob->hcol_;
  
  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  unsigned char *rowstat	= prob->rowstat_;
  double *acts	= prob->acts_;

  //presolvehlink *rlink = prob->rlink_;
  

  int i;
  int nactions = 0;
  for (i=0; i<nrows; i++)
    if (hinrow[i] == 0)
      nactions++;

  if (nactions == 0)
    return next;
  else {
    action *actions 	= new action[nactions];
    int * rowmapping = new int [nrows];

    nactions = 0;
    int nrows2=0;
    for (i=0; i<nrows; i++) {
      if (hinrow[i] == 0) {
	action &e = actions[nactions];
	nactions++;

#if	DEBUG_PRESOLVE
	if (nactions==1)
	  printf("unused rows:  ");

	printf("%d ", i);
#endif
	if (rlo[i] > 0.0 || rup[i] < 0.0) {
	  if (fabs(rlo[i])<=prob->feasibilityTolerance_ &&
	      fabs(rup[i])<=prob->feasibilityTolerance_) {
	    rlo[i]=0.0;
	    rup[i]=0.0;
	  } else {
	    prob->status_|= 1;
	  prob->originalModel_->messageHandler()->message(CLP_PRESOLVE_ROWINFEAS,
					     prob->originalModel_->messages())
					       <<i
					       <<rlo[i]
					       <<rup[i]
					       <<CoinMessageEol;
	    break;
	  }
	}
	e.row	= i;
	e.rlo	= rlo[i];
	e.rup	= rup[i];
	rowmapping[i]=-1;

      } else {
	// move down - we want to preserve order
	rlo[nrows2]=rlo[i];
	rup[nrows2]=rup[i];
	if (acts) {
	  acts[nrows2]=acts[i];
	  rowstat[nrows2]=rowstat[i];
	}
	rowmapping[i]=nrows2++;
      }
    }

    // remap matrix
    for (i=0;i<ncols;i++) {
      int j;
      for (j=mcstrt[i];j<mcstrt[i]+hincol[i];j++) 
	hrow[j] = rowmapping[hrow[j]];
    }
    delete [] rowmapping;

    prob->nrows_ = nrows2;

#if	DEBUG_PRESOLVE
    if (nactions)
      printf("\ndropped %d rows\n", nactions);
#endif
    return (new drop_empty_rows_action(nactions, actions, next));
  }
}

void drop_empty_rows_action::postsolve(PostsolveMatrix *prob) const
{
  const int nactions	= nactions_;
  const action *const actions = actions_;

  int ncols	= prob->ncols_;
  int *mcstrt	= prob->mcstrt_;
  int *hincol	= prob->hincol_;

  int *hrow	= prob->hrow_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;
  unsigned char *rowstat	= prob->rowstat_;
  double *rowduals = prob->rowduals_;
  double *acts	= prob->acts_;
  char *rdone	= prob->rdone_;

  int nrows0	= prob->nrows0_;
  int nrows	= prob->nrows_;

  int * rowmapping = new int [nrows0];

  int i, action_i;
  for (i=0; i<nrows0; i++)
    rowmapping[i] = 0;
  
  for (action_i = 0; action_i<nactions; action_i++) {
    const action *e = &actions[action_i];
    int hole = e->row;
    rowmapping[hole]=-1;
  }

  // move data
  for (i=nrows0-1; i>=0; i--) {
    if (!rowmapping[i]) {
      // not a hole
      nrows--;
      rlo[i]=rlo[nrows];
      rup[i]=rup[nrows];
      acts[i]=acts[nrows];
      rowduals[i]=rowduals[nrows];
      if (rowstat)
	rowstat[i] = rowstat[nrows];
    }
  }
  assert (!nrows);
  // set up mapping for matrix
  for (i=0;i<nrows0;i++) {
    if (!rowmapping[i])
      rowmapping[nrows++]=i;
  }

  for (int j=0; j<ncols; j++) {
    int start = mcstrt[j];
    int end   = start + hincol[j];
    
    for (int k=start; k<end; ++k) {
      hrow[k] = rowmapping[hrow[k]];
    }
  }


  delete [] rowmapping;

  for (action_i = 0; action_i < nactions; action_i++) {
    const action *e = &actions[action_i];
    int irow = e->row;

    // Now recreate irow
    rlo[irow] = e->rlo;
    rup[irow] = e->rup;

    if (rowstat)
      prob->setRowStatus(irow,PrePostsolveMatrix::basic);
    rowduals[irow] = 0.0;	// ???
    acts[irow] = 0.0;

    // hinrow is not maintained during postsolve
    //hinrow[irow] = 0;

    rdone[irow] = DROP_ROW;
  }

  prob->nrows_ = prob->nrows_+nactions;
}

