#include <stdio.h>
#include <math.h>
#include "PresolveMatrix.hpp"
#include "PresolveUseless.hpp"


// WHAT HAPPENS IF COLS ARE DROPPED AS A RESULT??
// should be like do_tighten.
// not really - one could fix costed variables to appropriate bound.
// ok, don't bother about it.  If it is costed, it will be checked
// when it is eliminated as an empty col; if it is costed in the
// wrong direction, the problem is unbounded, otherwise it is pegged
// at its bound.  no special action need be taken here.
const PresolveAction *useless_constraint_action::presolve(PresolveMatrix * prob,
								  const int *useless_rows,
								  int nuseless_rows,
				       const PresolveAction *next)
{
  // may be modified by useless constraint
  double *colels	= prob->colels_;

  // may be modified by useless constraint
        int *hrow	= prob->hrow_;

  const CoinBigIndex *mcstrt	= prob->mcstrt_;

  // may be modified by useless constraint
        int *hincol	= prob->hincol_;

	//  double *clo	= prob->clo_;
	//  double *cup	= prob->cup_;

  const double *rowels	= prob->rowels_;
  const int *hcol	= prob->hcol_;
  const CoinBigIndex *mrstrt	= prob->mrstrt_;

  // may be written by useless constraint
        int *hinrow	= prob->hinrow_;
	//  const int nrows	= prob->nrows_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  action *actions	= new action [nuseless_rows];

#if	PRESOLVE_SUMMARY
    printf("NUSELESS ROWS:  %d\n", nuseless_rows);
#endif

  for (int i=0; i<nuseless_rows; ++i) {
    int irow = useless_rows[i];
    CoinBigIndex krs = mrstrt[irow];
    CoinBigIndex kre = krs + hinrow[irow];

    action *f = &actions[i];

    f->row = irow;
    f->ninrow = hinrow[irow];
    f->rlo = rlo[irow];
    f->rup = rup[irow];
    f->rowcols = copyOfArray(&hcol[krs], hinrow[irow]);
    f->rowels  = copyOfArray(&rowels[krs], hinrow[irow]);

    for (CoinBigIndex k=krs; k<kre; k++)
      presolve_delete_from_row(hcol[k], irow, mcstrt, hincol, hrow, colels);
    hinrow[irow] = 0;

    // just to make things squeeky
    rlo[irow] = 0.0;
    rup[irow] = 0.0;
  }


  next = new useless_constraint_action(nuseless_rows, actions, next);

  return (next);
}

const char *useless_constraint_action::name() const
{
  return ("useless_constraint_action");
}

void useless_constraint_action::postsolve(PostsolveMatrix *prob) const
{
  const action *const actions = actions_;
  const int nactions = nactions_;

  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  CoinBigIndex *mcstrt		= prob->mcstrt_;
  int *link		= prob->link_;
  int *hincol		= prob->hincol_;
  
  //  double *rowduals	= prob->rowduals_;
  double *rowacts	= prob->acts_;
  const double *sol	= prob->sol_;


  CoinBigIndex free_list		= prob->free_list_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  for (const action *f = &actions[nactions-1]; actions<=f; f--) {

    int irow	= f->row;
    int ninrow	= f->ninrow;
    const int *rowcols	= f->rowcols;
    const double *rowels = f->rowels;
    double rowact = 0.0;

    rup[irow] = f->rup;
    rlo[irow] = f->rlo;

    for (CoinBigIndex k=0; k<ninrow; k++) {
      int jcol = rowcols[k];
      //      CoinBigIndex kk = mcstrt[jcol];

      // append deleted row element to each col
      {
	CoinBigIndex kk = free_list;
	free_list = link[free_list];

	check_free_list(free_list);

	hrow[kk] = irow;
	colels[kk] = rowels[k];
	link[kk] = mcstrt[jcol];
	mcstrt[jcol] = kk;
      }
      
      rowact += rowels[k] * sol[jcol];
      hincol[jcol]++;
    }
    
    // I don't know if this is always true
    PRESOLVEASSERT(prob->getRowStatus(irow)==PrePostsolveMatrix::basic);
    PRESOLVEASSERT(rowduals[irow] == 0.0);    
    // rcosts are unaffected since rowdual is 0

    rowacts[irow] = rowact;
  }
  prob->free_list_ = free_list;
}


