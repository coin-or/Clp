#include <stdio.h>
#include <math.h>
#include <strings.h>

#include "PresolveMatrix.hpp"
#include "PresolveIsolated.hpp"


// Rarely, there may a constraint whose variables only
// occur in that constraint.
// In this case it is a completely independent problem.
// We should be able to solve it right now.
// Since that is actually not trivial, I'm just going to ignore
// them and stick them back in at postsolve.
const PresolveAction *isolated_constraint_action::presolve(PresolveMatrix *prob,
							    int irow,
							    const PresolveAction *next)
{
  int *hincol	= prob->hincol_;
  const int *mcstrt	= prob->mcstrt_;
  int *hrow	= prob->hrow_;
  double *colels	= prob->colels_;

  double *clo	= prob->clo_;
  double *cup	= prob->cup_;

  const double *rowels	= prob->rowels_;
  const int *hcol	= prob->hcol_;
  const int *mrstrt	= prob->mrstrt_;

  // may be written by useless constraint
  int *hinrow	= prob->hinrow_;

  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  int krs = mrstrt[irow];
  int kre = krs + hinrow[irow];
  
  double *dcost	= prob->cost_;

#if	DEBUG_PRESOLVE
  printf("ISOLATED:  %d - ", irow);
  for (int k = krs; k<kre; ++k)
    printf("%d ", hcol[k]);
  printf("\n");
#endif

  if (rlo[irow] != 0.0 || rup[irow] != 0.0)
    DIE("can't handle non-trivial isolated constraints for now\n");

  for (int k = krs; k<kre; ++k) {
    int jcol = hcol[k];
    if (clo[jcol] != 0.0 && cup[jcol] != 0.0)
      DIE("can't handle non-trivial isolated constraints for now\n");
  }

  int nc = hinrow[irow];

#if 0
  double tableau = new double[nc];
  double sol = new double[nc];
  double clo = new double[nc];
  double cup = new double[nc];


  for (int i=0; i<nc; ++i) {
    int col = hcol[krs+1];
    tableau[i] = rowels[krs+i];
    clo[i] = prob->clo[krs+i];
    cup[i] = prob->cup[krs+i];

    sol[i] = clo[i];
  }
#endif

  // HACK - set costs to 0.0 so empty.cpp doesn't complain
  double *costs = new double[nc];
  for (int k = krs; k<kre; ++k) {
    costs[k-krs] = dcost[hcol[k]];
    dcost[hcol[k]] = 0.0;
  }
  
  next = new isolated_constraint_action(rlo[irow], rup[irow],
					irow, nc,
					copyOfArray(&hcol[krs], nc*sizeof(int)),
					copyOfArray(&rowels[krs], nc*sizeof(double)),
					costs,
					next);

  for (int k=krs; k<kre; k++)
    presolve_delete_from_row(hcol[k], irow, mcstrt, hincol, hrow, colels);
  hinrow[irow] = 0;

  // just to make things squeeky
  rlo[irow] = 0.0;
  rup[irow] = 0.0;

  return (next);
}

const char *isolated_constraint_action::name() const
{
  return ("isolated_constraint_action");
}

void isolated_constraint_action::postsolve(PostsolveMatrix *prob) const
{
  double *colels	= prob->colels_;
  int *hrow		= prob->hrow_;
  int *mcstrt		= prob->mcstrt_;
  int *link		= prob->link_;
  int *hincol		= prob->hincol_;
  
  double *rowduals	= prob->rowduals_;
  double *rowacts	= prob->acts_;
  double *sol		= prob->sol_;

  int free_list		= prob->free_list_;


  // hides fields
  double *rlo	= prob->rlo_;
  double *rup	= prob->rup_;

  double rowact = 0.0;

  int irow  = this->row_;

  rup[irow] = this->rup_;
  rlo[irow] = this->rlo_;

  for (int k=0; k<this->ninrow_; k++) {
    int jcol = this->rowcols_[k];

    sol[jcol] = 0.0;	// ONLY ACCEPTED SUCH CONSTRAINTS

    int kk = free_list;
    free_list = link[free_list];

    check_free_list(free_list);

    mcstrt[jcol] = kk;

    //rowact += rowels[k] * sol[jcol];

    colels[kk] = this->rowels_[k];
    hrow[kk]   = irow;

    hincol[jcol] = 1;
  }

  // ???
  prob->setRowStatus(irow,PrePostsolveMatrix::basic);
    rowduals[irow] = 0.0;

  rowacts[irow] = rowact;

  prob->free_list_ = free_list;
}

