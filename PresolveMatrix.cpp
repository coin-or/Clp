// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#define	CHECK_CONSISTENCY	1

#include <stdio.h>

#include <assert.h>
#include <iostream>

#include "CoinHelperFunctions.hpp"

#include "CoinPackedMatrix.hpp"
#include "ClpSimplex.hpp"

#include "PresolveMatrix.hpp"
#include "Presolve.hpp"

void PresolveAction::throwCoinError(const char *error,
				       const char *ps_routine)
{
  throw CoinError(error, ps_routine, "Presolve");
}


  class presolvehlink;

void presolve_make_memlists(CoinBigIndex *starts, int *lengths,
		       presolvehlink *link,
		       int n);

/*static*/ void presolve_prefix(const int *starts, int *diffs, int n, int limit);

/*static*/ void presolve_prefix(const CoinBigIndex *starts, int *diffs, int n, int limit)
{
  int i;

  for (i=0; i<n; i++) {
    diffs[i] = starts[i+1] - starts[i];
  }
}

// afterward, link[i].pre is the largest index less than i st lengths[i]!=0
//  (or -1 if all such lengths are 0).
// link[i].suc is the opposite.
// That is, it is the same thing as setting link[i].pre = i-1 and link[i].suc = i+1
// and then deleting all the empty elements.
// This list is maintained together with hrow/hcol so that as we relocate
// columns or rows, we can quickly determine what column/row precedes a given
// column/row in the memory region hrow/hcol.
// Note that n itself has a pre and a suc.
void presolve_make_memlists(CoinBigIndex *starts, int *lengths, presolvehlink *link, int n)
{
  int i;
  int pre = NO_LINK;
  
  for (i=0; i<n; i++) {
    if (lengths[i]) {
      link[i].pre = pre;
      if (pre != NO_LINK)
	link[pre].suc = i;
      pre = i;
    }
    else {
      link[i].pre = NO_LINK;
      link[i].suc = NO_LINK;
    }
  }
  if (pre != NO_LINK)
    link[pre].suc = n;

  // (1) Arbitrarily place the last non-empty entry in link[n].pre
  link[n].pre = pre;

  link[n].suc = NO_LINK;
}



double *presolve_duparray(const double *d, int n, int n2)
{
  double *d1 = new double[n2];
  memcpy(d1, d, n*sizeof(double));
  return (d1);
}
double *presolve_duparray(const double *d, int n)
{
  return presolve_duparray(d, n, n);
}

int *presolve_duparray(const int *d, int n, int n2)
{
  int *d1 = new int[n2];
  memcpy(d1, d, n*sizeof(int));
  return (d1);
}
int *presolve_duparray(const int *d, int n)
{
  return presolve_duparray(d, n, n);
}

char *presolve_duparray(const char *d, int n, int n2)
{
  char *d1 = new char[n2];
  memcpy(d1, d, n*sizeof(char));
  return (d1);
}
char *presolve_duparray(const char *d, int n)
{
  return presolve_duparray(d, n, n);
}

#if 0
double *presolve_duparray(const double *d, int n, char **end_mmapp)
{
  double *d1 = (double*)*end_mmapp;
  memcpy(d, d1, n*sizeof(double));
  *end_mmapp += ALIGN_DOUBLE(n*sizeof(double));
  return (d1);
}
int *presolve_duparray(const int *d, int n, char **end_mmapp)
{
  int *d1 = (int*)*end_mmapp;
  memcpy(d, d1, n*sizeof(int));
  *end_mmapp += ALIGN_DOUBLE(n*sizeof(int));
  return (d1);
}

double *presolve_duparray(const double *d, int n, int n2, char **end_mmapp)
{
  double *d1 = (double*)*end_mmapp;
  memcpy(d, d1, n*sizeof(double));
  *end_mmapp += ALIGN_DOUBLE(n2*sizeof(double));
  return (d1);
}
int *presolve_duparray(const int *d, int n, int n2, char **end_mmapp)
{
  int *d1 = (int*)*end_mmapp;
  memcpy(d, d1, n*sizeof(int));
  *end_mmapp += ALIGN_DOUBLE(n2*sizeof(int));
  return (d1);
}
#endif


int presolve_find_row(int row, CoinBigIndex kcs, CoinBigIndex kce, const int *hrow)
{
  CoinBigIndex k;
  for (k=kcs; k<kce; k++)
    if (hrow[k] == row)
      return (k);
  DIE("FIND_ROW");
  return (0);
}

int presolve_find_row1(int row, CoinBigIndex kcs, CoinBigIndex kce, const int *hrow)
{
  CoinBigIndex k;
  for (k=kcs; k<kce; k++)
    if (hrow[k] == row)
      return (k);
  return (kce);
}

int presolve_find_row2(int irow, CoinBigIndex ks, int nc, const int *hrow, const int *link)
{
  for (int i=0; i<nc; ++i) {
    if (hrow[ks] == irow)
      return (ks);
    ks = link[ks];
  }
  abort();
  return -1;
}

int presolve_find_row3(int irow, CoinBigIndex ks, int nc, const int *hrow, const int *link)
{
  for (int i=0; i<nc; ++i) {
    if (hrow[ks] == irow)
      return (ks);
    ks = link[ks];
  }
  return (-1);
}


// delete the entry for col from row
// this keeps the row loosely packed
// if we didn't want to maintain that property, we could just decrement hinrow[row].
void presolve_delete_from_row(int row, int col /* thing to delete */,
		     const CoinBigIndex *mrstrt,
		     int *hinrow, int *hcol, double *dels)
{
  CoinBigIndex krs = mrstrt[row];
  CoinBigIndex kre = krs + hinrow[row];

  CoinBigIndex kcol = presolve_find_row(col, krs, kre, hcol);

  swap(hcol[kcol], hcol[kre-1]);
  swap(dels[kcol], dels[kre-1]);
  hinrow[row]--;
}


void presolve_delete_from_row2(int row, int col /* thing to delete */,
		      CoinBigIndex *mrstrt,
		     int *hinrow, int *hcol, double *dels, int *link, 
			       CoinBigIndex *free_listp)
{
  CoinBigIndex k = mrstrt[row];

  if (hcol[k] == col) {
    mrstrt[row] = link[k];
    link[k] = *free_listp;
    *free_listp = k;
    check_free_list(k);
    hinrow[row]--;
  } else {
    int n = hinrow[row] - 1;
    CoinBigIndex k0 = k;
    k = link[k];
    for (int i=0; i<n; ++i) {
      if (hcol[k] == col) {
	link[k0] = link[k];
	link[k] = *free_listp;
	*free_listp = k;
	check_free_list(k);
	hinrow[row]--;
	return;
      }
      k0 = k;
      k = link[k];
    }
    abort();
  }
}

static inline double getTolerance(const ClpSimplex &si, ClpDblParam key)
{
  double tol;
  if (! si.getDblParam(key, tol)) {
    PresolveAction::throwCoinError("getDblParam failed",
				      "PrePostsolveMatrix::PrePostsolveMatrix");
  }
  return (tol);
}

#if 0
PresolveAction::~PresolveAction()
{
  cout << "OOPS\n";
}
#endif

// Assumptions:
// 1. nrows>=m.getNumRows()
// 2. ncols>=m.getNumCols()
//
// In presolve, these values are equal.
// In postsolve, they may be inequal, since the reduced problem
// may be smaller, but we need room for the large problem.
// ncols may be larger than si.getNumCols() in postsolve,
// this at that point si will be the reduced problem,
// but we need to reserve enough space for the original problem.
PrePostsolveMatrix::PrePostsolveMatrix(const ClpSimplex& si,
					     int ncols_in,
					     int nrows_in,
					     CoinBigIndex nelems_in) :
  ncols_(si.getNumCols()),
  ncols0_(ncols_in),
  nelems_(si.getNumElements()),

  mcstrt_(new CoinBigIndex[ncols_in+1]),
  hincol_(new int[ncols_in+1]),
  hrow_  (new int   [2*nelems_in]),
  colels_(new double[2*nelems_in]),

  cost_(new double[ncols_in]),
  clo_(new double[ncols_in]),
  cup_(new double[ncols_in]),
  rlo_(new double[nrows_in]),
  rup_(new double[nrows_in]),
  originalColumn_(new int[ncols_in]),

  ztolzb_(getTolerance(si, ClpPrimalTolerance)),
  ztoldj_(getTolerance(si, ClpDualTolerance)),

  maxmin_(si.getObjSense())

{
  si.getDblParam(ClpObjOffset,originalOffset_);
  int ncols = si.getNumCols();
  int nrows = si.getNumRows();

  ClpDisjointCopyN(si.getColLower(), ncols, clo_);
  ClpDisjointCopyN(si.getColUpper(), ncols, cup_);
  ClpDisjointCopyN(si.getObjCoefficients(), ncols, cost_);
  ClpDisjointCopyN(si.getRowLower(), nrows,  rlo_);
  ClpDisjointCopyN(si.getRowUpper(), nrows,  rup_);
  int i;
  for (i=0;i<ncols_in;i++) 
    originalColumn_[i]=i;
  sol_=NULL;
  rowduals_=NULL;
  acts_=NULL;

  rcosts_=NULL;
  colstat_=NULL;
  rowstat_=NULL;
}


PrePostsolveMatrix::~PrePostsolveMatrix()
{
  delete[]mcstrt_;
  delete[]hrow_;
  delete[]colels_;
  delete[]hincol_;

  delete[]cost_;
  delete[]clo_;
  delete[]cup_;
  delete[]rlo_;
  delete[]rup_;
  delete[]originalColumn_;
  delete[]rowduals_;

  delete[]rcosts_;
}

// Sets status (non -basic ) using value
void 
PrePostsolveMatrix::setRowStatusUsingValue(int iRow)
{
  double value = acts_[iRow];
  double lower = rlo_[iRow];
  double upper = rup_[iRow];
  if (lower<-1.0e20&&upper>1.0e20) {
    setRowStatus(iRow,isFree);
  } else if (fabs(lower-value)<=ztolzb_) {
    setRowStatus(iRow,atLowerBound);
  } else if (fabs(upper-value)<=ztolzb_) {
    setRowStatus(iRow,atUpperBound);
  } else {
    setRowStatus(iRow,superBasic);
  }
}
void 
PrePostsolveMatrix::setColumnStatusUsingValue(int iColumn)
{
  double value = sol_[iColumn];
  double lower = clo_[iColumn];
  double upper = cup_[iColumn];
  if (lower<-1.0e20&&upper>1.0e20) {
    setColumnStatus(iColumn,isFree);
  } else if (fabs(lower-value)<=ztolzb_) {
    setColumnStatus(iColumn,atLowerBound);
  } else if (fabs(upper-value)<=ztolzb_) {
    setColumnStatus(iColumn,atUpperBound);
  } else {
    setColumnStatus(iColumn,superBasic);
  }
}

// I am not familiar enough with CoinPackedMatrix to be confident
// that I will implement a row-ordered version of toColumnOrderedGapFree
// properly.
static bool isGapFree(const CoinPackedMatrix& matrix)
{
  const CoinBigIndex * start = matrix.getVectorStarts();
  const int * length = matrix.getVectorLengths();
  int i;
  for (i = matrix.getSizeVectorLengths() - 1; i >= 0; --i) {
    if (start[i+1] - start[i] != length[i])
      break;
  }
  return (! (i >= 0));
}


#if	DEBUG_PRESOLVE
static void matrix_bounds_ok(const double *lo, const double *up, int n)
{
  int i;
  for (i=0; i<n; i++) {
    PRESOLVEASSERT(lo[i] <= up[i]);
    PRESOLVEASSERT(lo[i] < PRESOLVE_INF);
    PRESOLVEASSERT(-PRESOLVE_INF < up[i]);
  }
}
#endif

PresolveMatrix::PresolveMatrix(int ncols0_in,
				     double maxmin_,
				     // end prepost members

				     ClpSimplex &si,

				     // rowrep
				     int nrows_in,
				     CoinBigIndex nelems_in,
			       bool doStatus) :

  PrePostsolveMatrix(si,
			ncols0_in, nrows_in, nelems_in),
  clink_(new presolvehlink[ncols0_in+1]),
  rlink_(new presolvehlink[nrows_in+1]),

  dobias_(0.0),

  nrows_(si.getNumRows()),

  // temporary init
  mrstrt_(new CoinBigIndex[nrows_in+1]),
  hinrow_(new int[nrows_in+1]),
  rowels_(new double[2*nelems_in]),
  hcol_(new int[2*nelems_in]),
  integerType_(new char[ncols0_in]),
  feasibilityTolerance_(0.0),
  status_(-1),
  rowsToDo_(new int [nrows_in]),
  numberRowsToDo_(0),
  nextRowsToDo_(new int[nrows_in]),
  numberNextRowsToDo_(0),
  colsToDo_(new int [ncols0_in]),
  numberColsToDo_(0),
  nextColsToDo_(new int[ncols0_in]),
  numberNextColsToDo_(0)

{
  const int bufsize = 2*nelems_in;

  // Set up change bits
  rowChanged_ = new unsigned int[(nrows_+31)>>5];
  memset(rowChanged_,0,((nrows_+31)>>5)*sizeof(unsigned int));
  colChanged_ = new unsigned int[(ncols_+31)>>5];
  memset(colChanged_,0,((ncols_+31)>>5)*sizeof(unsigned int));
  CoinPackedMatrix * m = si.matrix();

  // The coefficient matrix is a big hunk of stuff.
  // Do the copy here to try to avoid running out of memory.

  const CoinBigIndex * start = m->getVectorStarts();
  const int * length = m->getVectorLengths();
  const int * row = m->getIndices();
  const double * element = m->getElements();
  int icol,nel=0;
  mcstrt_[0]=0;
  for (icol=0;icol<ncols_;icol++) {
    int j;
    for (j=start[icol];j<start[icol]+length[icol];j++) {
      hrow_[nel]=row[j];
      colels_[nel++]=element[j];
    }
    mcstrt_[icol+1]=nel;
  }
  assert(mcstrt_[ncols_] == nelems_);
  ClpDisjointCopyN(m->getVectorLengths(),ncols_,  hincol_);

  // same thing for row rep
  m = new CoinPackedMatrix();
  m->reverseOrderedCopyOf(*si.matrix());
  m->removeGaps();


  ClpDisjointCopyN(m->getVectorStarts(),  nrows_,  mrstrt_);
  mrstrt_[nrows_] = nelems_;
  ClpDisjointCopyN(m->getVectorLengths(), nrows_,  hinrow_);
  ClpDisjointCopyN(m->getIndices(),       nelems_, hcol_);
  ClpDisjointCopyN(m->getElements(),      nelems_, rowels_);

  delete m;
  if (si.integerInformation()) {
    memcpy(integerType_,si.integerInformation(),ncols_*sizeof(char));
  } else {
    ClpFillN<char>(integerType_, ncols_, 0);
  }

  if (doStatus) {
    // allow for status and solution
    sol_ = new double[ncols_];
    memcpy(sol_,si.primalColumnSolution(),ncols_*sizeof(double));;
    acts_ = new double [nrows_];
    memcpy(acts_,si.primalRowSolution(),nrows_*sizeof(double));
    if (!si.statusArray())
      si.createStatus();
    colstat_ = new unsigned char [nrows_+ncols_];
    memcpy(colstat_,si.statusArray(),
	   (nrows_+ncols_)*sizeof(unsigned char));
    rowstat_ = colstat_+ncols_;
  }

  // the original model's fields are now unneeded - free them
  
  si.resize(0,0);

#if	DEBUG_PRESOLVE
  matrix_bounds_ok(rlo_, rup_, nrows_);
  matrix_bounds_ok(clo_, cup_, ncols_);
#endif

#if 0
  for (i=0; i<nrows; ++i)
    printf("NR: %6d\n", hinrow[i]);
  for (int i=0; i<ncols; ++i)
    printf("NC: %6d\n", hincol[i]);
#endif

  presolve_make_memlists(mcstrt_, hincol_, clink_, ncols_);
  presolve_make_memlists(mrstrt_, hinrow_, rlink_, nrows_);

  // this allows last col/row to expand up to bufsize-1 (22);
  // this must come after the calls to presolve_prefix
  mcstrt_[ncols_] = bufsize-1;
  mrstrt_[nrows_] = bufsize-1;

#if	CHECK_CONSISTENCY
  consistent(false);
#endif
}

PresolveMatrix::~PresolveMatrix()
{
  delete[]clink_;
  delete[]rlink_;
  
  delete[]mrstrt_;
  delete[]hinrow_;
  delete[]rowels_;
  delete[]hcol_;

  delete[]integerType_;
  delete[] rowChanged_;
  delete[] rowsToDo_;
  delete[] nextRowsToDo_;
  delete[] colChanged_;
  delete[] colsToDo_;
  delete[] nextColsToDo_;
}


void PresolveMatrix::update_model(ClpSimplex& si,
				     int nrows0,
				     int ncols0,
				     CoinBigIndex nelems0)
{
  si.loadProblem(ncols_, nrows_, mcstrt_, hrow_, colels_, hincol_,
		 clo_, cup_, cost_, rlo_, rup_);

  delete [] si.integerInformation();
  int numberIntegers=0;
  for (int i=0; i<ncols_; i++) {
    if (integerType_[i])
      numberIntegers++;
  }
  if (numberIntegers) 
    si.copyInIntegerInformation(integerType_);
  else
    si.copyInIntegerInformation(NULL);

#if	PRESOLVE_SUMMARY
  printf("NEW NCOL/NROW/NELS:  %d(-%d) %d(-%d) %d(-%d)\n",
	 ncols_, ncols0-ncols_,
	 nrows_, nrows0-nrows_,
	 si.getNumElements(), nelems0-si.getNumElements());
#endif
  si.setDblParam(ClpObjOffset,originalOffset_-dobias_);

}


// The matrix is represented redundantly in both row and column format,
// in what I call "loosely packed" format.
// "packed" format would be as in normal OSL:  a vector of column starts,
// together with two vectors that give the row indices and coefficients.
//
// This format is "loosely packed" because some of the elements may correspond
// to empty rows.  This is so that we can quickly delete rows without having
// to update the column rep and vice versa.
//
// Checks whether an entry is in the col rep iff it is also in the row rep,
// and also that their values are the same (if testvals is non-zero).
//
// Note that there may be entries in a row that correspond to empty columns
// and vice-versa.  --- HUH???
static void matrix_consistent(const CoinBigIndex *mrstrt, const int *hinrow, const int *hcol,
			      const CoinBigIndex *mcstrt, const int *hincol, const int *hrow,
			      const double *rowels,
			      const double *colels,
			      int nrows, int testvals,
			      const char *ROW, const char *COL)
{
#if	CHECK_CONSISTENCY
  for (int irow=0; irow<nrows; irow++) {
    if (hinrow[irow] > 0) {
      CoinBigIndex krs = mrstrt[irow];
      CoinBigIndex kre = krs + hinrow[irow];

      for (CoinBigIndex k=krs; k<kre; k++) {
	int jcol = hcol[k];
	CoinBigIndex kcs = mcstrt[jcol];
	CoinBigIndex kce = kcs + hincol[jcol];

	CoinBigIndex kk = presolve_find_row1(irow, kcs, kce, hrow);
	if (kk == kce) {
	  printf("MATRIX INCONSISTENT:  can't find %s %d in %s %d\n",
		 ROW, irow, COL, jcol);
	  fflush(stdout);
	  abort();
	}
	if (testvals && colels[kk] != rowels[k]) {
	  printf("MATRIX INCONSISTENT:  value for %s %d and %s %d\n",
		 ROW, irow, COL, jcol);
	  fflush(stdout);
	  abort();
	}
      }
    }
  }
#endif
}


void PresolveMatrix::consistent(bool testvals)
{
#if	CHECK_CONSISTENCY
  matrix_consistent(mrstrt_, hinrow_, hcol_,
		    mcstrt_, hincol_, hrow_,
		    rowels_, colels_,
		    nrows_, testvals,
		    "row", "col");
  matrix_consistent(mcstrt_, hincol_, hrow_,
		    mrstrt_, hinrow_, hcol_,
		    colels_, rowels_, 
		    ncols_, testvals,
		    "col", "row");
#endif
}











////////////////  POSTSOLVE

PostsolveMatrix::PostsolveMatrix(ClpSimplex& si,
				       int ncols0_in,
				       int nrows0_in,
				       CoinBigIndex nelems0,
				   
				       double maxmin_,
				       // end prepost members

				       double *sol_in,
				       double *acts_in,

				       unsigned char *colstat_in,
				       unsigned char *rowstat_in) :
  PrePostsolveMatrix(si,
			ncols0_in, nrows0_in, nelems0),

  free_list_(0),
  // link, free_list, maxlink
  maxlink_(2*nelems0),
  link_(new int[/*maxlink*/ 2*nelems0]),
      
  cdone_(new char[ncols0_]),
  rdone_(new char[nrows0_in]),

  nrows_(si.getNumRows()),
  nrows0_(nrows0_in)
{

  sol_=sol_in;
  rowduals_=NULL;
  acts_=acts_in;

  rcosts_=NULL;
  colstat_=colstat_in;
  rowstat_=rowstat_in;

  // this is the *reduced* model, which is probably smaller
  int ncols1 = si.getNumCols();
  int nrows1 = si.getNumRows();

  const CoinPackedMatrix * m = si.matrix();

  if (! isGapFree(*m)) {
    PresolveAction::throwCoinError("Matrix not gap free",
				      "PostsolveMatrix");
  }

  const CoinBigIndex nelemsr = m->getNumElements();

  ClpDisjointCopyN(m->getVectorStarts(), ncols1, mcstrt_);
  mcstrt_[ncols_] = nelems0;	// ??
  ClpDisjointCopyN(m->getVectorLengths(),ncols1,  hincol_);
  ClpDisjointCopyN(m->getIndices(),      nelemsr, hrow_);
  ClpDisjointCopyN(m->getElements(),     nelemsr, colels_);


#if	0 && DEBUG_PRESOLVE
  presolve_check_costs(model, &colcopy);
#endif

  // This determines the size of the data structure that contains
  // the matrix being postsolved.  Links are taken from the free_list
  // to recreate matrix entries that were presolved away,
  // and links are added to the free_list when entries created during
  // presolve are discarded.  There is never a need to gc this list.
  // Naturally, it should contain
  // exactly nelems0 entries "in use" when postsolving is done,
  // but I don't know whether the matrix could temporarily get
  // larger during postsolving.  Substitution into more than two
  // rows could do that, in principle.  I am being very conservative
  // here by reserving much more than the amount of space I probably need.
  // If this guess is wrong, check_free_list may be called.
  //  int bufsize = 2*nelems0;

  memset(cdone_, -1, ncols0_);
  memset(rdone_, -1, nrows0_);

  rowduals_ = new double[nrows0_];
  ClpDisjointCopyN(si.getRowPrice(), nrows1, rowduals_);

  rcosts_ = new double[ncols0_];
  ClpDisjointCopyN(si.getReducedCost(), ncols1, rcosts_);
  if (maxmin_<0.0) {
    // change so will look as if minimize
    int i;
    for (i=0;i<nrows1;i++)
      rowduals_[i] = - rowduals_[i];
    for (i=0;i<ncols1;i++) {
      rcosts_[i] = - rcosts_[i];
    }
  }

  //ClpDisjointCopyN(si.getRowUpper(), nrows1, rup_);
  //ClpDisjointCopyN(si.getRowLower(), nrows1, rlo_);

  ClpDisjointCopyN(si.getColSolution(), ncols1, sol_);
  si.setDblParam(ClpObjOffset,originalOffset_);

  for (int j=0; j<ncols1; j++) {
    CoinBigIndex kcs = mcstrt_[j];
    CoinBigIndex kce = kcs + hincol_[j];
    for (CoinBigIndex k=kcs; k<kce; ++k) {
      link_[k] = k+1;
    }
  }
  {
    int ml = maxlink_;
    for (CoinBigIndex k=nelemsr; k<ml; ++k)
      link_[k] = k+1;
    link_[ml-1] = NO_LINK;
  }
  free_list_ = nelemsr;
}

PostsolveMatrix::~PostsolveMatrix()
{
  delete[]link_;

  delete[]cdone_;
  delete[]rdone_;
}


void PostsolveMatrix::check_nbasic()
{
  int nbasic = 0;

  int i;
  for (i=0; i<ncols_; i++)
    if (columnIsBasic(i))
      nbasic++;

  for (i=0; i<nrows_; i++)
    if (rowIsBasic(i))
      nbasic++;

  if (nbasic != nrows_) {
    printf("WRONG NUMBER NBASIC:  is:  %d  should be:  %d\n",
	   nbasic, nrows_);
    fflush(stdout);
  }
}






