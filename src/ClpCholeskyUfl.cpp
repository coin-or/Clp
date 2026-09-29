// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpConfig.h"

#ifndef CLP_HAS_CHOLMOD
#ifndef CLP_HAS_AMD
#error "Need to have AMD or CHOLMOD to compile ClpCholeskyUfl."
#else
// Do NOT wrap this include in our own extern "C" block: amd.h (and the
// SuiteSparse_config.h it pulls in) already guards its own declarations with
// `#ifdef __cplusplus extern "C" { ... }`. Modern SuiteSparse versions'
// SuiteSparse_config.h also does `#include <complex>` under __cplusplus
// (for BLAS/LAPACK complex types) — if that happens while already inside an
// extra outer extern "C" opened here, the compiler correctly rejects the
// C++ templates pulled in by <complex> as having invalid C linkage.
#include "amd.h"
#endif
#else
#include "cholmod.h"
#endif

#include "CoinPragma.hpp"
#include "ClpCholeskyUfl.hpp"
#include "ClpMessage.hpp"
#include "ClpInterior.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpCholeskyUfl::ClpCholeskyUfl(int denseThreshold)
  : ClpCholeskyBase(denseThreshold)
{
  type_ = 14;
  L_ = NULL;
  c_ = NULL;
  B_ = NULL;
  X_ = NULL;
  Y_ = NULL;
  E_ = NULL;
  forceSimplicial_ = false;

#ifdef CLP_HAS_CHOLMOD
  c_ = (cholmod_common *)malloc(sizeof(cholmod_common));
  cholmod_start(c_);
  /* Let CHOLMOD choose between the simplicial and the supernodal method.
     A*D*A' + delta^2*I is positive definite, and the supernodal method is
     the only one that uses the (threaded) dense BLAS on the dense trailing
     blocks of the factor - which is exactly what the native ClpCholeskyBase
     does by hand ("going dense for last n rows").  Forcing simplicial here
     gave up that whole gain on every instance.  If a factorization ever does
     fail numerically we fall back to simplicial LDL' for good - see
     factorize(). */
  c_->supernodal = CHOLMOD_AUTO;
#endif
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpCholeskyUfl::ClpCholeskyUfl(const ClpCholeskyUfl &rhs)
  : ClpCholeskyBase(rhs)
{
  abort();
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpCholeskyUfl::~ClpCholeskyUfl()
{
#ifdef CLP_HAS_CHOLMOD
  cholmod_free_dense(&B_, c_);
  cholmod_free_dense(&X_, c_);
  cholmod_free_dense(&Y_, c_);
  cholmod_free_dense(&E_, c_);
  cholmod_free_factor(&L_, c_);
  cholmod_finish(c_);
  free(c_);
#endif
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpCholeskyUfl &
ClpCholeskyUfl::operator=(const ClpCholeskyUfl &rhs)
{
  if (this != &rhs) {
    ClpCholeskyBase::operator=(rhs);
    abort();
  }
  return *this;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase *ClpCholeskyUfl::clone() const
{
  return new ClpCholeskyUfl(*this);
}

#ifndef CLP_HAS_CHOLMOD
/* Orders rows and saves pointer to matrix.and model */
int ClpCholeskyUfl::order(ClpInterior *model)
{
  int iRow;
  model_ = model;
  if (preOrder(false, true, doKKT_))
    return -1;
  permuteInverse_ = new int[numberRows_];
  permute_ = new int[numberRows_];
  double Control[AMD_CONTROL];
  double Info[AMD_INFO];

  amd_defaults(Control);
  //amd_control(Control);

  int returnCode = amd_order(numberRows_, choleskyStart_, choleskyRow_,
    permute_, Control, Info);
  delete[] choleskyRow_;
  choleskyRow_ = NULL;
  delete[] choleskyStart_;
  choleskyStart_ = NULL;
  //amd_info(Info);

  if (returnCode != AMD_OK) {
    std::cout << "AMD ordering failed" << std::endl;
    return 1;
  }
  for (iRow = 0; iRow < numberRows_; iRow++) {
    permuteInverse_[permute_[iRow]] = iRow;
  }
  return 0;
}
#else
/* Orders rows and saves pointer to matrix.and model */
int ClpCholeskyUfl::order(ClpInterior *model)
{
  /* The CHOLMOD objects below are told the integer arrays are CHOLMOD_INT,
     so a 64 bit CoinBigIndex would be silently misread. */
  static_assert(sizeof(CoinBigIndex) == sizeof(int),
    "ClpCholeskyUfl needs a 32 bit CoinBigIndex for CHOLMOD_INT");
  numberRows_ = model->numberRows();
  if (doKKT_) {
    numberRows_ += numberRows_ + model->numberColumns();
    printf("finish coding UFL KKT!\n");
    abort();
  }
  /* This path never calls ClpCholeskyBase::preOrder, so the dense column
     machinery (whichDense_/denseColumn_/dense_) is never set up, and solve()
     below has no low rank correction for it either.  Make sure nothing else
     thinks dense columns are being handled. */
  denseThreshold_ = -1;
  delete[] rowsDropped_;
  rowsDropped_ = new char[numberRows_];
  memset(rowsDropped_, 0, numberRows_);
  numberRowsDropped_ = 0;
  model_ = model;
#if defined(CHOLMOD_VERSION) && CHOLMOD_VERSION >= CHOLMOD_VER_CODE(5, 0)
  /* Respect the model's thread budget instead of letting CHOLMOD grab every
     core it can see.  CHOLMOD defaults nthreads_max to the OpenMP maximum,
     so a barrier solve run inside a 32-way parallel harness - or inside one
     of Cbc's -threads N branch and bound threads, each with its own Clp -
     silently oversubscribes the machine by that factor.  What is actually
     paid for is mostly spinning: on a supernodal instance like seymour this
     was ~49% of all cycles in gomp_barrier_wait_end, 2.3x the CPU time for
     no wall clock gain at all. */
  int numberThreads = model_->numberThreads();
  c_->nthreads_max = (numberThreads > 0) ? numberThreads : 1;
#endif
  delete rowCopy_;
  rowCopy_ = model->clpMatrix()->reverseOrderedCopy();
  // Space for starts
  choleskyStart_ = new CoinBigIndex[numberRows_ + 1];
  const CoinBigIndex *columnStart = model_->clpMatrix()->getVectorStarts();
  const int *columnLength = model_->clpMatrix()->getVectorLengths();
  const int *row = model_->clpMatrix()->getIndices();
  const CoinBigIndex *rowStart = rowCopy_->getVectorStarts();
  const int *rowLength = rowCopy_->getVectorLengths();
  const int *column = rowCopy_->getIndices();
  // We need two arrays for counts
  CoinBigIndex *which = new CoinBigIndex[numberRows_];
  int *used = new int[numberRows_ + 1];
  CoinZeroN(used, numberRows_);
  int iRow;
  sizeFactor_ = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    int number = 1;
    // make sure diagonal exists
    which[0] = iRow;
    used[iRow] = 1;
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow = rowStart[iRow];
      CoinBigIndex endRow = rowStart[iRow] + rowLength[iRow];
      for (CoinBigIndex k = startRow; k < endRow; k++) {
        int iColumn = column[k];
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn] + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int jRow = row[j];
          if (jRow >= iRow && !rowsDropped_[jRow]) {
            if (!used[jRow]) {
              used[jRow] = 1;
              which[number++] = jRow;
            }
          }
        }
      }
      sizeFactor_ += number;
      int j;
      for (j = 0; j < number; j++)
        used[which[j]] = 0;
    }
  }
  delete[] which;
  // Now we have size - create arrays and fill in
  try {
    choleskyRow_ = new CoinBigIndex[sizeFactor_];
  } catch (...) {
    // no memory
    delete[] choleskyStart_;
    choleskyStart_ = NULL;
    return -1;
  }
  try {
    sparseFactor_ = new double[sizeFactor_];
  } catch (...) {
    // no memory
    delete[] choleskyRow_;
    choleskyRow_ = NULL;
    delete[] choleskyStart_;
    choleskyStart_ = NULL;
    return -1;
  }

  sizeFactor_ = 0;
  which = choleskyRow_;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    int number = 1;
    // make sure diagonal exists
    which[0] = iRow;
    used[iRow] = 1;
    choleskyStart_[iRow] = sizeFactor_;
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow = rowStart[iRow];
      CoinBigIndex endRow = rowStart[iRow] + rowLength[iRow];
      for (CoinBigIndex k = startRow; k < endRow; k++) {
        int iColumn = column[k];
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn] + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int jRow = row[j];
          if (jRow >= iRow && !rowsDropped_[jRow]) {
            if (!used[jRow]) {
              used[jRow] = 1;
              which[number++] = jRow;
            }
          }
        }
      }
      sizeFactor_ += number;
      int j;
      for (j = 0; j < number; j++)
        used[which[j]] = 0;
      // Sort
      std::sort(which, which + number);
      // move which on
      which += number;
    }
  }
  choleskyStart_[numberRows_] = sizeFactor_;
  delete[] used;
  permuteInverse_ = new CoinBigIndex[numberRows_];
  permute_ = new CoinBigIndex[numberRows_];
  cholmod_sparse A;
  A.nrow = numberRows_;
  A.ncol = numberRows_;
  A.nzmax = choleskyStart_[numberRows_];
  A.p = choleskyStart_;
  A.i = choleskyRow_;
  A.x = NULL;
  A.stype = -1;
  A.itype = CHOLMOD_INT;
  A.xtype = CHOLMOD_PATTERN;
  A.dtype = CHOLMOD_DOUBLE;
  A.sorted = 1;
  A.packed = 1;
  /* Use CHOLMOD's own ordering strategy: AMD, and METIS on top of it only
     when AMD's fill-in turns out to be poor.  The previous value (9) asked
     CHOLMOD to run *every* built-in method - AMD, COLAMD, METIS and four
     NESDIS nested dissection variants - and keep the best.  On the normal
     equations of a large LP the four NESDIS runs alone dominate the whole
     barrier solve, for an ordering that is rarely better than AMD/METIS. */
  c_->nmethods = 0;
  c_->postorder = true;
  //c_->dbound=1.0e-20;
  cholmod_free_factor(&L_, c_);
  cholmod_free_dense(&B_, c_);
  cholmod_free_dense(&X_, c_);
  cholmod_free_dense(&Y_, c_);
  cholmod_free_dense(&E_, c_);
  L_ = cholmod_analyze(&A, c_);
  if (c_->status || !L_) {
    COIN_DETAIL_PRINT(std::cout << "CHOLMOD ordering failed" << std::endl);
    return 1;
  } else if (c_->lnz > static_cast< double >(COIN_INT_MAX)) {
    // With 32-bit (CHOLMOD_INT) indices cholmod_factorize would fail with
    // CHOLMOD_TOO_LARGE and leave L_ symbolic - give up on barrier now.
    printf("CHOLMOD: factor too large for 32-bit indices (%g nonzeros predicted)\n",
      c_->lnz);
    return 1;
  } else {
    COIN_DETAIL_PRINT(printf("%g nonzeros, flop count %g\n", c_->lnz, c_->fl));
  }
  for (iRow = 0; iRow < numberRows_; iRow++) {
    permuteInverse_[iRow] = iRow;
    permute_[iRow] = iRow;
  }
  return 0;
}
#endif

/* Does Symbolic factorization given permutation.
   This is called immediately after order.  If user provides this then
   user must provide factorize and solve.  Otherwise the default factorization is used
   returns non-zero if not enough memory */
int ClpCholeskyUfl::symbolic()
{
#ifdef CLP_HAS_CHOLMOD
  return 0;
#else
  return ClpCholeskyBase::symbolic();
#endif
}

#ifdef CLP_HAS_CHOLMOD
/* Maximum number of rows dropped-and-retried before giving up on the
   supernodal method for the rest of the solve. */
#define CLP_UFL_MAX_DROP_RETRIES 2

/* Turn every dropped row into a unit row of the matrix being factorized.
   Clamping only its diagonal to 1e-10 while leaving its off-diagonals in
   place (what the old code did) hands CHOLMOD an almost singular matrix and
   the resulting direction is meaningless. */
void ClpCholeskyUfl::applyDroppedRows(const int *rowsDropped)
{
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    CoinBigIndex start = choleskyStart_[iRow];
    CoinBigIndex end = choleskyStart_[iRow + 1];
    if (rowsDropped[iRow]) {
      sparseFactor_[start] = 1.0;
      for (CoinBigIndex j = start + 1; j < end; j++)
        sparseFactor_[j] = 0.0;
    } else {
      for (CoinBigIndex j = start + 1; j < end; j++) {
        if (rowsDropped[choleskyRow_[j]])
          sparseFactor_[j] = 0.0;
      }
    }
  }
}

/* Factorize - filling in rowsDropped and returning number dropped */
int ClpCholeskyUfl::factorize(const double *diagonal, int *rowsDropped)
{
  const CoinBigIndex *columnStart = model_->clpMatrix()->getVectorStarts();
  const int *columnLength = model_->clpMatrix()->getVectorLengths();
  const int *row = model_->clpMatrix()->getIndices();
  const double *element = model_->clpMatrix()->getElements();
  const CoinBigIndex *rowStart = rowCopy_->getVectorStarts();
  const int *rowLength = rowCopy_->getVectorLengths();
  const int *column = rowCopy_->getIndices();
  const double *elementByRow = rowCopy_->getElements();
  int numberColumns = model_->clpMatrix()->getNumCols();
  int iRow;
  double *work = new double[numberRows_];
  CoinZeroN(work, numberRows_);
  const double *diagonalSlack = diagonal + numberColumns;
  int newDropped = 0;
  double largest;
  //double smallest;
  //perturbation
  double perturbation = model_->diagonalPerturbation() * model_->diagonalNorm();
  perturbation = 0.0;
  perturbation = perturbation * perturbation;
  if (perturbation > 1.0) {
#ifdef COIN_DEVELOP
    //if (model_->model()->logLevel()&4)
    std::cout << "large perturbation " << perturbation << std::endl;
#endif
    perturbation = sqrt(perturbation);
    ;
    perturbation = 1.0;
  }
  double delta2 = model_->delta(); // add delta*delta to diagonal
  delta2 *= delta2;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    double *put = sparseFactor_ + choleskyStart_[iRow];
    CoinBigIndex *which = choleskyRow_ + choleskyStart_[iRow];
    int number = choleskyStart_[iRow + 1] - choleskyStart_[iRow];
    if (!rowLength[iRow])
      rowsDropped_[iRow] = 1;
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow = rowStart[iRow];
      CoinBigIndex endRow = rowStart[iRow] + rowLength[iRow];
      work[iRow] = diagonalSlack[iRow] + delta2;
      for (CoinBigIndex k = startRow; k < endRow; k++) {
        int iColumn = column[k];
        if (!whichDense_ || !whichDense_[iColumn]) {
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = columnStart[iColumn] + columnLength[iColumn];
          double multiplier = diagonal[iColumn] * elementByRow[k];
          for (CoinBigIndex j = start; j < end; j++) {
            int jRow = row[j];
            if (jRow >= iRow && !rowsDropped_[jRow]) {
              double value = element[j] * multiplier;
              work[jRow] += value;
            }
          }
        }
      }
      int j;
      for (j = 0; j < number; j++) {
        int jRow = which[j];
        put[j] = work[jRow];
        work[jRow] = 0.0;
      }
    } else {
      // dropped
      int j;
      for (j = 1; j < number; j++) {
        put[j] = 0.0;
      }
      put[0] = 1.0;
    }
  }
  //check sizes
  double largest2 = maximumAbsElement(sparseFactor_, sizeFactor_);
  largest2 *= 1.0e-20;
  largest = std::min(largest2, 1.0e-11);
  int numberDroppedBefore = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    int dropped = rowsDropped_[iRow];
    // Move to int array
    rowsDropped[iRow] = dropped;
    if (!dropped) {
      CoinBigIndex start = choleskyStart_[iRow];
      double diagonal = sparseFactor_[start];
      if (diagonal > largest2) {
        sparseFactor_[start] = std::max(diagonal, 1.0e-10);
      } else {
        sparseFactor_[start] = std::max(diagonal, 1.0e-10);
        rowsDropped[iRow] = 2;
        numberDroppedBefore++;
      }
    }
  }
  delete[] work;
  /* rows dropped for a tiny diagonal have to be reported to the caller, and
     remembered, exactly as ClpCholeskyBase::factorize does - the old code
     computed numberDroppedBefore and then threw it away, so rowsDropped_ was
     never updated and the caller was always told that nothing was dropped. */
  newDropped = numberDroppedBefore;
  if (newDropped || numberRowsDropped_)
    applyDroppedRows(rowsDropped);
  cholmod_sparse A;
  A.nrow = numberRows_;
  A.ncol = numberRows_;
  A.nzmax = choleskyStart_[numberRows_];
  A.p = choleskyStart_;
  A.i = choleskyRow_;
  A.x = sparseFactor_;
  A.stype = -1;
  A.itype = CHOLMOD_INT;
  A.xtype = CHOLMOD_REAL;
  A.dtype = CHOLMOD_DOUBLE;
  A.sorted = 1;
  A.packed = 1;
  cholmod_factorize(&A, L_, c_); /* factorize */
  /* The supernodal method does no pivoting at all: it gives up at the first
     non positive pivot, reporting the offending column in L_->minor.  Near
     convergence the barrier diagonal becomes extreme and this happens on
     most instances, so simply switching to the simplicial method for the
     rest of the solve - what used to happen here - throws away BLAS3 for
     every remaining iteration.  On an instance with a dense-ish factor that
     is catastrophic (sorrell3 spends 99% of its time in simplicial rowfac).
     Instead drop the offending row, exactly as ClpCholeskyBase drops a row
     with a non positive pivot, and refactorize.  The sparsity pattern is
     unchanged, so the symbolic analysis - and with it the supernodal
     blocking - is reused; only the numeric factorization is repeated. */
  int numberRetries = 0;
  while (!forceSimplicial_ && L_ && c_->status >= CHOLMOD_OK
    && L_->xtype != CHOLMOD_PATTERN
    && L_->minor < static_cast< size_t >(numberRows_)
    && numberRetries < CLP_UFL_MAX_DROP_RETRIES) {
    /* L_->minor is in the permuted order used by the factor. */
    int badRow = static_cast< int >(L_->minor);
    if (L_->Perm)
      badRow = static_cast< int * >(L_->Perm)[L_->minor];
    if (badRow < 0 || badRow >= numberRows_ || rowsDropped[badRow])
      break;
    rowsDropped[badRow] = 2;
    newDropped++;
    applyDroppedRows(rowsDropped);
    numberRetries++;
    cholmod_factorize(&A, L_, c_);
  }
  if (!forceSimplicial_ && L_ && c_->status >= CHOLMOD_OK
    && L_->xtype != CHOLMOD_PATTERN
    && L_->minor < static_cast< size_t >(numberRows_)) {
    /* The supernodal method does no pivoting and gives up on the first non
       positive pivot.  Redo this factorization - and every later one - with
       the simplicial LDL' method, bounding tiny pivots instead. */
    forceSimplicial_ = true;
    c_->supernodal = CHOLMOD_SIMPLICIAL;
    c_->dbound = 1.0e-20;
    cholmod_free_factor(&L_, c_);
    L_ = cholmod_analyze(&A, c_);
    if (L_)
      cholmod_factorize(&A, L_, c_);
  }
  if (!L_ || c_->status < CHOLMOD_OK || L_->xtype == CHOLMOD_PATTERN) {
    // Hard failure (out of memory, too large, ...) - L_ is unusable.
    printf("CHOLMOD: factorization failed (status %d)\n", c_->status);
    return -1;
  }
  choleskyCondition_ = 1.0;
  bool cleanCholesky;
  if (model_->numberIterations() < 2000)
    cleanCholesky = true;
  else
    cleanCholesky = false;
  if (cleanCholesky) {
    //drop fresh makes some formADAT easier
    //int oldDropped=numberRowsDropped_;
    if (newDropped || numberRowsDropped_) {
      //std::cout <<"Rank "<<numberRows_-newDropped<<" ( "<<
      //  newDropped<<" dropped)";
      //if (newDropped>oldDropped)
      //std::cout<<" ( "<<newDropped-oldDropped<<" dropped this time)";
      //std::cout<<std::endl;
      newDropped = 0;
      for (int i = 0; i < numberRows_; i++) {
        int dropped = rowsDropped[i];
        rowsDropped_[i] = (char)dropped;
        if (dropped == 2) {
          //dropped this time
          rowsDropped[newDropped++] = i;
          rowsDropped_[i] = 0;
        }
      }
      numberRowsDropped_ = newDropped;
      newDropped = -(2 + newDropped);
    }
  } else {
    if (newDropped) {
      newDropped = 0;
      for (int i = 0; i < numberRows_; i++) {
        int dropped = rowsDropped[i];
        rowsDropped_[i] = (char)dropped;
        if (dropped == 2) {
          //dropped this time
          rowsDropped[newDropped++] = i;
          rowsDropped_[i] = 1;
        }
      }
    }
    numberRowsDropped_ += newDropped;
    if (numberRowsDropped_ && 0) {
      std::cout << "Rank " << numberRows_ - numberRowsDropped_ << " ( " << numberRowsDropped_ << " dropped)";
      if (newDropped) {
        std::cout << " ( " << newDropped << " dropped this time)";
      }
      std::cout << std::endl;
    }
  }
  status_ = 0;
  return newDropped;
}
#else
/* Factorize - filling in rowsDropped and returning number dropped */
int ClpCholeskyUfl::factorize(const double *diagonal, int *rowsDropped)
{
  return ClpCholeskyBase::factorize(diagonal, rowsDropped);
}
#endif

#ifdef CLP_HAS_CHOLMOD
/* Uses factorization to solve. */
void ClpCholeskyUfl::solve(double *region)
{
  /* cholmod_solve2 reuses B_/X_/Y_/E_ across calls; cholmod_solve used to
     allocate and free a dense n-vector (plus internal workspace) on every
     single solve, and the barrier does several of those per iteration. */
  if (!B_)
    B_ = cholmod_allocate_dense(numberRows_, 1, numberRows_, CHOLMOD_REAL, c_);
  CoinMemcpyN(region, numberRows_, (double *)B_->x);
  cholmod_solve2(CHOLMOD_A, L_, B_, NULL, &X_, NULL, &Y_, &E_, c_);
  CoinMemcpyN((double *)X_->x, numberRows_, region);
}
#else
void ClpCholeskyUfl::solve(double *region)
{
  ClpCholeskyBase::solve(region);
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
