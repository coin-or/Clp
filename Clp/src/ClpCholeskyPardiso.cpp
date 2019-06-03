/* $Id: ClpCholeskyPardiso.cpp 1723 2011-04-17 15:07:10Z forrest $ */
// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#ifdef PARDISO_BARRIER

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"

#include "ClpInterior.hpp"
#include "ClpCholeskyPardiso.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpCholeskyPardiso::ClpCholeskyPardiso(int denseThreshold)
  : ClpCholeskyBase(denseThreshold)
{
  type_ = 17;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpCholeskyPardiso::ClpCholeskyPardiso(const ClpCholeskyPardiso &rhs)
  : ClpCholeskyBase(rhs)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpCholeskyPardiso::~ClpCholeskyPardiso()
{
  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  int phase = -1; /* Release internal memory. */
  void *voidParameters = reinterpret_cast< void * >(doubleParameters_);
  MKL_INT mtype = -2; /* Real symmetric matrix */
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  memset(integerParameters_, 0, sizeof(integerParameters_));
  MKL_INT maxfct, mnum, error, msglvl = 0;
  /* Auxiliary variables. */
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
  PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
    &numberRows_, sparseFactor_,
    choleskyStart_, choleskyRow_, &idum,
    &nrhs, integerParameters_, &msglvl, &ddum, &ddum, &error);
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpCholeskyPardiso &
ClpCholeskyPardiso::operator=(const ClpCholeskyPardiso &rhs)
{
  if (this != &rhs) {
    ClpCholeskyBase::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase *ClpCholeskyPardiso::clone() const
{
  return new ClpCholeskyPardiso(*this);
}
/* Orders rows and saves pointer to matrix.and model */
int ClpCholeskyPardiso::order(ClpInterior *model)
{
  numberRows_ = model->numberRows();
  rowsDropped_ = new char[numberRows_];
  memset(rowsDropped_, 0, numberRows_);
  numberRowsDropped_ = 0;
  model_ = model;
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
  int *which = new int[numberRows_];
  int *used = new int[numberRows_ + 1];
  CoinZeroN(used, numberRows_);
  int iRow;
  sizeFactor_ = 0;
  int numberColumns = model->numberColumns();
  int numberDense = 0;
  denseThreshold_ = 0;
  if (denseThreshold_ > 0) {
    delete[] whichDense_;
    delete[] denseColumn_;
    delete dense_;
    whichDense_ = new char[numberColumns];
    int iColumn;
    used[numberRows_] = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int length = columnLength[iColumn];
      used[length] += 1;
    }
    int nLong = 0;
    int stop = CoinMax(denseThreshold_ / 2, 100);
    for (iRow = numberRows_; iRow >= stop; iRow--) {
      if (used[iRow])
        COIN_DETAIL_PRINT(printf("%d columns are of length %d\n", used[iRow], iRow));
      nLong += used[iRow];
      if (nLong > 50 || nLong > (numberColumns >> 2))
        break;
    }
    CoinZeroN(used, numberRows_);
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnLength[iColumn] < denseThreshold_) {
        whichDense_[iColumn] = 0;
      } else {
        whichDense_[iColumn] = 1;
        numberDense++;
      }
    }
    if (!numberDense || numberDense > 100) {
      // free
      delete[] whichDense_;
      whichDense_ = NULL;
      denseColumn_ = NULL;
      dense_ = NULL;
    } else {
      // space for dense columns
      denseColumn_ = new double[numberDense * numberRows_];
      // dense cholesky
      dense_ = new ClpCholeskyDense();
      dense_->reserveSpace(NULL, numberDense);
      COIN_DETAIL_PRINT(printf("Taking %d columns as dense\n", numberDense));
    }
  }
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
        if (!whichDense_ || !whichDense_[iColumn]) {
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
    choleskyRow_ = new int[sizeFactor_];
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
        if (!whichDense_ || !whichDense_[iColumn]) {
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
  permuteInverse_ = new int[numberRows_];
  permute_ = new int[numberRows_];
  // Assume void * same size as double
  void *voidParameters = reinterpret_cast< void * >(doubleParameters_);
  MKL_INT mtype = -2; /* Real symmetric matrix */
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  memset(integerParameters_, 0, sizeof(integerParameters_));
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  integerParameters_[0] = 1; /* No solver default */
  integerParameters_[1] = 2; /* Fill-in reordering from METIS */
  integerParameters_[3] = 0; /* No iterative-direct algorithm */
  integerParameters_[4] = 0; /* No user fill-in reducing permutation */
  integerParameters_[5] = 0; /* Write solution into x */
  integerParameters_[6] = 0; /* Number of refinements done */
  integerParameters_[7] = 2; /* Max numbers of iterative refinement steps */
  integerParameters_[8] = 0; /* Not in use */
  integerParameters_[9] = 13; /* Perturb the pivot elements with 1E-13 */
  integerParameters_[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  integerParameters_[11] = 0; /* Not in use */
  integerParameters_[12] = 1; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try integerParameters_[12] = 1 in case of inappropriate accuracy */
  integerParameters_[13] = 0; /* Output: Number of perturbed pivots */
  integerParameters_[14] = 0; /* Not in use */
  integerParameters_[15] = 0; /* Not in use */
  integerParameters_[16] = 0; /* Not in use */
  integerParameters_[17] = -1; /* Output: Number of nonzeros in the factor LU */
  integerParameters_[18] = -1; /* Output: Mflops for LU factorization */
  integerParameters_[19] = 0; /* Output: Numbers of CG Iterations */
  integerParameters_[19] = 1; // for 2x2 pivoting
  integerParameters_[34] = 1; /* PARDISO use C-style indexing for ia and ja arrays */
  integerParameters_[55] = 1; /* Pivoting control */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 0; /* Print statistical information in file */
  error = 0; /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  memset(voidParameters, 0, sizeof(doubleParameters_));
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  // temp - fill elements
  for (int i = 0; i < sizeFactor_; i++)
    sparseFactor_[i] = 1.0;
#define ALWAYS_ORDER
#ifdef ALWAYS_ORDER
  error = 0;
#else
  phase = 11;
  PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
    &numberRows_, sparseFactor_,
    choleskyStart_, choleskyRow_, &idum,
    &nrhs, integerParameters_, &msglvl, &ddum, &ddum, &error);
#endif
  if (error != 0) {
    printf("ERROR during symbolic factorization: %d\n", error);
    return 1;
  }
  printf("Number of nonzeros in factors = %d\n", integerParameters_[17]);
  printf("Number of factorization MFLOPS = %d\n", integerParameters_[18]);
  lastNumberDropped_ = 0;
  // Modify gamma and delta if 0.0
  if (!model->gamma() && !model->delta()) {
    model->setGamma(5.0e-5);
    model->setDelta(5.0e-5);
  }
  return 0;
}
/* Does Symbolic factorization given permutation.
   This is called immediately after order.  If user provides this then
   user must provide factorize and solve.  Otherwise the default factorization is used
   returns non-zero if not enough memory */
int ClpCholeskyPardiso::symbolic()
{
  return 0;
}
/* Factorize - filling in rowsDropped and returning number dropped */
int ClpCholeskyPardiso::factorize(const double *diagonal, int *rowsDropped)
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
  double largest;
  //int numberDense = 0;
  //if (dense_)
  //   numberDense = dense_->numberRows();
  //perturbation
  double perturbation = model_->diagonalPerturbation() * model_->diagonalNorm();
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
  if (whichDense_) {
    double *denseDiagonal = dense_->diagonal();
    double *dense = denseColumn_;
    int iDense = 0;
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (whichDense_[iColumn]) {
        denseDiagonal[iDense++] = 1.0 / diagonal[iColumn];
        CoinZeroN(dense, numberRows_);
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn] + columnLength[iColumn];
        for (CoinBigIndex j = start; j < end; j++) {
          int jRow = row[j];
          dense[jRow] = element[j];
        }
        dense += numberRows_;
      }
    }
  }
  double delta2 = model_->delta(); // add delta*delta to diagonal
  delta2 *= delta2;
  int *whichNon = new int[numberRows_];
  int nNon = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    //for (int i=0;i<numberRows_;i++)
    //assert(!work[i]);
    for (int i = 0; i < nNon; i++)
      work[whichNon[i]] = 0.0;
    nNon = 1;
    whichNon[0] = iRow;
    double *put = sparseFactor_ + choleskyStart_[iRow];
    int *which = choleskyRow_ + choleskyStart_[iRow];
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
              if (!work[jRow])
                whichNon[nNon++] = jRow;
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
      put[0] = 1.0e12;
    }
  }
  delete[] whichNon;
  //check sizes
  double largest2 = maximumAbsElement(sparseFactor_, sizeFactor_);
  largest2 *= 1.0e-20;
  largest = CoinMin(largest2, 1.0e-11);
  int numberDroppedBefore = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    int dropped = rowsDropped_[iRow];
    // Move to int array
    rowsDropped[iRow] = dropped;
    if (!dropped) {
      CoinBigIndex start = choleskyStart_[iRow];
      double diagonal = sparseFactor_[start];
      if (diagonal > largest2) {
        sparseFactor_[start] = diagonal + perturbation;
      } else {
        sparseFactor_[start] = diagonal + perturbation;
        rowsDropped[iRow] = 2;
        numberDroppedBefore++;
      }
    }
  }
  //integerParameters_[0] = 0;
  int maxfct = 1; /* Maximum number of numerical factorizations. */
  int mnum = 1; /* Which factorization to use. */
  int msglvl = 0; /* Print statistical information in file */
  int error = 0; /* Initialize error flag */
  int phase = 22;
#ifdef ALWAYS_ORDER
  {
    int phase = -1; /* Release internal memory. */
    void *voidParameters = reinterpret_cast< void * >(doubleParameters_);
    MKL_INT mtype = -2; /* Real symmetric matrix */
    MKL_INT nrhs = 1; /* Number of right hand sides. */
    memset(integerParameters_, 0, sizeof(integerParameters_));
    MKL_INT maxfct, mnum, error, msglvl = 0;
    /* Auxiliary variables. */
    double ddum; /* Double dummy */
    MKL_INT idum; /* Integer dummy. */
    PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
      &numberRows_, sparseFactor_,
      choleskyStart_, choleskyRow_, &idum,
      &nrhs, integerParameters_, &msglvl, &ddum, &ddum, &error);
    memset(integerParameters_, 0, sizeof(integerParameters_));
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    integerParameters_[0] = 1; /* No solver default */
    integerParameters_[1] = 2; /* Fill-in reordering from METIS */
    integerParameters_[3] = 0; /* No iterative-direct algorithm */
    integerParameters_[4] = 0; /* No user fill-in reducing permutation */
    integerParameters_[5] = 0; /* Write solution into x */
    integerParameters_[6] = 0; /* Number of refinements done */
    integerParameters_[7] = 2; /* Max numbers of iterative refinement steps */
    integerParameters_[8] = 0; /* Not in use */
    integerParameters_[9] = 13; /* Perturb the pivot elements with 1E-13 */
    integerParameters_[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
    integerParameters_[11] = 0; /* Not in use */
    integerParameters_[12] = 1; /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try integerParameters_[12] = 1 in case of inappropriate accuracy */
    integerParameters_[13] = 0; /* Output: Number of perturbed pivots */
    integerParameters_[14] = 0; /* Not in use */
    integerParameters_[15] = 0; /* Not in use */
    integerParameters_[16] = 0; /* Not in use */
    integerParameters_[17] = -1; /* Output: Number of nonzeros in the factor LU */
    integerParameters_[18] = -1; /* Output: Mflops for LU factorization */
    integerParameters_[19] = 0; /* Output: Numbers of CG Iterations */
    integerParameters_[19] = 1; // for 2x2 pivoting
    integerParameters_[34] = 1; /* PARDISO use C-style indexing for ia and ja arrays */
    integerParameters_[55] = 1; /* Pivoting control */
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum = 1; /* Which factorization to use. */
    msglvl = 0; /* Print statistical information in file */
    error = 0; /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    memset(voidParameters, 0, sizeof(doubleParameters_));
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    // pack down
    CoinBigIndex numberElements = 0;
    for (int iRow = 0; iRow < numberRows_; iRow++) {
      CoinBigIndex start = choleskyStart_[iRow];
      CoinBigIndex end = choleskyStart_[iRow + 1];
      choleskyStart_[iRow] = numberElements;
      for (CoinBigIndex j = start; j < end; j++) {
        if (sparseFactor_[j]) {
          choleskyRow_[numberElements] = choleskyRow_[j];
          sparseFactor_[numberElements++] = sparseFactor_[j];
        }
      }
    }
    choleskyStart_[numberRows_] = numberElements;
    phase = 11;
    PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
      &numberRows_, sparseFactor_,
      choleskyStart_, choleskyRow_, &idum,
      &nrhs, integerParameters_, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
      printf("ERROR during symbolic factorization: %d\n", error);
      return 1;
    }
    printf("Number of nonzeros in factors = %d\n", integerParameters_[17]);
    printf("Number of factorization MFLOPS = %d\n", integerParameters_[18]);
  }
#endif
  if (numberRowsDropped_ > lastNumberDropped_ && false) {
    phase = 12; //reorder
    CoinBigIndex numberElements = 0;
    for (int iRow = 0; iRow < numberRows_; iRow++) {
      CoinBigIndex start = choleskyStart_[iRow];
      CoinBigIndex end = choleskyStart_[iRow + 1];
      choleskyStart_[iRow] = numberElements;
      if (rowsDropped_[iRow] != 2) {
        for (CoinBigIndex j = start; j < end; j++) {
          choleskyRow_[numberElements] = choleskyRow_[j];
          sparseFactor_[numberElements++] = sparseFactor_[j];
        }
      } else {
        rowsDropped_[iRow] = 1;
        for (CoinBigIndex j = start; j < end; j++) {
          if (sparseFactor_[j]) {
            choleskyRow_[numberElements] = choleskyRow_[j];
            sparseFactor_[numberElements++] = sparseFactor_[j];
          }
        }
      }
    }
    choleskyStart_[numberRows_] = numberElements;
    lastNumberDropped_ = numberRowsDropped_;
  }
  MKL_INT mtype = -2; /* Real symmetric matrix */
  void *voidParameters = reinterpret_cast< void * >(doubleParameters_);
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  /* Auxiliary variables. */
  double ddum; /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
  PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
    &numberRows_, sparseFactor_,
    choleskyStart_, choleskyRow_, &idum,
    &nrhs, integerParameters_, &msglvl, &ddum, &ddum, &error);
  printf("%d positive, %d negative\n", integerParameters_[21],
    integerParameters_[22]);
  // use some other array?
  double *originalDiagonal = new double[numberRows_];
  pardiso_getdiag(voidParameters, work, originalDiagonal, &mnum, &error);
  delete[] originalDiagonal;
  largest = 1.0e-11;
  double smallest = COIN_DBL_MAX;
  int newDropped = 0;
  for (int i = 0; i < numberRows_; i++) {
    double value = work[i];
    if (value != 1.0e100) {
      //if (value>1.0e13)
      //printf("large diagonal %g at %d\n",value,i);
      largest = CoinMax(largest, value);
      smallest = CoinMin(smallest, value);
    } else if (!rowsDropped_[i]) {
      rowsDropped_[i] = 2;
      rowsDropped[i] = 2;
      printf("%d dropped\n", i);
      newDropped++;
      numberRowsDropped_++;
    }
  }
  if (smallest <= 0.0)
    printf("small negative diagonal %g\n", smallest);
  delete[] work;
  if (model_->messageHandler()->logLevel() > 1)
    std::cout << "Cholesky - largest " << largest << " smallest " << smallest << std::endl;
  choleskyCondition_ = largest / smallest;
  status_ = 0;
  return newDropped;
}
/* Uses factorization to solve. */
void ClpCholeskyPardiso::solve(double *region)
{
  /* -------------------------------------------------------------------- */
  /* .. Back substitution and iterative refinement. */
  /* -------------------------------------------------------------------- */
  int phase = 33;
  integerParameters_[7] = 2; /* Max numbers of iterative refinement steps. */
  double *x = new double[numberRows_];
  memcpy(x, region, numberRows_ * sizeof(double));
  MKL_INT mtype = -2; /* Real symmetric matrix */
  void *voidParameters = reinterpret_cast< void * >(doubleParameters_);
  MKL_INT nrhs = 1; /* Number of right hand sides. */
  MKL_INT maxfct = 1, mnum = 1, error, msglvl = 0;
  /* Auxiliary variables. */
  MKL_INT idum; /* Integer dummy. */
  PARDISO(voidParameters, &maxfct, &mnum, &mtype, &phase,
    &numberRows_, sparseFactor_,
    choleskyStart_, choleskyRow_, &idum,
    &nrhs, integerParameters_, &msglvl, x, region, &error);
  if (error != 0) {
    printf("ERROR during solution: %d\n", error);
    exit(3);
  }
  delete[] x;
}
/* -------------------------------------------------------------------- */
/* .. mkl_pardiso_pivot: This function replace appropriate function     */
/* ..                    in pardiso solver                              */
/* -------------------------------------------------------------------- */
// up rows dropped
extern "C" {
int mkl_pardiso_pivot(const double *aii, double *bii, const double *eps)
{
  //if ( (*bii > *eps) || ( *bii < -*eps ) )
  if (*bii > *eps)
    return 0;
  if (*bii > 0)
    *bii = *eps;
  else
    *bii = -*eps;
  *bii = 1.0e100;
  return 1;
}
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
