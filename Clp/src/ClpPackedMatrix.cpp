/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
//#define THREAD
//#define FAKE_CILK
#if ABOCA_LITE
// 1 is not owner of abcState_
#define ABCSTATE_LITE 1
#endif
#include "ClpSimplex.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpFactorization.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#endif
// at end to get min/max!
#include "ClpPackedMatrix.hpp"
#include "ClpMessage.hpp"
#ifdef INTEL_MKL
#include "mkl_spblas.h"
#endif
//#define DO_CHECK_FLAGS 1
//=============================================================================
#ifdef COIN_PREFETCH
#if 1
#define coin_prefetch(mem)              \
  __asm__ __volatile__("prefetchnta %0" \
                       :                \
                       : "m"(*(reinterpret_cast< char * >(mem))))
#define coin_prefetch_const(mem)        \
  __asm__ __volatile__("prefetchnta %0" \
                       :                \
                       : "m"(*(reinterpret_cast< const char * >(mem))))
#else
#define coin_prefetch(mem)           \
  __asm__ __volatile__("prefetch %0" \
                       :             \
                       : "m"(*(reinterpret_cast< char * >(mem))))
#define coin_prefetch_const(mem)     \
  __asm__ __volatile__("prefetch %0" \
                       :             \
                       : "m"(*(reinterpret_cast< const char * >(mem))))
#endif
#else
// dummy
#define coin_prefetch(mem)
#define coin_prefetch_const(mem)
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################
static void debug3(int iColumn, double &thisWeight, double pivotSquared, double devex,
  double pivot, double modification, double oldWeight)
{
  //double newWeight=oldWeight+pivotSquared * devex + pivot * modification;
  //if (/*iColumn==34845 && */newWeight!=thisWeight)
  //printf("YY %d %.30g %.30g %.30g %.30g %.30g %.30g\n",
  //	   iColumn,thisWeight,pivotSquared,devex,pivot,modification,oldWeight);
  thisWeight = oldWeight;
  thisWeight += pivotSquared * devex + pivot * modification;
}
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix()
  : ClpMatrixBase()
  , matrix_(NULL)
  , numberActiveColumns_(0)
  , flags_(2)
  , rowCopy_(NULL)
  , columnCopy_(NULL)
{
  setType(1);
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix(const ClpPackedMatrix &rhs)
  : ClpMatrixBase(rhs)
{
#ifdef DO_CHECK_FLAGS
  rhs.checkFlags(0);
#endif
#ifndef COIN_SPARSE_MATRIX
  // Guaranteed no gaps or small elements
  matrix_ = new CoinPackedMatrix(*(rhs.matrix_), -1, 0);
  flags_ = rhs.flags_ & (~0x02);
#else
  // Gaps & small elements preserved
  matrix_ = new CoinPackedMatrix(*(rhs.matrix_), 0, 0);
  flags_ = rhs.flags_;
  if (matrix_->hasGaps())
    flags_ |= 0x02;
#endif
  numberActiveColumns_ = rhs.numberActiveColumns_;
  int numberRows = matrix_->getNumRows();
  if (rhs.rhsOffset_ && numberRows) {
    rhsOffset_ = ClpCopyOfArray(rhs.rhsOffset_, numberRows);
  } else {
    rhsOffset_ = NULL;
  }
  if (rhs.rowCopy_) {
    assert((flags_ & 4) != 0);
    rowCopy_ = new ClpPackedMatrix2(*rhs.rowCopy_);
  } else {
    rowCopy_ = NULL;
  }
  if (rhs.columnCopy_) {
    assert((flags_ & (8 + 16)) == 8 + 16);
    columnCopy_ = new ClpPackedMatrix3(*rhs.columnCopy_);
  } else {
    columnCopy_ = NULL;
  }
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix(CoinPackedMatrix *rhs)
  : ClpMatrixBase()
{
  matrix_ = rhs;
  flags_ = ((matrix_->hasGaps()) ? 0x02 : 0);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  columnCopy_ = NULL;
  setType(1);
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}

ClpPackedMatrix::ClpPackedMatrix(const CoinPackedMatrix &rhs)
  : ClpMatrixBase()
{
#ifndef COIN_SPARSE_MATRIX
  matrix_ = new CoinPackedMatrix(rhs, -1, 0);
  flags_ = 0;
#else
  matrix_ = new CoinPackedMatrix(rhs, 0, 0);
  flags_ = ((matrix_->hasGaps()) ? 0x02 : 0);
#endif
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  columnCopy_ = NULL;
  setType(1);
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPackedMatrix::~ClpPackedMatrix()
{
  delete matrix_;
  delete rowCopy_;
  delete columnCopy_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPackedMatrix &
ClpPackedMatrix::operator=(const ClpPackedMatrix &rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete matrix_;
#ifndef COIN_SPARSE_MATRIX
    matrix_ = new CoinPackedMatrix(*(rhs.matrix_), -1, 0);
    flags_ = rhs.flags_ & (~0x02);
#else
    matrix_ = new CoinPackedMatrix(*(rhs.matrix_), 0, 0);
    flags_ = rhs.flags_;
    if (matrix_->hasGaps())
      flags_ |= 0x02;
#endif
    numberActiveColumns_ = rhs.numberActiveColumns_;
    delete rowCopy_;
    delete columnCopy_;
    if (rhs.rowCopy_) {
      assert((flags_ & 4) != 0);
      rowCopy_ = new ClpPackedMatrix2(*rhs.rowCopy_);
    } else {
      rowCopy_ = NULL;
    }
    if (rhs.columnCopy_) {
      assert((flags_ & (8 + 16)) == 8 + 16);
      columnCopy_ = new ClpPackedMatrix3(*rhs.columnCopy_);
    } else {
      columnCopy_ = NULL;
    }
#ifdef DO_CHECK_FLAGS
    checkFlags(0);
#endif
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase *ClpPackedMatrix::clone() const
{
  return new ClpPackedMatrix(*this);
}
// Copy contents - resizing if necessary - otherwise re-use memory
void ClpPackedMatrix::copy(const ClpPackedMatrix *rhs)
{
  //*this = *rhs;
  assert(numberActiveColumns_ == rhs->numberActiveColumns_);
  assert(matrix_->isColOrdered() == rhs->matrix_->isColOrdered());
  matrix_->copyReuseArrays(*rhs->matrix_);
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase *
ClpPackedMatrix::subsetClone(int numberRows, const int *whichRows,
  int numberColumns,
  const int *whichColumns) const
{
  return new ClpPackedMatrix(*this, numberRows, whichRows,
    numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpPackedMatrix::ClpPackedMatrix(
  const ClpPackedMatrix &rhs,
  int numberRows, const int *whichRows,
  int numberColumns, const int *whichColumns)
  : ClpMatrixBase(rhs)
{
  matrix_ = new CoinPackedMatrix(*(rhs.matrix_), numberRows, whichRows,
    numberColumns, whichColumns);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  flags_ = rhs.flags_ & (~0x02); // no gaps
  columnCopy_ = NULL;
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
ClpPackedMatrix::ClpPackedMatrix(
  const CoinPackedMatrix &rhs,
  int numberRows, const int *whichRows,
  int numberColumns, const int *whichColumns)
  : ClpMatrixBase()
{
  matrix_ = new CoinPackedMatrix(rhs, numberRows, whichRows,
    numberColumns, whichColumns);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  flags_ = 0; // no gaps
  columnCopy_ = NULL;
  setType(1);
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase *
ClpPackedMatrix::reverseOrderedCopy() const
{
  ClpPackedMatrix *copy = new ClpPackedMatrix();
  copy->matrix_ = new CoinPackedMatrix();
  copy->matrix_->setExtraGap(0.0);
  copy->matrix_->setExtraMajor(0.0);
  copy->matrix_->reverseOrderedCopyOf(*matrix_);
  //copy->matrix_->removeGaps();
  copy->numberActiveColumns_ = copy->matrix_->getNumCols();
  copy->flags_ = flags_ & (~0x02); // no gaps
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
  return copy;
}
//unscaled versions
void ClpPackedMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
  int iRow, iColumn;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  //memset(y,0,matrix_->getNumRows()*sizeof(double));
  assert(((flags_ & 0x02) != 0) == matrix_->hasGaps());
  if (!(flags_ & 2)) {
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      double value = x[iColumn];
      if (value) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn + 1];
        value *= scalar;
        for (j = start; j < end; j++) {
          iRow = row[j];
          y[iRow] += value * elementByColumn[j];
        }
      }
    }
  } else {
    const int *columnLength = matrix_->getVectorLengths();
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      double value = x[iColumn];
      if (value) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        value *= scalar;
        for (j = start; j < end; j++) {
          iRow = row[j];
          y[iRow] += value * elementByColumn[j];
        }
      }
    }
  }
}
#if ABOCA_LITE
static void
transposeTimesBit(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT y = info.spare;
  const double *COIN_RESTRICT x = info.work;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  CoinBigIndex start = columnStart[first];
  for (int iColumn = first; iColumn < last; iColumn++) {
    CoinBigIndex j;
    CoinBigIndex next = columnStart[iColumn + 1];
    double value = y[iColumn];
    for (j = start; j < next; j++) {
      int jRow = row[j];
      value -= x[jRow] * elementByColumn[j];
    }
    start = next;
    y[iColumn] = value;
  }
}
#endif
void ClpPackedMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
  int iColumn;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  if (!(flags_ & 2)) {
    if (scalar == -1.0) {
#if ABOCA_LITE
      int numberThreads = abcState();
      if (numberThreads) {
        clpTempInfo info[ABOCA_LITE];
        int chunk = (numberActiveColumns_ + numberThreads - 1) / numberThreads;
        int n = 0;
        for (int i = 0; i < numberThreads; i++) {
          info[i].spare = y;
          info[i].work = const_cast< double * >(x);
          info[i].startColumn = n;
          info[i].element = elementByColumn;
          info[i].start = columnStart;
          info[i].row = row;
          info[i].numberToDo = CoinMin(chunk, numberActiveColumns_ - n);
          n += chunk;
        }
        for (int i = 0; i < numberThreads; i++) {
          cilk_spawn transposeTimesBit(info[i]);
        }
        cilk_sync;
      } else {
#endif
        CoinBigIndex start = columnStart[0];
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          CoinBigIndex j;
          CoinBigIndex next = columnStart[iColumn + 1];
          double value = y[iColumn];
          for (j = start; j < next; j++) {
            int jRow = row[j];
            value -= x[jRow] * elementByColumn[j];
          }
          start = next;
          y[iColumn] = value;
        }
#if ABOCA_LITE
      }
#endif
    } else {
      CoinBigIndex start = columnStart[0];
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        CoinBigIndex j;
        CoinBigIndex next = columnStart[iColumn + 1];
        double value = 0.0;
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value += x[jRow] * elementByColumn[j];
        }
        start = next;
        y[iColumn] += value * scalar;
      }
    }
  } else {
    const int *columnLength = matrix_->getVectorLengths();
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      double value = 0.0;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      for (j = start; j < end; j++) {
        int jRow = row[j];
        value += x[jRow] * elementByColumn[j];
      }
      y[iColumn] += value * scalar;
    }
  }
}
void ClpPackedMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *COIN_RESTRICT rowScale,
  const double *COIN_RESTRICT columnScale) const
{
  if (rowScale) {
    int iRow, iColumn;
    // get matrix data pointers
    const int *COIN_RESTRICT row = matrix_->getIndices();
    const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
    const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
    if (!(flags_ & 2)) {
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        double value = x[iColumn];
        if (value) {
          // scaled
          value *= scalar * columnScale[iColumn];
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = columnStart[iColumn + 1];
          CoinBigIndex j;
          for (j = start; j < end; j++) {
            iRow = row[j];
            y[iRow] += value * elementByColumn[j] * rowScale[iRow];
          }
        }
      }
    } else {
      const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        double value = x[iColumn];
        if (value) {
          // scaled
          value *= scalar * columnScale[iColumn];
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex end = start + columnLength[iColumn];
          CoinBigIndex j;
          for (j = start; j < end; j++) {
            iRow = row[j];
            y[iRow] += value * elementByColumn[j] * rowScale[iRow];
          }
        }
      }
    }
  } else {
    times(scalar, x, y);
  }
}
void ClpPackedMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *COIN_RESTRICT rowScale,
  const double *COIN_RESTRICT columnScale,
  double *COIN_RESTRICT spare) const
{
  if (rowScale) {
    int iColumn;
    // get matrix data pointers
    const int *COIN_RESTRICT row = matrix_->getIndices();
    const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
    const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
    const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
    if (!spare) {
      if (!(flags_ & 2)) {
        CoinBigIndex start = columnStart[0];
        if (scalar == -1.0) {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            CoinBigIndex j;
            CoinBigIndex next = columnStart[iColumn + 1];
            double value = 0.0;
            // scaled
            for (j = start; j < next; j++) {
              int jRow = row[j];
              value += x[jRow] * elementByColumn[j] * rowScale[jRow];
            }
            start = next;
            y[iColumn] -= value * columnScale[iColumn];
          }
        } else {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            CoinBigIndex j;
            CoinBigIndex next = columnStart[iColumn + 1];
            double value = 0.0;
            // scaled
            for (j = start; j < next; j++) {
              int jRow = row[j];
              value += x[jRow] * elementByColumn[j] * rowScale[jRow];
            }
            start = next;
            y[iColumn] += value * scalar * columnScale[iColumn];
          }
        }
      } else {
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          CoinBigIndex j;
          double value = 0.0;
          // scaled
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int jRow = row[j];
            value += x[jRow] * elementByColumn[j] * rowScale[jRow];
          }
          y[iColumn] += value * scalar * columnScale[iColumn];
        }
      }
    } else {
      // can use spare region
      int iRow;
      int numberRows = matrix_->getNumRows();
      for (iRow = 0; iRow < numberRows; iRow++) {
        double value = x[iRow];
        if (value)
          spare[iRow] = value * rowScale[iRow];
        else
          spare[iRow] = 0.0;
      }
      if (!(flags_ & 2)) {
        CoinBigIndex start = columnStart[0];
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          CoinBigIndex j;
          CoinBigIndex next = columnStart[iColumn + 1];
          double value = 0.0;
          // scaled
          for (j = start; j < next; j++) {
            int jRow = row[j];
            value += spare[jRow] * elementByColumn[j];
          }
          start = next;
          y[iColumn] += value * scalar * columnScale[iColumn];
        }
      } else {
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          CoinBigIndex j;
          double value = 0.0;
          // scaled
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int jRow = row[j];
            value += spare[jRow] * elementByColumn[j];
          }
          y[iColumn] += value * scalar * columnScale[iColumn];
        }
      }
      // no need to zero out
      //for (iRow=0;iRow<numberRows;iRow++)
      //spare[iRow] = 0.0;
    }
  } else {
    transposeTimes(scalar, x, y);
  }
}
#if ABOCA_LITE
static void
transposeTimesSubsetBit(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT y = info.spare;
  const double *COIN_RESTRICT x = info.work;
  const int *COIN_RESTRICT which = info.which;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  for (int i = first; i < last; i++) {
    int iColumn = which[i];
    CoinBigIndex j;
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex next = columnStart[iColumn + 1];
    double value = 0.0;
    for (j = start; j < next; j++) {
      int jRow = row[j];
      value += x[jRow] * elementByColumn[j];
    }
    start = next;
    y[iColumn] -= value;
  }
}
#endif
void ClpPackedMatrix::transposeTimesSubset(int number,
  const int *which,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *COIN_RESTRICT rowScale,
  const double *COIN_RESTRICT columnScale,
  double *COIN_RESTRICT spare) const
{
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  if (!spare || !rowScale) {
    if (rowScale) {
      for (int jColumn = 0; jColumn < number; jColumn++) {
        int iColumn = which[jColumn];
        CoinBigIndex j;
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex next = columnStart[iColumn + 1];
        double value = 0.0;
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value += x[jRow] * elementByColumn[j] * rowScale[jRow];
        }
        y[iColumn] -= value * columnScale[iColumn];
      }
    } else {
#if ABOCA_LITE
      int numberThreads = abcState();
      if (numberThreads) {
        clpTempInfo info[ABOCA_LITE];
        int chunk = (number + numberThreads - 1) / numberThreads;
        int n = 0;
        for (int i = 0; i < numberThreads; i++) {
          info[i].spare = y;
          info[i].work = const_cast< double * >(x);
          info[i].which = const_cast< int * >(which);
          info[i].startColumn = n;
          info[i].element = elementByColumn;
          info[i].start = columnStart;
          info[i].row = row;
          info[i].numberToDo = CoinMin(chunk, number - n);
          n += chunk;
        }
        for (int i = 0; i < numberThreads; i++) {
          cilk_spawn transposeTimesSubsetBit(info[i]);
        }
        cilk_sync;
      } else {
#endif
        for (int jColumn = 0; jColumn < number; jColumn++) {
          int iColumn = which[jColumn];
          CoinBigIndex j;
          CoinBigIndex start = columnStart[iColumn];
          CoinBigIndex next = columnStart[iColumn + 1];
          double value = 0.0;
          for (j = start; j < next; j++) {
            int jRow = row[j];
            value += x[jRow] * elementByColumn[j];
          }
          y[iColumn] -= value;
        }
#if ABOCA_LITE
      }
#endif
    }
  } else {
    // can use spare region
    int iRow;
    int numberRows = matrix_->getNumRows();
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value = x[iRow];
      if (value)
        spare[iRow] = value * rowScale[iRow];
      else
        spare[iRow] = 0.0;
    }
    for (int jColumn = 0; jColumn < number; jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex next = columnStart[iColumn + 1];
      double value = 0.0;
      for (j = start; j < next; j++) {
        int jRow = row[j];
        value += spare[jRow] * elementByColumn[j];
      }
      y[iColumn] -= value * columnScale[iColumn];
    }
  }
}
/* Return <code>x * A in <code>z</code>.
	Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix::transposeTimes(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
#if 0 //def COIN_DEVELOP
     if (zeroTolerance != 1.0e-13) {
          printf("small element in matrix - zero tolerance %g\n", zeroTolerance);
     }
#endif
  int numberRows = model->numberRows();
  ClpPackedMatrix *rowCopy = static_cast< ClpPackedMatrix * >(model->rowCopy());
  bool packed = rowArray->packedMode();
  double factor = (numberRows < 100) ? 0.25 : 0.35;
  factor = 0.5;
  // We may not want to do by row if there may be cache problems
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberActiveColumns_ * sizeof(double) > 1000000) {
    if (numberRows * 10 < numberActiveColumns_)
      factor *= 0.333333333;
    else if (numberRows * 4 < numberActiveColumns_)
      factor *= 0.5;
    else if (numberRows * 2 < numberActiveColumns_)
      factor *= 0.66666666667;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  // if not packed then bias a bit more towards by column
  if (!packed)
    factor *= 0.9;
  assert(!y->getNumElements());
  double multiplierX = 0.8;
  double factor2 = factor * multiplierX;
  if (packed && rowCopy_ && numberInRowArray > 2 && numberInRowArray > factor2 * numberRows && numberInRowArray < 0.9 * numberRows && 0) {
    rowCopy_->transposeTimes(model, rowCopy->matrix_, rowArray, y, columnArray);
    return;
  }
  if (columnCopy_)
    factor *= 0.7;
  if (numberInRowArray > factor * numberRows || !rowCopy) {
    // do by column
    // If no gaps - can do a bit faster
    if (!(flags_ & 2) || columnCopy_) {
      transposeTimesByColumn(model, scalar,
        rowArray, y, columnArray);
      return;
    }
    int iColumn;
    // get matrix data pointers
    const int *COIN_RESTRICT row = matrix_->getIndices();
    const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
    const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
    const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
    const double *COIN_RESTRICT rowScale = model->rowScale();
#if 0
          ClpPackedMatrix * scaledMatrix = model->clpScaledMatrix();
          if (rowScale && scaledMatrix) {
               rowScale = NULL;
               // get matrix data pointers
               row = scaledMatrix->getIndices();
               columnStart = scaledMatrix->getVectorStarts();
               columnLength = scaledMatrix->getVectorLengths();
               elementByColumn = scaledMatrix->getElements();
          }
#endif
    if (packed) {
      // need to expand pi into y
      assert(y->capacity() >= numberRows);
      double *COIN_RESTRICT piOld = pi;
      pi = y->denseVector();
      const int *COIN_RESTRICT whichRow = rowArray->getIndices();
      int i;
      if (!rowScale) {
        // modify pi so can collapse to one loop
        if (scalar == -1.0) {
          for (i = 0; i < numberInRowArray; i++) {
            int iRow = whichRow[i];
            pi[iRow] = -piOld[i];
          }
        } else {
          for (i = 0; i < numberInRowArray; i++) {
            int iRow = whichRow[i];
            pi[iRow] = scalar * piOld[i];
          }
        }
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          double value = 0.0;
          CoinBigIndex j;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j];
          }
          if (fabs(value) > zeroTolerance) {
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      } else {
#ifdef CLP_INVESTIGATE
        if (model->clpScaledMatrix())
          printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
        // scaled
        // modify pi so can collapse to one loop
        if (scalar == -1.0) {
          for (i = 0; i < numberInRowArray; i++) {
            int iRow = whichRow[i];
            pi[iRow] = -piOld[i] * rowScale[iRow];
          }
        } else {
          for (i = 0; i < numberInRowArray; i++) {
            int iRow = whichRow[i];
            pi[iRow] = scalar * piOld[i] * rowScale[iRow];
          }
        }
        for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
          double value = 0.0;
          CoinBigIndex j;
          const double *columnScale = model->columnScale();
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j];
          }
          value *= columnScale[iColumn];
          if (fabs(value) > zeroTolerance) {
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      }
      // zero out
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = 0.0;
      }
    } else {
      if (!rowScale) {
        if (scalar == -1.0) {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            double value = 0.0;
            CoinBigIndex j;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              value += pi[iRow] * elementByColumn[j];
            }
            if (fabs(value) > zeroTolerance) {
              index[numberNonZero++] = iColumn;
              array[iColumn] = -value;
            }
          }
        } else {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            double value = 0.0;
            CoinBigIndex j;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              value += pi[iRow] * elementByColumn[j];
            }
            value *= scalar;
            if (fabs(value) > zeroTolerance) {
              index[numberNonZero++] = iColumn;
              array[iColumn] = value;
            }
          }
        }
      } else {
#ifdef CLP_INVESTIGATE
        if (model->clpScaledMatrix())
          printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
        // scaled
        if (scalar == -1.0) {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            double value = 0.0;
            CoinBigIndex j;
            const double *columnScale = model->columnScale();
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
            }
            value *= columnScale[iColumn];
            if (fabs(value) > zeroTolerance) {
              index[numberNonZero++] = iColumn;
              array[iColumn] = -value;
            }
          }
        } else {
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            double value = 0.0;
            CoinBigIndex j;
            const double *columnScale = model->columnScale();
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
            }
            value *= scalar * columnScale[iColumn];
            if (fabs(value) > zeroTolerance) {
              index[numberNonZero++] = iColumn;
              array[iColumn] = value;
            }
          }
        }
      }
    }
    columnArray->setNumElements(numberNonZero);
    y->setNumElements(0);
  } else {
    // do by row
    rowCopy->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  }
  if (packed)
    columnArray->setPackedMode(true);
  if (0) {
    columnArray->checkClean();
    int numberNonZero = columnArray->getNumElements();
    ;
    int *index = columnArray->getIndices();
    double *array = columnArray->denseVector();
    int i;
    for (i = 0; i < numberNonZero; i++) {
      int j = index[i];
      double value;
      if (packed)
        value = array[i];
      else
        value = array[j];
      printf("Ti %d %d %g\n", i, j, value);
    }
  }
}
//static int xA=0;
//static int xB=0;
//static int xC=0;
//static int xD=0;
//static double yA=0.0;
//static double yC=0.0;
/* Return <code>x * scalar * A in <code>z</code>.
   Note - If x packed mode - then z packed mode
   This does by column and knows no gaps
   Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix::transposeTimesByColumn(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
  bool packed = rowArray->packedMode();
  // do by column
  int iColumn;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  assert(!y->getNumElements());
  assert(numberActiveColumns_ > 0);
  const ClpPackedMatrix *thisMatrix = this;
#if 0
     ClpPackedMatrix * scaledMatrix = model->clpScaledMatrix();
     if (rowScale && scaledMatrix) {
          rowScale = NULL;
          // get matrix data pointers
          row = scaledMatrix->getIndices();
          columnStart = scaledMatrix->getVectorStarts();
          elementByColumn = scaledMatrix->getElements();
          thisMatrix = scaledMatrix;
          //printf("scaledMatrix\n");
     } else if (rowScale) {
          //printf("no scaledMatrix\n");
     } else {
          //printf("no rowScale\n");
     }
#endif
  if (packed) {
    // need to expand pi into y
    assert(y->capacity() >= model->numberRows());
    double *COIN_RESTRICT piOld = pi;
    pi = y->denseVector();
    const int *COIN_RESTRICT whichRow = rowArray->getIndices();
    int i;
    if (!rowScale) {
      // modify pi so can collapse to one loop
      if (scalar == -1.0) {
        //yA += numberInRowArray;
        for (i = 0; i < numberInRowArray; i++) {
          int iRow = whichRow[i];
          pi[iRow] = -piOld[i];
        }
      } else {
        for (i = 0; i < numberInRowArray; i++) {
          int iRow = whichRow[i];
          pi[iRow] = scalar * piOld[i];
        }
      }
      if (!columnCopy_) {
        if ((model->specialOptions(), 131072) != 0) {
          if (model->spareIntArray_[0] > 0) {
            CoinIndexedVector *spareArray = model->rowArray(3);
            // also do dualColumn stuff
            double *COIN_RESTRICT spare = spareArray->denseVector();
            int *COIN_RESTRICT spareIndex = spareArray->getIndices();
            const double *COIN_RESTRICT reducedCost = model->djRegion(0);
            double multiplier[] = { -1.0, 1.0 };
            double dualT = -model->currentDualTolerance();
            double acceptablePivot = model->spareDoubleArray_[0];
            // We can also see if infeasible or pivoting on free
            double tentativeTheta = 1.0e15;
            double upperTheta = 1.0e31;
            int addSequence = model->numberColumns();
            const unsigned char *COIN_RESTRICT statusArray = model->statusArray() + addSequence;
            int numberRemaining = 0;
            assert(scalar == -1.0);
            for (i = 0; i < numberInRowArray; i++) {
              int iSequence = whichRow[i];
              int iStatus = (statusArray[iSequence] & 3) - 1;
              if (iStatus) {
                double mult = multiplier[iStatus - 1];
                double alpha = piOld[i] * mult;
                double oldValue;
                double value;
                if (alpha > 0.0) {
                  oldValue = reducedCost[iSequence] * mult;
                  value = oldValue - tentativeTheta * alpha;
                  if (value < dualT) {
                    value = oldValue - upperTheta * alpha;
                    if (value < dualT && alpha >= acceptablePivot) {
                      upperTheta = (oldValue - dualT) / alpha;
                      //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
                    }
                    // add to list
                    spare[numberRemaining] = alpha * mult;
                    spareIndex[numberRemaining++] = iSequence + addSequence;
                  }
                }
              }
            }
            numberNonZero = thisMatrix->gutsOfTransposeTimesUnscaled(pi,
              columnArray->getIndices(),
              columnArray->denseVector(),
              model->statusArray(),
              spareIndex,
              spare,
              model->djRegion(1),
              upperTheta,
              acceptablePivot,
              model->currentDualTolerance(),
              numberRemaining,
              zeroTolerance);
            model->spareDoubleArray_[0] = upperTheta;
            spareArray->setNumElements(numberRemaining);
#if 0
                              columnArray->setNumElements(numberNonZero);
			      spareArray->setPackedMode(true);
			      columnArray->setPackedMode(true);
			      spareArray->checkClean();
			      columnArray->checkClean();
			      printf("col %d spare %d upper %g\n",
				     columnArray->getNumElements(),
				     spareArray->getNumElements(),upperTheta);
#endif
            // signal partially done
            model->spareIntArray_[0] = -2;
          } else {
            numberNonZero = thisMatrix->gutsOfTransposeTimesUnscaled(pi,
              columnArray->getIndices(),
              columnArray->denseVector(),
              model->statusArray(),
              zeroTolerance);
          }
        } else {
          numberNonZero = thisMatrix->gutsOfTransposeTimesUnscaled(pi,
            columnArray->getIndices(),
            columnArray->denseVector(),
            zeroTolerance);
        }
        columnArray->setNumElements(numberNonZero);
        //xA++;
      } else {
        if ((model->moreSpecialOptions() & 8) != 0 && model->algorithm() < 0) {
          columnCopy_->transposeTimes(model, pi, columnArray,
            model->rowArray(3), rowArray);
          model->spareIntArray_[0] = -2;
        } else {
          columnCopy_->transposeTimes(model, pi, columnArray);
        }
        numberNonZero = columnArray->getNumElements();
        //xB++;
      }
    } else {
#ifdef CLP_INVESTIGATE
      if (model->clpScaledMatrix())
        printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
      // scaled
      // modify pi so can collapse to one loop
      if (scalar == -1.0) {
        //yC += numberInRowArray;
        for (i = 0; i < numberInRowArray; i++) {
          int iRow = whichRow[i];
          pi[iRow] = -piOld[i] * rowScale[iRow];
        }
      } else {
        for (i = 0; i < numberInRowArray; i++) {
          int iRow = whichRow[i];
          pi[iRow] = scalar * piOld[i] * rowScale[iRow];
        }
      }
      const double *columnScale = model->columnScale();
      if (!columnCopy_) {
        if ((model->specialOptions(), 131072) != 0)
          numberNonZero = gutsOfTransposeTimesScaled(pi, columnScale,
            columnArray->getIndices(),
            columnArray->denseVector(),
            model->statusArray(),
            zeroTolerance);
        else
          numberNonZero = gutsOfTransposeTimesScaled(pi, columnScale,
            columnArray->getIndices(),
            columnArray->denseVector(),
            zeroTolerance);
        columnArray->setNumElements(numberNonZero);
        //xC++;
      } else {
        if ((model->moreSpecialOptions() & 8) != 0 && model->algorithm() < 0) {
          columnCopy_->transposeTimes(model, pi, columnArray,
            model->rowArray(3), rowArray);
          model->spareIntArray_[0] = -2;
        } else {
          columnCopy_->transposeTimes(model, pi, columnArray);
        }
        numberNonZero = columnArray->getNumElements();
        //xD++;
      }
    }
    // zero out
    int numberRows = model->numberRows();
    if (numberInRowArray * 4 < numberRows) {
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = 0.0;
      }
    } else {
      CoinZeroN(pi, numberRows);
    }
    //int kA=xA+xB;
    //int kC=xC+xD;
    //if ((kA+kC)%10000==0)
    //printf("AA %d %d %g, CC %d %d %g\n",
    //     xA,xB,kA ? yA/(double)(kA): 0.0,xC,xD,kC ? yC/(double) (kC) :0.0);
  } else {
    if (!rowScale) {
      if (scalar == -1.0) {
        double value = 0.0;
        CoinBigIndex j;
        CoinBigIndex end = columnStart[1];
        for (j = columnStart[0]; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j];
        }
        for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
          CoinBigIndex start = end;
          end = columnStart[iColumn + 2];
          if (fabs(value) > zeroTolerance) {
            array[iColumn] = -value;
            index[numberNonZero++] = iColumn;
          }
          value = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j];
          }
        }
        if (fabs(value) > zeroTolerance) {
          array[iColumn] = -value;
          index[numberNonZero++] = iColumn;
        }
      } else {
        double value = 0.0;
        CoinBigIndex j;
        CoinBigIndex end = columnStart[1];
        for (j = columnStart[0]; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j];
        }
        for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
          value *= scalar;
          CoinBigIndex start = end;
          end = columnStart[iColumn + 2];
          if (fabs(value) > zeroTolerance) {
            array[iColumn] = value;
            index[numberNonZero++] = iColumn;
          }
          value = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j];
          }
        }
        value *= scalar;
        if (fabs(value) > zeroTolerance) {
          array[iColumn] = value;
          index[numberNonZero++] = iColumn;
        }
      }
    } else {
#ifdef CLP_INVESTIGATE
      if (model->clpScaledMatrix())
        printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
      // scaled
      if (scalar == -1.0) {
        const double *COIN_RESTRICT columnScale = model->columnScale();
        double value = 0.0;
        double scale = columnScale[0];
        CoinBigIndex j;
        CoinBigIndex end = columnStart[1];
        for (j = columnStart[0]; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
        }
        for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
          value *= scale;
          CoinBigIndex start = end;
          end = columnStart[iColumn + 2];
          scale = columnScale[iColumn + 1];
          if (fabs(value) > zeroTolerance) {
            array[iColumn] = -value;
            index[numberNonZero++] = iColumn;
          }
          value = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
          }
        }
        value *= scale;
        if (fabs(value) > zeroTolerance) {
          array[iColumn] = -value;
          index[numberNonZero++] = iColumn;
        }
      } else {
        const double *COIN_RESTRICT columnScale = model->columnScale();
        double value = 0.0;
        double scale = columnScale[0] * scalar;
        CoinBigIndex j;
        CoinBigIndex end = columnStart[1];
        for (j = columnStart[0]; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
        }
        for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
          value *= scale;
          CoinBigIndex start = end;
          end = columnStart[iColumn + 2];
          scale = columnScale[iColumn + 1] * scalar;
          if (fabs(value) > zeroTolerance) {
            array[iColumn] = value;
            index[numberNonZero++] = iColumn;
          }
          value = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
          }
        }
        value *= scale;
        if (fabs(value) > zeroTolerance) {
          array[iColumn] = value;
          index[numberNonZero++] = iColumn;
        }
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
  if (packed)
    columnArray->setPackedMode(true);
}
/* Return <code>x * A in <code>z</code>.
	Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix::transposeTimesByRow(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = getVectorStarts();
  const double *COIN_RESTRICT element = getElements();
  const int *COIN_RESTRICT whichRow = rowArray->getIndices();
  bool packed = rowArray->packedMode();
  if (numberInRowArray > 2) {
    // do by rows
    // ** Row copy is already scaled
    int iRow;
    int i;
    int numberOriginal = 0;
    if (packed) {
      int *COIN_RESTRICT index = columnArray->getIndices();
      double *COIN_RESTRICT array = columnArray->denseVector();
#if 0
	       {
		 double  * COIN_RESTRICT array2 = y->denseVector();
		 int numberColumns = matrix_->getNumCols();
		 for (int i=0;i<numberColumns;i++) {
		   assert(!array[i]);
		   assert(!array2[i]);
		 }
	       }
#endif
#define COIN_SPARSE_MATRIX 2
      CoinBigIndex numberCovered = 0;
      int numberColumns = matrix_->getNumCols();
      bool sparse = true;
      for (int i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        numberCovered += rowStart[iRow + 1] - rowStart[iRow];
        // ? exact ratio
        if (numberCovered > numberColumns) {
          sparse = false;
          break;
        }
      }
      if (sparse) {
        assert(!y->getNumElements());
#if COIN_SPARSE_MATRIX != 2
        // and set up mark as char array
        char *COIN_RESTRICT marked = reinterpret_cast< char * >(index + columnArray->capacity());
        int *COIN_RESTRICT lookup = y->getIndices();
#ifndef NDEBUG
        //int numberColumns = matrix_->getNumCols();
        //for (int i=0;i<numberColumns;i++)
        //assert(!marked[i]);
#endif
        numberNonZero = gutsOfTransposeTimesByRowGE3a(rowArray, index, array,
          lookup, marked, zeroTolerance, scalar);
#else
        double *COIN_RESTRICT array2 = y->denseVector();
        numberNonZero = gutsOfTransposeTimesByRowGE3(rowArray, index, array,
          array2, zeroTolerance, scalar);
#endif
      } else {
        numberNonZero = gutsOfTransposeTimesByRowGEK(rowArray, index, array,
          numberColumns, zeroTolerance, scalar);
      }
      columnArray->setNumElements(numberNonZero);
    } else {
      double *COIN_RESTRICT markVector = y->denseVector();
      numberNonZero = 0;
      // and set up mark as char array
      char *COIN_RESTRICT marked = reinterpret_cast< char * >(markVector);
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        marked[iColumn] = 0;
      }

      for (i = 0; i < numberInRowArray; i++) {
        iRow = whichRow[i];
        double value = pi[iRow] * scalar;
        CoinBigIndex j;
        for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
          int iColumn = column[j];
          if (!marked[iColumn]) {
            marked[iColumn] = 1;
            index[numberNonZero++] = iColumn;
          }
          array[iColumn] += value * element[j];
        }
      }
      // get rid of tiny values and zero out marked
      numberOriginal = numberNonZero;
      numberNonZero = 0;
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        marked[iColumn] = 0;
        if (fabs(array[iColumn]) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
        } else {
          array[iColumn] = 0.0;
        }
      }
    }
  } else if (numberInRowArray == 2) {
    // do by rows when two rows
    int numberOriginal;
    int i;
    CoinBigIndex j;
    numberNonZero = 0;

    double value;
    if (packed) {
      gutsOfTransposeTimesByRowEQ2(rowArray, columnArray, y, zeroTolerance, scalar);
      numberNonZero = columnArray->getNumElements();
    } else {
      int iRow = whichRow[0];
      value = pi[iRow] * scalar;
      for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
        int iColumn = column[j];
        double value2 = value * element[j];
        index[numberNonZero++] = iColumn;
        array[iColumn] = value2;
      }
      iRow = whichRow[1];
      value = pi[iRow] * scalar;
      for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
        int iColumn = column[j];
        double value2 = value * element[j];
        // I am assuming no zeros in matrix
        if (array[iColumn])
          value2 += array[iColumn];
        else
          index[numberNonZero++] = iColumn;
        array[iColumn] = value2;
      }
      // get rid of tiny values and zero out marked
      numberOriginal = numberNonZero;
      numberNonZero = 0;
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        if (fabs(array[iColumn]) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
        } else {
          array[iColumn] = 0.0;
        }
      }
    }
  } else if (numberInRowArray == 1) {
    // Just one row
    int iRow = rowArray->getIndices()[0];
    numberNonZero = 0;
    CoinBigIndex j;
    if (packed) {
      gutsOfTransposeTimesByRowEQ1(rowArray, columnArray, zeroTolerance, scalar);
      numberNonZero = columnArray->getNumElements();
    } else {
      double value = pi[iRow] * scalar;
      for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
        int iColumn = column[j];
        double value2 = value * element[j];
        if (fabs(value2) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
          array[iColumn] = value2;
        }
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
}
// Meat of transposeTimes by column when not scaled
int ClpPackedMatrix::gutsOfTransposeTimesUnscaled(const double *COIN_RESTRICT pi,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT array,
  const double zeroTolerance) const
{
  int numberNonZero = 0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
#if 1 //ndef INTEL_MKL
  double value = 0.0;
  CoinBigIndex j;
  CoinBigIndex end = columnStart[1];
  for (j = columnStart[0]; j < end; j++) {
    int iRow = row[j];
    value += pi[iRow] * elementByColumn[j];
  }
  int iColumn;
  for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
    CoinBigIndex start = end;
    end = columnStart[iColumn + 2];
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = iColumn;
    }
    value = 0.0;
    for (j = start; j < end; j++) {
      int iRow = row[j];
      value += pi[iRow] * elementByColumn[j];
    }
  }
  if (fabs(value) > zeroTolerance) {
    array[numberNonZero] = value;
    index[numberNonZero++] = iColumn;
  }
#else
  char transA = 'N';
  //int numberRows = matrix_->getNumRows();
  mkl_cspblas_dcsrgemv(&transA, const_cast< int * >(&numberActiveColumns_),
    const_cast< double * >(elementByColumn),
    const_cast< int * >(columnStart),
    const_cast< int * >(row),
    const_cast< double * >(pi), array);
  int iColumn;
  for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
    double value = array[iColumn];
    if (value) {
      array[iColumn] = 0.0;
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = iColumn;
      }
    }
  }
#endif
  return numberNonZero;
}
// Meat of transposeTimes by column when scaled
int ClpPackedMatrix::gutsOfTransposeTimesScaled(const double *COIN_RESTRICT pi,
  const double *COIN_RESTRICT columnScale,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT array,
  const double zeroTolerance) const
{
  int numberNonZero = 0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  double value = 0.0;
  double scale = columnScale[0];
  CoinBigIndex j;
  CoinBigIndex end = columnStart[1];
  for (j = columnStart[0]; j < end; j++) {
    int iRow = row[j];
    value += pi[iRow] * elementByColumn[j];
  }
  int iColumn;
  for (iColumn = 0; iColumn < numberActiveColumns_ - 1; iColumn++) {
    value *= scale;
    CoinBigIndex start = end;
    scale = columnScale[iColumn + 1];
    end = columnStart[iColumn + 2];
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = iColumn;
    }
    value = 0.0;
    for (j = start; j < end; j++) {
      int iRow = row[j];
      value += pi[iRow] * elementByColumn[j];
    }
  }
  value *= scale;
  if (fabs(value) > zeroTolerance) {
    array[numberNonZero] = value;
    index[numberNonZero++] = iColumn;
  }
  return numberNonZero;
}
#if ABOCA_LITE
static void
transposeTimesUnscaledBit(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT array = info.infeas;
  int *COIN_RESTRICT index = info.which;
  const unsigned char *COIN_RESTRICT status = info.status;
  double zeroTolerance = info.tolerance;
  const double *COIN_RESTRICT pi = info.work;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  double value = 0.0;
  int jColumn = -1;
  int numberNonZero = 0;
  for (int iColumn = first; iColumn < last; iColumn++) {
    bool wanted = ((status[iColumn] & 3) != 1);
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = jColumn;
    }
    value = 0.0;
    if (wanted) {
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = columnStart[iColumn + 1];
      jColumn = iColumn;
      int n = static_cast< int >(end - start);
      bool odd = (n & 1) != 0;
      n = n >> 1;
      const int *COIN_RESTRICT rowThis = row + start;
      const double *COIN_RESTRICT elementThis = elementByColumn + start;
      for (; n; n--) {
        int iRow0 = *rowThis;
        int iRow1 = *(rowThis + 1);
        rowThis += 2;
        value += pi[iRow0] * (*elementThis);
        value += pi[iRow1] * (*(elementThis + 1));
        elementThis += 2;
      }
      if (odd) {
        int iRow = *rowThis;
        value += pi[iRow] * (*elementThis);
      }
    }
  }
  if (fabs(value) > zeroTolerance) {
    array[numberNonZero] = value;
    index[numberNonZero++] = jColumn;
  }
  info.numberAdded = numberNonZero;
}
#endif
// Meat of transposeTimes by column when not scaled
int ClpPackedMatrix::gutsOfTransposeTimesUnscaled(const double *COIN_RESTRICT pi,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT array,
  const unsigned char *COIN_RESTRICT status,
  const double zeroTolerance) const
{
  int numberNonZero = 0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
#if ABOCA_LITE
  int numberThreads = abcState();
  if (numberThreads) {
    clpTempInfo info[ABOCA_LITE];
    int chunk = (numberActiveColumns_ + numberThreads - 1) / numberThreads;
    int n = 0;
    for (int i = 0; i < numberThreads; i++) {
      info[i].which = index + n;
      info[i].infeas = array + n;
      info[i].work = const_cast< double * >(pi);
      info[i].status = const_cast< unsigned char * >(status);
      info[i].startColumn = n;
      info[i].element = elementByColumn;
      info[i].start = columnStart;
      info[i].row = row;
      info[i].numberToDo = CoinMin(chunk, numberActiveColumns_ - n);
      info[i].tolerance = zeroTolerance;
      n += chunk;
    }
    for (int i = 0; i < numberThreads; i++) {
      cilk_spawn transposeTimesUnscaledBit(info[i]);
    }
    cilk_sync;
    for (int i = 0; i < numberThreads; i++)
      numberNonZero += info[i].numberAdded;
    moveAndZero(info, 2, NULL);
  } else {
#endif
    double value = 0.0;
    int jColumn = -1;
    for (int iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      bool wanted = ((status[iColumn] & 3) != 1);
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = jColumn;
      }
      value = 0.0;
      if (wanted) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn + 1];
        jColumn = iColumn;
        int n = static_cast< int >(end - start);
        bool odd = (n & 1) != 0;
        n = n >> 1;
        const int *COIN_RESTRICT rowThis = row + start;
        const double *COIN_RESTRICT elementThis = elementByColumn + start;
        for (; n; n--) {
          int iRow0 = *rowThis;
          int iRow1 = *(rowThis + 1);
          rowThis += 2;
          value += pi[iRow0] * (*elementThis);
          value += pi[iRow1] * (*(elementThis + 1));
          elementThis += 2;
        }
        if (odd) {
          int iRow = *rowThis;
          value += pi[iRow] * (*elementThis);
        }
      }
    }
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = jColumn;
    }
#if ABOCA_LITE
  }
#endif
  return numberNonZero;
}
#if ABOCA_LITE
static void
transposeTimesUnscaledBit2(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT array = info.infeas;
  const double *COIN_RESTRICT reducedCost = info.reducedCost;
  int *COIN_RESTRICT index = info.which;
  double *COIN_RESTRICT spareArray = info.spare;
  int *COIN_RESTRICT spareIndex = info.index;
  const unsigned char *COIN_RESTRICT status = info.status;
  double zeroTolerance = info.tolerance;
  double dualTolerance = info.dualTolerance;
  double acceptablePivot = info.acceptablePivot;
  const double *COIN_RESTRICT pi = info.work;
  double tentativeTheta = 1.0e15;
  int numberRemaining = 0;
  double upperTheta = info.upperTheta;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  int numberNonZero = 0;
  double multiplier[] = { -1.0, 1.0 };
  double dualT = -dualTolerance;
  for (int iColumn = first; iColumn < last; iColumn++) {
    int wanted = (status[iColumn] & 3) - 1;
    if (wanted) {
      double value = 0.0;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = columnStart[iColumn + 1];
      int n = static_cast< int >(end - start);
      bool odd = (n & 1) != 0;
      n = n >> 1;
      const int *COIN_RESTRICT rowThis = row + start;
      const double *COIN_RESTRICT elementThis = elementByColumn + start;
      for (; n; n--) {
        int iRow0 = *rowThis;
        int iRow1 = *(rowThis + 1);
        rowThis += 2;
        value += pi[iRow0] * (*elementThis);
        value += pi[iRow1] * (*(elementThis + 1));
        elementThis += 2;
      }
      if (odd) {
        int iRow = *rowThis;
        value += pi[iRow] * (*elementThis);
      }
      if (fabs(value) > zeroTolerance) {
        double mult = multiplier[wanted - 1];
        double alpha = value * mult;
        array[numberNonZero] = value;
        index[numberNonZero++] = iColumn;
        if (alpha > 0.0) {
          double oldValue = reducedCost[iColumn] * mult;
          double value = oldValue - tentativeTheta * alpha;
          if (value < dualT) {
            value = oldValue - upperTheta * alpha;
            if (value < dualT && alpha >= acceptablePivot) {
              upperTheta = (oldValue - dualT) / alpha;
              //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
            }
            // add to list
            spareArray[numberRemaining] = alpha * mult;
            spareIndex[numberRemaining++] = iColumn;
          }
        }
      }
    }
  }
  info.numberAdded = numberNonZero;
  info.numberRemaining = numberRemaining;
  info.upperTheta = upperTheta;
}
#endif
/* Meat of transposeTimes by column when not scaled and skipping
   and doing part of dualColumn */
int ClpPackedMatrix::gutsOfTransposeTimesUnscaled(const double *COIN_RESTRICT pi,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT array,
  const unsigned char *COIN_RESTRICT status,
  int *COIN_RESTRICT spareIndex,
  double *COIN_RESTRICT spareArray,
  const double *COIN_RESTRICT reducedCost,
  double &upperThetaP,
  double acceptablePivot,
  double dualTolerance,
  int &numberRemainingP,
  const double zeroTolerance) const
{
  int numberRemaining = numberRemainingP;
  double upperTheta = upperThetaP;
  int numberNonZero = 0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
#if ABOCA_LITE
  int numberThreads = abcState();
  if (numberThreads) {
    clpTempInfo info[ABOCA_LITE];
    int chunk = (numberActiveColumns_ + numberThreads - 1) / numberThreads;
    int n = 0;
    for (int i = 0; i < numberThreads; i++) {
      info[i].upperTheta = upperThetaP;
      info[i].acceptablePivot = acceptablePivot;
      info[i].reducedCost = const_cast< double * >(reducedCost);
      info[i].index = spareIndex + n + numberRemaining;
      info[i].spare = spareArray + n + numberRemaining;
      info[i].which = index + n;
      info[i].infeas = array + n;
      info[i].work = const_cast< double * >(pi);
      info[i].status = const_cast< unsigned char * >(status);
      info[i].startColumn = n;
      info[i].element = elementByColumn;
      info[i].start = columnStart;
      info[i].row = row;
      info[i].numberToDo = CoinMin(chunk, numberActiveColumns_ - n);
      info[i].tolerance = zeroTolerance;
      info[i].dualTolerance = dualTolerance;
      n += chunk;
    }
    for (int i = 0; i < numberThreads; i++) {
      cilk_spawn transposeTimesUnscaledBit2(info[i]);
    }
    cilk_sync;
    for (int i = 0; i < numberThreads; i++) {
      numberNonZero += info[i].numberAdded;
      numberRemaining += info[i].numberRemaining;
      upperTheta = CoinMin(upperTheta, static_cast< double >(info[i].upperTheta));
    }
    moveAndZero(info, 1, NULL);
    moveAndZero(info, 2, NULL);
  } else {
#endif
    double tentativeTheta = 1.0e15;
    double multiplier[] = { -1.0, 1.0 };
    double dualT = -dualTolerance;
    for (int iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      int wanted = (status[iColumn] & 3) - 1;
      if (wanted) {
        double value = 0.0;
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn + 1];
        int n = static_cast< int >(end - start);
#if 1
        bool odd = (n & 1) != 0;
        n = n >> 1;
        const int *COIN_RESTRICT rowThis = row + start;
        const double *COIN_RESTRICT elementThis = elementByColumn + start;
        for (; n; n--) {
          int iRow0 = *rowThis;
          int iRow1 = *(rowThis + 1);
          rowThis += 2;
          value += pi[iRow0] * (*elementThis);
          value += pi[iRow1] * (*(elementThis + 1));
          elementThis += 2;
        }
        if (odd) {
          int iRow = *rowThis;
          value += pi[iRow] * (*elementThis);
        }
#else
      const int *COIN_RESTRICT rowThis = &row[end - 16];
      const double *COIN_RESTRICT elementThis = &elementByColumn[end - 16];
      bool odd = (n & 1) != 0;
      n = n >> 1;
      double value2 = 0.0;
      if (odd) {
        int iRow = row[start];
        value2 = pi[iRow] * elementByColumn[start];
      }
      switch (n) {
      default: {
        if (odd) {
          start++;
        }
        n -= 8;
        for (; n; n--) {
          int iRow0 = row[start];
          int iRow1 = row[start + 1];
          value += pi[iRow0] * elementByColumn[start];
          value2 += pi[iRow1] * elementByColumn[start + 1];
          start += 2;
        }
      case 8: {
        int iRow0 = rowThis[16 - 16];
        int iRow1 = rowThis[16 - 15];
        value += pi[iRow0] * elementThis[16 - 16];
        value2 += pi[iRow1] * elementThis[16 - 15];
      }
      case 7: {
        int iRow0 = rowThis[16 - 14];
        int iRow1 = rowThis[16 - 13];
        value += pi[iRow0] * elementThis[16 - 14];
        value2 += pi[iRow1] * elementThis[16 - 13];
      }
      case 6: {
        int iRow0 = rowThis[16 - 12];
        int iRow1 = rowThis[16 - 11];
        value += pi[iRow0] * elementThis[16 - 12];
        value2 += pi[iRow1] * elementThis[16 - 11];
      }
      case 5: {
        int iRow0 = rowThis[16 - 10];
        int iRow1 = rowThis[16 - 9];
        value += pi[iRow0] * elementThis[16 - 10];
        value2 += pi[iRow1] * elementThis[16 - 9];
      }
      case 4: {
        int iRow0 = rowThis[16 - 8];
        int iRow1 = rowThis[16 - 7];
        value += pi[iRow0] * elementThis[16 - 8];
        value2 += pi[iRow1] * elementThis[16 - 7];
      }
      case 3: {
        int iRow0 = rowThis[16 - 6];
        int iRow1 = rowThis[16 - 5];
        value += pi[iRow0] * elementThis[16 - 6];
        value2 += pi[iRow1] * elementThis[16 - 5];
      }
      case 2: {
        int iRow0 = rowThis[16 - 4];
        int iRow1 = rowThis[16 - 3];
        value += pi[iRow0] * elementThis[16 - 4];
        value2 += pi[iRow1] * elementThis[16 - 3];
      }
      case 1: {
        int iRow0 = rowThis[16 - 2];
        int iRow1 = rowThis[16 - 1];
        value += pi[iRow0] * elementThis[16 - 2];
        value2 += pi[iRow1] * elementThis[16 - 1];
      }
      case 0:;
      }
      }
      value += value2;
#endif
        if (fabs(value) > zeroTolerance) {
          double mult = multiplier[wanted - 1];
          double alpha = value * mult;
          array[numberNonZero] = value;
          index[numberNonZero++] = iColumn;
          if (alpha > 0.0) {
            double oldValue = reducedCost[iColumn] * mult;
            double value = oldValue - tentativeTheta * alpha;
            if (value < dualT) {
              value = oldValue - upperTheta * alpha;
              if (value < dualT && alpha >= acceptablePivot) {
                upperTheta = (oldValue - dualT) / alpha;
                //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
              }
              // add to list
              spareArray[numberRemaining] = alpha * mult;
              spareIndex[numberRemaining++] = iColumn;
            }
          }
        }
      }
    }
#if ABOCA_LITE
  }
#endif
  numberRemainingP = numberRemaining;
  upperThetaP = upperTheta;
  return numberNonZero;
}
// Meat of transposeTimes by column when scaled
int ClpPackedMatrix::gutsOfTransposeTimesScaled(const double *COIN_RESTRICT pi,
  const double *COIN_RESTRICT columnScale,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT array,
  const unsigned char *COIN_RESTRICT status, const double zeroTolerance) const
{
  int numberNonZero = 0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  double value = 0.0;
  int jColumn = -1;
  for (int iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
    bool wanted = ((status[iColumn] & 3) != 1);
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = jColumn;
    }
    value = 0.0;
    if (wanted) {
      double scale = columnScale[iColumn];
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = columnStart[iColumn + 1];
      jColumn = iColumn;
      for (CoinBigIndex j = start; j < end; j++) {
        int iRow = row[j];
        value += pi[iRow] * elementByColumn[j];
      }
      value *= scale;
    }
  }
  if (fabs(value) > zeroTolerance) {
    array[numberNonZero] = value;
    index[numberNonZero++] = jColumn;
  }
  return numberNonZero;
}
#if ABOCA_LITE
static void
packDownBit(clpTempInfo &info)
{
  double *COIN_RESTRICT output = info.infeas;
  int *COIN_RESTRICT index = info.which;
  double zeroTolerance = info.tolerance;
  int first = info.startColumn;
  int last = info.numberToDo;
  int numberNonZero = 0;
  for (int i = 0; i < last; i++) {
    double value = output[i];
    if (value) {
      output[i] = 0.0;
      if (fabs(value) > zeroTolerance) {
        output[numberNonZero] = value;
        index[numberNonZero++] = i + first;
      }
    }
  }
  info.numberAdded = numberNonZero;
}
#endif
// Meat of transposeTimes by row n > K if packed - returns number nonzero
int ClpPackedMatrix::gutsOfTransposeTimesByRowGEK(const CoinIndexedVector *COIN_RESTRICT piVector,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT output,
  int numberColumns,
  const double tolerance,
  const double scalar) const
{
  const double *COIN_RESTRICT pi = piVector->denseVector();
  int numberInRowArray = piVector->getNumElements();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT element = matrix_->getElements();
  const int *COIN_RESTRICT whichRow = piVector->getIndices();
  // ** Row copy is already scaled
  for (int i = 0; i < numberInRowArray; i++) {
    int iRow = whichRow[i];
    double value = pi[i] * scalar;
    CoinBigIndex start = rowStart[iRow];
    CoinBigIndex end = rowStart[iRow + 1];
    int n = static_cast< int >(end - start);
    const int *COIN_RESTRICT columnThis = column + start;
    const double *COIN_RESTRICT elementThis = element + start;

    // could do by twos
    for (; n; n--) {
      int iColumn = *columnThis;
      columnThis++;
      double elValue = *elementThis;
      elementThis++;
      elValue *= value;
      output[iColumn] += elValue;
    }
  }
  // get rid of tiny values and count
  int numberNonZero = 0;
#if ABOCA_LITE
  int numberThreads = abcState();
  if (numberThreads) {
    clpTempInfo info[ABOCA_LITE];
    int chunk = (numberColumns + numberThreads - 1) / numberThreads;
    int n = 0;
    for (int i = 0; i < numberThreads; i++) {
      info[i].which = index + n;
      info[i].infeas = output + n;
      info[i].startColumn = n;
      info[i].numberToDo = CoinMin(chunk, numberColumns - n);
      info[i].tolerance = tolerance;
      n += chunk;
    }
    for (int i = 0; i < numberThreads; i++) {
      cilk_spawn packDownBit(info[i]);
    }
    cilk_sync;
    for (int i = 0; i < numberThreads; i++) {
      numberNonZero += info[i].numberAdded;
    }
    moveAndZero(info, 2, NULL);
  } else {
#endif
    for (int i = 0; i < numberColumns; i++) {
      double value = output[i];
      if (value) {
        output[i] = 0.0;
        if (fabs(value) > tolerance) {
          output[numberNonZero] = value;
          index[numberNonZero++] = i;
        }
      }
    }
#if ABOCA_LITE
  }
#endif
#ifndef NDEBUG
  for (int i = numberNonZero; i < numberColumns; i++)
    assert(!output[i]);
#endif
  return numberNonZero;
}
// Meat of transposeTimes by row n == 2 if packed
void ClpPackedMatrix::gutsOfTransposeTimesByRowEQ2(const CoinIndexedVector *piVector, CoinIndexedVector *output,
  CoinIndexedVector *spareVector, const double tolerance, const double scalar) const
{
  double *COIN_RESTRICT pi = piVector->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = output->getIndices();
  double *COIN_RESTRICT array = output->denseVector();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT element = matrix_->getElements();
  const int *COIN_RESTRICT whichRow = piVector->getIndices();
  int iRow0 = whichRow[0];
  int iRow1 = whichRow[1];
  double pi0 = pi[0];
  double pi1 = pi[1];
  if (rowStart[iRow0 + 1] - rowStart[iRow0] > rowStart[iRow1 + 1] - rowStart[iRow1]) {
    // do one with fewer first
    iRow0 = iRow1;
    iRow1 = whichRow[0];
    pi0 = pi1;
    pi1 = pi[0];
  }
  // and set up mark as char array
  char *COIN_RESTRICT marked = reinterpret_cast< char * >(index + output->capacity());
  int *COIN_RESTRICT lookup = spareVector->getIndices();
  double value = pi0 * scalar;
  CoinBigIndex j;
  for (j = rowStart[iRow0]; j < rowStart[iRow0 + 1]; j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value * elValue;
    array[numberNonZero] = value2;
    marked[iColumn] = 1;
    lookup[iColumn] = numberNonZero;
    index[numberNonZero++] = iColumn;
  }
  value = pi1 * scalar;
  for (j = rowStart[iRow1]; j < rowStart[iRow1 + 1]; j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value * elValue;
    // I am assuming no zeros in matrix
    if (marked[iColumn]) {
      int iLookup = lookup[iColumn];
      array[iLookup] += value2;
    } else {
      if (fabs(value2) > tolerance) {
        array[numberNonZero] = value2;
        index[numberNonZero++] = iColumn;
      }
    }
  }
  // get rid of tiny values and zero out marked
  int i;
  int n = numberNonZero;
  numberNonZero = 0;
  for (i = 0; i < n; i++) {
    int iColumn = index[i];
    marked[iColumn] = 0;
    if (fabs(array[i]) > tolerance) {
      array[numberNonZero] = array[i];
      index[numberNonZero++] = index[i];
    }
  }
  memset(array + numberNonZero, 0, (n - numberNonZero) * sizeof(double));
  output->setNumElements(numberNonZero);
  spareVector->setNumElements(0);
}
// Meat of transposeTimes by row n == 1 if packed
void ClpPackedMatrix::gutsOfTransposeTimesByRowEQ1(const CoinIndexedVector *piVector, CoinIndexedVector *output,
  const double tolerance, const double scalar) const
{
  double *COIN_RESTRICT pi = piVector->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = output->getIndices();
  double *COIN_RESTRICT array = output->denseVector();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT element = matrix_->getElements();
  int iRow = piVector->getIndices()[0];
  numberNonZero = 0;
  CoinBigIndex j;
  double value = pi[0] * scalar;
  for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value * elValue;
    if (fabs(value2) > tolerance) {
      array[numberNonZero] = value2;
      index[numberNonZero++] = iColumn;
    }
  }
  output->setNumElements(numberNonZero);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y.
   Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix::subsetTransposeTimes(const ClpSimplex *model,
  const CoinIndexedVector *rowArray,
  const CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int jColumn;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  int numberToDo = y->getNumElements();
  const int *COIN_RESTRICT which = y->getIndices();
  assert(!rowArray->packedMode());
  columnArray->setPacked();
  ClpPackedMatrix *scaledMatrix = model->clpScaledMatrix();
  int flags = flags_;
  if (rowScale && scaledMatrix && !(scaledMatrix->flags() & 2)) {
    flags = 0;
    rowScale = NULL;
    // get matrix data pointers
    row = scaledMatrix->getIndices();
    columnStart = scaledMatrix->getVectorStarts();
    elementByColumn = scaledMatrix->getElements();
  }
  if (!(flags & 2) && numberToDo > 2) {
    // no gaps
    if (!rowScale) {
      int iColumn = which[0];
      double value = 0.0;
      CoinBigIndex j;
      int columnNext = which[1];
      CoinBigIndex startNext = columnStart[columnNext];
      //coin_prefetch_const(row+startNext);
      //coin_prefetch_const(elementByColumn+startNext);
      CoinBigIndex endNext = columnStart[columnNext + 1];
      for (j = columnStart[iColumn];
           j < columnStart[iColumn + 1]; j++) {
        int iRow = row[j];
        value += pi[iRow] * elementByColumn[j];
      }
      for (jColumn = 0; jColumn < numberToDo - 2; jColumn++) {
        CoinBigIndex start = startNext;
        CoinBigIndex end = endNext;
        columnNext = which[jColumn + 2];
        startNext = columnStart[columnNext];
        //coin_prefetch_const(row+startNext);
        //coin_prefetch_const(elementByColumn+startNext);
        endNext = columnStart[columnNext + 1];
        array[jColumn] = value;
        value = 0.0;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j];
        }
      }
      array[jColumn++] = value;
      value = 0.0;
      for (j = startNext; j < endNext; j++) {
        int iRow = row[j];
        value += pi[iRow] * elementByColumn[j];
      }
      array[jColumn] = value;
    } else {
#ifdef CLP_INVESTIGATE
      if (model->clpScaledMatrix())
        printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
      // scaled
      const double *columnScale = model->columnScale();
      int iColumn = which[0];
      double value = 0.0;
      double scale = columnScale[iColumn];
      CoinBigIndex j;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn + 1]; j++) {
        int iRow = row[j];
        value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
      }
      for (jColumn = 0; jColumn < numberToDo - 1; jColumn++) {
        int iColumn = which[jColumn + 1];
        value *= scale;
        scale = columnScale[iColumn];
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn + 1];
        array[jColumn] = value;
        value = 0.0;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
        }
      }
      value *= scale;
      array[jColumn] = value;
    }
  } else if (numberToDo) {
    // gaps
    if (!rowScale) {
      for (jColumn = 0; jColumn < numberToDo; jColumn++) {
        int iColumn = which[jColumn];
        double value = 0.0;
        CoinBigIndex j;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j];
        }
        array[jColumn] = value;
      }
    } else {
#ifdef CLP_INVESTIGATE
      if (model->clpScaledMatrix())
        printf("scaledMatrix_ at %d of ClpPackedMatrix - flags %d (%d) n %d\n",
          __LINE__, flags_, model->clpScaledMatrix()->flags(), numberToDo);
#endif
      // scaled
      const double *columnScale = model->columnScale();
      for (jColumn = 0; jColumn < numberToDo; jColumn++) {
        int iColumn = which[jColumn];
        double value = 0.0;
        CoinBigIndex j;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          value += pi[iRow] * elementByColumn[j] * rowScale[iRow];
        }
        value *= columnScale[iColumn];
        array[jColumn] = value;
      }
    }
  }
}
/* Returns true if can combine transposeTimes and subsetTransposeTimes
   and if it would be faster */
bool ClpPackedMatrix::canCombine(const ClpSimplex *model,
  const CoinIndexedVector *pi) const
{
  int numberInRowArray = pi->getNumElements();
  int numberRows = model->numberRows();
  bool packed = pi->packedMode();
  // factor should be smaller if doing both with two pi vectors
  double factor = 0.30;
  // We may not want to do by row if there may be cache problems
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberActiveColumns_ * sizeof(double) > 1000000) {
    if (numberRows * 10 < numberActiveColumns_)
      factor *= 0.333333333;
    else if (numberRows * 4 < numberActiveColumns_)
      factor *= 0.5;
    else if (numberRows * 2 < numberActiveColumns_)
      factor *= 0.66666666667;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  // if not packed then bias a bit more towards by column
  if (!packed)
    factor *= 0.9;
  // bias if columnCopy
  if (columnCopy_)
    factor *= 0.5;
  return ((numberInRowArray > factor * numberRows || !model->rowCopy()) && !(flags_ & 2));
}
#ifndef CLP_ALL_ONE_FILE
// These have to match ClpPrimalColumnSteepest version
#define reference(i) (((reference[i >> 5] >> (i & 31)) & 1) != 0)
#endif
#if ABOCA_LITE
static void
transposeTimes2UnscaledBit(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT array = info.infeas;
  double *COIN_RESTRICT weights = const_cast< double * >(info.lower);
  double *COIN_RESTRICT reducedCost = info.reducedCost;
  double *COIN_RESTRICT infeas = const_cast< double * >(info.upper);
  int *COIN_RESTRICT index = info.which;
  unsigned int *COIN_RESTRICT reference = reinterpret_cast< unsigned int * >(info.index);
  const unsigned char *COIN_RESTRICT status = info.status;
  double zeroTolerance = info.tolerance;
  double dualTolerance = info.dualTolerance;
  const double *COIN_RESTRICT pi = info.work;
  const double *COIN_RESTRICT piWeight = info.cost;
  double scaleFactor = info.primalRatio;
  double devex = info.changeObj;
  double referenceIn = info.theta;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  int killDjs = info.numberInfeasibilities;
  int numberNonZero = 0;
  double mmult[2] = { 1.0, -1.0 };
  for (int iColumn = first; iColumn < last; iColumn++) {
    int iStat = status[iColumn] & 7;
    if (iStat != 1 && iStat != 5) {
      double value = 0.0;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = columnStart[iColumn + 1];
      int n = static_cast< int >(end - start);
      const int *COIN_RESTRICT rowThis = row + start;
      const double *COIN_RESTRICT elementThis = elementByColumn + start;
      for (int i = 0; i < n; i++) {
        int iRow = rowThis[i];
        value -= pi[iRow] * elementThis[i];
      }
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (int i = 0; i < n; i++) {
          int iRow = rowThis[i];
          modification += piWeight[iRow] * elementThis[i];
        }
        double thisWeight = weights[iColumn];
        double oldWeight = thisWeight;
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        debug3(iColumn, thisWeight, pivotSquared, devex, pivot, modification, oldWeight);
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          value = reducedCost[iColumn] - value;
          reducedCost[iColumn] = value;
          if (iStat > 1 && iStat < 4) {
            value *= mmult[iStat - 2];
            if (value <= dualTolerance)
              value = 0.0;
          } else {
            // free or super basic
            if (fabs(value) > FREE_ACCEPT * dualTolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
            } else {
              value = 0.0;
            }
          }
          if (value) {
            value *= value;
            // store square in list
            if (infeas[iColumn]) {
              infeas[iColumn] = value; // already there
            } else {
              array[numberNonZero] = value;
              index[numberNonZero++] = iColumn;
            }
          } else {
            array[numberNonZero] = 0.0;
            index[numberNonZero++] = iColumn;
          }
        }
      }
    }
  }
  info.numberAdded = numberNonZero;
}
static void
transposeTimes2ScaledBit(clpTempInfo &info)
{
  const CoinBigIndex *COIN_RESTRICT columnStart = info.start;
  const int *COIN_RESTRICT row = info.row;
  const double *COIN_RESTRICT elementByColumn = info.element;
  double *COIN_RESTRICT array = info.infeas;
  double *COIN_RESTRICT weights = const_cast< double * >(info.lower);
  int *COIN_RESTRICT index = info.which;
  unsigned int *COIN_RESTRICT reference = reinterpret_cast< unsigned int * >(info.index);
  const unsigned char *COIN_RESTRICT status = info.status;
  double zeroTolerance = info.tolerance;
  double dualTolerance = info.dualTolerance;
  double *COIN_RESTRICT reducedCost = info.reducedCost;
  double *COIN_RESTRICT infeas = const_cast< double * >(info.upper);
  const double *COIN_RESTRICT pi = info.work;
  const double *COIN_RESTRICT piWeight = info.cost;
  const double *COIN_RESTRICT columnScale = info.solution;
  double scaleFactor = info.primalRatio;
  double devex = info.changeObj;
  double referenceIn = info.theta;
  int first = info.startColumn;
  int last = first + info.numberToDo;
  int killDjs = info.numberInfeasibilities;
  int numberNonZero = 0;
  double mmult[2] = { 1.0, -1.0 };
  for (int iColumn = first; iColumn < last; iColumn++) {
    int iStat = status[iColumn] & 7;
    if (iStat != 1 && iStat != 5) {
      double value = 0.0;
      double scale = columnScale[iColumn];
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = columnStart[iColumn + 1];
      int n = static_cast< int >(end - start);
      bool odd = (n & 1) != 0;
      n &= ~1;
      const int *COIN_RESTRICT rowThis = row + start;
      const double *COIN_RESTRICT elementThis = elementByColumn + start;
      for (int i = 0; i < n; i += 2) {
        int iRow0 = rowThis[i];
        int iRow1 = rowThis[i + 1];
        value -= pi[iRow0] * elementThis[i];
        value -= pi[iRow1] * elementThis[i + 1];
      }
      if (odd) {
        int iRow = rowThis[n];
        value -= pi[iRow] * elementThis[n];
      }
      value *= scale;
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (int i = 0; i < n; i += 2) {
          int iRow0 = rowThis[i];
          int iRow1 = rowThis[i + 1];
          modification += piWeight[iRow0] * elementThis[i];
          modification += piWeight[iRow1] * elementThis[i + 1];
        }
        if (odd) {
          int iRow = rowThis[n];
          modification += piWeight[iRow] * elementThis[n];
        }
        modification *= scale;
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          value = reducedCost[iColumn] - value;
          reducedCost[iColumn] = value;
          if (iStat > 1 && iStat < 4) {
            value *= mmult[iStat - 2];
            if (value <= dualTolerance)
              value = 0.0;
          } else {
            // free or super basic
            if (fabs(value) > FREE_ACCEPT * dualTolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
            } else {
              value = 0.0;
            }
          }
          if (value) {
            value *= value;
            // store square in list
            if (infeas[iColumn]) {
              infeas[iColumn] = value; // already there
            } else {
              array[numberNonZero] = value;
              index[numberNonZero++] = iColumn;
            }
          } else {
            array[numberNonZero] = 0.0;
            index[numberNonZero++] = iColumn;
          }
        }
      }
    }
  }
  info.numberAdded = numberNonZero;
}
#endif
/* Updates two arrays for steepest and does devex weights 
   Returns nonzero if updates reduced cost and infeas -
   new infeas in dj1 */
// Updates two arrays for steepest
int ClpPackedMatrix::transposeTimes2(const ClpSimplex *model,
  const CoinIndexedVector *pi1, CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2,
  CoinIndexedVector *spare,
  double *COIN_RESTRICT infeas,
  double *COIN_RESTRICT reducedCost,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  int returnCode = 0;
  // put row of tableau in dj1
  double *COIN_RESTRICT pi = pi1->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = dj1->getIndices();
  double *COIN_RESTRICT array = dj1->denseVector();
  int numberInRowArray = pi1->getNumElements();
  double zeroTolerance = model->zeroTolerance();
  double dualTolerance = model->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = CoinMin(1.0e-2, model->largestDualError());
  // allow tolerance at least slightly bigger than standard
  dualTolerance = dualTolerance + error;
  bool packed = pi1->packedMode();
  // do by column
  int iColumn;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  assert(!spare->getNumElements());
  assert(numberActiveColumns_ > 0);
  double *COIN_RESTRICT piWeight = pi2->denseVector();
  assert(!pi2->packedMode());
  bool killDjs = (scaleFactor == 0.0);
  if (!scaleFactor)
    scaleFactor = 1.0;
  if (packed) {
    // need to expand pi into y
    assert(spare->capacity() >= model->numberRows());
    double *COIN_RESTRICT piOld = pi;
    pi = spare->denseVector();
    const int *COIN_RESTRICT whichRow = pi1->getIndices();
    int i;
    ClpPackedMatrix *scaledMatrix = model->clpScaledMatrix();
    if (rowScale && scaledMatrix) {
      rowScale = NULL;
      // get matrix data pointers
      row = scaledMatrix->getIndices();
      columnStart = scaledMatrix->getVectorStarts();
      elementByColumn = scaledMatrix->getElements();
    }
    if (!rowScale) {
      // modify pi so can collapse to one loop
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = piOld[i];
      }
      if (!columnCopy_ || killDjs) {
        if (infeas)
          returnCode = 1;
#if ABOCA_LITE
        int numberThreads = abcState();
        if (numberThreads) {
          clpTempInfo info[ABOCA_LITE];
          int chunk = (numberActiveColumns_ + numberThreads - 1) / numberThreads;
          int n = 0;
          for (int i = 0; i < numberThreads; i++) {
            info[i].theta = referenceIn;
            info[i].changeObj = devex;
            info[i].primalRatio = scaleFactor;
            info[i].lower = weights;
            info[i].reducedCost = reducedCost;
            info[i].upper = infeas;
            info[i].which = index + n;
            info[i].infeas = array + n;
            info[i].index = reinterpret_cast< int * >(reference);
            info[i].work = const_cast< double * >(pi);
            info[i].cost = const_cast< double * >(piWeight);
            info[i].status = const_cast< unsigned char * >(model->statusArray());
            info[i].startColumn = n;
            info[i].element = elementByColumn;
            info[i].start = columnStart;
            info[i].row = row;
            info[i].numberToDo = CoinMin(chunk, numberActiveColumns_ - n);
            info[i].tolerance = zeroTolerance;
            info[i].dualTolerance = dualTolerance;
            info[i].numberInfeasibilities = killDjs ? 1 : 0;
            n += chunk;
          }
          for (int i = 0; i < numberThreads; i++) {
            cilk_spawn transposeTimes2UnscaledBit(info[i]);
          }
          cilk_sync;
          for (int i = 0; i < numberThreads; i++) {
            numberNonZero += info[i].numberAdded;
          }
          moveAndZero(info, 2, NULL);
        } else {
#endif
          CoinBigIndex j;
          CoinBigIndex end = columnStart[0];
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            CoinBigIndex start = end;
            end = columnStart[iColumn + 1];
            ClpSimplex::Status status = model->getStatus(iColumn);
            if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
              continue;
            double value = 0.0;
            for (j = start; j < end; j++) {
              int iRow = row[j];
              value -= pi[iRow] * elementByColumn[j];
            }
            if (fabs(value) > zeroTolerance) {
              // and do other array
              double modification = 0.0;
              for (j = start; j < end; j++) {
                int iRow = row[j];
                modification += piWeight[iRow] * elementByColumn[j];
              }
              double thisWeight = weights[iColumn];
              double oldWeight = thisWeight;
              double pivot = value * scaleFactor;
              double pivotSquared = pivot * pivot;
              thisWeight += pivotSquared * devex + pivot * modification;
              debug3(iColumn, thisWeight, pivotSquared, devex, pivot, modification, oldWeight);
              if (thisWeight < DEVEX_TRY_NORM) {
                if (referenceIn < 0.0) {
                  // steepest
                  thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
                } else {
                  // exact
                  thisWeight = referenceIn * pivotSquared;
                  if (reference(iColumn))
                    thisWeight += 1.0;
                  thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
                }
              }
              weights[iColumn] = thisWeight;
              if (!killDjs) {
                value = reducedCost[iColumn] - value;
                reducedCost[iColumn] = value;
                // simplify status
                ClpSimplex::Status status = model->getStatus(iColumn);

                switch (status) {

                case ClpSimplex::basic:
                case ClpSimplex::isFixed:
                  break;
                case ClpSimplex::isFree:
                case ClpSimplex::superBasic:
                  if (fabs(value) > FREE_ACCEPT * dualTolerance) {
                    // we are going to bias towards free (but only if reasonable)
                    value *= FREE_BIAS;
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                  break;
                case ClpSimplex::atUpperBound:
                  if (value > dualTolerance) {
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                  break;
                case ClpSimplex::atLowerBound:
                  if (value < -dualTolerance) {
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                }
              }
            }
          }
#if ABOCA_LITE
        }
#endif
      } else {
        // use special column copy
        // reset back
        assert(!killDjs);
        //if (killDjs)
        //   scaleFactor = 0.0;
        if (infeas)
          returnCode = 1;
        columnCopy_->transposeTimes2(model, pi, dj1, piWeight,
          infeas, reducedCost,
          referenceIn, devex,
          reference, weights, scaleFactor);
        numberNonZero = dj1->getNumElements();
      }
    } else {
      // scaled
      // modify pi so can collapse to one loop
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = piOld[i] * rowScale[iRow];
      }
      // can also scale piWeight as not used again
      int numberWeight = pi2->getNumElements();
      const int *indexWeight = pi2->getIndices();
      for (i = 0; i < numberWeight; i++) {
        int iRow = indexWeight[i];
        piWeight[iRow] *= rowScale[iRow];
      }
      if (!columnCopy_ || killDjs) {
        if (infeas)
          returnCode = 1;
        const double *COIN_RESTRICT columnScale = model->columnScale();
#if ABOCA_LITE
        int numberThreads = abcState();
        if (numberThreads) {
          clpTempInfo info[ABOCA_LITE];
          int chunk = (numberActiveColumns_ + numberThreads - 1) / numberThreads;
          int n = 0;
          for (int i = 0; i < numberThreads; i++) {
            info[i].theta = referenceIn;
            info[i].changeObj = devex;
            info[i].primalRatio = scaleFactor;
            info[i].lower = weights;
            info[i].reducedCost = reducedCost;
            info[i].upper = infeas;
            info[i].solution = const_cast< double * >(columnScale);
            info[i].which = index + n;
            info[i].infeas = array + n;
            info[i].index = reinterpret_cast< int * >(reference);
            info[i].work = const_cast< double * >(pi);
            info[i].cost = const_cast< double * >(piWeight);
            info[i].status = const_cast< unsigned char * >(model->statusArray());
            info[i].startColumn = n;
            info[i].element = elementByColumn;
            info[i].start = columnStart;
            info[i].row = row;
            info[i].numberToDo = CoinMin(chunk, numberActiveColumns_ - n);
            info[i].tolerance = zeroTolerance;
            info[i].dualTolerance = dualTolerance;
            info[i].numberInfeasibilities = killDjs ? 1 : 0;
            n += chunk;
          }
          for (int i = 0; i < numberThreads; i++) {
            cilk_spawn transposeTimes2ScaledBit(info[i]);
          }
          cilk_sync;
          for (int i = 0; i < numberThreads; i++) {
            numberNonZero += info[i].numberAdded;
          }
          moveAndZero(info, 2, NULL);
          if (infeas) {
            returnCode = 1;
            dj1->setNumElements(numberNonZero);
          }
        } else {
#endif
          CoinBigIndex j;
          CoinBigIndex end = columnStart[0];
          for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
            CoinBigIndex start = end;
            end = columnStart[iColumn + 1];
            ClpSimplex::Status status = model->getStatus(iColumn);
            if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
              continue;
            double scale = columnScale[iColumn];
            double value = 0.0;
            for (j = start; j < end; j++) {
              int iRow = row[j];
              value -= pi[iRow] * elementByColumn[j];
            }
            value *= scale;
            if (fabs(value) > zeroTolerance) {
              double modification = 0.0;
              for (j = start; j < end; j++) {
                int iRow = row[j];
                modification += piWeight[iRow] * elementByColumn[j];
              }
              modification *= scale;
              double thisWeight = weights[iColumn];
              double pivot = value * scaleFactor;
              double pivotSquared = pivot * pivot;
              thisWeight += pivotSquared * devex + pivot * modification;
              if (thisWeight < DEVEX_TRY_NORM) {
                if (referenceIn < 0.0) {
                  // steepest
                  thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
                } else {
                  // exact
                  thisWeight = referenceIn * pivotSquared;
                  if (reference(iColumn))
                    thisWeight += 1.0;
                  thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
                }
              }
              weights[iColumn] = thisWeight;
              if (!killDjs) {
                value = reducedCost[iColumn] - value;
                reducedCost[iColumn] = value;
                // simplify status
                ClpSimplex::Status status = model->getStatus(iColumn);

                switch (status) {

                case ClpSimplex::basic:
                case ClpSimplex::isFixed:
                  break;
                case ClpSimplex::isFree:
                case ClpSimplex::superBasic:
                  if (fabs(value) > FREE_ACCEPT * dualTolerance) {
                    // we are going to bias towards free (but only if reasonable)
                    value *= FREE_BIAS;
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                  break;
                case ClpSimplex::atUpperBound:
                  if (value > dualTolerance) {
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                  break;
                case ClpSimplex::atLowerBound:
                  if (value < -dualTolerance) {
                    value *= value;
                    // store square in list
                    if (infeas[iColumn]) {
                      infeas[iColumn] = value; // already there
                    } else {
                      array[numberNonZero] = value;
                      index[numberNonZero++] = iColumn;
                    }
                  } else {
                    array[numberNonZero] = 0.0;
                    index[numberNonZero++] = iColumn;
                  }
                }
              }
            }
          }
#if ABOCA_LITE
        }
#endif
        if (infeas && false) {
          double tolerance = model->currentDualTolerance();
          // we can't really trust infeasibilities if there is dual error
          // this coding has to mimic coding in checkDualSolution
          double error = CoinMin(1.0e-2, model->largestDualError());
          // allow tolerance at least slightly bigger than standard
          tolerance = tolerance + error;
          returnCode = 1;
          int number = 0;
          for (int j = 0; j < numberNonZero; j++) {
            int iSequence = index[j];
            double value = reducedCost[iSequence];
            double value2 = array[j];
            array[j] = 0.0;
            value -= value2;
            reducedCost[iSequence] = value;
            // simplify status
            ClpSimplex::Status status = model->getStatus(iSequence);

            switch (status) {

            case ClpSimplex::basic:
            case ClpSimplex::isFixed:
              break;
            case ClpSimplex::isFree:
            case ClpSimplex::superBasic:
              if (fabs(value) > FREE_ACCEPT * tolerance) {
                // we are going to bias towards free (but only if reasonable)
                value *= FREE_BIAS;
                value *= value;
                // store square in list
                if (infeas[iSequence]) {
                  infeas[iSequence] = value; // already there
                } else {
                  array[number] = value;
                  index[number++] = iSequence;
                }
              } else {
                array[number] = 0.0;
                index[number++] = iSequence;
              }
              break;
            case ClpSimplex::atUpperBound:
              if (value > tolerance) {
                value *= value;
                // store square in list
                if (infeas[iSequence]) {
                  infeas[iSequence] = value; // already there
                } else {
                  array[number] = value;
                  index[number++] = iSequence;
                }
              } else {
                array[number] = 0.0;
                index[number++] = iSequence;
              }
              break;
            case ClpSimplex::atLowerBound:
              if (value < -tolerance) {
                value *= value;
                // store square in list
                if (infeas[iSequence]) {
                  infeas[iSequence] = value; // already there
                } else {
                  array[number] = value;
                  index[number++] = iSequence;
                }
              } else {
                array[number] = 0.0;
                index[number++] = iSequence;
              }
            }
          }
          dj1->setNumElements(number);
        }
      } else {
        // use special column copy
        // reset back
        assert(!killDjs);
        //if (killDjs)
        //   scaleFactor = 0.0;
        if (infeas)
          returnCode = 1;
        columnCopy_->transposeTimes2(model, pi, dj1, piWeight,
          infeas, reducedCost,
          referenceIn, devex,
          reference, weights, scaleFactor);
        numberNonZero = dj1->getNumElements();
      }
    }
    // zero out
    int numberRows = model->numberRows();
    if (numberInRowArray * 4 < numberRows) {
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = 0.0;
      }
    } else {
      CoinZeroN(pi, numberRows);
    }
  } else {
    if (!rowScale) {
      CoinBigIndex j;
      CoinBigIndex end = columnStart[0];
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        CoinBigIndex start = end;
        end = columnStart[iColumn + 1];
        ClpSimplex::Status status = model->getStatus(iColumn);
        if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
          continue;
        double value = 0.0;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          value -= pi[iRow] * elementByColumn[j];
        }
        if (fabs(value) > zeroTolerance) {
          // and do other array
          double modification = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            modification += piWeight[iRow] * elementByColumn[j];
          }
          double thisWeight = weights[iColumn];
          double pivot = value * scaleFactor;
          double pivotSquared = pivot * pivot;
          thisWeight += pivotSquared * devex + pivot * modification;
          if (thisWeight < DEVEX_TRY_NORM) {
            if (referenceIn < 0.0) {
              // steepest
              thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
            } else {
              // exact
              thisWeight = referenceIn * pivotSquared;
              if (reference(iColumn))
                thisWeight += 1.0;
              thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
            }
          }
          weights[iColumn] = thisWeight;
          if (!killDjs) {
            array[iColumn] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      }
    } else {
#ifdef CLP_INVESTIGATE
      if (model->clpScaledMatrix())
        printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
      // scaled
      // can also scale piWeight as not used again
      int numberWeight = pi2->getNumElements();
      const int *COIN_RESTRICT indexWeight = pi2->getIndices();
      for (int i = 0; i < numberWeight; i++) {
        int iRow = indexWeight[i];
        piWeight[iRow] *= rowScale[iRow];
      }
      const double *COIN_RESTRICT columnScale = model->columnScale();
      CoinBigIndex j;
      CoinBigIndex end = columnStart[0];
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        CoinBigIndex start = end;
        end = columnStart[iColumn + 1];
        ClpSimplex::Status status = model->getStatus(iColumn);
        if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
          continue;
        double scale = columnScale[iColumn];
        double value = 0.0;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          value -= pi[iRow] * elementByColumn[j] * rowScale[iRow];
        }
        value *= scale;
        if (fabs(value) > zeroTolerance) {
          double modification = 0.0;
          for (j = start; j < end; j++) {
            int iRow = row[j];
            modification += piWeight[iRow] * elementByColumn[j];
          }
          modification *= scale;
          double thisWeight = weights[iColumn];
          double pivot = value * scaleFactor;
          double pivotSquared = pivot * pivot;
          thisWeight += pivotSquared * devex + pivot * modification;
          if (thisWeight < DEVEX_TRY_NORM) {
            if (referenceIn < 0.0) {
              // steepest
              thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
            } else {
              // exact
              thisWeight = referenceIn * pivotSquared;
              if (reference(iColumn))
                thisWeight += 1.0;
              thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
            }
          }
          weights[iColumn] = thisWeight;
          if (!killDjs) {
            array[iColumn] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      }
    }
  }
  dj1->setNumElements(numberNonZero);
  spare->setNumElements(0);
  if (packed)
    dj1->setPackedMode(true);
  return returnCode;
}
// Updates second array for steepest and does devex weights
void ClpPackedMatrix::subsetTimes2(const ClpSimplex *model,
  CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2, CoinIndexedVector *,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  int number = dj1->getNumElements();
  const int *COIN_RESTRICT index = dj1->getIndices();
  double *COIN_RESTRICT array = dj1->denseVector();
  assert(dj1->packedMode());

  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  double *COIN_RESTRICT piWeight = pi2->denseVector();
  bool killDjs = (scaleFactor == 0.0);
  if (!scaleFactor)
    scaleFactor = 1.0;
  if (!rowScale) {
    for (int k = 0; k < number; k++) {
      int iColumn = index[k];
      double pivot = array[k] * scaleFactor;
      if (killDjs)
        array[k] = 0.0;
      // and do other array
      double modification = 0.0;
      for (CoinBigIndex j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        modification += piWeight[iRow] * elementByColumn[j];
      }
      double thisWeight = weights[iColumn];
      double pivotSquared = pivot * pivot;
      thisWeight += pivotSquared * devex + pivot * modification;
      if (thisWeight < DEVEX_TRY_NORM) {
        if (referenceIn < 0.0) {
          // steepest
          thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn * pivotSquared;
          if (reference(iColumn))
            thisWeight += 1.0;
          thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
        }
      }
      weights[iColumn] = thisWeight;
    }
  } else {
#ifdef CLP_INVESTIGATE
    if (model->clpScaledMatrix())
      printf("scaledMatrix_ at %d of ClpPackedMatrix\n", __LINE__);
#endif
    // scaled
    const double *columnScale = model->columnScale();
    for (int k = 0; k < number; k++) {
      int iColumn = index[k];
      double pivot = array[k] * scaleFactor;
      double scale = columnScale[iColumn];
      if (killDjs)
        array[k] = 0.0;
      // and do other array
      double modification = 0.0;
      for (CoinBigIndex j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        modification += piWeight[iRow] * elementByColumn[j] * rowScale[iRow];
      }
      double thisWeight = weights[iColumn];
      modification *= scale;
      double pivotSquared = pivot * pivot;
      thisWeight += pivotSquared * devex + pivot * modification;
      if (thisWeight < DEVEX_TRY_NORM) {
        if (referenceIn < 0.0) {
          // steepest
          thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn * pivotSquared;
          if (reference(iColumn))
            thisWeight += 1.0;
          thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
        }
      }
      weights[iColumn] = thisWeight;
    }
  }
}
/// returns number of elements in column part of basis,
int ClpPackedMatrix::countBasis(const int *whichColumn,
  int &numberColumnBasic)
{
  const int *columnLength = matrix_->getVectorLengths();
  int i;
  CoinBigIndex numberElements = 0;
  // just count - can be over so ignore zero problem
  for (i = 0; i < numberColumnBasic; i++) {
    int iColumn = whichColumn[i];
    numberElements += columnLength[iColumn];
  }
#if COIN_BIG_INDEX
  if (numberElements > COIN_INT_MAX) {
    printf("Factorization too large\n");
    abort();
  }
#endif
  return static_cast< int >(numberElements);
}
void ClpPackedMatrix::fillBasis(ClpSimplex *model,
  const int *COIN_RESTRICT whichColumn,
  int &numberColumnBasic,
  int *COIN_RESTRICT indexRowU,
  int *COIN_RESTRICT start,
  int *COIN_RESTRICT rowCount,
  int *COIN_RESTRICT columnCount,
  CoinFactorizationDouble *COIN_RESTRICT elementU)
{
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  int i;
  CoinBigIndex numberElements = start[0];
  // fill
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  ClpPackedMatrix *scaledMatrix = model->clpScaledMatrix();
  if (scaledMatrix && true) {
    columnLength = scaledMatrix->matrix_->getVectorLengths();
    columnStart = scaledMatrix->matrix_->getVectorStarts();
    rowScale = NULL;
    row = scaledMatrix->matrix_->getIndices();
    elementByColumn = scaledMatrix->matrix_->getElements();
  }
  if ((flags_ & 1) == 0) {
    if (!rowScale) {
      // no scaling
      for (i = 0; i < numberColumnBasic; i++) {
        int iColumn = whichColumn[i];
        int length = columnLength[iColumn];
        CoinBigIndex startThis = columnStart[iColumn];
        columnCount[i] = length;
        CoinBigIndex endThis = startThis + length;
        for (CoinBigIndex j = startThis; j < endThis; j++) {
          int iRow = row[j];
          indexRowU[numberElements] = iRow;
          rowCount[iRow]++;
          assert(elementByColumn[j]);
          elementU[numberElements++] = elementByColumn[j];
        }
        start[i + 1] = static_cast< int >(numberElements);
      }
    } else {
      // scaling
      const double *COIN_RESTRICT columnScale = model->columnScale();
      for (i = 0; i < numberColumnBasic; i++) {
        int iColumn = whichColumn[i];
        double scale = columnScale[iColumn];
        int length = columnLength[iColumn];
        CoinBigIndex startThis = columnStart[iColumn];
        columnCount[i] = length;
        CoinBigIndex endThis = startThis + length;
        for (CoinBigIndex j = startThis; j < endThis; j++) {
          int iRow = row[j];
          indexRowU[numberElements] = iRow;
          rowCount[iRow]++;
          assert(elementByColumn[j]);
          elementU[numberElements++] = elementByColumn[j] * scale * rowScale[iRow];
        }
        start[i + 1] = static_cast< int >(numberElements);
      }
    }
  } else {
    // there are zero elements so need to look more closely
    if (!rowScale) {
      // no scaling
      for (i = 0; i < numberColumnBasic; i++) {
        int iColumn = whichColumn[i];
        CoinBigIndex j;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          double value = elementByColumn[j];
          if (value) {
            int iRow = row[j];
            indexRowU[numberElements] = iRow;
            rowCount[iRow]++;
            elementU[numberElements++] = value;
          }
        }
        start[i + 1] = static_cast< int >(numberElements);
        columnCount[i] = static_cast< int >(numberElements) - start[i];
      }
    } else {
      // scaling
      const double *columnScale = model->columnScale();
      for (i = 0; i < numberColumnBasic; i++) {
        int iColumn = whichColumn[i];
        CoinBigIndex j;
        double scale = columnScale[iColumn];
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[i]; j++) {
          double value = elementByColumn[j];
          if (value) {
            int iRow = row[j];
            indexRowU[numberElements] = iRow;
            rowCount[iRow]++;
            elementU[numberElements++] = value * scale * rowScale[iRow];
          }
        }
        start[i + 1] = static_cast< int >(numberElements);
        columnCount[i] = static_cast< int >(numberElements) - start[i];
      }
    }
  }
#if COIN_BIG_INDEX
  if (numberElements > COIN_INT_MAX)
    abort();
#endif
}
#if 0
int
ClpPackedMatrix::scale2(ClpModel * model) const
{
     ClpSimplex * baseModel = NULL;
#ifndef NDEBUG
     //checkFlags();
#endif
     int numberRows = model->numberRows();
     int numberColumns = matrix_->getNumCols();
     model->setClpScaledMatrix(NULL); // get rid of any scaled matrix
     // If empty - return as sanityCheck will trap
     if (!numberRows || !numberColumns) {
          model->setRowScale(NULL);
          model->setColumnScale(NULL);
          return 1;
     }
     ClpMatrixBase * rowCopyBase = model->rowCopy();
     double * rowScale;
     double * columnScale;
     //assert (!model->rowScale());
     bool arraysExist;
     double * inverseRowScale = NULL;
     double * inverseColumnScale = NULL;
     if (!model->rowScale()) {
          rowScale = new double [numberRows*2];
          columnScale = new double [numberColumns*2];
          inverseRowScale = rowScale + numberRows;
          inverseColumnScale = columnScale + numberColumns;
          arraysExist = false;
     } else {
          rowScale = model->mutableRowScale();
          columnScale = model->mutableColumnScale();
          inverseRowScale = model->mutableInverseRowScale();
          inverseColumnScale = model->mutableInverseColumnScale();
          arraysExist = true;
     }
     assert (inverseRowScale == rowScale + numberRows);
     assert (inverseColumnScale == columnScale + numberColumns);
     // we are going to mark bits we are interested in
     char * usefulRow = new char [numberRows];
     char * usefulColumn = new char [numberColumns];
     double * rowLower = model->rowLower();
     double * rowUpper = model->rowUpper();
     double * columnLower = model->columnLower();
     double * columnUpper = model->columnUpper();
     int iColumn, iRow;
     //#define LEAVE_FIXED
     // mark free rows
     for (iRow = 0; iRow < numberRows; iRow++) {
#if 0 //ndef LEAVE_FIXED
          if (rowUpper[iRow] < 1.0e20 ||
                    rowLower[iRow] > -1.0e20)
               usefulRow[iRow] = 1;
          else
               usefulRow[iRow] = 0;
#else
          usefulRow[iRow] = 1;
#endif
     }
     // mark empty and fixed columns
     // also see if worth scaling
     assert (model->scalingFlag() <= 4);
     //  scale_stats[model->scalingFlag()]++;
     double largest = 0.0;
     double smallest = 1.0e50;
     // get matrix data pointers
     int * row = matrix_->getMutableIndices();
     const CoinBigIndex * columnStart = matrix_->getVectorStarts();
     int * columnLength = matrix_->getMutableVectorLengths();
     double * elementByColumn = matrix_->getMutableElements();
     bool deletedElements = false;
     for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          CoinBigIndex j;
          char useful = 0;
          bool deleteSome = false;
          int start = columnStart[iColumn];
          int end = start + columnLength[iColumn];
#ifndef LEAVE_FIXED
          if (columnUpper[iColumn] >
                    columnLower[iColumn] + 1.0e-12) {
#endif
               for (j = start; j < end; j++) {
                    iRow = row[j];
                    double value = fabs(elementByColumn[j]);
                    if (value > 1.0e-20) {
                         if(usefulRow[iRow]) {
                              useful = 1;
                              largest = CoinMax(largest, fabs(elementByColumn[j]));
                              smallest = CoinMin(smallest, fabs(elementByColumn[j]));
                         }
                    } else {
                         // small
                         deleteSome = true;
                    }
               }
#ifndef LEAVE_FIXED
          } else {
               // just check values
               for (j = start; j < end; j++) {
                    double value = fabs(elementByColumn[j]);
                    if (value <= 1.0e-20) {
                         // small
                         deleteSome = true;
                    }
               }
          }
#endif
          usefulColumn[iColumn] = useful;
          if (deleteSome) {
               deletedElements = true;
               CoinBigIndex put = start;
               for (j = start; j < end; j++) {
                    double value = elementByColumn[j];
                    if (fabs(value) > 1.0e-20) {
                         row[put] = row[j];
                         elementByColumn[put++] = value;
                    }
               }
               columnLength[iColumn] = put - start;
          }
     }
     model->messageHandler()->message(CLP_PACKEDSCALE_INITIAL, *model->messagesPointer())
               << smallest << largest
               << CoinMessageEol;
     if (smallest >= 0.5 && largest <= 2.0 && !deletedElements) {
          // don't bother scaling
          model->messageHandler()->message(CLP_PACKEDSCALE_FORGET, *model->messagesPointer())
                    << CoinMessageEol;
          if (!arraysExist) {
               delete [] rowScale;
               delete [] columnScale;
          } else {
               model->setRowScale(NULL);
               model->setColumnScale(NULL);
          }
          delete [] usefulRow;
          delete [] usefulColumn;
          return 1;
     } else {
#ifdef CLP_INVESTIGATE
          if (deletedElements)
               printf("DEL_ELS\n");
#endif
          if (!rowCopyBase) {
               // temporary copy
               rowCopyBase = reverseOrderedCopy();
          } else if (deletedElements) {
               rowCopyBase = reverseOrderedCopy();
               model->setNewRowCopy(rowCopyBase);
          }
#ifndef NDEBUG
          ClpPackedMatrix* rowCopy =
               dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
          // Make sure it is really a ClpPackedMatrix
          assert (rowCopy != NULL);
#else
          ClpPackedMatrix* rowCopy =
               static_cast< ClpPackedMatrix*>(rowCopyBase);
#endif

          const int * column = rowCopy->getIndices();
          const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
          const double * element = rowCopy->getElements();
          // need to scale
          if (largest > 1.0e13 * smallest) {
               // safer to have smaller zero tolerance
               double ratio = smallest / largest;
               ClpSimplex * simplex = static_cast<ClpSimplex *> (model);
               double newTolerance = CoinMax(ratio * 0.5, 1.0e-18);
               if (simplex->zeroTolerance() > newTolerance)
                    simplex->setZeroTolerance(newTolerance);
          }
          int scalingMethod = model->scalingFlag();
          if (scalingMethod == 4) {
               // As auto
               scalingMethod = 3;
          }
          // and see if there any empty rows
          for (iRow = 0; iRow < numberRows; iRow++) {
               if (usefulRow[iRow]) {
                    CoinBigIndex j;
                    int useful = 0;
                    for (j = rowStart[iRow]; j < rowStart[iRow+1]; j++) {
                         int iColumn = column[j];
                         if (usefulColumn[iColumn]) {
                              useful = 1;
                              break;
                         }
                    }
                    usefulRow[iRow] = static_cast<char>(useful);
               }
          }
          double savedOverallRatio = 0.0;
          double tolerance = 5.0 * model->primalTolerance();
          double overallLargest = -1.0e-20;
          double overallSmallest = 1.0e20;
          bool finished = false;
          // if scalingMethod 3 then may change
          bool extraDetails = (model->logLevel() > 2);
          while (!finished) {
               int numberPass = 3;
               overallLargest = -1.0e-20;
               overallSmallest = 1.0e20;
               if (!baseModel) {
                    ClpFillN ( rowScale, numberRows, 1.0);
                    ClpFillN ( columnScale, numberColumns, 1.0);
               } else {
                    // Copy scales and do quick scale on extra rows
                    // Then just one? pass
                    assert (numberColumns == baseModel->numberColumns());
                    int numberRows2 = baseModel->numberRows();
                    assert (numberRows >= numberRows2);
                    assert (baseModel->rowScale());
                    CoinMemcpyN(baseModel->rowScale(), numberRows2, rowScale);
                    CoinMemcpyN(baseModel->columnScale(), numberColumns, columnScale);
                    if (numberRows > numberRows2) {
                         numberPass = 1;
                         // do some scaling
                         if (scalingMethod == 1 || scalingMethod == 3) {
                              // Maximum in each row
                              for (iRow = numberRows2; iRow < numberRows; iRow++) {
                                   if (usefulRow[iRow]) {
                                        CoinBigIndex j;
                                        largest = 1.0e-10;
                                        for (j = rowStart[iRow]; j < rowStart[iRow+1]; j++) {
                                             int iColumn = column[j];
                                             if (usefulColumn[iColumn]) {
                                                  double value = fabs(element[j] * columnScale[iColumn]);
                                                  largest = CoinMax(largest, value);
                                                  assert (largest < 1.0e40);
                                             }
                                        }
                                        rowScale[iRow] = 1.0 / largest;
#ifdef COIN_DEVELOP
                                        if (extraDetails) {
                                             overallLargest = CoinMax(overallLargest, largest);
                                             overallSmallest = CoinMin(overallSmallest, largest);
                                        }
#endif
                                   }
                              }
                         } else {
                              overallLargest = 0.0;
                              overallSmallest = 1.0e50;
                              // Geometric mean on row scales
                              for (iRow = 0; iRow < numberRows; iRow++) {
                                   if (usefulRow[iRow]) {
                                        CoinBigIndex j;
                                        largest = 1.0e-20;
                                        smallest = 1.0e50;
                                        for (j = rowStart[iRow]; j < rowStart[iRow+1]; j++) {
                                             int iColumn = column[j];
                                             if (usefulColumn[iColumn]) {
                                                  double value = fabs(element[j]);
                                                  value *= columnScale[iColumn];
                                                  largest = CoinMax(largest, value);
                                                  smallest = CoinMin(smallest, value);
                                             }
                                        }
                                        if (iRow >= numberRows2) {
                                             rowScale[iRow] = 1.0 / sqrt(smallest * largest);
                                             //rowScale[iRow]=CoinMax(1.0e-10,CoinMin(1.0e10,rowScale[iRow]));
                                        }
#ifdef COIN_DEVELOP
                                        if (extraDetails) {
                                             overallLargest = CoinMax(largest * rowScale[iRow], overallLargest);
                                             overallSmallest = CoinMin(smallest * rowScale[iRow], overallSmallest);
                                        }
#endif
                                   }
                              }
                         }
                    } else {
                         // just use
                         numberPass = 0;
                    }
               }
               if (!baseModel && (scalingMethod == 1 || scalingMethod == 3)) {
                    // Maximum in each row
                    for (iRow = 0; iRow < numberRows; iRow++) {
                         if (usefulRow[iRow]) {
                              CoinBigIndex j;
                              largest = 1.0e-10;
                              for (j = rowStart[iRow]; j < rowStart[iRow+1]; j++) {
                                   int iColumn = column[j];
                                   if (usefulColumn[iColumn]) {
                                        double value = fabs(element[j]);
                                        largest = CoinMax(largest, value);
                                        assert (largest < 1.0e40);
                                   }
                              }
                              rowScale[iRow] = 1.0 / largest;
#ifdef COIN_DEVELOP
                              if (extraDetails) {
                                   overallLargest = CoinMax(overallLargest, largest);
                                   overallSmallest = CoinMin(overallSmallest, largest);
                              }
#endif
                         }
                    }
               } else {
#ifdef USE_OBJECTIVE
                    // This will be used to help get scale factors
                    double * objective = new double[numberColumns];
                    CoinMemcpyN(model->costRegion(1), numberColumns, objective);
                    double objScale = 1.0;
#endif
                    while (numberPass) {
                         overallLargest = 0.0;
                         overallSmallest = 1.0e50;
                         numberPass--;
                         // Geometric mean on row scales
                         for (iRow = 0; iRow < numberRows; iRow++) {
                              if (usefulRow[iRow]) {
                                   CoinBigIndex j;
                                   largest = 1.0e-20;
                                   smallest = 1.0e50;
                                   for (j = rowStart[iRow]; j < rowStart[iRow+1]; j++) {
                                        int iColumn = column[j];
                                        if (usefulColumn[iColumn]) {
                                             double value = fabs(element[j]);
                                             value *= columnScale[iColumn];
                                             largest = CoinMax(largest, value);
                                             smallest = CoinMin(smallest, value);
                                        }
                                   }

                                   rowScale[iRow] = 1.0 / sqrt(smallest * largest);
                                   //rowScale[iRow]=CoinMax(1.0e-10,CoinMin(1.0e10,rowScale[iRow]));
                                   if (extraDetails) {
                                        overallLargest = CoinMax(largest * rowScale[iRow], overallLargest);
                                        overallSmallest = CoinMin(smallest * rowScale[iRow], overallSmallest);
                                   }
                              }
                         }
#ifdef USE_OBJECTIVE
                         largest = 1.0e-20;
                         smallest = 1.0e50;
                         for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                              if (usefulColumn[iColumn]) {
                                   double value = fabs(objective[iColumn]);
                                   value *= columnScale[iColumn];
                                   largest = CoinMax(largest, value);
                                   smallest = CoinMin(smallest, value);
                              }
                         }
                         objScale = 1.0 / sqrt(smallest * largest);
#endif
                         model->messageHandler()->message(CLP_PACKEDSCALE_WHILE, *model->messagesPointer())
                                   << overallSmallest
                                   << overallLargest
                                   << CoinMessageEol;
                         // skip last column round
                         if (numberPass == 1)
                              break;
                         // Geometric mean on column scales
                         for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                              if (usefulColumn[iColumn]) {
                                   CoinBigIndex j;
                                   largest = 1.0e-20;
                                   smallest = 1.0e50;
                                   for (j = columnStart[iColumn];
                                             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                                        iRow = row[j];
                                        double value = fabs(elementByColumn[j]);
                                        if (usefulRow[iRow]) {
                                             value *= rowScale[iRow];
                                             largest = CoinMax(largest, value);
                                             smallest = CoinMin(smallest, value);
                                        }
                                   }
#ifdef USE_OBJECTIVE
                                   if (fabs(objective[iColumn]) > 1.0e-20) {
                                        double value = fabs(objective[iColumn]) * objScale;
                                        largest = CoinMax(largest, value);
                                        smallest = CoinMin(smallest, value);
                                   }
#endif
                                   columnScale[iColumn] = 1.0 / sqrt(smallest * largest);
                                   //columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
                              }
                         }
                    }
#ifdef USE_OBJECTIVE
                    delete [] objective;
                    printf("obj scale %g - use it if you want\n", objScale);
#endif
               }
               // If ranges will make horrid then scale
               for (iRow = 0; iRow < numberRows; iRow++) {
                    if (usefulRow[iRow]) {
                         double difference = rowUpper[iRow] - rowLower[iRow];
                         double scaledDifference = difference * rowScale[iRow];
                         if (scaledDifference > tolerance && scaledDifference < 1.0e-4) {
                              // make gap larger
                              rowScale[iRow] *= 1.0e-4 / scaledDifference;
                              rowScale[iRow] = CoinMax(1.0e-10, CoinMin(1.0e10, rowScale[iRow]));
                              //printf("Row %d difference %g scaled diff %g => %g\n",iRow,difference,
                              // scaledDifference,difference*rowScale[iRow]);
                         }
                    }
               }
               // final pass to scale columns so largest is reasonable
               // See what smallest will be if largest is 1.0
               overallSmallest = 1.0e50;
               for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (usefulColumn[iColumn]) {
                         CoinBigIndex j;
                         largest = 1.0e-20;
                         smallest = 1.0e50;
                         for (j = columnStart[iColumn];
                                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                              iRow = row[j];
                              if(elementByColumn[j] && usefulRow[iRow]) {
                                   double value = fabs(elementByColumn[j] * rowScale[iRow]);
                                   largest = CoinMax(largest, value);
                                   smallest = CoinMin(smallest, value);
                              }
                         }
                         if (overallSmallest * largest > smallest)
                              overallSmallest = smallest / largest;
                    }
               }
               if (scalingMethod == 1 || scalingMethod == 2) {
                    finished = true;
               } else if (savedOverallRatio == 0.0 && scalingMethod != 4) {
                    savedOverallRatio = overallSmallest;
                    scalingMethod = 4;
               } else {
                    assert (scalingMethod == 4);
                    if (overallSmallest > 2.0 * savedOverallRatio) {
                         finished = true; // geometric was better
                         if (model->scalingFlag() == 4) {
                              // If in Branch and bound change
                              if ((model->specialOptions() & 1024) != 0) {
                                   model->scaling(2);
                              }
                         }
                    } else {
                         scalingMethod = 1; // redo equilibrium
                         if (model->scalingFlag() == 4) {
                              // If in Branch and bound change
                              if ((model->specialOptions() & 1024) != 0) {
                                   model->scaling(1);
                              }
                         }
                    }
#if 0
                    if (extraDetails) {
                         if (finished)
                              printf("equilibrium ratio %g, geometric ratio %g , geo chosen\n",
                                     savedOverallRatio, overallSmallest);
                         else
                              printf("equilibrium ratio %g, geometric ratio %g , equi chosen\n",
                                     savedOverallRatio, overallSmallest);
                    }
#endif
               }
          }
          //#define RANDOMIZE
#ifdef RANDOMIZE
          // randomize by up to 10%
          for (iRow = 0; iRow < numberRows; iRow++) {
               double value = 0.5 - randomNumberGenerator_.randomDouble(); //between -0.5 to + 0.5
               rowScale[iRow] *= (1.0 + 0.1 * value);
          }
#endif
          overallLargest = 1.0;
          if (overallSmallest < 1.0e-1)
               overallLargest = 1.0 / sqrt(overallSmallest);
          overallLargest = CoinMin(100.0, overallLargest);
          overallSmallest = 1.0e50;
          //printf("scaling %d\n",model->scalingFlag());
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
               if (columnUpper[iColumn] >
                         columnLower[iColumn] + 1.0e-12) {
                    //if (usefulColumn[iColumn]) {
                    CoinBigIndex j;
                    largest = 1.0e-20;
                    smallest = 1.0e50;
                    for (j = columnStart[iColumn];
                              j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                         iRow = row[j];
                         if(elementByColumn[j] && usefulRow[iRow]) {
                              double value = fabs(elementByColumn[j] * rowScale[iRow]);
                              largest = CoinMax(largest, value);
                              smallest = CoinMin(smallest, value);
                         }
                    }
                    columnScale[iColumn] = overallLargest / largest;
                    //columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
#ifdef RANDOMIZE
                    double value = 0.5 - randomNumberGenerator_.randomDouble(); //between -0.5 to + 0.5
                    columnScale[iColumn] *= (1.0 + 0.1 * value);
#endif
                    double difference = columnUpper[iColumn] - columnLower[iColumn];
                    if (difference < 1.0e-5 * columnScale[iColumn]) {
                         // make gap larger
                         columnScale[iColumn] = difference / 1.0e-5;
                         //printf("Column %d difference %g scaled diff %g => %g\n",iColumn,difference,
                         // scaledDifference,difference*columnScale[iColumn]);
                    }
                    double value = smallest * columnScale[iColumn];
                    if (overallSmallest > value)
                         overallSmallest = value;
                    //overallSmallest = CoinMin(overallSmallest,smallest*columnScale[iColumn]);
               }
          }
          model->messageHandler()->message(CLP_PACKEDSCALE_FINAL, *model->messagesPointer())
                    << overallSmallest
                    << overallLargest
                    << CoinMessageEol;
          if (overallSmallest < 1.0e-13) {
               // Change factorization zero tolerance
               double newTolerance = CoinMax(1.0e-15 * (overallSmallest / 1.0e-13),
                                             1.0e-18);
               ClpSimplex * simplex = static_cast<ClpSimplex *> (model);
               if (simplex->factorization()->zeroTolerance() > newTolerance)
                    simplex->factorization()->zeroTolerance(newTolerance);
               newTolerance = CoinMax(overallSmallest * 0.5, 1.0e-18);
               simplex->setZeroTolerance(newTolerance);
          }
          delete [] usefulRow;
          delete [] usefulColumn;
#ifndef SLIM_CLP
          // If quadratic then make symmetric
          ClpObjective * obj = model->objectiveAsObject();
#ifndef NO_RTTI
          ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(obj));
#else
          ClpQuadraticObjective * quadraticObj = NULL;
          if (obj->type() == 2)
               quadraticObj = (static_cast< ClpQuadraticObjective*>(obj));
#endif
          if (quadraticObj) {
               if (!rowCopyBase) {
                    // temporary copy
                    rowCopyBase = reverseOrderedCopy();
               }
#ifndef NDEBUG
               ClpPackedMatrix* rowCopy =
                    dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
               // Make sure it is really a ClpPackedMatrix
               assert (rowCopy != NULL);
#else
               ClpPackedMatrix* rowCopy =
                    static_cast< ClpPackedMatrix*>(rowCopyBase);
#endif
               const int * column = rowCopy->getIndices();
               const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
               CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
               int numberXColumns = quadratic->getNumCols();
               if (numberXColumns < numberColumns) {
                    // we assume symmetric
                    int numberQuadraticColumns = 0;
                    int i;
                    //const int * columnQuadratic = quadratic->getIndices();
                    //const int * columnQuadraticStart = quadratic->getVectorStarts();
                    const int * columnQuadraticLength = quadratic->getVectorLengths();
                    for (i = 0; i < numberXColumns; i++) {
                         int length = columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
                         length = 1;
#endif
                         if (length)
                              numberQuadraticColumns++;
                    }
                    int numberXRows = numberRows - numberQuadraticColumns;
                    numberQuadraticColumns = 0;
                    for (i = 0; i < numberXColumns; i++) {
                         int length = columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
                         length = 1;
#endif
                         if (length) {
                              rowScale[numberQuadraticColumns+numberXRows] = columnScale[i];
                              numberQuadraticColumns++;
                         }
                    }
                    int numberQuadraticRows = 0;
                    for (i = 0; i < numberXRows; i++) {
                         // See if any in row quadratic
                         CoinBigIndex j;
                         int numberQ = 0;
                         for (j = rowStart[i]; j < rowStart[i+1]; j++) {
                              int iColumn = column[j];
                              if (columnQuadraticLength[iColumn])
                                   numberQ++;
                         }
#ifndef CORRECT_ROW_COUNTS
                         numberQ = 1;
#endif
                         if (numberQ) {
                              columnScale[numberQuadraticRows+numberXColumns] = rowScale[i];
                              numberQuadraticRows++;
                         }
                    }
                    // and make sure Sj okay
                    for (iColumn = numberQuadraticRows + numberXColumns; iColumn < numberColumns; iColumn++) {
                         CoinBigIndex j = columnStart[iColumn];
                         assert(columnLength[iColumn] == 1);
                         int iRow = row[j];
                         double value = fabs(elementByColumn[j] * rowScale[iRow]);
                         columnScale[iColumn] = 1.0 / value;
                    }
               }
          }
#endif
          // make copy (could do faster by using previous values)
          // could just do partial
          for (iRow = 0; iRow < numberRows; iRow++)
               inverseRowScale[iRow] = 1.0 / rowScale[iRow] ;
          for (iColumn = 0; iColumn < numberColumns; iColumn++)
               inverseColumnScale[iColumn] = 1.0 / columnScale[iColumn] ;
          if (!arraysExist) {
               model->setRowScale(rowScale);
               model->setColumnScale(columnScale);
          }
          if (model->rowCopy()) {
               // need to replace row by row
               ClpPackedMatrix* rowCopy =
                    static_cast< ClpPackedMatrix*>(model->rowCopy());
               double * element = rowCopy->getMutableElements();
               const int * column = rowCopy->getIndices();
               const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
               // scale row copy
               for (iRow = 0; iRow < numberRows; iRow++) {
                    CoinBigIndex j;
                    double scale = rowScale[iRow];
                    double * elementsInThisRow = element + rowStart[iRow];
                    const int * columnsInThisRow = column + rowStart[iRow];
                    int number = rowStart[iRow+1] - rowStart[iRow];
                    assert (number <= numberColumns);
                    for (j = 0; j < number; j++) {
                         int iColumn = columnsInThisRow[j];
                         elementsInThisRow[j] *= scale * columnScale[iColumn];
                    }
               }
               if ((model->specialOptions() & 262144) != 0) {
                    //if ((model->specialOptions()&(COIN_CBC_USING_CLP|16384))!=0) {
                    //if (model->inCbcBranchAndBound()&&false) {
                    // copy without gaps
                    CoinPackedMatrix * scaledMatrix = new CoinPackedMatrix(*matrix_, 0, 0);
                    ClpPackedMatrix * scaled = new ClpPackedMatrix(scaledMatrix);
                    model->setClpScaledMatrix(scaled);
                    // get matrix data pointers
                    const int * row = scaledMatrix->getIndices();
                    const CoinBigIndex * columnStart = scaledMatrix->getVectorStarts();
#ifndef NDEBUG
                    const int * columnLength = scaledMatrix->getVectorLengths();
#endif
                    double * elementByColumn = scaledMatrix->getMutableElements();
                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                         CoinBigIndex j;
                         double scale = columnScale[iColumn];
                         assert (columnStart[iColumn+1] == columnStart[iColumn] + columnLength[iColumn]);
                         for (j = columnStart[iColumn];
                                   j < columnStart[iColumn+1]; j++) {
                              int iRow = row[j];
                              elementByColumn[j] *= scale * rowScale[iRow];
                         }
                    }
               } else {
                    //printf("not in b&b\n");
               }
          } else {
               // no row copy
               delete rowCopyBase;
          }
          return 0;
     }
}
#endif
//#define SQRT_ARRAY
#ifdef SQRT_ARRAY
static void doSqrts(double *array, int n)
{
  for (int i = 0; i < n; i++)
    array[i] = 1.0 / sqrt(array[i]);
}
#endif
//static int scale_stats[5]={0,0,0,0,0};
// Creates scales for column copy (rowCopy in model may be modified)
int ClpPackedMatrix::scale(ClpModel *model, ClpSimplex *simplex) const
{
  //const ClpSimplex * baseModel=NULL;
  //return scale2(model);
#if 0
     ClpMatrixBase * rowClone = NULL;
     if (model->rowCopy())
          rowClone = model->rowCopy()->clone();
     assert (!model->rowScale());
     assert (!model->columnScale());
     int returnCode = scale2(model);
     if (returnCode)
          return returnCode;
#endif
#ifndef NDEBUG
  //checkFlags();
#endif
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  model->setClpScaledMatrix(NULL); // get rid of any scaled matrix
  // If empty - return as sanityCheck will trap
  if (!numberRows || !numberColumns) {
    model->setRowScale(NULL);
    model->setColumnScale(NULL);
    return 1;
  }
#if 0
     // start fake
     double * rowScale2 = CoinCopyOfArray(model->rowScale(), numberRows);
     double * columnScale2 = CoinCopyOfArray(model->columnScale(), numberColumns);
     model->setRowScale(NULL);
     model->setColumnScale(NULL);
     model->setNewRowCopy(rowClone);
#endif
  ClpMatrixBase *rowCopyBase = model->rowCopy();
  double *COIN_RESTRICT rowScale;
  double *COIN_RESTRICT columnScale;
  //assert (!model->rowScale());
  bool arraysExist;
  double *COIN_RESTRICT inverseRowScale = NULL;
  double *COIN_RESTRICT inverseColumnScale = NULL;
  if (!model->rowScale()) {
    rowScale = new double[numberRows * 2];
    columnScale = new double[numberColumns * 2];
    inverseRowScale = rowScale + numberRows;
    inverseColumnScale = columnScale + numberColumns;
    arraysExist = false;
  } else {
    rowScale = model->mutableRowScale();
    columnScale = model->mutableColumnScale();
    inverseRowScale = model->mutableInverseRowScale();
    inverseColumnScale = model->mutableInverseColumnScale();
    arraysExist = true;
  }
  assert(inverseRowScale == rowScale + numberRows);
  assert(inverseColumnScale == columnScale + numberColumns);
  // we are going to mark bits we are interested in
  char *COIN_RESTRICT usefulColumn = new char[numberColumns];
  double *COIN_RESTRICT rowLower = model->rowLower();
  double *COIN_RESTRICT rowUpper = model->rowUpper();
  double *COIN_RESTRICT columnLower = model->columnLower();
  double *COIN_RESTRICT columnUpper = model->columnUpper();
  int iColumn, iRow;
  //#define LEAVE_FIXED
  // mark empty and fixed columns
  // also see if worth scaling
  assert(model->scalingFlag() <= 5);
  //  scale_stats[model->scalingFlag()]++;
  double largest = 0.0;
  double smallest = 1.0e50;
  // get matrix data pointers
  int *COIN_RESTRICT row = matrix_->getMutableIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  int *COIN_RESTRICT columnLength = matrix_->getMutableVectorLengths();
  double *COIN_RESTRICT elementByColumn = matrix_->getMutableElements();
  CoinBigIndex deletedElements = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    char useful = 0;
    bool deleteSome = false;
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
#ifndef LEAVE_FIXED
    if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-12 || (simplex && simplex->getColumnStatus(iColumn) == ClpSimplex::basic)) {
#endif
      for (j = start; j < end; j++) {
        iRow = row[j];
        double value = fabs(elementByColumn[j]);
        if (value > 1.0e-20) {
          useful = 1;
          largest = CoinMax(largest, fabs(elementByColumn[j]));
          smallest = CoinMin(smallest, fabs(elementByColumn[j]));
        } else {
          // small
          deleteSome = true;
        }
      }
#ifndef LEAVE_FIXED
    } else {
      // just check values
      for (j = start; j < end; j++) {
        double value = fabs(elementByColumn[j]);
        if (value <= 1.0e-20) {
          // small
          deleteSome = true;
        }
      }
    }
#endif
    usefulColumn[iColumn] = useful;
    if (deleteSome) {
      CoinBigIndex put = start;
      for (j = start; j < end; j++) {
        double value = elementByColumn[j];
        if (fabs(value) > 1.0e-20) {
          row[put] = row[j];
          elementByColumn[put++] = value;
        }
      }
      deletedElements += end - put;
      columnLength[iColumn] = static_cast< int >(put - start);
    }
  }
  // don't scale integers if option set
  if ((model->specialOptions() & 4194304) != 0 && model->integerInformation()) {
    const char *COIN_RESTRICT integer = model->integerInformation();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (integer[iColumn])
        usefulColumn[iColumn] = 0;
    }
  }
  if (deletedElements) {
    matrix_->setNumElements(matrix_->getNumElements() - deletedElements);
    flags_ |= 0x02;
  }
  model->messageHandler()->message(CLP_PACKEDSCALE_INITIAL, *model->messagesPointer())
    << smallest << largest
    << CoinMessageEol;
  if (smallest * 1.0e12 < largest && simplex) {
    // increase tolerances
    simplex->setCurrentDualTolerance(CoinMax(simplex->currentDualTolerance(), 5.0e-7));
    simplex->setCurrentPrimalTolerance(CoinMax(simplex->currentPrimalTolerance(), 5.0e-7));
  }
  if (smallest >= 0.5 && largest <= 2.0 && !deletedElements) {
    // don't bother scaling
    model->messageHandler()->message(CLP_PACKEDSCALE_FORGET, *model->messagesPointer())
      << CoinMessageEol;
    if (!arraysExist) {
      delete[] rowScale;
      delete[] columnScale;
    } else {
      model->setRowScale(NULL);
      model->setColumnScale(NULL);
    }
    delete[] usefulColumn;
    return 1;
  } else {
#ifdef CLP_INVESTIGATE
    if (deletedElements)
      printf("DEL_ELS\n");
#endif
    if (!rowCopyBase) {
      // temporary copy
      rowCopyBase = reverseOrderedCopy();
    } else if (deletedElements) {
      rowCopyBase = reverseOrderedCopy();
      model->setNewRowCopy(rowCopyBase);
    }
#ifndef NDEBUG
    ClpPackedMatrix *rowCopy = dynamic_cast< ClpPackedMatrix * >(rowCopyBase);
    // Make sure it is really a ClpPackedMatrix
    assert(rowCopy != NULL);
#else
    ClpPackedMatrix *rowCopy = static_cast< ClpPackedMatrix * >(rowCopyBase);
#endif

    const int *COIN_RESTRICT column = rowCopy->getIndices();
    const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
    const double *COIN_RESTRICT element = rowCopy->getElements();
    // need to scale
    if (largest > 1.0e13 * smallest) {
      // safer to have smaller zero tolerance
      double ratio = smallest / largest;
      ClpSimplex *simplex = static_cast< ClpSimplex * >(model);
      double newTolerance = CoinMax(ratio * 0.5, 1.0e-18);
      if (simplex->zeroTolerance() > newTolerance)
        simplex->setZeroTolerance(newTolerance);
    }
    int scalingMethod = model->scalingFlag();
    if (scalingMethod == 4) {
      // As auto
      scalingMethod = 3;
    } else if (scalingMethod == 5) {
      // As geometric
      scalingMethod = 2;
    }
    double savedOverallRatio = 0.0;
    double tolerance = 5.0 * model->primalTolerance();
    double overallLargest = -1.0e-20;
    double overallSmallest = 1.0e20;
    bool finished = false;
    // if scalingMethod 3 then may change
    bool extraDetails = (model->logLevel() > 2);
#if 0
	  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
	    if (columnUpper[iColumn] >
		columnLower[iColumn] + 1.0e-12 && columnLength[iColumn]) 
	      assert(usefulColumn[iColumn]!=0);
	    else
	      assert(usefulColumn[iColumn]==0);
	  }
#endif
    while (!finished) {
      int numberPass = 3;
      overallLargest = -1.0e-20;
      overallSmallest = 1.0e20;
      ClpFillN(rowScale, numberRows, 1.0);
      ClpFillN(columnScale, numberColumns, 1.0);
      if (scalingMethod == 1 || scalingMethod == 3) {
        // Maximum in each row
        for (iRow = 0; iRow < numberRows; iRow++) {
          CoinBigIndex j;
          largest = 1.0e-10;
          for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
            int iColumn = column[j];
            if (usefulColumn[iColumn]) {
              double value = fabs(element[j]);
              largest = CoinMax(largest, value);
              assert(largest < 1.0e40);
            }
          }
          rowScale[iRow] = 1.0 / largest;
#ifdef COIN_DEVELOP
          if (extraDetails) {
            overallLargest = CoinMax(overallLargest, largest);
            overallSmallest = CoinMin(overallSmallest, largest);
          }
#endif
        }
      } else {
#ifdef USE_OBJECTIVE
        // This will be used to help get scale factors
        double *COIN_RESTRICT objective = new double[numberColumns];
        CoinMemcpyN(model->costRegion(1), numberColumns, objective);
        double objScale = 1.0;
#endif
        while (numberPass) {
          overallLargest = 0.0;
          overallSmallest = 1.0e50;
          numberPass--;
          // Geometric mean on row scales
          for (iRow = 0; iRow < numberRows; iRow++) {
            CoinBigIndex j;
            largest = 1.0e-50;
            smallest = 1.0e50;
            for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
              int iColumn = column[j];
              if (usefulColumn[iColumn]) {
                double value = fabs(element[j]);
                value *= columnScale[iColumn];
                largest = CoinMax(largest, value);
                smallest = CoinMin(smallest, value);
              }
            }

#ifdef SQRT_ARRAY
            rowScale[iRow] = smallest * largest;
#else
            rowScale[iRow] = 1.0 / sqrt(smallest * largest);
#endif
            //rowScale[iRow]=CoinMax(1.0e-10,CoinMin(1.0e10,rowScale[iRow]));
            if (extraDetails) {
              overallLargest = CoinMax(largest * rowScale[iRow], overallLargest);
              overallSmallest = CoinMin(smallest * rowScale[iRow], overallSmallest);
            }
          }
          if (model->scalingFlag() == 5)
            break; // just scale rows
#ifdef SQRT_ARRAY
          doSqrts(rowScale, numberRows);
#endif
#ifdef USE_OBJECTIVE
          largest = 1.0e-20;
          smallest = 1.0e50;
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (usefulColumn[iColumn]) {
              double value = fabs(objective[iColumn]);
              value *= columnScale[iColumn];
              largest = CoinMax(largest, value);
              smallest = CoinMin(smallest, value);
            }
          }
          objScale = 1.0 / sqrt(smallest * largest);
#endif
          model->messageHandler()->message(CLP_PACKEDSCALE_WHILE, *model->messagesPointer())
            << overallSmallest
            << overallLargest
            << CoinMessageEol;
          // skip last column round
          if (numberPass == 1)
            break;
          // Geometric mean on column scales
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (usefulColumn[iColumn]) {
              CoinBigIndex j;
              largest = 1.0e-50;
              smallest = 1.0e50;
              for (j = columnStart[iColumn];
                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                iRow = row[j];
                double value = fabs(elementByColumn[j]);
                value *= rowScale[iRow];
                largest = CoinMax(largest, value);
                smallest = CoinMin(smallest, value);
              }
#ifdef USE_OBJECTIVE
              if (fabs(objective[iColumn]) > 1.0e-20) {
                double value = fabs(objective[iColumn]) * objScale;
                largest = CoinMax(largest, value);
                smallest = CoinMin(smallest, value);
              }
#endif
#ifdef SQRT_ARRAY
              columnScale[iColumn] = smallest * largest;
#else
              columnScale[iColumn] = 1.0 / sqrt(smallest * largest);
#endif
              //columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
            }
          }
#ifdef SQRT_ARRAY
          doSqrts(columnScale, numberColumns);
#endif
        }
#ifdef USE_OBJECTIVE
        delete[] objective;
        printf("obj scale %g - use it if you want\n", objScale);
#endif
      }
      // If ranges will make horrid then scale
      for (iRow = 0; iRow < numberRows; iRow++) {
        double difference = rowUpper[iRow] - rowLower[iRow];
        double scaledDifference = difference * rowScale[iRow];
        if (scaledDifference > tolerance && scaledDifference < 1.0e-4) {
          // make gap larger
          rowScale[iRow] *= 1.0e-4 / scaledDifference;
          rowScale[iRow] = CoinMax(1.0e-10, CoinMin(1.0e10, rowScale[iRow]));
          //printf("Row %d difference %g scaled diff %g => %g\n",iRow,difference,
          // scaledDifference,difference*rowScale[iRow]);
        }
      }
      // final pass to scale columns so largest is reasonable
      // See what smallest will be if largest is 1.0
      if (model->scalingFlag() != 5) {
        overallSmallest = 1.0e50;
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (usefulColumn[iColumn]) {
            CoinBigIndex j;
            largest = 1.0e-20;
            smallest = 1.0e50;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              iRow = row[j];
              double value = fabs(elementByColumn[j] * rowScale[iRow]);
              largest = CoinMax(largest, value);
              smallest = CoinMin(smallest, value);
            }
            if (overallSmallest * largest > smallest)
              overallSmallest = smallest / largest;
          }
        }
      }
      if (scalingMethod == 1 || scalingMethod == 2) {
        finished = true;
      } else if (savedOverallRatio == 0.0 && scalingMethod != 4) {
        savedOverallRatio = overallSmallest;
        scalingMethod = 4;
      } else {
        assert(scalingMethod == 4);
        if (overallSmallest > 2.0 * savedOverallRatio) {
          finished = true; // geometric was better
          if (model->scalingFlag() == 4) {
            // If in Branch and bound change
            if ((model->specialOptions() & 1024) != 0) {
              model->scaling(2);
            }
          }
        } else {
          scalingMethod = 1; // redo equilibrium
          if (model->scalingFlag() == 4) {
            // If in Branch and bound change
            if ((model->specialOptions() & 1024) != 0) {
              model->scaling(1);
            }
          }
        }
#if 0
                    if (extraDetails) {
                         if (finished)
                              printf("equilibrium ratio %g, geometric ratio %g , geo chosen\n",
                                     savedOverallRatio, overallSmallest);
                         else
                              printf("equilibrium ratio %g, geometric ratio %g , equi chosen\n",
                                     savedOverallRatio, overallSmallest);
                    }
#endif
      }
    }
    //#define RANDOMIZE
#ifdef RANDOMIZE
    // randomize by up to 10%
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value = 0.5 - randomNumberGenerator_.randomDouble(); //between -0.5 to + 0.5
      rowScale[iRow] *= (1.0 + 0.1 * value);
    }
#endif
    overallLargest = 1.0;
    if (overallSmallest < 1.0e-1)
      overallLargest = 1.0 / sqrt(overallSmallest);
    overallLargest = CoinMin(100.0, overallLargest);
    overallSmallest = 1.0e50;
    char *usedRow = reinterpret_cast< char * >(inverseRowScale);
    memset(usedRow, 0, numberRows);
    //printf("scaling %d\n",model->scalingFlag());
    if (model->scalingFlag() != 5) {
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (columnUpper[iColumn] > columnLower[iColumn] + 1.0e-12) {
          //if (usefulColumn[iColumn]) {
          CoinBigIndex j;
          largest = 1.0e-20;
          smallest = 1.0e50;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            iRow = row[j];
            usedRow[iRow] = 1;
            double value = fabs(elementByColumn[j] * rowScale[iRow]);
            largest = CoinMax(largest, value);
            smallest = CoinMin(smallest, value);
          }
          columnScale[iColumn] = overallLargest / largest;
          //columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
#ifdef RANDOMIZE
          double value = 0.5 - randomNumberGenerator_.randomDouble(); //between -0.5 to + 0.5
          columnScale[iColumn] *= (1.0 + 0.1 * value);
#endif
          double difference = columnUpper[iColumn] - columnLower[iColumn];
          if (difference < 1.0e-5 * columnScale[iColumn]) {
            // make gap larger
            columnScale[iColumn] = difference / 1.0e-5;
            //printf("Column %d difference %g scaled diff %g => %g\n",iColumn,difference,
            // scaledDifference,difference*columnScale[iColumn]);
          }
          double value = smallest * columnScale[iColumn];
          if (overallSmallest > value)
            overallSmallest = value;
          //overallSmallest = CoinMin(overallSmallest,smallest*columnScale[iColumn]);
        } else {
          //assert(columnScale[iColumn] == 1.0);
          columnScale[iColumn] = 1.0;
        }
      }
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (!usedRow[iRow]) {
          rowScale[iRow] = 1.0;
        }
      }
    }
    model->messageHandler()->message(CLP_PACKEDSCALE_FINAL, *model->messagesPointer())
      << overallSmallest
      << overallLargest
      << CoinMessageEol;
#if 0
          {
               for (iRow = 0; iRow < numberRows; iRow++) {
                    assert (rowScale[iRow] == rowScale2[iRow]);
               }
               delete [] rowScale2;
               for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    assert (columnScale[iColumn] == columnScale2[iColumn]);
               }
               delete [] columnScale2;
          }
#endif
    if (overallSmallest < 1.0e-13) {
      // Change factorization zero tolerance
      double newTolerance = CoinMax(1.0e-15 * (overallSmallest / 1.0e-13),
        1.0e-18);
      ClpSimplex *simplex = static_cast< ClpSimplex * >(model);
      if (simplex->factorization()->zeroTolerance() > newTolerance)
        simplex->factorization()->zeroTolerance(newTolerance);
      newTolerance = CoinMax(overallSmallest * 0.5, 1.0e-18);
      simplex->setZeroTolerance(newTolerance);
    }
    delete[] usefulColumn;
#ifndef SLIM_CLP
    // If quadratic then make symmetric
    ClpObjective *obj = model->objectiveAsObject();
#ifndef NO_RTTI
    ClpQuadraticObjective *quadraticObj = (dynamic_cast< ClpQuadraticObjective * >(obj));
#else
    ClpQuadraticObjective *quadraticObj = NULL;
    if (obj->type() == 2)
      quadraticObj = (static_cast< ClpQuadraticObjective * >(obj));
#endif
    if (quadraticObj) {
      if (!rowCopyBase) {
        // temporary copy
        rowCopyBase = reverseOrderedCopy();
      }
#ifndef NDEBUG
      ClpPackedMatrix *rowCopy = dynamic_cast< ClpPackedMatrix * >(rowCopyBase);
      // Make sure it is really a ClpPackedMatrix
      assert(rowCopy != NULL);
#else
      ClpPackedMatrix *rowCopy = static_cast< ClpPackedMatrix * >(rowCopyBase);
#endif
      const int *COIN_RESTRICT column = rowCopy->getIndices();
      const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
      CoinPackedMatrix *quadratic = quadraticObj->quadraticObjective();
      int numberXColumns = quadratic->getNumCols();
      if (numberXColumns < numberColumns) {
        // we assume symmetric
        int numberQuadraticColumns = 0;
        int i;
        //const int * columnQuadratic = quadratic->getIndices();
        //const int * columnQuadraticStart = quadratic->getVectorStarts();
        const int *COIN_RESTRICT columnQuadraticLength = quadratic->getVectorLengths();
        for (i = 0; i < numberXColumns; i++) {
          int length = columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
          length = 1;
#endif
          if (length)
            numberQuadraticColumns++;
        }
        int numberXRows = numberRows - numberQuadraticColumns;
        numberQuadraticColumns = 0;
        for (i = 0; i < numberXColumns; i++) {
          int length = columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
          length = 1;
#endif
          if (length) {
            rowScale[numberQuadraticColumns + numberXRows] = columnScale[i];
            numberQuadraticColumns++;
          }
        }
        int numberQuadraticRows = 0;
        for (i = 0; i < numberXRows; i++) {
          // See if any in row quadratic
          CoinBigIndex j;
          int numberQ = 0;
          for (j = rowStart[i]; j < rowStart[i + 1]; j++) {
            int iColumn = column[j];
            if (columnQuadraticLength[iColumn])
              numberQ++;
          }
#ifndef CORRECT_ROW_COUNTS
          numberQ = 1;
#endif
          if (numberQ) {
            columnScale[numberQuadraticRows + numberXColumns] = rowScale[i];
            numberQuadraticRows++;
          }
        }
        // and make sure Sj okay
        for (iColumn = numberQuadraticRows + numberXColumns; iColumn < numberColumns; iColumn++) {
          CoinBigIndex j = columnStart[iColumn];
          assert(columnLength[iColumn] == 1);
          int iRow = row[j];
          double value = fabs(elementByColumn[j] * rowScale[iRow]);
          columnScale[iColumn] = 1.0 / value;
        }
      }
    }
#endif
    // make copy (could do faster by using previous values)
    // could just do partial
    for (iRow = 0; iRow < numberRows; iRow++)
      inverseRowScale[iRow] = 1.0 / rowScale[iRow];
    for (iColumn = 0; iColumn < numberColumns; iColumn++)
      inverseColumnScale[iColumn] = 1.0 / columnScale[iColumn];
    if (!arraysExist) {
      model->setRowScale(rowScale);
      model->setColumnScale(columnScale);
    }
    if (model->rowCopy()) {
      // need to replace row by row
      ClpPackedMatrix *rowCopy = static_cast< ClpPackedMatrix * >(model->rowCopy());
      double *COIN_RESTRICT element = rowCopy->getMutableElements();
      const int *COIN_RESTRICT column = rowCopy->getIndices();
      const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
      // scale row copy
      for (iRow = 0; iRow < numberRows; iRow++) {
        CoinBigIndex j;
        double scale = rowScale[iRow];
        double *COIN_RESTRICT elementsInThisRow = element + rowStart[iRow];
        const int *COIN_RESTRICT columnsInThisRow = column + rowStart[iRow];
        int number = static_cast< int >(rowStart[iRow + 1] - rowStart[iRow]);
        assert(number <= numberColumns);
        for (j = 0; j < number; j++) {
          int iColumn = columnsInThisRow[j];
          elementsInThisRow[j] *= scale * columnScale[iColumn];
        }
      }
      if ((model->specialOptions() & 262144) != 0) {
        //if ((model->specialOptions()&(COIN_CBC_USING_CLP|16384))!=0) {
        //if (model->inCbcBranchAndBound()&&false) {
        // copy without gaps
        CoinPackedMatrix *scaledMatrix = new CoinPackedMatrix(*matrix_, 0, 0);
        ClpPackedMatrix *scaled = new ClpPackedMatrix(scaledMatrix);
        model->setClpScaledMatrix(scaled);
        // get matrix data pointers
        const int *COIN_RESTRICT row = scaledMatrix->getIndices();
        const CoinBigIndex *COIN_RESTRICT columnStart = scaledMatrix->getVectorStarts();
#ifndef NDEBUG
        const int *COIN_RESTRICT columnLength = scaledMatrix->getVectorLengths();
#endif
        double *COIN_RESTRICT elementByColumn = scaledMatrix->getMutableElements();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          CoinBigIndex j;
          double scale = columnScale[iColumn];
          assert(columnStart[iColumn + 1] == columnStart[iColumn] + columnLength[iColumn]);
          for (j = columnStart[iColumn];
               j < columnStart[iColumn + 1]; j++) {
            int iRow = row[j];
            elementByColumn[j] *= scale * rowScale[iRow];
          }
        }
      } else {
        //printf("not in b&b\n");
      }
    } else {
      // no row copy
      delete rowCopyBase;
    }
    return 0;
  }
}
// Creates scaled column copy if scales exist
void ClpPackedMatrix::createScaledMatrix(ClpSimplex *model) const
{
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  model->setClpScaledMatrix(NULL); // get rid of any scaled matrix
  // If empty - return as sanityCheck will trap
  if (!numberRows || !numberColumns) {
    model->setRowScale(NULL);
    model->setColumnScale(NULL);
    return;
  }
  if (!model->rowScale())
    return;
  double *COIN_RESTRICT rowScale = model->mutableRowScale();
  double *COIN_RESTRICT columnScale = model->mutableColumnScale();
  // copy without gaps
  CoinPackedMatrix *scaledMatrix = new CoinPackedMatrix(*matrix_, 0, 0);
  ClpPackedMatrix *scaled = new ClpPackedMatrix(scaledMatrix);
  model->setClpScaledMatrix(scaled);
  // get matrix data pointers
  const int *COIN_RESTRICT row = scaledMatrix->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = scaledMatrix->getVectorStarts();
#ifndef NDEBUG
  const int *COIN_RESTRICT columnLength = scaledMatrix->getVectorLengths();
#endif
  double *COIN_RESTRICT elementByColumn = scaledMatrix->getMutableElements();
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    double scale = columnScale[iColumn];
    assert(columnStart[iColumn + 1] == columnStart[iColumn] + columnLength[iColumn]);
    for (j = columnStart[iColumn];
         j < columnStart[iColumn + 1]; j++) {
      int iRow = row[j];
      elementByColumn[j] *= scale * rowScale[iRow];
    }
  }
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
/* Unpacks a column into an CoinIndexedvector
 */
void ClpPackedMatrix::unpack(const ClpSimplex *model, CoinIndexedVector *rowArray,
  int iColumn) const
{
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      rowArray->quickAdd(row[i], elementByColumn[i]);
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn];
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      rowArray->quickAdd(iRow, elementByColumn[i] * scale * rowScale[iRow]);
    }
  }
}
/* Unpacks a column into a CoinIndexedVector
** in packed format
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void ClpPackedMatrix::unpackPacked(ClpSimplex *model,
  CoinIndexedVector *rowArray,
  int iColumn) const
{
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  int *COIN_RESTRICT index = rowArray->getIndices();
  double *COIN_RESTRICT array = rowArray->denseVector();
  int number = 0;
  if (!rowScale) {
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      double value = elementByColumn[i];
      if (value) {
        array[number] = value;
        index[number++] = iRow;
      }
    }
    rowArray->setNumElements(number);
    rowArray->setPackedMode(true);
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn];
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      double value = elementByColumn[i] * scale * rowScale[iRow];
      if (value) {
        array[number] = value;
        index[number++] = iRow;
      }
    }
    rowArray->setNumElements(number);
    rowArray->setPackedMode(true);
  }
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void ClpPackedMatrix::add(const ClpSimplex *model, CoinIndexedVector *rowArray,
  int iColumn, double multiplier) const
{
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      rowArray->quickAdd(iRow, multiplier * elementByColumn[i]);
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn] * multiplier;
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      rowArray->quickAdd(iRow, elementByColumn[i] * scale * rowScale[iRow]);
    }
  }
}
/* Adds multiple of a column into an array */
void ClpPackedMatrix::add(const ClpSimplex *model, double *COIN_RESTRICT array,
  int iColumn, double multiplier) const
{
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      array[iRow] += multiplier * elementByColumn[i];
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn] * multiplier;
    for (i = columnStart[iColumn];
         i < columnStart[iColumn] + columnLength[iColumn]; i++) {
      int iRow = row[i];
      array[iRow] += elementByColumn[i] * scale * rowScale[iRow];
    }
  }
}
/* Checks if all elements are in valid range.  Can just
   return true if you are not paranoid.  For Clp I will
   probably expect no zeros.  Code can modify matrix to get rid of
   small elements.
   check bits (can be turned off to save time) :
   1 - check if matrix has gaps
   2 - check if zero elements
   4 - check and compress duplicates
   8 - report on large and small
*/
bool ClpPackedMatrix::allElementsInRange(ClpModel *model,
  double smallest, double largest,
  int check)
{
  int iColumn;
  // make sure matrix correct size
  assert(matrix_->getNumRows() <= model->numberRows());
  if (model->clpScaledMatrix())
    assert(model->clpScaledMatrix()->getNumElements() == matrix_->getNumElements());
  assert(matrix_->getNumRows() <= model->numberRows());
  matrix_->setDimensions(model->numberRows(), model->numberColumns());
  CoinBigIndex numberLarge = 0;
  ;
  CoinBigIndex numberSmall = 0;
  ;
  CoinBigIndex numberDuplicate = 0;
  ;
  int firstBadColumn = -1;
  int firstBadRow = -1;
  double firstBadElement = 0.0;
  // get matrix data pointers
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = matrix_->getVectorStarts();
  const int *COIN_RESTRICT columnLength = matrix_->getVectorLengths();
  const double *COIN_RESTRICT elementByColumn = matrix_->getElements();
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  // Say no gaps
  flags_ &= ~2;
  if (type_ >= 10)
    return true; // gub
  if (check == 14 || check == 10) {
    if (matrix_->getNumElements() < columnStart[numberColumns]) {
      // pack down!
#if 0
               matrix_->removeGaps();
#else
      checkGaps();
#ifdef DO_CHECK_FLAGS
      checkFlags(0);
#endif
#endif
#ifdef COIN_DEVELOP
      //printf("flags set to 2\n");
#endif
    } else if (numberColumns) {
      assert(columnStart[numberColumns] == columnStart[numberColumns - 1] + columnLength[numberColumns - 1]);
    }
    return true;
  }
  assert(check == 15 || check == 11);
  if (check == 15) {
    CoinBigIndex *COIN_RESTRICT mark = new CoinBigIndex[numberRows];
    int i;
    for (i = 0; i < numberRows; i++)
      mark[i] = -1;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (end != columnStart[iColumn + 1])
        flags_ |= 2;
      for (j = start; j < end; j++) {
        double value = fabs(elementByColumn[j]);
        int iRow = row[j];
        if (iRow < 0 || iRow >= numberRows) {
#ifndef COIN_BIG_INDEX
          printf("Out of range %d %d %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#elif COIN_BIG_INDEX == 0
          printf("Out of range %d %d %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#elif COIN_BIG_INDEX == 1
          printf("Out of range %d %ld %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#else
          printf("Out of range %d %lld %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#endif
          return false;
        }
        if (mark[iRow] == -1) {
          mark[iRow] = j;
        } else {
          // duplicate
          numberDuplicate++;
        }
        //printf("%d %d %d %g\n",iColumn,j,row[j],elementByColumn[j]);
        if (!value)
          flags_ |= 1; // there are zero elements
        if (value < smallest) {
          numberSmall++;
        } else if (!(value <= largest)) {
          numberLarge++;
          if (firstBadColumn < 0) {
            firstBadColumn = iColumn;
            firstBadRow = row[j];
            firstBadElement = elementByColumn[j];
          }
        }
      }
      //clear mark
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        mark[iRow] = -1;
      }
    }
    delete[] mark;
  } else {
    // just check for out of range - not for duplicates
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (end != columnStart[iColumn + 1])
        flags_ |= 2;
      for (j = start; j < end; j++) {
        double value = fabs(elementByColumn[j]);
        int iRow = row[j];
        if (iRow < 0 || iRow >= numberRows) {
#ifndef COIN_BIG_INDEX
          printf("Out of range %d %d %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#elif COIN_BIG_INDEX == 0
          printf("Out of range %d %d %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#elif COIN_BIG_INDEX == 1
          printf("Out of range %d %ld %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#else
          printf("Out of range %d %lld %d %g\n", iColumn, j, row[j], elementByColumn[j]);
#endif
          return false;
        }
        if (!value)
          flags_ |= 1; // there are zero elements
        if (value < smallest) {
          numberSmall++;
        } else if (!(value <= largest)) {
          numberLarge++;
          if (firstBadColumn < 0) {
            firstBadColumn = iColumn;
            firstBadRow = iRow;
            firstBadElement = value;
          }
        }
      }
    }
  }
  if (numberLarge) {
    model->messageHandler()->message(CLP_BAD_MATRIX, model->messages())
      << numberLarge
      << firstBadColumn << firstBadRow << firstBadElement
      << CoinMessageEol;
    return false;
  }
  if (numberSmall)
    model->messageHandler()->message(CLP_SMALLELEMENTS, model->messages())
      << numberSmall
      << CoinMessageEol;
  if (numberDuplicate)
    model->messageHandler()->message(CLP_DUPLICATEELEMENTS, model->messages())
      << numberDuplicate
      << CoinMessageEol;
  if (numberDuplicate)
    matrix_->eliminateDuplicates(smallest);
  else if (numberSmall)
    matrix_->compress(smallest);
  // If smallest >0.0 then there can't be zero elements
  if (smallest > 0.0)
    flags_ &= ~1;
  ;
  if (numberSmall || numberDuplicate)
    flags_ |= 2; // will have gaps
  return true;
}
int ClpPackedMatrix::gutsOfTransposeTimesByRowGE3a(const CoinIndexedVector *COIN_RESTRICT piVector,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT output,
  int *COIN_RESTRICT lookup,
  char *COIN_RESTRICT marked,
  const double tolerance,
  const double scalar) const
{
  const double *COIN_RESTRICT pi = piVector->denseVector();
  int numberNonZero = 0;
  int numberInRowArray = piVector->getNumElements();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT element = matrix_->getElements();
  const int *COIN_RESTRICT whichRow = piVector->getIndices();
  int *COIN_RESTRICT fakeRow = const_cast< int * >(whichRow);
  fakeRow[numberInRowArray] = 0; // so can touch
#ifndef NDEBUG
  int maxColumn = 0;
#endif
  // ** Row copy is already scaled
  int nextRow = whichRow[0];
  CoinBigIndex nextStart = rowStart[nextRow];
  CoinBigIndex nextEnd = rowStart[nextRow + 1];
  for (int i = 0; i < numberInRowArray; i++) {
    double value = pi[i] * scalar;
    CoinBigIndex start = nextStart;
    CoinBigIndex end = nextEnd;
    nextRow = whichRow[i + 1];
    nextStart = rowStart[nextRow];
    //coin_prefetch_const(column + nextStart);
    //coin_prefetch_const(element + nextStart);
    nextEnd = rowStart[nextRow + 1];
    CoinBigIndex j;
    for (j = start; j < end; j++) {
      int iColumn = column[j];
#ifndef NDEBUG
      maxColumn = CoinMax(maxColumn, iColumn);
#endif
      double elValue = element[j];
      elValue *= value;
      if (marked[iColumn]) {
        int k = lookup[iColumn];
        output[k] += elValue;
      } else {
        output[numberNonZero] = elValue;
        marked[iColumn] = 1;
        lookup[iColumn] = numberNonZero;
        index[numberNonZero++] = iColumn;
      }
    }
  }
#ifndef NDEBUG
  int saveN = numberNonZero;
#endif
  // get rid of tiny values and zero out lookup
  for (int i = 0; i < numberNonZero; i++) {
    int iColumn = index[i];
    marked[iColumn] = 0;
    double value = output[i];
    if (fabs(value) <= tolerance) {
      while (fabs(value) <= tolerance) {
        numberNonZero--;
        value = output[numberNonZero];
        iColumn = index[numberNonZero];
        marked[iColumn] = 0;
        if (i < numberNonZero) {
          output[numberNonZero] = 0.0;
          output[i] = value;
          index[i] = iColumn;
        } else {
          output[i] = 0.0;
          value = 1.0; // to force end of while
        }
      }
    }
  }
#ifndef NDEBUG
  for (int i = numberNonZero; i < saveN; i++)
    assert(!output[i]);
  for (int i = 0; i <= maxColumn; i++)
    assert(!marked[i]);
#endif
  return numberNonZero;
}
int ClpPackedMatrix::gutsOfTransposeTimesByRowGE3(const CoinIndexedVector *COIN_RESTRICT piVector,
  int *COIN_RESTRICT index,
  double *COIN_RESTRICT output,
  double *COIN_RESTRICT array,
  const double tolerance,
  const double scalar) const
{
  const double *COIN_RESTRICT pi = piVector->denseVector();
  int numberNonZero = 0;
  int numberInRowArray = piVector->getNumElements();
  const int *COIN_RESTRICT column = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = matrix_->getVectorStarts();
  const double *COIN_RESTRICT element = matrix_->getElements();
  const int *COIN_RESTRICT whichRow = piVector->getIndices();
  // ** Row copy is already scaled
  for (int i = 0; i < numberInRowArray; i++) {
    int iRow = whichRow[i];
    double value = pi[i] * scalar;
    CoinBigIndex j;
    for (j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
      int iColumn = column[j];
      double inValue = array[iColumn];
      double elValue = element[j];
      elValue *= value;
      if (inValue) {
        double outValue = inValue + elValue;
        if (!outValue)
          outValue = COIN_INDEXED_REALLY_TINY_ELEMENT;
        array[iColumn] = outValue;
      } else {
        array[iColumn] = elValue;
        assert(elValue);
        index[numberNonZero++] = iColumn;
      }
    }
  }
  int saveN = numberNonZero;
  // get rid of tiny values
  numberNonZero = 0;
  for (int i = 0; i < saveN; i++) {
    int iColumn = index[i];
    double value = array[iColumn];
    array[iColumn] = 0.0;
    if (fabs(value) > tolerance) {
      output[numberNonZero] = value;
      index[numberNonZero++] = iColumn;
    }
  }
  return numberNonZero;
}
/* Given positive integer weights for each row fills in sum of weights
   for each column (and slack).
   Returns weights vector
*/
CoinBigIndex *
ClpPackedMatrix::dubiousWeights(const ClpSimplex *model, int *inputWeights) const
{
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  int number = numberRows + numberColumns;
  CoinBigIndex *weights = new CoinBigIndex[number];
  // get matrix data pointers
  const int *row = matrix_->getIndices();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const int *columnLength = matrix_->getVectorLengths();
  int i;
  for (i = 0; i < numberColumns; i++) {
    CoinBigIndex j;
    CoinBigIndex count = 0;
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      int iRow = row[j];
      count += inputWeights[iRow];
    }
    weights[i] = count;
  }
  for (i = 0; i < numberRows; i++) {
    weights[i + numberColumns] = inputWeights[i];
  }
  return weights;
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void ClpPackedMatrix::rangeOfElements(double &smallestNegative, double &largestNegative,
  double &smallestPositive, double &largestPositive)
{
  smallestNegative = -COIN_DBL_MAX;
  largestNegative = 0.0;
  smallestPositive = COIN_DBL_MAX;
  largestPositive = 0.0;
  // get matrix data pointers
  const double *elementByColumn = matrix_->getElements();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const int *columnLength = matrix_->getVectorLengths();
  int numberColumns = matrix_->getNumCols();
  int i;
  for (i = 0; i < numberColumns; i++) {
    CoinBigIndex j;
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      double value = elementByColumn[j];
      if (value > 0.0) {
        smallestPositive = CoinMin(smallestPositive, value);
        largestPositive = CoinMax(largestPositive, value);
      } else if (value < 0.0) {
        smallestNegative = CoinMax(smallestNegative, value);
        largestNegative = CoinMin(largestNegative, value);
      }
    }
  }
}
// Says whether it can do partial pricing
bool ClpPackedMatrix::canDoPartialPricing() const
{
  return true;
}
// Partial pricing
void ClpPackedMatrix::partialPricing(ClpSimplex *model, double startFraction, double endFraction,
  int &bestSequence, int &numberWanted)
{
  numberWanted = currentWanted_;
  int start = static_cast< int >(startFraction * numberActiveColumns_);
  int end = CoinMin(static_cast< int >(endFraction * numberActiveColumns_ + 1), numberActiveColumns_);
  const double *COIN_RESTRICT element = matrix_->getElements();
  const int *COIN_RESTRICT row = matrix_->getIndices();
  const CoinBigIndex *COIN_RESTRICT startColumn = matrix_->getVectorStarts();
  const int *COIN_RESTRICT length = matrix_->getVectorLengths();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const double *COIN_RESTRICT columnScale = model->columnScale();
  int iSequence;
  CoinBigIndex j;
  double tolerance = model->currentDualTolerance();
  double *COIN_RESTRICT reducedCost = model->djRegion();
  const double *COIN_RESTRICT duals = model->dualRowSolution();
  const double *COIN_RESTRICT cost = model->costRegion();
  double bestDj;
  if (bestSequence >= 0)
    bestDj = fabs(model->clpMatrix()->reducedCost(model, bestSequence));
  else
    bestDj = tolerance;
  int sequenceOut = model->sequenceOut();
  int saveSequence = bestSequence;
  int lastScan = minimumObjectsScan_ < 0 ? end : start + minimumObjectsScan_;
  int minNeg = minimumGoodReducedCosts_ == -1 ? numberWanted : minimumGoodReducedCosts_;
  if (rowScale) {
    // scaled
    for (iSequence = start; iSequence < end; iSequence++) {
      if (iSequence != sequenceOut) {
        double value;
        ClpSimplex::Status status = model->getStatus(iSequence);

        switch (status) {

        case ClpSimplex::basic:
        case ClpSimplex::isFixed:
          break;
        case ClpSimplex::isFree:
        case ClpSimplex::superBasic:
          value = 0.0;
          // scaled
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j] * rowScale[jRow];
          }
          value = fabs(cost[iSequence] + value * columnScale[iSequence]);
          if (value > FREE_ACCEPT * tolerance) {
            numberWanted--;
            // we are going to bias towards free (but only if reasonable)
            value *= FREE_BIAS;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        case ClpSimplex::atUpperBound:
          value = 0.0;
          // scaled
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j] * rowScale[jRow];
          }
          value = cost[iSequence] + value * columnScale[iSequence];
          if (value > tolerance) {
            numberWanted--;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        case ClpSimplex::atLowerBound:
          value = 0.0;
          // scaled
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j] * rowScale[jRow];
          }
          value = -(cost[iSequence] + value * columnScale[iSequence]);
          if (value > tolerance) {
            numberWanted--;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        }
      }
      if (numberWanted + minNeg < originalWanted_ && iSequence > lastScan) {
        // give up
        break;
      }
      if (!numberWanted)
        break;
    }
    if (bestSequence != saveSequence) {
      // recompute dj
      double value = 0.0;
      // scaled
      for (j = startColumn[bestSequence];
           j < startColumn[bestSequence] + length[bestSequence]; j++) {
        int jRow = row[j];
        value -= duals[jRow] * element[j] * rowScale[jRow];
      }
      reducedCost[bestSequence] = cost[bestSequence] + value * columnScale[bestSequence];
      savedBestSequence_ = bestSequence;
      savedBestDj_ = reducedCost[savedBestSequence_];
    }
  } else {
    // not scaled
    for (iSequence = start; iSequence < end; iSequence++) {
      if (iSequence != sequenceOut) {
        double value;
        ClpSimplex::Status status = model->getStatus(iSequence);

        switch (status) {

        case ClpSimplex::basic:
        case ClpSimplex::isFixed:
          break;
        case ClpSimplex::isFree:
        case ClpSimplex::superBasic:
          value = cost[iSequence];
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j];
          }
          value = fabs(value);
          if (value > FREE_ACCEPT * tolerance) {
            numberWanted--;
            // we are going to bias towards free (but only if reasonable)
            value *= FREE_BIAS;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        case ClpSimplex::atUpperBound:
          value = cost[iSequence];
          // scaled
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j];
          }
          if (value > tolerance) {
            numberWanted--;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        case ClpSimplex::atLowerBound:
          value = cost[iSequence];
          for (j = startColumn[iSequence];
               j < startColumn[iSequence] + length[iSequence]; j++) {
            int jRow = row[j];
            value -= duals[jRow] * element[j];
          }
          value = -value;
          if (value > tolerance) {
            numberWanted--;
            if (value > bestDj) {
              // check flagged variable and correct dj
              if (!model->flagged(iSequence)) {
                bestDj = value;
                bestSequence = iSequence;
              } else {
                // just to make sure we don't exit before got something
                numberWanted++;
              }
            }
          }
          break;
        }
      }
      if (numberWanted + minNeg < originalWanted_ && iSequence > lastScan) {
        // give up
        break;
      }
      if (!numberWanted)
        break;
    }
    if (bestSequence != saveSequence) {
      // recompute dj
      double value = cost[bestSequence];
      for (j = startColumn[bestSequence];
           j < startColumn[bestSequence] + length[bestSequence]; j++) {
        int jRow = row[j];
        value -= duals[jRow] * element[j];
      }
      reducedCost[bestSequence] = value;
      savedBestSequence_ = bestSequence;
      savedBestDj_ = reducedCost[savedBestSequence_];
    }
  }
  currentWanted_ = numberWanted;
}
// Sets up an effective RHS
void ClpPackedMatrix::useEffectiveRhs(ClpSimplex *model)
{
  delete[] rhsOffset_;
  int numberRows = model->numberRows();
  rhsOffset_ = new double[numberRows];
  rhsOffset(model, true);
}
// Gets rid of special copies
void ClpPackedMatrix::clearCopies()
{
  delete rowCopy_;
  delete columnCopy_;
  rowCopy_ = NULL;
  columnCopy_ = NULL;
  flags_ &= ~(4 + 8);
  checkGaps();
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
// makes sure active columns correct
int ClpPackedMatrix::refresh(ClpSimplex *)
{
  numberActiveColumns_ = matrix_->getNumCols();
#if 0
     ClpMatrixBase * rowCopyBase = reverseOrderedCopy();
     ClpPackedMatrix* rowCopy =
          dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
     // Make sure it is really a ClpPackedMatrix
     assert (rowCopy != NULL);

     const int * column = rowCopy->matrix_->getIndices();
     const CoinBigIndex * rowStart = rowCopy->matrix_->getVectorStarts();
     const int * rowLength = rowCopy->matrix_->getVectorLengths();
     const double * element = rowCopy->matrix_->getElements();
     int numberRows = rowCopy->matrix_->getNumRows();
     for (int i = 0; i < numberRows; i++) {
          if (!rowLength[i])
               printf("zero row %d\n", i);
     }
     delete rowCopy;
#endif
  checkGaps();
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
  return 0;
}

/* Scales rowCopy if column copy scaled
   Only called if scales already exist */
void ClpPackedMatrix::scaleRowCopy(ClpModel *model) const
{
  if (model->rowCopy()) {
    // need to replace row by row
    int numberRows = model->numberRows();
#ifndef NDEBUG
    int numberColumns = matrix_->getNumCols();
#endif
    ClpMatrixBase *rowCopyBase = model->rowCopy();
#ifndef NDEBUG
    ClpPackedMatrix *rowCopy = dynamic_cast< ClpPackedMatrix * >(rowCopyBase);
    // Make sure it is really a ClpPackedMatrix
    assert(rowCopy != NULL);
#else
    ClpPackedMatrix *rowCopy = static_cast< ClpPackedMatrix * >(rowCopyBase);
#endif

    const int *COIN_RESTRICT column = rowCopy->getIndices();
    const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
    double *COIN_RESTRICT element = rowCopy->getMutableElements();
    const double *COIN_RESTRICT rowScale = model->rowScale();
    const double *COIN_RESTRICT columnScale = model->columnScale();
    // scale row copy
    for (int iRow = 0; iRow < numberRows; iRow++) {
      CoinBigIndex j;
      double scale = rowScale[iRow];
      double *COIN_RESTRICT elementsInThisRow = element + rowStart[iRow];
      const int *COIN_RESTRICT columnsInThisRow = column + rowStart[iRow];
      int number = static_cast< int >(rowStart[iRow + 1] - rowStart[iRow]);
      assert(number <= numberColumns);
      for (j = 0; j < number; j++) {
        int iColumn = columnsInThisRow[j];
        elementsInThisRow[j] *= scale * columnScale[iColumn];
      }
    }
  }
}
/* Realy really scales column copy
   Only called if scales already exist.
   Up to user ro delete */
ClpMatrixBase *
ClpPackedMatrix::scaledColumnCopy(ClpModel *model) const
{
  // need to replace column by column
#ifndef NDEBUG
  int numberRows = model->numberRows();
#endif
  int numberColumns = matrix_->getNumCols();
  ClpPackedMatrix *copy = new ClpPackedMatrix(*this);
  const int *COIN_RESTRICT row = copy->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = copy->getVectorStarts();
  const int *COIN_RESTRICT length = copy->getVectorLengths();
  double *COIN_RESTRICT element = copy->getMutableElements();
  const double *COIN_RESTRICT rowScale = model->rowScale();
  const double *COIN_RESTRICT columnScale = model->columnScale();
  // scale column copy
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    double scale = columnScale[iColumn];
    double *COIN_RESTRICT elementsInThisColumn = element + columnStart[iColumn];
    const int *COIN_RESTRICT rowsInThisColumn = row + columnStart[iColumn];
    int number = length[iColumn];
    assert(number <= numberRows);
    for (j = 0; j < number; j++) {
      int iRow = rowsInThisColumn[j];
      elementsInThisColumn[j] *= scale * rowScale[iRow];
    }
  }
  return copy;
}
// Really scale matrix
void ClpPackedMatrix::reallyScale(const double *rowScale, const double *columnScale)
{
  clearCopies();
  int numberColumns = matrix_->getNumCols();
  const int *row = matrix_->getIndices();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const int *length = matrix_->getVectorLengths();
  double *element = matrix_->getMutableElements();
  // scale
  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex j;
    double scale = columnScale[iColumn];
    for (j = columnStart[iColumn]; j < columnStart[iColumn] + length[iColumn]; j++) {
      int iRow = row[j];
      element[j] *= scale * rowScale[iRow];
    }
  }
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void ClpPackedMatrix::deleteCols(const int numDel, const int *indDel)
{
  if (matrix_->getNumCols())
    matrix_->deleteCols(numDel, indDel);
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  checkGaps();
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
  matrix_->setExtraGap(0.0);
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void ClpPackedMatrix::deleteRows(const int numDel, const int *indDel)
{
  if (matrix_->getNumRows())
    matrix_->deleteRows(numDel, indDel);
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  checkGaps();
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
  matrix_->setExtraGap(0.0);
}
#ifndef CLP_NO_VECTOR
// Append Columns
void ClpPackedMatrix::appendCols(int number, const CoinPackedVectorBase *const *columns)
{
  matrix_->appendCols(number, columns);
  numberActiveColumns_ = matrix_->getNumCols();
  clearCopies();
}
// Append Rows
void ClpPackedMatrix::appendRows(int number, const CoinPackedVectorBase *const *rows)
{
  matrix_->appendRows(number, rows);
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  checkGaps();
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
  clearCopies();
}
#endif
/* Set the dimensions of the matrix. In effect, append new empty
   columns/rows to the matrix. A negative number for either dimension
   means that that dimension doesn't change. Otherwise the new dimensions
   MUST be at least as large as the current ones otherwise an exception
   is thrown. */
void ClpPackedMatrix::setDimensions(int numrows, int numcols)
{
  matrix_->setDimensions(numrows, numcols);
#ifdef DO_CHECK_FLAGS
  checkFlags(0);
#endif
}
/* Append a set of rows/columns to the end of the matrix. Returns number of errors
   i.e. if any of the new rows/columns contain an index that's larger than the
   number of columns-1/rows-1 (if numberOther>0) or duplicates
   If 0 then rows, 1 if columns */
int ClpPackedMatrix::appendMatrix(int number, int type,
  const CoinBigIndex *starts, const int *index,
  const double *element, int numberOther)
{
  int numberErrors = 0;
  // make sure other dimension is big enough
  if (type == 0) {
    // rows
    if (matrix_->isColOrdered() && numberOther > matrix_->getNumCols())
      matrix_->setDimensions(-1, numberOther);
    if (!matrix_->isColOrdered() || numberOther >= 0 || matrix_->getExtraGap()) {
      numberErrors = matrix_->appendRows(number, starts, index, element, numberOther);
    } else {
      //CoinPackedMatrix mm(*matrix_);
      matrix_->appendMinorFast(number, starts, index, element);
      //mm.appendRows(number,starts,index,element,numberOther);
      //if (!mm.isEquivalent(*matrix_)) {
      //printf("bad append\n");
      //abort();
      //}
    }
  } else {
    // columns
    if (!matrix_->isColOrdered() && numberOther > matrix_->getNumRows())
      matrix_->setDimensions(numberOther, -1);
    if (element)
      numberErrors = matrix_->appendCols(number, starts, index, element, numberOther);
    else
      matrix_->setDimensions(-1, matrix_->getNumCols() + number); // resize
  }
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  return numberErrors;
}
void ClpPackedMatrix::specialRowCopy(ClpSimplex *model, const ClpMatrixBase *rowCopy)
{
  delete rowCopy_;
  rowCopy_ = new ClpPackedMatrix2(model, rowCopy->getPackedMatrix());
  // See if anything in it
  if (!rowCopy_->usefulInfo()) {
    delete rowCopy_;
    rowCopy_ = NULL;
    flags_ &= ~4;
  } else {
    flags_ |= 4;
  }
}
void ClpPackedMatrix::specialColumnCopy(ClpSimplex *model)
{
  delete columnCopy_;
  if (model->vectorMode() == 1) {
    flags_ |= 16;
    // go to exact devex (unless full steepest)
    ClpPrimalColumnSteepest *pricing = dynamic_cast< ClpPrimalColumnSteepest * >(model->primalColumnPivot());
    if (pricing && pricing->mode() > 1)
      pricing->setMode(0);
  }
  if ((flags_ & 16) != 0 && model->numberRows() > 200 && model->numberColumns() > 500) {
    columnCopy_ = new ClpPackedMatrix3(model, matrix_);
    flags_ |= 8;
  } else {
    columnCopy_ = NULL;
  }
}
// Say we don't want special column copy
void ClpPackedMatrix::releaseSpecialColumnCopy()
{
  flags_ &= ~(8 + 16);
  delete columnCopy_;
  columnCopy_ = NULL;
}
// Correct sequence in and out to give true value
void ClpPackedMatrix::correctSequence(const ClpSimplex *model, int &sequenceIn, int &sequenceOut)
{
  if (columnCopy_) {
    if (sequenceIn != -999) {
      columnCopy_->swapOne(model, this, sequenceIn);
      if (sequenceIn != sequenceOut)
        columnCopy_->swapOne(model, this, sequenceOut);
    } else {
      // do all
      columnCopy_->sortBlocks(model);
    }
#ifndef NDEBUG
    columnCopy_->checkBlocks(model);
#endif
  }
}
// Check validity
void ClpPackedMatrix::checkFlags(int type) const
{
  int iColumn;
  // get matrix data pointers
  //const int * row = matrix_->getIndices();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const int *columnLength = matrix_->getVectorLengths();
  const double *elementByColumn = matrix_->getElements();
  if (!zeros()) {
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      for (j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        if (!elementByColumn[j])
          abort();
      }
    }
  }
  if ((flags_ & 2) == 0) {
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      if (columnStart[iColumn + 1] != columnStart[iColumn] + columnLength[iColumn]) {
        abort();
      }
    }
  }
  if (type) {
    if ((flags_ & 2) != 0) {
      bool ok = true;
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        if (columnStart[iColumn + 1] != columnStart[iColumn] + columnLength[iColumn]) {
          ok = false;
          break;
        }
      }
      if (ok)
        COIN_DETAIL_PRINT(printf("flags_ could be 0\n"));
    }
  }
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2()
  : numberBlocks_(0)
  , numberRows_(0)
  , offset_(NULL)
  , count_(NULL)
  , rowStart_(NULL)
  , column_(NULL)
  , work_(NULL)
{
#ifdef THREAD
  threadId_ = NULL;
  info_ = NULL;
#endif
}
//-------------------------------------------------------------------
// Useful Constructor
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2(ClpSimplex *, const CoinPackedMatrix *rowCopy)
  : numberBlocks_(0)
  , numberRows_(0)
  , offset_(NULL)
  , count_(NULL)
  , rowStart_(NULL)
  , column_(NULL)
  , work_(NULL)
{
#ifdef THREAD
  threadId_ = NULL;
  info_ = NULL;
#endif
  numberRows_ = rowCopy->getNumRows();
  if (!numberRows_)
    return;
  int numberColumns = rowCopy->getNumCols();
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *length = rowCopy->getVectorLengths();
  const double *element = rowCopy->getElements();
  int chunk = 32768; // tune
  //chunk=100;
  // tune
#if 0
     int chunkY[7] = {1024, 2048, 4096, 8192, 16384, 32768, 65535};
     int its = model->maximumIterations();
     if (its >= 1000000 && its < 1000999) {
          its -= 1000000;
          its = its / 10;
          if (its >= 7) abort();
          chunk = chunkY[its];
          printf("chunk size %d\n", chunk);
#define cpuid(func, ax, bx, cx, dx)                             \
  __asm__ __volatile__("cpuid"                                  \
                       : "=a"(ax), "=b"(bx), "=c"(cx), "=d"(dx) \
                       : "a"(func));
          unsigned int a, b, c, d;
          int func = 0;
          cpuid(func, a, b, c, d);
          {
               int i;
               unsigned int value;
               value = b;
               for (i = 0; i < 4; i++) {
                    printf("%c", (value & 0xff));
                    value = value >> 8;
               }
               value = d;
               for (i = 0; i < 4; i++) {
                    printf("%c", (value & 0xff));
                    value = value >> 8;
               }
               value = c;
               for (i = 0; i < 4; i++) {
                    printf("%c", (value & 0xff));
                    value = value >> 8;
               }
               printf("\n");
               int maxfunc = a;
               if (maxfunc > 10) {
                    printf("not intel?\n");
                    abort();
               }
               for (func = 1; func <= maxfunc; func++) {
                    cpuid(func, a, b, c, d);
                    printf("func %d, %x %x %x %x\n", func, a, b, c, d);
               }
          }
#else
  if (numberColumns > 10000 || chunk == 100) {
#endif
}
else
{
  //printf("no chunk\n");
  return;
}
// Could also analyze matrix to get natural breaks
numberBlocks_ = (numberColumns + chunk - 1) / chunk;
#ifdef THREAD
// Get work areas
threadId_ = new pthread_t[numberBlocks_];
info_ = new dualColumn0Struct[numberBlocks_];
#endif
// Even out
chunk = (numberColumns + numberBlocks_ - 1) / numberBlocks_;
offset_ = new int[numberBlocks_ + 1];
offset_[numberBlocks_] = numberColumns;
int nRow = numberBlocks_ * numberRows_;
count_ = new unsigned short[nRow];
memset(count_, 0, nRow * sizeof(unsigned short));
rowStart_ = new CoinBigIndex[nRow + numberRows_ + 1];
CoinBigIndex nElement = rowStart[numberRows_];
rowStart_[nRow + numberRows_] = nElement;
column_ = new unsigned short[nElement];
// assumes int <= double
int sizeWork = 6 * numberBlocks_;
work_ = new double[sizeWork];
;
int iBlock;
//int nZero = 0;
for (iBlock = 0; iBlock < numberBlocks_; iBlock++) {
  int start = iBlock * chunk;
  offset_[iBlock] = start;
  int end = start + chunk;
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    if (rowStart[iRow + 1] != rowStart[iRow] + length[iRow]) {
      printf("not packed correctly - gaps\n");
      abort();
    }
    bool lastFound = false;
    int nFound = 0;
    for (CoinBigIndex j = rowStart[iRow];
         j < rowStart[iRow] + length[iRow]; j++) {
      int iColumn = column[j];
      if (iColumn >= start) {
        if (iColumn < end) {
          if (!element[j]) {
            printf("not packed correctly - zero element\n");
            abort();
          }
          column_[j] = static_cast< unsigned short >(iColumn - start);
          nFound++;
          if (lastFound) {
            printf("not packed correctly - out of order\n");
            abort();
          }
        } else {
          //can't find any more
          lastFound = true;
        }
      }
    }
    count_[iRow * numberBlocks_ + iBlock] = static_cast< unsigned short >(nFound);
    //if (!nFound)
    //  nZero++;
  }
}
//double fraction = ((double) nZero)/((double) (numberBlocks_*numberRows_));
//printf("%d empty blocks, %g%%\n",nZero,100.0*fraction);
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2(const ClpPackedMatrix2 &rhs)
  : numberBlocks_(rhs.numberBlocks_)
  , numberRows_(rhs.numberRows_)
{
  if (numberBlocks_) {
    offset_ = CoinCopyOfArray(rhs.offset_, numberBlocks_ + 1);
    int nRow = numberBlocks_ * numberRows_;
    count_ = CoinCopyOfArray(rhs.count_, nRow);
    rowStart_ = CoinCopyOfArray(rhs.rowStart_, nRow + numberRows_ + 1);
    CoinBigIndex nElement = rowStart_[nRow + numberRows_];
    column_ = CoinCopyOfArray(rhs.column_, nElement);
    int sizeWork = 6 * numberBlocks_;
    work_ = CoinCopyOfArray(rhs.work_, sizeWork);
#ifdef THREAD
    threadId_ = new pthread_t[numberBlocks_];
    info_ = new dualColumn0Struct[numberBlocks_];
#endif
  } else {
    offset_ = NULL;
    count_ = NULL;
    rowStart_ = NULL;
    column_ = NULL;
    work_ = NULL;
#ifdef THREAD
    threadId_ = NULL;
    info_ = NULL;
#endif
  }
}
//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPackedMatrix2::~ClpPackedMatrix2()
{
  delete[] offset_;
  delete[] count_;
  delete[] rowStart_;
  delete[] column_;
  delete[] work_;
#ifdef THREAD
  delete[] threadId_;
  delete[] info_;
#endif
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPackedMatrix2 &
ClpPackedMatrix2::operator=(const ClpPackedMatrix2 &rhs)
{
  if (this != &rhs) {
    numberBlocks_ = rhs.numberBlocks_;
    numberRows_ = rhs.numberRows_;
    delete[] offset_;
    delete[] count_;
    delete[] rowStart_;
    delete[] column_;
    delete[] work_;
#ifdef THREAD
    delete[] threadId_;
    delete[] info_;
#endif
    if (numberBlocks_) {
      offset_ = CoinCopyOfArray(rhs.offset_, numberBlocks_ + 1);
      int nRow = numberBlocks_ * numberRows_;
      count_ = CoinCopyOfArray(rhs.count_, nRow);
      rowStart_ = CoinCopyOfArray(rhs.rowStart_, nRow + numberRows_ + 1);
      CoinBigIndex nElement = rowStart_[nRow + numberRows_];
      column_ = CoinCopyOfArray(rhs.column_, nElement);
      int sizeWork = 6 * numberBlocks_;
      work_ = CoinCopyOfArray(rhs.work_, sizeWork);
#ifdef THREAD
      threadId_ = new pthread_t[numberBlocks_];
      info_ = new dualColumn0Struct[numberBlocks_];
#endif
    } else {
      offset_ = NULL;
      count_ = NULL;
      rowStart_ = NULL;
      column_ = NULL;
      work_ = NULL;
#ifdef THREAD
      threadId_ = NULL;
      info_ = NULL;
#endif
    }
  }
  return *this;
}
static int dualColumn0(const ClpSimplex *model, double *spare,
  int *spareIndex, const double *arrayTemp,
  const int *indexTemp, int numberIn,
  int offset, double acceptablePivot,
  double *upperThetaPtr, int *posFreePtr, double *freePivotPtr)
{
  // do dualColumn0
  int i;
  int numberRemaining = 0;
  double upperTheta = 1.0e31;
  double freePivot = acceptablePivot;
  int posFree = -1;
  const double *reducedCost = model->djRegion(1);
  double dualTolerance = model->dualTolerance();
  // We can also see if infeasible or pivoting on free
  double tentativeTheta = 1.0e25;
  for (i = 0; i < numberIn; i++) {
    double alpha = arrayTemp[i];
    int iSequence = indexTemp[i] + offset;
    double oldValue;
    double value;
    bool keep;

    switch (model->getStatus(iSequence)) {

    case ClpSimplex::basic:
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      oldValue = reducedCost[iSequence];
      // If free has to be very large - should come in via dualRow
      if (model->getStatus(iSequence) == ClpSimplex::isFree && fabs(alpha) < 1.0e-3)
        break;
      if (oldValue > dualTolerance) {
        keep = true;
      } else if (oldValue < -dualTolerance) {
        keep = true;
      } else {
        if (fabs(alpha) > CoinMax(10.0 * acceptablePivot, 1.0e-5))
          keep = true;
        else
          keep = false;
      }
      if (keep) {
        // free - choose largest
        if (fabs(alpha) > freePivot) {
          freePivot = fabs(alpha);
          posFree = i;
        }
      }
      break;
    case ClpSimplex::atUpperBound:
      oldValue = reducedCost[iSequence];
      value = oldValue - tentativeTheta * alpha;
      //assert (oldValue<=dualTolerance*1.0001);
      if (value > dualTolerance) {
        value = oldValue - upperTheta * alpha;
        if (value > dualTolerance && -alpha >= acceptablePivot)
          upperTheta = (oldValue - dualTolerance) / alpha;
        // add to list
        spare[numberRemaining] = alpha;
        spareIndex[numberRemaining++] = iSequence;
      }
      break;
    case ClpSimplex::atLowerBound:
      oldValue = reducedCost[iSequence];
      value = oldValue - tentativeTheta * alpha;
      //assert (oldValue>=-dualTolerance*1.0001);
      if (value < -dualTolerance) {
        value = oldValue - upperTheta * alpha;
        if (value < -dualTolerance && alpha >= acceptablePivot)
          upperTheta = (oldValue + dualTolerance) / alpha;
        // add to list
        spare[numberRemaining] = alpha;
        spareIndex[numberRemaining++] = iSequence;
      }
      break;
    }
  }
  *upperThetaPtr = upperTheta;
  *freePivotPtr = freePivot;
  *posFreePtr = posFree;
  return numberRemaining;
}
static int doOneBlock(double *array, int *index,
  const double *pi, const CoinBigIndex *rowStart, const double *element,
  const unsigned short *column, int numberInRowArray, int numberLook)
{
  int iWhich = 0;
  CoinBigIndex nextN = 0;
  CoinBigIndex nextStart = 0;
  double nextPi = 0.0;
  for (; iWhich < numberInRowArray; iWhich++) {
    nextStart = rowStart[0];
    nextN = rowStart[numberInRowArray] - nextStart;
    rowStart++;
    if (nextN) {
      nextPi = pi[iWhich];
      break;
    }
  }
  int i;
  while (iWhich < numberInRowArray) {
    double value = nextPi;
    CoinBigIndex j = nextStart;
    CoinBigIndex n = nextN;
    // get next
    iWhich++;
    for (; iWhich < numberInRowArray; iWhich++) {
      nextStart = rowStart[0];
      nextN = rowStart[numberInRowArray] - nextStart;
      rowStart++;
      if (nextN) {
        //coin_prefetch_const(element + nextStart);
        nextPi = pi[iWhich];
        break;
      }
    }
    CoinBigIndex end = j + n;
    //coin_prefetch_const(element+rowStart_[i+1]);
    //coin_prefetch_const(column_+rowStart_[i+1]);
    if (n < 100) {
      if ((n & 1) != 0) {
        unsigned int jColumn = column[j];
        array[jColumn] -= value * element[j];
        j++;
      }
      //coin_prefetch_const(column + nextStart);
      for (; j < end; j += 2) {
        unsigned int jColumn0 = column[j];
        double value0 = value * element[j];
        unsigned int jColumn1 = column[j + 1];
        double value1 = value * element[j + 1];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
      }
    } else {
      if ((n & 1) != 0) {
        unsigned int jColumn = column[j];
        array[jColumn] -= value * element[j];
        j++;
      }
      if ((n & 2) != 0) {
        unsigned int jColumn0 = column[j];
        double value0 = value * element[j];
        unsigned int jColumn1 = column[j + 1];
        double value1 = value * element[j + 1];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        j += 2;
      }
      if ((n & 4) != 0) {
        unsigned int jColumn0 = column[j];
        double value0 = value * element[j];
        unsigned int jColumn1 = column[j + 1];
        double value1 = value * element[j + 1];
        unsigned int jColumn2 = column[j + 2];
        double value2 = value * element[j + 2];
        unsigned int jColumn3 = column[j + 3];
        double value3 = value * element[j + 3];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
        j += 4;
      }
      //coin_prefetch_const(column+nextStart);
      for (; j < end; j += 8) {
        //coin_prefetch_const(element + j + 16);
        unsigned int jColumn0 = column[j];
        double value0 = value * element[j];
        unsigned int jColumn1 = column[j + 1];
        double value1 = value * element[j + 1];
        unsigned int jColumn2 = column[j + 2];
        double value2 = value * element[j + 2];
        unsigned int jColumn3 = column[j + 3];
        double value3 = value * element[j + 3];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
        //coin_prefetch_const(column + j + 16);
        jColumn0 = column[j + 4];
        value0 = value * element[j + 4];
        jColumn1 = column[j + 5];
        value1 = value * element[j + 5];
        jColumn2 = column[j + 6];
        value2 = value * element[j + 6];
        jColumn3 = column[j + 7];
        value3 = value * element[j + 7];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
      }
    }
  }
  // get rid of tiny values
  int nSmall = numberLook;
  int numberNonZero = 0;
  for (i = 0; i < nSmall; i++) {
    double value = array[i];
    array[i] = 0.0;
    if (fabs(value) > 1.0e-12) {
      array[numberNonZero] = value;
      index[numberNonZero++] = i;
    }
  }
  for (; i < numberLook; i += 4) {
    double value0 = array[i + 0];
    double value1 = array[i + 1];
    double value2 = array[i + 2];
    double value3 = array[i + 3];
    array[i + 0] = 0.0;
    array[i + 1] = 0.0;
    array[i + 2] = 0.0;
    array[i + 3] = 0.0;
    if (fabs(value0) > 1.0e-12) {
      array[numberNonZero] = value0;
      index[numberNonZero++] = i + 0;
    }
    if (fabs(value1) > 1.0e-12) {
      array[numberNonZero] = value1;
      index[numberNonZero++] = i + 1;
    }
    if (fabs(value2) > 1.0e-12) {
      array[numberNonZero] = value2;
      index[numberNonZero++] = i + 2;
    }
    if (fabs(value3) > 1.0e-12) {
      array[numberNonZero] = value3;
      index[numberNonZero++] = i + 3;
    }
  }
  return numberNonZero;
}
#ifdef THREAD
static void *doOneBlockThread(void *voidInfo)
{
  dualColumn0Struct *info = (dualColumn0Struct *)voidInfo;
  *(info->numberInPtr) = doOneBlock(info->arrayTemp, info->indexTemp, info->pi,
    info->rowStart, info->element, info->column,
    info->numberInRowArray, info->numberLook);
  return NULL;
}
static void *doOneBlockAnd0Thread(void *voidInfo)
{
  dualColumn0Struct *info = (dualColumn0Struct *)voidInfo;
  *(info->numberInPtr) = doOneBlock(info->arrayTemp, info->indexTemp, info->pi,
    info->rowStart, info->element, info->column,
    info->numberInRowArray, info->numberLook);
  *(info->numberOutPtr) = dualColumn0(info->model, info->spare,
    info->spareIndex, (const double *)info->arrayTemp,
    (const int *)info->indexTemp, *(info->numberInPtr),
    info->offset, info->acceptablePivot,
    info->upperThetaPtr, info->posFreePtr, info->freePivotPtr);
  return NULL;
}
#endif
/* Return <code>x * scalar * A in <code>z</code>.
   Note - x packed and z will be packed mode
   Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix2::transposeTimes(const ClpSimplex *model,
  const CoinPackedMatrix *rowCopy,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *spareArray,
  CoinIndexedVector *columnArray) const
{
  // See if dualColumn0 coding wanted
  bool dualColumn = model->spareIntArray_[0] == 1;
  double acceptablePivot = model->spareDoubleArray_[0];
  double upperTheta = 1.0e31;
  double freePivot = acceptablePivot;
  int posFree = -1;
  int numberRemaining = 0;
  //if (model->numberIterations()>=200000) {
  //printf("time %g\n",CoinCpuTime()-startTime);
  //exit(77);
  //}
  double *pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *index = columnArray->getIndices();
  double *array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  const int *whichRow = rowArray->getIndices();
  double *element = const_cast< double * >(rowCopy->getElements());
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  int i;
  CoinBigIndex *rowStart2 = rowStart_;
  if (!dualColumn) {
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      CoinBigIndex start = rowStart[iRow];
      *rowStart2 = start;
      unsigned short *count1 = count_ + iRow * numberBlocks_;
      int put = 0;
      for (int j = 0; j < numberBlocks_; j++) {
        put += numberInRowArray;
        start += count1[j];
        rowStart2[put] = start;
      }
      rowStart2++;
    }
  } else {
    // also do dualColumn stuff
    double *spare = spareArray->denseVector();
    int *spareIndex = spareArray->getIndices();
    const double *reducedCost = model->djRegion(0);
    double dualTolerance = model->dualTolerance();
    // We can also see if infeasible or pivoting on free
    double tentativeTheta = 1.0e25;
    int addSequence = model->numberColumns();
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      double alpha = pi[i];
      double oldValue;
      double value;
      bool keep;

      switch (model->getStatus(iRow + addSequence)) {

      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        oldValue = reducedCost[iRow];
        // If free has to be very large - should come in via dualRow
        if (model->getStatus(iRow + addSequence) == ClpSimplex::isFree && fabs(alpha) < 1.0e-3)
          break;
        if (oldValue > dualTolerance) {
          keep = true;
        } else if (oldValue < -dualTolerance) {
          keep = true;
        } else {
          if (fabs(alpha) > CoinMax(10.0 * acceptablePivot, 1.0e-5))
            keep = true;
          else
            keep = false;
        }
        if (keep) {
          // free - choose largest
          if (fabs(alpha) > freePivot) {
            freePivot = fabs(alpha);
            posFree = i + addSequence;
          }
        }
        break;
      case ClpSimplex::atUpperBound:
        oldValue = reducedCost[iRow];
        value = oldValue - tentativeTheta * alpha;
        //assert (oldValue<=dualTolerance*1.0001);
        if (value > dualTolerance) {
          value = oldValue - upperTheta * alpha;
          if (value > dualTolerance && -alpha >= acceptablePivot)
            upperTheta = (oldValue - dualTolerance) / alpha;
          // add to list
          spare[numberRemaining] = alpha;
          spareIndex[numberRemaining++] = iRow + addSequence;
        }
        break;
      case ClpSimplex::atLowerBound:
        oldValue = reducedCost[iRow];
        value = oldValue - tentativeTheta * alpha;
        //assert (oldValue>=-dualTolerance*1.0001);
        if (value < -dualTolerance) {
          value = oldValue - upperTheta * alpha;
          if (value < -dualTolerance && alpha >= acceptablePivot)
            upperTheta = (oldValue + dualTolerance) / alpha;
          // add to list
          spare[numberRemaining] = alpha;
          spareIndex[numberRemaining++] = iRow + addSequence;
        }
        break;
      }
      CoinBigIndex start = rowStart[iRow];
      *rowStart2 = start;
      unsigned short *count1 = count_ + iRow * numberBlocks_;
      int put = 0;
      for (int j = 0; j < numberBlocks_; j++) {
        put += numberInRowArray;
        start += count1[j];
        rowStart2[put] = start;
      }
      rowStart2++;
    }
  }
  double *spare = spareArray->denseVector();
  int *spareIndex = spareArray->getIndices();
  int saveNumberRemaining = numberRemaining;
  int iBlock;
  for (iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    double *dwork = work_ + 6 * iBlock;
    int *iwork = reinterpret_cast< int * >(dwork + 3);
    if (!dualColumn) {
#ifndef THREAD
      int offset = offset_[iBlock];
      int offset3 = offset;
      offset = numberNonZero;
      double *arrayTemp = array + offset;
      int *indexTemp = index + offset;
      iwork[0] = doOneBlock(arrayTemp, indexTemp, pi, rowStart_ + numberInRowArray * iBlock,
        element, column_, numberInRowArray, offset_[iBlock + 1] - offset);
      int number = iwork[0];
      for (i = 0; i < number; i++) {
        //double value = arrayTemp[i];
        //arrayTemp[i]=0.0;
        //array[numberNonZero]=value;
        index[numberNonZero++] = indexTemp[i] + offset3;
      }
#else
      int offset = offset_[iBlock];
      double *arrayTemp = array + offset;
      int *indexTemp = index + offset;
      dualColumn0Struct *infoPtr = info_ + iBlock;
      infoPtr->arrayTemp = arrayTemp;
      infoPtr->indexTemp = indexTemp;
      infoPtr->numberInPtr = &iwork[0];
      infoPtr->pi = pi;
      infoPtr->rowStart = rowStart_ + numberInRowArray * iBlock;
      infoPtr->element = element;
      infoPtr->column = column_;
      infoPtr->numberInRowArray = numberInRowArray;
      infoPtr->numberLook = offset_[iBlock + 1] - offset;
      pthread_create(&threadId_[iBlock], NULL, doOneBlockThread, infoPtr);
#endif
    } else {
#ifndef THREAD
      int offset = offset_[iBlock];
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      int offset3 = offset;
      offset = numberNonZero;
      offset2 = numberRemaining;
      double *arrayTemp = array + offset;
      int *indexTemp = index + offset;
      iwork[0] = doOneBlock(arrayTemp, indexTemp, pi, rowStart_ + numberInRowArray * iBlock,
        element, column_, numberInRowArray, offset_[iBlock + 1] - offset);
      iwork[1] = dualColumn0(model, spare + offset2,
        spareIndex + offset2,
        arrayTemp, indexTemp,
        iwork[0], offset3, acceptablePivot,
        &dwork[1], &iwork[2],
        &dwork[2]);
      int number = iwork[0];
      int numberLook = iwork[1];
#if 1
      numberRemaining += numberLook;
#else
      double *spareTemp = spare + offset2;
      const int *spareIndexTemp = spareIndex + offset2;
      for (i = 0; i < numberLook; i++) {
        double value = spareTemp[i];
        spareTemp[i] = 0.0;
        spare[numberRemaining] = value;
        spareIndex[numberRemaining++] = spareIndexTemp[i];
      }
#endif
      if (dwork[2] > freePivot) {
        freePivot = dwork[2];
        posFree = iwork[2] + numberNonZero;
      }
      upperTheta = CoinMin(dwork[1], upperTheta);
      for (i = 0; i < number; i++) {
        // double value = arrayTemp[i];
        //arrayTemp[i]=0.0;
        //array[numberNonZero]=value;
        index[numberNonZero++] = indexTemp[i] + offset3;
      }
#else
      int offset = offset_[iBlock];
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      double *arrayTemp = array + offset;
      int *indexTemp = index + offset;
      dualColumn0Struct *infoPtr = info_ + iBlock;
      infoPtr->model = model;
      infoPtr->spare = spare + offset2;
      infoPtr->spareIndex = spareIndex + offset2;
      infoPtr->arrayTemp = arrayTemp;
      infoPtr->indexTemp = indexTemp;
      infoPtr->numberInPtr = &iwork[0];
      infoPtr->offset = offset;
      infoPtr->acceptablePivot = acceptablePivot;
      infoPtr->upperThetaPtr = &dwork[1];
      infoPtr->posFreePtr = &iwork[2];
      infoPtr->freePivotPtr = &dwork[2];
      infoPtr->numberOutPtr = &iwork[1];
      infoPtr->pi = pi;
      infoPtr->rowStart = rowStart_ + numberInRowArray * iBlock;
      infoPtr->element = element;
      infoPtr->column = column_;
      infoPtr->numberInRowArray = numberInRowArray;
      infoPtr->numberLook = offset_[iBlock + 1] - offset;
      if (iBlock >= 2)
        pthread_join(threadId_[iBlock - 2], NULL);
      pthread_create(threadId_ + iBlock, NULL, doOneBlockAnd0Thread, infoPtr);
      //pthread_join(threadId_[iBlock],NULL);
#endif
    }
  }
  for (iBlock = CoinMax(0, numberBlocks_ - 2); iBlock < numberBlocks_; iBlock++) {
#ifdef THREAD
    pthread_join(threadId_[iBlock], NULL);
#endif
  }
#ifdef THREAD
  for (iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    //pthread_join(threadId_[iBlock],NULL);
    int offset = offset_[iBlock];
    double *dwork = work_ + 6 * iBlock;
    int *iwork = (int *)(dwork + 3);
    int number = iwork[0];
    if (dualColumn) {
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      int numberLook = iwork[1];
      double *spareTemp = spare + offset2;
      const int *spareIndexTemp = spareIndex + offset2;
      for (i = 0; i < numberLook; i++) {
        double value = spareTemp[i];
        spareTemp[i] = 0.0;
        spare[numberRemaining] = value;
        spareIndex[numberRemaining++] = spareIndexTemp[i];
      }
      if (dwork[2] > freePivot) {
        freePivot = dwork[2];
        posFree = iwork[2] + numberNonZero;
      }
      upperTheta = CoinMin(dwork[1], upperTheta);
    }
    double *arrayTemp = array + offset;
    const int *indexTemp = index + offset;
    for (i = 0; i < number; i++) {
      double value = arrayTemp[i];
      arrayTemp[i] = 0.0;
      array[numberNonZero] = value;
      index[numberNonZero++] = indexTemp[i] + offset;
    }
  }
#endif
  columnArray->setNumElements(numberNonZero);
  columnArray->setPackedMode(true);
  if (dualColumn) {
    model->spareDoubleArray_[0] = upperTheta;
    // and theta and alpha and sequence
    if (posFree < 0) {
      model->spareIntArray_[1] = -1;
    } else {
      const double *reducedCost = model->djRegion(0);
      double alpha;
      int numberColumns = model->numberColumns();
      if (posFree < numberColumns) {
        alpha = columnArray->denseVector()[posFree];
        posFree = columnArray->getIndices()[posFree];
      } else {
        alpha = rowArray->denseVector()[posFree - numberColumns];
        posFree = rowArray->getIndices()[posFree - numberColumns] + numberColumns;
      }
      model->spareDoubleArray_[2] = fabs(reducedCost[posFree] / alpha);
      ;
      model->spareDoubleArray_[3] = alpha;
      model->spareIntArray_[1] = posFree;
    }
    spareArray->setNumElements(numberRemaining);
    // signal done
    model->spareIntArray_[0] = -1;
  }
}
#define ALIGNMENT 32
static void *clp_align(void *memory)
{
  if (sizeof(int) == sizeof(void *) && ALIGNMENT) {
    CoinInt64 k = reinterpret_cast< CoinInt64 >(memory);
    if ((k & (ALIGNMENT - 1)) != 0) {
      k &= ~(ALIGNMENT - 1);
      k += ALIGNMENT;
      memory = reinterpret_cast< void * >(k);
    }
    return memory;
  } else {
    CoinInt64 k = reinterpret_cast< CoinInt64 >(memory);
    if ((k & (ALIGNMENT - 1)) != 0) {
      k &= ~(ALIGNMENT - 1);
      k += ALIGNMENT;
      memory = reinterpret_cast< void * >(k);
    }
    return memory;
  }
}
/* Default constructor. */
ClpPackedMatrix3::ClpPackedMatrix3()
  : numberBlocks_(0)
  , numberColumns_(0)
  , numberColumnsWithGaps_(0)
  ,
#if ABOCA_LITE
  numberChunks_(0)
  ,
#endif
  numberElements_(0)
  , maxBlockSize_(0)
  , column_(NULL)
  , start_(NULL)
  , row_(NULL)
  , element_(NULL)
  , temporary_(NULL)
  , block_(NULL)
  , ifActive_(0)
{
}
#ifdef _MSC_VER
#include <intrin.h>
#elif defined(__ARM_FEATURE_SIMD32) || defined(__ARM_NEON)
#include <arm_neon.h>
#else
//#include <immintrin.h> // deemed unnecessary in #126; if that was wrong, then try to do as suggested in #127
//#include <fmaintrin.h>
#endif
/* Constructor from copy. */
ClpPackedMatrix3::ClpPackedMatrix3(ClpSimplex *model, const CoinPackedMatrix *columnCopy)
  : numberBlocks_(0)
  , numberColumns_(0)
  , numberColumnsWithGaps_(0)
  ,
#if ABOCA_LITE
  numberChunks_(0)
  ,
#endif
  numberElements_(0)
  , maxBlockSize_(0)
  , column_(NULL)
  , start_(NULL)
  , row_(NULL)
  , element_(NULL)
  , temporary_(NULL)
  , block_(NULL)
  , ifActive_(0)
{
  //#undef COIN_AVX2
  //#define COIN_AVX2 8
  //#define NO_AVX_HARDWARE
#ifndef COIN_AVX2
#define COIN_AVX2 4
#else
#if COIN_AVX2 == 4
#ifndef NO_AVX_HARDWARE
#define HASWELL
#endif
#elif COIN_AVX2 == 8
#ifndef NO_AVX_HARDWARE
#define SKYLAKE
#endif
#else
  error
#endif
#endif
#if COIN_AVX2 == 1
#define COIN_AVX2_SHIFT 0
#elif COIN_AVX2 == 2
#define COIN_AVX2_SHIFT 1
#elif COIN_AVX2 == 4
#define COIN_AVX2_SHIFT 2
#elif COIN_AVX2 == 8
#define COIN_AVX2_SHIFT 3
#else
  error;
#endif
#define COIN_ALIGN 8 * COIN_AVX2 // later
#define COIN_ALIGN_DOUBLE COIN_AVX2
#define COIN_AVX2_CHUNK 128 // do this many at a time
#define MINBLOCK 6
#define MAXBLOCK 100
#define MAXUNROLL 10
#define roundUp(xx) ((xx + COIN_ALIGN_DOUBLE - 1) & (~(COIN_ALIGN_DOUBLE - 1)))
#define roundDown(xx) (xx & (~(COIN_ALIGN_DOUBLE - 1)))
  numberColumns_ = model->getNumCols();
  int numberColumns = columnCopy->getNumCols();
  assert(numberColumns_ >= numberColumns);
  int numberRows = columnCopy->getNumRows();
  int *counts = new int[numberRows + 1];
  CoinZeroN(counts, numberRows + 1);
  CoinBigIndex nels = 0;
  int iColumn;
  // get matrix data pointers
  const int *row = columnCopy->getIndices();
  const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
  const int *columnLength = columnCopy->getVectorLengths();
  const double *elementByColumn = columnCopy->getElements();
  unsigned char *status = model->statusArray();
  const double *lower = model->columnLower();
  const double *upper = model->columnUpper();
  int nFreeEls = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    int n = columnLength[iColumn];
    CoinBigIndex end = start + n;
    int kZero = 0;
    for (CoinBigIndex j = start; j < end; j++) {
      if (!elementByColumn[j])
        kZero++;
    }
    n -= kZero;
    nels += n;
    if ((lower[iColumn] == -COIN_DBL_MAX && upper[iColumn] == COIN_DBL_MAX) || (status[iColumn] & 3) == 0) {
      nFreeEls += n;
      n = 0;
      if ((status[iColumn] & 3) != 0) {
        status[iColumn] &= ~7;
        status[iColumn] |= 4;
      }
    }
    counts[n]++;
  }
  counts[0] += numberColumns_ - numberColumns;
  int nZeroColumns = counts[0]; // also free
  counts[0] = -1;
  int nOdd = nZeroColumns;
  CoinBigIndex nInOdd = nFreeEls;
  maxBlockSize_ = 0;
  int i;
  for (i = 1; i <= numberRows; i++) {
    int n = counts[i];
    if (n) {
      if (n < MINBLOCK || i > MAXBLOCK) {
        nOdd += n;
        counts[i] = -1;
        nInOdd += n * i;
      } else {
        numberBlocks_++;
        maxBlockSize_ = CoinMax(maxBlockSize_, n);
        ;
      }
    } else {
      counts[i] = -1;
    }
  }
  // align on boundaries
  nels = roundUp(nInOdd);
  numberColumnsWithGaps_ = nOdd;
  for (int i = 0; i <= CoinMin(MAXBLOCK, numberRows); i++) {
    if (counts[i] > 0) {
      int n = roundUp(counts[i]);
      nels += n * i;
      numberColumnsWithGaps_ += n;
    }
  }
  row_ = new int[nels + 15];
  element_ = new double[nels + 31];
  start_ = new CoinBigIndex[nOdd + 1];
  // also allow for rows!
  int numberColumnsWithGaps = roundUp(numberColumnsWithGaps_);
  numberColumnsWithGaps_ = roundUp(numberColumnsWithGaps + numberRows);
  column_ = new int[2 * numberColumnsWithGaps_];
  memset(row_, 0, nels * sizeof(int));
  memset(element_, 0, nels * sizeof(double));
  int *lookup = column_ + numberColumnsWithGaps_;
  for (int i = 0; i < numberColumnsWithGaps; i++) {
    column_[i] = -1;
    lookup[i] = -1;
  }
  for (int i = 0; i < numberRows; i++) {
    column_[i + numberColumnsWithGaps] = i + numberColumns;
    lookup[i + numberColumns] = i;
  }
  for (int i = numberRows + numberColumnsWithGaps;
       i < numberColumnsWithGaps_; i++) {
    column_[i] = -1;
    lookup[i] = -1;
  }
  // even if no blocks do a dummy one
  numberBlocks_ = CoinMax(numberBlocks_, 1);
  block_ = new blockStruct[numberBlocks_ + 1]; // +1 for slacks
  memset(block_, 0, (numberBlocks_ + 1) * sizeof(blockStruct));
  // Fill in what we can
  int nTotal = nOdd;
  // in case no blocks
  block_->startIndices_ = nTotal;
  // adjust nels so start on 8 double boundary
  double *elements2 = reinterpret_cast< double * >(clp_align(element_ + nInOdd));
  nels = elements2 - element_;
  int nBlock = 0;
  for (i = 0; i <= CoinMin(MAXBLOCK, numberRows); i++) {
    if (counts[i] > 0) {
      blockStruct *block = block_ + nBlock;
      int n = roundUp(counts[i]);
      counts[i] = nBlock; // backward pointer
      nBlock++;
      block->startIndices_ = nTotal;
      block->startElements_ = nels;
      block->numberElements_ = i;
      // up counts
      nTotal += n;
      nels += n * i;
    }
  }
  numberElements_ = nels;
  nBlock = CoinMax(nBlock, 1);
  // slacks
  block_[nBlock].numberElements_ = 0;
  block_[nBlock].numberInBlock_ = numberRows;
  block_[nBlock].startIndices_ = numberColumnsWithGaps;
  // fill
  start_[0] = 0;
  nOdd = 0;
  nInOdd = 0;
  const double *columnScale = model->columnScale();
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    int n = columnLength[iColumn];
    CoinBigIndex end = start + n;
    int kZero = 0;
    for (CoinBigIndex j = start; j < end; j++) {
      if (!elementByColumn[j])
        kZero++;
    }
    n -= kZero;
    if ((status[iColumn] & 3) == 0)
      n = 0;
    {
      int iBlock = counts[n];
      if (iBlock >= 0) {
        blockStruct *block = block_ + iBlock;
        int k = block->numberInBlock_;
        block->numberInBlock_++;
        column_[block->startIndices_ + k] = iColumn;
        lookup[iColumn] = k;
        CoinBigIndex put = block->startElements_
          + roundDown(k) * n + (k & (COIN_AVX2 - 1));
        for (CoinBigIndex j = start; j < end; j++) {
          double value = elementByColumn[j];
          if (value) {
            if (columnScale)
              value *= columnScale[iColumn];
            element_[put] = value;
            row_[put] = row[j];
            put += COIN_AVX2;
          }
        }
      } else {
        // odd ones
        for (CoinBigIndex j = start; j < end; j++) {
          double value = elementByColumn[j];
          if (value) {
            if (columnScale)
              value *= columnScale[iColumn];
            element_[nInOdd] = value;
            row_[nInOdd++] = row[j];
          }
        }
        column_[nOdd] = iColumn;
        lookup[iColumn] = -1;
        nOdd++;
        start_[nOdd] = nInOdd;
      }
    }
  }
  // This much temporary space
#if ABOCA_LITE
  temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK * 2 * ABOCA_LITE,
    -2 - COIN_ALIGN_DOUBLE);
#else
  temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK,
    -2 - COIN_ALIGN_DOUBLE);
#endif
#if ABOCA_LITE
  // balance work - maybe with 2*cpu chunks
#define COLUMN_WEIGHT 5
#define BLOCK_MULT 2
  // If do +1's then need two element weights
  double totalWork = 0.0;
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    totalWork += (COLUMN_WEIGHT + block_[iBlock].numberElements_) * block_[iBlock].numberInBlock_;
    if (model->logLevel() > 2)
      printf("block %d ncol %d nel %d total %g this %d\n",
        iBlock, block_[iBlock].numberInBlock_,
        block_[iBlock].numberElements_, totalWork,
        (COLUMN_WEIGHT + block_[iBlock].numberElements_) * block_[iBlock].numberInBlock_);
  }
  double eachWork = totalWork / (BLOCK_MULT * ABOCA_LITE);
  double thisWork = 0.0;
  numberChunks_ = 0;
  endChunk_[0] = 0;
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    thisWork += (COLUMN_WEIGHT + block_[iBlock].numberElements_) * block_[iBlock].numberInBlock_;
    if (thisWork >= eachWork) {
      totalWork -= thisWork;
      thisWork = 0.0;
      eachWork = totalWork / (BLOCK_MULT * ABOCA_LITE - numberChunks_ - 1);
      numberChunks_++;
      endChunk_[numberChunks_] = iBlock + 1;
      if (numberChunks_ == BLOCK_MULT * ABOCA_LITE)
        break;
    }
  }
  if (numberChunks_ == BLOCK_MULT * ABOCA_LITE)
    endChunk_[numberChunks_] = numberBlocks_;
  else if (endChunk_[numberChunks_] < numberBlocks_)
    endChunk_[++numberChunks_] = numberBlocks_;
  //printf("%d chunks for %d blocks\n",numberChunks_,numberBlocks_);
  assert(numberChunks_ <= 2 * ABOCA_LITE);
#endif
  delete[] counts;
}
/* Destructor */
ClpPackedMatrix3::~ClpPackedMatrix3()
{
  delete[] column_;
  delete[] start_;
  delete[] row_;
  delete[] element_;
  delete temporary_;
  delete[] block_;
}
/* The copy constructor. */
ClpPackedMatrix3::ClpPackedMatrix3(const ClpPackedMatrix3 &rhs)
  : numberBlocks_(rhs.numberBlocks_)
  , numberColumns_(rhs.numberColumns_)
  , numberColumnsWithGaps_(rhs.numberColumnsWithGaps_)
  ,
#if ABOCA_LITE
  numberChunks_(rhs.numberChunks_)
  ,
#endif
  numberElements_(rhs.numberElements_)
  , maxBlockSize_(rhs.maxBlockSize_)
  , column_(NULL)
  , start_(NULL)
  , row_(NULL)
  , element_(NULL)
  , temporary_(NULL)
  , block_(NULL)
  , ifActive_(rhs.ifActive_)
{
  if (rhs.numberBlocks_) {
    block_ = CoinCopyOfArray(rhs.block_, numberBlocks_);
    column_ = CoinCopyOfArray(rhs.column_, 2 * numberColumnsWithGaps_);
    int numberOdd = block_->startIndices_;
    start_ = CoinCopyOfArray(rhs.start_, numberOdd + 1);
    row_ = CoinCopyOfArray(rhs.row_, numberElements_);
    element_ = CoinCopyOfArray(rhs.element_, numberElements_ + 8);
    // This much temporary space
#if ABOCA_LITE
    temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK * 2 * ABOCA_LITE,
      -2 - COIN_ALIGN_DOUBLE);
#else
    temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK,
      -2 - COIN_ALIGN_DOUBLE);
#endif
#if ABOCA_LITE
    memcpy(endChunk_, rhs.endChunk_, sizeof(endChunk_));
#endif
  }
}
ClpPackedMatrix3 &
ClpPackedMatrix3::operator=(const ClpPackedMatrix3 &rhs)
{
  if (this != &rhs) {
    delete[] column_;
    delete[] start_;
    delete[] row_;
    delete[] element_;
    delete temporary_;
    temporary_ = NULL;
    delete[] block_;
    numberBlocks_ = rhs.numberBlocks_;
    numberColumns_ = rhs.numberColumns_;
    numberColumnsWithGaps_ = rhs.numberColumnsWithGaps_;
#if ABOCA_LITE
    numberChunks_ = rhs.numberChunks_;
    memcpy(endChunk_, rhs.endChunk_, sizeof(endChunk_));
#endif
    numberElements_ = rhs.numberElements_;
    maxBlockSize_ = rhs.maxBlockSize_;
    ifActive_ = rhs.ifActive_;
    if (rhs.numberBlocks_) {
      block_ = CoinCopyOfArray(rhs.block_, numberBlocks_);
      column_ = CoinCopyOfArray(rhs.column_, 2 * numberColumnsWithGaps_);
      int numberOdd = block_->startIndices_;
      start_ = CoinCopyOfArray(rhs.start_, numberOdd + 1);
      row_ = CoinCopyOfArray(rhs.row_, numberElements_);
      element_ = CoinCopyOfArray(rhs.element_, numberElements_ + 8);
      // This much temporary space
#if ABOCA_LITE
      temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK * 2 * ABOCA_LITE,
        -2 - COIN_ALIGN_DOUBLE);
#else
      temporary_ = new CoinDoubleArrayWithLength(2 * COIN_AVX2_CHUNK,
        -2 - COIN_ALIGN_DOUBLE);
#endif
    } else {
      column_ = NULL;
      start_ = NULL;
      row_ = NULL;
      element_ = NULL;
      block_ = NULL;
    }
  }
  return *this;
}
/* Sort blocks */
void ClpPackedMatrix3::sortBlocks(const ClpSimplex *model)
{
  ifActive_ = 1;
  int *lookup = column_ + numberColumnsWithGaps_;
  for (int iBlock = 0; iBlock < numberBlocks_ + 1; iBlock++) {
    blockStruct *block = block_ + iBlock;
    int numberInBlock = block->numberInBlock_;
    int nel = block->numberElements_;
    int *row = row_ + block->startElements_;
    double *element = element_ + block->startElements_;
    int *column = column_ + block->startIndices_;
    int lastPrice = 0;
    int firstNotPrice = numberInBlock - 1;
    while (lastPrice <= firstNotPrice) {
      // find first basic or fixed
      int iColumn = numberInBlock;
      for (; lastPrice <= firstNotPrice; lastPrice++) {
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) == ClpSimplex::basic || model->getColumnStatus(iColumn) == ClpSimplex::isFixed)
          break;
      }
      // find last non basic or fixed
      int jColumn = -1;
      for (; firstNotPrice > lastPrice; firstNotPrice--) {
        jColumn = column[firstNotPrice];
        if (model->getColumnStatus(jColumn) != ClpSimplex::basic && model->getColumnStatus(jColumn) != ClpSimplex::isFixed)
          break;
      }
      if (firstNotPrice > lastPrice) {
        assert(column[lastPrice] == iColumn);
        assert(column[firstNotPrice] == jColumn);
        // need to swap
        column[firstNotPrice] = iColumn;
        lookup[iColumn] = firstNotPrice;
        column[lastPrice] = jColumn;
        lookup[jColumn] = lastPrice;
        int startBit = roundDown(lastPrice);
        CoinBigIndex offset = startBit * nel + (lastPrice - startBit);
        double *elementA = element + offset;
        int *rowA = row + offset;
        startBit = roundDown(firstNotPrice);
        offset = startBit * nel + (firstNotPrice - startBit);
        double *elementB = element + offset;
        int *rowB = row + offset;
        for (int i = 0; i < nel * COIN_AVX2; i += COIN_AVX2) {
          int temp = rowA[i];
          double tempE = elementA[i];
          rowA[i] = rowB[i];
          elementA[i] = elementB[i];
          rowB[i] = temp;
          elementB[i] = tempE;
        }
        firstNotPrice--;
        lastPrice++;
      } else if (lastPrice == firstNotPrice) {
        // make sure correct side
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) != ClpSimplex::basic && model->getColumnStatus(iColumn) != ClpSimplex::isFixed)
          lastPrice++;
        break;
      }
    }
    block->firstBasic_ = lastPrice;
    // reuse ints even though names wrong
    // lower
    firstNotPrice = lastPrice - 1;
    lastPrice = 0;
    while (lastPrice <= firstNotPrice) {
      // find first upper
      int iColumn = numberInBlock;
      for (; lastPrice <= firstNotPrice; lastPrice++) {
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) == ClpSimplex::atUpperBound)
          break;
      }
      // find last upper
      int jColumn = -1;
      for (; firstNotPrice > lastPrice; firstNotPrice--) {
        jColumn = column[firstNotPrice];
        if (model->getColumnStatus(jColumn) != ClpSimplex::atUpperBound)
          break;
      }
      if (firstNotPrice > lastPrice) {
        assert(column[lastPrice] == iColumn);
        assert(column[firstNotPrice] == jColumn);
        // need to swap
        column[firstNotPrice] = iColumn;
        lookup[iColumn] = firstNotPrice;
        column[lastPrice] = jColumn;
        lookup[jColumn] = lastPrice;
        int startBit = roundDown(lastPrice);
        CoinBigIndex offset = startBit * nel + (lastPrice - startBit);
        double *elementA = element + offset;
        int *rowA = row + offset;
        startBit = roundDown(firstNotPrice);
        offset = startBit * nel + (firstNotPrice - startBit);
        double *elementB = element + offset;
        int *rowB = row + offset;
        for (int i = 0; i < nel * COIN_AVX2; i += COIN_AVX2) {
          int temp = rowA[i];
          double tempE = elementA[i];
          rowA[i] = rowB[i];
          elementA[i] = elementB[i];
          rowB[i] = temp;
          elementB[i] = tempE;
        }
        firstNotPrice--;
        lastPrice++;
      } else if (lastPrice == firstNotPrice) {
        // make sure correct side
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) != ClpSimplex::atUpperBound)
          lastPrice++;
        break;
      }
    }
    block->firstAtUpper_ = lastPrice;
    // now lower
    firstNotPrice = lastPrice - 1;
    lastPrice = 0;
    while (lastPrice <= firstNotPrice) {
      // find first basic or fixed
      int iColumn = numberInBlock;
      for (; lastPrice <= firstNotPrice; lastPrice++) {
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) == ClpSimplex::atLowerBound)
          break;
      }
      // find last non basic or fixed
      int jColumn = -1;
      for (; firstNotPrice > lastPrice; firstNotPrice--) {
        jColumn = column[firstNotPrice];
        if (model->getColumnStatus(jColumn) != ClpSimplex::atLowerBound)
          break;
      }
      if (firstNotPrice > lastPrice) {
        assert(column[lastPrice] == iColumn);
        assert(column[firstNotPrice] == jColumn);
        // need to swap
        column[firstNotPrice] = iColumn;
        lookup[iColumn] = firstNotPrice;
        column[lastPrice] = jColumn;
        lookup[jColumn] = lastPrice;
        int startBit = roundDown(lastPrice);
        CoinBigIndex offset = startBit * nel + (lastPrice - startBit);
        double *elementA = element + offset;
        int *rowA = row + offset;
        startBit = roundDown(firstNotPrice);
        offset = startBit * nel + (firstNotPrice - startBit);
        double *elementB = element + offset;
        int *rowB = row + offset;
        for (int i = 0; i < nel * COIN_AVX2; i += COIN_AVX2) {
          int temp = rowA[i];
          double tempE = elementA[i];
          rowA[i] = rowB[i];
          elementA[i] = elementB[i];
          rowB[i] = temp;
          elementB[i] = tempE;
        }
        firstNotPrice--;
        lastPrice++;
      } else if (lastPrice == firstNotPrice) {
        // make sure correct side
        iColumn = column[lastPrice];
        if (model->getColumnStatus(iColumn) != ClpSimplex::atLowerBound)
          lastPrice++;
        break;
      }
    }
    block->firstAtLower_ = lastPrice;
#ifndef NDEBUG
    // check
    int i;
    for (i = 0; i < block->firstBasic_; i++) {
      int iColumn = column[i];
      assert(model->getColumnStatus(iColumn) != ClpSimplex::basic && model->getColumnStatus(iColumn) != ClpSimplex::isFixed);
      assert(lookup[iColumn] == i);
      if (i < block->firstAtLower_) {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::isFree || model->getColumnStatus(iColumn) == ClpSimplex::superBasic);
      } else if (i < block->firstAtUpper_) {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::atLowerBound);
      } else {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::atUpperBound);
      }
    }
    for (; i < numberInBlock; i++) {
      int iColumn = column[i];
      assert(model->getColumnStatus(iColumn) == ClpSimplex::basic || model->getColumnStatus(iColumn) == ClpSimplex::isFixed);
      assert(lookup[iColumn] == i);
    }
#endif
  }
}
// Swap one variable - when found
void ClpPackedMatrix3::swapOne(int iBlock, int kA, int kB)
{
  int *lookup = column_ + numberColumnsWithGaps_;
  blockStruct *block = block_ + iBlock;
  int nel = block->numberElements_;
  int *row = row_ + block->startElements_;
  double *element = element_ + block->startElements_;
  int *column = column_ + block->startIndices_;
  int iColumn = column[kA];
  int jColumn = column[kB];
  column[kA] = jColumn;
  lookup[jColumn] = kA;
  column[kB] = iColumn;
  lookup[iColumn] = kB;
  int startBit = roundDown(kA);
  CoinBigIndex offset = startBit * nel + (kA - startBit);
  double *elementA = element + offset;
  int *rowA = row + offset;
  startBit = roundDown(kB);
  offset = startBit * nel + (kB - startBit);
  double *elementB = element + offset;
  int *rowB = row + offset;
  int i;
  for (i = 0; i < nel * COIN_AVX2; i += COIN_AVX2) {
    int temp = rowA[i];
    double tempE = elementA[i];
    rowA[i] = rowB[i];
    elementA[i] = elementB[i];
    rowB[i] = temp;
    elementB[i] = tempE;
  }
}
// Swap one variable
void ClpPackedMatrix3::swapOne(const ClpSimplex *model, const ClpPackedMatrix *matrix,
  int iColumn)
{
  if (!ifActive_)
    return;
  int *lookup = column_ + numberColumnsWithGaps_;
  // position in block
  int kA = lookup[iColumn];
  if (kA < 0)
    return; // odd one
  int iBlock = numberBlocks_;
  if (iColumn < model->numberColumns()) {
    // get matrix data pointers
    const CoinPackedMatrix *columnCopy = matrix->getPackedMatrix();
    //const int * row = columnCopy->getIndices();
    const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
    const int *columnLength = columnCopy->getVectorLengths();
    const double *elementByColumn = columnCopy->getElements();
    CoinBigIndex start = columnStart[iColumn];
    int n = columnLength[iColumn];
    if (matrix->zeros()) {
      CoinBigIndex end = start + n;
      for (CoinBigIndex j = start; j < end; j++) {
        if (!elementByColumn[j])
          n--;
      }
    }
    // find block - could do binary search
    iBlock = CoinMin(n, numberBlocks_) - 1;
    while (block_[iBlock].numberElements_ != n)
      iBlock--;
  }
  blockStruct *block = block_ + iBlock;
#ifndef NDEBUG
  int *column = column_ + block->startIndices_;
  assert(column[kA] == iColumn);
#endif
  int from, to;
  if (kA < block->firstBasic_) {
    if (kA >= block->firstAtUpper_) {
      from = 2;
    } else if (kA >= block->firstAtLower_) {
      from = 1;
    } else {
      from = 0;
    }
  } else {
    from = 3;
  }
  if (model->getColumnStatus(iColumn) == ClpSimplex::basic || model->getColumnStatus(iColumn) == ClpSimplex::isFixed) {
    to = 3;
  } else if (model->getColumnStatus(iColumn) == ClpSimplex::atUpperBound) {
    to = 2;
  } else if (model->getColumnStatus(iColumn) == ClpSimplex::atLowerBound) {
    to = 1;
  } else {
    to = 0;
  }
  int *statusCounts = (&block->firstAtLower_) - 1;
  if (from < to) {
    while (from < to) {
      int kB = statusCounts[from + 1] - 1;
      statusCounts[from + 1] = kB;
      swapOne(iBlock, kA, kB);
      kA = kB;
      from++;
    }
  } else if (from > to) {
    while (from > to) {
      int kB = statusCounts[from];
      statusCounts[from] = kB + 1;
      swapOne(iBlock, kA, kB);
      kA = kB;
      from--;
    }
  }
#ifndef NDEBUG
  // check
  int i;
  for (i = 0; i < block->firstBasic_; i++) {
    int iColumn = column[i];
    if (iColumn != model->sequenceIn() && iColumn != model->sequenceOut())
      assert(model->getColumnStatus(iColumn) != ClpSimplex::basic && model->getColumnStatus(iColumn) != ClpSimplex::isFixed);
    assert(lookup[iColumn] == i);
    if (model->algorithm() > 0) {
      if (i < block->firstAtLower_) {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::isFree || model->getColumnStatus(iColumn) == ClpSimplex::superBasic);
      } else if (i < block->firstAtUpper_) {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::atLowerBound);
      } else {
        assert(model->getColumnStatus(iColumn) == ClpSimplex::atUpperBound);
      }
    }
  }
  int numberInBlock = block->numberInBlock_;
  for (; i < numberInBlock; i++) {
    int iColumn = column[i];
    if (iColumn != model->sequenceIn() && iColumn != model->sequenceOut())
      assert(model->getColumnStatus(iColumn) == ClpSimplex::basic || model->getColumnStatus(iColumn) == ClpSimplex::isFixed);
    assert(lookup[iColumn] == i);
  }
#endif
}
#ifndef NDEBUG
/* Debug - check blocks */
void ClpPackedMatrix3::checkBlocks(const ClpSimplex *model)
{
  if (!ifActive_)
    return;
  for (int iBlock = 0; iBlock < numberBlocks_ + 1; iBlock++) {
    blockStruct *block = block_ + iBlock;
    int *column = column_ + block->startIndices_;
    for (int i = 0; i < block->firstAtLower_; i++) {
      int iSequence = column[i];
      assert(model->getColumnStatus(iSequence) == ClpSimplex::isFree || model->getColumnStatus(iSequence) == ClpSimplex::superBasic);
    }
    for (int i = block->firstAtLower_; i < block->firstAtUpper_; i++) {
      int iSequence = column[i];
      assert(model->getColumnStatus(iSequence) == ClpSimplex::atLowerBound);
    }
    for (int i = block->firstAtUpper_; i < block->firstBasic_; i++) {
      int iSequence = column[i];
      assert(model->getColumnStatus(iSequence) == ClpSimplex::atUpperBound);
    }
    for (int i = block->firstBasic_; i < block->numberInBlock_; i++) {
      int iSequence = column[i];
      assert(model->getColumnStatus(iSequence) == ClpSimplex::basic || model->getColumnStatus(iSequence) == ClpSimplex::isFixed);
    }
  }
}
#endif
/* Return <code>x * -1 * A in <code>z</code>.
   Note - x packed and z will be packed mode
   Squashes small elements and knows about ClpSimplex */
void ClpPackedMatrix3::transposeTimes(const ClpSimplex *model,
  const double *pi,
  CoinIndexedVector *output) const
{
  int numberNonZero = 0;
  int *index = output->getIndices();
  double *array = output->denseVector();
  double zeroTolerance = model->zeroTolerance();
  double value = 0.0;
  CoinBigIndex j;
  int numberOdd = block_->startIndices_;
  if (numberOdd) {
    // A) as probably long may be worth unrolling
    CoinBigIndex end = start_[1];
    for (j = start_[0]; j < end; j++) {
      int iRow = row_[j];
      value += pi[iRow] * element_[j];
    }
    int iColumn;
    // int jColumn=column_[0];

    for (iColumn = 0; iColumn < numberOdd - 1; iColumn++) {
      CoinBigIndex start = end;
      end = start_[iColumn + 2];
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = column_[iColumn];
        //index[numberNonZero++]=jColumn;
      }
      // jColumn = column_[iColumn+1];
      value = 0.0;
      //if (model->getColumnStatus(jColumn)!=ClpSimplex::basic) {
      for (j = start; j < end; j++) {
        int iRow = row_[j];
        value += pi[iRow] * element_[j];
      }
      //}
    }
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = column_[iColumn];
      //index[numberNonZero++]=jColumn;
    }
  }
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    // C) Can do two at a time (if so put odd one into start_)
    // D) can use switch
    blockStruct *block = block_ + iBlock;
    //int numberPrice = block->numberInBlock_;
    int numberPrice = block->firstBasic_;
    int nel = block->numberElements_;
    int *row = row_ + block->startElements_;
    double *element = element_ + block->startElements_;
    int *column = column_ + block->startIndices_;
    int nBlock = numberPrice >> COIN_AVX2_SHIFT;
    numberPrice -= roundDown(numberPrice);
#if defined(HASWELL) || defined(SKYLAKE)
    double *newValues = roundUpDouble((array + numberNonZero));
#ifdef HASWELL
    assert(COIN_AVX2 == 4);
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      __m256d arrayX = _mm256_setzero_pd();
      for (int j = 0; j < nel; j++) {
        __m128i rows = _mm_loadu_si128((const __m128i *)row);
        __m256d elements = _mm256_load_pd(element);
        __m256d pis = _mm256_i32gather_pd(pi, rows, 8);
        arrayX = _mm256_fmadd_pd(pis, elements, arrayX);
        row += 4;
        element += 4;
      }
      _mm256_store_pd(newValues, arrayX);
#else
    assert(COIN_AVX2 == 8);
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      __m512d arrayX = _mm512_setzero_pd();
      for (int j = 0; j < nel; j++) {
        __m256i rows = _mm256_load_si256((const __m256i *)row);
        __m512d elements = _mm512_load_pd(element);
        __m512d pis = _mm512_i32gather_pd(rows, pi, 8);
        arrayX = _mm512_fmadd_pd(pis, elements, arrayX);
        row += 4;
        element += 4;
      }
      _mm512_store_pd(newValues, arrayX);
#endif
      newValues += COIN_AVX2;
      //row += (nel-1)*COIN_AVX2;
      //element += (nel-1)*COIN_AVX2;
      assert(row == row_ + block->startElements_ + nel * COIN_AVX2 * (jBlock + 1));
    }
    int n = nBlock * COIN_AVX2;
    int nSave = static_cast< int >(newValues - array);
    newValues = roundUpDouble((array + numberNonZero));
    for (int j = 0; j < n; j++) {
      double value = newValues[j];
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = *column;
      }
      column++;
    }
    for (int j = numberNonZero; j < nSave; j++)
      array[j] = 0.0;
#else
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      for (int j = 0; j < COIN_AVX2; j++) {
        double value = 0.0;
        for (int i = 0; i < nel; i++) {
          int iRow = row[i * COIN_AVX2];
          value += pi[iRow] * element[i * COIN_AVX2];
        }
#if COIN_AVX2 > 1
        row++;
        element++;
#else
        row += nel;
        element += nel;
#endif
        if (fabs(value) > zeroTolerance) {
          array[numberNonZero] = value;
          index[numberNonZero++] = *column;
        }
        column++;
      }
#if COIN_AVX2 > 1
      row += (nel - 1) * COIN_AVX2;
      element += (nel - 1) * COIN_AVX2;
      assert(row == row_ + block->startElements_ + nel * COIN_AVX2 * (jBlock + 1));
#endif
    }
#endif
    // last lot
    for (int j = 0; j < numberPrice; j++) {
      double value = 0.0;
      for (int i = 0; i < nel; i++) {
        int iRow = row[i * COIN_AVX2];
        value += pi[iRow] * element[i * COIN_AVX2];
      }
#if COIN_AVX2 > 1
      row++;
      element++;
#else
      row += nel;
      element += nel;
#endif
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = *column;
      }
      column++;
    }
  }
  output->setNumElements(numberNonZero);
}
/* Return <code>x * -1 * A in <code>z</code>.
   Note - x packed and z will be packed mode
   Squashes small elements and knows about ClpSimplex
   - does dualColumn0 */
void ClpPackedMatrix3::transposeTimes(const ClpSimplex *model,
  const double *COIN_RESTRICT pi,
  CoinIndexedVector *output,
  CoinIndexedVector *candidate,
  const CoinIndexedVector *rowArray) const
{
  int numberNonZero = 0;
  int *index = output->getIndices();
  double *array = output->denseVector();
  double zeroTolerance = model->zeroTolerance();
  int numberColumns = model->numberColumns();
  const unsigned char *COIN_RESTRICT statusArray = model->statusArray() + numberColumns;
  int numberInRowArray = rowArray->getNumElements();
  const int *COIN_RESTRICT whichRow = rowArray->getIndices();
  const double *COIN_RESTRICT piOld = rowArray->denseVector();
  int *COIN_RESTRICT spareIndex = candidate->getIndices();
  double *COIN_RESTRICT spareArray = candidate->denseVector();
  const double *COIN_RESTRICT reducedCost = model->djRegion(0);
  double multiplier[] = { -1.0, 1.0 };
  double dualT = -model->currentDualTolerance();
  double acceptablePivot = model->spareDoubleArray_[0];
  double tentativeTheta = 1.0e15;
  double upperTheta = 1.0e31;
  int numberRemaining = 0;
  // dualColumn0 for slacks
  for (int i = 0; i < numberInRowArray; i++) {
    int iSequence = whichRow[i];
    int iStatus = (statusArray[iSequence] & 3) - 1;
    if (iStatus) {
      double mult = multiplier[iStatus - 1];
      double alpha = piOld[i] * mult;
      double oldValue;
      double value;
      if (alpha > 0.0) {
        oldValue = reducedCost[iSequence] * mult;
        value = oldValue - tentativeTheta * alpha;
        if (value < dualT) {
          value = oldValue - upperTheta * alpha;
          if (value < dualT && alpha >= acceptablePivot) {
            upperTheta = (oldValue - dualT) / alpha;
          }
          // add to list
          spareArray[numberRemaining] = alpha * mult;
          spareIndex[numberRemaining++] = iSequence + numberColumns;
        }
      }
    }
  }
  statusArray -= numberColumns;
  reducedCost -= numberColumns;
  double value = 0.0;
  CoinBigIndex j;
  int numberOdd = block_->startIndices_;
  if (numberOdd) {
    // A) as probably long may be worth unrolling
    CoinBigIndex end = start_[1];
    for (j = start_[0]; j < end; j++) {
      int iRow = row_[j];
      value += pi[iRow] * element_[j];
    }
    int iColumn;
    // int jColumn=column_[0];

    for (iColumn = 0; iColumn < numberOdd - 1; iColumn++) {
      CoinBigIndex start = end;
      end = start_[iColumn + 2];
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = column_[iColumn];
        //index[numberNonZero++]=jColumn;
      }
      // jColumn = column_[iColumn+1];
      value = 0.0;
      //if (model->getColumnStatus(jColumn)!=ClpSimplex::basic) {
      for (j = start; j < end; j++) {
        int iRow = row_[j];
        value += pi[iRow] * element_[j];
      }
      //}
    }
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = column_[iColumn];
      //index[numberNonZero++]=jColumn;
    }
  }
  // do odd ones
  for (int i = 0; i < numberNonZero; i++) {
    int iSequence = index[i];
    double alpha;
    double oldValue;
    double value;
    int iStatus = (statusArray[iSequence] & 3) - 1;
    if (iStatus) {
      double mult = multiplier[iStatus - 1];
      alpha = array[i] * mult;
      if (alpha > 0.0) {
        oldValue = reducedCost[iSequence] * mult;
        value = oldValue - tentativeTheta * alpha;
        if (value < dualT) {
          value = oldValue - upperTheta * alpha;
          if (value < dualT && alpha >= acceptablePivot) {
            upperTheta = (oldValue - dualT) / alpha;
          }
          // add to list
          spareArray[numberRemaining] = alpha * mult;
          spareIndex[numberRemaining++] = iSequence;
        }
      }
    }
  }
  int nMax = 0;
  for (int iBlock = 0; iBlock < numberBlocks_; iBlock++) {
    // C) Can do two at a time (if so put odd one into start_)
    // D) can use switch
    blockStruct *block = block_ + iBlock;
    //int numberPrice = block->numberInBlock_;
    int numberPrice = block->firstBasic_;
    int nel = block->numberElements_;
    const int *COIN_RESTRICT row = row_ + block->startElements_;
    const double *COIN_RESTRICT element = element_ + block->startElements_;
    const int *COIN_RESTRICT column = column_ + block->startIndices_;
    int nBlock = numberPrice >> COIN_AVX2_SHIFT;
    numberPrice -= roundDown(numberPrice);
#if defined(HASWELL) || defined(SKYLAKE)
    double *COIN_RESTRICT newValues = roundUpDouble((array + numberNonZero));
#ifdef HASWELL
    assert(COIN_AVX2 == 4);
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      __m256d arrayX = _mm256_setzero_pd();
      for (int j = 0; j < nel; j++) {
        __m128i rows = _mm_loadu_si128((const __m128i *)row);
        __m256d elements = _mm256_load_pd(element);
        __m256d pis = _mm256_i32gather_pd(pi, rows, 8);
        arrayX = _mm256_fmadd_pd(pis, elements, arrayX);
        row += 4;
        element += 4;
      }
      _mm256_store_pd(newValues, arrayX);
#else
    assert(COIN_AVX2 == 8);
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      __m512d arrayX = _mm512_setzero_pd();
      for (int j = 0; j < nel; j++) {
        __m256i rows = _mm256_loadu_si256((const __m256i *)row);
        __m512d elements = _mm512_load_pd(element);
        __m512d pis = _mm512_i32gather_pd(rows, pi, 8);
        arrayX = _mm512_fmadd_pd(pis, elements, arrayX);
        row += 4;
        element += 4;
      }
      _mm512_store_pd(newValues, arrayX);
#endif
      newValues += COIN_AVX2;
      //row += (nel-1)*COIN_AVX2;
      //element += (nel-1)*COIN_AVX2;
      assert(row == row_ + block->startElements_ + nel * COIN_AVX2 * (jBlock + 1));
    }
#else
    double *COIN_RESTRICT newValues = array + numberNonZero;
    for (int jBlock = 0; jBlock < nBlock; jBlock++) {
      for (int j = 0; j < COIN_AVX2; j++) {
        double value = 0.0;
        for (int i = 0; i < nel; i++) {
          int iRow = row[i * COIN_AVX2];
          value += pi[iRow] * element[i * COIN_AVX2];
        }
#if COIN_AVX2 > 1
        row++;
        element++;
#else
        row += nel;
        element += nel;
#endif
        *newValues = value;
        newValues++;
      }
#if COIN_AVX2 > 1
      row += (nel - 1) * COIN_AVX2;
      element += (nel - 1) * COIN_AVX2;
      assert(row == row_ + block->startElements_ + nel * COIN_AVX2 * (jBlock + 1));
#endif
    }
#endif
    for (int j = 0; j < numberPrice; j++) {
      double value = 0.0;
      for (int i = 0; i < nel; i++) {
        int iRow = row[i * COIN_AVX2];
        value += pi[iRow] * element[i * COIN_AVX2];
      }
#if COIN_AVX2 > 1
      row++;
      element++;
#else
      row += nel;
      element += nel;
#endif
      *newValues = value;
      newValues++;
    }
#ifdef HASWELL
    newValues = roundUpDouble((array + numberNonZero));
#else
    newValues = array + numberNonZero;
#endif
    int n = block->firstBasic_;
    int nL = block->firstAtUpper_;
    nMax = static_cast< int >(newValues - array) + n;
    for (int j = 0; j < nL; j++) {
      double value2 = newValues[j];
      if (fabs(value2) > zeroTolerance) {
        int iSequence = column[j];
        double alpha;
        double oldValue;
        double value;
        alpha = value2;
        if (alpha > 0.0) {
          oldValue = reducedCost[iSequence];
          value = oldValue - tentativeTheta * alpha;
          if (value < dualT) {
            value = oldValue - upperTheta * alpha;
            if (value < dualT && alpha >= acceptablePivot) {
              upperTheta = (oldValue - dualT) / alpha;
            }
            // add to list
            spareArray[numberRemaining] = value2;
            spareIndex[numberRemaining++] = iSequence;
          }
        }
        array[numberNonZero] = value2;
        index[numberNonZero++] = iSequence;
      }
    }
    for (int j = nL; j < n; j++) {
      double value2 = newValues[j];
      if (fabs(value2) > zeroTolerance) {
        int iSequence = column[j];
        double alpha;
        double oldValue;
        double value;
        alpha = -value2;
        if (alpha > 0.0) {
          oldValue = -reducedCost[iSequence];
          value = oldValue - tentativeTheta * alpha;
          if (value < dualT) {
            value = oldValue - upperTheta * alpha;
            if (value < dualT && alpha >= acceptablePivot) {
              upperTheta = (oldValue - dualT) / alpha;
            }
            // add to list
            spareArray[numberRemaining] = value2;
            spareIndex[numberRemaining++] = iSequence;
          }
        }
        array[numberNonZero] = value2;
        index[numberNonZero++] = iSequence;
      }
    }
  }
  for (int j = numberNonZero; j < nMax; j++)
    array[j] = 0.0;
  output->setNumElements(numberNonZero);
  candidate->setNumElements(numberRemaining);
  model->spareDoubleArray_[0] = upperTheta;
}
static void
transposeTimes3Bit2Odd(clpTempInfo &info)
{
  double zeroTolerance = info.tolerance;
  double dualTolerance = -info.dualTolerance; // sign swapped
  double *COIN_RESTRICT reducedCost = info.reducedCost;
  double *COIN_RESTRICT weights = info.solution;
  const unsigned int *reference = reinterpret_cast< const unsigned int * >(info.upper);
  const unsigned char *status = info.status;
  const CoinBigIndex *starts = info.start;
  const int *row = info.row;
  const int *column = info.which;
  const double *element = info.element;
  double scaleFactor = info.theta;
  const double *pi = info.cost;
  const double *piWeight = info.lower;
  double referenceIn = info.upperTheta;
  double devex = info.changeObj;
  assert(scaleFactor);
  int bestSequence = info.numberAdded;
  double bestRatio = info.bestPossible;
#define NO_CHANGE_MULTIPLIER 1
#define INVERSE_MULTIPLIER 1.0 / NO_CHANGE_MULTIPLIER
#if NO_CHANGE_MULTIPLIER != 1
  double bestRatio2 = NO_CHANGE_MULTIPLIER * bestRatio;
#else
#define bestRatio2 bestRatio
#endif
  CoinBigIndex end = starts[0];
  int numberOdd = info.numberToDo;
  for (int jColumn = 0; jColumn < numberOdd; jColumn++) {
    CoinBigIndex start = end;
    CoinBigIndex j;
    int iColumn = column[jColumn];
    end = starts[jColumn + 1];
    double value = 0.0;
    if ((status[iColumn] & 7) != 1) {
      for (j = start; j < end; j++) {
        int iRow = row[j];
        value -= pi[iRow] * element[j];
      }
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (j = start; j < end; j++) {
          int iRow = row[j];
          modification += piWeight[iRow] * element[j];
        }
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = CoinMax(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = CoinMax(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        value = reducedCost[iColumn] - value;
        reducedCost[iColumn] = value;
        unsigned char thisStatus = status[iColumn] & 7;
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > -FREE_ACCEPT * dualTolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < dualTolerance) {
          value *= value;
          if (value > bestRatio * weights[iColumn]) {
            bestSequence = iColumn;
            bestRatio = value / weights[iColumn];
#if NO_CHANGE_MULTIPLIER != 1
            bestRatio2 = bestRatio * NO_CHANGE_MULTIPLIER;
#endif
          }
        }
      } else {
        value = reducedCost[iColumn];
        unsigned char thisStatus = status[iColumn] & 7;
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > -FREE_ACCEPT * dualTolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < dualTolerance) {
          value *= value;
          if (value > bestRatio2 * weights[iColumn]) {
            bestSequence = iColumn;
            bestRatio2 = value / weights[iColumn];
#if NO_CHANGE_MULTIPLIER != 1
            bestRatio = bestRatio2 * INVERSE_MULTIPLIER;
#endif
          }
        }
      }
    }
  }
  info.numberAdded = bestSequence;
}

static void
transposeTimes3Bit2(clpTempInfo &info)
{
  double zeroTolerance = info.tolerance;
  double dualTolerance = -info.dualTolerance; // sign swapped
  double *COIN_RESTRICT reducedCost = info.reducedCost;
  double *COIN_RESTRICT weights = info.solution;
  double *COIN_RESTRICT work = info.work;
  double *COIN_RESTRICT work2 = work + COIN_AVX2_CHUNK;
  const unsigned int *reference = reinterpret_cast< const unsigned int * >(info.upper);
  const blockStruct *blocks = reinterpret_cast< const blockStruct * >(info.pivotVariable);
  const unsigned char *status = info.status;
  const int *rowBlock = info.row;
  const int *columnBlock = info.which;
  const double *elementBlock = info.element;
  double scaleFactor = info.theta;
  const double *pi = info.cost;
  const double *piWeight = info.lower;
  double referenceIn = info.upperTheta;
  double devex = info.changeObj;
  assert(scaleFactor);
  int bestSequence = info.numberAdded;
  double bestRatio = info.bestPossible;
#if NO_CHANGE_MULTIPLIER != 1
  double bestRatio2 = NO_CHANGE_MULTIPLIER * bestRatio;
#endif
#define INCLUDE_MATRIX3_PRICING 1
  int firstBlock = info.startColumn;
  int lastBlock = info.numberToDo;
  //double * COIN_RESTRICT tempArray =
  //const_cast<double *>(elementBlock+info.numberColumns);
  //double * COIN_RESTRICT tempArray2 = tempArray + roundUp(maxBlockSize);
  for (int iBlock = firstBlock; iBlock < lastBlock; iBlock++) {
    const blockStruct *block = blocks + iBlock;
    int numberPrice = block->firstBasic_;
    int nel = block->numberElements_;
    const int *row = rowBlock + block->startElements_;
    const double *element = elementBlock + block->startElements_;
    const int *column = columnBlock + block->startIndices_;
#if COIN_AVX2 > 1
    int numberDo = roundDown(numberPrice);
#if defined(HASWELL) || defined(SKYLAKE)
#ifdef HASWELL
    assert(COIN_AVX2 == 4);
    __m256d zero = _mm256_setzero_pd();
    for (int kColumn = 0; kColumn < numberDo; kColumn += COIN_AVX2_CHUNK) {
      int endBlock = CoinMin(COIN_AVX2_CHUNK, numberDo - kColumn);
      for (int i = 0; i < endBlock; i += COIN_AVX2) {
        __m256d arrayX = _mm256_setzero_pd();
        __m256d tempX = _mm256_setzero_pd();
        for (int j = 0; j < nel; j++) {
          __m128i rows = _mm_loadu_si128((const __m128i *)row);
          __m256d elements = _mm256_load_pd(element);
          __m256d pis = _mm256_i32gather_pd(pi, rows, 8);
          __m256d piWeights = _mm256_i32gather_pd(piWeight, rows, 8);
          // should be better way using sub
          arrayX = _mm256_fmadd_pd(pis, elements, arrayX);
          tempX = _mm256_fmadd_pd(piWeights, elements, tempX);
          row += 4;
          element += 4;
        }
        arrayX = _mm256_sub_pd(zero, arrayX);
        _mm256_store_pd(work + i, tempX);
        _mm256_store_pd(work2 + i, arrayX);
      }
#else
    assert(COIN_AVX2 == 8);
    __m512d zero = _mm512_setzero_pd();
    for (int kColumn = 0; kColumn < numberDo; kColumn += COIN_AVX2_CHUNK) {
      int endBlock = CoinMin(COIN_AVX2_CHUNK, numberDo - kColumn);
      for (int i = 0; i < endBlock; i += COIN_AVX2) {
        __m512d arrayX = _mm512_setzero_pd();
        __m512d tempX = _mm512_setzero_pd();
        for (int j = 0; j < nel; j++) {
          __m256i rows = _mm256_loadu_si256((const __m256i *)row);
          __m512d elements = _mm512_load_pd(element);
          __m512d pis = _mm512_i32gather_pd(rows, pi, 8);
          __m512d piWeights = _mm512_i32gather_pd(rows, piWeight, 8);
          // should be better way using sub
          arrayX = _mm512_fmadd_pd(pis, elements, arrayX);
          tempX = _mm512_fmadd_pd(piWeights, elements, tempX);
          row += 4;
          element += 4;
        }
        arrayX = _mm512_sub_pd(zero, arrayX);
        _mm512_store_pd(work + i, tempX);
        _mm512_store_pd(work2 + i, arrayX);
      }
#endif
      for (int i = 0; i < endBlock; i++) {
        double value = work2[i];
        double modification = work[i];
        // common coding
#include "ClpPackedMatrix.hpp"
      }
    }
#else
    for (int kColumn = 0; kColumn < numberDo; kColumn += COIN_AVX2_CHUNK) {
      int n = 0;
      int nBlock = CoinMin(COIN_AVX2_CHUNK, numberPrice - kColumn) >> COIN_AVX2_SHIFT;
      //int nBlock = numberPrice>>COIN_AVX2_SHIFT;
      for (int jBlock = 0; jBlock < nBlock; jBlock++) {
        for (int j = 0; j < COIN_AVX2; j++) {
          double value = 0.0;
          double modification = 0.0;
          for (int i = 0; i < nel; i++) {
            //assert (element+i*COIN_AVX2<tempArray);
            int iRow = row[i * COIN_AVX2];
            value -= pi[iRow] * element[i * COIN_AVX2];
            modification += piWeight[iRow] * element[i * COIN_AVX2];
          }
          work[n] = modification;
          work2[n++] = value;
          row++;
          element++;
        }
        row += (nel - 1) * COIN_AVX2;
        element += (nel - 1) * COIN_AVX2;
        //assert (row == rowBlock + block->startElements_ + nel*COIN_AVX2*(jBlock+1));
      }
      for (int j = 0; j < n; j++) {
        double value = work2[j];
        double modification = work[j];
        // common coding
#include "ClpPackedMatrix.hpp"
      }
    }
#endif
    // last lot
    for (int j = numberDo; j < numberPrice; j++) {
      double value = 0.0;
      double modification = 0.0;
      for (int i = 0; i < nel; i++) {
        int iRow = row[i * COIN_AVX2];
        value -= pi[iRow] * element[i * COIN_AVX2];
        modification += piWeight[iRow] * element[i * COIN_AVX2];
      }
// common coding
#include "ClpPackedMatrix.hpp"
      row++;
      element++;
    }
#else // COIN_AVX2 == 1
    for (int j = 0; j < numberPrice; j++) {
      double value = 0.0;
      double modification = 0.0;
      for (int i = 0; i < nel; i++) {
        //assert (element+i*COIN_AVX2<tempArray);
        int iRow = row[i * COIN_AVX2];
        value -= pi[iRow] * element[i * COIN_AVX2];
        modification += piWeight[iRow] * element[i * COIN_AVX2];
      }
// common coding
#include "ClpPackedMatrix.hpp"
      row += nel;
      element += nel;
    }
    int n = numberPrice;
#endif
  }
  info.numberAdded = bestSequence;
  info.bestPossible = bestRatio;
}
static void
transposeTimes3BitSlacks(clpTempInfo &info)
{
  double dualTolerance = info.dualTolerance;
  double *COIN_RESTRICT reducedCost = info.reducedCost;
  double *COIN_RESTRICT weights = info.solution;
  const int *columnBlock = info.which;
  const blockStruct *blocks = reinterpret_cast< const blockStruct * >(info.pivotVariable);
  //const unsigned char * status = info.status;
  int bestSequence = info.numberAdded;
  double bestRatio = info.bestPossible;
  int iBlock = info.startColumn;
  assert(info.numberToDo == iBlock + 1);
  const blockStruct *block = blocks + iBlock;
  int first = 0;
  int last = block->firstAtLower_;
  const int *column = columnBlock + block->startIndices_;
  for (int j = first; j < last; j++) {
    int iColumn = *column;
    column++;
    double value = reducedCost[iColumn];
    // free or superbasic
    if (fabs(value) > FREE_ACCEPT * dualTolerance) {
      // we are going to bias towards free (but only if reasonable)
      value = -fabs(value) * FREE_BIAS;
      value *= value;
      if (value > bestRatio * weights[iColumn]) {
        bestSequence = iColumn;
        bestRatio = value / weights[iColumn];
      }
    }
  }
  first = last;
  last = block->firstAtUpper_;
  dualTolerance = -dualTolerance; //flip sign
  for (int j = first; j < last; j++) {
    int iColumn = *column;
    column++;
    double value = reducedCost[iColumn];
    // at lower
    if (value < dualTolerance) {
      value *= value;
      if (value > bestRatio * weights[iColumn]) {
        bestSequence = iColumn;
        bestRatio = value / weights[iColumn];
      }
    }
  }
  first = last;
  last = block->firstBasic_;
  dualTolerance = -dualTolerance; //flip sign
  for (int j = first; j < last; j++) {
    int iColumn = *column;
    column++;
    double value = reducedCost[iColumn];
    // at upper
    if (value > dualTolerance) {
      value *= value;
      if (value > bestRatio * weights[iColumn]) {
        bestSequence = iColumn;
        bestRatio = value / weights[iColumn];
      }
    }
  }
  info.numberAdded = bestSequence;
  info.bestPossible = bestRatio;
}
// Updates two arrays for steepest
void ClpPackedMatrix3::transposeTimes2(const ClpSimplex *model,
  const double *COIN_RESTRICT pi,
  CoinIndexedVector *output,
  const double *COIN_RESTRICT piWeight,
  double *COIN_RESTRICT infeas,
  double *COIN_RESTRICT reducedCost,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  /*
    This is designed for meaty problems so probably should also go to exact devex
    also could do sequenceIn_ choice here
    could do ones at lb separately - ?? masking - ? masking for weights
   */
  double zeroTolerance = model->zeroTolerance();
  double dualTolerance = model->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = CoinMin(1.0e-2, model->largestDualError());
  // allow tolerance at least slightly bigger than standard
  dualTolerance = dualTolerance + error;
  int numberOdd = block_->startIndices_;
  assert(scaleFactor);
#if ABOCA_LITE
  //int nChunks=numberChunks_;
#define ODD_INFO 2 * ABOCA_LITE
#else
#define ODD_INFO 1
  //int nChunks=1;
#endif
  clpTempInfo info[ODD_INFO + 2] = {};
  unsigned char *status = model->statusArray();
  // fill in info
  for (int i = 0; i < ODD_INFO + 2; i++) {
    //memset(info+i,0,sizeof(clpTempInfo));
    info[i].tolerance = zeroTolerance;
    info[i].dualTolerance = dualTolerance;
    info[i].reducedCost = reducedCost;
    info[i].infeas = infeas;
    info[i].solution = weights;
    info[i].upper = reinterpret_cast< const double * >(reference);
    info[i].pivotVariable = reinterpret_cast< const int * >(block_);
    info[i].status = status;
    info[i].start = start_;
    info[i].row = row_;
    info[i].which = column_;
    info[i].element = element_;
    info[i].theta = scaleFactor;
    info[i].cost = pi;
    info[i].lower = piWeight;
    info[i].upperTheta = referenceIn;
    info[i].changeObj = devex;
    info[i].numberAdded = -1;
  }
  //assert (ODD_INFO==1);
  info[ODD_INFO].startColumn = 0;
  info[ODD_INFO].numberToDo = numberOdd;
  //
  double *temporary = temporary_->array();
#if ABOCA_LITE
  for (int iBlock = 0; iBlock < numberChunks_; iBlock++) {
    info[iBlock].startColumn = endChunk_[iBlock];
    info[iBlock].numberToDo = endChunk_[iBlock + 1];
    info[iBlock].work = temporary + 2 * COIN_AVX2_CHUNK * iBlock;
  }
#else
  info[0].startColumn = 0;
  info[0].numberToDo = numberBlocks_;
  info[0].work = temporary;
#endif
  info[ODD_INFO + 1].startColumn = numberBlocks_;
  info[ODD_INFO + 1].numberToDo = numberBlocks_ + 1;
#if ABOCA_LITE
  if (abcState() > 1) {
    cilk_spawn transposeTimes3Bit2Odd(info[ODD_INFO]);
    for (int iBlock = 0; iBlock < numberChunks_; iBlock++) {
      cilk_spawn transposeTimes3Bit2(info[iBlock]);
    }
    transposeTimes3BitSlacks(info[ODD_INFO + 1]);
    cilk_sync;
  } else {
#endif
    transposeTimes3Bit2Odd(info[ODD_INFO]);
#if ABOCA_LITE
    for (int iBlock = 0; iBlock < numberChunks_; iBlock++) {
      transposeTimes3Bit2(info[iBlock]);
    }
#else
  transposeTimes3Bit2(info[0]);
#endif
    transposeTimes3BitSlacks(info[ODD_INFO + 1]);
#if ABOCA_LITE
  }
#endif
  int sequenceIn = -1;
  double best = 0.0;
  for (int iBlock = 0; iBlock < ODD_INFO + 2; iBlock++) {
    if (info[iBlock].bestPossible > best) {
      best = info[iBlock].bestPossible;
      sequenceIn = info[iBlock].numberAdded;
    }
  }
  // Make sure sequenceOut not eligible
  int sequenceOut = model->sequenceOut();
  double saveOutDj = 0.0;
  if (sequenceOut >= 0) {
    saveOutDj = reducedCost[sequenceOut];
    unsigned char thisStatus = status[sequenceOut] & 7;
    assert(thisStatus != 0 && thisStatus != 4);
    if (thisStatus != 2)
      reducedCost[sequenceOut] = COIN_DBL_MAX;
    else
      reducedCost[sequenceOut] = -COIN_DBL_MAX;
  }
  if (sequenceIn >= 0) {
    // check not flagged
    if (model->flagged(sequenceIn) || sequenceIn == sequenceOut) {
      int number = model->numberRows() + model->numberColumns();
      const unsigned char *COIN_RESTRICT status = model->statusArray();
      sequenceIn = -2;
      double bestRatio = 0.0;
      for (int iSequence = 0; iSequence < number; iSequence++) {
        unsigned char thisStatus = status[iSequence] & 7;
        double value = reducedCost[iSequence];
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > FREE_ACCEPT * -dualTolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < dualTolerance) {
          value *= value;
          if (value > bestRatio * weights[iSequence]) {
            if (!model->flagged(iSequence)) {
              sequenceIn = iSequence;
              bestRatio = value / weights[iSequence];
            }
          }
        }
      }
    }
  }
  if (sequenceOut >= 0)
    reducedCost[sequenceOut] = saveOutDj;
  model->spareIntArray_[3] = sequenceIn;
}
/*
  type - 1 redo infeasible, 2 choose sequenceIn, 3 both
  returns sequenceIn (or -1) for type 2
 */
int ClpPackedMatrix3::redoInfeasibilities(const ClpSimplex *model,
  ClpPrimalColumnSteepest *pivotChoose,
  int type)
{
  CoinIndexedVector *infeasible = pivotChoose->infeasible();
  // change to use blocks
  double tolerance = model->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = CoinMin(1.0e-2, model->largestDualError());
  // allow tolerance at least slightly bigger than standard
  tolerance = tolerance + error;
  // reverse sign so test is cleaner
  tolerance = -tolerance;
  int number = model->numberRows() + model->numberColumns();

  const double *COIN_RESTRICT reducedCost = model->djRegion();
  const unsigned char *COIN_RESTRICT status = model->statusArray();
  const double *COIN_RESTRICT weights = pivotChoose->weights();
  int sequenceIn = -1;
  double bestRatio = 0.0;
  switch (type) {
  case 1:
    infeasible->clear();
    for (int iSequence = 0; iSequence < number; iSequence++) {
      unsigned char thisStatus = status[iSequence] & 7;
      double value = reducedCost[iSequence];
      if (thisStatus == 3) {
      } else if ((thisStatus & 1) != 0) {
        // basic or fixed
        value = 0.0;
      } else if (thisStatus == 2) {
        value = -value;
      } else {
        // free or superbasic
        if (fabs(value) > FREE_ACCEPT * -tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value = -fabs(value) * FREE_BIAS;
        } else {
          value = 0.0;
        }
      }
      if (value < tolerance) {
        // store square in list
        infeasible->quickAdd(iSequence, value * value);
      }
    }
    break;
  case 2:
    infeasible->clear(); // temp for debug
    for (int iSequence = 0; iSequence < number; iSequence++) {
      unsigned char thisStatus = status[iSequence] & 7;
      double value = reducedCost[iSequence];
      if (thisStatus == 3) {
      } else if ((thisStatus & 1) != 0) {
        // basic or fixed
        value = 0.0;
      } else if (thisStatus == 2) {
        value = -value;
      } else {
        // free or superbasic
        if (fabs(value) > FREE_ACCEPT * -tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value = -fabs(value) * FREE_BIAS;
        } else {
          value = 0.0;
        }
      }
      if (value < tolerance) {
        value *= value;
        if (value > bestRatio * weights[iSequence]) {
          sequenceIn = iSequence;
          bestRatio = value / weights[iSequence];
        }
      }
    }
    break;
  case 3:
    infeasible->clear();
    for (int iSequence = 0; iSequence < number; iSequence++) {
      unsigned char thisStatus = status[iSequence] & 7;
      double value = reducedCost[iSequence];
      if (thisStatus == 3) {
      } else if ((thisStatus & 1) != 0) {
        // basic or fixed
        value = 0.0;
      } else if (thisStatus == 2) {
        value = -value;
      } else {
        // free or superbasic
        if (fabs(value) > FREE_ACCEPT * -tolerance) {
          // we are going to bias towards free (but only if reasonable)
          value = -fabs(value) * FREE_BIAS;
        } else {
          value = 0.0;
        }
      }
      if (value < tolerance) {
        value *= value;
        // store square in list
        infeasible->quickAdd(iSequence, value);
        if (value > bestRatio * weights[iSequence]) {
          sequenceIn = iSequence;
          bestRatio = value / weights[iSequence];
        }
      }
    }
    break;
  }
  if (sequenceIn >= 0) {
    // check not flagged
    if (model->flagged(sequenceIn)) {
      sequenceIn = -1;
      for (int iSequence = 0; iSequence < number; iSequence++) {
        unsigned char thisStatus = status[iSequence] & 7;
        double value = reducedCost[iSequence];
        if (thisStatus == 3) {
        } else if ((thisStatus & 1) != 0) {
          // basic or fixed
          value = 0.0;
        } else if (thisStatus == 2) {
          value = -value;
        } else {
          // free or superbasic
          if (fabs(value) > FREE_ACCEPT * -tolerance) {
            // we are going to bias towards free (but only if reasonable)
            value = -fabs(value) * FREE_BIAS;
          } else {
            value = 0.0;
          }
        }
        if (value < tolerance) {
          value *= value;
          if (value > bestRatio * weights[iSequence]) {
            if (!model->flagged(iSequence)) {
              sequenceIn = iSequence;
              bestRatio = value / weights[iSequence];
            }
          }
        }
      }
    }
  }
  return sequenceIn;
}
#if COIN_LONG_WORK
// For long double versions
void ClpPackedMatrix::times(CoinWorkDouble scalar,
  const CoinWorkDouble *x, CoinWorkDouble *y) const
{
  int iRow, iColumn;
  // get matrix data pointers
  const int *row = matrix_->getIndices();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const double *elementByColumn = matrix_->getElements();
  //memset(y,0,matrix_->getNumRows()*sizeof(double));
  assert(((flags_ & 2) != 0) == matrix_->hasGaps());
  if (!(flags_ & 2)) {
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      CoinWorkDouble value = x[iColumn];
      if (value) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = columnStart[iColumn + 1];
        value *= scalar;
        for (j = start; j < end; j++) {
          iRow = row[j];
          y[iRow] += value * elementByColumn[j];
        }
      }
    }
  } else {
    const int *columnLength = matrix_->getVectorLengths();
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      CoinWorkDouble value = x[iColumn];
      if (value) {
        CoinBigIndex start = columnStart[iColumn];
        CoinBigIndex end = start + columnLength[iColumn];
        value *= scalar;
        for (j = start; j < end; j++) {
          iRow = row[j];
          y[iRow] += value * elementByColumn[j];
        }
      }
    }
  }
}
void ClpPackedMatrix::transposeTimes(CoinWorkDouble scalar,
  const CoinWorkDouble *x, CoinWorkDouble *y) const
{
  int iColumn;
  // get matrix data pointers
  const int *row = matrix_->getIndices();
  const CoinBigIndex *columnStart = matrix_->getVectorStarts();
  const double *elementByColumn = matrix_->getElements();
  if (!(flags_ & 2)) {
    if (scalar == -1.0) {
      CoinBigIndex start = columnStart[0];
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        CoinBigIndex j;
        CoinBigIndex next = columnStart[iColumn + 1];
        CoinWorkDouble value = y[iColumn];
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value -= x[jRow] * elementByColumn[j];
        }
        start = next;
        y[iColumn] = value;
      }
    } else {
      CoinBigIndex start = columnStart[0];
      for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
        CoinBigIndex j;
        CoinBigIndex next = columnStart[iColumn + 1];
        CoinWorkDouble value = 0.0;
        for (j = start; j < next; j++) {
          int jRow = row[j];
          value += x[jRow] * elementByColumn[j];
        }
        start = next;
        y[iColumn] += value * scalar;
      }
    }
  } else {
    const int *columnLength = matrix_->getVectorLengths();
    for (iColumn = 0; iColumn < numberActiveColumns_; iColumn++) {
      CoinBigIndex j;
      CoinWorkDouble value = 0.0;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      for (j = start; j < end; j++) {
        int jRow = row[j];
        value += x[jRow] * elementByColumn[j];
      }
      y[iColumn] += value * scalar;
    }
  }
}
#endif
#ifdef CLP_ALL_ONE_FILE
#undef reference
#endif
#if 0
#undef ABCSTATE_LITE
#if ABOCA_LITE
/* Meat of transposeTimes by column when not scaled and skipping
   and doing part of dualColumn */
static void
dualColumn00(clpTempInfo & info)
{
  const int * COIN_RESTRICT which = info.which;
  const double * COIN_RESTRICT work = info.work;
  int * COIN_RESTRICT index = info.index;
  double * COIN_RESTRICT spare = info.spare;
  const unsigned char * COIN_RESTRICT status = info.status;
  const double * COIN_RESTRICT reducedCost = info.reducedCost;
  double upperTheta = info.upperTheta;
  double acceptablePivot = info.acceptablePivot;
  double dualTolerance = info.tolerance;
  int numberToDo=info.numberToDo;
  double tentativeTheta = 1.0e15;
  int numberRemaining = 0;
  double multiplier[] = { -1.0, 1.0};
  double dualT = - dualTolerance;
  for (int i = 0; i < numberToDo; i++) {
    int iSequence = which[i];
    int wanted = (status[iSequence] & 3) - 1;
    if (wanted) {
      double mult = multiplier[wanted-1];
      double alpha = work[i] * mult;
      if (alpha > 0.0) {
	double oldValue = reducedCost[iSequence] * mult;
	double value = oldValue - tentativeTheta * alpha;
	if (value < dualT) {
	  value = oldValue - upperTheta * alpha;
	  if (value < dualT && alpha >= acceptablePivot) {
	    upperTheta = (oldValue - dualT) / alpha;
	  }
	  // add to list
	  spare[numberRemaining] = alpha * mult;
	  index[numberRemaining++] = iSequence;
	}
      }
    }
  }
  info.numberRemaining = numberRemaining;
  info.upperTheta = upperTheta;
}
#endif
int
ClpSimplexDual::dualColumn0(const CoinIndexedVector * rowArray,
                            const CoinIndexedVector * columnArray,
                            CoinIndexedVector * spareArray,
                            double acceptablePivot,
                            double & upperReturn, double &bestReturn, double & badFree)
{
     // do first pass to get possibles
     double * spare = spareArray->denseVector();
     int * index = spareArray->getIndices();
     const double * work;
     int number;
     const int * which;
     const double * reducedCost;
     // We can also see if infeasible or pivoting on free
     double tentativeTheta = 1.0e15;
     double upperTheta = 1.0e31;
     double freePivot = acceptablePivot;
     double bestPossible = 0.0;
     int numberRemaining = 0;
     int i;
     badFree = 0.0;
     if ((moreSpecialOptions_ & 8) != 0) {
          // No free or super basic
       // bestPossible will re recomputed if necessary
       bestPossible=1.0;
       //#define COIN_AVX7
#ifndef COIN_AVX2
          double multiplier[] = { -1.0, 1.0};
#else
	  double multiplier[4]={0.0,0.0,-1.0,1.0};
#endif
          double dualT = - dualTolerance_;
#if ABOCA_LITE == 0
	  int nSections=2;
#else
	  int numberThreads=abcState();
	  int nSections=numberThreads ? 1 : 2;
#endif
          for (int iSection = 0; iSection < nSections; iSection++) {

               int addSequence;
               unsigned char * statusArray;
               if (!iSection) {
                    work = rowArray->denseVector();
                    number = rowArray->getNumElements();
                    which = rowArray->getIndices();
                    reducedCost = rowReducedCost_;
                    addSequence = numberColumns_;
                    statusArray = status_ + numberColumns_;
               } else {
                    work = columnArray->denseVector();
                    number = columnArray->getNumElements();
                    which = columnArray->getIndices();
                    reducedCost = reducedCostWork_;
                    addSequence = 0;
                    statusArray = status_;
               }
#ifndef COIN_AVX2
               for (i = 0; i < number; i++) {
                    int iSequence = which[i];
                    double alpha;
                    double oldValue;
                    double value;

                    assert (getStatus(iSequence + addSequence) != isFree
                            && getStatus(iSequence + addSequence) != superBasic);
                    int iStatus = (statusArray[iSequence] & 3) - 1;
                    if (iStatus) {
                         double mult = multiplier[iStatus-1];
                         alpha = work[i] * mult;
                         if (alpha > 0.0) {
                              oldValue = reducedCost[iSequence] * mult;
                              value = oldValue - tentativeTheta * alpha;
                              if (value < dualT) {
                                   value = oldValue - upperTheta * alpha;
                                   if (value < dualT && alpha >= acceptablePivot) {
                                        upperTheta = (oldValue - dualT) / alpha;
                                        //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
                                   }
                                   // add to list
                                   spare[numberRemaining] = alpha * mult;
                                   index[numberRemaining++] = iSequence + addSequence;
                              }
                         }
                    }
               }
	       //
#else
	       //#define COIN_AVX2 4 // temp
#if COIN_AVX2 == 1
#define COIN_AVX2_SHIFT 0
#elif COIN_AVX2 == 2
#define COIN_AVX2_SHIFT 1
#elif COIN_AVX2 == 4
#define COIN_AVX2_SHIFT 2
#elif COIN_AVX2 == 8
#define COIN_AVX2_SHIFT 3
#else
  error;
#endif
  //#define COIN_ALIGN 8*COIN_AVX2 // later
  //#define COIN_ALIGN_DOUBLE COIN_AVX2
#define CHECK_CHUNK 4
	       // round up
	       int * whichX = const_cast<int *>(which);
	       int nBlocks = (number+CHECK_CHUNK-1)/CHECK_CHUNK;
	       int n=nBlocks*CHECK_CHUNK+1;
	       for (int i=number;i<n;i++)
		 whichX[i]=0; // alpha will be zero so not chosen
	       bool acceptableX[CHECK_CHUNK+1];
	       double oldValueX[CHECK_CHUNK+1];
	       double newValueX[CHECK_CHUNK+1];
	       double alphaX[CHECK_CHUNK+1];
	       newValueX[CHECK_CHUNK]=0.0;
#define USE_USE_AVX
	       //#define CHECK_H 1
#ifdef USE_USE_AVX
#define NEED_AVX
#elif CHECK_H
#define NEED_AVX
#endif
#ifdef NEED_AVX
	       double mult2[CHECK_CHUNK];
	       CoinInt64 acceptableY[CHECK_CHUNK];
	       CoinInt64 goodDj[CHECK_CHUNK];
	       double oldValueY[CHECK_CHUNK];
	       double alphaY[CHECK_CHUNK+1];
	       memset(acceptableY,0,sizeof(acceptableY));
	       memset(goodDj,0,sizeof(goodDj));
	       memset(oldValueY,0,sizeof(oldValueY));
	       memset(alphaY,0,sizeof(alphaY));
	       __m256d tentative2 = _mm256_set1_pd(-tentativeTheta);
	       __m256d dualT2 = _mm256_set1_pd(dualT);
	       __m256d acceptable2 = _mm256_set1_pd(acceptablePivot);
#endif
	       for (int iBlock=0;iBlock<nBlocks;iBlock++) {
		 bool store=false;
		 double alpha=0.0;
		 double oldValue=0.0;
		 double newValue=0.0;
		 double trueAlpha=0.0;
		 int jSequence=0;
#ifndef USE_USE_AVX
		 for (int i = 0; i < CHECK_CHUNK+1; i++) {
		   int iSequence = which[i];
		   int iStatus = (statusArray[iSequence] & 3);
		   double mult = multiplier[iStatus];
		   double newAlpha = work[i]*mult;
		   double oldDj = reducedCost[iSequence]*mult;
		   newValue = (oldDj - tentativeTheta * newAlpha)-dualT;
		   acceptableX[i]=newAlpha>= acceptablePivot;
		   oldValueX[i]=oldDj;
		   newValueX[i]=newValue;
		   alphaX[i]=newAlpha;
		 }
#endif
#ifdef NEED_AVX
		 __m128i columns = _mm_load_si128((const __m128i *)which);
		 // what do we get - this must be wrong
		 // probably only 1 and 2 - can we be clever
		 // fix
		 //__m128i status; // = _mm256_i32gather_ps(statusArray,columns,1);
		 //status.m128i_i32[0]=statusArray[columns.m128i_i32[0]];
		 for (int i=0;i<CHECK_CHUNK;i++) {
		   int iSequence = which[i];
		   int iStatus = (statusArray[iSequence] & 3);
		   mult2[i] = multiplier[iStatus];
		 }
		 //__m256d newAlpha2 = _mm256_i32gather_pd(multiplier,status,1); // mult here
		 __m256d newAlpha2 = _mm256_load_pd(mult2);
		 __m256d oldDj = _mm256_i32gather_pd(reducedCost,columns,8);
		 oldDj=_mm256_mul_pd(oldDj,newAlpha2); // remember newAlpha==mult
		 _mm256_store_pd(oldValueY,oldDj); // redo later
		 __m256d work2 = _mm256_load_pd(work);
		 newAlpha2=_mm256_mul_pd(newAlpha2,work2); // now really newAlpha
		 //__m256d newValue2 = _mm256_fmadd_pd(tentative2,newAlpha2,oldDj);
		 oldDj=_mm256_fmadd_pd(newAlpha2,tentative2,oldDj);
		 __v4df bitsDj=_mm256_cmp_pd(oldDj,dualT2,_CMP_LT_OS);
		 __v4df bitsAcceptable = _mm256_cmp_pd(newAlpha2,acceptable2,_CMP_GE_OS);
		 _mm256_store_pd(reinterpret_cast<double *>(goodDj),bitsDj);
		 _mm256_store_pd(reinterpret_cast<double *>(acceptableY),bitsAcceptable);
		 _mm256_store_pd(alphaY,newAlpha2);
#ifndef USE_USE_AVX
#undef NDEBUG
		 for (int i = 0; i < CHECK_CHUNK; i++) {
		   assert(newValueX[i]>0.0==(goodDj[i]));
		   //assert(acceptableX[i]==(acceptableY[i]));
		   assert(oldValueX[i]==oldValueY[i]);
		   assert(alphaX[i]==alphaY[i]);
		 }
		 for (int i = 0; i < CHECK_CHUNK; i++) {
		   bool g1=newValueX[i]<0.0;
		   bool g2=goodDj[i]!=0;
		   if(g1!=g2)abort();
		   //if(acceptableX[i]!=(acceptableY[i]))abort();
		   if(fabs(oldValueX[i]-oldValueY[i])>1.0e-5+
		      +(1.0e-10*fabs(oldValueX[i])))abort();
		   if(alphaX[i]!=alphaY[i])abort();
		 }
#endif
#endif
		 for (int i = 0; i < CHECK_CHUNK+1; i++) {
#ifndef USE_USE_AVX
		   double newValue=newValueX[i];
		   bool newStore = newValue < 0.0;
		   if (store) {
		     // add to list
		     bool acceptable = acceptableX[i-1];
		     spare[numberRemaining] = work[i-1];
		     index[numberRemaining++] = which[i-1] + addSequence;
		     double value = oldValueX[i-1] - upperTheta * alphaX[i-1];
		     if (value < dualT && acceptable) {
		       upperTheta = (oldValueX[i-1] - dualT) / alphaX[i-1];
		     }
		   }
#else
		   bool newStore = goodDj[i]!=0;
		   if (store) {
		     // add to list
		     bool acceptable = acceptableY[i-1];
		     spare[numberRemaining] = work[i-1];
		     index[numberRemaining++] = which[i-1] + addSequence;
		     double value = oldValueY[i-1] - upperTheta * alphaY[i-1];
		     if (value < dualT && acceptable) {
		       upperTheta = (oldValueY[i-1] - dualT) / alphaY[i-1];
		     }
		   }
#endif
		   store=newStore;
		 }
		 which += CHECK_CHUNK;
		 work += CHECK_CHUNK;
	       }
#endif
          }
#if ABOCA_LITE
	  if (numberThreads) {
	  work = columnArray->denseVector();
	  number = columnArray->getNumElements();
	  which = columnArray->getIndices();
	  reducedCost = reducedCostWork_;
	  unsigned char * statusArray = status_;
	  
	  clpTempInfo info[ABOCA_LITE];
	  int chunk = (number+numberThreads-1)/numberThreads;
	  int n=0;
	  int nR=numberRemaining;
	  for (int i=0;i<numberThreads;i++) {
	    info[i].which=const_cast<int *>(which+n);
	    info[i].work=const_cast<double *>(work+n);
	    info[i].numberToDo=CoinMin(chunk,number-n);
	    n += chunk;
	    info[i].index = index+nR;
	    info[i].spare = spare+nR;
	    nR += chunk;
	    info[i].reducedCost = const_cast<double *>(reducedCost);
	    info[i].upperTheta = upperTheta;
	    info[i].acceptablePivot = acceptablePivot;
	    info[i].status = statusArray;
	    info[i].tolerance=dualTolerance_;
	  }
	  for (int i=0;i<numberThreads;i++) {
	    cilk_spawn dualColumn00(info[i]);
	  }
	  cilk_sync;
	  moveAndZero(info,1,NULL);
	  for (int i=0;i<numberThreads;i++) {
	    numberRemaining += info[i].numberRemaining;
	    upperTheta = CoinMin(upperTheta,static_cast<double>(info[i].upperTheta));
	  }
	  }
#endif
     } else {
          // some free or super basic
          for (int iSection = 0; iSection < 2; iSection++) {

               int addSequence;

               if (!iSection) {
                    work = rowArray->denseVector();
                    number = rowArray->getNumElements();
                    which = rowArray->getIndices();
                    reducedCost = rowReducedCost_;
                    addSequence = numberColumns_;
               } else {
                    work = columnArray->denseVector();
                    number = columnArray->getNumElements();
                    which = columnArray->getIndices();
                    reducedCost = reducedCostWork_;
                    addSequence = 0;
               }

               for (i = 0; i < number; i++) {
                    int iSequence = which[i];
                    double alpha;
                    double oldValue;
                    double value;
                    bool keep;

                    switch(getStatus(iSequence + addSequence)) {

                    case basic:
                    case ClpSimplex::isFixed:
                         break;
                    case isFree:
                    case superBasic:
                         alpha = work[i];
                         bestPossible = CoinMax(bestPossible, fabs(alpha));
                         oldValue = reducedCost[iSequence];
                         // If free has to be very large - should come in via dualRow
                         //if (getStatus(iSequence+addSequence)==isFree&&fabs(alpha)<1.0e-3)
                         //break;
                         if (oldValue > dualTolerance_) {
                              keep = true;
                         } else if (oldValue < -dualTolerance_) {
                              keep = true;
                         } else {
                              if (fabs(alpha) > CoinMax(10.0 * acceptablePivot, 1.0e-5)) {
                                   keep = true;
                              } else {
                                   keep = false;
                                   badFree = CoinMax(badFree, fabs(alpha));
                              }
                         }
                         if (keep) {
                              // free - choose largest
                              if (fabs(alpha) > freePivot) {
                                   freePivot = fabs(alpha);
                                   sequenceIn_ = iSequence + addSequence;
                                   theta_ = oldValue / alpha;
                                   alpha_ = alpha;
                              }
			      // give fake bounds if possible
			      int jSequence=iSequence+addSequence;
			      if (2.0*fabs(solution_[jSequence])<
				  dualBound_) {
				FakeBound bound = getFakeBound(jSequence);
				assert (bound == ClpSimplexDual::noFake);
				setFakeBound(jSequence,ClpSimplexDual::bothFake);
				numberFake_++;
				value = oldValue - tentativeTheta * alpha;
				if (value > dualTolerance_) {
				  // pretend coming in from upper bound
				  upper_[jSequence] = solution_[jSequence];
				  lower_[jSequence] = upper_[jSequence] - dualBound_;
				  setColumnStatus(jSequence,ClpSimplex::atUpperBound);
				} else {
				  // pretend coming in from lower bound
				  lower_[jSequence] = solution_[jSequence];
				  upper_[jSequence] = lower_[jSequence] + dualBound_;
				  setColumnStatus(jSequence,ClpSimplex::atLowerBound);
				}
			      }
                         }
                         break;
                    case atUpperBound:
                         alpha = work[i];
                         oldValue = reducedCost[iSequence];
                         value = oldValue - tentativeTheta * alpha;
                         //assert (oldValue<=dualTolerance_*1.0001);
                         if (value > dualTolerance_) {
                              bestPossible = CoinMax(bestPossible, -alpha);
                              value = oldValue - upperTheta * alpha;
                              if (value > dualTolerance_ && -alpha >= acceptablePivot) {
                                   upperTheta = (oldValue - dualTolerance_) / alpha;
                                   //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
                              }
                              // add to list
                              spare[numberRemaining] = alpha;
                              index[numberRemaining++] = iSequence + addSequence;
                         }
                         break;
                    case atLowerBound:
                         alpha = work[i];
                         oldValue = reducedCost[iSequence];
                         value = oldValue - tentativeTheta * alpha;
                         //assert (oldValue>=-dualTolerance_*1.0001);
                         if (value < -dualTolerance_) {
                              bestPossible = CoinMax(bestPossible, alpha);
                              value = oldValue - upperTheta * alpha;
                              if (value < -dualTolerance_ && alpha >= acceptablePivot) {
                                   upperTheta = (oldValue + dualTolerance_) / alpha;
                                   //tentativeTheta = CoinMin(2.0*upperTheta,tentativeTheta);
                              }
                              // add to list
                              spare[numberRemaining] = alpha;
                              index[numberRemaining++] = iSequence + addSequence;
                         }
                         break;
                    }
               }
          }
     }
     upperReturn = upperTheta;
     bestReturn = bestPossible;
     return numberRemaining;
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
