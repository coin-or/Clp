/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "ClpPresolve.hpp"
#include "Idiot.hpp"
#include "CoinTime.hpp"
#include "CoinSort.hpp"
#include "CoinFactorization.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinHelperFunctions.hpp"
#ifdef CLP_HAS_ABC
#include "AbcCommon.hpp"
#endif
#include "ClpEventHandler.hpp"
// Redefine stuff for Clp
#ifndef OSI_IDIOT
#include "ClpMessage.hpp"
#define OsiObjOffset ClpObjOffset
#endif
/**** strategy 4 - drop, exitDrop and djTolerance all relative:
For first two major iterations these are small.  Then:

drop - exit a major iteration if drop over 5*checkFrequency < this is
used as info->drop*(10.0+fabs(last weighted objective))

exitDrop - exit idiot if feasible and drop < this is
used as info->exitDrop*(10.0+fabs(last objective))

djExit - exit a major iteration if largest dj (averaged over 5 checks)
drops below this - used as info->djTolerance*(10.0+fabs(last weighted objective)

djFlag - mostly skip variables with bad dj worse than this => 2*djExit

djTol - only look at variables with dj better than this => 0.01*djExit
****************/

#define IDIOT_FIX_TOLERANCE 1e-6
#define SMALL_IDIOT_FIX_TOLERANCE 1e-10
int Idiot::dropping(IdiotResult result,
  double tolerance,
  double small,
  int *nbad)
{
  if (result.infeas <= small) {
    double value = CoinMax(fabs(result.objval), fabs(result.dropThis)) + 1.0;
    if (result.dropThis > tolerance * value) {
      *nbad = 0;
      return 1;
    } else {
      (*nbad)++;
      if (*nbad > 4) {
        return 0;
      } else {
        return 1;
      }
    }
  } else {
    *nbad = 0;
    return 1;
  }
}
// Deals with whenUsed and slacks
int Idiot::cleanIteration(int iteration, int ordinaryStart, int ordinaryEnd,
  double *colsol, const double *lower, const double *upper,
  const double *rowLower, const double *rowUpper,
  const double *cost, const double *element, double fixTolerance,
  double &objValue, double &infValue, double &maxInfeasibility)
{
  int n = 0;
  if ((strategy_ & 16384) == 0) {
    for (int i = ordinaryStart; i < ordinaryEnd; i++) {
      if (colsol[i] > lower[i] + fixTolerance) {
        if (colsol[i] < upper[i] - fixTolerance) {
          n++;
        } else {
          colsol[i] = upper[i];
        }
        whenUsed_[i] = iteration;
      } else {
        colsol[i] = lower[i];
      }
    }
    return n;
  } else {
#ifdef COIN_DEVELOP
    printf("entering inf %g, obj %g\n", infValue, objValue);
#endif
    int nrows = model_->getNumRows();
    int ncols = model_->getNumCols();
    int *posSlack = whenUsed_ + ncols;
    int *negSlack = posSlack + nrows;
    int *nextSlack = negSlack + nrows;
    double *rowsol = reinterpret_cast< double * >(nextSlack + ncols);
    memset(rowsol, 0, nrows * sizeof(double));
#ifdef OSI_IDIOT
    const CoinPackedMatrix *matrix = model_->getMatrixByCol();
#else
    // safer for odd matrices
    const CoinPackedMatrix *matrix = model_->matrix();
    //ClpMatrixBase * matrix = model_->clpMatrix();
#endif
    const int *row = matrix->getIndices();
    const CoinBigIndex *columnStart = matrix->getVectorStarts();
    const int *columnLength = matrix->getVectorLengths();
    //const double * element = matrix->getElements();
    int i;
    objValue = 0.0;
    infValue = 0.0;
    maxInfeasibility = 0.0;
    for (i = 0; i < ncols; i++) {
      if (nextSlack[i] == -1) {
        // not a slack
        if (colsol[i] > lower[i] + fixTolerance) {
          if (colsol[i] < upper[i] - fixTolerance) {
            n++;
            whenUsed_[i] = iteration;
          } else {
            colsol[i] = upper[i];
          }
          whenUsed_[i] = iteration;
        } else {
          colsol[i] = lower[i];
        }
        double value = colsol[i];
        if (value) {
          objValue += cost[i] * value;
          CoinBigIndex j;
          for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            rowsol[iRow] += value * element[j];
          }
        }
      }
    }
    // temp fix for infinite lbs - just limit to -1000
    for (i = 0; i < nrows; i++) {
      double rowSave = rowsol[i];
      int iCol;
      iCol = posSlack[i];
      if (iCol >= 0) {
        // slide all slack down
        double rowValue = rowsol[i];
        CoinBigIndex j = columnStart[iCol];
        double lowerValue = CoinMax(CoinMin(colsol[iCol], 0.0) - 1000.0, lower[iCol]);
        rowSave += (colsol[iCol] - lowerValue) * element[j];
        colsol[iCol] = lowerValue;
        while (nextSlack[iCol] >= 0) {
          iCol = nextSlack[iCol];
          double lowerValue = CoinMax(CoinMin(colsol[iCol], 0.0) - 1000.0, lower[iCol]);
          j = columnStart[iCol];
          rowSave += (colsol[iCol] - lowerValue) * element[j];
          colsol[iCol] = lowerValue;
        }
        iCol = posSlack[i];
        while (rowValue < rowLower[i] && iCol >= 0) {
          // want to increase
          double distance = rowLower[i] - rowValue;
          double value = element[columnStart[iCol]];
          double thisCost = cost[iCol];
          if (distance <= value * (upper[iCol] - colsol[iCol])) {
            // can get there
            double movement = distance / value;
            objValue += movement * thisCost;
            rowValue = rowLower[i];
            colsol[iCol] += movement;
          } else {
            // can't get there
            double movement = upper[iCol] - colsol[iCol];
            objValue += movement * thisCost;
            rowValue += movement * value;
            colsol[iCol] = upper[iCol];
            iCol = nextSlack[iCol];
          }
        }
        if (iCol >= 0) {
          // may want to carry on - because of cost?
          while (iCol >= 0 && cost[iCol] < 0 && rowValue < rowUpper[i]) {
            // want to increase
            double distance = rowUpper[i] - rowValue;
            double value = element[columnStart[iCol]];
            double thisCost = cost[iCol];
            if (distance <= value * (upper[iCol] - colsol[iCol])) {
              // can get there
              double movement = distance / value;
              objValue += movement * thisCost;
              rowValue = rowUpper[i];
              colsol[iCol] += movement;
              iCol = -1;
            } else {
              // can't get there
              double movement = upper[iCol] - colsol[iCol];
              objValue += movement * thisCost;
              rowValue += movement * value;
              colsol[iCol] = upper[iCol];
              iCol = nextSlack[iCol];
            }
          }
          if (iCol >= 0 && colsol[iCol] > lower[iCol] + fixTolerance && colsol[iCol] < upper[iCol] - fixTolerance) {
            whenUsed_[i] = iteration;
            n++;
          }
        }
        rowsol[i] = rowValue;
      }
      iCol = negSlack[i];
      if (iCol >= 0) {
        // slide all slack down
        double rowValue = rowsol[i];
        CoinBigIndex j = columnStart[iCol];
        double lowerValue = CoinMax(CoinMin(colsol[iCol], 0.0) - 1000.0, lower[iCol]);
        rowSave += (colsol[iCol] - lowerValue) * element[j];
        colsol[iCol] = lowerValue;
        while (nextSlack[iCol] >= 0) {
          iCol = nextSlack[iCol];
          j = columnStart[iCol];
          double lowerValue = CoinMax(CoinMin(colsol[iCol], 0.0) - 1000.0, lower[iCol]);
          rowSave += (colsol[iCol] - lowerValue) * element[j];
          colsol[iCol] = lowerValue;
        }
        iCol = negSlack[i];
        while (rowValue > rowUpper[i] && iCol >= 0) {
          // want to increase
          double distance = -(rowUpper[i] - rowValue);
          double value = -element[columnStart[iCol]];
          double thisCost = cost[iCol];
          if (distance <= value * (upper[iCol] - lower[iCol])) {
            // can get there
            double movement = distance / value;
            objValue += movement * thisCost;
            rowValue = rowUpper[i];
            colsol[iCol] += movement;
          } else {
            // can't get there
            double movement = upper[iCol] - lower[iCol];
            objValue += movement * thisCost;
            rowValue -= movement * value;
            colsol[iCol] = upper[iCol];
            iCol = nextSlack[iCol];
          }
        }
        if (iCol >= 0) {
          // may want to carry on - because of cost?
          while (iCol >= 0 && cost[iCol] < 0 && rowValue > rowLower[i]) {
            // want to increase
            double distance = -(rowLower[i] - rowValue);
            double value = -element[columnStart[iCol]];
            double thisCost = cost[iCol];
            if (distance <= value * (upper[iCol] - colsol[iCol])) {
              // can get there
              double movement = distance / value;
              objValue += movement * thisCost;
              rowValue = rowLower[i];
              colsol[iCol] += movement;
              iCol = -1;
            } else {
              // can't get there
              double movement = upper[iCol] - colsol[iCol];
              objValue += movement * thisCost;
              rowValue -= movement * value;
              colsol[iCol] = upper[iCol];
              iCol = nextSlack[iCol];
            }
          }
          if (iCol >= 0 && colsol[iCol] > lower[iCol] + fixTolerance && colsol[iCol] < upper[iCol] - fixTolerance) {
            whenUsed_[i] = iteration;
            n++;
          }
        }
        rowsol[i] = rowValue;
      }
      double infeasibility = CoinMax(CoinMax(0.0, rowLower[i] - rowsol[i]), rowsol[i] - rowUpper[i]);
      infValue += infeasibility;
      maxInfeasibility = CoinMax(maxInfeasibility, infeasibility);
      // just change
      rowsol[i] -= rowSave;
    }
    return n;
  }
}

/* returns -1 if none or start of costed slacks or -2 if
   there are costed slacks but it is messy */
static int countCostedSlacks(OsiSolverInterface *model)
{
#ifdef OSI_IDIOT
  const CoinPackedMatrix *matrix = model->getMatrixByCol();
#else
  // safer for odd matrices (note really ClpSimplex not OsiSolverInterface)
  const CoinPackedMatrix *matrix = model->matrix();
  //ClpMatrixBase * matrix = model_->clpMatrix();
#endif
  const int *row = matrix->getIndices();
  const CoinBigIndex *columnStart = matrix->getVectorStarts();
  const int *columnLength = matrix->getVectorLengths();
  const double *element = matrix->getElements();
  const double *rowupper = model->getRowUpper();
  int nrows = model->getNumRows();
  int ncols = model->getNumCols();
  int slackStart = ncols - nrows;
  int nSlacks = nrows;
  int i;

  if (ncols <= nrows)
    return -1;
  while (1) {
    for (i = 0; i < nrows; i++) {
      int j = i + slackStart;
      CoinBigIndex k = columnStart[j];
      if (columnLength[j] == 1) {
        if (row[k] != i || element[k] != 1.0) {
          nSlacks = 0;
          break;
        }
      } else {
        nSlacks = 0;
        break;
      }
      if (rowupper[i] <= 0.0) {
        nSlacks = 0;
        break;
      }
    }
    if (nSlacks || !slackStart)
      break;
    slackStart = 0;
  }
  if (!nSlacks)
    slackStart = -1;
  return slackStart;
}
void Idiot::crash(int numberPass, CoinMessageHandler *handler,
  const CoinMessages *messages, bool doCrossover)
{
  // lightweight options
  int numberColumns = model_->getNumCols();
  const double *objective = model_->getObjCoefficients();
  int nnzero = 0;
  double sum = 0.0;
  int i;
  for (i = 0; i < numberColumns; i++) {
    if (objective[i]) {
      sum += fabs(objective[i]);
      nnzero++;
    }
  }
  sum /= static_cast< double >(nnzero + 1);
  if (maxIts_ == 5)
    maxIts_ = 2;
  if (numberPass <= 0)
    majorIterations_ = static_cast< int >(2 + log10(static_cast< double >(numberColumns + 1)));
  else
    majorIterations_ = numberPass;
  // If mu not changed then compute
  if (mu_ == 1e-4)
    mu_ = CoinMax(1.0e-3, sum * 1.0e-5);
  if (maxIts2_ == 100) {
    if (!lightWeight_) {
      maxIts2_ = 105;
    } else if (lightWeight_ == 1) {
      mu_ *= 1000.0;
      maxIts2_ = 23;
    } else if (lightWeight_ == 2) {
      maxIts2_ = 11;
    } else {
      maxIts2_ = 23;
    }
  }
  //printf("setting mu to %g and doing %d passes\n",mu_,majorIterations_);
  if (numberColumns)
    solve2(handler, messages);
#ifndef OSI_IDIOT
  if (doCrossover) {
    double averageInfeas = model_->sumPrimalInfeasibilities() / static_cast< double >(model_->numberRows());
    if ((averageInfeas < 0.01 && (strategy_ & 512) != 0) || (strategy_ & 8192) != 0)
      crossOver(16 + 1);
    else
      crossOver(majorIterations_ < 1000000 ? 3 : 2);
  }
#endif
}

void Idiot::solve()
{
  CoinMessages dummy;
  solve2(NULL, &dummy);
}
void Idiot::solve2(CoinMessageHandler *handler, const CoinMessages *messages)
{
  int strategy = 0;
  double d2;
  int i, n;
  int allOnes = 1;
  int iteration = 0;
  //int iterationTotal = 0;
  int nTry = 0; /* number of tries at same weight */
  double fixTolerance = IDIOT_FIX_TOLERANCE;
  int maxBigIts = maxBigIts_;
  int maxIts = maxIts_;
  int logLevel = logLevel_;
  int saveMajorIterations = majorIterations_;
  majorIterations_ = majorIterations_ % 1000000;
  if (handler) {
    if (handler->logLevel() > 0 && handler->logLevel() < 3)
      logLevel = 1;
    else if (!handler->logLevel())
      logLevel = 0;
    else
      logLevel = 7;
  }
  double djExit = djTolerance_;
  double djFlag = 1.0 + 100.0 * djExit;
  double djTol = 0.00001;
  double mu = mu_;
  double drop = drop_;
  int maxIts2 = maxIts2_;
  double factor = muFactor_;
  double smallInfeas = smallInfeas_;
  double reasonableInfeas = reasonableInfeas_;
  double stopMu = stopMu_;
  double maxmin, offset;
  double lastWeighted = 1.0e50;
  double exitDrop = exitDrop_;
  double fakeSmall = smallInfeas;
  double firstInfeas;
  int badIts = 0;
  int slackStart, ordStart, ordEnd;
  int checkIteration = 0;
  int lambdaIteration = 0;
  int belowReasonable = 0; /* set if ever gone below reasonable infeas */
  double bestWeighted = 1.0e60;
  double bestFeasible = 1.0e60; /* best solution while feasible */
  IdiotResult result, lastResult;
  int saveStrategy = strategy_;
  const int strategies[] = { 0, 2, 128 };
  double saveLambdaScale = 0.0;
  if ((saveStrategy & 128) != 0) {
    fixTolerance = SMALL_IDIOT_FIX_TOLERANCE;
  }
#ifdef OSI_IDIOT
  const CoinPackedMatrix *matrix = model_->getMatrixByCol();
#else
  // safer for odd matrices
  const CoinPackedMatrix *matrix = model_->matrix();
  //ClpMatrixBase * matrix = model_->clpMatrix();
#endif
  const int *row = matrix->getIndices();
  const CoinBigIndex *columnStart = matrix->getVectorStarts();
  const int *columnLength = matrix->getVectorLengths();
  const double *element = matrix->getElements();
  int nrows = model_->getNumRows();
  int ncols = model_->getNumCols();
  double *rowsol, *colsol;
  double *pi, *dj;
#ifndef OSI_IDIOT
  double *cost = model_->objective();
  double *lower = model_->columnLower();
  double *upper = model_->columnUpper();
#else
  double *cost = new double[ncols];
  CoinMemcpyN(model_->getObjCoefficients(), ncols, cost);
  const double *lower = model_->getColLower();
  const double *upper = model_->getColUpper();
#endif
  const double *elemXX;
  double *saveSol;
  double *rowupper = new double[nrows]; // not const as modified
  CoinMemcpyN(model_->getRowUpper(), nrows, rowupper);
  double *rowlower = new double[nrows]; // not const as modified
  CoinMemcpyN(model_->getRowLower(), nrows, rowlower);
  CoinThreadRandom *randomNumberGenerator = model_->randomNumberGenerator();
  int *whenUsed;
  double *lambda;
  saveSol = new double[ncols];
  lambda = new double[nrows];
  rowsol = new double[nrows];
  colsol = new double[ncols];
  CoinMemcpyN(model_->getColSolution(), ncols, colsol);
  pi = new double[nrows];
  dj = new double[ncols];
#ifndef OSI_IDIOT
  bool fixAfterSome = false; //(model_->specialOptions()&8388608)!=0;
  int exitAfter = 50; //(model_->specialOptions()&8388608)!=0 ? 50 : 1000000;
  {
    int numberColumns = model_->numberColumns();
    for (int i = 0; i < numberColumns; i++) {
      if (upper[i] == lower[i])
        model_->setColumnStatus(i, ClpSimplex::isFixed);
    }
  }
#endif
  delete[] whenUsed_;
  bool oddSlacks = false;
  // See if any costed slacks
  int numberSlacks = 0;
  for (i = 0; i < ncols; i++) {
    if (columnLength[i] == 1)
      numberSlacks++;
  }
  if (!numberSlacks || (strategy_ & 524288) != 0) {
    whenUsed_ = new int[ncols];
  } else {
#ifdef COIN_DEVELOP
    printf("%d slacks\n", numberSlacks);
#endif
    oddSlacks = true;
    int extra = static_cast< int >(nrows * sizeof(double) / sizeof(int));
    whenUsed_ = new int[2 * ncols + 2 * nrows + extra];
    int *posSlack = whenUsed_ + ncols;
    int *negSlack = posSlack + nrows;
    int *nextSlack = negSlack + nrows;
    for (i = 0; i < nrows; i++) {
      posSlack[i] = -1;
      negSlack[i] = -1;
    }
    for (i = 0; i < ncols; i++)
      nextSlack[i] = -1;
    for (i = 0; i < ncols; i++) {
      if (columnLength[i] == 1) {
        CoinBigIndex j = columnStart[i];
        int iRow = row[j];
        if (element[j] > 0.0) {
          if (posSlack[iRow] == -1) {
            posSlack[iRow] = i;
          } else {
            int iCol = posSlack[iRow];
            while (nextSlack[iCol] >= 0)
              iCol = nextSlack[iCol];
            nextSlack[iCol] = i;
          }
        } else {
          if (negSlack[iRow] == -1) {
            negSlack[iRow] = i;
          } else {
            int iCol = negSlack[iRow];
            while (nextSlack[iCol] >= 0)
              iCol = nextSlack[iCol];
            nextSlack[iCol] = i;
          }
        }
      }
    }
    // now sort
    for (i = 0; i < nrows; i++) {
      int iCol;
      iCol = posSlack[i];
      if (iCol >= 0) {
        CoinBigIndex j = columnStart[iCol];
#ifndef NDEBUG
        int iRow = row[j];
#endif
        assert(element[j] > 0.0);
        assert(iRow == i);
        dj[0] = cost[iCol] / element[j];
        whenUsed_[0] = iCol;
        int n = 1;
        while (nextSlack[iCol] >= 0) {
          iCol = nextSlack[iCol];
          CoinBigIndex j = columnStart[iCol];
#ifndef NDEBUG
          int iRow = row[j];
#endif
          assert(element[j] > 0.0);
          assert(iRow == i);
          dj[n] = cost[iCol] / element[j];
          whenUsed_[n++] = iCol;
        }
        for (j = 0; j < n; j++) {
          int jCol = whenUsed_[j];
          nextSlack[jCol] = -2;
        }
        CoinSort_2(dj, dj + n, whenUsed_);
        // put back
        iCol = whenUsed_[0];
        posSlack[i] = iCol;
        for (j = 1; j < n; j++) {
          int jCol = whenUsed_[j];
          nextSlack[iCol] = jCol;
          iCol = jCol;
        }
      }
      iCol = negSlack[i];
      if (iCol >= 0) {
        CoinBigIndex j = columnStart[iCol];
#ifndef NDEBUG
        int iRow = row[j];
#endif
        assert(element[j] < 0.0);
        assert(iRow == i);
        dj[0] = -cost[iCol] / element[j];
        whenUsed_[0] = iCol;
        int n = 1;
        while (nextSlack[iCol] >= 0) {
          iCol = nextSlack[iCol];
          CoinBigIndex j = columnStart[iCol];
#ifndef NDEBUG
          int iRow = row[j];
#endif
          assert(element[j] < 0.0);
          assert(iRow == i);
          dj[n] = -cost[iCol] / element[j];
          whenUsed_[n++] = iCol;
        }
        for (j = 0; j < n; j++) {
          int jCol = whenUsed_[j];
          nextSlack[jCol] = -2;
        }
        CoinSort_2(dj, dj + n, whenUsed_);
        // put back
        iCol = whenUsed_[0];
        negSlack[i] = iCol;
        for (j = 1; j < n; j++) {
          int jCol = whenUsed_[j];
          nextSlack[iCol] = jCol;
          iCol = jCol;
        }
      }
    }
  }
  whenUsed = whenUsed_;
  if (model_->getObjSense() == -1.0) {
    maxmin = -1.0;
  } else {
    maxmin = 1.0;
  }
  model_->getDblParam(OsiObjOffset, offset);
  if (!maxIts2)
    maxIts2 = maxIts;
  strategy = strategy_;
  strategy &= 3;
  memset(lambda, 0, nrows * sizeof(double));
  slackStart = countCostedSlacks(model_);
  // redo in case odd matrix
  row = matrix->getIndices();
  columnStart = matrix->getVectorStarts();
  columnLength = matrix->getVectorLengths();
  element = matrix->getElements();
  if (slackStart >= 0) {
    COIN_DETAIL_PRINT(printf("This model has costed slacks\n"));
    if (slackStart) {
      ordStart = 0;
      ordEnd = slackStart;
    } else {
      ordStart = nrows;
      ordEnd = ncols;
    }
  } else {
    ordStart = 0;
    ordEnd = ncols;
  }
  if (offset && logLevel > 2) {
    printf("** Objective offset is %g\n", offset);
  }
  /* compute reasonable solution cost */
  for (i = 0; i < nrows; i++) {
    rowsol[i] = 1.0e31;
  }
  for (i = 0; i < ncols; i++) {
    CoinBigIndex j;
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      if (element[j] != 1.0) {
        allOnes = 0;
        break;
      }
    }
  }
  if (allOnes) {
    elemXX = NULL;
  } else {
    elemXX = element;
  }
  // Do scaling if wanted
  bool scaled = false;
#ifndef OSI_IDIOT
  if ((strategy_ & 32) != 0 && !allOnes) {
    if (model_->scalingFlag() > 0)
      scaled = model_->clpMatrix()->scale(model_) == 0;
    if (scaled) {
#define IDIOT_SCALE 2
#ifndef IDIOT_SCALE
      const double *rowScale = model_->rowScale();
#else
      double *rowScale = model_->mutableRowScale();
#endif
      const double *columnScale = model_->columnScale();
      double *oldLower = lower;
      double *oldUpper = upper;
      double *oldCost = cost;
      lower = new double[ncols];
      upper = new double[ncols];
      cost = new double[ncols];
      CoinMemcpyN(oldLower, ncols, lower);
      CoinMemcpyN(oldUpper, ncols, upper);
      CoinMemcpyN(oldCost, ncols, cost);
      int icol, irow;
#if IDIOT_SCALE < 0
      for (irow = 0; irow < nrows; irow++) {
        rowlower[irow] = 1.0e100;
        rowupper[irow] = 1.0e-100;
      }
#endif
      for (icol = 0; icol < ncols; icol++) {
        double multiplier = 1.0 / columnScale[icol];
        if (lower[icol] > -1.0e50)
          lower[icol] *= multiplier;
        if (upper[icol] < 1.0e50)
          upper[icol] *= multiplier;
        colsol[icol] *= multiplier;
        cost[icol] *= columnScale[icol];
#if IDIOT_SCALE < 0
        CoinBigIndex j;
        double scale = columnScale[i];
        for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
          int jrow = row[j];
          double scaledValue = fabs(scale * element[j]);
          rowlower[jrow] = CoinMin(rowlower[jrow], scaledValue);
          rowupper[jrow] = CoinMax(rowupper[jrow], scaledValue);
        }
#endif
      }
#ifdef IDIOT_SCALE
#if IDIOT_SCALE > 1 || IDIOT_SCALE < -1
      const double *rowLower = model_->rowLower();
      const double *rowUpper = model_->rowUpper();
#endif
      for (irow = 0; irow < nrows; irow++) {
#if IDIOT_SCALE < 0
        double multiplier = 1.0 / sqrt(rowlower[irow] * rowupper[irow]);
#else
        double multiplier = rowScale[irow];
#endif
#if IDIOT_SCALE > 1 || IDIOT_SCALE < -1
#define EQUALITY_MULTIPLIER 2
        if (rowLower[irow] == rowUpper[irow])
          multiplier *= EQUALITY_MULTIPLIER;
#if IDIOT_SCALE > 2 || IDIOT_SCALE < -2
        if (rowLower[irow] == rowUpper[irow] && !rowlower[irow])
          multiplier *= EQUALITY_MULTIPLIER;
#endif
#endif
        rowScale[irow] = multiplier;
      }
      CoinMemcpyN(model_->rowUpper(), nrows, rowupper);
#endif
      CoinMemcpyN(model_->rowLower(), nrows, rowlower);
      for (irow = 0; irow < nrows; irow++) {
        double multiplier = rowScale[irow];
        if (rowlower[irow] > -1.0e50)
          rowlower[irow] *= multiplier;
        if (rowupper[irow] < 1.0e50)
          rowupper[irow] *= multiplier;
        rowsol[irow] *= multiplier;
      }
      CoinBigIndex length = columnStart[ncols - 1] + columnLength[ncols - 1];
      double *elemYY = new double[length];
      for (i = 0; i < ncols; i++) {
        CoinBigIndex j;
        double scale = columnScale[i];
        for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
          int irow = row[j];
          elemYY[j] = element[j] * scale * rowScale[irow];
        }
      }
      elemXX = elemYY;
    }
  }
#endif
  if ((strategy_ & 131072) != 0) {
    // going to mess with cost
    if (cost == model_->objective())
      cost = CoinCopyOfArray(cost, ncols);
  }
  for (i = 0; i < ncols; i++) {
    CoinBigIndex j;
    double dd = columnLength[i];
    dd = cost[i] / dd;
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      int irow = row[j];
      if (dd < rowsol[irow]) {
        rowsol[irow] = dd;
      }
    }
  }
  d2 = 0.0;
  for (i = 0; i < nrows; i++) {
    d2 += rowsol[i];
  }
  d2 *= 2.0; /* for luck */

  d2 = d2 / static_cast< double >(4 * nrows + 8000);
  d2 *= 0.5; /* halve with more flexible method */
  if (d2 < 5.0)
    d2 = 5.0;
  if (djExit == 0.0) {
    djExit = d2;
  }
  if ((saveStrategy & 4) != 0) {
    /* go to relative tolerances - first small */
    djExit = 1.0e-10;
    djFlag = 1.0e-5;
    drop = 1.0e-10;
  }
  memset(whenUsed, 0, ncols * sizeof(int));
  strategy = strategies[strategy];
  if ((saveStrategy & 8) != 0)
    strategy |= 64; /* don't allow large theta */
  CoinMemcpyN(colsol, ncols, saveSol);

  lastResult = IdiSolve(nrows, ncols, rowsol, colsol, pi,
    dj, cost, rowlower, rowupper,
    lower, upper, elemXX, row, columnStart, columnLength, lambda,
    0, mu, drop,
    maxmin, offset, strategy, djTol, djExit, djFlag, randomNumberGenerator);
  // update whenUsed_
  double maxInfeasibility = COIN_DBL_MAX;
  n = cleanIteration(iteration, ordStart, ordEnd,
    colsol, lower, upper,
    rowlower, rowupper,
    cost, elemXX, fixTolerance, lastResult.objval, lastResult.infeas, maxInfeasibility);
  if ((strategy_ & 16384) != 0) {
    int *posSlack = whenUsed_ + ncols;
    int *negSlack = posSlack + nrows;
    int *nextSlack = negSlack + nrows;
    double *rowsol2 = reinterpret_cast< double * >(nextSlack + ncols);
    for (i = 0; i < nrows; i++)
      rowsol[i] += rowsol2[i];
  }
  if ((logLevel_ & 1) != 0) {
#ifndef OSI_IDIOT
    if (!handler) {
#endif
      printf("Iteration %d infeasibility %g, objective %g - mu %g, its %d, %d interior\n",
        iteration, lastResult.infeas, lastResult.objval, mu, lastResult.iteration, n);
#ifndef OSI_IDIOT
    } else {
      handler->message(CLP_IDIOT_ITERATION, *messages)
        << iteration << lastResult.infeas << lastResult.objval << mu << lastResult.iteration << n
        << CoinMessageEol;
    }
#endif
  }
  int numberBaseTrys = 0; // for first time
  int numberAway = -1;
  //iterationTotal = lastResult.iteration;
  firstInfeas = lastResult.infeas;
  if ((strategy_ & 1024) != 0)
    reasonableInfeas = 0.5 * firstInfeas;
  if (lastResult.infeas < reasonableInfeas)
    lastResult.infeas = reasonableInfeas;
  double keepinfeas = 1.0e31;
  double lastInfeas = 1.0e31;
  double bestInfeas = 1.0e31;
  while ((mu > stopMu && lastResult.infeas > smallInfeas) || (lastResult.infeas <= smallInfeas && dropping(lastResult, exitDrop, smallInfeas, &badIts)) || checkIteration < 2 || lambdaIteration < lambdaIterations_) {
    if (lastResult.infeas <= exitFeasibility_)
      break;
    iteration++;
    //if (iteration>=exitAfter)
    //break;
    checkIteration++;
    if (lastResult.infeas <= smallInfeas && lastResult.objval < bestFeasible) {
      bestFeasible = lastResult.objval;
    }
    if (lastResult.infeas + mu * lastResult.objval < bestWeighted) {
      bestWeighted = lastResult.objval + mu * lastResult.objval;
    }
    if ((saveStrategy & 4096))
      strategy |= 256;
    if ((saveStrategy & 4) != 0 && iteration > 2) {
      /* go to relative tolerances */
      double weighted = 10.0 + fabs(lastWeighted);
      djExit = djTolerance_ * weighted;
      djFlag = 2.0 * djExit;
      drop = drop_ * weighted;
      djTol = 0.01 * djExit;
    }
    CoinMemcpyN(colsol, ncols, saveSol);
    result = IdiSolve(nrows, ncols, rowsol, colsol, pi, dj,
      cost, rowlower, rowupper,
      lower, upper, elemXX, row, columnStart, columnLength, lambda,
      maxIts, mu, drop,
      maxmin, offset, strategy, djTol, djExit, djFlag, randomNumberGenerator);
    n = cleanIteration(iteration, ordStart, ordEnd,
      colsol, lower, upper,
      rowlower, rowupper,
      cost, elemXX, fixTolerance, result.objval, result.infeas, maxInfeasibility);
    if ((strategy_ & 16384) != 0) {
      int *posSlack = whenUsed_ + ncols;
      int *negSlack = posSlack + nrows;
      int *nextSlack = negSlack + nrows;
      double *rowsol2 = reinterpret_cast< double * >(nextSlack + ncols);
      for (i = 0; i < nrows; i++)
        rowsol[i] += rowsol2[i];
    } else {
      maxInfeasibility = 0.0;
      for (i = 0; i < nrows; i++) {
        double infeasibility = CoinMax(CoinMax(0.0, rowlower[i] - rowsol[i]), rowsol[i] - rowupper[i]);
        maxInfeasibility = CoinMax(maxInfeasibility, infeasibility);
      }
    }
    if ((logLevel_ & 1) != 0) {
#ifndef OSI_IDIOT
      if (!handler) {
#endif
        printf("Iteration %d infeasibility %g, objective %g - mu %g, its %d, %d interior\n",
          iteration, result.infeas, result.objval, mu, result.iteration, n);
#ifndef OSI_IDIOT
      } else {
        handler->message(CLP_IDIOT_ITERATION, *messages)
          << iteration << result.infeas << result.objval << mu << result.iteration << n
          << CoinMessageEol;
      }
#endif
    }
#ifndef OSI_IDIOT
    if (fixAfterSome) {
      if (result.infeas < 0.01 * nrows && iteration > 10 && (3 * n > 2 * nrows || 4 * n > 2 * ncols)) {
        // flag
        int numberColumns = model_->numberColumns();
        printf("Flagging satisfied\n");
        fixAfterSome = false;
        for (int i = 0; i < numberColumns; i++) {
          if (colsol[i] > upper[i] - 1.0e-7 || colsol[i] < lower[i] + 1.0e-7)
            model_->setColumnStatus(i, ClpSimplex::isFixed);
        }
      }
    }
#endif
    if (iteration > exitAfter) {
      if ((result.infeas < 1.0e-4 && majorIterations_ < 200 && n == numberAway) || result.infeas < 1.0e-8 || maxInfeasibility < 1.0e-8) {
#ifdef CLP_INVESTIGATE
        printf("infeas small %g\n", result.infeas);
#endif
        if ((strategy_ & 131072) == 0) {
          break; // not much happening
        } else {
          int n = 0;
          for (int i = 0; i < ncols; i++) {
            if (cost[i])
              n++;
          }
          if (n * 10 < ncols) {
            // fix ones with costs
            exitAfter = 100;
            for (int i = 0; i < ncols; i++) {
              if (cost[i] || colsol[i] < lower[i] + 1.0e-9 || colsol[i] > upper[i] - 1.0e-9) {
                cost[i] = 0.0;
                model_->setColumnStatus(i, ClpSimplex::isFixed);
              } else {
                cost[i] = 0.0;
                double gap = upper[i] - lower[i];
                if (colsol[i] > upper[i] - 0.25 * gap) {
                  cost[i] = gap / (colsol[i] - upper[i]);
                } else if (colsol[i] < lower[i] + 0.25 * gap) {
                  cost[i] = gap / (colsol[i] - lower[i]);
                }
              }
            }
          }
        }
      }
    }
    if (lightWeight_ == 1 && iteration > 10 && result.infeas > 1.0 && maxIts != 7) {
      if (lastInfeas != bestInfeas && CoinMin(result.infeas, lastInfeas) > 0.95 * bestInfeas)
        majorIterations_ = CoinMin(majorIterations_, iteration); // not getting feasible
    }
    lastInfeas = result.infeas;
    numberAway = n;
    keepinfeas = result.infeas;
    lastWeighted = result.weighted;
    //iterationTotal += result.iteration;
    if (iteration == 1) {
      if ((strategy_ & 1024) != 0 && mu < 1.0e-10)
        result.infeas = firstInfeas * 0.8;
      if (majorIterations_ >= 50 || dropEnoughFeasibility_ <= 0.0)
        result.infeas *= 0.8;
      if (result.infeas > firstInfeas * 0.9
        && result.infeas > reasonableInfeas) {
        iteration--;
        if (majorIterations_ < 50)
          mu *= 1.0e-1;
        else
          mu *= 0.7;
        bestFeasible = 1.0e31;
        bestWeighted = 1.0e60;
        numberBaseTrys++;
        if (mu < 1.0e-30 || (numberBaseTrys > 10 && lightWeight_)) {
          // back to all slack basis
          lightWeight_ = 2;
          break;
        }
        CoinMemcpyN(saveSol, ncols, colsol);
      } else {
        maxIts = maxIts2;
        checkIteration = 0;
        if ((strategy_ & 1024) != 0)
          mu *= 1.0e-1;
      }
    } else {
    }
    bestInfeas = CoinMin(bestInfeas, result.infeas);
    if (majorIterations_ > 100 && majorIterations_ < 200) {
      if (iteration == majorIterations_ - 100) {
        // redo
        double muX = mu * 10.0;
        bestInfeas = 1.0e3;
        mu = muX;
        nTry = 0;
      }
    }
    if (iteration) {
      /* this code is in to force it to terminate sometime */
      double changeMu = factor;
      if ((saveStrategy & 64) != 0) {
        keepinfeas = 0.0; /* switch off ranga's increase */
        fakeSmall = smallInfeas;
      } else {
        fakeSmall = -1.0;
      }
      saveLambdaScale = 0.0;
      if (result.infeas > reasonableInfeas || (nTry + 1 == maxBigIts && result.infeas > fakeSmall)) {
        if (result.infeas > lastResult.infeas * (1.0 - dropEnoughFeasibility_) || nTry + 1 == maxBigIts || (result.infeas > lastResult.infeas * 0.9 && result.weighted > lastResult.weighted - dropEnoughWeighted_ * CoinMax(fabs(lastResult.weighted), fabs(result.weighted)))) {
          mu *= changeMu;
          if ((saveStrategy & 32) != 0 && result.infeas < reasonableInfeas && 0) {
            reasonableInfeas = CoinMax(smallInfeas, reasonableInfeas * sqrt(changeMu));
            COIN_DETAIL_PRINT(printf("reasonable infeas now %g\n", reasonableInfeas));
          }
          result.weighted = 1.0e60;
          nTry = 0;
          bestFeasible = 1.0e31;
          bestWeighted = 1.0e60;
          checkIteration = 0;
          lambdaIteration = 0;
#define LAMBDA
#ifdef LAMBDA
          if ((saveStrategy & 2048) == 0) {
            memset(lambda, 0, nrows * sizeof(double));
          }
#else
          memset(lambda, 0, nrows * sizeof(double));
#endif
        } else {
          nTry++;
        }
      } else if (lambdaIterations_ >= 0) {
        /* update lambda  */
        double scale = 1.0 / mu;
        int i, nnz = 0;
        saveLambdaScale = scale;
        lambdaIteration++;
        if ((saveStrategy & 4) == 0)
          drop = drop_ / 50.0;
        if (lambdaIteration > 4 && (((lambdaIteration % 10) == 0 && smallInfeas < keepinfeas) || ((lambdaIteration % 5) == 0 && 1.5 * smallInfeas < keepinfeas))) {
          //printf(" Increasing smallInfeas from %f to %f\n",smallInfeas,1.5*smallInfeas);
          smallInfeas *= 1.5;
        }
        if ((saveStrategy & 2048) == 0) {
          for (i = 0; i < nrows; i++) {
            if (lambda[i])
              nnz++;
            lambda[i] += scale * rowsol[i];
          }
        } else {
          nnz = 1;
#ifdef LAMBDA
          for (i = 0; i < nrows; i++) {
            lambda[i] += scale * rowsol[i];
          }
#else
          for (i = 0; i < nrows; i++) {
            lambda[i] = scale * rowsol[i];
          }
          for (i = 0; i < ncols; i++) {
            CoinBigIndex j;
            double value = cost[i] * maxmin;
            for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
              int irow = row[j];
              value += element[j] * lambda[irow];
            }
            cost[i] = value * maxmin;
          }
          for (i = 0; i < nrows; i++) {
            offset += lambda[i] * rowupper[i];
            lambda[i] = 0.0;
          }
#ifdef DEBUG
          printf("offset %g\n", offset);
#endif
          model_->setDblParam(OsiObjOffset, offset);
#endif
        }
        nTry++;
        if (!nnz) {
          bestFeasible = 1.0e32;
          bestWeighted = 1.0e60;
          checkIteration = 0;
          result.weighted = 1.0e31;
        }
#ifdef DEBUG
        double trueCost = 0.0;
        for (i = 0; i < ncols; i++) {
          int j;
          trueCost += cost[i] * colsol[i];
        }
        printf("True objective %g\n", trueCost - offset);
#endif
      } else {
        nTry++;
      }
      lastResult = result;
      if (result.infeas < reasonableInfeas && !belowReasonable) {
        belowReasonable = 1;
        bestFeasible = 1.0e32;
        bestWeighted = 1.0e60;
        checkIteration = 0;
        result.weighted = 1.0e31;
      }
    }
    if (iteration >= majorIterations_) {
      // If not feasible and crash then dive dive dive
      if (mu > 1.0e-12 && result.infeas > 1.0 && majorIterations_ < 40) {
        mu = 1.0e-30;
        majorIterations_ = iteration + 1;
        stopMu = 0.0;
      } else {
        if (logLevel > 2)
          printf("Exiting due to number of major iterations\n");
        break;
      }
    }
  }
  majorIterations_ = saveMajorIterations;
#ifndef OSI_IDIOT
  if (scaled) {
    // Scale solution and free arrays
    const double *rowScale = model_->rowScale();
    const double *columnScale = model_->columnScale();
    int icol, irow;
    for (icol = 0; icol < ncols; icol++) {
      colsol[icol] *= columnScale[icol];
      saveSol[icol] *= columnScale[icol];
      dj[icol] /= columnScale[icol];
    }
    for (irow = 0; irow < nrows; irow++) {
      rowsol[irow] /= rowScale[irow];
      pi[irow] *= rowScale[irow];
    }
    // Don't know why getting Microsoft problems
#if defined(_MSC_VER)
    delete[](double *) elemXX;
#else
    delete[] elemXX;
#endif
    model_->setRowScale(NULL);
    model_->setColumnScale(NULL);
    delete[] lower;
    delete[] upper;
    delete[] cost;
    lower = model_->columnLower();
    upper = model_->columnUpper();
    cost = model_->objective();
    //rowlower = model_->rowLower();
  } else if (cost != model_->objective()) {
    delete[] cost;
    cost = model_->objective();
  }
#endif
#define TRYTHIS
#ifdef TRYTHIS
  if ((saveStrategy & 2048) != 0) {
    double offset;
    model_->getDblParam(OsiObjOffset, offset);
    for (i = 0; i < ncols; i++) {
      CoinBigIndex j;
      double djval = cost[i] * maxmin;
      for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
        int irow = row[j];
        djval -= element[j] * lambda[irow];
      }
      cost[i] = djval;
    }
    for (i = 0; i < nrows; i++) {
      offset += lambda[i] * rowupper[i];
    }
    model_->setDblParam(OsiObjOffset, offset);
  }
#endif
  if (saveLambdaScale) {
    /* back off last update */
    for (i = 0; i < nrows; i++) {
      lambda[i] -= saveLambdaScale * rowsol[i];
    }
  }
  muAtExit_ = mu;
  // For last iteration make as feasible as possible
  if (oddSlacks)
    strategy_ |= 16384;
  // not scaled
  n = cleanIteration(iteration, ordStart, ordEnd,
    colsol, lower, upper,
    model_->rowLower(), model_->rowUpper(),
    cost, element, fixTolerance, lastResult.objval, lastResult.infeas, maxInfeasibility);
#if 0
     if ((logLevel & 1) == 0 || (strategy_ & 16384) != 0) {
          printf(
               "%d - mu %g, infeasibility %g, objective %g, %d interior\n",
               iteration, mu, lastResult.infeas, lastResult.objval, n);
     }
#endif
#ifndef OSI_IDIOT
  model_->setSumPrimalInfeasibilities(lastResult.infeas);
#endif
  // Put back more feasible solution
  double saveInfeas[] = { 0.0, 0.0 };
  for (int iSol = 0; iSol < 3; iSol++) {
    const double *solution = iSol ? colsol : saveSol;
    if (iSol == 2 && saveInfeas[0] < saveInfeas[1]) {
      // put back best solution
      CoinMemcpyN(saveSol, ncols, colsol);
    }
    double large = 0.0;
    int i;
    memset(rowsol, 0, nrows * sizeof(double));
    for (i = 0; i < ncols; i++) {
      CoinBigIndex j;
      double value = solution[i];
      for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
        int irow = row[j];
        rowsol[irow] += element[j] * value;
      }
    }
    for (i = 0; i < nrows; i++) {
      if (rowsol[i] > rowupper[i]) {
        double diff = rowsol[i] - rowupper[i];
        if (diff > large)
          large = diff;
      } else if (rowsol[i] < rowlower[i]) {
        double diff = rowlower[i] - rowsol[i];
        if (diff > large)
          large = diff;
      }
    }
    if (iSol < 2)
      saveInfeas[iSol] = large;
    if (logLevel > 2)
      printf("largest infeasibility is %g\n", large);
  }
  /* subtract out lambda */
  for (i = 0; i < nrows; i++) {
    pi[i] -= lambda[i];
  }
  for (i = 0; i < ncols; i++) {
    CoinBigIndex j;
    double djval = cost[i] * maxmin;
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      int irow = row[j];
      djval -= element[j] * pi[irow];
    }
    dj[i] = djval;
  }
  if ((strategy_ & 1024) != 0) {
    double ratio = static_cast< double >(ncols) / static_cast< double >(nrows);
    COIN_DETAIL_PRINT(printf("col/row ratio %g infeas ratio %g\n", ratio, lastResult.infeas / firstInfeas));
    if (lastResult.infeas > 0.01 * firstInfeas * ratio) {
      strategy_ &= (~1024);
      COIN_DETAIL_PRINT(printf(" - layer off\n"));
    } else {
      COIN_DETAIL_PRINT(printf(" - layer on\n"));
    }
  }
  delete[] saveSol;
  delete[] lambda;
  // save solution
  // duals not much use - but save anyway
#ifndef OSI_IDIOT
  CoinMemcpyN(rowsol, nrows, model_->primalRowSolution());
  CoinMemcpyN(colsol, ncols, model_->primalColumnSolution());
  CoinMemcpyN(pi, nrows, model_->dualRowSolution());
  CoinMemcpyN(dj, ncols, model_->dualColumnSolution());
#else
  model_->setColSolution(colsol);
  model_->setRowPrice(pi);
  delete[] cost;
#endif
  delete[] rowsol;
  delete[] colsol;
  delete[] pi;
  delete[] dj;
  delete[] rowlower;
  delete[] rowupper;
  return;
}
#ifndef OSI_IDIOT
void Idiot::crossOver(int mode)
{
  if (lightWeight_ == 2) {
    // total failure
    model_->allSlackBasis();
    return;
  }
  double fixTolerance = IDIOT_FIX_TOLERANCE;
#ifdef COIN_DEVELOP
  double startTime = CoinCpuTime();
#endif
  ClpSimplex *saveModel = NULL;
  ClpMatrixBase *matrix = model_->clpMatrix();
  const int *row = matrix->getIndices();
  const CoinBigIndex *columnStart = matrix->getVectorStarts();
  const int *columnLength = matrix->getVectorLengths();
  const double *element = matrix->getElements();
  const double *rowupper = model_->getRowUpper();
  model_->eventHandler()->event(ClpEventHandler::startOfCrossover);
  int nrows = model_->getNumRows();
  int ncols = model_->getNumCols();
  double *rowsol, *colsol;
  // different for Osi
  double *lower = model_->columnLower();
  double *upper = model_->columnUpper();
  const double *rowlower = model_->getRowLower();
  int *whenUsed = whenUsed_;
  rowsol = model_->primalRowSolution();
  colsol = model_->primalColumnSolution();
  ;
  double *cost = model_->objective();
  int slackEnd, ordStart, ordEnd;
  int slackStart = countCostedSlacks(model_);

  int addAll = mode & 7;
  int presolve = 0;

  double djTolerance = djTolerance_;
  if (djTolerance > 0.0 && djTolerance < 1.0)
    djTolerance = 1.0;
  int iteration;
  int i, n = 0;
  double ratio = 1.0;
#ifdef COIN_DETAIL
  double objValue = 0.0;
#endif
  if ((strategy_ & 128) != 0) {
    fixTolerance = SMALL_IDIOT_FIX_TOLERANCE;
  }
  if ((mode & 16) != 0 && addAll < 3)
    presolve = 1;
  double *saveUpper = NULL;
  double *saveLower = NULL;
  double *saveRowUpper = NULL;
  double *saveRowLower = NULL;
  bool allowInfeasible = ((strategy_ & 8192) != 0) || (majorIterations_ > 1000000);
  if (addAll < 3) {
    saveUpper = new double[ncols];
    saveLower = new double[ncols];
    CoinMemcpyN(upper, ncols, saveUpper);
    CoinMemcpyN(lower, ncols, saveLower);
    if (allowInfeasible) {
      saveRowUpper = new double[nrows];
      saveRowLower = new double[nrows];
      CoinMemcpyN(rowupper, nrows, saveRowUpper);
      CoinMemcpyN(rowlower, nrows, saveRowLower);
      double averageInfeas = model_->sumPrimalInfeasibilities() / static_cast< double >(model_->numberRows());
      fixTolerance = CoinMax(fixTolerance, 1.0e-5 * averageInfeas);
    }
  }
  if (slackStart >= 0) {
    slackEnd = slackStart + nrows;
    if (slackStart) {
      ordStart = 0;
      ordEnd = slackStart;
    } else {
      ordStart = nrows;
      ordEnd = ncols;
    }
  } else {
    slackEnd = slackStart;
    ordStart = 0;
    ordEnd = ncols;
  }
  /* get correct rowsol (without known slacks) */
  memset(rowsol, 0, nrows * sizeof(double));
  for (i = ordStart; i < ordEnd; i++) {
    CoinBigIndex j;
    double value = colsol[i];
    if (value < lower[i] + fixTolerance) {
      value = lower[i];
      colsol[i] = value;
    }
    for (j = columnStart[i]; j < columnStart[i] + columnLength[i]; j++) {
      int irow = row[j];
      rowsol[irow] += value * element[j];
    }
  }
  if (slackStart >= 0) {
    for (i = 0; i < nrows; i++) {
      if (ratio * rowsol[i] > rowlower[i] && rowsol[i] > 1.0e-8) {
        ratio = rowlower[i] / rowsol[i];
      }
    }
    for (i = 0; i < nrows; i++) {
      rowsol[i] *= ratio;
    }
    for (i = ordStart; i < ordEnd; i++) {
      double value = colsol[i] * ratio;
      colsol[i] = value;
#ifdef COIN_DETAIL
      objValue += value * cost[i];
#endif
    }
    for (i = 0; i < nrows; i++) {
      double value = rowlower[i] - rowsol[i];
      colsol[i + slackStart] = value;
#ifdef COIN_DETAIL
      objValue += value * cost[i + slackStart];
#endif
    }
    COIN_DETAIL_PRINT(printf("New objective after scaling %g\n", objValue));
  }
#if 0
     //maybe put back - but just get feasible ?
     // If not many fixed then just exit
     int numberFixed = 0;
     for (i = ordStart; i < ordEnd; i++) {
          if (colsol[i] < lower[i] + fixTolerance)
               numberFixed++;
          else if (colsol[i] > upper[i] - fixTolerance)
               numberFixed++;
     }
     if (numberFixed < ncols / 3) {
          addAll = 3;
          presolve = 0;
     }
#endif
#define FEB_TRY
#ifdef FEB_TRY
  int savePerturbation = model_->perturbation();
  int saveOptions = model_->specialOptions();
  model_->setSpecialOptions(saveOptions | 8192);
  //if (savePerturbation_ == 50)
  //   model_->setPerturbation(56);
#endif
  model_->createStatus();
  /* addAll
        0 - chosen,all used, all
        1 - chosen, all
        2 - all
        3 - do not do anything  - maybe basis
     */
  for (i = ordStart; i < ordEnd; i++) {
    if (addAll < 2) {
      if (colsol[i] < lower[i] + fixTolerance) {
        upper[i] = lower[i];
        colsol[i] = lower[i];
      } else if (colsol[i] > upper[i] - fixTolerance) {
        lower[i] = upper[i];
        colsol[i] = upper[i];
      }
    }
    model_->setColumnStatus(i, ClpSimplex::superBasic);
  }
  if ((strategy_ & 16384) != 0) {
    // put in basis
    int *posSlack = whenUsed_ + ncols;
    int *negSlack = posSlack + nrows;
    int *nextSlack = negSlack + nrows;
    /* Laci - try both ways - to see what works -
             you can change second part as much as you want */
#ifndef LACI_TRY // was #if 1
    // Array for sorting out slack values
    double *ratio = new double[ncols];
    int *which = new int[ncols];
    for (i = 0; i < nrows; i++) {
      if (posSlack[i] >= 0 || negSlack[i] >= 0) {
        int iCol;
        int nPlus = 0;
        int nMinus = 0;
        bool possible = true;
        // Get sum
        double sum = 0.0;
        iCol = posSlack[i];
        while (iCol >= 0) {
          double value = element[columnStart[iCol]];
          sum += value * colsol[iCol];
          if (lower[iCol]) {
            possible = false;
            break;
          } else {
            nPlus++;
          }
          iCol = nextSlack[iCol];
        }
        iCol = negSlack[i];
        while (iCol >= 0) {
          double value = -element[columnStart[iCol]];
          sum -= value * colsol[iCol];
          if (lower[iCol]) {
            possible = false;
            break;
          } else {
            nMinus++;
          }
          iCol = nextSlack[iCol];
        }
        //printf("%d plus, %d minus",nPlus,nMinus);
        //printf("\n");
        if ((rowsol[i] - rowlower[i] < 1.0e-7 || rowupper[i] - rowsol[i] < 1.0e-7) && nPlus + nMinus < 2)
          possible = false;
        if (possible) {
          // Amount contributed by other varaibles
          sum = rowsol[i] - sum;
          double lo = rowlower[i];
          if (lo > -1.0e20)
            lo -= sum;
          double up = rowupper[i];
          if (up < 1.0e20)
            up -= sum;
          //printf("row bounds %g %g\n",lo,up);
          if (0) {
            double sum = 0.0;
            double x = 0.0;
            for (int k = 0; k < ncols; k++) {
              CoinBigIndex j;
              double value = colsol[k];
              x += value * cost[k];
              for (j = columnStart[k]; j < columnStart[k] + columnLength[k]; j++) {
                int irow = row[j];
                if (irow == i)
                  sum += element[j] * value;
              }
            }
            printf("Before sum %g <= %g <= %g cost %.18g\n",
              rowlower[i], sum, rowupper[i], x);
          }
          // set all to zero
          iCol = posSlack[i];
          while (iCol >= 0) {
            colsol[iCol] = 0.0;
            iCol = nextSlack[iCol];
          }
          iCol = negSlack[i];
          while (iCol >= 0) {
            colsol[iCol] = 0.0;
            iCol = nextSlack[iCol];
          }
          {
            int iCol;
            iCol = posSlack[i];
            while (iCol >= 0) {
              //printf("col %d el %g sol %g bounds %g %g cost %g\n",
              //     iCol,element[columnStart[iCol]],
              //     colsol[iCol],lower[iCol],upper[iCol],cost[iCol]);
              iCol = nextSlack[iCol];
            }
            iCol = negSlack[i];
            while (iCol >= 0) {
              //printf("col %d el %g sol %g bounds %g %g cost %g\n",
              //     iCol,element[columnStart[iCol]],
              //     colsol[iCol],lower[iCol],upper[iCol],cost[iCol]);
              iCol = nextSlack[iCol];
            }
          }
          //printf("now what?\n");
          int n = 0;
          bool basic = false;
          if (lo > 0.0) {
            // Add in positive
            iCol = posSlack[i];
            while (iCol >= 0) {
              double value = element[columnStart[iCol]];
              ratio[n] = cost[iCol] / value;
              which[n++] = iCol;
              iCol = nextSlack[iCol];
            }
            CoinSort_2(ratio, ratio + n, which);
            for (int i = 0; i < n; i++) {
              iCol = which[i];
              double value = element[columnStart[iCol]];
              if (lo >= upper[iCol] * value) {
                value *= upper[iCol];
                sum += value;
                lo -= value;
                colsol[iCol] = upper[iCol];
              } else {
                value = lo / value;
                sum += lo;
                lo = 0.0;
                colsol[iCol] = value;
                model_->setColumnStatus(iCol, ClpSimplex::basic);
                basic = true;
              }
              if (lo < 1.0e-7)
                break;
            }
          } else if (up < 0.0) {
            // Use lo so coding is more similar
            lo = -up;
            // Add in negative
            iCol = negSlack[i];
            while (iCol >= 0) {
              double value = -element[columnStart[iCol]];
              ratio[n] = cost[iCol] / value;
              which[n++] = iCol;
              iCol = nextSlack[iCol];
            }
            CoinSort_2(ratio, ratio + n, which);
            for (int i = 0; i < n; i++) {
              iCol = which[i];
              double value = -element[columnStart[iCol]];
              if (lo >= upper[iCol] * value) {
                value *= upper[iCol];
                sum += value;
                lo -= value;
                colsol[iCol] = upper[iCol];
              } else {
                value = lo / value;
                sum += lo;
                lo = 0.0;
                colsol[iCol] = value;
                model_->setColumnStatus(iCol, ClpSimplex::basic);
                basic = true;
              }
              if (lo < 1.0e-7)
                break;
            }
          }
          if (0) {
            double sum2 = 0.0;
            double x = 0.0;
            for (int k = 0; k < ncols; k++) {
              CoinBigIndex j;
              double value = colsol[k];
              x += value * cost[k];
              for (j = columnStart[k]; j < columnStart[k] + columnLength[k]; j++) {
                int irow = row[j];
                if (irow == i)
                  sum2 += element[j] * value;
              }
            }
            printf("after sum %g <= %g <= %g cost %.18g (sum = %g)\n",
              rowlower[i], sum2, rowupper[i], x, sum);
          }
          rowsol[i] = sum;
          if (basic) {
            if (fabs(rowsol[i] - rowlower[i]) < fabs(rowsol[i] - rowupper[i]))
              model_->setRowStatus(i, ClpSimplex::atLowerBound);
            else
              model_->setRowStatus(i, ClpSimplex::atUpperBound);
          }
        } else {
          int n = 0;
          int iCol;
          iCol = posSlack[i];
          while (iCol >= 0) {
            if (colsol[iCol] > lower[iCol] + 1.0e-8 && colsol[iCol] < upper[iCol] - 1.0e-8) {
              model_->setColumnStatus(iCol, ClpSimplex::basic);
              n++;
            }
            iCol = nextSlack[iCol];
          }
          iCol = negSlack[i];
          while (iCol >= 0) {
            if (colsol[iCol] > lower[iCol] + 1.0e-8 && colsol[iCol] < upper[iCol] - 1.0e-8) {
              model_->setColumnStatus(iCol, ClpSimplex::basic);
              n++;
            }
            iCol = nextSlack[iCol];
          }
          if (n) {
            if (fabs(rowsol[i] - rowlower[i]) < fabs(rowsol[i] - rowupper[i]))
              model_->setRowStatus(i, ClpSimplex::atLowerBound);
            else
              model_->setRowStatus(i, ClpSimplex::atUpperBound);
#ifdef CLP_INVESTIGATE
            if (n > 1)
              printf("%d basic on row %d!\n", n, i);
#endif
          }
        }
      }
    }
    delete[] ratio;
    delete[] which;
#else
    for (i = 0; i < nrows; i++) {
      int n = 0;
      int iCol;
      iCol = posSlack[i];
      while (iCol >= 0) {
        if (colsol[iCol] > lower[iCol] + 1.0e-8 && colsol[iCol] < upper[iCol] - 1.0e-8) {
          model_->setColumnStatus(iCol, ClpSimplex::basic);
          n++;
        }
        iCol = nextSlack[iCol];
      }
      iCol = negSlack[i];
      while (iCol >= 0) {
        if (colsol[iCol] > lower[iCol] + 1.0e-8 && colsol[iCol] < upper[iCol] - 1.0e-8) {
          model_->setColumnStatus(iCol, ClpSimplex::basic);
          n++;
        }
        iCol = nextSlack[iCol];
      }
      if (n) {
        if (fabs(rowsol[i] - rowlower[i]) < fabs(rowsol[i] - rowupper[i]))
          model_->setRowStatus(i, ClpSimplex::atLowerBound);
        else
          model_->setRowStatus(i, ClpSimplex::atUpperBound);
#ifdef CLP_INVESTIGATE
        if (n > 1)
          printf("%d basic on row %d!\n", n, i);
#endif
      }
    }
#endif
  }
  double maxmin;
  if (model_->getObjSense() == -1.0) {
    maxmin = -1.0;
  } else {
    maxmin = 1.0;
  }
  bool justValuesPass = majorIterations_ > 1000000;
  if (slackStart >= 0) {
    for (i = 0; i < nrows; i++) {
      model_->setRowStatus(i, ClpSimplex::superBasic);
    }
    for (i = slackStart; i < slackEnd; i++) {
      model_->setColumnStatus(i, ClpSimplex::basic);
    }
  } else {
    /* still try and put singletons rather than artificials in basis */
    for (i = 0; i < nrows; i++) {
      model_->setRowStatus(i, ClpSimplex::basic);
    }
    int ninbas = 0;
    for (i = 0; i < ncols; i++) {
      if (columnLength[i] == 1 && upper[i] > lower[i] + 1.0e-5) {
        CoinBigIndex j = columnStart[i];
        double value = element[j];
        int irow = row[j];
        double rlo = rowlower[irow];
        double rup = rowupper[irow];
        double clo = lower[i];
        double cup = upper[i];
        double csol = colsol[i];
        /* adjust towards feasibility */
        double move = 0.0;
        if (rowsol[irow] > rup) {
          move = (rup - rowsol[irow]) / value;
          if (value > 0.0) {
            /* reduce */
            if (csol + move < clo)
              move = CoinMin(0.0, clo - csol);
          } else {
            /* increase */
            if (csol + move > cup)
              move = CoinMax(0.0, cup - csol);
          }
        } else if (rowsol[irow] < rlo) {
          move = (rlo - rowsol[irow]) / value;
          if (value > 0.0) {
            /* increase */
            if (csol + move > cup)
              move = CoinMax(0.0, cup - csol);
          } else {
            /* reduce */
            if (csol + move < clo)
              move = CoinMin(0.0, clo - csol);
          }
        } else {
          /* move to improve objective */
          if (cost[i] * maxmin > 0.0) {
            if (value > 0.0) {
              move = (rlo - rowsol[irow]) / value;
              /* reduce */
              if (csol + move < clo)
                move = CoinMin(0.0, clo - csol);
            } else {
              move = (rup - rowsol[irow]) / value;
              /* increase */
              if (csol + move > cup)
                move = CoinMax(0.0, cup - csol);
            }
          } else if (cost[i] * maxmin < 0.0) {
            if (value > 0.0) {
              move = (rup - rowsol[irow]) / value;
              /* increase */
              if (csol + move > cup)
                move = CoinMax(0.0, cup - csol);
            } else {
              move = (rlo - rowsol[irow]) / value;
              /* reduce */
              if (csol + move < clo)
                move = CoinMin(0.0, clo - csol);
            }
          }
        }
        rowsol[irow] += move * value;
        colsol[i] += move;
        /* put in basis if row was artificial */
        if (rup - rlo < 1.0e-7 && model_->getRowStatus(irow) == ClpSimplex::basic) {
          model_->setRowStatus(irow, ClpSimplex::superBasic);
          model_->setColumnStatus(i, ClpSimplex::basic);
          ninbas++;
        }
      }
    }
    /*printf("%d in basis\n",ninbas);*/
  }
  bool wantVector = false;
  if (dynamic_cast< ClpPackedMatrix * >(model_->clpMatrix())) {
    // See if original wanted vector
    ClpPackedMatrix *clpMatrixO = dynamic_cast< ClpPackedMatrix * >(model_->clpMatrix());
    wantVector = clpMatrixO->wantsSpecialColumnCopy();
  }
  if ((strategy_ & 32768) != 0)
    allowInfeasible = true;
  if ((strategy_ & 65536) != 0)
    justValuesPass = true;
  //double * saveBounds=NULL;
  if (addAll < 3) {
    ClpPresolve pinfo;
    if (presolve) {
      double *rhs = new double[nrows];
      double *saveBounds = new double[2 * ncols];
      char line[200];
      memcpy(saveBounds, lower, ncols * sizeof(double));
      memcpy(saveBounds + ncols, upper, ncols * sizeof(double));
      if (allowInfeasible) {
        // fix up so will be feasible
        const double *dual = model_->dualRowSolution();
        for (i = 0; i < nrows; i++)
          rhs[i] = fabs(dual[i]);
        std::sort(rhs, rhs + nrows);
        int nSmall = nrows;
        int nMedium = nrows;
        double largest = rhs[nrows - 1];
        double small = CoinMax(1.0e-4, 1.0e-5 * largest);
        small = CoinMin(small, 1.0e-2);
        double medium = small * 100.0;
        double *rowupper = model_->rowUpper();
        double *rowlower = model_->rowLower();
        // if tiny then drop row??
        for (i = 0; i < nrows; i++) {
          if (rhs[i] >= small) {
            nSmall = i - 1;
            break;
          }
        }
        for (; i < nrows; i++) {
          if (rhs[i] >= medium) {
            nMedium = i - 1;
            break;
          }
        }
        printf("%d < %g, %d < %g, %d <= %g\n",
          nSmall, small, nMedium - nSmall, medium, nrows - nMedium, largest);
        memset(rhs, 0, nrows * sizeof(double));
        int nFixed = 0;
        for (i = 0; i < ncols; i++) {
          if (colsol[i] < lower[i] + 1.0e-8) {
            upper[i] = lower[i];
            colsol[i] = lower[i];
            nFixed++;
          } else if (colsol[i] > upper[i] - 1.0e-8) {
            lower[i] = upper[i];
            colsol[i] = lower[i];
            nFixed++;
          }
        }
        model_->clpMatrix()->times(1.0, colsol, rhs);
        saveRowUpper = CoinCopyOfArray(rowupper, nrows);
        saveRowLower = CoinCopyOfArray(rowlower, nrows);
        double sum = 0.0;
        for (i = 0; i < nrows; i++) {
          if (rhs[i] > rowupper[i]) {
            sum += rhs[i] - rowupper[i];
          }
          if (rhs[i] < rowlower[i]) {
            sum += rowlower[i] - rhs[i];
          }
        }
        double averageInfeasibility = sum / nrows;
        double check = CoinMin(1.0e-3, 0.1 * averageInfeasibility);
        int nFixedRows = 0;
        int nFreed = 0;
#define MESS_UP 0
        for (i = 0; i < nrows; i++) {
          if (rowupper[i] > rowlower[i] + check) {
            // look at distance and sign of dual
            if (dual[i] < -medium && rowupper[i] - rhs[i] < check) {
              rowupper[i] = rhs[i];
              rowlower[i] = rowupper[i];
              nFixedRows++;
            } else if (dual[i] > medium && rhs[i] - rowlower[i] < check) {
              rowlower[i] = rhs[i];
              rowupper[i] = rowlower[i];
              nFixedRows++;
            } else if (fabs(dual[i]) < small && rhs[i] - rowlower[i] > check && rowupper[i] - rhs[i] > check) {
              nFreed++;
#if MESS_UP == 1 || MESS_UP == 2
              rowupper[i] = COIN_DBL_MAX;
              rowlower[i] = -COIN_DBL_MAX;
#endif
            }
          }
          if (rhs[i] > rowupper[i]) {
            rowupper[i] = rhs[i];
            // maybe make equality
#if MESS_UP == 2 || MESS_UP == 3
            rowlower[i] = rhs[i];
#endif
          }
          if (rhs[i] < rowlower[i]) {
            rowlower[i] = rhs[i];
            // maybe make equality
#if MESS_UP == 2 || MESS_UP == 3
            rowupper[i] = rhs[i];
#endif
          }
        }
        sprintf(line, "sum of infeasibilities %g - %d fixed rows, %d fixed columns - might free %d rows",
          sum, nFixedRows, nFixed, nFreed);
      } else {
        memset(rhs, 0, nrows * sizeof(double));
        int nFixed = 0;
        for (i = 0; i < ncols; i++) {
          if (colsol[i] < lower[i] + 1.0e-8) {
            upper[i] = lower[i];
            colsol[i] = lower[i];
            nFixed++;
          } else if (colsol[i] > upper[i] - 1.0e-8) {
            lower[i] = upper[i];
            colsol[i] = lower[i];
            nFixed++;
          }
        }
        model_->clpMatrix()->times(1.0, colsol, rhs);
        double sum = 0.0;
        for (i = 0; i < nrows; i++) {
          if (rhs[i] > rowupper[i]) {
            sum += rhs[i] - rowupper[i];
          }
          if (rhs[i] < rowlower[i]) {
            sum += rowlower[i] - rhs[i];
          }
        }

        double averageInfeasibility = sum / nrows;
        sprintf(line, "sum of infeasibilities %g - average %g, %d fixed columns",
          sum, averageInfeasibility, nFixed);
      }
      const CoinMessages *messages = model_->messagesPointer();
      model_->messageHandler()->message(CLP_GENERAL, *messages)
        << line
        << CoinMessageEol;
      delete[] rhs;
      saveModel = model_;
      pinfo.setPresolveActions(pinfo.presolveActions() | 16384);
      model_ = pinfo.presolvedModel(*model_, 1.0e-8, false, 5);
      if (saveBounds) {
        memcpy(saveModel->columnLower(), saveBounds, ncols * sizeof(double));
        memcpy(saveModel->columnUpper(), saveBounds + ncols, ncols * sizeof(double));
        delete[] saveBounds;
      }
      if (model_ && (strategy_ & 262144) != 0) {
        int nrows = model_->getNumRows();
        int ncols = model_->getNumCols();
        double *lower = model_->columnLower();
        double *upper = model_->columnUpper();
        const double *rowlower = model_->getRowLower();
        const double *rowupper = model_->getRowUpper();
        double *rowsol = model_->primalRowSolution();
        double *colsol = model_->primalColumnSolution();
        ;
        int ninbas = 0;
        int *which = new int[2 * ncols + nrows];
        double *dj = model_->dualColumnSolution();
        for (int i = 0; i < ncols; i++) {
          dj[i] = -CoinMin(upper[i] - colsol[i], colsol[i] - lower[i]);
          which[i] = i;
        }
        CoinSort_2(dj, dj + ncols, which);
        ninbas = CoinMin(ncols, nrows);
        int *columnIsBasic = which + ncols;
        int *rowIsBasic = columnIsBasic + ncols;
        for (int i = 0; i < nrows + ncols; i++)
          columnIsBasic[i] = -1;
        for (int i = 0; i < ninbas; i++) {
          int iColumn = which[i];
          columnIsBasic[iColumn] = i;
        }
        // factorize
        CoinFactorization factor;
        factor.pivotTolerance(0.1);
        factor.setDenseThreshold(0);
        int status = -1;
        // If initial is too dense - then all-slack may be better
        double areaFactor = 1.0; // was 4.0
        const CoinPackedMatrix *matrix = model_->matrix();
        while (status) {
          status = factor.factorize(*matrix, rowIsBasic, columnIsBasic, areaFactor);
          if (status == -99) {
            // put all slacks in
            for (int i = 0; i < nrows; i++)
              rowIsBasic[i] = i;
            for (int i = 0; i < ncols; i++)
              columnIsBasic[i] = -1;
            break;
          } else if (status == -1) {
            factor.pivotTolerance(0.99);
            // put all slacks in
            for (int i = 0; i < nrows; i++)
              rowIsBasic[i] = i;
            for (int i = 0; i < ncols; i++) {
              int iRow = columnIsBasic[i];
              if (iRow >= 0)
                rowIsBasic[iRow] = -1; // out
            }
          }
        }
        for (int i = 0; i < nrows; i++) {
          if (rowIsBasic[i] >= 0) {
            model_->setRowStatus(i, ClpSimplex::basic);
          } else if (rowlower[i] == rowupper[i]) {
            model_->setRowStatus(i, ClpSimplex::isFixed);
          } else if (rowsol[i] - rowlower[i] < rowupper[i] - rowsol[i]) {
            model_->setRowStatus(i, ClpSimplex::atLowerBound);
          } else {
            model_->setRowStatus(i, ClpSimplex::atUpperBound);
          }
        }
        for (int i = 0; i < ncols; i++) {
          if (colsol[i] > upper[i] - 1.0e-7 || colsol[i] < lower[i] + 1.0e-7) {
            model_->setColumnStatus(i, ClpSimplex::isFixed);
          } else if (columnIsBasic[i] >= 0) {
            model_->setColumnStatus(i, ClpSimplex::basic);
          } else {
            model_->setColumnStatus(i, ClpSimplex::superBasic);
          }
        }
        delete[] which;
      }
    }
    if (model_) {
      // See if we want to go all way
      int oldSize = 2 * saveModel->numberRows() + saveModel->numberColumns();
      int newSize = 2 * model_->numberRows() + model_->numberColumns();
      if (oldSize * 2 > newSize * 3)
        justValuesPass = false;
      if (!wantVector) {
        //#define TWO_GOES
#ifdef ABC_INHERIT
#ifndef TWO_GOES
        model_->dealWithAbc(1, justValuesPass ? 3 : 1);
#else
        model_->dealWithAbc(1, 1 + 11);
#endif
#else
#ifndef TWO_GOES
        model_->primal(justValuesPass ? 2 : 1);
#else
        model_->primal(1 + 11);
#endif
#endif
      } else {
        ClpMatrixBase *matrix = model_->clpMatrix();
        ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
        assert(clpMatrix);
        clpMatrix->makeSpecialColumnCopy();
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 1);
#else
        model_->primal(1);
#endif
        clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
        if (!justValuesPass)
          model_->primal(1);
        pinfo.postsolve(true);
        delete model_;
        model_ = saveModel;
        saveModel = NULL;
      }
    } else {
      // not feasible
      addAll = 1;
      presolve = 0;
      model_ = saveModel;
      saveModel = NULL;
      if (justValuesPass)
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 3);
#else
        model_->primal(2);
#endif
    }
    if (allowInfeasible) {
      CoinMemcpyN(saveRowUpper, nrows, model_->rowUpper());
      CoinMemcpyN(saveRowLower, nrows, model_->rowLower());
      delete[] saveRowUpper;
      delete[] saveRowLower;
      saveRowUpper = NULL;
      saveRowLower = NULL;
    }
    if (addAll < 2) {
      n = 0;
      if (!addAll) {
        /* could do scans to get a good number */
        iteration = 1;
        for (i = ordStart; i < ordEnd; i++) {
          if (whenUsed[i] >= iteration) {
            if (upper[i] - lower[i] < 1.0e-5 && saveUpper[i] - saveLower[i] > 1.0e-5) {
              n++;
              upper[i] = saveUpper[i];
              lower[i] = saveLower[i];
            }
          }
        }
      } else {
        for (i = ordStart; i < ordEnd; i++) {
          if (upper[i] - lower[i] < 1.0e-5 && saveUpper[i] - saveLower[i] > 1.0e-5) {
            n++;
            upper[i] = saveUpper[i];
            lower[i] = saveLower[i];
          }
        }
        delete[] saveUpper;
        delete[] saveLower;
        saveUpper = NULL;
        saveLower = NULL;
      }
#ifdef COIN_DEVELOP
      printf("Time so far %g, %d now added from previous iterations\n",
        CoinCpuTime() - startTime, n);
#endif
      if (justValuesPass)
        return;
      if (addAll)
        presolve = 0;
      if (presolve) {
        saveModel = model_;
        model_ = pinfo.presolvedModel(*model_, 1.0e-8, false, 5);
      } else {
        presolve = 0;
      }
      if (!wantVector) {
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 1);
#else
        model_->primal(1);
#endif
      } else {
        ClpMatrixBase *matrix = model_->clpMatrix();
        ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
        assert(clpMatrix);
        clpMatrix->makeSpecialColumnCopy();
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 1);
#else
        model_->primal(1);
#endif
        clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
        pinfo.postsolve(true);
        delete model_;
        model_ = saveModel;
        saveModel = NULL;
      }
      if (!addAll) {
        n = 0;
        for (i = ordStart; i < ordEnd; i++) {
          if (upper[i] - lower[i] < 1.0e-5 && saveUpper[i] - saveLower[i] > 1.0e-5) {
            n++;
            upper[i] = saveUpper[i];
            lower[i] = saveLower[i];
          }
        }
        delete[] saveUpper;
        delete[] saveLower;
        saveUpper = NULL;
        saveLower = NULL;
#ifdef COIN_DEVELOP
        printf("Time so far %g, %d now added from previous iterations\n",
          CoinCpuTime() - startTime, n);
#endif
      }
      if (presolve) {
        saveModel = model_;
        model_ = pinfo.presolvedModel(*model_, 1.0e-8, false, 5);
      } else {
        presolve = 0;
      }
      if (!wantVector) {
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 1);
#else
        model_->primal(1);
#endif
      } else {
        ClpMatrixBase *matrix = model_->clpMatrix();
        ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
        assert(clpMatrix);
        clpMatrix->makeSpecialColumnCopy();
#ifdef ABC_INHERIT
        model_->dealWithAbc(1, 1);
#else
        model_->primal(1);
#endif
        clpMatrix->releaseSpecialColumnCopy();
      }
      if (presolve) {
        pinfo.postsolve(true);
        delete model_;
        model_ = saveModel;
        saveModel = NULL;
      }
    }
#ifdef COIN_DEVELOP
    printf("Total time in crossover %g\n", CoinCpuTime() - startTime);
#endif
    delete[] saveUpper;
    delete[] saveLower;
  }
#ifdef FEB_TRY
  model_->setSpecialOptions(saveOptions);
  model_->setPerturbation(savePerturbation);
#endif
  return;
}
#endif
/*****************************************************************************/

// Default contructor
Idiot::Idiot()
{
  model_ = NULL;
  maxBigIts_ = 3;
  maxIts_ = 5;
  logLevel_ = 1;
  logFreq_ = 100;
  maxIts2_ = 100;
  djTolerance_ = 1e-1;
  mu_ = 1e-4;
  drop_ = 5.0;
  exitDrop_ = -1.0e20;
  muFactor_ = 0.3333;
  stopMu_ = 1e-12;
  smallInfeas_ = 1e-1;
  reasonableInfeas_ = 1e2;
  muAtExit_ = 1.0e31;
  strategy_ = 8;
  lambdaIterations_ = 0;
  checkFrequency_ = 100;
  whenUsed_ = NULL;
  majorIterations_ = 30;
  exitFeasibility_ = -1.0;
  dropEnoughFeasibility_ = 0.02;
  dropEnoughWeighted_ = 0.01;
  // adjust
  double nrows = 10000.0;
  int baseIts = static_cast< int >(sqrt(static_cast< double >(nrows)));
  baseIts = baseIts / 10;
  baseIts *= 10;
  maxIts2_ = 200 + baseIts + 5;
  maxIts2_ = 100;
  reasonableInfeas_ = static_cast< double >(nrows) * 0.05;
  lightWeight_ = 0;
}
// Constructor from model
Idiot::Idiot(OsiSolverInterface &model)
{
  model_ = &model;
  maxBigIts_ = 3;
  maxIts_ = 5;
  logLevel_ = 1;
  logFreq_ = 100;
  maxIts2_ = 100;
  djTolerance_ = 1e-1;
  mu_ = 1e-4;
  drop_ = 5.0;
  exitDrop_ = -1.0e20;
  muFactor_ = 0.3333;
  stopMu_ = 1e-12;
  smallInfeas_ = 1e-1;
  reasonableInfeas_ = 1e2;
  muAtExit_ = 1.0e31;
  strategy_ = 8;
  lambdaIterations_ = 0;
  checkFrequency_ = 100;
  whenUsed_ = NULL;
  majorIterations_ = 30;
  exitFeasibility_ = -1.0;
  dropEnoughFeasibility_ = 0.02;
  dropEnoughWeighted_ = 0.01;
  // adjust
  double nrows;
  if (model_)
    nrows = model_->getNumRows();
  else
    nrows = 10000.0;
  int baseIts = static_cast< int >(sqrt(static_cast< double >(nrows)));
  baseIts = baseIts / 10;
  baseIts *= 10;
  maxIts2_ = 200 + baseIts + 5;
  maxIts2_ = 100;
  reasonableInfeas_ = static_cast< double >(nrows) * 0.05;
  lightWeight_ = 0;
}
// Copy constructor.
Idiot::Idiot(const Idiot &rhs)
{
  model_ = rhs.model_;
  if (model_ && rhs.whenUsed_) {
    int numberColumns = model_->getNumCols();
    whenUsed_ = new int[numberColumns];
    CoinMemcpyN(rhs.whenUsed_, numberColumns, whenUsed_);
  } else {
    whenUsed_ = NULL;
  }
  djTolerance_ = rhs.djTolerance_;
  mu_ = rhs.mu_;
  drop_ = rhs.drop_;
  muFactor_ = rhs.muFactor_;
  stopMu_ = rhs.stopMu_;
  smallInfeas_ = rhs.smallInfeas_;
  reasonableInfeas_ = rhs.reasonableInfeas_;
  exitDrop_ = rhs.exitDrop_;
  muAtExit_ = rhs.muAtExit_;
  exitFeasibility_ = rhs.exitFeasibility_;
  dropEnoughFeasibility_ = rhs.dropEnoughFeasibility_;
  dropEnoughWeighted_ = rhs.dropEnoughWeighted_;
  maxBigIts_ = rhs.maxBigIts_;
  maxIts_ = rhs.maxIts_;
  majorIterations_ = rhs.majorIterations_;
  logLevel_ = rhs.logLevel_;
  logFreq_ = rhs.logFreq_;
  checkFrequency_ = rhs.checkFrequency_;
  lambdaIterations_ = rhs.lambdaIterations_;
  maxIts2_ = rhs.maxIts2_;
  strategy_ = rhs.strategy_;
  lightWeight_ = rhs.lightWeight_;
}
// Assignment operator. This copies the data
Idiot &
Idiot::operator=(const Idiot &rhs)
{
  if (this != &rhs) {
    delete[] whenUsed_;
    model_ = rhs.model_;
    if (model_ && rhs.whenUsed_) {
      int numberColumns = model_->getNumCols();
      whenUsed_ = new int[numberColumns];
      CoinMemcpyN(rhs.whenUsed_, numberColumns, whenUsed_);
    } else {
      whenUsed_ = NULL;
    }
    djTolerance_ = rhs.djTolerance_;
    mu_ = rhs.mu_;
    drop_ = rhs.drop_;
    muFactor_ = rhs.muFactor_;
    stopMu_ = rhs.stopMu_;
    smallInfeas_ = rhs.smallInfeas_;
    reasonableInfeas_ = rhs.reasonableInfeas_;
    exitDrop_ = rhs.exitDrop_;
    muAtExit_ = rhs.muAtExit_;
    exitFeasibility_ = rhs.exitFeasibility_;
    dropEnoughFeasibility_ = rhs.dropEnoughFeasibility_;
    dropEnoughWeighted_ = rhs.dropEnoughWeighted_;
    maxBigIts_ = rhs.maxBigIts_;
    maxIts_ = rhs.maxIts_;
    majorIterations_ = rhs.majorIterations_;
    logLevel_ = rhs.logLevel_;
    logFreq_ = rhs.logFreq_;
    checkFrequency_ = rhs.checkFrequency_;
    lambdaIterations_ = rhs.lambdaIterations_;
    maxIts2_ = rhs.maxIts2_;
    strategy_ = rhs.strategy_;
    lightWeight_ = rhs.lightWeight_;
  }
  return *this;
}
Idiot::~Idiot()
{
  delete[] whenUsed_;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
