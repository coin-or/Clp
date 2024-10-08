// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#if ABOCA_LITE
// 1 is not owner of abcState_
#define ABCSTATE_LITE 1
#endif
#include "ClpFactorization.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#endif
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpMatrixBase.hpp"
#ifndef SLIM_CLP
#include "ClpNetworkBasis.hpp"
#include "ClpNetworkMatrix.hpp"
//#define CHECK_NETWORK
#ifdef CHECK_NETWORK
const static bool doCheck = true;
#else
const static bool doCheck = false;
#endif
#endif
#ifndef CLP_FACTORIZATION_NEW_TIMING
#define CLP_FACTORIZATION_NEW_TIMING
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
#include "CoinTime.hpp"
double factorization_instrument(int type)
{
  static int times[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  static double startTime = 0.0;
  static double totalTimes[10] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  if (type < 0) {
    assert(!startTime);
    startTime = CoinCpuTime();
    return 0.0;
  } else if (type > 0) {
    times[type]++;
    double difference = CoinCpuTime() - startTime;
    totalTimes[type] += difference;
    startTime = 0.0;
    return difference;
  } else {
    // report
    const char *types[10] = {
      "", "fac=rhs_etc", "factorize", "replace", "update_FT",
      "update", "update_transpose", "gosparse", "getWeights!", "update2_FT"
    };
    double total = 0.0;
    for (int i = 1; i < 10; i++) {
      if (times[i]) {
        printf("%s was called %d times taking %g seconds\n",
          types[i], times[i], totalTimes[i]);
        total += totalTimes[i];
        times[i] = 0;
        totalTimes[i] = 0.0;
      }
    }
    return total;
  }
}
#endif
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################
#ifndef CLP_MULTIPLE_FACTORIZATIONS

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization()
  : CoinFactorization()
{
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization(const ClpFactorization &rhs,
  int dummyDenseIfSmaller)
  : CoinFactorization(rhs)
{
#ifndef SLIM_CLP
  if (rhs.networkBasis_)
    networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
  else
    networkBasis_ = NULL;
#endif
}

ClpFactorization::ClpFactorization(const CoinFactorization &rhs)
  : CoinFactorization(rhs)
{
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization()
{
#ifndef SLIM_CLP
  delete networkBasis_;
#endif
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization &rhs)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
  if (this != &rhs) {
    CoinFactorization::operator=(rhs);
#ifndef SLIM_CLP
    delete networkBasis_;
    if (rhs.networkBasis_)
      networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
    else
      networkBasis_ = NULL;
#endif
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(1);
#endif
  return *this;
}
#if 0
static unsigned int saveList[10000];
int numberSave = -1;
inline bool isDense(int i)
{
     return ((saveList[i>>5] >> (i & 31)) & 1) != 0;
}
inline void setDense(int i)
{
     unsigned int & value = saveList[i>>5];
     int bit = i & 31;
     value |= (1 << bit);
}
#endif
int ClpFactorization::factorize(ClpSimplex *model,
  int solveType, bool valuesPass)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
  ClpMatrixBase *matrix = model->clpMatrix();
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  if (!numberRows)
    return 0;
  // If too many compressions increase area
  if (numberPivots_ > 1 && numberCompressions_ * 10 > numberPivots_ + 10) {
    areaFactor_ *= 1.1;
  }
  //int numberPivots=numberPivots_;
#if 0
     if (model->algorithm() > 0)
          numberSave = -1;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_ || doCheck) {
#endif
    status_ = -99;
    int *COIN_RESTRICT pivotVariable = model->pivotVariable();
    //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
    while (status_ < -98) {

      int i;
      int numberBasic = 0;
      int numberRowBasic;
      // Move pivot variables across if they look good
      int *COIN_RESTRICT pivotTemp = model->rowArray(0)->getIndices();
      assert(!model->rowArray(0)->getNumElements());
      if (!matrix->rhsOffset(model)) {
#if 0
                    if (numberSave > 0) {
                         int nStill = 0;
                         int nAtBound = 0;
                         int nZeroDual = 0;
                         CoinIndexedVector * array = model->rowArray(3);
                         CoinIndexedVector * objArray = model->columnArray(1);
                         array->clear();
                         objArray->clear();
                         double *COIN_RESTRICT  cost = model->costRegion();
                         double tolerance = model->primalTolerance();
                         double offset = 0.0;
                         for (i = 0; i < numberRows; i++) {
                              int iPivot = pivotVariable[i];
                              if (iPivot < numberColumns && isDense(iPivot)) {
                                   if (model->getColumnStatus(iPivot) == ClpSimplex::basic) {
                                        nStill++;
                                        double value = model->solutionRegion()[iPivot];
                                        double dual = model->dualRowSolution()[i];
                                        double lower = model->lowerRegion()[iPivot];
                                        double upper = model->upperRegion()[iPivot];
                                        ClpSimplex::Status status;
                                        if (fabs(value - lower) < tolerance) {
                                             status = ClpSimplex::atLowerBound;
                                             nAtBound++;
                                        } else if (fabs(value - upper) < tolerance) {
                                             nAtBound++;
                                             status = ClpSimplex::atUpperBound;
                                        } else if (value > lower && value < upper) {
                                             status = ClpSimplex::superBasic;
                                        } else {
                                             status = ClpSimplex::basic;
                                        }
                                        if (status != ClpSimplex::basic) {
                                             if (model->getRowStatus(i) != ClpSimplex::basic) {
                                                  model->setColumnStatus(iPivot, ClpSimplex::atLowerBound);
                                                  model->setRowStatus(i, ClpSimplex::basic);
                                                  pivotVariable[i] = i + numberColumns;
                                                  model->dualRowSolution()[i] = 0.0;
                                                  model->djRegion(0)[i] = 0.0;
                                                  array->add(i, dual);
                                                  offset += dual * model->solutionRegion(0)[i];
                                             }
                                        }
                                        if (fabs(dual) < 1.0e-5)
                                             nZeroDual++;
                                   }
                              }
                         }
                         printf("out of %d dense, %d still in basis, %d at bound, %d with zero dual - offset %g\n",
                                numberSave, nStill, nAtBound, nZeroDual, offset);
                         if (array->getNumElements()) {
                              // modify costs
                              model->clpMatrix()->transposeTimes(model, 1.0, array, model->columnArray(0),
                                                                 objArray);
                              array->clear();
                              int n = objArray->getNumElements();
                              int *COIN_RESTRICT  indices = objArray->getIndices();
                              double *COIN_RESTRICT  elements = objArray->denseVector();
                              for (i = 0; i < n; i++) {
                                   int iColumn = indices[i];
                                   cost[iColumn] -= elements[iColumn];
                                   elements[iColumn] = 0.0;
                              }
                              objArray->setNumElements(0);
                         }
                    }
#endif
        // Seems to prefer things in order so quickest
        // way is to go though like this
        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic)
            pivotTemp[numberBasic++] = i;
        }
        numberRowBasic = numberBasic;
        /* Put column basic variables into pivotVariable
                       This is done by ClpMatrixBase to allow override for gub
                    */
        matrix->generalExpanded(model, 0, numberBasic);
      } else {
        // Long matrix - do a different way
        bool fullSearch = false;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot >= numberColumns) {
            pivotTemp[numberBasic++] = iPivot - numberColumns;
          }
        }
        numberRowBasic = numberBasic;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot < numberColumns) {
            if (iPivot >= 0) {
              pivotTemp[numberBasic++] = iPivot;
            } else {
              // not full basis
              fullSearch = true;
              break;
            }
          }
        }
        if (fullSearch) {
          // do slow way
          numberBasic = 0;
          for (i = 0; i < numberRows; i++) {
            if (model->getRowStatus(i) == ClpSimplex::basic)
              pivotTemp[numberBasic++] = i;
          }
          numberRowBasic = numberBasic;
          /* Put column basic variables into pivotVariable
                            This is done by ClpMatrixBase to allow override for gub
                         */
          matrix->generalExpanded(model, 0, numberBasic);
        }
      }
      if (numberBasic > model->maximumBasic()) {
#if 0 // ndef NDEBUG
                    printf("%d basic - should only be %d\n",
                           numberBasic, numberRows);
#endif
        // Take out some
        numberBasic = numberRowBasic;
        for (int i = 0; i < numberColumns; i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            if (numberBasic < numberRows)
              numberBasic++;
            else
              model->setColumnStatus(i, ClpSimplex::superBasic);
          }
        }
        numberBasic = numberRowBasic;
        matrix->generalExpanded(model, 0, numberBasic);
      }
#ifndef SLIM_CLP
      // see if matrix a network
#ifndef NO_RTTI
      ClpNetworkMatrix *networkMatrix = dynamic_cast< ClpNetworkMatrix * >(model->clpMatrix());
#else
      ClpNetworkMatrix *networkMatrix = NULL;
      if (model->clpMatrix()->type() == 11)
        networkMatrix = static_cast< ClpNetworkMatrix * >(model->clpMatrix());
#endif
      // If network - still allow ordinary factorization first time for laziness
      if (networkMatrix)
        biasLU_ = 0; // All to U if network
      //int saveMaximumPivots = maximumPivots();
      delete networkBasis_;
      networkBasis_ = NULL;
      if (networkMatrix && !doCheck)
        maximumPivots(1);
#endif
      //printf("L, U, R %d %d %d\n",numberElementsL(),numberElementsU(),numberElementsR());
      while (status_ == -99) {
        // maybe for speed will be better to leave as many regions as possible
        gutsOfDestructor();
        gutsOfInitialize(2);
        CoinBigIndex numberElements = numberRowBasic;

        // compute how much in basis

        int i;
        // can change for gub
        int numberColumnBasic = numberBasic - numberRowBasic;

        numberElements += matrix->countBasis(model,
          pivotTemp + numberRowBasic,
          numberRowBasic,
          numberColumnBasic);
        // and recompute as network side say different
        if (model->numberIterations())
          numberRowBasic = numberBasic - numberColumnBasic;
        numberElements = 3 * numberBasic + 3 * numberElements + 20000;
#if 0
                    // If iteration not zero then may be compressed
                    getAreas ( !model->numberIterations() ? numberRows : numberBasic,
                               numberRowBasic + numberColumnBasic, numberElements,
                               2 * numberElements );
#else
      getAreas(numberRows,
        numberRowBasic + numberColumnBasic, numberElements,
        2 * numberElements);
#endif
        //fill
        // Fill in counts so we can skip part of preProcess
        int *COIN_RESTRICT numberInRow = numberInRow_.array();
        int *COIN_RESTRICT numberInColumn = numberInColumn_.array();
        CoinZeroN(numberInRow, numberRows_ + 1);
        CoinZeroN(numberInColumn, maximumColumnsExtra_ + 1);
        double *COIN_RESTRICT elementU = elementU_.array();
        int *COIN_RESTRICT indexRowU = indexRowU_.array();
        CoinBigIndex *startColumnU = startColumnU_.array();
        for (i = 0; i < numberRowBasic; i++) {
          int iRow = pivotTemp[i];
          indexRowU[i] = iRow;
          startColumnU[i] = i;
          elementU[i] = slackValue_;
          numberInRow[iRow] = 1;
          numberInColumn[i] = 1;
        }
        startColumnU[numberRowBasic] = numberRowBasic;
        // can change for gub so redo
        numberColumnBasic = numberBasic - numberRowBasic;
        matrix->fillBasis(model,
          pivotTemp + numberRowBasic,
          numberColumnBasic,
          indexRowU,
          startColumnU + numberRowBasic,
          numberInRow,
          numberInColumn + numberRowBasic,
          elementU);
#if 0
                    {
                         printf("%d row basic, %d column basic\n", numberRowBasic, numberColumnBasic);
                         for (int i = 0; i < numberElements; i++)
                              printf("row %d col %d value %g\n", indexRowU_.array()[i], indexColumnU_[i],
                                     elementU_.array()[i]);
                    }
#endif
        // recompute number basic
        numberBasic = numberRowBasic + numberColumnBasic;
        if (numberBasic)
          numberElements = startColumnU[numberBasic - 1]
            + numberInColumn[numberBasic - 1];
        else
          numberElements = 0;
        lengthU_ = numberElements;
        //saveFactorization("dump.d");
        if (biasLU_ >= 3 || numberRows_ != numberColumns_)
          preProcess(2);
        else
          preProcess(3); // no row copy
        factor();
        if (status_ == -99) {
          // get more memory
          areaFactor(2.0 * areaFactor());
        } else if (status_ == -1 && model->numberIterations() == 0 && denseThreshold_) {
          // Round again without dense
          denseThreshold_ = 0;
          status_ = -99;
        }
      }
      // If we get here status is 0 or -1
      if (status_ == 0) {
        // We may need to tamper with order and redo - e.g. network with side
        int useNumberRows = numberRows;
        // **** we will also need to add test in dual steepest to do
        // as we do for network
        matrix->generalExpanded(model, 12, useNumberRows);
        const int *permuteBack = permuteBack_.array();
        const int *back = pivotColumnBack_.array();
        //int * pivotTemp = pivotColumn_.array();
        //ClpDisjointCopyN ( pivotVariable, numberRows , pivotTemp  );
        // Redo pivot order
        for (i = 0; i < numberRowBasic; i++) {
          int k = pivotTemp[i];
          // so rowIsBasic[k] would be permuteBack[back[i]]
          pivotVariable[permuteBack[back[i]]] = k + numberColumns;
        }
        for (; i < useNumberRows; i++) {
          int k = pivotTemp[i];
          // so rowIsBasic[k] would be permuteBack[back[i]]
          pivotVariable[permuteBack[back[i]]] = k;
        }
#if 0
                    if (numberSave >= 0) {
                         numberSave = numberDense_;
                         memset(saveList, 0, ((numberRows_ + 31) >> 5)*sizeof(int));
                         for (i = numberRows_ - numberSave; i < numberRows_; i++) {
                              int k = pivotTemp[pivotColumn_.array()[i]];
                              setDense(k);
                         }
                    }
#endif
        // Set up permutation vector
        // these arrays start off as copies of permute
        // (and we could use permute_ instead of pivotColumn (not back though))
        ClpDisjointCopyN(permute_.array(), useNumberRows, pivotColumn_.array());
        ClpDisjointCopyN(permuteBack_.array(), useNumberRows, pivotColumnBack_.array());
#ifndef SLIM_CLP
        if (networkMatrix) {
          maximumPivots(std::max(2000, maximumPivots()));
          // redo arrays
          for (int iRow = 0; iRow < 4; iRow++) {
            int length = model->numberRows() + maximumPivots();
            if (iRow == 3 || model->objectiveAsObject()->type() > 1)
              length += model->numberColumns();
            model->rowArray(iRow)->reserve(length);
          }
          // create network factorization
          if (doCheck)
            delete networkBasis_; // temp
          networkBasis_ = new ClpNetworkBasis(model, numberRows_,
            pivotRegion_.array(),
            permuteBack_.array(),
            startColumnU_.array(),
            numberInColumn_.array(),
            indexRowU_.array(),
            elementU_.array());
          // kill off arrays in ordinary factorization
          if (!doCheck) {
            gutsOfDestructor();
            // but make sure numberRows_ set
            numberRows_ = model->numberRows();
            status_ = 0;
#if 0
                              // but put back permute arrays so odd things will work
                              int numberRows = model->numberRows();
                              pivotColumnBack_ = new int [numberRows];
                              permute_ = new int [numberRows];
                              int i;
                              for (i = 0; i < numberRows; i++) {
                                   pivotColumnBack_[i] = i;
                                   permute_[i] = i;
                              }
#endif
          }
        } else {
#endif
          // See if worth going sparse and when
          if (numberFtranCounts_ > 100) {
            ftranCountInput_ = std::max(ftranCountInput_, 1.0);
            ftranAverageAfterL_ = std::max(ftranCountAfterL_ / ftranCountInput_, 1.0);
            ftranAverageAfterR_ = std::max(ftranCountAfterR_ / ftranCountAfterL_, 1.0);
            ftranAverageAfterU_ = std::max(ftranCountAfterU_ / ftranCountAfterR_, 1.0);
            if (btranCountInput_ && btranCountAfterU_ && btranCountAfterR_) {
              btranAverageAfterU_ = std::max(btranCountAfterU_ / btranCountInput_, 1.0);
              btranAverageAfterR_ = std::max(btranCountAfterR_ / btranCountAfterU_, 1.0);
              btranAverageAfterL_ = std::max(btranCountAfterL_ / btranCountAfterR_, 1.0);
            } else {
              // we have not done any useful btrans (values pass?)
              btranAverageAfterU_ = 1.0;
              btranAverageAfterR_ = 1.0;
              btranAverageAfterL_ = 1.0;
            }
          }
          // scale back

          ftranCountInput_ *= 0.8;
          ftranCountAfterL_ *= 0.8;
          ftranCountAfterR_ *= 0.8;
          ftranCountAfterU_ *= 0.8;
          btranCountInput_ *= 0.8;
          btranCountAfterU_ *= 0.8;
          btranCountAfterR_ *= 0.8;
          btranCountAfterL_ *= 0.8;
#ifndef SLIM_CLP
        }
#endif
      } else if (status_ == -1 && (solveType == 0 || solveType == 2)) {
        // This needs redoing as it was merged coding - does not need array
        int numberTotal = numberRows + numberColumns;
        int *COIN_RESTRICT isBasic = new int[numberTotal];
        int *COIN_RESTRICT rowIsBasic = isBasic + numberColumns;
        int *COIN_RESTRICT columnIsBasic = isBasic;
        for (i = 0; i < numberTotal; i++)
          isBasic[i] = -1;
        for (i = 0; i < numberRowBasic; i++) {
          int iRow = pivotTemp[i];
          rowIsBasic[iRow] = 1;
        }
        for (; i < numberBasic; i++) {
          int iColumn = pivotTemp[i];
          columnIsBasic[iColumn] = 1;
        }
        numberBasic = 0;
        for (i = 0; i < numberRows; i++)
          pivotVariable[i] = -1;
        // mark as basic or non basic
        const int *pivotColumn = pivotColumn_.array();
        for (i = 0; i < numberRows; i++) {
          if (rowIsBasic[i] >= 0) {
            if (pivotColumn[numberBasic] >= 0) {
              rowIsBasic[i] = pivotColumn[numberBasic];
            } else {
              rowIsBasic[i] = -1;
              model->setRowStatus(i, ClpSimplex::superBasic);
            }
            numberBasic++;
          }
        }
        for (i = 0; i < numberColumns; i++) {
          if (columnIsBasic[i] >= 0) {
            if (pivotColumn[numberBasic] >= 0)
              columnIsBasic[i] = pivotColumn[numberBasic];
            else
              columnIsBasic[i] = -1;
            numberBasic++;
          }
        }
        // leave pivotVariable in useful form for cleaning basis
        int *COIN_RESTRICT pivotVariable = model->pivotVariable();
        for (i = 0; i < numberRows; i++) {
          pivotVariable[i] = -1;
        }

        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic) {
            int iPivot = rowIsBasic[i];
            if (iPivot >= 0)
              pivotVariable[iPivot] = i + numberColumns;
          }
        }
        for (i = 0; i < numberColumns; i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            int iPivot = columnIsBasic[i];
            if (iPivot >= 0)
              pivotVariable[iPivot] = i;
          }
        }
        delete[] isBasic;
        double *COIN_RESTRICT columnLower = model->lowerRegion();
        double *COIN_RESTRICT columnUpper = model->upperRegion();
        double *COIN_RESTRICT columnActivity = model->solutionRegion();
        double *COIN_RESTRICT rowLower = model->lowerRegion(0);
        double *COIN_RESTRICT rowUpper = model->upperRegion(0);
        double *COIN_RESTRICT rowActivity = model->solutionRegion(0);
        //redo basis - first take ALL columns out
        int iColumn;
        double largeValue = model->largeValue();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) == ClpSimplex::basic) {
            // take out
            if (!valuesPass) {
              double lower = columnLower[iColumn];
              double upper = columnUpper[iColumn];
              double value = columnActivity[iColumn];
              if (lower > -largeValue || upper < largeValue) {
                if (fabs(value - lower) < fabs(value - upper)) {
                  model->setColumnStatus(iColumn, ClpSimplex::atLowerBound);
                  columnActivity[iColumn] = lower;
                } else {
                  model->setColumnStatus(iColumn, ClpSimplex::atUpperBound);
                  columnActivity[iColumn] = upper;
                }
              } else {
                model->setColumnStatus(iColumn, ClpSimplex::isFree);
              }
            } else {
              model->setColumnStatus(iColumn, ClpSimplex::superBasic);
            }
          }
        }
        int iRow;
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iSequence = pivotVariable[iRow];
          if (iSequence >= 0) {
            // basic
            if (iSequence >= numberColumns) {
              // slack in - leave
              //assert (iSequence-numberColumns==iRow);
            } else {
              assert(model->getRowStatus(iRow) != ClpSimplex::basic);
              // put back structural
              model->setColumnStatus(iSequence, ClpSimplex::basic);
            }
          } else {
            // put in slack
            model->setRowStatus(iRow, ClpSimplex::basic);
          }
        }
        // Put back any key variables for gub (status_ not touched)
        matrix->generalExpanded(model, 1, status_);
        // signal repeat
        status_ = -99;
        // set fixed if they are
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (model->getRowStatus(iRow) != ClpSimplex::basic) {
            if (rowLower[iRow] == rowUpper[iRow]) {
              rowActivity[iRow] = rowLower[iRow];
              model->setRowStatus(iRow, ClpSimplex::isFixed);
            }
          }
        }
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) != ClpSimplex::basic) {
            if (columnLower[iColumn] == columnUpper[iColumn]) {
              columnActivity[iColumn] = columnLower[iColumn];
              model->setColumnStatus(iColumn, ClpSimplex::isFixed);
            }
          }
        }
      }
    }
#ifndef SLIM_CLP
  } else {
    // network - fake factorization - do nothing
    status_ = 0;
    numberPivots_ = 0;
  }
#endif
#ifndef SLIM_CLP
  if (!status_) {
    // take out part if quadratic
    if (model->algorithm() == 2) {
      ClpObjective *obj = model->objectiveAsObject();
#ifndef NDEBUG
      ClpQuadraticObjective *quadraticObj = (dynamic_cast< ClpQuadraticObjective * >(obj));
      assert(quadraticObj);
#else
      ClpQuadraticObjective *quadraticObj = (static_cast< ClpQuadraticObjective * >(obj));
#endif
      CoinPackedMatrix *quadratic = quadraticObj->quadraticObjective();
      int numberXColumns = quadratic->getNumCols();
      assert(numberXColumns < numberColumns);
      int base = numberColumns - numberXColumns;
      int *COIN_RESTRICT which = new int[numberXColumns];
      int *COIN_RESTRICT pivotVariable = model->pivotVariable();
      int *COIN_RESTRICT permute = pivotColumn();
      int i;
      int n = 0;
      for (i = 0; i < numberRows; i++) {
        int iSj = pivotVariable[i] - base;
        if (iSj >= 0 && iSj < numberXColumns)
          which[n++] = permute[i];
      }
      if (n)
        emptyRows(n, which);
      delete[] which;
    }
  }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(2);
#endif
  return status_;
}
/* Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   If checkBeforeModifying is true will do all accuracy checks
   before modifying factorization.  Whether to set this depends on
   speed considerations.  You could just do this on first iteration
   after factorization and thereafter re-factorize
   partial update already in U */
int ClpFactorization::replaceColumn(const ClpSimplex *model,
  CoinIndexedVector *regionSparse,
  CoinIndexedVector *tableauColumn,
  int pivotRow,
  double pivotCheck,
  bool checkBeforeModifying,
  double acceptablePivot)
{
  int returnCode;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    // see if FT
    if (doForrestTomlin_) {
      returnCode = CoinFactorization::replaceColumn(regionSparse,
        pivotRow,
        pivotCheck,
        checkBeforeModifying,
        acceptablePivot);
    } else {
      returnCode = CoinFactorization::replaceColumnPFI(tableauColumn,
        pivotRow, pivotCheck); // Note array
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(3);
#endif

#ifndef SLIM_CLP
  } else {
    if (doCheck) {
      returnCode = CoinFactorization::replaceColumn(regionSparse,
        pivotRow,
        pivotCheck,
        checkBeforeModifying,
        acceptablePivot);
      networkBasis_->replaceColumn(regionSparse,
        pivotRow);
    } else {
      // increase number of pivots
      numberPivots_++;
      returnCode = networkBasis_->replaceColumn(regionSparse,
        pivotRow);
    }
  }
#endif
  return returnCode;
}

/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnFT(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2)
{
#ifdef CLP_DEBUG
  regionSparse->checkClear();
#endif
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    collectStatistics_ = true;
    int returnCode = CoinFactorization::updateColumnFT(regionSparse,
      regionSparse2);
    collectStatistics_ = false;
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(4);
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumnFT(regionSparse,
      regionSparse2);
    networkBasis_->updateColumn(regionSparse, save, -1);
    int i;
    double *array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, numberRows_ * sizeof(double));
    double *array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      for (i = 0; i < numberRows_; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse, regionSparse2, -1);
    return 1;
#endif
  }
#endif
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumn(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  bool noPermute) const
{
#ifdef CLP_DEBUG
  if (!noPermute)
    regionSparse->checkClear();
#endif
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    collectStatistics_ = true;
    int returnCode = CoinFactorization::updateColumn(regionSparse,
      regionSparse2,
      noPermute);
    collectStatistics_ = false;
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(5);
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumn(regionSparse,
      regionSparse2,
      noPermute);
    networkBasis_->updateColumn(regionSparse, save, -1);
    int i;
    double *array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, numberRows_ * sizeof(double));
    double *array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      for (i = 0; i < numberRows_; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse, regionSparse2, -1);
    return 1;
#endif
  }
#endif
}
/* Updates one column (FTRAN) from region2
   Tries to do FT update
   number returned is negative if no room.
   Also updates region3
   region1 starts as zero and is zero at end */
int ClpFactorization::updateTwoColumnsFT(CoinIndexedVector *regionSparse1,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3,
  bool noPermuteRegion3)
{
  int returnCode = updateColumnFT(regionSparse1, regionSparse2);
  updateColumn(regionSparse1, regionSparse3, noPermuteRegion3);
  return returnCode;
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnForDebug(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  bool noPermute) const
{
  if (!noPermute)
    regionSparse->checkClear();
  if (!numberRows_)
    return 0;
  collectStatistics_ = false;
  int returnCode = CoinFactorization::updateColumn(regionSparse,
    regionSparse2,
    noPermute);
  return returnCode;
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2) const
{
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    collectStatistics_ = true;
    int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
      regionSparse2);
    collectStatistics_ = false;
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(6);
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
      regionSparse2);
    networkBasis_->updateColumnTranspose(regionSparse, save);
    int i;
    double *array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, numberRows_ * sizeof(double));
    double *array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      for (i = 0; i < numberRows_; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    return networkBasis_->updateColumnTranspose(regionSparse, regionSparse2);
#endif
  }
#endif
}
/* Updates two columns (BTRAN) from regionSparse2 and 3
   regionSparse starts as zero and is zero at end 
   Note - if regionSparse2 packed on input - will be packed on output - same for 3
*/
void ClpFactorization::updateTwoColumnsTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3) const
{
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    collectStatistics_ = true;
    CoinFactorization::updateTwoColumnsTranspose(regionSparse,
      regionSparse2, regionSparse3);
    collectStatistics_ = false;
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(6);
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
    factorization__->updateColumnTranspose(regionSparse, regionSparse2);
    factorization__->updateColumnTranspose(regionSparse, regionSparse3);
  }
#endif
}
/* makes a row copy of L for speed and to allow very sparse problems */
void ClpFactorization::goSparse()
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  if (!networkBasis_)
#endif
    CoinFactorization::goSparse();
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(7);
#endif
}
// Cleans up i.e. gets rid of network basis
void ClpFactorization::cleanUp()
{
#ifndef SLIM_CLP
  delete networkBasis_;
  networkBasis_ = NULL;
#endif
  resetStatistics();
}
/// Says whether to redo pivot order
bool ClpFactorization::needToReorder() const
{
#ifdef CHECK_NETWORK
  return true;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_)
#endif
    return true;
#ifndef SLIM_CLP
  else
    return false;
#endif
}
// Get weighted row list
void ClpFactorization::getWeights(int *COIN_RESTRICT weights) const
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  if (networkBasis_) {
    // Network - just unit
    for (int i = 0; i < numberRows_; i++)
      weights[i] = 1;
    return;
  }
#endif
  int *COIN_RESTRICT numberInRow = numberInRow_.array();
  int *COIN_RESTRICT numberInColumn = numberInColumn_.array();
  int *COIN_RESTRICT permuteBack = pivotColumnBack_.array();
  int *COIN_RESTRICT indexRowU = indexRowU_.array();
  const CoinBigIndex *startColumnU = startColumnU_.array();
  const CoinBigIndex *startRowL = startRowL_.array();
  if (!startRowL || !numberInRow_.array()) {
    int *COIN_RESTRICT temp = new int[numberRows_];
    memset(temp, 0, numberRows_ * sizeof(int));
    int i;
    for (i = 0; i < numberRows_; i++) {
      // one for pivot
      temp[i]++;
      CoinBigIndex j;
      for (j = startColumnU[i]; j < startColumnU[i] + numberInColumn[i]; j++) {
        int iRow = indexRowU[j];
        temp[iRow]++;
      }
    }
    CoinBigIndex *startColumnL = startColumnL_.array();
    int *COIN_RESTRICT indexRowL = indexRowL_.array();
    for (i = baseL_; i < baseL_ + numberL_; i++) {
      CoinBigIndex j;
      for (j = startColumnL[i]; j < startColumnL[i + 1]; j++) {
        int iRow = indexRowL[j];
        temp[iRow]++;
      }
    }
    for (i = 0; i < numberRows_; i++) {
      int number = temp[i];
      int iPermute = permuteBack[i];
      weights[iPermute] = number;
    }
    delete[] temp;
  } else {
    int i;
    for (i = 0; i < numberRows_; i++) {
      int number = startRowL[i + 1] - startRowL[i] + numberInRow[i] + 1;
      //number = startRowL[i+1]-startRowL[i]+1;
      //number = numberInRow[i]+1;
      int iPermute = permuteBack[i];
      weights[iPermute] = number;
    }
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(8);
#endif
}
#else
#if COINUTILS_BIGINDEX_IS_INT
// This one allows multiple factorizations
#if CLP_MULTIPLE_FACTORIZATIONS == 1
typedef CoinDenseFactorization CoinOtherFactorization;
typedef CoinOslFactorization CoinOtherFactorization;
#elif CLP_MULTIPLE_FACTORIZATIONS == 2
#include "CoinSimpFactorization.hpp"
typedef CoinSimpFactorization CoinOtherFactorization;
typedef CoinOslFactorization CoinOtherFactorization;
#elif CLP_MULTIPLE_FACTORIZATIONS == 3
#include "CoinSimpFactorization.hpp"
#define CoinOslFactorization CoinDenseFactorization
#elif CLP_MULTIPLE_FACTORIZATIONS == 4
#include "CoinSimpFactorization.hpp"
//#define CoinOslFactorization CoinDenseFactorization
#include "CoinOslFactorization.hpp"
#endif
#endif

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization()
{
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
  //coinFactorizationA_ = NULL;
  coinFactorizationA_ = new CoinFactorization();
  coinFactorizationB_ = NULL;
  //coinFactorizationB_ = new CoinOtherFactorization();
  forceB_ = 0;
  goOslThreshold_ = -1;
  goDenseThreshold_ = -1;
  goSmallThreshold_ = -1;
  doStatistics_ = true;
  memset(&shortestAverage_, 0, 3 * (sizeof(double) + sizeof(int)));
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization(const ClpFactorization &rhs,
  int denseIfSmaller)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  if (rhs.networkBasis_)
    networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
  else
    networkBasis_ = NULL;
#endif
  forceB_ = rhs.forceB_;
  goOslThreshold_ = rhs.goOslThreshold_;
  goDenseThreshold_ = rhs.goDenseThreshold_;
  goSmallThreshold_ = rhs.goSmallThreshold_;
  doStatistics_ = rhs.doStatistics_;
  int goDense = 0;
#ifdef CLP_REUSE_ETAS
  model_ = rhs.model_;
#endif
#if COINUTILS_BIGINDEX_IS_INT
  if (denseIfSmaller > 0 && denseIfSmaller <= goDenseThreshold_) {
    CoinDenseFactorization *denseR = dynamic_cast< CoinDenseFactorization * >(rhs.coinFactorizationB_);
    if (!denseR)
      goDense = 1;
  }
  if (denseIfSmaller > 0 && !rhs.coinFactorizationB_) {
    if (denseIfSmaller <= goDenseThreshold_)
      goDense = 1;
    else if (denseIfSmaller <= goSmallThreshold_)
      goDense = 2;
    else if (denseIfSmaller <= goOslThreshold_)
      goDense = 3;
  } else if (denseIfSmaller < 0) {
    if (-denseIfSmaller <= goDenseThreshold_)
      goDense = 1;
    else if (-denseIfSmaller <= goSmallThreshold_)
      goDense = 2;
    else if (-denseIfSmaller <= goOslThreshold_)
      goDense = 3;
  }
#endif
  if (rhs.coinFactorizationA_ && !goDense)
    coinFactorizationA_ = new CoinFactorization(*(rhs.coinFactorizationA_));
  else
    coinFactorizationA_ = NULL;
  if (rhs.coinFactorizationB_ && (denseIfSmaller >= 0 || !goDense))
    coinFactorizationB_ = rhs.coinFactorizationB_->clone();
  else
    coinFactorizationB_ = NULL;
#if COINUTILS_BIGINDEX_IS_INT
  if (goDense) {
    delete coinFactorizationB_;
    if (goDense == 1)
      coinFactorizationB_ = new CoinDenseFactorization();
    else if (goDense == 2)
      coinFactorizationB_ = new CoinSimpFactorization();
    else
      coinFactorizationB_ = new CoinOslFactorization();
    if (rhs.coinFactorizationA_) {
      coinFactorizationB_->maximumPivots(rhs.coinFactorizationA_->maximumPivots());
      coinFactorizationB_->pivotTolerance(rhs.coinFactorizationA_->pivotTolerance());
      coinFactorizationB_->zeroTolerance(rhs.coinFactorizationA_->zeroTolerance());
    } else {
      assert(coinFactorizationB_);
      coinFactorizationB_->maximumPivots(rhs.coinFactorizationB_->maximumPivots());
      coinFactorizationB_->pivotTolerance(rhs.coinFactorizationB_->pivotTolerance());
      coinFactorizationB_->zeroTolerance(rhs.coinFactorizationB_->zeroTolerance());
    }
  }
#endif
  assert(!coinFactorizationA_ || !coinFactorizationB_);
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(1);
#endif
  memcpy(&shortestAverage_, &rhs.shortestAverage_, 3 * (sizeof(double) + sizeof(int)));
}

ClpFactorization::ClpFactorization(const CoinFactorization &rhs)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
  coinFactorizationA_ = new CoinFactorization(rhs);
  coinFactorizationB_ = NULL;
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(1);
#endif
  forceB_ = 0;
  goOslThreshold_ = -1;
  goDenseThreshold_ = -1;
  goSmallThreshold_ = -1;
  doStatistics_ = true;
  assert(!coinFactorizationA_ || !coinFactorizationB_);
  memset(&shortestAverage_, 0, 3 * (sizeof(double) + sizeof(int)));
}

ClpFactorization::ClpFactorization(const CoinOtherFactorization &rhs)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
  coinFactorizationA_ = NULL;
  coinFactorizationB_ = rhs.clone();
  //coinFactorizationB_ = new CoinOtherFactorization(rhs);
  forceB_ = 0;
  goOslThreshold_ = -1;
  goDenseThreshold_ = -1;
  goSmallThreshold_ = -1;
  doStatistics_ = true;
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(1);
#endif
  assert(!coinFactorizationA_ || !coinFactorizationB_);
  memset(&shortestAverage_, 0, 3 * (sizeof(double) + sizeof(int)));
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization()
{
#ifndef SLIM_CLP
  delete networkBasis_;
#endif
  delete coinFactorizationA_;
  delete coinFactorizationB_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization &rhs)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
  if (this != &rhs) {
#ifndef SLIM_CLP
    delete networkBasis_;
    if (rhs.networkBasis_)
      networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
    else
      networkBasis_ = NULL;
#endif
    forceB_ = rhs.forceB_;
#ifdef CLP_REUSE_ETAS
    model_ = rhs.model_;
#endif
    goOslThreshold_ = rhs.goOslThreshold_;
    goDenseThreshold_ = rhs.goDenseThreshold_;
    goSmallThreshold_ = rhs.goSmallThreshold_;
    doStatistics_ = rhs.doStatistics_;
    memcpy(&shortestAverage_, &rhs.shortestAverage_, 3 * (sizeof(double) + sizeof(int)));
    if (rhs.coinFactorizationA_) {
      if (coinFactorizationA_)
        *coinFactorizationA_ = *(rhs.coinFactorizationA_);
      else
        coinFactorizationA_ = new CoinFactorization(*rhs.coinFactorizationA_);
    } else {
      delete coinFactorizationA_;
      coinFactorizationA_ = NULL;
    }
#if COINUTILS_BIGINDEX_IS_INT
    if (rhs.coinFactorizationB_) {
      if (coinFactorizationB_) {
        CoinDenseFactorization *denseR = dynamic_cast< CoinDenseFactorization * >(rhs.coinFactorizationB_);
        CoinDenseFactorization *dense = dynamic_cast< CoinDenseFactorization * >(coinFactorizationB_);
        CoinOslFactorization *oslR = dynamic_cast< CoinOslFactorization * >(rhs.coinFactorizationB_);
        CoinOslFactorization *osl = dynamic_cast< CoinOslFactorization * >(coinFactorizationB_);
        CoinSimpFactorization *simpR = dynamic_cast< CoinSimpFactorization * >(rhs.coinFactorizationB_);
        CoinSimpFactorization *simp = dynamic_cast< CoinSimpFactorization * >(coinFactorizationB_);
        if (dense && denseR) {
          *dense = *denseR;
        } else if (osl && oslR) {
          *osl = *oslR;
        } else if (simp && simpR) {
          *simp = *simpR;
        } else {
          delete coinFactorizationB_;
          coinFactorizationB_ = rhs.coinFactorizationB_->clone();
        }
      } else {
        coinFactorizationB_ = rhs.coinFactorizationB_->clone();
      }
    } else {
      delete coinFactorizationB_;
      coinFactorizationB_ = NULL;
    }
#endif
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(1);
#endif
  assert(!coinFactorizationA_ || !coinFactorizationB_);
  return *this;
}
// Go over to dense code
void ClpFactorization::goDenseOrSmall(int numberRows)
{
  if (!forceB_) {
    if (numberRows <= goDenseThreshold_) {
      delete coinFactorizationA_;
      delete coinFactorizationB_;
      coinFactorizationA_ = NULL;
      coinFactorizationB_ = new CoinDenseFactorization();
      //printf("going dense\n");
#if COINUTILS_BIGINDEX_IS_INT
    } else if (numberRows <= goSmallThreshold_) {
      delete coinFactorizationA_;
      delete coinFactorizationB_;
      coinFactorizationA_ = NULL;
      coinFactorizationB_ = new CoinSimpFactorization();
      //printf("going small\n");
    } else if (numberRows <= goOslThreshold_) {
      delete coinFactorizationA_;
      delete coinFactorizationB_;
      coinFactorizationA_ = NULL;
      coinFactorizationB_ = new CoinOslFactorization();
      //printf("going small\n");
#endif
    }
  }
  assert(!coinFactorizationA_ || !coinFactorizationB_);
}
// If nonzero force use of 1,dense 2,small 3,osl
void ClpFactorization::forceOtherFactorization(int which)
{
  delete coinFactorizationB_;
  forceB_ = 0;
  coinFactorizationB_ = NULL;
  if (which > 0 && which < 4) {
    delete coinFactorizationA_;
    coinFactorizationA_ = NULL;
    forceB_ = which;
    switch (which) {
    case 1:
      coinFactorizationB_ = new CoinDenseFactorization();
      goDenseThreshold_ = COIN_INT_MAX;
      break;
#if COINUTILS_BIGINDEX_IS_INT
    case 2:
      coinFactorizationB_ = new CoinSimpFactorization();
      goSmallThreshold_ = COIN_INT_MAX;
      break;
    case 3:
      coinFactorizationB_ = new CoinOslFactorization();
      goOslThreshold_ = COIN_INT_MAX;
      break;
#endif
    }
  } else if (!coinFactorizationA_) {
    coinFactorizationA_ = new CoinFactorization();
    goOslThreshold_ = -1;
    goDenseThreshold_ = -1;
    goSmallThreshold_ = -1;
  }
}
#ifdef CLP_FACTORIZATION_NEW_TIMING
#ifdef CLP_FACTORIZATION_INSTRUMENT
extern double externalTimeStart;
extern double timeInFactorize;
extern double timeInUpdate;
extern double timeInFactorizeFake;
extern double timeInUpdateFake1;
extern double timeInUpdateFake2;
extern double timeInUpdateTranspose;
extern double timeInUpdateFT;
extern double timeInUpdateTwoFT;
extern double timeInReplace;
extern int numberUpdate;
extern int numberUpdateTranspose;
extern int numberUpdateFT;
extern int numberUpdateTwoFT;
extern int numberReplace;
extern int currentLengthR;
extern int currentLengthU;
extern int currentTakeoutU;
extern double averageLengthR;
extern double averageLengthL;
extern double averageLengthU;
extern double scaledLengthDense;
extern double scaledLengthDenseSquared;
extern double scaledLengthL;
extern double scaledLengthR;
extern double scaledLengthU;
extern int startLengthU;
extern int endLengthU;
extern int endLengthU2;
extern int numberAdded;
static int average[3];
static double shortest;
static int ifPrint;
#else
#ifdef CLP_STATIC
static int endLengthU_;
#endif
#endif
#ifdef CLP_STATIC
static double shortestAverage_;
static double totalInR_ = 0.0;
static double totalInIncreasingU_ = 0.0;
//static int lastR=0;
//static int lastU=0;
static int lastNumberPivots_ = 0;
static int effectiveStartNumberU_ = 0;
#else
//#define shortestAverage shortestAverage_
//#define totalInR totalInR_
//#define totalInIncreasingU totalInIncreasingU_
//#define lastNumberPivots lastNumberPivots_
//#define effectiveStartNumberU effectiveStartNumberU_
//#define endLengthU endLengthU_
#endif
#ifdef CLP_USEFUL_PRINTOUT
static bool readTwiddle = false;
static double weightIncU = 1.0;
static double weightR = 2.0;
static double weightRest = 1.0;
static double weightFactL = 30.0;
static double weightFactDense = 0.1;
static double weightNrows = 10.0;
static double increaseNeeded = 1.1;
static double constWeightIterate = 1.0;
static double weightNrowsIterate = 3.0;
#else
#define weightIncU 1.0
#define weightR 2.0
#define weightRest 1.0
#define weightFactL 30.0
#define weightFactDense 0.1
#define weightNrows 10.0
#define increaseNeeded 1.1
#define constWeightIterate 1.0
#define weightNrowsIterate 3.0
#endif
bool ClpFactorization::timeToRefactorize() const
{
  if (coinFactorizationA_) {
    bool reFactor = (coinFactorizationA_->pivots() * 3 > coinFactorizationA_->maximumPivots() * 2 && coinFactorizationA_->numberElementsR() * 3 > (coinFactorizationA_->numberElementsL() + coinFactorizationA_->numberElementsU()) * 2 + 1000 && !coinFactorizationA_->numberDense());
    reFactor = false;
    bool reFactor3 = false;
    int numberPivots = coinFactorizationA_->pivots();
    //if (coinFactorizationA_->pivots()<2)
    if (numberPivots > lastNumberPivots_) {
      if (!lastNumberPivots_) {
        //lastR=0;
        //lastU=endLengthU;
        totalInR_ = 0.0;
        totalInIncreasingU_ = 0.0;
        shortestAverage_ = COIN_DBL_MAX;
#ifdef CLP_USEFUL_PRINTOUT
        if (!readTwiddle) {
          readTwiddle = true;
          char *environ = getenv("CLP_TWIDDLE");
          if (environ) {
            sscanf(environ, "%lg %lg %lg %lg %lg %lg %lg %lg %lg",
              &weightIncU, &weightR, &weightRest, &weightFactL,
              &weightFactDense, &weightNrows, &increaseNeeded,
              &constWeightIterate, &weightNrowsIterate);
          }
          printf("weightIncU %g, weightR %g, weightRest %g, weightFactL %g, weightFactDense %g, weightNrows %g increaseNeeded %g constWeightIterate %g weightNrowsIterate %g\n",
            weightIncU, weightR, weightRest, weightFactL,
            weightFactDense, weightNrows, increaseNeeded,
            constWeightIterate, weightNrowsIterate);
        }
#endif
      }
      lastNumberPivots_ = numberPivots;
      int numberDense = coinFactorizationA_->numberDense();
      double nnd = numberDense * numberDense;
      int lengthL = coinFactorizationA_->numberElementsL();
      int lengthR = coinFactorizationA_->numberElementsR();
      int numberRows = coinFactorizationA_->numberRows();
      int lengthU = coinFactorizationA_->numberElementsU() - (numberRows - numberDense);
      totalInR_ += lengthR;
      int effectiveU = lengthU - effectiveStartNumberU_;
      totalInIncreasingU_ += effectiveU;
      //lastR=lengthR;
      //lastU=lengthU;
      double rest = lengthL + 0.05 * nnd;
      double constWeightFactor = weightFactL * lengthL + weightFactDense * nnd
        + weightNrows * numberRows;
      double constWeightIterateX = constWeightIterate * (lengthL + endLengthU_)
        + weightNrowsIterate * numberRows;
      double variableWeight = weightIncU * totalInIncreasingU_ + weightR * totalInR_ + weightRest * rest;
      double average = constWeightIterateX + (constWeightFactor + variableWeight) / static_cast< double >(numberPivots);
#if 0
      if ((numberPivots%20)==0&&!ifPrint3)
      printf("PIV %d nrow %d startU %d now %d L %d R %d dense %g average %g\n",
      numberPivots,numberRows,effectiveStartNumberU_,
      lengthU,lengthL,lengthR,nnd,average);
#endif
      shortestAverage_ = std::min(shortestAverage_, average);
      if (average > increaseNeeded * shortestAverage_ && coinFactorizationA_->pivots() > 30) {
        //printf("PIVX %d nrow %d startU %d now %d L %d R %d dense %g average %g\n",
        //numberPivots,numberRows,effectiveStartNumberU_,
        //lengthU,lengthL,lengthR,nnd,average);
        reFactor3 = true;
      }
    }
    if (reFactor || reFactor3) {
#if 0
	printf("%c%cftranCountInput_ %g  ,ftranCountAfterL_ %g  ,ftranCountAfterR_ %g  ,ftranCountAfterU_ %g  ,btranCountInput_ %g  ,btranCountAfterU_ %g  ,btranCountAfterR_ %g  ,btranCountAfterL_ %g  ,numberFtranCounts_ %d  ,numberBtranCounts_ %d  ,ftranAverageAfterL_ %g  ,ftranAverageAfterR_ %g  ,ftranAverageAfterU_ %g  ,btranAverageAfterU_ %g  ,btranAverageAfterR_ %g  ,btranAverageAfterL_ %g\n"
	       ,reFactor3 ? 'Y' : 'N'
	       ,reFactor ? 'Y' : 'N'
  ,coinFactorizationA_->ftranCountInput_
  ,coinFactorizationA_->ftranCountAfterL_
  ,coinFactorizationA_->ftranCountAfterR_
  ,coinFactorizationA_->ftranCountAfterU_
  ,coinFactorizationA_->btranCountInput_
  ,coinFactorizationA_->btranCountAfterU_
  ,coinFactorizationA_->btranCountAfterR_
  ,coinFactorizationA_->btranCountAfterL_
  ,coinFactorizationA_->numberFtranCounts_
  ,coinFactorizationA_->numberBtranCounts_
  ,coinFactorizationA_->ftranAverageAfterL_
  ,coinFactorizationA_->ftranAverageAfterR_
  ,coinFactorizationA_->ftranAverageAfterU_
  ,coinFactorizationA_->btranAverageAfterU_
  ,coinFactorizationA_->btranAverageAfterR_
	       ,coinFactorizationA_->btranAverageAfterL_);
#endif
      reFactor = true;
    }
    return reFactor;
  } else {
    return coinFactorizationB_->pivots() > coinFactorizationB_->numberRows() / 2.45 + 20;
  }
}
#if CLP_FACTORIZATION_NEW_TIMING > 1
void ClpFactorization::statsRefactor(char when) const
{
  int numberPivots = coinFactorizationA_->pivots();
  int numberDense = coinFactorizationA_->numberDense();
  double nnd = numberDense * numberDense;
  int lengthL = coinFactorizationA_->numberElementsL();
  int lengthR = coinFactorizationA_->numberElementsR();
  int numberRows = coinFactorizationA_->numberRows();
  int lengthU = coinFactorizationA_->numberElementsU() - (numberRows - numberDense);
  double rest = lengthL + 0.05 * nnd;
  double constWeightFactor = weightFactL * lengthL + weightFactDense * nnd
    + weightNrows * numberRows;
  double constWeightIterateX = constWeightIterate * (lengthL + endLengthU_)
    + weightNrowsIterate * numberRows;
  double variableWeight = weightIncU * totalInIncreasingU_ + weightR * totalInR_ + weightRest * rest;
  double average = constWeightIterateX + (constWeightFactor + variableWeight) / static_cast< double >(numberPivots);
  printf("PIV%c %d nrow %d startU %d now %d L %d R %d dense %g average %g - shortest %g\n",
    when, numberPivots, numberRows, effectiveStartNumberU_,
    lengthU, lengthL, lengthR, nnd, average, shortestAverage_);
}
#endif
#else
bool ClpFactorization::timeToRefactorize() const
{
  if (coinFactorizationA_) {
    return (coinFactorizationA_->pivots() * 3 > coinFactorizationA_->maximumPivots() * 2 && coinFactorizationA_->numberElementsR() * 3 > (coinFactorizationA_->numberElementsL() + coinFactorizationA_->numberElementsU()) * 2 + 1000 && !coinFactorizationA_->numberDense());
  } else {
    return coinFactorizationB_->pivots() > coinFactorizationB_->numberRows() / 2.45 + 20;
  }
#endif
int ClpFactorization::factorize(ClpSimplex *model,
  int solveType, bool valuesPass)
{
#ifdef CLP_REUSE_ETAS
  model_ = model;
#endif
  //if ((model->specialOptions()&16384))
  //printf("factor at %d iterations\n",model->numberIterations());
  ClpMatrixBase *matrix = model->clpMatrix();
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  if (!numberRows)
    return 0;
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
  if (!timeInUpdate) {
    printf("ZZZZ,timeInFactor,nRows,nrows_2,\
startU,endU,lengthL,dense,dense_2,nslacks\n");
    printf("YYYY,time,added,dense,dense_2,averageR,averageL,averageU,\
average2R,average2L,average2U,\
scaledDense,scaledDense_2,scaledL,scaledR,scaledU\n");
    printf("WWWW,time,size0,ratio1,ratio05,ratio02\n");
  }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
  externalTimeStart = CoinCpuTime();
  memset(average, 0, sizeof(average));
  shortest = COIN_DBL_MAX;
  ifPrint = 0;
  if (!model->numberIterations())
    printf("maxfactor %d\n", coinFactorizationA_->maximumPivots());
#endif
  bool anyChanged = false;
  if (coinFactorizationB_) {
    coinFactorizationB_->setStatus(-99);
    int *COIN_RESTRICT pivotVariable = model->pivotVariable();
    //returns 0 -okay, -1 singular, -2 too many in basis */
    // allow dense
    int solveMode = 2;
    if (model->numberIterations())
      solveMode += 8;
    if (valuesPass)
      solveMode += 4;
    coinFactorizationB_->setSolveMode(solveMode);
    while (status() < -98) {

      int i;
      int numberBasic = 0;
      int numberRowBasic;
      // Move pivot variables across if they look good
      int *COIN_RESTRICT pivotTemp = model->rowArray(0)->getIndices();
      assert(!model->rowArray(0)->getNumElements());
      if (!matrix->rhsOffset(model)) {
        // Seems to prefer things in order so quickest
        // way is to go though like this
        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic)
            pivotTemp[numberBasic++] = i;
        }
        numberRowBasic = numberBasic;
        /* Put column basic variables into pivotVariable
                       This is done by ClpMatrixBase to allow override for gub
                    */
        matrix->generalExpanded(model, 0, numberBasic);
      } else {
        // Long matrix - do a different way
        bool fullSearch = false;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot >= numberColumns) {
            pivotTemp[numberBasic++] = iPivot - numberColumns;
          }
        }
        numberRowBasic = numberBasic;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot < numberColumns) {
            if (iPivot >= 0) {
              pivotTemp[numberBasic++] = iPivot;
            } else {
              // not full basis
              fullSearch = true;
              break;
            }
          }
        }
        if (fullSearch) {
          // do slow way
          numberBasic = 0;
          for (i = 0; i < numberRows; i++) {
            if (model->getRowStatus(i) == ClpSimplex::basic)
              pivotTemp[numberBasic++] = i;
          }
          numberRowBasic = numberBasic;
          /* Put column basic variables into pivotVariable
                            This is done by ClpMatrixBase to allow override for gub
                         */
          matrix->generalExpanded(model, 0, numberBasic);
        }
      }
      if (numberBasic > model->maximumBasic()) {
        // Take out some
        numberBasic = numberRowBasic;
        for (int i = 0; i < numberColumns; i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            if (numberBasic < numberRows)
              numberBasic++;
            else
              model->setColumnStatus(i, ClpSimplex::superBasic);
          }
        }
        numberBasic = numberRowBasic;
        matrix->generalExpanded(model, 0, numberBasic);
      } else if (numberBasic < numberRows) {
        // add in some
        int needed = numberRows - numberBasic;
        // move up columns
        for (i = numberBasic - 1; i >= numberRowBasic; i--)
          pivotTemp[i + needed] = pivotTemp[i];
        numberRowBasic = 0;
        numberBasic = numberRows;
        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic) {
            pivotTemp[numberRowBasic++] = i;
          } else if (needed) {
            needed--;
            model->setRowStatus(i, ClpSimplex::basic);
            pivotTemp[numberRowBasic++] = i;
          }
        }
      }
      int numberElements = numberRowBasic;

      // compute how much in basis
      // can change for gub
      int numberColumnBasic = numberBasic - numberRowBasic;

      numberElements += matrix->countBasis(pivotTemp + numberRowBasic,
        numberColumnBasic);
#ifndef NDEBUG
//#define CHECK_CLEAN_BASIS
#ifdef CHECK_CLEAN_BASIS
      int saveNumberElements = numberElements;
#elif COIN_BIG_DOUBLE
      int saveNumberElements = numberElements;
#endif
#endif
      // Not needed for dense
      numberElements = 3 * numberBasic + 3 * numberElements + 20000;
      int numberIterations = model->numberIterations();
      coinFactorizationB_->setUsefulInformation(&numberIterations, 0);
      coinFactorizationB_->getAreas(numberRows,
        numberRowBasic + numberColumnBasic, numberElements,
        2 * numberElements);
      // Fill in counts so we can skip part of preProcess
      // This is NOT needed for dense but would be needed for later versions
#if COIN_BIG_DOUBLE == 0
      double *COIN_RESTRICT elementU = coinFactorizationB_->elements();
#else
      CoinFactorizationDouble *COIN_RESTRICT elementU =
	new CoinFactorizationDouble [saveNumberElements];
#endif
      int *COIN_RESTRICT indexRowU;
      int *COIN_RESTRICT startColumnU;
      int *COIN_RESTRICT numberInRow;
      int *COIN_RESTRICT numberInColumn;
      indexRowU = coinFactorizationB_->indices();
      startColumnU = coinFactorizationB_->starts();
#ifdef CHECK_CLEAN_BASIS
      for (int i = 0; i < saveNumberElements; i++) {
        elementU[i] = 0.0;
        indexRowU[i] = -1;
      }
      for (int i = 0; i < numberRows; i++)
        startColumnU[i] = -1;
#endif
#ifndef COIN_FAST_CODE
      double slackValue;
      slackValue = coinFactorizationB_->slackValue();
#else
#define slackValue -1.0
#endif
      numberInRow = coinFactorizationB_->numberInRow();
      numberInColumn = coinFactorizationB_->numberInColumn();
      CoinZeroN(numberInRow, numberRows);
      CoinZeroN(numberInColumn, numberRows);
      for (i = 0; i < numberRowBasic; i++) {
        int iRow = pivotTemp[i];
        // Change pivotTemp to correct sequence
        pivotTemp[i] = iRow + numberColumns;
        indexRowU[i] = iRow;
        startColumnU[i] = i;
        elementU[i] = slackValue;
        numberInRow[iRow] = 1;
        numberInColumn[i] = 1;
      }
      startColumnU[numberRowBasic] = numberRowBasic;
      // can change for gub so redo
      numberColumnBasic = numberBasic - numberRowBasic;
      matrix->fillBasis(model,
        pivotTemp + numberRowBasic,
        numberColumnBasic,
        indexRowU,
        startColumnU + numberRowBasic,
        numberInRow,
        numberInColumn + numberRowBasic,
        elementU);
      // recompute number basic
      numberBasic = numberRowBasic + numberColumnBasic;
      for (i = numberBasic; i < numberRows; i++)
        pivotTemp[i] = -1; // mark not there
      if (numberBasic)
        numberElements = startColumnU[numberBasic - 1]
          + numberInColumn[numberBasic - 1];
      else
        numberElements = 0;
#ifdef CHECK_CLEAN_BASIS
      assert(!startColumnU[0]);
      int lastStart = 0;
      for (int i = 0; i < numberRows; i++) {
        assert(startColumnU[i + 1] > lastStart);
        lastStart = startColumnU[i + 1];
      }
      assert(lastStart == saveNumberElements);
      for (int i = 0; i < saveNumberElements; i++) {
        assert(elementU[i]);
        assert(indexRowU[i] >= 0 && indexRowU[i] < numberRows);
      }
#endif
#if COIN_BIG_DOUBLE
      double *COIN_RESTRICT  elsU = coinFactorizationB_->elements();
      for (int i=0;i<saveNumberElements;i++)
	elsU[i] = elementU[i];;
#endif
      coinFactorizationB_->preProcess();
      coinFactorizationB_->factor();
      if (coinFactorizationB_->status() == -1 && (coinFactorizationB_->solveMode() % 3) != 0) {
        int solveMode = coinFactorizationB_->solveMode();
        solveMode -= solveMode % 3; // so bottom will be 0
        coinFactorizationB_->setSolveMode(solveMode);
        coinFactorizationB_->setStatus(-99);
      }
      if (coinFactorizationB_->status() == -99)
        continue;
      // If we get here status is 0 or -1
      if (coinFactorizationB_->status() == 0 && numberBasic == numberRows) {
        coinFactorizationB_->postProcess(pivotTemp, pivotVariable);
      } else if (solveType == 0 || solveType == 2 /*||solveType==1*/) {
        // Change pivotTemp to be correct list
        anyChanged = true;
        coinFactorizationB_->makeNonSingular(pivotTemp, numberColumns);
        double *COIN_RESTRICT columnLower = model->lowerRegion();
        double *COIN_RESTRICT columnUpper = model->upperRegion();
        double *COIN_RESTRICT columnActivity = model->solutionRegion();
        double *COIN_RESTRICT rowLower = model->lowerRegion(0);
        double *COIN_RESTRICT rowUpper = model->upperRegion(0);
        double *COIN_RESTRICT rowActivity = model->solutionRegion(0);
        //redo basis - first take ALL out
        int iColumn;
        double largeValue = model->largeValue();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) == ClpSimplex::basic) {
            // take out
            if (!valuesPass) {
              double lower = columnLower[iColumn];
              double upper = columnUpper[iColumn];
              double value = columnActivity[iColumn];
              if (lower > -largeValue || upper < largeValue) {
                if (fabs(value - lower) < fabs(value - upper)) {
                  model->setColumnStatus(iColumn, ClpSimplex::atLowerBound);
                  columnActivity[iColumn] = lower;
                } else {
                  model->setColumnStatus(iColumn, ClpSimplex::atUpperBound);
                  columnActivity[iColumn] = upper;
                }
              } else {
                model->setColumnStatus(iColumn, ClpSimplex::isFree);
              }
            } else {
              model->setColumnStatus(iColumn, ClpSimplex::superBasic);
            }
          }
        }
        int iRow;
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (model->getRowStatus(iRow) == ClpSimplex::basic) {
            // take out
            if (!valuesPass) {
              double lower = columnLower[iRow];
              double upper = columnUpper[iRow];
              double value = columnActivity[iRow];
              if (lower > -largeValue || upper < largeValue) {
                if (fabs(value - lower) < fabs(value - upper)) {
                  model->setRowStatus(iRow, ClpSimplex::atLowerBound);
                  columnActivity[iRow] = lower;
                } else {
                  model->setRowStatus(iRow, ClpSimplex::atUpperBound);
                  columnActivity[iRow] = upper;
                }
              } else {
                model->setRowStatus(iRow, ClpSimplex::isFree);
              }
            } else {
              model->setRowStatus(iRow, ClpSimplex::superBasic);
            }
          }
        }
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iSequence = pivotTemp[iRow];
          assert(iSequence >= 0);
          // basic
          model->setColumnStatus(iSequence, ClpSimplex::basic);
        }
        // signal repeat
        coinFactorizationB_->setStatus(-99);
        // set fixed if they are
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (model->getRowStatus(iRow) != ClpSimplex::basic) {
            if (rowLower[iRow] == rowUpper[iRow]) {
              rowActivity[iRow] = rowLower[iRow];
              model->setRowStatus(iRow, ClpSimplex::isFixed);
            }
          }
        }
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) != ClpSimplex::basic) {
            if (columnLower[iColumn] == columnUpper[iColumn]) {
              columnActivity[iColumn] = columnLower[iColumn];
              model->setColumnStatus(iColumn, ClpSimplex::isFixed);
            }
          }
        }
      }
    }
#ifdef CLP_DEBUG
    // check basic
    CoinIndexedVector region1(2 * numberRows);
    CoinIndexedVector region2B(2 * numberRows);
    int iPivot;
    double *COIN_RESTRICT arrayB = region2B.denseVector();
    int i;
    for (iPivot = 0; iPivot < numberRows; iPivot++) {
      int iSequence = pivotVariable[iPivot];
      model->unpack(&region2B, iSequence);
      coinFactorizationB_->updateColumn(&region1, &region2B);
      if (fabs(arrayB[iPivot] - 1.0) < 1.0e-4) {
        // OK?
        arrayB[iPivot] = 0.0;
      } else {
        assert(fabs(arrayB[iPivot]) < 1.0e-4);
        for (i = 0; i < numberRows; i++) {
          if (fabs(arrayB[i] - 1.0) < 1.0e-4)
            break;
        }
        assert(i < numberRows);
        printf("variable on row %d landed up on row %d\n", iPivot, i);
        arrayB[i] = 0.0;
      }
      for (i = 0; i < numberRows; i++)
        assert(fabs(arrayB[i]) < 1.0e-4);
      region2B.clear();
    }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(2);
#endif
    if (anyChanged && model->algorithm() < 0 && solveType > 0) {
      double dummyCost;
      static_cast< ClpSimplexDual * >(model)->changeBounds(3,
        NULL, dummyCost);
    }
    return coinFactorizationB_->status();
  }
  // If too many compressions increase area
  if (coinFactorizationA_->pivots() > 1 && coinFactorizationA_->numberCompressions() * 10 > coinFactorizationA_->pivots() + 10) {
    coinFactorizationA_->areaFactor(coinFactorizationA_->areaFactor() * 1.1);
  }
  //int numberPivots=coinFactorizationA_->pivots();
#if 0
     if (model->algorithm() > 0)
          numberSave = -1;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_ || doCheck) {
#endif
    coinFactorizationA_->setStatus(-99);
    int *COIN_RESTRICT pivotVariable = model->pivotVariable();
    int nTimesRound = 0;
    //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
    while (coinFactorizationA_->status() < -98) {
      nTimesRound++;

      int i;
      int numberBasic = 0;
      int numberRowBasic;
      // Move pivot variables across if they look good
      int *COIN_RESTRICT pivotTemp = model->rowArray(0)->getIndices();
      assert(!model->rowArray(0)->getNumElements());
      if (!matrix->rhsOffset(model)) {
#if 0
                    if (numberSave > 0) {
                         int nStill = 0;
                         int nAtBound = 0;
                         int nZeroDual = 0;
                         CoinIndexedVector * array = model->rowArray(3);
                         CoinIndexedVector * objArray = model->columnArray(1);
                         array->clear();
                         objArray->clear();
                         double *COIN_RESTRICT  cost = model->costRegion();
                         double tolerance = model->primalTolerance();
                         double offset = 0.0;
                         for (i = 0; i < numberRows; i++) {
                              int iPivot = pivotVariable[i];
                              if (iPivot < numberColumns && isDense(iPivot)) {
                                   if (model->getColumnStatus(iPivot) == ClpSimplex::basic) {
                                        nStill++;
                                        double value = model->solutionRegion()[iPivot];
                                        double dual = model->dualRowSolution()[i];
                                        double lower = model->lowerRegion()[iPivot];
                                        double upper = model->upperRegion()[iPivot];
                                        ClpSimplex::Status status;
                                        if (fabs(value - lower) < tolerance) {
                                             status = ClpSimplex::atLowerBound;
                                             nAtBound++;
                                        } else if (fabs(value - upper) < tolerance) {
                                             nAtBound++;
                                             status = ClpSimplex::atUpperBound;
                                        } else if (value > lower && value < upper) {
                                             status = ClpSimplex::superBasic;
                                        } else {
                                             status = ClpSimplex::basic;
                                        }
                                        if (status != ClpSimplex::basic) {
                                             if (model->getRowStatus(i) != ClpSimplex::basic) {
                                                  model->setColumnStatus(iPivot, ClpSimplex::atLowerBound);
                                                  model->setRowStatus(i, ClpSimplex::basic);
                                                  pivotVariable[i] = i + numberColumns;
                                                  model->dualRowSolution()[i] = 0.0;
                                                  model->djRegion(0)[i] = 0.0;
                                                  array->add(i, dual);
                                                  offset += dual * model->solutionRegion(0)[i];
                                             }
                                        }
                                        if (fabs(dual) < 1.0e-5)
                                             nZeroDual++;
                                   }
                              }
                         }
                         printf("out of %d dense, %d still in basis, %d at bound, %d with zero dual - offset %g\n",
                                numberSave, nStill, nAtBound, nZeroDual, offset);
                         if (array->getNumElements()) {
                              // modify costs
                              model->clpMatrix()->transposeTimes(model, 1.0, array, model->columnArray(0),
                                                                 objArray);
                              array->clear();
                              int n = objArray->getNumElements();
                              int *COIN_RESTRICT  indices = objArray->getIndices();
                              double *COIN_RESTRICT  elements = objArray->denseVector();
                              for (i = 0; i < n; i++) {
                                   int iColumn = indices[i];
                                   cost[iColumn] -= elements[iColumn];
                                   elements[iColumn] = 0.0;
                              }
                              objArray->setNumElements(0);
                         }
                    }
#endif
        // Seems to prefer things in order so quickest
        // way is to go though like this
        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic)
            pivotTemp[numberBasic++] = i;
        }
        numberRowBasic = numberBasic;
        /* Put column basic variables into pivotVariable
                       This is done by ClpMatrixBase to allow override for gub
                    */
        matrix->generalExpanded(model, 0, numberBasic);
      } else {
        // Long matrix - do a different way
        bool fullSearch = false;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot >= numberColumns) {
            pivotTemp[numberBasic++] = iPivot - numberColumns;
          }
        }
        numberRowBasic = numberBasic;
        for (i = 0; i < numberRows; i++) {
          int iPivot = pivotVariable[i];
          if (iPivot < numberColumns) {
            if (iPivot >= 0) {
              pivotTemp[numberBasic++] = iPivot;
            } else {
              // not full basis
              fullSearch = true;
              break;
            }
          }
        }
        if (fullSearch) {
          // do slow way
          numberBasic = 0;
          for (i = 0; i < numberRows; i++) {
            if (model->getRowStatus(i) == ClpSimplex::basic)
              pivotTemp[numberBasic++] = i;
          }
          numberRowBasic = numberBasic;
          /* Put column basic variables into pivotVariable
                            This is done by ClpMatrixBase to allow override for gub
                         */
          matrix->generalExpanded(model, 0, numberBasic);
        }
      }
      if (numberBasic > model->maximumBasic()) {
#if 0 // ndef NDEBUG
                    printf("%d basic - should only be %d\n",
                           numberBasic, numberRows);
#endif
        // Take out some
        numberBasic = numberRowBasic;
        for (int i = 0; i < numberColumns; i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            if (numberBasic < numberRows)
              numberBasic++;
            else
              model->setColumnStatus(i, ClpSimplex::superBasic);
          }
        }
        numberBasic = numberRowBasic;
        matrix->generalExpanded(model, 0, numberBasic);
      }
#ifndef SLIM_CLP
      // see if matrix a network
#ifndef NO_RTTI
      ClpNetworkMatrix *networkMatrix = dynamic_cast< ClpNetworkMatrix * >(model->clpMatrix());
#else
      ClpNetworkMatrix *networkMatrix = NULL;
      if (model->clpMatrix()->type() == 11)
        networkMatrix = static_cast< ClpNetworkMatrix * >(model->clpMatrix());
#endif
      // If network - still allow ordinary factorization first time for laziness
      if (networkMatrix)
        coinFactorizationA_->setBiasLU(0); // All to U if network
      //int saveMaximumPivots = maximumPivots();
      delete networkBasis_;
      networkBasis_ = NULL;
      if (networkMatrix && !doCheck)
        maximumPivots(1);
#endif
      //printf("L, U, R %d %d %d\n",numberElementsL(),numberElementsU(),numberElementsR());
      while (coinFactorizationA_->status() == -99) {
        // maybe for speed will be better to leave as many regions as possible
        coinFactorizationA_->gutsOfDestructor();
        coinFactorizationA_->gutsOfInitialize(2);
        int numberElements = numberRowBasic;

        // compute how much in basis

        int i;
        // can change for gub
        int numberColumnBasic = numberBasic - numberRowBasic;

        numberElements += matrix->countBasis(pivotTemp + numberRowBasic,
          numberColumnBasic);
        // and recompute as network side say different
        if (model->numberIterations())
          numberRowBasic = numberBasic - numberColumnBasic;
        numberElements = 3 * numberBasic + 3 * numberElements + 20000;
        coinFactorizationA_->getAreas(numberRows,
          numberRowBasic + numberColumnBasic, numberElements,
          2 * numberElements);
        //fill
        // Fill in counts so we can skip part of preProcess
        int *COIN_RESTRICT numberInRow = coinFactorizationA_->numberInRow();
        int *COIN_RESTRICT numberInColumn = coinFactorizationA_->numberInColumn();
        CoinZeroN(numberInRow, coinFactorizationA_->numberRows() + 1);
        CoinZeroN(numberInColumn, coinFactorizationA_->maximumColumnsExtra() + 1);
        CoinFactorizationDouble *COIN_RESTRICT elementU = coinFactorizationA_->elementU();
        int *COIN_RESTRICT indexRowU = coinFactorizationA_->indexRowU();
        int *COIN_RESTRICT startColumnU = coinFactorizationA_->startColumnU();
#ifndef COIN_FAST_CODE
        double slackValue = coinFactorizationA_->slackValue();
#endif
        for (i = 0; i < numberRowBasic; i++) {
          int iRow = pivotTemp[i];
          indexRowU[i] = iRow;
          startColumnU[i] = i;
          elementU[i] = slackValue;
          numberInRow[iRow] = 1;
          numberInColumn[i] = 1;
        }
        startColumnU[numberRowBasic] = numberRowBasic;
        // can change for gub so redo
        numberColumnBasic = numberBasic - numberRowBasic;
        matrix->fillBasis(model,
          pivotTemp + numberRowBasic,
          numberColumnBasic,
          indexRowU,
          startColumnU + numberRowBasic,
          numberInRow,
          numberInColumn + numberRowBasic,
          elementU);
#if 0
                    {
                         printf("%d row basic, %d column basic\n", numberRowBasic, numberColumnBasic);
                         for (int i = 0; i < numberElements; i++)
                              printf("row %d col %d value %g\n", indexRowU[i], indexColumnU_[i],
                                     elementU[i]);
                    }
#endif
        // recompute number basic
        numberBasic = numberRowBasic + numberColumnBasic;
        if (numberBasic)
          numberElements = startColumnU[numberBasic - 1]
            + numberInColumn[numberBasic - 1];
        else
          numberElements = 0;
#ifdef CLP_FACTORIZATION_NEW_TIMING
        lastNumberPivots_ = 0;
        effectiveStartNumberU_ = numberElements - numberRows;
        //printf("%d slacks,%d in U at beginning\n",
        //numberRowBasic,numberElements);
#endif
        coinFactorizationA_->setNumberElementsU(numberElements);
        //saveFactorization("dump.d");
        if (coinFactorizationA_->biasLU() >= 3 || coinFactorizationA_->numberRows() != coinFactorizationA_->numberColumns())
          coinFactorizationA_->preProcess(2);
        else
          coinFactorizationA_->preProcess(3); // no row copy
        coinFactorizationA_->factor();
#ifdef CLP_FACTORIZATION_NEW_TIMING
        endLengthU_ = coinFactorizationA_->numberElements() - coinFactorizationA_->numberDense() * coinFactorizationA_->numberDense()
          - coinFactorizationA_->numberElementsL();
#endif
        if (coinFactorizationA_->status() == -99) {
          // get more memory
          coinFactorizationA_->areaFactor(2.0 * coinFactorizationA_->areaFactor());
        } else if (coinFactorizationA_->status() == -1 && (model->numberIterations() == 0 || nTimesRound > 2) && coinFactorizationA_->denseThreshold()) {
          // Round again without dense
          coinFactorizationA_->setDenseThreshold(0);
          coinFactorizationA_->setStatus(-99);
        }
      }
      // If we get here status is 0 or -1
      if (coinFactorizationA_->status() == 0) {
        // We may need to tamper with order and redo - e.g. network with side
        int useNumberRows = numberRows;
        // **** we will also need to add test in dual steepest to do
        // as we do for network
        matrix->generalExpanded(model, 12, useNumberRows);
        const int *permuteBack = coinFactorizationA_->permuteBack();
        const int *back = coinFactorizationA_->pivotColumnBack();
        //int * pivotTemp = pivotColumn_.array();
        //ClpDisjointCopyN ( pivotVariable, numberRows , pivotTemp  );
#ifndef NDEBUG
        CoinFillN(pivotVariable, numberRows, -1);
#endif
        // Redo pivot order
        for (i = 0; i < numberRowBasic; i++) {
          int k = pivotTemp[i];
          // so rowIsBasic[k] would be permuteBack[back[i]]
          int j = permuteBack[back[i]];
          assert(pivotVariable[j] == -1);
          pivotVariable[j] = k + numberColumns;
        }
        for (; i < useNumberRows; i++) {
          int k = pivotTemp[i];
          // so rowIsBasic[k] would be permuteBack[back[i]]
          int j = permuteBack[back[i]];
          assert(pivotVariable[j] == -1);
          pivotVariable[j] = k;
        }
#if 0
                    if (numberSave >= 0) {
                         numberSave = numberDense_;
                         memset(saveList, 0, ((coinFactorizationA_->numberRows() + 31) >> 5)*sizeof(int));
                         for (i = coinFactorizationA_->numberRows() - numberSave; i < coinFactorizationA_->numberRows(); i++) {
                              int k = pivotTemp[pivotColumn_.array()[i]];
                              setDense(k);
                         }
                    }
#endif
        // Set up permutation vector
        // these arrays start off as copies of permute
        // (and we could use permute_ instead of pivotColumn (not back though))
        ClpDisjointCopyN(coinFactorizationA_->permute(), useNumberRows, coinFactorizationA_->pivotColumn());
        ClpDisjointCopyN(coinFactorizationA_->permuteBack(), useNumberRows, coinFactorizationA_->pivotColumnBack());
#ifndef SLIM_CLP
        if (networkMatrix) {
          maximumPivots(std::max(2000, maximumPivots()));
          // redo arrays
          for (int iRow = 0; iRow < 4; iRow++) {
            int length = model->numberRows() + maximumPivots();
            if (iRow == 3 || model->objectiveAsObject()->type() > 1)
              length += model->numberColumns();
            model->rowArray(iRow)->reserve(length);
          }
          // create network factorization
          if (doCheck)
            delete networkBasis_; // temp
          networkBasis_ = new ClpNetworkBasis(model, coinFactorizationA_->numberRows(),
            coinFactorizationA_->pivotRegion(),
            coinFactorizationA_->permuteBack(),
            coinFactorizationA_->startColumnU(),
            coinFactorizationA_->numberInColumn(),
            coinFactorizationA_->indexRowU(),
            coinFactorizationA_->elementU());
          // kill off arrays in ordinary factorization
          if (!doCheck) {
            coinFactorizationA_->gutsOfDestructor();
            // but make sure coinFactorizationA_->numberRows() set
            coinFactorizationA_->setNumberRows(model->numberRows());
            coinFactorizationA_->setStatus(0);
#if 0
                              // but put back permute arrays so odd things will work
                              int numberRows = model->numberRows();
                              pivotColumnBack_ = new int [numberRows];
                              permute_ = new int [numberRows];
                              int i;
                              for (i = 0; i < numberRows; i++) {
                                   pivotColumnBack_[i] = i;
                                   permute_[i] = i;
                              }
#endif
          }
        } else {
#endif
          // See if worth going sparse and when
          coinFactorizationA_->checkSparse();
#ifndef SLIM_CLP
        }
#endif
      } else if (coinFactorizationA_->status() == -1 && (solveType == 0 || solveType == 2)) {
        // This needs redoing as it was merged coding - does not need array
#if 1
        int numberTotal = numberRows + numberColumns;
        int *COIN_RESTRICT isBasic = new int[numberTotal];
        int *COIN_RESTRICT rowIsBasic = isBasic + numberColumns;
        int *COIN_RESTRICT columnIsBasic = isBasic;
        for (i = 0; i < numberTotal; i++)
          isBasic[i] = -1;
        for (i = 0; i < numberRowBasic; i++) {
          int iRow = pivotTemp[i];
          rowIsBasic[iRow] = 1;
        }
        for (; i < numberBasic; i++) {
          int iColumn = pivotTemp[i];
          columnIsBasic[iColumn] = 1;
        }
        numberBasic = 0;
        for (i = 0; i < numberRows; i++)
          pivotVariable[i] = -1;
        // mark as basic or non basic
        const int *pivotColumn = coinFactorizationA_->pivotColumn();
        for (i = 0; i < numberRows; i++) {
          if (rowIsBasic[i] >= 0) {
            if (pivotColumn[numberBasic] >= 0) {
              rowIsBasic[i] = pivotColumn[numberBasic];
            } else {
              rowIsBasic[i] = -1;
              model->setRowStatus(i, ClpSimplex::superBasic);
            }
            numberBasic++;
          }
        }
        for (i = 0; i < numberColumns; i++) {
          if (columnIsBasic[i] >= 0) {
            if (pivotColumn[numberBasic] >= 0)
              columnIsBasic[i] = pivotColumn[numberBasic];
            else
              columnIsBasic[i] = -1;
            numberBasic++;
          }
        }
        // leave pivotVariable in useful form for cleaning basis
        int *COIN_RESTRICT pivotVariable = model->pivotVariable();
        for (i = 0; i < numberRows; i++) {
          pivotVariable[i] = -1;
        }

        for (i = 0; i < numberRows; i++) {
          if (model->getRowStatus(i) == ClpSimplex::basic) {
            int iPivot = rowIsBasic[i];
            if (iPivot >= 0)
              pivotVariable[iPivot] = i + numberColumns;
          }
        }
        for (i = 0; i < numberColumns; i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            int iPivot = columnIsBasic[i];
            if (iPivot >= 0)
              pivotVariable[iPivot] = i;
          }
        }
        delete[] isBasic;
#else
        {
          //int * lastColumn = lastColumn_.array(); // -1 or pivot row
          int *COIN_RESTRICT lastRow = coinFactorizationA_->lastRow(); // -1 or pivot sequence (inside sequence)
          for (int i = 0; i < numberRows; i++) {
            int iSeq = lastRow[i];
            if (iSeq >= 0) {
              pivotVariable[i] = pivotTemp[iSeq];
              model->setRowStatus(i, ClpSimplex::superBasic);
            } else {
              pivotVariable[i] = i + numberColumns;
            }
          }
        }
#endif
        double *COIN_RESTRICT columnLower = model->lowerRegion();
        double *COIN_RESTRICT columnUpper = model->upperRegion();
        double *COIN_RESTRICT columnActivity = model->solutionRegion();
        double *COIN_RESTRICT rowLower = model->lowerRegion(0);
        double *COIN_RESTRICT rowUpper = model->upperRegion(0);
        double *COIN_RESTRICT rowActivity = model->solutionRegion(0);
        //redo basis - first take ALL columns out
        int iColumn;
        double largeValue = model->largeValue();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) == ClpSimplex::basic) {
            // take out
            if (!valuesPass) {
              double lower = columnLower[iColumn];
              double upper = columnUpper[iColumn];
              double value = columnActivity[iColumn];
              if (lower > -largeValue || upper < largeValue) {
                if (fabs(value - lower) < fabs(value - upper)) {
                  model->setColumnStatus(iColumn, ClpSimplex::atLowerBound);
                  columnActivity[iColumn] = lower;
                } else {
                  model->setColumnStatus(iColumn, ClpSimplex::atUpperBound);
                  columnActivity[iColumn] = upper;
                }
              } else {
                model->setColumnStatus(iColumn, ClpSimplex::isFree);
              }
            } else {
              model->setColumnStatus(iColumn, ClpSimplex::superBasic);
            }
          }
        }
        int iRow;
        for (iRow = 0; iRow < numberRows; iRow++) {
          int iSequence = pivotVariable[iRow];
          if (iSequence >= 0) {
            // basic
            if (iSequence >= numberColumns) {
              // slack in - leave
              //assert (iSequence-numberColumns==iRow);
            } else {
              assert(model->getRowStatus(iRow) != ClpSimplex::basic);
              // put back structural
              model->setColumnStatus(iSequence, ClpSimplex::basic);
            }
          } else {
            // put in slack
            model->setRowStatus(iRow, ClpSimplex::basic);
          }
        }
        // Put back any key variables for gub
        int dummy;
        matrix->generalExpanded(model, 1, dummy);
        // signal repeat
        coinFactorizationA_->setStatus(-99);
        // set fixed if they are
        for (iRow = 0; iRow < numberRows; iRow++) {
          if (model->getRowStatus(iRow) != ClpSimplex::basic) {
            if (rowLower[iRow] == rowUpper[iRow]) {
              rowActivity[iRow] = rowLower[iRow];
              model->setRowStatus(iRow, ClpSimplex::isFixed);
            }
          }
        }
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (model->getColumnStatus(iColumn) != ClpSimplex::basic) {
            if (columnLower[iColumn] == columnUpper[iColumn]) {
              columnActivity[iColumn] = columnLower[iColumn];
              model->setColumnStatus(iColumn, ClpSimplex::isFixed);
            }
          }
        }
      }
    }
#ifndef SLIM_CLP
  } else {
    // network - fake factorization - do nothing
    coinFactorizationA_->setStatus(0);
    coinFactorizationA_->setPivots(0);
  }
#endif
#ifndef SLIM_CLP
  if (!coinFactorizationA_->status()) {
    // take out part if quadratic
    if (model->algorithm() == 2) {
      ClpObjective *obj = model->objectiveAsObject();
#ifndef NDEBUG
      ClpQuadraticObjective *quadraticObj = (dynamic_cast< ClpQuadraticObjective * >(obj));
      assert(quadraticObj);
#else
      ClpQuadraticObjective *quadraticObj = (static_cast< ClpQuadraticObjective * >(obj));
#endif
      CoinPackedMatrix *quadratic = quadraticObj->quadraticObjective();
      int numberXColumns = quadratic->getNumCols();
      assert(numberXColumns < numberColumns);
      int base = numberColumns - numberXColumns;
      int *COIN_RESTRICT which = new int[numberXColumns];
      int *COIN_RESTRICT pivotVariable = model->pivotVariable();
      int *COIN_RESTRICT permute = pivotColumn();
      int i;
      int n = 0;
      for (i = 0; i < numberRows; i++) {
        int iSj = pivotVariable[i] - base;
        if (iSj >= 0 && iSj < numberXColumns)
          which[n++] = permute[i];
      }
      if (n)
        coinFactorizationA_->emptyRows(n, which);
      delete[] which;
    }
  }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(2);
#endif
  return coinFactorizationA_->status();
}
/* Replaces one Column in basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   If checkBeforeModifying is true will do all accuracy checks
   before modifying factorization.  Whether to set this depends on
   speed considerations.  You could just do this on first iteration
   after factorization and thereafter re-factorize
   partial update already in U */
int ClpFactorization::replaceColumn(const ClpSimplex *model,
  CoinIndexedVector *regionSparse,
  CoinIndexedVector *tableauColumn,
  int pivotRow,
  double pivotCheck,
  bool checkBeforeModifying,
  double acceptablePivot)
{
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_NEW_TIMING
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    int nOld = 0;
    int nNew = 0;
    int seq;
#if COINUTILS_BIGINDEX_IS_INT && CLP_POOL_MATRIX == 0
    const CoinPackedMatrix *matrix = model->matrix();
    const int *columnLength = matrix->getVectorLengths();
#else
    const int *columnLength = model->clpMatrix()->getVectorLengths();
#endif
    seq = model->sequenceIn();
    if (seq >= 0 && seq < model->numberColumns() + model->numberRows()) {
      if (seq < model->numberColumns()) {
        nNew = columnLength[seq];
      } else {
        nNew = 1;
      }
    }
    seq = model->sequenceOut();
    if (seq >= 0 && seq < model->numberColumns() + model->numberRows()) {
      if (seq < model->numberColumns()) {
        nOld = columnLength[seq];
      } else {
        nOld = 1;
      }
    }
    effectiveStartNumberU_ += nNew - nOld;
#endif
    int returnCode;
    // see if FT
    if (!coinFactorizationA_ || coinFactorizationA_->forrestTomlin()) {
      if (coinFactorizationA_) {
#if ABC_USE_COIN_FACTORIZATION < 2
        returnCode = coinFactorizationA_->replaceColumn(regionSparse,
          pivotRow,
          pivotCheck,
          checkBeforeModifying,
          acceptablePivot);
#else
        // fake btran alpha until I understand
        double btranAlpha = model->alpha();
        double ftAlpha = coinFactorizationA_->checkReplacePart1(regionSparse,
          pivotRow);
        returnCode = coinFactorizationA_->checkReplacePart2(pivotRow,
          btranAlpha,
          model->alpha(),
          ftAlpha,
          acceptablePivot);
        if (returnCode < 2)
          coinFactorizationA_->replaceColumnPart3(regionSparse,
            pivotRow,
            model->alpha());
#endif
      } else {
        bool tab = coinFactorizationB_->wantsTableauColumn();
#ifdef CLP_REUSE_ETAS
        int tempInfo[2];
        tempInfo[1] = model_->sequenceOut();
#else
        int tempInfo[1];
#endif
        tempInfo[0] = model->numberIterations();
        coinFactorizationB_->setUsefulInformation(tempInfo, 1);
        returnCode = coinFactorizationB_->replaceColumn(tab ? tableauColumn : regionSparse,
          pivotRow,
          pivotCheck,
          checkBeforeModifying,
          acceptablePivot);
#ifdef CLP_DEBUG
        // check basic
        int numberRows = coinFactorizationB_->numberRows();
        CoinIndexedVector region1(2 * numberRows);
        CoinIndexedVector region2A(2 * numberRows);
        CoinIndexedVector region2B(2 * numberRows);
        int iPivot;
        double *COIN_RESTRICT arrayB = region2B.denseVector();
        int *COIN_RESTRICT pivotVariable = model->pivotVariable();
        int i;
        for (iPivot = 0; iPivot < numberRows; iPivot++) {
          int iSequence = pivotVariable[iPivot];
          if (iPivot == pivotRow)
            iSequence = model->sequenceIn();
          model->unpack(&region2B, iSequence);
          coinFactorizationB_->updateColumn(&region1, &region2B);
          assert(fabs(arrayB[iPivot] - 1.0) < 1.0e-4);
          arrayB[iPivot] = 0.0;
          for (i = 0; i < numberRows; i++)
            assert(fabs(arrayB[i]) < 1.0e-4);
          region2B.clear();
        }
#endif
      }
    } else {
      returnCode = coinFactorizationA_->replaceColumnPFI(tableauColumn,
        pivotRow, pivotCheck); // Note array
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(3);
#endif
    return returnCode;

#ifndef SLIM_CLP
  } else {
    if (doCheck) {
      int returnCode = coinFactorizationA_->replaceColumn(regionSparse,
        pivotRow,
        pivotCheck,
        checkBeforeModifying,
        acceptablePivot);
      networkBasis_->replaceColumn(regionSparse,
        pivotRow);
      return returnCode;
    } else {
      // increase number of pivots
      coinFactorizationA_->setPivots(coinFactorizationA_->pivots() + 1);
      return networkBasis_->replaceColumn(regionSparse,
        pivotRow);
    }
  }
#endif
}

/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnFT(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2)
{
#ifdef CLP_DEBUG
  regionSparse->checkClear();
#endif
  if (!numberRows())
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    int returnCode;
    if (coinFactorizationA_) {
      coinFactorizationA_->setCollectStatistics(true);
      returnCode = coinFactorizationA_->updateColumnFT(regionSparse,
        regionSparse2);
      coinFactorizationA_->setCollectStatistics(false);
    } else {
#ifdef CLP_REUSE_ETAS
      int tempInfo[2];
      tempInfo[0] = model_->numberIterations();
      tempInfo[1] = model_->sequenceIn();
      coinFactorizationB_->setUsefulInformation(tempInfo, 2);
#endif
      returnCode = coinFactorizationB_->updateColumnFT(regionSparse,
        regionSparse2);
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(4);
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *COIN_RESTRICT check = new double[coinFactorizationA_->numberRows()];
    int returnCode = coinFactorizationA_->updateColumnFT(regionSparse,
      regionSparse2);
    networkBasis_->updateColumn(regionSparse, save, -1);
    int i;
    double *COIN_RESTRICT array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, coinFactorizationA_->numberRows() * sizeof(double));
    double *COIN_RESTRICT array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      int numberRows = coinFactorizationA_->numberRows();
      for (i = 0; i < numberRows; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse, regionSparse2, -1);
    return 1;
#endif
  }
#endif
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumn(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  bool noPermute) const
{
#ifdef CLP_DEBUG
  if (!noPermute)
    regionSparse->checkClear();
#endif
  if (!numberRows())
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    int returnCode;
    if (coinFactorizationA_) {
      coinFactorizationA_->setCollectStatistics(doStatistics_);
      returnCode = coinFactorizationA_->updateColumn(regionSparse,
        regionSparse2,
        noPermute);
      coinFactorizationA_->setCollectStatistics(false);
    } else {
      returnCode = coinFactorizationB_->updateColumn(regionSparse,
        regionSparse2,
        noPermute);
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(5);
#endif
    //#define PRINT_VECTOR
#ifdef PRINT_VECTOR
    printf("Update\n");
    regionSparse2->print();
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *COIN_RESTRICT check = new double[coinFactorizationA_->numberRows()];
    int returnCode = coinFactorizationA_->updateColumn(regionSparse,
      regionSparse2,
      noPermute);
    networkBasis_->updateColumn(regionSparse, save, -1);
    int i;
    double *COIN_RESTRICT array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, coinFactorizationA_->numberRows() * sizeof(double));
    double *COIN_RESTRICT array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      int numberRows = coinFactorizationA_->numberRows();
      for (i = 0; i < numberRows; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse, regionSparse2, -1);
    return 1;
#endif
  }
#endif
}
/* Updates one column (FTRAN) from region2
   Tries to do FT update
   number returned is negative if no room.
   Also updates region3
   region1 starts as zero and is zero at end */
int ClpFactorization::updateTwoColumnsFT(CoinIndexedVector *regionSparse1,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3,
  bool noPermuteRegion3)
{
#ifdef CLP_DEBUG
  regionSparse1->checkClear();
#endif
  if (!numberRows())
    return 0;
  int returnCode = 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    if (coinFactorizationA_) {
      coinFactorizationA_->setCollectStatistics(true);
      if (coinFactorizationA_->spaceForForrestTomlin()) {
        assert(regionSparse2->packedMode());
        assert(!regionSparse3->packedMode());
        returnCode = coinFactorizationA_->updateTwoColumnsFT(regionSparse1,
          regionSparse2,
          regionSparse3,
          noPermuteRegion3);
      } else {
        returnCode = coinFactorizationA_->updateColumnFT(regionSparse1,
          regionSparse2);
        coinFactorizationA_->updateColumn(regionSparse1,
          regionSparse3,
          noPermuteRegion3);
      }
      coinFactorizationA_->setCollectStatistics(false);
    } else {
#if 0
               CoinSimpFactorization * fact =
                    dynamic_cast< CoinSimpFactorization*>(coinFactorizationB_);
               if (!fact) {
                    returnCode = coinFactorizationB_->updateColumnFT(regionSparse1,
                                 regionSparse2);
                    coinFactorizationB_->updateColumn(regionSparse1,
                                                      regionSparse3,
                                                      noPermuteRegion3);
               } else {
                    returnCode = fact->updateTwoColumnsFT(regionSparse1,
                                                          regionSparse2,
                                                          regionSparse3,
                                                          noPermuteRegion3);
               }
#else
#ifdef CLP_REUSE_ETAS
      int tempInfo[2];
      tempInfo[0] = model_->numberIterations();
      tempInfo[1] = model_->sequenceIn();
      coinFactorizationB_->setUsefulInformation(tempInfo, 3);
#endif
      returnCode = coinFactorizationB_->updateTwoColumnsFT(
        regionSparse1,
        regionSparse2,
        regionSparse3,
        noPermuteRegion3);
#endif
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(9);
#endif
#ifdef PRINT_VECTOR
    printf("UpdateTwoFT\n");
    regionSparse2->print();
    regionSparse3->print();
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
    returnCode = updateColumnFT(regionSparse1, regionSparse2);
    updateColumn(regionSparse1, regionSparse3, noPermuteRegion3);
  }
#endif
  return returnCode;
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnForDebug(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  bool noPermute) const
{
  if (!noPermute)
    regionSparse->checkClear();
  if (!coinFactorizationA_->numberRows())
    return 0;
  coinFactorizationA_->setCollectStatistics(false);
  // above doesn't work any more - do by hand
  double save[15];
  memcpy(save, &coinFactorizationA_->ftranCountInput_, sizeof(save));
  int returnCode = coinFactorizationA_->updateColumn(regionSparse,
    regionSparse2,
    noPermute);
  memcpy(&coinFactorizationA_->ftranCountInput_, save, sizeof(save));
  return returnCode;
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int ClpFactorization::updateColumnTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2) const
{
  if (!numberRows())
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    int returnCode;

    if (coinFactorizationA_) {
      coinFactorizationA_->setCollectStatistics(doStatistics_);
      returnCode = coinFactorizationA_->updateColumnTranspose(regionSparse,
        regionSparse2);
      coinFactorizationA_->setCollectStatistics(false);
    } else {
      returnCode = coinFactorizationB_->updateColumnTranspose(regionSparse,
        regionSparse2);
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(6);
#endif
#ifdef PRINT_VECTOR
    printf("UpdateTranspose\n");
    regionSparse2->print();
#endif
    return returnCode;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector *save = new CoinIndexedVector(*regionSparse2);
    double *COIN_RESTRICT check = new double[coinFactorizationA_->numberRows()];
    int returnCode = coinFactorizationA_->updateColumnTranspose(regionSparse,
      regionSparse2);
    networkBasis_->updateColumnTranspose(regionSparse, save);
    int i;
    double *COIN_RESTRICT array = regionSparse2->denseVector();
    int *COIN_RESTRICT indices = regionSparse2->getIndices();
    int n = regionSparse2->getNumElements();
    memset(check, 0, coinFactorizationA_->numberRows() * sizeof(double));
    double *COIN_RESTRICT array2 = save->denseVector();
    int *COIN_RESTRICT indices2 = save->getIndices();
    int n2 = save->getNumElements();
    assert(n == n2);
    if (save->packedMode()) {
      for (i = 0; i < n; i++) {
        check[indices[i]] = array[i];
      }
      for (i = 0; i < n; i++) {
        double value2 = array2[i];
        assert(check[indices2[i]] == value2);
      }
    } else {
      int numberRows = coinFactorizationA_->numberRows();
      for (i = 0; i < numberRows; i++) {
        double value1 = array[i];
        double value2 = array2[i];
        assert(value1 == value2);
      }
    }
    delete save;
    delete[] check;
    return returnCode;
#else
    return networkBasis_->updateColumnTranspose(regionSparse, regionSparse2);
#endif
  }
#endif
}
/* Updates two columns (BTRAN) from regionSparse2 and 3
   regionSparse starts as zero and is zero at end 
   Note - if regionSparse2 packed on input - will be packed on output - same for 3
*/
void ClpFactorization::updateTwoColumnsTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3) const
{
  if (!numberRows())
    return;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(-1);
#endif
    if (coinFactorizationA_) {
      coinFactorizationA_->setCollectStatistics(doStatistics_);
#if ABOCA_LITE_FACTORIZATION
      coinFactorizationA_->updateTwoColumnsTranspose(regionSparse,
        regionSparse2, regionSparse3, abcState());
#else
      coinFactorizationA_->updateTwoColumnsTranspose(regionSparse,
        regionSparse2, regionSparse3, 0);
#endif
      coinFactorizationA_->setCollectStatistics(false);
    } else {
      coinFactorizationB_->updateColumnTranspose(regionSparse,
        regionSparse2);
      coinFactorizationB_->updateColumnTranspose(regionSparse,
        regionSparse3);
    }
#ifdef CLP_FACTORIZATION_INSTRUMENT
    factorization_instrument(6);
#endif
#ifndef SLIM_CLP
  } else {
    updateColumnTranspose(regionSparse, regionSparse2);
    updateColumnTranspose(regionSparse, regionSparse3);
  }
#endif
}
/* makes a row copy of L for speed and to allow very sparse problems */
void ClpFactorization::goSparse()
{
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
    if (coinFactorizationA_) {
#ifdef CLP_FACTORIZATION_INSTRUMENT
      factorization_instrument(-1);
#endif
      coinFactorizationA_->goSparse();
#ifdef CLP_FACTORIZATION_INSTRUMENT
      factorization_instrument(7);
#endif
    }
  }
}
// Cleans up i.e. gets rid of network basis
void ClpFactorization::cleanUp()
{
#ifndef SLIM_CLP
  delete networkBasis_;
  networkBasis_ = NULL;
#endif
  if (coinFactorizationA_)
    coinFactorizationA_->resetStatistics();
}
/// Says whether to redo pivot order
bool ClpFactorization::needToReorder() const
{
#ifdef CHECK_NETWORK
  return true;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_)
#endif
    return true;
#ifndef SLIM_CLP
  else
    return false;
#endif
}
// Get weighted row list
void ClpFactorization::getWeights(int *COIN_RESTRICT weights) const
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(-1);
#endif
#ifndef SLIM_CLP
  if (networkBasis_) {
    // Network - just unit
    int numberRows = coinFactorizationA_->numberRows();
    for (int i = 0; i < numberRows; i++)
      weights[i] = 1;
    return;
  }
#endif
  int *COIN_RESTRICT numberInRow = coinFactorizationA_->numberInRow();
  int *COIN_RESTRICT numberInColumn = coinFactorizationA_->numberInColumn();
  int *COIN_RESTRICT permuteBack = coinFactorizationA_->pivotColumnBack();
  int *COIN_RESTRICT indexRowU = coinFactorizationA_->indexRowU();
  const int *startColumnU = coinFactorizationA_->startColumnU();
  const int *startRowL = coinFactorizationA_->startRowL();
  int numberRows = coinFactorizationA_->numberRows();
  if (!startRowL || !coinFactorizationA_->numberInRow()) {
    int *COIN_RESTRICT temp = new int[numberRows];
    memset(temp, 0, numberRows * sizeof(int));
    int i;
    for (i = 0; i < numberRows; i++) {
      // one for pivot
      temp[i]++;
      int j;
      for (j = startColumnU[i]; j < startColumnU[i] + numberInColumn[i]; j++) {
        int iRow = indexRowU[j];
        temp[iRow]++;
      }
    }
    int *COIN_RESTRICT startColumnL = coinFactorizationA_->startColumnL();
    int *COIN_RESTRICT indexRowL = coinFactorizationA_->indexRowL();
    int numberL = coinFactorizationA_->numberL();
    int baseL = coinFactorizationA_->baseL();
    for (i = baseL; i < baseL + numberL; i++) {
      int j;
      for (j = startColumnL[i]; j < startColumnL[i + 1]; j++) {
        int iRow = indexRowL[j];
        temp[iRow]++;
      }
    }
    for (i = 0; i < numberRows; i++) {
      int number = temp[i];
      int iPermute = permuteBack[i];
      weights[iPermute] = number;
    }
    delete[] temp;
  } else {
    int i;
    for (i = 0; i < numberRows; i++) {
      int number = startRowL[i + 1] - startRowL[i] + numberInRow[i] + 1;
      //number = startRowL[i+1]-startRowL[i]+1;
      //number = numberInRow[i]+1;
      int iPermute = permuteBack[i];
      weights[iPermute] = number;
    }
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  factorization_instrument(8);
#endif
}
#if ABOCA_LITE_FACTORIZATION
// Does btranU part of replaceColumn (skipping entries)
void ClpFactorization::replaceColumn1(CoinIndexedVector *regionSparse,
  int pivotRow)
{
  if (coinFactorizationA_)
    coinFactorizationA_->replaceColumn1(regionSparse, pivotRow);
}
// Does replaceColumn - having already done btranU
int ClpFactorization::replaceColumn2(CoinIndexedVector *regionSparse,
  int pivotRow,
  double pivotCheck)
{
  if (coinFactorizationA_)
    return coinFactorizationA_->replaceColumn2(regionSparse, pivotRow, pivotCheck);
  else
    return 12345678;
}
#endif
// Set tolerances to safer of existing and given
void ClpFactorization::saferTolerances(double zeroValue,
  double pivotValue)
{
  double newValue;
  // better to have small tolerance even if slower
  if (zeroValue > 0.0)
    newValue = zeroValue;
  else
    newValue = -zeroTolerance() * zeroValue;
  zeroTolerance(std::min(zeroTolerance(), zeroValue));
  // better to have large tolerance even if slower
  if (pivotValue > 0.0)
    newValue = pivotValue;
  else
    newValue = -pivotTolerance() * pivotValue;
  pivotTolerance(std::min(std::max(pivotTolerance(), newValue), 0.999));
}
// Sets factorization
void ClpFactorization::setFactorization(ClpFactorization &rhs)
{
  ClpFactorization::operator=(rhs);
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
