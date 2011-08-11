/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpSimplexPrimal.hpp"
#include "ClpEventHandler.hpp"
#include "ClpHelperFunctions.hpp"
#include "ClpFactorization.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpNonLinearCost.hpp"
#include "ClpDynamicMatrix.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinBuild.hpp"
#include "CoinMpsIO.hpp"
#include "CoinFloatEqual.hpp"
#include "ClpMessage.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream> 
/* Dual ranging.
   This computes increase/decrease in cost for each given variable and corresponding
   sequence numbers which would change basis.  Sequence numbers are 0..numberColumns
   and numberColumns.. for artificials/slacks.
   For non-basic variables the sequence number will be that of the non-basic variables.

   Up to user to provide correct length arrays.

*/
void ClpSimplexOther::dualRanging(int numberCheck, const int * which,
                                  double * costIncreased, int * sequenceIncreased,
                                  double * costDecreased, int * sequenceDecreased,
                                  double * valueIncrease, double * valueDecrease)
{
     rowArray_[1]->clear();
     columnArray_[1]->clear();
     // long enough for rows+columns
     assert(rowArray_[3]->capacity() >= numberRows_ + numberColumns_);
     rowArray_[3]->clear();
     int * backPivot = rowArray_[3]->getIndices();
     int i;
     for ( i = 0; i < numberRows_ + numberColumns_; i++) {
          backPivot[i] = -1;
     }
     for (i = 0; i < numberRows_; i++) {
          int iSequence = pivotVariable_[i];
          backPivot[iSequence] = i;
     }
     // dualTolerance may be zero if from CBC.  In fact use that fact
     bool inCBC = !dualTolerance_;
     if (inCBC)
          assert (integerType_);
     dualTolerance_ = dblParam_[ClpDualTolerance];
     double * arrayX = rowArray_[0]->denseVector();
     for ( i = 0; i < numberCheck; i++) {
          rowArray_[0]->clear();
          //rowArray_[0]->checkClear();
          //rowArray_[1]->checkClear();
          //columnArray_[1]->checkClear();
          columnArray_[0]->clear();
          //columnArray_[0]->checkClear();
          int iSequence = which[i];
          if (iSequence < 0) {
               costIncreased[i] = 0.0;
               sequenceIncreased[i] = -1;
               costDecreased[i] = 0.0;
               sequenceDecreased[i] = -1;
               continue;
          }
          double costIncrease = COIN_DBL_MAX;
          double costDecrease = COIN_DBL_MAX;
          int sequenceIncrease = -1;
          int sequenceDecrease = -1;
          if (valueIncrease) {
               assert (valueDecrease);
               valueIncrease[i] = iSequence < numberColumns_ ? columnActivity_[iSequence] : rowActivity_[iSequence-numberColumns_];
               valueDecrease[i] = valueIncrease[i];
          }

          switch(getStatus(iSequence)) {

          case basic: {
               // non-trvial
               // Get pivot row
               int iRow = backPivot[iSequence];
               assert (iRow >= 0);
               double plusOne = 1.0;
               rowArray_[0]->createPacked(1, &iRow, &plusOne);
               factorization_->updateColumnTranspose(rowArray_[1], rowArray_[0]);
               // put row of tableau in rowArray[0] and columnArray[0]
               matrix_->transposeTimes(this, -1.0,
                                       rowArray_[0], columnArray_[1], columnArray_[0]);
               double alphaIncrease;
               double alphaDecrease;
               // do ratio test up and down
               checkDualRatios(rowArray_[0], columnArray_[0], costIncrease, sequenceIncrease, alphaIncrease,
                               costDecrease, sequenceDecrease, alphaDecrease);
               if (!inCBC) {
                    if (valueIncrease) {
                         if (sequenceIncrease >= 0)
                              valueIncrease[i] = primalRanging1(sequenceIncrease, iSequence);
                         if (sequenceDecrease >= 0)
                              valueDecrease[i] = primalRanging1(sequenceDecrease, iSequence);
                    }
               } else {
                    int number = rowArray_[0]->getNumElements();
                    double scale2 = 0.0;
                    int j;
                    for (j = 0; j < number; j++) {
                         scale2 += arrayX[j] * arrayX[j];
                    }
                    scale2 = 1.0 / sqrt(scale2);
                    //valueIncrease[i] = scale2;
                    if (sequenceIncrease >= 0) {
                         double djValue = dj_[sequenceIncrease];
                         if (fabs(djValue) > 10.0 * dualTolerance_) {
                              // we are going to use for cutoff so be exact
                              costIncrease = fabs(djValue / alphaIncrease);
                              /* Not sure this is good idea as I don't think correct e.g.
                                 suppose a continuous variable has dj slightly greater. */
                              if(false && sequenceIncrease < numberColumns_ && integerType_[sequenceIncrease]) {
                                   // can improve
                                   double movement = (columnScale_ == NULL) ? 1.0 :
                                                     rhsScale_ * inverseColumnScale_[sequenceIncrease];
                                   costIncrease = CoinMax(fabs(djValue * movement), costIncrease);
                              }
                         } else {
                              costIncrease = 0.0;
                         }
                    }
                    if (sequenceDecrease >= 0) {
                         double djValue = dj_[sequenceDecrease];
                         if (fabs(djValue) > 10.0 * dualTolerance_) {
                              // we are going to use for cutoff so be exact
                              costDecrease = fabs(djValue / alphaDecrease);
                              if(sequenceDecrease < numberColumns_ && integerType_[sequenceDecrease]) {
                                   // can improve
                                   double movement = (columnScale_ == NULL) ? 1.0 :
                                                     rhsScale_ * inverseColumnScale_[sequenceDecrease];
                                   costDecrease = CoinMax(fabs(djValue * movement), costDecrease);
                              }
                         } else {
                              costDecrease = 0.0;
                         }
                    }
                    costIncrease *= scale2;
                    costDecrease *= scale2;
               }
          }
          break;
          case isFixed:
               break;
          case isFree:
          case superBasic:
               costIncrease = 0.0;
               costDecrease = 0.0;
               sequenceIncrease = iSequence;
               sequenceDecrease = iSequence;
               break;
          case atUpperBound:
               costIncrease = CoinMax(0.0, -dj_[iSequence]);
               sequenceIncrease = iSequence;
               if (valueIncrease)
                    valueIncrease[i] = primalRanging1(iSequence, iSequence);
               break;
          case atLowerBound:
               costDecrease = CoinMax(0.0, dj_[iSequence]);
               sequenceDecrease = iSequence;
               if (valueIncrease)
                    valueDecrease[i] = primalRanging1(iSequence, iSequence);
               break;
          }
          double scaleFactor;
          if (rowScale_) {
               if (iSequence < numberColumns_)
                    scaleFactor = 1.0 / (objectiveScale_ * columnScale_[iSequence]);
               else
                    scaleFactor = rowScale_[iSequence-numberColumns_] / objectiveScale_;
          } else {
               scaleFactor = 1.0 / objectiveScale_;
          }
          if (costIncrease < 1.0e30)
               costIncrease *= scaleFactor;
          if (costDecrease < 1.0e30)
               costDecrease *= scaleFactor;
          if (optimizationDirection_ == 1.0) {
               costIncreased[i] = costIncrease;
               sequenceIncreased[i] = sequenceIncrease;
               costDecreased[i] = costDecrease;
               sequenceDecreased[i] = sequenceDecrease;
          } else if (optimizationDirection_ == -1.0) {
               costIncreased[i] = costDecrease;
               sequenceIncreased[i] = sequenceDecrease;
               costDecreased[i] = costIncrease;
               sequenceDecreased[i] = sequenceIncrease;
               if (valueIncrease) {
                    double temp = valueIncrease[i];
                    valueIncrease[i] = valueDecrease[i];
                    valueDecrease[i] = temp;
               }
          } else if (optimizationDirection_ == 0.0) {
               // !!!!!! ???
               costIncreased[i] = COIN_DBL_MAX;
               sequenceIncreased[i] = -1;
               costDecreased[i] = COIN_DBL_MAX;
               sequenceDecreased[i] = -1;
          } else {
               abort();
          }
     }
     rowArray_[0]->clear();
     //rowArray_[1]->clear();
     //columnArray_[1]->clear();
     columnArray_[0]->clear();
     //rowArray_[3]->clear();
     if (!optimizationDirection_)
          printf("*** ????? Ranging with zero optimization costs\n");
}
/*
   Row array has row part of pivot row
   Column array has column part.
   This is used in dual ranging
*/
void
ClpSimplexOther::checkDualRatios(CoinIndexedVector * rowArray,
                                 CoinIndexedVector * columnArray,
                                 double & costIncrease, int & sequenceIncrease, double & alphaIncrease,
                                 double & costDecrease, int & sequenceDecrease, double & alphaDecrease)
{
     double acceptablePivot = 1.0e-9;
     double * work;
     int number;
     int * which;
     int iSection;

     double thetaDown = 1.0e31;
     double thetaUp = 1.0e31;
     int sequenceDown = -1;
     int sequenceUp = -1;
     double alphaDown = 0.0;
     double alphaUp = 0.0;

     int addSequence;

     for (iSection = 0; iSection < 2; iSection++) {

          int i;
          if (!iSection) {
               work = rowArray->denseVector();
               number = rowArray->getNumElements();
               which = rowArray->getIndices();
               addSequence = numberColumns_;
          } else {
               work = columnArray->denseVector();
               number = columnArray->getNumElements();
               which = columnArray->getIndices();
               addSequence = 0;
          }

          for (i = 0; i < number; i++) {
               int iSequence = which[i];
               int iSequence2 = iSequence + addSequence;
               double alpha = work[i];
               if (fabs(alpha) < acceptablePivot)
                    continue;
               double oldValue = dj_[iSequence2];

               switch(getStatus(iSequence2)) {

               case basic:
                    break;
               case ClpSimplex::isFixed:
                    break;
               case isFree:
               case superBasic:
                    // treat dj as if zero
                    thetaDown = 0.0;
                    thetaUp = 0.0;
                    sequenceDown = iSequence2;
                    sequenceUp = iSequence2;
                    break;
               case atUpperBound:
                    if (alpha > 0.0) {
                         // test up
                         if (oldValue + thetaUp * alpha > dualTolerance_) {
                              thetaUp = (dualTolerance_ - oldValue) / alpha;
                              sequenceUp = iSequence2;
                              alphaUp = alpha;
                         }
                    } else {
                         // test down
                         if (oldValue - thetaDown * alpha > dualTolerance_) {
                              thetaDown = -(dualTolerance_ - oldValue) / alpha;
                              sequenceDown = iSequence2;
                              alphaDown = alpha;
                         }
                    }
                    break;
               case atLowerBound:
                    if (alpha < 0.0) {
                         // test up
                         if (oldValue + thetaUp * alpha < - dualTolerance_) {
                              thetaUp = -(dualTolerance_ + oldValue) / alpha;
                              sequenceUp = iSequence2;
                              alphaUp = alpha;
                         }
                    } else {
                         // test down
                         if (oldValue - thetaDown * alpha < -dualTolerance_) {
                              thetaDown = (dualTolerance_ + oldValue) / alpha;
                              sequenceDown = iSequence2;
                              alphaDown = alpha;
                         }
                    }
                    break;
               }
          }
     }
     if (sequenceUp >= 0) {
          costIncrease = thetaUp;
          sequenceIncrease = sequenceUp;
          alphaIncrease = alphaUp;
     }
     if (sequenceDown >= 0) {
          costDecrease = thetaDown;
          sequenceDecrease = sequenceDown;
          alphaDecrease = alphaDown;
     }
}
/** Primal ranging.
    This computes increase/decrease in value for each given variable and corresponding
    sequence numbers which would change basis.  Sequence numbers are 0..numberColumns
    and numberColumns.. for artificials/slacks.
    For basic variables the sequence number will be that of the basic variables.

    Up to user to provide correct length arrays.

    When here - guaranteed optimal
*/
void
ClpSimplexOther::primalRanging(int numberCheck, const int * which,
                               double * valueIncreased, int * sequenceIncreased,
                               double * valueDecreased, int * sequenceDecreased)
{
     rowArray_[0]->clear();
     rowArray_[1]->clear();
     lowerIn_ = -COIN_DBL_MAX;
     upperIn_ = COIN_DBL_MAX;
     valueIn_ = 0.0;
     for ( int i = 0; i < numberCheck; i++) {
          int iSequence = which[i];
          double valueIncrease = COIN_DBL_MAX;
          double valueDecrease = COIN_DBL_MAX;
          int sequenceIncrease = -1;
          int sequenceDecrease = -1;

          switch(getStatus(iSequence)) {

          case basic:
          case isFree:
          case superBasic:
               // Easy
               valueDecrease = CoinMax(0.0, upper_[iSequence] - solution_[iSequence]);
               valueIncrease = CoinMax(0.0, solution_[iSequence] - lower_[iSequence]);
               sequenceDecrease = iSequence;
               sequenceIncrease = iSequence;
               break;
          case isFixed:
          case atUpperBound:
          case atLowerBound: {
               // Non trivial
               // Other bound is ignored
               unpackPacked(rowArray_[1], iSequence);
               factorization_->updateColumn(rowArray_[2], rowArray_[1]);
               // Get extra rows
               matrix_->extendUpdated(this, rowArray_[1], 0);
               // do ratio test
               checkPrimalRatios(rowArray_[1], 1);
               if (pivotRow_ >= 0) {
                    valueIncrease = theta_;
                    sequenceIncrease = pivotVariable_[pivotRow_];
               }
               checkPrimalRatios(rowArray_[1], -1);
               if (pivotRow_ >= 0) {
                    valueDecrease = theta_;
                    sequenceDecrease = pivotVariable_[pivotRow_];
               }
               rowArray_[1]->clear();
          }
          break;
          }
          double scaleFactor;
          if (rowScale_) {
               if (iSequence < numberColumns_)
                    scaleFactor = columnScale_[iSequence] / rhsScale_;
               else
                    scaleFactor = 1.0 / (rowScale_[iSequence-numberColumns_] * rhsScale_);
          } else {
               scaleFactor = 1.0 / rhsScale_;
          }
          if (valueIncrease < 1.0e30)
               valueIncrease *= scaleFactor;
          else
               valueIncrease = COIN_DBL_MAX;
          if (valueDecrease < 1.0e30)
               valueDecrease *= scaleFactor;
          else
               valueDecrease = COIN_DBL_MAX;
          valueIncreased[i] = valueIncrease;
          sequenceIncreased[i] = sequenceIncrease;
          valueDecreased[i] = valueDecrease;
          sequenceDecreased[i] = sequenceDecrease;
     }
}
// Returns new value of whichOther when whichIn enters basis
double
ClpSimplexOther::primalRanging1(int whichIn, int whichOther)
{
     rowArray_[0]->clear();
     rowArray_[1]->clear();
     int iSequence = whichIn;
     double newValue = solution_[whichOther];
     double alphaOther = 0.0;
     Status status = getStatus(iSequence);
     assert (status == atLowerBound || status == atUpperBound);
     int wayIn = (status == atLowerBound) ? 1 : -1;

     switch(getStatus(iSequence)) {

     case basic:
     case isFree:
     case superBasic:
          assert (whichIn == whichOther);
          // Easy
          newValue = wayIn > 0 ? upper_[iSequence] : lower_[iSequence];
          break;
     case isFixed:
     case atUpperBound:
     case atLowerBound:
          // Non trivial
     {
          // Other bound is ignored
          unpackPacked(rowArray_[1], iSequence);
          factorization_->updateColumn(rowArray_[2], rowArray_[1]);
          // Get extra rows
          matrix_->extendUpdated(this, rowArray_[1], 0);
          // do ratio test
          double acceptablePivot = 1.0e-7;
          double * work = rowArray_[1]->denseVector();
          int number = rowArray_[1]->getNumElements();
          int * which = rowArray_[1]->getIndices();

          // we may need to swap sign
          double way = wayIn;
          double theta = 1.0e30;
          for (int iIndex = 0; iIndex < number; iIndex++) {

               int iRow = which[iIndex];
               double alpha = work[iIndex] * way;
               int iPivot = pivotVariable_[iRow];
               if (iPivot == whichOther) {
                    alphaOther = alpha;
                    continue;
               }
               double oldValue = solution_[iPivot];
               if (fabs(alpha) > acceptablePivot) {
                    if (alpha > 0.0) {
                         // basic variable going towards lower bound
                         double bound = lower_[iPivot];
                         oldValue -= bound;
                         if (oldValue - theta * alpha < 0.0) {
                              theta = CoinMax(0.0, oldValue / alpha);
                         }
                    } else {
                         // basic variable going towards upper bound
                         double bound = upper_[iPivot];
                         oldValue = oldValue - bound;
                         if (oldValue - theta * alpha > 0.0) {
                              theta = CoinMax(0.0, oldValue / alpha);
                         }
                    }
               }
          }
          if (whichIn != whichOther) {
               if (theta < 1.0e30)
                    newValue -= theta * alphaOther;
               else
                    newValue = alphaOther > 0.0 ? -1.0e30 : 1.0e30;
          } else {
               newValue += theta * wayIn;
          }
     }
     rowArray_[1]->clear();
     break;
     }
     double scaleFactor;
     if (rowScale_) {
          if (whichOther < numberColumns_)
               scaleFactor = columnScale_[whichOther] / rhsScale_;
          else
               scaleFactor = 1.0 / (rowScale_[whichOther-numberColumns_] * rhsScale_);
     } else {
          scaleFactor = 1.0 / rhsScale_;
     }
     if (newValue < 1.0e29)
          if (newValue > -1.0e29)
               newValue *= scaleFactor;
          else
               newValue = -COIN_DBL_MAX;
     else
          newValue = COIN_DBL_MAX;
     return newValue;
}
/*
   Row array has pivot column
   This is used in primal ranging
*/
void
ClpSimplexOther::checkPrimalRatios(CoinIndexedVector * rowArray,
                                   int direction)
{
     // sequence stays as row number until end
     pivotRow_ = -1;
     double acceptablePivot = 1.0e-7;
     double * work = rowArray->denseVector();
     int number = rowArray->getNumElements();
     int * which = rowArray->getIndices();

     // we need to swap sign if going down
     double way = direction;
     theta_ = 1.0e30;
     for (int iIndex = 0; iIndex < number; iIndex++) {

          int iRow = which[iIndex];
          double alpha = work[iIndex] * way;
          int iPivot = pivotVariable_[iRow];
          double oldValue = solution_[iPivot];
          if (fabs(alpha) > acceptablePivot) {
               if (alpha > 0.0) {
                    // basic variable going towards lower bound
                    double bound = lower_[iPivot];
                    oldValue -= bound;
                    if (oldValue - theta_ * alpha < 0.0) {
                         pivotRow_ = iRow;
                         theta_ = CoinMax(0.0, oldValue / alpha);
                    }
               } else {
                    // basic variable going towards upper bound
                    double bound = upper_[iPivot];
                    oldValue = oldValue - bound;
                    if (oldValue - theta_ * alpha > 0.0) {
                         pivotRow_ = iRow;
                         theta_ = CoinMax(0.0, oldValue / alpha);
                    }
               }
          }
     }
}
/* Write the basis in MPS format to the specified file.
   If writeValues true writes values of structurals
   (and adds VALUES to end of NAME card)

   Row and column names may be null.
   formatType is
   <ul>
   <li> 0 - normal
   <li> 1 - extra accuracy
   <li> 2 - IEEE hex (later)
   </ul>

   Returns non-zero on I/O error

   This is based on code contributed by Thorsten Koch
*/
int
ClpSimplexOther::writeBasis(const char *filename,
                            bool writeValues,
                            int formatType) const
{
     formatType = CoinMax(0, formatType);
     formatType = CoinMin(2, formatType);
     if (!writeValues)
          formatType = 0;
     // See if INTEL if IEEE
     if (formatType == 2) {
          // test intel here and add 1 if not intel
          double value = 1.0;
          char x[8];
          memcpy(x, &value, 8);
          if (x[0] == 63) {
               formatType ++; // not intel
          } else {
               assert (x[0] == 0);
          }
     }

     char number[20];
     FILE * fp = fopen(filename, "w");
     if (!fp)
          return -1;

     // NAME card

     if (strcmp(strParam_[ClpProbName].c_str(), "") == 0) {
          fprintf(fp, "NAME          BLANK      ");
     } else {
          fprintf(fp, "NAME          %s       ", strParam_[ClpProbName].c_str());
     }
     if (formatType >= 2)
          fprintf(fp, "FREEIEEE");
     else if (writeValues)
          fprintf(fp, "VALUES");
     // finish off name
     fprintf(fp, "\n");
     int iRow = 0;
     for(int iColumn = 0; iColumn < numberColumns_; iColumn++) {
          bool printit = false;
          if( getColumnStatus(iColumn) == ClpSimplex::basic) {
               printit = true;
               // Find non basic row
               for(; iRow < numberRows_; iRow++) {
                    if (getRowStatus(iRow) != ClpSimplex::basic)
                         break;
               }
               if (lengthNames_) {
                    if (iRow != numberRows_) {
                         fprintf(fp, " %s %-8s       %s",
                                 getRowStatus(iRow) == ClpSimplex::atUpperBound ? "XU" : "XL",
                                 columnNames_[iColumn].c_str(),
                                 rowNames_[iRow].c_str());
                         iRow++;
                    } else {
                         // Allow for too many basics!
                         fprintf(fp, " BS %-8s       ",
                                 columnNames_[iColumn].c_str());
                         // Dummy row name if values
                         if (writeValues)
                              fprintf(fp, "      _dummy_");
                    }
               } else {
                    // no names
                    if (iRow != numberRows_) {
                         fprintf(fp, " %s C%7.7d     R%7.7d",
                                 getRowStatus(iRow) == ClpSimplex::atUpperBound ? "XU" : "XL",
                                 iColumn, iRow);
                         iRow++;
                    } else {
                         // Allow for too many basics!
                         fprintf(fp, " BS C%7.7d", iColumn);
                         // Dummy row name if values
                         if (writeValues)
                              fprintf(fp, "      _dummy_");
                    }
               }
          } else  {
               if( getColumnStatus(iColumn) == ClpSimplex::atUpperBound) {
                    printit = true;
                    if (lengthNames_)
                         fprintf(fp, " UL %s", columnNames_[iColumn].c_str());
                    else
                         fprintf(fp, " UL C%7.7d", iColumn);
                    // Dummy row name if values
                    if (writeValues)
                         fprintf(fp, "      _dummy_");
               }
          }
          if (printit && writeValues) {
               // add value
               CoinConvertDouble(0, formatType, columnActivity_[iColumn], number);
               fprintf(fp, "     %s", number);
          }
          if (printit)
               fprintf(fp, "\n");
     }
     fprintf(fp, "ENDATA\n");
     fclose(fp);
     return 0;
}
// Read a basis from the given filename
int
ClpSimplexOther::readBasis(const char *fileName)
{
     int status = 0;
     bool canOpen = false;
     if (!strcmp(fileName, "-") || !strcmp(fileName, "stdin")) {
          // stdin
          canOpen = true;
     } else {
          FILE *fp = fopen(fileName, "r");
          if (fp) {
               // can open - lets go for it
               fclose(fp);
               canOpen = true;
          } else {
               handler_->message(CLP_UNABLE_OPEN, messages_)
                         << fileName << CoinMessageEol;
               return -1;
          }
     }
     CoinMpsIO m;
     m.passInMessageHandler(handler_);
     *m.messagesPointer() = coinMessages();
     bool savePrefix = m.messageHandler()->prefix();
     m.messageHandler()->setPrefix(handler_->prefix());
     status = m.readBasis(fileName, "", columnActivity_, status_ + numberColumns_,
                          status_,
                          columnNames_, numberColumns_,
                          rowNames_, numberRows_);
     m.messageHandler()->setPrefix(savePrefix);
     if (status >= 0) {
          if (!status) {
               // set values
               int iColumn, iRow;
               for (iRow = 0; iRow < numberRows_; iRow++) {
                    if (getRowStatus(iRow) == atLowerBound)
                         rowActivity_[iRow] = rowLower_[iRow];
                    else if (getRowStatus(iRow) == atUpperBound)
                         rowActivity_[iRow] = rowUpper_[iRow];
               }
               for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
                    if (getColumnStatus(iColumn) == atLowerBound)
                         columnActivity_[iColumn] = columnLower_[iColumn];
                    else if (getColumnStatus(iColumn) == atUpperBound)
                         columnActivity_[iColumn] = columnUpper_[iColumn];
               }
          } else {
               memset(rowActivity_, 0, numberRows_ * sizeof(double));
               matrix_->times(-1.0, columnActivity_, rowActivity_);
          }
     } else {
          // errors
          handler_->message(CLP_IMPORT_ERRORS, messages_)
                    << status << fileName << CoinMessageEol;
     }
     return status;
}
/* Creates dual of a problem if looks plausible
   (defaults will always create model)
   fractionRowRanges is fraction of rows allowed to have ranges
   fractionColumnRanges is fraction of columns allowed to have ranges
*/
ClpSimplex *
ClpSimplexOther::dualOfModel(double fractionRowRanges, double fractionColumnRanges) const
{
     const ClpSimplex * model2 = static_cast<const ClpSimplex *> (this);
     bool changed = false;
     int numberChanged = 0;
     int iColumn;
     // check if we need to change bounds to rows
     for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          if (columnUpper_[iColumn] < 1.0e20 &&
                    columnLower_[iColumn] > -1.0e20) {
               changed = true;
               numberChanged++;
          }
     }
     int iRow;
     int numberExtraRows = 0;
     if (numberChanged <= fractionColumnRanges * numberColumns_) {
          for (iRow = 0; iRow < numberRows_; iRow++) {
               if (rowLower_[iRow] > -1.0e20 &&
                         rowUpper_[iRow] < 1.0e20) {
                    if (rowUpper_[iRow] != rowLower_[iRow])
                         numberExtraRows++;
               }
          }
          if (numberExtraRows > fractionRowRanges * numberRows_)
               return NULL;
     } else {
          return NULL;
     }
     if (changed) {
          ClpSimplex * model3 = new ClpSimplex(*model2);
          CoinBuild build;
          double one = 1.0;
          int numberColumns = model3->numberColumns();
          const double * columnLower = model3->columnLower();
          const double * columnUpper = model3->columnUpper();
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
               if (columnUpper[iColumn] < 1.0e20 &&
                         columnLower[iColumn] > -1.0e20) {
                    if (fabs(columnLower[iColumn]) < fabs(columnUpper[iColumn])) {
                         double value = columnUpper[iColumn];
                         model3->setColumnUpper(iColumn, COIN_DBL_MAX);
                         build.addRow(1, &iColumn, &one, -COIN_DBL_MAX, value);
                    } else {
                         double value = columnLower[iColumn];
                         model3->setColumnLower(iColumn, -COIN_DBL_MAX);
                         build.addRow(1, &iColumn, &one, value, COIN_DBL_MAX);
                    }
               }
          }
          model3->addRows(build);
          model2 = model3;
     }
     int numberColumns = model2->numberColumns();
     const double * columnLower = model2->columnLower();
     const double * columnUpper = model2->columnUpper();
     int numberRows = model2->numberRows();
     double * rowLower = CoinCopyOfArray(model2->rowLower(), numberRows);
     double * rowUpper = CoinCopyOfArray(model2->rowUpper(), numberRows);

     const double * objective = model2->objective();
     CoinPackedMatrix * matrix = model2->matrix();
     // get transpose
     CoinPackedMatrix rowCopy = *matrix;
     const int * row = matrix->getIndices();
     const int * columnLength = matrix->getVectorLengths();
     const CoinBigIndex * columnStart = matrix->getVectorStarts();
     const double * elementByColumn = matrix->getElements();
     double objOffset = 0.0;
     for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          double offset = 0.0;
          double objValue = optimizationDirection_ * objective[iColumn];
          if (columnUpper[iColumn] > 1.0e20) {
               if (columnLower[iColumn] > -1.0e20)
                    offset = columnLower[iColumn];
          } else if (columnLower[iColumn] < -1.0e20) {
               offset = columnUpper[iColumn];
          } else {
               // taken care of before
               abort();
          }
          if (offset) {
               objOffset += offset * objValue;
               for (CoinBigIndex j = columnStart[iColumn];
                         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    if (rowLower[iRow] > -1.0e20)
                         rowLower[iRow] -= offset * elementByColumn[j];
                    if (rowUpper[iRow] < 1.0e20)
                         rowUpper[iRow] -= offset * elementByColumn[j];
               }
          }
     }
     int * which = new int[numberRows+numberExtraRows];
     rowCopy.reverseOrdering();
     rowCopy.transpose();
     double * fromRowsLower = new double[numberRows+numberExtraRows];
     double * fromRowsUpper = new double[numberRows+numberExtraRows];
     double * newObjective = new double[numberRows+numberExtraRows];
     double * fromColumnsLower = new double[numberColumns];
     double * fromColumnsUpper = new double[numberColumns];
     for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          double objValue = optimizationDirection_ * objective[iColumn];
          // Offset is already in
          if (columnUpper[iColumn] > 1.0e20) {
               if (columnLower[iColumn] > -1.0e20) {
                    fromColumnsLower[iColumn] = -COIN_DBL_MAX;
                    fromColumnsUpper[iColumn] = objValue;
               } else {
                    // free
                    fromColumnsLower[iColumn] = objValue;
                    fromColumnsUpper[iColumn] = objValue;
               }
          } else if (columnLower[iColumn] < -1.0e20) {
               fromColumnsLower[iColumn] = objValue;
               fromColumnsUpper[iColumn] = COIN_DBL_MAX;
          } else {
               abort();
          }
     }
     int kRow = 0;
     int kExtraRow = numberRows;
     for (iRow = 0; iRow < numberRows; iRow++) {
          if (rowLower[iRow] < -1.0e20) {
               assert (rowUpper[iRow] < 1.0e20);
               newObjective[kRow] = -rowUpper[iRow];
               fromRowsLower[kRow] = -COIN_DBL_MAX;
               fromRowsUpper[kRow] = 0.0;
               which[kRow] = iRow;
               kRow++;
          } else if (rowUpper[iRow] > 1.0e20) {
               newObjective[kRow] = -rowLower[iRow];
               fromRowsLower[kRow] = 0.0;
               fromRowsUpper[kRow] = COIN_DBL_MAX;
               which[kRow] = iRow;
               kRow++;
          } else {
               if (rowUpper[iRow] == rowLower[iRow]) {
                    newObjective[kRow] = -rowLower[iRow];
                    fromRowsLower[kRow] = -COIN_DBL_MAX;;
                    fromRowsUpper[kRow] = COIN_DBL_MAX;
                    which[kRow] = iRow;
                    kRow++;
               } else {
                    // range
                    newObjective[kRow] = -rowUpper[iRow];
                    fromRowsLower[kRow] = -COIN_DBL_MAX;
                    fromRowsUpper[kRow] = 0.0;
                    which[kRow] = iRow;
                    kRow++;
                    newObjective[kExtraRow] = -rowLower[iRow];
                    fromRowsLower[kExtraRow] = 0.0;
                    fromRowsUpper[kExtraRow] = COIN_DBL_MAX;
                    which[kExtraRow] = iRow;
                    kExtraRow++;
               }
          }
     }
     if (numberExtraRows) {
          CoinPackedMatrix newCopy;
          newCopy.setExtraGap(0.0);
          newCopy.setExtraMajor(0.0);
          newCopy.submatrixOfWithDuplicates(rowCopy, kExtraRow, which);
          rowCopy = newCopy;
     }
     ClpSimplex * modelDual = new ClpSimplex();
     modelDual->loadProblem(rowCopy, fromRowsLower, fromRowsUpper, newObjective,
                            fromColumnsLower, fromColumnsUpper);
     modelDual->setObjectiveOffset(objOffset);
     modelDual->setDualBound(model2->dualBound());
     modelDual->setInfeasibilityCost(model2->infeasibilityCost());
     modelDual->setDualTolerance(model2->dualTolerance());
     modelDual->setPrimalTolerance(model2->primalTolerance());
     modelDual->setPerturbation(model2->perturbation());
     modelDual->setSpecialOptions(model2->specialOptions());
     modelDual->setMoreSpecialOptions(model2->moreSpecialOptions());
     delete [] fromRowsLower;
     delete [] fromRowsUpper;
     delete [] fromColumnsLower;
     delete [] fromColumnsUpper;
     delete [] newObjective;
     delete [] which;
     delete [] rowLower;
     delete [] rowUpper;
     if (changed)
          delete model2;
     modelDual->createStatus();
     return modelDual;
}
// Restores solution from dualized problem
int
ClpSimplexOther::restoreFromDual(const ClpSimplex * dualProblem)
{
     int returnCode = 0;;
     createStatus();
     // Number of rows in dual problem was original number of columns
     assert (numberColumns_ == dualProblem->numberRows());
     // If slack on d-row basic then column at bound otherwise column basic
     // If d-column basic then rhs tight
     int numberBasic = 0;
     int iRow, iColumn = 0;
     // Get number of extra rows from ranges
     int numberExtraRows = 0;
     for (iRow = 0; iRow < numberRows_; iRow++) {
          if (rowLower_[iRow] > -1.0e20 &&
                    rowUpper_[iRow] < 1.0e20) {
               if (rowUpper_[iRow] != rowLower_[iRow])
                    numberExtraRows++;
          }
     }
     const double * objective = this->objective();
     const double * dualDual = dualProblem->dualRowSolution();
     const double * dualDj = dualProblem->dualColumnSolution();
     const double * dualSol = dualProblem->primalColumnSolution();
     const double * dualActs = dualProblem->primalRowSolution();
#if 0
     ClpSimplex thisCopy = *this;
     thisCopy.dual(); // for testing
     const double * primalDual = thisCopy.dualRowSolution();
     const double * primalDj = thisCopy.dualColumnSolution();
     const double * primalSol = thisCopy.primalColumnSolution();
     const double * primalActs = thisCopy.primalRowSolution();
     char ss[] = {'F', 'B', 'U', 'L', 'S', 'F'};
     printf ("Dual problem row info %d rows\n", dualProblem->numberRows());
     for (iRow = 0; iRow < dualProblem->numberRows(); iRow++)
          printf("%d at %c primal %g dual %g\n",
                 iRow, ss[dualProblem->getRowStatus(iRow)],
                 dualActs[iRow], dualDual[iRow]);
     printf ("Dual problem column info %d columns\n", dualProblem->numberColumns());
     for (iColumn = 0; iColumn < dualProblem->numberColumns(); iColumn++)
          printf("%d at %c primal %g dual %g\n",
                 iColumn, ss[dualProblem->getColumnStatus(iColumn)],
                 dualSol[iColumn], dualDj[iColumn]);
     printf ("Primal problem row info %d rows\n", thisCopy.numberRows());
     for (iRow = 0; iRow < thisCopy.numberRows(); iRow++)
          printf("%d at %c primal %g dual %g\n",
                 iRow, ss[thisCopy.getRowStatus(iRow)],
                 primalActs[iRow], primalDual[iRow]);
     printf ("Primal problem column info %d columns\n", thisCopy.numberColumns());
     for (iColumn = 0; iColumn < thisCopy.numberColumns(); iColumn++)
          printf("%d at %c primal %g dual %g\n",
                 iColumn, ss[thisCopy.getColumnStatus(iColumn)],
                 primalSol[iColumn], primalDj[iColumn]);
#endif
     // position at bound information
     int jColumn = numberRows_;
     for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double objValue = optimizationDirection_ * objective[iColumn];
          Status status = dualProblem->getRowStatus(iColumn);
          double otherValue = COIN_DBL_MAX;
          if (columnUpper_[iColumn] < 1.0e20 &&
                    columnLower_[iColumn] > -1.0e20) {
               if (fabs(columnLower_[iColumn]) < fabs(columnUpper_[iColumn])) {
                    otherValue = columnUpper_[iColumn] + dualDj[jColumn];
               } else {
                    otherValue = columnLower_[iColumn] + dualDj[jColumn];
               }
               jColumn++;
          }
          if (status == basic) {
               // column is at bound
               if (otherValue == COIN_DBL_MAX) {
                    reducedCost_[iColumn] = objValue - dualActs[iColumn];
                    if (columnUpper_[iColumn] > 1.0e20) {
                         if (columnLower_[iColumn] > -1.0e20) {
                              if (columnUpper_[iColumn] > columnLower_[iColumn])
                                   setColumnStatus(iColumn, atLowerBound);
                              else
                                   setColumnStatus(iColumn, isFixed);
                              columnActivity_[iColumn] = columnLower_[iColumn];
                         } else {
                              // free
                              setColumnStatus(iColumn, isFree);
                              columnActivity_[iColumn] = 0.0;
                         }
                    } else {
                         setColumnStatus(iColumn, atUpperBound);
                         columnActivity_[iColumn] = columnUpper_[iColumn];
                    }
               } else {
                    reducedCost_[iColumn] = objValue - dualActs[iColumn];
                    //printf("other dual sol %g\n",otherValue);
                    if (fabs(otherValue - columnLower_[iColumn]) < 1.0e-5) {
                         if (columnUpper_[iColumn] > columnLower_[iColumn])
                              setColumnStatus(iColumn, atLowerBound);
                         else
                              setColumnStatus(iColumn, isFixed);
                         columnActivity_[iColumn] = columnLower_[iColumn];
                    } else if (fabs(otherValue - columnUpper_[iColumn]) < 1.0e-5) {
                         if (columnUpper_[iColumn] > columnLower_[iColumn])
                              setColumnStatus(iColumn, atUpperBound);
                         else
                              setColumnStatus(iColumn, isFixed);
                         columnActivity_[iColumn] = columnUpper_[iColumn];
                    } else {
                         abort();
                    }
               }
          } else {
               if (otherValue == COIN_DBL_MAX) {
                    // column basic
                    setColumnStatus(iColumn, basic);
                    numberBasic++;
                    if (columnLower_[iColumn] > -1.0e20) {
                         columnActivity_[iColumn] = -dualDual[iColumn] + columnLower_[iColumn];
                    } else if (columnUpper_[iColumn] < 1.0e20) {
                         columnActivity_[iColumn] = -dualDual[iColumn] + columnUpper_[iColumn];
                    } else {
                         columnActivity_[iColumn] = -dualDual[iColumn];
                    }
                    reducedCost_[iColumn] = 0.0;
               } else {
                    // may be at other bound
                    //printf("xx %d %g jcol %d\n",iColumn,otherValue,jColumn-1);
                    if (dualProblem->getColumnStatus(jColumn - 1) != basic) {
                         // column basic
                         setColumnStatus(iColumn, basic);
                         numberBasic++;
                         //printf("Col %d otherV %g dualDual %g\n",iColumn,
                         // otherValue,dualDual[iColumn]);
                         columnActivity_[iColumn] = -dualDual[iColumn];
                         columnActivity_[iColumn] = otherValue;
                         reducedCost_[iColumn] = 0.0;
                    } else {
                         reducedCost_[iColumn] = objValue - dualActs[iColumn];
                         if (fabs(otherValue - columnLower_[iColumn]) < 1.0e-5) {
                              if (columnUpper_[iColumn] > columnLower_[iColumn])
                                   setColumnStatus(iColumn, atLowerBound);
                              else
                                   setColumnStatus(iColumn, isFixed);
                              columnActivity_[iColumn] = columnLower_[iColumn];
                         } else if (fabs(otherValue - columnUpper_[iColumn]) < 1.0e-5) {
                              if (columnUpper_[iColumn] > columnLower_[iColumn])
                                   setColumnStatus(iColumn, atUpperBound);
                              else
                                   setColumnStatus(iColumn, isFixed);
                              columnActivity_[iColumn] = columnUpper_[iColumn];
                         } else {
                              abort();
                         }
                    }
               }
          }
     }
     // now rows
     int kExtraRow = jColumn;
     int numberRanges = 0;
     for (iRow = 0; iRow < numberRows_; iRow++) {
          Status status = dualProblem->getColumnStatus(iRow);
          if (status == basic) {
               // row is at bound
               dual_[iRow] = dualSol[iRow];;
          } else {
               // row basic
               setRowStatus(iRow, basic);
               numberBasic++;
               dual_[iRow] = 0.0;
          }
          if (rowLower_[iRow] < -1.0e20) {
               if (status == basic) {
                    rowActivity_[iRow] = rowUpper_[iRow];
                    setRowStatus(iRow, atUpperBound);
               } else {
                    assert (dualDj[iRow] < 1.0e-5);
                    rowActivity_[iRow] = rowUpper_[iRow] + dualDj[iRow];
               }
          } else if (rowUpper_[iRow] > 1.0e20) {
               if (status == basic) {
                    rowActivity_[iRow] = rowLower_[iRow];
                    setRowStatus(iRow, atLowerBound);
               } else {
                    rowActivity_[iRow] = rowLower_[iRow] + dualDj[iRow];
                    assert (dualDj[iRow] > -1.0e-5);
               }
          } else {
               if (rowUpper_[iRow] == rowLower_[iRow]) {
                    rowActivity_[iRow] = rowLower_[iRow];
                    if (status == basic) {
                         setRowStatus(iRow, isFixed);
                    }
               } else {
                    // range
                    numberRanges++;
                    Status statusL = dualProblem->getColumnStatus(kExtraRow);
                    //printf("range row %d (%d), extra %d (%d) - dualSol %g,%g dualDj %g,%g\n",
                    //     iRow,status,kExtraRow,statusL, dualSol[iRow],
                    //     dualSol[kExtraRow],dualDj[iRow],dualDj[kExtraRow]);
                    if (status == basic) {
                         assert (statusL != basic);
                         rowActivity_[iRow] = rowUpper_[iRow];
                         setRowStatus(iRow, atUpperBound);
                    } else if (statusL == basic) {
                         numberBasic--; // already counted
                         rowActivity_[iRow] = rowLower_[iRow];
                         setRowStatus(iRow, atLowerBound);
                         dual_[iRow] = dualSol[kExtraRow];;
                    } else {
                         rowActivity_[iRow] = rowLower_[iRow] - dualDj[iRow];
                         assert (dualDj[iRow] < 1.0e-5);
                         // row basic
                         //setRowStatus(iRow,basic);
                         //numberBasic++;
                         dual_[iRow] = 0.0;
                    }
                    kExtraRow++;
               }
          }
     }
     if (numberBasic != numberRows_) {
          printf("Bad basis - ranges - coding needed\n");
          assert (numberRanges);
          abort();
     }
     if (optimizationDirection_ < 0.0) {
          for (iRow = 0; iRow < numberRows_; iRow++) {
               dual_[iRow] = -dual_[iRow];
          }
     }
     // redo row activities
     memset(rowActivity_, 0, numberRows_ * sizeof(double));
     matrix_->times(1.0, columnActivity_, rowActivity_);
     // redo reduced costs
     memcpy(reducedCost_, this->objective(), numberColumns_ * sizeof(double));
     matrix_->transposeTimes(-1.0, dual_, reducedCost_);
     checkSolutionInternal();
     if (sumDualInfeasibilities_ > 1.0e-5 || sumPrimalInfeasibilities_ > 1.0e-5) {
          returnCode = 1;
#ifdef CLP_INVESTIGATE
          printf("There are %d dual infeasibilities summing to %g ",
                 numberDualInfeasibilities_, sumDualInfeasibilities_);
          printf("and %d primal infeasibilities summing to %g\n",
                 numberPrimalInfeasibilities_, sumPrimalInfeasibilities_);
#endif
     }
     // Below will go to ..DEBUG later
#if 1 //ndef NDEBUG
     // Check if correct
     double * columnActivity = CoinCopyOfArray(columnActivity_, numberColumns_);
     double * rowActivity = CoinCopyOfArray(rowActivity_, numberRows_);
     double * reducedCost = CoinCopyOfArray(reducedCost_, numberColumns_);
     double * dual = CoinCopyOfArray(dual_, numberRows_);
     this->dual(); //primal();
     CoinRelFltEq eq(1.0e-5);
     for (iRow = 0; iRow < numberRows_; iRow++) {
          assert(eq(dual[iRow], dual_[iRow]));
     }
     for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          assert(eq(columnActivity[iColumn], columnActivity_[iColumn]));
     }
     for (iRow = 0; iRow < numberRows_; iRow++) {
          assert(eq(rowActivity[iRow], rowActivity_[iRow]));
     }
     for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          assert(eq(reducedCost[iColumn], reducedCost_[iColumn]));
     }
     delete [] columnActivity;
     delete [] rowActivity;
     delete [] reducedCost;
     delete [] dual;
#endif
     return returnCode;
}
/* Does very cursory presolve.
   rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns
*/
ClpSimplex *
ClpSimplexOther::crunch(double * rhs, int * whichRow, int * whichColumn,
                        int & nBound, bool moreBounds, bool tightenBounds)
{
     //#define CHECK_STATUS
#ifdef CHECK_STATUS
     {
          int n = 0;
          int i;
          for (i = 0; i < numberColumns_; i++)
               if (getColumnStatus(i) == ClpSimplex::basic)
                    n++;
          for (i = 0; i < numberRows_; i++)
               if (getRowStatus(i) == ClpSimplex::basic)
                    n++;
          assert (n == numberRows_);
     }
#endif

     const double * element = matrix_->getElements();
     const int * row = matrix_->getIndices();
     const CoinBigIndex * columnStart = matrix_->getVectorStarts();
     const int * columnLength = matrix_->getVectorLengths();

     CoinZeroN(rhs, numberRows_);
     int iColumn;
     int iRow;
     CoinZeroN(whichRow, numberRows_);
     int * backColumn = whichColumn + numberColumns_;
     int numberRows2 = 0;
     int numberColumns2 = 0;
     double offset = 0.0;
     const double * objective = this->objective();
     double * solution = columnActivity_;
     for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double lower = columnLower_[iColumn];
          double upper = columnUpper_[iColumn];
          if (upper > lower || getColumnStatus(iColumn) == ClpSimplex::basic) {
               backColumn[iColumn] = numberColumns2;
               whichColumn[numberColumns2++] = iColumn;
               for (CoinBigIndex j = columnStart[iColumn];
                         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    int n = whichRow[iRow];
                    if (n == 0 && element[j])
                         whichRow[iRow] = -iColumn - 1;
                    else if (n < 0)
                         whichRow[iRow] = 2;
               }
          } else {
               // fixed
               backColumn[iColumn] = -1;
               solution[iColumn] = upper;
               if (upper) {
                    offset += objective[iColumn] * upper;
                    for (CoinBigIndex j = columnStart[iColumn];
                              j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                         int iRow = row[j];
                         double value = element[j];
                         rhs[iRow] += upper * value;
                    }
               }
          }
     }
     int returnCode = 0;
     double tolerance = primalTolerance();
     nBound = 2 * numberRows_;
     for (iRow = 0; iRow < numberRows_; iRow++) {
          int n = whichRow[iRow];
          if (n > 0) {
               whichRow[numberRows2++] = iRow;
          } else if (n < 0) {
               //whichRow[numberRows2++]=iRow;
               //continue;
               // Can only do in certain circumstances as we don't know current value
               if (rowLower_[iRow] == rowUpper_[iRow] || getRowStatus(iRow) == ClpSimplex::basic) {
                    // save row and column for bound
                    whichRow[--nBound] = iRow;
                    whichRow[nBound+numberRows_] = -n - 1;
               } else if (moreBounds) {
                    // save row and column for bound
                    whichRow[--nBound] = iRow;
                    whichRow[nBound+numberRows_] = -n - 1;
               } else {
                    whichRow[numberRows2++] = iRow;
               }
          } else {
               // empty
               double rhsValue = rhs[iRow];
               if (rhsValue < rowLower_[iRow] - tolerance || rhsValue > rowUpper_[iRow] + tolerance) {
                    returnCode = 1; // infeasible
               }
          }
     }
     ClpSimplex * small = NULL;
     if (!returnCode) {
       //printf("CRUNCH from (%d,%d) to (%d,%d)\n",
       //     numberRows_,numberColumns_,numberRows2,numberColumns2);
          small = new ClpSimplex(this, numberRows2, whichRow,
                                 numberColumns2, whichColumn, true, false);
#if 0
          ClpPackedMatrix * rowCopy = dynamic_cast<ClpPackedMatrix *>(rowCopy_);
          if (rowCopy) {
               assert(!small->rowCopy());
               small->setNewRowCopy(new ClpPackedMatrix(*rowCopy, numberRows2, whichRow,
                                    numberColumns2, whichColumn));
          }
#endif
          // Set some stuff
          small->setDualBound(dualBound_);
          small->setInfeasibilityCost(infeasibilityCost_);
          small->setSpecialOptions(specialOptions_);
          small->setPerturbation(perturbation_);
          small->defaultFactorizationFrequency();
          small->setAlphaAccuracy(alphaAccuracy_);
          // If no rows left then no tightening!
          if (!numberRows2 || !numberColumns2)
               tightenBounds = false;

          int numberElements = getNumElements();
          int numberElements2 = small->getNumElements();
          small->setObjectiveOffset(objectiveOffset() - offset);
          handler_->message(CLP_CRUNCH_STATS, messages_)
                    << numberRows2 << -(numberRows_ - numberRows2)
                    << numberColumns2 << -(numberColumns_ - numberColumns2)
                    << numberElements2 << -(numberElements - numberElements2)
                    << CoinMessageEol;
          // And set objective value to match
          small->setObjectiveValue(this->objectiveValue());
          double * rowLower2 = small->rowLower();
          double * rowUpper2 = small->rowUpper();
          int jRow;
          for (jRow = 0; jRow < numberRows2; jRow++) {
               iRow = whichRow[jRow];
               if (rowLower2[jRow] > -1.0e20)
                    rowLower2[jRow] -= rhs[iRow];
               if (rowUpper2[jRow] < 1.0e20)
                    rowUpper2[jRow] -= rhs[iRow];
          }
          // and bounds
          double * columnLower2 = small->columnLower();
          double * columnUpper2 = small->columnUpper();
          const char * integerInformation = integerType_;
          for (jRow = nBound; jRow < 2 * numberRows_; jRow++) {
               iRow = whichRow[jRow];
               iColumn = whichRow[jRow+numberRows_];
               double lowerRow = rowLower_[iRow];
               if (lowerRow > -1.0e20)
                    lowerRow -= rhs[iRow];
               double upperRow = rowUpper_[iRow];
               if (upperRow < 1.0e20)
                    upperRow -= rhs[iRow];
               int jColumn = backColumn[iColumn];
               double lower = columnLower2[jColumn];
               double upper = columnUpper2[jColumn];
               double value = 0.0;
               for (CoinBigIndex j = columnStart[iColumn];
                         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    if (iRow == row[j]) {
                         value = element[j];
                         break;
                    }
               }
               assert (value);
               // convert rowLower and Upper to implied bounds on column
               double newLower = -COIN_DBL_MAX;
               double newUpper = COIN_DBL_MAX;
               if (value > 0.0) {
                    if (lowerRow > -1.0e20)
                         newLower = lowerRow / value;
                    if (upperRow < 1.0e20)
                         newUpper = upperRow / value;
               } else {
                    if (upperRow < 1.0e20)
                         newLower = upperRow / value;
                    if (lowerRow > -1.0e20)
                         newUpper = lowerRow / value;
               }
               if (integerInformation && integerInformation[iColumn]) {
                    if (newLower - floor(newLower) < 10.0 * tolerance)
                         newLower = floor(newLower);
                    else
                         newLower = ceil(newLower);
                    if (ceil(newUpper) - newUpper < 10.0 * tolerance)
                         newUpper = ceil(newUpper);
                    else
                         newUpper = floor(newUpper);
               }
               newLower = CoinMax(lower, newLower);
               newUpper = CoinMin(upper, newUpper);
               if (newLower > newUpper + tolerance) {
                    //printf("XXYY inf on bound\n");
                    returnCode = 1;
               }
               columnLower2[jColumn] = newLower;
               columnUpper2[jColumn] = CoinMax(newLower, newUpper);
               if (getRowStatus(iRow) != ClpSimplex::basic) {
                    if (getColumnStatus(iColumn) == ClpSimplex::basic) {
                         if (columnLower2[jColumn] == columnUpper2[jColumn]) {
                              // can only get here if will be fixed
                              small->setColumnStatus(jColumn, ClpSimplex::isFixed);
                         } else {
                              // solution is valid
                              if (fabs(columnActivity_[iColumn] - columnLower2[jColumn]) <
                                        fabs(columnActivity_[iColumn] - columnUpper2[jColumn]))
                                   small->setColumnStatus(jColumn, ClpSimplex::atLowerBound);
                              else
                                   small->setColumnStatus(jColumn, ClpSimplex::atUpperBound);
                         }
                    } else {
                         //printf("what now neither basic\n");
                    }
               }
          }
          if (returnCode) {
               delete small;
               small = NULL;
          } else if (tightenBounds && integerInformation) {
               // See if we can tighten any bounds
               // use rhs for upper and small duals for lower
               double * up = rhs;
               double * lo = small->dualRowSolution();
               const double * element = small->clpMatrix()->getElements();
               const int * row = small->clpMatrix()->getIndices();
               const CoinBigIndex * columnStart = small->clpMatrix()->getVectorStarts();
               //const int * columnLength = small->clpMatrix()->getVectorLengths();
               CoinZeroN(lo, numberRows2);
               CoinZeroN(up, numberRows2);
               for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
                    double upper = columnUpper2[iColumn];
                    double lower = columnLower2[iColumn];
                    //assert (columnLength[iColumn]==columnStart[iColumn+1]-columnStart[iColumn]);
                    for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn+1]; j++) {
                         int iRow = row[j];
                         double value = element[j];
                         if (value > 0.0) {
                              if (upper < 1.0e20)
                                   up[iRow] += upper * value;
                              else
                                   up[iRow] = COIN_DBL_MAX;
                              if (lower > -1.0e20)
                                   lo[iRow] += lower * value;
                              else
                                   lo[iRow] = -COIN_DBL_MAX;
                         } else {
                              if (upper < 1.0e20)
                                   lo[iRow] += upper * value;
                              else
                                   lo[iRow] = -COIN_DBL_MAX;
                              if (lower > -1.0e20)
                                   up[iRow] += lower * value;
                              else
                                   up[iRow] = COIN_DBL_MAX;
                         }
                    }
               }
               double * rowLower2 = small->rowLower();
               double * rowUpper2 = small->rowUpper();
               bool feasible = true;
               // make safer
               for (int iRow = 0; iRow < numberRows2; iRow++) {
                    double lower = lo[iRow];
                    if (lower > rowUpper2[iRow] + tolerance) {
                         feasible = false;
                         break;
                    } else {
                         lo[iRow] = CoinMin(lower - rowUpper2[iRow], 0.0) - tolerance;
                    }
                    double upper = up[iRow];
                    if (upper < rowLower2[iRow] - tolerance) {
                         feasible = false;
                         break;
                    } else {
                         up[iRow] = CoinMax(upper - rowLower2[iRow], 0.0) + tolerance;
                    }
               }
               if (!feasible) {
                    delete small;
                    small = NULL;
               } else {
                    // and tighten
                    for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
                         if (integerInformation[whichColumn[iColumn]]) {
                              double upper = columnUpper2[iColumn];
                              double lower = columnLower2[iColumn];
                              double newUpper = upper;
                              double newLower = lower;
                              double difference = upper - lower;
                              if (lower > -1000.0 && upper < 1000.0) {
                                   for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn+1]; j++) {
                                        int iRow = row[j];
                                        double value = element[j];
                                        if (value > 0.0) {
                                             double upWithOut = up[iRow] - value * difference;
                                             if (upWithOut < 0.0) {
                                                  newLower = CoinMax(newLower, lower - (upWithOut + tolerance) / value);
                                             }
                                             double lowWithOut = lo[iRow] + value * difference;
                                             if (lowWithOut > 0.0) {
                                                  newUpper = CoinMin(newUpper, upper - (lowWithOut - tolerance) / value);
                                             }
                                        } else {
                                             double upWithOut = up[iRow] + value * difference;
                                             if (upWithOut < 0.0) {
                                                  newUpper = CoinMin(newUpper, upper - (upWithOut + tolerance) / value);
                                             }
                                             double lowWithOut = lo[iRow] - value * difference;
                                             if (lowWithOut > 0.0) {
                                                  newLower = CoinMax(newLower, lower - (lowWithOut - tolerance) / value);
                                             }
                                        }
                                   }
                                   if (newLower > lower || newUpper < upper) {
                                        if (fabs(newUpper - floor(newUpper + 0.5)) > 1.0e-6)
                                             newUpper = floor(newUpper);
                                        else
                                             newUpper = floor(newUpper + 0.5);
                                        if (fabs(newLower - ceil(newLower - 0.5)) > 1.0e-6)
                                             newLower = ceil(newLower);
                                        else
                                             newLower = ceil(newLower - 0.5);
                                        // change may be too small - check
                                        if (newLower > lower || newUpper < upper) {
                                             if (newUpper >= newLower) {
                                                  // Could also tighten in this
                                                  //printf("%d bounds %g %g tightened to %g %g\n",
                                                  //     iColumn,columnLower2[iColumn],columnUpper2[iColumn],
                                                  //     newLower,newUpper);
#if 1
                                                  columnUpper2[iColumn] = newUpper;
                                                  columnLower2[iColumn] = newLower;
                                                  columnUpper_[whichColumn[iColumn]] = newUpper;
                                                  columnLower_[whichColumn[iColumn]] = newLower;
#endif
                                                  // and adjust bounds on rows
                                                  newUpper -= upper;
                                                  newLower -= lower;
                                                  for (CoinBigIndex j = columnStart[iColumn]; j < columnStart[iColumn+1]; j++) {
                                                       int iRow = row[j];
                                                       double value = element[j];
                                                       if (value > 0.0) {
                                                            up[iRow] += newUpper * value;
                                                            lo[iRow] += newLower * value;
                                                       } else {
                                                            lo[iRow] += newUpper * value;
                                                            up[iRow] += newLower * value;
                                                       }
                                                  }
                                             } else {
                                                  // infeasible
                                                  //printf("%d bounds infeasible %g %g tightened to %g %g\n",
                                                  //     iColumn,columnLower2[iColumn],columnUpper2[iColumn],
                                                  //     newLower,newUpper);
#if 1
                                                  delete small;
                                                  small = NULL;
                                                  break;
#endif
                                             }
                                        }
                                   }
                              }
                         }
                    }
               }
          }
     }
#if 0
     if (small) {
          static int which = 0;
          which++;
          char xxxx[20];
          sprintf(xxxx, "bad%d.mps", which);
          small->writeMps(xxxx, 0, 1);
          sprintf(xxxx, "largebad%d.mps", which);
          writeMps(xxxx, 0, 1);
          printf("bad%d %x old size %d %d new %d %d\n", which, small,
                 numberRows_, numberColumns_, small->numberRows(), small->numberColumns());
#if 0
          for (int i = 0; i < numberColumns_; i++)
               printf("Bound %d %g %g\n", i, columnLower_[i], columnUpper_[i]);
          for (int i = 0; i < numberRows_; i++)
               printf("Row bound %d %g %g\n", i, rowLower_[i], rowUpper_[i]);
#endif
     }
#endif
#ifdef CHECK_STATUS
     {
          int n = 0;
          int i;
          for (i = 0; i < small->numberColumns(); i++)
               if (small->getColumnStatus(i) == ClpSimplex::basic)
                    n++;
          for (i = 0; i < small->numberRows(); i++)
               if (small->getRowStatus(i) == ClpSimplex::basic)
                    n++;
          assert (n == small->numberRows());
     }
#endif
     return small;
}
/* After very cursory presolve.
   rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns.
*/
void
ClpSimplexOther::afterCrunch(const ClpSimplex & small,
                             const int * whichRow,
                             const int * whichColumn, int nBound)
{
#ifndef NDEBUG
     for (int i = 0; i < small.numberRows(); i++)
          assert (whichRow[i] >= 0 && whichRow[i] < numberRows_);
     for (int i = 0; i < small.numberColumns(); i++)
          assert (whichColumn[i] >= 0 && whichColumn[i] < numberColumns_);
#endif
     getbackSolution(small, whichRow, whichColumn);
     // and deal with status for bounds
     const double * element = matrix_->getElements();
     const int * row = matrix_->getIndices();
     const CoinBigIndex * columnStart = matrix_->getVectorStarts();
     const int * columnLength = matrix_->getVectorLengths();
     double tolerance = primalTolerance();
     double djTolerance = dualTolerance();
     for (int jRow = nBound; jRow < 2 * numberRows_; jRow++) {
          int iRow = whichRow[jRow];
          int iColumn = whichRow[jRow+numberRows_];
          if (getColumnStatus(iColumn) != ClpSimplex::basic) {
               double lower = columnLower_[iColumn];
               double upper = columnUpper_[iColumn];
               double value = columnActivity_[iColumn];
               double djValue = reducedCost_[iColumn];
               dual_[iRow] = 0.0;
               if (upper > lower) {
                    if (value < lower + tolerance && djValue > -djTolerance) {
                         setColumnStatus(iColumn, ClpSimplex::atLowerBound);
                         setRowStatus(iRow, ClpSimplex::basic);
                    } else if (value > upper - tolerance && djValue < djTolerance) {
                         setColumnStatus(iColumn, ClpSimplex::atUpperBound);
                         setRowStatus(iRow, ClpSimplex::basic);
                    } else {
                         // has to be basic
                         setColumnStatus(iColumn, ClpSimplex::basic);
                         reducedCost_[iColumn] = 0.0;
                         double value = 0.0;
                         for (CoinBigIndex j = columnStart[iColumn];
                                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                              if (iRow == row[j]) {
                                   value = element[j];
                                   break;
                              }
                         }
                         dual_[iRow] = djValue / value;
                         if (rowUpper_[iRow] > rowLower_[iRow]) {
                              if (fabs(rowActivity_[iRow] - rowLower_[iRow]) <
                                        fabs(rowActivity_[iRow] - rowUpper_[iRow]))
                                   setRowStatus(iRow, ClpSimplex::atLowerBound);
                              else
                                   setRowStatus(iRow, ClpSimplex::atUpperBound);
                         } else {
                              setRowStatus(iRow, ClpSimplex::isFixed);
                         }
                    }
               } else {
                    // row can always be basic
                    setRowStatus(iRow, ClpSimplex::basic);
               }
          } else {
               // row can always be basic
               setRowStatus(iRow, ClpSimplex::basic);
          }
     }
     //#ifndef NDEBUG
#if 0
     if  (small.status() == 0) {
          int n = 0;
          int i;
          for (i = 0; i < numberColumns; i++)
               if (getColumnStatus(i) == ClpSimplex::basic)
                    n++;
          for (i = 0; i < numberRows; i++)
               if (getRowStatus(i) == ClpSimplex::basic)
                    n++;
          assert (n == numberRows);
     }
#endif
}
/* Tightens integer bounds - returns number tightened or -1 if infeasible
 */
int
ClpSimplexOther::tightenIntegerBounds(double * rhsSpace)
{
     // See if we can tighten any bounds
     // use rhs for upper and small duals for lower
     double * up = rhsSpace;
     double * lo = dual_;
     const double * element = matrix_->getElements();
     const int * row = matrix_->getIndices();
     const CoinBigIndex * columnStart = matrix_->getVectorStarts();
     const int * columnLength = matrix_->getVectorLengths();
     CoinZeroN(lo, numberRows_);
     CoinZeroN(up, numberRows_);
     for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double upper = columnUpper_[iColumn];
          double lower = columnLower_[iColumn];
          //assert (columnLength[iColumn]==columnStart[iColumn+1]-columnStart[iColumn]);
          for (CoinBigIndex j = columnStart[iColumn];
                    j < columnStart[iColumn] + columnLength[iColumn]; j++) {
               int iRow = row[j];
               double value = element[j];
               if (value > 0.0) {
                    if (upper < 1.0e20)
                         up[iRow] += upper * value;
                    else
                         up[iRow] = COIN_DBL_MAX;
                    if (lower > -1.0e20)
                         lo[iRow] += lower * value;
                    else
                         lo[iRow] = -COIN_DBL_MAX;
               } else {
                    if (upper < 1.0e20)
                         lo[iRow] += upper * value;
                    else
                         lo[iRow] = -COIN_DBL_MAX;
                    if (lower > -1.0e20)
                         up[iRow] += lower * value;
                    else
                         up[iRow] = COIN_DBL_MAX;
               }
          }
     }
     bool feasible = true;
     // make safer
     double tolerance = primalTolerance();
     for (int iRow = 0; iRow < numberRows_; iRow++) {
          double lower = lo[iRow];
          if (lower > rowUpper_[iRow] + tolerance) {
               feasible = false;
               break;
          } else {
               lo[iRow] = CoinMin(lower - rowUpper_[iRow], 0.0) - tolerance;
          }
          double upper = up[iRow];
          if (upper < rowLower_[iRow] - tolerance) {
               feasible = false;
               break;
          } else {
               up[iRow] = CoinMax(upper - rowLower_[iRow], 0.0) + tolerance;
          }
     }
     int numberTightened = 0;
     if (!feasible) {
          return -1;
     } else if (integerType_) {
          // and tighten
          for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
               if (integerType_[iColumn]) {
                    double upper = columnUpper_[iColumn];
                    double lower = columnLower_[iColumn];
                    double newUpper = upper;
                    double newLower = lower;
                    double difference = upper - lower;
                    if (lower > -1000.0 && upper < 1000.0) {
                         for (CoinBigIndex j = columnStart[iColumn];
                                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                              int iRow = row[j];
                              double value = element[j];
                              if (value > 0.0) {
                                   double upWithOut = up[iRow] - value * difference;
                                   if (upWithOut < 0.0) {
                                        newLower = CoinMax(newLower, lower - (upWithOut + tolerance) / value);
                                   }
                                   double lowWithOut = lo[iRow] + value * difference;
                                   if (lowWithOut > 0.0) {
                                        newUpper = CoinMin(newUpper, upper - (lowWithOut - tolerance) / value);
                                   }
                              } else {
                                   double upWithOut = up[iRow] + value * difference;
                                   if (upWithOut < 0.0) {
                                        newUpper = CoinMin(newUpper, upper - (upWithOut + tolerance) / value);
                                   }
                                   double lowWithOut = lo[iRow] - value * difference;
                                   if (lowWithOut > 0.0) {
                                        newLower = CoinMax(newLower, lower - (lowWithOut - tolerance) / value);
                                   }
                              }
                         }
                         if (newLower > lower || newUpper < upper) {
                              if (fabs(newUpper - floor(newUpper + 0.5)) > 1.0e-6)
                                   newUpper = floor(newUpper);
                              else
                                   newUpper = floor(newUpper + 0.5);
                              if (fabs(newLower - ceil(newLower - 0.5)) > 1.0e-6)
                                   newLower = ceil(newLower);
                              else
                                   newLower = ceil(newLower - 0.5);
                              // change may be too small - check
                              if (newLower > lower || newUpper < upper) {
                                   if (newUpper >= newLower) {
                                        numberTightened++;
                                        //printf("%d bounds %g %g tightened to %g %g\n",
                                        //     iColumn,columnLower_[iColumn],columnUpper_[iColumn],
                                        //     newLower,newUpper);
                                        columnUpper_[iColumn] = newUpper;
                                        columnLower_[iColumn] = newLower;
                                        // and adjust bounds on rows
                                        newUpper -= upper;
                                        newLower -= lower;
                                        for (CoinBigIndex j = columnStart[iColumn];
                                                  j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                                             int iRow = row[j];
                                             double value = element[j];
                                             if (value > 0.0) {
                                                  up[iRow] += newUpper * value;
                                                  lo[iRow] += newLower * value;
                                             } else {
                                                  lo[iRow] += newUpper * value;
                                                  up[iRow] += newLower * value;
                                             }
                                        }
                                   } else {
                                        // infeasible
                                        //printf("%d bounds infeasible %g %g tightened to %g %g\n",
                                        //     iColumn,columnLower_[iColumn],columnUpper_[iColumn],
                                        //     newLower,newUpper);
                                        return -1;
                                   }
                              }
                         }
                    }
               }
          }
     }
     return numberTightened;
}
/* Parametrics
   This is an initial slow version.
   The code uses current bounds + theta * change (if change array not NULL)
   and similarly for objective.
   It starts at startingTheta and returns ending theta in endingTheta.
   If reportIncrement 0.0 it will report on any movement
   If reportIncrement >0.0 it will report at startingTheta+k*reportIncrement.
   If it can not reach input endingTheta return code will be 1 for infeasible,
   2 for unbounded, if error on ranges -1,  otherwise 0.
   Normal report is just theta and objective but
   if event handler exists it may do more
   On exit endingTheta is maximum reached (can be used for next startingTheta)
*/
int
ClpSimplexOther::parametrics(double startingTheta, double & endingTheta, double reportIncrement,
                             const double * lowerChangeBound, const double * upperChangeBound,
                             const double * lowerChangeRhs, const double * upperChangeRhs,
                             const double * changeObjective)
{
     bool needToDoSomething = true;
     bool canTryQuick = (reportIncrement) ? true : false;
     // Save copy of model
     ClpSimplex copyModel = *this;
     int savePerturbation = perturbation_;
     perturbation_ = 102; // switch off
     while (needToDoSomething) {
          needToDoSomething = false;
          algorithm_ = -1;

          // save data
          ClpDataSave data = saveData();
	  // Dantzig 
	  ClpDualRowPivot * savePivot = dualRowPivot_;
	  dualRowPivot_ = new ClpDualRowDantzig();
	  dualRowPivot_->setModel(this);
          int returnCode = reinterpret_cast<ClpSimplexDual *> (this)->startupSolve(0, NULL, 0);
          int iRow, iColumn;
          double * chgUpper = NULL;
          double * chgLower = NULL;
          double * chgObjective = NULL;


          if (!returnCode) {
               // Find theta when bounds will cross over and create arrays
               int numberTotal = numberRows_ + numberColumns_;
               chgLower = new double[numberTotal];
               memset(chgLower, 0, numberTotal * sizeof(double));
               chgUpper = new double[numberTotal];
               memset(chgUpper, 0, numberTotal * sizeof(double));
               chgObjective = new double[numberTotal];
               memset(chgObjective, 0, numberTotal * sizeof(double));
               assert (!rowScale_);
               double maxTheta = 1.0e50;
               if (lowerChangeRhs || upperChangeRhs) {
                    for (iRow = 0; iRow < numberRows_; iRow++) {
                         double lower = rowLower_[iRow];
                         double upper = rowUpper_[iRow];
                         if (lower > upper) {
                              maxTheta = -1.0;
                              break;
                         }
                         double lowerChange = (lowerChangeRhs) ? lowerChangeRhs[iRow] : 0.0;
                         double upperChange = (upperChangeRhs) ? upperChangeRhs[iRow] : 0.0;
                         if (lower > -1.0e20 && upper < 1.0e20) {
                              if (lower + maxTheta * lowerChange > upper + maxTheta * upperChange) {
                                   maxTheta = (upper - lower) / (lowerChange - upperChange);
                              }
                         }
                         if (lower > -1.0e20) {
                              lower_[numberColumns_+iRow] += startingTheta * lowerChange;
                              chgLower[numberColumns_+iRow] = lowerChange;
                         }
                         if (upper < 1.0e20) {
                              upper_[numberColumns_+iRow] += startingTheta * upperChange;
                              chgUpper[numberColumns_+iRow] = upperChange;
                         }
                    }
               }
               if (maxTheta > 0.0) {
                    if (lowerChangeBound || upperChangeBound) {
                         for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
                              double lower = columnLower_[iColumn];
                              double upper = columnUpper_[iColumn];
                              if (lower > upper) {
                                   maxTheta = -1.0;
                                   break;
                              }
                              double lowerChange = (lowerChangeBound) ? lowerChangeBound[iColumn] : 0.0;
                              double upperChange = (upperChangeBound) ? upperChangeBound[iColumn] : 0.0;
                              if (lower > -1.0e20 && upper < 1.0e20) {
                                   if (lower + maxTheta * lowerChange > upper + maxTheta * upperChange) {
                                        maxTheta = (upper - lower) / (lowerChange - upperChange);
                                   }
                              }
                              if (lower > -1.0e20) {
                                   lower_[iColumn] += startingTheta * lowerChange;
                                   chgLower[iColumn] = lowerChange;
                              }
                              if (upper < 1.0e20) {
                                   upper_[iColumn] += startingTheta * upperChange;
                                   chgUpper[iColumn] = upperChange;
                              }
                         }
                    }
                    if (maxTheta == 1.0e50)
                         maxTheta = COIN_DBL_MAX;
               }
               if (maxTheta < 0.0) {
                    // bad ranges or initial
                    returnCode = -1;
               }
	       if (maxTheta < endingTheta) {
		 char line[100];
		 sprintf(line,"Crossover considerations reduce ending  theta from %g to %g\n", 
			 endingTheta,maxTheta);
		 handler_->message(CLP_GENERAL,messages_)
		   << line << CoinMessageEol;
		 endingTheta = maxTheta;
	       }
               if (endingTheta < startingTheta) {
                    // bad initial
                    returnCode = -2;
               }
          }
          double saveEndingTheta = endingTheta;
          if (!returnCode) {
               if (changeObjective) {
                    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
                         chgObjective[iColumn] = changeObjective[iColumn];
                         cost_[iColumn] += startingTheta * changeObjective[iColumn];
                    }
               }
               double * saveDuals = NULL;
               reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(0, saveDuals, -1, data);
               assert (!problemStatus_);
	       for (int i=0;i<numberRows_+numberColumns_;i++)
		 setFakeBound(i, noFake);
               // Now do parametrics
	       handler_->message(CLP_PARAMETRICS_STATS, messages_)
		 << startingTheta << objectiveValue() << CoinMessageEol;
               while (!returnCode) {
		    //assert (reportIncrement);
                    returnCode = parametricsLoop(startingTheta, endingTheta, reportIncrement,
                                                 chgLower, chgUpper, chgObjective, data,
                                                 canTryQuick);
                    if (!returnCode) {
                         //double change = endingTheta-startingTheta;
                         startingTheta = endingTheta;
                         endingTheta = saveEndingTheta;
                         //for (int i=0;i<numberTotal;i++) {
                         //lower_[i] += change*chgLower[i];
                         //upper_[i] += change*chgUpper[i];
                         //cost_[i] += change*chgObjective[i];
                         //}
			 handler_->message(CLP_PARAMETRICS_STATS, messages_)
			   << startingTheta << objectiveValue() << CoinMessageEol;
                         if (startingTheta >= endingTheta)
                              break;
                    } else if (returnCode == -1) {
                         // trouble - do external solve
                         needToDoSomething = true;
		    } else if (problemStatus_==1) {
		      // can't move any further
		      if (!canTryQuick) {
			handler_->message(CLP_PARAMETRICS_STATS, messages_)
			  << endingTheta << objectiveValue() << CoinMessageEol;
			problemStatus_=0;
		      }
                    } else {
                         abort();
                    }
               }
          }
          reinterpret_cast<ClpSimplexDual *> (this)->finishSolve(0);

          delete dualRowPivot_;
          dualRowPivot_ = savePivot;
          // Restore any saved stuff
          restoreData(data);
          if (needToDoSomething) {
               double saveStartingTheta = startingTheta; // known to be feasible
               int cleanedUp = 1;
               while (cleanedUp) {
                    // tweak
                    if (cleanedUp == 1) {
                         if (!reportIncrement)
                              startingTheta = CoinMin(startingTheta + 1.0e-5, saveEndingTheta);
                         else
                              startingTheta = CoinMin(startingTheta + reportIncrement, saveEndingTheta);
                    } else {
                         // restoring to go slowly
                         startingTheta = saveStartingTheta;
                    }
                    // only works if not scaled
                    int i;
                    const double * obj1 = objective();
                    double * obj2 = copyModel.objective();
                    const double * lower1 = columnLower_;
                    double * lower2 = copyModel.columnLower();
                    const double * upper1 = columnUpper_;
                    double * upper2 = copyModel.columnUpper();
                    for (i = 0; i < numberColumns_; i++) {
                         obj2[i] = obj1[i] + startingTheta * chgObjective[i];
                         lower2[i] = lower1[i] + startingTheta * chgLower[i];
                         upper2[i] = upper1[i] + startingTheta * chgUpper[i];
                    }
                    lower1 = rowLower_;
                    lower2 = copyModel.rowLower();
                    upper1 = rowUpper_;
                    upper2 = copyModel.rowUpper();
                    for (i = 0; i < numberRows_; i++) {
                         lower2[i] = lower1[i] + startingTheta * chgLower[i+numberColumns_];
                         upper2[i] = upper1[i] + startingTheta * chgUpper[i+numberColumns_];
                    }
                    copyModel.dual();
                    if (copyModel.problemStatus()) {
		      char line[100];
		      sprintf(line,"Can not get to theta of %g\n", startingTheta);
		      handler_->message(CLP_GENERAL,messages_)
			<< line << CoinMessageEol;
                         canTryQuick = false; // do slowly to get exact amount
                         // back to last known good
                         if (cleanedUp == 1)
                              cleanedUp = 2;
                         else
                              abort();
                    } else {
                         // and move stuff back
                         int numberTotal = numberRows_ + numberColumns_;
                         CoinMemcpyN(copyModel.statusArray(), numberTotal, status_);
                         CoinMemcpyN(copyModel.primalColumnSolution(), numberColumns_, columnActivity_);
                         CoinMemcpyN(copyModel.primalRowSolution(), numberRows_, rowActivity_);
                         cleanedUp = 0;
                    }
               }
          }
          delete [] chgLower;
          delete [] chgUpper;
          delete [] chgObjective;
     }
     perturbation_ = savePerturbation;
     char line[100];
     sprintf(line,"Ending theta %g\n", endingTheta);
     handler_->message(CLP_GENERAL,messages_)
       << line << CoinMessageEol;
     return problemStatus_;
}
/* Version of parametrics which reads from file
   See CbcClpParam.cpp for details of format
   Returns -2 if unable to open file */
int 
ClpSimplexOther::parametrics(const char * dataFile)
{
  int returnCode=-2;
  FILE *fp = fopen(dataFile, "r");
  char line[200];
  if (!fp) {
    handler_->message(CLP_UNABLE_OPEN, messages_)
      << dataFile << CoinMessageEol;
    return -2;
  }

  if (!fgets(line, 200, fp)) {
    sprintf(line,"Empty parametrics file %s?",dataFile);
    handler_->message(CLP_GENERAL,messages_)
      << line << CoinMessageEol;
    fclose(fp);
    return -2;
  }
  char * pos = line;
  char * put = line;
  while (*pos >= ' ' && *pos != '\n') {
    if (*pos != ' ' && *pos != '\t') {
      *put = static_cast<char>(tolower(*pos));
      put++;
    }
    pos++;
  }
  *put = '\0';
  pos = line;
  double startTheta=0.0;
  double endTheta=0.0;
  double intervalTheta=COIN_DBL_MAX;
  int detail=0;
  bool good = true;
  while (good) {
    good=false;
    // check ROWS
    char * comma = strchr(pos, ',');
    if (!comma)
      break;
    *comma = '\0';
    if (strcmp(pos,"rows"))
      break;
    *comma = ',';
    pos = comma+1;
    // check lower theta
    comma = strchr(pos, ',');
    if (!comma)
      break;
    *comma = '\0';
    startTheta = atof(pos);
    *comma = ',';
    pos = comma+1;
    // check upper theta
    comma = strchr(pos, ',');
    good=true;
    if (comma)
      *comma = '\0';
    endTheta = atof(pos);
    if (comma) {
      *comma = ',';
      pos = comma+1;
      comma = strchr(pos, ',');
      if (comma)
	*comma = '\0';
      intervalTheta = atof(pos);
      if (comma) {
	*comma = ',';
	pos = comma+1;
	comma = strchr(pos, ',');
	if (comma)
	  *comma = '\0';
	detail = atoi(pos);
	if (comma) 
	*comma = ',';
      }
    }
    break;
  }
  if (good) {
    if (startTheta<0.0||
	startTheta>endTheta||
	intervalTheta<0.0)
      good=false;
    if (detail<0||detail>1)
      good=false;
  }
  if (intervalTheta>=endTheta)
    intervalTheta=0.0;
  if (!good) {
    sprintf(line,"Odd first line %s on file %s?",line,dataFile);
    handler_->message(CLP_GENERAL,messages_)
      << line << CoinMessageEol;
    fclose(fp);
    return -2;
  }
  if (!fgets(line, 200, fp)) {
    sprintf(line,"Not enough records on parametrics file %s?",dataFile);
    handler_->message(CLP_GENERAL,messages_)
      << line << CoinMessageEol;
    fclose(fp);
    return -2;
  }
  double * lowerRowMove = NULL;
  double * upperRowMove = NULL;
  double * lowerColumnMove = NULL;
  double * upperColumnMove = NULL;
  double * objectiveMove = NULL;
  char saveLine[200];
  saveLine[0]='\0';
  std::string headingsRow[] = {"name", "number", "lower", "upper", "rhs"};
  int gotRow[] = { -1, -1, -1, -1, -1};
  int orderRow[5];
  assert(sizeof(gotRow) == sizeof(orderRow));
  int nAcross = 0;
  pos = line;
  put = line;
  while (*pos >= ' ' && *pos != '\n') {
    if (*pos != ' ' && *pos != '\t') {
      *put = static_cast<char>(tolower(*pos));
      put++;
    }
    pos++;
  }
  *put = '\0';
  pos = line;
  int i;
  good = true;
  if (strncmp(line,"column",6)) {
    while (pos) {
      char * comma = strchr(pos, ',');
      if (comma)
	*comma = '\0';
      for (i = 0; i < static_cast<int> (sizeof(gotRow) / sizeof(int)); i++) {
	if (headingsRow[i] == pos) {
	  if (gotRow[i] < 0) {
	    orderRow[nAcross] = i;
	    gotRow[i] = nAcross++;
	  } else {
	    // duplicate
	    good = false;
	  }
	  break;
	}
      }
      if (i == static_cast<int> (sizeof(gotRow) / sizeof(int)))
	good = false;
      if (comma) {
	*comma = ',';
	pos = comma + 1;
      } else {
	break;
      }
    }
    if (gotRow[0] < 0 && gotRow[1] < 0)
      good = false;
    if (gotRow[0] >= 0 && gotRow[1] >= 0)
      good = false;
    if (gotRow[0] >= 0 && !lengthNames())
      good = false;
    if (gotRow[4]<0) {
      if (gotRow[2] < 0 && gotRow[3] >= 0)
	good = false;
      else if (gotRow[3] < 0 && gotRow[2] >= 0)
	good = false;
    } else if (gotRow[2]>=0||gotRow[3]>=0) {
      good = false;
    }
    if (good) {
      char ** rowNames = new char * [numberRows_];
      int iRow;
      for (iRow = 0; iRow < numberRows_; iRow++) {
	rowNames[iRow] =
	  CoinStrdup(rowName(iRow).c_str());
      }
      lowerRowMove = new double [numberRows_];
      memset(lowerRowMove,0,numberRows_*sizeof(double));
      upperRowMove = new double [numberRows_];
      memset(upperRowMove,0,numberRows_*sizeof(double));
      int nLine = 0;
      int nBadLine = 0;
      int nBadName = 0;
      bool goodLine=false;
      while (fgets(line, 200, fp)) {
	goodLine=true;
	if (!strncmp(line, "ENDATA", 6)||
	    !strncmp(line, "COLUMN",6))
	  break;
	goodLine=false;
	nLine++;
	iRow = -1;
	double upper = 0.0;
	double lower = 0.0;
	char * pos = line;
	char * put = line;
	while (*pos >= ' ' && *pos != '\n') {
	  if (*pos != ' ' && *pos != '\t') {
	    *put = *pos;
	    put++;
	  }
	  pos++;
	}
	*put = '\0';
	pos = line;
	for (int i = 0; i < nAcross; i++) {
	  char * comma = strchr(pos, ',');
	  if (comma) {
	    *comma = '\0';
	  } else if (i < nAcross - 1) {
	    nBadLine++;
	    break;
	  }
	  switch (orderRow[i]) {
	    // name
	  case 0:
	    // For large problems this could be slow
	    for (iRow = 0; iRow < numberRows_; iRow++) {
	      if (!strcmp(rowNames[iRow], pos))
		break;
	    }
	    if (iRow == numberRows_)
	      iRow = -1;
	    break;
	    // number
	  case 1:
	    iRow = atoi(pos);
	    if (iRow < 0 || iRow >= numberRows_)
	      iRow = -1;
	    break;
	    // lower
	  case 2:
	    upper = atof(pos);
	    break;
	    // upper
	  case 3:
	    lower = atof(pos);
	    break;
	    // rhs
	  case 4:
	    lower = atof(pos);
	    upper = lower;
	    break;
	  }
	  if (comma) {
	    *comma = ',';
	    pos = comma + 1;
	  }
	}
	if (iRow >= 0) {
	  if (rowLower_[iRow]>-1.0e20)
	    lowerRowMove[iRow] = lower;
	  else
	    lowerRowMove[iRow]=0.0;
	  if (rowUpper_[iRow]<1.0e20)
	    upperRowMove[iRow] = upper;
	  else
	    upperRowMove[iRow] = lower;
	} else {
	  nBadName++;
	  if(saveLine[0]=='\0')
	    strcpy(saveLine,line);
	}
      }
      sprintf(line,"%d Row fields and %d records", nAcross, nLine);
      handler_->message(CLP_GENERAL,messages_)
	<< line << CoinMessageEol;
      if (nBadName) {
	sprintf(line," ** %d records did not match on name/sequence, first bad %s", nBadName,saveLine);
	handler_->message(CLP_GENERAL,messages_)
	  << line << CoinMessageEol;
	returnCode=-1;
	good=false;
      }
      for (iRow = 0; iRow < numberRows_; iRow++) {
	free(rowNames[iRow]);
      }
      delete [] rowNames;
    } else {
      sprintf(line,"Duplicate or unknown keyword - or name/number fields wrong");
      handler_->message(CLP_GENERAL,messages_)
	<< line << CoinMessageEol;
      returnCode=-1;
      good=false;
    }
  }
  if (good&&(!strncmp(line, "COLUMN",6)||!strncmp(line, "column",6))) {
    if (!fgets(line, 200, fp)) {
      sprintf(line,"Not enough records on parametrics file %s after COLUMNS?",dataFile);
      handler_->message(CLP_GENERAL,messages_)
	<< line << CoinMessageEol;
      fclose(fp);
      return -2;
    }
    std::string headingsColumn[] = {"name", "number", "lower", "upper", "objective"};
    saveLine[0]='\0';
    int gotColumn[] = { -1, -1, -1, -1, -1};
    int orderColumn[5];
    assert(sizeof(gotColumn) == sizeof(orderColumn));
    nAcross = 0;
    pos = line;
    put = line;
    while (*pos >= ' ' && *pos != '\n') {
      if (*pos != ' ' && *pos != '\t') {
	*put = static_cast<char>(tolower(*pos));
	put++;
      }
      pos++;
    }
    *put = '\0';
    pos = line;
    int i;
    if (strncmp(line,"endata",6)&&good) {
      while (pos) {
	char * comma = strchr(pos, ',');
	if (comma)
	  *comma = '\0';
	for (i = 0; i < static_cast<int> (sizeof(gotColumn) / sizeof(int)); i++) {
	  if (headingsColumn[i] == pos) {
	    if (gotColumn[i] < 0) {
	      orderColumn[nAcross] = i;
	      gotColumn[i] = nAcross++;
	    } else {
	      // duplicate
	      good = false;
	    }
	    break;
	  }
	}
	if (i == static_cast<int> (sizeof(gotColumn) / sizeof(int)))
	  good = false;
	if (comma) {
	  *comma = ',';
	  pos = comma + 1;
	} else {
	  break;
	}
      }
      if (gotColumn[0] < 0 && gotColumn[1] < 0)
	good = false;
      if (gotColumn[0] >= 0 && gotColumn[1] >= 0)
	good = false;
      if (gotColumn[0] >= 0 && !lengthNames())
	good = false;
      if (good) {
	char ** columnNames = new char * [numberColumns_];
	int iColumn;
	for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
	  columnNames[iColumn] =
	    CoinStrdup(columnName(iColumn).c_str());
	}
	lowerColumnMove = reinterpret_cast<double *> (malloc(numberColumns_ * sizeof(double)));
	memset(lowerColumnMove,0,numberColumns_*sizeof(double));
	upperColumnMove = reinterpret_cast<double *> (malloc(numberColumns_ * sizeof(double)));
	memset(upperColumnMove,0,numberColumns_*sizeof(double));
	objectiveMove = reinterpret_cast<double *> (malloc(numberColumns_ * sizeof(double)));
	memset(objectiveMove,0,numberColumns_*sizeof(double));
	int nLine = 0;
	int nBadLine = 0;
	int nBadName = 0;
	bool goodLine=false;
	while (fgets(line, 200, fp)) {
	  goodLine=true;
	  if (!strncmp(line, "ENDATA", 6))
	    break;
	  goodLine=false;
	  nLine++;
	  iColumn = -1;
	  double upper = 0.0;
	  double lower = 0.0;
	  double obj =0.0;
	  char * pos = line;
	  char * put = line;
	  while (*pos >= ' ' && *pos != '\n') {
	    if (*pos != ' ' && *pos != '\t') {
	      *put = *pos;
	      put++;
	    }
	    pos++;
	  }
	  *put = '\0';
	  pos = line;
	  for (int i = 0; i < nAcross; i++) {
	    char * comma = strchr(pos, ',');
	    if (comma) {
	      *comma = '\0';
	    } else if (i < nAcross - 1) {
	      nBadLine++;
	      break;
	    }
	    switch (orderColumn[i]) {
	      // name
	    case 0:
	      // For large problems this could be slow
	      for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
		if (!strcmp(columnNames[iColumn], pos))
		  break;
	      }
	      if (iColumn == numberColumns_)
		iColumn = -1;
	      break;
	      // number
	    case 1:
	      iColumn = atoi(pos);
	      if (iColumn < 0 || iColumn >= numberColumns_)
		iColumn = -1;
	      break;
	      // lower
	    case 2:
	      upper = atof(pos);
	      break;
	      // upper
	    case 3:
	      lower = atof(pos);
	      break;
	      // objective
	    case 4:
	      obj = atof(pos);
	      upper = lower;
	      break;
	    }
	    if (comma) {
	      *comma = ',';
	      pos = comma + 1;
	    }
	  }
	  if (iColumn >= 0) {
	    if (columnLower_[iColumn]>-1.0e20)
	      lowerColumnMove[iColumn] = lower;
	    else
	      lowerColumnMove[iColumn]=0.0;
	    if (columnUpper_[iColumn]<1.0e20)
	      upperColumnMove[iColumn] = upper;
	    else
	      upperColumnMove[iColumn] = lower;
	    objectiveMove[iColumn] = obj;
	  } else {
	    nBadName++;
	    if(saveLine[0]=='\0')
	      strcpy(saveLine,line);
	  }
	}
	sprintf(line,"%d Column fields and %d records", nAcross, nLine);
	handler_->message(CLP_GENERAL,messages_)
	  << line << CoinMessageEol;
	if (nBadName) {
	  sprintf(line," ** %d records did not match on name/sequence, first bad %s", nBadName,saveLine);
	  handler_->message(CLP_GENERAL,messages_)
	    << line << CoinMessageEol;
	  returnCode=-1;
	  good=false;
	}
	for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
	  free(columnNames[iColumn]);
	}
	delete [] columnNames;
      } else {
	sprintf(line,"Duplicate or unknown keyword - or name/number fields wrong");
	handler_->message(CLP_GENERAL,messages_)
	  << line << CoinMessageEol;
	returnCode=-1;
	good=false;
      }
    }
  }
  returnCode=-1;
  if (good) {
    // clean arrays
    if (lowerRowMove) {
      bool empty=true;
      for (int i=0;i<numberRows_;i++) {
	if (lowerRowMove[i]) {
	  empty=false;
	break;
	}
      }
      if (empty) {
	delete [] lowerRowMove;
	lowerRowMove=NULL;
      }
    }
    if (upperRowMove) {
      bool empty=true;
      for (int i=0;i<numberRows_;i++) {
	if (upperRowMove[i]) {
	  empty=false;
	break;
	}
      }
      if (empty) {
	delete [] upperRowMove;
	upperRowMove=NULL;
      }
    }
    if (lowerColumnMove) {
      bool empty=true;
      for (int i=0;i<numberColumns_;i++) {
	if (lowerColumnMove[i]) {
	  empty=false;
	break;
	}
      }
      if (empty) {
	delete [] lowerColumnMove;
	lowerColumnMove=NULL;
      }
    }
    if (upperColumnMove) {
      bool empty=true;
      for (int i=0;i<numberColumns_;i++) {
	if (upperColumnMove[i]) {
	  empty=false;
	break;
	}
      }
      if (empty) {
	delete [] upperColumnMove;
	upperColumnMove=NULL;
      }
    }
    if (objectiveMove) {
      bool empty=true;
      for (int i=0;i<numberColumns_;i++) {
	if (objectiveMove[i]) {
	  empty=false;
	break;
	}
      }
      if (empty) {
	delete [] objectiveMove;
	objectiveMove=NULL;
      }
    }
    int saveScaling = scalingFlag_;
    scalingFlag_ = 0;
    int saveLogLevel = handler_->logLevel();
    if (detail>0&&!intervalTheta)
      handler_->setLogLevel(3);
    else
      handler_->setLogLevel(1);
    returnCode = parametrics(startTheta,endTheta,intervalTheta,
			     lowerColumnMove,upperColumnMove,
			     lowerRowMove,upperRowMove,
			     objectiveMove);
    scalingFlag_ = saveScaling;
    handler_->setLogLevel(saveLogLevel);
  }
  delete [] lowerRowMove;
  delete [] upperRowMove;
  delete [] lowerColumnMove;
  delete [] upperColumnMove;
  delete [] objectiveMove;
  fclose(fp);
  return returnCode;
}
int
ClpSimplexOther::parametricsLoop(double startingTheta, double & endingTheta, double reportIncrement,
                                 const double * lowerChange, const double * upperChange,
                                 const double * changeObjective, ClpDataSave & data,
                                 bool canTryQuick)
{
     // stuff is already at starting
     // For this crude version just try and go to end
     double change = 0.0;
     if (reportIncrement && canTryQuick) {
          endingTheta = CoinMin(endingTheta, startingTheta + reportIncrement);
          change = endingTheta - startingTheta;
     }
     int numberTotal = numberRows_ + numberColumns_;
     int i;
     for ( i = 0; i < numberTotal; i++) {
          lower_[i] += change * lowerChange[i];
          upper_[i] += change * upperChange[i];
          switch(getStatus(i)) {

          case basic:
          case isFree:
          case superBasic:
               break;
          case isFixed:
          case atUpperBound:
               solution_[i] = upper_[i];
               break;
          case atLowerBound:
               solution_[i] = lower_[i];
               break;
          }
          cost_[i] += change * changeObjective[i];
     }
     problemStatus_ = -1;

     // This says whether to restore things etc
     // startup will have factorized so can skip
     int factorType = 0;
     // Start check for cycles
     progress_.startCheck();
     // Say change made on first iteration
     changeMade_ = 1;
     /*
       Status of problem:
       0 - optimal
       1 - infeasible
       2 - unbounded
       -1 - iterating
       -2 - factorization wanted
       -3 - redo checking without factorization
       -4 - looks infeasible
     */
     while (problemStatus_ < 0) {
          int iRow, iColumn;
          // clear
          for (iRow = 0; iRow < 4; iRow++) {
               rowArray_[iRow]->clear();
          }

          for (iColumn = 0; iColumn < 2; iColumn++) {
               columnArray_[iColumn]->clear();
          }

          // give matrix (and model costs and bounds a chance to be
          // refreshed (normally null)
          matrix_->refresh(this);
          // may factorize, checks if problem finished
          statusOfProblemInParametrics(factorType, data);
          // Say good factorization
          factorType = 1;
          if (data.sparseThreshold_) {
               // use default at present
               factorization_->sparseThreshold(0);
               factorization_->goSparse();
          }

          // exit if victory declared
          if (problemStatus_ >= 0 && 
	      (canTryQuick || startingTheta>=endingTheta-1.0e-7) )
               break;

          // test for maximum iterations
          if (hitMaximumIterations()) {
               problemStatus_ = 3;
               break;
          }
          // Check event
          {
               int status = eventHandler_->event(ClpEventHandler::endOfFactorization);
               if (status >= 0) {
                    problemStatus_ = 5;
                    secondaryStatus_ = ClpEventHandler::endOfFactorization;
                    break;
               }
          }
          // Do iterations
	  problemStatus_=-1;
          if (canTryQuick) {
               double * saveDuals = NULL;
               reinterpret_cast<ClpSimplexDual *> (this)->whileIterating(saveDuals, 0);
          } else {
               whileIterating(startingTheta,  endingTheta, reportIncrement,
                              lowerChange, upperChange,
                              changeObjective);
	       startingTheta = endingTheta;
          }
     }
     if (!problemStatus_) {
          theta_ = change + startingTheta;
          eventHandler_->event(ClpEventHandler::theta);
          return 0;
     } else if (problemStatus_ == 10) {
          return -1;
     } else {
          return problemStatus_;
     }
}
/* Parametrics
   The code uses current bounds + theta * change (if change array not NULL)
   It starts at startingTheta and returns ending theta in endingTheta.
   If it can not reach input endingTheta return code will be 1 for infeasible,
   2 for unbounded, if error on ranges -1,  otherwise 0.
   Event handler may do more
   On exit endingTheta is maximum reached (can be used for next startingTheta)
*/
int
ClpSimplexOther::parametrics(double startingTheta, double & endingTheta,
                             const double * lowerChangeBound, const double * upperChangeBound,
                             const double * lowerChangeRhs, const double * upperChangeRhs)
{
  int savePerturbation = perturbation_;
  perturbation_ = 102; // switch off
  algorithm_ = -1;
  // extra region
  int maximumPivots = factorization_->maximumPivots();
  int numberDense = factorization_->numberDense();
  int length = numberRows_ + numberDense + maximumPivots;
  assert (!rowArray_[4]);
  rowArray_[4]=new CoinIndexedVector(length);
  assert (!rowArray_[5]);
  rowArray_[5]=new CoinIndexedVector(length);

  // save data
  ClpDataSave data = saveData();
  int numberTotal = numberRows_ + numberColumns_;
  int ratio = (2*sizeof(int))/sizeof(double);
  assert (ratio==1||ratio==2);
  // allow for unscaled - even if not needed
  int lengthArrays = 4*numberTotal+(numberTotal+2)*ratio;
  /*
    Save information and modify
  */
  double * saveLower = new double [lengthArrays];
  double * saveUpper = new double [lengthArrays];
  double * lowerCopy = saveLower+2*numberTotal;
  double * upperCopy = saveUpper+2*numberTotal;
  double * lowerChange = saveLower+numberTotal;
  double * upperChange = saveUpper+numberTotal;
  int * lowerList = (reinterpret_cast<int *>(saveLower+4*numberTotal))+2;
  int * upperList = (reinterpret_cast<int *>(saveUpper+4*numberTotal))+2;
  // To mark as odd
  char * markDone = reinterpret_cast<char *>(lowerList+numberTotal);
  memset(markDone,0,numberTotal);
  // Find theta when bounds will cross over and create arrays
  memset(lowerChange, 0, numberTotal * sizeof(double));
  memset(upperChange, 0, numberTotal * sizeof(double));
  if (lowerChangeBound)
    memcpy(lowerChange,lowerChangeBound,numberColumns_*sizeof(double));
  if (upperChangeBound)
    memcpy(upperChange,upperChangeBound,numberColumns_*sizeof(double));
  if (lowerChangeRhs)
    memcpy(lowerChange+numberColumns_,
	   lowerChangeRhs,numberRows_*sizeof(double));
  if (upperChangeRhs)
    memcpy(upperChange+numberColumns_,
	   upperChangeRhs,numberRows_*sizeof(double));
  int nLowerChange=0;
  int nUpperChange=0;
  for (int i=0;i<numberColumns_;i++) {
    if (lowerChange[i]) { 
      lowerList[nLowerChange++]=i;
    }
    if (upperChange[i]) { 
      upperList[nUpperChange++]=i;
    }
  }
  lowerList[-2]=nLowerChange;
  upperList[-2]=nUpperChange;
  for (int i=numberColumns_;i<numberTotal;i++) {
    if (lowerChange[i]) { 
      lowerList[nLowerChange++]=i;
    }
    if (upperChange[i]) { 
      upperList[nUpperChange++]=i;
    }
  }
  lowerList[-1]=nLowerChange;
  upperList[-1]=nUpperChange;
  memcpy(lowerCopy,columnLower_,numberColumns_*sizeof(double));
  memcpy(upperCopy,columnUpper_,numberColumns_*sizeof(double));
  memcpy(lowerCopy+numberColumns_,
	 rowLower_,numberRows_*sizeof(double));
  memcpy(upperCopy+numberColumns_,
	 rowUpper_,numberRows_*sizeof(double));
  {
    //  extra for unscaled
    double * unscaledCopy;
    unscaledCopy = lowerCopy + numberTotal; 
    memcpy(unscaledCopy,columnLower_,numberColumns_*sizeof(double));
    memcpy(unscaledCopy+numberColumns_,
	   rowLower_,numberRows_*sizeof(double));
    unscaledCopy = upperCopy + numberTotal; 
    memcpy(unscaledCopy,columnUpper_,numberColumns_*sizeof(double));
    memcpy(unscaledCopy+numberColumns_,
	   rowUpper_,numberRows_*sizeof(double));
  }
  double maxTheta = 1.0e50;
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    double lower = rowLower_[iRow];
    double upper = rowUpper_[iRow];
    if (lower<-1.0e30)
      lowerChange[numberColumns_+iRow]=0.0;
    double chgLower = lowerChange[numberColumns_+iRow];
    if (upper>1.0e30)
      upperChange[numberColumns_+iRow]=0.0;
    double chgUpper = upperChange[numberColumns_+iRow];
    if (lower > -1.0e30 && upper < 1.0e30) {
      if (lower + maxTheta * chgLower > upper + maxTheta * chgUpper) {
	maxTheta = (upper - lower) / (chgLower - chgUpper);
      }
    }
    lower+=startingTheta*chgLower;
    upper+=startingTheta*chgUpper;
    if (lower > upper) {
      maxTheta = -1.0;
      break;
    }
    rowLower_[iRow]=lower;
    rowUpper_[iRow]=upper;
  }
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    double lower = columnLower_[iColumn];
    double upper = columnUpper_[iColumn];
    if (lower<-1.0e30)
      lowerChange[iColumn]=0.0;
    double chgLower = lowerChange[iColumn];
    if (upper>1.0e30)
      upperChange[iColumn]=0.0;
    double chgUpper = upperChange[iColumn];
    if (lower > -1.0e30 && upper < 1.0e30) {
      if (lower + maxTheta * chgLower > upper + maxTheta * chgUpper) {
	maxTheta = (upper - lower) / (chgLower - chgUpper);
      }
    }
    lower+=startingTheta*chgLower;
    upper+=startingTheta*chgUpper;
    if (lower > upper) {
      maxTheta = -1.0;
      break;
    }
    columnLower_[iColumn]=lower;
    columnUpper_[iColumn]=upper;
  }
  if (maxTheta == 1.0e50)
    maxTheta = COIN_DBL_MAX;
  int returnCode=0;
  if (maxTheta < 0.0) {
    // bad ranges or initial
    returnCode = -1;
  }
  if (maxTheta < endingTheta) {
    char line[100];
    sprintf(line,"Crossover considerations reduce ending  theta from %g to %g\n", 
	    endingTheta,maxTheta);
    handler_->message(CLP_GENERAL,messages_)
      << line << CoinMessageEol;
    endingTheta = maxTheta;
  }
  if (endingTheta < startingTheta) {
    // bad initial
    returnCode = -2;
  }
  bool swapped=false;
  // Dantzig 
#define ALL_DANTZIG
#ifdef ALL_DANTZIG
  ClpDualRowPivot * savePivot = dualRowPivot_;
  dualRowPivot_ = new ClpDualRowDantzig();
  dualRowPivot_->setModel(this);
#else
  ClpDualRowPivot * savePivot = NULL;
#endif
  if (!returnCode) {
    returnCode = reinterpret_cast<ClpSimplexDual *> (this)->startupSolve(0, NULL, 0);
    if (!returnCode) {
      swapped=true;
      double * temp;
      memcpy(saveLower,lower_,numberTotal*sizeof(double));
      temp=saveLower;
      saveLower=lower_;
      lower_=temp;
      //columnLowerWork_ = lower_;
      //rowLowerWork_ = lower_ + numberColumns_;
      memcpy(saveUpper,upper_,numberTotal*sizeof(double));
      temp=saveUpper;
      saveUpper=upper_;
      upper_=temp;
      //columnUpperWork_ = upper_;
      //rowUpperWork_ = upper_ + numberColumns_;
      if (rowScale_) {
	// scale saved and change arrays
	double * lowerChange = lower_+numberTotal;
	double * upperChange = upper_+numberTotal;
	double * lowerSave = lowerChange+numberTotal;
	double * upperSave = upperChange+numberTotal;
	for (int i=0;i<numberColumns_;i++) {
	  double multiplier = inverseColumnScale_[i];
	  if (lowerSave[i]>-1.0e20)
	    lowerSave[i] *= multiplier;
	  if (upperSave[i]<1.0e20)
	    upperSave[i] *= multiplier;
	  lowerChange[i] *= multiplier;
	  upperChange[i] *= multiplier;
	}
	lowerChange += numberColumns_;
	upperChange += numberColumns_;
	lowerSave += numberColumns_;
	upperSave += numberColumns_;
	for (int i=0;i<numberRows_;i++) {
	  double multiplier = rowScale_[i];
	  if (lowerSave[i]>-1.0e20)
	    lowerSave[i] *= multiplier;
	  if (upperSave[i]<1.0e20)
	    upperSave[i] *= multiplier;
	  lowerChange[i] *= multiplier;
	  upperChange[i] *= multiplier;
	}
      }
      //double saveEndingTheta = endingTheta;
      double * saveDuals = NULL;
      reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(0, saveDuals, -1, data);
      if (numberPrimalInfeasibilities_&&sumPrimalInfeasibilities_<1.0e-4) {
	// probably can get rid of this if we adjust every change in theta
	//printf("INFEAS_A %d %g\n",numberPrimalInfeasibilities_,
	//   sumPrimalInfeasibilities_);
	int pass=100;
	while(sumPrimalInfeasibilities_) {
	  pass--;
	  if (!pass)
	    break;
	  problemStatus_=-1;
	  for (int iSequence=numberColumns_;iSequence<numberTotal;iSequence++) {
	    double value=solution_[iSequence];
	    // remember scaling
	    if (value<lower_[iSequence]-1.0e-9) {
	      lower_[iSequence]=value;
	      lowerCopy[iSequence]=value;
	    } else if (value>upper_[iSequence]+1.0e-9) {
	      upper_[iSequence]=value;
	      upperCopy[iSequence]=value;
	    }
	  }
	  reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(1, saveDuals, -1, data);
	}
      }
      assert (!problemStatus_);
      if (nLowerChange||nUpperChange) {
#ifndef ALL_DANTZIG
	// Dantzig 
	savePivot = dualRowPivot_;
	dualRowPivot_ = new ClpDualRowDantzig();
	dualRowPivot_->setModel(this);
#endif
	for (int i=0;i<numberRows_+numberColumns_;i++)
	  setFakeBound(i, noFake);
	// Now do parametrics
	handler_->message(CLP_PARAMETRICS_STATS, messages_)
	  << startingTheta << objectiveValue() << CoinMessageEol;
	bool canSkipFactorization=true;
	while (!returnCode) {
	  returnCode = parametricsLoop(startingTheta, endingTheta,
				       data,canSkipFactorization);
	  canSkipFactorization=false;
	  if (!returnCode) {
	    //startingTheta = endingTheta;
	    //endingTheta = saveEndingTheta;
	    handler_->message(CLP_PARAMETRICS_STATS, messages_)
	      << startingTheta << objectiveValue() << CoinMessageEol;
	    if (startingTheta >= endingTheta-primalTolerance_
		||problemStatus_==2)
	      break;
	  } else if (returnCode == -1) {
	    // trouble - do external solve
	    abort(); //needToDoSomething = true;
	  } else if (problemStatus_==1) {
	    // can't move any further
	    handler_->message(CLP_PARAMETRICS_STATS, messages_)
	      << endingTheta << objectiveValue() << CoinMessageEol;
	    problemStatus_=0;
	  }
	}
      }
      //reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(0, saveDuals, -1, data);
    }
    if (problemStatus_==2) {
      delete [] ray_;
      ray_ = new double [numberColumns_];
    }
    reinterpret_cast<ClpSimplexDual *> (this)->finishSolve(1);
    if (swapped&&lower_) {
      double * temp=saveLower;
      saveLower=lower_;
      lower_=temp;
      temp=saveUpper;
      saveUpper=upper_;
      upper_=temp;
    }
  }    
  if (!scalingFlag_) {
    memcpy(columnLower_,lowerCopy,numberColumns_*sizeof(double));
    memcpy(columnUpper_,upperCopy,numberColumns_*sizeof(double));
    memcpy(rowLower_,lowerCopy+numberColumns_,
	   numberRows_*sizeof(double));
    memcpy(rowUpper_,upperCopy+numberColumns_,
	   numberRows_*sizeof(double));
  } else {
    //  extra for unscaled
    double * unscaledCopy;
    unscaledCopy = lowerCopy + numberTotal; 
    memcpy(columnLower_,unscaledCopy,numberColumns_*sizeof(double));
    memcpy(rowLower_,unscaledCopy+numberColumns_,
	   numberRows_*sizeof(double));
    unscaledCopy = upperCopy + numberTotal; 
    memcpy(columnUpper_,unscaledCopy,numberColumns_*sizeof(double));
    memcpy(rowUpper_,unscaledCopy+numberColumns_,
	   numberRows_*sizeof(double));
  }
  delete [] saveLower;
  delete [] saveUpper;
#ifdef ALL_DANTZIG
  if (savePivot) {
#endif
    delete dualRowPivot_;
    dualRowPivot_ = savePivot;
#ifdef ALL_DANTZIG
  }
#endif
  // Restore any saved stuff
  restoreData(data);
  perturbation_ = savePerturbation;
  delete rowArray_[4];
  rowArray_[4]=NULL;
  delete rowArray_[5];
  rowArray_[5]=NULL;
  char line[100];
  sprintf(line,"Ending theta %g\n", endingTheta);
  handler_->message(CLP_GENERAL,messages_)
    << line << CoinMessageEol;
  return problemStatus_;
}
int
ClpSimplexOther::parametricsLoop(double & startingTheta, double & endingTheta,
                                 ClpDataSave & data,bool canSkipFactorization)
{
  int numberTotal = numberRows_+numberColumns_;
  const double * lowerChange = lower_+numberTotal;
  const double * upperChange = upper_+numberTotal;
  // stuff is already at starting
  int * lowerList = (reinterpret_cast<int *>(lower_+4*numberTotal))+2;
  int * upperList = (reinterpret_cast<int *>(upper_+4*numberTotal))+2;
  problemStatus_ = -1;
  //double saveEndingTheta=endingTheta;

  // This says whether to restore things etc
  // startup will have factorized so can skip
  int factorType = 0;
  // Start check for cycles
  progress_.startCheck();
  // Say change made on first iteration
  changeMade_ = 1;
  /*
    Status of problem:
    0 - optimal
    1 - infeasible
    2 - unbounded
    -1 - iterating
    -2 - factorization wanted
    -3 - redo checking without factorization
    -4 - looks infeasible
  */
  while (problemStatus_ < 0) {
    int iRow, iColumn;
    // clear
    for (iRow = 0; iRow < 6; iRow++) {
      rowArray_[iRow]->clear();
    }
    
    for (iColumn = 0; iColumn < 2; iColumn++) {
      columnArray_[iColumn]->clear();
    }
    
    // give matrix (and model costs and bounds a chance to be
    // refreshed (normally null)
    matrix_->refresh(this);
    // may factorize, checks if problem finished
    if (!canSkipFactorization)
      statusOfProblemInParametrics(factorType, data);
    canSkipFactorization=false;
    if (numberPrimalInfeasibilities_) {
      // probably can get rid of this if we adjust every change in theta
      //printf("INFEAS %d %g\n",numberPrimalInfeasibilities_,
      //     sumPrimalInfeasibilities_);
      const double * lowerChange = lower_+numberTotal;
      const double * upperChange = upper_+numberTotal;
      const double * startLower = lowerChange+numberTotal;
      const double * startUpper = upperChange+numberTotal;
      //startingTheta -= 1.0e-7;
      int nLowerChange = lowerList[-1];
      for (int i = 0; i < nLowerChange; i++) {
	int iSequence = lowerList[i];
	lower_[iSequence] = startLower[iSequence] + startingTheta * lowerChange[iSequence];
      }
      int nUpperChange = upperList[-1];
      for (int i = 0; i < nUpperChange; i++) {
	int iSequence = upperList[i];
	upper_[iSequence] = startUpper[iSequence] + startingTheta * upperChange[iSequence];
      }
      // adjust rhs in case dual uses
      memcpy(columnLower_,lower_,numberColumns_*sizeof(double));
      memcpy(rowLower_,lower_+numberColumns_,numberRows_*sizeof(double));
      memcpy(columnUpper_,upper_,numberColumns_*sizeof(double));
      memcpy(rowUpper_,upper_+numberColumns_,numberRows_*sizeof(double));
      if (rowScale_) {
	for (int i=0;i<numberColumns_;i++) {
	  double multiplier = columnScale_[i];
	  if (columnLower_[i]>-1.0e20)
	    columnLower_[i] *= multiplier;
	  if (columnUpper_[i]<1.0e20)
	    columnUpper_[i] *= multiplier;
	}
	for (int i=0;i<numberRows_;i++) {
	  double multiplier = inverseRowScale_[i];
	  if (rowLower_[i]>-1.0e20)
	    rowLower_[i] *= multiplier;
	  if (rowUpper_[i]<1.0e20)
	    rowUpper_[i] *= multiplier;
	}
      }
      double * saveDuals = NULL;
      problemStatus_=-1;
      reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(0, saveDuals, -1, data);
      int pass=100;
      double moved=0.0;
      while(sumPrimalInfeasibilities_) {
	pass--;
	if (!pass)
	  break;
	problemStatus_=-1;
	for (int iSequence=numberColumns_;iSequence<numberTotal;iSequence++) {
	  double value=solution_[iSequence];
	  if (value<lower_[iSequence]-1.0e-9) {
	    moved += lower_[iSequence]-value;
	    lower_[iSequence]=value;
	  } else if (value>upper_[iSequence]+1.0e-9) {
	    moved = upper_[iSequence]-value;
	    upper_[iSequence]=value;
	  }
	}
	reinterpret_cast<ClpSimplexDual *> (this)->gutsOfDual(1, saveDuals, -1, data);
      }
      // adjust
      //printf("Should adjust - moved %g\n",moved);
    }
    // Say good factorization
    factorType = 1;
    if (data.sparseThreshold_) {
      // use default at present
      factorization_->sparseThreshold(0);
      factorization_->goSparse();
    }
    
    // exit if victory declared
    if (problemStatus_ >= 0 && startingTheta>=endingTheta-1.0e-7 )
      break;
    
    // test for maximum iterations
    if (hitMaximumIterations()) {
      problemStatus_ = 3;
      break;
    }
#ifdef CLP_USER_DRIVEN
    // Check event
    {
      int status = eventHandler_->event(ClpEventHandler::endOfFactorization);
      if (status >= 0) {
	problemStatus_ = 5;
	secondaryStatus_ = ClpEventHandler::endOfFactorization;
	break;
      }
    }
#endif
    // Do iterations
    problemStatus_=-1;
    whileIterating(startingTheta,  endingTheta, 0.0,
		   lowerChange, upperChange,
		   NULL);
    //startingTheta = endingTheta;
    //endingTheta = saveEndingTheta;
  }
  if (!problemStatus_/*||problemStatus_==2*/) {
    theta_ = endingTheta;
#ifdef CLP_USER_DRIVEN
    {
      double saveTheta=theta_;
      theta_ = endingTheta;
      int status=eventHandler_->event(ClpEventHandler::theta);
      if (status>=0&&status<10) {
	endingTheta=theta_;
	theta_=saveTheta;
	problemStatus_=-1;
      } else {
	if (status>=10) {
	  problemStatus_=status-10;
	  startingTheta=endingTheta;
	}
	theta_=saveTheta;
      }
    }
#endif
    return 0;
  } else if (problemStatus_ == 10) {
    return -1;
  } else {
    return problemStatus_;
  }
}
/* Checks if finished.  Updates status */
void
ClpSimplexOther::statusOfProblemInParametrics(int type, ClpDataSave & saveData)
{
     if (type == 2) {
          // trouble - go to recovery
          problemStatus_ = 10;
          return;
     }
     if (problemStatus_ > -3 || factorization_->pivots()) {
          // factorize
          // later on we will need to recover from singularities
          // also we could skip if first time
          if (type) {
               // is factorization okay?
               if (internalFactorize(1)) {
                    // trouble - go to recovery
                    problemStatus_ = 10;
                    return;
               }
          }
          if (problemStatus_ != -4 || factorization_->pivots() > 10)
               problemStatus_ = -3;
     }
     // at this stage status is -3 or -4 if looks infeasible
     // get primal and dual solutions
     gutsOfSolution(NULL, NULL);
     double realDualInfeasibilities = sumDualInfeasibilities_;
     // If bad accuracy treat as singular
     if ((largestPrimalError_ > 1.0e15 || largestDualError_ > 1.0e15) && numberIterations_) {
          // trouble - go to recovery
          problemStatus_ = 10;
          return;
     } else if (largestPrimalError_ < 1.0e-7 && largestDualError_ < 1.0e-7) {
          // Can reduce tolerance
          double newTolerance = CoinMax(0.99 * factorization_->pivotTolerance(), saveData.pivotTolerance_);
          factorization_->pivotTolerance(newTolerance);
     }
     // Check if looping
     int loop;
     if (type != 2)
          loop = progress_.looping();
     else
          loop = -1;
     if (loop >= 0) {
          problemStatus_ = loop; //exit if in loop
          if (!problemStatus_) {
               // declaring victory
               numberPrimalInfeasibilities_ = 0;
               sumPrimalInfeasibilities_ = 0.0;
          } else {
               problemStatus_ = 10; // instead - try other algorithm
          }
          return;
     } else if (loop < -1) {
          // something may have changed
          gutsOfSolution(NULL, NULL);
     }
     progressFlag_ = 0; //reset progress flag
     if (handler_->detail(CLP_SIMPLEX_STATUS, messages_) < 100) {
          handler_->message(CLP_SIMPLEX_STATUS, messages_)
                    << numberIterations_ << objectiveValue();
          handler_->printing(sumPrimalInfeasibilities_ > 0.0)
                    << sumPrimalInfeasibilities_ << numberPrimalInfeasibilities_;
          handler_->printing(sumDualInfeasibilities_ > 0.0)
                    << sumDualInfeasibilities_ << numberDualInfeasibilities_;
          handler_->printing(numberDualInfeasibilitiesWithoutFree_
                             < numberDualInfeasibilities_)
                    << numberDualInfeasibilitiesWithoutFree_;
          handler_->message() << CoinMessageEol;
     }
     /* If we are primal feasible and any dual infeasibilities are on
        free variables then it is better to go to primal */
     if (!numberPrimalInfeasibilities_ && !numberDualInfeasibilitiesWithoutFree_ &&
               numberDualInfeasibilities_) {
          problemStatus_ = 10;
          return;
     }

     // check optimal
     // give code benefit of doubt
     if (sumOfRelaxedDualInfeasibilities_ == 0.0 &&
               sumOfRelaxedPrimalInfeasibilities_ == 0.0) {
          // say optimal (with these bounds etc)
          numberDualInfeasibilities_ = 0;
          sumDualInfeasibilities_ = 0.0;
          numberPrimalInfeasibilities_ = 0;
          sumPrimalInfeasibilities_ = 0.0;
     }
     if (dualFeasible() || problemStatus_ == -4) {
          progress_.modifyObjective(objectiveValue_
                                    - sumDualInfeasibilities_ * dualBound_);
     }
     if (numberPrimalInfeasibilities_) {
          if (problemStatus_ == -4 || problemStatus_ == -5) {
               problemStatus_ = 1; // infeasible
          }
     } else if (numberDualInfeasibilities_) {
          // clean up
          problemStatus_ = 10;
     } else {
          problemStatus_ = 0;
     }
     lastGoodIteration_ = numberIterations_;
     if (problemStatus_ < 0) {
          sumDualInfeasibilities_ = realDualInfeasibilities; // back to say be careful
          if (sumDualInfeasibilities_)
               numberDualInfeasibilities_ = 1;
     }
     // Allow matrices to be sorted etc
     int fake = -999; // signal sort
     matrix_->correctSequence(this, fake, fake);
}
//static double lastThetaX=0.0;
/* This has the flow between re-factorizations
   Reasons to come out:
   -1 iterations etc
   -2 inaccuracy
   -3 slight inaccuracy (and done iterations)
   +0 looks optimal (might be unbounded - but we will investigate)
   +1 looks infeasible
   +3 max iterations
   +4 accuracy problems
*/
int
ClpSimplexOther::whileIterating(double & startingTheta, double & endingTheta, double /*reportIncrement*/,
                                const double * lowerChange, const double * upperChange,
                                const double * /*changeObjective*/)
{
  int numberTotal = numberColumns_ + numberRows_;
  const int * lowerList = (reinterpret_cast<int *>(lower_+4*numberTotal))+2;
  const int * upperList = (reinterpret_cast<int *>(upper_+4*numberTotal))+2;
     {
          int i;
          for (i = 0; i < 4; i++) {
               rowArray_[i]->clear();
          }
          for (i = 0; i < 2; i++) {
               columnArray_[i]->clear();
          }
     }
     // if can't trust much and long way from optimal then relax
     if (largestPrimalError_ > 10.0)
          factorization_->relaxAccuracyCheck(CoinMin(1.0e2, largestPrimalError_ / 10.0));
     else
          factorization_->relaxAccuracyCheck(1.0);
     // status stays at -1 while iterating, >=0 finished, -2 to invert
     // status -3 to go to top without an invert
     int returnCode = -1;
     double lastTheta = startingTheta;
     double useTheta = startingTheta;
     while (problemStatus_ == -1) {
          double increaseTheta = CoinMin(endingTheta - lastTheta, 1.0e50);
          // Get theta for bounds - we know can't crossover
          int pivotType = nextTheta(1, increaseTheta, 
                                    lowerChange, upperChange, NULL);
	  useTheta += theta_;
	  double change = useTheta - lastTheta;
	  if (change>1.0e-14) {
	    int n;
	    n=lowerList[-1];
	    for (int i=0;i<n;i++) {
	      int iSequence = lowerList[i];
	      lower_[iSequence] += change * lowerChange[iSequence];
	      if(getStatus(iSequence)==atLowerBound) {
		solution_[iSequence] = lower_[iSequence];
	      }
	    }
	    n=upperList[-1];
	    for (int i=0;i<n;i++) {
	      int iSequence = upperList[i];
	      upper_[iSequence] += change * upperChange[iSequence];
	      if(getStatus(iSequence)==atUpperBound||
		 getStatus(iSequence)==isFixed) {
		solution_[iSequence] = upper_[iSequence];
	      }
	    }
	  }
	  sequenceIn_=-1;
          if (pivotType) {
	    if (useTheta>lastTheta+1.0e-9) {
	      handler_->message(CLP_PARAMETRICS_STATS, messages_)
		<< useTheta << objectiveValue() << CoinMessageEol;
	      lastTheta = useTheta;
	    }
	    problemStatus_ = -2;
	    if (!factorization_->pivots()&&pivotRow_<0)
	      problemStatus_=2;
#ifdef CLP_USER_DRIVEN
	    {
	      double saveTheta=theta_;
	      theta_ = endingTheta;
	      int status=eventHandler_->event(ClpEventHandler::theta);
	      if (status>=0&&status<10) {
		endingTheta=theta_;
		theta_=saveTheta;
		problemStatus_=-1;
		continue;
	      } else {
		if (status>=10)
		  problemStatus_=status-10;
		if (status<0)
		  startingTheta = useTheta;
		theta_=saveTheta;
	      }
	    }
#else
	    startingTheta = useTheta;
#endif
	    return 4;
	  }
          // choose row to go out
          //reinterpret_cast<ClpSimplexDual *> ( this)->dualRow(-1);
          if (pivotRow_ >= 0) {
               // we found a pivot row
               if (handler_->detail(CLP_SIMPLEX_PIVOTROW, messages_) < 100) {
                    handler_->message(CLP_SIMPLEX_PIVOTROW, messages_)
                              << pivotRow_
                              << CoinMessageEol;
               }
               // check accuracy of weights
               dualRowPivot_->checkAccuracy();
               // do ratio test for normal iteration
               double bestPossiblePivot = bestPivot();
               if (sequenceIn_ >= 0) {
                    // normal iteration
                    // update the incoming column
                    double btranAlpha = -alpha_ * directionOut_; // for check
                    unpackPacked(rowArray_[1]);
                    // and update dual weights (can do in parallel - with extra array)
		    rowArray_[2]->clear();
                    alpha_ = dualRowPivot_->updateWeights(rowArray_[0],
                                                          rowArray_[2],
                                                          rowArray_[3],
                                                          rowArray_[1]);
                    // see if update stable
#ifdef CLP_DEBUG
                    if ((handler_->logLevel() & 32))
                         printf("btran alpha %g, ftran alpha %g\n", btranAlpha, alpha_);
#endif
                    double checkValue = 1.0e-7;
                    // if can't trust much and long way from optimal then relax
                    if (largestPrimalError_ > 10.0)
                         checkValue = CoinMin(1.0e-4, 1.0e-8 * largestPrimalError_);
                    if (fabs(btranAlpha) < 1.0e-12 || fabs(alpha_) < 1.0e-12 ||
                              fabs(btranAlpha - alpha_) > checkValue*(1.0 + fabs(alpha_))) {
                         handler_->message(CLP_DUAL_CHECK, messages_)
                                   << btranAlpha
                                   << alpha_
                                   << CoinMessageEol;
			 // clear arrays
			 rowArray_[4]->clear();
                         if (factorization_->pivots()) {
                              dualRowPivot_->unrollWeights();
                              problemStatus_ = -2; // factorize now
                              rowArray_[0]->clear();
                              rowArray_[1]->clear();
                              columnArray_[0]->clear();
                              returnCode = -2;
                              break;
                         } else {
                              // take on more relaxed criterion
                              double test;
                              if (fabs(btranAlpha) < 1.0e-8 || fabs(alpha_) < 1.0e-8)
                                   test = 1.0e-1 * fabs(alpha_);
                              else
                                   test = 1.0e-4 * (1.0 + fabs(alpha_));
                              if (fabs(btranAlpha) < 1.0e-12 || fabs(alpha_) < 1.0e-12 ||
                                        fabs(btranAlpha - alpha_) > test) {
                                   dualRowPivot_->unrollWeights();
                                   // need to reject something
                                   char x = isColumn(sequenceOut_) ? 'C' : 'R';
                                   handler_->message(CLP_SIMPLEX_FLAG, messages_)
                                             << x << sequenceWithin(sequenceOut_)
                                             << CoinMessageEol;
                                   setFlagged(sequenceOut_);
                                   progress_.clearBadTimes();
                                   lastBadIteration_ = numberIterations_; // say be more cautious
                                   rowArray_[0]->clear();
                                   rowArray_[1]->clear();
                                   columnArray_[0]->clear();
                                   if (fabs(alpha_) < 1.0e-10 && fabs(btranAlpha) < 1.0e-8 && numberIterations_ > 100) {
                                        //printf("I think should declare infeasible\n");
                                        problemStatus_ = 1;
                                        returnCode = 1;
                                        break;
                                   }
                                   continue;
                              }
                         }
                    }
                    // update duals BEFORE replaceColumn so can do updateColumn
                    double objectiveChange = 0.0;
                    // do duals first as variables may flip bounds
                    // rowArray_[0] and columnArray_[0] may have flips
                    // so use rowArray_[3] for work array from here on
                    int nswapped = 0;
                    //rowArray_[0]->cleanAndPackSafe(1.0e-60);
                    //columnArray_[0]->cleanAndPackSafe(1.0e-60);
                    nswapped = reinterpret_cast<ClpSimplexDual *> ( this)->updateDualsInDual(rowArray_[0], columnArray_[0],
                               rowArray_[2], theta_,
                               objectiveChange, false);
		    assert (!nswapped);
                    // which will change basic solution
                    if (nswapped) {
                         factorization_->updateColumn(rowArray_[3], rowArray_[2]);
                         dualRowPivot_->updatePrimalSolution(rowArray_[2],
                                                             1.0, objectiveChange);
                         // recompute dualOut_
                         valueOut_ = solution_[sequenceOut_];
                         if (directionOut_ < 0) {
                              dualOut_ = valueOut_ - upperOut_;
                         } else {
                              dualOut_ = lowerOut_ - valueOut_;
                         }
                    }
                    // amount primal will move
                    double movement = -dualOut_ * directionOut_ / alpha_;
                    // so objective should increase by fabs(dj)*movement
                    // but we already have objective change - so check will be good
                    if (objectiveChange + fabs(movement * dualIn_) < -1.0e-5) {
#ifdef CLP_DEBUG
                         if (handler_->logLevel() & 32)
                              printf("movement %g, swap change %g, rest %g  * %g\n",
                                     objectiveChange + fabs(movement * dualIn_),
                                     objectiveChange, movement, dualIn_);
#endif
			 assert (objectiveChange + fabs(movement * dualIn_) >= -1.0e-5);
                         if(factorization_->pivots()) {
                              // going backwards - factorize
                              dualRowPivot_->unrollWeights();
                              problemStatus_ = -2; // factorize now
                              returnCode = -2;
                              break;
                         }
                    }
                    CoinAssert(fabs(dualOut_) < 1.0e50);
                    // if stable replace in basis
                    int updateStatus = factorization_->replaceColumn(this,
                                       rowArray_[2],
                                       rowArray_[1],
                                       pivotRow_,
                                       alpha_);
                    // if no pivots, bad update but reasonable alpha - take and invert
                    if (updateStatus == 2 &&
                              !factorization_->pivots() && fabs(alpha_) > 1.0e-5)
                         updateStatus = 4;
                    if (updateStatus == 1 || updateStatus == 4) {
                         // slight error
                         if (factorization_->pivots() > 5 || updateStatus == 4) {
                              problemStatus_ = -2; // factorize now
                              returnCode = -3;
                         }
                    } else if (updateStatus == 2) {
                         // major error
                         dualRowPivot_->unrollWeights();
                         // later we may need to unwind more e.g. fake bounds
                         if (factorization_->pivots()) {
                              problemStatus_ = -2; // factorize now
                              returnCode = -2;
                              break;
                         } else {
                              // need to reject something
                              char x = isColumn(sequenceOut_) ? 'C' : 'R';
                              handler_->message(CLP_SIMPLEX_FLAG, messages_)
                                        << x << sequenceWithin(sequenceOut_)
                                        << CoinMessageEol;
                              setFlagged(sequenceOut_);
                              progress_.clearBadTimes();
                              lastBadIteration_ = numberIterations_; // say be more cautious
                              rowArray_[0]->clear();
                              rowArray_[1]->clear();
                              columnArray_[0]->clear();
                              // make sure dual feasible
                              // look at all rows and columns
                              double objectiveChange = 0.0;
                              reinterpret_cast<ClpSimplexDual *> ( this)->updateDualsInDual(rowArray_[0], columnArray_[0], rowArray_[1],
                                        0.0, objectiveChange, true);
                              continue;
                         }
                    } else if (updateStatus == 3) {
                         // out of memory
                         // increase space if not many iterations
                         if (factorization_->pivots() <
                                   0.5 * factorization_->maximumPivots() &&
                                   factorization_->pivots() < 200)
                              factorization_->areaFactor(
                                   factorization_->areaFactor() * 1.1);
                         problemStatus_ = -2; // factorize now
                    } else if (updateStatus == 5) {
                         problemStatus_ = -2; // factorize now
                    }
		    // update change vector
		    {
		      double * work = rowArray_[1]->denseVector();
		      int number = rowArray_[1]->getNumElements();
		      int * which = rowArray_[1]->getIndices();
		      assert (!rowArray_[4]->packedMode());
		      assert (rowArray_[1]->packedMode());
		      double pivotValue = rowArray_[4]->denseVector()[pivotRow_];
		      double multiplier = -pivotValue/alpha_;
		      if (multiplier) {
			for (int i = 0; i < number; i++) {
			  int iRow = which[i];
			  rowArray_[4]->quickAddNonZero(iRow,multiplier*work[i]);
			}
		      }
		      pivotValue = rowArray_[4]->denseVector()[pivotRow_];
		      // we want pivot to be -multiplier
		      rowArray_[4]->quickAdd(pivotRow_,-multiplier-pivotValue);
		    }
                    // update primal solution
                    if (theta_ < 0.0) {
#ifdef CLP_DEBUG
                         if (handler_->logLevel() & 32)
                              printf("negative theta %g\n", theta_);
#endif
                         theta_ = 0.0;
                    }
                    // do actual flips
                    reinterpret_cast<ClpSimplexDual *> ( this)->flipBounds(rowArray_[0], columnArray_[0]);
                    //rowArray_[1]->expand();
                    dualRowPivot_->updatePrimalSolution(rowArray_[1],
                                                        movement,
                                                        objectiveChange);
                    // modify dualout
                    dualOut_ /= alpha_;
                    dualOut_ *= -directionOut_;
                    //setStatus(sequenceIn_,basic);
                    dj_[sequenceIn_] = 0.0;
                    //double oldValue = valueIn_;
                    if (directionIn_ == -1) {
                         // as if from upper bound
                         valueIn_ = upperIn_ + dualOut_;
                    } else {
                         // as if from lower bound
                         valueIn_ = lowerIn_ + dualOut_;
                    }
		    objectiveChange = 0.0;
		    for (int i=0;i<numberTotal;i++)
		      objectiveChange += solution_[i]*cost_[i];
                    objectiveChange -= objectiveValue_;
                    // outgoing
                    originalBound(sequenceOut_,useTheta,lowerChange,upperChange);
		    lowerOut_=lower_[sequenceOut_];
		    upperOut_=upper_[sequenceOut_];
                    // set dj to zero unless values pass
                    if (directionOut_ > 0) {
                         valueOut_ = lowerOut_;
                         dj_[sequenceOut_] = theta_;
                    } else {
                         valueOut_ = upperOut_;
                         dj_[sequenceOut_] = -theta_;
                    }
                    solution_[sequenceOut_] = valueOut_;
                    int whatNext = housekeeping(objectiveChange);
		    {
		      char in[200],out[200];
		      int iSequence=sequenceIn_;
		      if (iSequence<numberColumns_) {
			if (lengthNames_) 
			  strcpy(in,columnNames_[iSequence].c_str());
			 else 
			  sprintf(in,"C%7.7d",iSequence);
		      } else {
			iSequence -= numberColumns_;
			if (lengthNames_) 
			  strcpy(in,rowNames_[iSequence].c_str());
			 else 
			  sprintf(in,"R%7.7d",iSequence);
		      }
		      iSequence=sequenceOut_;
		      if (iSequence<numberColumns_) {
			if (lengthNames_) 
			  strcpy(out,columnNames_[iSequence].c_str());
			 else 
			  sprintf(out,"C%7.7d",iSequence);
		      } else {
			iSequence -= numberColumns_;
			if (lengthNames_) 
			  strcpy(out,rowNames_[iSequence].c_str());
			 else 
			  sprintf(out,"R%7.7d",iSequence);
		      }
		      handler_->message(CLP_PARAMETRICS_STATS2, messages_)
			<< useTheta << objectiveValue() 
			<< in << out << CoinMessageEol;
		    }
		    if (useTheta>lastTheta+1.0e-9) {
		      handler_->message(CLP_PARAMETRICS_STATS, messages_)
			<< useTheta << objectiveValue() << CoinMessageEol;
		      lastTheta = useTheta;
		    }
                    // and set bounds correctly
                    originalBound(sequenceIn_,useTheta,lowerChange,upperChange);
                    reinterpret_cast<ClpSimplexDual *> ( this)->changeBound(sequenceOut_);
                    if (whatNext == 1) {
                         problemStatus_ = -2; // refactorize
                    } else if (whatNext == 2) {
                         // maximum iterations or equivalent
                         problemStatus_ = 3;
                         returnCode = 3;
                         break;
                    }
#ifdef CLP_USER_DRIVEN
                    // Check event
                    {
                         int status = eventHandler_->event(ClpEventHandler::endOfIteration);
                         if (status >= 0) {
                              problemStatus_ = 5;
                              secondaryStatus_ = ClpEventHandler::endOfIteration;
                              returnCode = 4;
                              break;
                         }
                    }
#endif
               } else {
                    // no incoming column is valid
#ifdef CLP_USER_DRIVEN
		 rowArray_[0]->clear();
		 columnArray_[0]->clear();
		 theta_ = useTheta;
		 lastTheta = useTheta;
		 int action = eventHandler_->event(ClpEventHandler::noTheta);
		 if (action>=0) {
		   endingTheta=theta_;
		   theta_ = 0.0;
		   //adjust [4] from handler - but
		   rowArray_[4]->clear(); // temp
		   if (action>=0&&action<10)
		     problemStatus_=-1; // carry on
		   else if (action==15)
		     problemStatus_ =5; // say stopped
		   returnCode = 1;
		   if (action==0||action>=10) 
		     break;
		   else
		     continue;
		 } else {
		 theta_ = 0.0;
		 }
#endif
                    pivotRow_ = -1;
#ifdef CLP_DEBUG
                    if (handler_->logLevel() & 32)
                         printf("** no column pivot\n");
#endif
                    if (factorization_->pivots() < 10) { 
                         // If we have just factorized and infeasibility reasonable say infeas
                         if (((specialOptions_ & 4096) != 0 || bestPossiblePivot < 1.0e-11) && dualBound_ > 1.0e8) {
                              if (valueOut_ > upperOut_ + 1.0e-3 || valueOut_ < lowerOut_ - 1.0e-3
                                        || (specialOptions_ & 64) == 0) {
                                   // say infeasible
                                   problemStatus_ = 1;
                                   // unless primal feasible!!!!
                                   //printf("%d %g %d %g\n",numberPrimalInfeasibilities_,sumPrimalInfeasibilities_,
                                   //     numberDualInfeasibilities_,sumDualInfeasibilities_);
                                   if (numberDualInfeasibilities_)
                                        problemStatus_ = 10;
                                   rowArray_[0]->clear();
                                   columnArray_[0]->clear();
                              }
                         }
                         // If special option set - put off as long as possible
                         if ((specialOptions_ & 64) == 0) {
                              problemStatus_ = -4; //say looks infeasible
                         } else {
                              // flag
                              char x = isColumn(sequenceOut_) ? 'C' : 'R';
                              handler_->message(CLP_SIMPLEX_FLAG, messages_)
                                        << x << sequenceWithin(sequenceOut_)
                                        << CoinMessageEol;
                              setFlagged(sequenceOut_);
                              if (!factorization_->pivots()) {
                                   rowArray_[0]->clear();
                                   columnArray_[0]->clear();
                                   continue;
                              }
                         }
                    }
                    rowArray_[0]->clear();
                    columnArray_[0]->clear();
                    returnCode = 1;
                    break;
               }
          } else {
               // no pivot row
#ifdef CLP_USER_DRIVEN
	    {
	      double saveTheta=theta_;
	      theta_ = endingTheta;
	      int status=eventHandler_->event(ClpEventHandler::theta);
	      if (status>=0&&status<10) {
		endingTheta=theta_;
		theta_=saveTheta;
		continue;
	      } else {
		theta_=saveTheta;
	      }
	    }
#endif
#ifdef CLP_DEBUG
               if (handler_->logLevel() & 32)
                    printf("** no row pivot\n");
#endif
               int numberPivots = factorization_->pivots();
               bool specialCase;
               int useNumberFake;
               returnCode = 0;
               if (numberPivots < 20 &&
                         (specialOptions_ & 2048) != 0 && !numberChanged_ && perturbation_ >= 100
                         && dualBound_ > 1.0e8) {
                    specialCase = true;
                    // as dual bound high - should be okay
                    useNumberFake = 0;
               } else {
                    specialCase = false;
                    useNumberFake = numberFake_;
               }
               if (!numberPivots || specialCase) {
                    // may have crept through - so may be optimal
                    // check any flagged variables
                    int iRow;
                    for (iRow = 0; iRow < numberRows_; iRow++) {
                         int iPivot = pivotVariable_[iRow];
                         if (flagged(iPivot))
                              break;
                    }
                    if (iRow < numberRows_ && numberPivots) {
                         // try factorization
                         returnCode = -2;
                    }

                    if (useNumberFake || numberDualInfeasibilities_) {
                         // may be dual infeasible
                         problemStatus_ = -5;
                    } else {
                         if (iRow < numberRows_) {
                              problemStatus_ = -5;
                         } else {
                              if (numberPivots) {
                                   // objective may be wrong
                                   objectiveValue_ = innerProduct(cost_,
                                                                  numberColumns_ + numberRows_,
                                                                  solution_);
                                   objectiveValue_ += objective_->nonlinearOffset();
                                   objectiveValue_ /= (objectiveScale_ * rhsScale_);
                                   if ((specialOptions_ & 16384) == 0) {
                                        // and dual_ may be wrong (i.e. for fixed or basic)
                                        CoinIndexedVector * arrayVector = rowArray_[1];
                                        arrayVector->clear();
                                        int iRow;
                                        double * array = arrayVector->denseVector();
                                        /* Use dual_ instead of array
                                           Even though dual_ is only numberRows_ long this is
                                           okay as gets permuted to longer rowArray_[2]
                                        */
                                        arrayVector->setDenseVector(dual_);
                                        int * index = arrayVector->getIndices();
                                        int number = 0;
                                        for (iRow = 0; iRow < numberRows_; iRow++) {
                                             int iPivot = pivotVariable_[iRow];
                                             double value = cost_[iPivot];
                                             dual_[iRow] = value;
                                             if (value) {
                                                  index[number++] = iRow;
                                             }
                                        }
                                        arrayVector->setNumElements(number);
                                        // Extended duals before "updateTranspose"
                                        matrix_->dualExpanded(this, arrayVector, NULL, 0);
                                        // Btran basic costs
                                        rowArray_[2]->clear();
                                        factorization_->updateColumnTranspose(rowArray_[2], arrayVector);
                                        // and return vector
                                        arrayVector->setDenseVector(array);
                                   }
                              }
                              problemStatus_ = 0;
                              sumPrimalInfeasibilities_ = 0.0;
                              if ((specialOptions_&(1024 + 16384)) != 0) {
                                   CoinIndexedVector * arrayVector = rowArray_[1];
                                   arrayVector->clear();
                                   double * rhs = arrayVector->denseVector();
                                   times(1.0, solution_, rhs);
                                   bool bad2 = false;
                                   int i;
                                   for ( i = 0; i < numberRows_; i++) {
                                        if (rhs[i] < rowLowerWork_[i] - primalTolerance_ ||
                                                  rhs[i] > rowUpperWork_[i] + primalTolerance_) {
                                             bad2 = true;
                                        } else if (fabs(rhs[i] - rowActivityWork_[i]) > 1.0e-3) {
                                        }
                                        rhs[i] = 0.0;
                                   }
                                   for ( i = 0; i < numberColumns_; i++) {
                                        if (solution_[i] < columnLowerWork_[i] - primalTolerance_ ||
                                                  solution_[i] > columnUpperWork_[i] + primalTolerance_) {
                                             bad2 = true;
                                        }
                                   }
                                   if (bad2) {
                                        problemStatus_ = -3;
                                        returnCode = -2;
                                        // Force to re-factorize early next time
                                        int numberPivots = factorization_->pivots();
                                        forceFactorization_ = CoinMin(forceFactorization_, (numberPivots + 1) >> 1);
                                   }
                              }
                         }
                    }
               } else {
                    problemStatus_ = -3;
                    returnCode = -2;
                    // Force to re-factorize early next time
                    int numberPivots = factorization_->pivots();
                    forceFactorization_ = CoinMin(forceFactorization_, (numberPivots + 1) >> 1);
               }
               break;
          }
     }
     startingTheta = lastTheta+theta_;
     return returnCode;
}
// Finds best possible pivot
double 
ClpSimplexOther::bestPivot(bool justColumns)
{
  // Get good size for pivot
  // Allow first few iterations to take tiny
  double acceptablePivot = 1.0e-9;
  if (numberIterations_ > 100)
    acceptablePivot = 1.0e-8;
  if (factorization_->pivots() > 10 ||
      (factorization_->pivots() && sumDualInfeasibilities_))
    acceptablePivot = 1.0e-5; // if we have iterated be more strict
  else if (factorization_->pivots() > 5)
    acceptablePivot = 1.0e-6; // if we have iterated be slightly more strict
  else if (factorization_->pivots())
    acceptablePivot = 1.0e-8; // relax
  double bestPossiblePivot = 1.0;
  // get sign for finding row of tableau
  // normal iteration
  // create as packed
  double direction = directionOut_;
  rowArray_[0]->createPacked(1, &pivotRow_, &direction);
  factorization_->updateColumnTranspose(rowArray_[1], rowArray_[0]);
  // put row of tableau in rowArray[0] and columnArray[0]
  matrix_->transposeTimes(this, -1.0,
			  rowArray_[0], rowArray_[3], columnArray_[0]);
  sequenceIn_=-1;
  if (justColumns)
    rowArray_[0]->clear();
  // do ratio test for normal iteration
  bestPossiblePivot = 
    reinterpret_cast<ClpSimplexDual *> 
    ( this)->dualColumn(rowArray_[0],
			columnArray_[0], columnArray_[1],
			rowArray_[3], acceptablePivot, NULL);
  return bestPossiblePivot;
}
// Computes next theta and says if objective or bounds (0= bounds, 1 objective, -1 none)
int
ClpSimplexOther::nextTheta(int type, double maxTheta, 
                           const double * lowerChange, const double * upperChange,
                           const double * /*changeObjective*/)
{
  int numberTotal = numberColumns_ + numberRows_;
  const int * lowerList = (reinterpret_cast<int *>(lower_+4*numberTotal))+2;
  const int * upperList = (reinterpret_cast<int *>(upper_+4*numberTotal))+2;
  int iSequence;
  theta_ = maxTheta;
  bool toLower = false;
  assert (type==1);
  // may need to decide based on model?
  bool needFullUpdate = rowArray_[4]->getNumElements()==0;
  double * array = rowArray_[4]->denseVector();
  const int * row = matrix_->getIndices();
  const int * columnLength = matrix_->getVectorLengths();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * elementByColumn = matrix_->getElements();
  if (!factorization_->pivots()||needFullUpdate) {
    rowArray_[4]->clear();
    // get change
    if (!rowScale_) {
      int n;
      n=lowerList[-2];
      int i;
      for (i=0;i<n;i++) {
	int iSequence = lowerList[i];
	assert (iSequence<numberColumns_);
	if (getColumnStatus(iSequence)==atLowerBound) {
	  double value=lowerChange[iSequence];
	  for (CoinBigIndex j = columnStart[iSequence];
	       j < columnStart[iSequence] + columnLength[iSequence]; j++) {
	    rowArray_[4]->quickAddNonZero(row[j], elementByColumn[j]*value);
	  }
	}
      }
      n=lowerList[-1];
      const double * change = lowerChange+numberColumns_;
      for (;i<n;i++) {
	int iSequence = lowerList[i]-numberColumns_;
	assert (iSequence>=0);
	if (getRowStatus(iSequence)==atLowerBound) {
	  double value=change[iSequence];
	  rowArray_[4]->quickAddNonZero(iSequence, -value);
	}
      }
      n=upperList[-2];
      for (i=0;i<n;i++) {
	int iSequence = upperList[i];
	assert (iSequence<numberColumns_);
	if (getColumnStatus(iSequence)==atUpperBound) {
	  double value=upperChange[iSequence];
	  for (CoinBigIndex j = columnStart[iSequence];
	       j < columnStart[iSequence] + columnLength[iSequence]; j++) {
	    rowArray_[4]->quickAddNonZero(row[j], elementByColumn[j]*value);
	  }
	}
      }
      n=upperList[-1];
      change = upperChange+numberColumns_;
      for (;i<n;i++) {
	int iSequence = upperList[i]-numberColumns_;
	assert (iSequence>=0);
	if (getRowStatus(iSequence)==atUpperBound) {
	  double value=change[iSequence];
	  rowArray_[4]->quickAddNonZero(iSequence, -value);
	}
      }
    } else {
      int n;
      n=lowerList[-2];
      int i;
      for (i=0;i<n;i++) {
	int iSequence = lowerList[i];
	assert (iSequence<numberColumns_);
	if (getColumnStatus(iSequence)==atLowerBound) {
	  double value=lowerChange[iSequence];
	  // apply scaling
	  double scale = columnScale_[iSequence];
	  for (CoinBigIndex j = columnStart[iSequence];
	       j < columnStart[iSequence] + columnLength[iSequence]; j++) {
	    int iRow = row[j];
	    rowArray_[4]->quickAddNonZero(iRow, elementByColumn[j]*scale * rowScale_[iRow]*value);
	  }
	}
      }
      n=lowerList[-1];
      const double * change = lowerChange+numberColumns_;
      for (;i<n;i++) {
	int iSequence = lowerList[i]-numberColumns_;
	assert (iSequence>=0);
	if (getRowStatus(iSequence)==atLowerBound) {
	  double value=change[iSequence];
	  rowArray_[4]->quickAddNonZero(iSequence, -value);
	}
      }
      n=upperList[-2];
      for (i=0;i<n;i++) {
	int iSequence = upperList[i];
	assert (iSequence<numberColumns_);
	if (getColumnStatus(iSequence)==atUpperBound) {
	  double value=upperChange[iSequence];
	  // apply scaling
	  double scale = columnScale_[iSequence];
	  for (CoinBigIndex j = columnStart[iSequence];
	       j < columnStart[iSequence] + columnLength[iSequence]; j++) {
	    int iRow = row[j];
	    rowArray_[4]->quickAddNonZero(iRow, elementByColumn[j]*scale * rowScale_[iRow]*value);
	  }
	}
      }
      n=upperList[-1];
      change = upperChange+numberColumns_;
      for (;i<n;i++) {
	int iSequence = upperList[i]-numberColumns_;
	assert (iSequence>=0);
	if (getRowStatus(iSequence)==atUpperBound) {
	  double value=change[iSequence];
	  rowArray_[4]->quickAddNonZero(iSequence, -value);
	}
      }
    }
    // ftran it
    factorization_->updateColumn(rowArray_[0], rowArray_[4]);
  } else {
    assert (sequenceIn_>=0);
    assert (sequenceOut_>=0);
    assert (sequenceIn_!=sequenceOut_);
    double change = (directionIn_>0) ? -lowerChange[sequenceIn_] : -upperChange[sequenceIn_];
    int needed=0;
    assert (!rowArray_[5]->getNumElements());
    if (change) {
      if (sequenceIn_<numberColumns_) {
	if (!rowScale_) {
	  for (CoinBigIndex i = columnStart[sequenceIn_];
	       i < columnStart[sequenceIn_] + columnLength[sequenceIn_]; i++) {
	    rowArray_[5]->quickAddNonZero(row[i], elementByColumn[i]*change);
	  }
	} else {
	  // apply scaling
	  double scale = columnScale_[sequenceIn_];
	  for (CoinBigIndex i = columnStart[sequenceIn_];
	       i < columnStart[sequenceIn_] + columnLength[sequenceIn_]; i++) {
	    int iRow = row[i];
	    rowArray_[5]->quickAddNonZero(iRow, elementByColumn[i]*scale * rowScale_[iRow]*change);
	  }
	}
      } else {
	rowArray_[5]->insert(sequenceIn_-numberColumns_,-change);
      }
      needed++;
    }
    if (getStatus(sequenceOut_)==atLowerBound)
      change=lowerChange[sequenceOut_];
    else
      change=upperChange[sequenceOut_];
    if (change) {
      if (sequenceOut_<numberColumns_) {
	if (!rowScale_) {
	  for (CoinBigIndex i = columnStart[sequenceOut_];
	       i < columnStart[sequenceOut_] + columnLength[sequenceOut_]; i++) {
	    rowArray_[5]->quickAddNonZero(row[i], elementByColumn[i]*change);
	  }
	} else {
	  // apply scaling
	  double scale = columnScale_[sequenceOut_];
	  for (CoinBigIndex i = columnStart[sequenceOut_];
	       i < columnStart[sequenceOut_] + columnLength[sequenceOut_]; i++) {
	    int iRow = row[i];
	    rowArray_[5]->quickAddNonZero(iRow, elementByColumn[i]*scale * rowScale_[iRow]*change);
	  }
	}
      } else {
	rowArray_[5]->quickAddNonZero(sequenceOut_-numberColumns_,-change);
      }
      needed++;
    }
    //printf("seqin %d seqout %d needed %d\n",
    //	   sequenceIn_,sequenceOut_,needed);
    if (needed) {
      // ftran it
      factorization_->updateColumn(rowArray_[0], rowArray_[5]);
      // add
      double * array5 = rowArray_[5]->denseVector();
      int * index5 = rowArray_[5]->getIndices();
      int number5 = rowArray_[5]->getNumElements();
      for (int i = 0; i < number5; i++) {
	int iPivot = index5[i];
	rowArray_[4]->quickAddNonZero(iPivot,array5[iPivot]);
	array5[iPivot]=0.0;
      }
      rowArray_[5]->setNumElements(0);
    }
  }
  pivotRow_ = -1;
  const int * index = rowArray_[4]->getIndices();
  int number = rowArray_[4]->getNumElements();
  int * lowerList2 = (reinterpret_cast<int *>(lower_+4*numberTotal))+2;
  char * markDone = reinterpret_cast<char *>(lowerList2+numberTotal);
#if 0
  for (int iPivot = 0; iPivot < numberRows_; iPivot++) {
    //int iPivot = index[iRow];
    iSequence = pivotVariable_[iPivot];
    // solution value will be sol - theta*alpha
    // bounds will be bounds + change *theta
    double currentSolution = solution_[iSequence];
    double currentLower = lower_[iSequence];
    double currentUpper = upper_[iSequence];
    double alpha = array[iPivot];
    assert (currentSolution >= currentLower - 100.0*primalTolerance_);
    assert (currentSolution <= currentUpper + 100.0*primalTolerance_);
    double thetaCoefficient;
    double hitsLower = COIN_DBL_MAX;
    thetaCoefficient = lowerChange[iSequence] + alpha;
    if (thetaCoefficient > 1.0e-8)
      hitsLower = (currentSolution - currentLower) / thetaCoefficient;
    //if (hitsLower < 0.0) {
    // does not hit - but should we check further
    //   hitsLower = COIN_DBL_MAX;
    //}
    double hitsUpper = COIN_DBL_MAX;
    thetaCoefficient = upperChange[iSequence] + alpha;
    if (thetaCoefficient < -1.0e-8)
      hitsUpper = (currentSolution - currentUpper) / thetaCoefficient;
    //if (hitsUpper < 0.0) {
    // does not hit - but should we check further
    //   hitsUpper = COIN_DBL_MAX;
    //}
    if (CoinMin(hitsLower, hitsUpper) < theta_) {
      theta_ = CoinMin(hitsLower, hitsUpper);
      toLower = hitsLower < hitsUpper;
      pivotRow_ = iPivot;
    }
  }
#else
  // first ones with alpha
  for (int i=0;i<number;i++) {
    int iPivot=index[i];
    iSequence = pivotVariable_[iPivot];
    assert(!markDone[iSequence]);
    markDone[iSequence]=1;
    // solution value will be sol - theta*alpha
    // bounds will be bounds + change *theta
    double currentSolution = solution_[iSequence];
    double currentLower = lower_[iSequence];
    double currentUpper = upper_[iSequence];
    double alpha = array[iPivot];
    assert (currentSolution >= currentLower - 100.0*primalTolerance_);
    assert (currentSolution <= currentUpper + 100.0*primalTolerance_);
    double thetaCoefficient;
    thetaCoefficient = lowerChange[iSequence] + alpha;
    if (thetaCoefficient > 1.0e-8) {
      double gap=currentSolution-currentLower;
      if (thetaCoefficient*theta_>gap) {
	theta_ = gap/thetaCoefficient;
	toLower=true;
	pivotRow_=iPivot;
      }
    }
    thetaCoefficient = upperChange[iSequence] + alpha;
    if (thetaCoefficient < -1.0e-8) {
      double gap=currentSolution-currentUpper; //negative
      if (thetaCoefficient*theta_<gap) {
	theta_ = gap/thetaCoefficient;
	toLower=false;
	pivotRow_=iPivot;
      }
    }
  }
  // now others
  int nLook=lowerList[-1];
  for (int i=0;i<nLook;i++) {
    int iSequence = lowerList[i];
    if (getColumnStatus(iSequence)==basic&&!markDone[iSequence]) {
      double currentSolution = solution_[iSequence];
      double currentLower = lower_[iSequence];
      assert (currentSolution >= currentLower - 100.0*primalTolerance_);
      double thetaCoefficient = lowerChange[iSequence];
      if (thetaCoefficient > 0.0) {
	double gap=currentSolution-currentLower;
	if (thetaCoefficient*theta_>gap) {
	  theta_ = gap/thetaCoefficient;
	  toLower=true;
	  pivotRow_ = -2-iSequence;
	}
      }
    }
  }
  nLook=upperList[-1];
  for (int i=0;i<nLook;i++) {
    int iSequence = upperList[i];
    if (getColumnStatus(iSequence)==basic&&!markDone[iSequence]) {
      double currentSolution = solution_[iSequence];
      double currentUpper = upper_[iSequence];
      assert (currentSolution <= currentUpper + 100.0*primalTolerance_);
      double thetaCoefficient = upperChange[iSequence];
      if (thetaCoefficient < 0) {
	double gap=currentSolution-currentUpper; //negative
	if (thetaCoefficient*theta_<gap) {
	  theta_ = gap/thetaCoefficient;
	  toLower=false;
	  pivotRow_ = -2-iSequence;
	}
      }
    }
  }
  if (pivotRow_<-1) {
    // find
    int iSequence = -pivotRow_-2;
    for (int iPivot = 0; iPivot < numberRows_; iPivot++) {
      if (iSequence == pivotVariable_[iPivot]) {
	pivotRow_=iPivot;
	break;
      }
    }
    assert (pivotRow_>=0);
  }
#endif
  theta_ = CoinMax(theta_,0.0);
  if (theta_>1.0e-15) {
    // update solution
    for (int iRow = 0; iRow < number; iRow++) {
      int iPivot = index[iRow];
      iSequence = pivotVariable_[iPivot];
      markDone[iSequence]=0;
      // solution value will be sol - theta*alpha
      double alpha = array[iPivot];
      solution_[iSequence] -= theta_ * alpha;
    }
  } else {
    for (int iRow = 0; iRow < number; iRow++) {
      int iPivot = index[iRow];
      iSequence = pivotVariable_[iPivot];
      markDone[iSequence]=0;
    }
  }
#if 0
  for (int i=0;i<numberTotal;i++)
    assert(!markDone[i]);
#endif
  if (pivotRow_ >= 0) {
    sequenceOut_ = pivotVariable_[pivotRow_];
    valueOut_ = solution_[sequenceOut_];
    lowerOut_ = lower_[sequenceOut_]+theta_*lowerChange[sequenceOut_];
    upperOut_ = upper_[sequenceOut_]+theta_*upperChange[sequenceOut_];
    if (!toLower) {
      directionOut_ = -1;
      dualOut_ = valueOut_ - upperOut_;
    } else {
      directionOut_ = 1;
      dualOut_ = lowerOut_ - valueOut_;
    }
    return 0;
  } else {
    //theta_=0.0;
    return -1;
  }
}
// Restores bound to original bound
void 
ClpSimplexOther::originalBound(int iSequence, double theta, 
			       const double * lowerChange,
			       const double * upperChange)
{
     if (getFakeBound(iSequence) != noFake) {
          numberFake_--;
          setFakeBound(iSequence, noFake);
          if (iSequence >= numberColumns_) {
               // rows
               int iRow = iSequence - numberColumns_;
               rowLowerWork_[iRow] = rowLower_[iRow]+theta*lowerChange[iSequence];
               rowUpperWork_[iRow] = rowUpper_[iRow]+theta*upperChange[iSequence];
               if (rowScale_) {
                    if (rowLowerWork_[iRow] > -1.0e50)
                         rowLowerWork_[iRow] *= rowScale_[iRow] * rhsScale_;
                    if (rowUpperWork_[iRow] < 1.0e50)
                         rowUpperWork_[iRow] *= rowScale_[iRow] * rhsScale_;
               } else if (rhsScale_ != 1.0) {
                    if (rowLowerWork_[iRow] > -1.0e50)
                         rowLowerWork_[iRow] *= rhsScale_;
                    if (rowUpperWork_[iRow] < 1.0e50)
                         rowUpperWork_[iRow] *= rhsScale_;
               }
          } else {
               // columns
               columnLowerWork_[iSequence] = columnLower_[iSequence]+theta*lowerChange[iSequence];
               columnUpperWork_[iSequence] = columnUpper_[iSequence]+theta*upperChange[iSequence];
               if (rowScale_) {
                    double multiplier = 1.0 * inverseColumnScale_[iSequence];
                    if (columnLowerWork_[iSequence] > -1.0e50)
                         columnLowerWork_[iSequence] *= multiplier * rhsScale_;
                    if (columnUpperWork_[iSequence] < 1.0e50)
                         columnUpperWork_[iSequence] *= multiplier * rhsScale_;
               } else if (rhsScale_ != 1.0) {
                    if (columnLowerWork_[iSequence] > -1.0e50)
                         columnLowerWork_[iSequence] *= rhsScale_;
                    if (columnUpperWork_[iSequence] < 1.0e50)
                         columnUpperWork_[iSequence] *= rhsScale_;
               }
          }
     }
}
/* Expands out all possible combinations for a knapsack
   If buildObj NULL then just computes space needed - returns number elements
   On entry numberOutput is maximum allowed, on exit it is number needed or
   -1 (as will be number elements) if maximum exceeded.  numberOutput will have at
   least space to return values which reconstruct input.
   Rows returned will be original rows but no entries will be returned for
   any rows all of whose entries are in knapsack.  So up to user to allow for this.
   If reConstruct >=0 then returns number of entrie which make up item "reConstruct"
   in expanded knapsack.  Values in buildRow and buildElement;
*/
int
ClpSimplexOther::expandKnapsack(int knapsackRow, int & numberOutput,
                                double * buildObj, CoinBigIndex * buildStart,
                                int * buildRow, double * buildElement, int reConstruct) const
{
     int iRow;
     int iColumn;
     // Get column copy
     CoinPackedMatrix * columnCopy = matrix();
     // Get a row copy in standard format
     CoinPackedMatrix matrixByRow;
     matrixByRow.reverseOrderedCopyOf(*columnCopy);
     const double * elementByRow = matrixByRow.getElements();
     const int * column = matrixByRow.getIndices();
     const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
     const int * rowLength = matrixByRow.getVectorLengths();
     CoinBigIndex j;
     int * whichColumn = new int [numberColumns_];
     int * whichRow = new int [numberRows_];
     int numJ = 0;
     // Get what other columns can compensate for
     double * lo = new double [numberRows_];
     double * high = new double [numberRows_];
     {
          // Use to get tight column bounds
          ClpSimplex tempModel(*this);
          tempModel.tightenPrimalBounds(0.0, 0, true);
          // Now another model without knapsacks
          int nCol = 0;
          for (iRow = 0; iRow < numberRows_; iRow++) {
               whichRow[iRow] = iRow;
          }
          for (iColumn = 0; iColumn < numberColumns_; iColumn++)
               whichColumn[iColumn] = -1;
          for (j = rowStart[knapsackRow]; j < rowStart[knapsackRow] + rowLength[knapsackRow]; j++) {
               int iColumn = column[j];
               if (columnUpper_[iColumn] > columnLower_[iColumn]) {
                    whichColumn[iColumn] = 0;
               } else {
                    assert (!columnLower_[iColumn]); // fix later
               }
          }
          for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
               if (whichColumn[iColumn] < 0)
                    whichColumn[nCol++] = iColumn;
          }
          ClpSimplex tempModel2(&tempModel, numberRows_, whichRow, nCol, whichColumn, false, false, false);
          // Row copy
          CoinPackedMatrix matrixByRow;
          matrixByRow.reverseOrderedCopyOf(*tempModel2.matrix());
          const double * elementByRow = matrixByRow.getElements();
          const int * column = matrixByRow.getIndices();
          const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
          const int * rowLength = matrixByRow.getVectorLengths();
          const double * columnLower = tempModel2.getColLower();
          const double * columnUpper = tempModel2.getColUpper();
          for (iRow = 0; iRow < numberRows_; iRow++) {
               lo[iRow] = -COIN_DBL_MAX;
               high[iRow] = COIN_DBL_MAX;
               if (rowLower_[iRow] > -1.0e20 || rowUpper_[iRow] < 1.0e20) {

                    // possible row
                    int infiniteUpper = 0;
                    int infiniteLower = 0;
                    double maximumUp = 0.0;
                    double maximumDown = 0.0;
                    CoinBigIndex rStart = rowStart[iRow];
                    CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
                    CoinBigIndex j;
                    // Compute possible lower and upper ranges

                    for (j = rStart; j < rEnd; ++j) {
                         double value = elementByRow[j];
                         iColumn = column[j];
                         if (value > 0.0) {
                              if (columnUpper[iColumn] >= 1.0e20) {
                                   ++infiniteUpper;
                              } else {
                                   maximumUp += columnUpper[iColumn] * value;
                              }
                              if (columnLower[iColumn] <= -1.0e20) {
                                   ++infiniteLower;
                              } else {
                                   maximumDown += columnLower[iColumn] * value;
                              }
                         } else if (value < 0.0) {
                              if (columnUpper[iColumn] >= 1.0e20) {
                                   ++infiniteLower;
                              } else {
                                   maximumDown += columnUpper[iColumn] * value;
                              }
                              if (columnLower[iColumn] <= -1.0e20) {
                                   ++infiniteUpper;
                              } else {
                                   maximumUp += columnLower[iColumn] * value;
                              }
                         }
                    }
                    // Build in a margin of error
                    maximumUp += 1.0e-8 * fabs(maximumUp) + 1.0e-7;
                    maximumDown -= 1.0e-8 * fabs(maximumDown) + 1.0e-7;
                    // we want to save effective rhs
                    double up = (infiniteUpper) ? COIN_DBL_MAX : maximumUp;
                    double down = (infiniteLower) ? -COIN_DBL_MAX : maximumDown;
                    if (up == COIN_DBL_MAX || rowLower_[iRow] == -COIN_DBL_MAX) {
                         // However low we go it doesn't matter
                         lo[iRow] = -COIN_DBL_MAX;
                    } else {
                         // If we go below this then can not be feasible
                         lo[iRow] = rowLower_[iRow] - up;
                    }
                    if (down == -COIN_DBL_MAX || rowUpper_[iRow] == COIN_DBL_MAX) {
                         // However high we go it doesn't matter
                         high[iRow] = COIN_DBL_MAX;
                    } else {
                         // If we go above this then can not be feasible
                         high[iRow] = rowUpper_[iRow] - down;
                    }
               }
          }
     }
     numJ = 0;
     for (iColumn = 0; iColumn < numberColumns_; iColumn++)
          whichColumn[iColumn] = -1;
     int * markRow = new int [numberRows_];
     for (iRow = 0; iRow < numberRows_; iRow++)
          markRow[iRow] = 1;
     for (j = rowStart[knapsackRow]; j < rowStart[knapsackRow] + rowLength[knapsackRow]; j++) {
          int iColumn = column[j];
          if (columnUpper_[iColumn] > columnLower_[iColumn]) {
               whichColumn[iColumn] = numJ;
               numJ++;
          }
     }
     /* mark rows
        -n in knapsack and n other variables
        1 no entries
        n+1000 not involved in knapsack but n entries
        0 only in knapsack
     */
     for (iRow = 0; iRow < numberRows_; iRow++) {
          int type = 1;
          for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
               int iColumn = column[j];
               if (whichColumn[iColumn] >= 0) {
                    if (type == 1) {
                         type = 0;
                    } else if (type > 0) {
                         assert (type > 1000);
                         type = -(type - 1000);
                    }
               } else if (type == 1) {
                    type = 1001;
               } else if (type < 0) {
                    type --;
               } else if (type == 0) {
                    type = -1;
               } else {
                    assert (type > 1000);
                    type++;
               }
          }
          markRow[iRow] = type;
          if (type < 0 && type > -30 && false)
               printf("markrow on row %d is %d\n", iRow, markRow[iRow]);
     }
     int * bound = new int [numberColumns_+1];
     int * stack = new int [numberColumns_+1];
     int * flip = new int [numberColumns_+1];
     double * offset = new double[numberColumns_+1];
     double * size = new double [numberColumns_+1];
     double * rhsOffset = new double[numberRows_];
     int * build = new int[numberColumns_];
     int maxNumber = numberOutput;
     numJ = 0;
     double minSize = rowLower_[knapsackRow];
     double maxSize = rowUpper_[knapsackRow];
     double knapsackOffset = 0.0;
     for (j = rowStart[knapsackRow]; j < rowStart[knapsackRow] + rowLength[knapsackRow]; j++) {
          int iColumn = column[j];
          double lowerColumn = columnLower_[iColumn];
          double upperColumn = columnUpper_[iColumn];
          if (lowerColumn == upperColumn)
               continue;
          double gap = upperColumn - lowerColumn;
          if (gap > 1.0e8)
               gap = 1.0e8;
          assert (fabs(floor(gap + 0.5) - gap) < 1.0e-5);
          whichColumn[numJ] = iColumn;
          bound[numJ] = static_cast<int> (gap);
          if (elementByRow[j] > 0.0) {
               flip[numJ] = 1;
               offset[numJ] = lowerColumn;
               size[numJ++] = elementByRow[j];
          } else {
               flip[numJ] = -1;
               offset[numJ] = upperColumn;
               size[numJ++] = -elementByRow[j];
               lowerColumn = upperColumn;
          }
          knapsackOffset += elementByRow[j] * lowerColumn;
     }
     int jRow;
     for (iRow = 0; iRow < numberRows_; iRow++)
          whichRow[iRow] = iRow;
     ClpSimplex smallModel(this, numberRows_, whichRow, numJ, whichColumn, true, true, true);
     // modify rhs to allow for nonzero lower bounds
     //double * rowLower = smallModel.rowLower();
     //double * rowUpper = smallModel.rowUpper();
     //const double * columnLower = smallModel.columnLower();
     //const double * columnUpper = smallModel.columnUpper();
     const CoinPackedMatrix * matrix = smallModel.matrix();
     const double * element = matrix->getElements();
     const int * row = matrix->getIndices();
     const CoinBigIndex * columnStart = matrix->getVectorStarts();
     const int * columnLength = matrix->getVectorLengths();
     const double * objective = smallModel.objective();
     //double objectiveOffset=0.0;
     // would use for fixed?
     CoinZeroN(rhsOffset, numberRows_);
     double * rowActivity = smallModel.primalRowSolution();
     CoinZeroN(rowActivity, numberRows_);
     maxSize -= knapsackOffset;
     minSize -= knapsackOffset;
     // now generate
     int i;
     int iStack = numJ;
     for (i = 0; i < numJ; i++) {
          stack[i] = 0;
     }
     double tooMuch = 10.0 * maxSize + 10000;
     stack[numJ] = 1;
     size[numJ] = tooMuch;
     bound[numJ] = 0;
     double sum = tooMuch;
     // allow for all zero being OK
     stack[numJ-1] = -1;
     sum -= size[numJ-1];
     numberOutput = 0;
     int nelCreate = 0;
     /* typeRun is - 0 for initial sizes
                     1 for build
     	  2 for reconstruct
     */
     int typeRun = buildObj ? 1 : 0;
     if (reConstruct >= 0) {
          assert (buildRow && buildElement);
          typeRun = 2;
     }
     if (typeRun == 1)
          buildStart[0] = 0;
     while (iStack >= 0) {
          if (sum >= minSize && sum <= maxSize) {
               double checkSize = 0.0;
               bool good = true;
               int nRow = 0;
               double obj = 0.0;
               CoinZeroN(rowActivity, numberRows_);
               for (iColumn = 0; iColumn < numJ; iColumn++) {
                    int iValue = stack[iColumn];
                    if (iValue > bound[iColumn]) {
                         good = false;
                         break;
                    } else {
                         double realValue = offset[iColumn] + flip[iColumn] * iValue;
                         if (realValue) {
                              obj += objective[iColumn] * realValue;
                              for (CoinBigIndex j = columnStart[iColumn];
                                        j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                                   double value = element[j] * realValue;
                                   int kRow = row[j];
                                   if (rowActivity[kRow]) {
                                        rowActivity[kRow] += value;
                                        if (!rowActivity[kRow])
                                             rowActivity[kRow] = 1.0e-100;
                                   } else {
                                        build[nRow++] = kRow;
                                        rowActivity[kRow] = value;
                                   }
                              }
                         }
                    }
               }
               if (good) {
                    for (jRow = 0; jRow < nRow; jRow++) {
                         int kRow = build[jRow];
                         double value = rowActivity[kRow];
                         if (value > high[kRow] || value < lo[kRow]) {
                              good = false;
                              break;
                         }
                    }
               }
               if (good) {
                    if (typeRun == 1) {
                         buildObj[numberOutput] = obj;
                         for (jRow = 0; jRow < nRow; jRow++) {
                              int kRow = build[jRow];
                              double value = rowActivity[kRow];
                              if (markRow[kRow] < 0 && fabs(value) > 1.0e-13) {
                                   buildElement[nelCreate] = value;
                                   buildRow[nelCreate++] = kRow;
                              }
                         }
                         buildStart[numberOutput+1] = nelCreate;
                    } else if (!typeRun) {
                         for (jRow = 0; jRow < nRow; jRow++) {
                              int kRow = build[jRow];
                              double value = rowActivity[kRow];
                              if (markRow[kRow] < 0 && fabs(value) > 1.0e-13) {
                                   nelCreate++;
                              }
                         }
                    }
                    if (typeRun == 2 && reConstruct == numberOutput) {
                         // build and exit
                         nelCreate = 0;
                         for (iColumn = 0; iColumn < numJ; iColumn++) {
                              int iValue = stack[iColumn];
                              double realValue = offset[iColumn] + flip[iColumn] * iValue;
                              if (realValue) {
                                   buildRow[nelCreate] = whichColumn[iColumn];
                                   buildElement[nelCreate++] = realValue;
                              }
                         }
                         numberOutput = 1;
                         for (i = 0; i < numJ; i++) {
                              bound[i] = 0;
                         }
                         break;
                    }
                    numberOutput++;
                    if (numberOutput > maxNumber) {
                         nelCreate = -numberOutput;
                         numberOutput = -1;
                         for (i = 0; i < numJ; i++) {
                              bound[i] = 0;
                         }
                         break;
                    } else if (typeRun == 1 && numberOutput == maxNumber) {
                         // On second run
                         for (i = 0; i < numJ; i++) {
                              bound[i] = 0;
                         }
                         break;
                    }
                    for (int j = 0; j < numJ; j++) {
                         checkSize += stack[j] * size[j];
                    }
                    assert (fabs(sum - checkSize) < 1.0e-3);
               }
               for (jRow = 0; jRow < nRow; jRow++) {
                    int kRow = build[jRow];
                    rowActivity[kRow] = 0.0;
               }
          }
          if (sum > maxSize || stack[iStack] > bound[iStack]) {
               sum -= size[iStack] * stack[iStack];
               stack[iStack--] = 0;
               if (iStack >= 0) {
                    stack[iStack] ++;
                    sum += size[iStack];
               }
          } else {
               // must be less
               // add to last possible
               iStack = numJ - 1;
               sum += size[iStack];
               stack[iStack]++;
          }
     }
     //printf("%d will be created\n",numberOutput);
     delete [] whichColumn;
     delete [] whichRow;
     delete [] bound;
     delete [] stack;
     delete [] flip;
     delete [] size;
     delete [] offset;
     delete [] rhsOffset;
     delete [] build;
     delete [] markRow;
     delete [] lo;
     delete [] high;
     return nelCreate;
}
// Quick try at cleaning up duals if postsolve gets wrong
void
ClpSimplexOther::cleanupAfterPostsolve()
{
     // First mark singleton equality rows
     char * mark = new char [ numberRows_];
     memset(mark, 0, numberRows_);
     const int * row = matrix_->getIndices();
     const CoinBigIndex * columnStart = matrix_->getVectorStarts();
     const int * columnLength = matrix_->getVectorLengths();
     const double * element = matrix_->getElements();
     for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
          for (CoinBigIndex j = columnStart[iColumn];
                    j < columnStart[iColumn] + columnLength[iColumn]; j++) {
               int iRow = row[j];
               if (mark[iRow])
                    mark[iRow] = 2;
               else
                    mark[iRow] = 1;
          }
     }
     // for now just == rows
     for (int iRow = 0; iRow < numberRows_; iRow++) {
          if (rowUpper_[iRow] > rowLower_[iRow])
               mark[iRow] = 3;
     }
     double dualTolerance = dblParam_[ClpDualTolerance];
     double primalTolerance = dblParam_[ClpPrimalTolerance];
     int numberCleaned = 0;
     double maxmin = optimizationDirection_;
     for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double dualValue = reducedCost_[iColumn] * maxmin;
          double primalValue = columnActivity_[iColumn];
          double lower = columnLower_[iColumn];
          double upper = columnUpper_[iColumn];
          int way = 0;
          switch(getColumnStatus(iColumn)) {

          case basic:
               // dual should be zero
               if (dualValue > dualTolerance) {
                    way = -1;
               } else if (dualValue < -dualTolerance) {
                    way = 1;
               }
               break;
          case ClpSimplex::isFixed:
               break;
          case atUpperBound:
               // dual should not be positive
               if (dualValue > dualTolerance) {
                    way = -1;
               }
               break;
          case atLowerBound:
               // dual should not be negative
               if (dualValue < -dualTolerance) {
                    way = 1;
               }
               break;
          case superBasic:
          case isFree:
               if (primalValue < upper - primalTolerance) {
                    // dual should not be negative
                    if (dualValue < -dualTolerance) {
                         way = 1;
                    }
               }
               if (primalValue > lower + primalTolerance) {
                    // dual should not be positive
                    if (dualValue > dualTolerance) {
                         way = -1;
                    }
               }
               break;
          }
          if (way) {
               // see if can find singleton row
               for (CoinBigIndex j = columnStart[iColumn];
                         j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                    int iRow = row[j];
                    if (mark[iRow] == 1) {
                         double value = element[j];
                         // dj - addDual*value == 0.0
                         double addDual = dualValue / value;
                         dual_[iRow] += addDual;
                         reducedCost_[iColumn] = 0.0;
                         numberCleaned++;
                         break;
                    }
               }
          }
     }
     delete [] mark;
#ifdef CLP_INVESTIGATE
     printf("cleanupAfterPostsolve cleaned up %d columns\n", numberCleaned);
#endif
     // Redo
     memcpy(reducedCost_, this->objective(), numberColumns_ * sizeof(double));
     matrix_->transposeTimes(-1.0, dual_, reducedCost_);
     checkSolutionInternal();
}
// Returns gub version of model or NULL
ClpSimplex * 
ClpSimplexOther::gubVersion(int * whichRows, int * whichColumns,
			    int neededGub,
			    int factorizationFrequency)
{
  // find gub
  int numberRows = this->numberRows();
  int numberColumns = this->numberColumns();
  int iRow, iColumn;
  int * columnIsGub = new int [numberColumns];
  const double * columnLower = this->columnLower();
  const double * columnUpper = this->columnUpper();
  int numberFixed=0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    if (columnUpper[iColumn] == columnLower[iColumn]) {
      columnIsGub[iColumn]=-2;
      numberFixed++;
    } else if (columnLower[iColumn]>=0) {
      columnIsGub[iColumn]=-1;
    } else {
      columnIsGub[iColumn]=-3;
    }
  }
  CoinPackedMatrix * matrix = this->matrix();
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int * column = rowCopy.getIndices();
  const int * rowLength = rowCopy.getVectorLengths();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const double * element = rowCopy.getElements();
  int numberNonGub = 0;
  int numberEmpty = numberRows;
  int * rowIsGub = new int [numberRows];
  int smallestGubRow=-1;
  int count=numberColumns+1;
  double * rowLower = this->rowLower();
  double * rowUpper = this->rowUpper();
  // make sure we can get rid of upper bounds
  double * fixedRow = new double [numberRows];
  for (iRow = 0 ; iRow < numberRows ; iRow++) {
    double sumFixed=0.0;
    for (int j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      double value = columnLower[iColumn];
      if (value) 
	sumFixed += element[j] * value;
    }
    fixedRow[iRow]=rowUpper[iRow]-sumFixed;
  }
  for (iRow = numberRows - 1; iRow >= 0; iRow--) {
    bool gubRow = true;
    int numberInRow=0;
    double sumFixed=0.0;
    double gap = fixedRow[iRow]-1.0e-12;
    for (int j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      if (columnIsGub[iColumn]!=-2) {
	if (element[j] != 1.0||columnIsGub[iColumn]==-3||
	    columnUpper[iColumn]-columnLower[iColumn]<gap) {
	  gubRow = false;
	  break;
	} else {
	  numberInRow++;
	  if (columnIsGub[iColumn] >= 0) {
	    gubRow = false;
	    break;
	  }
	}
      } else {
	sumFixed += columnLower[iColumn]*element[j];
      }
    }
    if (!gubRow) {
      whichRows[numberNonGub++] = iRow;
      rowIsGub[iRow] = -1;
    } else if (numberInRow) {
      if (numberInRow<count) {
	count = numberInRow;
	smallestGubRow=iRow;
      }
      for (int j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
	int iColumn = column[j];
	if (columnIsGub[iColumn]!=-2) 
	  columnIsGub[iColumn] = iRow;
      }
      rowIsGub[iRow] = 0;
    } else {
      // empty row!
      whichRows[--numberEmpty] = iRow;
      rowIsGub[iRow] = -2;
      if (sumFixed>rowUpper[iRow]+1.0e-4||
	  sumFixed<rowLower[iRow]-1.0e-4) {
	fprintf(stderr,"******** No infeasible empty rows - please!\n");
	abort();
      }
    }
  }
  delete [] fixedRow;
  char message[100];
  int numberGub = numberEmpty - numberNonGub;
  if (numberGub >= neededGub) { 
    sprintf(message,"%d gub rows", numberGub);
    handler_->message(CLP_GENERAL2, messages_)
      << message << CoinMessageEol;
    int numberNormal = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnIsGub[iColumn] < 0  && columnIsGub[iColumn] !=-2) {
	whichColumns[numberNormal++] = iColumn;
      }
    }
    if (!numberNormal) {
      sprintf(message,"Putting back one gub row to make non-empty");
      handler_->message(CLP_GENERAL2, messages_)
	<< message << CoinMessageEol;
      rowIsGub[smallestGubRow]=-1;
      whichRows[numberNonGub++] = smallestGubRow;
      for (int j = rowStart[smallestGubRow]; 
	   j < rowStart[smallestGubRow] + rowLength[smallestGubRow]; j++) {
	int iColumn = column[j];
	if (columnIsGub[iColumn]>=0) {
	  columnIsGub[iColumn]=-4;
	  whichColumns[numberNormal++] = iColumn;
	}
      }
    }
    std::sort(whichRows,whichRows+numberNonGub);
    std::sort(whichColumns,whichColumns+numberNormal);
    double * lower = CoinCopyOfArray(this->rowLower(),numberRows);
    double * upper = CoinCopyOfArray(this->rowUpper(),numberRows);
    // leave empty rows at end
    numberEmpty = numberRows-numberEmpty;
    const int * row = matrix->getIndices();
    const int * columnLength = matrix->getVectorLengths();
    const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const double * elementByColumn = matrix->getElements();
    // Fixed at end
    int put2 = numberColumns-numberFixed;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnIsGub[iColumn] ==-2) {
	whichColumns[put2++] = iColumn;
	double value = columnLower[iColumn];
	for (int j = columnStart[iColumn]; 
	     j < columnStart[iColumn] + columnLength[iColumn]; j++) {
	  int iRow = row[j];
	  if (lower[iRow]>-1.0e20)
	    lower[iRow] -= value*element[j];
	  if (upper[iRow]<1.0e20)
	    upper[iRow] -= value*element[j];
	}
      }
    }
    int put = numberNormal;
    ClpSimplex * model2 =
      new ClpSimplex(this, numberNonGub, whichRows , numberNormal, whichColumns);
    // scale
    double * scaleArray = new double [numberRows];
    for (int i=0;i<numberRows;i++) {
      scaleArray[i]=1.0;
      if (rowIsGub[i]==-1) {
	double largest = 1.0e-30;
	double smallest = 1.0e30;
	for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
	  int iColumn = column[j];
	  if (columnIsGub[iColumn]!=-2) {
	    double value =fabs(element[j]);
	    largest = CoinMax(value,largest);
	    smallest = CoinMin(value,smallest);
	  }
	}
	double scale = CoinMax(0.001,1.0/sqrt(largest*smallest));
	scaleArray[i]=scale;
	if (lower[i]>-1.0e30)
	  lower[i] *= scale;
	if (upper[i]<1.0e30)
	  upper[i] *= scale;
      }
    }
    // scale partial matrix
    {
      CoinPackedMatrix * matrix = model2->matrix();
      const int * row = matrix->getIndices();
      const int * columnLength = matrix->getVectorLengths();
      const CoinBigIndex * columnStart = matrix->getVectorStarts();
      double * element = matrix->getMutableElements();
      for (int i=0;i<numberNormal;i++) {
	for (int j = columnStart[i]; 
	     j < columnStart[i] + columnLength[i]; j++) {
	  int iRow = row[j];
	  iRow = whichRows[iRow];
	  double scaleBy = scaleArray[iRow];
	  element[j] *= scaleBy;
	}
      }
    }
    // adjust rhs
    double * rowLower = model2->rowLower();
    double * rowUpper = model2->rowUpper();
    for (int i=0;i<numberNonGub;i++) {
      int iRow = whichRows[i];
      rowLower[i] = lower[iRow];
      rowUpper[i] = upper[iRow];
    }
    int numberGubColumns = numberColumns - put - numberFixed;
    CoinBigIndex numberElements=0;
    int * temp1 = new int [numberRows+1];
    // get counts
    memset(temp1,0,numberRows*sizeof(int));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int iGub = columnIsGub[iColumn];
      if (iGub>=0) {
	numberElements += columnLength[iColumn]-1;
	temp1[iGub]++;
      }
    }
    /* Optional but means can eventually simplify coding
       could even add in fixed slacks to deal with
       singularities - but should not be necessary */
    int numberSlacks=0;
    for (int i = 0; i < numberRows; i++) {
      if (rowIsGub[i]>=0) {
	if (lower[i]<upper[i]) {
	  numberSlacks++;
	  temp1[i]++;
	}
      }
    }
    int * gubStart = new int [numberGub+1];
    numberGub=0;
    gubStart[0]=0;
    for (int i = 0; i < numberRows; i++) {
      if (rowIsGub[i]>=0) {
	rowIsGub[i]=numberGub;
	gubStart[numberGub+1]=gubStart[numberGub]+temp1[i];
	temp1[numberGub]=0;
	lower[numberGub]=lower[i];
	upper[numberGub]=upper[i];
	whichRows[numberNonGub+numberGub]=i;
	numberGub++;
      }
    }
    int numberGubColumnsPlus = numberGubColumns + numberSlacks;
    double * lowerColumn2 = new double [numberGubColumnsPlus];
    CoinFillN(lowerColumn2, numberGubColumnsPlus, 0.0);
    double * upperColumn2 = new double [numberGubColumnsPlus];
    CoinFillN(upperColumn2, numberGubColumnsPlus, COIN_DBL_MAX);
    int * start2 = new int[numberGubColumnsPlus+1];
    int * row2 = new int[numberElements];
    double * element2 = new double[numberElements];
    double * cost2 = new double [numberGubColumnsPlus];
    CoinFillN(cost2, numberGubColumnsPlus, 0.0);
    const double * cost = this->objective();
    put = numberNormal;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int iGub = columnIsGub[iColumn];
      if (iGub>=0) {
	// TEMP
	//this->setColUpper(iColumn,COIN_DBL_MAX);
	iGub = rowIsGub[iGub];
	assert (iGub>=0);
	int kPut = put+gubStart[iGub]+temp1[iGub];
	temp1[iGub]++;
	whichColumns[kPut]=iColumn;
      }
    }
    for (int i = 0; i < numberRows; i++) {
      if (rowIsGub[i]>=0) {
	int iGub = rowIsGub[i];
	if (lower[iGub]<upper[iGub]) {
	  int kPut = put+gubStart[iGub]+temp1[iGub];
	  temp1[iGub]++;
	  whichColumns[kPut]=iGub+numberColumns;
	}
      }
    }
    //this->primal(1); // TEMP
    // redo rowIsGub to give lookup
    for (int i=0;i<numberRows;i++)
      rowIsGub[i]=-1;
    for (int i=0;i<numberNonGub;i++) 
      rowIsGub[whichRows[i]]=i;
    start2[0]=0;
    numberElements = 0;
    for (int i=0;i<numberGubColumnsPlus;i++) {
      int iColumn = whichColumns[put++];
      if (iColumn<numberColumns) {
	cost2[i] = cost[iColumn];
	lowerColumn2[i] = columnLower[iColumn];
	upperColumn2[i] = columnUpper[iColumn];
	upperColumn2[i] = COIN_DBL_MAX; 
	for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
	  int iRow = row[j];
	  double scaleBy = scaleArray[iRow];
	  iRow = rowIsGub[iRow];
	  if (iRow >= 0) {
	    row2[numberElements] = iRow;
	    element2[numberElements++] = elementByColumn[j]*scaleBy;
	  }
	}
      } else {
	// slack
	int iGub = iColumn-numberColumns;
	double slack = upper[iGub]-lower[iGub];
	assert (upper[iGub]<1.0e20);
	lower[iGub]=upper[iGub];
	cost2[i] = 0;
	lowerColumn2[i] = 0;
	upperColumn2[i] = slack;
	upperColumn2[i] = COIN_DBL_MAX;
      }
      start2[i+1] = numberElements;
    }
    // clean up bounds on variables
    for (int iSet=0;iSet<numberGub;iSet++) {
      double lowerValue=0.0;
      for (int i=gubStart[iSet];i<gubStart[iSet+1];i++) {
	lowerValue += lowerColumn2[i];
      }
      assert (lowerValue<upper[iSet]+1.0e-6);
      double gap = CoinMax(0.0,upper[iSet]-lowerValue);
      for (int i=gubStart[iSet];i<gubStart[iSet+1];i++) {
	if (upperColumn2[i]<1.0e30) {
	  upperColumn2[i] = CoinMin(upperColumn2[i],
				    lowerColumn2[i]+gap);
	}
      }
    }
    sprintf(message,"** Before adding matrix there are %d rows and %d columns",
	   model2->numberRows(), model2->numberColumns());
    handler_->message(CLP_GENERAL2, messages_)
      << message << CoinMessageEol;
    delete [] scaleArray;
    delete [] temp1;
    model2->setFactorizationFrequency(factorizationFrequency);
    ClpDynamicMatrix * newMatrix = 
      new ClpDynamicMatrix(model2, numberGub,
				  numberGubColumnsPlus, gubStart,
				  lower, upper,
				  start2, row2, element2, cost2,
				  lowerColumn2, upperColumn2);
    delete [] gubStart;
    delete [] lowerColumn2;
    delete [] upperColumn2;
    delete [] start2;
    delete [] row2;
    delete [] element2;
    delete [] cost2;
    delete [] lower;
    delete [] upper;
    model2->replaceMatrix(newMatrix,true);
#ifdef EVERY_ITERATION
    {
      ClpDynamicMatrix * gubMatrix =
	dynamic_cast< ClpDynamicMatrix*>(model2->clpMatrix());
      assert(gubMatrix);
      gubMatrix->writeMps("gub.mps");
    }
#endif
    delete [] columnIsGub;
    delete [] rowIsGub;
    newMatrix->switchOffCheck();
#ifdef EVERY_ITERATION
    newMatrix->setRefreshFrequency(1/*000*/);
#else
    newMatrix->setRefreshFrequency(1000);
#endif
    sprintf(message,
	    "** While after adding matrix there are %d rows and %d columns",
	    model2->numberRows(), model2->numberColumns());
    handler_->message(CLP_GENERAL2, messages_)
      << message << CoinMessageEol;
    model2->setSpecialOptions(4);    // exactly to bound
    // Scaling off (done by hand)
    model2->scaling(0);
    return model2;
  } else {
    delete [] columnIsGub;
    delete [] rowIsGub;
    return NULL;
  }
}
// Sets basis from original
void 
ClpSimplexOther::setGubBasis(ClpSimplex &original,const int * whichRows,
			     const int * whichColumns)
{
  ClpDynamicMatrix * gubMatrix =
    dynamic_cast< ClpDynamicMatrix*>(clpMatrix());
  assert(gubMatrix);
  int numberGubColumns = gubMatrix->numberGubColumns();
  int numberNormal = gubMatrix->firstDynamic();
  //int lastOdd = gubMatrix->firstAvailable();
  //int numberTotalColumns = numberNormal + numberGubColumns;
  //assert (numberTotalColumns==numberColumns+numberSlacks);
  int numberRows = original.numberRows();
  int numberColumns = original.numberColumns();
  int * columnIsGub = new int [numberColumns];
  int numberNonGub = gubMatrix->numberStaticRows();
  //assert (firstOdd==numberNormal);
  double * solution = primalColumnSolution();
  double * originalSolution = original.primalColumnSolution();
  const double * upperSet = gubMatrix->upperSet();
  // Column copy of GUB part
  int numberSets = gubMatrix->numberSets();
  const int * startSet = gubMatrix->startSets();
  const CoinBigIndex * columnStart = gubMatrix->startColumn();
  const double * columnLower = gubMatrix->columnLower();
#ifdef TRY_IMPROVE
  const double * columnUpper = gubMatrix->columnUpper();
  const double * lowerSet = gubMatrix->lowerSet();
  const double * element = gubMatrix->element();
  const int * row = gubMatrix->row();
  bool allPositive=true;
  double * rowActivity = new double[numberNonGub];
  memset(rowActivity, 0, numberNonGub*sizeof(double));
  {
    // Non gub contribution
    const double * element = matrix_->getElements();
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths();
    for (int i=0;i<numberNormal;i++) {
      int iColumn = whichColumns[i];
      double value = originalSolution[iColumn];
      if (value) {
	for (CoinBigIndex j = columnStart[i];
	     j < columnStart[i] + columnLength[i]; j++) {
	  int iRow = row[j];
	  rowActivity[iRow] += value*element[j];
	}
      }
    }
  }
  double * newSolution = new double [numberGubColumns];
  int * slacks = new int [numberSets];
  for (int i=0;i<numberSets;i++) {
    double sum=0.0;
    int iSlack=-1;
    for (int j=startSet[i];j<startSet[i+1];j++) {
      gubMatrix->setDynamicStatus(j,ClpDynamicMatrix::atLowerBound);
      int iColumn = whichColumns[j+numberNormal];
      if (iColumn<numberColumns) {
	columnIsGub[iColumn] = whichRows[numberNonGub+i];
	double value = originalSolution[iColumn];
	sum += value;
	newSolution[j]=value;
	for (CoinBigIndex k = columnStart[j]; k < columnStart[j+1] ; k++) {
	  int iRow = row[k];
	  rowActivity[iRow] += value*element[k];
	  if (element[k] < 0.0)
	    allPositive=false;
	}
	if (columnStart[j]==columnStart[j+1])
	  iSlack=j;
      } else {
	newSolution[j]=0.0;
	iSlack=j;
	allPositive=false; // for now
      }
    }
    slacks[i]=iSlack;
    if (sum>upperSet[i]+1.0e-8) {
      double gap = sum-upperSet[i];
      if (iSlack>=0) {
	double value=newSolution[iSlack];
	if (value>0.0) {
	  double down = CoinMin(gap,value);
	  gap -= down;
	  sum -= down;
	  newSolution[iSlack] = value-down;
	}
      }
      if (gap>1.0e-8) {
	for (int j=startSet[i];j<startSet[i+1];j++) {
	  int iColumn = whichColumns[j+numberNormal];
	  if (newSolution[j]>0.0&&iColumn<numberColumns) {
	    double value = newSolution[j];
	    double down = CoinMin(gap,value);
	    gap -= down;
	    sum -= down;
	    newSolution[iSlack] = value-down;
	    for (CoinBigIndex k = columnStart[j]; k < columnStart[j+1] ; k++) {
	      int iRow = row[k];
	      rowActivity[iRow] -= down*element[k];
	    }
	  }
	}
      }
      assert (gap<1.0e-8);
    } else if (sum<lowerSet[i]-1.0e-8) {
      double gap = lowerSet[i]-sum;
      if (iSlack>=0) {
	double value=newSolution[iSlack];
	if (value<columnUpper[iSlack]) {
	  double up = CoinMin(gap,columnUpper[iSlack]-value);
	  gap -= up;
	  sum += up;
	  newSolution[iSlack] = value+up;
	}
      }
      if (gap>1.0e-8) {
	for (int j=startSet[i];j<startSet[i+1];j++) {
	  int iColumn = whichColumns[j+numberNormal];
	  if (newSolution[j]<columnUpper[j]&&iColumn<numberColumns) {
	    double value = newSolution[j];
	    double up = CoinMin(gap,columnUpper[j]-value);
	    gap -= up;
	    sum += up;
	    newSolution[iSlack] = value+up;
	    for (CoinBigIndex k = columnStart[j]; k < columnStart[j+1] ; k++) {
	      int iRow = row[k];
	      rowActivity[iRow] += up*element[k];
	    }
	  }
	}
      }
      assert (gap<1.0e-8);
    } 
    if (fabs(sum-upperSet[i])>1.0e-7)
      printf("Sum for set %d is %g - lower %g, upper %g\n",i,
	     sum,lowerSet[i],upperSet[i]);
  }
  if (allPositive) {
    // See if we can improve solution
    // first reduce if over
    double * gaps = new double [numberNonGub];
    double direction = optimizationDirection_; 
    const double * cost = gubMatrix->cost();
    bool over=false;
    for (int i=0;i<numberNonGub;i++) {
      double activity = rowActivity[i];
      gaps[i]=0.0;
      if (activity>rowUpper_[i]+1.0e-6) {
	gaps[i]=activity-rowUpper_[i];
	over=true;
      }
    }
    double * weights = new double [numberGubColumns];
    int * which = new int [numberGubColumns];
    int * whichSet = new int [numberGubColumns];
    if (over) {
      int n=0;
      for (int i=0;i<numberSets;i++) {
	int iSlack = slacks[i];
	if (iSlack<0||newSolution[iSlack]>upperSet[i]-1.0e-8)
	  continue;
	double slackCost = cost[iSlack]*direction;
	for (int j=startSet[i];j<startSet[i+1];j++) {
	  whichSet[j]=i;
	  double value = newSolution[j];
	  double thisCost = cost[j]*direction;
	  if (value>columnLower[j]&&j!=iSlack) {
	    if(thisCost<slackCost) {
	      double sum = 1.0e-30;
	      for (CoinBigIndex k = columnStart[j]; 
		   k < columnStart[j+1] ; k++) {
		int iRow = row[k];
		sum += gaps[iRow]*element[k];
	      }
	      which[n]=j;
	      // big drop and small difference in cost better
	      weights[n++]=(slackCost-thisCost)/sum;
	    } else {
	      // slack better anyway
	      double move = value-columnLower[j];
	      newSolution[iSlack]=CoinMin(upperSet[i],
					  newSolution[iSlack]+move);
	      newSolution[j]=columnLower[j];
	      for (CoinBigIndex k = columnStart[j]; 
		   k < columnStart[j+1] ; k++) {
		int iRow = row[k];
		rowActivity[iRow] -= move*element[k];
	      }
	    }
	  }
	}
      }
      // sort
      CoinSort_2(weights,weights+n,which);
      for (int i=0;i<n;i++) {
	int j= which[i];
	int iSet = whichSet[j];
	int iSlack = slacks[iSet];
	assert (iSlack>=0);
	double move = 0.0;
	for (CoinBigIndex k = columnStart[j]; 
	     k < columnStart[j+1] ; k++) {
	  int iRow = row[k];
	  if(rowActivity[iRow]-rowUpper_[iRow]>move*element[k]) {
	    move = (rowActivity[iRow]-rowUpper_[iRow])/element[k];
	  }
	}
	move=CoinMin(move,newSolution[j]-columnLower[j]);
	if (move) {
	  newSolution[j] -= move;
	  newSolution[iSlack] += move;
	  for (CoinBigIndex k = columnStart[j]; 
	       k < columnStart[j+1] ; k++) {
	    int iRow = row[k];
	    rowActivity[iRow] -= move*element[k];
	  }
	}
      }
    }
    delete [] whichSet;
    delete [] which;
    delete [] weights;
    delete [] gaps;
    // redo original status!
    for (int i=0;i<numberSets;i++) {
      int numberBasic=0;
      int numberNewBasic=0;
      int j1=-1;
      int j2=-1;
      for (int j=startSet[i];j<startSet[i+1];j++) {
	if (newSolution[j]>columnLower[j]) {
	  numberNewBasic++;
	  j2=j;
	}
	int iOrig = whichColumns[j+numberNormal];
	if (iOrig<numberColumns) {
	  if (original.getColumnStatus(iOrig)!=ClpSimplex::atLowerBound) {
	    numberBasic++;
	    j1=j;
	  }
	} else {
	  int iSet = iOrig - numberColumns;
	  int iRow = whichRows[iSet+numberNonGub];
	  if (original.getRowStatus(iRow)==ClpSimplex::basic) {
	    numberBasic++;
	    j1=j;
	    abort();
	  }
	}
      }
      if (numberBasic==1&&numberNewBasic==1&&
	  j1!=j2) {
	int iOrig1=whichColumns[j1+numberNormal];
	int iOrig2=whichColumns[j2+numberNormal];
	ClpSimplex::Status status1 = original.getColumnStatus(iOrig1);
	ClpSimplex::Status status2 = original.getColumnStatus(iOrig2);
	originalSolution[iOrig1] = newSolution[j1];
	originalSolution[iOrig2] = newSolution[j2];
	original.setColumnStatus(iOrig1,status2);
	original.setColumnStatus(iOrig2,status1);
      }
    }
  }
  delete [] newSolution;
  delete [] slacks;
  delete [] rowActivity;
#else
  for (int i=0;i<numberSets;i++) {
    for (int j=startSet[i];j<startSet[i+1];j++) {
      gubMatrix->setDynamicStatus(j,ClpDynamicMatrix::atLowerBound);
      int iColumn = whichColumns[j+numberNormal];
      if (iColumn<numberColumns) {
	columnIsGub[iColumn] = whichRows[numberNonGub+i];
      }
    }
  }
#endif
  int * numberKey = new int [numberRows];
  memset(numberKey,0,numberRows*sizeof(int));
  for (int i=0;i<numberGubColumns;i++) {
    int iOrig = whichColumns[i+numberNormal];
    if (iOrig<numberColumns) {
      if (original.getColumnStatus(iOrig)==ClpSimplex::basic) {
	int iRow = columnIsGub[iOrig];
	assert (iRow>=0);
	numberKey[iRow]++;
      }
    } else {
      // Set slack
      int iSet = iOrig - numberColumns;
      int iRow = whichRows[iSet+numberNonGub];
      if (original.getRowStatus(iRow)==ClpSimplex::basic) 
	numberKey[iRow]++;
    }
  }
  /* Before going into cleanMatrix we need
     gub status set (inSmall just means basic and active)
     row status set
  */
  for (int i = 0; i < numberSets; i++) {
    gubMatrix->setStatus(i,ClpSimplex::isFixed);
  }
  for (int i = 0; i < numberGubColumns; i++) {
    int iOrig = whichColumns[i+numberNormal];
    if (iOrig<numberColumns) {
      ClpSimplex::Status status = original.getColumnStatus(iOrig);
      if (status==ClpSimplex::atUpperBound) {
	gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::atUpperBound);
      } else if (status==ClpSimplex::atLowerBound) {
	gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::atLowerBound);
      } else if (status==ClpSimplex::basic) {
	int iRow = columnIsGub[iOrig];
	assert (iRow>=0);
	assert(numberKey[iRow]);
	if (numberKey[iRow]==1)
	  gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::soloKey);
	else
	  gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::inSmall);
      }
    } else {
      // slack
      int iSet = iOrig - numberColumns;
      int iRow = whichRows[iSet+numberNonGub];
      if (original.getRowStatus(iRow)==ClpSimplex::basic
#ifdef TRY_IMPROVE
	  ||newSolution[i]>columnLower[i]+1.0e-8
#endif
	  ) {
	assert(numberKey[iRow]);
	if (numberKey[iRow]==1)
	  gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::soloKey);
	else
	  gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::inSmall);
      } else {
	gubMatrix->setDynamicStatus(i,ClpDynamicMatrix::atLowerBound);
      }
    }
  }
  // deal with sets without key
  for (int i = 0; i < numberSets; i++) {
    int iRow = whichRows[numberNonGub+i];
    if (!numberKey[iRow]) {
      double upper = upperSet[i]-1.0e-7;
      if (original.getRowStatus(iRow)==ClpSimplex::basic) 
	gubMatrix->setStatus(i,ClpSimplex::basic);
      // If not at lb make key otherwise one with smallest number els
      double largest=0.0;
      int fewest=numberRows+1;
      int chosen=-1;
      for (int j=startSet[i];j<startSet[i+1];j++) {
	int length=columnStart[j+1]-columnStart[j];
	int iOrig = whichColumns[j+numberNormal];
	double value;
	if (iOrig<numberColumns) {
#ifdef TRY_IMPROVE
	  value=newSolution[j]-columnLower[j];
#else
	  value = originalSolution[iOrig]-columnLower[j];
#endif
	  if (value>upper)
	    gubMatrix->setStatus(i,ClpSimplex::atLowerBound);
	} else {
	  // slack - take value as 0.0 as will win on length
	  value=0.0;
	}
	if (value>largest+1.0e-8) {
	  largest=value;
	  fewest=length;
	  chosen=j;
	} else if (fabs(value-largest)<=1.0e-8&&length<fewest) {
	  largest=value;
	  fewest=length;
	  chosen=j;
	}
      }
      assert(chosen>=0);
      if (gubMatrix->getStatus(i)!=ClpSimplex::basic) {
	// set as key
	for (int j=startSet[i];j<startSet[i+1];j++) {
	  if (j!=chosen)
	    gubMatrix->setDynamicStatus(j,ClpDynamicMatrix::atLowerBound);
	  else
	    gubMatrix->setDynamicStatus(j,ClpDynamicMatrix::soloKey);
	}
      }
    }
  }
  for (int i = 0; i < numberNormal; i++) {
    int iOrig = whichColumns[i];
    setColumnStatus(i,original.getColumnStatus(iOrig));
    solution[i]=originalSolution[iOrig];
  }
  for (int i = 0; i < numberNonGub; i++) {
    int iOrig = whichRows[i];
    setRowStatus(i,original.getRowStatus(iOrig));
  }
  // Fill in current matrix
  gubMatrix->initialProblem();
  delete [] numberKey;
  delete [] columnIsGub;
}
// Restores basis to original
void 
ClpSimplexOther::getGubBasis(ClpSimplex &original,const int * whichRows,
			     const int * whichColumns) const
{
  ClpDynamicMatrix * gubMatrix =
    dynamic_cast< ClpDynamicMatrix*>(clpMatrix());
  assert(gubMatrix);
  int numberGubColumns = gubMatrix->numberGubColumns();
  int numberNormal = gubMatrix->firstDynamic();
  //int lastOdd = gubMatrix->firstAvailable();
  //int numberRows = original.numberRows();
  int numberColumns = original.numberColumns();
  int numberNonGub = gubMatrix->numberStaticRows();
  //assert (firstOdd==numberNormal);
  double * solution = primalColumnSolution();
  double * originalSolution = original.primalColumnSolution();
  int numberSets = gubMatrix->numberSets();
  const double * cost = original.objective();
  int lastOdd = gubMatrix->firstAvailable();
  //assert (numberTotalColumns==numberColumns+numberSlacks);
  int numberRows = original.numberRows();
  //int numberStaticRows = gubMatrix->numberStaticRows();
  const int * startSet = gubMatrix->startSets();
  unsigned char * status = original.statusArray();
  unsigned char * rowStatus = status+numberColumns;
  //assert (firstOdd==numberNormal);
  for (int i=0;i<numberSets;i++) {
    int iRow = whichRows[i+numberNonGub];
    original.setRowStatus(iRow,ClpSimplex::atLowerBound);
  }
  const int * id = gubMatrix->id();
  const double * columnLower = gubMatrix->columnLower();
  const double * columnUpper = gubMatrix->columnUpper();
  for (int i = 0; i < numberGubColumns; i++) {
    int iOrig = whichColumns[i+numberNormal];
    if (iOrig<numberColumns) {
      if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::atUpperBound) {
	originalSolution[iOrig] = columnUpper[i];
	status[iOrig] = 2;
      } else if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::atLowerBound && columnLower) {
	originalSolution[iOrig] = columnLower[i];
	status[iOrig] = 3;
      } else if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::soloKey) {
	int iSet = gubMatrix->whichSet(i);
	originalSolution[iOrig] = gubMatrix->keyValue(iSet);
	status[iOrig] = 1;
      } else {
	originalSolution[iOrig] = 0.0;
	status[iOrig] = 4;
      }
    } else {
      // slack
      int iSet = iOrig - numberColumns;
      int iRow = whichRows[iSet+numberNonGub];
      if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::atUpperBound) {
	original.setRowStatus(iRow,ClpSimplex::atLowerBound);
      } else if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::atLowerBound) {
	original.setRowStatus(iRow,ClpSimplex::atUpperBound);
      } else if (gubMatrix->getDynamicStatus(i) == ClpDynamicMatrix::soloKey) {
	original.setRowStatus(iRow,ClpSimplex::basic);
      }
    }
  }
  for (int i = 0; i < numberNormal; i++) {
    int iOrig = whichColumns[i];
    ClpSimplex::Status thisStatus = getStatus(i);
    if (thisStatus == ClpSimplex::basic)
      status[iOrig] = 1;
    else if (thisStatus == ClpSimplex::atLowerBound)
      status[iOrig] = 3;
    else if (thisStatus == ClpSimplex::atUpperBound)
      status[iOrig] = 2;
    else if (thisStatus == ClpSimplex::isFixed)
      status[iOrig] = 5;
    else
      abort();
    originalSolution[iOrig] = solution[i];
  }
  for (int i = numberNormal; i < lastOdd; i++) {
    int iOrig = whichColumns[id[i-numberNormal] + numberNormal];
    if (iOrig<numberColumns) {
      ClpSimplex::Status thisStatus = getStatus(i);
      if (thisStatus == ClpSimplex::basic)
	status[iOrig] = 1;
      else if (thisStatus == ClpSimplex::atLowerBound)
	status[iOrig] = 3;
      else if (thisStatus == ClpSimplex::atUpperBound)
	status[iOrig] = 2;
      else if (thisStatus == ClpSimplex::isFixed)
	status[iOrig] = 5;
      else
	abort();
      originalSolution[iOrig] = solution[i];
    } else {
      // slack (basic probably)
      int iSet = iOrig - numberColumns;
      int iRow = whichRows[iSet+numberNonGub];
      ClpSimplex::Status thisStatus = getStatus(i);
      if (thisStatus == ClpSimplex::atLowerBound)
	thisStatus = ClpSimplex::atUpperBound;
      else if (thisStatus == ClpSimplex::atUpperBound)
	thisStatus = ClpSimplex::atLowerBound;
      original.setRowStatus(iRow,thisStatus);
    }
  }
  for (int i = 0; i < numberNonGub; i++) {
    int iOrig = whichRows[i];
    ClpSimplex::Status thisStatus = getRowStatus(i);
    if (thisStatus == ClpSimplex::basic)
      rowStatus[iOrig] = 1;
    else if (thisStatus == ClpSimplex::atLowerBound)
      rowStatus[iOrig] = 3;
    else if (thisStatus == ClpSimplex::atUpperBound)
      rowStatus[iOrig] = 2;
    else if (thisStatus == ClpSimplex::isFixed)
      rowStatus[iOrig] = 5;
    else
      abort();
  }
  int * numberKey = new int [numberRows];
  memset(numberKey,0,numberRows*sizeof(int));
  for (int i=0;i<numberSets;i++) {
    int iRow = whichRows[i+numberNonGub];
    for (int j=startSet[i];j<startSet[i+1];j++) {
      int iOrig = whichColumns[j+numberNormal];
      if (iOrig<numberColumns) {
	if (original.getColumnStatus(iOrig)==ClpSimplex::basic) {
	  numberKey[iRow]++;
	}
      } else {
	// slack
	if (original.getRowStatus(iRow)==ClpSimplex::basic) 
	  numberKey[iRow]++;
      }
    }
  }
  for (int i=0;i<numberSets;i++) {
    int iRow = whichRows[i+numberNonGub];
    if (!numberKey[iRow]) {
      original.setRowStatus(iRow,ClpSimplex::basic);
    }
  }
  delete [] numberKey;
  double objValue = 0.0;
  for (int i = 0; i < numberColumns; i++)
    objValue += cost[i] * originalSolution[i];
  //printf("objective value is %g\n", objValue);
}
/*
  Modifies coefficients etc and if necessary pivots in and out.
  All at same status will be done (basis may go singular).
  User can tell which others have been done (i.e. if status matches).
  If called from outside will change status and return 0 (-1 if not ClpPacked)
  If called from event handler returns non-zero if user has to take action.
  -1 - not ClpPackedMatrix
  0 - no pivots
  1 - pivoted
  2 - pivoted and optimal
  3 - refactorize
  4 - worse than that
  indices>=numberColumns are slacks (obviously no coefficients)
  status array is (char) Status enum
*/
int 
ClpSimplex::modifyCoefficientsAndPivot(int number,
				       const int * which,
				       const CoinBigIndex * start,
				       const int * row,
				       const double * newCoefficient,
				       const unsigned char * newStatus,
				       const double * newLower,
				       const double * newUpper,
				       const double * newObjective)
{
  ClpPackedMatrix* clpMatrix =
    dynamic_cast< ClpPackedMatrix*>(matrix_);
  bool canPivot = lower_!=NULL && factorization_!=NULL;
  int returnCode=0;
  if (!clpMatrix) {
    canPivot=false;
    returnCode=-1;
    // very slow
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	for (CoinBigIndex j=start[i];j<start[i+1];j++) {
 	  matrix_->modifyCoefficient(row[j],iSequence,newCoefficient[j]);
	}
      } else {
	assert (start[i]==start[i+1]);
      }
    }
  } else {
#if 0
    // when in stable
    CoinPackedMatrix * matrix = clpMatrix->getPackedMatrix();
    matrix->modifyCoefficients(number,which,start,
			       row,newCoefficient);
#else
    // Copy and sort which
    int * which2 = new int [2*number+2];
    int * sort = which2+number+1;
    int n=0;
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	which2[n]=iSequence;
	sort[n++]=i;
      } else {
	assert (start[i]==start[i+1]);
      }
    }
    if (n) {
      CoinIndexedVector * rowVector=NULL;
      for (int i=0;i<4;i++) {
	if (rowArray_[i]&&!rowArray_[i]->getNumElements()) {
	  rowVector = rowArray_[i];
	  break;
	}
      }
      bool tempVector=false;
      if (!rowVector) {
	tempVector=true;
	rowVector=new CoinIndexedVector(numberRows_);
      }
      CoinSort_2(which2,which2+n,sort);
      // Stop at end
      which2[n]=numberColumns_;
      sort[n]=n;
      CoinPackedMatrix * matrix = clpMatrix->getPackedMatrix();
      int * rowNow = matrix->getMutableIndices();
      CoinBigIndex * columnStart = matrix->getMutableVectorStarts();
      int * columnLength = matrix->getMutableVectorLengths();
      double * elementByColumn = matrix->getMutableElements();
      double * array = rowVector->denseVector();
      //int * index = rowVector->getIndices();
      int needed=0;
      bool moveUp=false;
      for (int i=0;i<n;i++) {
	int inWhich=sort[i];
	int iSequence=which2[inWhich];
	int nZeroNew=0;
	int nZeroOld=0;
	for (CoinBigIndex j=start[inWhich];j<start[inWhich+1];j++) {
	  int iRow=row[j];
	  double newValue=newCoefficient[j];
	  if (!newValue) {
	    newValue = COIN_INDEXED_REALLY_TINY_ELEMENT;
	    nZeroNew++;
	  }
	  array[iRow]=newValue;
	}
	for (CoinBigIndex j=columnStart[iSequence];
	     j<columnStart[iSequence]+columnLength[iSequence];j++) {
	  int iRow=rowNow[j];
	  double oldValue=elementByColumn[j];
	  if (fabs(oldValue)>COIN_INDEXED_REALLY_TINY_ELEMENT) {
	    double newValue=array[iRow];
	    if (oldValue!=newValue) {
	      if (newValue) {
		array[iRow]=0.0;
		if (newValue==COIN_INDEXED_REALLY_TINY_ELEMENT) {
		  needed--;
		}
	      }
	    }
	  } else {
	    nZeroOld++;
	  }
	}
	assert (!nZeroOld);
	for (CoinBigIndex j=start[inWhich];j<start[inWhich+1];j++) {
	  int iRow=row[j];
	  double newValue=array[iRow];
	  if (newValue) {
	    array[iRow]=0.0;
	    needed++;
	    if (needed>0)
	      moveUp=true;
	  }
	}
      }
      int numberElements = matrix->getNumElements();
      assert (numberElements==columnStart[numberColumns_]);
      if (needed>0) {
	// need more space
	matrix->reserve(numberColumns_, numberElements+needed);
	rowNow = matrix->getMutableIndices();
	elementByColumn = matrix->getMutableElements();
      }
      if (moveUp) {
	// move up from top
	CoinBigIndex top = numberElements+needed;
	for (int iColumn=numberColumns_-1;iColumn>=0;iColumn--) {
	  CoinBigIndex end = columnStart[iColumn+1];
	  columnStart[iColumn+1]=top;
	  CoinBigIndex startThis = columnStart[iColumn];
	  for (CoinBigIndex j=end-1;j>=startThis;j--) {
	    if (elementByColumn[j]) {
	      top--;
	      elementByColumn[top]=elementByColumn[j];
	      rowNow[top]=rowNow[j];
	    }
	  }
	}
	columnStart[0]=top;
      }
      // now move down and insert
      CoinBigIndex put=0;
      int iColumn=0;
      for (int i=0;i<n+1;i++) {
	int inWhich=sort[i];
	int nextMod=which2[inWhich];
	for (;iColumn<nextMod;iColumn++) {
	  CoinBigIndex startThis = columnStart[iColumn];
	  columnStart[iColumn]=put;
	  for (CoinBigIndex j=startThis;
	       j<columnStart[iColumn+1];j++) {
	    int iRow=rowNow[j];
	    double oldValue=elementByColumn[j];
	    if (oldValue) {
	      rowNow[put]=iRow;
	      elementByColumn[put++]=oldValue;
	    }
	  }
	}
	if (i==n) {
	  columnStart[iColumn]=put;
	  break;
	}
	// Now 
	for (CoinBigIndex j=start[inWhich];j<start[inWhich+1];j++) {
	  int iRow=row[j];
	  double newValue=newCoefficient[j];
	  if (!newValue) {
	    newValue = COIN_INDEXED_REALLY_TINY_ELEMENT;
	  }
	  array[iRow]=newValue;
	}
	CoinBigIndex startThis = columnStart[iColumn];
	columnStart[iColumn]=put;
	for (CoinBigIndex j=startThis;
	     j<columnStart[iColumn+1];j++) {
	  int iRow=rowNow[j];
	  double oldValue=elementByColumn[j];
	  if (array[iRow]) {
	    oldValue=array[iRow];
	    if (oldValue==COIN_INDEXED_REALLY_TINY_ELEMENT) 
	      oldValue=0.0;
	    array[iRow]=0.0;
	  }
	  if (fabs(oldValue)>COIN_INDEXED_REALLY_TINY_ELEMENT) {
	    rowNow[put]=iRow;
	    elementByColumn[put++]=oldValue;
	  }
	}
	for (CoinBigIndex j=start[inWhich];j<start[inWhich+1];j++) {
	  int iRow=row[j];
	  double newValue=array[iRow];
	  if (newValue) {
	    array[iRow]=0.0;
	    rowNow[put]=iRow;
	    elementByColumn[put++]=newValue;
	  }
	}
	iColumn++;
      }
      matrix->setNumElements(put);
      if (tempVector)
	delete rowVector;
      for (int i=0;i<numberColumns_;i++) {
	columnLength[i]=columnStart[i+1]-columnStart[i];
      }
    }
#endif
    if (canPivot) {
      // ? faster to modify row copy??
      if (rowCopy_&&start[number]) {
	delete rowCopy_;
	rowCopy_ = clpMatrix->reverseOrderedCopy();
      }
      assert (!newStatus); // do later
      int numberPivots = factorization_->pivots();
      int needed=0;
      for (int i=0;i<number;i++) {
	int iSequence=which[i];
	if (start[i+1]>start[i]&&getStatus(iSequence)==basic)
	  needed++;
      }
      if (needed&&numberPivots+needed<20&&needed<-2) {
	// update factorization
	int saveIn = sequenceIn_;
	int savePivot = pivotRow_;
	int nArray=0;
	CoinIndexedVector * vec[2];
	for (int i=0;i<4;i++) {
	  if (!rowArray_[i]->getNumElements()) {
	    vec[nArray++]=rowArray_[i];
	    if (nArray==2)
	      break;
	  }
	}
	assert (nArray==2); // could use temp array
	for (int i=0;i<number;i++) {
	  int sequenceIn_=which[i];
	  if (start[i+1]>start[i]&&getStatus(sequenceIn_)==basic) {
	    // find pivot row
	    for (pivotRow_=0;pivotRow_<numberRows_;pivotRow_++) {
	      if (pivotVariable_[pivotRow_]==sequenceIn_) 
		break;
	    }
	    assert(pivotRow_<numberRows_);
	    // unpack column
	    assert(!vec[0]->getNumElements());
	    unpackPacked(vec[0]);
	    // update
	    assert(!vec[1]->getNumElements());
	    factorization_->updateColumnFT(vec[1], vec[0]);
	    // Find alpha
	    const int * indices = vec[0]->getIndices();
	    const double * array = vec[0]->denseVector();
	    int n=vec[0]->getNumElements();
	    alpha_=0.0;
	    for (int i=0;i<n;i++) {
	      if (indices[i]==pivotRow_) {
		alpha_ = array[i];
		break;
	      }
	    }
	    int updateStatus=2;
	    if (fabs(alpha_)>1.0e-7) 
	      updateStatus = factorization_->replaceColumn(this,
							   vec[1],
							   vec[0],
							   pivotRow_,
							   alpha_);
	    assert(!vec[1]->getNumElements());
	    vec[0]->clear();
	    if (updateStatus) {
	      returnCode=3;
	      break;
	    }
	  }
	}
	sequenceIn_=saveIn;
	pivotRow_ = savePivot;
	if (!returnCode)
	  returnCode=100; // say can do more
      } else if (needed) {
	returnCode=3; // refactorize
      }
    }
  }
  if (newStatus) {
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      status_[iSequence]=newStatus[i];
    }
  }
  if (newLower) {
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	if (columnLower_[iSequence]!=newLower[i]) {
	  columnLower_[iSequence]=newLower[i];
	}
      } else {
	iSequence -= numberColumns_;
	if (rowLower_[iSequence]!=newLower[i]) {
	  rowLower_[iSequence]=newLower[i];
	}
      }
    }
  }
  if (newLower) {
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	if (columnLower_[iSequence]!=newLower[i]) {
	  columnLower_[iSequence]=newLower[i];
	}
      } else {
	iSequence -= numberColumns_;
	if (rowLower_[iSequence]!=newLower[i]) {
	  rowLower_[iSequence]=newLower[i];
	}
      }
    }
  }
  if (newUpper) {
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	if (columnUpper_[iSequence]!=newUpper[i]) {
	  columnUpper_[iSequence]=newUpper[i];
	}
      } else {
	iSequence -= numberColumns_;
	if (rowUpper_[iSequence]!=newUpper[i]) {
	  rowUpper_[iSequence]=newUpper[i];
	}
      }
    }
  }
  if (newObjective) {
    double * obj = objective();
    for (int i=0;i<number;i++) {
      int iSequence=which[i];
      if (iSequence<numberColumns_) {
	if (obj[iSequence]!=newObjective[i]) {
	  obj[iSequence]=newObjective[i];
	}
      } else {
	  assert (!newObjective[i]);
      }
    }
  }
  if (canPivot) {
    // update lower, upper, objective and nonLinearCost
    assert (!rowScale_); // for now
    memcpy(lower_,columnLower_,numberColumns_*sizeof(double));
    memcpy(lower_+numberColumns_,rowLower_,numberRows_*sizeof(double));
    memcpy(upper_,columnUpper_,numberColumns_*sizeof(double));
    memcpy(upper_+numberColumns_,rowUpper_,numberRows_*sizeof(double));
    memcpy(cost_,objective(),numberColumns_*sizeof(double));
    memset(cost_+numberColumns_,0,numberRows_*sizeof(double));
    // ? parameter to say no gutsOfSolution needed
    // make sure slacks have correct value
    // see if still optimal
    if (returnCode==100) {
      // is this needed
      if (nonLinearCost_) {
	// speed up later
	//nonLinearCost_->checkInfeasibilities(oldTolerance);
	delete nonLinearCost_;
	nonLinearCost_ = new ClpNonLinearCost(this);
      }
      gutsOfSolution(NULL,NULL,false);
      assert (!newStatus);
      printf("%d primal %d dual\n",numberPrimalInfeasibilities_,
	     numberDualInfeasibilities_);
      returnCode=3;
    } else {
      // is this needed
      if (nonLinearCost_) {
	// speed up later
#if 1
	for (int i=0;i<number;i++) {
	  int iSequence=which[i];
	  nonLinearCost_->setOne(iSequence,solution_[iSequence],
				 lower_[iSequence],upper_[iSequence],
				 cost_[iSequence]);
	}
#else	
	//nonLinearCost_->checkInfeasibilities(oldTolerance);
	delete nonLinearCost_;
	nonLinearCost_ = new ClpNonLinearCost(this);
	//nonLinearCost_->checkInfeasibilities(0.0);
#endif
	//gutsOfSolution(NULL,NULL,false);
	assert (!newStatus);
      }
    }
  }
  return returnCode;
}

/* Pivot out a variable and choose an incoing one.  Assumes dual
   feasible - will not go through a reduced cost.
   Returns step length in theta
   Return codes as before but -1 means no acceptable pivot
*/
int
ClpSimplex::dualPivotResultPart1()
{
     return static_cast<ClpSimplexDual *> (this)->pivotResultPart1();
}
/* Do actual pivot
   state is 1,3 if got tableau column in rowArray_[1]
   2,3 if got tableau row in rowArray_[0] and columnArray_[0]
*/
int 
ClpSimplex::pivotResultPart2(int algorithm,int state)
{
  if (!(state&1)) {
    // update the incoming column
    unpackPacked(rowArray_[1]);
    factorization_->updateColumnFT(rowArray_[2], rowArray_[1]);
  }
#define CHECK_TABLEAU 0
  if (!(state&2)||CHECK_TABLEAU) {
    // get tableau row
    // create as packed
    double direction = directionOut_;
    assert (!rowArray_[2]->getNumElements());
    assert (!columnArray_[1]->getNumElements());
#if CHECK_TABLEAU
    printf("rowArray0 old\n");
    rowArray_[0]->print();
    rowArray_[0]->clear();
    printf("columnArray0 old\n");
    columnArray_[0]->print();
    columnArray_[0]->clear();
#else
    assert (!columnArray_[0]->getNumElements());
    assert (!rowArray_[0]->getNumElements());
#endif
    rowArray_[0]->createPacked(1, &pivotRow_, &direction);
    factorization_->updateColumnTranspose(rowArray_[2], rowArray_[0]);
    rowArray_[3]->clear();
    // put row of tableau in rowArray[0] and columnArray[0]
    assert (!rowArray_[2]->getNumElements());
    matrix_->transposeTimes(this, -1.0,
			    rowArray_[0], rowArray_[2], columnArray_[0]);
#if CHECK_TABLEAU
    printf("rowArray0 new\n");
    rowArray_[0]->print();
    printf("columnArray0 new\n");
    columnArray_[0]->print();
#endif
  }
  assert (pivotRow_>=0);
  assert (sequenceIn_>=0);
  assert (sequenceOut_>=0);
  int returnCode=-1;
  if (algorithm>0) {
    // replace in basis
    int updateStatus = factorization_->replaceColumn(this,
						     rowArray_[2],
						     rowArray_[1],
						     pivotRow_,
						     alpha_);
    if (!updateStatus) {
      dualIn_ = cost_[sequenceIn_];
      double * work = rowArray_[1]->denseVector();
      int number = rowArray_[1]->getNumElements();
      int * which = rowArray_[1]->getIndices();
      for (int i = 0; i < number; i++) {
	int iRow = which[i];
	double alpha = work[i];
	int iPivot = pivotVariable_[iRow];
	dualIn_ -= alpha * cost_[iPivot];
      }
      returnCode=0;
      double multiplier = dualIn_ / alpha_;
      // update column djs
      int i;
      int * index = columnArray_[0]->getIndices();
      number = columnArray_[0]->getNumElements();
      double * element = columnArray_[0]->denseVector();
      assert (columnArray_[0]->packedMode());
      for (i = 0; i < number; i++) {
	int iSequence = index[i];
	dj_[iSequence] += multiplier*element[i];
	reducedCost_[iSequence] = dj_[iSequence];
	element[i] = 0.0;
      }
      columnArray_[0]->setNumElements(0);
      // and row djs
      index = rowArray_[0]->getIndices();
      number = rowArray_[0]->getNumElements();
      element = rowArray_[0]->denseVector();
      assert (rowArray_[0]->packedMode());
      for (i = 0; i < number; i++) {
	int iSequence = index[i];
	dj_[iSequence+numberColumns_] += multiplier*element[i];
	dual_[iSequence] = dj_[iSequence+numberColumns_];
	element[i] = 0.0;
      }
      rowArray_[0]->setNumElements(0);
      double oldCost = cost_[sequenceOut_];
      // update primal solution
      
      double objectiveChange = 0.0;
      // after this rowArray_[1] is not empty - used to update djs
      static_cast<ClpSimplexPrimal *>(this)->updatePrimalsInPrimal(rowArray_[1], theta_, objectiveChange, 0);
      double oldValue = valueIn_;
      if (directionIn_ == -1) {
	// as if from upper bound
	if (sequenceIn_ != sequenceOut_) {
	  // variable becoming basic
	  valueIn_ -= fabs(theta_);
	} else {
	  valueIn_ = lowerIn_;
	}
      } else {
	// as if from lower bound
	if (sequenceIn_ != sequenceOut_) {
	  // variable becoming basic
	  valueIn_ += fabs(theta_);
	} else {
	  valueIn_ = upperIn_;
	}
      }
      objectiveChange += dualIn_ * (valueIn_ - oldValue);
      // outgoing
      if (sequenceIn_ != sequenceOut_) {
	if (directionOut_ > 0) {
	  valueOut_ = lowerOut_;
	} else {
	  valueOut_ = upperOut_;
	}
	if(valueOut_ < lower_[sequenceOut_] - primalTolerance_)
	  valueOut_ = lower_[sequenceOut_] - 0.9 * primalTolerance_;
	else if (valueOut_ > upper_[sequenceOut_] + primalTolerance_)
	  valueOut_ = upper_[sequenceOut_] + 0.9 * primalTolerance_;
	// may not be exactly at bound and bounds may have changed
	// Make sure outgoing looks feasible
	directionOut_ = nonLinearCost_->setOneOutgoing(sequenceOut_, valueOut_);
	// May have got inaccurate
	//if (oldCost!=cost_[sequenceOut_])
	//printf("costchange on %d from %g to %g\n",sequenceOut_,
	//       oldCost,cost_[sequenceOut_]);
	dj_[sequenceOut_] = cost_[sequenceOut_] - oldCost; // normally updated next iteration
	solution_[sequenceOut_] = valueOut_;
      }
      // change cost and bounds on incoming if primal
      nonLinearCost_->setOne(sequenceIn_, valueIn_);
      progress_.startCheck(); // make sure won't worry about cycling
      int whatNext = housekeeping(objectiveChange);
      if (whatNext == 1) {
	returnCode = -2; // refactorize
      } else if (whatNext == 2) {
	// maximum iterations or equivalent
	returnCode = 3;
      } else if(numberIterations_ == lastGoodIteration_
		+ 2 * factorization_->maximumPivots()) {
	// done a lot of flips - be safe
	returnCode = -2; // refactorize
      }
    } else {
      // ?
      abort();
    }
  } else {
    // dual
    // recompute dualOut_
    if (directionOut_ < 0) {
      dualOut_ = valueOut_ - upperOut_;
    } else {
      dualOut_ = lowerOut_ - valueOut_;
    }
    // update the incoming column
    double btranAlpha = -alpha_ * directionOut_; // for check
    rowArray_[1]->clear();
    unpackPacked(rowArray_[1]);
    // moved into updateWeights - factorization_->updateColumnFT(rowArray_[2],rowArray_[1]);
    // and update dual weights (can do in parallel - with extra array)
    alpha_ = dualRowPivot_->updateWeights(rowArray_[0],
					  rowArray_[2],
					  rowArray_[3],
					  rowArray_[1]);
    // see if update stable
#ifdef CLP_DEBUG
    if ((handler_->logLevel() & 32))
      printf("btran alpha %g, ftran alpha %g\n", btranAlpha, alpha_);
#endif
    double checkValue = 1.0e-7;
    // if can't trust much and long way from optimal then relax
    if (largestPrimalError_ > 10.0)
      checkValue = CoinMin(1.0e-4, 1.0e-8 * largestPrimalError_);
    if (fabs(btranAlpha) < 1.0e-12 || fabs(alpha_) < 1.0e-12 ||
	fabs(btranAlpha - alpha_) > checkValue*(1.0 + fabs(alpha_))) {
      handler_->message(CLP_DUAL_CHECK, messages_)
	<< btranAlpha
	<< alpha_
	<< CoinMessageEol;
      if (factorization_->pivots()) {
	dualRowPivot_->unrollWeights();
	problemStatus_ = -2; // factorize now
	rowArray_[0]->clear();
	rowArray_[1]->clear();
	columnArray_[0]->clear();
	returnCode = -2;
	abort();
	return returnCode;
      } else {
	// take on more relaxed criterion
	double test;
	if (fabs(btranAlpha) < 1.0e-8 || fabs(alpha_) < 1.0e-8)
	  test = 1.0e-1 * fabs(alpha_);
	else
	  test = 1.0e-4 * (1.0 + fabs(alpha_));
	if (fabs(btranAlpha) < 1.0e-12 || fabs(alpha_) < 1.0e-12 ||
	    fabs(btranAlpha - alpha_) > test) {
	  abort();
	}
      }
    }
    // update duals BEFORE replaceColumn so can do updateColumn
    double objectiveChange = 0.0;
    // do duals first as variables may flip bounds
    // rowArray_[0] and columnArray_[0] may have flips
    // so use rowArray_[3] for work array from here on
    int nswapped = 0;
    //rowArray_[0]->cleanAndPackSafe(1.0e-60);
    //columnArray_[0]->cleanAndPackSafe(1.0e-60);
    // make sure incoming doesn't count
    Status saveStatus = getStatus(sequenceIn_);
    setStatus(sequenceIn_, basic);
    nswapped = 
      static_cast<ClpSimplexDual *>(this)->updateDualsInDual(rowArray_[0], columnArray_[0],
				 rowArray_[2], theta_,
				 objectiveChange, false);
    assert (!nswapped);
    setStatus(sequenceIn_, saveStatus);
    double oldDualOut = dualOut_;
    // which will change basic solution
    if (nswapped) {
      if (rowArray_[2]->getNumElements()) {
	factorization_->updateColumn(rowArray_[3], rowArray_[2]);
	dualRowPivot_->updatePrimalSolution(rowArray_[2],
					    1.0, objectiveChange);
      }
      // recompute dualOut_
      valueOut_ = solution_[sequenceOut_];
      if (directionOut_ < 0) {
	dualOut_ = valueOut_ - upperOut_;
      } else {
	dualOut_ = lowerOut_ - valueOut_;
      }
    }
    // amount primal will move
    double movement = -dualOut_ * directionOut_ / alpha_;
    double movementOld = oldDualOut * directionOut_ / alpha_;
    // so objective should increase by fabs(dj)*movement
    // but we already have objective change - so check will be good
    if (objectiveChange + fabs(movementOld * dualIn_) < -CoinMax(1.0e-5, 1.0e-12 * fabs(objectiveValue_))) {
      if (handler_->logLevel() & 32)
	printf("movement %g, swap change %g, rest %g  * %g\n",
	       objectiveChange + fabs(movement * dualIn_),
	       objectiveChange, movement, dualIn_);
    }
    // if stable replace in basis
    int updateStatus = factorization_->replaceColumn(this,
						     rowArray_[2],
						     rowArray_[1],
						     pivotRow_,
						     alpha_);
    // If looks like bad pivot - refactorize
    if (fabs(dualOut_) > 1.0e50)
      updateStatus = 2;
    // if no pivots, bad update but reasonable alpha - take and invert
    if (updateStatus == 2 &&
	!factorization_->pivots() && fabs(alpha_) > 1.0e-5)
      updateStatus = 4;
    if (updateStatus == 1 || updateStatus == 4) {
      // slight error
      if (factorization_->pivots() > 5 || updateStatus == 4) {
	problemStatus_ = -2; // factorize now
	returnCode = -3;
      }
    } else if (updateStatus == 2) {
      // major error
      dualRowPivot_->unrollWeights();
      // later we may need to unwind more e.g. fake bounds
      if (factorization_->pivots() &&
	  ((moreSpecialOptions_ & 16) == 0 || factorization_->pivots() > 4)) {
	problemStatus_ = -2; // factorize now
	returnCode = -2;
	moreSpecialOptions_ |= 16;
	return returnCode;
      } else {
	// need to reject something
	abort();
      }
    } else if (updateStatus == 3) {
      // out of memory
      // increase space if not many iterations
      if (factorization_->pivots() <
	  0.5 * factorization_->maximumPivots() &&
	  factorization_->pivots() < 200)
	factorization_->areaFactor(
				   factorization_->areaFactor() * 1.1);
      problemStatus_ = -2; // factorize now
    } else if (updateStatus == 5) {
      problemStatus_ = -2; // factorize now
    }
    // update primal solution
    if (theta_ < 0.0) {
      if (handler_->logLevel() & 32)
	printf("negative theta %g\n", theta_);
      theta_ = 0.0;
    }
    // do actual flips (should not be any?)
    static_cast<ClpSimplexDual *>(this)->flipBounds(rowArray_[0], columnArray_[0]);
      //rowArray_[1]->expand();
    dualRowPivot_->updatePrimalSolution(rowArray_[1],
					movement,
					objectiveChange);
    // modify dualout
    dualOut_ /= alpha_;
    dualOut_ *= -directionOut_;
    //setStatus(sequenceIn_,basic);
    dj_[sequenceIn_] = 0.0;
    double oldValue = valueIn_;
    if (directionIn_ == -1) {
      // as if from upper bound
      valueIn_ = upperIn_ + dualOut_;
    } else {
      // as if from lower bound
      valueIn_ = lowerIn_ + dualOut_;
    }
    objectiveChange += cost_[sequenceIn_] * (valueIn_ - oldValue);
    // outgoing
    // set dj to zero unless values pass
    if (directionOut_ > 0) {
      valueOut_ = lowerOut_;
      dj_[sequenceOut_] = theta_;
    } else {
      valueOut_ = upperOut_;
      dj_[sequenceOut_] = -theta_;
    }
    solution_[sequenceOut_] = valueOut_;
    int whatNext = housekeeping(objectiveChange);
    // and set bounds correctly
    static_cast<ClpSimplexDual *>(this)->originalBound(sequenceIn_);
    static_cast<ClpSimplexDual *>(this)->changeBound(sequenceOut_);
    if (whatNext == 1) {
      problemStatus_ = -2; // refactorize
    } else if (whatNext == 2) {
      // maximum iterations or equivalent
      problemStatus_ = 3;
      returnCode = 3;
      abort();
    }
  }
  // Check event
  {
    int status = eventHandler_->event(ClpEventHandler::endOfIteration);
    if (status >= 0) {
      problemStatus_ = 5;
      secondaryStatus_ = ClpEventHandler::endOfIteration;
      returnCode = 3;
    }
  }
  // need to be able to refactoriza
  //printf("return code %d problem status %d\n",
  //	 returnCode,problemStatus_);
  return returnCode;
}
#include "CoinPresolveMatrix.hpp"
/* Take out duplicate rows (includes scaled rows and intersections).
   On exit whichRows has rows to delete - return code is number can be deleted 
   or -1 if would be infeasible.
   If tolerance is -1.0 use primalTolerance for equality rows and infeasibility
   If cleanUp not zero then spend more time trying to leave more stable row
   and make row bounds exact multiple of cleanUp if close enough
   moves status to try and delete a basic slack
*/
int 
ClpSimplex::outDuplicateRows(int numberLook,int * whichRows, double tolerance,
			  double cleanUp) 
{
  double * weights = new double [numberLook+numberColumns_];
  double * columnWeights = weights+numberLook;
  coin_init_random_vec(columnWeights,numberColumns_);
  // get row copy
  CoinPackedMatrix rowCopy = *matrix();
  rowCopy.reverseOrdering();
  int * column = rowCopy.getMutableIndices();
  CoinBigIndex * rowStart = rowCopy.getMutableVectorStarts();
  int * rowLength = rowCopy.getMutableVectorLengths();
  double * element = rowCopy.getMutableElements();
  for (int i=0;i<numberLook;i++) {
    int iRow=whichRows[i];
    double value = 0.0;
    CoinBigIndex start=rowStart[iRow];
    CoinBigIndex end = start+rowLength[iRow];
    // sort (probably in order anyway)
    CoinSort_2(column+start,column+end,element+start);
    for (CoinBigIndex j=start;j<end;j++) {
      int iColumn = column[j];
      value += columnWeights[iColumn];
    }
    weights[iRow]=value;
    whichRows[iRow]=iRow;
  }
  CoinSort_2(weights,weights+numberLook,whichRows);
  if (tolerance<0.0)
    tolerance = primalTolerance_;
  int nPossible=0;
  int nDelete=0;
  double value = weights[0];
  double inverseCleanup = (cleanUp>0.0) ? 1.0/cleanUp : 0.0;
  int firstSame=-1;
  int lastSame=-1;
  for (int iLook = 1; iLook < numberLook; iLook++) {
    if (weights[iLook]==value) {
      int iLast=whichRows[iLook-1];
      if (firstSame<0) {
	/* see how many same - if >2 but < ? may be 
	   worth looking at all combinations
	*/
	firstSame=iLook-1;
	for (lastSame=iLook+1;lastSame<numberLook;lastSame++) {
	  if (weights[lastSame]!=value)
	    break;
	}
	//#define PRINT_DUP
#ifdef PRINT_DUP
	if (lastSame-firstSame>2)
	  printf("dupsame %d rows have same weight\n",lastSame-firstSame);
#endif
      }
      int iThis=whichRows[iLook];
      CoinBigIndex start = rowStart[iThis];
      CoinBigIndex end = start + rowLength[iThis];
      if (rowLength[iThis] == rowLength[iLast]) {
	nPossible++;
#ifdef PRINT_DUP
	char line[520],temp[50];
#endif
	int ishift = rowStart[iLast] - start;
	CoinBigIndex k;
#ifdef PRINT_DUP
	sprintf(line,"dupj %d,%d %d els ",
		iThis,iLast,rowLength[iLast]);
	int n=strlen(line);
#endif
	bool bad=false;
	double multiplier=0.0;
	for (k=start;k<end;k++) {
	  if (column[k] != column[k+ishift]) {
	    bad=true;
	    break;
	  } else {
#ifdef PRINT_DUP
	    sprintf(temp,"(%g,%g) ",element[k],element[k+ishift]);
#endif
	    if (!multiplier) {
	      multiplier = element[k]/element[k+ishift];
	    } else if(fabs(element[k+ishift]*multiplier-element[k])>1.0e-8) {
	      bad=true;
	    }
#ifdef PRINT_DUP
	    int n2=strlen(temp);
	    if (n+n2<500) {
	      strcat(line,temp);
	      n += n2;
	    } else {
	      strcat(line,"...");
	      break;
	    }
#endif
	  }
	}
	if (!bad) {
#ifdef PRINT_DUP
	  printf("%s lo (%g,%g) up (%g,%g) - multiplier %g\n",line,rowLower_[iThis],rowUpper_[iThis],
		 rowLower_[iLast],rowUpper_[iLast],multiplier);
#endif
	  double rlo1=rowLower_[iLast];
	  double rup1=rowUpper_[iLast];
	  double rlo2=rowLower_[iThis];
	  double rup2=rowUpper_[iThis];
	  // scale
	  rlo1 *= multiplier;
	  rup1 *= multiplier;
	  //swap bounds if neg
	  if (multiplier<0.0) {
	    double temp = rup1;
	    rup1=rlo1;
	    rlo1=temp;
	  }
	  /* now check rhs to see what is what */
#ifdef PRINT_DUP
	  printf("duplicate row %g %g, %g %g\n",
		 rlo1,rup1,rlo2,rup2);
#endif
	  /* we always keep this and delete last 
	     later maybe keep better formed one */
	  rlo2 = CoinMax(rlo1,rlo2);
	  if (rlo2<-1.0e30)
	    rlo2=-COIN_DBL_MAX;
	  rup2 = CoinMin(rup1,rup2);
	  if (rup2>1.0e30)
	    rup2=COIN_DBL_MAX;
#ifdef PRINT_DUP
	  printf("pre_duprow %dR %dR keep this\n",iLast,iThis);
#endif
	  if (rup2<rlo2-tolerance) {
	    // infeasible
	    nDelete=-1;
	    break;
	  } else if (fabs(rup2-rlo2)<=tolerance) {
	    // equal - choose closer to zero
	    if (fabs(rup2)<fabs(rlo2))
	      rlo2=rup2;
	    else
	      rup2=rlo2;
#ifdef PRINT_DUP
	    printf("Row %d has %d elements == %g\n",
		   iThis,end-start,rlo2);
#endif
	  }
	  if (cleanUp>0.0) {
	    /* see if close to multiple
	       always allow integer values */
	    if (rlo2>-1.0e30) {
	      double value = rlo2;
	      double value2 = floor(value+0.5);
	      if (fabs(value-value2)<1.0e-9) {
		rlo2=value2;
	      } else {
		value = rlo2*inverseCleanup;
		value2 = floor(value+0.5);
		if (fabs(value-value2)<1.0e-9)
		  rlo2=value2*cleanUp;
	      }
	    }
	    if (rup2<1.0e30) {
	      double value = rup2;
	      double value2 = floor(value+0.5);
	      if (fabs(value-value2)<1.0e-9) {
		rup2=value2;
	      } else {
		value = rup2*inverseCleanup;
		value2 = floor(value+0.5);
		if (fabs(value-value2)<1.0e-9)
		  rup2=value2*cleanUp;
	      }
	    }
	  }
	  rowLower_[iThis]=rlo2;
	  rowUpper_[iThis]=rup2;
	  whichRows[nDelete++]=iLast;
	  if (getRowStatus(iLast)!=basic) {
	    if (getRowStatus(iThis)==basic) {
	      setRowStatus(iThis,superBasic);
	      setRowStatus(iLast,basic);
	    }
	  }
	} else {
#ifdef PRINT_DUP
	  printf("%s lo (%g,%g) up (%g,%g) - ODD\n",line,rowLower_[iThis],rowUpper_[iThis],
		 rowLower_[iLast],rowUpper_[iLast]);
#endif
	}
      }
    } else {
      // say no match
      firstSame=-1;
      value=weights[iLook];
    }
  }
#ifdef PRINT_DUP
  printf("%d possible duplicate rows - deleting %d\n",
	 nPossible,nDelete);
#endif
  delete [] weights;
  return nDelete;
}
/* Try simple crash like techniques to get closer to primal feasibility
   returns final sum of infeasibilities */
double 
ClpSimplex::moveTowardsPrimalFeasible()
{
  memset (rowActivity_,0,numberRows_*sizeof(double));
  matrix()->times(columnActivity_,rowActivity_);
  double sum=0.0;
  int * which = new int[numberRows_];
  int numberLook=0;
  for (int iRow=0;iRow<numberRows_;iRow++) {
    double value = rowActivity_[iRow];
    double infeasibility = 0.0;
    if (value<rowLower_[iRow]-primalTolerance_) 
      infeasibility = rowLower_[iRow]-value;
    else if (value>rowUpper_[iRow]+primalTolerance_) 
      infeasibility = value-rowUpper_[iRow];
    if (infeasibility) {
      sum += infeasibility;
      which[numberLook++]=iRow;
    }
  }
  if (numberLook) {
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths();
    const double * element = matrix_->getElements();
    // get row copy
    CoinPackedMatrix rowCopy = *matrix();
    rowCopy.reverseOrdering();
    const int * column = rowCopy.getIndices();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    const int * rowLength = rowCopy.getVectorLengths();
    const double * elementByRow = rowCopy.getElements();
    double lastSum=COIN_DBL_MAX;
    while (sum>primalTolerance_&&numberLook) {
      sum =0.0;
      double worst=primalTolerance_;
      int iWorst=-1;
      int n=numberLook;
      numberLook=0;
      for (int iLook=0;iLook<n;iLook++) {
	int iRow=which[iLook];
	double value = rowActivity_[iRow];
	double infeasibility = 0.0;
	if (value<rowLower_[iRow]-primalTolerance_) 
	  infeasibility = rowLower_[iRow]-value;
	else if (value>rowUpper_[iRow]+primalTolerance_) 
	  infeasibility = value-rowUpper_[iRow];
	if (infeasibility) {
	  sum += infeasibility;
	  which[numberLook++]=iRow;
	  if (infeasibility>worst) {
	    worst = infeasibility;
	    iWorst=iRow;
	  }
	}
      }
      if (sum==0.0||sum>=lastSum)
	break;
      lastSum=sum;
      double direction;
      if (rowActivity_[iWorst]<rowLower_[iWorst])
	direction=1.0; // increase
      else
	direction=-1.0;
      for (CoinBigIndex k=rowStart[iWorst];
	   k<rowStart[iWorst]+rowLength[iWorst];k++) {
	if (worst<primalTolerance_) 
	  break;
	int iColumn = column[k];
	double value=elementByRow[k]*direction;
	double distance=worst;
	double multiplier = (value>0.0) ? 1.0 : -1.0;
	// but allow for column bounds
	double currentValue = columnActivity_[iColumn];
	if (multiplier>0.0)
	  distance = CoinMin(worst,columnUpper_[iColumn]-currentValue);
	else
	  distance = CoinMin(worst,currentValue-columnLower_[iColumn]);
	distance /= fabs(value);
	for (CoinBigIndex i=columnStart[iColumn];
	     i<columnStart[iColumn]+columnLength[iColumn];i++) {
	  int iRow=row[i];
	  if (iRow!=iWorst) {
	    double value2=element[i]*multiplier;
	    if (value2>0.0) {
	      double distance2 = rowUpper_[iRow]-rowActivity_[iRow]; 
	      if (value2*distance>distance2)
		distance = distance2/value2;
	    } else {
	      double distance2 = rowLower_[iRow]-rowActivity_[iRow]; 
	      if (value2*distance<distance2)
		distance = distance2/value2;
	    }
	  }
	}
	if (distance>1.0e-12) {
	  worst-=distance*fabs(value);
	  distance *= multiplier;
	  columnActivity_[iColumn] = currentValue+distance;
	  for (CoinBigIndex i=columnStart[iColumn];
	       i<columnStart[iColumn]+columnLength[iColumn];i++) {
	    int iRow=row[i];
	    rowActivity_[iRow] += distance*element[i];
	  }
	}
      }
    }
  }
  delete [] which;
  return sum;
}
/* Try simple crash like techniques to remove super basic slacks
   but only if > threshold */
void 
ClpSimplex::removeSuperBasicSlacks(int threshold)
{
  // could try going both ways - for first attempt to nearer bound
  memset (rowActivity_,0,numberRows_*sizeof(double));
  matrix()->times(columnActivity_,rowActivity_);
  double * distance = new double [numberRows_];
  int * whichRows = new int [numberRows_];
  int numberLook=0;
  for (int iRow=0;iRow<numberRows_;iRow++) {
    if (getRowStatus(iRow)!=basic) {
      double value = rowActivity_[iRow];
      if (value>rowLower_[iRow]+primalTolerance_&&
	  value<rowUpper_[iRow]-primalTolerance_) {
	setRowStatus(iRow,superBasic);
	distance[numberLook]=CoinMin(value-rowLower_[iRow],
				      rowUpper_[iRow]-value);
	whichRows[numberLook++]=iRow;
      }
    }
  }
  if (numberLook>threshold) {
    CoinSort_2(distance,distance+numberLook,whichRows);
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths();
    const double * element = matrix_->getElements();
    // get row copy
    CoinPackedMatrix rowCopy = *matrix();
    rowCopy.reverseOrdering();
    const int * column = rowCopy.getIndices();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    const int * rowLength = rowCopy.getVectorLengths();
    const double * elementByRow = rowCopy.getElements();
    int nMoved=0;
    for (int iLook=0;iLook<numberLook;iLook++) {
      int kRow = whichRows[iLook];
      double direction;
      double needed;
      if (rowUpper_[kRow]-rowActivity_[kRow]<rowActivity_[kRow]-rowLower_[kRow]) {
	direction=1.0; // increase
	needed = rowUpper_[kRow]-rowActivity_[kRow];
      } else {
	direction=-1.0;
	needed = rowActivity_[kRow]-rowLower_[kRow];
      }
      for (CoinBigIndex k=rowStart[kRow];
	 k<rowStart[kRow]+rowLength[kRow];k++) {
	if (needed<primalTolerance_) 
	  break;
	int iColumn = column[k];
	if (getColumnStatus(iColumn)!=basic)
	  continue;
	double value=elementByRow[k]*direction;
	double distance;
	double multiplier = (value>0.0) ? 1.0 : -1.0;
	// but allow for column bounds
	double currentValue = columnActivity_[iColumn];
	if (multiplier>0.0)
	  distance = columnUpper_[iColumn]-currentValue;
	else
	  distance = currentValue-columnLower_[iColumn];
	for (CoinBigIndex i=columnStart[iColumn];
	     i<columnStart[iColumn]+columnLength[iColumn];i++) {
	  int iRow=row[i];
	  double value2=element[i]*multiplier;
	  if (value2>0.0) {
	    double distance2 = rowUpper_[iRow]-rowActivity_[iRow]; 
	    if (value2*distance>distance2)
	      distance = distance2/value2;
	  } else {
	    double distance2 = rowLower_[iRow]-rowActivity_[iRow]; 
	    if (value2*distance<distance2)
	      distance = distance2/value2;
	  }
	}
	if (distance>1.0e-12) {
	  distance *= multiplier;
	  columnActivity_[iColumn] = currentValue+distance;
	  for (CoinBigIndex i=columnStart[iColumn];
	       i<columnStart[iColumn]+columnLength[iColumn];i++) {
	    int iRow=row[i];
	    rowActivity_[iRow] += distance*element[i];
	  }
	  if (direction>0.0) {
	    needed = rowUpper_[kRow]-rowActivity_[kRow];
	  } else {
	    needed = rowActivity_[kRow]-rowLower_[kRow];
	  }
	}
      }
      if (needed<primalTolerance_) {
	nMoved++;
	if (rowUpper_[kRow]-rowActivity_[kRow]<primalTolerance_) 
	  setRowStatus(kRow,atUpperBound);
	else if (rowActivity_[kRow]-rowLower_[kRow]<primalTolerance_)
	  setRowStatus(kRow,atLowerBound);
	else
	  assert (rowUpper_[kRow]-rowActivity_[kRow]<primalTolerance_||
		  rowActivity_[kRow]-rowLower_[kRow]<primalTolerance_);
      }
    }
    char line[100];
    sprintf(line,"Threshold %d found %d fixed %d",threshold,numberLook,nMoved);
    handler_->message(CLP_GENERAL,messages_)
      << line << CoinMessageEol;
  }
  delete [] distance;
  delete [] whichRows;
}
