// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpFactorization.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinBuild.hpp"
#include "CoinMpsIO.hpp"
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
void ClpSimplexOther::dualRanging(int numberCheck,const int * which,
			    double * costIncreased, int * sequenceIncreased,
			    double * costDecreased, int * sequenceDecreased)
{
  rowArray_[0]->clear();
  rowArray_[1]->clear();
  rowArray_[3]->clear();
  columnArray_[0]->clear();
  for ( int i=0;i<numberCheck;i++) {
    int iSequence = which[i];
    double costIncrease=COIN_DBL_MAX;
    double costDecrease=COIN_DBL_MAX;
    int sequenceIncrease=-1;
    int sequenceDecrease=-1;
    
    switch(getStatus(iSequence)) {
      
    case basic:
      {
	// non-trvial
	// Find pivot row (could be faster)
	int iRow=-1;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  if (iSequence == pivotVariable_[iRow]) 
	    break;
	}
	assert (iRow>=0);
	double plusOne=1.0;
        rowArray_[0]->createPacked(1,&iRow,&plusOne);
	factorization_->updateColumnTranspose(rowArray_[1],rowArray_[0]);
	// put row of tableau in rowArray[0] and columnArray[0]
	matrix_->transposeTimes(this,-1.0,
				rowArray_[0],rowArray_[3],columnArray_[0]);
	// do ratio test up and down
	checkDualRatios(rowArray_[0],columnArray_[0],costIncrease,sequenceIncrease,
		    costDecrease,sequenceDecrease);
      }
      break;
    case isFixed:
      break;
    case isFree:
    case superBasic:
      costIncrease=0.0;
      costDecrease=0.0;
      sequenceIncrease=iSequence;
      sequenceDecrease=iSequence;
      break;
    case atUpperBound:
      costIncrease = CoinMax(0.0,-dj_[iSequence]);
      sequenceIncrease = iSequence;
      break;
    case atLowerBound:
      costDecrease = CoinMax(0.0,dj_[iSequence]);
      sequenceDecrease = iSequence;
      break;
    }
    double scaleFactor;
    if (rowScale_) {
      if (iSequence<numberColumns_) 
	scaleFactor = 1.0/(objectiveScale_*columnScale_[iSequence]);
      else
	scaleFactor = rowScale_[iSequence-numberColumns_]/objectiveScale_;
    } else {
      scaleFactor = 1.0/objectiveScale_;
    }
    if (costIncrease<1.0e30)
      costIncrease *= scaleFactor;
    if (costDecrease<1.0e30)
      costDecrease *= scaleFactor;
    if (optimizationDirection_==1.0) {
      costIncreased[i] = costIncrease;
      sequenceIncreased[i] = sequenceIncrease;
      costDecreased[i] = costDecrease;
      sequenceDecreased[i] = sequenceDecrease;
    } else if (optimizationDirection_==-1.0) {
      costIncreased[i] = costDecrease;
      sequenceIncreased[i] = sequenceDecrease;
      costDecreased[i] = costIncrease;
      sequenceDecreased[i] = sequenceIncrease;
    } else if (optimizationDirection_==0.0) {
      // !!!!!! ???
      costIncreased[i] = COIN_DBL_MAX;
      sequenceIncreased[i] = -1;
      costDecreased[i] = COIN_DBL_MAX;
      sequenceDecreased[i] = -1;
    } else {
      abort();
    }
  }
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
			     double & costIncrease, int & sequenceIncrease,
			     double & costDecrease, int & sequenceDecrease)
{
  double acceptablePivot = 1.0e-7;
  double * work;
  int number;
  int * which;
  int iSection;

  double thetaDown = 1.0e31;
  double thetaUp = 1.0e31;
  int sequenceDown =-1;
  int sequenceUp = -1;

  int addSequence;

  for (iSection=0;iSection<2;iSection++) {

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
    
    for (i=0;i<number;i++) {
      int iSequence = which[i];
      int iSequence2 = iSequence + addSequence;
      double alpha=work[i];
      if (fabs(alpha)<acceptablePivot)
	continue;
      double oldValue=dj_[iSequence2];

      switch(getStatus(iSequence2)) {
	  
      case basic: 
	break;
      case ClpSimplex::isFixed:
	break;
      case isFree:
      case superBasic:
	// treat dj as if zero
	thetaDown=0.0;
	thetaUp=0.0;
	sequenceDown=iSequence2;
	sequenceUp=iSequence2;
	break;
      case atUpperBound:
	if (alpha>0.0) {
	  // test up
	  if (oldValue + thetaUp*alpha > dualTolerance_) {
	    thetaUp = (dualTolerance_-oldValue)/alpha;
	    sequenceUp = iSequence2;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha > dualTolerance_) {
	    thetaDown = -(dualTolerance_-oldValue)/alpha;
	    sequenceDown = iSequence2;
	  }
	}
	break;
      case atLowerBound:
	if (alpha<0.0) {
	  // test up
	  if (oldValue + thetaUp*alpha <- dualTolerance_) {
	    thetaUp = -(dualTolerance_+oldValue)/alpha;
	    sequenceUp = iSequence2;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha < -dualTolerance_) {
	    thetaDown = (dualTolerance_+oldValue)/alpha;
	    sequenceDown = iSequence2;
	  }
	}
	break;
      }
    }
  }
  if (sequenceUp>=0) {
    costIncrease = thetaUp;
    sequenceIncrease = sequenceUp;
  }
  if (sequenceDown>=0) {
    costDecrease = thetaDown;
    sequenceDecrease = sequenceDown;
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
ClpSimplexOther::primalRanging(int numberCheck,const int * which,
		  double * valueIncreased, int * sequenceIncreased,
		  double * valueDecreased, int * sequenceDecreased)
{
  rowArray_[0]->clear();
  rowArray_[1]->clear();
  lowerIn_=-COIN_DBL_MAX;
  upperIn_=COIN_DBL_MAX;
  valueIn_ = 0.0;
  for ( int i=0;i<numberCheck;i++) {
    int iSequence = which[i];
    double valueIncrease=COIN_DBL_MAX;
    double valueDecrease=COIN_DBL_MAX;
    int sequenceIncrease=-1;
    int sequenceDecrease=-1;
    
    switch(getStatus(iSequence)) {
      
    case basic:
    case isFree:
    case superBasic:
      // Easy
      valueIncrease=CoinMax(0.0,upper_[iSequence]-solution_[iSequence]);
      valueDecrease=CoinMax(0.0,solution_[iSequence]-lower_[iSequence]);
      sequenceIncrease=iSequence;
      sequenceDecrease=iSequence;
      break;
    case isFixed:
    case atUpperBound:
    case atLowerBound:
      {
	// Non trivial
	// Other bound is ignored
	unpackPacked(rowArray_[1],iSequence);
	factorization_->updateColumn(rowArray_[2],rowArray_[1]);
	// Get extra rows
	matrix_->extendUpdated(this,rowArray_[1],0);
	// do ratio test
	checkPrimalRatios(rowArray_[1],-1);
	if (pivotRow_>=0) {
	  valueIncrease = theta_;
	  sequenceIncrease=pivotVariable_[pivotRow_];
	}
	directionIn_=-1; // down
	checkPrimalRatios(rowArray_[1],1);
	if (pivotRow_>=0) {
	  valueDecrease = theta_;
	  sequenceDecrease=pivotVariable_[pivotRow_];
	}
	rowArray_[1]->clear();
      }
      break;
    }
    double scaleFactor;
    if (rowScale_) {
      if (iSequence<numberColumns_) 
	scaleFactor = columnScale_[iSequence]/rhsScale_;
      else
	scaleFactor = 1.0/(rowScale_[iSequence-numberColumns_]*rhsScale_);
    } else {
      scaleFactor = 1.0/rhsScale_;
    }
    if (valueIncrease<1.0e30)
      valueIncrease *= scaleFactor;
    else
      valueIncrease = COIN_DBL_MAX;
    if (valueDecrease<1.0e30)
      valueDecrease *= scaleFactor;
    else
      valueDecrease = COIN_DBL_MAX;
    valueIncreased[i] = valueIncrease;
    sequenceIncreased[i] = sequenceIncrease;
    valueDecreased[i] = valueDecrease;
    sequenceDecreased[i] = sequenceDecrease;
  }
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
  pivotRow_=-1;
  double acceptablePivot=1.0e-7;
  double * work=rowArray->denseVector();
  int number=rowArray->getNumElements();
  int * which=rowArray->getIndices();

  // we need to swap sign if going down
  double way = direction;
  theta_ = 1.0e30;
  for (int iIndex=0;iIndex<number;iIndex++) {
    
    int iRow = which[iIndex];
    double alpha = work[iIndex]*way;
    int iPivot=pivotVariable_[iRow];
    double oldValue = solution_[iPivot];
    if (fabs(alpha)>acceptablePivot) {
      if (alpha>0.0) {
	// basic variable going towards lower bound
	double bound = lower_[iPivot];
	oldValue -= bound;
	if (oldValue-theta_*alpha<0.0) {
	  pivotRow_ = iRow;
	  theta_ = CoinMax(0.0,oldValue/alpha);
	}
      } else {
	// basic variable going towards upper bound
	double bound = upper_[iPivot];
	oldValue = oldValue-bound;
	if (oldValue-theta_*alpha>0.0) {
	  pivotRow_ = iRow;
	  theta_ = CoinMax(0.0,oldValue/alpha);
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
  formatType=CoinMax(0,formatType);
  formatType=CoinMin(2,formatType);
  if (!writeValues)
    formatType=0;
  // See if INTEL if IEEE
  if (formatType==2) {
    // test intel here and add 1 if not intel
    double value=1.0;
    char x[8];
    memcpy(x,&value,8);
    if (x[0]==63) {
      formatType ++; // not intel
    } else {
      assert (x[0]==0);
    }
  }
  
  char number[20];
  FILE * fp = fopen(filename,"w");
  if (!fp)
    return -1;
   
  // NAME card

  if (strcmp(strParam_[ClpProbName].c_str(),"")==0) {
    fprintf(fp, "NAME          BLANK      ");
  } else {
    fprintf(fp, "NAME          %s       ",strParam_[ClpProbName].c_str());
  }
  if (formatType>=2)
    fprintf(fp,"IEEE");
  else if (writeValues)
    fprintf(fp,"VALUES");
  // finish off name 
  fprintf(fp,"\n");
  int iRow=0;
  for(int iColumn =0; iColumn < numberColumns_; iColumn++) {
    bool printit=false;
    if( getColumnStatus(iColumn) == ClpSimplex::basic) {
      printit=true;
      // Find non basic row
      for(; iRow < numberRows_; iRow++) {
	if (getRowStatus(iRow) != ClpSimplex::basic) 
	  break;
      }
      if (lengthNames_) {
	if (iRow!=numberRows_) {
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
	    fprintf(fp,"      _dummy_");
	}
      } else {
	// no names
	if (iRow!=numberRows_) {
	  fprintf(fp, " %s C%7.7d     R%7.7d",
		  getRowStatus(iRow) == ClpSimplex::atUpperBound ? "XU" : "XL",
		  iColumn,iRow);
	  iRow++;
	} else {
	  // Allow for too many basics!
	  fprintf(fp, " BS C%7.7d",iColumn);
	  // Dummy row name if values
	  if (writeValues)
	    fprintf(fp,"      _dummy_");
	}
      }
    } else  {
      if( getColumnStatus(iColumn) == ClpSimplex::atUpperBound) {
	printit=true;
	if (lengthNames_) 
	  fprintf(fp, " UL %s", columnNames_[iColumn].c_str());
	else 
	  fprintf(fp, " UL C%7.7d", iColumn);
	// Dummy row name if values
	if (writeValues)
	  fprintf(fp,"      _dummy_");
      }
    }
    if (printit&&writeValues) {
      // add value
      CoinConvertDouble(0,formatType,columnActivity_[iColumn],number);
      fprintf(fp,"     %s",number);
    }
    if (printit)
      fprintf(fp,"\n");
  }
  fprintf(fp, "ENDATA\n");
  fclose(fp);
  return 0;
}
// Read a basis from the given filename
int 
ClpSimplexOther::readBasis(const char *fileName)
{
  int status=0;
  bool canOpen=false;
  if (!strcmp(fileName,"-")||!strcmp(fileName,"stdin")) {
    // stdin
    canOpen=true;
  } else {
    FILE *fp=fopen(fileName,"r");
    if (fp) {
      // can open - lets go for it
      fclose(fp);
      canOpen=true;
    } else {
      handler_->message(CLP_UNABLE_OPEN,messages_)
	<<fileName<<CoinMessageEol;
      return -1;
    }
  }
  CoinMpsIO m;
  m.passInMessageHandler(handler_);
  bool savePrefix =m.messageHandler()->prefix();
  m.messageHandler()->setPrefix(handler_->prefix());
  status=m.readBasis(fileName,"",columnActivity_,status_+numberColumns_,
		     status_,
		     columnNames_,numberColumns_,
		     rowNames_,numberRows_);
  m.messageHandler()->setPrefix(savePrefix);
  if (status>=0) {
    if (!status) {
      // set values
      int iColumn,iRow;
      for (iRow=0;iRow<numberRows_;iRow++) {
	if (getRowStatus(iRow)==atLowerBound)
	  rowActivity_[iRow]=rowLower_[iRow];
	else if (getRowStatus(iRow)==atUpperBound)
	  rowActivity_[iRow]=rowUpper_[iRow];
      }
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	if (getColumnStatus(iColumn)==atLowerBound)
	  columnActivity_[iColumn]=columnLower_[iColumn];
	else if (getColumnStatus(iColumn)==atUpperBound)
	  columnActivity_[iColumn]=columnUpper_[iColumn];
      }
    } else {
      memset(rowActivity_,0,numberRows_*sizeof(double));
      matrix_->times(-1.0,columnActivity_,rowActivity_);
    }
  } else {
    // errors
    handler_->message(CLP_IMPORT_ERRORS,messages_)
      <<status<<fileName<<CoinMessageEol;
  }
  return status;
}
// Creates dual of a problem
ClpSimplex * 
ClpSimplexOther::dualOfModel() const
{
  const ClpSimplex * model2 = (const ClpSimplex *) this;
  bool changed=false;
  int iColumn;
  // check if we need to change bounds to rows
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (columnUpper_[iColumn]<1.0e20&&
        columnLower_[iColumn]>-1.0e20) {
      changed=true;
      break;
    }
  }
  if (changed) {
    ClpSimplex * model3 = new ClpSimplex(*model2);
    CoinBuild build;
    double one=1.0;
    int numberColumns = model3->numberColumns();
    const double * columnLower = model3->columnLower();
    const double * columnUpper = model3->columnUpper();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnUpper[iColumn]<1.0e20&&
          columnLower[iColumn]>-1.0e20) {
        if (fabs(columnLower[iColumn])<fabs(columnUpper[iColumn])) {
          double value = columnUpper[iColumn];
          model3->setColumnUpper(iColumn,COIN_DBL_MAX);
          build.addRow(1,&iColumn,&one,-COIN_DBL_MAX,value);
        } else {
          double value = columnLower[iColumn];
          model3->setColumnLower(iColumn,-COIN_DBL_MAX);
          build.addRow(1,&iColumn,&one,value,COIN_DBL_MAX);
        }
      }
    }
    model3->addRows(build);
    model2=model3;
  }
  int numberColumns = model2->numberColumns();
  const double * columnLower = model2->columnLower();
  const double * columnUpper = model2->columnUpper();
  int numberRows = model2->numberRows();
  double * rowLower = CoinCopyOfArray(model2->rowLower(),numberRows);
  double * rowUpper = CoinCopyOfArray(model2->rowUpper(),numberRows);

  const double * objective = model2->objective();
  CoinPackedMatrix * matrix = model2->matrix();
  // get transpose
  CoinPackedMatrix rowCopy = *matrix;
  int iRow;
  int numberExtraRows=0;
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]>-1.0e20&&
        rowUpper[iRow]<1.0e20) {
      if (rowUpper[iRow]!=rowLower[iRow])
         numberExtraRows++;
    }
  }
  const int * row = matrix->getIndices();
  const int * columnLength = matrix->getVectorLengths();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const double * elementByColumn = matrix->getElements();
  double objOffset=0.0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double offset=0.0;
    double objValue =optimizationDirection_*objective[iColumn];
    if (columnUpper[iColumn]>1.0e20) {
      if (columnLower[iColumn]>-1.0e20)
        offset=columnLower[iColumn];
    } else if (columnLower[iColumn]<-1.0e20) {
      offset=columnUpper[iColumn];
    } else {
      // taken care of before
      abort();
    }
    if (offset) {
      objOffset += offset*objValue;
      for (CoinBigIndex j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow = row[j];
        if (rowLower[iRow]>-1.0e20)
          rowLower[iRow] -= offset*elementByColumn[j];
        if (rowUpper[iRow]<1.0e20)
          rowUpper[iRow] -= offset*elementByColumn[j];
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
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double objValue =optimizationDirection_*objective[iColumn];
    // Offset is already in
    if (columnUpper[iColumn]>1.0e20) {
      if (columnLower[iColumn]>-1.0e20) {
        fromColumnsLower[iColumn]=-COIN_DBL_MAX;
        fromColumnsUpper[iColumn]=objValue;
      } else {
        // free
        fromColumnsLower[iColumn]=objValue;
        fromColumnsUpper[iColumn]=objValue;
      }
    } else if (columnLower[iColumn]<-1.0e20) {
      fromColumnsLower[iColumn]=objValue;
      fromColumnsUpper[iColumn]=COIN_DBL_MAX;
    } else {
      abort();
    }
  }
  int kRow=0;
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]<-1.0e20) {
      assert (rowUpper[iRow]<1.0e20);
      newObjective[kRow]=-rowUpper[iRow];
      fromRowsLower[kRow]=-COIN_DBL_MAX;
      fromRowsUpper[kRow]=0.0;
      which[kRow]=iRow;
      kRow++;
    } else if (rowUpper[iRow]>1.0e20) {
      newObjective[kRow]=-rowLower[iRow];
      fromRowsLower[kRow]=0.0;
      fromRowsUpper[kRow]=COIN_DBL_MAX;
      which[kRow]=iRow;
      kRow++;
    } else {
      if (rowUpper[iRow]==rowLower[iRow]) {
        newObjective[kRow]=-rowLower[iRow];
        fromRowsLower[kRow]=-COIN_DBL_MAX;;
        fromRowsUpper[kRow]=COIN_DBL_MAX;
        which[kRow]=iRow;
        kRow++;
      } else {
        // range
        newObjective[kRow]=-rowUpper[iRow];
        fromRowsLower[kRow]=-COIN_DBL_MAX;
        fromRowsUpper[kRow]=0.0;
        which[kRow]=iRow;
        kRow++;
        newObjective[kRow]=-rowLower[iRow];
        fromRowsLower[kRow]=0.0;
        fromRowsUpper[kRow]=COIN_DBL_MAX;
        which[kRow]=iRow;
        kRow++;
      }
    }
  }
  if (numberExtraRows) {
    CoinPackedMatrix newCopy;
    newCopy.submatrixOfWithDuplicates(rowCopy,kRow,which);
    rowCopy=newCopy;
  }
  ClpSimplex * modelDual = new ClpSimplex();
  modelDual->loadProblem(rowCopy,fromRowsLower,fromRowsUpper,newObjective,
                        fromColumnsLower,fromColumnsUpper);
  modelDual->setObjectiveOffset(objOffset);
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
void 
ClpSimplexOther::restoreFromDual(const ClpSimplex * dualProblem)
{
  createStatus();
  // Number of rows in dual problem was original number of columns
  assert (numberColumns_==dualProblem->numberRows());
  // If slack on d-row basic then column at bound otherwise column basic
  // If d-column basic then rhs tight
  int numberBasic=0;
  int iRow,iColumn=0;
  // Get number of extra rows from ranges
  int numberExtraRows=0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (rowLower_[iRow]>-1.0e20&&
        rowUpper_[iRow]<1.0e20) {
      if (rowUpper_[iRow]!=rowLower_[iRow])
         numberExtraRows++;
    }
  }
  const double * objective = this->objective();
  const double * dualDual = dualProblem->dualRowSolution();
  const double * dualDj = dualProblem->dualColumnSolution();
  const double * dualSol = dualProblem->primalColumnSolution();
  const double * dualActs = dualProblem->primalRowSolution();
#if 0
  const double * primalDual = this->dualRowSolution();
  const double * primalDj = this->dualColumnSolution();
  const double * primalSol = this->primalColumnSolution();
  const double * primalActs = this->primalRowSolution();
  char ss[]={'F','B','U','L','S','F'};
  dual(); // for testing
  printf ("Dual problem row info %d rows\n",dualProblem->numberRows());
  for (iRow=0;iRow<dualProblem->numberRows();iRow++)
    printf("%d at %c primal %g dual %g\n",
           iRow,ss[dualProblem->getRowStatus(iRow)],
           dualActs[iRow],dualDual[iRow]);
  printf ("Dual problem column info %d columns\n",dualProblem->numberColumns());
  for (iColumn=0;iColumn<dualProblem->numberColumns();iColumn++)
    printf("%d at %c primal %g dual %g\n",
           iColumn,ss[dualProblem->getColumnStatus(iColumn)],
           dualSol[iColumn],dualDj[iColumn]);
  printf ("Primal problem row info %d rows\n",this->numberRows());
  for (iRow=0;iRow<this->numberRows();iRow++)
    printf("%d at %c primal %g dual %g\n",
           iRow,ss[this->getRowStatus(iRow)],
           primalActs[iRow],primalDual[iRow]);
  printf ("Primal problem column info %d columns\n",this->numberColumns());
  for (iColumn=0;iColumn<this->numberColumns();iColumn++)
    printf("%d at %c primal %g dual %g\n",
           iColumn,ss[this->getColumnStatus(iColumn)],
           primalSol[iColumn],primalDj[iColumn]);
#endif
  // position at bound information
  int jColumn=numberRows_+numberExtraRows;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double objValue =optimizationDirection_*objective[iColumn];
    Status status = dualProblem->getRowStatus(iColumn);
    double otherValue = COIN_DBL_MAX;
    if (columnUpper_[iColumn]<1.0e20&&
        columnLower_[iColumn]>-1.0e20) {
        if (fabs(columnLower_[iColumn])<fabs(columnUpper_[iColumn])) {
          otherValue = columnUpper_[iColumn]+dualDj[jColumn];
        } else {
          otherValue = columnLower_[iColumn]+dualDj[jColumn];
        }
        jColumn++;
    }
    if (status==basic) {
      // column is at bound
      if (otherValue==COIN_DBL_MAX) {
        reducedCost_[iColumn]=objValue-dualActs[iColumn];
        if (columnUpper_[iColumn]>1.0e20) {
          setColumnStatus(iColumn,atLowerBound);
          columnActivity_[iColumn]=columnLower_[iColumn];
        } else {
          setColumnStatus(iColumn,atUpperBound);
          columnActivity_[iColumn]=columnUpper_[iColumn];
        }
      } else {
        reducedCost_[iColumn]=objValue-dualActs[iColumn];
        //printf("other dual sol %g\n",otherValue);
        if (fabs(otherValue-columnLower_[iColumn])<1.0e-5) {
          setColumnStatus(iColumn,atLowerBound);
          columnActivity_[iColumn]=columnLower_[iColumn];
        } else if (fabs(otherValue-columnUpper_[iColumn])<1.0e-5) {
          setColumnStatus(iColumn,atUpperBound);
          columnActivity_[iColumn]=columnUpper_[iColumn];
        } else {
          abort();
        }
      }
    } else {
      if (otherValue==COIN_DBL_MAX) {
        // column basic
        setColumnStatus(iColumn,basic);
        numberBasic++;
        columnActivity_[iColumn]=-dualDual[iColumn];
        reducedCost_[iColumn]=0.0;
      } else {
        // may be at other bound
        //printf("xx %d %g jcol %d\n",iColumn,otherValue,jColumn-1);
        if (dualProblem->getColumnStatus(jColumn-1)!=basic) {
          // column basic
          setColumnStatus(iColumn,basic);
          numberBasic++;
          columnActivity_[iColumn]=-dualDual[iColumn];
          reducedCost_[iColumn]=0.0;
        } else {
          reducedCost_[iColumn]=objValue-dualActs[iColumn];
          if (fabs(otherValue-columnLower_[iColumn])<1.0e-5) {
            setColumnStatus(iColumn,atLowerBound);
            columnActivity_[iColumn]=columnLower_[iColumn];
          } else if (fabs(otherValue-columnUpper_[iColumn])<1.0e-5) {
            setColumnStatus(iColumn,atUpperBound);
            columnActivity_[iColumn]=columnUpper_[iColumn];
          } else {
            abort();
          }
        }
      }
    }
  }
  // now rows
  int kRow=0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    Status status = dualProblem->getColumnStatus(kRow);
    if (status==basic) {
      // row is at bound
      dual_[iRow]=dualSol[kRow];;
    } else {
      // row basic
      setRowStatus(iRow,basic);
      numberBasic++;
      dual_[iRow]=0.0;
    }
    if (rowLower_[iRow]<-1.0e20) {
      if (status==basic) {
        rowActivity_[iRow]=rowUpper_[iRow];
        setRowStatus(iRow,atUpperBound);
      } else {
        rowActivity_[iRow]=rowUpper_[iRow]+dualSol[kRow];
      }        
      kRow++;
    } else if (rowUpper_[iRow]>1.0e20) {
      if (status==basic) {
        rowActivity_[iRow]=rowLower_[iRow];
        setRowStatus(iRow,atLowerBound);
      } else {
        rowActivity_[iRow]=rowLower_[iRow]+dualSol[kRow];
      }        
      kRow++;
    } else {
      if (rowUpper_[iRow]==rowLower_[iRow]) {
        rowActivity_[iRow]=rowLower_[iRow];
        if (status==basic) {
          setRowStatus(iRow,atLowerBound);
        }        
        kRow++;
      } else {
        // range
        abort();
        kRow++;
        kRow++;
      }
    }
  }
  assert (numberBasic==numberRows_);
  if (optimizationDirection_<0.0) {
    for (iRow=0;iRow<numberRows_;iRow++) {
      dual_[iRow]=-dual_[iRow];
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      reducedCost_[iColumn]=-reducedCost_[iColumn];
    }
  }
  checkSolutionInternal();
}
/* Does very cursory presolve.
   rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns
*/
ClpSimplex * 
ClpSimplexOther::crunch(double * rhs, int * whichRow, int * whichColumn,
                        int & nBound, bool moreBounds)
{
#if 0
  //#ifndef NDEBUG
  {
    int n=0;
    int i;
    for (i=0;i<numberColumns_;i++)
      if (getColumnStatus(i)==ClpSimplex::basic)
        n++;
    for (i=0;i<numberRows_;i++)
      if (getRowStatus(i)==ClpSimplex::basic)
        n++;
    assert (n==numberRows_);
  }
#endif
  
  const double * element = matrix_->getElements();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths();

  CoinZeroN(rhs,numberRows_);
  int iColumn;
  int iRow;
  CoinZeroN(whichRow,numberRows_);
  int * backColumn = whichColumn+numberColumns_;
  int numberRows2=0;
  int numberColumns2=0;
  double offset=0.0;
  const double * objective = this->objective();
  double * solution = columnActivity_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double lower = columnLower_[iColumn];
    double upper = columnUpper_[iColumn];
    if (upper>lower||getColumnStatus(iColumn)==ClpSimplex::basic) {
      backColumn[iColumn]=numberColumns2;
      whichColumn[numberColumns2++]=iColumn;
      for (CoinBigIndex j = columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow = row[j];
        int n=whichRow[iRow];
        if (n==0)
          whichRow[iRow]=-iColumn-1;
        else if (n<0)
          whichRow[iRow]=2;
      }
    } else {
      // fixed
      backColumn[iColumn]=-1;
      solution[iColumn]=upper;
      if (upper) {
        offset += objective[iColumn]*upper;
        for (CoinBigIndex j = columnStart[iColumn];
             j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          double value = element[j];
          rhs[iRow] += upper*value;
        }
      }
    }
  }
  int returnCode=0;
  double tolerance = primalTolerance();
  nBound=2*numberRows_;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int n=whichRow[iRow];
    if (n>0) {
      whichRow[numberRows2++]=iRow;
    } else if (n<0) {
      //whichRow[numberRows2++]=iRow;
      //continue;
      // Can only do in certain circumstances as we don't know current value
      if (rowLower_[iRow]==rowUpper_[iRow]||getRowStatus(iRow)==ClpSimplex::basic) {
        // save row and column for bound
        whichRow[--nBound]=iRow;
        whichRow[nBound+numberRows_]=-n-1;
      } else if (moreBounds) {
        // save row and column for bound
        whichRow[--nBound]=iRow;
        whichRow[nBound+numberRows_]=-n-1;
      } else {
        whichRow[numberRows2++]=iRow;
      }
    } else {
      // empty
      double rhsValue = rhs[iRow];
      if (rhsValue<rowLower_[iRow]-tolerance||rhsValue>rowUpper_[iRow]+tolerance) {
        returnCode=1; // infeasible
      }
    }
  }
  ClpSimplex * small=NULL;
  if (!returnCode) {
    small = new ClpSimplex(this,numberRows2,whichRow,
                     numberColumns2,whichColumn);
    int numberElements=getNumElements();
    int numberElements2=small->getNumElements();
    small->setObjectiveOffset(objectiveOffset()-offset);
    handler_->message(CLP_CRUNCH_STATS,messages_)
      <<numberRows2<< -(numberRows_ - numberRows2)
      <<numberColumns2<< -(numberColumns_ - numberColumns2)
      <<numberElements2<< -(numberElements - numberElements2)
      <<CoinMessageEol;
    // And set objective value to match
    small->setObjectiveValue(this->objectiveValue());
    double * rowLower2 = small->rowLower();
    double * rowUpper2 = small->rowUpper();
    int jRow;
    for (jRow=0;jRow<numberRows2;jRow++) {
      iRow = whichRow[jRow];
      if (rowLower2[jRow]>-1.0e20)
        rowLower2[jRow] -= rhs[iRow];
      if (rowUpper2[jRow]<1.0e20)
        rowUpper2[jRow] -= rhs[iRow];
    }
    // and bounds
    double * columnLower2 = small->columnLower();
    double * columnUpper2 = small->columnUpper();
    for (jRow=nBound;jRow<2*numberRows_;jRow++) {
      iRow = whichRow[jRow];
      iColumn = whichRow[jRow+numberRows_];
      double lowerRow = rowLower_[iRow];
      if (lowerRow>-1.0e20)
        lowerRow -= rhs[iRow];
      double upperRow = rowUpper_[iRow];
      if (upperRow<1.0e20)
        upperRow -= rhs[iRow];
      int jColumn = backColumn[iColumn];
      double lower = columnLower2[jColumn];
      double upper = columnUpper2[jColumn];
      double value=0.0;
      for (CoinBigIndex j = columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        if (iRow==row[j]) {
          value=element[j];
          break;
        }
      }
      assert (value);
      // convert rowLower and Upper to implied bounds on column
      double newLower=-COIN_DBL_MAX;
      double newUpper=COIN_DBL_MAX;
      if (value>0.0) {
        if (lowerRow>-1.0e20)
          newLower = lowerRow/value;
        if (upperRow<1.0e20)
          newUpper = upperRow/value;
      } else {
        if (upperRow<1.0e20)
          newLower = upperRow/value;
        if (lowerRow>-1.0e20)
          newUpper = lowerRow/value;
      }
      newLower = CoinMax(lower,newLower);
      newUpper = CoinMin(upper,newUpper);
      if (newLower>newUpper+tolerance) {
        //printf("XXYY inf on bound\n");
        returnCode=1;
      }
      columnLower2[jColumn]=newLower;
      columnUpper2[jColumn]=CoinMax(newLower,newUpper);
      if (getRowStatus(iRow)!=ClpSimplex::basic) {
        if (getColumnStatus(iColumn)==ClpSimplex::basic) {
          if (columnLower2[jColumn]==columnUpper2[jColumn]) {
            // can only get here if will be fixed
            small->setColumnStatus(jColumn,ClpSimplex::isFixed);
          } else {
            // solution is valid
            if (fabs(columnActivity_[iColumn]-columnLower2[jColumn])<
                fabs(columnActivity_[iColumn]-columnUpper2[jColumn]))
              small->setColumnStatus(jColumn,ClpSimplex::atLowerBound);
            else
              small->setColumnStatus(jColumn,ClpSimplex::atUpperBound);
          }
        } else {
          printf("what now neither basic\n");
        }
      }
    }
    if (returnCode) {
      delete small;
      small=NULL;
    }
  }
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
  getbackSolution(small,whichRow,whichColumn);
  // and deal with status for bounds
  const double * element = matrix_->getElements();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths();
  double tolerance = primalTolerance();
  double djTolerance = dualTolerance();
  for (int jRow=nBound;jRow<2*numberRows_;jRow++) {
    int iRow = whichRow[jRow];
    int iColumn = whichRow[jRow+numberRows_];
    if (getColumnStatus(iColumn)!=ClpSimplex::basic) {
      double lower = columnLower_[iColumn];
      double upper = columnUpper_[iColumn];
      double value = columnActivity_[iColumn];
      double djValue = reducedCost_[iColumn];
      dual_[iRow]=0.0;
      if (upper>lower) {
        if (value<lower+tolerance&&djValue>-djTolerance) {
          setColumnStatus(iColumn,ClpSimplex::atLowerBound);
          setRowStatus(iRow,ClpSimplex::basic);
        } else if (value>upper-tolerance&&djValue<djTolerance) {
          setColumnStatus(iColumn,ClpSimplex::atUpperBound);
          setRowStatus(iRow,ClpSimplex::basic);
        } else {
          // has to be basic
          setColumnStatus(iColumn,ClpSimplex::basic);
          double value=0.0;
          for (CoinBigIndex j = columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            if (iRow==row[j]) {
              value=element[j];
              break;
            }
          }
          dual_[iRow]=djValue/value;
          if (rowUpper_[iRow]>rowLower_[iRow]) {
            if (fabs(rowActivity_[iRow]-rowLower_[iRow])<
                fabs(rowActivity_[iRow]-rowUpper_[iRow]))
              setRowStatus(iRow,ClpSimplex::atLowerBound);
            else
              setRowStatus(iRow,ClpSimplex::atUpperBound);
          } else {
            setRowStatus(iRow,ClpSimplex::isFixed);
          }
        }
      } else {
        // row can always be basic
        setRowStatus(iRow,ClpSimplex::basic);
      }
    } else {
      // row can always be basic
      setRowStatus(iRow,ClpSimplex::basic);
    }
  }
  //#ifndef NDEBUG
#if 0
  if  (small.status()==0) {
    int n=0;
    int i;
    for (i=0;i<numberColumns;i++)
      if (getColumnStatus(i)==ClpSimplex::basic)
        n++;
    for (i=0;i<numberRows;i++)
      if (getRowStatus(i)==ClpSimplex::basic)
        n++;
    assert (n==numberRows);
  }
#endif
}
