// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpFactorization.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
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
      CoinConvertDouble(formatType,columnActivity_[iColumn],number);
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
