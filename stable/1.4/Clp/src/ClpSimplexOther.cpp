// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

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
				  double * costDecreased, int * sequenceDecreased,
				  double * valueIncrease, double * valueDecrease)
{
  rowArray_[1]->clear();
  columnArray_[1]->clear();
  // long enough for rows+columns
  assert(rowArray_[3]->capacity()>=numberRows_+numberColumns_);
  rowArray_[3]->clear();
  int * backPivot = rowArray_[3]->getIndices();
  int i;
  for ( i=0;i<numberRows_+numberColumns_;i++) {
    backPivot[i]=-1;
  }
  for (i=0;i<numberRows_;i++) {
    int iSequence = pivotVariable_[i];
    backPivot[iSequence]=i;
  }
  // dualTolerance may be zero if from CBC.  In fact use that fact
  bool inCBC = !dualTolerance_;
  if (inCBC)
    assert (integerType_);
  dualTolerance_ = dblParam_[ClpDualTolerance];
  for ( i=0;i<numberCheck;i++) {
    rowArray_[0]->clear();
    //rowArray_[0]->checkClear();
    //rowArray_[1]->checkClear();
    //columnArray_[1]->checkClear();
    columnArray_[0]->clear();
    //columnArray_[0]->checkClear();
    int iSequence = which[i];
    double costIncrease=COIN_DBL_MAX;
    double costDecrease=COIN_DBL_MAX;
    int sequenceIncrease=-1;
    int sequenceDecrease=-1;
    if (valueIncrease) {
      assert (valueDecrease);
      valueIncrease[i]=iSequence<numberColumns_ ? columnActivity_[iSequence] : rowActivity_[iSequence-numberColumns_];
      valueDecrease[i]=valueIncrease[i];
    }
    
    switch(getStatus(iSequence)) {
      
    case basic:
      {
	// non-trvial
	// Get pivot row
	int iRow=backPivot[iSequence];
	assert (iRow>=0);
	double plusOne=1.0;
        rowArray_[0]->createPacked(1,&iRow,&plusOne);
	factorization_->updateColumnTranspose(rowArray_[1],rowArray_[0]);
	// put row of tableau in rowArray[0] and columnArray[0]
	matrix_->transposeTimes(this,-1.0,
				rowArray_[0],columnArray_[1],columnArray_[0]);
        double alphaIncrease;
        double alphaDecrease;
	// do ratio test up and down
	checkDualRatios(rowArray_[0],columnArray_[0],costIncrease,sequenceIncrease,alphaIncrease,
		    costDecrease,sequenceDecrease,alphaDecrease);
	if (valueIncrease) {
	  if (sequenceIncrease>=0)
	    valueIncrease[i] = primalRanging1(sequenceIncrease,iSequence);
	  if (sequenceDecrease>=0)
	    valueDecrease[i] = primalRanging1(sequenceDecrease,iSequence);
	}
        if (inCBC) { 
          if (sequenceIncrease>=0) {
            double djValue = dj_[sequenceIncrease];
            if (fabs(djValue)>10.0*dualTolerance_) {
              // we are going to use for cutoff so be exact
              costIncrease = fabs(djValue/alphaIncrease); 
              /* Not sure this is good idea as I don't think correct e.g.
                 suppose a continuous variable has dj slightly greater. */
              if(false&&sequenceIncrease<numberColumns_&&integerType_[sequenceIncrease]) {
                // can improve
                double movement = (columnScale_==NULL) ? 1.0 : 
                  rhsScale_/columnScale_[sequenceIncrease];
                costIncrease = CoinMax(fabs(djValue*movement),costIncrease);
              }
            } else {
              costIncrease=0.0;
            }
          }
          if (sequenceDecrease>=0) {
            double djValue = dj_[sequenceDecrease];
            if (fabs(djValue)>10.0*dualTolerance_) {
              // we are going to use for cutoff so be exact
              costDecrease = fabs(djValue/alphaDecrease); 
              if(sequenceDecrease<numberColumns_&&integerType_[sequenceDecrease]) {
                // can improve
                double movement = (columnScale_==NULL) ? 1.0 : 
                  rhsScale_/columnScale_[sequenceDecrease];
                costDecrease = CoinMax(fabs(djValue*movement),costDecrease);
              }
            } else {
              costDecrease=0.0;
            }
          }
        }
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
      if (valueIncrease) 
	valueIncrease[i] = primalRanging1(iSequence,iSequence);
      break;
    case atLowerBound:
      costDecrease = CoinMax(0.0,dj_[iSequence]);
      sequenceDecrease = iSequence;
      if (valueIncrease) 
	valueDecrease[i] = primalRanging1(iSequence,iSequence);
      break;
    }
    double scaleFactor;
    if (!auxiliaryModel_) {
      if (rowScale_) {
        if (iSequence<numberColumns_) 
          scaleFactor = 1.0/(objectiveScale_*columnScale_[iSequence]);
        else
          scaleFactor = rowScale_[iSequence-numberColumns_]/objectiveScale_;
      } else {
        scaleFactor = 1.0/objectiveScale_;
      }
    } else {
      if (auxiliaryModel_->rowScale()) {
        if (iSequence<numberColumns_) 
          scaleFactor = 1.0/(objectiveScale_*auxiliaryModel_->columnScale()[iSequence]);
        else
          scaleFactor = auxiliaryModel_->rowScale()[iSequence-numberColumns_]/objectiveScale_;
      } else {
        scaleFactor = 1.0/objectiveScale_;
      }
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
      if (valueIncrease) {
	double temp = valueIncrease[i];
	valueIncrease[i]=valueDecrease[i];
	valueDecrease[i]=temp;
      }
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
  //rowArray_[0]->clear();
  //rowArray_[1]->clear();
  //columnArray_[1]->clear();
  //columnArray_[0]->clear();
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
  int sequenceDown =-1;
  int sequenceUp = -1;
  double alphaDown=0.0;
  double alphaUp=0.0;

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
            alphaUp=alpha;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha > dualTolerance_) {
	    thetaDown = -(dualTolerance_-oldValue)/alpha;
	    sequenceDown = iSequence2;
            alphaDown=alpha;
	  }
	}
	break;
      case atLowerBound:
	if (alpha<0.0) {
	  // test up
	  if (oldValue + thetaUp*alpha <- dualTolerance_) {
	    thetaUp = -(dualTolerance_+oldValue)/alpha;
	    sequenceUp = iSequence2;
            alphaUp=alpha;
	  }
	} else {
	  // test down
	  if (oldValue - thetaDown*alpha < -dualTolerance_) {
	    thetaDown = (dualTolerance_+oldValue)/alpha;
	    sequenceDown = iSequence2;
            alphaDown=alpha;
	  }
	}
	break;
      }
    }
  }
  if (sequenceUp>=0) {
    costIncrease = thetaUp;
    sequenceIncrease = sequenceUp;
    alphaIncrease = alphaUp;
  }
  if (sequenceDown>=0) {
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
	checkPrimalRatios(rowArray_[1],1);
	if (pivotRow_>=0) {
	  valueIncrease = theta_;
	  sequenceIncrease=pivotVariable_[pivotRow_];
	}
	checkPrimalRatios(rowArray_[1],-1);
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
// Returns new value of whichOther when whichIn enters basis
double 
ClpSimplexOther::primalRanging1(int whichIn, int whichOther)
{
  rowArray_[0]->clear();
  rowArray_[1]->clear();
  int iSequence = whichIn;
  double newValue=solution_[whichOther];
  double alphaOther=0.0;
  Status status = getStatus(iSequence);
  assert (status==atLowerBound||status==atUpperBound);
  int wayIn = (status==atLowerBound) ? 1 : -1;
  
  switch(getStatus(iSequence)) {
    
  case basic:
  case isFree:
  case superBasic:
    assert (whichIn==whichOther);
    // Easy
    newValue = wayIn>0 ? upper_[iSequence] : lower_[iSequence];
    break;
  case isFixed:
  case atUpperBound:
  case atLowerBound:
    // Non trivial
    {
      // Other bound is ignored
      unpackPacked(rowArray_[1],iSequence);
      factorization_->updateColumn(rowArray_[2],rowArray_[1]);
      // Get extra rows
      matrix_->extendUpdated(this,rowArray_[1],0);
      // do ratio test
      double acceptablePivot=1.0e-7;
      double * work=rowArray_[1]->denseVector();
      int number=rowArray_[1]->getNumElements();
      int * which=rowArray_[1]->getIndices();
      
      // we may need to swap sign
      double way = wayIn;
      double theta = 1.0e30;
      for (int iIndex=0;iIndex<number;iIndex++) {
	
	int iRow = which[iIndex];
	double alpha = work[iIndex]*way;
	int iPivot=pivotVariable_[iRow];
	if (iPivot==whichOther) {
	  alphaOther=alpha;
	  continue;
	}
	double oldValue = solution_[iPivot];
	if (fabs(alpha)>acceptablePivot) {
	  if (alpha>0.0) {
	    // basic variable going towards lower bound
	    double bound = lower_[iPivot];
	    oldValue -= bound;
	    if (oldValue-theta*alpha<0.0) {
	      theta = CoinMax(0.0,oldValue/alpha);
	    }
	  } else {
	    // basic variable going towards upper bound
	    double bound = upper_[iPivot];
	    oldValue = oldValue-bound;
	    if (oldValue-theta*alpha>0.0) {
	      theta = CoinMax(0.0,oldValue/alpha);
	    }
	  }
	}
      }
      if (whichIn!=whichOther) {
	if (theta<1.0e30)
	  newValue -= theta*alphaOther;
	else
	  newValue = alphaOther>0.0 ? -1.0e30 : 1.0e30;
      } else {
	newValue += theta*wayIn;
      }
    }
    rowArray_[1]->clear();
    break;
  }
  double scaleFactor;
  if (rowScale_) {
    if (whichOther<numberColumns_) 
      scaleFactor = columnScale_[whichOther]/rhsScale_;
    else
      scaleFactor = 1.0/(rowScale_[whichOther-numberColumns_]*rhsScale_);
  } else {
    scaleFactor = 1.0/rhsScale_;
  }
  if (newValue<1.0e29)
    if (newValue>-1.0e29)
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
  *m.messagesPointer()=coinMessages();
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
  int returnCode=0;;
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
	if (columnLower_[iColumn]>-1.0e20) {
	  columnActivity_[iColumn]=-dualDual[iColumn] + columnLower_[iColumn];
	} else if (columnUpper_[iColumn]<1.0e20) {
	  columnActivity_[iColumn]=-dualDual[iColumn] + columnUpper_[iColumn];
	}
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
  int numberRanges=0;
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
        numberRanges++;
	Status statusL = dualProblem->getColumnStatus(kRow+1);
	if (status==basic) {
	  assert (statusL!=basic);
	  rowActivity_[iRow]=rowUpper_[iRow];
	  setRowStatus(iRow,atUpperBound);
	} else if (statusL==basic) {
	  rowActivity_[iRow]=rowLower_[iRow];
	  setRowStatus(iRow,atLowerBound);
	} else {
	  rowActivity_[iRow]=rowLower_[iRow]+dualSol[kRow];
	  // row basic
	  setRowStatus(iRow,basic);
	  numberBasic++;
	  dual_[iRow]=0.0;
	}
        kRow++;
        kRow++;
      }
    }
  }
  if (numberRanges) {
    printf("%d ranges - coding needed\n",numberRanges);
    returnCode=1;
  }
  if (numberBasic!=numberRows_) {
    printf("Bad basis - ranges?\n");
    assert (numberRanges);
  }
  if (optimizationDirection_<0.0) {
    for (iRow=0;iRow<numberRows_;iRow++) {
      dual_[iRow]=-dual_[iRow];
    }
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      reducedCost_[iColumn]=-reducedCost_[iColumn];
    }
  }
  // redo row activities
  memset(rowActivity_,0,numberRows_*sizeof(double));
  matrix_->times(-1.0,columnActivity_,rowActivity_);
  checkSolutionInternal();
  //return returnCode;
}
/* Does very cursory presolve.
   rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns
*/
ClpSimplex * 
ClpSimplexOther::crunch(double * rhs, int * whichRow, int * whichColumn,
                        int & nBound, bool moreBounds, bool tightenBounds)
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
        if (n==0&&element[j])
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
                     numberColumns2,whichColumn,true,false);
    // Set some stuff
    small->setDualBound(dualBound_);
    small->setInfeasibilityCost(infeasibilityCost_);
    small->setSpecialOptions(specialOptions_);
    small->setPerturbation(perturbation_);
    small->defaultFactorizationFrequency();
    // If no rows left then no tightening!
    if (!numberRows2||!numberColumns2) 
      tightenBounds=false;

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
    const char * integerInformation = integerType_;
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
      if (integerInformation&&integerInformation[iColumn]) {
        if (newLower-floor(newLower)<10.0*tolerance) 
          newLower=floor(newLower);
        else
          newLower=ceil(newLower);
        if (ceil(newUpper)-newUpper<10.0*tolerance) 
          newUpper=ceil(newUpper);
        else
          newUpper=floor(newUpper);
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
          //printf("what now neither basic\n");
        }
      }
    }
    if (returnCode) {
      delete small;
      small=NULL;
    } else if (tightenBounds&&integerInformation) {
      // See if we can tighten any bounds
      // use rhs for upper and small duals for lower
      double * up = rhs;
      double * lo = small->dualRowSolution();
      const double * element = small->clpMatrix()->getElements();
      const int * row = small->clpMatrix()->getIndices();
      const CoinBigIndex * columnStart = small->clpMatrix()->getVectorStarts();
      //const int * columnLength = small->clpMatrix()->getVectorLengths();
      CoinZeroN(lo,numberRows2);
      CoinZeroN(up,numberRows2);
      for (int iColumn=0;iColumn<numberColumns2;iColumn++) {
        double upper=columnUpper2[iColumn];
        double lower=columnLower2[iColumn];
        //assert (columnLength[iColumn]==columnStart[iColumn+1]-columnStart[iColumn]);
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn+1];j++) {
          int iRow=row[j];
          double value = element[j];
          if (value>0.0) {
            if (upper<1.0e20)
              up[iRow] += upper*value;
            else
              up[iRow] = COIN_DBL_MAX;
            if (lower>-1.0e20)
              lo[iRow] += lower*value;
            else
              lo[iRow] = -COIN_DBL_MAX;
          } else {
            if (upper<1.0e20)
              lo[iRow] += upper*value;
            else
              lo[iRow] = -COIN_DBL_MAX;
            if (lower>-1.0e20)
              up[iRow] += lower*value;
            else
              up[iRow] = COIN_DBL_MAX;
          }
        }
      }
      double * rowLower2 = small->rowLower();
      double * rowUpper2 = small->rowUpper();
      bool feasible=true;
      // make safer
      for (int iRow=0;iRow<numberRows2;iRow++) {
        double lower = lo[iRow];
        if (lower>rowUpper2[iRow]+tolerance) {
          feasible=false;
          break;
        } else {
          lo[iRow] = CoinMin(lower-rowUpper2[iRow],0.0)-tolerance;
        }
        double upper = up[iRow];
        if (upper<rowLower2[iRow]-tolerance) {
          feasible=false;
          break;
        } else {
          up[iRow] = CoinMax(upper-rowLower2[iRow],0.0)+tolerance;
        }
      }
      if (!feasible) {
        delete small;
        small=NULL;
      } else {
        // and tighten
        for (int iColumn=0;iColumn<numberColumns2;iColumn++) {
          if (integerInformation[whichColumn[iColumn]]) {
            double upper=columnUpper2[iColumn];
            double lower=columnLower2[iColumn];
            double newUpper = upper;
            double newLower = lower;
            double difference = upper-lower;
            if (lower>-1000.0&&upper<1000.0) {
              for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn+1];j++) {
                int iRow=row[j];
                double value = element[j];
                if (value>0.0) {
                  double upWithOut = up[iRow] - value*difference;
                  if (upWithOut<0.0) {
                    newLower = CoinMax(newLower,lower-(upWithOut+tolerance)/value);
                  }
                  double lowWithOut = lo[iRow] + value*difference;
                  if (lowWithOut>0.0) {
                    newUpper = CoinMin(newUpper,upper-(lowWithOut-tolerance)/value);
                  }
                } else {
                  double upWithOut = up[iRow] + value*difference;
                  if (upWithOut<0.0) {
                    newUpper = CoinMin(newUpper,upper-(upWithOut+tolerance)/value);
                  }
                  double lowWithOut = lo[iRow] - value*difference;
                  if (lowWithOut>0.0) {
                    newLower = CoinMax(newLower,lower-(lowWithOut-tolerance)/value);
                  }
                }
              }
              if (newLower>lower||newUpper<upper) {
                if (fabs(newUpper-floor(newUpper+0.5))>1.0e-6)
                  newUpper = floor(newUpper);
                else
                  newUpper = floor(newUpper+0.5);
                if (fabs(newLower-ceil(newLower-0.5))>1.0e-6)
                  newLower = ceil(newLower);
                else
                  newLower = ceil(newLower-0.5);
                // change may be too small - check
                if (newLower>lower||newUpper<upper) {
                  if (newUpper>=newLower) {
                    // Could also tighten in this
                    //printf("%d bounds %g %g tightened to %g %g\n",
                    //     iColumn,columnLower2[iColumn],columnUpper2[iColumn],
                    //     newLower,newUpper);
#if 1
                    columnUpper2[iColumn]=newUpper;
                    columnLower2[iColumn]=newLower;
                    columnUpper_[whichColumn[iColumn]]=newUpper;
                    columnLower_[whichColumn[iColumn]]=newLower;
#endif
                    // and adjust bounds on rows
                    newUpper -= upper;
                    newLower -= lower;
                    for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn+1];j++) {
                      int iRow=row[j];
                      double value = element[j];
                      if (value>0.0) {
                        up[iRow] += newUpper*value;
                        lo[iRow] += newLower*value;
                      } else {
                        lo[iRow] += newUpper*value;
                        up[iRow] += newLower*value;
                      }
                    }
                  } else {
                    // infeasible
                    //printf("%d bounds infeasible %g %g tightened to %g %g\n",
                    //     iColumn,columnLower2[iColumn],columnUpper2[iColumn],
                    //     newLower,newUpper);
#if 1
                    delete small;
                    small=NULL;
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
          reducedCost_[iColumn] = 0.0;
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
  CoinZeroN(lo,numberRows_);
  CoinZeroN(up,numberRows_);
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    double upper=columnUpper_[iColumn];
    double lower=columnLower_[iColumn];
    //assert (columnLength[iColumn]==columnStart[iColumn+1]-columnStart[iColumn]);
    for (CoinBigIndex j=columnStart[iColumn];
         j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow=row[j];
      double value = element[j];
      if (value>0.0) {
        if (upper<1.0e20)
          up[iRow] += upper*value;
        else
          up[iRow] = COIN_DBL_MAX;
        if (lower>-1.0e20)
          lo[iRow] += lower*value;
        else
          lo[iRow] = -COIN_DBL_MAX;
      } else {
        if (upper<1.0e20)
          lo[iRow] += upper*value;
        else
          lo[iRow] = -COIN_DBL_MAX;
        if (lower>-1.0e20)
          up[iRow] += lower*value;
        else
          up[iRow] = COIN_DBL_MAX;
      }
    }
  }
  bool feasible=true;
  // make safer
  double tolerance = primalTolerance();
  for (int iRow=0;iRow<numberRows_;iRow++) {
    double lower = lo[iRow];
    if (lower>rowUpper_[iRow]+tolerance) {
      feasible=false;
      break;
    } else {
      lo[iRow] = CoinMin(lower-rowUpper_[iRow],0.0)-tolerance;
    }
    double upper = up[iRow];
    if (upper<rowLower_[iRow]-tolerance) {
      feasible=false;
      break;
    } else {
      up[iRow] = CoinMax(upper-rowLower_[iRow],0.0)+tolerance;
    }
  }
  int numberTightened=0;
  if (!feasible) {
    return -1;
  } else {
    // and tighten
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (integerType_[iColumn]) {
        double upper=columnUpper_[iColumn];
        double lower=columnLower_[iColumn];
        double newUpper = upper;
        double newLower = lower;
        double difference = upper-lower;
        if (lower>-1000.0&&upper<1000.0) {
          for (CoinBigIndex j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            double value = element[j];
            if (value>0.0) {
              double upWithOut = up[iRow] - value*difference;
              if (upWithOut<0.0) {
                newLower = CoinMax(newLower,lower-(upWithOut+tolerance)/value);
              }
              double lowWithOut = lo[iRow] + value*difference;
              if (lowWithOut>0.0) {
                newUpper = CoinMin(newUpper,upper-(lowWithOut-tolerance)/value);
              }
            } else {
              double upWithOut = up[iRow] + value*difference;
              if (upWithOut<0.0) {
                newUpper = CoinMin(newUpper,upper-(upWithOut+tolerance)/value);
              }
              double lowWithOut = lo[iRow] - value*difference;
              if (lowWithOut>0.0) {
                newLower = CoinMax(newLower,lower-(lowWithOut-tolerance)/value);
              }
            }
          }
          if (newLower>lower||newUpper<upper) {
            if (fabs(newUpper-floor(newUpper+0.5))>1.0e-6)
              newUpper = floor(newUpper);
            else
              newUpper = floor(newUpper+0.5);
            if (fabs(newLower-ceil(newLower-0.5))>1.0e-6)
              newLower = ceil(newLower);
            else
              newLower = ceil(newLower-0.5);
            // change may be too small - check
            if (newLower>lower||newUpper<upper) {
              if (newUpper>=newLower) {
                numberTightened++;
                //printf("%d bounds %g %g tightened to %g %g\n",
                //     iColumn,columnLower_[iColumn],columnUpper_[iColumn],
                //     newLower,newUpper);
                columnUpper_[iColumn]=newUpper;
                columnLower_[iColumn]=newLower;
                // and adjust bounds on rows
                newUpper -= upper;
                newLower -= lower;
                for (CoinBigIndex j=columnStart[iColumn];
                     j<columnStart[iColumn]+columnLength[iColumn];j++) {
                  int iRow=row[j];
                  double value = element[j];
                  if (value>0.0) {
                    up[iRow] += newUpper*value;
                    lo[iRow] += newLower*value;
                  } else {
                    lo[iRow] += newUpper*value;
                    up[iRow] += newLower*value;
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
ClpSimplexOther::parametrics(double startingTheta, double & endingTheta,double reportIncrement,
                             const double * changeLowerBound, const double * changeUpperBound,
                             const double * changeLowerRhs, const double * changeUpperRhs,
                             const double * changeObjective)
{
  bool needToDoSomething=true;
  bool canTryQuick = (reportIncrement) ? true : false;
  // Save copy of model
  ClpSimplex copyModel = *this;
  int savePerturbation = perturbation_;
  perturbation_=102; // switch off
  while (needToDoSomething) {
    needToDoSomething=false;
    algorithm_ = -1;
    
    // save data
    ClpDataSave data = saveData();
    int returnCode = ((ClpSimplexDual *) this)->startupSolve(0,NULL,0);
    int iRow,iColumn;
    double * chgUpper=NULL;
    double * chgLower=NULL;
    double * chgObjective=NULL;
    
    // Dantzig (as will not be used) (out later)
    ClpDualRowPivot * savePivot = dualRowPivot_;
    dualRowPivot_ = new ClpDualRowDantzig();
    
    if (!returnCode) {
      // Find theta when bounds will cross over and create arrays
      int numberTotal = numberRows_+numberColumns_;
      chgLower = new double[numberTotal];
      memset(chgLower,0,numberTotal*sizeof(double));
      chgUpper = new double[numberTotal];
      memset(chgUpper,0,numberTotal*sizeof(double));
      chgObjective = new double[numberTotal];
      memset(chgObjective,0,numberTotal*sizeof(double));
      assert (!rowScale_);
      double maxTheta=1.0e50;
      if (changeLowerRhs||changeUpperRhs) {
        for (iRow=0;iRow<numberRows_;iRow++) {
          double lower = rowLower_[iRow];
          double upper = rowUpper_[iRow];
          if (lower>upper) {
            maxTheta=-1.0;
            break;
          }
          double changeLower = (changeLowerRhs) ? changeLowerRhs[iRow] : 0.0;
          double changeUpper = (changeUpperRhs) ? changeUpperRhs[iRow] : 0.0;
          if (lower>-1.0e20&&upper<1.0e20) {
            if (lower+maxTheta*changeLower>upper+maxTheta*changeUpper) {
              maxTheta = (upper-lower)/(changeLower-changeUpper);
            }
          }
          if (lower>-1.0e20) {
            lower_[numberColumns_+iRow] += startingTheta*changeLower;
            chgLower[numberColumns_+iRow]=changeLower;
          }
          if (upper<1.0e20) {
            upper_[numberColumns_+iRow] += startingTheta*changeUpper;
            chgUpper[numberColumns_+iRow]=changeUpper;
          }
        }
      }
      if (maxTheta>0.0) {
        if (changeLowerBound||changeUpperBound) {
          for (iColumn=0;iColumn<numberColumns_;iColumn++) {
            double lower = columnLower_[iColumn];
            double upper = columnUpper_[iColumn];
            if (lower>upper) {
              maxTheta=-1.0;
              break;
            }
            double changeLower = (changeLowerBound) ? changeLowerBound[iColumn] : 0.0;
            double changeUpper = (changeUpperBound) ? changeUpperBound[iColumn] : 0.0;
            if (lower>-1.0e20&&upper<1.0e20) {
              if (lower+maxTheta*changeLower>upper+maxTheta*changeUpper) {
                maxTheta = (upper-lower)/(changeLower-changeUpper);
              }
            }
            if (lower>-1.0e20) {
              lower_[iColumn] += startingTheta*changeLower;
              chgLower[iColumn]=changeLower;
            }
            if (upper<1.0e20) {
              upper_[iColumn] += startingTheta*changeUpper;
              chgUpper[iColumn]=changeUpper;
            }
          }
        }
        if (maxTheta==1.0e50)
          maxTheta = COIN_DBL_MAX;
      }
      if (maxTheta<0.0) {
        // bad ranges or initial
        returnCode = -1;
      }
      endingTheta = CoinMin(endingTheta,maxTheta);
      if (endingTheta<startingTheta) {
        // bad initial
        returnCode = -2;
      }
    }
    double saveEndingTheta=endingTheta;
    if (!returnCode) {
      if (changeObjective) {
        for (iColumn=0;iColumn<numberColumns_;iColumn++) {
          chgObjective[iColumn] = changeObjective[iColumn];
          cost_[iColumn] += startingTheta*changeObjective[iColumn];
        }
      }
      double * saveDuals=NULL;
      ((ClpSimplexDual *) this)->gutsOfDual(0,saveDuals,-1,data);
      assert (!problemStatus_);
      // Now do parametrics
      printf("at starting theta of %g objective value is %g\n",startingTheta,
             objectiveValue());
      while (!returnCode) {
        assert (reportIncrement);
        returnCode = parametricsLoop(startingTheta,endingTheta,reportIncrement,
                                     chgLower,chgUpper,chgObjective,data,
                                     canTryQuick);
        if (!returnCode) {
          //double change = endingTheta-startingTheta;
          startingTheta=endingTheta;
          endingTheta = saveEndingTheta;
          //for (int i=0;i<numberTotal;i++) {
          //lower_[i] += change*chgLower[i];
          //upper_[i] += change*chgUpper[i];
          //cost_[i] += change*chgObjective[i];
          //}
          printf("at theta of %g objective value is %g\n",startingTheta,
                 objectiveValue());
          if (startingTheta>=endingTheta)
            break;
        } else if (returnCode==-1) {
          // trouble - do external solve
          needToDoSomething=true;
        } else {
          abort();
        }
      }
    }
    ((ClpSimplexDual *) this)->finishSolve(0);
    
    delete dualRowPivot_;
    dualRowPivot_ = savePivot;
    // Restore any saved stuff
    restoreData(data);
    if (needToDoSomething) {
      double saveStartingTheta=startingTheta; // known to be feasible
      int cleanedUp=1;
      while (cleanedUp) {
        // tweak
        if (cleanedUp==1) {
          if (!reportIncrement)
            startingTheta = CoinMin(startingTheta+1.0e-5,saveEndingTheta);
          else
            startingTheta = CoinMin(startingTheta+reportIncrement,saveEndingTheta);
        } else {
          // restoring to go slowly
          startingTheta=saveStartingTheta;
        }
        // only works if not scaled
        int i;
        const double * obj1 = objective();
        double * obj2 = copyModel.objective();
        const double * lower1 = columnLower_;
        double * lower2 = copyModel.columnLower();
        const double * upper1 = columnUpper_;
        double * upper2 = copyModel.columnUpper();
        for (i=0;i<numberColumns_;i++) {
          obj2[i] = obj1[i] + startingTheta*chgObjective[i];
          lower2[i] = lower1[i] + startingTheta*chgLower[i];
          upper2[i] = upper1[i] + startingTheta*chgUpper[i];
        }
        lower1 = rowLower_;
        lower2 = copyModel.rowLower();
        upper1 = rowUpper_;
        upper2 = copyModel.rowUpper();
        for (i=0;i<numberRows_;i++) {
          lower2[i] = lower1[i] + startingTheta*chgLower[i+numberColumns_];
          upper2[i] = upper1[i] + startingTheta*chgUpper[i+numberColumns_];
        }
        copyModel.dual();
        if (copyModel.problemStatus()) {
          printf("Can not get to theta of %g\n",startingTheta);
          canTryQuick=false; // do slowly to get exact amount
          // back to last known good
          if (cleanedUp==1)
            cleanedUp=2;
          else
            abort();
        } else {
          // and move stuff back
          int numberTotal = numberRows_+numberColumns_;
          memcpy(status_,copyModel.statusArray(),numberTotal);
          memcpy(columnActivity_,copyModel.primalColumnSolution(),numberColumns_*sizeof(double));
          memcpy(rowActivity_,copyModel.primalRowSolution(),numberRows_*sizeof(double));
          cleanedUp=0;
        }
      }
    }
    delete [] chgLower;
    delete [] chgUpper;
    delete [] chgObjective;
  }
  perturbation_ = savePerturbation;
  return problemStatus_;
}
int 
ClpSimplexOther::parametricsLoop(double startingTheta, double & endingTheta,double reportIncrement,
                                 const double * changeLower, const double * changeUpper,
                                 const double * changeObjective, ClpDataSave & data,
                                 bool canTryQuick)
{
  // stuff is already at starting
  // For this crude version just try and go to end
  double change=0.0;
  if (reportIncrement&&canTryQuick) { 
    endingTheta = CoinMin(endingTheta,startingTheta+reportIncrement);
    change = endingTheta-startingTheta;
  }
  int numberTotal = numberRows_+numberColumns_;
  int i;
  for ( i=0;i<numberTotal;i++) {
    lower_[i] += change*changeLower[i];
    upper_[i] += change*changeUpper[i];
    switch(getStatus(i)) {
      
    case basic:
    case isFree:
    case superBasic:
      break;
    case isFixed:
    case atUpperBound:
      solution_[i]=upper_[i];
      break;
    case atLowerBound:
      solution_[i]=lower_[i];
      break;
    }
    cost_[i] += change*changeObjective[i];
  }
  problemStatus_=-1;
  
  // This says whether to restore things etc
  // startup will have factorized so can skip
  int factorType=0;
  // Start check for cycles
  progress_->startCheck();
  // Say change made on first iteration
  changeMade_=1;
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
  while (problemStatus_<0) {
    int iRow, iColumn;
    // clear
    for (iRow=0;iRow<4;iRow++) {
      rowArray_[iRow]->clear();
    }    
    
    for (iColumn=0;iColumn<2;iColumn++) {
      columnArray_[iColumn]->clear();
    }    
    
    // give matrix (and model costs and bounds a chance to be
    // refreshed (normally null)
    matrix_->refresh(this);
    // may factorize, checks if problem finished
    statusOfProblemInParametrics(factorType,data);
    // Say good factorization
    factorType=1;
    if (data.sparseThreshold_) {
      // use default at present
      factorization_->sparseThreshold(0);
      factorization_->goSparse();
    }
    
    // exit if victory declared
    if (problemStatus_>=0)
      break;
    
    // test for maximum iterations
    if (hitMaximumIterations()) {
      problemStatus_=3;
      break;
    }
    // Check event
    {
      int status = eventHandler_->event(ClpEventHandler::endOfFactorization);
      if (status>=0) {
        problemStatus_=5;
        secondaryStatus_=ClpEventHandler::endOfFactorization;
        break;
      }
    }
    // Do iterations
    if (canTryQuick) {
      double * saveDuals=NULL;
      ((ClpSimplexDual *)this)->whileIterating(saveDuals,0);
    } else {
      whileIterating(startingTheta,  endingTheta, reportIncrement,
                     changeLower, changeUpper,
                     changeObjective);
    }
  }
  if (!problemStatus_) {
    theta_=change+startingTheta;
    eventHandler_->event(ClpEventHandler::theta);
    return 0;
  } else if (problemStatus_==10) {
    return -1;
  } else {
    return problemStatus_;
  }
}
/* Checks if finished.  Updates status */
void 
ClpSimplexOther::statusOfProblemInParametrics(int type, ClpDataSave & saveData)
{
  if (type==2) {
    // trouble - go to recovery
    problemStatus_=10;
    return;
  }
  if (problemStatus_>-3||factorization_->pivots()) {
    // factorize
    // later on we will need to recover from singularities
    // also we could skip if first time
    if (type) {
      // is factorization okay?
      if (internalFactorize(1)) {
        // trouble - go to recovery
        problemStatus_=10;
        return;
      }
    }
    if (problemStatus_!=-4||factorization_->pivots()>10)
      problemStatus_=-3;
  }
  // at this stage status is -3 or -4 if looks infeasible
  // get primal and dual solutions
  gutsOfSolution(NULL,NULL);
  double realDualInfeasibilities=sumDualInfeasibilities_;
  // If bad accuracy treat as singular
  if ((largestPrimalError_>1.0e15||largestDualError_>1.0e15)&&numberIterations_) {
    // trouble - go to recovery
    problemStatus_=10;
    return;
  } else if (largestPrimalError_<1.0e-7&&largestDualError_<1.0e-7) {
    // Can reduce tolerance
    double newTolerance = CoinMax(0.99*factorization_->pivotTolerance(),saveData.pivotTolerance_);
    factorization_->pivotTolerance(newTolerance);
  } 
  // Check if looping
  int loop;
  if (type!=2) 
    loop = progress_->looping();
  else
    loop=-1;
  if (loop>=0) {
    problemStatus_ = loop; //exit if in loop
    if (!problemStatus_) {
      // declaring victory
      numberPrimalInfeasibilities_ = 0;
      sumPrimalInfeasibilities_ = 0.0;
    } else {
      problemStatus_ = 10; // instead - try other algorithm
    }
    return;
  } else if (loop<-1) {
    // something may have changed
    gutsOfSolution(NULL,NULL);
  }
  progressFlag_ = 0; //reset progress flag
  if (handler_->detail(CLP_SIMPLEX_STATUS,messages_)<100) {
    handler_->message(CLP_SIMPLEX_STATUS,messages_)
      <<numberIterations_<<objectiveValue();
    handler_->printing(sumPrimalInfeasibilities_>0.0)
      <<sumPrimalInfeasibilities_<<numberPrimalInfeasibilities_;
    handler_->printing(sumDualInfeasibilities_>0.0)
      <<sumDualInfeasibilities_<<numberDualInfeasibilities_;
    handler_->printing(numberDualInfeasibilitiesWithoutFree_
                       <numberDualInfeasibilities_)
                         <<numberDualInfeasibilitiesWithoutFree_;
    handler_->message()<<CoinMessageEol;
  }
  /* If we are primal feasible and any dual infeasibilities are on
     free variables then it is better to go to primal */
  if (!numberPrimalInfeasibilities_&&!numberDualInfeasibilitiesWithoutFree_&&
      numberDualInfeasibilities_) {
    problemStatus_=10;
    return;
  }
  
  // check optimal
  // give code benefit of doubt
  if (sumOfRelaxedDualInfeasibilities_ == 0.0&&
      sumOfRelaxedPrimalInfeasibilities_ == 0.0) {
    // say optimal (with these bounds etc)
    numberDualInfeasibilities_ = 0;
    sumDualInfeasibilities_ = 0.0;
    numberPrimalInfeasibilities_ = 0;
    sumPrimalInfeasibilities_ = 0.0;
  }
  if (dualFeasible()||problemStatus_==-4) {
    progress_->modifyObjective(objectiveValue_
			       -sumDualInfeasibilities_*dualBound_);
  }
  if (numberPrimalInfeasibilities_) {
    if (problemStatus_==-4||problemStatus_==-5) {
      problemStatus_=1; // infeasible
    }
  } else if (numberDualInfeasibilities_) {
    // clean up
    problemStatus_=10;
  } else {
    problemStatus_=0;
  }
  lastGoodIteration_ = numberIterations_;
  if (problemStatus_<0) {
    sumDualInfeasibilities_=realDualInfeasibilities; // back to say be careful
    if (sumDualInfeasibilities_)
      numberDualInfeasibilities_=1;
  }
  // Allow matrices to be sorted etc
  int fake=-999; // signal sort
  matrix_->correctSequence(fake,fake);
}
/* This has the flow between re-factorizations
   Reasons to come out:
   -1 iterations etc
   -2 inaccuracy 
   -3 slight inaccuracy (and done iterations)
   +0 looks optimal (might be unbounded - but we will investigate)
   +1 looks infeasible
   +3 max iterations 
*/
int 
ClpSimplexOther::whileIterating(double startingTheta, double & endingTheta,double reportIncrement,
                                const double * changeLower, const double * changeUpper,
                                const double * changeObjective)
{
  {
    int i;
    for (i=0;i<4;i++) {
      rowArray_[i]->clear();
    }    
    for (i=0;i<2;i++) {
      columnArray_[i]->clear();
    }    
  }      
  // if can't trust much and long way from optimal then relax
  if (largestPrimalError_>10.0)
    factorization_->relaxAccuracyCheck(CoinMin(1.0e2,largestPrimalError_/10.0));
  else
    factorization_->relaxAccuracyCheck(1.0);
  // status stays at -1 while iterating, >=0 finished, -2 to invert
  // status -3 to go to top without an invert
  int returnCode = -1;
  double saveSumDual = sumDualInfeasibilities_; // so we know to be careful
  double useTheta = startingTheta;
  double * primalChange = new double[numberRows_];
  double * dualChange = new double[numberColumns_];
  int numberTotal = numberColumns_+numberRows_;
  int iSequence;
  // See if bounds
  int type=0;
  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    if (changeLower[iSequence]||changeUpper[iSequence]) {
      type=1;
      break;
    }
  }
  // See if objective
  for (iSequence=0;iSequence<numberTotal;iSequence++) {
    if (changeObjective[iSequence]) {
      type |= 2;
      break;
    }
  }
  assert (type);
  while (problemStatus_==-1) {
    double increaseTheta = CoinMin(endingTheta-useTheta,1.0e50);
    
    // Get theta for bounds - we know can't crossover
    int pivotType = nextTheta(type,increaseTheta,primalChange,dualChange,
                              changeLower,changeUpper,changeObjective);
    if (pivotType)
      abort();
    // choose row to go out
    // dualRow will go to virtual row pivot choice algorithm
    ((ClpSimplexDual *) this)->dualRow(-1);
    if (pivotRow_>=0) {
      // we found a pivot row
      if (handler_->detail(CLP_SIMPLEX_PIVOTROW,messages_)<100) {
        handler_->message(CLP_SIMPLEX_PIVOTROW,messages_)
          <<pivotRow_
          <<CoinMessageEol;
      }
      // check accuracy of weights
      dualRowPivot_->checkAccuracy();
      // Get good size for pivot
      // Allow first few iterations to take tiny
      double acceptablePivot=1.0e-9;
      if (numberIterations_>100)
        acceptablePivot=1.0e-8;
      if (factorization_->pivots()>10||
	  (factorization_->pivots()&&saveSumDual))
	acceptablePivot=1.0e-5; // if we have iterated be more strict
      else if (factorization_->pivots()>5)
	acceptablePivot=1.0e-6; // if we have iterated be slightly more strict
      else if (factorization_->pivots())
        acceptablePivot=1.0e-8; // relax
      double bestPossiblePivot=1.0;
      // get sign for finding row of tableau
      // normal iteration
      // create as packed
      double direction=directionOut_;
      rowArray_[0]->createPacked(1,&pivotRow_,&direction);
      factorization_->updateColumnTranspose(rowArray_[1],rowArray_[0]);
      // put row of tableau in rowArray[0] and columnArray[0]
      matrix_->transposeTimes(this,-1.0,
			      rowArray_[0],rowArray_[3],columnArray_[0]);
      // do ratio test for normal iteration
      bestPossiblePivot = ((ClpSimplexDual *) this)->dualColumn(rowArray_[0],
                                                                columnArray_[0],columnArray_[1],
                                                                rowArray_[3],acceptablePivot,NULL);
      if (sequenceIn_>=0) {
	// normal iteration
	// update the incoming column
	double btranAlpha = -alpha_*directionOut_; // for check
	unpackPacked(rowArray_[1]);
	factorization_->updateColumnFT(rowArray_[2],rowArray_[1]);
	// and update dual weights (can do in parallel - with extra array)
	alpha_ = dualRowPivot_->updateWeights(rowArray_[0],
					      rowArray_[2],
					      rowArray_[1]);
	// see if update stable
#ifdef CLP_DEBUG
	if ((handler_->logLevel()&32))
	  printf("btran alpha %g, ftran alpha %g\n",btranAlpha,alpha_);
#endif
	double checkValue=1.0e-7;
	// if can't trust much and long way from optimal then relax
	if (largestPrimalError_>10.0)
	  checkValue = CoinMin(1.0e-4,1.0e-8*largestPrimalError_);
	if (fabs(btranAlpha)<1.0e-12||fabs(alpha_)<1.0e-12||
	    fabs(btranAlpha-alpha_)>checkValue*(1.0+fabs(alpha_))) {
	  handler_->message(CLP_DUAL_CHECK,messages_)
	    <<btranAlpha
	    <<alpha_
	    <<CoinMessageEol;
	  if (factorization_->pivots()) {
	    dualRowPivot_->unrollWeights();
	    problemStatus_=-2; // factorize now
	    rowArray_[0]->clear();
	    rowArray_[1]->clear();
	    columnArray_[0]->clear();
	    returnCode=-2;
	    break;
	  } else {
	    // take on more relaxed criterion
            double test;
	    if (fabs(btranAlpha)<1.0e-8||fabs(alpha_)<1.0e-8)
              test = 1.0e-1*fabs(alpha_);
            else
              test = 1.0e-4*(1.0+fabs(alpha_));
	    if (fabs(btranAlpha)<1.0e-12||fabs(alpha_)<1.0e-12||
		fabs(btranAlpha-alpha_)>test) {
	      dualRowPivot_->unrollWeights();
	      // need to reject something
	      char x = isColumn(sequenceOut_) ? 'C' :'R';
	      handler_->message(CLP_SIMPLEX_FLAG,messages_)
		<<x<<sequenceWithin(sequenceOut_)
		<<CoinMessageEol;
	      setFlagged(sequenceOut_);
	      progress_->clearBadTimes();
	      lastBadIteration_ = numberIterations_; // say be more cautious
	      rowArray_[0]->clear();
	      rowArray_[1]->clear();
	      columnArray_[0]->clear();
              if (fabs(alpha_)<1.0e-10&&fabs(btranAlpha)<1.0e-8&&numberIterations_>100) {
                //printf("I think should declare infeasible\n");
                problemStatus_=1;
                returnCode=1;
                break;
              }
	      continue;
	    }
	  }
	}
	// update duals BEFORE replaceColumn so can do updateColumn
	double objectiveChange=0.0;
	// do duals first as variables may flip bounds
	// rowArray_[0] and columnArray_[0] may have flips
	// so use rowArray_[3] for work array from here on
	int nswapped = 0;
	//rowArray_[0]->cleanAndPackSafe(1.0e-60);
	//columnArray_[0]->cleanAndPackSafe(1.0e-60);
        nswapped = ((ClpSimplexDual *) this)->updateDualsInDual(rowArray_[0],columnArray_[0],
                                     rowArray_[2],theta_,
                                     objectiveChange,false);

	// which will change basic solution
	if (nswapped) {
	  factorization_->updateColumn(rowArray_[3],rowArray_[2]);
	  dualRowPivot_->updatePrimalSolution(rowArray_[2],
					      1.0,objectiveChange);
	  // recompute dualOut_
	  valueOut_ = solution_[sequenceOut_];
	  if (directionOut_<0) {
	    dualOut_ = valueOut_ - upperOut_;
	  } else {
	    dualOut_ = lowerOut_ - valueOut_;
	  }
	}
	// amount primal will move
	double movement = -dualOut_*directionOut_/alpha_;
	// so objective should increase by fabs(dj)*movement
	// but we already have objective change - so check will be good
	if (objectiveChange+fabs(movement*dualIn_)<-1.0e-5) {
#ifdef CLP_DEBUG
	  if (handler_->logLevel()&32)
	    printf("movement %g, swap change %g, rest %g  * %g\n",
		   objectiveChange+fabs(movement*dualIn_),
		   objectiveChange,movement,dualIn_);
#endif
	  if(factorization_->pivots()) {
	    // going backwards - factorize
	    dualRowPivot_->unrollWeights();
	    problemStatus_=-2; // factorize now
	    returnCode=-2;
	    break;
	  }
	}
	CoinAssert(fabs(dualOut_)<1.0e50);
	// if stable replace in basis
	int updateStatus = factorization_->replaceColumn(this,
							 rowArray_[2],
							 rowArray_[1],
							 pivotRow_,
							 alpha_);
	// if no pivots, bad update but reasonable alpha - take and invert
	if (updateStatus==2&&
		   !factorization_->pivots()&&fabs(alpha_)>1.0e-5)
	  updateStatus=4;
	if (updateStatus==1||updateStatus==4) {
	  // slight error
	  if (factorization_->pivots()>5||updateStatus==4) {
	    problemStatus_=-2; // factorize now
	    returnCode=-3;
	  }
	} else if (updateStatus==2) {
	  // major error
	  dualRowPivot_->unrollWeights();
	  // later we may need to unwind more e.g. fake bounds
	  if (factorization_->pivots()) {
	    problemStatus_=-2; // factorize now
	    returnCode=-2;
	    break;
	  } else {
	    // need to reject something
	    char x = isColumn(sequenceOut_) ? 'C' :'R';
	    handler_->message(CLP_SIMPLEX_FLAG,messages_)
	      <<x<<sequenceWithin(sequenceOut_)
	      <<CoinMessageEol;
	    setFlagged(sequenceOut_);
	    progress_->clearBadTimes();
	    lastBadIteration_ = numberIterations_; // say be more cautious
	    rowArray_[0]->clear();
	    rowArray_[1]->clear();
	    columnArray_[0]->clear();
	    // make sure dual feasible
	    // look at all rows and columns
	    double objectiveChange=0.0;
	    ((ClpSimplexDual *) this)->updateDualsInDual(rowArray_[0],columnArray_[0],rowArray_[1],
			      0.0,objectiveChange,true);
	    continue;
	  }
	} else if (updateStatus==3) {
	  // out of memory
	  // increase space if not many iterations
	  if (factorization_->pivots()<
	      0.5*factorization_->maximumPivots()&&
	      factorization_->pivots()<200)
	    factorization_->areaFactor(
				       factorization_->areaFactor() * 1.1);
	  problemStatus_=-2; // factorize now
	} else if (updateStatus==5) {
	  problemStatus_=-2; // factorize now
	}
	// update primal solution
	if (theta_<0.0) {
#ifdef CLP_DEBUG
	  if (handler_->logLevel()&32)
	    printf("negative theta %g\n",theta_);
#endif
	  theta_=0.0;
	}
	// do actual flips
	((ClpSimplexDual *) this)->flipBounds(rowArray_[0],columnArray_[0],theta_);
	//rowArray_[1]->expand();
	dualRowPivot_->updatePrimalSolution(rowArray_[1],
					    movement,
					    objectiveChange);
	// modify dualout
	dualOut_ /= alpha_;
	dualOut_ *= -directionOut_;
	//setStatus(sequenceIn_,basic);
	dj_[sequenceIn_]=0.0;
	double oldValue=valueIn_;
	if (directionIn_==-1) {
	  // as if from upper bound
	  valueIn_ = upperIn_+dualOut_;
	} else {
	  // as if from lower bound
	  valueIn_ = lowerIn_+dualOut_;
	}
	objectiveChange += cost_[sequenceIn_]*(valueIn_-oldValue);
	// outgoing
	// set dj to zero unless values pass
	if (directionOut_>0) {
	  valueOut_ = lowerOut_;
          dj_[sequenceOut_] = theta_;
	} else {
	  valueOut_ = upperOut_;
          dj_[sequenceOut_] = -theta_;
	}
	solution_[sequenceOut_]=valueOut_;
	int whatNext=housekeeping(objectiveChange);
	// and set bounds correctly
	((ClpSimplexDual *) this)->originalBound(sequenceIn_); 
	((ClpSimplexDual *) this)->changeBound(sequenceOut_);
	if (whatNext==1) {
	  problemStatus_ =-2; // refactorize
	} else if (whatNext==2) {
	  // maximum iterations or equivalent
	  problemStatus_= 3;
	  returnCode=3;
	  break;
	}
	// Check event
	{
	  int status = eventHandler_->event(ClpEventHandler::endOfIteration);
	  if (status>=0) {
	    problemStatus_=5;
	    secondaryStatus_=ClpEventHandler::endOfIteration;
	    returnCode=4;
	    break;
	  }
	}
      } else {
	// no incoming column is valid
        pivotRow_=-1;
#ifdef CLP_DEBUG
	if (handler_->logLevel()&32)
	  printf("** no column pivot\n");
#endif
	if (factorization_->pivots()<5) {
          // If not in branch and bound etc save ray
          if ((specialOptions_&(1024|4096))==0) {
	    // create ray anyway
	    delete [] ray_;
	    ray_ = new double [ numberRows_];
	    rowArray_[0]->expand(); // in case packed
	    ClpDisjointCopyN(rowArray_[0]->denseVector(),numberRows_,ray_);
          }
	  // If we have just factorized and infeasibility reasonable say infeas
	  if (((specialOptions_&4096)!=0||bestPossiblePivot<1.0e-11)&&dualBound_>1.0e8) {
	    if (valueOut_>upperOut_+1.0e-3||valueOut_<lowerOut_-1.0e-3
		|| (specialOptions_&64)==0) {
	      // say infeasible
	      problemStatus_=1;
              // unless primal feasible!!!!
              //printf("%d %g %d %g\n",numberPrimalInfeasibilities_,sumPrimalInfeasibilities_,
              //     numberDualInfeasibilities_,sumDualInfeasibilities_);
              if (numberDualInfeasibilities_)
                problemStatus_=10;
	      rowArray_[0]->clear();
	      columnArray_[0]->clear();
	      returnCode=1;
	      break;
	    }
	  }
	  // If special option set - put off as long as possible
	  if ((specialOptions_&64)==0) {
	    problemStatus_=-4; //say looks infeasible
	  } else {
	    // flag
	    char x = isColumn(sequenceOut_) ? 'C' :'R';
	    handler_->message(CLP_SIMPLEX_FLAG,messages_)
	      <<x<<sequenceWithin(sequenceOut_)
	      <<CoinMessageEol;
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
	returnCode=1;
	break;
      }
    } else {
      // no pivot row
#ifdef CLP_DEBUG
      if (handler_->logLevel()&32)
	printf("** no row pivot\n");
#endif
      int numberPivots = factorization_->pivots();
      bool specialCase;
      int useNumberFake;
      returnCode=0;
      if (numberPivots<20&&
	  (specialOptions_&2048)!=0&&!numberChanged_&&perturbation_>=100
	  &&dualBound_>1.0e8) {
	specialCase=true;
	// as dual bound high - should be okay
	useNumberFake=0;
      } else {
	specialCase=false;
	useNumberFake=numberFake_;
      }
      if (!numberPivots||specialCase) {
	// may have crept through - so may be optimal
	// check any flagged variables
	int iRow;
	for (iRow=0;iRow<numberRows_;iRow++) {
	  int iPivot=pivotVariable_[iRow];
	  if (flagged(iPivot))
	    break;
	}
	if (iRow<numberRows_&&numberPivots) {
	  // try factorization
	  returnCode=-2;
	}
        
	if (useNumberFake||numberDualInfeasibilities_) {
	  // may be dual infeasible
	  problemStatus_=-5;
	} else {
	  if (iRow<numberRows_) {
	    problemStatus_=-5;
	  } else {
            if (numberPivots) {
              // objective may be wrong
              objectiveValue_ = innerProduct(cost_,
                                                                        numberColumns_+numberRows_,
                                                                        solution_);
              objectiveValue_ += objective_->nonlinearOffset();
              objectiveValue_ /= (objectiveScale_*rhsScale_);
              if ((specialOptions_&16384)==0) {
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
                int number=0;
                for (iRow=0;iRow<numberRows_;iRow++) {
                  int iPivot=pivotVariable_[iRow];
                  double value = cost_[iPivot];
                  dual_[iRow]=value;
                  if (value) {
                    index[number++]=iRow;
                  }
                }
                arrayVector->setNumElements(number);
                // Extended duals before "updateTranspose"
                matrix_->dualExpanded(this,arrayVector,NULL,0);
                // Btran basic costs
                rowArray_[2]->clear();
                factorization_->updateColumnTranspose(rowArray_[2],arrayVector);
                // and return vector
                arrayVector->setDenseVector(array);
              }
            }
	    problemStatus_=0;
	    sumPrimalInfeasibilities_=0.0;
            if ((specialOptions_&(1024+16384))!=0) {
              CoinIndexedVector * arrayVector = rowArray_[1];
              arrayVector->clear();
              double * rhs = arrayVector->denseVector();
              times(1.0,solution_,rhs);
              bool bad2=false;
              int i;
              for ( i=0;i<numberRows_;i++) {
                if (rhs[i]<rowLowerWork_[i]-primalTolerance_||
                    rhs[i]>rowUpperWork_[i]+primalTolerance_) {
                  bad2=true;
                } else if (fabs(rhs[i]-rowActivityWork_[i])>1.0e-3) {
                }
                rhs[i]=0.0;
              }
              for ( i=0;i<numberColumns_;i++) {
                if (solution_[i]<columnLowerWork_[i]-primalTolerance_||
                    solution_[i]>columnUpperWork_[i]+primalTolerance_) {
                  bad2=true;
                }
              }
              if (bad2) {
                problemStatus_=-3;
                returnCode=-2;
                // Force to re-factorize early next time
                int numberPivots = factorization_->pivots();
                forceFactorization_=CoinMin(forceFactorization_,(numberPivots+1)>>1);
              }
            }
	  }
	}
      } else {
	problemStatus_=-3;
        returnCode=-2;
	// Force to re-factorize early next time
	int numberPivots = factorization_->pivots();
	forceFactorization_=CoinMin(forceFactorization_,(numberPivots+1)>>1);
      }
      break;
    }
  }
  delete [] primalChange;
  delete [] dualChange;
  return returnCode;
}
// Computes next theta and says if objective or bounds (0= bounds, 1 objective, -1 none)
int 
ClpSimplexOther::nextTheta(int type, double maxTheta, double * primalChange, double * dualChange,
                           const double * changeLower, const double * changeUpper,
                           const double * changeObjective)
{
  int numberTotal = numberColumns_+numberRows_;
  int iSequence;
  int iRow;
  theta_=maxTheta;
  bool toLower=false;
  if ((type&1)!=0) {
    // get change 
    for (iSequence=0;iSequence<numberTotal;iSequence++) {
      primalChange[iSequence]=0.0;
      switch(getStatus(iSequence)) {
        
      case basic:
      case isFree:
      case superBasic:
        break;
      case isFixed:
      case atUpperBound:
        primalChange[iSequence]=changeUpper[iSequence];
        break;
      case atLowerBound:
        primalChange[iSequence]=changeLower[iSequence];
        break;
      }
    }
    // use array
    double * array = rowArray_[1]->denseVector();
    times(1.0,primalChange,array);
    int * index = rowArray_[1]->getIndices();
    int number=0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      double value = array[iRow];
      if (value) {
	array[iRow]=value;
	index[number++]=iRow;
      }
    }
    // ftran it
    rowArray_[1]->setNumElements(number);
    factorization_->updateColumn(rowArray_[0],rowArray_[1]);
    number=rowArray_[1]->getNumElements();
    pivotRow_=-1;
    for (iRow=0;iRow<number;iRow++) {
      int iPivot = index[iRow];
      iSequence = pivotVariable_[iPivot];
      // solution value will be sol - theta*alpha
      // bounds will be bounds + change *theta
      double currentSolution = solution_[iSequence];
      double currentLower = lower_[iSequence];
      double currentUpper = upper_[iSequence];
      double alpha = array[iPivot];
      assert (currentSolution>=currentLower-primalTolerance_);
      assert (currentSolution<=currentUpper+primalTolerance_);
      double thetaCoefficient;
      double hitsLower = COIN_DBL_MAX;
      thetaCoefficient = changeLower[iSequence]+alpha;
      if (fabs(thetaCoefficient)>1.0e-8)
        hitsLower = (currentSolution-currentLower)/thetaCoefficient;
      if (hitsLower<0.0) {
        // does not hit - but should we check further
        hitsLower=COIN_DBL_MAX;
      }
      double hitsUpper = COIN_DBL_MAX;
      thetaCoefficient = changeUpper[iSequence]+alpha;
      if (fabs(thetaCoefficient)>1.0e-8)
        hitsUpper = (currentSolution-currentUpper)/thetaCoefficient;
      if (hitsUpper<0.0) {
        // does not hit - but should we check further
        hitsUpper=COIN_DBL_MAX;
      }
      if (CoinMin(hitsLower,hitsUpper)<theta_) {
        theta_ = CoinMin(hitsLower,hitsUpper);
        toLower = hitsLower<hitsUpper;
        pivotRow_=iPivot;
      }
    }
  }
  if ((type&2)!=0) {
    abort();
  }
  if (pivotRow_>=0) {
    sequenceOut_ = pivotVariable_[pivotRow_];
    valueOut_ = solution_[sequenceOut_];
    lowerOut_ = lower_[sequenceOut_];
    upperOut_ = upper_[sequenceOut_];
    if (!toLower) {
      directionOut_ = -1;
      dualOut_ = valueOut_ - upperOut_;
    } else if (valueOut_<lowerOut_) {
      directionOut_ = 1;
      dualOut_ = lowerOut_ - valueOut_;
    }
    return 0;
  } else {
    return -1;
  }
}
