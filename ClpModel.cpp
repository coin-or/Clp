// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cmath>
#include <cassert>
#include <cfloat>
#include <string>
#include <cstdio>
#include <iostream>


#include "CoinPragma.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "ClpModel.hpp"
#include "ClpEventHandler.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinMpsIO.hpp"
#include "ClpMessage.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpQuadraticObjective.hpp"

//#############################################################################

ClpModel::ClpModel () :

  optimizationDirection_(1),
  objectiveValue_(0.0),
  smallElement_(1.0e-20),
  objectiveScale_(1.0),
  rhsScale_(1.0),
  numberRows_(0),
  numberColumns_(0),
  rowActivity_(NULL),
  columnActivity_(NULL),
  dual_(NULL),
  reducedCost_(NULL),
  rowLower_(NULL),
  rowUpper_(NULL),
  objective_(NULL),
  rowObjective_(NULL),
  columnLower_(NULL),
  columnUpper_(NULL),
  matrix_(NULL),
  rowCopy_(NULL),
  ray_(NULL),
  rowScale_(NULL),
  columnScale_(NULL),
  scalingFlag_(3),
  status_(NULL),
  integerType_(NULL),
  userPointer_(NULL),
  numberIterations_(0),
  solveType_(0),
  problemStatus_(-1),
  secondaryStatus_(0),
  lengthNames_(0),
  defaultHandler_(true),
  rowNames_(),
  columnNames_()
{
  intParam_[ClpMaxNumIteration] = 99999999;
  intParam_[ClpMaxNumIterationHotStart] = 9999999;

  dblParam_[ClpDualObjectiveLimit] = COIN_DBL_MAX;
  dblParam_[ClpPrimalObjectiveLimit] = COIN_DBL_MAX;
  dblParam_[ClpDualTolerance] = 1e-7;
  dblParam_[ClpPrimalTolerance] = 1e-7;
  dblParam_[ClpObjOffset] = 0.0;
  dblParam_[ClpMaxSeconds] = -1.0;

  strParam_[ClpProbName] = "ClpDefaultName";
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(1);
  eventHandler_ = new ClpEventHandler();
  messages_ = ClpMessage();
  CoinSeedRandom(1234567);
}

//-----------------------------------------------------------------------------

ClpModel::~ClpModel ()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  gutsOfDelete();
}
void ClpModel::gutsOfDelete()
{
  delete [] rowActivity_;
  rowActivity_=NULL;
  delete [] columnActivity_;
  columnActivity_=NULL;
  delete [] dual_;
  dual_=NULL;
  delete [] reducedCost_;
  reducedCost_=NULL;
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] rowObjective_;
  rowLower_=NULL;
  rowUpper_=NULL;
  rowObjective_=NULL;
  delete [] columnLower_;
  delete [] columnUpper_;
  delete objective_;
  columnLower_=NULL;
  columnUpper_=NULL;
  objective_=NULL;
  delete matrix_;
  matrix_=NULL;
  delete rowCopy_;
  rowCopy_=NULL;
  delete [] ray_;
  ray_ = NULL;
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
  delete [] integerType_;
  integerType_ = NULL;
  delete [] status_;
  status_=NULL;
  delete eventHandler_;
  eventHandler_=NULL;
}
//#############################################################################
void ClpModel::setPrimalTolerance( double value) 
{
  if (value>0.0&&value<1.0e10)
    dblParam_[ClpPrimalTolerance]=value;
}
void ClpModel::setDualTolerance( double value) 
{
  if (value>0.0&&value<1.0e10)
    dblParam_[ClpDualTolerance]=value;
}
void ClpModel::setOptimizationDirection( double value) 
{
  optimizationDirection_=value;
}
void
ClpModel::gutsOfLoadModel (int numberRows, int numberColumns, 
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
				const double * rowObjective)
{
  // save event handler in case already set
  ClpEventHandler * handler = eventHandler_->clone();
  gutsOfDelete();
  eventHandler_ = handler;
  numberRows_=numberRows;
  numberColumns_=numberColumns;
  rowActivity_=new double[numberRows_];
  columnActivity_=new double[numberColumns_];
  dual_=new double[numberRows_];
  reducedCost_=new double[numberColumns_];

  ClpFillN(dual_,numberRows_,0.0);
  ClpFillN(reducedCost_,numberColumns_,0.0);
  int iRow,iColumn;

  rowLower_=ClpCopyOfArray(rowlb,numberRows_,-COIN_DBL_MAX);
  rowUpper_=ClpCopyOfArray(rowub,numberRows_,COIN_DBL_MAX);
  double * objective=ClpCopyOfArray(obj,numberColumns_,0.0);
  objective_ = new ClpLinearObjective(objective,numberColumns_);
  delete [] objective;
  rowObjective_=ClpCopyOfArray(rowObjective,numberRows_);
  columnLower_=ClpCopyOfArray(collb,numberColumns_,0.0);
  columnUpper_=ClpCopyOfArray(colub,numberColumns_,COIN_DBL_MAX);
  // set default solution
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (rowLower_[iRow]>0.0) {
      rowActivity_[iRow]=rowLower_[iRow];
    } else if (rowUpper_[iRow]<0.0) {
      rowActivity_[iRow]=rowUpper_[iRow];
    } else {
      rowActivity_[iRow]=0.0;
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (columnLower_[iColumn]>0.0) {
      columnActivity_[iColumn]=columnLower_[iColumn];
    } else if (columnUpper_[iColumn]<0.0) {
      columnActivity_[iColumn]=columnUpper_[iColumn];
    } else {
      columnActivity_[iColumn]=0.0;
    }
  }
}
// This just loads up a row objective
void ClpModel::setRowObjective(const double * rowObjective)
{
  delete [] rowObjective_;
  rowObjective_=ClpCopyOfArray(rowObjective,numberRows_);
}
void
ClpModel::loadProblem (  const ClpMatrixBase& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
				const double * rowObjective)
{
  gutsOfLoadModel(matrix.getNumRows(),matrix.getNumCols(),
		  collb, colub, obj, rowlb, rowub, rowObjective);
  if (matrix.isColOrdered()) {
    matrix_=matrix.clone();
  } else {
    // later may want to keep as unknown class
    CoinPackedMatrix matrix2;
    matrix2.reverseOrderedCopyOf(*matrix.getPackedMatrix());
    matrix.releasePackedMatrix();
    matrix_=new ClpPackedMatrix(matrix2);
  }    
}
void
ClpModel::loadProblem (  const CoinPackedMatrix& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
				const double * rowObjective)
{
  gutsOfLoadModel(matrix.getNumRows(),matrix.getNumCols(),
		  collb, colub, obj, rowlb, rowub, rowObjective);
  if (matrix.isColOrdered()) {
    matrix_=new ClpPackedMatrix(matrix);
  } else {
    CoinPackedMatrix matrix2;
    matrix2.reverseOrderedCopyOf(matrix);
    matrix_=new ClpPackedMatrix(matrix2);
  }    
}
void
ClpModel::loadProblem ( 
			      const int numcols, const int numrows,
			      const CoinBigIndex* start, const int* index,
			      const double* value,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub,
			      const double * rowObjective)
{
  gutsOfLoadModel(numrows, numcols,
		  collb, colub, obj, rowlb, rowub, rowObjective);
  CoinPackedMatrix matrix(true,numrows,numcols,start[numcols],
			      value,index,start,NULL);
  matrix_ = new ClpPackedMatrix(matrix);
}
void
ClpModel::loadProblem ( 
			      const int numcols, const int numrows,
			      const CoinBigIndex* start, const int* index,
			      const double* value,const int* length,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub,
			      const double * rowObjective)
{
  gutsOfLoadModel(numrows, numcols,
		  collb, colub, obj, rowlb, rowub, rowObjective);
  // Compute number of elements
  int numberElements = 0;
  int i;
  for (i=0;i<numcols;i++) 
    numberElements += length[i];
  CoinPackedMatrix matrix(true,numrows,numcols,numberElements,
			      value,index,start,length);
  matrix_ = new ClpPackedMatrix(matrix);
}
void
ClpModel::getRowBound(int iRow, double& lower, double& upper) const
{
  lower=-COIN_DBL_MAX;
  upper=COIN_DBL_MAX;
  if (rowUpper_)
    upper=rowUpper_[iRow];
  if (rowLower_)
    lower=rowLower_[iRow];
}
//#############################################################################
// Copy constructor. 
ClpModel::ClpModel(const ClpModel &rhs, int scalingMode) :
  optimizationDirection_(rhs.optimizationDirection_),
  numberRows_(rhs.numberRows_),
  numberColumns_(rhs.numberColumns_)
{
  gutsOfCopy(rhs);
  if (scalingMode>=0&&matrix_&&
      matrix_->allElementsInRange(this,smallElement_,1.0e20)) {
    // really do scaling
    scalingFlag_=scalingMode;
    delete [] rowScale_;
    rowScale_ = NULL;
    delete [] columnScale_;
    columnScale_ = NULL;
    if (!matrix_->scale(this)) {
      // scaling worked - now apply
      gutsOfScaling();
      // pretend not scaled
      scalingFlag_ = -scalingFlag_;
    } else {
      // not scaled
      scalingFlag_=0;
    }
  }
  CoinSeedRandom(1234567);
}
// Assignment operator. This copies the data
ClpModel & 
ClpModel::operator=(const ClpModel & rhs)
{
  if (this != &rhs) {
    if (defaultHandler_) {
      delete handler_;
      handler_ = NULL;
    }
    gutsOfDelete();
    optimizationDirection_ = rhs.optimizationDirection_;
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    gutsOfCopy(rhs);
  }
  return *this;
}
// Does most of copying
void 
ClpModel::gutsOfCopy(const ClpModel & rhs, bool trueCopy)
{
  defaultHandler_ = rhs.defaultHandler_;
  if (defaultHandler_) 
    handler_ = new CoinMessageHandler(*rhs.handler_);
   else 
    handler_ = rhs.handler_;
  eventHandler_ = rhs.eventHandler_->clone();
  messages_ = rhs.messages_;
  intParam_[ClpMaxNumIteration] = rhs.intParam_[ClpMaxNumIteration];
  intParam_[ClpMaxNumIterationHotStart] = 
    rhs.intParam_[ClpMaxNumIterationHotStart];

  dblParam_[ClpDualObjectiveLimit] = rhs.dblParam_[ClpDualObjectiveLimit];
  dblParam_[ClpPrimalObjectiveLimit] = rhs.dblParam_[ClpPrimalObjectiveLimit];
  dblParam_[ClpDualTolerance] = rhs.dblParam_[ClpDualTolerance];
  dblParam_[ClpPrimalTolerance] = rhs.dblParam_[ClpPrimalTolerance];
  dblParam_[ClpObjOffset] = rhs.dblParam_[ClpObjOffset];
  dblParam_[ClpMaxSeconds] = rhs.dblParam_[ClpMaxSeconds];

  strParam_[ClpProbName] = rhs.strParam_[ClpProbName];

  optimizationDirection_ = rhs.optimizationDirection_;
  objectiveValue_=rhs.objectiveValue_;
  smallElement_ = rhs.smallElement_;
  objectiveScale_ = rhs.objectiveScale_;
  rhsScale_ = rhs.rhsScale_;
  numberIterations_ = rhs.numberIterations_;
  solveType_ = rhs.solveType_;
  problemStatus_ = rhs.problemStatus_;
  secondaryStatus_ = rhs.secondaryStatus_;
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
  userPointer_ = rhs.userPointer_;
  scalingFlag_ = rhs.scalingFlag_;
  if (trueCopy) {
    lengthNames_ = rhs.lengthNames_;
    rowNames_ = rhs.rowNames_;
    columnNames_ = rhs.columnNames_;
    if (rhs.integerType_) {
      integerType_ = new char[numberColumns_];
      memcpy(integerType_,rhs.integerType_,numberColumns_*sizeof(char));
    } else {
      integerType_ = NULL;
    }
    if (rhs.rowActivity_) {
      rowActivity_=new double[numberRows_];
      columnActivity_=new double[numberColumns_];
      dual_=new double[numberRows_];
      reducedCost_=new double[numberColumns_];
      ClpDisjointCopyN ( rhs.rowActivity_, numberRows_ ,
			  rowActivity_);
      ClpDisjointCopyN ( rhs.columnActivity_, numberColumns_ ,
			  columnActivity_);
      ClpDisjointCopyN ( rhs.dual_, numberRows_ , 
			  dual_);
      ClpDisjointCopyN ( rhs.reducedCost_, numberColumns_ ,
			  reducedCost_);
    } else {
      rowActivity_=NULL;
      columnActivity_=NULL;
      dual_=NULL;
      reducedCost_=NULL;
    }
    rowLower_ = ClpCopyOfArray ( rhs.rowLower_, numberRows_ );
    rowUpper_ = ClpCopyOfArray ( rhs.rowUpper_, numberRows_ );
    columnLower_ = ClpCopyOfArray ( rhs.columnLower_, numberColumns_ );
    columnUpper_ = ClpCopyOfArray ( rhs.columnUpper_, numberColumns_ );
    rowScale_ = ClpCopyOfArray(rhs.rowScale_,numberRows_);
    columnScale_ = ClpCopyOfArray(rhs.columnScale_,numberColumns_);
    if (rhs.objective_)
      objective_  = rhs.objective_->clone();
    else
      objective_ = NULL;
    rowObjective_ = ClpCopyOfArray ( rhs.rowObjective_, numberRows_ );
    status_ = ClpCopyOfArray( rhs.status_,numberColumns_+numberRows_);
    ray_ = NULL;
    if (problemStatus_==1&&!secondaryStatus_)
      ray_ = ClpCopyOfArray (rhs.ray_,numberRows_);
    else if (problemStatus_==2)
      ray_ = ClpCopyOfArray (rhs.ray_,numberColumns_);
    if (rhs.rowCopy_) {
      rowCopy_ = rhs.rowCopy_->clone();
    } else {
      rowCopy_=NULL;
    }
    matrix_=NULL;
    if (rhs.matrix_) {
      matrix_ = rhs.matrix_->clone();
    }
  } else {
    rowActivity_ = rhs.rowActivity_;
    columnActivity_ = rhs.columnActivity_;
    dual_ = rhs.dual_;
    reducedCost_ = rhs.reducedCost_;
    rowLower_ = rhs.rowLower_;
    rowUpper_ = rhs.rowUpper_;
    objective_ = rhs.objective_;
    rowObjective_ = rhs.rowObjective_;
    columnLower_ = rhs.columnLower_;
    columnUpper_ = rhs.columnUpper_;
    matrix_ = rhs.matrix_;
    rowCopy_ = NULL;
    ray_ = rhs.ray_;
    //rowScale_ = rhs.rowScale_;
    //columnScale_ = rhs.columnScale_;
    lengthNames_ = 0;
    rowNames_ = std::vector<std::string> ();
    columnNames_ = std::vector<std::string> ();
    integerType_ = NULL;
    status_ = rhs.status_;
  }
}
/* Borrow model.  This is so we dont have to copy large amounts
   of data around.  It assumes a derived class wants to overwrite
   an empty model with a real one - while it does an algorithm */
void 
ClpModel::borrowModel(ClpModel & rhs)
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  gutsOfDelete();
  optimizationDirection_ = rhs.optimizationDirection_;
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
  delete [] rhs.ray_;
  rhs.ray_=NULL;
  gutsOfCopy(rhs,false);
}
// Return model - nulls all arrays so can be deleted safely
void 
ClpModel::returnModel(ClpModel & otherModel)
{
  otherModel.objectiveValue_=objectiveValue_;
  otherModel.numberIterations_ = numberIterations_;
  otherModel.problemStatus_ = problemStatus_;
  otherModel.secondaryStatus_ = secondaryStatus_;
  rowActivity_ = NULL;
  columnActivity_ = NULL;
  dual_ = NULL;
  reducedCost_ = NULL;
  rowLower_ = NULL;
  rowUpper_ = NULL;
  objective_ = NULL;
  rowObjective_ = NULL;
  columnLower_ = NULL;
  columnUpper_ = NULL;
  matrix_ = NULL;
  rowCopy_ = NULL;
  delete [] otherModel.ray_;
  otherModel.ray_ = ray_;
  ray_ = NULL;
  //rowScale_=NULL;
  //columnScale_=NULL;
  // do status
  if (otherModel.status_!=status_) {
    delete [] otherModel.status_;
    otherModel.status_ = status_;
  }
  status_ = NULL;
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
}
//#############################################################################
// Parameter related methods
//#############################################################################

bool
ClpModel::setIntParam(ClpIntParam key, int value)
{
  switch (key) {
  case ClpMaxNumIteration:
    if (value < 0)
      return false;
    break;
  case ClpMaxNumIterationHotStart:
    if (value < 0)
      return false;
    break;
  case ClpLastIntParam:
    return false;
  }
  intParam_[key] = value;
  return true;
}

//-----------------------------------------------------------------------------

bool
ClpModel::setDblParam(ClpDblParam key, double value)
{

  switch (key) {
  case ClpDualObjectiveLimit:
    break;

  case ClpPrimalObjectiveLimit:
    break;

  case ClpDualTolerance: 
    if (value<=0.0||value>1.0e10)
      return false;
    break;
    
  case ClpPrimalTolerance: 
    if (value<=0.0||value>1.0e10)
      return false;
    break;
    
  case ClpObjOffset: 
    break;

  case ClpMaxSeconds: 
    if(value>=0)
      value += CoinCpuTime();
    else
      value = -1.0;
    break;

  case ClpLastDblParam:
    return false;
  }
  dblParam_[key] = value;
  return true;
}

//-----------------------------------------------------------------------------

bool
ClpModel::setStrParam(ClpStrParam key, const std::string & value)
{

  switch (key) {
  case ClpProbName:
    break;

  case ClpLastStrParam:
    return false;
  }
  strParam_[key] = value;
  return true;
}
// Useful routines
// Returns resized array and deletes incoming
double * resizeDouble(double * array , int size, int newSize, double fill,
		      bool createArray)
{
  if ((array||createArray)&&size!=newSize) {
    int i;
    double * newArray = new double[newSize];
    if (array)
      memcpy(newArray,array,CoinMin(newSize,size)*sizeof(double));
    delete [] array;
    array = newArray;
    for (i=size;i<newSize;i++) 
      array[i]=fill;
  } 
  return array;
}
// Returns resized array and updates size
double * deleteDouble(double * array , int size, 
		      int number, const int * which,int & newSize)
{
  if (array) {
    int i ;
    char * deleted = new char[size];
    int numberDeleted=0;
    memset(deleted,0,size*sizeof(char));
    for (i=0;i<number;i++) {
      int j = which[i];
      if (j>=0&&j<size&&!deleted[j]) {
	numberDeleted++;
	deleted[j]=1;
      }
    }
    newSize = size-numberDeleted;
    double * newArray = new double[newSize];
    int put=0;
    for (i=0;i<size;i++) {
      if (!deleted[i]) {
	newArray[put++]=array[i];
      }
    }
    delete [] array;
    array = newArray;
    delete [] deleted;
  }
  return array;
}
char * deleteChar(char * array , int size, 
		  int number, const int * which,int & newSize,
		  bool ifDelete)
{
  if (array) {
    int i ;
    char * deleted = new char[size];
    int numberDeleted=0;
    memset(deleted,0,size*sizeof(char));
    for (i=0;i<number;i++) {
      int j = which[i];
      if (j>=0&&j<size&&!deleted[j]) {
	numberDeleted++;
	deleted[j]=1;
      }
    }
    newSize = size-numberDeleted;
    char * newArray = new char[newSize];
    int put=0;
    for (i=0;i<size;i++) {
      if (!deleted[i]) {
	newArray[put++]=array[i];
      }
    }
    if (ifDelete)
      delete [] array;
    array = newArray;
    delete [] deleted;
  }
  return array;
}
// Create empty ClpPackedMatrix
void 
ClpModel::createEmptyMatrix()
{
  delete matrix_;
  CoinPackedMatrix matrix2;
  matrix_=new ClpPackedMatrix(matrix2);
}
// Resizes 
void 
ClpModel::resize (int newNumberRows, int newNumberColumns)
{
  rowActivity_ = resizeDouble(rowActivity_,numberRows_,
			      newNumberRows,0.0,true);
  dual_ = resizeDouble(dual_,numberRows_,
		       newNumberRows,0.0,true);
  rowObjective_ = resizeDouble(rowObjective_,numberRows_,
			       newNumberRows,0.0,false);
  rowLower_ = resizeDouble(rowLower_,numberRows_,
			   newNumberRows,-COIN_DBL_MAX,true);
  rowUpper_ = resizeDouble(rowUpper_,numberRows_,
			   newNumberRows,COIN_DBL_MAX,true);
  columnActivity_ = resizeDouble(columnActivity_,numberColumns_,
				 newNumberColumns,0.0,true);
  reducedCost_ = resizeDouble(reducedCost_,numberColumns_,
			      newNumberColumns,0.0,true);
  if (objective_)
    objective_->resize(newNumberColumns);
  else 
    objective_ = new ClpLinearObjective(NULL,newNumberColumns);
  columnLower_ = resizeDouble(columnLower_,numberColumns_,
			      newNumberColumns,0.0,true);
  columnUpper_ = resizeDouble(columnUpper_,numberColumns_,
			      newNumberColumns,COIN_DBL_MAX,true);
  if (newNumberRows<numberRows_) {
    int * which = new int[numberRows_-newNumberRows];
    int i;
    for (i=newNumberRows;i<numberRows_;i++) 
      which[i-newNumberRows]=i;
    matrix_->deleteRows(numberRows_-newNumberRows,which);
    delete [] which;
  }
  if (numberRows_!=newNumberRows||numberColumns_!=newNumberColumns) {
    // set state back to unknown
    problemStatus_ = -1;
    secondaryStatus_ = 0;
    delete [] ray_;
    ray_ = NULL;
  }
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
  if (status_) {
    unsigned char * tempC = new unsigned char [newNumberColumns+newNumberRows];
    unsigned char * tempR = tempC + newNumberColumns;
    memset(tempC,0,(newNumberColumns+newNumberRows)*sizeof(unsigned char));
    memcpy(tempC,status_,CoinMin(newNumberColumns,numberColumns_)*sizeof(unsigned char));
    memcpy(tempR,status_+numberColumns_,CoinMin(newNumberRows,numberRows_)*sizeof(unsigned char));
    delete [] status_;
    status_ = tempC;
  }
  numberRows_ = newNumberRows;
  if (newNumberColumns<numberColumns_) {
    int * which = new int[numberColumns_-newNumberColumns];
    int i;
    for (i=newNumberColumns;i<numberColumns_;i++) 
      which[i-newNumberColumns]=i;
    matrix_->deleteCols(numberColumns_-newNumberColumns,which);
    delete [] which;
  }
  if (integerType_) {
    char * temp = new char [newNumberColumns];
    memset(temp,0,newNumberColumns*sizeof(char));
    memcpy(temp,integerType_,
	   CoinMin(newNumberColumns,numberColumns_)*sizeof(char));
    delete [] integerType_;
    integerType_ = temp;
  }
  numberColumns_ = newNumberColumns;
  // for now gets rid of names
  lengthNames_ = 0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
}
// Deletes rows
void 
ClpModel::deleteRows(int number, const int * which)
{
  int newSize=0;
  rowActivity_ = deleteDouble(rowActivity_,numberRows_,
			      number, which, newSize);
  dual_ = deleteDouble(dual_,numberRows_,
			      number, which, newSize);
  rowObjective_ = deleteDouble(rowObjective_,numberRows_,
			      number, which, newSize);
  rowLower_ = deleteDouble(rowLower_,numberRows_,
			      number, which, newSize);
  rowUpper_ = deleteDouble(rowUpper_,numberRows_,
			      number, which, newSize);
  matrix_->deleteRows(number,which);
  //matrix_->removeGaps();
  // status
  if (status_) {
    unsigned char * tempR  = (unsigned char *) deleteChar((char *)status_+numberColumns_,
					numberRows_,
					number, which, newSize,false);
    unsigned char * tempC = new unsigned char [numberColumns_+newSize];
    memcpy(tempC,status_,numberColumns_*sizeof(unsigned char));
    memcpy(tempC+numberColumns_,tempR,newSize*sizeof(unsigned char));
    delete [] tempR;
    delete [] status_;
    status_ = tempC;
  }
#if 1
  if (lengthNames_) {
    int i, j, k;
    for (k = 0, j = 0, i = 0; j < number && i < numberRows_; ++i) {
      if (which[j] == i) {
	++j;
      } else {
	rowNames_[k++] = rowNames_[i];
      }
    }
    for ( ; i < numberRows_; ++i) {
      rowNames_[k++] = rowNames_[i];
    }
    rowNames_.erase(rowNames_.begin()+k, rowNames_.end());
  }
#else
  // for now gets rid of names
  lengthNames_ = 0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
#endif
  numberRows_=newSize;
  // set state back to unknown
  problemStatus_ = -1;
  secondaryStatus_ = 0;
  delete [] ray_;
  ray_ = NULL;
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// Deletes columns
void 
ClpModel::deleteColumns(int number, const int * which)
{
  int newSize=0;
  columnActivity_ = deleteDouble(columnActivity_,numberColumns_,
			      number, which, newSize);
  reducedCost_ = deleteDouble(reducedCost_,numberColumns_,
			      number, which, newSize);
  objective_->deleteSome(number, which);
  columnLower_ = deleteDouble(columnLower_,numberColumns_,
			      number, which, newSize);
  columnUpper_ = deleteDouble(columnUpper_,numberColumns_,
			      number, which, newSize);
  // possible matrix is not full
  if (matrix_->getNumCols()<numberColumns_) {
    int * which2 = new int [number];
    int n=0;
    int nMatrix = matrix_->getNumCols();
    for (int i=0;i<number;i++) {
      if (which[i]<nMatrix)
	which2[n++]=which[i];
    }
    matrix_->deleteCols(n,which2);
    delete [] which2;
  } else {
    matrix_->deleteCols(number,which);
  }
  //matrix_->removeGaps();
  // status
  if (status_) {
    unsigned char * tempC  = (unsigned char *) deleteChar((char *)status_,
					numberColumns_,
					number, which, newSize,false);
    unsigned char * temp = new unsigned char [numberRows_+newSize];
    memcpy(temp,tempC,newSize*sizeof(unsigned char));
    memcpy(temp+newSize,status_+numberColumns_,
	   numberRows_*sizeof(unsigned char));
    delete [] tempC;
    delete [] status_;
    status_ = temp;
  }
  integerType_ = deleteChar(integerType_,numberColumns_,
			    number, which, newSize,true);
#if 1
  if (lengthNames_) {
    int i, j, k;
    for (k = 0, j = 0, i = 0; j < number && i < numberColumns_; ++i) {
      if (which[j] == i) {
	++j;
      } else {
	columnNames_[k++] = columnNames_[i];
      }
    }
    for ( ; i < numberColumns_; ++i) {
      columnNames_[k++] = columnNames_[i];
    }
    columnNames_.erase(columnNames_.begin()+k, columnNames_.end());
  }
#else
  // for now gets rid of names
  lengthNames_ = 0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
#endif
  numberColumns_=newSize;
  // set state back to unknown
  problemStatus_ = -1;
  secondaryStatus_ = 0;
  delete [] ray_;
  ray_ = NULL;
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// Add rows
void 
ClpModel::addRows(int number, const double * rowLower, 
		  const double * rowUpper,
		  const int * rowStarts, const int * columns,
		  const double * elements)
{
  // Create a list of CoinPackedVectors
  if (number) {
    CoinPackedVectorBase ** rows=
      new CoinPackedVectorBase * [number];
    int iRow;
    for (iRow=0;iRow<number;iRow++) {
      int iStart = rowStarts[iRow];
      rows[iRow] = 
	new CoinPackedVector(rowStarts[iRow+1]-iStart,
			     columns+iStart,elements+iStart);
    }
    addRows(number, rowLower, rowUpper,
	    rows);
    for (iRow=0;iRow<number;iRow++) 
      delete rows[iRow];
    delete [] rows;
  }
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// Add rows
void 
ClpModel::addRows(int number, const double * rowLower, 
		  const double * rowUpper,
		  const int * rowStarts, 
		  const int * rowLengths, const int * columns,
		  const double * elements)
{
  // Create a list of CoinPackedVectors
  if (number) {
    CoinPackedVectorBase ** rows=
      new CoinPackedVectorBase * [number];
    int iRow;
    for (iRow=0;iRow<number;iRow++) {
      int iStart = rowStarts[iRow];
      rows[iRow] = 
	new CoinPackedVector(rowLengths[iRow],
			     columns+iStart,elements+iStart);
    }
    addRows(number, rowLower, rowUpper,
	    rows);
    for (iRow=0;iRow<number;iRow++) 
      delete rows[iRow];
    delete [] rows;
  }
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
void 
ClpModel::addRows(int number, const double * rowLower, 
		  const double * rowUpper,
		  const CoinPackedVectorBase * const * rows)
{
  if (!number)
    return;
  int numberRowsNow = numberRows_;
  resize(numberRowsNow+number,numberColumns_);
  double * lower = rowLower_+numberRowsNow;
  double * upper = rowUpper_+numberRowsNow;
  int iRow;
  if (rowLower) {
    for (iRow = 0; iRow < number; iRow++) {
      double value = rowLower[iRow];
      if (value<-1.0e20)
	value = -COIN_DBL_MAX;
      lower[iRow]= value;
    }
  } else {
    for (iRow = 0; iRow < number; iRow++) {
      lower[iRow]= -COIN_DBL_MAX;
    }
  }
  if (rowUpper) {
    for (iRow = 0; iRow < number; iRow++) {
      double value = rowUpper[iRow];
      if (value>1.0e20)
	value = COIN_DBL_MAX;
      upper[iRow]= value;
    }
  } else {
    for (iRow = 0; iRow < number; iRow++) {
      upper[iRow]= COIN_DBL_MAX;
    }
  }
  // Deal with matrix

  delete rowCopy_;
  rowCopy_=NULL;
  if (!matrix_)
    createEmptyMatrix();
  matrix_->appendRows(number,rows);
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// Add columns
void 
ClpModel::addColumns(int number, const double * columnLower, 
		     const double * columnUpper,
		     const double * objIn,
		     const int * columnStarts, const int * rows,
		     const double * elements)
{
  // Create a list of CoinPackedVectors
  if (number) {
    CoinPackedVectorBase ** columns=
      new CoinPackedVectorBase * [number];
    int iColumn;
    for (iColumn=0;iColumn<number;iColumn++) {
      int iStart = columnStarts[iColumn];
      columns[iColumn] = 
	new CoinPackedVector(columnStarts[iColumn+1]-iStart,
			     rows+iStart,elements+iStart);
    }
    addColumns(number, columnLower, columnUpper,
	       objIn, columns);
    for (iColumn=0;iColumn<number;iColumn++) 
      delete columns[iColumn];
    delete [] columns;

  }
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// Add columns
void 
ClpModel::addColumns(int number, const double * columnLower, 
		     const double * columnUpper,
		     const double * objIn,
		     const int * columnStarts, 
		     const int * columnLengths, const int * rows,
		     const double * elements)
{
  // Create a list of CoinPackedVectors
  if (number) {
    CoinPackedVectorBase ** columns=
      new CoinPackedVectorBase * [number];
    int iColumn;
    for (iColumn=0;iColumn<number;iColumn++) {
      int iStart = columnStarts[iColumn];
      columns[iColumn] = 
	new CoinPackedVector(columnLengths[iColumn],
			     rows+iStart,elements+iStart);
    }
    addColumns(number, columnLower, columnUpper,
	       objIn, columns);
    for (iColumn=0;iColumn<number;iColumn++) 
      delete columns[iColumn];
    delete [] columns;

  }
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
void 
ClpModel::addColumns(int number, const double * columnLower, 
		     const double * columnUpper,
		     const double * objIn,
		     const CoinPackedVectorBase * const * columns)
{
  if (!number)
    return;
  int numberColumnsNow = numberColumns_;
  resize(numberRows_,numberColumnsNow+number);
  double * lower = columnLower_+numberColumnsNow;
  double * upper = columnUpper_+numberColumnsNow;
  double * obj = objective()+numberColumnsNow;
  int iColumn;
  if (columnLower) {
    for (iColumn = 0; iColumn < number; iColumn++) {
      double value = columnLower[iColumn];
      if (value<-1.0e20)
	value = -COIN_DBL_MAX;
      lower[iColumn]= value;
    }
  } else {
    for (iColumn = 0; iColumn < number; iColumn++) {
      lower[iColumn]= 0.0;
    }
  }
  if (columnUpper) {
    for (iColumn = 0; iColumn < number; iColumn++) {
      double value = columnUpper[iColumn];
      if (value>1.0e20)
	value = COIN_DBL_MAX;
      upper[iColumn]= value;
    }
  } else {
    for (iColumn = 0; iColumn < number; iColumn++) {
      upper[iColumn]= COIN_DBL_MAX;
    }
  }
  if (objIn) {
    for (iColumn = 0; iColumn < number; iColumn++) {
      obj[iColumn] = objIn[iColumn];
    }
  } else {
    for (iColumn = 0; iColumn < number; iColumn++) {
      obj[iColumn]= 0.0;
    }
  }
  // Deal with matrix

  delete rowCopy_;
  rowCopy_=NULL;
  if (!matrix_)
    createEmptyMatrix();
  matrix_->appendCols(number,columns);
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
// chgRowLower
void 
ClpModel::chgRowLower(const double * rowLower) 
{
  int numberRows = numberRows_;
  int iRow;
  if (rowLower) {
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value = rowLower[iRow];
      if (value<-1.0e20)
		 value = -COIN_DBL_MAX;
      rowLower_[iRow]= value;
    }
  } else {
    for (iRow = 0; iRow < numberRows; iRow++) {
      rowLower_[iRow]= -COIN_DBL_MAX;
    }
  }
}
// chgRowUpper
void 
ClpModel::chgRowUpper(const double * rowUpper) 
{
  int numberRows = numberRows_;
  int iRow;
  if (rowUpper) {
    for (iRow = 0; iRow < numberRows; iRow++) {
      double value = rowUpper[iRow];
      if (value>1.0e20)
		 value = COIN_DBL_MAX;
      rowUpper_[iRow]= value;
    }
  } else {
    for (iRow = 0; iRow < numberRows; iRow++) {
      rowUpper_[iRow]= COIN_DBL_MAX;;
    }
  }
}
// chgColumnLower
void 
ClpModel::chgColumnLower(const double * columnLower) 
{
  int numberColumns = numberColumns_;
  int iColumn;
  if (columnLower) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = columnLower[iColumn];
      if (value<-1.0e20)
		 value = -COIN_DBL_MAX;
      columnLower_[iColumn]= value;
    }
  } else {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      columnLower_[iColumn]= 0.0;
    }
  }
}
// chgColumnUpper
void 
ClpModel::chgColumnUpper(const double * columnUpper) 
{
  int numberColumns = numberColumns_;
  int iColumn;
  if (columnUpper) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = columnUpper[iColumn];
      if (value>1.0e20)
		 value = COIN_DBL_MAX;
      columnUpper_[iColumn]= value;
    }
  } else {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      columnUpper_[iColumn]= COIN_DBL_MAX;;
    }
  }
}
// chgObjCoefficients
void 
ClpModel::chgObjCoefficients(const double * objIn) 
{
  double * obj = objective();
  int numberColumns = numberColumns_;
  int iColumn;
  if (objIn) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      obj[iColumn] = objIn[iColumn];
    }
  } else {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      obj[iColumn]= 0.0;
    }
  }
}
// Infeasibility/unbounded ray (NULL returned if none/wrong)
double * 
ClpModel::infeasibilityRay() const
{
  double * array = NULL;
  if (problemStatus_==1&&!secondaryStatus_) 
    array = ClpCopyOfArray(ray_,numberRows_);
  return array;
}
double * 
ClpModel::unboundedRay() const
{
  double * array = NULL;
  if (problemStatus_==2) 
    array = ClpCopyOfArray(ray_,numberColumns_);
  return array;
}
void 
ClpModel::setMaximumIterations(int value)
{
  if(value>=0)
    intParam_[ClpMaxNumIteration]=value;
}
void 
ClpModel::setMaximumSeconds(double value)
{
  if(value>=0)
    dblParam_[ClpMaxSeconds]=value+CoinCpuTime();
  else
    dblParam_[ClpMaxSeconds]=-1.0;
}
// Returns true if hit maximum iterations (or time)
bool 
ClpModel::hitMaximumIterations() const
{
  bool hitMax= (numberIterations_>=maximumIterations());
  if (dblParam_[ClpMaxSeconds]>=0.0&&!hitMax)
    hitMax = (CoinCpuTime()>=dblParam_[ClpMaxSeconds]);
  return hitMax;
}
// Pass in Message handler (not deleted at end)
void 
ClpModel::passInMessageHandler(CoinMessageHandler * handler)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_=false;
  handler_=handler;
}
// Pass in Message handler (not deleted at end) and return current
CoinMessageHandler *
ClpModel::pushMessageHandler(CoinMessageHandler * handler,
			     bool & oldDefault)
{
  CoinMessageHandler * returnValue = handler_;
  oldDefault = defaultHandler_;
  defaultHandler_=false;
  handler_=handler;
  return returnValue;
}
// back to previous message handler
void
ClpModel::popMessageHandler(CoinMessageHandler * oldHandler,bool oldDefault)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_=oldDefault;
  handler_=oldHandler;
}
// Set language
void 
ClpModel::newLanguage(CoinMessages::Language language)
{
  messages_ = ClpMessage(language);
}
// Read an mps file from the given filename
int 
ClpModel::readMps(const char *fileName,
		  bool keepNames,
		  bool ignoreErrors)
{
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
  double time1 = CoinCpuTime(),time2;
  int status=m.readMps(fileName,"");
  m.messageHandler()->setPrefix(savePrefix);
  if (!status||ignoreErrors) {
    loadProblem(*m.getMatrixByCol(),
		m.getColLower(),m.getColUpper(),
		m.getObjCoefficients(),
		m.getRowLower(),m.getRowUpper());
    if (m.integerColumns()) {
      integerType_ = new char[numberColumns_];
      memcpy(integerType_,m.integerColumns(),numberColumns_*sizeof(char));
    } else {
      integerType_ = NULL;
    }
    // get quadratic part
    if (m.reader()->whichSection (  ) == COIN_QUAD_SECTION ) {
      int * start=NULL;
      int * column = NULL;
      double * element = NULL;
      status=m.readQuadraticMps(NULL,start,column,element,2);
      if (!status||ignoreErrors) 
	loadQuadraticObjective(numberColumns_,start,column,element);
      delete [] start;
      delete [] column;
      delete [] element;
    }
    // set problem name
    setStrParam(ClpProbName,m.getProblemName());
    // do names
    if (keepNames) {
      unsigned int maxLength=0;
      int iRow;
      rowNames_ = std::vector<std::string> ();
      columnNames_ = std::vector<std::string> ();
      rowNames_.reserve(numberRows_);
      for (iRow=0;iRow<numberRows_;iRow++) {
	const char * name = m.rowName(iRow);
	maxLength = CoinMax(maxLength,(unsigned int) strlen(name));
	  rowNames_.push_back(name);
      }
      
      int iColumn;
      columnNames_.reserve(numberColumns_);
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	const char * name = m.columnName(iColumn);
	maxLength = CoinMax(maxLength,(unsigned int) strlen(name));
	columnNames_.push_back(name);
      }
      lengthNames_=(int) maxLength;
    } else {
      lengthNames_=0;
    }
    setDblParam(ClpObjOffset,m.objectiveOffset());
    time2 = CoinCpuTime();
    handler_->message(CLP_IMPORT_RESULT,messages_)
      <<fileName
      <<time2-time1<<CoinMessageEol;
  } else {
    // errors
    handler_->message(CLP_IMPORT_ERRORS,messages_)
      <<status<<fileName<<CoinMessageEol;
  }

  return status;
}
bool ClpModel::isPrimalObjectiveLimitReached() const
{
  double limit = 0.0;
  getDblParam(ClpPrimalObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = objectiveValue();
  const double maxmin = optimizationDirection();

  if (problemStatus_ == 0) // optimal
    return maxmin > 0 ? (obj < limit) /*minim*/ : (-obj < limit) /*maxim*/;
  else if (problemStatus_==2)
    return true;
  else
    return false;
}

bool ClpModel::isDualObjectiveLimitReached() const
{

  double limit = 0.0;
  getDblParam(ClpDualObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = objectiveValue();
  const double maxmin = optimizationDirection();

  if (problemStatus_ == 0) // optimal
    return maxmin > 0 ? (obj > limit) /*minim*/ : (-obj > limit) /*maxim*/;
  else if (problemStatus_==1)
    return true;
  else
    return false;

}
void 
ClpModel::copyInIntegerInformation(const char * information)
{
  delete [] integerType_;
  if (information) {
    integerType_ = new char[numberColumns_];
    memcpy(integerType_,information,numberColumns_*sizeof(char));
  } else {
    integerType_ = NULL;
  }
}
// Drops names - makes lengthnames 0 and names empty
void 
ClpModel::dropNames()
{
  lengthNames_=0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
}
// Drop integer informations
void 
ClpModel::deleteIntegerInformation()
{
  delete [] integerType_;
  integerType_ = NULL;
}
/* Return copy of status array (char[numberRows+numberColumns]),
   use delete [] */
unsigned char *  
ClpModel::statusCopy() const
{
  return ClpCopyOfArray(status_,numberRows_+numberColumns_);
}
// Copy in status vector
void 
ClpModel::copyinStatus(const unsigned char * statusArray)
{
  delete [] status_;
  if (statusArray) {
    status_ = new unsigned char [numberRows_+numberColumns_];
    memcpy(status_,statusArray,(numberRows_+numberColumns_)*sizeof(unsigned char));
  } else {
    status_=NULL;
  }
}

// Load up quadratic objective 
void 
ClpModel::loadQuadraticObjective(const int numberColumns, const CoinBigIndex * start,
			      const int * column, const double * element)
{
  assert (numberColumns==numberColumns_);
  assert ((dynamic_cast< ClpLinearObjective*>(objective_)));
  double offset;
  ClpObjective * obj = new ClpQuadraticObjective(objective_->gradient(NULL,NULL,offset,false),
						 numberColumns,
						 start,column,element);
  delete objective_;
  objective_ = obj;

}
void 
ClpModel::loadQuadraticObjective (  const CoinPackedMatrix& matrix)
{
  assert (matrix.getNumCols()==numberColumns_);
  assert ((dynamic_cast< ClpLinearObjective*>(objective_)));
  double offset;
  ClpQuadraticObjective * obj = 
    new ClpQuadraticObjective(objective_->gradient(NULL,NULL,offset,false),
			      numberColumns_,
			      NULL,NULL,NULL);
  delete objective_;
  objective_ = obj;
  obj->loadQuadraticObjective(matrix);
}
// Get rid of quadratic objective
void 
ClpModel::deleteQuadraticObjective()
{
  ClpQuadraticObjective * obj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
  if (obj)
    obj->deleteQuadraticObjective();
}
void 
ClpModel::setObjective(ClpObjective * objective)
{
  delete objective_;
  objective_=objective->clone();
}
// Returns resized array and updates size
double * whichDouble(double * array , int number, const int * which)
{
  double * newArray=NULL;
  if (array&&number) {
    int i ;
    newArray = new double[number];
    for (i=0;i<number;i++) 
      newArray[i]=array[which[i]];
  }
  return newArray;
}
char * whichChar(char * array , int number, const int * which)
{
  char * newArray=NULL;
  if (array&&number) {
    int i ;
    newArray = new char[number];
    for (i=0;i<number;i++) 
      newArray[i]=array[which[i]];
  }
  return newArray;
}
unsigned char * whichUnsignedChar(unsigned char * array , 
				  int number, const int * which)
{
  unsigned char * newArray=NULL;
  if (array&&number) {
    int i ;
    newArray = new unsigned char[number];
    for (i=0;i<number;i++) 
      newArray[i]=array[which[i]];
  }
  return newArray;
}
// Replace Clp Matrix (current is not deleted)
void 
ClpModel::replaceMatrix( ClpMatrixBase * matrix)
{
  matrix_=matrix;
}
// Subproblem constructor
ClpModel::ClpModel ( const ClpModel * rhs,
		     int numberRows, const int * whichRow,
		     int numberColumns, const int * whichColumn,
		     bool dropNames, bool dropIntegers)
{
  defaultHandler_ = rhs->defaultHandler_;
  if (defaultHandler_) 
    handler_ = new CoinMessageHandler(*rhs->handler_);
   else 
    handler_ = rhs->handler_;
  eventHandler_ = rhs->eventHandler_->clone();
  messages_ = rhs->messages_;
  intParam_[ClpMaxNumIteration] = rhs->intParam_[ClpMaxNumIteration];
  intParam_[ClpMaxNumIterationHotStart] = 
    rhs->intParam_[ClpMaxNumIterationHotStart];

  dblParam_[ClpDualObjectiveLimit] = rhs->dblParam_[ClpDualObjectiveLimit];
  dblParam_[ClpPrimalObjectiveLimit] = rhs->dblParam_[ClpPrimalObjectiveLimit];
  dblParam_[ClpDualTolerance] = rhs->dblParam_[ClpDualTolerance];
  dblParam_[ClpPrimalTolerance] = rhs->dblParam_[ClpPrimalTolerance];
  dblParam_[ClpObjOffset] = rhs->dblParam_[ClpObjOffset];
  dblParam_[ClpMaxSeconds] = rhs->dblParam_[ClpMaxSeconds];

  strParam_[ClpProbName] = rhs->strParam_[ClpProbName];

  optimizationDirection_ = rhs->optimizationDirection_;
  objectiveValue_=rhs->objectiveValue_;
  smallElement_ = rhs->smallElement_;
  objectiveScale_ = rhs->objectiveScale_;
  rhsScale_ = rhs->rhsScale_;
  numberIterations_ = rhs->numberIterations_;
  solveType_ = rhs->solveType_;
  problemStatus_ = rhs->problemStatus_;
  secondaryStatus_ = rhs->secondaryStatus_;
  // check valid lists
  int numberBad=0;
  int i;
  for (i=0;i<numberRows;i++)
    if (whichRow[i]<0||whichRow[i]>=rhs->numberRows_)
      numberBad++;
  if (numberBad)
    throw CoinError("bad row list", "subproblem constructor", "ClpModel");
  numberBad=0;
  for (i=0;i<numberColumns;i++)
    if (whichColumn[i]<0||whichColumn[i]>=rhs->numberColumns_)
      numberBad++;
  if (numberBad)
    throw CoinError("bad column list", "subproblem constructor", "ClpModel");
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  userPointer_ = rhs->userPointer_;
  if (!dropNames) {
    unsigned int maxLength=0;
    int iRow;
    rowNames_ = std::vector<std::string> ();
    columnNames_ = std::vector<std::string> ();
    rowNames_.reserve(numberRows_);
    for (iRow=0;iRow<numberRows_;iRow++) {
      rowNames_.push_back(rhs->rowNames_[whichRow[iRow]]);
      maxLength = CoinMax(maxLength,(unsigned int) strlen(rowNames_[iRow].c_str()));
    }
    int iColumn;
    columnNames_.reserve(numberColumns_);
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      columnNames_.push_back(rhs->columnNames_[whichColumn[iColumn]]);
      maxLength = CoinMax(maxLength,(unsigned int) strlen(columnNames_[iColumn].c_str()));
    }
    lengthNames_=(int) maxLength;
  } else {
    lengthNames_ = 0;
    rowNames_ = std::vector<std::string> ();
    columnNames_ = std::vector<std::string> ();
  }
  if (rhs->integerType_&&!dropIntegers) {
    integerType_ = whichChar(rhs->integerType_,numberColumns,whichColumn);
  } else {
    integerType_ = NULL;
  }
  if (rhs->rowActivity_) {
    rowActivity_=whichDouble(rhs->rowActivity_,numberRows,whichRow);
    dual_=whichDouble(rhs->dual_,numberRows,whichRow);
    columnActivity_=whichDouble(rhs->columnActivity_,numberColumns,
				whichColumn);
    reducedCost_=whichDouble(rhs->reducedCost_,numberColumns,
				whichColumn);
  } else {
    rowActivity_=NULL;
    columnActivity_=NULL;
    dual_=NULL;
    reducedCost_=NULL;
  }
  rowLower_=whichDouble(rhs->rowLower_,numberRows,whichRow);
  rowUpper_=whichDouble(rhs->rowUpper_,numberRows,whichRow);
  columnLower_=whichDouble(rhs->columnLower_,numberColumns,whichColumn);
  columnUpper_=whichDouble(rhs->columnUpper_,numberColumns,whichColumn);
  if (rhs->objective_)
    objective_  = rhs->objective_->subsetClone(numberColumns,whichColumn);
  else
    objective_ = NULL;
  rowObjective_=whichDouble(rhs->rowObjective_,numberRows,whichRow);
  // status has to be done in two stages
  status_ = new unsigned char[numberColumns_+numberRows_];
  unsigned char * rowStatus = whichUnsignedChar(rhs->status_+rhs->numberColumns_,
						numberRows_,whichRow);
  unsigned char * columnStatus = whichUnsignedChar(rhs->status_,
						numberColumns_,whichColumn);
  memcpy(status_+numberColumns_,rowStatus,numberRows_);
  delete [] rowStatus;
  memcpy(status_,columnStatus,numberColumns_);
  delete [] columnStatus;
  ray_ = NULL;
  if (problemStatus_==1&&!secondaryStatus_)
    ray_ = whichDouble (rhs->ray_,numberRows,whichRow);
  else if (problemStatus_==2)
    ray_ = whichDouble (rhs->ray_,numberColumns,whichColumn);
  rowScale_ = NULL;
  columnScale_ = NULL;
  scalingFlag_ = rhs->scalingFlag_;
  if (rhs->rowCopy_) {
    rowCopy_ = rhs->rowCopy_->subsetClone(numberRows,whichRow,
					  numberColumns,whichColumn);
  } else {
    rowCopy_=NULL;
  }
  matrix_=NULL;
  if (rhs->matrix_) {
    matrix_ = rhs->matrix_->subsetClone(numberRows,whichRow,
					numberColumns,whichColumn);
  }
  CoinSeedRandom(1234567);
}
// Copies in names
void 
ClpModel::copyNames(std::vector<std::string> & rowNames,
		 std::vector<std::string> & columnNames)
{
  unsigned int maxLength=0;
  int iRow;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
  rowNames_.reserve(numberRows_);
  for (iRow=0;iRow<numberRows_;iRow++) {
    rowNames_.push_back(rowNames[iRow]);
    maxLength = CoinMax(maxLength,(unsigned int) strlen(rowNames_[iRow].c_str()));
  }
  int iColumn;
  columnNames_.reserve(numberColumns_);
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    columnNames_.push_back(columnNames[iColumn]);
    maxLength = CoinMax(maxLength,(unsigned int) strlen(columnNames_[iColumn].c_str()));
  }
  lengthNames_=(int) maxLength;
}
// Dual objective limit
void 
ClpModel::setDualObjectiveLimit(double value)
{
  dblParam_[ClpDualObjectiveLimit]=value;
}
// Objective offset
void 
ClpModel::setObjectiveOffset(double value)
{
  dblParam_[ClpObjOffset]=value;
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpDataSave::ClpDataSave () 
{
  dualBound_ = 0.0;
  infeasibilityCost_ = 0.0;
  sparseThreshold_ = 0;
  perturbation_ = 0;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpDataSave::ClpDataSave (const ClpDataSave & rhs) 
{  
  dualBound_ = rhs.dualBound_;
  infeasibilityCost_ = rhs.infeasibilityCost_;
  sparseThreshold_ = rhs.sparseThreshold_;
  perturbation_ = rhs.perturbation_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpDataSave::~ClpDataSave ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpDataSave &
ClpDataSave::operator=(const ClpDataSave& rhs)
{
  if (this != &rhs) {
    dualBound_ = rhs.dualBound_;
    infeasibilityCost_ = rhs.infeasibilityCost_;
    sparseThreshold_ = rhs.sparseThreshold_;
    perturbation_ = rhs.perturbation_;
  }
  return *this;
}
// Solve a problem with no elements - return status
int ClpModel::emptyProblem(int * infeasNumber, double * infeasSum,bool printMessage)
{
  if (printMessage)
    handler_->message(CLP_EMPTY_PROBLEM,messages_)
      <<numberRows_
      <<numberColumns_
      <<0
      <<CoinMessageEol;
  int returnCode=0;
  if (numberRows_||numberColumns_) {
    if (!status_) {
      status_ = new unsigned char[numberRows_+numberColumns_];
      memset(status_,0,numberRows_+numberColumns_);
    }
  }
  // status is set directly (as can be used by Interior methods)
  // check feasible
  int numberPrimalInfeasibilities=0;
  double sumPrimalInfeasibilities=0.0;
  int numberDualInfeasibilities=0;
  double sumDualInfeasibilities=0.0;
  if (numberRows_) {
    for (int i=0;i<numberRows_;i++) {
      dual_[i]=0.0;
      if (rowLower_[i]<=rowUpper_[i]) {
	if (rowLower_[i]>-1.0e30||rowUpper_[i]<1.0e30) {
	  if (fabs(rowLower_[i])<fabs(rowUpper_[i]))
	    rowActivity_[i]=rowLower_[i];
	  else
	    rowActivity_[i]=rowUpper_[i];
	} else {
	  rowActivity_[i]=0.0;
	}
      } else {
	rowActivity_[i]=0.0;
	numberPrimalInfeasibilities++;
	sumPrimalInfeasibilities += rowLower_[i]-rowUpper_[i];
	returnCode=1;
      }
      status_[i+numberColumns_]=1;
    }
  }
  objectiveValue_=0.0;
  if (numberColumns_) {
    const double * cost = objective();
    for (int i=0;i<numberColumns_;i++) {
      reducedCost_[i]=cost[i];
      double objValue = cost[i]*optimizationDirection_;
      if (columnLower_[i]<=columnUpper_[i]) {
	if (columnLower_[i]>-1.0e30||columnUpper_[i]<1.0e30) {
	  if (!objValue) {
	    if (fabs(columnLower_[i])<fabs(columnUpper_[i])) {
	      columnActivity_[i]=columnLower_[i];
	      status_[i]=3;
	    } else {
	      columnActivity_[i]=columnUpper_[i];
	      status_[i]=2;
	    }
	  } else if (objValue>0.0) {
	    if (columnLower_[i]>-1.0e30) {
	      columnActivity_[i]=columnLower_[i];
	      status_[i]=3;
	    } else {
	      columnActivity_[i]=columnUpper_[i];
	      status_[i]=2;
	      numberDualInfeasibilities++;;
	      sumDualInfeasibilities += fabs(objValue);
	      returnCode |= 2;
	    }
	    objectiveValue_ += columnActivity_[i]*objValue;
	  } else {
	    if (columnUpper_[i]<1.0e30) {
	      columnActivity_[i]=columnUpper_[i];
	      status_[i]=2;
	    } else {
	      columnActivity_[i]=columnLower_[i];
	      status_[i]=3;
	      numberDualInfeasibilities++;;
	      sumDualInfeasibilities += fabs(objValue);
	      returnCode |= 2;
	    }
	    objectiveValue_ += columnActivity_[i]*objValue;
	  }
	} else {
	  columnActivity_[i]=0.0;
	  if (objValue) {
	    numberDualInfeasibilities++;;
	    sumDualInfeasibilities += fabs(objValue);
	    returnCode |= 2;
	  }
	  status_[i]=0;
	}
      } else {
	if (fabs(columnLower_[i])<fabs(columnUpper_[i])) {
	  columnActivity_[i]=columnLower_[i];
	  status_[i]=3;
	} else {
	  columnActivity_[i]=columnUpper_[i];
	  status_[i]=2;
	}
	numberPrimalInfeasibilities++;
	sumPrimalInfeasibilities += columnLower_[i]-columnUpper_[i];
	returnCode |= 1;
      }
    }
  }
  objectiveValue_ *= optimizationDirection_;
  if (infeasNumber) {
    infeasNumber[0]=numberDualInfeasibilities;
    infeasSum[0]=sumDualInfeasibilities;
    infeasNumber[1]=numberPrimalInfeasibilities;
    infeasSum[1]=sumPrimalInfeasibilities;
  }
  if (returnCode==3) 
    returnCode=4;
  return returnCode;
}
/* Write the problem in MPS format to the specified file.
   
Row and column names may be null.
formatType is
<ul>
<li> 0 - normal
<li> 1 - extra accuracy 
<li> 2 - IEEE hex (later)
</ul>

Returns non-zero on I/O error
*/
int 
ClpModel::writeMps(const char *filename, 
		   int formatType,int numberAcross,
		   double objSense) const 
{
  
  // Get multiplier for objective function - default 1.0
  double * objective = new double[numberColumns_];
  memcpy(objective,getObjCoefficients(),numberColumns_*sizeof(double));
  if (objSense*getObjSense()<0.0) {
    for (int i = 0; i < numberColumns_; ++i) 
      objective [i] = - objective[i];
  }
  
  char ** rowNames = NULL;
  char ** columnNames = NULL;
  if (lengthNames()) {
    rowNames = new char * [numberRows_];
    for (int iRow=0;iRow<numberRows_;iRow++) {
      rowNames[iRow] = 
	strdup(rowName(iRow).c_str());
#ifdef STRIPBLANKS
      char * xx = rowNames[iRow];
      int i;
      int length = strlen(xx);
      int n=0;
      for (i=0;i<length;i++) {
	if (xx[i]!=' ')
	  xx[n++]=xx[i];
      }
      xx[n]='\0';
#endif
    }
    
    columnNames = new char * [numberColumns_];
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      columnNames[iColumn] = 
	strdup(columnName(iColumn).c_str());
#ifdef STRIPBLANKS
      char * xx = columnNames[iColumn];
      int i;
      int length = strlen(xx);
      int n=0;
      for (i=0;i<length;i++) {
	if (xx[i]!=' ')
	  xx[n++]=xx[i];
      }
      xx[n]='\0';
#endif
    }
  }
  CoinMpsIO writer;
  writer.passInMessageHandler(handler_);
  writer.setMpsData(*(matrix_->getPackedMatrix()), COIN_DBL_MAX,
		    getColLower(), getColUpper(),
		    objective,
		    (const char*) 0 /*integrality*/,
		    getRowLower(), getRowUpper(),
		    columnNames, rowNames);
  // Pass in array saying if each variable integer
  writer.copyInIntegerInformation(integerInformation());
  writer.setObjectiveOffset(objectiveOffset());
  delete [] objective;
  // allow for quadratic objective
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  CoinPackedMatrix * quadratic=NULL;
  if (quadraticObj) 
    quadratic = quadraticObj->quadraticObjective();
  return writer.writeMps(filename, 0 /* do not gzip it*/, formatType, numberAcross,
			 quadratic);
  if (rowNames) {
    for (int iRow=0;iRow<numberRows_;iRow++) {
      free(rowNames[iRow]);
    }
    delete [] rowNames;
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      free(columnNames[iColumn]);
    }
    delete [] columnNames;
  }
}
// Pass in Event handler (cloned and deleted at end)
void 
ClpModel::passInEventHandler(const ClpEventHandler * eventHandler)
{
  delete eventHandler_;
  eventHandler_ = eventHandler->clone();
}
// Sets or unsets scaling, 0 -off, 1 on, 2 dynamic(later)
void 
ClpModel::scaling(int mode)
{
  if (mode>0&&mode<4) {
    scalingFlag_=mode;
  } else if (!mode) {
    scalingFlag_=0;
    delete [] rowScale_;
    rowScale_ = NULL;
    delete [] columnScale_;
    columnScale_ = NULL;
  }
}
void 
ClpModel::times(double scalar,
		  const double * x, double * y) const
{
  if (rowScale_)
    matrix_->times(scalar,x,y,rowScale_,columnScale_);
  else
    matrix_->times(scalar,x,y);
}
void 
ClpModel::transposeTimes(double scalar,
			   const double * x, double * y) const 
{
  if (rowScale_)
    matrix_->transposeTimes(scalar,x,y,rowScale_,columnScale_);
  else
    matrix_->transposeTimes(scalar,x,y);
}
// Does much of scaling
void 
ClpModel::gutsOfScaling()
{
  int i;
  if (rowObjective_) {
    for (i=0;i<numberRows_;i++) 
      rowObjective_[i] /= rowScale_[i];
  }
  for (i=0;i<numberRows_;i++) {
    double multiplier = rowScale_[i];
    double inverseMultiplier = 1.0/multiplier;
    rowActivity_[i] *= multiplier;
    dual_[i] *= inverseMultiplier;
    if (rowLower_[i]>-1.0e30)
      rowLower_[i] *= multiplier;
    else
      rowLower_[i] = -COIN_DBL_MAX;
    if (rowUpper_[i]<1.0e30)
      rowUpper_[i] *= multiplier;
    else
      rowUpper_[i] = COIN_DBL_MAX;
  }
  for (i=0;i<numberColumns_;i++) {
    double multiplier = 1.0/columnScale_[i];
    columnActivity_[i] *= multiplier;
    reducedCost_[i] *= columnScale_[i];
    if (columnLower_[i]>-1.0e30)
      columnLower_[i] *= multiplier;
    else
      columnLower_[i] = -COIN_DBL_MAX;
    if (columnUpper_[i]<1.0e30)
      columnUpper_[i] *= multiplier;
    else
      columnUpper_[i] = COIN_DBL_MAX;
    
  }
  //now replace matrix
  //and objective
  matrix_->reallyScale(rowScale_,columnScale_);
  objective_->reallyScale(columnScale_);
}
/* If we constructed a "really" scaled model then this reverses the operation.
      Quantities may not be exactly as they were before due to rounding errors */
void 
ClpModel::unscale()
{
  if (rowScale_) {
    int i;
    // reverse scaling
    for (i=0;i<numberRows_;i++) 
      rowScale_[i] = 1.0/rowScale_[i];
    for (i=0;i<numberColumns_;i++) 
      columnScale_[i] = 1.0/columnScale_[i];
    gutsOfScaling();
  }
  
  scalingFlag_=0;
  delete [] rowScale_;
  rowScale_ = NULL;
  delete [] columnScale_;
  columnScale_ = NULL;
}
