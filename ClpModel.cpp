// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cmath>
#include <cassert>
#include <cfloat>
#include <string>
#include <cstdio>
#include <iostream>

#include <time.h>

#include "CoinPragma.hpp"
#ifndef _MSC_VER
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include "CoinHelperFunctions.hpp"
#include "ClpModel.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinMpsIO.hpp"
#include "ClpMessage.hpp"
#include "ClpLinearObjective.hpp"

static double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow;        /* clock_t is same as int */
  
  ticksnow = (unsigned int)clock();
  
  cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}
//#############################################################################

ClpModel::ClpModel () :

  optimizationDirection_(1),
  objectiveValue_(0.0),
  smallElement_(1.0e-20),
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
  status_(NULL),
  integerType_(NULL),
  numberIterations_(0),
  solveType_(0),
  problemStatus_(-1),
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

  strParam_[ClpProbName] = "ClpDefaultName";
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = ClpMessage();
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
  delete [] integerType_;
  integerType_ = NULL;
  delete [] status_;
  status_=NULL;
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
void ClpModel::setOptimizationDirection( int value) 
{
  if (value>=-1&&value<=1)
    optimizationDirection_=value;
}
void
ClpModel::gutsOfLoadModel (int numberRows, int numberColumns, 
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
				const double * rowObjective)
{
  gutsOfDelete();
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
ClpModel::ClpModel(const ClpModel &rhs) :
  optimizationDirection_(rhs.optimizationDirection_),
  numberRows_(rhs.numberRows_),
  numberColumns_(rhs.numberColumns_)
{
  gutsOfCopy(rhs);
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
  messages_ = rhs.messages_;
  intParam_[ClpMaxNumIteration] = rhs.intParam_[ClpMaxNumIteration];
  intParam_[ClpMaxNumIterationHotStart] = 
    rhs.intParam_[ClpMaxNumIterationHotStart];

  dblParam_[ClpDualObjectiveLimit] = rhs.dblParam_[ClpDualObjectiveLimit];
  dblParam_[ClpPrimalObjectiveLimit] = rhs.dblParam_[ClpPrimalObjectiveLimit];
  dblParam_[ClpDualTolerance] = rhs.dblParam_[ClpDualTolerance];
  dblParam_[ClpPrimalTolerance] = rhs.dblParam_[ClpPrimalTolerance];
  dblParam_[ClpObjOffset] = rhs.dblParam_[ClpObjOffset];

  strParam_[ClpProbName] = rhs.strParam_[ClpProbName];

  objectiveValue_=rhs.objectiveValue_;
  smallElement_ = rhs.smallElement_;
  numberIterations_ = rhs.numberIterations_;
  solveType_ = rhs.solveType_;
  problemStatus_ = rhs.problemStatus_;
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
    if (rhs.objective_)
      objective_  = rhs.objective_->clone();
    else
      objective_ = NULL;
    rowObjective_ = ClpCopyOfArray ( rhs.rowObjective_, numberRows_ );
    status_ = ClpCopyOfArray( rhs.status_,numberColumns_+numberRows_);
    ray_ = NULL;
    if (problemStatus_==1)
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
    rowCopy_ = rhs.rowCopy_;
    ray_ = rhs.ray_;
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
      memcpy(newArray,array,min(newSize,size)*sizeof(double));
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
    delete [] ray_;
    ray_ = NULL;
  }
  if (status_) {
    unsigned char * tempC = new unsigned char [newNumberColumns+newNumberRows];
    unsigned char * tempR = tempC + newNumberColumns;
    memset(tempC,0,(newNumberColumns+newNumberRows)*sizeof(unsigned char));
    memcpy(tempC,status_,min(newNumberColumns,numberColumns_)*sizeof(unsigned char));
    memcpy(tempR,status_+numberColumns_,min(newNumberRows,numberRows_)*sizeof(unsigned char));
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
	   min(newNumberColumns,numberColumns_)*sizeof(char));
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
  numberRows_=newSize;
  // set state back to unknown
  problemStatus_ = -1;
  delete [] ray_;
  ray_ = NULL;
  // for now gets rid of names
  lengthNames_ = 0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
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
  matrix_->deleteCols(number,which);
  // status
  if (status_) {
    unsigned char * tempC  = (unsigned char *) deleteChar((char *)status_,
					numberColumns_,
					number, which, newSize,false);
    unsigned char * temp = new unsigned char [numberRows_+newSize];
    memcpy(temp,status_,newSize*sizeof(unsigned char));
    memcpy(temp+newSize,status_+numberColumns_,
	   numberRows_*sizeof(unsigned char));
    delete [] tempC;
    delete [] status_;
    status_ = temp;
  }
  numberColumns_=newSize;
  // set state back to unknown
  problemStatus_ = -1;
  delete [] ray_;
  ray_ = NULL;
  // for now gets rid of names
  lengthNames_ = 0;
  rowNames_ = std::vector<std::string> ();
  columnNames_ = std::vector<std::string> ();
  integerType_ = deleteChar(integerType_,numberColumns_,
			    number, which, newSize,true);
}
// Infeasibility/unbounded ray (NULL returned if none/wrong)
double * 
ClpModel::infeasibilityRay() const
{
  double * array = NULL;
  if (problemStatus_==1) 
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
// Pass in Message handler (not deleted at end)
void 
ClpModel::passInMessageHandler(CoinMessageHandler * handler)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_=false;
  handler_=handler;
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
  if (fileName=="-") {
    // stdin
    canOpen=true;
    fileName = "-";
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
  double time1 = cpuTime(),time2;
  int status=m.readMps(fileName,"");
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
    // set problem name
    setStrParam(ClpProbName,m.getProblemName());
    // do names
    if (keepNames) {
      unsigned int maxLength=0;
      int iRow;
      rowNames_.reserve(numberRows_);
      for (iRow=0;iRow<numberRows_;iRow++) {
	const char * name = m.rowName(iRow);
	maxLength = max(maxLength,(unsigned int) strlen(name));
	  rowNames_.push_back(name);
      }
      
      int iColumn;
      columnNames_.reserve(numberColumns_);
      for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	const char * name = m.columnName(iColumn);
	maxLength = max(maxLength,(unsigned int) strlen(name));
	columnNames_.push_back(name);
      }
      lengthNames_=(int) maxLength;
    } else {
      lengthNames_=0;
    }
    setDblParam(ClpObjOffset,m.objectiveOffset());
    time2 = cpuTime();
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
  const int maxmin = optimizationDirection();

  if (problemStatus_ == 0) // optimal
    return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
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
  const int maxmin = optimizationDirection();

  if (problemStatus_ == 0) // optimal
    return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
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
