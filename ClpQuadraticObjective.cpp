// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpModel.hpp"
#include "ClpQuadraticObjective.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpQuadraticObjective::ClpQuadraticObjective () 
: ClpObjective()
{
  type_=2;
  objective_=NULL;
  quadraticObjective_=NULL;
  gradient_ = NULL;
  numberColumns_=0;
  numberExtendedColumns_=0;
}

//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
ClpQuadraticObjective::ClpQuadraticObjective (const double * objective , 
					      int numberColumns,
					      const CoinBigIndex * start,
					      const int * column, const double * element,
					      int numberExtendedColumns)
  : ClpObjective()
{
  type_=2;
  numberColumns_ = numberColumns;
  if (numberExtendedColumns>=0)
    numberExtendedColumns_= max(numberColumns_,numberExtendedColumns);
  else
    numberExtendedColumns_= numberColumns_;
  if (objective) {
    objective_ = new double [numberExtendedColumns_];
    memcpy(objective_,objective,numberColumns_*sizeof(double));
    memset(objective_+numberColumns_,0,(numberExtendedColumns_-numberColumns_)*sizeof(double));
  } else {
    objective_ = new double [numberExtendedColumns_];
    memset(objective_,0,numberExtendedColumns_*sizeof(double));
  }
  if (start) 
    quadraticObjective_ = new CoinPackedMatrix(true,numberColumns,numberColumns,
					     start[numberColumns],element,column,start,NULL);
  else
  quadraticObjective_=NULL;
  gradient_ = NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpQuadraticObjective::ClpQuadraticObjective (const ClpQuadraticObjective & rhs) 
: ClpObjective(rhs)
{  
  numberColumns_=rhs.numberColumns_;
  numberExtendedColumns_=rhs.numberExtendedColumns_;
  if (rhs.objective_) {
    objective_ = new double [numberExtendedColumns_];
    memcpy(objective_,rhs.objective_,numberExtendedColumns_*sizeof(double));
  } else {
    objective_=NULL;
  }
  if (rhs.gradient_) {
    gradient_ = new double [numberExtendedColumns_];
    memcpy(gradient_,rhs.gradient_,numberExtendedColumns_*sizeof(double));
  } else {
    gradient_=NULL;
  }
  if (rhs.quadraticObjective_) 
    quadraticObjective_ = new CoinPackedMatrix(*rhs.quadraticObjective_);
  else 
    quadraticObjective_=NULL;
}
/* Subset constructor.  Duplicates are allowed
   and order is as given.
*/
ClpQuadraticObjective::ClpQuadraticObjective (const ClpQuadraticObjective &rhs,
					int numberColumns, 
					const int * whichColumn) 
: ClpObjective(rhs)
{
  objective_=NULL;
  int extra = rhs.numberExtendedColumns_-rhs.numberColumns_;
  numberColumns_=0;
  if (numberColumns>0) {
    // check valid lists
    int numberBad=0;
    int i;
    for (i=0;i<numberColumns;i++)
      if (whichColumn[i]<0||whichColumn[i]>=rhs.numberColumns_)
	numberBad++;
    if (numberBad)
      throw CoinError("bad column list", "subset constructor", 
		      "ClpQuadraticObjective");
    numberColumns_ = numberColumns;
    numberExtendedColumns_ = numberColumns+extra;
    objective_ = new double[numberExtendedColumns_];
    for (i=0;i<numberColumns_;i++) 
      objective_[i]=rhs.objective_[whichColumn[i]];
    memcpy(objective_+numberColumns_,rhs.objective_+rhs.numberColumns_,
	   (numberExtendedColumns_-numberColumns_)*sizeof(double));
    if (rhs.gradient_) {
      gradient_ = new double[numberExtendedColumns_];
      for (i=0;i<numberColumns_;i++) 
	gradient_[i]=rhs.gradient_[whichColumn[i]];
      memcpy(gradient_+numberColumns_,rhs.gradient_+rhs.numberColumns_,
	     (numberExtendedColumns_-numberColumns_)*sizeof(double));
    }
  }
  if (rhs.quadraticObjective_) {
    quadraticObjective_ = new CoinPackedMatrix(*rhs.quadraticObjective_,
					       numberColumns,whichColumn,
					       numberColumns,whichColumn);
  } else {
    quadraticObjective_=NULL;
  }
}
  

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpQuadraticObjective::~ClpQuadraticObjective ()
{
  delete [] objective_;
  delete [] gradient_;
  delete quadraticObjective_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpQuadraticObjective &
ClpQuadraticObjective::operator=(const ClpQuadraticObjective& rhs)
{
  if (this != &rhs) {
    delete quadraticObjective_;
    quadraticObjective_ = NULL;
    ClpObjective::operator=(rhs);
    numberColumns_=rhs.numberColumns_;
    numberExtendedColumns_=rhs.numberExtendedColumns_;
    if (rhs.objective_) {
      objective_ = new double [numberExtendedColumns_];
      memcpy(objective_,rhs.objective_,numberExtendedColumns_*sizeof(double));
    } else {
      objective_=NULL;
    }
    if (rhs.gradient_) {
      gradient_ = new double [numberExtendedColumns_];
      memcpy(gradient_,rhs.gradient_,numberExtendedColumns_*sizeof(double));
    } else {
      gradient_=NULL;
    }
    if (rhs.quadraticObjective_) {
      quadraticObjective_ = new CoinPackedMatrix(*rhs.quadraticObjective_);
    } else {
      quadraticObjective_=NULL;
    }
  }
  return *this;
}

// Returns gradient
double *  
ClpQuadraticObjective::gradient(const double * solution, double & offset)
{
  offset=0.0;
  if (!quadraticObjective_||!solution) {
    return objective_;
  } else {
    if (!gradient_) 
      gradient_ = new double[numberExtendedColumns_];
    const int * columnQuadratic = quadraticObjective_->getIndices();
    const int * columnQuadraticStart = quadraticObjective_->getVectorStarts();
    const int * columnQuadraticLength = quadraticObjective_->getVectorLengths();
    const double * quadraticElement = quadraticObjective_->getElements();
    double offset=0.0;
    memcpy(gradient_,objective_,numberExtendedColumns_*sizeof(double));
    int iColumn;
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      double valueI = solution[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution[jColumn];
	double elementValue = quadraticElement[j];
	elementValue *= 0.5;
	offset += valueI*valueJ*elementValue;
	double gradientI = valueJ*elementValue;
	double gradientJ = valueI*elementValue;
	offset -= gradientI*valueI;
	gradient_[iColumn] += gradientI;
	offset -= gradientJ*valueJ;
	gradient_[jColumn] += gradientJ;
      }
    }
    return gradient_;
  }
}
  
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpObjective * ClpQuadraticObjective::clone() const
{
  return new ClpQuadraticObjective(*this);
}
/* Subset clone.  Duplicates are allowed
   and order is as given.
*/
ClpObjective * 
ClpQuadraticObjective::subsetClone (int numberColumns, 
			   const int * whichColumns) const
{
  return new ClpQuadraticObjective(*this, numberColumns, whichColumns);
}
// Resize objective
void 
ClpQuadraticObjective::resize(int newNumberColumns)
{
  if (numberColumns_!=newNumberColumns) {
    int newExtended = newNumberColumns + (numberExtendedColumns_-numberColumns_);
    int i;
    double * newArray = new double[newExtended];
    if (objective_)
      memcpy(newArray,objective_,
	     min(newExtended,numberExtendedColumns_)*sizeof(double));
    delete [] objective_;
    objective_ = newArray;
    for (i=numberColumns_;i<newNumberColumns;i++) 
      objective_[i]=0.0;
    if (gradient_) {
      newArray = new double[newExtended];
      if (gradient_)
	memcpy(newArray,gradient_,
	       min(newExtended,numberExtendedColumns_)*sizeof(double));
      delete [] gradient_;
      gradient_ = newArray;
      for (i=numberColumns_;i<newNumberColumns;i++) 
	gradient_[i]=0.0;
    }
    if (quadraticObjective_) {
      if (newNumberColumns<numberColumns_) {
	int * which = new int[numberColumns_-newNumberColumns];
	int i;
	for (i=newNumberColumns;i<numberColumns_;i++) 
	  which[i-newNumberColumns]=i;
	quadraticObjective_->deleteCols(numberColumns_-newNumberColumns,which);
	quadraticObjective_->deleteRows(numberColumns_-newNumberColumns,which);
	delete [] which;
      } else {
	quadraticObjective_->setDimensions(newNumberColumns,newNumberColumns);
      }
    }
    numberColumns_ = newNumberColumns;
    numberExtendedColumns_ = newExtended;
  } 
  
}
// Delete columns in  objective
void 
ClpQuadraticObjective::deleteSome(int numberToDelete, const int * which) 
{
  int newNumberColumns = numberColumns_-numberToDelete;
  int newExtended = numberExtendedColumns_ - numberToDelete;
  if (objective_) {
    int i ;
    char * deleted = new char[numberColumns_];
    int numberDeleted=0;
    memset(deleted,0,numberColumns_*sizeof(char));
    for (i=0;i<numberToDelete;i++) {
      int j = which[i];
      if (j>=0&&j<numberColumns_&&!deleted[j]) {
	numberDeleted++;
	deleted[j]=1;
      }
    }
    newNumberColumns = numberColumns_-numberDeleted;
    newExtended = numberExtendedColumns_ - numberDeleted;
    double * newArray = new double[newExtended];
    int put=0;
    for (i=0;i<numberColumns_;i++) {
      if (!deleted[i]) {
	newArray[put++]=objective_[i];
      }
    }
    delete [] objective_;
    objective_ = newArray;
    delete [] deleted;
    memcpy(objective_+newNumberColumns,objective_+numberColumns_,
	   (numberExtendedColumns_-numberColumns_)*sizeof(double));
  }
  if (gradient_) {
    int i ;
    char * deleted = new char[numberColumns_];
    int numberDeleted=0;
    memset(deleted,0,numberColumns_*sizeof(char));
    for (i=0;i<numberToDelete;i++) {
      int j = which[i];
      if (j>=0&&j<numberColumns_&&!deleted[j]) {
	numberDeleted++;
	deleted[j]=1;
      }
    }
    newNumberColumns = numberColumns_-numberDeleted;
    newExtended = numberExtendedColumns_ - numberDeleted;
    double * newArray = new double[newExtended];
    int put=0;
    for (i=0;i<numberColumns_;i++) {
      if (!deleted[i]) {
	newArray[put++]=gradient_[i];
      }
    }
    delete [] gradient_;
    gradient_ = newArray;
    delete [] deleted;
    memcpy(gradient_+newNumberColumns,gradient_+numberColumns_,
	   (numberExtendedColumns_-numberColumns_)*sizeof(double));
  }
  numberColumns_ = newNumberColumns;
  numberExtendedColumns_ = newExtended;
  if (quadraticObjective_) {
    quadraticObjective_->deleteCols(numberToDelete,which);
    quadraticObjective_->deleteRows(numberToDelete,which);
  }
}

// Load up quadratic objective 
void 
ClpQuadraticObjective::loadQuadraticObjective(const int numberColumns, const CoinBigIndex * start,
			      const int * column, const double * element,int numberExtended)
{
  delete quadraticObjective_;
  quadraticObjective_ = new CoinPackedMatrix(true,numberColumns,numberColumns,
					     start[numberColumns],element,column,start,NULL);
  numberColumns_=numberColumns;
  if (numberExtended>numberExtendedColumns_) {
    if (objective_) {
      // make correct size
      double * newArray = new double[numberExtended];
      memcpy(newArray,objective_,numberColumns_*sizeof(double));
      delete [] objective_;
      objective_ = newArray;
      memset(objective_+numberColumns_,0,(numberExtended-numberColumns_)*sizeof(double));
    }
    if (gradient_) {
      // make correct size
      double * newArray = new double[numberExtended];
      memcpy(newArray,gradient_,numberColumns_*sizeof(double));
      delete [] gradient_;
      gradient_ = newArray;
      memset(gradient_+numberColumns_,0,(numberExtended-numberColumns_)*sizeof(double));
    }
    numberExtendedColumns_ = numberExtended;
  } else {
    numberExtendedColumns_ = numberColumns_;
  }
}
void 
ClpQuadraticObjective::loadQuadraticObjective (  const CoinPackedMatrix& matrix)
{
  delete quadraticObjective_;
  quadraticObjective_ = new CoinPackedMatrix(matrix);
}
// Get rid of quadratic objective
void 
ClpQuadraticObjective::deleteQuadraticObjective()
{
  delete quadraticObjective_;
  quadraticObjective_ = NULL;
}
