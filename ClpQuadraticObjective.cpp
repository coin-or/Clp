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
}

//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
ClpQuadraticObjective::ClpQuadraticObjective (const double * objective , 
					      int numberColumns,
					      const CoinBigIndex * start,
					      const int * column, const double * element)
  : ClpObjective()
{
  type_=2;
  numberColumns_ = numberColumns;
  if (objective) {
    objective_ = new double [numberColumns_];
    memcpy(objective_,objective,numberColumns_*sizeof(double));
  } else {
    objective_ = new double [numberColumns_];
    memset(objective_,0,numberColumns_*sizeof(double));
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
  if (rhs.objective_) {
    objective_ = new double [numberColumns_];
    memcpy(objective_,rhs.objective_,numberColumns_*sizeof(double));
  } else {
    objective_=NULL;
  }
  if (rhs.gradient_) {
    gradient_ = new double [numberColumns_];
    memcpy(gradient_,rhs.gradient_,numberColumns_*sizeof(double));
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
    objective_ = new double[numberColumns_];
    for (i=0;i<numberColumns_;i++) 
      objective_[i]=rhs.objective_[whichColumn[i]];
    if (rhs.gradient_) {
      gradient_ = new double[numberColumns_];
      for (i=0;i<numberColumns_;i++) 
	gradient_[i]=rhs.gradient_[whichColumn[i]];
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
    if (rhs.objective_) {
      objective_ = new double [numberColumns_];
      memcpy(objective_,rhs.objective_,numberColumns_*sizeof(double));
    } else {
      objective_=NULL;
    }
    if (rhs.gradient_) {
      gradient_ = new double [numberColumns_];
      memcpy(gradient_,rhs.gradient_,numberColumns_*sizeof(double));
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
      gradient_ = new double[numberColumns_];
    const int * columnQuadratic = quadraticObjective_->getIndices();
    const int * columnQuadraticStart = quadraticObjective_->getVectorStarts();
    const int * columnQuadraticLength = quadraticObjective_->getVectorLengths();
    const double * quadraticElement = quadraticObjective_->getElements();
    double offset=0.0;
    memcpy(gradient_,objective_,numberColumns_*sizeof(double));
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
    int i;
    double * newArray = new double[newNumberColumns];
    if (objective_)
      memcpy(newArray,objective_,
	     min(newNumberColumns,numberColumns_)*sizeof(double));
    delete [] objective_;
    objective_ = newArray;
    for (i=numberColumns_;i<newNumberColumns;i++) 
      objective_[i]=0.0;
    if (gradient_) {
      newArray = new double[newNumberColumns];
      if (gradient_)
	memcpy(newArray,gradient_,
	       min(newNumberColumns,numberColumns_)*sizeof(double));
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
  } 
  
}
// Delete columns in  objective
void 
ClpQuadraticObjective::deleteSome(int numberToDelete, const int * which) 
{
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
    int newNumberColumns = numberColumns_-numberDeleted;
    double * newArray = new double[newNumberColumns];
    int put=0;
    for (i=0;i<numberColumns_;i++) {
      if (!deleted[i]) {
	newArray[put++]=objective_[i];
      }
    }
    delete [] objective_;
    objective_ = newArray;
    delete [] deleted;
    numberColumns_ = newNumberColumns;
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
    int newNumberColumns = numberColumns_-numberDeleted;
    double * newArray = new double[newNumberColumns];
    int put=0;
    for (i=0;i<numberColumns_;i++) {
      if (!deleted[i]) {
	newArray[put++]=gradient_[i];
      }
    }
    delete [] gradient_;
    gradient_ = newArray;
    delete [] deleted;
    numberColumns_ = newNumberColumns;
  }
  if (quadraticObjective_) {
    quadraticObjective_->deleteCols(numberToDelete,which);
    quadraticObjective_->deleteRows(numberToDelete,which);
  }
}

// Load up quadratic objective 
void 
ClpQuadraticObjective::loadQuadraticObjective(const int numberColumns, const CoinBigIndex * start,
			      const int * column, const double * element)
{
  delete quadraticObjective_;
  quadraticObjective_ = new CoinPackedMatrix(true,numberColumns,numberColumns,
					     start[numberColumns],element,column,start,NULL);
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
