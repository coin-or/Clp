// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpModel.hpp"
#include "ClpLinearObjective.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpLinearObjective::ClpLinearObjective () 
: ClpObjective()
{
  type_=1;
  objective_=NULL;
  numberColumns_=0;
}

//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
ClpLinearObjective::ClpLinearObjective (const double * objective , 
					int numberColumns) 
  : ClpObjective()
{
  type_=1;
  numberColumns_=numberColumns;
  if (objective) {
    objective_ = new double [numberColumns_];
    memcpy(objective_,objective,numberColumns_*sizeof(double));
  } else {
    objective_ = new double [numberColumns_];
    memset(objective_,0,numberColumns_*sizeof(double));
  }
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpLinearObjective::ClpLinearObjective (const ClpLinearObjective & rhs) 
: ClpObjective(rhs)
{  
  numberColumns_=rhs.numberColumns_;
  if (rhs.objective_) {
    objective_ = new double [numberColumns_];
    memcpy(objective_,rhs.objective_,numberColumns_*sizeof(double));
  } else {
    objective_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpLinearObjective::~ClpLinearObjective ()
{
  delete [] objective_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpLinearObjective &
ClpLinearObjective::operator=(const ClpLinearObjective& rhs)
{
  if (this != &rhs) {
    ClpObjective::operator=(rhs);
    numberColumns_=rhs.numberColumns_;
    if (rhs.objective_) {
      objective_ = new double [numberColumns_];
      memcpy(objective_,rhs.objective_,numberColumns_*sizeof(double));
    } else {
      objective_=NULL;
    }
  }
  return *this;
}

// Returns gradient
double *  
ClpLinearObjective::gradient(const double * solution, double & offset)
{
  offset=0.0;
  return objective_;
}
  
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpObjective * ClpLinearObjective::clone() const
{
  return new ClpLinearObjective(*this);
}
// Resize objective
void 
ClpLinearObjective::resize(int newNumberColumns)
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
    numberColumns_ = newNumberColumns;
  } 
  
}
// Delete columns in  objective
void 
ClpLinearObjective::deleteSome(int numberToDelete, const int * which) 
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
}

