// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpQuadraticObjective.hpp"
// at end to get min/max!
#include "ClpGubMatrix.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpGubMatrix::ClpGubMatrix () 
  : ClpPackedMatrix(),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberSets_(0)
{
  setType(11);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpGubMatrix::ClpGubMatrix (const ClpGubMatrix & rhs) 
: ClpPackedMatrix(rhs)
{  
  numberSets_ = rhs.numberSets_;
  start_ = ClpCopyOfArray(rhs.start_,numberSets_);
  end_ = ClpCopyOfArray(rhs.end_,numberSets_);
  lower_ = ClpCopyOfArray(rhs.lower_,numberSets_);
  upper_ = ClpCopyOfArray(rhs.upper_,numberSets_);
  status_ = ClpCopyOfArray(rhs.status_,numberSets_);
  int numberColumns = getNumCols();
  backward_ = ClpCopyOfArray(rhs.backward_,numberColumns);
  keyVariable_ = ClpCopyOfArray(rhs.keyVariable_,numberSets_);
  next_ = ClpCopyOfArray(rhs.next_,numberColumns+numberSets_);
}

//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpGubMatrix::ClpGubMatrix (CoinPackedMatrix * rhs) 
  : ClpPackedMatrix(rhs),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberSets_(0)
{  
  setType(11);
}

/* This takes over ownership (for space reasons) and is the
   real constructor*/
ClpGubMatrix::ClpGubMatrix(ClpPackedMatrix * matrix, int numberSets,
			   const int * start, const int * end,
			   const double * lower, const double * upper,
			   const unsigned char * status)
  : ClpPackedMatrix(matrix->matrix())
{
  numberSets_ = numberSets;
  start_ = ClpCopyOfArray(start,numberSets_);
  end_ = ClpCopyOfArray(end,numberSets_);
  lower_ = ClpCopyOfArray(lower,numberSets_);
  upper_ = ClpCopyOfArray(upper,numberSets_);
  // Check valid and ordered
  int last=-1;
  int numberColumns = matrix_->getNumCols();
  backward_ = new int[numberColumns];
  keyVariable_ = new int[numberSets_];
  // signal to need new ordering
  delete [] next_;
  next_ = NULL;
  for (int iColumn=0;iColumn<numberColumns;iColumn++) 
    backward_[iColumn]=-1;

  for (int iSet=0;iSet<numberSets_;iSet++) {
    if (start_[iSet]<0||start_[iSet]>=numberColumns)
      throw CoinError("Index out of range","constructor","ClpGubMatrix");
    if (end_[iSet]<0||end_[iSet]>numberColumns)
      throw CoinError("Index out of range","constructor","ClpGubMatrix");
    if (end_[iSet]<=start_[iSet])
      throw CoinError("Empty or negative set","constructor","ClpGubMatrix");
    if (start_[iSet]<last)
      throw CoinError("overlapping or non-monotonic sets","constructor","ClpGubMatrix");
    last=end_[iSet];
    int j;
    for (j=start_[iSet];j<end_[iSet];j++)
      backward_[j]=iSet;
  }
  if (status) {
    status_ = ClpCopyOfArray(status,numberSets_);
  } else {
    status_= new unsigned char [numberSets_];
    int i;
    for (i=0;i<numberSets_;i++) {
      // make slack key
      setStatus(i,ClpSimplex::basic);
    }
  }
}

ClpGubMatrix::ClpGubMatrix (const CoinPackedMatrix & rhs) 
  : ClpPackedMatrix(rhs),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberSets_(0)
{  
  setType(11);
  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpGubMatrix::~ClpGubMatrix ()
{
  delete [] start_;
  delete [] end_;
  delete [] lower_;
  delete [] upper_;
  delete [] status_;
  delete [] backward_;
  delete [] keyVariable_;
  delete [] next_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpGubMatrix &
ClpGubMatrix::operator=(const ClpGubMatrix& rhs)
{
  if (this != &rhs) {
    ClpPackedMatrix::operator=(rhs);
    delete [] start_;
    delete [] end_;
    delete [] lower_;
    delete [] upper_;
    delete [] status_;
    delete [] backward_;
    delete [] keyVariable_;
    delete [] next_;
    start_ = ClpCopyOfArray(rhs.start_,numberSets_);
    end_ = ClpCopyOfArray(rhs.end_,numberSets_);
    lower_ = ClpCopyOfArray(rhs.lower_,numberSets_);
    upper_ = ClpCopyOfArray(rhs.upper_,numberSets_);
    status_ = ClpCopyOfArray(rhs.status_,numberSets_);
    int numberColumns = getNumCols();
    backward_ = ClpCopyOfArray(rhs.backward_,numberColumns);
    keyVariable_ = ClpCopyOfArray(rhs.keyVariable_,numberSets_);
    next_ = ClpCopyOfArray(rhs.next_,numberColumns+numberSets_);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase * ClpGubMatrix::clone() const
{
  return new ClpGubMatrix(*this);
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase * 
ClpGubMatrix::subsetClone (int numberRows, const int * whichRows,
			   int numberColumns, 
			   const int * whichColumns) const 
{
  return new ClpGubMatrix(*this, numberRows, whichRows,
				   numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpGubMatrix::ClpGubMatrix (
			    const ClpGubMatrix & rhs,
			    int numberRows, const int * whichRows,
			    int numberColumns, const int * whichColumns)
  : ClpPackedMatrix(rhs, numberRows, whichRows, numberColumns, whichColumns)
{
  // Assuming no gub rows deleted
  // We also assume all sets in same order
  // Get array with backward pointers
  int numberColumnsOld = rhs.matrix_->getNumCols();
  int * array = new int [ numberColumnsOld];
  int i;
  for (i=0;i<numberColumnsOld;i++)
    array[i]=-1;
  for (int iSet=0;iSet<numberSets_;iSet++) {
    for (int j=start_[iSet];j<end_[iSet];j++)
      array[j]=iSet;
  }
  numberSets_=-1;
  int lastSet=-1;
  bool inSet=false;
  for (i=0;i<numberColumns;i++) {
    int iColumn = whichColumns[i];
    int iSet=array[iColumn];
    if (iSet<0) {
      inSet=false;
    } else {
      if (!inSet) {
	// start of new set but check okay
	if (iSet<=lastSet)
	  throw CoinError("overlapping or non-monotonic sets","subset constructor","ClpGubMatrix");
	lastSet = iSet;
	numberSets_++;
	start_[numberSets_]=i;
	end_[numberSets_]=i+1;
	lower_[numberSets_]=lower_[iSet];
	upper_[numberSets_]=upper_[iSet];
	inSet=true;
      } else {
	if (iSet<lastSet) {
	  throw CoinError("overlapping or non-monotonic sets","subset constructor","ClpGubMatrix");
	} else if (iSet==lastSet) {
	  end_[numberSets_]=i+1;
	} else {
	  // new set
	  lastSet = iSet;
	  numberSets_++;
	  start_[numberSets_]=i;
	  end_[numberSets_]=i+1;
	  lower_[numberSets_]=lower_[iSet];
	  upper_[numberSets_]=upper_[iSet];
	}
      }
    }
  }
  numberSets_++; // adjust
}
ClpGubMatrix::ClpGubMatrix (
		       const CoinPackedMatrix & rhs,
		       int numberRows, const int * whichRows,
		       int numberColumns, const int * whichColumns)
  : ClpPackedMatrix(rhs, numberRows, whichRows, numberColumns, whichColumns),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberSets_(0)
{
  setType(11);
}
#if 0
/* Returns a new matrix in reverse order without gaps */
ClpPackedMatrix * 
ClpGubMatrix::reverseOrderedCopy() const
{
  ClpPackedMatrix * copy = new ClpPackedMatrix();
  copy->matrix_= new CoinPackedMatrix();
  copy->matrix_->reverseOrderedCopyOf(*matrix_);
  copy->matrix_->removeGaps();
  // Row copy does not need stuff???
  return new ClpGubMatrix(copy,numberSets_,NULL,NULL,NULL,NULL);
}
#endif
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpGubMatrix::transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  // Do packed part
  ClpPackedMatrix::transposeTimes(model, scalar, rowArray, y, columnArray);
  if (numberSets_) {
    abort();
  }
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpGubMatrix::transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  // Do packed part
  ClpPackedMatrix::transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  if (numberSets_) {
    abort();
  }
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y.
   Squashes small elements and knows about ClpSimplex */
void 
ClpGubMatrix::subsetTransposeTimes(const ClpSimplex * model,
			      const CoinIndexedVector * rowArray,
			      const CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  // Do packed part
  ClpPackedMatrix::subsetTransposeTimes(model, rowArray, y, columnArray);
  if (numberSets_) {
    abort();
  }
}
/* If element NULL returns number of elements in column part of basis,
   If not NULL fills in as well */
CoinBigIndex 
ClpGubMatrix::fillBasis(const ClpSimplex * model,
			   const int * whichColumn, 
			   int numberBasic,
			   int numberColumnBasic,
			   int * indexRowU, int * indexColumnU,
			   double * elementU) const 
{
  int i;
  int numberColumns = getNumCols();
  const int * columnLength = matrix_->getVectorLengths(); 
  int numberRows = getNumRows();
  assert (next_ ||!elementU) ;
  if (!next_ ) {
    // do ordering
    next_ = new int[numberColumns+numberSets_];
    char * mark = new char[numberColumns];
    memset(mark,0,numberColumns);
    for (int iColumn=0;iColumn<numberColumns;iColumn++) 
      next_[iColumn]=INT_MAX;
    for (i=0;i<numberColumnBasic;i++) 
      mark[whichColumn[i]]=1;
    for (i=0;i<numberSets_;i++) {
      if (getStatus(i)!=ClpSimplex::basic) {
	// slack not key - choose one with smallest length
	int j;
	int smallest=numberRows+1;
	int key=-1;
	for (j=start_[i];j<end_[i];j++) {
	  if (mark[j]&&columnLength[j]<smallest) {
	    key=j;
	    smallest=columnLength[j];
	  }
	}
	if (key>=0) {
	  keyVariable_[i]=key;
	  int lastMarker = -(key+1);
	  next_[key]=lastMarker;
	  int last = key;
	  int j;
	  for (j=start_[i];j<end_[i];j++) {
	    if (mark[i]&&j!=key) {
	      next_[last]=j;
	      next_[j]=lastMarker;
	    }
	  }
	} else {
	  // nothing basic - make slack key
	  //((ClpGubMatrix *)this)->setStatus(i,ClpSimplex::basic);
	  // fudge to avoid const problem
	  status_[i]=1;
	}
      }
      if (getStatus(i)==ClpSimplex::basic) {
	// slack key
	keyVariable_[i]=numberColumns+i;
	int lastMarker = -(i+numberColumns+1);
	next_[numberColumns+i]=lastMarker;
	int last = numberColumns+i;
	int j;
	for (j=start_[i];j<end_[i];j++) {
	  if (mark[j]) {
	    next_[last]=j;
	    next_[j]=lastMarker;
	  }
	}
      }
    }
    delete [] mark;
  }
  CoinBigIndex numberElements=0;
  int lastSet=-1;
  int key=-1;
  int keyLength=-1;
  double * work = new double[numberRows];
  CoinZeroN(work,numberRows);
  char * mark = new char[numberRows];
  CoinZeroN(mark,numberRows);
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * row = matrix_->getIndices();
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  if (elementU!=NULL) {
    // fill
    if (!rowScale) {
    } else {
      // scaling
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	int iSet = backward_[iColumn];
	int length = columnLength[iColumn];
	CoinBigIndex j;
	if (iSet<0||keyVariable_[iSet]>=numberColumns) {
	  double scale = columnScale[iColumn];
	  for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    double value = elementByColumn[j];
	    if (fabs(value)>1.0e-20) {
	      int iRow = row[j];
	      indexRowU[numberElements]=iRow;
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value*scale*rowScale[iRow];
	    }
	  }
	  numberBasic++;
	} else {
	  // in gub set
	  if (iColumn!=keyVariable_[iSet]) {
	    double scale = columnScale[iColumn];
	    // not key 
	    if (lastSet<iSet) {
	      // erase work
	      if (key>=0) {
		for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
		  int iRow=row[j];
		  work[iRow]=0.0;
		  mark[iRow]=0;
		}
	      }
	      key=keyVariable_[iSet];
	      lastSet=iSet;
	      keyLength = columnLength[key];
	      double scale = columnScale[key];
	      for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
		int iRow=row[j];
		work[iRow]=elementByColumn[j]*scale*rowScale[iRow];
		mark[iRow]=1;
	      }
	    }
	    for (j=columnStart[iColumn];j<columnStart[iColumn]+length;j++) {
	      int iRow = row[j];
	      double value=elementByColumn[j]*scale*rowScale[iRow];
	      if (mark[iRow]) {
		mark[iRow]=0;
		double keyValue = work[iRow];
		value -= keyValue;
	      }
	      if (fabs(value)>1.0e-20) {
		indexRowU[numberElements]=iRow;
		indexColumnU[numberElements]=numberBasic;
		elementU[numberElements++]=value*scale*rowScale[iRow];
	      }
	    }
	    for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
	      int iRow = row[j];
	      if (mark[iRow]) {
		double value = -work[iRow];
		if (fabs(value)>1.0e-20) {
		  indexRowU[numberElements]=iRow;
		  indexColumnU[numberElements]=numberBasic;
		  elementU[numberElements++]=value*scale*rowScale[iRow];
		}
	      } else {
		// just put back mark
		mark[iRow]=1;
	      }
	    }
	    numberBasic++;
	  }
	}
      }
    }
  } else {
    // just count 
    if (!rowScale) {
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	int iSet = backward_[iColumn];
	int length = columnLength[iColumn];
	if (iSet<0||keyVariable_[iSet]>=numberColumns) {
	  numberElements += length;
	} else {
	  // in gub set
	  if (iColumn!=keyVariable_[iSet]) {
	    CoinBigIndex j;
	    // not key 
	    if (lastSet<iSet) {
	      // erase work
	      if (key>=0) {
		for (j=columnStart[key];j<columnStart[key]+keyLength;j++)
		  work[row[j]]=0.0;
	      }
	      key=keyVariable_[iSet];
	      lastSet=iSet;
	      keyLength = columnLength[key];
	      for (j=columnStart[key];j<columnStart[key]+keyLength;j++)
		work[row[j]]=elementByColumn[j];
	    }
	    int extra=keyLength;
	    for (j=columnStart[iColumn];j<columnStart[iColumn]+length;j++) {
	      int iRow = row[j];
	      double keyValue = work[iRow];
	      double value=elementByColumn[j];
	      if (!keyValue) {
		if (fabs(value)>1.0e-20)
		  extra++;
	      } else {
		value -= keyValue;
		if (fabs(value)<=1.0e-20)
		  extra--;
	      }
	    }
	    numberElements+=extra;
	  }
	}
      }
    } else {
      // scaled
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	int iSet = backward_[iColumn];
	int length = columnLength[iColumn];
	if (iSet<0||keyVariable_[iSet]>=numberColumns) {
	  numberElements += length;
	} else {
	  // in gub set
	  if (iColumn!=keyVariable_[iSet]) {
	    CoinBigIndex j;
	    double scale = columnScale[iColumn];
	    // not key 
	    if (lastSet<iSet) {
	      // erase work
	      if (key>=0) {
		for (j=columnStart[key];j<columnStart[key]+keyLength;j++)
		  work[row[j]]=0.0;
	      }
	      key=keyVariable_[iSet];
	      lastSet=iSet;
	      keyLength = columnLength[key];
	      double scale = columnScale[key];
	      for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
		int iRow = row[j];
		work[iRow]=elementByColumn[j]*scale*rowScale[iRow];
	      }
	    }
	    int extra=keyLength;
	    for (j=columnStart[iColumn];j<columnStart[iColumn]+length;j++) {
	      int iRow = row[j];
	      double keyValue = work[iRow];
	      double value=elementByColumn[j]*scale*rowScale[iRow];
	      if (!keyValue) {
		if (fabs(value)>1.0e-20)
		  extra++;
	      } else {
		value -= keyValue;
		if (fabs(value)<=1.0e-20)
		  extra--;
	      }
	    }
	    numberElements+=extra;
	  }
	}
      }
    }
  }
  delete [] work;
  delete [] mark;
  return numberElements;
}
/* Unpacks a column into an CoinIndexedvector
 */
void 
ClpGubMatrix::unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn) const 
{
  // Do packed part
  ClpPackedMatrix::unpack(model,rowArray,iColumn);
  if (numberSets_) {
    abort();
  }
}
/* Unpacks a column into a CoinIndexedVector
** in packed format
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void 
ClpGubMatrix::unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int iColumn) const
{
  // Do packed part
  ClpPackedMatrix::unpackPacked(model,rowArray,iColumn);
  if (numberSets_) {
    abort();
  }
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void 
ClpGubMatrix::add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn, double multiplier) const 
{
  // Do packed part
  ClpPackedMatrix::add(model,rowArray,iColumn,multiplier);
  if (numberSets_) {
    abort();
  }
}
// Partial pricing 
void 
ClpGubMatrix::partialPricing(ClpSimplex * model, int start, int end,
			      int & bestSequence, int & numberWanted)
{
  // Do packed part
  ClpPackedMatrix::partialPricing(model,start,end,bestSequence,numberWanted);
  if (numberSets_) {
    abort();
  }
}
/* expands an updated column to allow for extra rows which the main
   solver does not know about and returns number added.  If the arrays are NULL 
   then returns number of extra entries needed.
   
   This will normally be a no-op - it is in for GUB!
*/
int 
ClpGubMatrix::extendUpdated(CoinIndexedVector * update, double * lower,
			     double * solution, double * upper)
{
  assert(!update);
  return getNumRows();
}


