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
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    sumOfRelaxedDualInfeasibilities_(0.0),
    sumOfRelaxedPrimalInfeasibilities_(0.0),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberDualInfeasibilities_(0),
    numberPrimalInfeasibilities_(0),
    numberSets_(0),
    firstGub_(0),
    lastGub_(0),
    gubType_(0)
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
  sumDualInfeasibilities_ = rhs. sumDualInfeasibilities_;
  sumPrimalInfeasibilities_ = rhs.sumPrimalInfeasibilities_;
  sumOfRelaxedDualInfeasibilities_ = rhs.sumOfRelaxedDualInfeasibilities_;
  sumOfRelaxedPrimalInfeasibilities_ = rhs.sumOfRelaxedPrimalInfeasibilities_;
  numberDualInfeasibilities_ = rhs.numberDualInfeasibilities_;
  numberPrimalInfeasibilities_ = rhs.numberPrimalInfeasibilities_;
  firstGub_ = rhs.firstGub_;
  lastGub_ = rhs.lastGub_;
  gubType_ = rhs.gubType_;
}

//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpGubMatrix::ClpGubMatrix (CoinPackedMatrix * rhs) 
  : ClpPackedMatrix(rhs),
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    sumOfRelaxedDualInfeasibilities_(0.0),
    sumOfRelaxedPrimalInfeasibilities_(0.0),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberDualInfeasibilities_(0),
    numberPrimalInfeasibilities_(0),
    numberSets_(0),
    firstGub_(0),
    lastGub_(0),
    gubType_(0)
{  
  setType(11);
}

/* This takes over ownership (for space reasons) and is the
   real constructor*/
ClpGubMatrix::ClpGubMatrix(ClpPackedMatrix * matrix, int numberSets,
			   const int * start, const int * end,
			   const double * lower, const double * upper,
			   const unsigned char * status)
  : ClpPackedMatrix(matrix->matrix()),
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    sumOfRelaxedDualInfeasibilities_(0.0),
    sumOfRelaxedPrimalInfeasibilities_(0.0),
    numberDualInfeasibilities_(0),
    numberPrimalInfeasibilities_(0)
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
  next_ = NULL;
  for (int iColumn=0;iColumn<numberColumns;iColumn++) 
    backward_[iColumn]=-1;

  int iSet;
  for (iSet=0;iSet<numberSets_;iSet++) {
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
  // Find type of gub
  firstGub_=numberColumns+1;
  lastGub_=-1;
  int i;
  for (i=0;i<numberColumns;i++) {
    if (backward_[i]>=0) {
      firstGub_ = min(firstGub_,i);
      lastGub_ = max(lastGub_,i);
    }
  }
  gubType_=0;
  for (i=firstGub_;i<lastGub_;i++) {
    if (backward_[i]<0) {
      gubType_=1;
      break;
    }
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
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    sumOfRelaxedDualInfeasibilities_(0.0),
    sumOfRelaxedPrimalInfeasibilities_(0.0),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    status_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberDualInfeasibilities_(0),
    numberPrimalInfeasibilities_(0),
    numberSets_(0),
    firstGub_(0),
    lastGub_(0),
    gubType_(0)
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
    sumDualInfeasibilities_ = rhs. sumDualInfeasibilities_;
    sumPrimalInfeasibilities_ = rhs.sumPrimalInfeasibilities_;
    sumOfRelaxedDualInfeasibilities_ = rhs.sumOfRelaxedDualInfeasibilities_;
    sumOfRelaxedPrimalInfeasibilities_ = rhs.sumOfRelaxedPrimalInfeasibilities_;
    numberDualInfeasibilities_ = rhs.numberDualInfeasibilities_;
    numberPrimalInfeasibilities_ = rhs.numberPrimalInfeasibilities_;
    firstGub_ = rhs.firstGub_;
    lastGub_ = rhs.lastGub_;
    gubType_ = rhs.gubType_;
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
  // Find type of gub
  firstGub_=numberColumns+1;
  lastGub_=-1;
  for (i=0;i<numberColumns;i++) {
    if (backward_[i]>=0) {
      firstGub_ = min(firstGub_,i);
      lastGub_ = max(lastGub_,i);
    }
  }
  gubType_=0;
  for (i=firstGub_;i<lastGub_;i++) {
    if (backward_[i]<0) {
      gubType_=1;
      break;
    }
  }

  // Make sure key is feasible if only key in set
}
ClpGubMatrix::ClpGubMatrix (
		       const CoinPackedMatrix & rhs,
		       int numberRows, const int * whichRows,
		       int numberColumns, const int * whichColumns)
  : ClpPackedMatrix(rhs, numberRows, whichRows, numberColumns, whichColumns),
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    sumOfRelaxedDualInfeasibilities_(0.0),
    sumOfRelaxedPrimalInfeasibilities_(0.0),
    start_(NULL),
    end_(NULL),
    lower_(NULL),
    upper_(NULL),
    backward_(NULL),
    keyVariable_(NULL),
    next_(NULL),
    numberDualInfeasibilities_(0),
    numberPrimalInfeasibilities_(0),
    numberSets_(0),
    firstGub_(0),
    lastGub_(0),
    gubType_(0)
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
ClpGubMatrix::fillBasis(ClpSimplex * model,
			   const int * whichColumn, 
			   int numberBasic,
			   int numberColumnBasic,
			   int * indexRowU, int * indexColumnU,
			   double * elementU) 
{
  int i;
  int numberColumns = getNumCols();
  const int * columnLength = matrix_->getVectorLengths(); 
  int numberRows = getNumRows();
  assert (next_ ||!elementU) ;
  if (!next_ ) {
    // do ordering
    assert (!effectiveRhs_);
    // create and do gub crash
    useEffectiveRhs(model,false);
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
      // no scaling
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	int iSet = backward_[iColumn];
	int length = columnLength[iColumn];
	CoinBigIndex j;
	if (iSet<0||keyVariable_[iSet]>=numberColumns) {
	  for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    double value = elementByColumn[j];
	    if (fabs(value)>1.0e-20) {
	      int iRow = row[j];
	      indexRowU[numberElements]=iRow;
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value;
	    }
	  }
	  numberBasic++;
	} else {
	  // in gub set
	  if (iColumn!=keyVariable_[iSet]) {
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
	      for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
		int iRow=row[j];
		work[iRow]=elementByColumn[j];
		mark[iRow]=1;
	      }
	    }
	    for (j=columnStart[iColumn];j<columnStart[iColumn]+length;j++) {
	      int iRow = row[j];
	      double value=elementByColumn[j];
	      if (mark[iRow]) {
		mark[iRow]=0;
		double keyValue = work[iRow];
		value -= keyValue;
	      }
	      if (fabs(value)>1.0e-20) {
		indexRowU[numberElements]=iRow;
		indexColumnU[numberElements]=numberBasic;
		elementU[numberElements++]=value;
	      }
	    }
	    for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
	      int iRow = row[j];
	      if (mark[iRow]) {
		double value = -work[iRow];
		if (fabs(value)>1.0e-20) {
		  indexRowU[numberElements]=iRow;
		  indexColumnU[numberElements]=numberBasic;
		  elementU[numberElements++]=value;
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
	    int iRow = row[j];
	    double value = elementByColumn[j]*scale*rowScale[iRow];
	    if (fabs(value)>1.0e-20) {
	      indexRowU[numberElements]=iRow;
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value;
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
		elementU[numberElements++]=value;
	      }
	    }
	    for (j=columnStart[key];j<columnStart[key]+keyLength;j++) {
	      int iRow = row[j];
	      if (mark[iRow]) {
		double value = -work[iRow];
		if (fabs(value)>1.0e-20) {
		  indexRowU[numberElements]=iRow;
		  indexColumnU[numberElements]=numberBasic;
		  elementU[numberElements++]=value;
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
  int iSet = backward_[iColumn];
  if (iSet>=0) {
    int iBasic = keyVariable_[iSet];
    if (iBasic <model->numberColumns()) {
      add(model,rowArray,iBasic,-1.0);
    }
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
  int iSet = backward_[iColumn];
  if (iSet>=0) {
    // column in order
    abort();
    int iBasic = keyVariable_[iSet];
    if (iBasic <model->numberColumns()) {
      add(model,rowArray,iBasic,-1.0);
    }
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
  if (numberSets_) 
    assert(!gubType_);
  // Do packed part before gub
  ClpPackedMatrix::partialPricing(model,start,firstGub_,bestSequence,numberWanted);
  if (numberWanted) {
    // do gub
    const double * element =matrix_->getElements();
    const int * row = matrix_->getIndices();
    const CoinBigIndex * startColumn = matrix_->getVectorStarts();
    const int * length = matrix_->getVectorLengths();
    const double * rowScale = model->rowScale();
    const double * columnScale = model->columnScale();
    int iSequence;
    CoinBigIndex j;
    double tolerance=model->currentDualTolerance();
    double * reducedCost = model->djRegion();
    const double * duals = model->dualRowSolution();
    const double * cost = model->costRegion();
    double bestDj;
    int numberColumns = model->numberColumns();
    if (bestSequence>=0)
      bestDj = fabs(reducedCost[bestSequence]);
    else
      bestDj=tolerance;
    int sequenceOut = model->sequenceOut();
    int saveSequence = bestSequence;
    start = max (start,firstGub_);
    end = min(lastGub_,end);
    int iSet = -1;
    double djMod=0.0;
    if (rowScale) {
      // scaled
      for (iSequence=start;iSequence<end;iSequence++) {
	if (backward_[iSequence]!=iSet) {
	  // get pi on gub row
	  iSet = backward_[iSequence];
	  int iBasic = keyVariable_[iSet];
	  if (iBasic>numberColumns) {
	    djMod = 0.0;
	  } else {
	    // get dj without 
	    assert (model->getStatus(iBasic)==ClpSimplex::basic);
	    djMod=0.0;
	    // scaled
	    for (j=startColumn[iBasic];
		 j<startColumn[iBasic]+length[iBasic];j++) {
	      int jRow=row[j];
	      djMod -= duals[jRow]*element[j]*rowScale[jRow];
	    }
	    // allow for scaling
	    djMod +=  cost[iBasic]/columnScale[iBasic];
	  }
	}
	if (iSequence!=sequenceOut) {
	  double value;
	  ClpSimplex::Status status = model->getStatus(iSequence);
	  
	  switch(status) {
	    
	  case ClpSimplex::basic:
	  case ClpSimplex::isFixed:
	    break;
	  case ClpSimplex::isFree:
	  case ClpSimplex::superBasic:
	    value=djMod;
	    // scaled
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j]*rowScale[jRow];
	    }
	    value = fabs(cost[iSequence] +value*columnScale[iSequence]);
	    if (value>FREE_ACCEPT*tolerance) {
	      numberWanted--;
	      // we are going to bias towards free (but only if reasonable)
	      value *= FREE_BIAS;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  case ClpSimplex::atUpperBound:
	    value=djMod;
	    // scaled
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j]*rowScale[jRow];
	    }
	    value = cost[iSequence] +value*columnScale[iSequence];
	    if (value>tolerance) {
	      numberWanted--;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  case ClpSimplex::atLowerBound:
	    value=djMod;
	    // scaled
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j]*rowScale[jRow];
	    }
	    value = -(cost[iSequence] +value*columnScale[iSequence]);
	    if (value>tolerance) {
	      numberWanted--;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  }
	}
	if (!numberWanted)
	  break;
      }
      if (bestSequence!=saveSequence) {
	// recompute dj
	double value=0.0;
	// scaled
	for (j=startColumn[bestSequence];
	     j<startColumn[bestSequence]+length[bestSequence];j++) {
	  int jRow=row[j];
	  value -= duals[jRow]*element[j]*rowScale[jRow];
	}
	reducedCost[bestSequence] = cost[bestSequence] +value*columnScale[bestSequence];
      }
    } else {
      // not scaled
      for (iSequence=start;iSequence<end;iSequence++) {
	if (backward_[iSequence]!=iSet) {
	  // get pi on gub row
	  iSet = backward_[iSequence];
	  int iBasic = keyVariable_[iSet];
	  if (iBasic>numberColumns) {
	    djMod = 0.0;
	  } else {
	    // get dj without 
	    assert (model->getStatus(iBasic)==ClpSimplex::basic);
	    djMod=0.0;
	    // scaled
	    for (j=startColumn[iBasic];
		 j<startColumn[iBasic]+length[iBasic];j++) {
	      int jRow=row[j];
	      djMod -= duals[jRow]*element[j];
	    }
	    djMod += cost[iBasic];
	  }
	}
	if (iSequence!=sequenceOut) {
	  double value;
	  ClpSimplex::Status status = model->getStatus(iSequence);
	  
	  switch(status) {
	    
	  case ClpSimplex::basic:
	  case ClpSimplex::isFixed:
	    break;
	  case ClpSimplex::isFree:
	  case ClpSimplex::superBasic:
	    value=cost[iSequence]-djMod;
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j];
	    }
	    value = fabs(value);
	    if (value>FREE_ACCEPT*tolerance) {
	      numberWanted--;
	      // we are going to bias towards free (but only if reasonable)
	      value *= FREE_BIAS;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  case ClpSimplex::atUpperBound:
	    value=cost[iSequence]-djMod;
	    // scaled
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j];
	    }
	    if (value>tolerance) {
	      numberWanted--;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  case ClpSimplex::atLowerBound:
	    value=cost[iSequence]-djMod;
	    for (j=startColumn[iSequence];
		 j<startColumn[iSequence]+length[iSequence];j++) {
	      int jRow=row[j];
	      value -= duals[jRow]*element[j];
	    }
	    value = -value;
	    if (value>tolerance) {
	      numberWanted--;
	      if (value>bestDj) {
		// check flagged variable and correct dj
		if (!model->flagged(iSequence)) {
		  bestDj=value;
		  bestSequence = iSequence;
		} else {
		  // just to make sure we don't exit before got something
		  numberWanted++;
		}
	      }
	    }
	    break;
	  }
	}
	if (!numberWanted)
	  break;
      }
      if (bestSequence!=saveSequence) {
	// recompute dj
	double value=cost[bestSequence];
	for (j=startColumn[bestSequence];
	     j<startColumn[bestSequence]+length[bestSequence];j++) {
	  int jRow=row[j];
	  value -= duals[jRow]*element[j];
	}
	reducedCost[bestSequence] = value;
      }
    }
    if (numberWanted) {
      // Do packed part after gub
      ClpPackedMatrix::partialPricing(model,lastGub_,end,bestSequence,numberWanted);
    }
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
/*
     utility primal function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
void 
ClpGubMatrix::primalExpanded(ClpSimplex * model,int mode)
{
  int numberColumns = model->numberColumns();
  switch (mode) {
    // If key variable then slot in gub rhs so will get correct contribution
  case 0:
    {
      int i;
      double * solution = model->solutionRegion();
      ClpSimplex::Status iStatus;
      for (i=0;i<numberSets_;i++) {
	int iColumn = keyVariable_[i];
	if (iColumn<numberColumns) {
	  // key is structural - where is slack
	  iStatus = getStatus(i);
	  assert (iStatus!=ClpSimplex::basic);
	  if (iStatus==ClpSimplex::atLowerBound)
	    solution[iColumn]=lower_[i];
	  else
	    solution[iColumn]=upper_[i];
	}
      }
    }
    break;
    // Compute values of key variables
  case 1:
    {
      int i;
      double * solution = model->solutionRegion();
      ClpSimplex::Status iStatus;
      //const int * columnLength = matrix_->getVectorLengths(); 
      //const CoinBigIndex * columnStart = matrix_->getVectorStarts();
      //const int * row = matrix_->getIndices();
      //const double * elementByColumn = matrix_->getElements();
      //int * pivotVariable = model->pivotVariable();
      sumPrimalInfeasibilities_=0.0;
      numberPrimalInfeasibilities_=0;
      double primalTolerance = model->primalTolerance();
      double relaxedTolerance=primalTolerance;
      // we can't really trust infeasibilities if there is primal error
      double error = min(1.0e-3,model->largestPrimalError());
      // allow tolerance at least slightly bigger than standard
      relaxedTolerance = relaxedTolerance +  error;
      // but we will be using difference
      relaxedTolerance -= primalTolerance;
      sumOfRelaxedPrimalInfeasibilities_ = 0.0;
      for (i=0;i<numberSets_;i++) {
	int kColumn = keyVariable_[i];
	int iColumn =next_[kColumn];
	// sum all non-key variables
	double value=0.0;
	while(iColumn>=0) {
	  value+=solution[iColumn];
	  iColumn=next_[iColumn];
	}
	if (kColumn<numberColumns) {
	  // feasibility will be done later
	  solution[kColumn] -= value;
	} else {
	  // slack is key
	  iStatus = getStatus(i);
	  assert (iStatus==ClpSimplex::basic);
	  double infeasibility=0.0;
	  if (value>upper_[i]+primalTolerance) {
	    infeasibility=value-upper_[i]-primalTolerance;
	    setAbove(i);
	  } else if (value<lower_[i]-primalTolerance) {
	    infeasibility=lower_[i]-value-primalTolerance ;
	    setBelow(i);
	  } else {
	    setFeasible(i);
	  }
	  if (infeasibility>0.0) {
	    sumPrimalInfeasibilities_ += infeasibility;
	    if (infeasibility>relaxedTolerance) 
	      sumOfRelaxedPrimalInfeasibilities_ += infeasibility;
	    numberPrimalInfeasibilities_ ++;
	  }
	}
      }
    }
    break;
    // Report on infeasibilities of key variables
  case 2:
    {
      model->setSumPrimalInfeasibilities(model->sumPrimalInfeasibilities()+
					 sumPrimalInfeasibilities_);
      model->setNumberPrimalInfeasibilities(model->numberPrimalInfeasibilities()+
					 numberPrimalInfeasibilities_);
      model->setSumOfRelaxedPrimalInfeasibilities(model->sumOfRelaxedPrimalInfeasibilities()+
					 sumOfRelaxedPrimalInfeasibilities_);
    }
    break;
  }
}
/*
     utility dual function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
void 
ClpGubMatrix::dualExpanded(ClpSimplex * model,
			    CoinIndexedVector * array,
			    double * other,int mode)
{
  switch (mode) {
    // modify costs before transposeUpdate
  case 0:
    {
      int i;
      double * cost = model->costRegion();
      ClpSimplex::Status iStatus;
      // not dual values yet
      assert (!other);
      double * work = array->denseVector();
      double infeasibilityCost = model->infeasibilityCost();
      int * pivotVariable = model->pivotVariable();
      int numberRows = model->numberRows();
      int numberColumns = model->numberColumns();
      for (i=0;i<numberRows;i++) {
	int iPivot = pivotVariable[i];
	if (iPivot<numberColumns) {
	  int iSet = backward_[iPivot];
	  if (iSet>=0) {
	    int kColumn = keyVariable_[iSet];
	    double costValue;
	    if (kColumn<numberColumns) {
	      // structural has cost
	      costValue = cost[kColumn];
	    } else {
	      // slack is key
	      iStatus = getStatus(iSet);
	      assert (iStatus==ClpSimplex::basic);
	      costValue=weight(iSet)*infeasibilityCost;
	    }
	    array->add(i,work[i]-costValue);
	  }
	}
      }
    }
    break;
    // create duals for key variables (without check on dual infeasible)
  case 1:
    {
      // If key slack then dual 0.0
      // dj for key is zero so that defines dual on set
      int i;
      double * dj = model->djRegion();
      double * dual = model->dualRowSolution();
      double * cost = model->costRegion();
      const int * columnLength = matrix_->getVectorLengths(); 
      const CoinBigIndex * columnStart = matrix_->getVectorStarts();
      const int * row = matrix_->getIndices();
      const double * elementByColumn = matrix_->getElements();
      int numberColumns = model->numberColumns();
      for (i=0;i<numberSets_;i++) {
	int kColumn = keyVariable_[i];
	// If slack key we need not do anything
	if (kColumn<numberColumns) {
	  // dj without set
	  double value = cost[kColumn];
	  for (CoinBigIndex j=columnStart[kColumn];
	       j<columnStart[kColumn]+columnLength[kColumn];j++) {
	    int iRow = row[j];
	    value -= dual[iRow]*elementByColumn[j];
	  }
	  // Now subtract out from all 
	  dj[kColumn] -= value;
	  int iColumn =next_[kColumn];
	  // modify all non-key variables
	  while(iColumn>=0) {
	  dj[iColumn]-=value;
	  iColumn=next_[iColumn];
	  }
	}
      }
    }
    // as 1 but check slacks
  case 2:
    {
      // If key slack then dual 0.0
      // If not then slack could be dual infeasible
      // dj for key is zero so that defines dual on set
      int i;
      double * dj = model->djRegion();
      double * dual = model->dualRowSolution();
      double * cost = model->costRegion();
      ClpSimplex::Status iStatus;
      const int * columnLength = matrix_->getVectorLengths(); 
      const CoinBigIndex * columnStart = matrix_->getVectorStarts();
      const int * row = matrix_->getIndices();
      const double * elementByColumn = matrix_->getElements();
      int numberColumns = model->numberColumns();
      sumDualInfeasibilities_=0.0;
      numberDualInfeasibilities_=0;
      double dualTolerance = model->dualTolerance();
      double relaxedTolerance=dualTolerance;
      // we can't really trust infeasibilities if there is dual error
      double error = min(1.0e-3,model->largestDualError());
      // allow tolerance at least slightly bigger than standard
      relaxedTolerance = relaxedTolerance +  error;
      // but we will be using difference
      relaxedTolerance -= dualTolerance;
      sumOfRelaxedDualInfeasibilities_ = 0.0;
      for (i=0;i<numberSets_;i++) {
	int kColumn = keyVariable_[i];
	// If slack key we need not do anything
	if (kColumn<numberColumns) {
	  // dj without set
	  double value = cost[kColumn];
	  for (CoinBigIndex j=columnStart[kColumn];
	       j<columnStart[kColumn]+columnLength[kColumn];j++) {
	    int iRow = row[j];
	    value -= dual[iRow]*elementByColumn[j];
	  }
	  // Now subtract out from all 
	  dj[kColumn] -= value;
	  int iColumn =next_[kColumn];
	  // modify all non-key variables
	  while(iColumn>=0) {
	  dj[iColumn]-=value;
	  iColumn=next_[iColumn];
	  }
	  // check slack
	  iStatus = getStatus(i);
	  assert (iStatus!=ClpSimplex::basic);
	  double infeasibility=0.0;
	  // dj of slack is -(-1.0)value
	  if (iStatus==ClpSimplex::atLowerBound) {
	    if (value>-dualTolerance) 
	      infeasibility=-value-dualTolerance;
	  } else {
	    // at upper bound
	    if (value>dualTolerance) 
	      infeasibility=value-dualTolerance;
	  }
	  if (infeasibility>0.0) {
	    sumDualInfeasibilities_ += infeasibility;
	    if (infeasibility>relaxedTolerance) 
	      sumOfRelaxedDualInfeasibilities_ += infeasibility;
	    numberDualInfeasibilities_ ++;
	  }
	}
      }
    }
    // Report on infeasibilities of key variables
  case 3:
    {
      model->setSumDualInfeasibilities(model->sumDualInfeasibilities()+
					 sumDualInfeasibilities_);
      model->setNumberDualInfeasibilities(model->numberDualInfeasibilities()+
					 numberDualInfeasibilities_);
      model->setSumOfRelaxedDualInfeasibilities(model->sumOfRelaxedDualInfeasibilities()+
					 sumOfRelaxedDualInfeasibilities_);
    }
    break;
  }
}
/*
     general utility function for dealing with dynamic constraints
     mode=n see ClpGubMatrix.hpp for definition
     Remember to update here when settled down
*/
int
ClpGubMatrix::generalExpanded(ClpSimplex * model,int mode,int &number)
{
  int returnCode=0;
  int numberColumns = model->numberColumns();
  switch (mode) {
    // Fill in pivotVariable but not for key variables
  case 0:
    {
      int i;
      int numberBasic=number;
      // Use different array so can build from true pivotVariable_
      //int * pivotVariable = model->pivotVariable();
      int * pivotVariable = model->rowArray(0)->getIndices();
      for (i=0;i<numberColumns;i++) {
	if (model->getColumnStatus(i) == ClpSimplex::basic) {
	  int iSet = backward_[i];
	  if (iSet<0||i!=keyVariable_[iSet])
	    pivotVariable[numberBasic++]=i;
	}
      }
      number = numberBasic;
    }
    break;
    // Make all key variables basic
  case 1:
    {
      int i;
      for (i=0;i<numberSets_;i++) {
	int iColumn = keyVariable_[i];
	if (iColumn<numberColumns)
	  model->setColumnStatus(iColumn,ClpSimplex::basic);
      }
    }
    break;
  }
  return returnCode;
}
// Sets up an effective RHS and does gub crash if needed
void 
ClpGubMatrix::useEffectiveRhs(ClpSimplex * model, bool cheapest)
{
  // Do basis - cheapest or slack if feasible (unless cheapest set)
  int longestSet=0;
  int iSet;
  for (iSet=0;iSet<numberSets_;iSet++) 
    longestSet = max(longestSet,end_[iSet]-start_[iSet]);
    
  double * upper = new double[longestSet+1];
  double * cost = new double[longestSet+1];
  double * lower = new double[longestSet+1];
  double * solution = new double[longestSet+1];
  assert (!next_);
  int numberColumns = getNumCols();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * columnLower = model->lowerRegion();
  const double * columnUpper = model->upperRegion();
  double * columnSolution = model->solutionRegion();
  const double * objective = model->costRegion();
  int numberRows = getNumRows();
  next_ = new int[numberColumns+numberSets_];
  for (int iColumn=0;iColumn<numberColumns;iColumn++) 
    next_[iColumn]=INT_MAX;
  double tolerance = model->primalTolerance();
  for (iSet=0;iSet<numberSets_;iSet++) {
    int j;
    int numberBasic=0;
    int iBasic=-1;
    int iStart = start_[iSet];
    int iEnd=end_[iSet];
    // find one with smallest length
    int smallest=numberRows+1;
    double value=0.0;
    for (j=iStart;j<iEnd;j++) {
      if (model->getStatus(j)== ClpSimplex::basic) {
	if (columnLength[j]<smallest) {
	  smallest=columnLength[j];
	  iBasic=j;
	}
	numberBasic++;
      }
      value += columnSolution[j];
    }
    bool done=false;
    if (numberBasic>1||(numberBasic==1&&getStatus(iSet)==ClpSimplex::basic)) {
      if (getStatus(iSet)==ClpSimplex::basic) 
	iBasic = iSet+numberColumns;// slack key - use
      done=true;
    } else if (numberBasic==1) {
      // see if can be key
      double thisSolution = columnSolution[iBasic];
      if (thisSolution>columnUpper[iBasic]) {
	value -= thisSolution-columnUpper[iBasic];
	thisSolution = columnUpper[iBasic];
	solution[iBasic]=thisSolution;
      }
      if (thisSolution<columnLower[iBasic]) {
	value -= thisSolution-columnLower[iBasic];
	thisSolution = columnLower[iBasic];
	solution[iBasic]=thisSolution;
      }
      // try setting slack to a bound
      assert (upper_[iSet]<1.0e20||lower_[iSet]>-1.0e20);
      double cost1 = COIN_DBL_MAX;
      int whichBound=-1;
      if (upper_[iSet]<1.0e20) {
	// try slack at ub
	double newBasic = thisSolution +upper_[iSet]-value;
	if (newBasic>=columnLower[iBasic]-tolerance&&
	    newBasic<=columnUpper[iBasic]+tolerance) {
	  // can go
	  whichBound=1;
	  cost1 = newBasic*objective[iBasic];
	  // But if exact then may be good solution
	  if (fabs(upper_[iSet]-value)<tolerance)
	    cost1=-COIN_DBL_MAX;
	}
      }
      if (lower_[iSet]>-1.0e20) {
	// try slack at lb
	double newBasic = thisSolution +lower_[iSet]-value;
	if (newBasic>=columnLower[iBasic]-tolerance&&
	    newBasic<=columnUpper[iBasic]+tolerance) {
	  // can go but is it cheaper
	  double cost2 = newBasic*objective[iBasic];
	  // But if exact then may be good solution
	  if (fabs(lower_[iSet]-value)<tolerance)
	    cost2=-COIN_DBL_MAX;
	  if (cost2<cost1)
	    whichBound=0;
	}
      }
      if (whichBound!=-1) {
	// key
	done=true;
	if (whichBound) {
	  // slack to upper
	  columnSolution[iBasic]=thisSolution + upper_[iSet]-value;
	  setStatus(iSet,ClpSimplex::atUpperBound);
	} else {
	  // slack to lower
	  columnSolution[iBasic]=thisSolution + lower_[iSet]-value;
	  setStatus(iSet,ClpSimplex::atLowerBound);
	}
      }
    }
    if (!done) {
      if (!cheapest) {
	// see if slack can be key
	if (value>=lower_[iSet]-tolerance&&value<=upper_[iSet]+tolerance) {
	  done=true;
	  setStatus(iSet,ClpSimplex::basic);
	  iBasic=iSet+numberColumns;
	}
      }
      if (!done) {
	// find cheapest
	int numberInSet = iEnd-iStart;
	CoinMemcpyN(columnLower+iStart,numberInSet,lower);
	CoinMemcpyN(columnUpper+iStart,numberInSet,upper);
	CoinMemcpyN(columnSolution+iStart,numberInSet,solution);
	// and slack
	iBasic=numberInSet;
	solution[iBasic]=-value;
	lower[iBasic]=-upper_[iSet];
	upper[iBasic]=-lower_[iSet];
	int kphase;
	if (value>=lower_[iSet]-tolerance&&value<=upper_[iSet]+tolerance) {
	  // feasible
	  kphase=1;
	  cost[iBasic]=0.0;
	  CoinMemcpyN(objective+iStart,numberInSet,cost);
	} else {
	  // infeasible
	  kphase=0;
	  // remember bounds are flipped so opposite to natural
	  if (value<lower_[iSet]-tolerance)
	    cost[iBasic]=1.0;
	  else
	    cost[iBasic]=-1.0;
	  CoinZeroN(cost,numberInSet);
	}
	double dualTolerance =model->dualTolerance();
	for (int iphase =kphase;iphase<2;iphase++) {
	  if (iphase) {
	    cost[iBasic]=0.0;
	    CoinMemcpyN(objective+iStart,numberInSet,cost);
	  }
	  // now do one row lp
	  bool improve=true;
	  while (improve) {
	    improve=false;
	    double dual = cost[iBasic];
	    int chosen =-1;
	    double best=dualTolerance;
	    int way=0;
	    for (int i=0;i<=numberInSet;i++) {
	      double dj = cost[i]-dual;
	      double improvement =0.0;
	      double distance=0.0;
	      if (iphase||i<numberInSet)
		assert (solution[i]>=lower[i]&&solution[i]<=upper[i]);
	      if (dj>dualTolerance)
		improvement = dj*(solution[i]-lower[i]);
	      else if (dj<-dualTolerance)
		improvement = dj*(solution[i]-upper[i]);
	      if (improvement>best) {
		best=improvement;
		chosen=i;
		if (dj<0.0) {
		  way = 1;
		  distance = upper[i]-solution[i];
		} else {
		  way = -1;
		  distance = solution[i]-lower[i];
		}
	      }
	    }
	    if (chosen>=0) {
	      // now see how far
	      if (way>0) {
		// incoming increasing so basic decreasing
		// if phase 0 then go to nearest bound
		double distance=upper[chosen]-solution[chosen];
		double basicDistance;
		if (!iphase) {
		  assert (iBasic==numberInSet);
		  assert (solution[iBasic]>upper[iBasic]);
		  basicDistance = solution[iBasic]-upper[iBasic];
		} else {
		  basicDistance = solution[iBasic]-lower[iBasic];
		}
		// need extra coding for unbounded
		assert (min(distance,basicDistance)<1.0e20);
		if (distance>basicDistance) {
		  // incoming becomes basic
		  solution[chosen] += basicDistance;
		  if (!iphase) 
		    solution[iBasic]=upper[iBasic];
		  else 
		    solution[iBasic]=lower[iBasic];
		  iBasic = chosen;
		} else {
		  // flip
		  solution[chosen]=upper[chosen];
		  solution[iBasic] -= distance;
		}
	      } else {
		// incoming decreasing so basic increasing
		// if phase 0 then go to nearest bound
		double distance=solution[chosen]-lower[chosen];
		double basicDistance;
		if (!iphase) {
		  assert (iBasic==numberInSet);
		  assert (solution[iBasic]<lower[iBasic]);
		  basicDistance = lower[iBasic]-solution[iBasic];
		} else {
		  basicDistance = upper[iBasic]-solution[iBasic];
		}
		// need extra coding for unbounded
		assert (min(distance,basicDistance)<1.0e20);
		if (distance>basicDistance) {
		  // incoming becomes basic
		  solution[chosen] -= basicDistance;
		  if (!iphase) 
		    solution[iBasic]=lower[iBasic];
		  else 
		    solution[iBasic]=upper[iBasic];
		  iBasic = chosen;
		} else {
		  // flip
		  solution[chosen]=lower[chosen];
		  solution[iBasic] += distance;
		}
	      }
	    }
	  }
	}
      }
    } 
    keyVariable_[iSet]=iBasic;
    int lastMarker = -(iBasic+1);
    next_[iBasic]=lastMarker;
    int last = iBasic;
    for (j=iStart;j<iEnd;j++) {
      if (model->getStatus(j)==ClpSimplex::basic&&j!=iBasic) {
	next_[last]=j;
	next_[j]=lastMarker;
      }
    }
  }
  delete [] lower;
  delete [] solution;
  delete [] upper;
  delete [] cost;
  // make sure matrix is in good shape
  matrix_->orderMatrix();
  // create effective rhs
}


