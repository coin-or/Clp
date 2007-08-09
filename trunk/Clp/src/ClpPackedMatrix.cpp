// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
//#define THREAD

#include "ClpSimplex.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpFactorization.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#endif
// at end to get min/max!
#include "ClpPackedMatrix.hpp"
#include "ClpMessage.hpp"
//#define COIN_PREFETCH
#ifdef COIN_PREFETCH
#if 1
#define coin_prefetch(mem) \
         __asm__ __volatile__ ("prefetchnta %0" : : "m" (*((char *)(mem))))
#else
#define coin_prefetch(mem) \
         __asm__ __volatile__ ("prefetch %0" : : "m" (*((char *)(mem))))
#endif
#else
// dummy
#define coin_prefetch(mem)
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix () 
  : ClpMatrixBase(),
    matrix_(NULL),
    numberActiveColumns_(0),
    flags_(2),
    rowCopy_(NULL),
    columnCopy_(NULL)
{
  setType(1);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix (const ClpPackedMatrix & rhs) 
: ClpMatrixBase(rhs)
{  
  matrix_ = new CoinPackedMatrix(*(rhs.matrix_));
  numberActiveColumns_ = rhs.numberActiveColumns_;
  flags_ = rhs.flags_;
  int numberRows = getNumRows();
  if (rhs.rhsOffset_&&numberRows) {
    rhsOffset_ = ClpCopyOfArray(rhs.rhsOffset_,numberRows);
  } else {
    rhsOffset_=NULL;
  }
  if (rhs.rowCopy_) {
    assert ((flags_&4)!=0);
    rowCopy_ = new ClpPackedMatrix2(*rhs.rowCopy_);
  } else {
    rowCopy_ = NULL;
  }
  if (rhs.columnCopy_) {
    assert ((flags_&8)!=0);
    columnCopy_ = new ClpPackedMatrix3(*rhs.columnCopy_);
  } else {
    columnCopy_ = NULL;
  }
}

//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix (CoinPackedMatrix * rhs) 
: ClpMatrixBase()
{  
  matrix_ = rhs;
  flags_ = 2;
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  columnCopy_=NULL;
  setType(1);
  
}

ClpPackedMatrix::ClpPackedMatrix (const CoinPackedMatrix & rhs) 
: ClpMatrixBase()
{  
  matrix_ = new CoinPackedMatrix(rhs);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  flags_ = 2;
  columnCopy_=NULL;
  setType(1);
  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPackedMatrix::~ClpPackedMatrix ()
{
  delete matrix_;
  delete rowCopy_;
  delete columnCopy_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPackedMatrix &
ClpPackedMatrix::operator=(const ClpPackedMatrix& rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete matrix_;
    matrix_ = new CoinPackedMatrix(*(rhs.matrix_));
    numberActiveColumns_ = rhs.numberActiveColumns_;
    flags_ = rhs.flags_;
    delete rowCopy_;
    delete columnCopy_;
    if (rhs.rowCopy_) {
      assert ((flags_&4)!=0);
      rowCopy_ = new ClpPackedMatrix2(*rhs.rowCopy_);
    } else {
      rowCopy_ = NULL;
    }
    if (rhs.columnCopy_) {
      assert ((flags_&8)!=0);
      columnCopy_ = new ClpPackedMatrix3(*rhs.columnCopy_);
    } else {
      columnCopy_ = NULL;
    }
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase * ClpPackedMatrix::clone() const
{
  return new ClpPackedMatrix(*this);
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase * 
ClpPackedMatrix::subsetClone (int numberRows, const int * whichRows,
			      int numberColumns, 
			      const int * whichColumns) const 
{
  return new ClpPackedMatrix(*this, numberRows, whichRows,
				   numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpPackedMatrix::ClpPackedMatrix (
		       const ClpPackedMatrix & rhs,
		       int numberRows, const int * whichRows,
		       int numberColumns, const int * whichColumns)
: ClpMatrixBase(rhs)
{
  matrix_ = new CoinPackedMatrix(*(rhs.matrix_),numberRows,whichRows,
				 numberColumns,whichColumns);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  flags_ = rhs.flags_;
  columnCopy_=NULL;
}
ClpPackedMatrix::ClpPackedMatrix (
		       const CoinPackedMatrix & rhs,
		       int numberRows, const int * whichRows,
		       int numberColumns, const int * whichColumns)
: ClpMatrixBase()
{
  matrix_ = new CoinPackedMatrix(rhs,numberRows,whichRows,
				 numberColumns,whichColumns);
  numberActiveColumns_ = matrix_->getNumCols();
  rowCopy_ = NULL;
  flags_ = 2;
  columnCopy_=NULL;
  setType(1);
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase * 
ClpPackedMatrix::reverseOrderedCopy() const
{
  ClpPackedMatrix * copy = new ClpPackedMatrix();
  copy->matrix_= new CoinPackedMatrix();
  copy->matrix_->setExtraGap(0.0);
  copy->matrix_->setExtraMajor(0.0);
  copy->matrix_->reverseOrderedCopyOf(*matrix_);
  //copy->matrix_->removeGaps();
  copy->numberActiveColumns_ = copy->matrix_->getNumCols();
  copy->flags_ = flags_&(~2); // no gaps
  return copy;
}
//unscaled versions
void 
ClpPackedMatrix::times(double scalar,
		   const double * x, double * y) const
{
  int iRow,iColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  //memset(y,0,matrix_->getNumRows()*sizeof(double));
  for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
    CoinBigIndex j;
    double value = scalar*x[iColumn];
    if (value) {
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow=row[j];
	y[iRow] += value*elementByColumn[j];
      }
    }
  }
}
void 
ClpPackedMatrix::transposeTimes(double scalar,
				const double * x, double * y) const
{
  int iColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  if (!(flags_&2)) {
    if (scalar==1.0) {
      CoinBigIndex start=columnStart[0];
      for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	CoinBigIndex j;
	CoinBigIndex next=columnStart[iColumn+1];
	double value=y[iColumn];
	for (j=start;j<next;j++) {
	  int jRow=row[j];
	  value += x[jRow]*elementByColumn[j];
	}
	start=next;
	y[iColumn] = value;
      }
    } else if (scalar==-1.0) {
      CoinBigIndex start=columnStart[0];
      for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	CoinBigIndex j;
	CoinBigIndex next=columnStart[iColumn+1];
	double value=y[iColumn];
	for (j=start;j<next;j++) {
	  int jRow=row[j];
	  value -= x[jRow]*elementByColumn[j];
	}
	start=next;
	y[iColumn] = value;
      }
    } else {
      CoinBigIndex start=columnStart[0];
      for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	CoinBigIndex j;
	CoinBigIndex next=columnStart[iColumn+1];
	double value=0.0;
	for (j=start;j<next;j++) {
	  int jRow=row[j];
	  value += x[jRow]*elementByColumn[j];
	}
	start=next;
	y[iColumn] += value*scalar;
      }
    }
  } else {
    for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
      CoinBigIndex j;
      double value=0.0;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int jRow=row[j];
	value += x[jRow]*elementByColumn[j];
      }
      y[iColumn] += value*scalar;
    }
  }
}
void 
ClpPackedMatrix::times(double scalar,
		       const double * x, double * y,
		       const double * rowScale, 
		       const double * columnScale) const
{
  if (rowScale) {
    int iRow,iColumn;
    // get matrix data pointers
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths(); 
    const double * elementByColumn = matrix_->getElements();
    for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
      CoinBigIndex j;
      double value = x[iColumn];
      if (value) {
	// scaled
	value *= scalar*columnScale[iColumn];
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  iRow=row[j];
	  y[iRow] += value*elementByColumn[j]*rowScale[iRow];
	}
      }
    }
  } else {
    times(scalar,x,y);
  }
}
void 
ClpPackedMatrix::transposeTimes( double scalar,
				 const double * x, double * y,
				 const double * rowScale, 
				 const double * columnScale,
				 double * spare) const
{
  if (rowScale) {
    int iColumn;
    // get matrix data pointers
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths(); 
    const double * elementByColumn = matrix_->getElements();
    if (!spare) {
      if (!(flags_&2)) {
	CoinBigIndex start=columnStart[0];
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex j;
	  CoinBigIndex next=columnStart[iColumn+1];
	  double value=0.0;
	  // scaled
	  for (j=start;j<next;j++) {
	    int jRow=row[j];
	    value += x[jRow]*elementByColumn[j]*rowScale[jRow];
	  }
	  start=next;
	  y[iColumn] += value*scalar*columnScale[iColumn];
	}
      } else {
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex j;
	  double value=0.0;
	  // scaled
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int jRow=row[j];
	    value += x[jRow]*elementByColumn[j]*rowScale[jRow];
	  }
	  y[iColumn] += value*scalar*columnScale[iColumn];
	}
      }
    } else {
      // can use spare region
      int iRow;
      int numberRows = getNumRows();
      for (iRow=0;iRow<numberRows;iRow++)
	spare[iRow] = x[iRow]*rowScale[iRow];
      if (!(flags_&2)) {
	CoinBigIndex start=columnStart[0];
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex j;
	  CoinBigIndex next=columnStart[iColumn+1];
	  double value=0.0;
	  // scaled
	  for (j=start;j<next;j++) {
	    int jRow=row[j];
	    value += spare[jRow]*elementByColumn[j];
	  }
	  start=next;
	  y[iColumn] += value*scalar*columnScale[iColumn];
	}
      } else {
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex j;
	  double value=0.0;
	  // scaled
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int jRow=row[j];
	    value += spare[jRow]*elementByColumn[j];
	  }
	  y[iColumn] += value*scalar*columnScale[iColumn];
	}
      }
      // no need to zero out
      //for (iRow=0;iRow<numberRows;iRow++)
      //spare[iRow] = 0.0;
    }
  } else {
    transposeTimes(scalar,x,y);
  }
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix::transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  int numberRows = model->numberRows();
#ifndef NO_RTTI
  ClpPackedMatrix* rowCopy =
    dynamic_cast< ClpPackedMatrix*>(model->rowCopy());
#else
  ClpPackedMatrix* rowCopy =
    static_cast< ClpPackedMatrix*>(model->rowCopy());
#endif
  bool packed = rowArray->packedMode();
  double factor = 0.30;
  // We may not want to do by row if there may be cache problems
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberActiveColumns_*sizeof(double)>1000000) {
    if (numberRows*10<numberActiveColumns_)
      factor *= 0.333333333;
    else if (numberRows*4<numberActiveColumns_)
      factor *= 0.5;
    else if (numberRows*2<numberActiveColumns_)
      factor *= 0.66666666667;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  // if not packed then bias a bit more towards by column
  if (!packed)
    factor *= 0.9;
  assert (!y->getNumElements());
  double multiplierX=0.8;
  double factor2 = factor*multiplierX;
  if (packed&&rowCopy_&&numberInRowArray>2&&numberInRowArray>factor2*numberRows&&
      numberInRowArray<0.9*numberRows&&0) {
    rowCopy_->transposeTimes(model,rowCopy->getPackedMatrix(),rowArray,y,columnArray);
    return;
  }
  if (numberInRowArray>factor*numberRows||!rowCopy) {
    // do by column
    // If no gaps - can do a bit faster
    if (!(flags_&2)||columnCopy_) {
      transposeTimesByColumn( model,  scalar,
			      rowArray, y, columnArray);
      return;
    }
    int iColumn;
    // get matrix data pointers
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths(); 
    const double * elementByColumn = matrix_->getElements();
    const double * rowScale = model->rowScale();
    if (packed) {
      // need to expand pi into y
      assert(y->capacity()>=numberRows);
      double * piOld = pi;
      pi = y->denseVector();
      const int * whichRow = rowArray->getIndices();
      int i;
      if (!rowScale) {
	// modify pi so can collapse to one loop
        if (scalar==-1.0) {
          for (i=0;i<numberInRowArray;i++) {
            int iRow = whichRow[i];
            pi[iRow]=-piOld[i];
          }
        } else {
          for (i=0;i<numberInRowArray;i++) {
            int iRow = whichRow[i];
            pi[iRow]=scalar*piOld[i];
          }
        }
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  double value = 0.0;
	  CoinBigIndex j;
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j];
	  }
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      } else {
	// scaled
	// modify pi so can collapse to one loop
        if (scalar==-1.0) {
          for (i=0;i<numberInRowArray;i++) {
            int iRow = whichRow[i];
            pi[iRow]=-piOld[i]*rowScale[iRow];
          }
        } else {
          for (i=0;i<numberInRowArray;i++) {
            int iRow = whichRow[i];
            pi[iRow]=scalar*piOld[i]*rowScale[iRow];
          }
        }
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  double value = 0.0;
	  CoinBigIndex j;
	  const double * columnScale = model->columnScale();
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j];
	  }
	  value *= columnScale[iColumn];
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
      // zero out
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=0.0;
      }
    } else {
      if (!rowScale) {
	if (scalar==-1.0) {
	  for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	    double value = 0.0;
	    CoinBigIndex j;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      value += pi[iRow]*elementByColumn[j];
	    }
	    if (fabs(value)>zeroTolerance) {
	      index[numberNonZero++]=iColumn;
	      array[iColumn]=-value;
	    }
	  }
	} else {
	  for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	    double value = 0.0;
	    CoinBigIndex j;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      value += pi[iRow]*elementByColumn[j];
	    }
	    value *= scalar;
	    if (fabs(value)>zeroTolerance) {
	      index[numberNonZero++]=iColumn;
	      array[iColumn]=value;
	    }
	  }
	}
      } else {
	// scaled
	if (scalar==-1.0) {
	  for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	    double value = 0.0;
	    CoinBigIndex j;
	    const double * columnScale = model->columnScale();
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	    }
	    value *= columnScale[iColumn];
	    if (fabs(value)>zeroTolerance) {
	      index[numberNonZero++]=iColumn;
	      array[iColumn]=-value;
	    }
	  }
	} else {
	  for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	    double value = 0.0;
	    CoinBigIndex j;
	    const double * columnScale = model->columnScale();
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	    }
	    value *= scalar*columnScale[iColumn];
	    if (fabs(value)>zeroTolerance) {
	      index[numberNonZero++]=iColumn;
	      array[iColumn]=value;
	    }
	  }
	}
      }
    }
    columnArray->setNumElements(numberNonZero);
    y->setNumElements(0);
  } else {
    // do by row
    rowCopy->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  }
  if (packed)
    columnArray->setPackedMode(true);
  if (0) {
    columnArray->checkClean();
    int numberNonZero=columnArray->getNumElements();;
    int * index = columnArray->getIndices();
    double * array = columnArray->denseVector();
    int i;
    for (i=0;i<numberNonZero;i++) {
      int j=index[i];
      double value;
      if (packed)
	value=array[i];
      else
	value=array[j];
      printf("Ti %d %d %g\n",i,j,value);
    }
  }
}
/* Return <code>x * scalar * A + y</code> in <code>z</code>. 
   Note - If x packed mode - then z packed mode
   This does by column and knows no gaps
   Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix::transposeTimesByColumn(const ClpSimplex * model, double scalar,
					const CoinIndexedVector * rowArray,
					CoinIndexedVector * y,
					CoinIndexedVector * columnArray) const
{
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  bool packed = rowArray->packedMode();
  // do by column
  int iColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  assert (!y->getNumElements());
  assert (numberActiveColumns_>0);
  if (packed) {
    // need to expand pi into y
    assert(y->capacity()>=model->numberRows());
    double * piOld = pi;
    pi = y->denseVector();
    const int * whichRow = rowArray->getIndices();
    int i;
    if (!rowScale) {
      // modify pi so can collapse to one loop
      if (scalar==-1.0) {
        for (i=0;i<numberInRowArray;i++) {
          int iRow = whichRow[i];
          pi[iRow]=-piOld[i];
        }
      } else {
        for (i=0;i<numberInRowArray;i++) {
          int iRow = whichRow[i];
          pi[iRow]=scalar*piOld[i];
        }
      }
#if 0
      double value = 0.0;
      CoinBigIndex j;
      CoinBigIndex end = columnStart[1];
      for (j=columnStart[0]; j<end;j++) {
	int iRow = row[j];
	value += pi[iRow]*elementByColumn[j];
      }
      for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	CoinBigIndex start = end;
	end = columnStart[iColumn+2];
	if (fabs(value)>zeroTolerance) {
	  array[numberNonZero]=value;
	  index[numberNonZero++]=iColumn;
	}
	value = 0.0;
	for (j=start; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j];
	}
      }
      if (fabs(value)>zeroTolerance) {
	array[numberNonZero]=value;
	index[numberNonZero++]=iColumn;
      }
#else
      if (!columnCopy_)
	gutsOfTransposeTimesUnscaled(pi,columnArray,zeroTolerance);
      else
	columnCopy_->transposeTimes(model,pi,columnArray);
      numberNonZero = columnArray->getNumElements();
#endif
    } else {
      // scaled
      // modify pi so can collapse to one loop
      if (scalar==-1.0) {
        for (i=0;i<numberInRowArray;i++) {
          int iRow = whichRow[i];
          pi[iRow]=-piOld[i]*rowScale[iRow];
        }
      } else {
        for (i=0;i<numberInRowArray;i++) {
          int iRow = whichRow[i];
          pi[iRow]=scalar*piOld[i]*rowScale[iRow];
        }
      }
      const double * columnScale = model->columnScale();
#if 0
      double value = 0.0;
      double scale=columnScale[0];
      CoinBigIndex j;
      CoinBigIndex end = columnStart[1];
      for (j=columnStart[0]; j<end;j++) {
	int iRow = row[j];
	value += pi[iRow]*elementByColumn[j];
      }
      for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	value *= scale;
	CoinBigIndex start = end;
	scale = columnScale[iColumn+1];
	end = columnStart[iColumn+2];
	if (fabs(value)>zeroTolerance) {
	  array[numberNonZero]=value;
	  index[numberNonZero++]=iColumn;
	}
	value = 0.0;
	if (model->getColumnStatus(iColumn+1)!=ClpSimplex::basic) {
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j];
	  }
	}
      }
      value *= scale;
      if (fabs(value)>zeroTolerance) {
	array[numberNonZero]=value;
	index[numberNonZero++]=iColumn;
      }
#else
      if (!columnCopy_)
	gutsOfTransposeTimesScaled(pi,columnScale,columnArray,zeroTolerance);
      else
	columnCopy_->transposeTimes(model,pi,columnArray);
      numberNonZero = columnArray->getNumElements();
#endif
    }
    // zero out
    int numberRows = model->numberRows();
    if (numberInRowArray*4<numberRows) {
      for (i=0;i<numberInRowArray;i++) {
        int iRow = whichRow[i];
        pi[iRow]=0.0;
      }
    } else {
      CoinZeroN(pi,numberRows);
    }
  } else {
    if (!rowScale) {
      if (scalar==-1.0) {
	double value = 0.0;
	CoinBigIndex j;
	CoinBigIndex end = columnStart[1];
	for (j=columnStart[0]; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j];
	}
	for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+2];
	  if (fabs(value)>zeroTolerance) {
	    array[iColumn]=-value;
	    index[numberNonZero++]=iColumn;
	  }
	  value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j];
	  }
	}
	if (fabs(value)>zeroTolerance) {
	  array[iColumn]=-value;
	  index[numberNonZero++]=iColumn;
	}
      } else {
	double value = 0.0;
	CoinBigIndex j;
	CoinBigIndex end = columnStart[1];
	for (j=columnStart[0]; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j];
	}
	for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	  value *= scalar;
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+2];
	  if (fabs(value)>zeroTolerance) {
	    array[iColumn]=value;
	    index[numberNonZero++]=iColumn;
	  }
	  value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j];
	  }
	}
	value *= scalar;
	if (fabs(value)>zeroTolerance) {
	  array[iColumn]=value;
	  index[numberNonZero++]=iColumn;
	}
      }
    } else {
      // scaled
      if (scalar==-1.0) {
	const double * columnScale = model->columnScale();
	double value = 0.0;
	double scale=columnScale[0];
	CoinBigIndex j;
	CoinBigIndex end = columnStart[1];
	for (j=columnStart[0]; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	}
	for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	  value *= scale;
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+2];
	  scale = columnScale[iColumn+1];
	  if (fabs(value)>zeroTolerance) {
	    array[iColumn]=-value;
	    index[numberNonZero++]=iColumn;
	  }
	  value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	  }
	}
	value *= scale;
	if (fabs(value)>zeroTolerance) {
	  array[iColumn]=-value;
	  index[numberNonZero++]=iColumn;
	}
      } else {
	const double * columnScale = model->columnScale();
	double value = 0.0;
	double scale=columnScale[0]*scalar;
	CoinBigIndex j;
	CoinBigIndex end = columnStart[1];
	for (j=columnStart[0]; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	}
	for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
	  value *= scale;
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+2];
	  scale = columnScale[iColumn+1]*scalar;
	  if (fabs(value)>zeroTolerance) {
	    array[iColumn]=value;
	    index[numberNonZero++]=iColumn;
	  }
	  value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	  }
	}
	value *= scale;
	if (fabs(value)>zeroTolerance) {
	  array[iColumn]=value;
	  index[numberNonZero++]=iColumn;
	}
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
  if (packed)
    columnArray->setPackedMode(true);
}
/* Return <code>x * A + y</code> in <code>z</code>. 
	Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix::transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * rowArray,
			      CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  const int * column = getIndices();
  const CoinBigIndex * rowStart = getVectorStarts();
  const double * element = getElements();
  const int * whichRow = rowArray->getIndices();
  bool packed = rowArray->packedMode();
  if (numberInRowArray>2) {
    // do by rows
    // ** Row copy is already scaled
    int iRow;
    int i;
    int numberOriginal = 0;
    if (packed) {
#if 0
      numberNonZero=0;
      // and set up mark as char array
      char * marked = (char *) (index+columnArray->capacity());
      double * array2 = y->denseVector();
#ifdef CLP_DEBUG
      for (i=0;i<numberActiveColumns_;i++) {
	assert(!marked[i]);
	assert(!array2[i]);
      }
#endif      
      for (i=0;i<numberInRowArray;i++) {
	iRow = whichRow[i]; 
	double value = pi[i]*scalar;
	CoinBigIndex j;
	for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array2[iColumn] += value*element[j];
	}
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	if (marked[iColumn]) {
	  double value = array2[iColumn];
	  array2[iColumn]=0.0; 
	  marked[iColumn]=0;
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
#else
      gutsOfTransposeTimesByRowGE3(rowArray,columnArray,y,zeroTolerance,scalar);
      numberNonZero = columnArray->getNumElements();
#endif
    } else {
      double * markVector = y->denseVector();
      numberNonZero=0;
      // and set up mark as char array
      char * marked = (char *) markVector;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	marked[iColumn]=0;
      }
      
      for (i=0;i<numberInRowArray;i++) {
	iRow = whichRow[i]; 
	double value = pi[iRow]*scalar;
	CoinBigIndex j;
	for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	  int iColumn = column[j];
	  if (!marked[iColumn]) {
	    marked[iColumn]=1;
	    index[numberNonZero++]=iColumn;
	  }
	  array[iColumn] += value*element[j];
	}
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	marked[iColumn]=0;
	if (fabs(array[iColumn])>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	} else {
	  array[iColumn]=0.0;
	}
      }
    }
  } else if (numberInRowArray==2) {
    // do by rows when two rows
    int numberOriginal;
    int i;
    CoinBigIndex j;
    numberNonZero=0;

    double value;
    if (packed) {
#if 0
      int iRow0 = whichRow[0]; 
      int iRow1 = whichRow[1]; 
      double pi0 = pi[0];
      double pi1 = pi[1];
      if (rowStart[iRow0+1]-rowStart[iRow0]>
	  rowStart[iRow1+1]-rowStart[iRow1]) {
	// do one with fewer first
	iRow0=iRow1;
	iRow1=whichRow[0];
	pi0=pi1;
	pi1=pi[0];
      }
      // and set up mark as char array
      char * marked = (char *) (index+columnArray->capacity());
      int * lookup = y->getIndices();
      value = pi0*scalar;
      for (j=rowStart[iRow0];j<rowStart[iRow0+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	array[numberNonZero] = value2;
	marked[iColumn]=1;
	lookup[iColumn]=numberNonZero;
	index[numberNonZero++]=iColumn;
      }
      numberOriginal = numberNonZero;
      value = pi1*scalar;
      for (j=rowStart[iRow1];j<rowStart[iRow1+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	// I am assuming no zeros in matrix
	if (marked[iColumn]) {
	  int iLookup = lookup[iColumn];
	  array[iLookup] += value2;
	} else {
	  if (fabs(value2)>zeroTolerance) {
	    array[numberNonZero] = value2;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
      // get rid of tiny values and zero out marked
      int nDelete=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	marked[iColumn]=0;
	if (fabs(array[i])<=zeroTolerance) 
	  nDelete++;
      }
      if (nDelete) {
	numberOriginal=numberNonZero;
	numberNonZero=0;
	for (i=0;i<numberOriginal;i++) {
	  int iColumn = index[i];
	  double value = array[i];
	  array[i]=0.0;
	  if (fabs(value)>zeroTolerance) {
	    array[numberNonZero]=value;
	    index[numberNonZero++]=iColumn;
	  }
	}
      }
#else
      gutsOfTransposeTimesByRowEQ2(rowArray,columnArray,y,zeroTolerance,scalar);
      numberNonZero = columnArray->getNumElements();
#endif
    } else {
      int iRow = whichRow[0]; 
      value = pi[iRow]*scalar;
      for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	index[numberNonZero++]=iColumn;
	array[iColumn] = value2;
      }
      iRow = whichRow[1]; 
      value = pi[iRow]*scalar;
      for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	// I am assuming no zeros in matrix
	if (array[iColumn])
	  value2 += array[iColumn];
	else
	  index[numberNonZero++]=iColumn;
	array[iColumn] = value2;
      }
      // get rid of tiny values and zero out marked
      numberOriginal=numberNonZero;
      numberNonZero=0;
      for (i=0;i<numberOriginal;i++) {
	int iColumn = index[i];
	if (fabs(array[iColumn])>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	} else {
	  array[iColumn]=0.0;
	}
      }
    }
  } else if (numberInRowArray==1) {
    // Just one row
    int iRow=rowArray->getIndices()[0];
    numberNonZero=0;
    CoinBigIndex j;
    if (packed) {
#if 0
      double value = pi[0]*scalar;
      for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	if (fabs(value2)>zeroTolerance) {
	  array[numberNonZero] = value2;
	  index[numberNonZero++]=iColumn;
	}
      }
#else
      gutsOfTransposeTimesByRowEQ1(rowArray,columnArray,zeroTolerance,scalar);
      numberNonZero = columnArray->getNumElements();
#endif
    } else {
      double value = pi[iRow]*scalar;
      for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	if (fabs(value2)>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	  array[iColumn] = value2;
	}
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
}
// Meat of transposeTimes by column when not scaled
void 
ClpPackedMatrix::gutsOfTransposeTimesUnscaled(const double * pi,CoinIndexedVector * output, const double zeroTolerance) const
{
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * elementByColumn = matrix_->getElements();
  double value = 0.0;
  CoinBigIndex j;
  CoinBigIndex end = columnStart[1];
  for (j=columnStart[0]; j<end;j++) {
    int iRow = row[j];
    value += pi[iRow]*elementByColumn[j];
  }
  int iColumn;
  for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
    CoinBigIndex start = end;
    end = columnStart[iColumn+2];
    if (fabs(value)>zeroTolerance) {
      array[numberNonZero]=value;
      index[numberNonZero++]=iColumn;
    }
    value = 0.0;
    for (j=start; j<end;j++) {
      int iRow = row[j];
      value += pi[iRow]*elementByColumn[j];
    }
  }
  if (fabs(value)>zeroTolerance) {
    array[numberNonZero]=value;
    index[numberNonZero++]=iColumn;
  }
  output->setNumElements(numberNonZero);
}
// Meat of transposeTimes by column when scaled
void 
ClpPackedMatrix::gutsOfTransposeTimesScaled(const double * pi,const double * columnScale, 
					   CoinIndexedVector * output, const double zeroTolerance) const
{
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * elementByColumn = matrix_->getElements();
  double value = 0.0;
  double scale=columnScale[0];
  CoinBigIndex j;
  CoinBigIndex end = columnStart[1];
  for (j=columnStart[0]; j<end;j++) {
    int iRow = row[j];
    value += pi[iRow]*elementByColumn[j];
  }
  int iColumn;
  for (iColumn=0;iColumn<numberActiveColumns_-1;iColumn++) {
    value *= scale;
    CoinBigIndex start = end;
    scale = columnScale[iColumn+1];
    end = columnStart[iColumn+2];
    if (fabs(value)>zeroTolerance) {
      array[numberNonZero]=value;
      index[numberNonZero++]=iColumn;
    }
    value = 0.0;
    for (j=start; j<end;j++) {
      int iRow = row[j];
      value += pi[iRow]*elementByColumn[j];
    }
  }
  value *= scale;
  if (fabs(value)>zeroTolerance) {
    array[numberNonZero]=value;
    index[numberNonZero++]=iColumn;
  }
  output->setNumElements(numberNonZero);
}
// Meat of transposeTimes by row n > 2 if packed
void 
ClpPackedMatrix::gutsOfTransposeTimesByRowGE3(const CoinIndexedVector * piVector, CoinIndexedVector * output,
					     CoinIndexedVector * spareVector, const double tolerance, const double scalar) const
{
  double * pi = piVector->denseVector();
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  int numberInRowArray = piVector->getNumElements();
  const int * column = getIndices();
  const CoinBigIndex * rowStart = getVectorStarts();
  const double * element = getElements();
  const int * whichRow = piVector->getIndices();
  // ** Row copy is already scaled
  int iRow;
  int i;
  // and set up mark as char array
  char * marked = (char *) (index+output->capacity());
  int * lookup = spareVector->getIndices();
#if 1
  for (i=0;i<numberInRowArray;i++) {
    iRow = whichRow[i]; 
    double value = pi[i]*scalar;
    CoinBigIndex j;
    for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
      int iColumn = column[j];
      double elValue = element[j];
      if (!marked[iColumn]) {
	marked[iColumn]=1;
	lookup[iColumn]=numberNonZero;
	array[numberNonZero] = value*elValue;
	index[numberNonZero++]=iColumn;
      } else {
	int k = lookup[iColumn];
	array[k] += value*elValue;
      }
    }
  }
#else
  int nextRow = whichRow[0];
  CoinBigIndex nextStart = rowStart[nextRow];
  CoinBigIndex nextEnd = rowStart[nextRow+1];
  whichRow[numberInRowArray]=0; // for electricfence etc
  for (i=0;i<numberInRowArray;i++) {
    iRow = nextRow; 
    nextRow = whichRow[i+1];
    CoinBigIndex start = nextStart;
    CoinBigIndex end = nextEnd;
    double value = pi[i]*scalar;
    CoinBigIndex j;
    nextRow = whichRow[i+1];
    nextStart = rowStart[nextRow];
    nextEnd = rowStart[nextRow+1];
    for (j=start;j<end;j++) {
      int iColumn = column[j];
      double elValue = element[j];
      if (!marked[iColumn]) {
	marked[iColumn]=1;
	lookup[iColumn]=numberNonZero;
	array[numberNonZero] = value*elValue;
	index[numberNonZero++]=iColumn;
      } else {
	int k = lookup[iColumn];
	array[k] += value*elValue;
      }
    }
  }
#endif
  // get rid of tiny values and zero out marked
  for (i=0;i<numberNonZero;i++) {
    int iColumn = index[i];
    marked[iColumn]=0;
    double value = array[i];
    while (fabs(value)<=tolerance) {
      numberNonZero--;
      value = array[numberNonZero];
      iColumn = index[numberNonZero];
      marked[iColumn]=0;
      if (i<numberNonZero) {
	array[numberNonZero]=0.0;
	array[i] = value;
	index[i] = iColumn;
      } else {
	array[i]=0.0;
	value =1.0; // to force end of while
      }
    }
  }
  output->setNumElements(numberNonZero);
  spareVector->setNumElements(0);
}
// Meat of transposeTimes by row n == 2 if packed
void 
ClpPackedMatrix::gutsOfTransposeTimesByRowEQ2(const CoinIndexedVector * piVector, CoinIndexedVector * output,
				   CoinIndexedVector * spareVector, const double tolerance, const double scalar) const
{
  double * pi = piVector->denseVector();
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  const int * column = getIndices();
  const CoinBigIndex * rowStart = getVectorStarts();
  const double * element = getElements();
  const int * whichRow = piVector->getIndices();
  int iRow0 = whichRow[0]; 
  int iRow1 = whichRow[1]; 
  double pi0 = pi[0];
  double pi1 = pi[1];
  if (rowStart[iRow0+1]-rowStart[iRow0]>
      rowStart[iRow1+1]-rowStart[iRow1]) {
    // do one with fewer first
    iRow0=iRow1;
    iRow1=whichRow[0];
    pi0=pi1;
    pi1=pi[0];
  }
  // and set up mark as char array
  char * marked = (char *) (index+output->capacity());
  int * lookup = spareVector->getIndices();
  double value = pi0*scalar;
  CoinBigIndex j;
  for (j=rowStart[iRow0];j<rowStart[iRow0+1];j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value*elValue;
    array[numberNonZero] = value2;
    marked[iColumn]=1;
    lookup[iColumn]=numberNonZero;
    index[numberNonZero++]=iColumn;
  }
  int numberOriginal = numberNonZero;
  value = pi1*scalar;
  for (j=rowStart[iRow1];j<rowStart[iRow1+1];j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value*elValue;
    // I am assuming no zeros in matrix
    if (marked[iColumn]) {
      int iLookup = lookup[iColumn];
      array[iLookup] += value2;
    } else {
      if (fabs(value2)>tolerance) {
	array[numberNonZero] = value2;
	index[numberNonZero++]=iColumn;
      }
    }
  }
  // get rid of tiny values and zero out marked
  int i;
  int iFirst=numberNonZero;
  for (i=0;i<numberOriginal;i++) {
    int iColumn = index[i];
    marked[iColumn]=0;
    if (fabs(array[i])<=tolerance) {
      if (numberNonZero>numberOriginal) {
	numberNonZero--;
	double value = array[numberNonZero];
	array[numberNonZero]=0.0;
	array[i]=value;
	index[i]=index[numberNonZero];
      } else {
	iFirst = i;
      }
    }
  }
  
  if (iFirst<numberNonZero) {
    int n=iFirst;
    for (i=n;i<numberOriginal;i++) {
      int iColumn = index[i];
      double value = array[i];
      array[i]=0.0;
      if (fabs(value)>tolerance) {
	array[n]=value;
	index[n++]=iColumn;
      }
    }
    for (;i<numberNonZero;i++) {
      int iColumn = index[i];
      double value = array[i];
      array[i]=0.0;
      array[n]=value;
      index[n++]=iColumn;
    }
    numberNonZero=n;
  }
  output->setNumElements(numberNonZero);
  spareVector->setNumElements(0);
}
// Meat of transposeTimes by row n == 1 if packed
void 
ClpPackedMatrix::gutsOfTransposeTimesByRowEQ1(const CoinIndexedVector * piVector, CoinIndexedVector * output,
				   const double tolerance, const double scalar) const
{
  double * pi = piVector->denseVector();
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  const int * column = getIndices();
  const CoinBigIndex * rowStart = getVectorStarts();
  const double * element = getElements();
  int iRow=piVector->getIndices()[0];
  numberNonZero=0;
  CoinBigIndex j;
  double value = pi[0]*scalar;
  for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
    int iColumn = column[j];
    double elValue = element[j];
    double value2 = value*elValue;
    if (fabs(value2)>tolerance) {
      array[numberNonZero] = value2;
      index[numberNonZero++]=iColumn;
    }
  }
  output->setNumElements(numberNonZero);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y.
   Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix::subsetTransposeTimes(const ClpSimplex * model,
			      const CoinIndexedVector * rowArray,
			      const CoinIndexedVector * y,
			      CoinIndexedVector * columnArray) const
{
  columnArray->clear();
  double * pi = rowArray->denseVector();
  double * array = columnArray->denseVector();
  int jColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  int numberToDo = y->getNumElements();
  const int * which = y->getIndices();
  assert (!rowArray->packedMode());
  columnArray->setPacked();
  if (!(flags_&2)&&numberToDo>5) {
    // no gaps and a reasonable number
    if (!rowScale) {
      int iColumn = which[0];
      double value = 0.0;
      CoinBigIndex j;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn+1];j++) {
	int iRow = row[j];
	value += pi[iRow]*elementByColumn[j];
      }
      for (jColumn=0;jColumn<numberToDo-1;jColumn++) {
	int iColumn = which[jColumn+1];
	CoinBigIndex start = columnStart[iColumn];
	CoinBigIndex end = columnStart[iColumn+1];
	array[jColumn]=value;
	value = 0.0;
	for (j=start; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j];
	}
      }
      array[jColumn]=value;
    } else {
      // scaled
      const double * columnScale = model->columnScale();
      int iColumn = which[0];
      double value = 0.0;
      double scale=columnScale[iColumn];
      CoinBigIndex j;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn+1];j++) {
	int iRow = row[j];
	value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
      }
      for (jColumn=0;jColumn<numberToDo-1;jColumn++) {
	int iColumn = which[jColumn+1];
	value *= scale;
	scale = columnScale[iColumn];
	CoinBigIndex start = columnStart[iColumn];
	CoinBigIndex end = columnStart[iColumn+1];
	array[jColumn]=value;
	value = 0.0;
	for (j=start; j<end;j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	}
      }
      value *= scale;
      array[jColumn]=value;
    }
  } else {
    // gaps
    if (!rowScale) {
      for (jColumn=0;jColumn<numberToDo;jColumn++) {
	int iColumn = which[jColumn];
	double value = 0.0;
	CoinBigIndex j;
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j];
	}
	array[jColumn]=value;
      }
    } else {
      // scaled
      const double * columnScale = model->columnScale();
      for (jColumn=0;jColumn<numberToDo;jColumn++) {
	int iColumn = which[jColumn];
	double value = 0.0;
	CoinBigIndex j;
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	}
	value *= columnScale[iColumn];
	array[jColumn]=value;
      }
    }
  }
}
/* Returns true if can combine transposeTimes and subsetTransposeTimes
   and if it would be faster */
bool 
ClpPackedMatrix::canCombine(const ClpSimplex * model,
                            const CoinIndexedVector * pi) const
{
  int numberInRowArray = pi->getNumElements();
  int numberRows = model->numberRows();
  bool packed = pi->packedMode();
  // factor should be smaller if doing both with two pi vectors 
  double factor = 0.30;
  // We may not want to do by row if there may be cache problems
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberActiveColumns_*sizeof(double)>1000000) {
    if (numberRows*10<numberActiveColumns_)
      factor *= 0.333333333;
    else if (numberRows*4<numberActiveColumns_)
      factor *= 0.5;
    else if (numberRows*2<numberActiveColumns_)
      factor *= 0.66666666667;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  // if not packed then bias a bit more towards by column
  if (!packed)
    factor *= 0.9;
  return ((numberInRowArray>factor*numberRows||!model->rowCopy())&&!(flags_&2));
}
#ifndef CLP_ALL_ONE_FILE
// These have to match ClpPrimalColumnSteepest version
#define reference(i)  (((reference[i>>5]>>(i&31))&1)!=0)
#endif
// Updates two arrays for steepest 
void 
ClpPackedMatrix::transposeTimes2(const ClpSimplex * model,
                                 const CoinIndexedVector * pi1, CoinIndexedVector * dj1,
                                 const CoinIndexedVector * pi2, CoinIndexedVector * dj2,
                                 CoinIndexedVector * spare,
                                 double referenceIn, double devex,
                                 // Array for exact devex to say what is in reference framework
                                 unsigned int * reference,
                                 double * weights, double scaleFactor)
{
  // put row of tableau in dj1
  double * pi = pi1->denseVector();
  int numberNonZero=0;
  int * index = dj1->getIndices();
  double * array = dj1->denseVector();
  int numberInRowArray = pi1->getNumElements();
  double zeroTolerance = model->factorization()->zeroTolerance();
  bool packed = pi1->packedMode();
  // do by column
  int iColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  assert (!spare->getNumElements());
  assert (numberActiveColumns_>0);
  double * piWeight = pi2->denseVector();
  assert (!pi2->packedMode());
  bool killDjs = (scaleFactor==0.0);
  if (!scaleFactor)
    scaleFactor=1.0;
  if (packed) {
    // need to expand pi into y
    assert(spare->capacity()>=model->numberRows());
    double * piOld = pi;
    pi = spare->denseVector();
    const int * whichRow = pi1->getIndices();
    int i;
    if (!rowScale) {
      // modify pi so can collapse to one loop
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=piOld[i];
      }
      if (!columnCopy_) {
	CoinBigIndex j;
	CoinBigIndex end = columnStart[0];
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+1];
	  ClpSimplex::Status status = model->getStatus(iColumn);
	  if (status==ClpSimplex::basic||status==ClpSimplex::isFixed) continue;
	  double value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value -= pi[iRow]*elementByColumn[j];
	  }
	  if (fabs(value)>zeroTolerance) {
	    // and do other array
	    double modification = 0.0;
	    for (j=start; j<end;j++) {
	      int iRow = row[j];
	      modification += piWeight[iRow]*elementByColumn[j];
	    }
	    double thisWeight = weights[iColumn];
	    double pivot = value*scaleFactor;
	    double pivotSquared = pivot * pivot;
	    thisWeight += pivotSquared * devex + pivot * modification;
	    if (thisWeight<DEVEX_TRY_NORM) {
	      if (referenceIn<0.0) {
		// steepest
		thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
	      } else {
		// exact
		thisWeight = referenceIn*pivotSquared;
		if (reference(iColumn))
		  thisWeight += 1.0;
		thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
	      }
	    }
	    weights[iColumn] = thisWeight;
	    if (!killDjs) {
	      array[numberNonZero]=value;
	      index[numberNonZero++]=iColumn;
	    }
	  }
	}
      } else {
	// use special column copy
	// reset back
	if (killDjs)
	  scaleFactor=0.0;
	columnCopy_->transposeTimes2(model,pi,dj1,piWeight,referenceIn,devex,
				     reference,weights,scaleFactor);
	numberNonZero = dj1->getNumElements();
      }
    } else {
      // scaled
      // modify pi so can collapse to one loop
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=piOld[i]*rowScale[iRow];
      }
      // can also scale piWeight as not used again
      int numberWeight = pi2->getNumElements();
      const int * indexWeight = pi2->getIndices();
      for (i=0;i<numberWeight;i++) {
	int iRow = indexWeight[i];
	piWeight[iRow] *= rowScale[iRow];
      }
      if (!columnCopy_) {
	const double * columnScale = model->columnScale();
	CoinBigIndex j;
	CoinBigIndex end = columnStart[0];
	for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
	  CoinBigIndex start = end;
	  end = columnStart[iColumn+1];
	  ClpSimplex::Status status = model->getStatus(iColumn);
	  if (status==ClpSimplex::basic||status==ClpSimplex::isFixed) continue;
	  double scale=columnScale[iColumn];
	  double value = 0.0;
	  for (j=start; j<end;j++) {
	    int iRow = row[j];
	    value -= pi[iRow]*elementByColumn[j];
	  }
	  value *= scale;
	  if (fabs(value)>zeroTolerance) {
	    double modification = 0.0;
	    for (j=start; j<end;j++) {
	      int iRow = row[j];
	      modification += piWeight[iRow]*elementByColumn[j];
	    }
	    modification *= scale;
	    double thisWeight = weights[iColumn];
	    double pivot = value*scaleFactor;
	    double pivotSquared = pivot * pivot;
	    thisWeight += pivotSquared * devex + pivot * modification;
	    if (thisWeight<DEVEX_TRY_NORM) {
	      if (referenceIn<0.0) {
		// steepest
		thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
	      } else {
		// exact
		thisWeight = referenceIn*pivotSquared;
		if (reference(iColumn))
		  thisWeight += 1.0;
		thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
	      }
	    }
	    weights[iColumn] = thisWeight;
	    if (!killDjs) {
	      array[numberNonZero]=value;
	      index[numberNonZero++]=iColumn;
	    }
	  }
	}
      } else {
	// use special column copy
	// reset back
	if (killDjs)
	  scaleFactor=0.0;
	columnCopy_->transposeTimes2(model,pi,dj1,piWeight,referenceIn,devex,
				     reference,weights,scaleFactor);
	numberNonZero = dj1->getNumElements();
      }
    }
    // zero out
    int numberRows = model->numberRows();
    if (numberInRowArray*4<numberRows) {
      for (i=0;i<numberInRowArray;i++) {
        int iRow = whichRow[i];
        pi[iRow]=0.0;
      }
    } else {
      CoinZeroN(pi,numberRows);
    }
  } else {
    if (!rowScale) {
      CoinBigIndex j;
      CoinBigIndex end = columnStart[0];
      for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
        CoinBigIndex start = end;
        end = columnStart[iColumn+1];
        ClpSimplex::Status status = model->getStatus(iColumn);
        if (status==ClpSimplex::basic||status==ClpSimplex::isFixed) continue;
        double value = 0.0;
        for (j=start; j<end;j++) {
          int iRow = row[j];
          value -= pi[iRow]*elementByColumn[j];
        }
        if (fabs(value)>zeroTolerance) {
          // and do other array
          double modification = 0.0;
          for (j=start; j<end;j++) {
            int iRow = row[j];
            modification += piWeight[iRow]*elementByColumn[j];
          }
          double thisWeight = weights[iColumn];
          double pivot = value*scaleFactor;
          double pivotSquared = pivot * pivot;
          thisWeight += pivotSquared * devex + pivot * modification;
          if (thisWeight<DEVEX_TRY_NORM) {
            if (referenceIn<0.0) {
              // steepest
              thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
            } else {
              // exact
              thisWeight = referenceIn*pivotSquared;
              if (reference(iColumn))
                thisWeight += 1.0;
              thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
            }
          }
          weights[iColumn] = thisWeight;
          if (!killDjs) {
            array[iColumn]=value;
            index[numberNonZero++]=iColumn;
          }
        }
      }
    } else {
      // scaled
      // can also scale piWeight as not used again
      int numberWeight = pi2->getNumElements();
      const int * indexWeight = pi2->getIndices();
      for (int i=0;i<numberWeight;i++) {
	int iRow = indexWeight[i];
	piWeight[iRow] *= rowScale[iRow];
      }
      const double * columnScale = model->columnScale();
      CoinBigIndex j;
      CoinBigIndex end = columnStart[0];
      for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
        CoinBigIndex start = end;
        end = columnStart[iColumn+1];
        ClpSimplex::Status status = model->getStatus(iColumn);
        if (status==ClpSimplex::basic||status==ClpSimplex::isFixed) continue;
        double scale=columnScale[iColumn];
        double value = 0.0;
        for (j=start; j<end;j++) {
          int iRow = row[j];
          value -= pi[iRow]*elementByColumn[j]*rowScale[iRow];
        }
        value *= scale;
        if (fabs(value)>zeroTolerance) {
          double modification = 0.0;
          for (j=start; j<end;j++) {
            int iRow = row[j];
            modification += piWeight[iRow]*elementByColumn[j];
          }
          modification *= scale;
          double thisWeight = weights[iColumn];
          double pivot = value*scaleFactor;
          double pivotSquared = pivot * pivot;
          thisWeight += pivotSquared * devex + pivot * modification;
          if (thisWeight<DEVEX_TRY_NORM) {
            if (referenceIn<0.0) {
              // steepest
              thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
            } else {
              // exact
              thisWeight = referenceIn*pivotSquared;
              if (reference(iColumn))
                thisWeight += 1.0;
              thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
            }
          }
          weights[iColumn] = thisWeight;
          if (!killDjs) {
            array[iColumn]=value;
            index[numberNonZero++]=iColumn;
          }
        }
      }
    }
  }
  dj1->setNumElements(numberNonZero);
  spare->setNumElements(0);
  if (packed)
    dj1->setPackedMode(true);
}
// Updates second array for steepest and does devex weights
void 
ClpPackedMatrix::subsetTimes2(const ClpSimplex * model,
                              CoinIndexedVector * dj1,
                            const CoinIndexedVector * pi2, CoinIndexedVector * dj2,
                            double referenceIn, double devex,
                            // Array for exact devex to say what is in reference framework
                            unsigned int * reference,
                            double * weights, double scaleFactor)
{
  int number = dj1->getNumElements();
  const int * index = dj1->getIndices();
  double * array = dj1->denseVector();
  assert( dj1->packedMode());

  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  double * piWeight = pi2->denseVector();
  bool killDjs = (scaleFactor==0.0);
  if (!scaleFactor)
    scaleFactor=1.0;
  if (!rowScale) {
    for (int k=0;k<number;k++) {
      int iColumn = index[k];
      double pivot = array[k]*scaleFactor;
      if (killDjs)
        array[k]=0.0;
      // and do other array
      double modification = 0.0;
      for (CoinBigIndex j=columnStart[iColumn]; 
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow = row[j];
        modification += piWeight[iRow]*elementByColumn[j];
      }
      double thisWeight = weights[iColumn];
      double pivotSquared = pivot * pivot;
      thisWeight += pivotSquared * devex + pivot * modification;
      if (thisWeight<DEVEX_TRY_NORM) {
        if (referenceIn<0.0) {
          // steepest
          thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn*pivotSquared;
          if (reference(iColumn))
            thisWeight += 1.0;
          thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
        }
      }
      weights[iColumn] = thisWeight;
    }
  } else {
    // scaled
    const double * columnScale = model->columnScale();
    for (int k=0;k<number;k++) {
      int iColumn = index[k];
      double pivot = array[k]*scaleFactor;
      double scale=columnScale[iColumn];
      if (killDjs)
        array[k]=0.0;
      // and do other array
      double modification = 0.0;
      for (CoinBigIndex j=columnStart[iColumn]; 
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow = row[j];
        modification += piWeight[iRow]*elementByColumn[j]*rowScale[iRow];
      }
      double thisWeight = weights[iColumn];
      modification *= scale;
      double pivotSquared = pivot * pivot;
      thisWeight += pivotSquared * devex + pivot * modification;
      if (thisWeight<DEVEX_TRY_NORM) {
        if (referenceIn<0.0) {
          // steepest
          thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
        } else {
          // exact
          thisWeight = referenceIn*pivotSquared;
          if (reference(iColumn))
            thisWeight += 1.0;
          thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
        }
      }
      weights[iColumn] = thisWeight;
    }
  }
}
/// returns number of elements in column part of basis,
CoinBigIndex 
ClpPackedMatrix::countBasis(ClpSimplex * model,
			   const int * whichColumn, 
			   int numberBasic,
			    int & numberColumnBasic)
{
  const int * columnLength = matrix_->getVectorLengths(); 
  int i;
  CoinBigIndex numberElements=0;
  // just count - can be over so ignore zero problem
  for (i=0;i<numberColumnBasic;i++) {
    int iColumn = whichColumn[i];
    numberElements += columnLength[iColumn];
  }
  return numberElements;
}
void
ClpPackedMatrix::fillBasis(ClpSimplex * model,
			 const int * whichColumn, 
			 int & numberColumnBasic,
			 int * indexRowU, int * start,
			 int * rowCount, int * columnCount,
			 double * elementU)
{
  const int * columnLength = matrix_->getVectorLengths(); 
  int i;
  CoinBigIndex numberElements=start[0];
  // fill
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const double * elementByColumn = matrix_->getElements();
  if ((flags_&1)==0) {
    if (!rowScale) {
      // no scaling
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	CoinBigIndex j;
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow=row[j];
	  indexRowU[numberElements]=iRow;
	  rowCount[iRow]++;
	  elementU[numberElements++]=elementByColumn[j];
	}
	start[i+1]=numberElements;
	columnCount[i]=columnLength[iColumn];
      }
    } else {
      // scaling
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	CoinBigIndex j;
	double scale = columnScale[iColumn];
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  indexRowU[numberElements]=iRow;
	  rowCount[iRow]++;
	  elementU[numberElements++]=
	    elementByColumn[j]*scale*rowScale[iRow];
	}
	start[i+1]=numberElements;
	columnCount[i]=columnLength[iColumn];
      }
    }
  } else {
    // there are zero elements so need to look more closely
    if (!rowScale) {
      // no scaling
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	CoinBigIndex j;
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  double value = elementByColumn[j];
	  if (value) {
	    int iRow = row[j];
	    indexRowU[numberElements]=iRow;
	    rowCount[iRow]++;
	    elementU[numberElements++]=value;
	  }
	}
	start[i+1]=numberElements;
	columnCount[i]=numberElements-start[i];
      }
    } else {
      // scaling
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumnBasic;i++) {
	int iColumn = whichColumn[i];
	CoinBigIndex j;
	double scale = columnScale[iColumn];
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[i];j++) {
	  double value = elementByColumn[j];
	  if (value) {
	    int iRow = row[j];
	    indexRowU[numberElements]=iRow;
	    rowCount[iRow]++;
	    elementU[numberElements++]=value*scale*rowScale[iRow];
	  }
	}
	start[i+1]=numberElements;
	columnCount[i]=numberElements-start[i];
      }
    }
  }
}
// Creates scales for column copy (rowCopy in model may be modified)
int 
ClpPackedMatrix::scale(ClpModel * model) const 
{
#ifndef NDEBUG
  //checkFlags();
#endif
  int numberRows = model->numberRows(); 
  int numberColumns = matrix_->getNumCols();
  // If empty - return as sanityCheck will trap
  if (!numberRows||!numberColumns)
    return 1;
  ClpMatrixBase * rowCopyBase=model->rowCopy();
  if (!rowCopyBase) {
    // temporary copy
    rowCopyBase = reverseOrderedCopy();
  }
#ifndef NO_RTTI
  ClpPackedMatrix* rowCopy =
    dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
  // Make sure it is really a ClpPackedMatrix
  assert (rowCopy!=NULL);
#else
  ClpPackedMatrix* rowCopy =
    static_cast< ClpPackedMatrix*>(rowCopyBase);
#endif

  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const double * element = rowCopy->getElements();
  int scaleLength = ((model->specialOptions()&131072)==0) ? 1 : 2;
  double * rowScale = new double [numberRows*scaleLength];
  double * columnScale = new double [numberColumns*scaleLength];
  // we are going to mark bits we are interested in
  char * usefulRow = new char [numberRows];
  char * usefulColumn = new char [numberColumns];
  double * rowLower = model->rowLower();
  double * rowUpper = model->rowUpper();
  double * columnLower = model->columnLower();
  double * columnUpper = model->columnUpper();
  int iColumn, iRow;
  //#define LEAVE_FIXED
  // mark free rows
  for (iRow=0;iRow<numberRows;iRow++) {
#ifndef LEAVE_FIXED
    usefulRow[iRow]=0;
    if (rowUpper[iRow]<1.0e20||
	rowLower[iRow]>-1.0e20)
#endif
      usefulRow[iRow]=1;
  }
  // mark empty and fixed columns
  // also see if worth scaling
  assert (model->scalingFlag()<4); // dynamic not implemented
  double largest=0.0;
  double smallest=1.0e50;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    char useful=0;
#ifndef LEAVE_FIXED
    if (columnUpper[iColumn]>
	columnLower[iColumn]+1.0e-12) {
#endif
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow=row[j];
	if(elementByColumn[j]&&usefulRow[iRow]) {
	  useful=1;
	  largest = CoinMax(largest,fabs(elementByColumn[j]));
	  smallest = CoinMin(smallest,fabs(elementByColumn[j]));
	}
      }
#ifndef LEAVE_FIXED
    }
#endif
    usefulColumn[iColumn]=useful;
  }
  model->messageHandler()->message(CLP_PACKEDSCALE_INITIAL,*model->messagesPointer())
    <<smallest<<largest
    <<CoinMessageEol;
  if (smallest>=0.5&&largest<=2.0) {
    // don't bother scaling
    model->messageHandler()->message(CLP_PACKEDSCALE_FORGET,*model->messagesPointer())
      <<CoinMessageEol;
    delete [] rowScale;
    delete [] usefulRow;
    delete [] columnScale;
    delete [] usefulColumn;
    if (!model->rowCopy()) 
      delete rowCopyBase; // get rid of temporary row copy
    return 1;
  } else {
      // need to scale 
    int scalingMethod = model->scalingFlag();
    // and see if there any empty rows
    for (iRow=0;iRow<numberRows;iRow++) {
      if (usefulRow[iRow]) {
	CoinBigIndex j;
	int useful=0;
	for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	  int iColumn = column[j];
	  if (usefulColumn[iColumn]) {
	    useful=1;
	    break;
	  }
	}
	usefulRow[iRow]=useful;
      }
    }
    double savedOverallRatio=0.0;
    double tolerance = 5.0*model->primalTolerance();
    double overallLargest=-1.0e-20;
    double overallSmallest=1.0e20;
    bool finished=false;
    // if scalingMethod 3 then may change
    while (!finished) {
      ClpFillN ( rowScale, numberRows,1.0);
      ClpFillN ( columnScale, numberColumns,1.0);
      overallLargest=-1.0e-20;
      overallSmallest=1.0e20;
      if (scalingMethod==1||scalingMethod==3) {
        // Maximum in each row
        for (iRow=0;iRow<numberRows;iRow++) {
          if (usefulRow[iRow]) {
            CoinBigIndex j;
            largest=1.0e-10;
            for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
              int iColumn = column[j];
              if (usefulColumn[iColumn]) {
                double value = fabs(element[j]);
                largest = CoinMax(largest,value);
                assert (largest<1.0e40);
              }
            }
            rowScale[iRow]=1.0/largest;
            overallLargest = CoinMax(overallLargest,largest);
            overallSmallest = CoinMin(overallSmallest,largest);
          }
        }
      } else {
        int numberPass=3;
#ifdef USE_OBJECTIVE
        // This will be used to help get scale factors
        double * objective = new double[numberColumns];
        memcpy(objective,model->costRegion(1),numberColumns*sizeof(double));
        double objScale=1.0;
#endif
        while (numberPass) {
          overallLargest=0.0;
          overallSmallest=1.0e50;
          numberPass--;
          // Geometric mean on row scales
          for (iRow=0;iRow<numberRows;iRow++) {
            if (usefulRow[iRow]) {
              CoinBigIndex j;
              largest=1.0e-20;
              smallest=1.0e50;
              for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
                int iColumn = column[j];
                if (usefulColumn[iColumn]) {
                  double value = fabs(element[j]);
                  // Don't bother with tiny elements
                  if (value>1.0e-20) {
                    value *= columnScale[iColumn];
                    largest = CoinMax(largest,value);
                    smallest = CoinMin(smallest,value);
                  }
                }
              }
              rowScale[iRow]=1.0/sqrt(smallest*largest);
              rowScale[iRow]=CoinMax(1.0e-10,CoinMin(1.0e10,rowScale[iRow]));
              overallLargest = CoinMax(largest*rowScale[iRow],overallLargest);
              overallSmallest = CoinMin(smallest*rowScale[iRow],overallSmallest);
            }
          }
#ifdef USE_OBJECTIVE
          largest=1.0e-20;
          smallest=1.0e50;
          for (iColumn=0;iColumn<numberColumns;iColumn++) {
            if (usefulColumn[iColumn]) {
              double value = fabs(objective[iColumn]);
              // Don't bother with tiny elements
              if (value>1.0e-20) {
                value *= columnScale[iColumn];
                largest = CoinMax(largest,value);
                smallest = CoinMin(smallest,value);
              }
            }
          }
          objScale=1.0/sqrt(smallest*largest);
#endif
          model->messageHandler()->message(CLP_PACKEDSCALE_WHILE,*model->messagesPointer())
            <<overallSmallest
            <<overallLargest
            <<CoinMessageEol;
          // skip last column round
          if (numberPass==1)
            break;
          // Geometric mean on column scales
          for (iColumn=0;iColumn<numberColumns;iColumn++) {
            if (usefulColumn[iColumn]) {
              CoinBigIndex j;
              largest=1.0e-20;
              smallest=1.0e50;
              for (j=columnStart[iColumn];
                   j<columnStart[iColumn]+columnLength[iColumn];j++) {
                iRow=row[j];
                double value = fabs(elementByColumn[j]);
                // Don't bother with tiny elements
                if (value>1.0e-20&&usefulRow[iRow]) {
                  value *= rowScale[iRow];
                  largest = CoinMax(largest,value);
                  smallest = CoinMin(smallest,value);
                }
              }
#ifdef USE_OBJECTIVE
              if (fabs(objective[iColumn])>1.0e-20) {
                double value = fabs(objective[iColumn])*objScale;
                largest = CoinMax(largest,value);
                smallest = CoinMin(smallest,value);
              }
#endif
              columnScale[iColumn]=1.0/sqrt(smallest*largest);
              columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
            }
          }
        }
#ifdef USE_OBJECTIVE
        delete [] objective;
        printf("obj scale %g - use it if you want\n",objScale);
#endif
      }
      // If ranges will make horrid then scale
      for (iRow=0;iRow<numberRows;iRow++) {
        if (usefulRow[iRow]) {
          double difference = rowUpper[iRow]-rowLower[iRow];
          double scaledDifference = difference*rowScale[iRow];
          if (scaledDifference>tolerance&&scaledDifference<1.0e-4) {
            // make gap larger
            rowScale[iRow] *= 1.0e-4/scaledDifference;
            rowScale[iRow]=CoinMax(1.0e-10,CoinMin(1.0e10,rowScale[iRow]));
            //printf("Row %d difference %g scaled diff %g => %g\n",iRow,difference,
            // scaledDifference,difference*rowScale[iRow]);
          }
        }
      }
      // final pass to scale columns so largest is reasonable
      // See what smallest will be if largest is 1.0
      overallSmallest=1.0e50;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
        if (usefulColumn[iColumn]) {
          CoinBigIndex j;
          largest=1.0e-20;
          smallest=1.0e50;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            iRow=row[j];
            if(elementByColumn[j]&&usefulRow[iRow]) {
              double value = fabs(elementByColumn[j]*rowScale[iRow]);
              largest = CoinMax(largest,value);
              smallest = CoinMin(smallest,value);
            }
          }
          if (overallSmallest*largest>smallest)
            overallSmallest = smallest/largest;
        }
      }
      if (scalingMethod==1||scalingMethod==2) {
        finished=true;
      } else if (savedOverallRatio==0.0&&scalingMethod!=4) {
        savedOverallRatio=overallSmallest;
        scalingMethod=4;
      } else {
        assert (scalingMethod==4);
        if (overallSmallest>2.0*savedOverallRatio)
          finished=true; // geometric was better
        else
          scalingMethod=1; // redo equilibrium
#if 0
        if (model->logLevel()>2) {
          if (finished)
            printf("equilibrium ratio %g, geometric ratio %g , geo chosen\n",
                   savedOverallRatio,overallSmallest);
          else
            printf("equilibrium ratio %g, geometric ratio %g , equi chosen\n",
                   savedOverallRatio,overallSmallest);
        }
#endif
      }
    }
    //#define RANDOMIZE
#ifdef RANDOMIZE
    // randomize by up to 10%
    for (iRow=0;iRow<numberRows;iRow++) {
      double value = 0.5-CoinDrand48();//between -0.5 to + 0.5
      rowScale[iRow] *= (1.0+0.1*value);
    }
#endif
    overallLargest=1.0;
    if (overallSmallest<1.0e-1)
      overallLargest = 1.0/sqrt(overallSmallest);
    overallLargest = CoinMin(100.0,overallLargest);
    overallSmallest=1.0e50;
    //printf("scaling %d\n",model->scalingFlag());
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnUpper[iColumn]>
	columnLower[iColumn]+1.0e-12) {
        //if (usefulColumn[iColumn]) {
	CoinBigIndex j;
	largest=1.0e-20;
	smallest=1.0e50;
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  iRow=row[j];
	  if(elementByColumn[j]&&usefulRow[iRow]) {
	    double value = fabs(elementByColumn[j]*rowScale[iRow]);
	    largest = CoinMax(largest,value);
	    smallest = CoinMin(smallest,value);
	  }
	}
	columnScale[iColumn]=overallLargest/largest;
        columnScale[iColumn]=CoinMax(1.0e-10,CoinMin(1.0e10,columnScale[iColumn]));
#ifdef RANDOMIZE
	double value = 0.5-CoinDrand48();//between -0.5 to + 0.5
	columnScale[iColumn] *= (1.0+0.1*value);
#endif
	double difference = columnUpper[iColumn]-columnLower[iColumn];
	if (difference<1.0e-5*columnScale[iColumn]) {
	  // make gap larger
	  columnScale[iColumn] = difference/1.0e-5;
	  //printf("Column %d difference %g scaled diff %g => %g\n",iColumn,difference,
	  // scaledDifference,difference*columnScale[iColumn]);
	}
	overallSmallest = CoinMin(overallSmallest,smallest*columnScale[iColumn]);
      }
    }
    model->messageHandler()->message(CLP_PACKEDSCALE_FINAL,*model->messagesPointer())
      <<overallSmallest
      <<overallLargest
      <<CoinMessageEol;
    delete [] usefulRow;
    delete [] usefulColumn;
#ifndef SLIM_CLP
    // If quadratic then make symmetric
    ClpObjective * obj = model->objectiveAsObject();
#ifndef NO_RTTI
    ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(obj));
#else
    ClpQuadraticObjective * quadraticObj = NULL;
    if (obj->type()==2)
      quadraticObj = (static_cast< ClpQuadraticObjective*>(obj));
#endif
    if (quadraticObj) {
      CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
      int numberXColumns = quadratic->getNumCols();
      if (numberXColumns<numberColumns) {
	// we assume symmetric
	int numberQuadraticColumns=0;
	int i;
	//const int * columnQuadratic = quadratic->getIndices();
	//const int * columnQuadraticStart = quadratic->getVectorStarts();
	const int * columnQuadraticLength = quadratic->getVectorLengths();
	for (i=0;i<numberXColumns;i++) {
	  int length=columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
	  length=1;
#endif
	  if (length)
	    numberQuadraticColumns++;
	}
	int numberXRows = numberRows-numberQuadraticColumns;
	numberQuadraticColumns=0;
	for (i=0;i<numberXColumns;i++) { 
	  int length=columnQuadraticLength[i];
#ifndef CORRECT_COLUMN_COUNTS
	  length=1;
#endif
	  if (length) {
	    rowScale[numberQuadraticColumns+numberXRows] = columnScale[i];
	    numberQuadraticColumns++;
	  }
	}    
	int numberQuadraticRows=0;
	for (i=0;i<numberXRows;i++) {
	  // See if any in row quadratic
	  CoinBigIndex j;
	  int numberQ=0;
	  for (j=rowStart[i];j<rowStart[i+1];j++) {
	    int iColumn = column[j];
	    if (columnQuadraticLength[iColumn])
	      numberQ++;
	  }
#ifndef CORRECT_ROW_COUNTS
	  numberQ=1;
#endif
	  if (numberQ) {
	    columnScale[numberQuadraticRows+numberXColumns] = rowScale[i];
	    numberQuadraticRows++;
	  }
	}
	// and make sure Sj okay
	for (iColumn=numberQuadraticRows+numberXColumns;iColumn<numberColumns;iColumn++) {
	  CoinBigIndex j=columnStart[iColumn];
	  assert(columnLength[iColumn]==1);
	  int iRow=row[j];
	  double value = fabs(elementByColumn[j]*rowScale[iRow]);
	  columnScale[iColumn]=1.0/value;
	}
      }
    }
#endif
    if (scaleLength>1) {
      // make copy
      double * inverseScale = rowScale + numberRows;
      for (iRow=0;iRow<numberRows;iRow++) 
	inverseScale[iRow] = 1.0/rowScale[iRow] ;
      inverseScale = columnScale + numberColumns;
      for (iColumn=0;iColumn<numberColumns;iColumn++) 
	inverseScale[iColumn] = 1.0/columnScale[iColumn] ;
    }
    model->setRowScale(rowScale);
    model->setColumnScale(columnScale);
    if (model->rowCopy()) {
      // need to replace row by row
      double * newElement = new double[numberColumns];
      // scale row copy
      for (iRow=0;iRow<numberRows;iRow++) {
	CoinBigIndex j;
	double scale = rowScale[iRow];
	const double * elementsInThisRow = element + rowStart[iRow];
	const int * columnsInThisRow = column + rowStart[iRow];
	int number = rowStart[iRow+1]-rowStart[iRow];
	assert (number<=numberColumns);
	for (j=0;j<number;j++) {
	  int iColumn = columnsInThisRow[j];
	  newElement[j] = elementsInThisRow[j]*scale*columnScale[iColumn];
	}
	rowCopy->replaceVector(iRow,number,newElement);
      }
      delete [] newElement;
    } else {
      // no row copy
      delete rowCopyBase;
    }
    return 0;
  }
}
/* Unpacks a column into an CoinIndexedvector
 */
void 
ClpPackedMatrix::unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn) const 
{
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      rowArray->add(row[i],elementByColumn[i]);
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn];
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      rowArray->add(iRow,elementByColumn[i]*scale*rowScale[iRow]);
    }
  }
}
/* Unpacks a column into a CoinIndexedVector
** in packed format
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void 
ClpPackedMatrix::unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int iColumn) const
{
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  int * index = rowArray->getIndices();
  double * array = rowArray->denseVector();
  int number = 0;
  if (!rowScale) {
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      double value = elementByColumn[i];
      if (value) {
	array[number]=value;
	index[number++]=iRow;
      }
    }
    rowArray->setNumElements(number);
    rowArray->setPackedMode(true);
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn];
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      double value = elementByColumn[i]*scale*rowScale[iRow];
      if (value) {
	array[number]=value;
	index[number++]=iRow;
      }
    }
    rowArray->setNumElements(number);
    rowArray->setPackedMode(true);
  }
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void 
ClpPackedMatrix::add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int iColumn, double multiplier) const 
{
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      rowArray->quickAdd(iRow,multiplier*elementByColumn[i]);
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn]*multiplier;
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      rowArray->quickAdd(iRow,elementByColumn[i]*scale*rowScale[iRow]);
    }
  }
}
/* Adds multiple of a column into an array */
void 
ClpPackedMatrix::add(const ClpSimplex * model,double * array,
		    int iColumn, double multiplier) const
{
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  CoinBigIndex i;
  if (!rowScale) {
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      array[iRow] += multiplier*elementByColumn[i];
    }
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn]*multiplier;
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      array[iRow] += elementByColumn[i]*scale*rowScale[iRow];
    }
  }
}
/* Checks if all elements are in valid range.  Can just
   return true if you are not paranoid.  For Clp I will
   probably expect no zeros.  Code can modify matrix to get rid of
   small elements.
   check bits (can be turned off to save time) :
   1 - check if matrix has gaps
   2 - check if zero elements
   4 - check and compress duplicates
   8 - report on large and small
*/
bool 
ClpPackedMatrix::allElementsInRange(ClpModel * model,
				    double smallest, double largest,
				    int check)
{
  int iColumn;
  // make sure matrix correct size
  matrix_->setDimensions(model->numberRows(),model->numberColumns());
  CoinBigIndex numberLarge=0;;
  CoinBigIndex numberSmall=0;;
  CoinBigIndex numberDuplicate=0;;
  int firstBadColumn=-1;
  int firstBadRow=-1;
  double firstBadElement=0.0;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  int numberRows = matrix_->getNumRows();
  int numberColumns = matrix_->getNumCols();
  // Say no gaps
  flags_ &= ~2;
  if (check==14||check==10) {
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (end!=columnStart[iColumn+1]) {
	flags_ |= 2;
	break;
      }
    }
    return true;
  }
  assert (check==15||check==11);
  if (check==15) {
    int * mark = new int [numberRows];
    int i;
    for (i=0;i<numberRows;i++)
      mark[i]=-1;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (end!=columnStart[iColumn+1])
	flags_ |= 2;
      for (j=start;j<end;j++) {
	double value = fabs(elementByColumn[j]);
	int iRow = row[j];
	if (iRow<0||iRow>=numberRows) {
#ifndef COIN_BIG_INDEX
	  printf("Out of range %d %d %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#elif COIN_BIG_INDEX==1
	  printf("Out of range %d %ld %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#else
	  printf("Out of range %d %lld %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#endif
	  return false;
	}
	if (mark[iRow]==-1) {
	  mark[iRow]=j;
	} else {
	  // duplicate
	  numberDuplicate++;
	}
	//printf("%d %d %d %g\n",iColumn,j,row[j],elementByColumn[j]);
	if (!value)
	  flags_ |= 1; // there are zero elements
	if (value<smallest) {
	  numberSmall++;
	} else if (!(value<=largest)) {
	  numberLarge++;
	  if (firstBadColumn<0) {
	    firstBadColumn=iColumn;
	    firstBadRow=row[j];
	    firstBadElement=elementByColumn[j];
	  }
	}
      }
      //clear mark
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int iRow = row[j];
	mark[iRow]=-1;
      }
    }
    delete [] mark;
  } else {
    // just check for out of range - not for duplicates
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (end!=columnStart[iColumn+1])
	flags_ |= 2;
      for (j=start;j<end;j++) {
	double value = fabs(elementByColumn[j]);
	int iRow = row[j];
	if (iRow<0||iRow>=numberRows) {
#ifndef COIN_BIG_INDEX
	  printf("Out of range %d %d %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#elif COIN_BIG_INDEX==1
	  printf("Out of range %d %ld %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#else
	  printf("Out of range %d %lld %d %g\n",iColumn,j,row[j],elementByColumn[j]);
#endif
	  return false;
	}
	if (!value)
	  flags_ |= 1; // there are zero elements
	if (value<smallest) {
	  numberSmall++;
	} else if (!(value<=largest)) {
	  numberLarge++;
	  if (firstBadColumn<0) {
	    firstBadColumn=iColumn;
	    firstBadRow=iRow;
	    firstBadElement=value;
	  }
	}
      }
    }
  }
  if (numberLarge) {
    model->messageHandler()->message(CLP_BAD_MATRIX,model->messages())
      <<numberLarge
      <<firstBadColumn<<firstBadRow<<firstBadElement
      <<CoinMessageEol;
    return false;
  }
  if (numberSmall) 
    model->messageHandler()->message(CLP_SMALLELEMENTS,model->messages())
      <<numberSmall
      <<CoinMessageEol;
  if (numberDuplicate) 
    model->messageHandler()->message(CLP_DUPLICATEELEMENTS,model->messages())
      <<numberDuplicate
      <<CoinMessageEol;
  if (numberDuplicate) 
    matrix_->eliminateDuplicates(smallest);
  else if (numberSmall) 
    matrix_->compress(smallest);
  // If smallest >0.0 then there can't be zero elements
  if (smallest>0.0)
    flags_ &= ~1;;
  return true;
}
/* Given positive integer weights for each row fills in sum of weights
   for each column (and slack).
   Returns weights vector
*/
CoinBigIndex * 
ClpPackedMatrix::dubiousWeights(const ClpSimplex * model,int * inputWeights) const
{
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  int number = numberRows+numberColumns;
  CoinBigIndex * weights = new CoinBigIndex[number];
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  int i;
  for (i=0;i<numberColumns;i++) {
    CoinBigIndex j;
    CoinBigIndex count=0;
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      int iRow=row[j];
      count += inputWeights[iRow];
    }
    weights[i]=count;
  }
  for (i=0;i<numberRows;i++) {
    weights[i+numberColumns]=inputWeights[i];
  }
  return weights;
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void 
ClpPackedMatrix::rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive)
{
  smallestNegative=-COIN_DBL_MAX;
  largestNegative=0.0;
  smallestPositive=COIN_DBL_MAX;
  largestPositive=0.0;
  // get matrix data pointers
  const double * elementByColumn = matrix_->getElements();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  int numberColumns = matrix_->getNumCols();
  int i;
  for (i=0;i<numberColumns;i++) {
    CoinBigIndex j;
    for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
      double value = elementByColumn[j];
      if (value>0.0) {
	smallestPositive = CoinMin(smallestPositive,value);
	largestPositive = CoinMax(largestPositive,value);
      } else if (value<0.0) {
	smallestNegative = CoinMax(smallestNegative,value);
	largestNegative = CoinMin(largestNegative,value);
      }
    }
  }
}
// Says whether it can do partial pricing
bool 
ClpPackedMatrix::canDoPartialPricing() const
{
  return true;
}
// Partial pricing 
void 
ClpPackedMatrix::partialPricing(ClpSimplex * model, double startFraction, double endFraction,
			      int & bestSequence, int & numberWanted)
{
  numberWanted=currentWanted_;
  int start = (int) (startFraction*numberActiveColumns_);
  int end = CoinMin((int) (endFraction*numberActiveColumns_+1),numberActiveColumns_);
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
  if (bestSequence>=0)
    bestDj = fabs(model->clpMatrix()->reducedCost(model,bestSequence));
  else
    bestDj=tolerance;
  int sequenceOut = model->sequenceOut();
  int saveSequence = bestSequence;
  int lastScan = minimumObjectsScan_<0 ? end : start+minimumObjectsScan_; 
  int minNeg = minimumGoodReducedCosts_==-1 ? numberWanted : minimumGoodReducedCosts_;
  if (rowScale) {
    // scaled
    for (iSequence=start;iSequence<end;iSequence++) {
      if (iSequence!=sequenceOut) {
	double value;
	ClpSimplex::Status status = model->getStatus(iSequence);
	
	switch(status) {
	  
	case ClpSimplex::basic:
	case ClpSimplex::isFixed:
	  break;
	case ClpSimplex::isFree:
	case ClpSimplex::superBasic:
	  value=0.0;
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
	  value=0.0;
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
	  value=0.0;
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
      if (numberWanted+minNeg<originalWanted_&&iSequence>lastScan) {
	// give up
	break;
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
      savedBestSequence_ = bestSequence;
      savedBestDj_ = reducedCost[savedBestSequence_];
    }
  } else {
    // not scaled
    for (iSequence=start;iSequence<end;iSequence++) {
      if (iSequence!=sequenceOut) {
	double value;
	ClpSimplex::Status status = model->getStatus(iSequence);
	
	switch(status) {
	  
	case ClpSimplex::basic:
	case ClpSimplex::isFixed:
	  break;
	case ClpSimplex::isFree:
	case ClpSimplex::superBasic:
	  value=cost[iSequence];
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
	  value=cost[iSequence];
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
	  value=cost[iSequence];
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
      if (numberWanted+minNeg<originalWanted_&&iSequence>lastScan) {
	// give up
	break;
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
      savedBestSequence_ = bestSequence;
      savedBestDj_ = reducedCost[savedBestSequence_];
    }
  }
  currentWanted_=numberWanted;
}
// Sets up an effective RHS 
void 
ClpPackedMatrix::useEffectiveRhs(ClpSimplex * model)
{
  delete [] rhsOffset_;
  int numberRows = model->numberRows();
  rhsOffset_ = new double[numberRows];
  rhsOffset(model,true);
}
// Gets rid of special copies
void 
ClpPackedMatrix::clearCopies()
{
  delete rowCopy_;
  delete columnCopy_;
  rowCopy_=NULL;
  columnCopy_=NULL;
  flags_ &= ~(4+8);
 }
// makes sure active columns correct
int 
ClpPackedMatrix::refresh(ClpSimplex * model)
{
  numberActiveColumns_ = matrix_->getNumCols();
#if 0
  ClpMatrixBase * rowCopyBase = reverseOrderedCopy();
  ClpPackedMatrix* rowCopy =
    dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
  // Make sure it is really a ClpPackedMatrix
  assert (rowCopy!=NULL);

  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths();
  const double * element = rowCopy->getElements();
  for (int i=0;i<rowCopy->getNumRows();i++) {
    if (!rowLength[i])
      printf("zero row %d\n",i);
  }
  delete rowCopy;
#endif
  return 0;
}

/* Scales rowCopy if column copy scaled
   Only called if scales already exist */
void 
ClpPackedMatrix::scaleRowCopy(ClpModel * model) const 
{
  if (model->rowCopy()) {
    // need to replace row by row
    int numberRows = model->numberRows();
    int numberColumns = matrix_->getNumCols();
    double * newElement = new double[numberColumns];
    ClpMatrixBase * rowCopyBase=model->rowCopy();
#ifndef NO_RTTI
    ClpPackedMatrix* rowCopy =
      dynamic_cast< ClpPackedMatrix*>(rowCopyBase);
    // Make sure it is really a ClpPackedMatrix
    assert (rowCopy!=NULL);
#else
    ClpPackedMatrix* rowCopy =
      static_cast< ClpPackedMatrix*>(rowCopyBase);
#endif

    const int * column = rowCopy->getIndices();
    const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
    const double * element = rowCopy->getElements();
    const double * rowScale = model->rowScale();
    const double * columnScale = model->columnScale();
    // scale row copy
    for (int iRow=0;iRow<numberRows;iRow++) {
      CoinBigIndex j;
      double scale = rowScale[iRow];
      const double * elementsInThisRow = element + rowStart[iRow];
      const int * columnsInThisRow = column + rowStart[iRow];
      int number = rowStart[iRow+1]-rowStart[iRow];
      assert (number<=numberColumns);
      for (j=0;j<number;j++) {
	int iColumn = columnsInThisRow[j];
	newElement[j] = elementsInThisRow[j]*scale*columnScale[iColumn];
      }
      rowCopy->replaceVector(iRow,number,newElement);
    }
    delete [] newElement;
  }
}
/* Realy really scales column copy 
   Only called if scales already exist.
   Up to user ro delete */
ClpMatrixBase * 
ClpPackedMatrix::scaledColumnCopy(ClpModel * model) const 
{
  // need to replace column by column
  int numberRows = model->numberRows();
  int numberColumns = matrix_->getNumCols();
  double * newElement = new double[numberRows];
  ClpPackedMatrix * copy = new ClpPackedMatrix(*this);
  const int * row = copy->getIndices();
  const CoinBigIndex * columnStart = copy->getVectorStarts();
  const int * length = copy->getVectorLengths();
  const double * element = copy->getElements();
  const double * rowScale = model->rowScale();
  const double * columnScale = model->columnScale();
  // scale column copy
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double scale = columnScale[iColumn];
    const double * elementsInThisColumn = element + columnStart[iColumn];
    const int * rowsInThisColumn = row + columnStart[iColumn];
    int number = length[iColumn];
    assert (number<=numberRows);
    for (j=0;j<number;j++) {
      int iRow = rowsInThisColumn[j];
      newElement[j] = elementsInThisColumn[j]*scale*rowScale[iRow];
    }
    copy->replaceVector(iColumn,number,newElement);
  }
  delete [] newElement;
  return copy;
}
// Really scale matrix
void 
ClpPackedMatrix::reallyScale(const double * rowScale, const double * columnScale)
{
  clearCopies();
  int numberColumns = matrix_->getNumCols();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * length = matrix_->getVectorLengths();
  double * element = matrix_->getMutableElements();
  // scale 
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double scale = columnScale[iColumn];
    for (j=columnStart[iColumn];j<columnStart[iColumn]+length[iColumn];j++) {
      int iRow = row[j];
      element[j] *= scale*rowScale[iRow];
    }
  }
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void 
ClpPackedMatrix::deleteCols(const int numDel, const int * indDel)
{ 
  if (matrix_->getNumCols())
    matrix_->deleteCols(numDel,indDel);
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  flags_ |= 2;
  matrix_->setExtraGap(1.0e-50);
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void 
ClpPackedMatrix::deleteRows(const int numDel, const int * indDel)
{ 
  if (matrix_->getNumRows())
    matrix_->deleteRows(numDel,indDel);
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  flags_ |= 2;
  matrix_->setExtraGap(1.0e-50);
}
#ifndef CLP_NO_VECTOR
// Append Columns
void 
ClpPackedMatrix::appendCols(int number, const CoinPackedVectorBase * const * columns)
{ 
  matrix_->appendCols(number,columns);
  numberActiveColumns_ = matrix_->getNumCols();
  clearCopies();
}
// Append Rows
void 
ClpPackedMatrix::appendRows(int number, const CoinPackedVectorBase * const * rows)
{ 
  matrix_->appendRows(number,rows);
  numberActiveColumns_ = matrix_->getNumCols();
  // may now have gaps
  flags_ |= 2;
  clearCopies();
}
#endif
/* Set the dimensions of the matrix. In effect, append new empty
   columns/rows to the matrix. A negative number for either dimension
   means that that dimension doesn't change. Otherwise the new dimensions
   MUST be at least as large as the current ones otherwise an exception
   is thrown. */
void 
ClpPackedMatrix::setDimensions(int numrows, int numcols)
{
  matrix_->setDimensions(numrows,numcols);
}
/* Append a set of rows/columns to the end of the matrix. Returns number of errors
   i.e. if any of the new rows/columns contain an index that's larger than the
   number of columns-1/rows-1 (if numberOther>0) or duplicates
   If 0 then rows, 1 if columns */
int 
ClpPackedMatrix::appendMatrix(int number, int type,
                            const CoinBigIndex * starts, const int * index,
                            const double * element, int numberOther)
{
  int numberErrors=0;
  // make sure other dimension is big enough
  if (type==0) {
    // rows
    if (matrix_->isColOrdered()&&numberOther>matrix_->getNumCols())
      matrix_->setDimensions(-1,numberOther);
    numberErrors=matrix_->appendRows(number,starts,index,element,numberOther);
  } else {
    // columns
    if (!matrix_->isColOrdered()&&numberOther>matrix_->getNumRows())
      matrix_->setDimensions(numberOther,-1);
    numberErrors=matrix_->appendCols(number,starts,index,element,numberOther);
  }
  clearCopies();
  numberActiveColumns_ = matrix_->getNumCols();
  return numberErrors;
}
void
ClpPackedMatrix::specialRowCopy(ClpSimplex * model,const ClpMatrixBase * rowCopy)
{
  delete rowCopy_;
  rowCopy_ = new ClpPackedMatrix2(model,rowCopy->getPackedMatrix());
  // See if anything in it
  if (!rowCopy_->usefulInfo()) {
    delete rowCopy_;
    rowCopy_=NULL;
    flags_ &= ~4;
  } else {
    flags_ |= 4;
  }
}
void
ClpPackedMatrix::specialColumnCopy(ClpSimplex * model)
{
  delete columnCopy_;
  if ((flags_&8)!=0)
    columnCopy_ = new ClpPackedMatrix3(model,matrix_);
  else
    columnCopy_=NULL;
}
// Say we don't want special column copy
void 
ClpPackedMatrix::releaseSpecialColumnCopy()
{ 
  flags_ &= ~8;
  delete columnCopy_;
  columnCopy_=NULL;
}
// Correct sequence in and out to give true value
void 
ClpPackedMatrix::correctSequence(const ClpSimplex * model,int & sequenceIn, int & sequenceOut) 
{
  if (columnCopy_) {
    if (sequenceIn!=-999) {
      if (sequenceIn!=sequenceOut) {
	if (sequenceIn<numberActiveColumns_)
	  columnCopy_->swapOne(model,this,sequenceIn);
	if (sequenceOut<numberActiveColumns_)
	  columnCopy_->swapOne(model,this,sequenceOut);
      }
    } else {
      // do all
      columnCopy_->sortBlocks(model);
    }
  }
}
// Check validity
void 
ClpPackedMatrix::checkFlags() const
{
  int iColumn;
  // get matrix data pointers
  //const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  if (!zeros()) {
    for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
      CoinBigIndex j;
      for (j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	if (!elementByColumn[j])
	  abort();
      }
    }
  }
  if ((flags_&2)==0) {
    for (iColumn=0;iColumn<numberActiveColumns_;iColumn++) {
      if (columnStart[iColumn+1]!=columnStart[iColumn]+columnLength[iColumn]) {
	abort();
      }
    }
  }
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2 () 
  : numberBlocks_(0),
    numberRows_(0),
    offset_(NULL),
    count_(NULL),
    rowStart_(NULL),
    column_(NULL),
    work_(NULL)
{
#ifdef THREAD
  threadId_ = NULL;
  info_ = NULL;
#endif
}
//-------------------------------------------------------------------
// Useful Constructor 
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2 (ClpSimplex * model,const CoinPackedMatrix * rowCopy) 
  : numberBlocks_(0),
    numberRows_(0),
    offset_(NULL),
    count_(NULL),
    rowStart_(NULL),
    column_(NULL),
    work_(NULL)
{
#ifdef THREAD
  threadId_ = NULL;
  info_ = NULL;
#endif
  numberRows_ = rowCopy->getNumRows();
  if (!numberRows_)
    return;
  int numberColumns = rowCopy->getNumCols();
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * length = rowCopy->getVectorLengths();
  const double * element = rowCopy->getElements();
  int chunk = 32768; // tune
  //chunk=100;
  // tune
#if 0
  int chunkY[7]={1024,2048,4096,8192,16384,32768,65535};
  int its = model->maximumIterations();
  if (its>=1000000&&its<1000999) {
    its -= 1000000;
    its = its/10;
    if (its>=7) abort();
    chunk=chunkY[its];
    printf("chunk size %d\n",chunk);
#define cpuid(func,ax,bx,cx,dx)\
    __asm__ __volatile__ ("cpuid":\
    "=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func));
    unsigned int a,b,c,d;
    int func=0;
    cpuid(func,a,b,c,d);
    {
      int i;
      unsigned int value;
      value=b;
      for (i=0;i<4;i++) {
        printf("%c",(value&0xff));
        value = value >>8;
      }
      value=d;
      for (i=0;i<4;i++) {
        printf("%c",(value&0xff));
        value = value >>8;
      }
      value=c;
      for (i=0;i<4;i++) {
        printf("%c",(value&0xff));
        value = value >>8;
      }
      printf("\n");
      int maxfunc = a;
      if (maxfunc>10) {
        printf("not intel?\n");
        abort();
      }
      for (func=1;func<=maxfunc;func++) {
        cpuid(func,a,b,c,d);
        printf("func %d, %x %x %x %x\n",func,a,b,c,d);
      }
    }
#else
    if (numberColumns>10000||chunk==100) {
#endif
  } else {
    //printf("no chunk\n");
    return;
  }
  // Could also analyze matrix to get natural breaks
  numberBlocks_ = (numberColumns+chunk-1)/chunk;
#ifdef THREAD
  // Get work areas
  threadId_ = new pthread_t [numberBlocks_];
  info_ = new dualColumn0Struct[numberBlocks_];
#endif
  // Even out
  chunk = (numberColumns+numberBlocks_-1)/numberBlocks_;
  offset_ = new int[numberBlocks_+1];
  offset_[numberBlocks_]=numberColumns;
  int nRow = numberBlocks_*numberRows_;
  count_ = new unsigned short[nRow];
  memset(count_,0,nRow*sizeof(unsigned short));
  rowStart_ = new CoinBigIndex[nRow+numberRows_+1];
  CoinBigIndex nElement = rowStart[numberRows_];
  rowStart_[nRow+numberRows_]=nElement;
  column_ = new unsigned short [nElement];
  // assumes int <= double
  int sizeWork = 6*numberBlocks_;
  work_ = new double[sizeWork];;
  int iBlock;
  int nZero=0;
  for (iBlock=0;iBlock<numberBlocks_;iBlock++) {
    int start = iBlock*chunk;
    offset_[iBlock]=start;
    int end = start+chunk;
    for (int iRow=0;iRow<numberRows_;iRow++) {
      if (rowStart[iRow+1]!=rowStart[iRow]+length[iRow]) {
        printf("not packed correctly - gaps\n");
        abort();
      }
      bool lastFound=false;
      int nFound=0;
      for (CoinBigIndex j = rowStart[iRow];
           j<rowStart[iRow]+length[iRow];j++) {
        int iColumn = column[j];
        if (iColumn>=start) {
          if (iColumn<end) {
            if (!element[j]) {
              printf("not packed correctly - zero element\n");
              abort();
            }
            column_[j]=iColumn-start;
            nFound++;
            if (lastFound) {
              printf("not packed correctly - out of order\n");
              abort();
            }
          } else {
            //can't find any more
            lastFound=true;
          }
        } 
      }
      count_[iRow*numberBlocks_+iBlock]=nFound;
      if (!nFound)
        nZero++;
    }
  }
  //double fraction = ((double) nZero)/((double) (numberBlocks_*numberRows_));
  //printf("%d empty blocks, %g%%\n",nZero,100.0*fraction);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpPackedMatrix2::ClpPackedMatrix2 (const ClpPackedMatrix2 & rhs) 
  : numberBlocks_(rhs.numberBlocks_),
    numberRows_(rhs.numberRows_)
{  
  if (numberBlocks_) {
    offset_ = CoinCopyOfArray(rhs.offset_,numberBlocks_+1);
    int nRow = numberBlocks_*numberRows_;
    count_ = CoinCopyOfArray(rhs.count_,nRow);
    rowStart_ = CoinCopyOfArray(rhs.rowStart_,nRow+numberRows_+1);
    CoinBigIndex nElement = rowStart_[nRow+numberRows_];
    column_ = CoinCopyOfArray(rhs.column_,nElement);
    int sizeWork = 6*numberBlocks_;
    work_ = CoinCopyOfArray(rhs.work_,sizeWork);
#ifdef THREAD
    threadId_ = new pthread_t [numberBlocks_];
    info_ = new dualColumn0Struct[numberBlocks_];
#endif
  } else {
    offset_ = NULL;
    count_ = NULL;
    rowStart_ = NULL;
    column_ = NULL;
    work_ = NULL;
#ifdef THREAD
    threadId_ = NULL;
    info_ = NULL;
#endif
  }
}
//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPackedMatrix2::~ClpPackedMatrix2 ()
{
  delete [] offset_;
  delete [] count_;
  delete [] rowStart_;
  delete [] column_;
  delete [] work_;
#ifdef THREAD
  delete [] threadId_;
  delete [] info_;
#endif
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpPackedMatrix2 &
ClpPackedMatrix2::operator=(const ClpPackedMatrix2& rhs)
{
  if (this != &rhs) {
    numberBlocks_ = rhs.numberBlocks_;
    numberRows_ = rhs.numberRows_;
    delete [] offset_;
    delete [] count_;
    delete [] rowStart_;
    delete [] column_;
    delete [] work_;
#ifdef THREAD
    delete [] threadId_;
    delete [] info_;
#endif
    if (numberBlocks_) {
      offset_ = CoinCopyOfArray(rhs.offset_,numberBlocks_+1);
      int nRow = numberBlocks_*numberRows_;
      count_ = CoinCopyOfArray(rhs.count_,nRow);
      rowStart_ = CoinCopyOfArray(rhs.rowStart_,nRow+numberRows_+1);
      CoinBigIndex nElement = rowStart_[nRow+numberRows_];
      column_ = CoinCopyOfArray(rhs.column_,nElement);
      int sizeWork = 6*numberBlocks_;
      work_ = CoinCopyOfArray(rhs.work_,sizeWork);
#ifdef THREAD
      threadId_ = new pthread_t [numberBlocks_];
      info_ = new dualColumn0Struct[numberBlocks_];
#endif
    } else {
      offset_ = NULL;
      count_ = NULL;
      rowStart_ = NULL;
      column_ = NULL;
      work_ = NULL;
#ifdef THREAD
      threadId_ = NULL;
      info_ = NULL;
#endif
    }
  }
  return *this;
}
static int dualColumn0(const ClpSimplex * model,double * spare,
                       int * spareIndex,const double * arrayTemp,
                       const int * indexTemp,int numberIn,
                       int offset, double acceptablePivot,double * bestPossiblePtr,
                       double * upperThetaPtr,int * posFreePtr, double * freePivotPtr)
{
  // do dualColumn0
  int i;
  int numberRemaining=0;
  double bestPossible=0.0;
  double upperTheta=1.0e31;
  double freePivot = acceptablePivot;
  int posFree=-1;
  const double * reducedCost=model->djRegion(1);
  double dualTolerance = model->dualTolerance();
  // We can also see if infeasible or pivoting on free
  double tentativeTheta = 1.0e25;
  for (i=0;i<numberIn;i++) {
    double alpha = arrayTemp[i];
    int iSequence = indexTemp[i]+offset;
    double oldValue;
    double value;
    bool keep;
    
    switch(model->getStatus(iSequence)) {
      
    case ClpSimplex::basic:
    case ClpSimplex::isFixed:
      break;
    case ClpSimplex::isFree:
    case ClpSimplex::superBasic:
      bestPossible = CoinMax(bestPossible,fabs(alpha));
      oldValue = reducedCost[iSequence];
      // If free has to be very large - should come in via dualRow
      if (model->getStatus(iSequence)==ClpSimplex::isFree&&fabs(alpha)<1.0e-3)
        break;
      if (oldValue>dualTolerance) {
        keep = true;
      } else if (oldValue<-dualTolerance) {
        keep = true;
      } else {
        if (fabs(alpha)>CoinMax(10.0*acceptablePivot,1.0e-5)) 
          keep = true;
        else
          keep=false;
      }
      if (keep) {
        // free - choose largest
        if (fabs(alpha)>freePivot) {
          freePivot=fabs(alpha);
          posFree = i;
        }
      }
      break;
    case ClpSimplex::atUpperBound:
      oldValue = reducedCost[iSequence];
      value = oldValue-tentativeTheta*alpha;
      //assert (oldValue<=dualTolerance*1.0001);
      if (value>dualTolerance) {
        bestPossible = CoinMax(bestPossible,-alpha);
        value = oldValue-upperTheta*alpha;
        if (value>dualTolerance && -alpha>=acceptablePivot) 
          upperTheta = (oldValue-dualTolerance)/alpha;
        // add to list
        spare[numberRemaining]=alpha;
        spareIndex[numberRemaining++]=iSequence;
      }
      break;
    case ClpSimplex::atLowerBound:
      oldValue = reducedCost[iSequence];
      value = oldValue-tentativeTheta*alpha;
      //assert (oldValue>=-dualTolerance*1.0001);
      if (value<-dualTolerance) {
        bestPossible = CoinMax(bestPossible,alpha);
        value = oldValue-upperTheta*alpha;
        if (value<-dualTolerance && alpha>=acceptablePivot) 
          upperTheta = (oldValue+dualTolerance)/alpha;
        // add to list
        spare[numberRemaining]=alpha;
        spareIndex[numberRemaining++]=iSequence;
      }
      break;
    }
  }
  *bestPossiblePtr = bestPossible;
  *upperThetaPtr = upperTheta;
  *freePivotPtr = freePivot;
  *posFreePtr = posFree;
  return numberRemaining;
}
static int doOneBlock(double * array, int * index, 
                      const double * pi,const CoinBigIndex * rowStart, const double * element,
                      const unsigned short * column, int numberInRowArray, int numberLook)
{
  int iWhich=0;
  int nextN=0;
  CoinBigIndex nextStart=0;
  double nextPi=0.0;
  for (;iWhich<numberInRowArray;iWhich++) {
    nextStart = rowStart[0];
    nextN = rowStart[numberInRowArray]-nextStart;
    rowStart++;
    if (nextN) {
      nextPi = pi[iWhich];
      break;
    }
  }
  int i;
  while (iWhich<numberInRowArray) {
    double value = nextPi;
    CoinBigIndex j=  nextStart;
    int n = nextN;
    // get next
    iWhich++;
    for (;iWhich<numberInRowArray;iWhich++) {
      nextStart = rowStart[0];
      nextN = rowStart[numberInRowArray]-nextStart;
      rowStart++;
      if (nextN) {
        coin_prefetch(element+nextStart);
        nextPi = pi[iWhich];
        break;
      }
    }
    CoinBigIndex end =j+n;
    //coin_prefetch(element+rowStart_[i+1]);
    //coin_prefetch(column_+rowStart_[i+1]);
    if (n<100) {
      if ((n&1)!=0) {
        unsigned int jColumn = column[j];
        array[jColumn] -= value*element[j];
        j++;
      }      
      coin_prefetch(column+nextStart);
      for (;j<end;j+=2) {
        unsigned int jColumn0 = column[j];
        double value0 = value*element[j];
        unsigned int jColumn1 = column[j+1];
        double value1 = value*element[j+1];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
      }
    } else {
      if ((n&1)!=0) {
        unsigned int jColumn = column[j];
        array[jColumn] -= value*element[j];
        j++;
      }      
      if ((n&2)!=0) {
        unsigned int jColumn0 = column[j];
        double value0 = value*element[j];
        unsigned int jColumn1 = column[j+1];
        double value1 = value*element[j+1];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        j+=2;
      }      
      if ((n&4)!=0) {
        unsigned int jColumn0 = column[j];
        double value0 = value*element[j];
        unsigned int jColumn1 = column[j+1];
        double value1 = value*element[j+1];
        unsigned int jColumn2 = column[j+2];
        double value2 = value*element[j+2];
        unsigned int jColumn3 = column[j+3];
        double value3 = value*element[j+3];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
        j+=4;
      }
      //coin_prefetch(column+nextStart);
      for (;j<end;j+=8) {
        coin_prefetch(element+j+16);
        unsigned int jColumn0 = column[j];
        double value0 = value*element[j];
        unsigned int jColumn1 = column[j+1];
        double value1 = value*element[j+1];
        unsigned int jColumn2 = column[j+2];
        double value2 = value*element[j+2];
        unsigned int jColumn3 = column[j+3];
        double value3 = value*element[j+3];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
        coin_prefetch(column+j+16);
        jColumn0 = column[j+4];
        value0 = value*element[j+4];
        jColumn1 = column[j+5];
        value1 = value*element[j+5];
        jColumn2 = column[j+6];
        value2 = value*element[j+6];
        jColumn3 = column[j+7];
        value3 = value*element[j+7];
        array[jColumn0] -= value0;
        array[jColumn1] -= value1;
        array[jColumn2] -= value2;
        array[jColumn3] -= value3;
      }
    }
  }
  // get rid of tiny values
  int nSmall = numberLook;
  int numberNonZero=0;
  for (i=0;i<nSmall;i++) {
    double value = array[i];
    array[i]=0.0; 
    if (fabs(value)>1.0e-12) {
      array[numberNonZero]=value;
      index[numberNonZero++]=i;
    }
  }
  for (;i<numberLook;i+=4) {
    double value0 = array[i+0];
    double value1 = array[i+1];
    double value2 = array[i+2];
    double value3 = array[i+3];
    array[i+0]=0.0; 
    array[i+1]=0.0; 
    array[i+2]=0.0; 
    array[i+3]=0.0; 
    if (fabs(value0)>1.0e-12) {
      array[numberNonZero]=value0;
      index[numberNonZero++]=i+0;
    }
    if (fabs(value1)>1.0e-12) {
      array[numberNonZero]=value1;
      index[numberNonZero++]=i+1;
    }
    if (fabs(value2)>1.0e-12) {
      array[numberNonZero]=value2;
      index[numberNonZero++]=i+2;
    }
    if (fabs(value3)>1.0e-12) {
      array[numberNonZero]=value3;
      index[numberNonZero++]=i+3;
    }
  }
  return numberNonZero;
}
#ifdef THREAD
static void * doOneBlockThread(void * voidInfo)
{
  dualColumn0Struct * info = (dualColumn0Struct *) voidInfo;
  *(info->numberInPtr) =  doOneBlock(info->arrayTemp,info->indexTemp,info->pi,
                                  info->rowStart,info->element,info->column,
                                  info->numberInRowArray,info->numberLook);
  return NULL;
}
static void * doOneBlockAnd0Thread(void * voidInfo)
{
  dualColumn0Struct * info = (dualColumn0Struct *) voidInfo;
  *(info->numberInPtr) =  doOneBlock(info->arrayTemp,info->indexTemp,info->pi,
                                  info->rowStart,info->element,info->column,
                                  info->numberInRowArray,info->numberLook);
  *(info->numberOutPtr) =  dualColumn0(info->model,info->spare,
                                    info->spareIndex,(const double *)info->arrayTemp,
                                    (const int *) info->indexTemp,*(info->numberInPtr),
                                    info->offset, info->acceptablePivot,info->bestPossiblePtr,
                                    info->upperThetaPtr,info->posFreePtr, info->freePivotPtr);
  return NULL;
}
#endif
/* Return <code>x * scalar * A in <code>z</code>. 
   Note - x packed and z will be packed mode
   Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix2::transposeTimes(const ClpSimplex * model, 
                                 const CoinPackedMatrix * rowCopy,
                                 const CoinIndexedVector * rowArray,
                                 CoinIndexedVector * spareArray,
                                 CoinIndexedVector * columnArray) const
{
  // See if dualColumn0 coding wanted
  bool dualColumn = model->spareIntArray_[0]==1;
  double acceptablePivot = model->spareDoubleArray_[0];
  double bestPossible=0.0;
  double upperTheta=1.0e31;
  double freePivot = acceptablePivot;
  int posFree=-1;
  int numberRemaining=0;
  //if (model->numberIterations()>=200000) {
  //printf("time %g\n",CoinCpuTime()-startTime);
  //exit(77);
  //}
  double * pi = rowArray->denseVector();
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  const int * whichRow = rowArray->getIndices();
  double * element = const_cast<double *>(rowCopy->getElements());
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  int i;
  CoinBigIndex * rowStart2 = rowStart_;
  if (!dualColumn) {
    for (i=0;i<numberInRowArray;i++) {
      int iRow = whichRow[i]; 
      CoinBigIndex start = rowStart[iRow];
      *rowStart2 = start;
      unsigned short * count1 = count_+iRow*numberBlocks_;
      int put=0;
      for (int j=0;j<numberBlocks_;j++) {
        put += numberInRowArray;
        start += count1[j];
        rowStart2[put]=start;
      }
      rowStart2 ++;
    }
  } else {
    // also do dualColumn stuff
    double * spare = spareArray->denseVector();
    int * spareIndex = spareArray->getIndices();
    const double * reducedCost=model->djRegion(0);
    double dualTolerance = model->dualTolerance();
    // We can also see if infeasible or pivoting on free
    double tentativeTheta = 1.0e25;
    int addSequence=model->numberColumns();
    for (i=0;i<numberInRowArray;i++) {
      int iRow = whichRow[i]; 
      double alpha=pi[i];
      double oldValue;
      double value;
      bool keep;

      switch(model->getStatus(iRow+addSequence)) {
	  
      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
	break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        bestPossible = CoinMax(bestPossible,fabs(alpha));
	oldValue = reducedCost[iRow];
        // If free has to be very large - should come in via dualRow
        if (model->getStatus(iRow+addSequence)==ClpSimplex::isFree&&fabs(alpha)<1.0e-3)
          break;
	if (oldValue>dualTolerance) {
	  keep = true;
	} else if (oldValue<-dualTolerance) {
	  keep = true;
	} else {
	  if (fabs(alpha)>CoinMax(10.0*acceptablePivot,1.0e-5)) 
	    keep = true;
	  else
	    keep=false;
	}
	if (keep) {
	  // free - choose largest
	  if (fabs(alpha)>freePivot) {
	    freePivot=fabs(alpha);
	    posFree = i + addSequence;
	  }
	}
	break;
      case ClpSimplex::atUpperBound:
	oldValue = reducedCost[iRow];
	value = oldValue-tentativeTheta*alpha;
	//assert (oldValue<=dualTolerance*1.0001);
	if (value>dualTolerance) {
          bestPossible = CoinMax(bestPossible,-alpha);
	  value = oldValue-upperTheta*alpha;
	  if (value>dualTolerance && -alpha>=acceptablePivot) 
	    upperTheta = (oldValue-dualTolerance)/alpha;
	  // add to list
	  spare[numberRemaining]=alpha;
	  spareIndex[numberRemaining++]=iRow+addSequence;
	}
	break;
      case ClpSimplex::atLowerBound:
	oldValue = reducedCost[iRow];
	value = oldValue-tentativeTheta*alpha;
	//assert (oldValue>=-dualTolerance*1.0001);
	if (value<-dualTolerance) {
          bestPossible = CoinMax(bestPossible,alpha);
	  value = oldValue-upperTheta*alpha;
	  if (value<-dualTolerance && alpha>=acceptablePivot) 
	    upperTheta = (oldValue+dualTolerance)/alpha;
	  // add to list
	  spare[numberRemaining]=alpha;
	  spareIndex[numberRemaining++]=iRow+addSequence;
	}
	break;
      }
      CoinBigIndex start = rowStart[iRow];
      *rowStart2 = start;
      unsigned short * count1 = count_+iRow*numberBlocks_;
      int put=0;
      for (int j=0;j<numberBlocks_;j++) {
        put += numberInRowArray;
        start += count1[j];
        rowStart2[put]=start;
      }
      rowStart2 ++;
    }
  }
  double * spare = spareArray->denseVector();
  int * spareIndex = spareArray->getIndices();
  int saveNumberRemaining=numberRemaining;
  int iBlock;
  for (iBlock=0;iBlock<numberBlocks_;iBlock++) {
    double * dwork = work_+6*iBlock;
    int * iwork = (int *) (dwork+3);
    if (!dualColumn) {
#ifndef THREAD
      int offset = offset_[iBlock];
      int offset3 = offset;
      offset = numberNonZero;
      double * arrayTemp = array+offset;
      int * indexTemp = index+offset; 
      iwork[0] = doOneBlock(arrayTemp,indexTemp,pi,rowStart_+numberInRowArray*iBlock,
                            element,column_,numberInRowArray,offset_[iBlock+1]-offset);
      int number = iwork[0];
      for (i=0;i<number;i++) {
        //double value = arrayTemp[i];
        //arrayTemp[i]=0.0; 
        //array[numberNonZero]=value;
        index[numberNonZero++]=indexTemp[i]+offset3;
      }
#else
      int offset = offset_[iBlock];
      double * arrayTemp = array+offset;
      int * indexTemp = index+offset; 
      dualColumn0Struct * infoPtr = info_+iBlock;
      infoPtr->arrayTemp=arrayTemp;
      infoPtr->indexTemp=indexTemp;
      infoPtr->numberInPtr=&iwork[0];
      infoPtr->pi=pi;
      infoPtr->rowStart=rowStart_+numberInRowArray*iBlock;
      infoPtr->element=element;
      infoPtr->column=column_;
      infoPtr->numberInRowArray=numberInRowArray;
      infoPtr->numberLook=offset_[iBlock+1]-offset;
      pthread_create(&threadId_[iBlock],NULL,doOneBlockThread,infoPtr);
#endif
    } else {
#ifndef THREAD
      int offset = offset_[iBlock];
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      int offset3 = offset;
      offset = numberNonZero;
      offset2 = numberRemaining;
      double * arrayTemp = array+offset;
      int * indexTemp = index+offset; 
      iwork[0] = doOneBlock(arrayTemp,indexTemp,pi,rowStart_+numberInRowArray*iBlock,
                            element,column_,numberInRowArray,offset_[iBlock+1]-offset);
      iwork[1] = dualColumn0(model,spare+offset2,
                                   spareIndex+offset2,
                                   arrayTemp,indexTemp,
                                   iwork[0],offset3,acceptablePivot,
                                   &dwork[0],&dwork[1],&iwork[2],
                                   &dwork[2]);
      int number = iwork[0];
      int numberLook = iwork[1];
#if 1
      numberRemaining += numberLook;
#else
      double * spareTemp = spare+offset2;
      const int * spareIndexTemp = spareIndex+offset2; 
      for (i=0;i<numberLook;i++) {
        double value = spareTemp[i];
        spareTemp[i]=0.0;
        spare[numberRemaining] = value;
        spareIndex[numberRemaining++] = spareIndexTemp[i];
      }
#endif
      if (dwork[2]>freePivot) {
        freePivot=dwork[2];
        posFree = iwork[2]+numberNonZero;
      } 
      upperTheta =  CoinMin(dwork[1],upperTheta);
      bestPossible = CoinMax(dwork[0],bestPossible);
      for (i=0;i<number;i++) {
        // double value = arrayTemp[i];
        //arrayTemp[i]=0.0; 
        //array[numberNonZero]=value;
        index[numberNonZero++]=indexTemp[i]+offset3;
      }
#else
      int offset = offset_[iBlock];
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      double * arrayTemp = array+offset;
      int * indexTemp = index+offset; 
      dualColumn0Struct * infoPtr = info_+iBlock;
      infoPtr->model=model;
      infoPtr->spare=spare+offset2;
      infoPtr->spareIndex=spareIndex+offset2;
      infoPtr->arrayTemp=arrayTemp;
      infoPtr->indexTemp=indexTemp;
      infoPtr->numberInPtr=&iwork[0];
      infoPtr->offset=offset;
      infoPtr->acceptablePivot=acceptablePivot;
      infoPtr->bestPossiblePtr=&dwork[0];
      infoPtr->upperThetaPtr=&dwork[1];
      infoPtr->posFreePtr=&iwork[2];
      infoPtr->freePivotPtr=&dwork[2];
      infoPtr->numberOutPtr=&iwork[1];
      infoPtr->pi=pi;
      infoPtr->rowStart=rowStart_+numberInRowArray*iBlock;
      infoPtr->element=element;
      infoPtr->column=column_;
      infoPtr->numberInRowArray=numberInRowArray;
      infoPtr->numberLook=offset_[iBlock+1]-offset;
      if (iBlock>=2)
        pthread_join(threadId_[iBlock-2],NULL);
      pthread_create(threadId_+iBlock,NULL,doOneBlockAnd0Thread,infoPtr);
      //pthread_join(threadId_[iBlock],NULL);
#endif
    }
  }
  for ( iBlock=CoinMax(0,numberBlocks_-2);iBlock<numberBlocks_;iBlock++) {
#ifdef THREAD
    pthread_join(threadId_[iBlock],NULL);
#endif
  }
#ifdef THREAD
  for ( iBlock=0;iBlock<numberBlocks_;iBlock++) {
    //pthread_join(threadId_[iBlock],NULL);
    int offset = offset_[iBlock];
    double * dwork = work_+6*iBlock;
    int * iwork = (int *) (dwork+3);
    int number = iwork[0];
    if (dualColumn) {
      // allow for already saved
      int offset2 = offset + saveNumberRemaining;
      int numberLook = iwork[1];
      double * spareTemp = spare+offset2;
      const int * spareIndexTemp = spareIndex+offset2; 
      for (i=0;i<numberLook;i++) {
        double value = spareTemp[i];
        spareTemp[i]=0.0;
        spare[numberRemaining] = value;
        spareIndex[numberRemaining++] = spareIndexTemp[i];
      }
      if (dwork[2]>freePivot) {
        freePivot=dwork[2];
        posFree = iwork[2]+numberNonZero;
      } 
      upperTheta =  CoinMin(dwork[1],upperTheta);
      bestPossible = CoinMax(dwork[0],bestPossible);
    }
    double * arrayTemp = array+offset;
    const int * indexTemp = index+offset; 
    for (i=0;i<number;i++) {
      double value = arrayTemp[i];
      arrayTemp[i]=0.0; 
      array[numberNonZero]=value;
      index[numberNonZero++]=indexTemp[i]+offset;
    }
  }
#endif
  columnArray->setNumElements(numberNonZero);
  columnArray->setPackedMode(true);
  if (dualColumn) {
    model->spareDoubleArray_[0]=upperTheta;
    model->spareDoubleArray_[1]=bestPossible;
    // and theta and alpha and sequence
    if (posFree<0) {
      model->spareIntArray_[1]=-1;
    } else {
      const double * reducedCost=model->djRegion(0);
      double alpha;
      int numberColumns=model->numberColumns();
      if (posFree<numberColumns) {
        alpha = columnArray->denseVector()[posFree];
        posFree = columnArray->getIndices()[posFree];
      } else {
        alpha = rowArray->denseVector()[posFree-numberColumns];
        posFree = rowArray->getIndices()[posFree-numberColumns]+numberColumns;
      }
      model->spareDoubleArray_[2]=fabs(reducedCost[posFree]/alpha);;
      model->spareDoubleArray_[3]=alpha;
      model->spareIntArray_[1]=posFree;
    }
    spareArray->setNumElements(numberRemaining);
    // signal done
    model->spareIntArray_[0]=-1;
  }
}
/* Default constructor. */
ClpPackedMatrix3::ClpPackedMatrix3()
  : numberBlocks_(0),
    numberColumns_(0),
    column_(NULL),
    start_(NULL),
    row_(NULL),
    element_(NULL),
    block_(NULL)
{
}
/* Constructor from copy. */
ClpPackedMatrix3::ClpPackedMatrix3(ClpSimplex * model,const CoinPackedMatrix * columnCopy)
  : numberBlocks_(0),
    numberColumns_(0),
    column_(NULL),
    start_(NULL),
    row_(NULL),
    element_(NULL),
    block_(NULL)
{
#define MINBLOCK 6
#define MAXBLOCK 100
#define MAXUNROLL 10
  numberColumns_ = columnCopy->getNumCols();
  int numberRows = columnCopy->getNumRows();
  int * counts = new int[numberRows+1];
  CoinZeroN(counts,numberRows+1);
  CoinBigIndex nels=0;
  CoinBigIndex nZeroEl=0;
  int iColumn;
  // get matrix data pointers
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * elementByColumn = columnCopy->getElements();
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    int n = columnLength[iColumn];
    CoinBigIndex end = start + n;
    int kZero=0;
    for (CoinBigIndex j=start;j<end;j++) {
      if(!elementByColumn[j])
	kZero++;
    }
    n -= kZero;
    nZeroEl += kZero;
    nels += n;
    counts[n]++;
  }
  int nZeroColumns = counts[0];
  counts[0]=-1;
  numberColumns_ -= nZeroColumns;
  column_ = new int [2*numberColumns_];
  int * lookup = column_ + numberColumns_;
  row_ = new int[nels];
  element_ = new double[nels];
  int nOdd=0;
  CoinBigIndex nInOdd=0;
  int i;
  for (i=1;i<=numberRows;i++) {
    int n = counts[i];
    if (n) {
      if (n<MINBLOCK||i>MAXBLOCK) {
	nOdd += n;
	counts[i]=-1;
	nInOdd += n*i;
      } else {
	numberBlocks_++;
      }
    } else {
      counts[i]=-1;
    }
  }
  start_ = new CoinBigIndex [nOdd+1];
  // even if no blocks do a dummy one
  numberBlocks_ = CoinMax(numberBlocks_,1);
  block_ = (blockStruct *) new char [numberBlocks_*sizeof(blockStruct)];
  memset(block_,0,numberBlocks_*sizeof(blockStruct));
  // Fill in what we can
  int nTotal=nOdd;
  // in case no blocks
  block_->startIndices_=nTotal;
  nels=nInOdd;
  int nBlock=0;
  for (i=0;i<=CoinMin(MAXBLOCK,numberRows);i++) {
    if (counts[i]>0) {
      blockStruct * block = block_ + nBlock;
      int n=counts[i];
      counts[i]= nBlock; // backward pointer
      nBlock++;
      block->startIndices_ = nTotal;
      block->startElements_ = nels;
      block->numberElements_ = i;
      // up counts
      nTotal += n;
      nels += n*i;
    }
  }
  // fill
  start_[0]=0;
  nOdd=0;
  nInOdd=0;
  const double * columnScale = model->columnScale();
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    int n = columnLength[iColumn];
    CoinBigIndex end = start + n;
    int kZero=0;
    for (CoinBigIndex j=start;j<end;j++) {
      if(!elementByColumn[j])
	kZero++;
    }
    n -= kZero;
    if (n) {
      int iBlock = counts[n];
      if (iBlock>=0) {
	blockStruct * block = block_ + iBlock;
	int k = block->numberInBlock_;
	block->numberInBlock_ ++;
	column_[block->startIndices_+k]=iColumn;
	lookup[iColumn]=k;
	CoinBigIndex put = block->startElements_+k*n;
	for (CoinBigIndex j=start;j<end;j++) {
	  double value = elementByColumn[j];
	  if(value) {
	    if (columnScale)
	      value *= columnScale[iColumn];
	    element_[put] = value;
	    row_[put++]=row[j];
	  }
	}
      } else {
	// odd ones
	for (CoinBigIndex j=start;j<end;j++) {
	  double value = elementByColumn[j];
	  if(value) {
	    if (columnScale)
	      value *= columnScale[iColumn];
	    element_[nInOdd] = value;
	    row_[nInOdd++]=row[j];
	  }
	}
	column_[nOdd]=iColumn;
	lookup[iColumn]=-1;
	nOdd++;
	start_[nOdd]=nInOdd;
      }
    }
  }
  delete [] counts;
}
/* Destructor */
ClpPackedMatrix3::~ClpPackedMatrix3()
{
  delete [] column_;
  delete [] start_;
  delete [] row_;
  delete [] element_;
  delete [] block_;
}
/* The copy constructor. */
ClpPackedMatrix3::ClpPackedMatrix3(const ClpPackedMatrix3 & rhs)
  : numberBlocks_(rhs.numberBlocks_),
    numberColumns_(rhs.numberColumns_),
    column_(NULL),
    start_(NULL),
    row_(NULL),
    element_(NULL),
    block_(NULL)
{
  if (rhs.numberBlocks_) {
    block_ = CoinCopyOfArray(rhs.block_,numberBlocks_);
    column_ = CoinCopyOfArray(rhs.column_,2*numberColumns_);
    int numberOdd = block_->startIndices_;
    start_ = CoinCopyOfArray(rhs.start_,numberOdd+1);
    blockStruct * lastBlock = block_ + (numberBlocks_-1);
    CoinBigIndex numberElements = lastBlock->startElements_+
      lastBlock->numberInBlock_*lastBlock->numberElements_;
    row_ = CoinCopyOfArray(rhs.row_,numberElements);
    element_ = CoinCopyOfArray(rhs.element_,numberElements);
  }
}
ClpPackedMatrix3& 
ClpPackedMatrix3::operator=(const ClpPackedMatrix3 & rhs)
{
  if (this != &rhs) {
    delete [] column_;
    delete [] start_;
    delete [] row_;
    delete [] element_;
    delete [] block_;
    numberBlocks_ = rhs.numberBlocks_;
    numberColumns_ = rhs.numberColumns_;
    if (rhs.numberBlocks_) {
      block_ = CoinCopyOfArray(rhs.block_,numberBlocks_);
      column_ = CoinCopyOfArray(rhs.column_,2*numberColumns_);
      int numberOdd = block_->startIndices_;
      start_ = CoinCopyOfArray(rhs.start_,numberOdd+1);
      blockStruct * lastBlock = block_ + (numberBlocks_-1);
      CoinBigIndex numberElements = lastBlock->startElements_+
	lastBlock->numberInBlock_*lastBlock->numberElements_;
      row_ = CoinCopyOfArray(rhs.row_,numberElements);
      element_ = CoinCopyOfArray(rhs.element_,numberElements);
    } else {
      column_ = NULL;
      start_ = NULL;
      row_ = NULL;
      element_ = NULL;
      block_ = NULL;
    }
  }
  return *this;
}
/* Sort blocks */
void 
ClpPackedMatrix3::sortBlocks(const ClpSimplex * model)
{
  int * lookup = column_ + numberColumns_;
  for (int iBlock=0;iBlock<numberBlocks_;iBlock++) {
    blockStruct * block = block_ + iBlock;
    int numberInBlock = block->numberInBlock_; 
    int nel = block->numberElements_;
    int * row = row_ + block->startElements_;
    double * element = element_ + block->startElements_;
    int * column = column_ + block->startIndices_;
    int lastPrice=0;
    int firstNotPrice=numberInBlock-1;
    while (lastPrice<=firstNotPrice) {
      // find first basic or fixed
      int iColumn=numberInBlock;
      for (;lastPrice<=firstNotPrice;lastPrice++) {
	iColumn = column[lastPrice];
	if (model->getColumnStatus(iColumn)==ClpSimplex::basic||
	    model->getColumnStatus(iColumn)==ClpSimplex::isFixed)
	  break;
      }
      // find last non basic or fixed
      int jColumn=-1;
      for (;firstNotPrice>lastPrice;firstNotPrice--) {
	jColumn = column[firstNotPrice];
	if (model->getColumnStatus(jColumn)!=ClpSimplex::basic&&
	    model->getColumnStatus(jColumn)!=ClpSimplex::isFixed)
	  break;
      }
      if (firstNotPrice>lastPrice) {
	assert (column[lastPrice]==iColumn);
	assert (column[firstNotPrice]==jColumn);
	// need to swap
	column[firstNotPrice]=iColumn;
	lookup[iColumn]=firstNotPrice;
	column[lastPrice]=jColumn;
	lookup[jColumn]=lastPrice;
	double * elementA = element + lastPrice*nel;
	int * rowA = row + lastPrice*nel;
	double * elementB = element + firstNotPrice*nel;
	int * rowB = row + firstNotPrice*nel;
	for (int i=0;i<nel;i++) {
	  int temp = rowA[i];
	  double tempE = elementA[i];
	  rowA[i]=rowB[i];
	  elementA[i]=elementB[i];
	  rowB[i]=temp;
	  elementB[i]=tempE;
	}
	firstNotPrice--;
	lastPrice++;
      } else if (lastPrice==firstNotPrice) {
	// make sure correct side
	iColumn = column[lastPrice];
	if (model->getColumnStatus(iColumn)!=ClpSimplex::basic&&
	    model->getColumnStatus(iColumn)!=ClpSimplex::isFixed)
	  lastPrice++;
	break;
      }
    }
    block->numberPrice_=lastPrice;
#ifndef NDEBUG
    // check
    int i;
    for (i=0;i<lastPrice;i++) {
      int iColumn = column[i];
      assert (model->getColumnStatus(iColumn)!=ClpSimplex::basic&&
	      model->getColumnStatus(iColumn)!=ClpSimplex::isFixed);
      assert (lookup[iColumn]==i);
    }
    for (;i<numberInBlock;i++) {
      int iColumn = column[i];
      assert (model->getColumnStatus(iColumn)==ClpSimplex::basic||
	      model->getColumnStatus(iColumn)==ClpSimplex::isFixed);
      assert (lookup[iColumn]==i);
    }
#endif
  }
}
// Swap one variable
void 
ClpPackedMatrix3::swapOne(const ClpSimplex * model,const ClpPackedMatrix * matrix,
			  int iColumn)
{
  int * lookup = column_ + numberColumns_;
  // position in block
  int kA = lookup[iColumn];
  if (kA<0)
    return; // odd one
  // get matrix data pointers
  const CoinPackedMatrix * columnCopy = matrix->getPackedMatrix();
  //const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * elementByColumn = columnCopy->getElements();
  CoinBigIndex start = columnStart[iColumn];
  int n = columnLength[iColumn];
  if (matrix->zeros()) {
    CoinBigIndex end = start + n;
    for (CoinBigIndex j=start;j<end;j++) {
      if(!elementByColumn[j])
	n--;
    }
  }
  // find block - could do binary search
  int iBlock=CoinMin(n,numberBlocks_)-1;
  while (block_[iBlock].numberElements_!=n)
    iBlock--;
  blockStruct * block = block_ + iBlock;
  int nel = block->numberElements_;
  int * row = row_ + block->startElements_;
  double * element = element_ + block->startElements_;
  int * column = column_ + block->startIndices_;
  assert (column[kA]==iColumn);
  bool moveUp = (model->getColumnStatus(iColumn)==ClpSimplex::basic||
		 model->getColumnStatus(iColumn)==ClpSimplex::isFixed);
  int lastPrice=block->numberPrice_;
  int kB;
  if (moveUp) {
    kB=lastPrice-1;
    block->numberPrice_--;
  } else {
    kB=lastPrice;
    block->numberPrice_++;
  }
  int jColumn = column[kB];
  column[kA]=jColumn;
  lookup[jColumn]=kA;
  column[kB]=iColumn;
  lookup[iColumn]=kB;
  double * elementA = element + kB*nel;
  int * rowA = row + kB*nel;
  double * elementB = element + kA*nel;
  int * rowB = row + kA*nel;
  int i;
  for (i=0;i<nel;i++) {
    int temp = rowA[i];
    double tempE = elementA[i];
    rowA[i]=rowB[i];
    elementA[i]=elementB[i];
    rowB[i]=temp;
    elementB[i]=tempE;
  }
#if 1
#ifndef NDEBUG
  // check
  for (i=0;i<block->numberPrice_;i++) {
    int iColumn = column[i];
    if (iColumn!=model->sequenceIn()&&iColumn!=model->sequenceOut())
      assert (model->getColumnStatus(iColumn)!=ClpSimplex::basic&&
	      model->getColumnStatus(iColumn)!=ClpSimplex::isFixed);
    assert (lookup[iColumn]==i);
  }
  int numberInBlock = block->numberInBlock_; 
  for (;i<numberInBlock;i++) {
    int iColumn = column[i];
    if (iColumn!=model->sequenceIn()&&iColumn!=model->sequenceOut())
      assert (model->getColumnStatus(iColumn)==ClpSimplex::basic||
	      model->getColumnStatus(iColumn)==ClpSimplex::isFixed);
    assert (lookup[iColumn]==i);
  }
#endif
#endif
}
/* Return <code>x * -1 * A in <code>z</code>. 
   Note - x packed and z will be packed mode
   Squashes small elements and knows about ClpSimplex */
void 
ClpPackedMatrix3::transposeTimes(const ClpSimplex * model,
				 const double * pi,
				 CoinIndexedVector * output) const
{
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  double zeroTolerance = model->factorization()->zeroTolerance();
  double value = 0.0;
  CoinBigIndex j;
  int numberOdd = block_->startIndices_;
  if (numberOdd) {
    // A) as probably long may be worth unrolling
    CoinBigIndex end = start_[1];
    for (j=start_[0]; j<end;j++) {
      int iRow = row_[j];
      value += pi[iRow]*element_[j];
    }
    int iColumn;
    // int jColumn=column_[0];
    
    for (iColumn=0;iColumn<numberOdd-1;iColumn++) {
      CoinBigIndex start = end;
      end = start_[iColumn+2];
      if (fabs(value)>zeroTolerance) {
	array[numberNonZero]=value;
	index[numberNonZero++]=column_[iColumn];
	//index[numberNonZero++]=jColumn;
      }
      // jColumn = column_[iColumn+1];
      value = 0.0;
      //if (model->getColumnStatus(jColumn)!=ClpSimplex::basic) {
	for (j=start; j<end;j++) {
	  int iRow = row_[j];
	  value += pi[iRow]*element_[j];
	}
	//}
    }
    if (fabs(value)>zeroTolerance) {
      array[numberNonZero]=value;
      index[numberNonZero++]=column_[iColumn];
      //index[numberNonZero++]=jColumn;
    }
  }
  for (int iBlock=0;iBlock<numberBlocks_;iBlock++) {
    // B) Can sort so just do nonbasic (and nonfixed)
    // C) Can do two at a time (if so put odd one into start_)
    // D) can use switch
    blockStruct * block = block_ + iBlock;
    //int numberPrice = block->numberInBlock_; 
    int numberPrice = block->numberPrice_; 
    int nel = block->numberElements_;
    int * row = row_ + block->startElements_;
    double * element = element_ + block->startElements_;
    int * column = column_ + block->startIndices_;
#if 0
    // two at a time
    if ((numberPrice&1)!=0) {
      double value = 0.0;
      int nel2=nel;
      for (; nel2;nel2--) {
	int iRow = *row++;
	value += pi[iRow]*(*element++);
      }
      if (fabs(value)>zeroTolerance) {
	array[numberNonZero]=value;
	index[numberNonZero++]=*column;
      }
      column++;
    }
    numberPrice = numberPrice>>1;
    switch ((nel%2)) {
      // 2 k +0
    case 0:
      for (;numberPrice;numberPrice--) {
	double value0 = 0.0;
	double value1 = 0.0;
	int nel2=nel;
	for (; nel2;nel2--) {
	  int iRow0 = *row;
	  int iRow1 = *(row+nel);
	  row++;
	  double element0 = *element;
	  double element1 = *(element+nel);
	  element++;
	  value0 += pi[iRow0]*element0;
	  value1 += pi[iRow1]*element1;
	}
	row += nel;
	element += nel;
	if (fabs(value0)>zeroTolerance) {
	  array[numberNonZero]=value0;
	  index[numberNonZero++]=*column;
	}
	column++;
	if (fabs(value1)>zeroTolerance) {
	  array[numberNonZero]=value1;
	  index[numberNonZero++]=*column;
	}
	column++;
      }  
      break;
      // 2 k +1
    case 1:
      for (;numberPrice;numberPrice--) {
	double value0;
	double value1;
	int nel2=nel-1;
	{
	  int iRow0 = row[0];
	  int iRow1 = row[nel];
	  double pi0 = pi[iRow0];
	  double pi1 = pi[iRow1];
	  value0 = pi0*element[0];
	  value1 = pi1*element[nel];
	  row++;
	  element++;
	}
	for (; nel2;nel2--) {
	  int iRow0 = *row;
	  int iRow1 = *(row+nel);
	  row++;
	  double element0 = *element;
	  double element1 = *(element+nel);
	  element++;
	  value0 += pi[iRow0]*element0;
	  value1 += pi[iRow1]*element1;
	}
	row += nel;
	element += nel;
	if (fabs(value0)>zeroTolerance) {
	  array[numberNonZero]=value0;
	  index[numberNonZero++]=*column;
	}
	column++;
	if (fabs(value1)>zeroTolerance) {
	  array[numberNonZero]=value1;
	  index[numberNonZero++]=*column;
	}
	column++;
      }
      break;
    }
#else
    for (;numberPrice;numberPrice--) {
      double value = 0.0;
      int nel2=nel;
      for (; nel2;nel2--) {
	int iRow = *row++;
	value += pi[iRow]*(*element++);
      }
      if (fabs(value)>zeroTolerance) {
	array[numberNonZero]=value;
	index[numberNonZero++]=*column;
      }
      column++;
    }
#endif
  }
  output->setNumElements(numberNonZero);
}
// Updates two arrays for steepest 
void 
ClpPackedMatrix3::transposeTimes2(const ClpSimplex * model,
				  const double * pi, CoinIndexedVector * output,
				  const double * piWeight,
				  double referenceIn, double devex,
				  // Array for exact devex to say what is in reference framework
				  unsigned int * reference,
				  double * weights, double scaleFactor)
{
  int numberNonZero=0;
  int * index = output->getIndices();
  double * array = output->denseVector();
  double zeroTolerance = model->factorization()->zeroTolerance();
  double value = 0.0;
  bool killDjs = (scaleFactor==0.0);
  if (!scaleFactor)
    scaleFactor=1.0;
  int numberOdd = block_->startIndices_;
  int iColumn;
  CoinBigIndex end = start_[0];
  for (iColumn=0;iColumn<numberOdd;iColumn++) {
    CoinBigIndex start = end;
    CoinBigIndex j;
    int jColumn = column_[iColumn];
    end = start_[iColumn+1];
    value = 0.0;
    if (model->getColumnStatus(jColumn)!=ClpSimplex::basic) {
      for (j=start; j<end;j++) {
	int iRow = row_[j];
	value -= pi[iRow]*element_[j];
      }
      if (fabs(value)>zeroTolerance) {
	// and do other array
	double modification = 0.0;
	for (j=start; j<end;j++) {
	  int iRow = row_[j];
	  modification += piWeight[iRow]*element_[j];
	}
	double thisWeight = weights[jColumn];
	double pivot = value*scaleFactor;
	double pivotSquared = pivot * pivot;
	thisWeight += pivotSquared * devex + pivot * modification;
	if (thisWeight<DEVEX_TRY_NORM) {
	  if (referenceIn<0.0) {
	    // steepest
	    thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
	  } else {
	    // exact
	    thisWeight = referenceIn*pivotSquared;
	    if (reference(jColumn))
	      thisWeight += 1.0;
	    thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
	  }
	}
	weights[jColumn] = thisWeight;
	if (!killDjs) {
	  array[numberNonZero]=value;
	  index[numberNonZero++]=jColumn;
	}
      }
    }
  }
  for (int iBlock=0;iBlock<numberBlocks_;iBlock++) {
    // B) Can sort so just do nonbasic (and nonfixed)
    // C) Can do two at a time (if so put odd one into start_)
    // D) can use switch
    blockStruct * block = block_ + iBlock;
    //int numberPrice = block->numberInBlock_; 
    int numberPrice = block->numberPrice_; 
    int nel = block->numberElements_;
    int * row = row_ + block->startElements_;
    double * element = element_ + block->startElements_;
    int * column = column_ + block->startIndices_;
    for (;numberPrice;numberPrice--) {
      double value = 0.0;
      int nel2=nel;
      for (; nel2;nel2--) {
	int iRow = *row++;
	value -= pi[iRow]*(*element++);
      }
      if (fabs(value)>zeroTolerance) {
	int jColumn = *column;
	// back to beginning
	row -= nel;
	element -= nel;
	// and do other array
	double modification = 0.0;
	nel2=nel;
	for (; nel2;nel2--) {
	  int iRow = *row++;
	  modification += piWeight[iRow]*(*element++);
	}
	double thisWeight = weights[jColumn];
	double pivot = value*scaleFactor;
	double pivotSquared = pivot * pivot;
	thisWeight += pivotSquared * devex + pivot * modification;
	if (thisWeight<DEVEX_TRY_NORM) {
	  if (referenceIn<0.0) {
	    // steepest
	    thisWeight = CoinMax(DEVEX_TRY_NORM,DEVEX_ADD_ONE+pivotSquared);
	  } else {
	    // exact
	    thisWeight = referenceIn*pivotSquared;
	    if (reference(jColumn))
	      thisWeight += 1.0;
	    thisWeight = CoinMax(thisWeight,DEVEX_TRY_NORM);
	  }
	}
	weights[jColumn] = thisWeight;
	if (!killDjs) {
	  array[numberNonZero]=value;
	  index[numberNonZero++]=jColumn;
	}
      }
      column++;
    }
  }
  output->setNumElements(numberNonZero);
  output->setPackedMode(true);
}
#ifdef CLP_ALL_ONE_FILE
#undef reference
#endif
