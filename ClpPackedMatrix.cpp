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
#include "ClpPackedMatrix.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix () 
  : ClpMatrixBase(),
    matrix_(NULL),
    zeroElements_(false)
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
  zeroElements_ = rhs.zeroElements_;
  
}

//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix (CoinPackedMatrix * rhs) 
: ClpMatrixBase()
{  
  matrix_ = rhs;
  zeroElements_ = false;
  setType(1);
  
}

ClpPackedMatrix::ClpPackedMatrix (const CoinPackedMatrix & rhs) 
: ClpMatrixBase()
{  
  matrix_ = new CoinPackedMatrix(rhs);
  zeroElements_ = false;
  setType(1);
  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpPackedMatrix::~ClpPackedMatrix ()
{
  delete matrix_;
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
    zeroElements_ = rhs.zeroElements_;
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
  zeroElements_ = rhs.zeroElements_;
}
ClpPackedMatrix::ClpPackedMatrix (
		       const CoinPackedMatrix & rhs,
		       int numberRows, const int * whichRows,
		       int numberColumns, const int * whichColumns)
: ClpMatrixBase()
{
  matrix_ = new CoinPackedMatrix(rhs,numberRows,whichRows,
				 numberColumns,whichColumns);
  zeroElements_ = false;
  setType(1);
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase * 
ClpPackedMatrix::reverseOrderedCopy() const
{
  ClpPackedMatrix * copy = new ClpPackedMatrix();
  copy->matrix_= new CoinPackedMatrix();
  copy->matrix_->reverseOrderedCopyOf(*matrix_);
  copy->matrix_->removeGaps();
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
  int numberColumns = matrix_->getNumCols();
  //memset(y,0,matrix_->getNumRows()*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
  int numberColumns = matrix_->getNumCols();
  //memset(y,0,numberColumns*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
    int numberColumns = matrix_->getNumCols();
    //memset(y,0,matrix_->getNumRows()*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      double value = x[iColumn];
      if (value) {
	// scaled
	value *= scalar*columnScale[iColumn];
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  iRow=row[j];
	  y[iRow] += value*elementByColumn[j];
	}
      }
    }
    int numberRows = getNumRows();
    for (iRow=0;iRow<numberRows;iRow++) {
      y[iRow] *= rowScale[iRow];
    }
  } else {
    times(scalar,x,y);
  }
}
void 
ClpPackedMatrix::transposeTimes( double scalar,
				 const double * x, double * y,
				 const double * rowScale, 
				 const double * columnScale) const
{
  if (rowScale) {
    int iColumn;
    // get matrix data pointers
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths(); 
    const double * elementByColumn = matrix_->getElements();
    int numberColumns = matrix_->getNumCols();
    //memset(y,0,numberColumns*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
  ClpPackedMatrix* rowCopy =
    dynamic_cast< ClpPackedMatrix*>(model->rowCopy());
  bool packed = rowArray->packedMode();
  double factor = 0.3;
  // We may not want to do by row if there may be cache problems
  int numberColumns = model->numberColumns();
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberColumns*sizeof(double)>1000000) {
    if (numberRows*10<numberColumns)
      factor=0.1;
    else if (numberRows*4<numberColumns)
      factor=0.15;
    else if (numberRows*2<numberColumns)
      factor=0.2;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  if (numberInRowArray>factor*numberRows||!rowCopy) {
    // do by column
    int iColumn;
    // get matrix data pointers
    const int * row = matrix_->getIndices();
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths(); 
    const double * elementByColumn = matrix_->getElements();
    const double * rowScale = model->rowScale();
    int numberColumns = model->numberColumns();
    if (!y->getNumElements()) {
      if (packed) {
	// need to expand pi into y
	assert(y->capacity()>=numberRows);
	double * piOld = pi;
	pi = y->denseVector();
	const int * whichRow = rowArray->getIndices();
	int i;
	if (!rowScale) {
	  // modify pi so can collapse to one loop
	  for (i=0;i<numberInRowArray;i++) {
	    int iRow = whichRow[i];
	    pi[iRow]=scalar*piOld[i];
	  }
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
	  for (i=0;i<numberInRowArray;i++) {
	    int iRow = whichRow[i];
	    pi[iRow]=scalar*piOld[i]*rowScale[iRow];
	  }
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
	  } else if (scalar==1.0) {
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
	      double value = 0.0;
	      CoinBigIndex j;
	      for (j=columnStart[iColumn];
		   j<columnStart[iColumn]+columnLength[iColumn];j++) {
		int iRow = row[j];
		value += pi[iRow]*elementByColumn[j];
	      }
	      if (fabs(value)>zeroTolerance) {
		index[numberNonZero++]=iColumn;
		array[iColumn]=value;
	      }
	    }
	  } else {
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
	  } else if (scalar==1.0) {
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
		array[iColumn]=value;
	      }
	    }
	  } else {
	    for (iColumn=0;iColumn<numberColumns;iColumn++) {
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
    } else {
      assert(!packed);
      double * markVector = y->denseVector(); // not empty
      if (!rowScale) {
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  double value = markVector[iColumn];
	  markVector[iColumn]=0.0;
	  double value2 = 0.0;
	  CoinBigIndex j;
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    value2 += pi[iRow]*elementByColumn[j];
	  }
	  value += value2*scalar;
	  if (fabs(value)>zeroTolerance) {
	    index[numberNonZero++]=iColumn;
	    array[iColumn]=value;
	  }
	}
      } else {
	// scaled
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  double value = markVector[iColumn];
	  markVector[iColumn]=0.0;
	  CoinBigIndex j;
	  const double * columnScale = model->columnScale();
	  double value2 = 0.0;
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    int iRow = row[j];
	    value2 += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	  }
	  value += value2*scalar*columnScale[iColumn];
	  if (fabs(value)>zeroTolerance) {
	    index[numberNonZero++]=iColumn;
	    array[iColumn]=value;
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
  if (numberInRowArray>2||y->getNumElements()) {
    // do by rows
    // ** Row copy is already scaled
    int iRow;
    int * mark = y->getIndices();
    int numberOriginal=y->getNumElements();
    int i;
    if (packed) {
      assert(!numberOriginal);
      numberNonZero=0;
      // and set up mark as char array
      char * marked = (char *) (index+columnArray->capacity());
      double * array2 = y->denseVector();
#ifdef CLP_DEBUG
      int numberColumns = model->numberColumns();
      for (i=0;i<numberColumns;i++) {
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
    } else {
      double * markVector = y->denseVector(); // probably empty .. but
      for (i=0;i<numberOriginal;i++) {
	int iColumn = mark[i];
	index[i]=iColumn;
	array[iColumn]=markVector[iColumn];
	markVector[iColumn]=0.0;
      }
      numberNonZero=numberOriginal;
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
      double value = pi[0]*scalar;
      for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	int iColumn = column[j];
	double value2 = value*element[j];
	if (fabs(value2)>zeroTolerance) {
	  array[numberNonZero] = value2;
	  index[numberNonZero++]=iColumn;
	}
      }
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
  int numberNonZero=0;
  int * index = columnArray->getIndices();
  double * array = columnArray->denseVector();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->factorization()->zeroTolerance();
  int jColumn;
  // get matrix data pointers
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  const double * rowScale = model->rowScale();
  int numberToDo = y->getNumElements();
  const int * which = y->getIndices();
  bool packed = rowArray->packedMode();
  if (packed) {
    // need to expand pi into y
    int numberInRowArray = rowArray->getNumElements();
    assert(y->capacity()>=model->numberRows());
    double * piOld = pi;
    pi = y->denseVector();
    const int * whichRow = rowArray->getIndices();
    int i;
    // Do NOT squash small elements - must line up with  y
    if (!rowScale) {
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=piOld[i];
      }
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
      for (i=0;i<numberInRowArray;i++) {
	int iRow = whichRow[i];
	pi[iRow]=rowScale[iRow]*piOld[i];
      }
      for (jColumn=0;jColumn<numberToDo;jColumn++) {
	int iColumn = which[jColumn];
	double value = 0.0;
	CoinBigIndex j;
	const double * columnScale = model->columnScale();
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  value += pi[iRow]*elementByColumn[j]*rowScale[iRow];
	}
	value *= columnScale[iColumn];
	array[jColumn]=value;
      }
    }
    for (i=0;i<numberInRowArray;i++) {
      int iRow = whichRow[i];
      pi[iRow]=0.0;
    }
  } else {
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
	if (fabs(value)>zeroTolerance) {
	  index[numberNonZero++]=iColumn;
	  array[iColumn]=value;
	}
      }
    } else {
      // scaled
      for (jColumn=0;jColumn<numberToDo;jColumn++) {
	int iColumn = which[jColumn];
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
	  array[iColumn]=value;
	}
      }
    }
  }
}
/* Returns number of elements in basis
   column is basic if entry >=0 */
CoinBigIndex 
ClpPackedMatrix::numberInBasis(const int * columnIsBasic) const 
{
  int i;
  int numberColumns = getNumCols();
  const int * columnLength = matrix_->getVectorLengths(); 
  CoinBigIndex numberElements=0;
  if (!zeroElements_) {
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	numberElements += columnLength[i];
      }
    }
  } else {
    // there are zero elements so need to look more closely
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const double * elementByColumn = matrix_->getElements();
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	CoinBigIndex j;
	for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	  if (elementByColumn[j])
	    numberElements++;
	}
      }
    }
  }
  return numberElements;
}
// Fills in basis (Returns number of elements and updates numberBasic)
CoinBigIndex 
ClpPackedMatrix::fillBasis(const ClpSimplex * model,
				const int * columnIsBasic, int & numberBasic,
				int * indexRowU, int * indexColumnU,
				double * elementU) const 
{
  const double * rowScale = model->rowScale();
  const int * row = matrix_->getIndices();
  const CoinBigIndex * columnStart = matrix_->getVectorStarts();
  const int * columnLength = matrix_->getVectorLengths(); 
  const double * elementByColumn = matrix_->getElements();
  int i;
  int numberColumns = getNumCols();
  CoinBigIndex numberElements=0;
  if (!zeroElements_) {
    if (!rowScale) {
      // no scaling
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  CoinBigIndex j;
	  for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	    indexRowU[numberElements]=row[j];
	    indexColumnU[numberElements]=numberBasic;
	    elementU[numberElements++]=elementByColumn[j];
	  }
	  numberBasic++;
	}
      }
    } else {
      // scaling
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  CoinBigIndex j;
	  double scale = columnScale[i];
	  for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	    int iRow = row[j];
	    indexRowU[numberElements]=iRow;
	    indexColumnU[numberElements]=numberBasic;
	    elementU[numberElements++]=elementByColumn[j]*scale*rowScale[iRow];
	  }
	  numberBasic++;
	}
      }
    }
  } else {
    // there are zero elements so need to look more closely
    if (!rowScale) {
      // no scaling
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  CoinBigIndex j;
	  for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	    double value = elementByColumn[j];
	    if (value) {
	      indexRowU[numberElements]=row[j];
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value;
	    }
	  }
	  numberBasic++;
	}
      }
    } else {
      // scaling
      const double * columnScale = model->columnScale();
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  CoinBigIndex j;
	  double scale = columnScale[i];
	  for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	    double value = elementByColumn[j];
	    if (value) {
	      int iRow = row[j];
	      indexRowU[numberElements]=iRow;
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value*scale*rowScale[iRow];
	    }
	  }
	  numberBasic++;
	}
      }
    }
  }
  return numberElements;
}
/* If element NULL returns number of elements in column part of basis,
   If not NULL fills in as well */
CoinBigIndex 
ClpPackedMatrix::fillBasis(const ClpSimplex * model,
			   const int * whichColumn, 
			   int numberBasic,
			   int numberColumnBasic,
			   int * indexRowU, int * indexColumnU,
			   double * elementU) const 
{
  const int * columnLength = matrix_->getVectorLengths(); 
  int i;
  CoinBigIndex numberElements=0;
  if (elementU!=NULL) {
    // fill
    const CoinBigIndex * columnStart = matrix_->getVectorStarts();
    const double * rowScale = model->rowScale();
    const int * row = matrix_->getIndices();
    const double * elementByColumn = matrix_->getElements();
    if (!zeroElements_) {
      if (!rowScale) {
	// no scaling
	for (i=0;i<numberColumnBasic;i++) {
	  int iColumn = whichColumn[i];
	  CoinBigIndex j;
	  for (j=columnStart[iColumn];
	       j<columnStart[iColumn]+columnLength[iColumn];j++) {
	    indexRowU[numberElements]=row[j];
	    indexColumnU[numberElements]=numberBasic;
	    elementU[numberElements++]=elementByColumn[j];
	  }
	  numberBasic++;
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
	    indexColumnU[numberElements]=numberBasic;
	    elementU[numberElements++]=
	      elementByColumn[j]*scale*rowScale[iRow];
	  }
	  numberBasic++;
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
	      indexRowU[numberElements]=row[j];
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value;
	    }
	  }
	  numberBasic++;
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
	      indexColumnU[numberElements]=numberBasic;
	      elementU[numberElements++]=value*scale*rowScale[iRow];
	    }
	  }
	  numberBasic++;
	}
      }
    }
  } else {
    // just count - can be over so ignore zero problem
    for (i=0;i<numberColumnBasic;i++) {
      int iColumn = whichColumn[i];
      numberElements += columnLength[iColumn];
    }
  }
  return numberElements;
}
// Creates scales for column copy (rowCopy in model may be modified)
int 
ClpPackedMatrix::scale(ClpSimplex * model) const 
{
  ClpMatrixBase * rowCopyBase=model->rowCopy();
  if (!rowCopyBase) {
    // temporary copy
    rowCopyBase = reverseOrderedCopy();
  }
  ClpPackedMatrix* rowCopy =
    dynamic_cast< ClpPackedMatrix*>(rowCopyBase);

  // Make sure it is really a ClpPackedMatrix
  assert (rowCopy!=NULL);
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const double * element = rowCopy->getElements();
  double * rowScale = new double [numberRows];
  double * columnScale = new double [numberColumns];
  // we are going to mark bits we are interested in
  char * usefulRow = new char [numberRows];
  char * usefulColumn = new char [numberColumns];
  double * rowLower = model->rowLower();
  double * rowUpper = model->rowUpper();
  double * columnLower = model->columnLower();
  double * columnUpper = model->columnUpper();
  int iColumn, iRow;
  // mark free rows
  for (iRow=0;iRow<numberRows;iRow++) {
    usefulRow[iRow]=0;
    if (rowUpper[iRow]<1.0e20||
	rowLower[iRow]>-1.0e20)
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
    if (columnUpper[iColumn]>
	columnLower[iColumn]+1.0e-9) {
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	iRow=row[j];
	if(elementByColumn[j]&&usefulRow[iRow]) {
	  useful=1;
	  largest = max(largest,fabs(elementByColumn[j]));
	  smallest = min(smallest,fabs(elementByColumn[j]));
	}
      }
    }
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
    return 1;
  } else {
    int scalingMethod = model->scalingFlag();
    if (scalingMethod==3) {
      // Choose between 1 and 2
      if (smallest<1.0e-5||smallest*largest<1.0)
	scalingMethod=1;
      else
	scalingMethod=2;
    }
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
    ClpFillN ( rowScale, numberRows,1.0);
    ClpFillN ( columnScale, numberColumns,1.0);
    double overallLargest=-1.0e-30;
    double overallSmallest=1.0e30;
    if (scalingMethod==1) {
      // Maximum in each row
      for (iRow=0;iRow<numberRows;iRow++) {
	if (usefulRow[iRow]) {
	  CoinBigIndex j;
	  largest=1.0e-20;
	  for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	    int iColumn = column[j];
	    if (usefulColumn[iColumn]) {
	      double value = fabs(element[j]);
	      largest = max(largest,value);
	    }
	  }
	  rowScale[iRow]=1.0/largest;
	  overallLargest = max(overallLargest,largest);
	  overallSmallest = min(overallSmallest,largest);
	}
      }
    } else {
      assert(scalingMethod==2);
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
		if (value>1.0e-30) {
		  value *= columnScale[iColumn];
		  largest = max(largest,value);
		  smallest = min(smallest,value);
		}
	      }
	    }
	    rowScale[iRow]=1.0/sqrt(smallest*largest);
	    overallLargest = max(largest*rowScale[iRow],overallLargest);
	    overallSmallest = min(smallest*rowScale[iRow],overallSmallest);
	  }
	}
#ifdef USE_OBJECTIVE
	largest=1.0e-20;
	smallest=1.0e50;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (usefulColumn[iColumn]) {
	    double value = fabs(objective[iColumn]);
	    // Don't bother with tiny elements
	    if (value>1.0e-30) {
	      value *= columnScale[iColumn];
	      largest = max(largest,value);
	      smallest = min(smallest,value);
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
	      if (value>1.0e-30&&usefulRow[iRow]) {
		value *= rowScale[iRow];
		largest = max(largest,value);
		smallest = min(smallest,value);
	      }
	    }
#ifdef USE_OBJECTIVE
	    if (fabs(objective[iColumn])>1.0e-30) {
	      double value = fabs(objective[iColumn])*objScale;
	      largest = max(largest,value);
	      smallest = min(smallest,value);
	    }
#endif
	    columnScale[iColumn]=1.0/sqrt(smallest*largest);
	  }
	}
      }
#ifdef USE_OBJECTIVE
      delete [] objective;
      printf("obj scale %g - use it if you want\n",objScale);
#endif
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
	    largest = max(largest,value);
	    smallest = min(smallest,value);
	  }
	}
	if (overallSmallest*largest>smallest)
	  overallSmallest = smallest/largest;
      }
    }
#define RANDOMIZE
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
    overallLargest = min(1000.0,overallLargest);
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
	    largest = max(largest,value);
	    smallest = min(smallest,value);
	  }
	}
	columnScale[iColumn]=overallLargest/largest;
#ifdef RANDOMIZE
	double value = 0.5-CoinDrand48();//between -0.5 to + 0.5
	columnScale[iColumn] *= (1.0+0.1*value);
#endif
	overallSmallest = min(overallSmallest,smallest*columnScale[iColumn]);
      }
    }
    model->messageHandler()->message(CLP_PACKEDSCALE_FINAL,*model->messagesPointer())
      <<overallSmallest
      <<overallLargest
      <<CoinMessageEol;
    delete [] usefulRow;
    delete [] usefulColumn;
    // If quadratic then make symmetric
    ClpObjective * obj = model->objectiveAsObject();
    ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(obj));
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
	  int j;
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
    model->setRowScale(rowScale);
    model->setColumnScale(columnScale);
    if (model->rowCopy()) {
      // need to replace row by row
      double * newElement = new double[numberColumns];
      // scale row copy
      for (iRow=0;iRow<numberRows;iRow++) {
	int j;
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
  if (!rowScale) {
    int j=columnStart[iColumn];
    rowArray->createPacked(columnLength[iColumn],
			   row+j,elementByColumn+j);
  } else {
    // apply scaling
    double scale = model->columnScale()[iColumn];
    int * index = rowArray->getIndices();
    double * array = rowArray->denseVector();
    int number = 0;
    for (i=columnStart[iColumn];
	 i<columnStart[iColumn]+columnLength[iColumn];i++) {
      int iRow = row[i];
      array[number]=elementByColumn[i]*scale*rowScale[iRow];
      index[number++]=iRow;
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
/* Checks if all elements are in valid range.  Can just
   return true if you are not paranoid.  For Clp I will
   probably expect no zeros.  Code can modify matrix to get rid of
   small elements.
*/
bool 
ClpPackedMatrix::allElementsInRange(ClpModel * model,
				    double smallest, double largest)
{
  int iColumn;
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
  int numberColumns = matrix_->getNumCols();
  int numberRows = matrix_->getNumRows();
  int * mark = new int [numberRows];
  int i;
  for (i=0;i<numberRows;i++)
    mark[i]=-1;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    for (j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];j++) {
      double value = fabs(elementByColumn[j]);
      int iRow = row[j];
      if (iRow<0||iRow>=numberRows) {
	printf("Out of range %d %d %d %g\n",iColumn,j,row[j],elementByColumn[j]);
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
	zeroElements_ = true; // there are zero elements
      if (value<smallest) {
	numberSmall++;
      } else if (value>largest) {
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
    zeroElements_=false;
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
  int numberColumns =model->numberColumns();
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


