// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.


#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
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
    matrix_(NULL)
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
  
}

//-------------------------------------------------------------------
// assign matrix (for space reasons)
//-------------------------------------------------------------------
ClpPackedMatrix::ClpPackedMatrix (CoinPackedMatrix * rhs) 
: ClpMatrixBase()
{  
  matrix_ = rhs;
  setType(1);
  
}

ClpPackedMatrix::ClpPackedMatrix (const CoinPackedMatrix & rhs) 
: ClpMatrixBase()
{  
  matrix_ = new CoinPackedMatrix(rhs);
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
}
void 
ClpPackedMatrix::transposeTimes( double scalar,
				 const double * x, double * y,
				 const double * rowScale, 
				 const double * columnScale) const
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
    // scaled
    for (j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int jRow=row[j];
      value += x[jRow]*elementByColumn[j]*rowScale[jRow];
    }
    y[iColumn] += value*scalar*columnScale[iColumn];
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
  if (numberInRowArray>0.333*numberRows||!rowCopy) {
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
    } else {
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
  if (numberInRowArray>2||y->getNumElements()) {
    // do by rows
    // ** Row copy is already scaled
    int iRow;
    double * markVector = y->denseVector(); // probably empty .. but
    int * mark = y->getIndices();
    int numberOriginal=y->getNumElements();
    int i;
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
  } else if (numberInRowArray==2) {
    // do by rows when two rows
    int iRow;
    int numberOriginal;
    int i;
    numberNonZero=0;

    double value;
    iRow = whichRow[0]; 
    value = pi[iRow]*scalar;
    CoinBigIndex j;
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
  } else if (numberInRowArray==1) {
    // Just one row
    int iRow=rowArray->getIndices()[0];
    numberNonZero=0;
    double value = pi[iRow]*scalar;
    CoinBigIndex j;
    for (j=rowStart[iRow];j<rowStart[iRow+1];j++) {
      int iColumn = column[j];
      double value2 = value*element[j];
      if (fabs(value2)>zeroTolerance) {
	index[numberNonZero++]=iColumn;
	array[iColumn] = value2;
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
/* Returns number of elements in basis
   column is basic if entry >=0 */
CoinBigIndex 
ClpPackedMatrix::numberInBasis(const int * columnIsBasic) const 
{
  int i;
  int numberColumns = getNumCols();
  const int * columnLength = matrix_->getVectorLengths(); 
  CoinBigIndex numberElements=0;
  for (i=0;i<numberColumns;i++) {
    if (columnIsBasic[i]>=0) {
      numberElements += columnLength[i];
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
  assert (model->scalingFlag()==1); // dynamic not implemented
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
    int numberPass=3;
    double overallLargest;
    double overallSmallest;
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
	      double value = fabs(element[j]*columnScale[iColumn]);
	      largest = max(largest,value);
	      smallest = min(smallest,value);
	    }
	  }
	  rowScale[iRow]=1.0/sqrt(smallest*largest);
	  overallLargest = max(largest*rowScale[iRow],overallLargest);
	  overallSmallest = min(smallest*rowScale[iRow],overallSmallest);
	}
      }
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
	    if(elementByColumn[j]&&usefulRow[iRow]) {
	      double value = fabs(elementByColumn[j]*rowScale[iRow]);
	      largest = max(largest,value);
	      smallest = min(smallest,value);
	    }
	  }
	  columnScale[iColumn]=1.0/sqrt(smallest*largest);
	}
      }
    }
    // final pass to scale columns so largest is 1.0
    overallLargest=0.0;
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
	columnScale[iColumn]=1.0/largest;
	overallLargest = max(overallLargest,largest*columnScale[iColumn]);
	overallSmallest = min(overallSmallest,smallest*columnScale[iColumn]);
      }
    }
    model->messageHandler()->message(CLP_PACKEDSCALE_FINAL,*model->messagesPointer())
      <<overallSmallest
      <<overallLargest
      <<CoinMessageEol;
    
    delete [] usefulRow;
    delete [] usefulColumn;
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
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation */
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
ClpPackedMatrix::allElementsInRange(ClpSimplex * model,
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
  return true;
}


