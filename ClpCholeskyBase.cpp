// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <iostream>

#include "ClpCholeskyBase.hpp"
#include "ClpInterior.hpp"
#include "ClpHelperFunctions.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpMessage.hpp"
#include "ClpQuadraticObjective.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpCholeskyBase::ClpCholeskyBase (int denseThreshold) :
  type_(0),
  doKKT_(false),
  goDense_(0.7),
  choleskyCondition_(0.0),
  model_(NULL),
  numberTrials_(),
  numberRows_(0),
  status_(0),
  rowsDropped_(NULL),
  permuteInverse_(NULL),
  permute_(NULL),
  numberRowsDropped_(0),
  sparseFactor_(NULL),
  choleskyStart_(NULL),
  choleskyRow_(NULL),
  indexStart_(NULL),
  diagonal_(NULL),
  workDouble_(NULL),
  link_(NULL),
  workInteger_(NULL),
  clique_(NULL),
  sizeFactor_(0),
  sizeIndex_(0),
  firstDense_(0),
  rowCopy_(NULL),
  whichDense_(NULL),
  denseColumn_(NULL),
  dense_(NULL),
  denseThreshold_(denseThreshold)
{
  memset(integerParameters_,0,64*sizeof(int));
  memset(doubleParameters_,0,64*sizeof(double));
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpCholeskyBase::ClpCholeskyBase (const ClpCholeskyBase & rhs) :
  type_(rhs.type_),
  doKKT_(rhs.doKKT_),
  goDense_(rhs.goDense_),
  choleskyCondition_(rhs.choleskyCondition_),
  model_(rhs.model_),
  numberTrials_(rhs.numberTrials_),
  numberRows_(rhs.numberRows_),
  status_(rhs.status_),
  numberRowsDropped_(rhs.numberRowsDropped_)
{  
  rowsDropped_ = ClpCopyOfArray(rhs.rowsDropped_,numberRows_);
  permuteInverse_ = ClpCopyOfArray(rhs.permuteInverse_,numberRows_);
  permute_ = ClpCopyOfArray(rhs.permute_,numberRows_);
  sizeFactor_=rhs.sizeFactor_;
  sizeIndex_ = rhs.sizeIndex_;
  firstDense_ = rhs.firstDense_;
  sparseFactor_ = ClpCopyOfArray(rhs.sparseFactor_,rhs.sizeFactor_);
  choleskyStart_ = ClpCopyOfArray(rhs.choleskyStart_,numberRows_+1);
  indexStart_ = ClpCopyOfArray(rhs.indexStart_,numberRows_);
  choleskyRow_ = ClpCopyOfArray(rhs.choleskyRow_,sizeIndex_);
  diagonal_ = ClpCopyOfArray(rhs.diagonal_,numberRows_);
#if CLP_LONG_CHOLESKY!=1
  workDouble_ = ClpCopyOfArray(rhs.workDouble_,numberRows_);
#else
  // actually long double
  workDouble_ = (double *) ClpCopyOfArray((longWork *) rhs.workDouble_,numberRows_);
#endif
  link_ = ClpCopyOfArray(rhs.link_,numberRows_);
  workInteger_ = ClpCopyOfArray(rhs.workInteger_,numberRows_);
  clique_ = ClpCopyOfArray(rhs.clique_,numberRows_);
  memcpy(integerParameters_,rhs.integerParameters_,64*sizeof(int));
  memcpy(doubleParameters_,rhs.doubleParameters_,64*sizeof(double));
  rowCopy_ = rhs.rowCopy_->clone();
  whichDense_ = NULL;
  denseColumn_=NULL;
  dense_=NULL;
  denseThreshold_ = rhs.denseThreshold_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpCholeskyBase::~ClpCholeskyBase ()
{
  delete [] rowsDropped_;
  delete [] permuteInverse_;
  delete [] permute_;
  delete [] sparseFactor_;
  delete [] choleskyStart_;
  delete [] choleskyRow_;
  delete [] indexStart_;
  delete [] diagonal_;
  delete [] workDouble_;
  delete [] link_;
  delete [] workInteger_;
  delete [] clique_;
  delete rowCopy_;
  delete [] whichDense_;
  delete [] denseColumn_;
  delete dense_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpCholeskyBase &
ClpCholeskyBase::operator=(const ClpCholeskyBase& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    doKKT_ = rhs.doKKT_;
    goDense_ = rhs.goDense_;
    choleskyCondition_ = rhs.choleskyCondition_;
    model_ = rhs.model_;
    numberTrials_ = rhs.numberTrials_;
    numberRows_ = rhs.numberRows_;
    status_ = rhs.status_;
    numberRowsDropped_ = rhs.numberRowsDropped_;
    delete [] rowsDropped_;
    delete [] permuteInverse_;
    delete [] permute_;
    delete [] sparseFactor_;
    delete [] choleskyStart_;
    delete [] choleskyRow_;
    delete [] indexStart_;
    delete [] diagonal_;
    delete [] workDouble_;
    delete [] link_;
    delete [] workInteger_;
    delete [] clique_;
    delete rowCopy_;
    delete [] whichDense_;
    delete [] denseColumn_;
    delete dense_;
    rowsDropped_ = ClpCopyOfArray(rhs.rowsDropped_,numberRows_);
    permuteInverse_ = ClpCopyOfArray(rhs.permuteInverse_,numberRows_);
    permute_ = ClpCopyOfArray(rhs.permute_,numberRows_);
    sizeFactor_=rhs.sizeFactor_;
    sizeIndex_ = rhs.sizeIndex_;
    firstDense_ = rhs.firstDense_;
    sparseFactor_ = ClpCopyOfArray(rhs.sparseFactor_,rhs.sizeFactor_);
    choleskyStart_ = ClpCopyOfArray(rhs.choleskyStart_,numberRows_+1);
    choleskyRow_ = ClpCopyOfArray(rhs.choleskyRow_,rhs.sizeFactor_);
    indexStart_ = ClpCopyOfArray(rhs.indexStart_,numberRows_);
    choleskyRow_ = ClpCopyOfArray(rhs.choleskyRow_,sizeIndex_);
    diagonal_ = ClpCopyOfArray(rhs.diagonal_,numberRows_);
#if CLP_LONG_CHOLESKY!=1
    workDouble_ = ClpCopyOfArray(rhs.workDouble_,numberRows_);
#else
    // actually long double
    workDouble_ = (double *) ClpCopyOfArray((longWork *) rhs.workDouble_,numberRows_);
#endif
    link_ = ClpCopyOfArray(rhs.link_,numberRows_);
    workInteger_ = ClpCopyOfArray(rhs.workInteger_,numberRows_);
    clique_ = ClpCopyOfArray(rhs.clique_,numberRows_);
    delete rowCopy_;
    rowCopy_ = rhs.rowCopy_->clone();
    whichDense_ = NULL;
    denseColumn_=NULL;
    dense_=NULL;
    denseThreshold_ = rhs.denseThreshold_;
  }
  return *this;
}
// reset numberRowsDropped and rowsDropped.
void 
ClpCholeskyBase::resetRowsDropped()
{
  numberRowsDropped_=0;
  memset(rowsDropped_,0,numberRows_);
}
/* Uses factorization to solve. - given as if KKT.
   region1 is rows+columns, region2 is rows */
void 
ClpCholeskyBase::solveKKT (double * region1, double * region2, const double * diagonal,
			   double diagonalScaleFactor)
{
  if (!doKKT_) {
    int iColumn;
    int numberColumns = model_->numberColumns();
    int numberTotal = numberRows_+numberColumns;
    double * region1Save = new double[numberTotal];
    for (iColumn=0;iColumn<numberTotal;iColumn++) {
      region1[iColumn] *= diagonal[iColumn];
      region1Save[iColumn]=region1[iColumn];
    }
    multiplyAdd(region1+numberColumns,numberRows_,-1.0,region2,1.0);
    model_->clpMatrix()->times(1.0,region1,region2);
    double maximumRHS = maximumAbsElement(region2,numberRows_);
    double scale=1.0;
    double unscale=1.0;
    if (maximumRHS>1.0e-30) {
      if (maximumRHS<=0.5) {
	double factor=2.0;
	while (maximumRHS<=0.5) {
	  maximumRHS*=factor;
	  scale*=factor;
	} /* endwhile */
      } else if (maximumRHS>=2.0&&maximumRHS<=COIN_DBL_MAX) {
	double factor=0.5;
	while (maximumRHS>=2.0) {
	  maximumRHS*=factor;
	  scale*=factor;
	} /* endwhile */
      } 
      unscale=diagonalScaleFactor/scale;
    } else {
      //effectively zero
      scale=0.0;
      unscale=0.0;
    } 
    multiplyAdd(NULL,numberRows_,0.0,region2,scale);
    solve(region2);
    multiplyAdd(NULL,numberRows_,0.0,region2,unscale);
    multiplyAdd(region2,numberRows_,-1.0,region1+numberColumns,0.0);
    CoinZeroN(region1,numberColumns);
    model_->clpMatrix()->transposeTimes(1.0,region2,region1);
    for (iColumn=0;iColumn<numberTotal;iColumn++)
      region1[iColumn] = region1[iColumn]*diagonal[iColumn]-region1Save[iColumn];
    delete [] region1Save;
  } else {
    // KKT
    int numberRowsModel = model_->numberRows();
    int numberColumns = model_->numberColumns();
    int numberTotal = numberColumns + numberRowsModel;
    double * array = new double [numberRows_];
    CoinMemcpyN(region1,numberTotal,array);
    CoinMemcpyN(region2,numberRowsModel,array+numberTotal);
    solve(array);
    int iRow;
    for (iRow=0;iRow<numberTotal;iRow++) { 
      if (rowsDropped_[iRow]&&fabs(array[iRow])>1.0e-8) {
	printf("row region1 %d dropped %g\n",iRow,array[iRow]);
      }
    }
    for (;iRow<numberRows_;iRow++) {
      if (rowsDropped_[iRow]&&fabs(array[iRow])>1.0e-8) {
	printf("row region2 %d dropped %g\n",iRow,array[iRow]);
      }
    }
    CoinMemcpyN(array+numberTotal,numberRowsModel,region2);
    CoinMemcpyN(array,numberTotal,region1);
    delete [] array;
  }
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase * ClpCholeskyBase::clone() const
{
  return new ClpCholeskyBase(*this);
}
// Forms ADAT - returns nonzero if not enough memory
int 
ClpCholeskyBase::preOrder(bool lowerTriangular, bool includeDiagonal, bool doKKT)
{
  delete rowCopy_;
  rowCopy_ = model_->clpMatrix()->reverseOrderedCopy();
  if (!doKKT) {
    numberRows_ = model_->numberRows();
    rowsDropped_ = new char [numberRows_];
    memset(rowsDropped_,0,numberRows_);
    numberRowsDropped_=0;
    // Space for starts
    choleskyStart_ = new CoinBigIndex[numberRows_+1];
    const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
    const int * columnLength = model_->clpMatrix()->getVectorLengths();
    const int * row = model_->clpMatrix()->getIndices();
    const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
    const int * rowLength = rowCopy_->getVectorLengths();
    const int * column = rowCopy_->getIndices();
    // We need two arrays for counts
    int * which = new int [numberRows_];
    int * used = new int[numberRows_+1];
    CoinZeroN(used,numberRows_);
    int iRow;
    sizeFactor_=0;
    int numberColumns = model_->numberColumns();
    int numberDense=0;
    //denseThreshold_=3;
    if (denseThreshold_>0) {
      delete [] whichDense_;
      delete [] denseColumn_;
      delete dense_;
      whichDense_ = new char[numberColumns];
      int iColumn;
      used[numberRows_]=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	int length = columnLength[iColumn];
	used[length] += 1;
      }
      int nLong=0;
      int stop = CoinMax(denseThreshold_/2,100);
      for (iRow=numberRows_;iRow>=stop;iRow--) {
	if (used[iRow]) 
	  printf("%d columns are of length %d\n",used[iRow],iRow);
	nLong += used[iRow];
	if (nLong>50||nLong>(numberColumns>>2))
	  break;
      }
      CoinZeroN(used,numberRows_);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (columnLength[iColumn]<denseThreshold_) {
	  whichDense_[iColumn]=0;
	} else {
	  whichDense_[iColumn]=1;
	  numberDense++;
	}
      }
      if (!numberDense||numberDense>100) {
	// free
	delete [] whichDense_;
	whichDense_=NULL;
	denseColumn_=NULL;
	dense_=NULL;
      } else {
	// space for dense columns
	denseColumn_ = new longDouble [numberDense*numberRows_];
	// dense cholesky
	dense_ = new ClpCholeskyDense();
	dense_->reserveSpace(NULL,numberDense);
	printf("Taking %d columns as dense\n",numberDense);
      }
    }
    int offset = includeDiagonal ? 0 : 1;
    if (lowerTriangular)
      offset=-offset;
    for (iRow=0;iRow<numberRows_;iRow++) {
      int number=0;
      // make sure diagonal exists if includeDiagonal
      if (!offset) {
	which[0]=iRow;
	used[iRow]=1;
	number=1;
      }
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      if (lowerTriangular) {
	for (CoinBigIndex k=startRow;k<endRow;k++) {
	  int iColumn=column[k];
	  if (!whichDense_||!whichDense_[iColumn]) {
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      if (jRow<=iRow+offset) {
		if (!used[jRow]) {
		  used[jRow]=1;
		  which[number++]=jRow;
		}
	      }
	    }
	  }
	}
      } else {
	for (CoinBigIndex k=startRow;k<endRow;k++) {
	  int iColumn=column[k];
	  if (!whichDense_||!whichDense_[iColumn]) {
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      if (jRow>=iRow+offset) {
		if (!used[jRow]) {
		  used[jRow]=1;
		  which[number++]=jRow;
		}
	      }
	    }
	  }
	}
      }
      sizeFactor_ += number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
    }
    delete [] which;
    // Now we have size - create arrays and fill in
    try { 
      choleskyRow_ = new int [sizeFactor_];
    }
    catch (...) {
      // no memory
      delete [] choleskyStart_;
      choleskyStart_=NULL;
      return -1;
    }
    sizeFactor_=0;
    which = choleskyRow_;
    for (iRow=0;iRow<numberRows_;iRow++) {
      int number=0;
      // make sure diagonal exists if includeDiagonal
      if (!offset) {
	which[0]=iRow;
	used[iRow]=1;
	number=1;
      }
      choleskyStart_[iRow]=sizeFactor_;
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      if (lowerTriangular) {
	for (CoinBigIndex k=startRow;k<endRow;k++) {
	  int iColumn=column[k];
	  if (!whichDense_||!whichDense_[iColumn]) {
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      if (jRow<=iRow+offset) {
		if (!used[jRow]) {
		  used[jRow]=1;
		  which[number++]=jRow;
		}
	      }
	    }
	  }
	}
      } else {
	for (CoinBigIndex k=startRow;k<endRow;k++) {
	  int iColumn=column[k];
	  if (!whichDense_||!whichDense_[iColumn]) {
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      if (jRow>=iRow+offset) {
		if (!used[jRow]) {
		  used[jRow]=1;
		  which[number++]=jRow;
		}
	      }
	    }
	  }
	}
      }
      sizeFactor_ += number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
      // Sort
      std::sort(which,which+number);
      // move which on
      which += number;
    }
    choleskyStart_[numberRows_]=sizeFactor_;
    delete [] used;
    return 0;
  } else {
    int numberRowsModel = model_->numberRows();
    int numberColumns = model_->numberColumns();
    int numberTotal = numberColumns + numberRowsModel;
    numberRows_ = 2*numberRowsModel+numberColumns;
    rowsDropped_ = new char [numberRows_];
    memset(rowsDropped_,0,numberRows_);
    numberRowsDropped_=0;
    CoinPackedMatrix * quadratic = NULL;
    ClpQuadraticObjective * quadraticObj = 
      (dynamic_cast< ClpQuadraticObjective*>(model_->objectiveAsObject()));
    if (quadraticObj) 
      quadratic = quadraticObj->quadraticObjective();
    int numberElements = model_->clpMatrix()->getNumElements();
    numberElements = numberElements+2*numberRowsModel+numberTotal;
    if (quadratic)
      numberElements += quadratic->getNumElements(); 
    // Space for starts
    choleskyStart_ = new CoinBigIndex[numberRows_+1];
    const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
    const int * columnLength = model_->clpMatrix()->getVectorLengths();
    const int * row = model_->clpMatrix()->getIndices();
    //const double * element = model_->clpMatrix()->getElements();
    // Now we have size - create arrays and fill in
    try { 
      choleskyRow_ = new int [numberElements];
    }
    catch (...) {
      // no memory
      delete [] choleskyStart_;
      choleskyStart_=NULL;
      return -1;
    }
    int iRow,iColumn;
  
    sizeFactor_=0;
    // matrix
    if (lowerTriangular) {
      if (!quadratic) {
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  choleskyStart_[iColumn]=sizeFactor_;
	  choleskyRow_[sizeFactor_++]=iColumn;
	  CoinBigIndex start=columnStart[iColumn];
	  CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	  if (!includeDiagonal)
	    start++;
	  for (CoinBigIndex j=start;j<end;j++) {
	    choleskyRow_[sizeFactor_++]=row[j]+numberTotal;
	  }
	}
      } else {
	// Quadratic
	const int * columnQuadratic = quadratic->getIndices();
	const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	const int * columnQuadraticLength = quadratic->getVectorLengths();
	//const double * quadraticElement = quadratic->getElements();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  choleskyStart_[iColumn]=sizeFactor_;
	  if (includeDiagonal) 
	    choleskyRow_[sizeFactor_++]=iColumn;
	  for (CoinBigIndex j=columnQuadraticStart[iColumn];
	       j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	    int jColumn = columnQuadratic[j];
	    if (jColumn>iColumn)
	      choleskyRow_[sizeFactor_++]=jColumn;
	  }
	  CoinBigIndex start=columnStart[iColumn];
	  CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	  for (CoinBigIndex j=start;j<end;j++) {
	    choleskyRow_[sizeFactor_++]=row[j]+numberTotal;
	  }
	}
      }
      // slacks
      for (;iColumn<numberTotal;iColumn++) {
	choleskyStart_[iColumn]=sizeFactor_;
	if (includeDiagonal) 
	  choleskyRow_[sizeFactor_++]=iColumn;
	choleskyRow_[sizeFactor_++]=iColumn-numberColumns+numberTotal;
      }
      // Transpose - nonzero diagonal (may regularize)
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	choleskyStart_[iRow+numberTotal]=sizeFactor_;
	// diagonal
	if (includeDiagonal) 
	  choleskyRow_[sizeFactor_++]=iRow+numberTotal;
      }
      choleskyStart_[numberRows_]=sizeFactor_;
    } else {
      // transpose
      ClpMatrixBase * rowCopy = model_->clpMatrix()->reverseOrderedCopy();
      const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
      const int * rowLength = rowCopy->getVectorLengths();
      const int * column = rowCopy->getIndices();
      if (!quadratic) {
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  choleskyStart_[iColumn]=sizeFactor_;
	  if (includeDiagonal) 
	    choleskyRow_[sizeFactor_++]=iColumn;
	}
      } else {
	// Quadratic
	// transpose
	CoinPackedMatrix quadraticT;
	quadraticT.reverseOrderedCopyOf(*quadratic);
	const int * columnQuadratic = quadraticT.getIndices();
	const CoinBigIndex * columnQuadraticStart = quadraticT.getVectorStarts();
	const int * columnQuadraticLength = quadraticT.getVectorLengths();
	//const double * quadraticElement = quadraticT.getElements();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  choleskyStart_[iColumn]=sizeFactor_;
	  for (CoinBigIndex j=columnQuadraticStart[iColumn];
	       j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	    int jColumn = columnQuadratic[j];
	    if (jColumn<iColumn)
	      choleskyRow_[sizeFactor_++]=jColumn;
	  }
	  if (includeDiagonal) 
	    choleskyRow_[sizeFactor_++]=iColumn;
	}
      }
      int iRow;
      // slacks
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	choleskyStart_[iRow+numberColumns]=sizeFactor_;
	if (includeDiagonal) 
	  choleskyRow_[sizeFactor_++]=iRow+numberColumns;
      }
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	choleskyStart_[iRow+numberTotal]=sizeFactor_;
	CoinBigIndex start=rowStart[iRow];
	CoinBigIndex end=rowStart[iRow]+rowLength[iRow];
	for (CoinBigIndex j=start;j<end;j++) {
	  choleskyRow_[sizeFactor_++]=column[j];
	}
	// slack
	choleskyRow_[sizeFactor_++]=numberColumns+iRow;
	if (includeDiagonal)
	  choleskyRow_[sizeFactor_++]=iRow+numberTotal;
      }
      choleskyStart_[numberRows_]=sizeFactor_;
    }
  }
  return 0;
}
/* Orders rows and saves pointer to matrix.and model */
int 
ClpCholeskyBase::order(ClpInterior * model) 
{
  model_=model;
  int numberRowsModel = model_->numberRows();
  int numberColumns = model_->numberColumns();
  int numberTotal = numberColumns + numberRowsModel;
  CoinPackedMatrix * quadratic = NULL;
  ClpQuadraticObjective * quadraticObj = 
    (dynamic_cast< ClpQuadraticObjective*>(model_->objectiveAsObject()));
  if (quadraticObj) 
    quadratic = quadraticObj->quadraticObjective();
  if (!doKKT_) {
    numberRows_ = model->numberRows();
  } else {
    numberRows_ = 2*numberRowsModel+numberColumns;
  }
  rowsDropped_ = new char [numberRows_];
  numberRowsDropped_=0;
  rowCopy_ = model->clpMatrix()->reverseOrderedCopy();
  const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
  const int * columnLength = model_->clpMatrix()->getVectorLengths();
  const int * row = model_->clpMatrix()->getIndices();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  // We need arrays for counts
  int * which = new int [numberRows_];
  int * used = new int[numberRows_+1];
  int * count = new int[numberRows_];
  CoinZeroN(count,numberRows_);
  CoinZeroN(used,numberRows_);
  int iRow;
  sizeFactor_=0;
  permute_ = new int[numberRows_];
  for (iRow=0;iRow<numberRows_;iRow++) 
    permute_[iRow]=iRow; 
  if (!doKKT_) {
    int numberDense=0;
    if (denseThreshold_>0) {
      delete [] whichDense_;
      delete [] denseColumn_;
      delete dense_;
      whichDense_ = new char[numberColumns];
      int iColumn;
      used[numberRows_]=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	int length = columnLength[iColumn];
	used[length] += 1;
      }
      int nLong=0;
      int stop = CoinMax(denseThreshold_/2,100);
      for (iRow=numberRows_;iRow>=stop;iRow--) {
	if (used[iRow]) 
	  printf("%d columns are of length %d\n",used[iRow],iRow);
	nLong += used[iRow];
	if (nLong>50||nLong>(numberColumns>>2))
	  break;
      }
      CoinZeroN(used,numberRows_);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (columnLength[iColumn]<denseThreshold_) {
	  whichDense_[iColumn]=0;
	} else {
	  whichDense_[iColumn]=1;
	  numberDense++;
	}
      }
      if (!numberDense||numberDense>100) {
	// free
	delete [] whichDense_;
	whichDense_=NULL;
	denseColumn_=NULL;
	dense_=NULL;
      } else {
	// space for dense columns
	denseColumn_ = new longDouble [numberDense*numberRows_];
	// dense cholesky
	dense_ = new ClpCholeskyDense();
	dense_->reserveSpace(NULL,numberDense);
	printf("Taking %d columns as dense\n",numberDense);
      }
    }
    /* 
       Get row counts and size
    */
    for (iRow=0;iRow<numberRows_;iRow++) {
      int number=1;
      // make sure diagonal exists
      which[0]=iRow;
      used[iRow]=1;
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	if (!whichDense_||!whichDense_[iColumn]) {
	  CoinBigIndex start=columnStart[iColumn];
	  CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	  for (CoinBigIndex j=start;j<end;j++) {
	    int jRow=row[j];
	    if (jRow<iRow) {
	      if (!used[jRow]) {
		used[jRow]=1;
		which[number++]=jRow;
		count[jRow]++;
	      }
	    }
	  }
	}
      }
      sizeFactor_ += number;
      count[iRow]+=number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
    }
    CoinSort_2(count,count+numberRows_,permute_);
  } else {
    // KKT
    int numberElements = model_->clpMatrix()->getNumElements();
    numberElements = numberElements+2*numberRowsModel+numberTotal;
    if (quadratic)
      numberElements += quadratic->getNumElements(); 
    // off diagonal
    numberElements -= numberRows_;
    sizeFactor_=numberElements;
    // If we sort we need to redo symbolic etc
  }
  delete [] which;
  delete [] used;
  delete [] count;
  permuteInverse_ = new int [numberRows_];
  memset(rowsDropped_,0,numberRows_);
  for (iRow=0;iRow<numberRows_;iRow++) {
    //permute_[iRow]=iRow; // force no permute
    //permute_[iRow]=numberRows_-1-iRow; // force odd permute
    //permute_[iRow]=(iRow+1)%numberRows_; // force odd permute
    permuteInverse_[permute_[iRow]]=iRow;
  }
  return 0;
}
/* Does Symbolic factorization given permutation.
   This is called immediately after order.  If user provides this then
   user must provide factorize and solve.  Otherwise the default factorization is used
   returns non-zero if not enough memory */
int 
ClpCholeskyBase::symbolic()
{
  const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
  const int * columnLength = model_->clpMatrix()->getVectorLengths();
  const int * row = model_->clpMatrix()->getIndices();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  int numberRowsModel = model_->numberRows();
  int numberColumns = model_->numberColumns();
  int numberTotal = numberColumns + numberRowsModel;
  CoinPackedMatrix * quadratic = NULL;
  ClpQuadraticObjective * quadraticObj = 
    (dynamic_cast< ClpQuadraticObjective*>(model_->objectiveAsObject()));
  if (quadraticObj) 
    quadratic = quadraticObj->quadraticObjective();
  // We need an array for counts
  int * used = new int[numberRows_+1];
  // If KKT then re-order so negative first
  if (doKKT_) {
    int nn=0;
    int np=0;
    int iRow;
    for (iRow=0;iRow<numberRows_;iRow++) {
      int originalRow = permute_[iRow];
      if (originalRow<numberTotal)
	permute_[nn++]=originalRow;
      else
	used[np++]=originalRow;
    }
    memcpy(permute_+nn,used,np*sizeof(int));
    for (iRow=0;iRow<numberRows_;iRow++) 
      permuteInverse_[permute_[iRow]]=iRow;
  }
  CoinZeroN(used,numberRows_);
  int iRow;
  int iColumn;
  bool noMemory=false;
  CoinBigIndex * Astart = new CoinBigIndex[numberRows_+1];
  int * Arow=NULL;
  try { 
    Arow = new int [sizeFactor_];
  }
  catch (...) {
    // no memory
    delete [] Astart;
    return -1;
  }
  choleskyStart_ = new int[numberRows_+1];
  link_ = new int[numberRows_];
  workInteger_ = new CoinBigIndex[numberRows_];
  indexStart_ = new CoinBigIndex[numberRows_];
  clique_ = new int[numberRows_];
  // Redo so permuted upper triangular
  sizeFactor_=0;
  int * which = Arow;
  if (!doKKT_) {
    for (iRow=0;iRow<numberRows_;iRow++) {
      int number=0;
      int iOriginalRow = permute_[iRow];
      Astart[iRow]=sizeFactor_;
      CoinBigIndex startRow=rowStart[iOriginalRow];
      CoinBigIndex endRow=rowStart[iOriginalRow]+rowLength[iOriginalRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	if (!whichDense_||!whichDense_[iColumn]) {
	  CoinBigIndex start=columnStart[iColumn];
	  CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	  for (CoinBigIndex j=start;j<end;j++) {
	    int jRow=row[j];
	    int jNewRow = permuteInverse_[jRow];
	    if (jNewRow<iRow) {
	      if (!used[jNewRow]) {
		used[jNewRow]=1;
		which[number++]=jNewRow;
	      }
	    }
	  }
	}
      }
      sizeFactor_ += number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
      // Sort
      std::sort(which,which+number);
      // move which on
      which += number;
    }
  } else {
    // KKT
    // transpose
    ClpMatrixBase * rowCopy = model_->clpMatrix()->reverseOrderedCopy();
    const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
    const int * rowLength = rowCopy->getVectorLengths();
    const int * column = rowCopy->getIndices();
    // temp
    bool permuted=false;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (permute_[iRow]!=iRow) {
	permuted=true;
	break;
      }
    }
    if (permuted) {
      // Need to permute - ugly
      if (!quadratic) {
	for (iRow=0;iRow<numberRows_;iRow++) {
	  Astart[iRow]=sizeFactor_;
	  int iOriginalRow = permute_[iRow];
	  if (iOriginalRow<numberColumns) {
	    // A may be upper triangular by mistake
	    iColumn=iOriginalRow;
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int kRow = row[j]+numberTotal;
	      kRow=permuteInverse_[kRow];
	      if (kRow<iRow)
		Arow[sizeFactor_++]=kRow;
	    }
	  } else if (iOriginalRow<numberTotal) {
	    int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	    if (kRow<iRow)
	      Arow[sizeFactor_++]=kRow;
	  } else {
	    int kRow = iOriginalRow-numberTotal;
	    CoinBigIndex start=rowStart[kRow];
	    CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow = column[j];
	      int jNewRow = permuteInverse_[jRow];
	      if (jNewRow<iRow)
		Arow[sizeFactor_++]=jNewRow;
	    }
	    // slack - should it be permute
	    kRow = permuteInverse_[kRow+numberColumns];
	    if (kRow<iRow)
	      Arow[sizeFactor_++]=kRow;
	  }
	  // Sort
	  std::sort(Arow+Astart[iRow],Arow+sizeFactor_);
	}
      } else {
	// quadratic
	// transpose
	CoinPackedMatrix quadraticT;
	quadraticT.reverseOrderedCopyOf(*quadratic);
	const int * columnQuadratic = quadraticT.getIndices();
	const CoinBigIndex * columnQuadraticStart = quadraticT.getVectorStarts();
	const int * columnQuadraticLength = quadraticT.getVectorLengths();
	for (iRow=0;iRow<numberRows_;iRow++) {
	  Astart[iRow]=sizeFactor_;
	  int iOriginalRow = permute_[iRow];
	  if (iOriginalRow<numberColumns) {
	    // Quadratic bit
	    CoinBigIndex j;
	    for ( j=columnQuadraticStart[iOriginalRow];
		  j<columnQuadraticStart[iOriginalRow]+columnQuadraticLength[iOriginalRow];j++) {
	      int jColumn = columnQuadratic[j];
	      int jNewColumn = permuteInverse_[jColumn];
	      if (jNewColumn<iRow)
		Arow[sizeFactor_++]=jNewColumn;
	    }
	    // A may be upper triangular by mistake
	    iColumn=iOriginalRow;
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (j=start;j<end;j++) {
	      int kRow = row[j]+numberTotal;
	      kRow=permuteInverse_[kRow];
	      if (kRow<iRow)
		Arow[sizeFactor_++]=kRow;
	    }
	  } else if (iOriginalRow<numberTotal) {
	    int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	    if (kRow<iRow)
	      Arow[sizeFactor_++]=kRow;
	  } else {
	    int kRow = iOriginalRow-numberTotal;
	    CoinBigIndex start=rowStart[kRow];
	    CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow = column[j];
	      int jNewRow = permuteInverse_[jRow];
	      if (jNewRow<iRow)
		Arow[sizeFactor_++]=jNewRow;
	    }
	    // slack - should it be permute
	    kRow = permuteInverse_[kRow+numberColumns];
	    if (kRow<iRow)
	      Arow[sizeFactor_++]=kRow;
	  }
	  // Sort
	  std::sort(Arow+Astart[iRow],Arow+sizeFactor_);
	}
      }
    } else {
      if (!quadratic) {
	for (iRow=0;iRow<numberRows_;iRow++) {
	  Astart[iRow]=sizeFactor_;
	}
      } else {
	// Quadratic
	// transpose
	CoinPackedMatrix quadraticT;
	quadraticT.reverseOrderedCopyOf(*quadratic);
	const int * columnQuadratic = quadraticT.getIndices();
	const CoinBigIndex * columnQuadraticStart = quadraticT.getVectorStarts();
	const int * columnQuadraticLength = quadraticT.getVectorLengths();
	//const double * quadraticElement = quadraticT.getElements();
	for (iRow=0;iRow<numberColumns;iRow++) {
	  int iOriginalRow = permute_[iRow];
	  Astart[iRow]=sizeFactor_;
	  for (CoinBigIndex j=columnQuadraticStart[iOriginalRow];
	       j<columnQuadraticStart[iOriginalRow]+columnQuadraticLength[iOriginalRow];j++) {
	    int jColumn = columnQuadratic[j];
	    int jNewColumn = permuteInverse_[jColumn];
	    if (jNewColumn<iRow)
	      Arow[sizeFactor_++]=jNewColumn;
	  }
	}
      }
      int iRow;
      // slacks
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	Astart[iRow+numberColumns]=sizeFactor_;
      }
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	Astart[iRow+numberTotal]=sizeFactor_;
	CoinBigIndex start=rowStart[iRow];
	CoinBigIndex end=rowStart[iRow]+rowLength[iRow];
	for (CoinBigIndex j=start;j<end;j++) {
	  Arow[sizeFactor_++]=column[j];
	}
	// slack
	Arow[sizeFactor_++]=numberColumns+iRow;
      }
    }
  }
  Astart[numberRows_]=sizeFactor_;
  firstDense_=numberRows_;
  symbolic1(Astart,Arow);
  // Now fill in indices
  try { 
    // too big
    choleskyRow_ = new int[sizeFactor_];
  }
  catch (...) {
    // no memory
    noMemory=true;
  } 
  double sizeFactor=sizeFactor_;
  if (!noMemory) {
    // Do lower triangular
    sizeFactor_=0;
    int * which = Arow;
    if (!doKKT_) {
      for (iRow=0;iRow<numberRows_;iRow++) {
	int number=0;
	int iOriginalRow = permute_[iRow];
	Astart[iRow]=sizeFactor_;
	if (!rowsDropped_[iOriginalRow]) {
	  CoinBigIndex startRow=rowStart[iOriginalRow];
	  CoinBigIndex endRow=rowStart[iOriginalRow]+rowLength[iOriginalRow];
	  for (CoinBigIndex k=startRow;k<endRow;k++) {
	    int iColumn=column[k];
	    if (!whichDense_||!whichDense_[iColumn]) {
	      CoinBigIndex start=columnStart[iColumn];
	      CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	      for (CoinBigIndex j=start;j<end;j++) {
		int jRow=row[j];
		int jNewRow = permuteInverse_[jRow];
		if (jNewRow>iRow&&!rowsDropped_[jRow]) {
		  if (!used[jNewRow]) {
		    used[jNewRow]=1;
		    which[number++]=jNewRow;
		  }
		}
	      }
	    }
	  }
	  sizeFactor_ += number;
	  int j;
	  for (j=0;j<number;j++)
	    used[which[j]]=0;
	  // Sort
	  std::sort(which,which+number);
	  // move which on
	  which += number;
	}
      }
    } else {
      // KKT
      // temp
      bool permuted=false;
      for (iRow=0;iRow<numberRows_;iRow++) {
	if (permute_[iRow]!=iRow) {
	  permuted=true;
	  break;
	}
      }
      // but fake it
      for (iRow=0;iRow<numberRows_;iRow++) {
	//permute_[iRow]=iRow; // force no permute
	//permuteInverse_[permute_[iRow]]=iRow;
      }
      if (permuted) {
	// Need to permute - ugly
	if (!quadratic) {
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    Astart[iRow]=sizeFactor_;
	    int iOriginalRow = permute_[iRow];
	    if (iOriginalRow<numberColumns) {
	      iColumn=iOriginalRow;
	      CoinBigIndex start=columnStart[iColumn];
	      CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	      for (CoinBigIndex j=start;j<end;j++) {
		int kRow = row[j]+numberTotal;
		kRow=permuteInverse_[kRow];
		if (kRow>iRow)
		  Arow[sizeFactor_++]=kRow;
	      }
	    } else if (iOriginalRow<numberTotal) {
	      int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	      if (kRow>iRow)
		Arow[sizeFactor_++]=kRow;
	    } else {
	      int kRow = iOriginalRow-numberTotal;
	      CoinBigIndex start=rowStart[kRow];
	      CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	      for (CoinBigIndex j=start;j<end;j++) {
		int jRow = column[j];
		int jNewRow = permuteInverse_[jRow];
		if (jNewRow>iRow)
		  Arow[sizeFactor_++]=jNewRow;
	      }
	      // slack - should it be permute
	      kRow = permuteInverse_[kRow+numberColumns];
	      if (kRow>iRow)
		Arow[sizeFactor_++]=kRow;
	    }
	    // Sort
	    std::sort(Arow+Astart[iRow],Arow+sizeFactor_);
	  }
	} else {
	  // quadratic
	  const int * columnQuadratic = quadratic->getIndices();
	  const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	  const int * columnQuadraticLength = quadratic->getVectorLengths();
	  for (iRow=0;iRow<numberRows_;iRow++) {
	    Astart[iRow]=sizeFactor_;
	    int iOriginalRow = permute_[iRow];
	    if (iOriginalRow<numberColumns) {
	      // Quadratic bit
	      CoinBigIndex j;
	      for ( j=columnQuadraticStart[iOriginalRow];
		    j<columnQuadraticStart[iOriginalRow]+columnQuadraticLength[iOriginalRow];j++) {
		int jColumn = columnQuadratic[j];
		int jNewColumn = permuteInverse_[jColumn];
		if (jNewColumn>iRow)
		  Arow[sizeFactor_++]=jNewColumn;
	      }
	      iColumn=iOriginalRow;
	      CoinBigIndex start=columnStart[iColumn];
	      CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	      for (j=start;j<end;j++) {
		int kRow = row[j]+numberTotal;
		kRow=permuteInverse_[kRow];
		if (kRow>iRow)
		  Arow[sizeFactor_++]=kRow;
	      }
	    } else if (iOriginalRow<numberTotal) {
	      int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	      if (kRow>iRow)
		Arow[sizeFactor_++]=kRow;
	    } else {
	      int kRow = iOriginalRow-numberTotal;
	      CoinBigIndex start=rowStart[kRow];
	      CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	      for (CoinBigIndex j=start;j<end;j++) {
		int jRow = column[j];
		int jNewRow = permuteInverse_[jRow];
		if (jNewRow>iRow)
		  Arow[sizeFactor_++]=jNewRow;
	      }
	      // slack - should it be permute
	      kRow = permuteInverse_[kRow+numberColumns];
	      if (kRow>iRow)
		Arow[sizeFactor_++]=kRow;
	    }
	    // Sort
	    std::sort(Arow+Astart[iRow],Arow+sizeFactor_);
	  }
	}
      } else {
	if (!quadratic) {
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
	    Astart[iColumn]=sizeFactor_;
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      Arow[sizeFactor_++]=row[j]+numberTotal;
	    }
	  }
	} else {
	  // Quadratic
	  const int * columnQuadratic = quadratic->getIndices();
	  const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	  const int * columnQuadraticLength = quadratic->getVectorLengths();
	  //const double * quadraticElement = quadratic->getElements();
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
	    Astart[iColumn]=sizeFactor_;
	    for (CoinBigIndex j=columnQuadraticStart[iColumn];
		 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	      int jColumn = columnQuadratic[j];
	      if (jColumn>iColumn)
		Arow[sizeFactor_++]=jColumn;
	    }
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      Arow[sizeFactor_++]=row[j]+numberTotal;
	    }
	  }
	}
	// slacks
	for (iRow=0;iRow<numberRowsModel;iRow++) {
	  Astart[iRow+numberColumns]=sizeFactor_;
	  Arow[sizeFactor_++]=iRow+numberTotal;
	}
	// Transpose - nonzero diagonal (may regularize)
	for (iRow=0;iRow<numberRowsModel;iRow++) {
	  Astart[iRow+numberTotal]=sizeFactor_;
	}
      }
      Astart[numberRows_]=sizeFactor_;
    }
    symbolic2(Astart,Arow);
    if (sizeIndex_<sizeFactor_) {
      int * indices=NULL;
      try { 
	indices = new int[sizeIndex_];
      }
      catch (...) {
	// no memory
	noMemory=true;
      } 
      if (!noMemory)  {
	memcpy(indices,choleskyRow_,sizeIndex_*sizeof(int));
	delete [] choleskyRow_;
	choleskyRow_=indices;
      }
    }
  }
  delete [] used;
  // Use cholesky regions
  delete [] Astart;
  delete [] Arow;
  double flops=0.0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int length = choleskyStart_[iRow+1]-choleskyStart_[iRow];
    flops += ((double) length) * (length + 2.0);
  }
  if (model_->messageHandler()->logLevel()>0) 
    std::cout<<sizeFactor<<" elements in sparse Cholesky, flop count "<<flops<<std::endl;
  try { 
    sparseFactor_ = new longDouble [sizeFactor_];
#if CLP_LONG_CHOLESKY!=1
    workDouble_ = new longDouble[numberRows_];
#else
    // actually long double
    workDouble_ = (double *) (new longWork[numberRows_]);
#endif
    diagonal_ = new longDouble[numberRows_];
  }
  catch (...) {
    // no memory
    noMemory=true;
  }
  if (noMemory) {
    delete [] choleskyRow_;
    choleskyRow_=NULL;
    delete [] choleskyStart_;
    choleskyStart_=NULL;
    delete [] permuteInverse_;
    permuteInverse_ = NULL;
    delete [] permute_;
    permute_ = NULL;
    delete [] choleskyStart_;
    choleskyStart_ = NULL;
    delete [] indexStart_;
    indexStart_ = NULL;
    delete [] link_;
    link_ = NULL;
    delete [] workInteger_;
    workInteger_ = NULL;
    delete [] sparseFactor_;
    sparseFactor_ = NULL;
    delete [] workDouble_;
    workDouble_ = NULL;
    delete [] diagonal_;
    diagonal_ = NULL;
    delete [] clique_;
    clique_ = NULL;
    return -1;
  }
  return 0;
}
int
ClpCholeskyBase::symbolic1(const CoinBigIndex * Astart, const int * Arow)
{
  int * marked = (int *) workInteger_;
  int iRow;
  // may not need to do this here but makes debugging easier
  for (iRow=0;iRow<numberRows_;iRow++) {
    marked[iRow]=-1;
    link_[iRow]=-1;
    choleskyStart_[iRow]=0; // counts
  }
  for (iRow=0;iRow<numberRows_;iRow++) {
    marked[iRow]=iRow;
    for (CoinBigIndex j=Astart[iRow];j<Astart[iRow+1];j++) {
      int kRow=Arow[j];
      while (marked[kRow] != iRow ) {
	if (link_[kRow] <0 )
	  link_[kRow]=iRow;
	choleskyStart_[kRow]++;
	marked[kRow]=iRow;
	kRow = link_[kRow];
      }
    }
  }
  sizeFactor_=0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int number = choleskyStart_[iRow];
    choleskyStart_[iRow]=sizeFactor_;
    sizeFactor_ += number;
  }
  choleskyStart_[numberRows_]=sizeFactor_;
  return sizeFactor_;;
}
void
ClpCholeskyBase::symbolic2(const CoinBigIndex * Astart, const int * Arow)
{
  int * mergeLink = clique_;
  int * marker = (int *) workInteger_;
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    marker[iRow]=-1;
    mergeLink[iRow]=-1;
    link_[iRow]=-1; // not needed but makes debugging easier
  }
  int start=0;
  int end=0;
  choleskyStart_[0]=0;
    
  for (iRow=0;iRow<numberRows_;iRow++) {
    int nz=0;
    int merge = mergeLink[iRow];
    bool marked=false;
    if (merge<0)
      marker[iRow]=iRow;
    else
      marker[iRow]=merge;
    start = end;
    int startSub=start;
    link_[iRow]=numberRows_;
    for (CoinBigIndex j=Astart[iRow];j<Astart[iRow+1];j++) {
      int kRow=Arow[j];
      int k=iRow;
      int linked = link_[iRow];
      while (linked<=kRow) {
	k=linked;
	linked = link_[k];
      }
      nz++;
      link_[k]=kRow;
      link_[kRow]=linked;
      if (marker[kRow] != marker[iRow])
	marked=true;
    }
    bool reuse=false;
    // Check if we can re-use indices
    if (!marked && merge>=0&&mergeLink[merge]<0) {
      // can re-use all
      startSub = indexStart_[merge]+1;
      nz=choleskyStart_[merge+1]-(choleskyStart_[merge]+1);
      reuse=true;
    } else {
      // See if we can re-use any
      int k=mergeLink[iRow];
      int maxLength=0;
      while (k>=0) {
	int length = choleskyStart_[k+1]-(choleskyStart_[k]+1);
	int start = indexStart_[k]+1;
	int stop = start+length;
	if (length>maxLength) {
	  maxLength = length;
	  startSub = start;
	}
	int linked = iRow;
	for (CoinBigIndex j=start;j<stop;j++) {
	  int kRow=choleskyRow_[j];
	  int kk=linked;
	  linked = link_[kk];
	  while (linked<kRow) {
	    kk=linked;
	    linked = link_[kk];
	  }
	  if (linked!=kRow) {
	    nz++;
	    link_[kk]=kRow;
	    link_[kRow]=linked;
	    linked=kRow;
	  }
	}
	k=mergeLink[k];
      }
      if (nz== maxLength) 
	reuse=true; // can re-use
    }
    //reuse=false; //temp
    if (!reuse) {
      end += nz;
      startSub=start;
      int kRow = iRow;
      for (int j=start;j<end;j++) {
	kRow=link_[kRow];
	choleskyRow_[j]=kRow;
	assert (kRow<numberRows_);
	marker[kRow]=iRow;
      }
      marker[iRow]=iRow;
    }
    indexStart_[iRow]=startSub;
    choleskyStart_[iRow+1]=choleskyStart_[iRow]+nz;
    if (nz>1) {
      int kRow = choleskyRow_[startSub];
      mergeLink[iRow]=mergeLink[kRow];
      mergeLink[kRow]=iRow;
    }
    // should not be needed
    //std::sort(choleskyRow_+indexStart_[iRow]
    //      ,choleskyRow_+indexStart_[iRow]+nz);
#define CLP_DEBUG
#ifdef CLP_DEBUG
    int last=-1;
    for (CoinBigIndex j=indexStart_[iRow];j<indexStart_[iRow]+nz;j++) {
      int kRow=choleskyRow_[j];
      assert (kRow>last);
      last=kRow;
    }
#endif    
  }
  sizeFactor_ = choleskyStart_[numberRows_];
  sizeIndex_ = start;
  // find dense segment here
  int numberleft=numberRows_;
  for (iRow=0;iRow<numberRows_;iRow++) {
    CoinBigIndex left=sizeFactor_-choleskyStart_[iRow];
    double n=numberleft;
    double threshold = n*(n-1.0)*0.5*goDense_;
    if ((double) left >= threshold) 
      break;
    numberleft--;
  }
  //iRow=numberRows_;
  int nDense = numberRows_-iRow;
#define DENSE_THRESHOLD 8
  // don't do if dense columns
  if (nDense>=DENSE_THRESHOLD&&!dense_) {
    printf("Going dense for last %d rows\n",nDense);
    // make sure we don't disturb any indices
    CoinBigIndex k=0;
    for (int jRow=0;jRow<iRow;jRow++) {
      int nz=choleskyStart_[jRow+1]-choleskyStart_[jRow];
      k=CoinMax(k,indexStart_[jRow]+nz);
    }
    indexStart_[iRow]=k;
    int j;
    for (j=iRow+1;j<numberRows_;j++) {
      choleskyRow_[k++]=j;
      indexStart_[j]=k;
    }
    sizeIndex_=k;
    assert (k<=sizeFactor_); // can't happen with any reasonable defaults
    k=choleskyStart_[iRow];
    for (j=iRow+1;j<=numberRows_;j++) {
      k += numberRows_-j;
      choleskyStart_[j]=k;
    }
    // allow for blocked dense
    ClpCholeskyDense dense;
    sizeFactor_=choleskyStart_[iRow]+dense.space(nDense);
    firstDense_=iRow;
    if (doKKT_) {
      // redo permute so negative ones first
      int putN=firstDense_;
      int putP=0;
      int numberRowsModel = model_->numberRows();
      int numberColumns = model_->numberColumns();
      int numberTotal = numberColumns + numberRowsModel;
      for (iRow=firstDense_;iRow<numberRows_;iRow++) {
	int originalRow=permute_[iRow];
	if (originalRow<numberTotal)
	  permute_[putN++]=originalRow;
	else
	  permuteInverse_[putP++]=originalRow;
      }
      for (iRow=putN;iRow<numberRows_;iRow++) {
	permute_[iRow]=permuteInverse_[iRow-putN];
      }
      for (iRow=0;iRow<numberRows_;iRow++) {
	permuteInverse_[permute_[iRow]]=iRow;
      }
    }
  }
  // Clean up clique info
  for (iRow=0;iRow<numberRows_;iRow++)
    clique_[iRow]=0;
  int lastClique=-1;
  bool inClique=false;
  for (iRow=1;iRow<firstDense_;iRow++) {
    int sizeLast = choleskyStart_[iRow]-choleskyStart_[iRow-1];
    int sizeThis = choleskyStart_[iRow+1]-choleskyStart_[iRow];
    if (indexStart_[iRow]==indexStart_[iRow-1]+1&&
	sizeThis==sizeLast-1&&
	sizeThis) {
      // in clique
      if (!inClique) {
	inClique=true;
	lastClique=iRow-1;
      }
    } else if (inClique) {
      int sizeClique=iRow-lastClique;
      for (int i=lastClique;i<iRow;i++) {
	clique_[i]=sizeClique;
	sizeClique--;
      }
      inClique=false;
    }
  }
  if (inClique) {
    int sizeClique=iRow-lastClique;
    for (int i=lastClique;i<iRow;i++) {
      clique_[i]=sizeClique;
      sizeClique--;
    }
  }
  //for (iRow=0;iRow<numberRows_;iRow++)
  //clique_[iRow]=0;
}
/* Factorize - filling in rowsDropped and returning number dropped */
int 
ClpCholeskyBase::factorize(const double * diagonal, int * rowsDropped) 
{
  const CoinBigIndex * columnStart = model_->clpMatrix()->getVectorStarts();
  const int * columnLength = model_->clpMatrix()->getVectorLengths();
  const int * row = model_->clpMatrix()->getIndices();
  const double * element = model_->clpMatrix()->getElements();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  const double * elementByRow = rowCopy_->getElements();
  int numberColumns=model_->clpMatrix()->getNumCols();
  //perturbation
  longDouble perturbation=model_->diagonalPerturbation()*model_->diagonalNorm();
  //perturbation=perturbation*perturbation*100000000.0;
  if (perturbation>1.0) {
    //if (model_->model()->logLevel()&4) 
    //std::cout <<"large perturbation "<<perturbation<<std::endl;
    perturbation=sqrt(perturbation);;
    perturbation=1.0;
  }
  int iRow;
  int iColumn;
  longDouble * work = workDouble_;
  CoinZeroN(work,numberRows_);
  int newDropped=0;
  double largest=1.0;
  double smallest=COIN_DBL_MAX;
  int numberDense=0;
  if (!doKKT_) {
    const double * diagonalSlack = diagonal + numberColumns;
    if (dense_)
      numberDense=dense_->numberRows();
    if (whichDense_) {
      longDouble * denseDiagonal = dense_->diagonal();
      longDouble * dense = denseColumn_;
      int iDense=0;
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (whichDense_[iColumn]) {
	  CoinZeroN(dense,numberRows_);
	  CoinBigIndex start=columnStart[iColumn];
	  CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	  if (diagonal[iColumn]) {
	    denseDiagonal[iDense++]=1.0/diagonal[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      dense[jRow] = element[j];
	    }
	  } else {
	    denseDiagonal[iDense++]=1.0;
	  }
	  dense += numberRows_;
	}
      }
    }
    longDouble delta2 = model_->delta(); // add delta*delta to diagonal
    delta2 *= delta2;
    // largest in initial matrix
    double largest2=1.0e-20;
    for (iRow=0;iRow<numberRows_;iRow++) {
      longDouble * put = sparseFactor_+choleskyStart_[iRow];
      int * which = choleskyRow_+indexStart_[iRow];
      int iOriginalRow = permute_[iRow];
      int number = choleskyStart_[iRow+1]-choleskyStart_[iRow];
      if (!rowLength[iOriginalRow])
	rowsDropped_[iOriginalRow]=1;
      if (!rowsDropped_[iOriginalRow]) {
	CoinBigIndex startRow=rowStart[iOriginalRow];
	CoinBigIndex endRow=rowStart[iOriginalRow]+rowLength[iOriginalRow];
	work[iRow] = diagonalSlack[iOriginalRow]+delta2;
	for (CoinBigIndex k=startRow;k<endRow;k++) {
	  int iColumn=column[k];
	  if (!whichDense_||!whichDense_[iColumn]) {
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    longDouble multiplier = diagonal[iColumn]*elementByRow[k];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow=row[j];
	      int jNewRow = permuteInverse_[jRow];
	      if (jNewRow>=iRow&&!rowsDropped_[jRow]) {
		longDouble value=element[j]*multiplier;
		work[jNewRow] += value;
	      }
	    }
	  }
	}
	diagonal_[iRow]=work[iRow];
	largest2 = CoinMax(largest2,fabs(work[iRow]));
	work[iRow]=0.0;
	int j;
	for (j=0;j<number;j++) {
	  int jRow = which[j];
	  put[j]=work[jRow];
	  largest2 = CoinMax(largest2,fabs(work[jRow]));
	  work[jRow]=0.0;
	}
      } else {
	// dropped
	diagonal_[iRow]=1.0;
	int j;
	for (j=1;j<number;j++) {
	  put[j]=0.0;
	}
      }
    }
    //check sizes
    largest2*=1.0e-20;
    largest = CoinMin(largest2,1.0e-11);
    int numberDroppedBefore=0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      int dropped=rowsDropped_[iRow];
      // Move to int array
      rowsDropped[iRow]=dropped;
      if (!dropped) {
	longDouble diagonal = diagonal_[iRow];
	if (diagonal>largest2) {
	  diagonal_[iRow]=diagonal+perturbation;
	} else {
	  diagonal_[iRow]=diagonal+perturbation;
	  rowsDropped[iRow]=2;
	  numberDroppedBefore++;
	  //printf("dropped - small diagonal %g\n",diagonal);
	} 
      } 
    }
    doubleParameters_[10]=CoinMax(1.0e-20,largest);
    integerParameters_[20]=0;
    doubleParameters_[3]=0.0;
    doubleParameters_[4]=COIN_DBL_MAX;
    integerParameters_[34]=0; // say all must be positive
    factorizePart2(rowsDropped);
    newDropped=integerParameters_[20]+numberDroppedBefore;
    largest=doubleParameters_[3];
    smallest=doubleParameters_[4];
    if (model_->messageHandler()->logLevel()>1) 
      std::cout<<"Cholesky - largest "<<largest<<" smallest "<<smallest<<std::endl;
    choleskyCondition_=largest/smallest;
    if (whichDense_) {
      for (int i=0;i<numberRows_;i++) {
	assert (diagonal_[i]>=0.0);
	diagonal_[i]=sqrt(diagonal_[i]);
      }
      // Update dense columns (just L)
      // Zero out dropped rows
      int i;
      for (i=0;i<numberDense;i++) {
	longDouble * a = denseColumn_+i*numberRows_;
	for (int j=0;j<numberRows_;j++) {
	  if (rowsDropped[j])
	    a[j]=0.0;
	}
	solveLong(a,1);
      }
      dense_->resetRowsDropped();
      longDouble * denseBlob = dense_->aMatrix();
      longDouble * denseDiagonal = dense_->diagonal();
      // Update dense matrix
      for (i=0;i<numberDense;i++) {
	const longDouble * a = denseColumn_+i*numberRows_;
	// do diagonal
	longDouble value = denseDiagonal[i];
	const longDouble * b = denseColumn_+i*numberRows_;
	for (int k=0;k<numberRows_;k++) 
	  value += a[k]*b[k];
	denseDiagonal[i]=value;
	for (int j=i+1;j<numberDense;j++) {
	  longDouble value = 0.0;
	  const longDouble * b = denseColumn_+j*numberRows_;
	  for (int k=0;k<numberRows_;k++) 
	    value += a[k]*b[k];
	  *denseBlob=value;
	  denseBlob++;
	}
      }
      // dense cholesky (? long double)
      int * dropped = new int [numberDense];
      dense_->factorizePart2(dropped);
      delete [] dropped;
    }
    // try allowing all every time
    //printf("trying ?\n");
    //for (iRow=0;iRow<numberRows_;iRow++) {
    //rowsDropped[iRow]=0;
    //rowsDropped_[iRow]=0;
    //}
    bool cleanCholesky;
    //if (model_->numberIterations()<20||(model_->numberIterations()&1)==0) 
    if (model_->numberIterations()<2000) 
      cleanCholesky=true;
    else 
      cleanCholesky=false;
    if (cleanCholesky) {
      //drop fresh makes some formADAT easier
      if (newDropped||numberRowsDropped_) {
	newDropped=0;
	for (int i=0;i<numberRows_;i++) {
	  char dropped = rowsDropped[i];
	  rowsDropped_[i]=dropped;
	  rowsDropped_[i]=0;
	  if (dropped==2) {
	    //dropped this time
	    rowsDropped[newDropped++]=i;
	    rowsDropped_[i]=0;
	  } 
	} 
	numberRowsDropped_=newDropped;
	newDropped=-(2+newDropped);
      } 
    } else {
      if (newDropped) {
	newDropped=0;
	for (int i=0;i<numberRows_;i++) {
	  char dropped = rowsDropped[i];
	  rowsDropped_[i]=dropped;
	  if (dropped==2) {
	    //dropped this time
	    rowsDropped[newDropped++]=i;
	    rowsDropped_[i]=1;
	  } 
	} 
      } 
      numberRowsDropped_+=newDropped;
      if (numberRowsDropped_&&0) {
	std::cout <<"Rank "<<numberRows_-numberRowsDropped_<<" ( "<<
          numberRowsDropped_<<" dropped)";
	if (newDropped) {
	  std::cout<<" ( "<<newDropped<<" dropped this time)";
	} 
	std::cout<<std::endl;
      } 
    }
  } else {
    //KKT
    CoinPackedMatrix * quadratic = NULL;
    ClpQuadraticObjective * quadraticObj = 
      (dynamic_cast< ClpQuadraticObjective*>(model_->objectiveAsObject()));
    if (quadraticObj) 
      quadratic = quadraticObj->quadraticObjective();
    // matrix
    int numberRowsModel = model_->numberRows();
    int numberColumns = model_->numberColumns();
    int numberTotal = numberColumns + numberRowsModel;
    // temp
    bool permuted=false;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (permute_[iRow]!=iRow) {
	permuted=true;
	break;
      }
    }
    // but fake it
    for (iRow=0;iRow<numberRows_;iRow++) {
      //permute_[iRow]=iRow; // force no permute
      //permuteInverse_[permute_[iRow]]=iRow;
    }
    if (permuted) {
      longDouble delta2 = model_->delta(); // add delta*delta to bottom
      delta2 *= delta2;
      // Need to permute - ugly
      if (!quadratic) {
	for (iRow=0;iRow<numberRows_;iRow++) {
	  longDouble * put = sparseFactor_+choleskyStart_[iRow];
	  int * which = choleskyRow_+indexStart_[iRow];
	  int iOriginalRow = permute_[iRow];
	  if (iOriginalRow<numberColumns) {
	    iColumn=iOriginalRow;
	    longDouble value = diagonal[iColumn];
	    if (fabs(value)>1.0e-100) {
	      value = 1.0/value;
	      largest = CoinMax(largest,fabs(value));
	      diagonal_[iRow] = -value;
	      CoinBigIndex start=columnStart[iColumn];
	      CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	      for (CoinBigIndex j=start;j<end;j++) {
		int kRow = row[j]+numberTotal;
		kRow=permuteInverse_[kRow];
		if (kRow>iRow) {
		  work[kRow]=element[j];
		  largest = CoinMax(largest,fabs(element[j]));
		}
	      }
	    } else {
	      diagonal_[iRow] = -value;
	    }
	  } else if (iOriginalRow<numberTotal) {
	    longDouble value = diagonal[iOriginalRow];
	    if (fabs(value)>1.0e-100) {
	      value = 1.0/value;
	      largest = CoinMax(largest,fabs(value));
	    } else {
	      value = 1.0e100;
	    }
	    diagonal_[iRow] = -value;
	    int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	    if (kRow>iRow) 
	      work[kRow]=-1.0;
	  } else {
	    diagonal_[iRow]=delta2;
	    int kRow = iOriginalRow-numberTotal;
	    CoinBigIndex start=rowStart[kRow];
	    CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow = column[j];
	      int jNewRow = permuteInverse_[jRow];
	      if (jNewRow>iRow) {
		work[jNewRow]=elementByRow[j];
		largest = CoinMax(largest,fabs(elementByRow[j]));
	      }
	    }
	    // slack - should it be permute
	    kRow = permuteInverse_[kRow+numberColumns];
	    if (kRow>iRow)
	      work[kRow]=-1.0;
	  }
	  CoinBigIndex j;
	  int number = choleskyStart_[iRow+1]-choleskyStart_[iRow];
	  for (j=0;j<number;j++) {
	    int jRow = which[j];
	    put[j]=work[jRow];
	    work[jRow]=0.0;
	  }
	}
      } else {
	// quadratic
	const int * columnQuadratic = quadratic->getIndices();
	const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	const int * columnQuadraticLength = quadratic->getVectorLengths();
	const double * quadraticElement = quadratic->getElements();
	for (iRow=0;iRow<numberRows_;iRow++) {
	  longDouble * put = sparseFactor_+choleskyStart_[iRow];
	  int * which = choleskyRow_+indexStart_[iRow];
	  int iOriginalRow = permute_[iRow];
	  if (iOriginalRow<numberColumns) {
	    CoinBigIndex j;
	    iColumn=iOriginalRow;
	    longDouble value = diagonal[iColumn];
	    if (fabs(value)>1.0e-100) {
	      value = 1.0/value;
	      for (j=columnQuadraticStart[iColumn];
		   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
		int jColumn = columnQuadratic[j];
		int jNewColumn = permuteInverse_[jColumn];
		if (jNewColumn>iRow) {
		  work[jNewColumn]=-quadraticElement[j];
		} else if (iColumn==jColumn) {
		  value += quadraticElement[j];
		}
	      }
	      largest = CoinMax(largest,fabs(value));
	      diagonal_[iRow] = -value;
	      CoinBigIndex start=columnStart[iColumn];
	      CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	      for (j=start;j<end;j++) {
		int kRow = row[j]+numberTotal;
		kRow=permuteInverse_[kRow];
		if (kRow>iRow) {
		  work[kRow]=element[j];
		  largest = CoinMax(largest,fabs(element[j]));
		}
	      }
	    } else {
	      diagonal_[iRow] = -value;
	    }
	  } else if (iOriginalRow<numberTotal) {
	    longDouble value = diagonal[iOriginalRow];
	    if (fabs(value)>1.0e-100) {
	      value = 1.0/value;
	      largest = CoinMax(largest,fabs(value));
	    } else {
	      value = 1.0e100;
	    }
	    diagonal_[iRow] = -value;
	    int kRow = permuteInverse_[iOriginalRow+numberRowsModel];
	    if (kRow>iRow) 
	      work[kRow]=-1.0;
	  } else {
	    diagonal_[iRow]=delta2;
	    int kRow = iOriginalRow-numberTotal;
	    CoinBigIndex start=rowStart[kRow];
	    CoinBigIndex end=rowStart[kRow]+rowLength[kRow];
	    for (CoinBigIndex j=start;j<end;j++) {
	      int jRow = column[j];
	      int jNewRow = permuteInverse_[jRow];
	      if (jNewRow>iRow) {
		work[jNewRow]=elementByRow[j];
		largest = CoinMax(largest,fabs(elementByRow[j]));
	      }
	    }
	    // slack - should it be permute
	    kRow = permuteInverse_[kRow+numberColumns];
	    if (kRow>iRow)
	      work[kRow]=-1.0;
	  }
	  CoinBigIndex j;
	  int number = choleskyStart_[iRow+1]-choleskyStart_[iRow];
	  for (j=0;j<number;j++) {
	    int jRow = which[j];
	    put[j]=work[jRow];
	    work[jRow]=0.0;
	  }
	  for (j=0;j<numberRows_;j++)
	    assert (!work[j]);
	}
      }
    } else {
      if (!quadratic) {
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  longDouble * put = sparseFactor_+choleskyStart_[iColumn];
	  int * which = choleskyRow_+indexStart_[iColumn];
	  longDouble value = diagonal[iColumn];
	  if (fabs(value)>1.0e-100) {
	    value = 1.0/value;
	    largest = CoinMax(largest,fabs(value));
	    diagonal_[iColumn] = -value;
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (CoinBigIndex j=start;j<end;j++) {
	      //choleskyRow_[numberElements]=row[j]+numberTotal;
	      //sparseFactor_[numberElements++]=element[j];
	      work[row[j]+numberTotal]=element[j];
	      largest = CoinMax(largest,fabs(element[j]));
	    }
	  } else {
	    diagonal_[iColumn] = -value;
	  }
	  CoinBigIndex j;
	  int number = choleskyStart_[iColumn+1]-choleskyStart_[iColumn];
	  for (j=0;j<number;j++) {
	    int jRow = which[j];
	    put[j]=work[jRow];
	    work[jRow]=0.0;
	  }
	}
      } else {
	// Quadratic
	const int * columnQuadratic = quadratic->getIndices();
	const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
	const int * columnQuadraticLength = quadratic->getVectorLengths();
	const double * quadraticElement = quadratic->getElements();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  longDouble * put = sparseFactor_+choleskyStart_[iColumn];
	  int * which = choleskyRow_+indexStart_[iColumn];
	  int number = choleskyStart_[iColumn+1]-choleskyStart_[iColumn];
	  longDouble value = diagonal[iColumn];
	  CoinBigIndex j;
	  if (fabs(value)>1.0e-100) {
	    value = 1.0/value;
	    for (j=columnQuadraticStart[iColumn];
		 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	      int jColumn = columnQuadratic[j];
	      if (jColumn>iColumn) {
		work[jColumn]=-quadraticElement[j];
	      } else if (iColumn==jColumn) {
		value += quadraticElement[j];
	      }
	    }
	    largest = CoinMax(largest,fabs(value));
	    diagonal_[iColumn] = -value;
	    CoinBigIndex start=columnStart[iColumn];
	    CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	    for (j=start;j<end;j++) {
	      work[row[j]+numberTotal]=element[j];
	      largest = CoinMax(largest,fabs(element[j]));
	    }
	    for (j=0;j<number;j++) {
	      int jRow = which[j];
	      put[j]=work[jRow];
	      work[jRow]=0.0;
	    }
	  } else {
	    value = 1.0e100;
	    diagonal_[iColumn] = -value;
	    for (j=0;j<number;j++) {
	      int jRow = which[j];
	      put[j]=work[jRow];
	    }
	  }
	}
      }
      // slacks
      for (iColumn=numberColumns;iColumn<numberTotal;iColumn++) {
	longDouble * put = sparseFactor_+choleskyStart_[iColumn];
	int * which = choleskyRow_+indexStart_[iColumn];
	longDouble value = diagonal[iColumn];
	if (fabs(value)>1.0e-100) {
	  value = 1.0/value;
	  largest = CoinMax(largest,fabs(value));
	} else {
	  value = 1.0e100;
	}
	diagonal_[iColumn] = -value;
	work[iColumn-numberColumns+numberTotal]=-1.0;
	CoinBigIndex j;
	int number = choleskyStart_[iColumn+1]-choleskyStart_[iColumn];
	for (j=0;j<number;j++) {
	  int jRow = which[j];
	  put[j]=work[jRow];
	  work[jRow]=0.0;
	}
      }
      // Finish diagonal
      longDouble delta2 = model_->delta(); // add delta*delta to bottom
      delta2 *= delta2;
      for (iRow=0;iRow<numberRowsModel;iRow++) {
	longDouble * put = sparseFactor_+choleskyStart_[iRow+numberTotal];
	diagonal_[iRow+numberTotal]=delta2;
	CoinBigIndex j;
	int number = choleskyStart_[iRow+numberTotal+1]-choleskyStart_[iRow+numberTotal];
	for (j=0;j<number;j++) {
	  put[j]=0.0;
	}
      }
    }
    //check sizes
    largest*=1.0e-20;
    largest = CoinMin(largest,1.0e-11);
    doubleParameters_[10]=CoinMax(1.0e-20,largest);
    integerParameters_[20]=0;
    doubleParameters_[3]=0.0;
    doubleParameters_[4]=COIN_DBL_MAX;
    // Set up LDL cutoff
    integerParameters_[34]=numberTotal;
    // KKT
    int * rowsDropped2 = new int[numberRows_];
    CoinZeroN(rowsDropped2,numberRows_);
    factorizePart2(rowsDropped2);
    newDropped=integerParameters_[20];
    largest=doubleParameters_[3];
    smallest=doubleParameters_[4];
    if (model_->messageHandler()->logLevel()>1) 
      std::cout<<"Cholesky - largest "<<largest<<" smallest "<<smallest<<std::endl;
    choleskyCondition_=largest/smallest;
    // Should save adjustments in ..R_
    int n1=0,n2=0;
    double * primalR = model_->primalR();
    double * dualR = model_->dualR();
    for (iRow=0;iRow<numberTotal;iRow++) { 
      if (rowsDropped2[iRow]) {
	n1++;
	//printf("row region1 %d dropped\n",iRow);
	//rowsDropped_[iRow]=1;
	rowsDropped_[iRow]=0;
	primalR[iRow]=doubleParameters_[20];
      } else {
	rowsDropped_[iRow]=0;
	primalR[iRow]=0.0;
      }
    }
    for (;iRow<numberRows_;iRow++) {
      if (rowsDropped2[iRow]) {
	n2++;
	//printf("row region2 %d dropped\n",iRow);
	//rowsDropped_[iRow]=1;
	rowsDropped_[iRow]=0;
	dualR[iRow-numberTotal]=doubleParameters_[34];
      } else {
	rowsDropped_[iRow]=0;
	dualR[iRow-numberTotal]=0.0;
      }
    }
  }
  status_=0;
  return newDropped;
}
/* Factorize - filling in rowsDropped and returning number dropped
   in integerParam.
*/
void 
ClpCholeskyBase::factorizePart2(int * rowsDropped) 
{
  longDouble largest=doubleParameters_[3];
  longDouble smallest=doubleParameters_[4];
  int numberDropped=integerParameters_[20];
  // probably done before
  largest=0.0;
  smallest=COIN_DBL_MAX;
  numberDropped=0;
  double dropValue = doubleParameters_[10];
  int firstPositive=integerParameters_[34];
  longDouble * d = ClpCopyOfArray(diagonal_,numberRows_);
  int iRow;
  // minimum size before clique done
  //#define MINCLIQUE INT_MAX
#define MINCLIQUE 3
  longDouble * work = workDouble_;
  CoinBigIndex * first = workInteger_;
  
  for (iRow=0;iRow<numberRows_;iRow++) {
    link_[iRow]=-1;
    work[iRow]=0.0;
    first[iRow]=choleskyStart_[iRow];
  }

  int lastClique=-1;
  bool inClique=false;
  bool newClique=false;
  bool endClique=false;
  int lastRow=0;
  int cliqueSize=0;
  CoinBigIndex cliquePointer=0;
  int nextRow2=-1;
  
  for (iRow=0;iRow<firstDense_+1;iRow++) {
    if (iRow<firstDense_) {
      endClique=false;
      if (clique_[iRow]>0) {
	// this is in a clique
	inClique=true;
	if (clique_[iRow]>lastClique) {
	  // new Clique
	  newClique=true;
	  // If we have clique going then signal to do old one
	  endClique=(lastClique>0);
	} else {
	  // Still in clique
	  newClique=false;
	}
      } else {
	// not in clique
	inClique=false;
	newClique=false;
	// If we have clique going then signal to do old one
	endClique=(lastClique>0);
      }
      lastClique=clique_[iRow];
    } else if (inClique) {
      // Finish off
      endClique=true;
    } else {
      break;
    }
    if (endClique) {
      // We have just finished updating a clique - do block pivot and clean up
      int jRow;
      for ( jRow=lastRow;jRow<iRow;jRow++) {
	int jCount = jRow-lastRow;
	longDouble diagonalValue = diagonal_[jRow];
	CoinBigIndex start=choleskyStart_[jRow];
	CoinBigIndex end=choleskyStart_[jRow+1];
	for (int kRow=lastRow;kRow<jRow;kRow++) {
	  jCount--;
	  CoinBigIndex get = choleskyStart_[kRow]+jCount;
	  longDouble a_jk = sparseFactor_[get];
	  longDouble value1 = d[kRow]*a_jk;
	  diagonalValue -= a_jk*value1;
	  for (CoinBigIndex j=start;j<end;j++)
	    sparseFactor_[j] -= value1*sparseFactor_[++get];
	}
	// check
	int originalRow = permute_[jRow];
	if (originalRow<firstPositive) {
	  // must be negative
	  if (diagonalValue<=-dropValue) {
	    smallest = CoinMin(smallest,-diagonalValue);
	    largest = CoinMax(largest,-diagonalValue);
	    d[jRow]=diagonalValue;
	    diagonalValue = 1.0/diagonalValue;
	  } else {
	    rowsDropped[originalRow]=2;
	    d[jRow]=-1.0e100;
	    diagonalValue=0.0;
	    integerParameters_[20]++;
	  }
	} else {
	  // must be positive
	  if (diagonalValue>=dropValue) {
	    smallest = CoinMin(smallest,diagonalValue);
	    largest = CoinMax(largest,diagonalValue);
	    d[jRow]=diagonalValue;
	    diagonalValue = 1.0/diagonalValue;
	  } else {
	    rowsDropped[originalRow]=2;
	    d[jRow]=1.0e100;
	    diagonalValue=0.0;
	    integerParameters_[20]++;
	  }
	}
	diagonal_[jRow]=diagonalValue;
	for (CoinBigIndex j=start;j<end;j++) {
	  sparseFactor_[j] *= diagonalValue;
	}
      }
      if (nextRow2>=0) {
	for ( jRow=lastRow;jRow<iRow-1;jRow++) {
	  link_[jRow]=jRow+1;
	}
	link_[iRow-1]=link_[nextRow2];
	link_[nextRow2]=lastRow;
      }
    }
    if (iRow==firstDense_)
      break; // we were just cleaning up
    if (newClique) {
      // initialize new clique
      lastRow=iRow;
      cliquePointer=choleskyStart_[iRow];
      cliqueSize=choleskyStart_[iRow+1]-cliquePointer+1;
    }
    // for each column L[*,kRow] that affects L[*,iRow] 
    longDouble diagonalValue=diagonal_[iRow];
    int nextRow = link_[iRow];
    int kRow=0;
    while (1) {
      kRow=nextRow;
      if (kRow<0)
	break; // out of loop
      nextRow=link_[kRow];
      // Modify by outer product of L[*,irow] by L[*,krow] from first
      CoinBigIndex k=first[kRow];
      CoinBigIndex end=choleskyStart_[kRow+1];
      assert(k<end);
      longDouble a_ik=sparseFactor_[k++];
      longDouble value1=d[kRow]*a_ik;
      // update first
      first[kRow]=k;
      diagonalValue -= value1*a_ik;
      CoinBigIndex offset = indexStart_[kRow]-choleskyStart_[kRow];
      int jRow = choleskyRow_[k+offset];
      if (k<end) {
	if (clique_[kRow]<MINCLIQUE) {
	  link_[kRow]=link_[jRow];
	  link_[jRow]=kRow;
	  for (;k<end;k++) {
	    int jRow = choleskyRow_[k+offset];
	    work[jRow] += sparseFactor_[k]*value1;
	  }
	} else {
	  // Clique
	  CoinBigIndex currentIndex = k+offset;
	  int linkSave=link_[jRow];
	  link_[jRow]=kRow;
	  work[kRow]=value1; // ? or a_jk
	  int last = kRow+clique_[kRow];
	  for (int kkRow=kRow+1;kkRow<last;kkRow++) {
	    CoinBigIndex j=first[kkRow];
	    //int iiRow = choleskyRow_[j+indexStart_[kkRow]-choleskyStart_[kkRow]];
	    longDouble a = sparseFactor_[j];
	    longDouble dValue = d[kkRow]*a;
	    diagonalValue -= a*dValue;
	    work[kkRow]=dValue;
	    first[kkRow]++;
	    link_[kkRow-1]=kkRow;
	  }
	  nextRow = link_[last-1];
	  link_[last-1]=linkSave;
	  int length = end-k;
	  for (int i=0;i<length;i++) {
	    int lRow = choleskyRow_[currentIndex++];
	    longDouble t0 = work[lRow];
	    for (int kkRow=kRow;kkRow<last;kkRow++) {
	      CoinBigIndex j = first[kkRow]+i;
	      t0 += work[kkRow]*sparseFactor_[j];
	    }
	    work[lRow]=t0;
	  }
	}
      }
    }
    // Now apply
    if (inClique) {
      // in clique
      diagonal_[iRow]=diagonalValue;
      CoinBigIndex start=choleskyStart_[iRow];
      CoinBigIndex end=choleskyStart_[iRow+1];
      CoinBigIndex currentIndex = indexStart_[iRow];
      nextRow2=-1;
      CoinBigIndex get=start+clique_[iRow]-1;
      if (get<end) {
	nextRow2 = choleskyRow_[currentIndex+get-start];
	first[iRow]=get;
      }
      for (CoinBigIndex j=start;j<end;j++) {
	int kRow = choleskyRow_[currentIndex++];
	sparseFactor_[j] -= work[kRow]; // times?
	work[kRow]=0.0;
      }
    } else {
      // not in clique
      int originalRow = permute_[iRow];
      if (originalRow<firstPositive) {
	// must be negative
	if (diagonalValue<=-dropValue) {
	  smallest = CoinMin(smallest,-diagonalValue);
	  largest = CoinMax(largest,-diagonalValue);
	  d[iRow]=diagonalValue;
	  diagonalValue = 1.0/diagonalValue;
	} else {
	  rowsDropped[originalRow]=2;
	  d[iRow]=-1.0e100;
	  diagonalValue=0.0;
	  integerParameters_[20]++;
	}
      } else {
	// must be positive
	if (diagonalValue>=dropValue) {
	  smallest = CoinMin(smallest,diagonalValue);
	  largest = CoinMax(largest,diagonalValue);
	  d[iRow]=diagonalValue;
	  diagonalValue = 1.0/diagonalValue;
	} else {
	  rowsDropped[originalRow]=2;
	  d[iRow]=1.0e100;
	  diagonalValue=0.0;
	  integerParameters_[20]++;
	}
      }
      diagonal_[iRow]=diagonalValue;
      CoinBigIndex offset = indexStart_[iRow]-choleskyStart_[iRow];
      CoinBigIndex start = choleskyStart_[iRow];
      CoinBigIndex end = choleskyStart_[iRow+1];
      assert (first[iRow]==start);
      if (start<end) {
	int nextRow = choleskyRow_[start+offset];
	link_[iRow]=link_[nextRow];
	link_[nextRow]=iRow;
	for (CoinBigIndex j=start;j<end;j++) {
	  int jRow = choleskyRow_[j+offset];
	  longDouble value = sparseFactor_[j]-work[jRow];
	  work[jRow]=0.0;
	  sparseFactor_[j]= diagonalValue*value;
	}
      }
    }
  }
  if (firstDense_<numberRows_) {
    // do dense
    // update dense part
    updateDense(d,work,first);
    ClpCholeskyDense dense;
    // just borrow space
    int nDense = numberRows_-firstDense_;
    if (doKKT_) {
      for (iRow=firstDense_;iRow<numberRows_;iRow++) {
	int originalRow=permute_[iRow];
	if (originalRow>=firstPositive) {
	  firstPositive=iRow-firstDense_;
	  break;
	}
      }
    }
    dense.reserveSpace(this,nDense);
    int * dropped = new int[nDense];
    memset(dropped,0,nDense*sizeof(int));
    dense.setDoubleParameter(3,largest);
    dense.setDoubleParameter(4,smallest);
    dense.setDoubleParameter(10,dropValue);
    dense.setIntegerParameter(20,0);
    dense.setIntegerParameter(34,firstPositive);
    dense.factorizePart2(dropped);
    largest=dense.getDoubleParameter(3);
    smallest=dense.getDoubleParameter(4);
    integerParameters_[20]+=dense.getIntegerParameter(20);
    for (iRow=firstDense_;iRow<numberRows_;iRow++) {
      int originalRow=permute_[iRow];
      rowsDropped[originalRow]=dropped[iRow-firstDense_];
    }
    delete [] dropped;
  }
  delete [] d;
  doubleParameters_[3]=largest;
  doubleParameters_[4]=smallest;
  return;
}
// Updates dense part (broken out for profiling)
void ClpCholeskyBase::updateDense(longDouble * d, longDouble * work, int * first)
{
  for (int iRow=0;iRow<firstDense_;iRow++) {
    CoinBigIndex start=first[iRow];
    CoinBigIndex end=choleskyStart_[iRow+1];
    if (start<end) {
      CoinBigIndex offset = indexStart_[iRow]-choleskyStart_[iRow];
      if (clique_[iRow]<2) {
	longDouble dValue=d[iRow];
	for (CoinBigIndex k=start;k<end;k++) {
	  int kRow = choleskyRow_[k+offset]; 
	  assert(kRow>=firstDense_);
	  longDouble a_ik=sparseFactor_[k];
	  longDouble value1=dValue*a_ik;
	  diagonal_[kRow] -= value1*a_ik;
	  CoinBigIndex base = choleskyStart_[kRow]-kRow-1;
	  for (CoinBigIndex j=k+1;j<end;j++) {
	    int jRow=choleskyRow_[j+offset];
	    longDouble a_jk = sparseFactor_[j];
	    sparseFactor_[base+jRow] -= a_jk*value1;
	  }
	}
      } else if (clique_[iRow]<3) {
	// do as pair
	longDouble dValue0=d[iRow];
	longDouble dValue1=d[iRow+1];
	int offset1 = first[iRow+1]-start;
	// skip row
	iRow++;
	for (CoinBigIndex k=start;k<end;k++) {
	  int kRow = choleskyRow_[k+offset]; 
	  assert(kRow>=firstDense_);
	  longDouble a_ik0=sparseFactor_[k];
	  longDouble value0=dValue0*a_ik0;
	  longDouble a_ik1=sparseFactor_[k+offset1];
	  longDouble value1=dValue1*a_ik1;
	  diagonal_[kRow] -= value0*a_ik0 + value1*a_ik1;
	  CoinBigIndex base = choleskyStart_[kRow]-kRow-1;
	  for (CoinBigIndex j=k+1;j<end;j++) {
	    int jRow=choleskyRow_[j+offset];
	    longDouble a_jk0 = sparseFactor_[j];
	    longDouble a_jk1 = sparseFactor_[j+offset1];
	    sparseFactor_[base+jRow] -= a_jk0*value0+a_jk1*value1;
	  }
	}
#define MANY_REGISTERS
#ifdef MANY_REGISTERS
      } else if (clique_[iRow]==3) {
#else
      } else {
#endif
	// do as clique
	// maybe later get fancy on big cliques and do transpose copy
	// seems only worth going to 3 on Intel
	longDouble dValue0=d[iRow];
	longDouble dValue1=d[iRow+1];
	longDouble dValue2=d[iRow+2];
	// get offsets and skip rows
	int offset1 = first[++iRow]-start;
	int offset2 = first[++iRow]-start;
	for (CoinBigIndex k=start;k<end;k++) {
	  int kRow = choleskyRow_[k+offset]; 
	  assert(kRow>=firstDense_);
	  double diagonalValue=diagonal_[kRow];
	  longDouble a_ik0=sparseFactor_[k];
	  longDouble value0=dValue0*a_ik0;
	  longDouble a_ik1=sparseFactor_[k+offset1];
	  longDouble value1=dValue1*a_ik1;
	  longDouble a_ik2=sparseFactor_[k+offset2];
	  longDouble value2=dValue2*a_ik2;
	  CoinBigIndex base = choleskyStart_[kRow]-kRow-1;
	  diagonal_[kRow] = diagonalValue - value0*a_ik0 - value1*a_ik1 - value2*a_ik2;
	  for (CoinBigIndex j=k+1;j<end;j++) {
	    int jRow=choleskyRow_[j+offset];
	    longDouble a_jk0 = sparseFactor_[j];
	    longDouble a_jk1 = sparseFactor_[j+offset1];
	    longDouble a_jk2 = sparseFactor_[j+offset2];
	    sparseFactor_[base+jRow] -= a_jk0*value0+a_jk1*value1+a_jk2*value2;
	  }
	}
#ifdef MANY_REGISTERS
      } else {
	// do as clique
	// maybe later get fancy on big cliques and do transpose copy
	// maybe only worth going to 3 on Intel (but may have hidden registers)
	longDouble dValue0=d[iRow];
	longDouble dValue1=d[iRow+1];
	longDouble dValue2=d[iRow+2];
	longDouble dValue3=d[iRow+3];
	// get offsets and skip rows
	int offset1 = first[++iRow]-start;
	int offset2 = first[++iRow]-start;
	int offset3 = first[++iRow]-start;
	for (CoinBigIndex k=start;k<end;k++) {
	  int kRow = choleskyRow_[k+offset]; 
	  assert(kRow>=firstDense_);
	  double diagonalValue=diagonal_[kRow];
	  longDouble a_ik0=sparseFactor_[k];
	  longDouble value0=dValue0*a_ik0;
	  longDouble a_ik1=sparseFactor_[k+offset1];
	  longDouble value1=dValue1*a_ik1;
	  longDouble a_ik2=sparseFactor_[k+offset2];
	  longDouble value2=dValue2*a_ik2;
	  longDouble a_ik3=sparseFactor_[k+offset3];
	  longDouble value3=dValue3*a_ik3;
	  CoinBigIndex base = choleskyStart_[kRow]-kRow-1;
	  diagonal_[kRow] = diagonalValue - (value0*a_ik0 + value1*a_ik1 + value2*a_ik2+value3*a_ik3);
	  for (CoinBigIndex j=k+1;j<end;j++) {
	    int jRow=choleskyRow_[j+offset];
	    longDouble a_jk0 = sparseFactor_[j];
	    longDouble a_jk1 = sparseFactor_[j+offset1];
	    longDouble a_jk2 = sparseFactor_[j+offset2];
	    longDouble a_jk3 = sparseFactor_[j+offset3];
	    sparseFactor_[base+jRow] -= a_jk0*value0+a_jk1*value1+a_jk2*value2+a_jk3*value3;
	  }
	}
#endif
      }
    }
  }
}
/* Uses factorization to solve. */
void 
ClpCholeskyBase::solve (double * region) 
{
  if (!whichDense_) {
    solve(region,3);
  } else {
    // dense columns
    int i;
    solve(region,1);
    // do change;
    int numberDense = dense_->numberRows();
    double * change = new double[numberDense];
    for (i=0;i<numberDense;i++) {
      const longDouble * a = denseColumn_+i*numberRows_;
      longDouble value =0.0;
      for (int iRow=0;iRow<numberRows_;iRow++) 
	value += a[iRow]*region[iRow];
      change[i]=value;
    }
    // solve
    dense_->solve(change);
    for (i=0;i<numberDense;i++) {
      const longDouble * a = denseColumn_+i*numberRows_;
      longDouble value = change[i];
      for (int iRow=0;iRow<numberRows_;iRow++) 
	region[iRow] -= value*a[iRow];
    }
    delete [] change;
    // and finish off
    solve(region,2);
  }
}
/* solve - 1 just first half, 2 just second half - 3 both.
   If 1 and 2 then diagonal has sqrt of inverse otherwise inverse
*/
void 
ClpCholeskyBase::solve(double * region, int type)
{
#ifdef CLP_DEBUG
  double * regionX=NULL;
  if (sizeof(longWork)!=sizeof(double)&&type==3) {
    regionX=ClpCopyOfArray(region,numberRows_);
  }
#endif
  longWork * work = (longWork *) workDouble_;
  int i;
  CoinBigIndex j;
  for (i=0;i<numberRows_;i++) {
    int iRow = permute_[i];
    work[i] = region[iRow];
  }
  switch (type) {
  case 1:
    for (i=0;i<numberRows_;i++) {
      longDouble value=work[i];
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	work[iRow] -= sparseFactor_[j]*value;
      }
    }
    for (i=0;i<numberRows_;i++) {
      int iRow = permute_[i];
      region[iRow]=work[i]*diagonal_[i];
    }
    break;
  case 2:
    for (i=numberRows_-1;i>=0;i--) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=work[i]*diagonal_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	value -= sparseFactor_[j]*work[iRow];
      }
      work[i]=value;
      int iRow = permute_[i];
      region[iRow]=value;
    }
    break;
  case 3:
    for (i=0;i<firstDense_;i++) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=work[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	work[iRow] -= sparseFactor_[j]*value;
      }
    }
    if (firstDense_<numberRows_) {
      // do dense
      ClpCholeskyDense dense;
      // just borrow space
      int nDense = numberRows_-firstDense_;
      dense.reserveSpace(this,nDense);
#if CLP_LONG_CHOLESKY!=1
      dense.solveLong(work+firstDense_);
#else
      dense.solveLongWork(work+firstDense_);
#endif
      for (i=numberRows_-1;i>=firstDense_;i--) {
	longDouble value=work[i];
	int iRow = permute_[i];
	region[iRow]=value;
      }
    }
    for (i=firstDense_-1;i>=0;i--) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=work[i]*diagonal_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	value -= sparseFactor_[j]*work[iRow];
      }
      work[i]=value;
      int iRow = permute_[i];
      region[iRow]=value;
    }
    break;
  }
#ifdef CLP_DEBUG
  if (regionX) {
    longDouble * work = workDouble_;
    int i;
    CoinBigIndex j;
    double largestO=0.0;
    for (i=0;i<numberRows_;i++) {
      largestO = CoinMax(largestO,fabs(regionX[i]));
    }
    for (i=0;i<numberRows_;i++) {
      int iRow = permute_[i];
      work[i] = regionX[iRow];
    }
    for (i=0;i<firstDense_;i++) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=work[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	work[iRow] -= sparseFactor_[j]*value;
      }
    }
    if (firstDense_<numberRows_) {
      // do dense
      ClpCholeskyDense dense;
      // just borrow space
      int nDense = numberRows_-firstDense_;
      dense.reserveSpace(this,nDense);
      dense.solveLong(work+firstDense_);
      for (i=numberRows_-1;i>=firstDense_;i--) {
	longDouble value=work[i];
	int iRow = permute_[i];
	regionX[iRow]=value;
      }
    }
    for (i=firstDense_-1;i>=0;i--) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=work[i]*diagonal_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	value -= sparseFactor_[j]*work[iRow];
      }
      work[i]=value;
      int iRow = permute_[i];
      regionX[iRow]=value;
    }
    double largest=0.0;
    double largestV=0.0;
    for (i=0;i<numberRows_;i++) {
      largest = CoinMax(largest,fabs(region[i]-regionX[i]));
      largestV = CoinMax(largestV,fabs(region[i]));
    }
    printf("largest difference %g, largest %g, largest original %g\n",
	   largest,largestV,largestO);
    delete [] regionX;
  }
#endif
}
/* solve - 1 just first half, 2 just second half - 3 both.
   If 1 and 2 then diagonal has sqrt of inverse otherwise inverse
*/
void 
ClpCholeskyBase::solveLong(longDouble * region, int type)
{
  int i;
  CoinBigIndex j;
  for (i=0;i<numberRows_;i++) {
    int iRow = permute_[i];
    workDouble_[i] = region[iRow];
  }
  switch (type) {
  case 1:
    for (i=0;i<numberRows_;i++) {
      longDouble value=workDouble_[i];
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	workDouble_[iRow] -= sparseFactor_[j]*value;
      }
    }
    for (i=0;i<numberRows_;i++) {
      int iRow = permute_[i];
      region[iRow]=workDouble_[i]*diagonal_[i];
    }
    break;
  case 2:
    for (i=numberRows_-1;i>=0;i--) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=workDouble_[i]*diagonal_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	value -= sparseFactor_[j]*workDouble_[iRow];
      }
      workDouble_[i]=value;
      int iRow = permute_[i];
      region[iRow]=value;
    }
    break;
  case 3:
    for (i=0;i<firstDense_;i++) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=workDouble_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	workDouble_[iRow] -= sparseFactor_[j]*value;
      }
    }
    if (firstDense_<numberRows_) {
      // do dense
      ClpCholeskyDense dense;
      // just borrow space
      int nDense = numberRows_-firstDense_;
      dense.reserveSpace(this,nDense);
      dense.solveLong(workDouble_+firstDense_);
      for (i=numberRows_-1;i>=firstDense_;i--) {
	longDouble value=workDouble_[i];
	int iRow = permute_[i];
	region[iRow]=value;
      }
    }
    for (i=firstDense_-1;i>=0;i--) {
      CoinBigIndex offset = indexStart_[i]-choleskyStart_[i];
      longDouble value=workDouble_[i]*diagonal_[i];
      for (j=choleskyStart_[i];j<choleskyStart_[i+1];j++) {
	int iRow = choleskyRow_[j+offset];
	value -= sparseFactor_[j]*workDouble_[iRow];
      }
      workDouble_[i]=value;
      int iRow = permute_[i];
      region[iRow]=value;
    }
    break;
  }
}

