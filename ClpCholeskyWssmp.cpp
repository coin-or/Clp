// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.



#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"

#include "ClpInterior.hpp"
#include "ClpCholeskyWssmp.hpp"
#include "ClpMessage.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpCholeskyWssmp::ClpCholeskyWssmp () 
  : ClpCholeskyBase(),
    sparseFactor_(NULL),
    choleskyStart_(NULL),
    choleskyRow_(NULL),
    sizeFactor_(0),
    rowCopy_(NULL)
{
  type_=12;
  memset(integerParameters_,0,64*sizeof(int));
  memset(doubleParameters_,0,64*sizeof(double));
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpCholeskyWssmp::ClpCholeskyWssmp (const ClpCholeskyWssmp & rhs) 
: ClpCholeskyBase(rhs)
{
  type_=rhs.type_;
  sparseFactor_ = ClpCopyOfArray(rhs.sparseFactor_,rhs.sizeFactor_);
  choleskyStart_ = ClpCopyOfArray(rhs.choleskyStart_,numberRows_+1);
  choleskyRow_ = ClpCopyOfArray(rhs.choleskyRow_,rhs.sizeFactor_);
  sizeFactor_=rhs.sizeFactor_;
  memcpy(integerParameters_,rhs.integerParameters_,64*sizeof(int));
  memcpy(doubleParameters_,rhs.doubleParameters_,64*sizeof(double));
  rowCopy_ = rhs.rowCopy_->clone();
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpCholeskyWssmp::~ClpCholeskyWssmp ()
{
  delete [] sparseFactor_;
  delete [] choleskyStart_;
  delete [] choleskyRow_;
  delete rowCopy_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpCholeskyWssmp &
ClpCholeskyWssmp::operator=(const ClpCholeskyWssmp& rhs)
{
  if (this != &rhs) {
    ClpCholeskyBase::operator=(rhs);
    delete [] sparseFactor_;
    delete [] choleskyStart_;
    delete [] choleskyRow_;
    sparseFactor_ = ClpCopyOfArray(rhs.sparseFactor_,rhs.sizeFactor_);
    choleskyStart_ = ClpCopyOfArray(rhs.choleskyStart_,numberRows_+1);
    choleskyRow_ = ClpCopyOfArray(rhs.choleskyRow_,rhs.sizeFactor_);
    sizeFactor_=rhs.sizeFactor_;
    delete rowCopy_;
    rowCopy_ = rhs.rowCopy_->clone();
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpCholeskyBase * ClpCholeskyWssmp::clone() const
{
  return new ClpCholeskyWssmp(*this);
}
  extern "C" void wssmp(int * n,
                        int * columnStart , int * rowIndex , double * element,
                        double * diagonal , int * perm , int * invp ,
                        double * rhs , int * ldb , int * nrhs ,
                        double * aux , int * naux ,
                        int   * mrp , int * iparm , double * dparm);
/* Orders rows and saves pointer to matrix.and model */
int 
ClpCholeskyWssmp::order(ClpInterior * model) 
{
  numberRows_ = model->numberRows();
  rowsDropped_ = new char [numberRows_];
  memset(rowsDropped_,0,numberRows_);
  numberRowsDropped_=0;
  model_=model;
  rowCopy_ = model->clpMatrix()->reverseOrderedCopy();
  // Space for starts
  choleskyStart_ = new CoinBigIndex[numberRows_+1];
  const CoinBigIndex * columnStart = model_->matrix()->getVectorStarts();
  const int * columnLength = model_->matrix()->getVectorLengths();
  const int * row = model_->matrix()->getIndices();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  // We need two arrays for counts
  int * which = new int [numberRows_];
  int * used = new int[numberRows_];
  CoinZeroN(used,numberRows_);
  int iRow;
  sizeFactor_=0;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int number=0;
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	CoinBigIndex start=columnStart[iColumn];
	CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	for (CoinBigIndex j=start;j<end;j++) {
	  int jRow=row[j];
	  if (jRow>=iRow&&!rowsDropped_[jRow]) {
	    if (!used[jRow]) {
	      used[jRow]=1;
	      which[number++]=jRow;
	    }
	  }
	}
      }
      sizeFactor_ += number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
    }
  }
  delete [] which;
  // Now we have size - create arrays and fill in
  choleskyRow_ = new int [sizeFactor_];
  sparseFactor_ = new double[sizeFactor_];
  sizeFactor_=0;
  which = choleskyRow_;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int number=0;
    choleskyStart_[iRow]=sizeFactor_;
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	CoinBigIndex start=columnStart[iColumn];
	CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	for (CoinBigIndex j=start;j<end;j++) {
	  int jRow=row[j];
	  if (jRow>=iRow&&!rowsDropped_[jRow]) {
	    if (!used[jRow]) {
	      used[jRow]=1;
	      which[number++]=jRow;
	    }
	  }
	}
      }
      sizeFactor_ += number;
      int j;
      for (j=0;j<number;j++)
	used[which[j]]=0;
      // move which on
      which += number;
    }
  }
  choleskyStart_[numberRows_]=sizeFactor_;
  delete [] used;
  permuteIn_ = new int [numberRows_];
  permuteOut_ = new int[numberRows_];
  integerParameters_[0]=0;
  int i0=0;
  int i1=1;
  wssmp(&numberRows_,choleskyStart_,choleskyRow_,sparseFactor_,
         NULL,permuteOut_,permuteIn_,0,&numberRows_,&i1,
         NULL,&i0,NULL,integerParameters_,doubleParameters_);
  integerParameters_[1]=1;//order and symbolic
  integerParameters_[2]=2;
  integerParameters_[3]=0;//CSR
  integerParameters_[4]=0;//C style
  integerParameters_[13]=1;//reuse initial factorization space
  integerParameters_[15+0]=1;//ordering
  integerParameters_[15+1]=0;
  integerParameters_[15+2]=1;
  integerParameters_[15+3]=0;
  integerParameters_[15+4]=1;
  doubleParameters_[10]=1.0e-20;
  doubleParameters_[11]=1.0e-15;
  wssmp(&numberRows_,choleskyStart_,choleskyRow_,sparseFactor_,
         NULL,permuteOut_,permuteIn_,NULL,&numberRows_,&i1,
         NULL,&i0,NULL,integerParameters_,doubleParameters_);
  std::cout<<"Ordering and symbolic factorization took "<<doubleParameters_[0]<<std::endl;
  std::cout<<integerParameters_[23]<<" elements in sparse Cholesky"<<std::endl;
  return 0;
}
/* Factorize - filling in rowsDropped and returning number dropped */
int 
ClpCholeskyWssmp::factorize(const double * diagonal, int * rowsDropped) 
{
  const CoinBigIndex * columnStart = model_->matrix()->getVectorStarts();
  const int * columnLength = model_->matrix()->getVectorLengths();
  const int * row = model_->matrix()->getIndices();
  const double * element = model_->matrix()->getElements();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths();
  const int * column = rowCopy_->getIndices();
  const double * elementByRow = rowCopy_->getElements();
  int numberColumns=model_->matrix()->getNumCols();
  int iRow;
  double * work = new double[numberRows_];
  CoinZeroN(work,numberRows_);
  const double * diagonalSlack = diagonal + numberColumns;
  int newDropped=0;
  double largest;
  double smallest;
  //perturbation
  double perturbation=model_->diagonalPerturbation()*model_->diagonalNorm();
  perturbation=perturbation*perturbation;
  if (perturbation>1.0) {
    //if (model_->model()->logLevel()&4) 
      std::cout <<"large perturbation "<<perturbation<<std::endl;
    perturbation=sqrt(perturbation);;
    perturbation=1.0;
  } 
  for (iRow=0;iRow<numberRows_;iRow++) {
    double * put = sparseFactor_+choleskyStart_[iRow];
    int * which = choleskyRow_+choleskyStart_[iRow];
    int number = choleskyStart_[iRow+1]-choleskyStart_[iRow];
    if (!rowsDropped_[iRow]) {
      CoinBigIndex startRow=rowStart[iRow];
      CoinBigIndex endRow=rowStart[iRow]+rowLength[iRow];
      work[iRow] = diagonalSlack[iRow];
      for (CoinBigIndex k=startRow;k<endRow;k++) {
	int iColumn=column[k];
	CoinBigIndex start=columnStart[iColumn];
	CoinBigIndex end=columnStart[iColumn]+columnLength[iColumn];
	double multiplier = diagonal[iColumn]*elementByRow[k];
	for (CoinBigIndex j=start;j<end;j++) {
	  int jRow=row[j];
	  if (jRow>=iRow&&!rowsDropped_[jRow]) {
	    double value=element[j]*multiplier;
	    work[jRow] += value;
	  }
	}
      }
      int j;
      for (j=0;j<number;j++) {
	int jRow = which[j];
	put[j]=work[jRow];
	work[jRow]=0.0;
      }
    }
  }
  //check sizes
  double largest2=maximumAbsElement(sparseFactor_,sizeFactor_);
  largest2*=1.0e-19;
  largest = min (largest2,1.0e-11);
  int numberDroppedBefore=0;
  for (int iRow=0;iRow<numberRows_;iRow++) {
    int dropped=rowsDropped_[iRow];
    // Move to int array
    rowsDropped[iRow]=dropped;
    if (!dropped) {
      CoinBigIndex start = choleskyStart_[iRow];
      double diagonal = sparseFactor_[start];
      if (diagonal>largest2) {
	sparseFactor_[start]=diagonal+perturbation;
      } else {
	sparseFactor_[start]=diagonal+perturbation;
	rowsDropped[iRow]=2;
	numberDroppedBefore++;
      } 
    } 
  } 
  int i1=1;
  int i0=0;
  integerParameters_[1]=3;
  integerParameters_[2]=3;
  integerParameters_[10]=2;
  //integerParameters_[11]=1;
  integerParameters_[12]=2;
  wssmp(&numberRows_,choleskyStart_,choleskyRow_,sparseFactor_,
	NULL,permuteOut_,permuteIn_,NULL,&numberRows_,&i1,
	NULL,&i0,rowsDropped,integerParameters_,doubleParameters_);
  //    NULL,&i0,(int *) diagonal,integerParameters_,doubleParameters_);
  std::cout<<"factorization took "<<doubleParameters_[0]<<std::endl;
  if (integerParameters_[9]) {
    std::cout<<"scaling applied"<<std::endl;
  } 
  newDropped=integerParameters_[20]+numberDroppedBefore;
  if (integerParameters_[20]) 
    std::cout<<integerParameters_[20]<<" rows dropped"<<std::endl;
  largest=doubleParameters_[3];
  smallest=doubleParameters_[4];
  delete [] work;
  //if (model_->model()->logLevel()&1) 
    std::cout<<"Cholesky - largest "<<largest<<" smallest "<<smallest<<std::endl;
  choleskyCondition_=largest/smallest;
  bool cleanCholesky;
  if (model_->numberIterations()<10000) 
    cleanCholesky=true;
  else 
    cleanCholesky=false;
  if (cleanCholesky) {
    //drop fresh makes some formADAT easier
    int oldDropped=numberRowsDropped_;
    if (newDropped||numberRowsDropped_) {
      std::cout <<"Rank "<<numberRows_-newDropped<<" ( "<<
          newDropped<<" dropped)";
      if (newDropped>oldDropped) 
        std::cout<<" ( "<<newDropped-oldDropped<<" dropped this time)";
      std::cout<<std::endl;
      newDropped=0;
      for (int i=0;i<numberRows_;i++) {
	char dropped = rowsDropped[i];
	rowsDropped_[i]=dropped;
        if (dropped==2) {
          //dropped this time
          rowsDropped[newDropped++]=i;
          rowsDropped_[i]=0;
        } 
      } 
      numberRowsDropped_=newDropped;
      newDropped=-(1+newDropped);
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
    if (numberRowsDropped_) {
      std::cout <<"Rank "<<numberRows_-numberRowsDropped_<<" ( "<<
          numberRowsDropped_<<" dropped)";
      if (newDropped) {
        std::cout<<" ( "<<newDropped<<" dropped this time)";
      } 
      std::cout<<std::endl;
    } 
  } 
  status_=0;
  return newDropped;
}
/* Uses factorization to solve. */
void 
ClpCholeskyWssmp::solve (double * region) 
{
  int i1=1;
  int i0=0;
  integerParameters_[1]=4;
  integerParameters_[2]=4;
#if 0
  integerParameters_[5]=3;
  doubleParameters_[5]=1.0e-10;
  integerParameters_[6]=6;
#endif
  wssmp(&numberRows_,choleskyStart_,choleskyRow_,sparseFactor_,
       NULL,permuteOut_,permuteIn_,region,&numberRows_,&i1,
       NULL,&i0,NULL,integerParameters_,doubleParameters_);
#if 0
  if (integerParameters_[5]) {
    std::cout<<integerParameters_[5]<<" refinements ";
  } 
  std::cout<<doubleParameters_[6]<<std::endl;
#endif
}
