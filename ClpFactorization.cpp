// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpFactorization.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpNetworkBasis.hpp"
#include "ClpNetworkMatrix.hpp"
//#define CHECK_NETWORK
#ifdef CHECK_NETWORK
const static bool doCheck=true;
#else
const static bool doCheck=false;
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization () :
   CoinFactorization() 
{
  networkBasis_ = NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization (const ClpFactorization & rhs) :
   CoinFactorization(rhs) 
{
  if (rhs.networkBasis_)
    networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
  else
    networkBasis_=NULL;
}

ClpFactorization::ClpFactorization (const CoinFactorization & rhs) :
   CoinFactorization(rhs) 
{
  networkBasis_=NULL;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization () 
{
  delete networkBasis_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization& rhs)
{
  if (this != &rhs) {
    CoinFactorization::operator=(rhs);
    delete networkBasis_;
    if (rhs.networkBasis_)
      networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
    else
      networkBasis_=NULL;
  }
  return *this;
}
int 
ClpFactorization::factorize ( const ClpSimplex * model,
			      const ClpMatrixBase * matrix, 
			      int numberRows, int numberColumns,
			      int rowIsBasic[], int columnIsBasic[] , 
			      double areaFactor )
{
  if (!networkBasis_||doCheck) {
    // see if matrix a network
    ClpNetworkMatrix* networkMatrix =
      dynamic_cast< ClpNetworkMatrix*>(model->clpMatrix());
    // If network - still allow ordinary factorization first time for laziness
    int saveMaximumPivots = maximumPivots();
    delete networkBasis_;
    networkBasis_ = NULL;
    if (networkMatrix&&!doCheck)
      maximumPivots(1);
    // maybe for speed will be better to leave as many regions as possible
    gutsOfDestructor();
    gutsOfInitialize(2);
    if (areaFactor)
      areaFactor_ = areaFactor;
    int numberBasic = 0;
    CoinBigIndex numberElements=0;
    int numberRowBasic=0;
    
    // compute how much in basis
    
    int i;
    
    for (i=0;i<numberRows;i++) {
      if (rowIsBasic[i]>=0)
	numberRowBasic++;
    }
    
    numberBasic = numberRowBasic;
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	numberBasic++;
      }
    }
    numberElements += matrix->numberInBasis(columnIsBasic);
    if ( numberBasic > numberRows ) {
    return -2; // say too many in basis
    }
    numberElements = 3 * numberBasic + 3 * numberElements + 10000;
    getAreas ( numberRows, numberBasic, numberElements,
	       2 * numberElements );
    //fill
    //copy
    numberBasic=0;
    numberElements=0;
    for (i=0;i<numberRows;i++) {
      if (rowIsBasic[i]>=0) {
	indexRowU_[numberElements]=i;
	indexColumnU_[numberElements]=numberBasic;
	elementU_[numberElements++]=slackValue_;
	numberBasic++;
      }
    }
    numberElements +=matrix->fillBasis(model, columnIsBasic, numberBasic, 
				       indexRowU_+numberElements, 
				       indexColumnU_+numberElements,
				       elementU_+numberElements);
    lengthU_ = numberElements;
    
    preProcess ( 0 );
    factor (  );
    numberBasic=0;
    if (status_ == 0) {
      int * permuteBack = permuteBack_;
      int * back = pivotColumnBack_;
      for (i=0;i<numberRows;i++) {
	if (rowIsBasic[i]>=0) {
	  rowIsBasic[i]=permuteBack[back[numberBasic++]];
	}
      }
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  columnIsBasic[i]=permuteBack[back[numberBasic++]];
	}
      }
      if (increasingRows_>1) {
	// Set up permutation vector
	if (increasingRows_<3) {
	  // these arrays start off as copies of permute
	  // (and we could use permute_ instead of pivotColumn (not back though))
	  ClpDisjointCopyN ( permute_, numberRows_ , pivotColumn_  );
	  ClpDisjointCopyN ( permuteBack_, numberRows_ , pivotColumnBack_  );
	}
      } else {
	// Set up permutation vector
	// (we could use permute_ instead of pivotColumn (not back though))
	for (i=0;i<numberRows_;i++) {
	  int k=pivotColumn_[i];
	  pivotColumn_[i]=pivotColumnBack_[i];
	  pivotColumnBack_[i]=k;
	}
      }
    } else if (status_ == -1) {
      // mark as basic or non basic
      for (i=0;i<numberRows_;i++) {
	if (rowIsBasic[i]>=0) {
	  if (pivotColumn_[numberBasic]>=0) 
	    rowIsBasic[i]=pivotColumn_[numberBasic];
	  else
	    rowIsBasic[i]=-1;
	  numberBasic++;
	}
      }
      for (i=0;i<numberColumns;i++) {
	if (columnIsBasic[i]>=0) {
	  if (pivotColumn_[numberBasic]>=0) 
	    columnIsBasic[i]=pivotColumn_[numberBasic];
	  else
	    columnIsBasic[i]=-1;
	  numberBasic++;
	}
      }
    }
    if (networkMatrix) {
      maximumPivots(saveMaximumPivots);
      if (!status_) {
	// create network factorization
	if (doCheck)
	  delete networkBasis_; // temp
        networkBasis_ = new ClpNetworkBasis(model,numberRows_,
					    pivotRegion_,
					    permuteBack_,
					    startColumnU_,
					    numberInColumn_,
					    indexRowU_,
					    elementU_);
	// kill off arrays in ordinary factorization
	if (!doCheck)
	  gutsOfDestructor();
      }
    }
    if (!status_&&!networkBasis_) {
      // See if worth going sparse and when
      if (numberFtranCounts_>100) {
	ftranAverageAfterL_ = max(ftranCountAfterL_/ftranCountInput_,1.0);
	ftranAverageAfterR_ = max(ftranCountAfterR_/ftranCountAfterL_,1.0);
	ftranAverageAfterU_ = max(ftranCountAfterU_/ftranCountAfterR_,1.0);
        if (btranCountInput_) {
	  btranAverageAfterU_ = max(btranCountAfterU_/btranCountInput_,1.0);
	  btranAverageAfterR_ = max(btranCountAfterR_/btranCountAfterU_,1.0);
	  btranAverageAfterL_ = max(btranCountAfterL_/btranCountAfterR_,1.0);
	} else {
	  // odd - we have not done any btrans
	  btranAverageAfterU_ = 1.0;
	  btranAverageAfterR_ = 1.0;
	  btranAverageAfterL_ = 1.0;
	}
      }
      // scale back

      ftranCountInput_ *= 0.8;
      ftranCountAfterL_ *= 0.8;
      ftranCountAfterR_ *= 0.8;
      ftranCountAfterU_ *= 0.8;
      btranCountInput_ *= 0.8;
      btranCountAfterU_ *= 0.8;
      btranCountAfterR_ *= 0.8;
      btranCountAfterL_ *= 0.8;
    }
  } else {
    // network - fake factorization - do nothing
    status_=0;
  }

  return status_;
}
/* Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   If checkBeforeModifying is true will do all accuracy checks
   before modifying factorization.  Whether to set this depends on
   speed considerations.  You could just do this on first iteration
   after factorization and thereafter re-factorize
   partial update already in U */
int 
ClpFactorization::replaceColumn ( CoinIndexedVector * regionSparse,
		      int pivotRow,
		      double pivotCheck ,
		      bool checkBeforeModifying)
{
  if (!networkBasis_) {
    return CoinFactorization::replaceColumn(regionSparse,
					    pivotRow,
					    pivotCheck,
					    checkBeforeModifying);
  } else {
    if (doCheck) {
      int returnCode = CoinFactorization::replaceColumn(regionSparse,
							pivotRow,
							pivotCheck,
							checkBeforeModifying);
      networkBasis_->replaceColumn(regionSparse,
				   pivotRow);
      return returnCode;
    } else {
      // increase number of pivots
      numberPivots_++;
      return networkBasis_->replaceColumn(regionSparse,
				   pivotRow);
    }
  }
}

/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumn ( CoinIndexedVector * regionSparse,
				 CoinIndexedVector * regionSparse2,
				 bool FTUpdate) 
{
#ifdef CLP_DEBUG
  regionSparse->checkClear();
#endif
  if (!networkBasis_) {
    collectStatistics_ = true;
    return CoinFactorization::updateColumn(regionSparse,
					   regionSparse2,
					   FTUpdate);
    collectStatistics_ = false;
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    int returnCode = CoinFactorization::updateColumn(regionSparse,
						     regionSparse2, FTUpdate);
    networkBasis_->updateColumn(regionSparse,save);
    int i;
    double * array = regionSparse2->denseVector();
    double * array2 = save->denseVector();
    for (i=0;i<numberRows_;i++) {
      double value1 = array[i];
      double value2 = array2[i];
      assert (value1==value2);
    }
    delete save;
    return returnCode;
#else
    return networkBasis_->updateColumn(regionSparse,regionSparse2);
#endif
  }
}
/* Updates one column (FTRAN) to/from array 
   number returned is negative if no room
   ** For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. rhs)
   region starts as zero and is zero at end */
int 
ClpFactorization::updateColumn ( CoinIndexedVector * regionSparse,
			double array[] ) const
{
  if (!networkBasis_) {
    return CoinFactorization::updateColumn(regionSparse,
					   array);
  } else {
#ifdef CHECK_NETWORK
    double * save = new double [numberRows_+1];
    memcpy(save,array,(numberRows_+1)*sizeof(double));
    int returnCode = CoinFactorization::updateColumn(regionSparse,
						     array);
    networkBasis_->updateColumn(regionSparse, save);
    int i;
    for (i=0;i<numberRows_;i++)
      assert (fabs(save[i]-array[i])<1.0e-8*(1.0+fabs(array[i])));
    delete [] save;
    return returnCode;
#else
    return networkBasis_->updateColumn(regionSparse, array);
#endif
  }
}
/* Updates one column transpose (BTRAN)
   For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. dense objective)
   returns number of nonzeros */
int 
ClpFactorization::updateColumnTranspose ( CoinIndexedVector * regionSparse,
					  double array[] ) const
{
  if (!networkBasis_) {
    return CoinFactorization::updateColumnTranspose(regionSparse,
						    array);
  } else {
#ifdef CHECK_NETWORK
    double * save = new double [numberRows_+1];
    memcpy(save,array,(numberRows_+1)*sizeof(double));
    int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
							      array);
    networkBasis_->updateColumnTranspose(regionSparse, save);
    int i;
    for (i=0;i<numberRows_;i++)
      assert (fabs(save[i]-array[i])<1.0e-8*(1.0+fabs(array[i])));
    delete [] save;
    return returnCode;
#else
    return networkBasis_->updateColumnTranspose(regionSparse, array);
#endif
  }
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnTranspose ( CoinIndexedVector * regionSparse,
    			  CoinIndexedVector * regionSparse2) const
{
  if (!networkBasis_) {
    collectStatistics_ = true;
    return CoinFactorization::updateColumnTranspose(regionSparse,
						    regionSparse2);
    collectStatistics_ = false;
  } else {
#ifdef CHECK_NETWORK
      CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
      int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
								regionSparse2);
      networkBasis_->updateColumnTranspose(regionSparse,save);
      int i;
      double * array = regionSparse2->denseVector();
      double * array2 = save->denseVector();
      for (i=0;i<numberRows_;i++) {
	double value1 = array[i];
	double value2 = array2[i];
	assert (value1==value2);
      }
      delete save;
      return returnCode;
#else
      return networkBasis_->updateColumnTranspose(regionSparse,regionSparse2);
#endif
  }
}
/* makes a row copy of L for speed and to allow very sparse problems */
void 
ClpFactorization::goSparse()
{
  if (!networkBasis_) 
    CoinFactorization::goSparse();
}
// Cleans up i.e. gets rid of network basis 
void 
ClpFactorization::cleanUp()
{
  delete networkBasis_;
  networkBasis_=NULL;
  resetStatistics();
}
/// Says whether to redo pivot order
bool 
ClpFactorization::needToReorder() const
{
#ifdef CHECK_NETWORK
  return true;
#endif
  if (!networkBasis_)
    return true;
  else
    return false;
}
