// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpFactorization.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpNetworkBasis.hpp"


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
    return networkBasis_->replaceColumn(regionSparse,
					pivotRow);
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
  if (!networkBasis_) {
    return CoinFactorization::updateColumn(regionSparse,
					   regionSparse2,
					   FTUpdate);
  } else {
    return networkBasis_->updateColumn(regionSparse,
				       regionSparse2);
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
    return networkBasis_->updateColumn(regionSparse,
						    array);
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
    return networkBasis_->updateColumnTranspose(regionSparse,
						    array);
  }
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnTranspose ( CoinIndexedVector * regionSparse,
    			  CoinIndexedVector * regionSparse2) const
{
  if (!networkBasis_) {
    return CoinFactorization::updateColumnTranspose(regionSparse,
						    regionSparse2);
  } else {
    return networkBasis_->updateColumnTranspose(regionSparse,
						    regionSparse2);
  }
}
