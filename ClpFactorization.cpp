// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpFactorization.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization () :
   CoinFactorization() {}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization (const ClpFactorization & rhs) :
   CoinFactorization(rhs) {}

ClpFactorization::ClpFactorization (const CoinFactorization & rhs) :
   CoinFactorization(rhs) {}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization () {}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization& rhs)
{
  if (this != &rhs) {
    CoinFactorization::operator=(rhs);
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
