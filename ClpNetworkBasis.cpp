// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpNetworkBasis.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"
#include "CoinIndexedVector.hpp"


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNetworkBasis::ClpNetworkBasis () 
{
  slackValue_=-1.0;
  numberRows_=0;
  numberColumns_=0;
  root_ = -1;
  leaf_ = -1;
  parent_ = NULL;
  descendant_ = NULL;
  pivot_ = NULL;
  rightSibling_ = NULL;
  leftSibling_ = NULL;
  sign_ = NULL;
  stack_ = NULL;
  toLeaf_ = NULL;
  toRoot_ = NULL;
  mark_ = NULL;
  model_=NULL;
}
// Constructor from CoinFactorization
ClpNetworkBasis::ClpNetworkBasis(const ClpSimplex * model,
				 int numberRows, const double * pivotRegion,
				 const int * permuteBack,
				 const int * startColumn, 
				 const int * numberInColumn,
				 const int * indexRow, const double * element)
{
  slackValue_=-1.0;
  numberRows_=numberRows;
  numberColumns_=numberRows;
  parent_ = new int [ numberRows_+1];
  descendant_ = new int [ numberRows_+1];
  pivot_ = new int [ numberRows_+1];
  rightSibling_ = new int [ numberRows_+1];
  leftSibling_ = new int [ numberRows_+1];
  sign_ = new double [ numberRows_+1];
  stack_ = new int [ numberRows_+1];
  toLeaf_ = new int [numberRows_+1];
  toRoot_ = new int [numberRows_+1];
  mark_ = new char[numberRows_+1];
  int i;
  for (i=0;i<numberRows_+1;i++) {
    parent_[i]=-1;
    descendant_[i]=-1;
    pivot_[i]=-1;
    rightSibling_[i]=-1;
    leftSibling_[i]=-1;
    sign_[i]=-1.0;
    stack_[i]=-1;
    mark_[i]=0;
  }
  // pivotColumnBack gives order of pivoting into basis
  // so pivotColumnback[0] is first slack in basis and
  // it pivots on row permuteBack[0]
  // a known root is given by permuteBack[numberRows_-1]
  root_ = numberRows_;
  int lastPivot=numberRows_;
  for (i=0;i<numberRows_;i++) {
    int iPivot = permuteBack[i];
    toRoot_[iPivot] = lastPivot;
    toLeaf_[lastPivot]=iPivot;
    lastPivot=iPivot;
    double sign;
    if (pivotRegion[i]>0.0)
      sign = 1.0;
    else
      sign =-1.0;
    int other;
    if (numberInColumn[i]>0) {
      int iRow = indexRow[startColumn[i]];
      other = permuteBack[iRow];
      assert (parent_[other]!=-1);
    } else {
      other = numberRows_;
    }
    sign_[iPivot] = sign;
    int iParent = other;
    parent_[iPivot] = other;
    if (descendant_[iParent]>=0) {
      // we have a sibling
      int iRight = descendant_[iParent];
      rightSibling_[iPivot]=iRight;
      leftSibling_[iRight]=iPivot;
    } else {
      rightSibling_[iPivot]=-1;
    }	    
    descendant_[iParent] = iPivot;
    leftSibling_[iPivot]=-1;
  }
  toLeaf_[lastPivot]=numberRows_;
  leaf_ = lastPivot;
  model_=model;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNetworkBasis::ClpNetworkBasis (const ClpNetworkBasis & rhs) 
{
  slackValue_=rhs.slackValue_;
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  root_ = rhs.root_;
  leaf_ = rhs.leaf_;
  if (rhs.parent_) {
    parent_ = new int [numberRows_+1];
    memcpy(parent_,rhs.parent_,(numberRows_+1)*sizeof(int));
  } else {
    parent_ = NULL;
  }
  if (rhs.descendant_) {
    descendant_ = new int [numberRows_+1];
    memcpy(descendant_,rhs.descendant_,(numberRows_+1)*sizeof(int));
  } else {
    descendant_ = NULL;
  }
  if (rhs.pivot_) {
    pivot_ = new int [numberRows_+1];
    memcpy(pivot_,rhs.pivot_,(numberRows_+1)*sizeof(int));
  } else {
    pivot_ = NULL;
  }
  if (rhs.rightSibling_) {
    rightSibling_ = new int [numberRows_+1];
    memcpy(rightSibling_,rhs.rightSibling_,(numberRows_+1)*sizeof(int));
  } else {
    rightSibling_ = NULL;
  }
  if (rhs.leftSibling_) {
    leftSibling_ = new int [numberRows_+1];
    memcpy(leftSibling_,rhs.leftSibling_,(numberRows_+1)*sizeof(int));
  } else {
    leftSibling_ = NULL;
  }
  if (rhs.sign_) {
    sign_ = new double [numberRows_+1];
    memcpy(sign_,rhs.sign_,(numberRows_+1)*sizeof(double));
  } else {
    sign_ = NULL;
  }
  if (rhs.stack_) {
    stack_ = new int [numberRows_+1];
    memcpy(stack_,rhs.stack_,(numberRows_+1)*sizeof(int));
  } else {
    stack_ = NULL;
  }
  if (rhs.toLeaf_) {
    toLeaf_ = new int [numberRows_+1];
    memcpy(toLeaf_,rhs.toLeaf_,(numberRows_+1)*sizeof(int));
  } else {
    toLeaf_ = NULL;
  }
  if (rhs.toRoot_) {
    toRoot_ = new int [numberRows_+1];
    memcpy(toRoot_,rhs.toRoot_,(numberRows_+1)*sizeof(int));
  } else {
    toRoot_ = NULL;
  }
  if (rhs.mark_) {
    mark_ = new char [numberRows_+1];
    memcpy(mark_,rhs.mark_,(numberRows_+1)*sizeof(char));
  } else {
    mark_ = NULL;
  }
  model_=rhs.model_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNetworkBasis::~ClpNetworkBasis () 
{
  delete [] parent_;
  delete [] descendant_;
  delete [] pivot_;
  delete [] rightSibling_;
  delete [] leftSibling_;
  delete [] sign_;
  delete [] stack_;
  delete [] toLeaf_;
  delete [] toRoot_;
  delete [] mark_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNetworkBasis &
ClpNetworkBasis::operator=(const ClpNetworkBasis& rhs)
{
  if (this != &rhs) {
    delete [] parent_;
    delete [] descendant_;
    delete [] pivot_;
    delete [] rightSibling_;
    delete [] leftSibling_;
    delete [] sign_;
    delete [] stack_;
    delete [] toLeaf_;
    delete [] toRoot_;
    delete [] mark_;
    slackValue_=rhs.slackValue_;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    root_ = rhs.root_;
    leaf_ = rhs.leaf_;
    if (rhs.parent_) {
      parent_ = new int [numberRows_+1];
      memcpy(parent_,rhs.parent_,(numberRows_+1)*sizeof(int));
    } else {
      parent_ = NULL;
    }
    if (rhs.descendant_) {
      descendant_ = new int [numberRows_+1];
      memcpy(descendant_,rhs.descendant_,(numberRows_+1)*sizeof(int));
    } else {
      descendant_ = NULL;
    }
    if (rhs.pivot_) {
      pivot_ = new int [numberRows_+1];
      memcpy(pivot_,rhs.pivot_,(numberRows_+1)*sizeof(int));
    } else {
      pivot_ = NULL;
    }
    if (rhs.rightSibling_) {
      rightSibling_ = new int [numberRows_+1];
      memcpy(rightSibling_,rhs.rightSibling_,(numberRows_+1)*sizeof(int));
    } else {
      rightSibling_ = NULL;
    }
    if (rhs.leftSibling_) {
      leftSibling_ = new int [numberRows_+1];
      memcpy(leftSibling_,rhs.leftSibling_,(numberRows_+1)*sizeof(int));
    } else {
      leftSibling_ = NULL;
    }
    if (rhs.sign_) {
      sign_ = new double [numberRows_+1];
      memcpy(sign_,rhs.sign_,(numberRows_+1)*sizeof(double));
    } else {
      sign_ = NULL;
    }
    if (rhs.stack_) {
      stack_ = new int [numberRows_+1];
      memcpy(stack_,rhs.stack_,(numberRows_+1)*sizeof(int));
    } else {
      stack_ = NULL;
    }
    if (rhs.toLeaf_) {
      toLeaf_ = new int [numberRows_+1];
      memcpy(toLeaf_,rhs.toLeaf_,(numberRows_+1)*sizeof(int));
    } else {
      toLeaf_ = NULL;
    }
    if (rhs.toRoot_) {
      toRoot_ = new int [numberRows_+1];
      memcpy(toRoot_,rhs.toRoot_,(numberRows_+1)*sizeof(int));
    } else {
      toRoot_ = NULL;
    }
    if (rhs.mark_) {
      mark_ = new char [numberRows_+1];
      memcpy(mark_,rhs.mark_,(numberRows_+1)*sizeof(char));
    } else {
    mark_ = NULL;
    }
    model_=rhs.model_;
  }
  return *this;
}
/* Replaces one Column to basis,
   returns 0=OK
*/
int 
ClpNetworkBasis::replaceColumn ( CoinIndexedVector * regionSparse,
				 int pivotRow)
{
  // regionSparse is empty
  assert (!regionSparse->getNumElements());
  model_->unpack(regionSparse, model_->sequenceIn());
  // arc given by pivotRow is leaving basis
  int iParent = parent_[pivotRow];
  // arc coming in has these two nodes
  int * indices = regionSparse->getIndices();
  int iRow0 = indices[0];
  int iRow1;
  if (regionSparse->getNumElements()==2)
    iRow1 = indices[1];
  else
    iRow1 = numberRows_;
  double sign = -regionSparse->denseVector()[iRow0];
  printf("In %d (%g) %d pivoting on %d\n",
	 iRow1, sign, iRow0,pivotRow);
  regionSparse->clear();
  // take out of tree
  int iLeft = leftSibling_[pivotRow];
  int iRight = rightSibling_[pivotRow];
  if (iLeft>=0) {
    rightSibling_[iLeft] = iRight;
    if (iRight>=0) 
      leftSibling_[iRight]=iLeft;
  } else if (iRight>=0) {
    leftSibling_[iRight]=-1;
    descendant_[iParent]=iRight;
  } else {
    descendant_[iParent]=-1;;
  }
  // move to other end of chain
  int descendant = descendant_[pivotRow];
  if (descendant>=0) {
    // make this descendant of that
    if (descendant_[descendant]>=0) {
      // we have a sibling
      int iRight = descendant_[descendant];
      rightSibling_[pivotRow]=iRight;
      leftSibling_[iRight]=pivotRow;
    } else {
      rightSibling_[pivotRow]=-1;
    }	    
    descendant_[descendant] = pivotRow;
    leftSibling_[pivotRow]=-1;
  }
  // now insert new one
  descendant = descendant_[iRow1];
  parent_[iRow0] = iRow1;
  sign_[iRow1]= sign;
  if (descendant>=0) {
    // we have a sibling
    int iRight = descendant;
    rightSibling_[iRow1]=iRight;
    leftSibling_[iRight]=iRow1;
  } else {
    rightSibling_[iRow1]=-1;
  }	    
  descendant_[descendant] = iRow1;
  leftSibling_[iRow1]=-1;
  return 0;
}

/* Updates one column (FTRAN) from region2 */
int 
ClpNetworkBasis::updateColumn ( CoinIndexedVector * regionSparse2)
{
  int iPivot=leaf_;
  int numberNonZero=0;
  double * array = regionSparse2->denseVector();
  while (iPivot!=numberRows_) {
    double pivotValue = array[iPivot];
    if (pivotValue) {
      numberNonZero++;
      int otherRow = parent_[iPivot];
      if (sign_[iPivot]<0) {
	array[iPivot] = - pivotValue;
	array[otherRow] += pivotValue;
      } else {
	array[otherRow] += pivotValue;
      }
    }
    iPivot = toRoot_[iPivot];
  }
  array[numberRows_]=0.0;
  return numberNonZero;
}
/* Updates one column (FTRAN) to/from array 
    ** For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. rhs) */
int 
ClpNetworkBasis::updateColumn ( double array[] ) const
{
  int iPivot=leaf_;
  int numberNonZero=0;
  while (iPivot!=numberRows_) {
    double pivotValue = array[iPivot];
    if (pivotValue) {
      numberNonZero++;
      int otherRow = parent_[iPivot];
      if (sign_[iPivot]<0) {
	array[iPivot] = - pivotValue;
	array[otherRow] += pivotValue;
      } else {
	array[otherRow] += pivotValue;
      }
    }
    iPivot = toRoot_[iPivot];
  }
  array[numberRows_]=0.0;
  return numberNonZero;
}
/* Updates one column transpose (BTRAN)
   For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. dense objective)
   returns number of nonzeros */
int 
ClpNetworkBasis::updateColumnTranspose ( double array[] ) const
{
  abort();
  return 1;
}
/* Updates one column (BTRAN) from region2 */
int 
ClpNetworkBasis::updateColumnTranspose ( 
					CoinIndexedVector * regionSparse2) const
{
  abort();
  return 1;
}
