// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpNetworkBasis_H
#define ClpNetworkBasis_H

class ClpMatrixBase;
class CoinIndexedVector;
class ClpSimplex;

/** This deals with Factorization and Updates for network structures
 */


class ClpNetworkBasis {

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    ClpNetworkBasis (  );
  /// Constructor from CoinFactorization
  ClpNetworkBasis(const ClpSimplex * model,
		  int numberRows, const double * pivotRegion,
		  const int * permuteBack,const int * startColumn,
		  const int * numberInColumn,
		  const int * indexRow, const double * element);
  /// Copy constructor 
  ClpNetworkBasis ( const ClpNetworkBasis &other);

  /// Destructor
   ~ClpNetworkBasis (  );
  /// = copy
    ClpNetworkBasis & operator = ( const ClpNetworkBasis & other );
  //@}

  /**@name Do factorization */
  //@{
  /** When part of LP - given by basic variables.
  Actually does factorization.
  Arrays passed in have non negative value to say basic.
  If status is okay, basic variables have pivot row - this is only needed
  if increasingRows_ >1.
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -2 too many in basis */
  int factorize ( const ClpMatrixBase * matrix, 
		  int rowIsBasic[], int columnIsBasic[]);
  //@}

  /**@name rank one updates which do exist */
  //@{

  /** Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular!!
  */
  int replaceColumn ( CoinIndexedVector * column,
		      int pivotRow);
  //@}

  /**@name various uses of factorization (return code number elements) 
   which user may want to know about */
  //@{
  /** Updates one column (FTRAN) from region */
  int updateColumn ( CoinIndexedVector * regionSparse2);
  /** Updates one column (FTRAN) to/from array 
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. rhs) */
  int updateColumn ( double array[] ) const;
  /** Updates one column transpose (BTRAN)
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. dense objective)
      returns number of nonzeros */
  int updateColumnTranspose ( double array[] ) const;
  /** Updates one column (BTRAN) from region2 */
  int updateColumnTranspose ( CoinIndexedVector * regionSparse2) const;
  //@}
////////////////// data //////////////////
private:

  /**@name data */
  //@{
  /// Whether slack value is  +1 or -1
  double slackValue_;
  /// Number of Rows in factorization
  int numberRows_;
  /// Number of Columns in factorization
  int numberColumns_;
  /// Index of root
  int root_;
  /// Index of extreme leaf
  int leaf_;
  /// model
  const ClpSimplex * model_; 
  /// Parent for each column
  int * parent_;
  /// Descendant
  int * descendant_;
  /// Pivot row
  int * pivot_;
  /// Right sibling
  int * rightSibling_;
  /// Left sibling
  int * leftSibling_;
  /// Sign of pivot
  double * sign_;
  /// Stack
  int * stack_;
  /// Next one towards leaf
  int * toLeaf_;
  /// Next one towards root
  int * toRoot_;
  /// To mark rows
  char * mark_;
  //@}
};
#endif
