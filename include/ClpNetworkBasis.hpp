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
/** This deals with Factorization and Updates for network structures
 */


class ClpNetworkBasis {

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    ClpNetworkBasis (  );
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
  /** Updates one column (FTRAN) from region2
      number returned is negative if no room
      region1 starts as zero and is zero at end */
  int updateColumn ( CoinIndexedVector * regionSparse,
		     CoinIndexedVector * regionSparse2);
  /** Updates one column (FTRAN) to/from array 
      number returned is negative if no room
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. rhs)
      region starts as zero and is zero at end */
  int updateColumn ( CoinIndexedVector * regionSparse,
			double array[] ) const;
  /** Updates one column transpose (BTRAN)
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. dense objective)
      returns number of nonzeros */
  int updateColumnTranspose ( CoinIndexedVector * regionSparse,
				 double array[] ) const;
  /** Updates one column (BTRAN) from region2
      region1 starts as zero and is zero at end */
  int updateColumnTranspose ( CoinIndexedVector * regionSparse,
			      CoinIndexedVector * regionSparse2) const;
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
  /// Maximum number of pivots before factorization
  int maximumPivots_;
  /// Number pivots since last factorization
  int numberPivots_;
  /// Pivot order for each Column
  int *pivotColumn_;
  /// Permutation vector for pivot row order
  int *permute_;
  /// DePermutation vector for pivot row order
  int *permuteBack_;
  //@}
};
#endif
