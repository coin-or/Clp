// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpFactorization_H
#define ClpFactorization_H


#include "CoinPragma.hpp"

#include "CoinFactorization.hpp"

/** This just implements CoinFactorization when an ClpMatrixBase object
    is passed.  If a network then has a dummy CoinFactorization and
    a genuine ClpNetworkBasis object
*/
class ClpMatrixBase;
class ClpSimplex;
class ClpNetworkBasis;

class ClpFactorization : public CoinFactorization {
  
public:
  /**@name factorization */
   //@{
  /** When part of LP - given by basic variables.
  Actually does factorization.
  Arrays passed in have non negative value to say basic.
  If status is okay, basic variables have pivot row - this is only needed
  if increasingRows_ >1.
  Allows scaling
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  int factorize (ClpSimplex * model,int solveType, bool valuesPass);
   //@}


  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpFactorization();
   /** Destructor */
   ~ClpFactorization();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpFactorization(const ClpFactorization&);
   /** The copy constructor from an CoinFactorization. */
   ClpFactorization(const CoinFactorization&);

   ClpFactorization& operator=(const ClpFactorization&);
   //@}
   
  /*  **** below here is so can use networkish basis */
  /**@name rank one updates which do exist */
  //@{

  /** Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
      If checkBeforeModifying is true will do all accuracy checks
      before modifying factorization.  Whether to set this depends on
      speed considerations.  You could just do this on first iteration
      after factorization and thereafter re-factorize
   partial update already in U */
  int replaceColumn ( CoinIndexedVector * regionSparse,
		      int pivotRow,
		      double pivotCheck ,
		      bool checkBeforeModifying=false);
  //@}

  /**@name various uses of factorization (return code number elements) 
   which user may want to know about */
  //@{
  /** Updates one column (FTRAN) from region2
      Tries to do FT update
      number returned is negative if no room
      region1 starts as zero and is zero at end */
  int updateColumnFT ( CoinIndexedVector * regionSparse,
		       CoinIndexedVector * regionSparse2);
  int updateColumn ( CoinIndexedVector * regionSparse,
		     CoinIndexedVector * regionSparse2,
		     bool noPermute=false) const;
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
    
  /**@name other stuff */
  //@{ 
  /** makes a row copy of L for speed and to allow very sparse problems */
  void goSparse();
  /// Cleans up i.e. gets rid of network basis 
  void cleanUp();
  /// Says whether to redo pivot order
  bool needToReorder() const;
  /// Says if a network basis
  bool inline networkBasis() const
  { return (networkBasis_!=NULL);};
  //@}

////////////////// data //////////////////
private:

  /**@name data */
  //@{
  /// Pointer to network basis
  ClpNetworkBasis * networkBasis_;
  //@}
};

#endif
