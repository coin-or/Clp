// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpFactorization_H
#define ClpFactorization_H


#include "CoinPragma.hpp"

#include "CoinFactorization.hpp"

/** This just implements CoinFactorization when an ClpMatrixBase object
    is passed.  It has no data.
*/
class ClpMatrixBase;
class ClpSimplex;

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
  int factorize (const ClpSimplex * model, 
		 const ClpMatrixBase * matrix, 
		  int numberRows, int numberColumns,
		  int rowIsBasic[], int columnIsBasic[] , 
		 double areaFactor = 0.0);
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
   
    
};

#endif
