// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpNetworkBasis.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"


//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNetworkBasis::ClpNetworkBasis () 
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNetworkBasis::ClpNetworkBasis (const ClpNetworkBasis & rhs) 
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNetworkBasis::~ClpNetworkBasis () 
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNetworkBasis &
ClpNetworkBasis::operator=(const ClpNetworkBasis& rhs)
{
  if (this != &rhs) {
  }
  return *this;
}
/* Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular
*/
int 
ClpNetworkBasis::replaceColumn ( CoinIndexedVector * regionSparse,
				 int pivotRow)
{
  abort();
  return 1;
}

/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int 
ClpNetworkBasis::updateColumn ( CoinIndexedVector * regionSparse,
				CoinIndexedVector * regionSparse2)
{
  abort();
  return 1;
}
/* Updates one column (FTRAN) to/from array 
   number returned is negative if no room
   ** For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. rhs)
   region starts as zero and is zero at end */
int 
ClpNetworkBasis::updateColumn ( CoinIndexedVector * regionSparse,
			double array[] ) const
{
  abort();
  return 1;
}
/* Updates one column transpose (BTRAN)
   For large problems you should ALWAYS know where the nonzeros
   are, so please try and migrate to previous method after you
   have got code working using this simple method - thank you!
   (the only exception is if you know input is dense e.g. dense objective)
   returns number of nonzeros */
int 
ClpNetworkBasis::updateColumnTranspose ( CoinIndexedVector * regionSparse,
					  double array[] ) const
{
  abort();
  return 1;
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int 
ClpNetworkBasis::updateColumnTranspose ( CoinIndexedVector * regionSparse,
    			  CoinIndexedVector * regionSparse2) const
{
  abort();
  return 1;
}
