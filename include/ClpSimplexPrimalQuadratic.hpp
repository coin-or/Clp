// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpSimplexPrimalQuadratic_H
#define ClpSimplexPrimalQuadratic_H

#include "ClpSimplexPrimal.hpp"

/** This solves LPs using the primal simplex method

    It inherits from ClpSimplexPrimal.  It has no data of its own and 
    is never created - only cast from a ClpSimplexPrimal object at algorithm time. 
    If needed create new class and pass around

*/

class ClpSimplexPrimalQuadratic : public ClpSimplexPrimal {

public:

  /**@name Description of algorithm */
  //@{
  /** Primal algorithms for quadratic
      At present we have two algorithms:

      a) Beale's algorithm - this is in for sentimental reasons
         not because I think it is best one.
      b) Using a semi-trust region approach as for pooling problem
         This is in because I have it lying around

  */
  /// A sequential LP method
  int primalSLP(int numberPasses, double deltaTolerance);
  /// Wolfe's method (actually a mixture with Jensen and King)
  int primalWolfe();
  /// This is done after first pass
  int primalWolfe2 (const ClpSimplexPrimalQuadratic * originalModel);
  /// Main part
  int whileIterating (const ClpSimplexPrimalQuadratic * originalModel,
		      int & sequenceIn,
		      int & crucialSj);
  /** 
      Row array has pivot column
      This chooses pivot row.
      Rhs array is used for distance to next bound (for speed)
      For speed, we may need to go to a bucket approach when many
      variables go through bounds
      On exit rhsArray will have changes in costs of basic variables
      Initially no go thru
      Returns 0 - can do normal iteration
      1 - losing complementarity
  */
  int primalRow(CoinIndexedVector * rowArray,
		CoinIndexedVector * rhsArray,
		CoinIndexedVector * spareArray,
		CoinIndexedVector * spareArray2,
		const ClpSimplexPrimalQuadratic * originalModel,
		int crucialSj);
  //@}

};
#endif

