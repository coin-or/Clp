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

    It inherits from ClpSimplex.  It has no data of its own and 
    is never created - only cast from a ClpSimplex object at algorithm time. 

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
  /// Beale's method
  int primalBeale();
  //@}

};
#endif

