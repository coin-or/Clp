// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpSimplexOther_H
#define ClpSimplexOther_H

#include "ClpSimplex.hpp"

/** This is for Simplex stuff which is neither dual nor primal

    It inherits from ClpSimplex.  It has no data of its own and 
    is never created - only cast from a ClpSimplex object at algorithm time. 

*/

class ClpSimplexOther : public ClpSimplex {

public:

  /**@name Methods */
  //@{
  /** Dual ranging.
      This computes increase/decrease in cost for each given variable and corresponding
      sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
      and numberColumns.. for artificials/slacks.
      For non-basic variables the sequence number will be that of the non-basic variables.

      Up to user to provide correct length arrays.

      When here - guaranteed optimal
  */
  void dualRanging(int numberCheck,const int * which,
		  double * costIncrease, int * sequenceIncrease,
		  double * costDecrease, int * sequenceDecrease);
  /** 
      Row array has row part of pivot row
      Column array has column part.
      This is used in dual ranging
  */
  void checkRatios(CoinIndexedVector * rowArray,
		   CoinIndexedVector * columnArray,
		   double & costIncrease, int & sequenceIncrease,
		   double & costDecrease, int & sequencedecrease);
  //@}
};
#endif
