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
  /** Primal ranging.
      This computes increase/decrease in value for each given variable and corresponding
      sequence numbers which would change basis.  Sequence numbers are 0..numberColumns 
      and numberColumns.. for artificials/slacks.
      For basic variables the sequence number will be that of the basic variables.

      Up to user to providen correct length arrays.

      When here - guaranteed optimal
  */
  void primalRanging(int numberCheck,const int * which,
		  double * valueIncrease, int * sequenceIncrease,
		  double * valueDecrease, int * sequenceDecrease);
  /** 
      Row array has row part of pivot row
      Column array has column part.
      This is used in dual ranging
  */
  void checkDualRatios(CoinIndexedVector * rowArray,
		   CoinIndexedVector * columnArray,
		   double & costIncrease, int & sequenceIncrease,
		   double & costDecrease, int & sequencedecrease);
  /** 
      Row array has pivot column
      This is used in primal ranging
  */
  void checkPrimalRatios(CoinIndexedVector * rowArray,
			 int direction);
    /** Write the basis in MPS format to the specified file.
	If writeValues true writes values of structurals
	(and adds VALUES to end of NAME card)

	Row and column names may be null.
	formatType is
	<ul>
	  <li> 0 - normal
	  <li> 1 - extra accuracy 
	  <li> 2 - IEEE hex (later)
	</ul>

	Returns non-zero on I/O error
    */
    int writeBasis(const char *filename,
		 bool writeValues=false,
		 int formatType=0) const;
    /// Read a basis from the given filename
    int readBasis(const char *filename);
  //@}
};
#endif
