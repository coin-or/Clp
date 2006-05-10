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
      For non-basic variables the information is trivial to compute and the change in cost is just minus the 
      reduced cost and the sequence number will be that of the non-basic variables.
      For basic variables a ratio test is between the reduced costs for non-basic variables
      and the row of the tableau corresponding to the basic variable.
      The increase/decrease value is always >= 0.0

      Up to user to provide correct length arrays where each array is of length numberCheck.
      which contains list of variables for which information is desired.  All other
      arrays will be filled in by function.  If fifth entry in which is variable 7 then fifth entry in output arrays
      will information for variable 7.

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

      Up to user to provide correct length arrays where each array is of length numberCheck.
      which contains list of variables for which information is desired.  All other
      arrays will be filled in by function.  If fifth entry in which is variable 7 then fifth entry in output arrays
      will information for variable 7.

      When here - guaranteed optimal
  */
  void primalRanging(int numberCheck,const int * which,
		  double * valueIncrease, int * sequenceIncrease,
		  double * valueDecrease, int * sequenceDecrease);
  /** Parametrics
      This is an initial slow version.
      The code uses current bounds + theta * change (if change array not NULL)
      and similarly for objective.
      It starts at startingTheta and returns ending theta in endingTheta.
      If reportIncrement 0.0 it will report on any movement
      If reportIncrement >0.0 it will report at startingTheta+k*reportIncrement.
      If it can not reach input endingTheta return code will be 1 for infeasible,
      2 for unbounded, if error on ranges -1,  otherwise 0.
      Normal report is just theta and objective but
      if event handler exists it may do more
      On exit endingTheta is maximum reached (can be used for next startingTheta)
  */
  int parametrics(double startingTheta, double & endingTheta,double reportIncrement,
                  const double * changeLowerBound, const double * changeUpperBound,
                  const double * changeLowerRhs, const double * changeUpperRhs,
                  const double * changeObjective);
  /** Parametrics - inner loop
      This first attempt is when reportIncrement non zero and may
      not report endingTheta correctly
      If it can not reach input endingTheta return code will be 1 for infeasible,
      2 for unbounded,  otherwise 0.
      Normal report is just theta and objective but
      if event handler exists it may do more
  */
  int parametricsLoop(double startingTheta, double & endingTheta,double reportIncrement,
                      const double * changeLower, const double * changeUpper,
                      const double * changeObjective, ClpDataSave & data,
                      bool canTryQuick);
  /**  Refactorizes if necessary 
       Checks if finished.  Updates status.

       type - 0 initial so set up save arrays etc
            - 1 normal -if good update save
	    - 2 restoring from saved 
  */
  void statusOfProblemInParametrics(int type,ClpDataSave & saveData);
  /** This has the flow between re-factorizations

      Reasons to come out:
      -1 iterations etc
      -2 inaccuracy 
      -3 slight inaccuracy (and done iterations)
      +0 looks optimal (might be unbounded - but we will investigate)
      +1 looks infeasible
      +3 max iterations 
   */
  int whileIterating(double startingTheta, double & endingTheta,double reportIncrement,
                      const double * changeLower, const double * changeUpper,
                      const double * changeObjective);
  /** Computes next theta and says if objective or bounds (0= bounds, 1 objective, -1 none).
      theta is in theta_.
      type 1 bounds, 2 objective, 3 both.
  */
  int nextTheta(int type, double maxTheta, double * primalChange, double * dualChange,
                      const double * changeLower, const double * changeUpper,
                      const double * changeObjective);
  /** 
      Row array has row part of pivot row
      Column array has column part.
      This is used in dual ranging
  */
  void checkDualRatios(CoinIndexedVector * rowArray,
		   CoinIndexedVector * columnArray,
                       double & costIncrease, int & sequenceIncrease, double & alphaIncrease,
                       double & costDecrease, int & sequenceDecrease, double & alphaDecrease);
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
  /// Creates dual of a problem
  ClpSimplex * dualOfModel() const;
  /// Restores solution from dualized problem
  void restoreFromDual(const ClpSimplex * dualProblem);
  /** Does very cursory presolve.
      rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns.
  */
  ClpSimplex * crunch(double * rhs, int * whichRows, int * whichColumns,
                      int & nBound, bool moreBounds=false, bool tightenBounds=false);
  /** After very cursory presolve.
      rhs is numberRows, whichRows is 3*numberRows and whichColumns is 2*numberColumns.
  */
  void afterCrunch(const ClpSimplex & small,
                   const int * whichRows, const int * whichColumns,
                   int nBound);
  /** Tightens integer bounds - returns number tightened or -1 if infeasible
  */
  int tightenIntegerBounds(double * rhsSpace); 
  //@}
};
#endif
