/* Copyright (C) 2002, 2003 International Business Machines
   Corporation and others.  All Rights Reserved.*/
#ifndef ClpSimplexC_H
#define ClpSimplexC_H

/* include all defines and ugly stuff */
#include "Clp_C_defines.h"

/** This is a first "C" interface to Clp.
    It has similarities to the OSL V3 interface
    and only has most common functions
*/

#ifdef __cplusplus
extern "C"{
#endif
  
  /**@name Constructors and destructor 
     These do not have an exact analogue in C++.
     The user does not need to know structure of Clp_Simplex.
     
     For all functions outside this group there is an exact C++
     analogue created by taking the first parameter out, removing the Clp_
     from name and applying the method to an object of type ClpSimplex.
  */
  /*@{*/

  /** Default constructor */
  CLPLIBAPI Clp_Simplex * CLPLINKAGE Clp_newModel(void);
  /** Destructor */
  CLPLIBAPI void CLPLINKAGE Clp_deleteModel(Clp_Simplex * model);
  /*@}*/

  /**@name Load model - loads some stuff and initializes others */
  /*@{*/
  /** Loads a problem (the constraints on the
      rows are given by lower and upper bounds). If a pointer is NULL then the
      following values are the default:
      <ul>
      <li> <code>colub</code>: all columns have upper bound infinity
      <li> <code>collb</code>: all columns have lower bound 0 
      <li> <code>rowub</code>: all rows have upper bound infinity
      <li> <code>rowlb</code>: all rows have lower bound -infinity
      <li> <code>obj</code>: all variables have 0 objective coefficient
      </ul>
  */
  /** Just like the other loadProblem() method except that the matrix is
	given in a standard column major ordered format (without gaps). */
  CLPLIBAPI void CLPLINKAGE Clp_loadProblem (Clp_Simplex * model,  const int numcols, const int numrows,
		      const CoinBigIndex * start, const int* index,
		      const double* value,
		      const double* collb, const double* colub,   
		      const double* obj,
		      const double* rowlb, const double* rowub);
  /** Read an mps file from the given filename */
  CLPLIBAPI int CLPLINKAGE Clp_readMps(Clp_Simplex * model,const char *filename,
	      int keepNames,
	      int ignoreErrors);
  /** Copy in integer informations */
  CLPLIBAPI void CLPLINKAGE Clp_copyInIntegerInformation(Clp_Simplex * model,const char * information);
  /** Drop integer informations */
  CLPLIBAPI void CLPLINKAGE Clp_deleteIntegerInformation(Clp_Simplex * model);
  /** Resizes rim part of model  */
  CLPLIBAPI void CLPLINKAGE Clp_resize (Clp_Simplex * model, int newNumberRows, int newNumberColumns);
  /** Deletes rows */
  CLPLIBAPI void CLPLINKAGE Clp_deleteRows(Clp_Simplex * model, int number, const int * which);
  /** Add rows */
  CLPLIBAPI void CLPLINKAGE Clp_addRows(Clp_Simplex * model, int number, const double * rowLower, 
	       const double * rowUpper,
	       const int * rowStarts, const int * columns,
	       const double * elements);

  /** Deletes columns */
  CLPLIBAPI void CLPLINKAGE Clp_deleteColumns(Clp_Simplex * model, int number, const int * which);
  /** Add columns */
  CLPLIBAPI void CLPLINKAGE Clp_addColumns(Clp_Simplex * model, int number, const double * columnLower, 
		  const double * columnUpper,
		  const double * objective,
		  const int * columnStarts, const int * rows,
		  const double * elements);
  /** Drops names - makes lengthnames 0 and names empty */
  CLPLIBAPI void CLPLINKAGE Clp_dropNames(Clp_Simplex * model);
  /** Copies in names */
  CLPLIBAPI void CLPLINKAGE Clp_copyNames(Clp_Simplex * model, const char * const * rowNames,
		 const char * const * columnNames);
  
  /*@}*/
  /**@name gets and sets - you will find some synonyms at the end of this file */
  /*@{*/ 
  /** Number of rows */
  CLPLIBAPI int CLPLINKAGE Clp_numberRows(Clp_Simplex * model);
  /** Number of columns */
  CLPLIBAPI int CLPLINKAGE Clp_numberColumns(Clp_Simplex * model);
  /** Primal tolerance to use */
  CLPLIBAPI double CLPLINKAGE Clp_primalTolerance(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setPrimalTolerance(Clp_Simplex * model,  double value) ;
  /** Dual tolerance to use */
  CLPLIBAPI double CLPLINKAGE Clp_dualTolerance(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setDualTolerance(Clp_Simplex * model,  double value) ;
  /** Dual objective limit */
  CLPLIBAPI double CLPLINKAGE Clp_dualObjectiveLimit(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setDualObjectiveLimit(Clp_Simplex * model, double value);
  /** Objective offset */
  CLPLIBAPI double CLPLINKAGE Clp_objectiveOffset(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setObjectiveOffset(Clp_Simplex * model, double value);
  /** Fills in array with problem name  */
  CLPLIBAPI void CLPLINKAGE Clp_problemName(Clp_Simplex * model, int maxNumberCharacters, char * array);
  /** Number of iterations */
  CLPLIBAPI int CLPLINKAGE Clp_numberIterations(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setNumberIterations(Clp_Simplex * model, int numberIterations);
  /** Maximum number of iterations */
  CLPLIBAPI int maximumIterations(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setMaximumIterations(Clp_Simplex * model, int value);
  /** Maximum time in seconds (from when set called) */
  CLPLIBAPI double CLPLINKAGE Clp_maximumSeconds(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setMaximumSeconds(Clp_Simplex * model, double value);
  /** Returns true if hit maximum iterations (or time) */
  CLPLIBAPI int CLPLINKAGE Clp_hitMaximumIterations(Clp_Simplex * model);
  /** Status of problem:
      0 - optimal
      1 - primal infeasible
      2 - dual infeasible
      3 - stopped on iterations etc
      4 - stopped due to errors
  */
  CLPLIBAPI int CLPLINKAGE Clp_status(Clp_Simplex * model);
  /** Set problem status */
  CLPLIBAPI void CLPLINKAGE Clp_setProblemStatus(Clp_Simplex * model, int problemStatus);
  /** Secondary status of problem - may get extended
      0 - none
      1 - primal infeasible because dual limit reached
      2 - scaled problem optimal - unscaled has primal infeasibilities
      3 - scaled problem optimal - unscaled has dual infeasibilities
      4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
  */
  CLPLIBAPI int CLPLINKAGE Clp_secondaryStatus(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setSecondaryStatus(Clp_Simplex * model, int status);
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  CLPLIBAPI double CLPLINKAGE Clp_optimizationDirection(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setOptimizationDirection(Clp_Simplex * model, double value);
  /** Primal row solution */
  CLPLIBAPI double * CLPLINKAGE Clp_primalRowSolution(Clp_Simplex * model);
  /** Primal column solution */
  CLPLIBAPI double * CLPLINKAGE Clp_primalColumnSolution(Clp_Simplex * model);
  /** Dual row solution */
  CLPLIBAPI double * CLPLINKAGE Clp_dualRowSolution(Clp_Simplex * model);
  /** Reduced costs */
  CLPLIBAPI double * CLPLINKAGE Clp_dualColumnSolution(Clp_Simplex * model);
  /** Row lower */
  CLPLIBAPI double* CLPLINKAGE Clp_rowLower(Clp_Simplex * model);
  /** Row upper  */
  CLPLIBAPI double* CLPLINKAGE Clp_rowUpper(Clp_Simplex * model);
  /** Objective */
  CLPLIBAPI double * CLPLINKAGE Clp_objective(Clp_Simplex * model);            
  /** Column Lower */
  CLPLIBAPI double * CLPLINKAGE Clp_columnLower(Clp_Simplex * model);
  /** Column Upper */
  CLPLIBAPI double * CLPLINKAGE Clp_columnUpper(Clp_Simplex * model);
  /** Number of elements in matrix */
  CLPLIBAPI int CLPLINKAGE Clp_getNumElements(Clp_Simplex * model); 
  /** Objective value */
  CLPLIBAPI double CLPLINKAGE Clp_objectiveValue(Clp_Simplex * model);
  /** Integer information */
  CLPLIBAPI char * CLPLINKAGE Clp_integerInformation(Clp_Simplex * model);
  /** Infeasibility/unbounded ray (NULL returned if none/wrong)
      Up to user to use delete [] on these arrays.  */
  CLPLIBAPI double * CLPLINKAGE Clp_infeasibilityRay(Clp_Simplex * model);
  CLPLIBAPI double * CLPLINKAGE Clp_unboundedRay(Clp_Simplex * model);
  /** See if status array exists (partly for OsiClp) */
  CLPLIBAPI int CLPLINKAGE Clp_statusExists(Clp_Simplex * model);
  /** Return address of status array (char[numberRows+numberColumns]) */
  CLPLIBAPI unsigned char *  CLPLINKAGE Clp_statusArray(Clp_Simplex * model);
  /** Copy in status vector */
  CLPLIBAPI void CLPLINKAGE Clp_copyinStatus(Clp_Simplex * model, const unsigned char * statusArray);
  /* Get variable basis info */
  CLPLIBAPI const int CLPLINKAGE Clp_getColumnStatus(Clp_Simplex * model,int sequence);
  /* Get row basis info */
  CLPLIBAPI const int CLPLINKAGE Clp_getRowStatus(Clp_Simplex * model,int sequence);
  
  /** User pointer for whatever reason */
  CLPLIBAPI void CLPLINKAGE Clp_setUserPointer (Clp_Simplex * model, void * pointer);
  CLPLIBAPI void * CLPLINKAGE Clp_getUserPointer (Clp_Simplex * model);
  /*@}*/
  /**@name Message handling.  Call backs are handled by ONE function */
  /*@{*/
  /** Pass in Callback function.
   Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
  CLPLIBAPI void CLPLINKAGE Clp_registerCallBack(Clp_Simplex * model, 
						   clp_callback userCallBack);
  /** Unset Callback function */
  CLPLIBAPI void CLPLINKAGE Clp_clearCallBack(Clp_Simplex * model);
  /** Amount of print out:
      0 - none
      1 - just final
      2 - just factorizations
      3 - as 2 plus a bit more
      4 - verbose
      above that 8,16,32 etc just for selective debug
  */
  CLPLIBAPI void CLPLINKAGE Clp_setLogLevel(Clp_Simplex * model, int value);
  CLPLIBAPI int CLPLINKAGE Clp_logLevel(Clp_Simplex * model);
  /** length of names (0 means no names0 */
  CLPLIBAPI int CLPLINKAGE Clp_lengthNames(Clp_Simplex * model);
  /** Fill in array (at least lengthNames+1 long) with a row name */
  CLPLIBAPI void CLPLINKAGE Clp_rowName(Clp_Simplex * model, int iRow, char * name);
  /** Fill in array (at least lengthNames+1 long) with a column name */
  CLPLIBAPI void CLPLINKAGE Clp_columnName(Clp_Simplex * model, int iColumn, char * name);

  /*@}*/


  /**@name Functions most useful to user */
  /*@{*/
  /** General solve algorithm which can do presolve.
      See  ClpSolve.hpp for options
   */
  CLPLIBAPI int CLPLINKAGE Clp_initialSolve(Clp_Simplex * model);
  /** Dual initial solve */
  CLPLIBAPI int CLPLINKAGE Clp_initialDualSolve(Clp_Simplex * model);
  /** Primal initial solve */
  CLPLIBAPI int CLPLINKAGE Clp_initialPrimalSolve(Clp_Simplex * model);
  /** Dual algorithm - see ClpSimplexDual.hpp for method */
  CLPLIBAPI int CLPLINKAGE Clp_dual(Clp_Simplex * model, int ifValuesPass);
  /** Primal algorithm - see ClpSimplexPrimal.hpp for method */
  CLPLIBAPI int CLPLINKAGE Clp_primal(Clp_Simplex * model, int ifValuesPass);
  /** Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
  CLPLIBAPI void CLPLINKAGE Clp_scaling(Clp_Simplex * model, int mode);
  /** Gets scalingFlag */
  CLPLIBAPI int CLPLINKAGE Clp_scalingFlag(Clp_Simplex * model);
  /** Crash - at present just aimed at dual, returns
      -2 if dual preferred and crash basis created
      -1 if dual preferred and all slack basis preferred
       0 if basis going in was not all slack
       1 if primal preferred and all slack basis preferred
       2 if primal preferred and crash basis created.
       
       if gap between bounds <="gap" variables can be flipped

       If "pivot" is
       0 No pivoting (so will just be choice of algorithm)
       1 Simple pivoting e.g. gub
       2 Mini iterations
  */
  CLPLIBAPI int CLPLINKAGE Clp_crash(Clp_Simplex * model, double gap,int pivot);
  /*@}*/


  /**@name most useful gets and sets */
  /*@{*/ 
  /** If problem is primal feasible */
  CLPLIBAPI int CLPLINKAGE Clp_primalFeasible(Clp_Simplex * model);
  /** If problem is dual feasible */
  CLPLIBAPI int CLPLINKAGE Clp_dualFeasible(Clp_Simplex * model);
  /** Dual bound */
  CLPLIBAPI double CLPLINKAGE Clp_dualBound(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setDualBound(Clp_Simplex * model, double value);
  /** Infeasibility cost */
  CLPLIBAPI double CLPLINKAGE Clp_infeasibilityCost(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setInfeasibilityCost(Clp_Simplex * model, double value);
  /** Perturbation:
      50  - switch on perturbation
      100 - auto perturb if takes too long (1.0e-6 largest nonzero)
      101 - we are perturbed
      102 - don't try perturbing again
      default is 100
      others are for playing
  */
  CLPLIBAPI int CLPLINKAGE Clp_perturbation(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setPerturbation(Clp_Simplex * model, int value);
  /** Current (or last) algorithm */
  CLPLIBAPI int CLPLINKAGE Clp_algorithm(Clp_Simplex * model); 
  /** Set algorithm */
  CLPLIBAPI void CLPLINKAGE Clp_setAlgorithm(Clp_Simplex * model, int value);
  /** Sum of dual infeasibilities */
  CLPLIBAPI double CLPLINKAGE Clp_sumDualInfeasibilities(Clp_Simplex * model); 
  /** Number of dual infeasibilities */
  CLPLIBAPI int CLPLINKAGE Clp_numberDualInfeasibilities(Clp_Simplex * model); 
  /** Sum of primal infeasibilities */
  CLPLIBAPI double CLPLINKAGE Clp_sumPrimalInfeasibilities(Clp_Simplex * model); 
  /** Number of primal infeasibilities */
  CLPLIBAPI int CLPLINKAGE Clp_numberPrimalInfeasibilities(Clp_Simplex * model); 
  /** Save model to file, returns 0 if success.  This is designed for
      use outside algorithms so does not save iterating arrays etc.
  It does not save any messaging information. 
  Does not save scaling values.
  It does not know about all types of virtual functions.
  */
  CLPLIBAPI int CLPLINKAGE Clp_saveModel(Clp_Simplex * model, const char * fileName);
  /** Restore model from file, returns 0 if success,
      deletes current model */
  CLPLIBAPI int CLPLINKAGE Clp_restoreModel(Clp_Simplex * model, const char * fileName);
  
  /** Just check solution (for external use) - sets sum of
      infeasibilities etc */
  CLPLIBAPI void CLPLINKAGE Clp_checkSolution(Clp_Simplex * model);
  /*@}*/

  /******************** End of most useful part **************/
  /**@name gets and sets - some synonyms */
  /*@{*/ 
  /** Number of rows */
  CLPLIBAPI int CLPLINKAGE Clp_getNumRows(Clp_Simplex * model);
  /** Number of columns */
  CLPLIBAPI int CLPLINKAGE Clp_getNumCols(Clp_Simplex * model);
  /** Number of iterations */
  CLPLIBAPI int CLPLINKAGE Clp_getIterationCount(Clp_Simplex * model);
  /** Are there a numerical difficulties? */
  CLPLIBAPI int CLPLINKAGE Clp_isAbandoned(Clp_Simplex * model);
  /** Is optimality proven? */
  CLPLIBAPI int CLPLINKAGE Clp_isProvenOptimal(Clp_Simplex * model);
  /** Is primal infeasiblity proven? */
  CLPLIBAPI int CLPLINKAGE Clp_isProvenPrimalInfeasible(Clp_Simplex * model);
  /** Is dual infeasiblity proven? */
  CLPLIBAPI int CLPLINKAGE Clp_isProvenDualInfeasible(Clp_Simplex * model);
  /** Is the given primal objective limit reached? */
  CLPLIBAPI int CLPLINKAGE Clp_isPrimalObjectiveLimitReached(Clp_Simplex * model) ;
  /** Is the given dual objective limit reached? */
  CLPLIBAPI int CLPLINKAGE Clp_isDualObjectiveLimitReached(Clp_Simplex * model) ;
  /** Iteration limit reached? */
  CLPLIBAPI int CLPLINKAGE Clp_isIterationLimitReached(Clp_Simplex * model);
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  CLPLIBAPI double CLPLINKAGE Clp_getObjSense(Clp_Simplex * model);
  /** Primal row solution */
  CLPLIBAPI const double * CLPLINKAGE Clp_getRowActivity(Clp_Simplex * model);
  /** Primal column solution */
  CLPLIBAPI const double * CLPLINKAGE Clp_getColSolution(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setColSolution(Clp_Simplex * model, const double * input);
  /** Dual row solution */
  CLPLIBAPI const double * CLPLINKAGE Clp_getRowPrice(Clp_Simplex * model);
  /** Reduced costs */
  CLPLIBAPI const double * CLPLINKAGE Clp_getReducedCost(Clp_Simplex * model);
  /** Row lower */
  CLPLIBAPI const double* CLPLINKAGE Clp_getRowLower(Clp_Simplex * model);
  /** Row upper  */
  CLPLIBAPI const double* CLPLINKAGE Clp_getRowUpper(Clp_Simplex * model);
  /** Objective */
  CLPLIBAPI const double * CLPLINKAGE Clp_getObjCoefficients(Clp_Simplex * model); 
  /** Column Lower */
  CLPLIBAPI const double * CLPLINKAGE Clp_getColLower(Clp_Simplex * model);
  /** Column Upper */
  CLPLIBAPI const double * CLPLINKAGE Clp_getColUpper(Clp_Simplex * model);
  /** Objective value */
  CLPLIBAPI double CLPLINKAGE Clp_getObjValue(Clp_Simplex * model);
  /* Small element value - elements less than this set to zero,
     default is 1.0e-20 */
  CLPLIBAPI double CLPLINKAGE Clp_getSmallElementValue(Clp_Simplex * model);
  CLPLIBAPI void CLPLINKAGE Clp_setSmallElementValue(Clp_Simplex * model,double value);
  /*@}*/
#ifdef __cplusplus
          }
#endif
#endif
