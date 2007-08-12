/* Copyright (C) 2002, 2003 International Business Machines
   Corporation and others.  All Rights Reserved.*/
#ifndef ClpSimplexC_H
#define ClpSimplexC_H

/* include all defines and ugly stuff */
#include "Coin_C_defines.h"

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
  COINLIBAPI Clp_Simplex * COINLINKAGE Clp_newModel(void);
  /** Destructor */
  COINLIBAPI void COINLINKAGE Clp_deleteModel(Clp_Simplex * model);
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
  COINLIBAPI void COINLINKAGE Clp_loadProblem (Clp_Simplex * model,  const int numcols, const int numrows,
		      const CoinBigIndex * start, const int* index,
		      const double* value,
		      const double* collb, const double* colub,   
		      const double* obj,
		      const double* rowlb, const double* rowub);

/* read quadratic part of the objective (the matrix part) */
COINLIBAPI void COINLINKAGE 
Clp_loadQuadraticObjective(Clp_Simplex * model,  
                           const int numberColumns, 
                           const CoinBigIndex * start,
                           const int * column, 
                           const double * element);
  /** Read an mps file from the given filename */
  COINLIBAPI int COINLINKAGE Clp_readMps(Clp_Simplex * model,const char *filename,
	      int keepNames,
	      int ignoreErrors);
  /** Copy in integer informations */
  COINLIBAPI void COINLINKAGE Clp_copyInIntegerInformation(Clp_Simplex * model,const char * information);
  /** Drop integer informations */
  COINLIBAPI void COINLINKAGE Clp_deleteIntegerInformation(Clp_Simplex * model);
  /** Resizes rim part of model  */
  COINLIBAPI void COINLINKAGE Clp_resize (Clp_Simplex * model, int newNumberRows, int newNumberColumns);
  /** Deletes rows */
  COINLIBAPI void COINLINKAGE Clp_deleteRows(Clp_Simplex * model, int number, const int * which);
  /** Add rows */
  COINLIBAPI void COINLINKAGE Clp_addRows(Clp_Simplex * model, int number, const double * rowLower, 
	       const double * rowUpper,
	       const int * rowStarts, const int * columns,
	       const double * elements);

  /** Deletes columns */
  COINLIBAPI void COINLINKAGE Clp_deleteColumns(Clp_Simplex * model, int number, const int * which);
  /** Add columns */
  COINLIBAPI void COINLINKAGE Clp_addColumns(Clp_Simplex * model, int number, const double * columnLower, 
		  const double * columnUpper,
		  const double * objective,
		  const int * columnStarts, const int * rows,
		  const double * elements);
  /** Change row lower bounds */
  COINLIBAPI void COINLINKAGE Clp_chgRowLower(Clp_Simplex * model, const double * rowLower);
  /** Change row upper bounds */
  COINLIBAPI void COINLINKAGE Clp_chgRowUpper(Clp_Simplex * model, const double * rowUpper); 
  /** Change column lower bounds */
  COINLIBAPI void COINLINKAGE Clp_chgColumnLower(Clp_Simplex * model, const double * columnLower); 
  /** Change column upper bounds */
  COINLIBAPI void COINLINKAGE Clp_chgColumnUpper(Clp_Simplex * model, const double * columnUpper); 
  /** Change objective coefficients */
  COINLIBAPI void COINLINKAGE Clp_chgObjCoefficients(Clp_Simplex * model, const double * objIn); 
  /** Drops names - makes lengthnames 0 and names empty */
  COINLIBAPI void COINLINKAGE Clp_dropNames(Clp_Simplex * model);
  /** Copies in names */
  COINLIBAPI void COINLINKAGE Clp_copyNames(Clp_Simplex * model, const char * const * rowNames,
		 const char * const * columnNames);
  
  /*@}*/
  /**@name gets and sets - you will find some synonyms at the end of this file */
  /*@{*/ 
  /** Number of rows */
  COINLIBAPI int COINLINKAGE Clp_numberRows(Clp_Simplex * model);
  /** Number of columns */
  COINLIBAPI int COINLINKAGE Clp_numberColumns(Clp_Simplex * model);
  /** Primal tolerance to use */
  COINLIBAPI double COINLINKAGE Clp_primalTolerance(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setPrimalTolerance(Clp_Simplex * model,  double value) ;
  /** Dual tolerance to use */
  COINLIBAPI double COINLINKAGE Clp_dualTolerance(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setDualTolerance(Clp_Simplex * model,  double value) ;
  /** Dual objective limit */
  COINLIBAPI double COINLINKAGE Clp_dualObjectiveLimit(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setDualObjectiveLimit(Clp_Simplex * model, double value);
  /** Objective offset */
  COINLIBAPI double COINLINKAGE Clp_objectiveOffset(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setObjectiveOffset(Clp_Simplex * model, double value);
  /** Fills in array with problem name  */
  COINLIBAPI void COINLINKAGE Clp_problemName(Clp_Simplex * model, int maxNumberCharacters, char * array);
  /* Sets problem name.  Must have \0 at end.  */
  COINLIBAPI int COINLINKAGE 
    Clp_setProblemName(Clp_Simplex * model, int maxNumberCharacters, char * array);
  /** Number of iterations */
  COINLIBAPI int COINLINKAGE Clp_numberIterations(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setNumberIterations(Clp_Simplex * model, int numberIterations);
  /** Maximum number of iterations */
  COINLIBAPI int maximumIterations(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setMaximumIterations(Clp_Simplex * model, int value);
  /** Maximum time in seconds (from when set called) */
  COINLIBAPI double COINLINKAGE Clp_maximumSeconds(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setMaximumSeconds(Clp_Simplex * model, double value);
  /** Returns true if hit maximum iterations (or time) */
  COINLIBAPI int COINLINKAGE Clp_hitMaximumIterations(Clp_Simplex * model);
  /** Status of problem:
      0 - optimal
      1 - primal infeasible
      2 - dual infeasible
      3 - stopped on iterations etc
      4 - stopped due to errors
  */
  COINLIBAPI int COINLINKAGE Clp_status(Clp_Simplex * model);
  /** Set problem status */
  COINLIBAPI void COINLINKAGE Clp_setProblemStatus(Clp_Simplex * model, int problemStatus);
  /** Secondary status of problem - may get extended
      0 - none
      1 - primal infeasible because dual limit reached
      2 - scaled problem optimal - unscaled has primal infeasibilities
      3 - scaled problem optimal - unscaled has dual infeasibilities
      4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
  */
  COINLIBAPI int COINLINKAGE Clp_secondaryStatus(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setSecondaryStatus(Clp_Simplex * model, int status);
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  COINLIBAPI double COINLINKAGE Clp_optimizationDirection(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setOptimizationDirection(Clp_Simplex * model, double value);
  /** Primal row solution */
  COINLIBAPI double * COINLINKAGE Clp_primalRowSolution(Clp_Simplex * model);
  /** Primal column solution */
  COINLIBAPI double * COINLINKAGE Clp_primalColumnSolution(Clp_Simplex * model);
  /** Dual row solution */
  COINLIBAPI double * COINLINKAGE Clp_dualRowSolution(Clp_Simplex * model);
  /** Reduced costs */
  COINLIBAPI double * COINLINKAGE Clp_dualColumnSolution(Clp_Simplex * model);
  /** Row lower */
  COINLIBAPI double* COINLINKAGE Clp_rowLower(Clp_Simplex * model);
  /** Row upper  */
  COINLIBAPI double* COINLINKAGE Clp_rowUpper(Clp_Simplex * model);
  /** Objective */
  COINLIBAPI double * COINLINKAGE Clp_objective(Clp_Simplex * model);            
  /** Column Lower */
  COINLIBAPI double * COINLINKAGE Clp_columnLower(Clp_Simplex * model);
  /** Column Upper */
  COINLIBAPI double * COINLINKAGE Clp_columnUpper(Clp_Simplex * model);
  /** Number of elements in matrix */
  COINLIBAPI int COINLINKAGE Clp_getNumElements(Clp_Simplex * model); 
  /* Column starts in matrix */
  COINLIBAPI const CoinBigIndex * COINLINKAGE Clp_getVectorStarts(Clp_Simplex * model); 
  /* Row indices in matrix */
  COINLIBAPI const int * COINLINKAGE Clp_getIndices(Clp_Simplex * model); 
  /* Column vector lengths in matrix */
  COINLIBAPI const int * COINLINKAGE Clp_getVectorLengths(Clp_Simplex * model); 
  /* Element values in matrix */
  COINLIBAPI const double * COINLINKAGE Clp_getElements(Clp_Simplex * model); 
  /** Objective value */
  COINLIBAPI double COINLINKAGE Clp_objectiveValue(Clp_Simplex * model);
  /** Integer information */
  COINLIBAPI char * COINLINKAGE Clp_integerInformation(Clp_Simplex * model);
  /** Infeasibility/unbounded ray (NULL returned if none/wrong)
      Up to user to use delete [] on these arrays.  */
  COINLIBAPI double * COINLINKAGE Clp_infeasibilityRay(Clp_Simplex * model);
  COINLIBAPI double * COINLINKAGE Clp_unboundedRay(Clp_Simplex * model);
  /** See if status array exists (partly for OsiClp) */
  COINLIBAPI int COINLINKAGE Clp_statusExists(Clp_Simplex * model);
  /** Return address of status array (char[numberRows+numberColumns]) */
  COINLIBAPI unsigned char *  COINLINKAGE Clp_statusArray(Clp_Simplex * model);
  /** Copy in status vector */
  COINLIBAPI void COINLINKAGE Clp_copyinStatus(Clp_Simplex * model, const unsigned char * statusArray);
  /* status values are as in ClpSimplex.hpp i.e. 0 - free, 1 basic, 2 at upper,
     3 at lower, 4 superbasic, (5 fixed) */
  /* Get variable basis info */
  COINLIBAPI int COINLINKAGE Clp_getColumnStatus(Clp_Simplex * model,int sequence);
  /* Get row basis info */
  COINLIBAPI int COINLINKAGE Clp_getRowStatus(Clp_Simplex * model,int sequence);
  /* Set variable basis info (and value if at bound) */
  COINLIBAPI void COINLINKAGE Clp_setColumnStatus(Clp_Simplex * model,
						  int sequence, int value);
  /* Set row basis info (and value if at bound) */
  COINLIBAPI void COINLINKAGE Clp_setRowStatus(Clp_Simplex * model,
					       int sequence, int value);
  
  /** User pointer for whatever reason */
  COINLIBAPI void COINLINKAGE Clp_setUserPointer (Clp_Simplex * model, void * pointer);
  COINLIBAPI void * COINLINKAGE Clp_getUserPointer (Clp_Simplex * model);
  /*@}*/
  /**@name Message handling.  Call backs are handled by ONE function */
  /*@{*/
  /** Pass in Callback function.
   Message numbers up to 1000000 are Clp, Coin ones have 1000000 added */
  COINLIBAPI void COINLINKAGE Clp_registerCallBack(Clp_Simplex * model, 
						   clp_callback userCallBack);
  /** Unset Callback function */
  COINLIBAPI void COINLINKAGE Clp_clearCallBack(Clp_Simplex * model);
  /** Amount of print out:
      0 - none
      1 - just final
      2 - just factorizations
      3 - as 2 plus a bit more
      4 - verbose
      above that 8,16,32 etc just for selective debug
  */
  COINLIBAPI void COINLINKAGE Clp_setLogLevel(Clp_Simplex * model, int value);
  COINLIBAPI int COINLINKAGE Clp_logLevel(Clp_Simplex * model);
  /** length of names (0 means no names0 */
  COINLIBAPI int COINLINKAGE Clp_lengthNames(Clp_Simplex * model);
  /** Fill in array (at least lengthNames+1 long) with a row name */
  COINLIBAPI void COINLINKAGE Clp_rowName(Clp_Simplex * model, int iRow, char * name);
  /** Fill in array (at least lengthNames+1 long) with a column name */
  COINLIBAPI void COINLINKAGE Clp_columnName(Clp_Simplex * model, int iColumn, char * name);

  /*@}*/


  /**@name Functions most useful to user */
  /*@{*/
  /** General solve algorithm which can do presolve.
      See  ClpSolve.hpp for options
   */
  COINLIBAPI int COINLINKAGE Clp_initialSolve(Clp_Simplex * model);
  /** Dual initial solve */
  COINLIBAPI int COINLINKAGE Clp_initialDualSolve(Clp_Simplex * model);
  /** Primal initial solve */
  COINLIBAPI int COINLINKAGE Clp_initialPrimalSolve(Clp_Simplex * model);
  /** Barrier initial solve */
  COINLIBAPI int COINLINKAGE Clp_initialBarrierSolve(Clp_Simplex * model);
  /** Barrier initial solve, no crossover */
  COINLIBAPI int COINLINKAGE Clp_initialBarrierNoCrossSolve(Clp_Simplex * model);
  /** Dual algorithm - see ClpSimplexDual.hpp for method */
  COINLIBAPI int COINLINKAGE Clp_dual(Clp_Simplex * model, int ifValuesPass);
  /** Primal algorithm - see ClpSimplexPrimal.hpp for method */
  COINLIBAPI int COINLINKAGE Clp_primal(Clp_Simplex * model, int ifValuesPass);
  /** Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
  COINLIBAPI void COINLINKAGE Clp_scaling(Clp_Simplex * model, int mode);
  /** Gets scalingFlag */
  COINLIBAPI int COINLINKAGE Clp_scalingFlag(Clp_Simplex * model);
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
  COINLIBAPI int COINLINKAGE Clp_crash(Clp_Simplex * model, double gap,int pivot);
  /*@}*/


  /**@name most useful gets and sets */
  /*@{*/ 
  /** If problem is primal feasible */
  COINLIBAPI int COINLINKAGE Clp_primalFeasible(Clp_Simplex * model);
  /** If problem is dual feasible */
  COINLIBAPI int COINLINKAGE Clp_dualFeasible(Clp_Simplex * model);
  /** Dual bound */
  COINLIBAPI double COINLINKAGE Clp_dualBound(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setDualBound(Clp_Simplex * model, double value);
  /** Infeasibility cost */
  COINLIBAPI double COINLINKAGE Clp_infeasibilityCost(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setInfeasibilityCost(Clp_Simplex * model, double value);
  /** Perturbation:
      50  - switch on perturbation
      100 - auto perturb if takes too long (1.0e-6 largest nonzero)
      101 - we are perturbed
      102 - don't try perturbing again
      default is 100
      others are for playing
  */
  COINLIBAPI int COINLINKAGE Clp_perturbation(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setPerturbation(Clp_Simplex * model, int value);
  /** Current (or last) algorithm */
  COINLIBAPI int COINLINKAGE Clp_algorithm(Clp_Simplex * model); 
  /** Set algorithm */
  COINLIBAPI void COINLINKAGE Clp_setAlgorithm(Clp_Simplex * model, int value);
  /** Sum of dual infeasibilities */
  COINLIBAPI double COINLINKAGE Clp_sumDualInfeasibilities(Clp_Simplex * model); 
  /** Number of dual infeasibilities */
  COINLIBAPI int COINLINKAGE Clp_numberDualInfeasibilities(Clp_Simplex * model); 
  /** Sum of primal infeasibilities */
  COINLIBAPI double COINLINKAGE Clp_sumPrimalInfeasibilities(Clp_Simplex * model); 
  /** Number of primal infeasibilities */
  COINLIBAPI int COINLINKAGE Clp_numberPrimalInfeasibilities(Clp_Simplex * model); 
  /** Save model to file, returns 0 if success.  This is designed for
      use outside algorithms so does not save iterating arrays etc.
  It does not save any messaging information. 
  Does not save scaling values.
  It does not know about all types of virtual functions.
  */
  COINLIBAPI int COINLINKAGE Clp_saveModel(Clp_Simplex * model, const char * fileName);
  /** Restore model from file, returns 0 if success,
      deletes current model */
  COINLIBAPI int COINLINKAGE Clp_restoreModel(Clp_Simplex * model, const char * fileName);
  
  /** Just check solution (for external use) - sets sum of
      infeasibilities etc */
  COINLIBAPI void COINLINKAGE Clp_checkSolution(Clp_Simplex * model);
  /*@}*/

  /******************** End of most useful part **************/
  /**@name gets and sets - some synonyms */
  /*@{*/ 
  /** Number of rows */
  COINLIBAPI int COINLINKAGE Clp_getNumRows(Clp_Simplex * model);
  /** Number of columns */
  COINLIBAPI int COINLINKAGE Clp_getNumCols(Clp_Simplex * model);
  /** Number of iterations */
  COINLIBAPI int COINLINKAGE Clp_getIterationCount(Clp_Simplex * model);
  /** Are there a numerical difficulties? */
  COINLIBAPI int COINLINKAGE Clp_isAbandoned(Clp_Simplex * model);
  /** Is optimality proven? */
  COINLIBAPI int COINLINKAGE Clp_isProvenOptimal(Clp_Simplex * model);
  /** Is primal infeasiblity proven? */
  COINLIBAPI int COINLINKAGE Clp_isProvenPrimalInfeasible(Clp_Simplex * model);
  /** Is dual infeasiblity proven? */
  COINLIBAPI int COINLINKAGE Clp_isProvenDualInfeasible(Clp_Simplex * model);
  /** Is the given primal objective limit reached? */
  COINLIBAPI int COINLINKAGE Clp_isPrimalObjectiveLimitReached(Clp_Simplex * model) ;
  /** Is the given dual objective limit reached? */
  COINLIBAPI int COINLINKAGE Clp_isDualObjectiveLimitReached(Clp_Simplex * model) ;
  /** Iteration limit reached? */
  COINLIBAPI int COINLINKAGE Clp_isIterationLimitReached(Clp_Simplex * model);
  /** Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
  COINLIBAPI double COINLINKAGE Clp_getObjSense(Clp_Simplex * model);
  /** Primal row solution */
  COINLIBAPI const double * COINLINKAGE Clp_getRowActivity(Clp_Simplex * model);
  /** Primal column solution */
  COINLIBAPI const double * COINLINKAGE Clp_getColSolution(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setColSolution(Clp_Simplex * model, const double * input);
  /** Dual row solution */
  COINLIBAPI const double * COINLINKAGE Clp_getRowPrice(Clp_Simplex * model);
  /** Reduced costs */
  COINLIBAPI const double * COINLINKAGE Clp_getReducedCost(Clp_Simplex * model);
  /** Row lower */
  COINLIBAPI const double* COINLINKAGE Clp_getRowLower(Clp_Simplex * model);
  /** Row upper  */
  COINLIBAPI const double* COINLINKAGE Clp_getRowUpper(Clp_Simplex * model);
  /** Objective */
  COINLIBAPI const double * COINLINKAGE Clp_getObjCoefficients(Clp_Simplex * model); 
  /** Column Lower */
  COINLIBAPI const double * COINLINKAGE Clp_getColLower(Clp_Simplex * model);
  /** Column Upper */
  COINLIBAPI const double * COINLINKAGE Clp_getColUpper(Clp_Simplex * model);
  /** Objective value */
  COINLIBAPI double COINLINKAGE Clp_getObjValue(Clp_Simplex * model);
  /** Print model for debugging purposes */
  COINLIBAPI void COINLINKAGE Clp_printModel(Clp_Simplex * model, const char * prefix);
  /* Small element value - elements less than this set to zero,
     default is 1.0e-20 */
  COINLIBAPI double COINLINKAGE Clp_getSmallElementValue(Clp_Simplex * model);
  COINLIBAPI void COINLINKAGE Clp_setSmallElementValue(Clp_Simplex * model,double value);
  /*@}*/
#ifdef __cplusplus
          }
#endif
#endif
