// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.




#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include <cfloat>
// This sections needs to match Clp_C_defines.h but with extern C
#define ClpSimplexCDefine_H

/** This has #defines etc for the "C" interface to Clp.

*/

/* Plus infinity */
#ifndef COIN_DBL_MAX
#define COIN_DBL_MAX DBL_MAX
#endif

/* We need to allow for Microsoft */
#ifndef CLPLIBAPI

#if defined (CLPMSDLL)
#   define CLPLIBAPI __declspec(dllexport) extern "C"
#   define CLPLINKAGE  __stdcall
#   define CLPLINKAGE_CB  __cdecl
#else
#   define CLPLIBAPI extern "C"
#   define CLPLINKAGE
#   define CLPLINKAGE_CB 
#endif

#endif
class CMessageHandler;
// Real typedef for structure
typedef struct {
  ClpSimplex * model_;
  CMessageHandler * handler_;
} Clp_Simplex;
/** typedef for user call back */
typedef  void (CLPLINKAGE_CB *clp_callback) (Clp_Simplex * model,int  msgno, int ndouble,
                            const double * dvec, int nint, const int * ivec,
                            int nchar,  char ** cvec);

// To allow call backs
class CMessageHandler : public CoinMessageHandler {
  
public:
  /**@name Overrides */
  //@{
  virtual int print();
  //@}
  /**@name set and get */
  //@{
  /// Model
  const Clp_Simplex * model() const;
  void setModel(Clp_Simplex * model);
  /// Call back
  void setCallBack(clp_callback callback);
  //@}

  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  CMessageHandler();
  /// Constructor with pointer to model
  CMessageHandler(Clp_Simplex * model,
			   FILE * userPointer=NULL);
  /** Destructor */
  virtual ~CMessageHandler();
  //@}

  /**@name Copy method */
  //@{
  /** The copy constructor. */
  CMessageHandler(const CMessageHandler&);
  /** The copy constructor from an CoinSimplexMessageHandler. */
  CMessageHandler(const CoinMessageHandler&);
  
  CMessageHandler& operator=(const CMessageHandler&);
  /// Clone
  virtual CoinMessageHandler * clone() const ;
  //@}
   
    
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer back to model
  Clp_Simplex * model_;
  /// call back
  clp_callback callback_;
  //@}
};


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CMessageHandler::CMessageHandler () 
  : CoinMessageHandler(),
    model_(NULL),
    callback_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CMessageHandler::CMessageHandler (const CMessageHandler & rhs) 
: CoinMessageHandler(rhs),
    model_(rhs.model_),
    callback_(rhs.callback_)
{  
}

CMessageHandler::CMessageHandler (const CoinMessageHandler & rhs) 
  : CoinMessageHandler(),
    model_(NULL),
    callback_(NULL)
{  
}

// Constructor with pointer to model
CMessageHandler::CMessageHandler(Clp_Simplex * model,
               FILE * userPointer)
  : CoinMessageHandler(),
    model_(model),
    callback_(NULL)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CMessageHandler::~CMessageHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CMessageHandler &
CMessageHandler::operator=(const CMessageHandler& rhs)
{
  if (this != &rhs) {
    CoinMessageHandler::operator=(rhs);
    model_ = rhs.model_;
    callback_ = rhs.callback_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CoinMessageHandler * CMessageHandler::clone() const
{
  return new CMessageHandler(*this);
}

int 
CMessageHandler::print()
{
  if (callback_) {
    int messageNumber = currentMessage().externalNumber();
    if (currentSource()!="Clp")
      messageNumber += 1000000;
    int i;
    int nDouble=numberDoubleFields();
    assert (nDouble<=10);
    double vDouble[10];
    for (i=0;i<nDouble;i++)
      vDouble[i]=doubleValue(i);
    int nInt=numberIntFields();
    assert (nInt<=10);
    int vInt[10];
    for (i=0;i<nInt;i++)
      vInt[i]=intValue(i);
    int nString=numberStringFields();
    assert (nString<=10);
    char * vString[10];
    for (i=0;i<nString;i++) {
      std::string value = stringValue(i);
      vString[i]=strdup(value.c_str());
    }
    callback_(model_,messageNumber,
	      nDouble,vDouble,
	      nInt,vInt,
	      nString,vString);
    for (i=0;i<nString;i++) 
      free(vString[i]);
    
  }
  return CoinMessageHandler::print();
}
const Clp_Simplex *
CMessageHandler::model() const
{
  return model_;
}
void 
CMessageHandler::setModel(Clp_Simplex * model)
{
  model_ = model;
}
// Call back
void 
CMessageHandler::setCallBack(clp_callback callback)
{
  callback_ = callback;
}

#include "Clp_C_Interface.h"
#include <string>
#include <stdio.h>
#include <iostream>

/* Default constructor */
CLPLIBAPI Clp_Simplex *  CLPLINKAGE 
Clp_newModel()
{
  Clp_Simplex * model = new Clp_Simplex;
  model->model_ = new ClpSimplex();
  model->handler_=NULL;
  return model;
}
/* Destructor */
CLPLIBAPI void CLPLINKAGE 
Clp_deleteModel(Clp_Simplex * model)
{
  delete model->model_;
  delete model->handler_;
  delete model;
}

/* Loads a problem (the constraints on the
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
/* Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). */
CLPLIBAPI void CLPLINKAGE 
Clp_loadProblem (Clp_Simplex * model,  const int numcols, const int numrows,
		 const CoinBigIndex * start, const int* index,
		 const double* value,
		 const double* collb, const double* colub,   
		 const double* obj,
		 const double* rowlb, const double* rowub)
{
  model->model_->loadProblem(numcols,numrows,start,index,value,
			     collb,colub,obj,rowlb,rowub);
}
/* Read an mps file from the given filename */
CLPLIBAPI int CLPLINKAGE 
Clp_readMps(Clp_Simplex * model,const char *filename,
	    int keepNames,
	    int ignoreErrors)
{
  return model->model_->readMps(filename,keepNames!=0,ignoreErrors!=0);
}
/* Copy in integer informations */
CLPLIBAPI void CLPLINKAGE 
Clp_copyInIntegerInformation(Clp_Simplex * model,const char * information)
{
  model->model_->copyInIntegerInformation(information);
}
/* Drop integer informations */
CLPLIBAPI void CLPLINKAGE 
Clp_deleteIntegerInformation(Clp_Simplex * model)
{
  model->model_->deleteIntegerInformation();
}
/* Resizes rim part of model  */
CLPLIBAPI void CLPLINKAGE 
Clp_resize (Clp_Simplex * model, int newNumberRows, int newNumberColumns)
{
  model->model_->resize(newNumberRows,newNumberColumns);
}
/* Deletes rows */
CLPLIBAPI void CLPLINKAGE 
Clp_deleteRows(Clp_Simplex * model, int number, const int * which)
{
  model->model_->deleteRows(number,which);
}
/* Add rows */
CLPLIBAPI void CLPLINKAGE 
Clp_addRows(Clp_Simplex * model, int number, const double * rowLower, 
	       const double * rowUpper,
	       const int * rowStarts, const int * columns,
	       const double * elements)
{
  model->model_->addRows(number,rowLower,rowUpper,rowStarts,columns,elements);
}

/* Deletes columns */
CLPLIBAPI void CLPLINKAGE 
Clp_deleteColumns(Clp_Simplex * model, int number, const int * which)
{
  model->model_->deleteColumns(number,which);
}
/* Add columns */
CLPLIBAPI void CLPLINKAGE 
Clp_addColumns(Clp_Simplex * model, int number, const double * columnLower, 
		  const double * columnUpper,
		  const double * objective,
		  const int * columnStarts, const int * rows,
		  const double * elements)
{
  model->model_->addColumns(number,columnLower,columnUpper,objective,
			    columnStarts,rows,elements);
}
/* Change row lower bounds */
CLPLIBAPI void CLPLINKAGE 
Clp_chgRowLower(Clp_Simplex * model, const double * rowLower) 
{
  model->model_->chgRowLower(rowLower);
}
/* Change row upper bounds */
CLPLIBAPI void CLPLINKAGE 
Clp_chgRowUpper(Clp_Simplex * model, const double * rowUpper) 
{
  model->model_->chgRowUpper(rowUpper);
}
/* Change column lower bounds */
CLPLIBAPI void CLPLINKAGE 
Clp_chgColumnLower(Clp_Simplex * model, const double * columnLower) 
{
  model->model_->chgColumnLower(columnLower);
}
/* Change column upper bounds */
CLPLIBAPI void CLPLINKAGE 
Clp_chgColumnUpper(Clp_Simplex * model, const double * columnUpper) 
{
  model->model_->chgColumnUpper(columnUpper);
}
/* Change objective coefficients */
CLPLIBAPI void CLPLINKAGE 
Clp_chgObjCoefficients(Clp_Simplex * model, const double * objIn) 
{
  model->model_->chgObjCoefficients(objIn);
}
/* Drops names - makes lengthnames 0 and names empty */
CLPLIBAPI void CLPLINKAGE 
Clp_dropNames(Clp_Simplex * model)
{
  model->model_->dropNames();
}
/* Copies in names */
CLPLIBAPI void CLPLINKAGE 
Clp_copyNames(Clp_Simplex * model, const char * const * rowNamesIn,
	      const char * const * columnNamesIn)
{
  int iRow;
  std::vector<std::string> rowNames;
  int numberRows = model->model_->numberRows();
  rowNames.reserve(numberRows);
  for (iRow=0;iRow<numberRows;iRow++) {
    rowNames.push_back(rowNamesIn[iRow]);
  }
  
  int iColumn;
  std::vector<std::string> columnNames;
  int numberColumns = model->model_->numberColumns();
  columnNames.reserve(numberColumns);
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    columnNames.push_back(columnNamesIn[iColumn]);
  }
  model->model_->copyNames(rowNames,columnNames);
}

/* Number of rows */
CLPLIBAPI int CLPLINKAGE 
Clp_numberRows(Clp_Simplex * model)
{
  return model->model_->numberRows();
}
/* Number of columns */
CLPLIBAPI int CLPLINKAGE 
Clp_numberColumns(Clp_Simplex * model)
{
  return model->model_->numberColumns();
}
/* Primal tolerance to use */
CLPLIBAPI double CLPLINKAGE 
Clp_primalTolerance(Clp_Simplex * model)
{
  return model->model_->primalTolerance();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setPrimalTolerance(Clp_Simplex * model,  double value) 
{
  model->model_->setPrimalTolerance(value);
}
/* Dual tolerance to use */
CLPLIBAPI double CLPLINKAGE 
Clp_dualTolerance(Clp_Simplex * model)
{
  return model->model_->dualTolerance();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setDualTolerance(Clp_Simplex * model,  double value) 
{
  model->model_->setDualTolerance(value);
}
/* Dual objective limit */
CLPLIBAPI double CLPLINKAGE 
Clp_dualObjectiveLimit(Clp_Simplex * model)
{
  return model->model_->dualObjectiveLimit();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setDualObjectiveLimit(Clp_Simplex * model, double value)
{
  model->model_->setDualObjectiveLimit(value);
}
/* Objective offset */
CLPLIBAPI double CLPLINKAGE 
Clp_objectiveOffset(Clp_Simplex * model)
{
  return model->model_->objectiveOffset();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setObjectiveOffset(Clp_Simplex * model, double value)
{
  model->model_->setObjectiveOffset(value);
}
/* Fills in array with problem name  */
CLPLIBAPI void CLPLINKAGE 
Clp_problemName(Clp_Simplex * model, int maxNumberCharacters, char * array)
{
  std::string name = model->model_->problemName();
  maxNumberCharacters = CoinMin(maxNumberCharacters,(int)strlen(name.c_str()));
  strncpy(array,name.c_str(),maxNumberCharacters-1);
  array[maxNumberCharacters-1]='\0';
}
/* Number of iterations */
CLPLIBAPI int CLPLINKAGE 
Clp_numberIterations(Clp_Simplex * model)
{
  return model->model_->numberIterations();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setNumberIterations(Clp_Simplex * model, int numberIterations)
{
  model->model_->setNumberIterations(numberIterations);
}
/* Maximum number of iterations */
CLPLIBAPI int maximumIterations(Clp_Simplex * model)
{
  return model->model_->maximumIterations();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setMaximumIterations(Clp_Simplex * model, int value)
{
  model->model_->setMaximumIterations(value);
}
/* Maximum time in seconds (from when set called) */
CLPLIBAPI double CLPLINKAGE 
Clp_maximumSeconds(Clp_Simplex * model)
{
  return model->model_->maximumSeconds();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setMaximumSeconds(Clp_Simplex * model, double value)
{
  model->model_->setMaximumSeconds(value);
}
/* Returns true if hit maximum iteratio`ns (or time) */
CLPLIBAPI int CLPLINKAGE 
Clp_hitMaximumIterations(Clp_Simplex * model)
{
  return model->model_->hitMaximumIterations() ? 1 : 0;
}
/* Status of problem:
   0 - optimal
   1 - primal infeasible
   2 - dual infeasible
   3 - stopped on iterations etc
   4 - stopped due to errors
*/
CLPLIBAPI int CLPLINKAGE 
Clp_status(Clp_Simplex * model)
{
  return model->model_->status();
}
/* Set problem status */
CLPLIBAPI void CLPLINKAGE 
Clp_setProblemStatus(Clp_Simplex * model, int problemStatus)
{
  model->model_->setProblemStatus(problemStatus);
}
/* Secondary status of problem - may get extended
   0 - none
   1 - primal infeasible because dual limit reached
   2 - scaled problem optimal - unscaled has primal infeasibilities
   3 - scaled problem optimal - unscaled has dual infeasibilities
   4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
*/
CLPLIBAPI int CLPLINKAGE 
Clp_secondaryStatus(Clp_Simplex * model)
{
  return model->model_->secondaryStatus();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setSecondaryStatus(Clp_Simplex * model, int status)
{
  model->model_->setSecondaryStatus(status);
}
/* Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
CLPLIBAPI double CLPLINKAGE 
Clp_optimizationDirection(Clp_Simplex * model)
{
  return model->model_->optimizationDirection();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setOptimizationDirection(Clp_Simplex * model, double value)
{
  model->model_->setOptimizationDirection(value);
}
/* Primal row solution */
CLPLIBAPI double * CLPLINKAGE 
Clp_primalRowSolution(Clp_Simplex * model)
{
  return model->model_->primalRowSolution();
}
/* Primal column solution */
CLPLIBAPI double * CLPLINKAGE 
Clp_primalColumnSolution(Clp_Simplex * model)
{
  return model->model_->primalColumnSolution();
}
/* Dual row solution */
CLPLIBAPI double * CLPLINKAGE 
Clp_dualRowSolution(Clp_Simplex * model)
{
  return model->model_->dualRowSolution();
}
/* Reduced costs */
CLPLIBAPI double * CLPLINKAGE 
Clp_dualColumnSolution(Clp_Simplex * model)
{
  return model->model_->dualColumnSolution();
}
/* Row lower */
CLPLIBAPI double* CLPLINKAGE 
Clp_rowLower(Clp_Simplex * model)
{
  return model->model_->rowLower();
}
/* Row upper  */
CLPLIBAPI double* CLPLINKAGE 
Clp_rowUpper(Clp_Simplex * model)
{
  return model->model_->rowUpper();
}
/* Objective */
CLPLIBAPI double * CLPLINKAGE 
Clp_objective(Clp_Simplex * model)
{
  return model->model_->objective();
}
/* Column Lower */
CLPLIBAPI double * CLPLINKAGE 
Clp_columnLower(Clp_Simplex * model)
{
  return model->model_->columnLower();
}
/* Column Upper */
CLPLIBAPI double * CLPLINKAGE 
Clp_columnUpper(Clp_Simplex * model)
{
  return model->model_->columnUpper();
}
/* Number of elements in matrix */
CLPLIBAPI int CLPLINKAGE 
Clp_getNumElements(Clp_Simplex * model)
{
  return model->model_->getNumElements();
}
/* Objective value */
CLPLIBAPI double CLPLINKAGE 
Clp_objectiveValue(Clp_Simplex * model)
{
  return model->model_->objectiveValue();
}
/* Integer information */
CLPLIBAPI char * CLPLINKAGE 
Clp_integerInformation(Clp_Simplex * model)
{
  return model->model_->integerInformation();
}
/* Infeasibility/unbounded ray (NULL returned if none/wrong)
   Up to user to use delete [] on these arrays.  */
CLPLIBAPI double * CLPLINKAGE 
Clp_infeasibilityRay(Clp_Simplex * model)
{
  return model->model_->infeasibilityRay();
}
CLPLIBAPI double * CLPLINKAGE 
Clp_unboundedRay(Clp_Simplex * model)
{
  return model->model_->unboundedRay();
}
/* See if status array exists (partly for OsiClp) */
CLPLIBAPI int CLPLINKAGE 
Clp_statusExists(Clp_Simplex * model)
{
  return model->model_->statusExists() ? 1 : 0;
}
/* Return address of status array (char[numberRows+numberColumns]) */
CLPLIBAPI unsigned char *  CLPLINKAGE 
Clp_statusArray(Clp_Simplex * model)
{
  return model->model_->statusArray();
}
/* Copy in status vector */
CLPLIBAPI void CLPLINKAGE 
Clp_copyinStatus(Clp_Simplex * model, const unsigned char * statusArray)
{
  model->model_->copyinStatus(statusArray);
}

/* User pointer for whatever reason */
CLPLIBAPI void CLPLINKAGE 
Clp_setUserPointer (Clp_Simplex * model, void * pointer)
{
  model->model_->setUserPointer(pointer);
}
CLPLIBAPI void * CLPLINKAGE 
Clp_getUserPointer (Clp_Simplex * model)
{
  return model->model_->getUserPointer();
}
/* Pass in Callback function */
CLPLIBAPI void CLPLINKAGE 
Clp_registerCallBack(Clp_Simplex * model, 
		     clp_callback userCallBack)
{
  // Will be copy of users one
  delete model->handler_;
  model->handler_ = new CMessageHandler(*(model->model_->messageHandler()));
  model->handler_->setCallBack(userCallBack);
  model->handler_->setModel(model);
  model->model_->passInMessageHandler(model->handler_);
}
/* Unset Callback function */
CLPLIBAPI void CLPLINKAGE 
Clp_clearCallBack(Clp_Simplex * model)
{
  delete model->handler_;
  model->handler_=NULL;
}
/* Amount of print out:
   0 - none
   1 - just final
   2 - just factorizations
   3 - as 2 plus a bit more
   4 - verbose
   above that 8,16,32 etc just for selective debug
*/
CLPLIBAPI void CLPLINKAGE 
Clp_setLogLevel(Clp_Simplex * model, int value)
{
  model->model_->setLogLevel(value);
}
CLPLIBAPI int CLPLINKAGE 
Clp_logLevel(Clp_Simplex * model)
{
  return model->model_->logLevel();
}
/* length of names (0 means no names0 */
CLPLIBAPI int CLPLINKAGE 
Clp_lengthNames(Clp_Simplex * model)
{
  return model->model_->lengthNames();
}
/* Fill in array (at least lengthNames+1 long) with a row name */
CLPLIBAPI void CLPLINKAGE 
Clp_rowName(Clp_Simplex * model, int iRow, char * name)
{
  std::string rowName=model->model_->rowName(iRow);
  strcpy(name,rowName.c_str());
}
/* Fill in array (at least lengthNames+1 long) with a column name */
CLPLIBAPI void CLPLINKAGE 
Clp_columnName(Clp_Simplex * model, int iColumn, char * name)
{
  std::string columnName= model->model_->columnName(iColumn);
  strcpy(name,columnName.c_str());
}

/* General solve algorithm which can do presolve.
   See  ClpSolve.hpp for options
*/
CLPLIBAPI int CLPLINKAGE 
Clp_initialSolve(Clp_Simplex * model)
{
  return model->model_->initialSolve();
}
/* Dual initial solve */
CLPLIBAPI int CLPLINKAGE 
Clp_initialDualSolve(Clp_Simplex * model)
{
  return model->model_->initialDualSolve();
}
/* Primal initial solve */
CLPLIBAPI int CLPLINKAGE 
Clp_initialPrimalSolve(Clp_Simplex * model)
{
  return model->model_->initialPrimalSolve();
}
/* Dual algorithm - see ClpSimplexDual.hpp for method */
CLPLIBAPI int CLPLINKAGE 
Clp_dual(Clp_Simplex * model, int ifValuesPass)
{
  return model->model_->dual(ifValuesPass);
}
/* Primal algorithm - see ClpSimplexPrimal.hpp for method */
CLPLIBAPI int CLPLINKAGE 
Clp_primal(Clp_Simplex * model, int ifValuesPass)
{
  return model->model_->primal(ifValuesPass);
}
/* Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) */
CLPLIBAPI void CLPLINKAGE 
Clp_scaling(Clp_Simplex * model, int mode)
{
  model->model_->scaling(mode);
}
/* Gets scalingFlag */
CLPLIBAPI int CLPLINKAGE 
Clp_scalingFlag(Clp_Simplex * model)
{
  return model->model_->scalingFlag();
}
/* Crash - at present just aimed at dual, returns
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
CLPLIBAPI int CLPLINKAGE 
Clp_crash(Clp_Simplex * model, double gap,int pivot)
{
  return model->model_->crash(gap,pivot);
}
/* If problem is primal feasible */
CLPLIBAPI int CLPLINKAGE 
Clp_primalFeasible(Clp_Simplex * model)
{
  return model->model_->primalFeasible() ? 1 : 0;
}
/* If problem is dual feasible */
CLPLIBAPI int CLPLINKAGE 
Clp_dualFeasible(Clp_Simplex * model)
{
  return model->model_->dualFeasible() ? 1 : 0;
}
/* Dual bound */
CLPLIBAPI double CLPLINKAGE 
Clp_dualBound(Clp_Simplex * model)
{
  return model->model_->dualBound();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setDualBound(Clp_Simplex * model, double value)
{
  model->model_->setDualBound(value);
}
/* Infeasibility cost */
CLPLIBAPI double CLPLINKAGE 
Clp_infeasibilityCost(Clp_Simplex * model)
{
  return model->model_->infeasibilityCost();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setInfeasibilityCost(Clp_Simplex * model, double value)
{
  model->model_->setInfeasibilityCost(value);
}
/* Perturbation:
   50  - switch on perturbation
   100 - auto perturb if takes too long (1.0e-6 largest nonzero)
   101 - we are perturbed
   102 - don't try perturbing again
   default is 100
   others are for playing
*/
CLPLIBAPI int CLPLINKAGE 
Clp_perturbation(Clp_Simplex * model)
{
  return model->model_->perturbation();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setPerturbation(Clp_Simplex * model, int value)
{
  model->model_->setPerturbation(value);
}
/* Current (or last) algorithm */
CLPLIBAPI int CLPLINKAGE 
Clp_algorithm(Clp_Simplex * model)
{
  return model->model_->algorithm();
}
/* Set algorithm */
CLPLIBAPI void CLPLINKAGE 
Clp_setAlgorithm(Clp_Simplex * model, int value)
{
  model->model_->setAlgorithm(value);
}
/* Sum of dual infeasibilities */
CLPLIBAPI double CLPLINKAGE 
Clp_sumDualInfeasibilities(Clp_Simplex * model)
{
  return model->model_->sumDualInfeasibilities();
}
/* Number of dual infeasibilities */
CLPLIBAPI int CLPLINKAGE 
Clp_numberDualInfeasibilities(Clp_Simplex * model)
{
  return model->model_->numberDualInfeasibilities();
}
/* Sum of primal infeasibilities */
CLPLIBAPI double CLPLINKAGE 
Clp_sumPrimalInfeasibilities(Clp_Simplex * model)
{
  return model->model_->sumPrimalInfeasibilities();
}
/* Number of primal infeasibilities */
CLPLIBAPI int CLPLINKAGE 
Clp_numberPrimalInfeasibilities(Clp_Simplex * model)
{
  return model->model_->numberPrimalInfeasibilities();
}
/* Save model to file, returns 0 if success.  This is designed for
   use outside algorithms so does not save iterating arrays etc.
   It does not save any messaging information. 
   Does not save scaling values.
   It does not know about all types of virtual functions.
*/
CLPLIBAPI int CLPLINKAGE 
Clp_saveModel(Clp_Simplex * model, const char * fileName)
{
  return model->model_->saveModel(fileName);
}
/* Restore model from file, returns 0 if success,
   deletes current model */
CLPLIBAPI int CLPLINKAGE 
Clp_restoreModel(Clp_Simplex * model, const char * fileName)
{
  return model->model_->restoreModel(fileName);
}
  
/* Just check solution (for external use) - sets sum of
   infeasibilities etc */
CLPLIBAPI void CLPLINKAGE 
Clp_checkSolution(Clp_Simplex * model)
{
  model->model_->checkSolution();
}
/* Number of rows */
CLPLIBAPI int CLPLINKAGE 
Clp_getNumRows(Clp_Simplex * model)
{
  return model->model_->getNumRows();
}
/* Number of columns */
CLPLIBAPI int CLPLINKAGE 
Clp_getNumCols(Clp_Simplex * model)
{
  return model->model_->getNumCols();
}
/* Number of iterations */
CLPLIBAPI int CLPLINKAGE 
Clp_getIterationCount(Clp_Simplex * model)
{
  return model->model_->getIterationCount();
}
/* Are there a numerical difficulties? */
CLPLIBAPI int CLPLINKAGE 
Clp_isAbandoned(Clp_Simplex * model)
{
  return model->model_->isAbandoned() ? 1 : 0;
}
/* Is optimality proven? */
CLPLIBAPI int CLPLINKAGE 
Clp_isProvenOptimal(Clp_Simplex * model)
{
  return model->model_->isProvenOptimal() ? 1 : 0;
}
/* Is primal infeasiblity proven? */
CLPLIBAPI int CLPLINKAGE 
Clp_isProvenPrimalInfeasible(Clp_Simplex * model)
{
  return model->model_->isProvenPrimalInfeasible() ? 1 : 0;
}
/* Is dual infeasiblity proven? */
CLPLIBAPI int CLPLINKAGE 
Clp_isProvenDualInfeasible(Clp_Simplex * model)
{
  return model->model_->isProvenDualInfeasible() ? 1 : 0;
}
/* Is the given primal objective limit reached? */
CLPLIBAPI int CLPLINKAGE 
Clp_isPrimalObjectiveLimitReached(Clp_Simplex * model) 
{
  return model->model_->isPrimalObjectiveLimitReached() ? 1 : 0;
}
/* Is the given dual objective limit reached? */
CLPLIBAPI int CLPLINKAGE 
Clp_isDualObjectiveLimitReached(Clp_Simplex * model) 
{
  return model->model_->isDualObjectiveLimitReached() ? 1 : 0;
}
/* Iteration limit reached? */
CLPLIBAPI int CLPLINKAGE 
Clp_isIterationLimitReached(Clp_Simplex * model)
{
  return model->model_->isIterationLimitReached() ? 1 : 0;
}
/* Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore */
CLPLIBAPI double CLPLINKAGE 
Clp_getObjSense(Clp_Simplex * model)
{
  return model->model_->getObjSense();
}
/* Primal row solution */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getRowActivity(Clp_Simplex * model)
{
  return model->model_->getRowActivity();
}
/* Primal column solution */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getColSolution(Clp_Simplex * model)
{
  return model->model_->getColSolution();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setColSolution(Clp_Simplex * model, const double * input)
{
  model->model_->setColSolution(input);
}
/* Dual row solution */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getRowPrice(Clp_Simplex * model)
{
  return model->model_->getRowPrice();
}
/* Reduced costs */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getReducedCost(Clp_Simplex * model)
{
  return model->model_->getReducedCost();
}
/* Row lower */
CLPLIBAPI const double* CLPLINKAGE 
Clp_getRowLower(Clp_Simplex * model)
{
  return model->model_->getRowLower();
}
/* Row upper  */
CLPLIBAPI const double* CLPLINKAGE 
Clp_getRowUpper(Clp_Simplex * model)
{
  return model->model_->getRowUpper();
}
/* Objective */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getObjCoefficients(Clp_Simplex * model)
{
  return model->model_->getObjCoefficients();
}
/* Column Lower */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getColLower(Clp_Simplex * model)
{
  return model->model_->getColLower();
}
/* Column Upper */
CLPLIBAPI const double * CLPLINKAGE 
Clp_getColUpper(Clp_Simplex * model)
{
  return model->model_->getColUpper();
}
/* Objective value */
CLPLIBAPI double CLPLINKAGE 
Clp_getObjValue(Clp_Simplex * model)
{
  return model->model_->getObjValue();
}
/* Get variable basis info */
CLPLIBAPI const int CLPLINKAGE
Clp_getColumnStatus(Clp_Simplex * model,int sequence)
{
  return (int) model->model_->getColumnStatus(sequence);
}
/* Get row basis info */
CLPLIBAPI const int CLPLINKAGE
Clp_getRowStatus(Clp_Simplex * model,int sequence)
{
  return (int) model->model_->getRowStatus(sequence);
}
/* Small element value - elements less than this set to zero,
   default is 1.0e-20 */
CLPLIBAPI double CLPLINKAGE 
Clp_getSmallElementValue(Clp_Simplex * model)
{
  return model->model_->getSmallElementValue();
}
CLPLIBAPI void CLPLINKAGE 
Clp_setSmallElementValue(Clp_Simplex * model,double value)
{
  model->model_->setSmallElementValue(value);
}

