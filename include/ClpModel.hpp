// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpModel_H
#define ClpModel_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>

#include "ClpMatrixBase.hpp"
#include "CoinMessageHandler.hpp"
#include "ClpParameters.hpp"
#include "ClpObjective.hpp"
class ClpEventHandler;

// Plus infinity
#ifndef COIN_DBL_MAX
#define COIN_DBL_MAX DBL_MAX
#endif

/** This is the base class for Linear and quadratic Models
    This knows nothing about the algorithm, but it seems to
    have a reasonable amount of information

    I would welcome suggestions for what should be in this and
    how it relates to OsiSolverInterface.  Some methods look
    very similar.

*/

class ClpModel {

public:

  /**@name Constructors and destructor 
     Note - copy methods copy ALL data so can chew up memory
     until other copy is freed
   */
  //@{
    /// Default constructor
    ClpModel (  );

    /// Copy constructor. 
    ClpModel(const ClpModel &);
    /// Assignment operator. This copies the data
    ClpModel & operator=(const ClpModel & rhs);
  /** Subproblem constructor.  A subset of whole model is created from the 
      row and column lists given.  The new order is given by list order and
      duplicates are allowed.  Name and integer information can be dropped
  */
    ClpModel (const ClpModel * wholeModel,
      int numberRows, const int * whichRows,
      int numberColumns, const int * whichColumns,
	      bool dropNames=true, bool dropIntegers=true);
    /// Destructor
    ~ClpModel (  );
  //@}

  /**@name Load model - loads some stuff and initializes others */
  //@{
    /** Loads a problem (the constraints on the
        rows are given by lower and upper bounds). If a pointer is 0 then the
        following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>
    */
  void loadProblem (  const ClpMatrixBase& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  void loadProblem (  const CoinPackedMatrix& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);

  /** Just like the other loadProblem() method except that the matrix is
	given in a standard column major ordered format (without gaps). */
  void loadProblem (  const int numcols, const int numrows,
		     const CoinBigIndex* start, const int* index,
		     const double* value,
		     const double* collb, const double* colub,   
		     const double* obj,
		      const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  /// This one is for after presolve to save memory
  void loadProblem (  const int numcols, const int numrows,
		     const CoinBigIndex* start, const int* index,
		      const double* value,const int * length,
		     const double* collb, const double* colub,   
		     const double* obj,
		      const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  /** Load up quadratic objective.  This is stored as a CoinPackedMatrix */
  void loadQuadraticObjective(const int numberColumns, 
			      const CoinBigIndex * start,
			      const int * column, const double * element);
  void loadQuadraticObjective (  const CoinPackedMatrix& matrix);
  /// Get rid of quadratic objective
  void deleteQuadraticObjective();
  /// This just loads up a row objective
  void setRowObjective(const double * rowObjective);
  /// Read an mps file from the given filename
  int readMps(const char *filename,
	      bool keepNames=false,
	      bool ignoreErrors = false);
  /// Copy in integer informations
  void copyInIntegerInformation(const char * information);
  /// Drop integer informations
  void deleteIntegerInformation();
  /// Resizes rim part of model 
  void resize (int newNumberRows, int newNumberColumns);
  /// Deletes rows
  void deleteRows(int number, const int * which);
  /// Add rows
  void addRows(int number, const double * rowLower, 
	       const double * rowUpper,
	       const int * rowStarts, const int * columns,
	       const double * elements);
  /// Add rows
  void addRows(int number, const double * rowLower, 
	       const double * rowUpper,
	       const int * rowStarts, const int * rowLengths,
	       const int * columns,
	       const double * elements);
  void addRows(int number, const double * rowLower, 
	       const double * rowUpper,
	       const CoinPackedVectorBase * const * rows);

  /// Deletes columns
  void deleteColumns(int number, const int * which);
  /// Add columns
  void addColumns(int number, const double * columnLower, 
		  const double * columnUpper,
		  const double * objective,
		  const int * columnStarts, const int * rows,
		  const double * elements);
  void addColumns(int number, const double * columnLower, 
		  const double * columnUpper,
		  const double * objective,
		  const int * columnStarts, const int * columnLengths,
		  const int * rows,
		  const double * elements);
  void addColumns(int number, const double * columnLower, 
	       const double * columnUpper,
		  const double * objective,
	       const CoinPackedVectorBase * const * columns);
  /** Borrow model.  This is so we don't have to copy large amounts
      of data around.  It assumes a derived class wants to overwrite
      an empty model with a real one - while it does an algorithm */
  void borrowModel(ClpModel & otherModel);
  /** Return model - nulls all arrays so can be deleted safely
      also updates any scalars */
  void returnModel(ClpModel & otherModel);

  /// Create empty ClpPackedMatrix
  void createEmptyMatrix();
  /// Drops names - makes lengthnames 0 and names empty
  void dropNames();
  /// Copies in names
  void copyNames(std::vector<std::string> & rowNames,
		 std::vector<std::string> & columnNames);
  
    /** Write the problem in MPS format to the specified file.

	Row and column names may be null.
	formatType is
	<ul>
	  <li> 0 - normal
	  <li> 1 - extra accuracy 
	  <li> 2 - IEEE hex (later)
	</ul>

	Returns non-zero on I/O error
    */
    int writeMps(const char *filename, 
		  int formatType=0,int numberAcross=2,
		 double objSense=0.0) const ;
  //@}
  /**@name gets and sets */
  //@{ 
   /// Number of rows
   inline int numberRows() const {
      return numberRows_;
   }
   inline int getNumRows() const {
      return numberRows_;
   }
   /// Number of columns
   inline int getNumCols() const {
      return numberColumns_;
   }
   inline int numberColumns() const {
      return numberColumns_;
   }
   /// Primal tolerance to use
   inline double primalTolerance() const {
      return dblParam_[ClpPrimalTolerance];
   }
   void setPrimalTolerance( double value) ;
   /// Dual tolerance to use
   inline double dualTolerance() const  { return dblParam_[ClpDualTolerance]; }
   void setDualTolerance( double value) ;
  /// Dual objective limit
  inline double dualObjectiveLimit() const { return dblParam_[ClpDualObjectiveLimit];}
  void setDualObjectiveLimit(double value);
  /// Objective offset
  inline double objectiveOffset() const { return dblParam_[ClpObjOffset];}
  void setObjectiveOffset(double value);
  inline std::string problemName() const { return strParam_[ClpProbName]; };
   /// Number of iterations
   inline int numberIterations() const  { return numberIterations_; }
   inline int getIterationCount() const { return numberIterations_; }
  inline void setNumberIterations(int numberIterations)
  { numberIterations_ = numberIterations;};
  /** Solve type - 1 simplex, 2 simplex interface, 3 Interior.*/
  inline int solveType() const
  { return solveType_;};
  inline void setSolveType(int type)
  { solveType_=type;};
   /// Maximum number of iterations
   inline int maximumIterations() const { return intParam_[ClpMaxNumIteration]; }
   void setMaximumIterations(int value);
  /// Maximum time in seconds (from when set called)
   inline double maximumSeconds() const { return dblParam_[ClpMaxSeconds]; }
   void setMaximumSeconds(double value);
  /// Returns true if hit maximum iterations (or time)
  bool hitMaximumIterations() const;
   /** Status of problem:
       0 - optimal
       1 - primal infeasible
       2 - dual infeasible
       3 - stopped on iterations or time
       4 - stopped due to errors
       5 - stopped by event handler (virtual int ClpEventHandler::event())
   */
   inline int status() const            { return problemStatus_; }
  /// Set problem status
  inline void setProblemStatus(int problemStatus)
  { problemStatus_ = problemStatus;};
   /** Secondary status of problem - may get extended
       0 - none
       1 - primal infeasible because dual limit reached
       2 - scaled problem optimal - unscaled problem has primal infeasibilities
       3 - scaled problem optimal - unscaled problem has dual infeasibilities
       4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
       100 up - translation of enum from ClpEventHandler
   */
   inline int secondaryStatus() const            { return secondaryStatus_; }
  inline void setSecondaryStatus(int status)
  { secondaryStatus_ = status;};
   /// Are there a numerical difficulties?
   bool isAbandoned() const             { return problemStatus_==4; }
   /// Is optimality proven?
   bool isProvenOptimal() const         { return problemStatus_==0; }
   /// Is primal infeasiblity proven?
   bool isProvenPrimalInfeasible() const { return problemStatus_==1; }
   /// Is dual infeasiblity proven?
   bool isProvenDualInfeasible() const  { return problemStatus_==2; }
   /// Is the given primal objective limit reached?
   bool isPrimalObjectiveLimitReached() const ;
   /// Is the given dual objective limit reached?
   bool isDualObjectiveLimitReached() const ;
   /// Iteration limit reached?
   bool isIterationLimitReached() const { return problemStatus_==3; }
   /// Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore
   inline double optimizationDirection() const {
      return  optimizationDirection_;
   }
   inline double getObjSense() const    { return optimizationDirection_; }
   void setOptimizationDirection(double value);
   /// Primal row solution
   inline double * primalRowSolution() const    { return rowActivity_; }
   inline const double * getRowActivity() const { return rowActivity_; }
   /// Primal column solution
   inline double * primalColumnSolution() const { return columnActivity_; }
   inline const double * getColSolution() const { return columnActivity_; }
   inline void setColSolution(const double * input)
   { memcpy(columnActivity_,input,numberColumns_*sizeof(double));};
   /// Dual row solution
   inline double * dualRowSolution() const      { return dual_; }
   inline const double * getRowPrice() const    { return dual_; }
   /// Reduced costs
   inline double * dualColumnSolution() const   { return reducedCost_; }
   inline const double * getReducedCost() const { return reducedCost_; }
   /// Row lower
   inline double* rowLower() const              { return rowLower_; }
   inline const double* getRowLower() const     { return rowLower_; }
   /// Row upper 
   inline double* rowUpper() const              { return rowUpper_; }
   inline const double* getRowUpper() const     { return rowUpper_; }
   /// Objective
   inline double * objective() const            
  {
    if (objective_) {
      double offset; 
      return objective_->gradient(NULL,offset,false);
    } else {
      return NULL;
    }
  }
   inline double * objective(const double * solution, double & offset,bool refresh=true) const            
  {
    offset=0.0;
    if (objective_) {
      return objective_->gradient(solution,offset,refresh);
    } else {
      return NULL;
    }
  }
   inline const double * getObjCoefficients() const 
  { 
    if (objective_) {
      double offset; 
      return objective_->gradient(NULL,offset,false);
    } else {
      return NULL;
    }
  }
   /// Row Objective
   inline double * rowObjective() const         { return rowObjective_; }
   inline const double * getRowObjCoefficients() const {
      return rowObjective_;
   }
   /// Column Lower
   inline double * columnLower() const          { return columnLower_; }
   inline const double * getColLower() const    { return columnLower_; }
   /// Column Upper
   inline double * columnUpper() const          { return columnUpper_; }
   inline const double * getColUpper() const    { return columnUpper_; }
   /// Matrix (if not ClpPackedmatrix be careful about memory leak
   inline CoinPackedMatrix * matrix() const {
     if ( matrix_ == NULL ) return NULL;
     else return matrix_->getPackedMatrix();
   }
   /// Number of elements in matrix
   inline int getNumElements() const 
     { return matrix_->getNumElements();};
   /** Small element value - elements less than this set to zero,
      default is 1.0e-20 */
   inline double getSmallElementValue() const
  { return smallElement_;}; 
  inline void setSmallElementValue(double value)
  { smallElement_=value;}; 
   /// Row Matrix 
   inline ClpMatrixBase * rowCopy() const       { return rowCopy_; }
   /// Clp Matrix 
   inline ClpMatrixBase * clpMatrix() const     { return matrix_; }
  /** Replace Clp Matrix (current is not deleted and new is used)
      So up to user to delete one
  */
   void replaceMatrix(ClpMatrixBase * matrix);
   /// Objective value
   inline double objectiveValue() const {
      return objectiveValue_*optimizationDirection_ - dblParam_[ClpObjOffset];
   }
   inline double getObjValue() const {
      return objectiveValue_*optimizationDirection_ - dblParam_[ClpObjOffset];
   }
   /// Integer information
   inline char * integerInformation() const     { return integerType_; }
   /** Infeasibility/unbounded ray (NULL returned if none/wrong)
       Up to user to use delete [] on these arrays.  */
   double * infeasibilityRay() const;
   double * unboundedRay() const;
  /// See if status array exists (partly for OsiClp)
  inline bool statusExists() const
  { return (status_!=NULL);};
  /// Return address of status array (char[numberRows+numberColumns])
  inline unsigned char *  statusArray() const
  { return status_;};
  /** Return copy of status array (char[numberRows+numberColumns]),
      use delete [] */
  unsigned char *  statusCopy() const;
  /// Copy in status vector
  void copyinStatus(const unsigned char * statusArray);

  /// User pointer for whatever reason
  inline void setUserPointer (void * pointer)
  { userPointer_=pointer;};
  inline void * getUserPointer () const
  { return userPointer_;};
  //@}
  /**@name Message handling */
  //@{
   /// Pass in Message handler (not deleted at end)
   void passInMessageHandler(CoinMessageHandler * handler);
  /// Pass in Message handler (not deleted at end) and return current
  CoinMessageHandler * pushMessageHandler(CoinMessageHandler * handler,
					  bool & oldDefault);
  /// back to previous message handler
  void popMessageHandler(CoinMessageHandler * oldHandler,bool oldDefault);
   /// Set language
   void newLanguage(CoinMessages::Language language);
   void setLanguage(CoinMessages::Language language) { newLanguage(language); }
   /// Return handler
   CoinMessageHandler * messageHandler() const       { return handler_; }
   /// Return messages
   CoinMessages messages() const                     { return messages_; }
   /// Return pointer to messages
   CoinMessages * messagesPointer()                  { return & messages_; }
  /** Amount of print out:
      0 - none
      1 - just final
      2 - just factorizations
      3 - as 2 plus a bit more
      4 - verbose
      above that 8,16,32 etc just for selective debug
  */
   void setLogLevel(int value)    { handler_->setLogLevel(value); }
   int logLevel() const           { return handler_->logLevel(); }
   /// Pass in Event handler (cloned and deleted at end)
   void passInEventHandler(const ClpEventHandler * eventHandler);
   /// length of names (0 means no names0
   inline int lengthNames() const { return lengthNames_; }
   /// Row names
   const std::vector<std::string> * rowNames() const {
      return &rowNames_;
   }
   const std::string& rowName(int iRow) const {
      return rowNames_[iRow];
   }
   /// Column names
   const std::vector<std::string> * columnNames() const {
      return &columnNames_;
   }
   const std::string& columnName(int iColumn) const {
      return columnNames_[iColumn];
   }
  /// Objective methods
  inline ClpObjective * objectiveAsObject() const
  { return objective_;};
  void setObjective(ClpObjective * objective);
  void setObjectivePointer(ClpObjective * objective)
  { objective_ = objective;};
  /** Solve a problem with no elements - return status and
      dual and primal infeasibilites */
  int emptyProblem(int * infeasNumber=NULL, double * infeasSum=NULL,bool printMessage=true);
  
  //@}


  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the volume algorithm), the parameter is not yet implemented
     for the solver or simply the value of the parameter is out of the range
     the solver accepts. If a parameter setting call returns false check the
     details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.

     ** once it has been decided where solver sits this may be redone
  */
  //@{
    /// Set an integer parameter
    bool setIntParam(ClpIntParam key, int value) ;
    /// Set an double parameter
    bool setDblParam(ClpDblParam key, double value) ;
    /// Set an string parameter
    bool setStrParam(ClpStrParam key, const std::string & value);
    // Get an integer parameter
    inline bool getIntParam(ClpIntParam key, int& value) const {
      if (key!=ClpLastIntParam) {
	value = intParam_[key];
	return true;
      } else {
	return false;
      }
    }
    // Get an double parameter
    inline bool getDblParam(ClpDblParam key, double& value) const {
      if (key!=ClpLastDblParam) {
	value = dblParam_[key];
	return true;
      } else {
	return false;
      }
    }
    // Get a string parameter
    inline bool getStrParam(ClpStrParam key, std::string& value) const {
      if (key!=ClpLastStrParam) {
	value = strParam_[key];
	return true;
      } else {
	return false;
      }
    }
  //@}

  /**@name private or protected methods */
  //@{
protected:
  /// Does most of deletion
  void gutsOfDelete();
  /** Does most of copying
      If trueCopy false then just points to arrays */
  void gutsOfCopy(const ClpModel & rhs, bool trueCopy=true);
  /// gets lower and upper bounds on rows
  void getRowBound(int iRow, double& lower, double& upper) const;
  /// puts in format I like - 4 array matrix - may make row copy 
  void gutsOfLoadModel ( int numberRows, int numberColumns,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  //@}


////////////////// data //////////////////
protected:

  /**@name data */
  //@{
  /// Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore
  double optimizationDirection_;
  /// Array of double parameters
  double dblParam_[ClpLastDblParam];
  /// Objective value
  double objectiveValue_;
  /// Small element value
  double smallElement_;
  /// Number of rows
  int numberRows_;
  /// Number of columns
  int numberColumns_;
  /// Row activities
  double * rowActivity_;
  /// Column activities
  double * columnActivity_;
  /// Duals
  double * dual_;
  /// Reduced costs
  double * reducedCost_;
  /// Row lower 
  double* rowLower_;
  /// Row upper 
  double* rowUpper_;
  /// Objective
  ClpObjective * objective_;
  /// Row Objective (? sign)  - may be NULL
  double * rowObjective_;
  /// Column Lower
  double * columnLower_;
  /// Column Upper
  double * columnUpper_;
  /// Packed matrix
  ClpMatrixBase * matrix_;
  /// Row copy if wanted
  ClpMatrixBase * rowCopy_;
  /// Infeasible/unbounded ray
  double * ray_;
  /** Status Region.  I know that not all algorithms need a status
      array, but it made sense for things like crossover and put
      all permanent stuff in one place.  No assumption is made
      about what is in status array (although it might be good to reserve
      bottom 3 bits (i.e. 0-7 numeric) for classic status).  This
      is number of columns + number of rows long (in that order).
  */
  unsigned char * status_;
  /// Integer information
  char * integerType_;
  /// User pointer for whatever reason
  void * userPointer_;
  /// Array of integer parameters
  int intParam_[ClpLastIntParam];
  /// Number of iterations
  int numberIterations_;
  /** Solve type - 1 simplex, 2 simplex interface, 3 Interior.*/
  int solveType_;
  /// Status of problem
  int problemStatus_;
  /// Secondary status of problem
  int secondaryStatus_;
  /// length of names (0 means no names)
  int lengthNames_;
  /// Message handler
  CoinMessageHandler * handler_;
  /// Flag to say if default handler (so delete)
  bool defaultHandler_;
  /// Event handler
  ClpEventHandler * eventHandler_;
  /// Row names
  std::vector<std::string> rowNames_;
  /// Column names
  std::vector<std::string> columnNames_;
  /// Messages
  CoinMessages messages_;
  /// Array of string parameters
  std::string strParam_[ClpLastStrParam];
  //@}
};
/** This is a tiny class where data can be saved round calls */
class ClpDataSave {

public:
  /**@name Constructors and destructor 
   */
  //@{
    /// Default constructor
    ClpDataSave (  );

    /// Copy constructor. 
    ClpDataSave(const ClpDataSave &);
    /// Assignment operator. This copies the data
    ClpDataSave & operator=(const ClpDataSave & rhs);
    /// Destructor
    ~ClpDataSave (  );

  //@}

////////////////// data //////////////////
public:

  /**@name data - with same names as in other classes*/
  //@{
  double dualBound_;
  double infeasibilityCost_;
  int sparseThreshold_;
  int perturbation_;

  //@}
};

#endif
