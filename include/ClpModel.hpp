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

#define CLP_INFINITY 1e30

/** This is the base class for Linear Models
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
  /// Deletes columns
  
  void deleteColumns(int number, const int * which);
  /** Borrow model.  This is so we dont have to copy large amounts
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
   /// Number of iterations
   inline int numberIterations() const  { return numberIterations_; }
   inline int getIterationCount() const { return numberIterations_; }
  inline void setNumberIterations(int numberIterations)
  { numberIterations_ = numberIterations;};
  /// Solve type - 1 simplex, 2 simplex interface
  inline int solveType() const
  { return solveType_;};
  inline void setSolveType(int type)
  { solveType_=type;};
   /// Maximum number of iterations
   inline int maximumIterations() const { return intParam_[ClpMaxNumIteration]; }
   void setMaximumIterations(int value);
   /** Status of problem:
       0 - optimal
       1 - primal infeasible
       2 - dual infeasible
       3 - stopped on iterations etc
       4 - stopped due to errors
   */
   inline int status() const            { return problemStatus_; }
  /// Set problem status
  inline void setProblemStatus(int problemStatus)
  { problemStatus_ = problemStatus;};
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
   inline int optimizationDirection() const {
      return (int) optimizationDirection_;
   }
   inline double getObjSense() const    { return optimizationDirection_; }
   void setOptimizationDirection(int value);
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
  { double offset; return objective_->gradient(NULL,offset); }
   inline const double * getObjCoefficients() const 
  { double offset; return objective_->gradient(NULL,offset); }
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

  //@}
  /**@name Message handling */
  //@{
   /// Pass in Message handler (not deleted at end)
   void passInMessageHandler(CoinMessageHandler * handler);
   /// Set language
   void newLanguage(CoinMessages::Language language);
   void setLanguage(CoinMessages::Language language) { newLanguage(language); }
   /// Return handler
   CoinMessageHandler * messageHandler() const       { return handler_; }
   /// Return messages
   CoinMessages messages() const                     { return messages_; }
   /// Return pointer to messages
   CoinMessages * messagesPointer()                  { return & messages_; }
   /// Log level
   void setLogLevel(int value)    { handler_->setLogLevel(value); }
   int logLevel() const           { return handler_->logLevel(); }
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
  /// Array of integer parameters
  int intParam_[ClpLastIntParam];
  /// Number of iterations
  int numberIterations_;
  /// Solve type - 1 simplex, 2 simplex interface
  int solveType_;
  /// Status of problem
  int problemStatus_;
  /// length of names (0 means no names)
  int lengthNames_;
  /// Message handler
  CoinMessageHandler * handler_;
  /// Flag to say if default handler (so delete)
  bool defaultHandler_;
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
