// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyBase_H
#define ClpCholeskyBase_H

#include "CoinPragma.hpp"

class ClpInterior;

/** Abstract Base class for Clp Cholesky factorization

    Derived classes will be using sophisticated methods apart from ClpCholeskyDense
*/

class ClpCholeskyBase  {
  
public:
   /**@name Virtual methods that the derived classes must provide  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) = 0;
  /** Factorize - filling in rowsDropped and returning number dropped.
      If return code negative then out of memory */
  virtual int factorize(const double * diagonal, int * rowsDropped) =0;
  /** Uses factorization to solve. */
  virtual void solve (double * region) = 0;
  /** Uses factorization to solve. - given as if KKT.
   region1 is rows+columns, region2 is rows */
  virtual void solveKKT (double * region1, double * region2, const double * diagonal,
			 double diagonalScaleFactor);
  //@}

  /**@name Gets */
  //@{
  /// status.  Returns status
  inline int status() const 
  {return status_;};
  /// numberRowsDropped.  Number of rows gone
  inline int numberRowsDropped() const 
  {return numberRowsDropped_;};
  /// reset numberRowsDropped and rowsDropped.
  void resetRowsDropped();
  /// rowsDropped - which rows are gone
  inline char * rowsDropped() const 
  {return rowsDropped_;};
  /// choleskyCondition.
  inline double choleskyCondition() const 
  {return choleskyCondition_;};
  /// rank.  Returns rank
  inline int rank() const 
  {return numberRows_-numberRowsDropped_;};
  /// Return number of rows
  inline int numberRows() const 
  {return numberRows_;};
   //@}
  
  
protected:

   /**@name Constructors, destructor<br>
      <strong>NOTE</strong>: All constructors are protected. There's no need
      to expose them, after all, this is an abstract class. */
   //@{
   /** Default constructor. */
   ClpCholeskyBase();
   /** Destructor (has to be public) */
public:
   virtual ~ClpCholeskyBase();
protected:
  // Copy
   ClpCholeskyBase(const ClpCholeskyBase&);
  // Assignment
   ClpCholeskyBase& operator=(const ClpCholeskyBase&);
   //@}
  //@{
  ///@name Other
  /// Clone
public:
  virtual ClpCholeskyBase * clone() const = 0;
protected:
 
  /// Returns type
  inline int type() const
  { return type_;};
  /// Sets type
  void setType(int type) {type_=type;};
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
   /// type (may be useful)
   int type_;
  /// pivotTolerance.
  double pivotTolerance_;
  /// zeroTolerance.
  double zeroTolerance_;
  /// choleskyCondition.
  double choleskyCondition_;
  /// model.
  ClpInterior * model_;
  /// numberTrials.  Number of trials before rejection
  int numberTrials_;
  /// numberRows.  Number of Rows in factorization
  int numberRows_;
  /// status.  Status of factorization
  int status_;
  /// rowsDropped
  char * rowsDropped_;
  /// permuteIn.
  int * permuteIn_;
  /// permuteOut.
  int * permuteOut_;
  /// numberRowsDropped.  Number of rows gone
  int numberRowsDropped_;
  //@}
};

#endif
