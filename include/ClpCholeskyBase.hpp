// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyBase_H
#define ClpCholeskyBase_H

#include "CoinPragma.hpp"
#include "CoinFinite.hpp"

class ClpInterior;

/** Base class for Clp Cholesky factorization
    Will do better factorization.  very crude ordering

    Derived classes will be using sophisticated methods apart from ClpCholeskyDense
*/

class ClpCholeskyDense;
class ClpMatrixBase;
class ClpCholeskyBase  {
  
public:
   /**@name Virtual methods that the derived classes may provide  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   returns non-zero if not enough memory.
   You can use preOrder to set up ADAT
   If using default symbolic etc then must set sizeFactor_ to
   size of input matrix to order (and to symbolic).
   Also just permute_ and permuteInverse_ should be created */
  virtual int order(ClpInterior * model);
  /** Does Symbolic factorization given permutation.
      This is called immediately after order.  If user provides this then
      user must provide factorize and solve.  Otherwise the default factorization is used
      returns non-zero if not enough memory */
  virtual int symbolic();
  /** Factorize - filling in rowsDropped and returning number dropped.
      If return code negative then out of memory */
  virtual int factorize(const double * diagonal, int * rowsDropped) ;
  /** Uses factorization to solve. */
  virtual void solve (double * region) ;
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
  /// Return size
  inline CoinBigIndex size() const
  { return sizeFactor_;};
  /// Return sparseFactor
  inline double * sparseFactor() const
  { return sparseFactor_;};
  /// Return diagonal
  inline double * diagonal() const
  { return diagonal_;};
  /// Return workDouble
  inline double * workDouble() const
  { return workDouble_;};
  /// If KKT on
  inline bool kkt() const
  { return doKKT_;};
  /// Set KKT
  inline void setKKT(bool yesNo)
  { doKKT_ = yesNo;};
   //@}
  
  
public:

   /**@name Constructors, destructor
    */
   //@{
  /** Constructor which has dense columns activated.
      Default is off. */
  ClpCholeskyBase(int denseThreshold=-1);
   /** Destructor (has to be public) */
   virtual ~ClpCholeskyBase();
  // Copy
   ClpCholeskyBase(const ClpCholeskyBase&);
  // Assignment
   ClpCholeskyBase& operator=(const ClpCholeskyBase&);
   //@}
  //@{
  ///@name Other
  /// Clone
  virtual ClpCholeskyBase * clone() const;
 
  /// Returns type
  inline int type() const
  { if (doKKT_) return 100; else return type_;};
protected:
  /// Sets type
  void setType(int type) {type_=type;};
   //@}
   
  /**@name Symbolic, factor and solve */
  //@{
  /** Symbolic1  - works out size without clever stuff.
      Uses upper triangular as much easier.
      Returns size
   */
  int symbolic1(const CoinBigIndex * Astart, const int * Arow);
  /** Symbolic2  - Fills in indices
      Uses lower triangular so can do cliques etc
   */
  void symbolic2(const CoinBigIndex * Astart, const int * Arow);
  /** Factorize - filling in rowsDropped and returning number dropped
      in integerParam.
   */
  void factorizePart2(int * rowsDropped) ;
  /** solve - 1 just first half, 2 just second half - 3 both.
  If 1 and 2 then diagonal has sqrt of inverse otherwise inverse
  */
  void solve(double * region, int type);
  /// Forms ADAT - returns nonzero if not enough memory
  int preOrder(bool lowerTriangular, bool includeDiagonal, bool doKKT);
  //@}
    
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// type (may be useful) if > 20 do KKT
  int type_;
  /// Doing full KKT (only used if default symbolic and factorization)
  bool doKKT_;
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
  /// permute inverse.
  int * permuteInverse_;
  /// main permute.
  int * permute_;
  /// numberRowsDropped.  Number of rows gone
  int numberRowsDropped_;
  /// sparseFactor.
  double * sparseFactor_;
  /// choleskyStart - element starts
  CoinBigIndex * choleskyStart_;
  /// choleskyRow (can be shorter than sparsefactor)
  int * choleskyRow_;
  /// Index starts
  CoinBigIndex * indexStart_;
  /// Diagonal
  double * diagonal_;
  /// double work array
  double * workDouble_;
  /// link array
  int * link_;
  // Integer work array
  CoinBigIndex * workInteger_;
  // Clique information
  int * clique_;
  /// sizeFactor.
  CoinBigIndex sizeFactor_;
  /// Size of index array
  CoinBigIndex sizeIndex_;
  /// First dense row
  int firstDense_;
  /// integerParameters
  int integerParameters_[64];
  /// doubleParameters;
  double doubleParameters_[64];
  /// Row copy of matrix
  ClpMatrixBase * rowCopy_;
  /// Dense indicators
  char * whichDense_;
  /// Dense columns (updated)
  double * denseColumn_;
  /// Dense cholesky
  ClpCholeskyDense * dense_;
  /// Dense threshold
  int denseThreshold_;
  //@}
};

#endif
