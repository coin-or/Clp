// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyWssmp_H
#define ClpCholeskyWssmp_H

#include "ClpCholeskyBase.hpp"


/** Wssmp class for Clp Cholesky factorization

*/
class ClpMatrixBase;
class ClpCholeskyDense;
class ClpCholeskyWssmp : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
  /** Factorize - filling in rowsDropped and returning number dropped */
  virtual int factorize(const double * diagonal, int * rowsDropped) ;
  /** Uses factorization to solve. */
  virtual void solve (double * region) ;
  //@}


  /**@name Constructors, destructor */
  //@{
  /** Constructor which has dense columns activated.
      Default is off. */
  ClpCholeskyWssmp(int denseThreshold=-1);
  /** Destructor  */
  virtual ~ClpCholeskyWssmp();
  // Copy
  ClpCholeskyWssmp(const ClpCholeskyWssmp&);
  // Assignment
  ClpCholeskyWssmp& operator=(const ClpCholeskyWssmp&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
  /**@name Data members */
   //@{
  /// sparseFactor.
  double * sparseFactor_;
  /// choleskyStart
  CoinBigIndex * choleskyStart_;
  /// choleskyRow
  int * choleskyRow_;
  /// sizeFactor.
  CoinBigIndex sizeFactor_;
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
  /// Dense square cholesky
  ClpCholeskyDense * dense_;
  /// Dense threshold
  int denseThreshold_;
  //@}
};

#endif
