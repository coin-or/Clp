// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyDense_H
#define ClpCholeskyDense_H

#include "ClpCholeskyBase.hpp"

typedef double longDouble;

/** Dense class for Clp Cholesky factorization

*/
class ClpMatrixBase;
class ClpCholeskyDense : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
  /** Factorize - filling in rowsDropped and returning number dropped.
      If return code negative then out of memory */
  virtual int factorize(const double * diagonal, int * rowsDropped) ;
  /** Uses factorization to solve. */
  virtual void solve (double * region) ;
  //@}

   /**@name Non virtual methods for ClpCholeskyDense  */
   //@{
  /** Reserves space..
   Returns non-zero if not enough memory */
   int reserveSpace(int numberRows) ;
  /** part 2 of Factorize - filling in rowsDropped and returning number dropped */
  int factorizePart2(int * rowsDropped) ;
  /// A
  inline double * aMatrix() const
  { return work_;}
  //@}


  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  ClpCholeskyDense();
  /** Destructor  */
  virtual ~ClpCholeskyDense();
  // Copy
  ClpCholeskyDense(const ClpCholeskyDense&);
  // Assignment
  ClpCholeskyDense& operator=(const ClpCholeskyDense&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
  /**@name Data members */
   //@{
  /// ADAT stored in full
  double * work_;
  /// Row copy of matrix
  ClpMatrixBase * rowCopy_;
  //@}
};

#endif
