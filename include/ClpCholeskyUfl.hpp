// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyUfl_H
#define ClpCholeskyUfl_H

#include "ClpCholeskyBase.hpp"
#ifdef __cplusplus
extern "C"{
#endif
#include "amd.h"
#ifdef __cplusplus
          }
#endif
class ClpMatrixBase;
class ClpCholeskyDense;

/** Ufl class for Clp Cholesky factorization

If  you wish to use AMD code from University of Florida see

    http://www.cise.ufl.edu/research/sparse/amd

for terms of use

*/
class ClpCholeskyUfl : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
  //@}


  /**@name Constructors, destructor */
  //@{
  /** Constructor which has dense columns activated.
      Default is off. */
  ClpCholeskyUfl(int denseThreshold=-1);
  /** Destructor  */
  virtual ~ClpCholeskyUfl();
  // Copy
  ClpCholeskyUfl(const ClpCholeskyUfl&);
  // Assignment
  ClpCholeskyUfl& operator=(const ClpCholeskyUfl&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
};

#endif
