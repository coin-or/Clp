// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyUfl_H
#define ClpCholeskyUfl_H
#include "ClpCholeskyBase.hpp"
#ifdef __cplusplus
extern "C"{
#endif
#ifndef CLP_USE_CHOLMOD
#include "amd.h"
#else
#include "cholmod.h"
#endif  
#ifdef __cplusplus
          }
#endif
class ClpMatrixBase;
class ClpCholeskyDense;

/** Ufl class for Clp Cholesky factorization

If  you wish to use AMD code from University of Florida see

    http://www.cise.ufl.edu/research/sparse/amd

for terms of use

If  you wish to use CHOLMOD code from University of Florida see

    http://www.cise.ufl.edu/research/sparse/cholmod

for terms of use

*/
class ClpCholeskyUfl : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
#ifdef CLP_USE_CHOLMOD
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
#endif
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
#ifdef CLP_USE_CHOLMOD
    cholmod_factor * L_ ;
    cholmod_common c_ ;
#endif
};

#endif
