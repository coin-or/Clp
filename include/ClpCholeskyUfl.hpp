// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyUfl_H
#define ClpCholeskyUfl_H

#include "ClpCholeskyBase.hpp"
#ifdef __cplusplus
extern "C"{
#endif
#include "/home/forrest/ordering/AMD/Include/amd.h"
#ifdef __cplusplus
          }
#endif
class ClpMatrixBase;
class ClpCholeskyDense;

/** Ufl class for Clp Cholesky factorization


-------------------------------------------------------------------------------

AMD Version 1.1 (Jan. 21, 2004),  Copyright (c) 2004 by Timothy A.
Davis, Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.

AMD License:

    Your use or distribution of AMD or any modified version of
    AMD implies that you agree to this License.

    THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
    EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Permission is hereby granted to use or copy this program, provided
    that the Copyright, this License, and the Availability of the original
    version is retained on all copies.  User documentation of any code that
    uses AMD or any modified version of AMD code must cite the
    Copyright, this License, the Availability note, and "Used by permission."
    Permission to modify the code and to distribute modified code is granted,
    provided the Copyright, this License, and the Availability note are
    retained, and a notice that the code was modified is included.  This
    software was developed with support from the National Science Foundation,
    and is provided to you free of charge.

Availability:

    http://www.cise.ufl.edu/research/sparse/amd

-------------------------------------------------------------------------------

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
