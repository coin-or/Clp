/* $Id: ClpCholeskyMumps.hpp 1370 2009-06-04 09:37:13Z forrest $ */
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyMumps_H
#define ClpCholeskyMumps_H
#include "ClpCholeskyBase.hpp"
#ifdef __cplusplus
extern "C"{
#endif
#include "amd.h"
#include "dmumps_c.h"
  //#include "mpi.h"
  //#include "/home/jjforre/cbc-trunk/ThirdParty/Mumps/MUMPS/libseq/mpi.h"
#ifdef __cplusplus
          }
#endif
class ClpMatrixBase;
class ClpCholeskyDense;

/** Mumps class for Clp Cholesky factorization

*/
class ClpCholeskyMumps : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
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
  //@}


  /**@name Constructors, destructor */
  //@{
  /** Constructor which has dense columns activated.
      Default is off. */
  ClpCholeskyMumps(int denseThreshold=-1);
  /** Destructor  */
  virtual ~ClpCholeskyMumps();
  // Copy
  ClpCholeskyMumps(const ClpCholeskyMumps&);
  // Assignment
  ClpCholeskyMumps& operator=(const ClpCholeskyMumps&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
  // Mumps structure
  DMUMPS_STRUC_C mumps_;
};

#endif
