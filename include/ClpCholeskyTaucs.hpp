// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyTaucs_H
#define ClpCholeskyTaucs_H
#include "taucs.h"
#include "ClpCholeskyBase.hpp"


/** Taucs class for Clp Cholesky factorization

Taucs comes with no warranty whatsoever and is distributed under the GNU LGPL (Library or 
Lesser GNU Public Library). The license is available in www.gnu.org. Alternatively, you can also 
elect to use Taucs under the following UMFPACK-style license, which is simpler to understand 
than the LGPL: 
Taucs. Version 1.0, November 29, 2001. Copyright (c) 2001 by Sivan Toledo, Tel-Aviv
Univesity, stoledo@tau.ac.il. All Rights Reserved.
TAUCS License:
Your use or distribution of TAUCS or any derivative code implies that you agree to
this License OR to the GNU LGPL.
THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
Permission is hereby granted to use or copy this program, provided that the Copyright,
this License, and the Availability of the original version is retained on all copies. User
documentation of any code that uses this code or any derivative code must cite the
Copyright, this License, the Availability note, and"Used by permission."If this code or
any derivative code is accessible from within M....., then typing "help taucs" must
cite the Copyright, and "type taucs" must also cite this License and the Availability
note. Permission to modify the code and to distribute modified code is granted,
provided the Copyright, this License, and the Availability note are retained, and a
notice that the code was modified is included. This software is provided to you free
of charge.

The taucs.h file was modified to put 

#ifdef __cplusplus
extern "C"{
#endif
               after line 440 (#endif) and
#ifdef __cplusplus
          }
#endif
               at end

I also modified LAPACK dpotf2.f (two places) to change the GO TO 30 on AJJ.Lt.0.0

to

            IF( AJJ.LE.1.0e-20 ) THEN
               AJJ = 1.0e100;
            ELSE
               AJJ = SQRT( AJJ )
            END IF

*/
class ClpMatrixBase;
class ClpCholeskyTaucs : public ClpCholeskyBase {
  
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
  /** Default constructor. */
  ClpCholeskyTaucs();
  /** Destructor  */
  virtual ~ClpCholeskyTaucs();
  // Copy
  ClpCholeskyTaucs(const ClpCholeskyTaucs&);
  // Assignment
  ClpCholeskyTaucs& operator=(const ClpCholeskyTaucs&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
  /**@name Data members */
   //@{
  /// Taucs matrix (== sparseFactor etc)
  taucs_ccs_matrix * matrix_;
  /// Taucs factot
  void * factorization_;
  /// sparseFactor.
  double * sparseFactor_;
  /// choleskyStart
  CoinBigIndex * choleskyStart_;
  /// choleskyRow
  int * choleskyRow_;
  /// sizeFactor.
  CoinBigIndex sizeFactor_;
  /// Row copy of matrix
  ClpMatrixBase * rowCopy_;
  //@}
};

#endif
