// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Tomlin

 */
#ifndef ClpPdco_H
#define ClpPdco_H

#include "ClpInterior.hpp"

/** This solves problems in Primal Dual Convex Optimization

    It inherits from ClpInterior.  It has no data of its own and 
    is never created - only cast from a ClpInterior object at algorithm time. 

*/

class ClpPdco : public ClpInterior {

public:

  /**@name Description of algorithm */
  //@{
  /** Pdco algorithm

      Method


  */

  int pdco();
  // ** Temporary version
  int  pdco( Lsqr *lsqr, Options &options, Info &info, Outfo &outfo);

  //@}

  /**@name Functions used in dual */
  //@{
  /// LSQR
  void lsqr();
  //@}
};
#endif
