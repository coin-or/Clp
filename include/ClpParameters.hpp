// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef _ClpParameters_H
#define _ClpParameters_H

enum ClpIntParam {
   /** The maximum number of iterations Clp can execute in the simplex methods
    */
  ClpMaxNumIteration = 0,
  /** The maximum number of iterations Clp can execute in hotstart before
      terminating */
  ClpMaxNumIterationHotStart,
  /** Just a marker, so that we can allocate a static sized array to store
      parameters. */
  ClpLastIntParam
};

enum ClpDblParam {
  /** Set Dual objective limit. This is to be used as a termination criteria
      in methods where the dual objective monotonically changes (dual
      simplex). */
  ClpDualObjectiveLimit,
  /** Primal objective limit. This is to be used as a termination
      criteria in methods where the primal objective monotonically changes
      (e.g., primal simplex) */
  ClpPrimalObjectiveLimit,
  /** The maximum amount the dual constraints can be violated and still be
      considered feasible. */
  ClpDualTolerance,
  /** The maximum amount the primal constraints can be violated and still be
      considered feasible. */
  ClpPrimalTolerance,
  /** Objective function constant. This the value of the constant term in
      the objective function. */
  ClpObjOffset,
  /** Just a marker, so that we can allocate a static sized array to store
      parameters. */
  ClpLastDblParam
};


enum ClpStrParam {
  /** Name of the problem. This is the found on the Name card of
      an mps file. */
  ClpProbName = 0,
  /** Just a marker, so that we can allocate a static sized array to store
      parameters. */
  ClpLastStrParam
};

#endif
