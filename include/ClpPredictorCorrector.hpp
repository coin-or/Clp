// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpPredictorCorrector_H
#define ClpPredictorCorrector_H

#include "ClpInterior.hpp"

/** This solves LPs using the predictor-corrector method.

    It is rather basic as Interior point is not my speciality

    It inherits from ClpInterior.  It has no data of its own and 
    is never created - only cast from a ClpInterior object at algorithm time. 

*/

class ClpPredictorCorrector : public ClpInterior {

public:

  /**@name Description of algorithm */
  //@{
  /** Primal Dual Predictor Corrector algorithm

      Method

      Big TODO
  */

  int solve();
  //@}

  /**@name Functions used in algorithm */
  //@{
  /// findStepLength.
  //phase  - 0 predictor
  //         1 corrector
  //         2 primal dual
  double findStepLength( int phase);
  /// findDirectionVector.
  double findDirectionVector(const int phase);
  /// createSolution.  Creates solution from scratch (- code if no memory)
  int createSolution();
  /// complementarityGap.  Computes gap
  //phase 0=as is , 1 = after predictor , 2 after corrector
  double complementarityGap(int & numberComplementarityPairs,int & numberComplementarityItems,
			    const int phase);
  /// setupForSolve.
  //phase 0=affine , 1 = corrector , 2 = primal-dual
  void setupForSolve(const int phase);
  /** Does solve. region1 is for deltaX (columns+rows), region2 for deltaPi (rows) */
  void solveSystem(double * region1, double * region2,
		   const double * region1In, const double * region2In,
		   const double * saveRegion1, const double * saveRegion2,
		   bool gentleRefine);
  /// sees if looks plausible change in complementarity
  bool checkGoodMove(const bool doCorrector,double & bestNextGap,
		     bool allowIncreasingGap);
  ///:  checks for one step size
  bool checkGoodMove2(double move,double & bestNextGap,
		      bool allowIncreasingGap);
  /// updateSolution.  Updates solution at end of iteration
  //returns number fixed
  int updateSolution(double nextGap);
  ///  Save info on products of affine deltaT*deltaW and deltaS*deltaZ
  double affineProduct();
  ///See exactly what would happen given current deltas
  void debugMove(int phase,double primalStep, double dualStep);
  //@}

};
#endif
