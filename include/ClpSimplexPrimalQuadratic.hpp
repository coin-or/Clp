// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpSimplexPrimalQuadratic_H
#define ClpSimplexPrimalQuadratic_H

class ClpQuadraticInfo;

#include "ClpSimplexPrimal.hpp"

/** This solves Quadratic LPs using the primal simplex method

    It inherits from ClpSimplexPrimal.  It has no data of its own and 
    is never created - only cast from a ClpSimplexPrimal object at algorithm time. 
    If needed create new class and pass around

*/

class ClpSimplexPrimalQuadratic : public ClpSimplexPrimal {

public:

  /**@name Description of algorithm */
  //@{
  /** Primal algorithms for quadratic
      At present we have two algorithms:

      a) Dantzig's algorithm 
      b) Using a semi-trust region approach as for pooling problem
         This is in because I have it lying around

  */
  /// A sequential LP method
  int primalSLP(int numberPasses, double deltaTolerance);
  /** Dantzig's method (actually a mixture with Jensen and King)
      phase - 0 normal, 1 getting complementary solution,
      2 getting basic solution. 
      Returns 0 if okay, 1 if LP infeasible.
  */
  int primalQuadratic(int phase=2);
  /** This is done after first pass
      phase - 0 normal, 1 getting complementary solution,
      2 getting basic solution. */
  int primalQuadratic2 (ClpQuadraticInfo * info,
			int phase=2);
  /** This creates the large version of QP and
      fills in quadratic information.
      Returns NULL if no quadratic information
  */
  ClpSimplexPrimalQuadratic * makeQuadratic(ClpQuadraticInfo & info);

  /// This moves solution back
  int endQuadratic(ClpSimplexPrimalQuadratic * quadraticModel,
		   ClpQuadraticInfo & info);
  /// Just for debug
  int checkComplementary (const ClpQuadraticInfo * info);
  /** Main part.
      phase - 0 normal, 1 getting complementary solution,
      2 getting basic solution. */
  int whileIterating (int & sequenceIn,
		      ClpQuadraticInfo * info,
		      int phase);
  /** 
      Row array has pivot column
      This chooses pivot row.
      Rhs array is used for distance to next bound (for speed)
      For speed, we may need to go to a bucket approach when many
      variables go through bounds
      On exit rhsArray will have changes in costs of basic variables
      Initially no go thru
      Returns 0 - can do normal iteration
      1 - losing complementarity
  */
  int primalRow(CoinIndexedVector * rowArray,
		CoinIndexedVector * rhsArray,
		CoinIndexedVector * spareArray,
		CoinIndexedVector * spareArray2,
		ClpQuadraticInfo * info,
		bool cleanupIteration);
  //@}

};

/** Trivial class to keep quadratic iterating info around

*/

class ClpQuadraticInfo  {
  
public:
  
public:

  /**@name Constructors, destructor */
  //@{
  /// Default constructor. 
  ClpQuadraticInfo();
  /** Constructor from original model
  */
  ClpQuadraticInfo(const ClpSimplex * model);
  /// Destructor
  ~ClpQuadraticInfo();
  // Copy
  ClpQuadraticInfo(const ClpQuadraticInfo&);
  // Assignment
  ClpQuadraticInfo& operator=(const ClpQuadraticInfo&);
  //@}
     

  /**@name Gets and sets */
  //@{
  /// Number of Original columns
  inline int numberXColumns() const
  {return numberXColumns_;};
  /// Number of Quadratic columns
  inline int numberQuadraticColumns() const
  {return numberQuadraticColumns_;};
  /// Number of Original rows
  inline int numberXRows() const
  {return numberXRows_;};
  /// Sequence number of binding Sj
  inline int crucialSj() const
  {return crucialSj_;};
  inline void setCrucialSj(int sequence) 
  {crucialSj_=sequence;};
  /// Original objective
  inline const double * originalObjective() const
  {return originalModel_->objective();};
  /// Quadratic sequence or -1 if linear
  inline const int * quadraticSequence() const
  {return quadraticSequence_;};
  /// Which ones are quadratic
  inline const int * backSequence() const
  {return backSequence_;};
  /// Returns pointer to original model
  inline const ClpSimplex * originalModel() const
  { return originalModel_;};
  //@}
    
private:
  /**@name Data members */
  //@{
  /// Model
  const ClpSimplex * originalModel_;
  /// Quadratic sequence
  int * quadraticSequence_;
  /// Which ones are quadratic
  int * backSequence_;
  /// Crucial Sj
  int crucialSj_;
  /// Number of original rows
  int numberXRows_;
  /// Number of original columns 
  int numberXColumns_;
  /// Number of quadratic columns 
  int numberQuadraticColumns_;
  //@}
};
#endif

