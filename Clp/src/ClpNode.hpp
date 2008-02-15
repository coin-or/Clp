// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpNode_H
#define ClpNode_H

#include "CoinPragma.hpp"

// This implements all stuff for Clp fathom
/** This contains what is in a Clp "node"
    
 */

class ClpFactorization;
class ClpDualRowSteepest;
class ClpNodeStuff;
class ClpNode {
  
public:
  /**@name Useful methods */
  //@{
  /// Applies node to model
  void applyNode(ClpSimplex * model, bool justBounds );
  /// Fix on reduced costs
  int fixOnReducedCosts(ClpSimplex * model);
  //@}
  
  /**@name Gets and sets */
  //@{
  /// Initial value of integer variable
  inline double branchingValue() const
  { return branchingValue_;}
  /** Way for integer variable -1 down , +1 up */
  int way() const;
  /// Return true if branch exhausted
  bool fathomed() const;
  /// Change state of variable i.e. go other way
  void changeState();
  /// Sequence number of integer variable (-1 if none)
  inline int sequence() const
  { return sequence_;}
  //@}
  
  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  ClpNode();
  /// Constructor from model
  ClpNode (const ClpSimplex * model, const ClpNodeStuff * stuff);
  /// Does work of constructor (partly so gdb will work)
  void gutsOfConstructor(const ClpSimplex * model, const ClpNodeStuff * stuff);
  /** Destructor */
  virtual ~ClpNode();
  //@}
  
  /**@name Copy methods (at present illegal - will abort) */
  //@{
  /** The copy constructor. */
  ClpNode(const ClpNode&);
  /// Operator =
  ClpNode& operator=(const ClpNode&);
  //@}
  
protected:
// For state of branch
typedef struct {
  unsigned int firstBranch:1; //  nonzero if first branch on variable is up
  unsigned int branch:2; //  0 means do first branch next, 1 second, 2 finished
  unsigned int spare:29;
} branchState;
  /**@name Data */
  //@{
  /// Initial value of integer variable
  double branchingValue_;
  /// Factorization
  ClpFactorization * factorization_;
  /// Steepest edge weights
  ClpDualRowSteepest * weights_;
  /// Status vector
  unsigned char * status_;
  /// Primal solution
  double * primalSolution_;
  /// Dual solution
  double * dualSolution_;
  /// Pivot variables for factorization
  int * pivotVariables_;
  /// Variables fixed by reduced costs (at end of branch) 0x10000000 added if fixed to UB
  int * fixed_;
  /// State of branch
  branchState branchState_;
  /// Sequence number of integer variable (-1 if none)
  int sequence_;
  /// Number fixed by reduced cost
  int numberFixed_;
  //@}
};
class ClpNodeStuff {
  
public:
  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  ClpNodeStuff();
  /** Destructor */
  virtual ~ClpNodeStuff();
  //@}
  
  /**@name Copy methods (at present illegal - will abort) */
  //@{
  /** The copy constructor. */
  ClpNodeStuff(const ClpNodeStuff&);
  /// Operator =
  ClpNodeStuff& operator=(const ClpNodeStuff&);
  //@}
  
public:
  /**@name Data */
  //@{
  /// Integer tolerance
  double integerTolerance_;
  /// Integer increment
  double integerIncrement_;
  //@}
};
#endif
