// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPrimalcolumnPivot_H
#define ClpPrimalcolumnPivot_H

class ClpSimplex;
class OsiIndexedVector;

//#############################################################################

/** Primal Column Pivot Abstract Base Class

Abstract Base Class for describing an interface to an algorithm
to choose column pivot in primal simplex algorithm.  For some algorithms
e.g. Dantzig choice then some functions may be null.

*/

class ClpPrimalColumnPivot  {
  
public:
  
  ///@name Algorithmic methods 
  //@{
  
  /** Returns pivot column, -1 if none
      updateArray has cost updates (also use pivotRow_ from last iteration)
      We can use other arrays to help updates
  */
  virtual int pivotColumn(OsiIndexedVector * updates,
			  OsiIndexedVector * spareRow1,
			  OsiIndexedVector * spareRow2,
			  OsiIndexedVector * spareColumn1,
			  OsiIndexedVector * spareColumn2) = 0;
  
  /// Updates weights - part 1 (may be empty)
  virtual void updateWeights(OsiIndexedVector * input);
  
  /** Saves any weights round factorization as pivot rows may change
      Will be empty unless steepest edge (will save model)
      May also recompute infeasibility stuff
      1) before factorization
      2) after good factorization (if weights empty may initialize)
      3) after something happened but no factorization 
         (e.g. check for infeasible)
      4) as 2 but restore weights from previous snapshot
      5) forces some initialization e.g. weights
      Also sets model
  */
  virtual void saveWeights(ClpSimplex * model,int mode)=0;
  /** Signals pivot row choice:
      -2 (default) - use normal pivot row choice
      -1 to numberRows-1 - use this (will be checked)
      way should be -1 to go to lower bound, +1 to upper bound
  */
  virtual int pivotRow(double & way)
  {way=0;return -2;};
  //@}
  
  
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  ClpPrimalColumnPivot(); 
  
  /// Copy constructor 
  ClpPrimalColumnPivot(const ClpPrimalColumnPivot &);
  
  /// Assignment operator 
  ClpPrimalColumnPivot & operator=(const ClpPrimalColumnPivot& rhs);
  
  /// Destructor 
  virtual ~ClpPrimalColumnPivot ();

  /// Clone
  virtual ClpPrimalColumnPivot * clone(bool copyData = true) const = 0;
 
  //@}

  ///@name Other
  //@{
  /// Returns model
  inline ClpSimplex * model()
  { return model_;};
  
  /// Returns type
  inline int type()
  { return type_;};
  
  //@}

  //---------------------------------------------------------------------------
  
protected:
  ///@name Protected member data 
  //@{
  /// Pointer to model
  ClpSimplex * model_;
  /// Type of column pivot algorithm 
  int type_;
  //@}
};

#endif
