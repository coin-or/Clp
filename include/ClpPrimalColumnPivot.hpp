// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPrimalcolumnPivot_H
#define ClpPrimalcolumnPivot_H

class ClpSimplex;
class CoinIndexedVector;

//#############################################################################

/** Primal Column Pivot Abstract Base Class

Abstract Base Class for describing an interface to an algorithm
to choose column pivot in primal simplex algorithm.  For some algorithms
e.g. Dantzig choice then some functions may be null.  For Dantzig
the only one of any importance is pivotColumn.

If  you wish to inherit from this look at ClpPrimalColumnDantzig.cpp
as that is simplest version.
*/

class ClpPrimalColumnPivot  {
  
public:
  
  ///@name Algorithmic methods 
  //@{
  
  /** Returns pivot column, -1 if none

      Normally updates reduced costs using result of last iteration
      before selecting incoming column.

      The CoinIndexedVector updates has cost updates - for normal LP
      that is just +-weight where a feasibility changed.

      Inside pivotColumn the pivotRow_ and reduced cost from last iteration
      are also used.

      So in the simplest case i.e. feasible we compute the row of the
      tableau corresponding to last pivot and add a multiple of this
      to current reduced costs.

      We can use other arrays to help updates
  */
  virtual int pivotColumn(CoinIndexedVector * updates,
			  CoinIndexedVector * spareRow1,
			  CoinIndexedVector * spareRow2,
			  CoinIndexedVector * spareColumn1,
			  CoinIndexedVector * spareColumn2) = 0;
  
  /// Updates weights - part 1 (may be empty)
  virtual void updateWeights(CoinIndexedVector * input);
  
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
  /// Gets rid of all arrays (may be empty)
  virtual void clearArrays();
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
  
  /// Returns type (above 63 is extra information)
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
