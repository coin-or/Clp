// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPrimalColumnSteepest_H
#define ClpPrimalColumnSteepest_H

#include "ClpPrimalColumnPivot.hpp"
#include <bitset>

//#############################################################################


/** Primal Column Pivot Steepest Edge Algorithm Class

See Forrest-Goldfarb paper for algorithm

*/

class CoinIndexedVector;

class ClpPrimalColumnSteepest : public ClpPrimalColumnPivot {
  
public:
  
  ///@name Algorithmic methods 
  //@{
  
  /** Returns pivot column, -1 if none.
      updateArray has cost updates (also use pivotRow_ from last iteration).
      Parts of operation split out into seperate functions for
      profiling and speed
  */
  virtual int pivotColumn(CoinIndexedVector * updates,
			  CoinIndexedVector * spareRow1,
			  CoinIndexedVector * spareRow2,
			  CoinIndexedVector * spareColumn1,
			  CoinIndexedVector * spareColumn2);
  /// For quadratic or funny nonlinearities
  int pivotColumnOldMethod(CoinIndexedVector * updates,
			  CoinIndexedVector * spareRow1,
			  CoinIndexedVector * spareRow2,
			  CoinIndexedVector * spareColumn1,
			  CoinIndexedVector * spareColumn2);
  /// Just update djs
  void justDjs(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update djs, weights for Devex using djs
  void djsAndDevex(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update djs, weights for Steepest using djs
  void djsAndSteepest(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update djs, weights for Devex using pivot row
  void djsAndDevex2(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update djs, weights for Steepest using pivot row
  void djsAndSteepest2(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update weights for Devex
  void justDevex(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);
  /// Update weights for Steepest
  void justSteepest(CoinIndexedVector * updates,
	       CoinIndexedVector * spareRow1,
	       CoinIndexedVector * spareRow2,
	       CoinIndexedVector * spareColumn1,
	       CoinIndexedVector * spareColumn2);

  /// Updates weights - part 1 - also checks accuracy
  virtual void updateWeights(CoinIndexedVector * input);

  /// Checks accuracy - just for debug
  void checkAccuracy(int sequence,double relativeTolerance,
		     CoinIndexedVector * rowArray1,
		     CoinIndexedVector * rowArray2);

  /// Initialize weights
  void initializeWeights();

  /// Save weights
  virtual void saveWeights(ClpSimplex * model,int mode);
  /// Gets rid of last update
  virtual void unrollWeights();
  /// Gets rid of all arrays
  virtual void clearArrays();
  /// Returns true if would not find any column
  virtual bool looksOptimal() const;
  //@}
  
  /**@name gets and sets */
  //@{ 
  /// Mode
  inline int mode() const
    { return mode_;};
  /** Returns number of extra columns for sprint algorithm - 0 means off.
      Also number of iterations before recompute
  */
  virtual int numberSprintColumns(int & numberIterations) const;
  /// Switch off sprint idea
  virtual void switchOffSprint();
  
 //@}

  
  ///@name Constructors and destructors
  //@{
  /** Default Constructor 
      0 is exact devex, 1 full steepest, 2 is partial exact devex
      3 switches between 0 and 2 depending on factorization
      4 starts as partial dantzig/devex but then may switch between 0 and 2.
      By partial exact devex is meant that the weights are updated as normal
      but only part of the nonbasic variables are scanned.  
      This can be faster on very easy problems.
  */
  ClpPrimalColumnSteepest(int mode=3); 
  
  /// Copy constructor 
  ClpPrimalColumnSteepest(const ClpPrimalColumnSteepest &);
  
  /// Assignment operator 
  ClpPrimalColumnSteepest & operator=(const ClpPrimalColumnSteepest& rhs);
  
  /// Destructor 
  virtual ~ClpPrimalColumnSteepest ();

  /// Clone
  virtual ClpPrimalColumnPivot * clone(bool copyData = true) const;
 
  //@}

  ///@name Private functions to deal with devex 
  /** reference would be faster using ClpSimplex's status_,
      but I prefer to keep modularity.
  */
  inline bool reference(int i) const {
    return ((reference_[i>>5]>>(i&31))&1)!=0;
  }
  inline void setReference(int i,bool trueFalse) {
    unsigned int & value = reference_[i>>5];
    int bit = i&31;
    if (trueFalse)
      value |= (1<<bit);
    else
      value &= ~(1<<bit);
  }
  //@}
  //---------------------------------------------------------------------------
  
private:
  ///@name Private member data 
  // Update weight
  double devex_;
  /// weight array 
  double * weights_;
  /// square of infeasibility array (just for infeasible columns)
  CoinIndexedVector * infeasible_;
  /// alternate weight array (so we can unroll)
  CoinIndexedVector * alternateWeights_;
  /// save weight array (so we can use checkpoint)
  double * savedWeights_;
  // Array for exact devex to say what is in reference framework
  unsigned int * reference_;
  /** Status
      0) Normal 
      -1) Needs initialization
      1) Weights are stored by sequence number
  */
  int state_;
  /**
      0 is exact devex, 1 full steepest, 2 is partial exact devex
      3 switches between 0 and 2 depending on factorization
      4 starts as partial dantzig/devex but then may switch between 0 and 2.
      By partial exact devex is meant that the weights are updated as normal
      but only part of the nonbasic variables are scanned.  
      This can be faster on very easy problems.

      New dubious option is >=10 which does mini-sprint
  */
  int mode_;
  /// Number of times switched from partial dantzig to 0/2
  int numberSwitched_;
  // This is pivot row (or pivot sequence round re-factorization)
  int pivotSequence_;  
  // This is saved pivot sequence
  int savedPivotSequence_;  
  // This is saved outgoing variable
  int savedSequenceOut_;  
  //@}
};

#endif
