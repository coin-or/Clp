// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPrimalColumnSteepest_H
#define ClpPrimalColumnSteepest_H

#include "ClpPrimalColumnPivot.hpp"

//#############################################################################


/** Primal Column Pivot Steepest Edge Algorithm Class

See Forrest-Goldfarb paper for algorithm

*/

class OsiIndexedVector;

class ClpPrimalColumnSteepest : public ClpPrimalColumnPivot {
  
public:
  
  ///@name Algorithmic methods 
  //@{
  
  /** Returns pivot column, -1 if none.
      updateArray has cost updates (also use pivotRow_ from last iteration)
  */
  virtual int pivotColumn(OsiIndexedVector * updates,
			  OsiIndexedVector * spareRow1,
			  OsiIndexedVector * spareRow2,
			  OsiIndexedVector * spareColumn1,
			  OsiIndexedVector * spareColumn2);

  /// Updates weights - part 1 - also checks accuracy
  virtual void updateWeights(OsiIndexedVector * input);

  /// Checks accuracy - just for debug
  void checkAccuracy(int sequence,double relativeTolerance,
		     OsiIndexedVector * rowArray1,
		     OsiIndexedVector * rowArray2);

  /// Initialize weights
  void initializeWeights();

  /// Save weights
  virtual void saveWeights(ClpSimplex * model,int mode);
  /// Gets rid of last update
  virtual void unrollWeights();
  //@}
  
  
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  ClpPrimalColumnSteepest(int mode=0); 
  
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
    return (reference_[i>>5]>>(i&31))!=0;
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
  /** Status
      0) Normal 
      -1) Needs initialization
      1) Weights are stored by sequence number
  */
  int state_;
  /// If 0 then we are using exact devex, 1 then full
  int mode_;
  /// weight array 
  double * weights_;
  /// square of infeasibility array (just for infeasible columns)
  OsiIndexedVector * infeasible_;
  /// alternate weight array (so we can unroll)
  OsiIndexedVector * alternateWeights_;
  /// save weight array (so we can use checkpoint)
  double * savedWeights_;
  // This is pivot row (or pivot sequence round re-factorization)
  int pivotSequence_;  
  // This is saved pivot sequence
  int savedPivotSequence_;  
  // This is saved outgoing variable
  int savedSequenceOut_;  
  // Array for exact devex to say what is in reference framework
  unsigned int * reference_;
  // Update weight
  double devex_;
  //@}
};

#endif
