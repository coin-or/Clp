// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef ClpPresolve_H
#define ClpPresolve_H
#include "ClpSimplex.hpp"

class CoinPresolveAction;
#include "CoinPresolveMatrix.hpp"
class ClpPresolve {
public:
  /// Default constructor
  ClpPresolve();

  /// Virtual destructor
  virtual ~ClpPresolve();
  //@}
  /**@name presolve - presolves a model, transforming the model
   * and saving information in the ClpPresolve object needed for postsolving.
   * This underlying (protected) method is virtual; the idea is that in the future,
   * one could override this method to customize how the various
   * presolve techniques are applied.

   This version of presolve returns a pointer to a new presolved 
      model.  NULL if infeasible or unbounded.  
      This should be paired with postsolve
      below.  The advantage of going back to original model is that it
      will be exactly as it was i.e. 0.0 will not become 1.0e-19.
      If keepIntegers is true then bounds may be tightened in
      original.  Bounds will be moved by up to feasibilityTolerance
      to try and stay feasible.
      Names will be dropped in presolved model if asked
  */
  ClpSimplex * presolvedModel(ClpSimplex & si,
				      double feasibilityTolerance=0.0,
				      bool keepIntegers=true,
				      int numberPasses=5,
				      bool dropNames=false);
  /** This version saves data in a file.  The passed in model
      is updated to be presolved model.  names are always dropped.
      Returns non-zero if infeasible*/
  int presolvedModel(ClpSimplex &si,std::string fileName,
			      double feasibilityTolerance=0.0,
			      bool keepIntegers=true,
			      int numberPasses=5);
  /** Return pointer to presolved model,
      Up to user to destroy */
  ClpSimplex * model() const;
  /// Return pointer to original model
  ClpSimplex * originalModel() const;
  /// Set pointer to original model
  void setOriginalModel(ClpSimplex * model);
    
  /// return pointer to original columns
  const int * originalColumns() const;
  /// return pointer to original rows
  const int * originalRows() const;
  /** "Magic" number. If this is non-zero then any elements with this value
      may change and so presolve is very limited in what can be done
      to the row and column.  This is for non-linear problems.
  */
  inline void setNonLinearValue(double value)
  { nonLinearValue_ = value;};
  inline double nonLinearValue() const
    { return nonLinearValue_;};

  /**@name postsolve - postsolve the problem.  If the problem 
    has not been solved to optimality, there are no guarantees.
   If you are using an algorithm like simplex that has a concept
   of "basic" rows/cols, then set updateStatus
  
   Note that if you modified the original problem after presolving,
   then you must ``undo'' these modifications before calling postsolve.
  This version updates original*/
  virtual void postsolve(bool updateStatus=true);

  /**@name private or protected data */
private:
  /// Original model - must not be destroyed before postsolve
  ClpSimplex * originalModel_;

  /// ClpPresolved model - up to user to destroy by deleteClpPresolvedModel
  ClpSimplex * presolvedModel_;
  /** "Magic" number. If this is non-zero then any elements with this value
      may change and so presolve is very limited in what can be done
      to the row and column.  This is for non-linear problems.
      One could also allow for cases where sign of coefficient is known.
  */
  double nonLinearValue_;
  /// Original column numbers
  int * originalColumn_;
  /// Original row numbers
  int * originalRow_;
  /// The list of transformations applied.
  const CoinPresolveAction *paction_;

  /// The postsolved problem will expand back to its former size
  /// as postsolve transformations are applied.
  /// It is efficient to allocate data structures for the final size
  /// of the problem rather than expand them as needed.
  /// These fields give the size of the original problem.
  int ncols_;
  int nrows_;
  CoinBigIndex nelems_;
  /// Number of major passes
  int numberPasses_;
  /// Name of saved model file
  std::string saveFile_;

protected:
  /// If you want to apply the individual presolve routines differently,
  /// or perhaps add your own to the mix,
  /// define a derived class and override this method
  virtual const CoinPresolveAction *presolve(CoinPresolveMatrix *prob);

  /// Postsolving is pretty generic; just apply the transformations
  /// in reverse order.
  /// You will probably only be interested in overriding this method
  /// if you want to add code to test for consistency
  /// while debugging new presolve techniques.
  virtual void postsolve(CoinPostsolveMatrix &prob);
  /// Gets rid of presolve actions (e.g.when infeasible)
  void gutsOfDestroy();
  /** This is main part of Presolve */
  virtual ClpSimplex * gutsOfPresolvedModel(ClpSimplex * originalModel
					    ,double feasibilityTolerance,
					    bool keepIntegers,
					    int numberPasses,
					    bool dropNames);
};
#endif
