// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpNonLinearCost_H
#define ClpNonLinearCost_H


#include "CoinPragma.hpp"

class ClpSimplex;
class CoinIndexedVector;

/** Trivial class to deal with non linear costs

    I don't make any explicit assumptions about convexity but I am
    sure I do make implicit ones.

    One interesting idea for normal LP's will be to allow non-basic
    variables to come into basis as infeasible i.e. if variable at
    lower bound has very large positive reduced cost (when problem
    is infeasible) could it reduce overall problem infeasibility more
    by bringing it into basis below its lower bound.

    Another feature would be to automatically discover when problems
    are convex piecewise linear and re-formulate to use non-linear.
    I did some work on this many years ago on "grade" problems, but
    while it improved primal interior point algorithms were much better
    for that particular problem.
*/

class ClpNonLinearCost  {
  
public:
  
public:

  /**@name Constructors, destructor */
  //@{
  /// Default constructor. 
  ClpNonLinearCost();
  /** Constructor from simplex.
      This will just set up wasteful arrays for linear, but
      later may do dual analysis and even finding duplicate columns .
  */
  ClpNonLinearCost(ClpSimplex * model);
  /** Constructor from simplex and list of non-linearities (columns only)
      First lower of each column has to match real lower
      Last lower has to be <= upper (if == then cost ignored)
      This could obviously be changed to make more user friendly
  */
  ClpNonLinearCost(ClpSimplex * model,const int * starts,
		   const double * lower, const double * cost);
  /// Destructor
  ~ClpNonLinearCost();
  // Copy
  ClpNonLinearCost(const ClpNonLinearCost&);
  // Assignment
  ClpNonLinearCost& operator=(const ClpNonLinearCost&);
  //@}
     

  /**@name Actual work in primal */
  //@{
  /** Changes infeasible costs and computes number and cost of infeas
      Puts all non-basic (non free) variables to bounds
      and all free variables to zero if oldTolerance is non-zero
      - but does not move those <= oldTolerance away*/
  void checkInfeasibilities(double oldTolerance=0.0);
  /** Changes infeasible costs for each variable
      The indices are row indices and need converting to sequences
  */
  void checkInfeasibilities(int numberInArray, const int * index);
  /** Puts back correct infeasible costs for each variable
      The input indices are row indices and need converting to sequences
      for costs.
      On input array is empty (but indices exist).  On exit just
      changed costs will be stored as normal CoinIndexedVector
  */
  void checkChanged(int numberInArray, CoinIndexedVector * update);
  /** Goes through one bound for each variable.
      If multiplier*work[iRow]>0 goes down, otherwise up.
      The indices are row indices and need converting to sequences
      Temporary offsets may be set
      Rhs entries are increased
  */
  void goThru(int numberInArray, double multiplier,
	      const int * index, const double * work,
	      double * rhs);
  /** Takes off last iteration (i.e. offsets closer to 0)
  */
  void goBack(int numberInArray, const int * index, 
	      double * rhs);
  /** Puts back correct infeasible costs for each variable
      The input indices are row indices and need converting to sequences
      for costs.
      At the end of this all temporary offsets are zero
  */
  void goBackAll(const CoinIndexedVector * update);
  /** Sets bounds and cost for one variable 
      Returns change in cost
   May need to be inline for speed */
  double setOne(int sequence, double solutionValue);
  /** Sets bounds and cost for outgoing variable 
      may change value
      Returns direction */
  int setOneOutgoing(int sequence, double &solutionValue);
  /// Returns nearest bound
  double nearest(int sequence, double solutionValue);
  /** Returns change in cost - one down if alpha >0.0, up if <0.0
      Value is current - new
   */
  inline double changeInCost(int sequence, double alpha) const
  {
    int iRange = whichRange_[sequence]+offset_[sequence];
    if (alpha>0.0)
      return cost_[iRange]-cost_[iRange-1];
    else
      return cost_[iRange]-cost_[iRange+1];
  }
  inline double changeUpInCost(int sequence) const
  {
    int iRange = whichRange_[sequence]+offset_[sequence];
    if (iRange+1!=start_[sequence+1]&&!infeasible(iRange+1))
      return cost_[iRange]-cost_[iRange+1];
    else
      return -1.0e100;
  }
  inline double changeDownInCost(int sequence) const
  {
    int iRange = whichRange_[sequence]+offset_[sequence];
    if (iRange!=start_[sequence]&&!infeasible(iRange-1))
      return cost_[iRange]-cost_[iRange-1];
    else
      return 1.0e100;
  }
  /// Returns current lower bound
  inline double lower(int sequence) const
  { return lower_[whichRange_[sequence]+offset_[sequence]];};
  /// Returns current upper bound
  inline double upper(int sequence) const
  { return lower_[whichRange_[sequence]+offset_[sequence]+1];};
  /// Returns current cost
  inline double cost(int sequence) const
  { return cost_[whichRange_[sequence]+offset_[sequence]];};
  //@}


  /**@name Gets and sets */
  //@{
  /// Number of infeasibilities
  inline int numberInfeasibilities() const
  {return numberInfeasibilities_;};
  /// Change in cost
  inline double changeInCost() const
  {return changeCost_;};
  /// Sum of infeasibilities
  inline double sumInfeasibilities() const
  {return sumInfeasibilities_;};
  /// Largest infeasibility
  inline double largestInfeasibility() const
  {return largestInfeasibility_;};
  inline void setChangeInCost(double value) 
  {changeCost_ = value;};
  /// See if may want to look both ways
  inline bool lookBothWays() const
  { return bothWays_;};
  //@}
  ///@name Private functions to deal with infeasible regions 
  inline bool infeasible(int i) const {
    return ((infeasible_[i>>5]>>(i&31))&1)!=0;
  }
  inline void setInfeasible(int i,bool trueFalse) {
    unsigned int & value = infeasible_[i>>5];
    int bit = i&31;
    if (trueFalse)
      value |= (1<<bit);
    else
      value &= ~(1<<bit);
  }
  //@}
    
private:
  /**@name Data members */
  //@{
  /// Change in cost because of infeasibilities
  double changeCost_;
  /// Largest infeasibility
  double largestInfeasibility_;
  /// Sum of infeasibilities
  double sumInfeasibilities_;
  /// Number of rows (mainly for checking and copy)
  int numberRows_;
  /// Number of columns (mainly for checking and copy)
  int numberColumns_;
  /// Starts for each entry (columns then rows)
  int * start_;
  /// Range for each entry (columns then rows)
  int * whichRange_;
  /// Temporary range offset for each entry (columns then rows)
  int * offset_;
  /** Lower bound for each range (upper bound is next lower).
      For various reasons there is always an infeasible range
      at bottom - even if lower bound is - infinity */
  double * lower_;
  /// Cost for each range
  double * cost_;
  /// Model
  ClpSimplex * model_;
  // Array to say which regions are infeasible
  unsigned int * infeasible_;
  /// Number of infeasibilities found
  int numberInfeasibilities_;
  /// If all non-linear costs convex
  bool convex_;
  /// If we should look both ways for djs
  bool bothWays_;
  //@}
};

#endif
