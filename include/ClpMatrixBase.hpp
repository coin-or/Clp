// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpMatrixBase_H
#define ClpMatrixBase_H

#include "CoinPragma.hpp"

#include "CoinPackedMatrix.hpp"
class CoinIndexedVector;
class ClpSimplex;
class ClpModel;

/** Abstract base class for Clp Matrices
    
Since this class is abstract, no object of this type can be created.

If a derived class provides all methods then all Clp algorithms
should work.  Some can be very inefficient e.g. getElements etc is
only used for tightening bounds for dual and the copies are
deleted.  Many methods can just be dummy i.e. abort(); if not
all features are being used.  So if column generation was being done
then it makes no sense to do steepest edge so there would be
no point providing subsetTransposeTimes.
*/

class ClpMatrixBase  {
  
public:
  /**@name Virtual methods that the derived classes must provide */
  //@{
  /// Return a complete CoinPackedMatrix
  virtual CoinPackedMatrix * getPackedMatrix() const = 0;
  /** Whether the packed matrix is column major ordered or not. */
  virtual bool isColOrdered() const = 0;
  /** Number of entries in the packed matrix. */
  virtual CoinBigIndex getNumElements() const = 0;
  /** Number of columns. */
  virtual int getNumCols() const = 0;
  /** Number of rows. */
  virtual int getNumRows() const = 0;
  
  /** A vector containing the elements in the packed matrix. Note that there
      might be gaps in this list, entries that do not belong to any
      major-dimension vector. To get the actual elements one should look at
      this vector together with vectorStarts and vectorLengths. */
  virtual const double * getElements() const = 0;
  /** A vector containing the minor indices of the elements in the packed
      matrix. Note that there might be gaps in this list, entries that do not
      belong to any major-dimension vector. To get the actual elements one
      should look at this vector together with vectorStarts and
      vectorLengths. */
  virtual const int * getIndices() const = 0;
  
  virtual const CoinBigIndex * getVectorStarts() const = 0;
  /** The lengths of the major-dimension vectors. */
  virtual const int * getVectorLengths() const = 0 ;
  /** Delete the columns whose indices are listed in <code>indDel</code>. */
  virtual void deleteCols(const int numDel, const int * indDel) = 0;
  /** Delete the rows whose indices are listed in <code>indDel</code>. */
  virtual void deleteRows(const int numDel, const int * indDel) = 0;
  /// Append Columns
  virtual void appendCols(int number, const CoinPackedVectorBase * const * columns);
  /// Append Rows
  virtual void appendRows(int number, const CoinPackedVectorBase * const * rows);
  
  /** Returns a new matrix in reverse order without gaps
      Is allowed to return NULL if doesn't want to have row copy */
  virtual ClpMatrixBase * reverseOrderedCopy() const {return NULL;};
  
  /** If element NULL returns number of elements in column part of basis,
      If not NULL fills in as well */
  virtual CoinBigIndex fillBasis(ClpSimplex * model,
				 const int * whichColumn, 
				 int numberRowBasic,
				 int & numberColumnBasic,
				 int * row, int * column,
				 double * element)  = 0;
  /** Creates scales for column copy (rowCopy in model may be modified)
      default does not allow scaling
      returns non-zero if no scaling done */
  virtual int scale(ClpSimplex * model) const 
  { return 1;};
  
  /** Checks if all elements are in valid range.  Can just
      return true if you are not paranoid.  For Clp I will
      probably expect no zeros.  Code can modify matrix to get rid of
      small elements.
  */
  virtual bool allElementsInRange(ClpModel * model,
				  double smallest, double largest)
  { return true;};
  /** Returns largest and smallest elements of both signs.
      Largest refers to largest absolute value.
      If returns zeros then can't tell anything */
  virtual void rangeOfElements(double & smallestNegative, double & largestNegative,
			       double & smallestPositive, double & largestPositive);
  
  /** Unpacks a column into an CoinIndexedvector
   */
  virtual void unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		      int column) const =0;
  /** Unpacks a column into an CoinIndexedvector
   ** in packed foramt
   Note that model is NOT const.  Bounds and objective could
   be modified if doing column generation (just for this variable) */
  virtual void unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int column) const =0;
  /** Purely for column generation and similar ideas.  Allows
      matrix and any bounds or costs to be updated (sensibly).
      Returns non-zero if any changes.
  */
  int refresh(ClpSimplex * model)
  { return 0;};
  /** Given positive integer weights for each row fills in sum of weights
      for each column (and slack).
      Returns weights vector
      Default returns vector of ones
  */
  virtual CoinBigIndex * dubiousWeights(const ClpSimplex * model,int * inputWeights) const;
  /** Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
  virtual void add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column, double multiplier) const =0;
  /** Adds multiple of a column into an array */
  virtual void add(const ClpSimplex * model,double * array,
		   int column, double multiplier) const =0;
  /// Allow any parts of a created CoinPackedMatrix to be deleted
  virtual void releasePackedMatrix() const =0;
  /// Says whether it can do partial pricing
  virtual bool canDoPartialPricing() const;
  /// Returns number of hidden rows e.g. gub
  virtual int hiddenRows() const;
  /// Partial pricing 
  virtual void partialPricing(ClpSimplex * model, int start, int end,
			      int & bestSequence, int & numberWanted);
  /** expands an updated column to allow for extra rows which the main
      solver does not know about and returns number added. 
      
      This will normally be a no-op - it is in for GUB but may get extended to
      general non-overlapping and embedded networks.
      
      mode 0 - extend
      mode 1 - delete etc
  */
  virtual int extendUpdated(ClpSimplex * model,CoinIndexedVector * update,int mode);
  /**
     utility primal function for dealing with dynamic constraints
     mode=0  - Set up before "update" and "times" for primal solution using extended rows
     mode=1  - Cleanup primal solution after "times" using extended rows.
     mode=2  - Check (or report on) primal infeasibilities
  */
  virtual void primalExpanded(ClpSimplex * model,int mode);
  /** 
      utility dual function for dealing with dynamic constraints
      mode=0  - Set up before "updateTranspose" and "transposeTimes" for duals using extended
                updates array (and may use other if dual values pass)
      mode=1  - Update dual solution after "transposeTimes" using extended rows.
      mode=2  - Compute all djs and compute key dual infeasibilities
      mode=3  - Report on key dual infeasibilities
      mode=4  - Modify before updateTranspose in partial pricing
  */
  virtual void dualExpanded(ClpSimplex * model,CoinIndexedVector * array,
			    double * other,int mode);
  /** 
      general utility function for dealing with dynamic constraints
      mode=0  - Create list of non-key basics in pivotVariable_ using
                number as numberBasic in and out
      mode=1  - Set all key variables as basic
      mode=2  - return number extra rows needed, number gives maximum number basic
      mode=3  - before replaceColumn
      mode=4  - return 1 if can do primal, 2 if dual, 3 if both
      mode=5  - save any status stuff (when in good state)
      mode=6  - restore status stuff
  */
  virtual int generalExpanded(ClpSimplex * model,int mode,int & number);
  /** 
     update information for a pivot (and effective rhs)
  */
  virtual int updatePivot(ClpSimplex * model,double oldInValue, double oldOutValue);
  //@}
  
  //---------------------------------------------------------------------------
  /**@name Matrix times vector methods 
     They can be faster if scalar is +- 1
     Also for simplex I am not using basic/non-basic split */
  //@{
  /** Return <code>y + A * x * scalar</code> in <code>y</code>.
      @pre <code>x</code> must be of size <code>numColumns()</code>
      @pre <code>y</code> must be of size <code>numRows()</code> */
  virtual void times(double scalar,
		     const double * x, double * y) const=0;
  /** And for scaling - default aborts for when scaling not supported
      (unless pointers NULL when as normal)
  */
  virtual void times(double scalar,
		     const double * x, double * y,
		     const double * rowScale, 
		     const double * columnScale) const;
  /** Return <code>y + x * scalar * A</code> in <code>y</code>.
      @pre <code>x</code> must be of size <code>numRows()</code>
      @pre <code>y</code> must be of size <code>numColumns()</code> */
  virtual void transposeTimes(double scalar,
			      const double * x, double * y) const = 0;
  /** And for scaling - default aborts for when scaling not supported
      (unless pointers NULL when as normal)
  */
  virtual void transposeTimes(double scalar,
			      const double * x, double * y,
			      const double * rowScale, 
			      const double * columnScale) const;
  /** Return <code>x * scalar *A + y</code> in <code>z</code>. 
      Can use y as temporary array (will be empty at end)
      Note - If x packed mode - then z packed mode
      Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const = 0;
  /** Return <code>x *A</code> in <code>z</code> but
      just for indices in y.
      This is only needed for primal steepest edge.
      Note - z always packed mode
      Squashes small elements and knows about ClpSimplex */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
				    const CoinIndexedVector * x,
				    const CoinIndexedVector * y,
				    CoinIndexedVector * z) const = 0;
  //@}
  //@{
  ///@name Other
  /// Clone
  virtual ClpMatrixBase * clone() const = 0;
  /** Subset clone (without gaps).  Duplicates are allowed
      and order is as given.
      Derived classes need not provide this as it may not always make
      sense */
  virtual ClpMatrixBase * subsetClone (
				       int numberRows, const int * whichRows,
				       int numberColumns, const int * whichColumns) const;
  
  /// Returns type
  inline int type() const
  { return type_;};
  /// Sets type
  void setType(int type) {type_=type;};
  /// Sets up an effective RHS
  void useEffectiveRhs(ClpSimplex * model);
  /** Returns effective RHS if it is being used.  This is used for long problems
      or big gub or anywhere where going through full columns is
      expensive.  This may re-compute */
  virtual double * effectiveRhs(ClpSimplex * model,bool forceRefresh=false,
				bool check=false);
  /// If effectiveRhs used this is iteration last refreshed
  inline int lastRefresh() const
  { return lastRefresh_;};
  /// If effectiveRhs used this is refresh frequency (0==off)
  inline int refreshFrequency() const
  { return refreshFrequency_;};
  inline void setRefreshFrequency(int value)
  { refreshFrequency_=value;};
  /// whether to skip dual checks most of time
  inline bool skipDualCheck() const
  { return skipDualCheck_;};
  inline void setSkipDualCheck(bool yes)
  { skipDualCheck_=yes;};
  //@}
  
  
protected:
  
  /**@name Constructors, destructor<br>
     <strong>NOTE</strong>: All constructors are protected. There's no need
     to expose them, after all, this is an abstract class. */
  //@{
  /** Default constructor. */
  ClpMatrixBase();
  /** Destructor (has to be public) */
public:
  virtual ~ClpMatrixBase();
protected:
  // Copy
  ClpMatrixBase(const ClpMatrixBase&);
  // Assignment
  ClpMatrixBase& operator=(const ClpMatrixBase&);
  //@}
  
  
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /** Effective RHS if it is being used.  This is used for long problems
      or big gub or anywhere where going through full columns is
      expensive */
  double * effectiveRhs_;
  /// type (may be useful)
  int type_;
  /// If effectiveRhs used this is iteration last refreshed
  int lastRefresh_;
  /// If effectiveRhs used this is refresh frequency (0==off)
  int refreshFrequency_;
  /// whether to skip dual checks most of time
  bool skipDualCheck_;
  //@}
};
// bias for free variables
#define FREE_BIAS 1.0e1
// Acceptance criteria for free variables
#define FREE_ACCEPT 1.0e2

#endif
