// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpMatrixBase_H
#define ClpMatrixBase_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
using std::min;
using std::max;

#include "CoinPackedMatrix.hpp"
class CoinIndexedVector;
class ClpSimplex;
typedef int ClpBigIndex;

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
  virtual int getNumElements() const = 0;
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

   virtual const int * getVectorStarts() const = 0;
   /** The lengths of the major-dimension vectors. */
   virtual const int * getVectorLengths() const = 0 ;
    /** Delete the columns whose indices are listed in <code>indDel</code>. */
    virtual void deleteCols(const int numDel, const int * indDel) = 0;
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
    virtual void deleteRows(const int numDel, const int * indDel) = 0;

  /** Returns a new matrix in reverse order without gaps
      Is allowed to return NULL if doesn't want to have row copy */
  virtual ClpMatrixBase * reverseOrderedCopy() const {return NULL;};

  /** Returns number of elements in basis
      column is basic if entry >=0 */
  virtual ClpBigIndex numberInBasis(const int * columnIsBasic) const = 0;
  /// Fills in basis (Returns number of elements and updates numberBasic)
  virtual ClpBigIndex fillBasis(const ClpSimplex * model,
				const int * columnIsBasic, int & numberBasic,
				int * row, int * column,
				double * element) const = 0;
  /** Creates scales for column copy (rowCopy in model may be modified)
      default does not allow scaling
      returns non-zero if no scaling done */
  virtual int scale(ClpSimplex * model) const 
  { return 1;};
  /// Creates row copy and scales if necessary
  virtual ClpMatrixBase * scaledRowCopy(ClpSimplex * model) const
  { return reverseOrderedCopy();};

  /** Checks if all elements are in valid range.  Can just
      return true if you are not paranoid.  For Clp I will
      probably expect no zeros.  Code can modify matrix to get rid of
      small elements.
  */
  virtual bool allElementsInRange(ClpSimplex * model,
				  double smallest, double largest)
  { return true;};

  /** Unpacks a column into an CoinIndexedvector
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation (just for this variable) */
  virtual void unpack(ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column) const =0;
  /** Purely for column generation and similar ideas.  Allows
      matrix and any bounds or costs to be updated (sensibly).
      Returns non-zero if any changes.
  */
  int refresh(ClpSimplex * model)
    { return 0;};

  /** Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
  virtual void add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column, double multiplier) const =0;
   /// Allow any parts of a created CoinPackedMatrix to be deleted
   virtual void releasePackedMatrix() const {};
   //@}

  //---------------------------------------------------------------------------
  /**@name Matrix times vector methods 
     They can be faster if scalar is +- 1
     Also for simplex I am not using basic/non-basic split */
  //@{
    /** Return <code>y + A * x * scalar</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numColumns()</code>
        @precond <code>y<code> must be of size <code>numRows()</code> */
  virtual void times(double scalar,
		       const double * x, double * y) const=0;
  /// And for scaling - default aborts for when scaling not supported
  virtual void times(double scalar,
		     const double * x, double * y,
		     const double * rowScale, 
		     const double * columnScale) const;
    /** Return <code>y + x * scalar * A</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numRows()</code>
        @precond <code>y<code> must be of size <code>numColumns()</code> */
    virtual void transposeTimes(double scalar,
				const double * x, double * y) const = 0;
  /// And for scaling - default aborts for when scaling not supported
    virtual void transposeTimes(double scalar,
				const double * x, double * y,
				const double * rowScale, 
				const double * columnScale) const;
    /** Return <code>x * scalar *A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const = 0;
    /** Return <code>x *A in <code>z</code> but
	just for indices in y.
	This is only needed for primal steepest edge.
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
 
  /// Returns type
  inline int type() const
  { return type_;};
  /// Sets type
  void setType(int type) {type_=type;};
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
   /// type (may be useful)
   int type_;
   //@}
};

#endif
