// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPackedMatrix_H
#define ClpPackedMatrix_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "ClpMatrixBase.hpp"

/** This implements OsiPackedMatrix as derived from ClpMatrixBase.

    It adds a few methods that know about model as well as matrix

    For details see OsiPackedMatrix */

class ClpPackedMatrix : public ClpMatrixBase {
  
public:
  /**@name Useful methods */
   //@{
   /// Return a complete OsiPackedMatrix
  virtual OsiPackedMatrix * getPackedMatrix() const { return matrix_;};
    /** Whether the packed matrix is column major ordered or not. */
    virtual bool isColOrdered() const { return matrix_->isColOrdered(); }
   /** Number of entries in the packed matrix. */
  virtual  int getNumElements() const { return matrix_->getNumElements(); }
   /** Number of columns. */
   virtual int getNumCols() const { return matrix_->getNumCols(); }
   /** Number of rows. */
  virtual int getNumRows() const { return matrix_->getNumRows(); };

   /** A vector containing the elements in the packed matrix. Note that there
	might be gaps in this list, entries that do not belong to any
	major-dimension vector. To get the actual elements one should look at
	this vector together with vectorStarts and vectorLengths. */
   virtual const double * getElements() const 
  { return matrix_->getElements();};
   /** A vector containing the minor indices of the elements in the packed
        matrix. Note that there might be gaps in this list, entries that do not
        belong to any major-dimension vector. To get the actual elements one
        should look at this vector together with vectorStarts and
        vectorLengths. */
   virtual const int * getIndices() const 
  { return matrix_->getIndices();};

   virtual const int * getVectorStarts() const 
  { return matrix_->getVectorStarts();};
   /** The lengths of the major-dimension vectors. */
   virtual const int * getVectorLengths() const 
  { return matrix_->getVectorLengths();} ;

    /** Delete the columns whose indices are listed in <code>indDel</code>. */
    virtual void deleteCols(const int numDel, const int * indDel)
  { matrix_->deleteCols(numDel,indDel);};
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
    virtual void deleteRows(const int numDel, const int * indDel)
  { matrix_->deleteRows(numDel,indDel);};
  /** Replace the elements of a vector.  The indices remain the same.
      This is only needed if scaling and a row copy is used.
      At most the number specified will be replaced.
      The index is between 0 and major dimension of matrix */
  void replaceVector(const int index,
		       const int numReplace, const double * newElements)
      {matrix_->replaceVector(index,numReplace,newElements);};
  /** Returns a new matrix in reverse order without gaps */
  virtual ClpMatrixBase * reverseOrderedCopy() const;
  /** Returns number of elements in basis
      column is basic if entry >=0 */
  virtual OsiBigIndex numberInBasis(const int * columnIsBasic) const ;
  /// Fills in basis (Returns number of elements and updates numberBasic)
  virtual OsiBigIndex fillBasis(const ClpSimplex * model,
				const int * columnIsBasic, int & numberBasic,
				int * row, int * column,
				double * element) const ;
  /** Creates scales for column copy (rowCopy in model may be modified)
      returns non-zero if no scaling done */
  virtual int scale(ClpSimplex * model) const ;
  /** Unpacks a column into an OsiIndexedvector
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation */
  virtual void unpack(ClpSimplex * model,OsiIndexedVector * rowArray,
		   int column) const ;
  /** Adds multiple of a column into an OsiIndexedvector
      You can use quickAdd to add to vector */
  virtual void add(const ClpSimplex * model,OsiIndexedVector * rowArray,
		   int column, double multiplier) const ;
   /// Allow any parts of a created OsiPackedMatrix to be deleted
   virtual void releasePackedMatrix() const { };
   //@}

  /**@name Matrix times vector methods */
  //@{
    /** Return <code>y + A * scalar *x</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numColumns()</code>
        @precond <code>y<code> must be of size <code>numRows()</code> */
  virtual void times(double scalar,
		       const double * x, double * y) const;
  /// And for scaling
  virtual void times(double scalar,
		     const double * x, double * y,
		     const double * rowScale, 
		     const double * columnScale) const;
    /** Return <code>y + x * scalar * A</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numRows()</code>
        @precond <code>y<code> must be of size <code>numColumns()</code> */
    virtual void transposeTimes(double scalar,
				const double * x, double * y) const;
  /// And for scaling 
    virtual void transposeTimes(double scalar,
				const double * x, double * y,
				const double * rowScale, 
				const double * columnScale) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const OsiIndexedVector * x,
			      OsiIndexedVector * y,
			      OsiIndexedVector * z) const;
    /** Return <code>x *A in <code>z</code> but
	just for indices in y.
	Squashes small elements and knows about ClpSimplex */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
			      const OsiIndexedVector * x,
			      const OsiIndexedVector * y,
			      OsiIndexedVector * z) const;
  //@}

  /**@name Other */
   //@{
  /// Returns OsiPackedMatrix (non const)
  inline OsiPackedMatrix * matrix() const { return matrix_;};
   //@}


  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpPackedMatrix();
   /** Destructor */
   virtual ~ClpPackedMatrix();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpPackedMatrix(const ClpPackedMatrix&);
   /** The copy constructor from an OsiPackedMatrix. */
   ClpPackedMatrix(const OsiPackedMatrix&);

   ClpPackedMatrix& operator=(const ClpPackedMatrix&);
  /// Clone
  virtual ClpMatrixBase * clone() const ;
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// Data
  OsiPackedMatrix * matrix_;
   //@}
};

#endif
