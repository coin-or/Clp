// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPackedMatrix_H
#define ClpPackedMatrix_H


#include "CoinPragma.hpp"

#include "ClpMatrixBase.hpp"

/** This implements OsiPackedMatrix as derived from ClpMatrixBase.

    It adds a few methods that know about model as well as matrix

    For details see OsiPackedMatrix */

class ClpPackedMatrix : public ClpMatrixBase {
  
public:
  /**@name Useful methods */
   //@{
   /// Return a complete CoinPackedMatrix
  virtual CoinPackedMatrix * getPackedMatrix() const { return matrix_;};
    /** Whether the packed matrix is column major ordered or not. */
    virtual bool isColOrdered() const { return matrix_->isColOrdered(); }
   /** Number of entries in the packed matrix. */
  virtual  CoinBigIndex getNumElements() const 
  { return matrix_->getNumElements(); }
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

   virtual const CoinBigIndex * getVectorStarts() const 
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
  /// Append Columns
  virtual void appendCols(int number, const CoinPackedVectorBase * const * columns)
  { matrix_->appendCols(number,columns);};
  /// Append Rows
  virtual void appendRows(int number, const CoinPackedVectorBase * const * rows)
  { matrix_->appendRows(number,rows);};
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
  virtual CoinBigIndex numberInBasis(const int * columnIsBasic) const ;
  /// Fills in basis (Returns number of elements and updates numberBasic)
  virtual CoinBigIndex fillBasis(const ClpSimplex * model,
				const int * columnIsBasic, int & numberBasic,
				int * row, int * column,
				double * element) const ;
  /** If element NULL returns number of elements in column part of basis,
      If not NULL fills in as well */
  virtual CoinBigIndex fillBasis(const ClpSimplex * model,
				 const int * whichColumn, 
				 int numberRowBasic,
				 int numberColumnBasic,
				 int * row, int * column,
				 double * element) const ;
  /** Creates scales for column copy (rowCopy in model may be modified)
      returns non-zero if no scaling done */
  virtual int scale(ClpSimplex * model) const ;
  /** Checks if all elements are in valid range.  Can just
      return true if you are not paranoid.  For Clp I will
      probably expect no zeros.  Code can modify matrix to get rid of
      small elements.
  */
  virtual bool allElementsInRange(ClpModel * model,
				  double smallest, double largest);

  /** Unpacks a column into an CoinIndexedvector
   */
  virtual void unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column) const ;
  /** Unpacks a column into an CoinIndexedvector
   ** in packed foramt
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation (just for this variable) */
  virtual void unpackPacked(ClpSimplex * model,
			    CoinIndexedVector * rowArray,
			    int column) const;
  /** Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
  virtual void add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column, double multiplier) const ;
   /// Allow any parts of a created CoinPackedMatrix to be deleted
   virtual void releasePackedMatrix() const { };
  /** Given positive integer weights for each row fills in sum of weights
      for each column (and slack).
      Returns weights vector
  */
  virtual CoinBigIndex * dubiousWeights(const ClpSimplex * model,int * inputWeights) const;
   //@}

  /**@name Matrix times vector methods */
  //@{
    /** Return <code>y + A * scalar *x</code> in <code>y</code>.
        @pre <code>x</code> must be of size <code>numColumns()</code>
        @pre <code>y</code> must be of size <code>numRows()</code> */
  virtual void times(double scalar,
		       const double * x, double * y) const;
  /// And for scaling
  virtual void times(double scalar,
		     const double * x, double * y,
		     const double * rowScale, 
		     const double * columnScale) const;
    /** Return <code>y + x * scalar * A</code> in <code>y</code>.
        @pre <code>x</code> must be of size <code>numRows()</code>
        @pre <code>y</code> must be of size <code>numColumns()</code> */
    virtual void transposeTimes(double scalar,
				const double * x, double * y) const;
  /// And for scaling 
    virtual void transposeTimes(double scalar,
				const double * x, double * y,
				const double * rowScale, 
				const double * columnScale) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Note - If x packed mode - then z packed mode
	Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Note - If x packed mode - then z packed mode
	Squashes small elements and knows about ClpSimplex.
    This version uses row copy*/
  virtual void transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
    /** Return <code>x *A</code> in <code>z</code> but
	just for indices in y.
	Note - If x packed mode - then z packed mode
	Squashes small elements and knows about ClpSimplex */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
				    const CoinIndexedVector * x,
				    const CoinIndexedVector * y,
				    CoinIndexedVector * z) const;
  //@}

  /**@name Other */
   //@{
  /// Returns CoinPackedMatrix (non const)
  inline CoinPackedMatrix * matrix() const { return matrix_;};
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
   /** The copy constructor from an CoinPackedMatrix. */
   ClpPackedMatrix(const CoinPackedMatrix&);
  /** Subset constructor (without gaps).  Duplicates are allowed
      and order is as given */
  ClpPackedMatrix (const ClpPackedMatrix & wholeModel,
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns);
  ClpPackedMatrix (const CoinPackedMatrix & wholeModel,
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns);

  /** This takes over ownership (for space reasons) */
   ClpPackedMatrix(CoinPackedMatrix * matrix);

   ClpPackedMatrix& operator=(const ClpPackedMatrix&);
  /// Clone
  virtual ClpMatrixBase * clone() const ;
  /** Subset clone (without gaps).  Duplicates are allowed
      and order is as given */
  virtual ClpMatrixBase * subsetClone (
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns) const ;
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// Data
  CoinPackedMatrix * matrix_;
  /// Zero element flag - set true if any zero elements
  bool zeroElements_;
   //@}
};

#endif
