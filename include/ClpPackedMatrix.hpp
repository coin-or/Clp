// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPackedMatrix_H
#define ClpPackedMatrix_H


#include "CoinPragma.hpp"

#include "ClpMatrixBase.hpp"

/** This implements CoinPackedMatrix as derived from ClpMatrixBase.

    It adds a few methods that know about model as well as matrix

    For details see CoinPackedMatrix */

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
  virtual void deleteCols(const int numDel, const int * indDel);
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
  virtual void deleteRows(const int numDel, const int * indDel);
  /// Append Columns
  virtual void appendCols(int number, const CoinPackedVectorBase * const * columns);
  /// Append Rows
  virtual void appendRows(int number, const CoinPackedVectorBase * const * rows);
  /** Replace the elements of a vector.  The indices remain the same.
      This is only needed if scaling and a row copy is used.
      At most the number specified will be replaced.
      The index is between 0 and major dimension of matrix */
  virtual void replaceVector(const int index,
		       const int numReplace, const double * newElements)
      {matrix_->replaceVector(index,numReplace,newElements);};
  /** Returns a new matrix in reverse order without gaps */
  virtual ClpMatrixBase * reverseOrderedCopy() const;
  /// Returns number of elements in column part of basis 
  virtual CoinBigIndex countBasis(ClpSimplex * model,
				 const int * whichColumn, 
				 int numberRowBasic,
				  int & numberColumnBasic);
  /// Fills in column part of basis
  virtual void fillBasis(ClpSimplex * model,
				 const int * whichColumn, 
				 int & numberColumnBasic,
				 int * row, int * start,
				 int * rowCount, int * columnCount,
				 double * element);
  /** Creates scales for column copy (rowCopy in model may be modified)
      returns non-zero if no scaling done */
  virtual int scale(ClpModel * model) const ;
  /** Scales rowCopy if column copy scaled
      Only called if scales already exist */
  virtual void scaleRowCopy(ClpModel * model) const ;
  /** Realy really scales column copy 
      Only called if scales already exist.
      Up to user ro delete */
  virtual ClpMatrixBase * scaledColumnCopy(ClpModel * model) const ;
  /** Checks if all elements are in valid range.  Can just
      return true if you are not paranoid.  For Clp I will
      probably expect no zeros.  Code can modify matrix to get rid of
      small elements.
      check bits (can be turned off to save time) :
      1 - check if matrix has gaps
      2 - check if zero elements
      4 - check and compress duplicates
      8 - report on large and small
  */
  virtual bool allElementsInRange(ClpModel * model,
				  double smallest, double largest,
				  int check=15);
  /** Returns largest and smallest elements of both signs.
      Largest refers to largest absolute value.
  */
  virtual void rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive);

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
  /** Adds multiple of a column into an array */
  virtual void add(const ClpSimplex * model,double * array,
		   int column, double multiplier) const;
   /// Allow any parts of a created CoinPackedMatrix to be deleted
   virtual void releasePackedMatrix() const { };
  /** Given positive integer weights for each row fills in sum of weights
      for each column (and slack).
      Returns weights vector
  */
  virtual CoinBigIndex * dubiousWeights(const ClpSimplex * model,int * inputWeights) const;
  /// Says whether it can do partial pricing
  virtual bool canDoPartialPricing() const;
  /// Partial pricing 
  virtual void partialPricing(ClpSimplex * model, double start, double end,
		      int & bestSequence, int & numberWanted);
  /// makes sure active columns correct
  virtual int refresh(ClpSimplex * model);
  // Really scale matrix
  virtual void reallyScale(const double * rowScale, const double * columnScale);
  /** Set the dimensions of the matrix. In effect, append new empty
      columns/rows to the matrix. A negative number for either dimension
      means that that dimension doesn't change. Otherwise the new dimensions
      MUST be at least as large as the current ones otherwise an exception
      is thrown. */
  virtual void setDimensions(int numrows, int numcols) throw(CoinError);
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
				const double * columnScale,
				double * spare=NULL) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Note - If x packed mode - then z packed mode
	Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Note - If x packed mode - then z packed mode
	This does by column and knows no gaps
	Squashes small elements and knows about ClpSimplex */
  void transposeTimesByColumn(const ClpSimplex * model, double scalar,
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
	Note - z always packed mode */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
				    const CoinIndexedVector * x,
				    const CoinIndexedVector * y,
				    CoinIndexedVector * z) const;
  /** Returns true if can combine transposeTimes and subsetTransposeTimes
      and if it would be faster */
  virtual bool canCombine(const ClpSimplex * model,
                          const CoinIndexedVector * pi) const;
  /// Updates two arrays for steepest 
  virtual void transposeTimes2(const ClpSimplex * model,
                               const CoinIndexedVector * pi1, CoinIndexedVector * dj1,
                               const CoinIndexedVector * pi2, CoinIndexedVector * dj2,
                               CoinIndexedVector * spare,
                               double referenceIn, double devex,
                               // Array for exact devex to say what is in reference framework
                               unsigned int * reference,
                               double * weights, double scaleFactor);
  /// Updates second array for steepest and does devex weights 
  virtual void subsetTimes2(const ClpSimplex * model,
                                CoinIndexedVector * dj1,
                               const CoinIndexedVector * pi2, CoinIndexedVector * dj2,
                               double referenceIn, double devex,
                               // Array for exact devex to say what is in reference framework
                               unsigned int * reference,
                               double * weights, double scaleFactor);
  /// Sets up an effective RHS
  void useEffectiveRhs(ClpSimplex * model);
//@}

  /**@name Other */
   //@{
  /// Returns CoinPackedMatrix (non const)
  inline CoinPackedMatrix * matrix() const { return matrix_;};
  /** Just sets matrix_ to NULL so it can be used elsewhere.
      used in GUB
  */
  inline void setMatrixNull()
  { matrix_=NULL;};
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
  /// number of active columns (normally same as number of columns)
  int numberActiveColumns_;
  /// Zero element flag - set true if any zero elements
  bool zeroElements_;
  /// Gaps flag - set true if column start and length don't say contiguous
  bool hasGaps_;
   //@}
};

#endif
