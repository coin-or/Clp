// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpSmallMatrix_H
#define ClpSmallMatrix_H


#include "CoinPragma.hpp"

#include "ClpMatrixBase.hpp"
// Change this to what you want
typedef float Double;
typedef unsigned short Int;
//typedef double Double;
//typedef int Int;

/** This implements a small matrix as derived from ClpMatrixBase.
*/


class ClpSmallMatrix : public ClpMatrixBase {
  
public:
  /**@name Useful methods */
   //@{
   /// Return a complete CoinPackedMatrix
   virtual CoinPackedMatrix * getPackedMatrix() const;
    /** Whether the packed matrix is column major ordered or not. */
   virtual bool isColOrdered() const { return true; }
   /** Number of entries in the packed matrix. */
  virtual  CoinBigIndex getNumElements() const 
  { return columnStart_[numberColumns_]; }
   /** Number of columns. */
   virtual int getNumCols() const { return numberColumns_; }
   /** Number of rows. */
  virtual int getNumRows() const { return numberRows_; };

   /** A vector containing the elements in the packed matrix. Note that there
	might be gaps in this list, entries that do not belong to any
	major-dimension vector. To get the actual elements one should look at
	this vector together with vectorStarts and vectorLengths. */
  virtual const double * getElements() const; 
   /** A vector containing the minor indices of the elements in the packed
        matrix. Note that there might be gaps in this list, entries that do not
        belong to any major-dimension vector. To get the actual elements one
        should look at this vector together with vectorStarts and
        vectorLengths. */
  virtual const int * getIndices() const;

  virtual const CoinBigIndex * getVectorStarts() const;
   /** The lengths of the major-dimension vectors. */
  virtual const int * getVectorLengths() const;
  /** The length of a single major-dimension vector. */
  virtual int getVectorLength(int index) const 
  { return columnStart_[index+1]-columnStart_[index];};
  /** Returns largest and smallest elements of both signs.
      Largest refers to largest absolute value.
  */
  virtual void rangeOfElements(double & smallestNegative, double & largestNegative,
		       double & smallestPositive, double & largestPositive);

    /** Delete the columns whose indices are listed in <code>indDel</code>. */
  virtual void deleteCols(const int numDel, const int * indDel);
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
  virtual void deleteRows(const int numDel, const int * indDel);
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
   virtual void releasePackedMatrix() const {};
  /// Returns true if can create row copy
  virtual bool canGetRowCopy() const
  { return true;};
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
	Note - If x packed mode - then z packed mode */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
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
  void transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
  //@}

  /**@name Other */
   //@{
   //@}


  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpSmallMatrix();
  /// Constructor with data (takes ownership)
  ClpSmallMatrix(int numberColumns, int numberRows,
		 int * starts, Int * index, Double * element);
   /** Destructor */
   virtual ~ClpSmallMatrix();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpSmallMatrix(const ClpSmallMatrix &rhs);
   /** The copy constructor from an CoinPackedMatrix. */
   ClpSmallMatrix(const CoinPackedMatrix &rhs);
  /// = operator
   ClpSmallMatrix& operator=(const ClpSmallMatrix & rhs);
  /// Clone
  virtual ClpMatrixBase * clone() const ;
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// Number of rows
  int numberRows_;
  /// Number of columns
  int numberColumns_;
  /// Column starts
  int * columnStart_;
  /// Rows
  Int * row_;
  /// Elements
  Double * element_;
   //@}
};

#endif
