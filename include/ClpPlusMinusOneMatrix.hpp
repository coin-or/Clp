// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpPlusMinusOneMatrix_H
#define ClpPlusMinusOneMatrix_H


#include "CoinPragma.hpp"

#include "ClpMatrixBase.hpp"

/** This implements a simple +- one matrix as derived from ClpMatrixBase.

*/

class ClpPlusMinusOneMatrix : public ClpMatrixBase {
  
public:
  /**@name Useful methods */
   //@{
   /// Return a complete CoinPackedMatrix
   virtual CoinPackedMatrix * getPackedMatrix() const;
    /** Whether the packed matrix is column major ordered or not. */
  virtual bool isColOrdered() const ;
   /** Number of entries in the packed matrix. */
  virtual  CoinBigIndex getNumElements() const; 
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
  virtual const int * getIndices() const
  { return indices_;};

  virtual const CoinBigIndex * getVectorStarts() const;
   /** The lengths of the major-dimension vectors. */
  virtual const int * getVectorLengths() const;

    /** Delete the columns whose indices are listed in <code>indDel</code>. */
  virtual void deleteCols(const int numDel, const int * indDel);
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
  virtual void deleteRows(const int numDel, const int * indDel);
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
  /** Unpacks a column into an CoinIndexedvector
      Note that model is NOT const.  Bounds and objective could
      be modified if doing column generation */
  virtual void unpack(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column) const ;
  /** Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
  virtual void add(const ClpSimplex * model,CoinIndexedVector * rowArray,
		   int column, double multiplier) const ;
   /// Allow any parts of a created CoinMatrix to be deleted
   virtual void releasePlusMinusOneMatrix() const { };
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
	Squashes small elements and knows about ClpSimplex */
  virtual void transposeTimes(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
    /** Return <code>x * scalar * A + y</code> in <code>z</code>. 
	Can use y as temporary array (will be empty at end)
	Squashes small elements and knows about ClpSimplex.
    This version uses row copy*/
  virtual void transposeTimesByRow(const ClpSimplex * model, double scalar,
			      const CoinIndexedVector * x,
			      CoinIndexedVector * y,
			      CoinIndexedVector * z) const;
    /** Return <code>x *A</code> in <code>z</code> but
	just for indices in y.
	Squashes small elements and knows about ClpSimplex */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
				    const CoinIndexedVector * x,
				    const CoinIndexedVector * y,
				    CoinIndexedVector * z) const;
  //@}

  /**@name Other */
   //@{
  /// Return starts of +1s
  inline int * startPositive() const
  { return startPositive_;};
  /// Return starts of -1s
  inline int * startNegative() const
  { return startNegative_;};
   //@}


  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpPlusMinusOneMatrix();
   /** Destructor */
   virtual ~ClpPlusMinusOneMatrix();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpPlusMinusOneMatrix(const ClpPlusMinusOneMatrix&);
   /** The copy constructor from an CoinPlusMinusOneMatrix. */
   ClpPlusMinusOneMatrix(const CoinPackedMatrix&);

   ClpPlusMinusOneMatrix& operator=(const ClpPlusMinusOneMatrix&);
  /// Clone
  virtual ClpMatrixBase * clone() const ;
  /// pass in copy (object takes ownership)
  void passInCopy(int numberRows, int numberColumns,
		  bool columnOrdered, int * indices,
		  int * startPositive, int * startNegative);
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// For fake CoinPackedMatrix
  mutable double * elements_;
  mutable int * lengths_;
  /// Start of +1's for each
  int * startPositive_;
  /// Start of -1's for each
  int * startNegative_;
  /// Data -1, then +1 rows in pairs (row==-1 if one entry)
  int * indices_;
  /// Number of rows
  int numberRows_;
  /// Number of columns
  int numberColumns_;
  /// True if column ordered
  bool columnOrdered_;
  
   //@}
};

#endif
