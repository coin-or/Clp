// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpGubMatrix_H
#define ClpGubMatrix_H


#include "CoinPragma.hpp"

#include "ClpPackedMatrix.hpp"
class ClpSimplex;
/** This implements Gub rows plus a ClpPackedMatrix.

    There will be a version using ClpPlusMinusOne matrix but
    there is no point doing one with ClpNetworkMatrix (although
    an embedded network is attractive).

*/

class ClpGubMatrix : public ClpPackedMatrix {
  
public:
  /** Returns a new matrix in reverse order without gaps */
  //virtual ClpMatrixBase * reverseOrderedCopy() const;
  /** If element NULL returns number of elements in column part of basis,
      If not NULL fills in as well */
  virtual CoinBigIndex fillBasis(const ClpSimplex * model,
				 const int * whichColumn, 
				 int numberRowBasic,
				 int numberColumnBasic,
				 int * row, int * column,
				 double * element) const ;
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
  /// Partial pricing 
  virtual void partialPricing(ClpSimplex * model, int start, int end,
		      int & bestSequence, int & numberWanted);
   //@}

  /**@name Matrix times vector methods */
  //@{
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
	Note - z always packed mode
	Squashes small elements and knows about ClpSimplex */
  virtual void subsetTransposeTimes(const ClpSimplex * model,
				    const CoinIndexedVector * x,
				    const CoinIndexedVector * y,
				    CoinIndexedVector * z) const;
  /** expands an updated column to allow for extra rows which the main
      solver does not know about and returns number added.  If the arrays are NULL 
      then returns number of extra entries needed.

      This active in Gub
  */
  virtual int extendUpdated(CoinIndexedVector * update, double * lower,
			    double * solution, double * upper);
  //@}



  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpGubMatrix();
   /** Destructor */
   virtual ~ClpGubMatrix();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpGubMatrix(const ClpGubMatrix&);
   /** The copy constructor from an CoinPackedMatrix. */
   ClpGubMatrix(const CoinPackedMatrix&);
  /** Subset constructor (without gaps).  Duplicates are allowed
      and order is as given */
  ClpGubMatrix (const ClpGubMatrix & wholeModel,
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns);
  ClpGubMatrix (const CoinPackedMatrix & wholeModel,
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns);

  /** This takes over ownership (for space reasons) */
   ClpGubMatrix(CoinPackedMatrix * matrix);

  /** This takes over ownership (for space reasons) and is the
      real constructor*/
   ClpGubMatrix(ClpPackedMatrix * matrix, int numberSets,
		const int * start, const int * end,
		const double * lower, const double * upper,
		const unsigned char * status=NULL);

   ClpGubMatrix& operator=(const ClpGubMatrix&);
  /// Clone
  virtual ClpMatrixBase * clone() const ;
  /** Subset clone (without gaps).  Duplicates are allowed
      and order is as given */
  virtual ClpMatrixBase * subsetClone (
		    int numberRows, const int * whichRows,
		    int numberColumns, const int * whichColumns) const ;
  //@}
  /**@name gets and sets */
  //@{
  /// Status
  inline ClpSimplex::Status getStatus(int sequence) const
  {return static_cast<ClpSimplex::Status> (status_[sequence]&7);};
  inline void setStatus(int sequence, ClpSimplex::Status status)
  {
    unsigned char & st_byte = status_[sequence];
    st_byte &= ~7;
    st_byte |= status;
  };
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// Starts
  int * start_;
  /// End
  int * end_;
  /// Lower bounds on sets
  double * lower_;
  /// Upper bounds on sets
  double * upper_;
  /// Status of slacks
  mutable unsigned char * status_;
  /// Backward pointer to set number
  int * backward_;
  /// Key variable of set
  mutable int * keyVariable_;
  /** Next basic variable in set - starts at key and end with -(set+1) */
  mutable int * next_;
  /// Number of sets (gub rows)
  int numberSets_;
   //@}
};

#endif
