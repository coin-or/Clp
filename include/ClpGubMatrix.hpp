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
  virtual CoinBigIndex fillBasis(ClpSimplex * model,
				 const int * whichColumn, 
				 int numberRowBasic,
				 int numberColumnBasic,
				 int * row, int * column,
				 double * element)  ;
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
  /**
     mode=0  - Set up before "update" and "times" for primal solution using extended rows
     mode=1  - Cleanup primal solution after "times" using extended rows.
     mode=2  - Check (or report on) primal infeasibilities
  */
  virtual void primalExpanded(ClpSimplex * model,int mode);
  /** 
      mode=0  - Set up before "updateTranspose" and "transposeTimes" for duals using extended
                updates array (and may use other if dual values pass)
      mode=1  - Update dual solution after "transposeTimes" using extended rows.
      mode=2  - Check (or report on) dual infeasibilities
  */
  virtual void dualExpanded(ClpSimplex * model,CoinIndexedVector * array,
			    double * other,int mode);
  /** 
      mode=0  - Create list of non-key basics in pivotVariable_ using
                number as numberBasic in and out
      mode=1  - Set all key variables as basic
  */
  virtual int generalExpanded(ClpSimplex * model,int mode,int & number);
  /// Sets up an effective RHS and does gub crash if needed
  void useEffectiveRhs(ClpSimplex * model,bool cheapest=true);
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
  /// To say key is above ub
  inline void setAbove( int sequence)
  {
    unsigned char iStat = status_[sequence];
    iStat &= ~24;
    status_[sequence] = iStat|16;
  };
  /// To say key is feasible
  inline void setFeasible( int sequence)
  {
    unsigned char iStat = status_[sequence];
    iStat &= ~24;
    status_[sequence] = iStat|8;
  };
  /// To say key is below lb
  inline void setBelow( int sequence)
  {
    unsigned char iStat = status_[sequence];
    iStat &= ~24;
    status_[sequence] = iStat;
  };
  inline double weight( int sequence) const
  {
    int iStat = status_[sequence]&31;
    iStat = iStat>>3;
    return (double) (iStat-1);
  };
   //@}
   
    
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
  /// Sum of dual infeasibilities
  double sumDualInfeasibilities_;
  /// Sum of primal infeasibilities
  double sumPrimalInfeasibilities_;
  /// Sum of Dual infeasibilities using tolerance based on error in duals
  double sumOfRelaxedDualInfeasibilities_;
  /// Sum of Primal infeasibilities using tolerance based on error in primals
  double sumOfRelaxedPrimalInfeasibilities_;
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
  /// Number of dual infeasibilities
  int numberDualInfeasibilities_;
  /// Number of primal infeasibilities
  int numberPrimalInfeasibilities_;
  /// Number of sets (gub rows)
  int numberSets_;
  /// First gub variables (same as start_[0] at present)
  int firstGub_;
  /// last gub variable (same as end_[numberSets_-1] at present)
  int lastGub_;
  /// type of gub - 0 not contiguous, 1 contiguous
  int gubType_;
   //@}
};

#endif
