// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Tomlin (with some help from John Forrest)

 */
#ifndef ClpInterior_H
#define ClpInterior_H

#include <iostream>
#include <cfloat>
#include "ClpModel.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpSolve.hpp"
class ClpDualRowPivot;
class ClpPrimalColumnPivot;
class ClpFactorization;
class CoinIndexedVector;
class ClpNonLinearCost;
class ClpInteriorProgress;
// ******** DATA to be moved into protected section of ClpInterior
typedef struct{
  double  atolmin;
  double  r3norm;
  double  LSdamp;
  double* deltay;
} Info;

typedef struct{
  double  atolold;
  double  atolnew;
  double  r3ratio;
  int   istop;
  int   itncg;
} Outfo;
  
typedef struct{
double  gamma;
double  delta;
int MaxIter;
double  FeaTol;
double  OptTol;
double  StepTol;
double  x0min;
double  z0min;
double  mu0;
int   LSmethod;   // 1=Cholesky    2=QR    3=LSQR
int   LSproblem;  // See below
int LSQRMaxIter;
double  LSQRatol1; // Initial  atol
double  LSQRatol2; // Smallest atol (unless atol1 is smaller)
double  LSQRconlim;
int  wait;
} Options;
class Lsqr;
// ***** END
/** This solves LPs using interior point methods

    It inherits from ClpModel and all its arrays are created at
    algorithm time.

*/

class ClpInterior : public ClpModel {
   friend void ClpInteriorUnitTest(const std::string & mpsDir,
				  const std::string & netlibDir);

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    ClpInterior (  );

  /// Copy constructor. 
  ClpInterior(const ClpInterior &);
  /// Copy constructor from model. 
  ClpInterior(const ClpModel &);
  /** Subproblem constructor.  A subset of whole model is created from the 
      row and column lists given.  The new order is given by list order and
      duplicates are allowed.  Name and integer information can be dropped
  */
  ClpInterior (const ClpModel * wholeModel,
	      int numberRows, const int * whichRows,
	      int numberColumns, const int * whichColumns,
	      bool dropNames=true, bool dropIntegers=true);
  /// Assignment operator. This copies the data
    ClpInterior & operator=(const ClpInterior & rhs);
  /// Destructor
   ~ClpInterior (  );
  // Ones below are just ClpModel with some changes
  /** Loads a problem (the constraints on the
        rows are given by lower and upper bounds). If a pointer is 0 then the
        following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>
    */
  void loadProblem (  const ClpMatrixBase& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  void loadProblem (  const CoinPackedMatrix& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);

  /** Just like the other loadProblem() method except that the matrix is
	given in a standard column major ordered format (without gaps). */
  void loadProblem (  const int numcols, const int numrows,
		     const CoinBigIndex* start, const int* index,
		     const double* value,
		     const double* collb, const double* colub,   
		     const double* obj,
		      const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  /// This one is for after presolve to save memory
  void loadProblem (  const int numcols, const int numrows,
		     const CoinBigIndex* start, const int* index,
		      const double* value,const int * length,
		     const double* collb, const double* colub,   
		     const double* obj,
		      const double* rowlb, const double* rowub,
		      const double * rowObjective=NULL);
  /// Read an mps file from the given filename
  int readMps(const char *filename,
	      bool keepNames=false,
	      bool ignoreErrors = false);
  //@}

  /**@name Functions most useful to user */
  //@{
  /** Pdco algorithm - see ClpPdco.hpp for method */
  int pdco();
  // ** Temporary version
  int  pdco( Lsqr *lsqr, Options &options, Info &info, Outfo &outfo);
  //@}

  /**@name most useful gets and sets */
  //@{ 
  /// If problem is primal feasible
  inline bool primalFeasible() const
         { return (sumPrimalInfeasibilities_<=1.0e-5);};
  /// If problem is dual feasible
  inline bool dualFeasible() const
         { return (sumDualInfeasibilities_<=1.0e-5);};
  /// Current (or last) algorithm
  inline int algorithm() const 
  {return algorithm_; } ;
  /// Set algorithm
  inline void setAlgorithm(int value)
  {algorithm_=value; } ;
  /// Sum of dual infeasibilities
  inline double sumDualInfeasibilities() const 
          { return sumDualInfeasibilities_;} ;
  /// Sum of primal infeasibilities
  inline double sumPrimalInfeasibilities() const 
          { return sumPrimalInfeasibilities_;} ;
  //@}

  /******************** End of most useful part **************/
  /**@name Functions less likely to be useful to casual user */
  //@{
  //@}
  /**@name Matrix times vector methods 
     They can be faster if scalar is +- 1
     These are covers so user need not worry about scaling
     Also for simplex I am not using basic/non-basic split */
  //@{
    /** Return <code>y + A * x * scalar</code> in <code>y</code>.
        @pre <code>x</code> must be of size <code>numColumns()</code>
        @pre <code>y</code> must be of size <code>numRows()</code> */
   void times(double scalar,
		       const double * x, double * y) const;
    /** Return <code>y + x * scalar * A</code> in <code>y</code>.
        @pre <code>x</code> must be of size <code>numRows()</code>
        @pre <code>y</code> must be of size <code>numColumns()</code> */
    void transposeTimes(double scalar,
				const double * x, double * y) const ;
  //@}

  /**@name most useful gets and sets */
  //@{ 
  /// Largest error on Ax-b
  inline double largestPrimalError() const
          { return largestPrimalError_;} ;
  /// Largest error on basic duals
  inline double largestDualError() const
          { return largestDualError_;} ;
  //@}

  protected:
  /**@name protected methods */
  //@{
  /// Does most of deletion
  void gutsOfDelete();
  /// Does most of copying
  void gutsOfCopy(const ClpInterior & rhs);
  /// Returns true if data looks okay, false if not
  bool createWorkingData();
  void deleteWorkingData();
  /// Sanity check on input rim data
  bool sanityCheck();
  ///  This does housekeeping
  int housekeeping();
  //@}
  public:
  /**@name public methods */
  //@{
  /// Raw objective value (so always minimize)
  inline double rawObjectiveValue() const
  { return objectiveValue_;};
  /// Returns 1 if sequence indicates column
  inline int isColumn(int sequence) const
  { return sequence<numberColumns_ ? 1 : 0;};
  /// Returns sequence number within section
  inline int sequenceWithin(int sequence) const
  { return sequence<numberColumns_ ? sequence : sequence-numberColumns_;};
  //@}

////////////////// data //////////////////
protected:

  /**@name data.  Many arrays have a row part and a column part.
   There is a single array with both - columns then rows and
   then normally two arrays pointing to rows and columns.  The
   single array is the owner of memory 
  */
  //@{
  /// Largest error on Ax-b
  double largestPrimalError_;
  /// Largest error on basic duals
  double largestDualError_;
  /// Sum of dual infeasibilities
  double sumDualInfeasibilities_;
  /// Sum of primal infeasibilities
  double sumPrimalInfeasibilities_;
  /// Working copy of lower bounds (Owner of arrays below)
  double * lower_;
  /// Row lower bounds - working copy
  double * rowLowerWork_;
  /// Column lower bounds - working copy
  double * columnLowerWork_;
  /// Working copy of upper bounds (Owner of arrays below)
  double * upper_;
  /// Row upper bounds - working copy
  double * rowUpperWork_;
  /// Column upper bounds - working copy
  double * columnUpperWork_;
  /// Working copy of objective 
  double * cost_;
  /// Which algorithm being used
  int algorithm_;
  //@}
};
//#############################################################################
/** A function that tests the methods in the ClpInterior class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging.

    It also does some testing of ClpFactorization class
 */
void
ClpInteriorUnitTest(const std::string & mpsDir,
		   const std::string & netlibDir);


#endif
