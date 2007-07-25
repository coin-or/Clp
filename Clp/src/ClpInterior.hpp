// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Tomlin (pdco)
   John Forrest (standard predictor-corrector)

   Note JJF has added arrays - this takes more memory but makes
   flow easier to understand and hopefully easier to extend

 */
#ifndef ClpInterior_H
#define ClpInterior_H

#include <iostream>
#include <cfloat>
#include "ClpModel.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpSolve.hpp"
class ClpLsqr;
class ClpPdcoBase;
/// ******** DATA to be moved into protected section of ClpInterior
typedef struct{
  double  atolmin;
  double  r3norm;
  double  LSdamp;
  double* deltay;
} Info;
/// ******** DATA to be moved into protected section of ClpInterior

typedef struct{
  double  atolold;
  double  atolnew;
  double  r3ratio;
  int   istop;
  int   itncg;
} Outfo;
/// ******** DATA to be moved into protected section of ClpInterior
  
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
class ClpCholeskyBase;
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
  /** Borrow model.  This is so we dont have to copy large amounts
      of data around.  It assumes a derived class wants to overwrite
      an empty model with a real one - while it does an algorithm.
      This is same as ClpModel one. */
  void borrowModel(ClpModel & otherModel);
  /** Return model - updates any scalars */
  void returnModel(ClpModel & otherModel);
  //@}

  /**@name Functions most useful to user */
  //@{
  /** Pdco algorithm - see ClpPdco.hpp for method */
  int pdco();
  // ** Temporary version
  int  pdco( ClpPdcoBase * stuff, Options &options, Info &info, Outfo &outfo);
  /// Primal-Dual Predictor-Corrector barrier
  int primalDual();
  //@}

  /**@name most useful gets and sets */
  //@{ 
  /// If problem is primal feasible
  inline bool primalFeasible() const
         { return (sumPrimalInfeasibilities_<=1.0e-5);}
  /// If problem is dual feasible
  inline bool dualFeasible() const
         { return (sumDualInfeasibilities_<=1.0e-5);}
  /// Current (or last) algorithm
  inline int algorithm() const 
  {return algorithm_; } 
  /// Set algorithm
  inline void setAlgorithm(int value)
  {algorithm_=value; } 
  /// Sum of dual infeasibilities
  inline double sumDualInfeasibilities() const 
          { return sumDualInfeasibilities_;} 
  /// Sum of primal infeasibilities
  inline double sumPrimalInfeasibilities() const 
          { return sumPrimalInfeasibilities_;} 
  /// dualObjective.
  inline double dualObjective() const
  { return dualObjective_;}
  /// primalObjective.
  inline double primalObjective() const
  { return primalObjective_;}
  /// diagonalNorm
  inline double diagonalNorm() const
  { return diagonalNorm_;}
  /// linearPerturbation
  inline double linearPerturbation() const
  { return linearPerturbation_;}
  inline void setLinearPerturbation(double value)
  { linearPerturbation_=value;}
  /// diagonalPerturbation
  inline double diagonalPerturbation() const
  { return diagonalPerturbation_;}
  inline void setDiagonalPerturbation(double value)
  { diagonalPerturbation_=value;}
  /// gamma
  inline double gamma() const
  { return gamma_;}
  inline void setGamma(double value)
  { gamma_=value;}
  /// delta
  inline double delta() const
  { return delta_;}
  inline void setDelta(double value)
  { delta_=value;}
  /// ComplementarityGap
  inline double complementarityGap() const 
          { return complementarityGap_;} 
  //@}

  /**@name most useful gets and sets */
  //@{ 
  /// Largest error on Ax-b
  inline double largestPrimalError() const
          { return largestPrimalError_;} 
  /// Largest error on basic duals
  inline double largestDualError() const
          { return largestDualError_;} 
  /// Maximum iterations
  inline int maximumBarrierIterations() const
  { return maximumBarrierIterations_;}
  inline void setMaximumBarrierIterations(int value)
  { maximumBarrierIterations_=value;}
  /// Set cholesky (and delete present one)
  void setCholesky(ClpCholeskyBase * cholesky);
  /// Return number fixed to see if worth presolving
  int numberFixed() const;
  /** fix variables interior says should be.  If reallyFix false then just
      set values to exact bounds */
  void fixFixed(bool reallyFix=true);
  /// Primal erturbation vector
  inline double * primalR() const
  { return primalR_;}
  /// Dual erturbation vector
  inline double * dualR() const
  { return dualR_;}
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
  { return objectiveValue_;}
  /// Returns 1 if sequence indicates column
  inline int isColumn(int sequence) const
  { return sequence<numberColumns_ ? 1 : 0;}
  /// Returns sequence number within section
  inline int sequenceWithin(int sequence) const
  { return sequence<numberColumns_ ? sequence : sequence-numberColumns_;}
  /// Checks solution
  void checkSolution();
  /** Modifies djs to allow for quadratic.
      returns quadratic offset */
  double quadraticDjs(double * djRegion, const double * solution,
		      double scaleFactor);

  /// To say a variable is fixed
  inline void setFixed( int sequence)
  {
    status_[sequence] |= 1;
  }
  inline void clearFixed( int sequence)
  {
    status_[sequence] &= ~1;
  }
  inline bool fixed(int sequence) const
  {return ((status_[sequence]&1)!=0);}

  /// To flag a variable
  inline void setFlagged( int sequence)
  {
    status_[sequence] |= 2;
  }
  inline void clearFlagged( int sequence)
  {
    status_[sequence] &= ~2;
  }
  inline bool flagged(int sequence) const
  {return ((status_[sequence]&2)!=0);}

  /// To say a variable is fixed OR free
  inline void setFixedOrFree( int sequence)
  {
    status_[sequence] |= 4;
  }
  inline void clearFixedOrFree( int sequence)
  {
    status_[sequence] &= ~4;
  }
  inline bool fixedOrFree(int sequence) const
  {return ((status_[sequence]&4)!=0);}

  /// To say a variable has lower bound
  inline void setLowerBound( int sequence)
  {
    status_[sequence] |= 8;
  }
  inline void clearLowerBound( int sequence)
  {
    status_[sequence] &= ~8;
  }
  inline bool lowerBound(int sequence) const
  {return ((status_[sequence]&8)!=0);}

  /// To say a variable has upper bound
  inline void setUpperBound( int sequence)
  {
    status_[sequence] |= 16;
  }
  inline void clearUpperBound( int sequence)
  {
    status_[sequence] &= ~16;
  }
  inline bool upperBound(int sequence) const
  {return ((status_[sequence]&16)!=0);}

  /// To say a variable has fake lower bound
  inline void setFakeLower( int sequence)
  {
    status_[sequence] |= 32;
  }
  inline void clearFakeLower( int sequence)
  {
    status_[sequence] &= ~32;
  }
  inline bool fakeLower(int sequence) const
  {return ((status_[sequence]&32)!=0);}

  /// To say a variable has fake upper bound
  inline void setFakeUpper( int sequence)
  {
    status_[sequence] |= 64;
  }
  inline void clearFakeUpper( int sequence)
  {
    status_[sequence] &= ~64;
  }
  inline bool fakeUpper(int sequence) const
  {return ((status_[sequence]&64)!=0);}
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
  /// Worst complementarity
  double worstComplementarity_;
  /// 
public:
  double xsize_;
  double zsize_;
protected:
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
public:
  /// Rhs
  double * rhs_;
  double * x_;
  double * y_;
  double * dj_;
protected:
  /// Pointer to Lsqr object
  ClpLsqr * lsqrObject_;
  /// Pointer to stuff
  ClpPdcoBase * pdcoStuff_;
  /// Below here is standard barrier stuff
  /// mu.
  double mu_;
  /// objectiveNorm.
  double objectiveNorm_;
  /// rhsNorm.
  double rhsNorm_;
  /// solutionNorm.
  double solutionNorm_;
  /// dualObjective.
  double dualObjective_;
  /// primalObjective.
  double primalObjective_;
  /// diagonalNorm.
  double diagonalNorm_;
  /// stepLength
  double stepLength_;
  /// linearPerturbation
  double linearPerturbation_;
  /// diagonalPerturbation
  double diagonalPerturbation_;
  // gamma from Saunders and Tomlin regularized
  double gamma_;
  // delta from Saunders and Tomlin regularized
  double delta_;
  /// targetGap
  double targetGap_;
  /// projectionTolerance
  double projectionTolerance_;
  /// maximumRHSError.  maximum Ax
  double maximumRHSError_;
  /// maximumBoundInfeasibility.
  double maximumBoundInfeasibility_;
  /// maximumDualError.
  double maximumDualError_;
  /// diagonalScaleFactor.
  double diagonalScaleFactor_;
  /// scaleFactor.  For scaling objective
  double scaleFactor_;
  /// actualPrimalStep
  double actualPrimalStep_;
  /// actualDualStep
  double actualDualStep_;
  /// smallestInfeasibility
  double smallestInfeasibility_;
  /// historyInfeasibility.
#define LENGTH_HISTORY 5
  double historyInfeasibility_[LENGTH_HISTORY];
  /// complementarityGap.
  double complementarityGap_;
  /// baseObjectiveNorm
  double baseObjectiveNorm_;
  /// worstDirectionAccuracy
  double worstDirectionAccuracy_;
  /// maximumRHSChange
  double maximumRHSChange_;
  /// errorRegion. i.e. Ax
  double * errorRegion_;
  /// rhsFixRegion.
  double * rhsFixRegion_;
  /// upperSlack
  double * upperSlack_;
  /// lowerSlack
  double * lowerSlack_;
  /// diagonal
  double * diagonal_;
  /// solution
  double * solution_;
  /// work array
  double * workArray_;
  /// delta X
  double * deltaX_;
  /// delta Y
  double * deltaY_;
  /// deltaZ.
  double * deltaZ_;
  /// deltaW.
  double * deltaW_;
  /// deltaS.
  double * deltaSU_;
  double * deltaSL_;
  /// Primal regularization array
  double * primalR_;
  /// Dual regularization array
  double * dualR_;
  /// rhs B
  double * rhsB_;
  /// rhsU.
  double * rhsU_;
  /// rhsL.
  double * rhsL_;
  /// rhsZ.
  double * rhsZ_;
  /// rhsW.
  double * rhsW_;
  /// rhs C
  double * rhsC_;
  /// zVec
  double * zVec_;
  /// wVec
  double * wVec_;
  /// cholesky.
  ClpCholeskyBase * cholesky_;
  /// numberComplementarityPairs i.e. ones with lower and/or upper bounds (not fixed)
  int numberComplementarityPairs_;
  /// numberComplementarityItems_ i.e. number of active bounds
  int numberComplementarityItems_;
  /// Maximum iterations
  int maximumBarrierIterations_;
  /// gonePrimalFeasible.
  bool gonePrimalFeasible_;
  /// goneDualFeasible.
  bool goneDualFeasible_;
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
