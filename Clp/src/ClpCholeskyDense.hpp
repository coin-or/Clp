// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpCholeskyDense_H
#define ClpCholeskyDense_H

#include "ClpCholeskyBase.hpp"
class ClpMatrixBase;

/** Dense class for Clp Cholesky factorization

*/
class ClpCholeskyDense : public ClpCholeskyBase {
  
public:
   /**@name Virtual methods that the derived classes provides  */
   //@{
  /** Orders rows and saves pointer to matrix.and model.
   Returns non-zero if not enough memory */
  virtual int order(ClpInterior * model) ;
  /** Does Symbolic factorization given permutation.
      This is called immediately after order.  If user provides this then
      user must provide factorize and solve.  Otherwise the default factorization is used
      returns non-zero if not enough memory */
  virtual int symbolic();
  /** Factorize - filling in rowsDropped and returning number dropped.
      If return code negative then out of memory */
  virtual int factorize(const double * diagonal, int * rowsDropped) ;
  /** Uses factorization to solve. */
  virtual void solve (double * region) ;
  //@}

   /**@name Non virtual methods for ClpCholeskyDense  */
   //@{
  /** Reserves space.
      If factor not NULL then just uses passed space
   Returns non-zero if not enough memory */
   int reserveSpace(const ClpCholeskyBase * factor, int numberRows) ;
  /** Returns space needed */
   CoinBigIndex space( int numberRows) const;
  /** part 2 of Factorize - filling in rowsDropped */
  void factorizePart2(int * rowsDropped) ;
  /** part 2 of Factorize - filling in rowsDropped - blocked */
  void factorizePart3(int * rowsDropped) ;
  /// Non leaf recursive factor
  void factor(longDouble * a, int n, int numberBlocks,
	      longDouble * diagonal, longDouble * work, int * rowsDropped);
  /// Non leaf recursive triangle rectangle update
  void triRec(longDouble * aTri, int nThis, longDouble * aUnder, longDouble * diagonal, longDouble * work,
	      int nLeft, int iBlock, int jBlock,
	      int numberBlocks);
  /// Non leaf recursive rectangle triangle update
  void recTri(longDouble * aUnder, int nTri, int nDo,
	      int iBlock, int jBlock,longDouble * aTri,
	      longDouble * diagonal, longDouble * work, 
	      int numberBlocks);
  /** Non leaf recursive rectangle rectangle update,
      nUnder is number of rows in iBlock,
      nUnderK is number of rows in kBlock
  */
  void recRec(longDouble * above, int nUnder, int nUnderK,
	      int nDo, longDouble * aUnder, longDouble *aOther,
	      longDouble * work,
	      int iBlock, int jBlock,
	      int numberBlocks);
  /// Leaf recursive factor
  void factorLeaf(longDouble * a, int n, 
	      longDouble * diagonal, longDouble * work, int * rowsDropped);
  /// Leaf recursive triangle rectangle update
  void triRecLeaf(longDouble * aTri, longDouble * aUnder,
		  longDouble * diagonal, longDouble * work,
		  int nUnder);
  /// Leaf recursive rectangle triangle update
  void recTriLeaf(longDouble * aUnder, longDouble * aTri, 
		  longDouble * diagonal, longDouble * work, int nUnder);
  /** Leaf recursive rectangle rectangle update,
      nUnder is number of rows in iBlock,
      nUnderK is number of rows in kBlock
  */
  void recRecLeaf(const longDouble * COIN_RESTRICT above, 
		  const longDouble * COIN_RESTRICT aUnder, 
		  longDouble * COIN_RESTRICT aOther, 
		  const longDouble * COIN_RESTRICT work,
		  int nUnder);
  /// Forward part of solve
  void solveF1(longDouble * a,int n,double * region);
  void solveF2(longDouble * a,int n,double * region,double * region2);
  /// Backward part of solve
  void solveB1(longDouble * a,int n,double * region);
  void solveB2(longDouble * a,int n,double * region,double * region2);
  /** Uses factorization to solve. */
  void solveLong (longDouble * region) ;
  /// Forward part of solve
  void solveF1Long(longDouble * a,int n,longDouble * region);
  void solveF2Long(longDouble * a,int n,longDouble * region,longDouble * region2);
  /// Backward part of solve
  void solveB1Long(longDouble * a,int n,longDouble * region);
  void solveB2Long(longDouble * a,int n,longDouble * region,longDouble * region2);
  /** Uses factorization to solve. */
  void solveLongWork (longWork * region) ;
  /// Forward part of solve
  void solveF1LongWork(longDouble * a,int n,longWork * region);
  void solveF2LongWork(longDouble * a,int n,longWork * region,longWork * region2);
  /// Backward part of solve
  void solveB1LongWork(longDouble * a,int n,longWork * region);
  void solveB2LongWork(longDouble * a,int n,longWork * region,longWork * region2);
  int bNumber(const longDouble * array,int &, int&);
  /// A
  inline longDouble * aMatrix() const
  { return sparseFactor_;}
  /// Diagonal
  inline longDouble * diagonal() const
  { return diagonal_;}
  //@}


  /**@name Constructors, destructor */
  //@{
  /** Default constructor. */
  ClpCholeskyDense();
  /** Destructor  */
  virtual ~ClpCholeskyDense();
  // Copy
  ClpCholeskyDense(const ClpCholeskyDense&);
  // Assignment
  ClpCholeskyDense& operator=(const ClpCholeskyDense&);
  /// Clone
  virtual ClpCholeskyBase * clone() const ;
  //@}
   
    
private:
  /**@name Data members */
  //@{
  /// Just borrowing space
  bool borrowSpace_;
  //@}
};

#endif
