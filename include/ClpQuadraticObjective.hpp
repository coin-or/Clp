// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpQuadraticObjective_H
#define ClpQuadraticObjective_H

#include "ClpObjective.hpp"
#include "CoinPackedMatrix.hpp"

//#############################################################################

/** Quadratic Objective Class

*/

class ClpQuadraticObjective : public ClpObjective {
  
public:
  
  ///@name Stuff
  //@{
  
  /** Returns gradient.  If Quadratic then solution may be NULL,
      also returns an offset (to be added to current one)
      If refresh is false then uses last solution
  */
  virtual double * gradient(const double * solution, double & offset,bool refresh);
  /// Resize objective
  virtual void resize(int newNumberColumns) ; 
  /// Delete columns in  objective
  virtual void deleteSome(int numberToDelete, const int * which) ; 
  
  //@}
  
  
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  ClpQuadraticObjective(); 
  
  /// Constructor from objective
  ClpQuadraticObjective(const double * linearObjective, int numberColumns,
			const CoinBigIndex * start,
			const int * column, const double * element,
			int numberExtendedColumns_=-1);
  
  /** Copy constructor .
      If type is -1 then make sure half symmetric,
      if +1 then make sure full
  */
  ClpQuadraticObjective(const ClpQuadraticObjective & rhs,int type=0);
  /** Subset constructor.  Duplicates are allowed
      and order is as given.
  */
  ClpQuadraticObjective (const ClpQuadraticObjective &rhs,int numberColumns, 
				      const int * whichColumns) ;
  
  /// Assignment operator 
  ClpQuadraticObjective & operator=(const ClpQuadraticObjective& rhs);
  
  /// Destructor 
  virtual ~ClpQuadraticObjective ();

  /// Clone
  virtual ClpObjective * clone() const;
  /** Subset clone.  Duplicates are allowed
      and order is as given.
  */
  virtual ClpObjective * subsetClone (int numberColumns, 
				      const int * whichColumns) const;
 
  /** Load up quadratic objective.  This is stored as a CoinPackedMatrix */
  void loadQuadraticObjective(const int numberColumns, 
			      const CoinBigIndex * start,
			      const int * column, const double * element,
			      int numberExtendedColumns=-1);
  void loadQuadraticObjective (  const CoinPackedMatrix& matrix);
  /// Get rid of quadratic objective
  void deleteQuadraticObjective();
  //@}
  ///@name Gets and sets
  //@{
   /// Quadratic objective
   inline CoinPackedMatrix * quadraticObjective() const     { return quadraticObjective_; }
   /// Linear objective
   inline double * linearObjective() const     { return objective_; }
  /// Length of linear objective which could be bigger
  inline int numberExtendedColumns() const
  {return numberExtendedColumns_;};
  /// Number of columns in quadratic objective
  inline int numberColumns() const
  {return numberColumns_;};
  //@}

  //---------------------------------------------------------------------------
  
private:
  ///@name Private member data 
  /// Quadratic objective
  CoinPackedMatrix * quadraticObjective_;
  /// Objective
  double * objective_;
  /// Gradient
  double * gradient_;
  /// Useful to have number of columns about
  int numberColumns_;
  /// Also length of linear objective which could be bigger
  int numberExtendedColumns_;
  //@}
};

#endif
