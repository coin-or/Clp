// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpLinearObjective_H
#define ClpLinearObjective_H

#include "ClpObjective.hpp"

//#############################################################################

/** Linear Objective Class

*/

class ClpLinearObjective : public ClpObjective {
  
public:
  
  ///@name Stuff
  //@{
  
  /** Returns gradient.  If Linear then solution may be NULL,
      also returns an offset (to be added to current one)
  */
  virtual double * gradient(const double * solution, double & offset);
  /// Resize objective
  virtual void resize(int newNumberColumns) ; 
  /// Delete columns in  objective
  virtual void deleteSome(int numberToDelete, const int * which) ; 
  
  //@}
  
  
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  ClpLinearObjective(); 
  
  /// Constructor from objective
  ClpLinearObjective(const double * objective, int numberColumns); 
  
  /// Copy constructor 
  ClpLinearObjective(const ClpLinearObjective &);
  
  /// Assignment operator 
  ClpLinearObjective & operator=(const ClpLinearObjective& rhs);
  
  /// Destructor 
  virtual ~ClpLinearObjective ();

  /// Clone
  virtual ClpObjective * clone() const;
 
  //@}

  //---------------------------------------------------------------------------
  
private:
  ///@name Private member data 
  /// Objective
  double * objective_;
  /// number of columns
  int numberColumns_;
  //@}
};

#endif
