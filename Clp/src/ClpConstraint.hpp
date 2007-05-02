// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpConstraint_H
#define ClpConstraint_H


//#############################################################################
class ClpSimplex;
class ClpModel;

/** Constraint Abstract Base Class

Abstract Base Class for describing a constraint or objective function

*/
class ClpConstraint  {
  
public:
  
  ///@name Stuff
  //@{
  
  /** Fills gradient.  If Linear then solution may be NULL,
      also returns value of function
      If refresh is false then uses last solution
      Uses model for scaling
      Returns non-zero if gradient udefined at current solution
  */
  virtual int gradient(const ClpSimplex * model,
		       const double * solution,
		       double * gradient,
		       double & functionValue ,
		       bool refresh=true) const =0;
  /// Resize constraint
  virtual void resize(int newNumberColumns) = 0; 
  /// Delete columns in  constraint
  virtual void deleteSome(int numberToDelete, const int * which) = 0; 
  /// Scale constraint 
  virtual void reallyScale(const double * columnScale) =0;
  /** Given a zeroed array sets nonlinear columns to 1.
      Returns number of nonlinear columns
   */
  virtual int markNonlinear(char * which) const = 0;
  //@}
  
  
  ///@name Constructors and destructors
  //@{
  /// Default Constructor
  ClpConstraint(); 
  
  /// Copy constructor 
  ClpConstraint(const ClpConstraint &);
  
  /// Assignment operator 
  ClpConstraint & operator=(const ClpConstraint& rhs);
  
  /// Destructor 
  virtual ~ClpConstraint ();

  /// Clone
  virtual ClpConstraint * clone() const = 0;
 
  //@}

  ///@name Other
  //@{
  /// Returns type, 0 linear, 1 nonlinear
  inline int type()
  { return type_;};
  /// Row number (-1 is objective)
  inline int rowNumber() const
  {return rowNumber_;};
  
  /// Constraint function value
  inline double functionValue () const
  { return functionValue_;};
  //@}

  //---------------------------------------------------------------------------
  
protected:
  ///@name Protected member data 
  //@{
  /// Gradient at last evaluation
  mutable double * lastGradient_;
  /// Value of non-linear part of constraint
  mutable double functionValue_;
  /// Type of constraint - linear is 1
  int type_;
  /// Row number (-1 is objective)
  int rowNumber_;
  //@}
};

#endif
