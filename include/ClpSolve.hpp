// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef ClpSolve_H
#define ClpSolve_H

/** 
    This is a very simple class to guide algorithms.  It is used to tidy up 
    passing parameters to initialSolve and maybe for output from that

*/

class ClpSolve  {

public:

  /** enums for solve function */
  enum SolveType {
    useDual=0,
    usePrimal,
    usePrimalorSprint,
    automatic
  };
  enum PresolveType {
    presolveOn=0,
    presolveOff,
    presolveNumber
  };

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    ClpSolve (  );

  /// Copy constructor. 
  ClpSolve(const ClpSolve &);
  /// Assignment operator. This copies the data
    ClpSolve & operator=(const ClpSolve & rhs);
  /// Destructor
   ~ClpSolve (  );
  //@}

  /**@name Functions most useful to user */
  //@{
  /** Special options - bits
0      4 - use crash (default allslack in dual, idiot in primal)
      8 - all slack basis in primal
2      16 - switch off interrupt handling
3      32 - do not try and make plus minus one matrix
      64 - do not use sprint even if problem looks good
   */
  /** which translation is:
      which:
      0 - startup in Dual  (nothing if basis exists).:
                   0 - no basis, 1 crash
      1 - startup in Primal (nothing if basis exists):
                   0 - use initiative
		   1 - use crash
		   2 - use idiot and look at further info
		   3 - use sprint and look at further info
		   4 - use all slack
		   5 - use initiative but no idiot
		   6 - use initiative but no sprint
		   7 - use initiative but no crash
                   8 - do allslack or idiot
                   9 - do allslack or sprint
      2 - interrupt handling - 0 yes, 1 no (for threadsafe)
      3 - whether to make +- 1matrix - 0 yes, 1 no
  */
  void setSpecialOption(int which,int value,int extraInfo=-1);
  int getSpecialOption(int which) const;

  /// Solve types
  void setSolveType(SolveType method, int extraInfo=-1);
  SolveType getSolveType();

  // Presolve types
  void setPresolveType(PresolveType amount, int extraInfo=-1);
  PresolveType getPresolveType();
  int getPresolvePasses() const;
  /// Extra info for idiot (or sprint)
  int getExtraInfo(int which) const;
  //@}

////////////////// data //////////////////
private:

  /**@name data.
  */
  //@{
  /// Solve type
  SolveType method_;
  /// Presolve type
  PresolveType presolveType_;
  /// Amount of presolve
  int numberPasses_;
  /// Options
  int options_[4];
  /// Extra information
  int extraInfo_[4];
  //@}
};
#endif
