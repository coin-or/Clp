// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpBuild_H
#define ClpBuild_H


#include "CoinPragma.hpp"
#include "CoinFinite.hpp"


/** 
    In many cases it is natural to build a model by adding one row at a time.  In Clp this
    is inefficient so this class gives some help.  An instance of ClpBuild can be built up
    more efficiently and then added to the ClpModel in one go.

    It may be more efficient to have fewer arrays and re-allocate them but this should
    give a large gain over addRow.

    I may extend it to columns if asked.

*/

class ClpBuild {
  
public:
  /**@name Useful methods */
   //@{
   /// add a row
   void addRow(int numberInRow, const int * columns,
	       const double * elements, double rowLower=-COIN_DBL_MAX, 
              double rowUpper=COIN_DBL_MAX);
   /// Return number of rows
  inline int numberRows() const
  { return numberRows_;};
   /// Return maximum number of columns found so far
  inline int numberColumns() const
  { return numberColumns_;};
   /// Return number of elements
  inline CoinBigIndex numberElements() const
  { return numberElements_;};
  /**  Returns number of elements in a row and information in row
   */
  int row(int whichRow, double & rowLower, double & rowUpper,
          int * & indices, double * & elements);
  /**  Returns number of elements in current row and information in row
       Used as rows may be stored in a chain
   */
  int currentRow(double & rowLower, double & rowUpper,
          int * & indices, double * & elements);
  /// Set current row
  void setCurrentRow(int whichRow);
  /// Returns current row number
  int currentRow() const;
   //@}


  /**@name Constructors, destructor */
   //@{
   /** Default constructor. */
   ClpBuild();
   /** Destructor */
   ~ClpBuild();
   //@}

   /**@name Copy method */
   //@{
   /** The copy constructor. */
   ClpBuild(const ClpBuild&);
  /// =
   ClpBuild& operator=(const ClpBuild&);
   //@}
   
    
private:
  /**@name Data members */
   //@{
  /// Current number of rows
  int numberRows_;
  /// Current number of Columns (i.e. max)
  int numberColumns_;
  /// Current number of elements
  CoinBigIndex numberElements_;
  /// Current row pointer
  double * currentRow_;
  /// First row pointer
  double * firstRow_;
  /// Last row pointer
  double * lastRow_;
  /// Type of build (just 0 now for row)
  int type_;
   //@}
};

#endif
