// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.



#include <cmath>
#include <cassert>
#include <cfloat>
#include <string>
#include <cstdio>
#include <iostream>


#include "CoinPragma.hpp"

#include "CoinHelperFunctions.hpp"

#include "ClpBuild.hpp"

/*
  Format of each row is a bit sleazy.
  First we have pointer to next row
  Then we have two ints giving row number and number of elements
  Then we have two double for lower and upper
  Then we have elements
  Then indices
*/
struct buildFormat {
  buildFormat * next;
  int rowNumber;
  int numberElements;
  double lower;
  double upper;
  double restDouble[1]; 
  int restInt[1]; // just to make correct size 
} ;

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpBuild::ClpBuild () 
  : numberRows_(0),
    numberColumns_(0),
    numberElements_(0),
    currentRow_(NULL),
    firstRow_(NULL),
    lastRow_(NULL),
    type_(0)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpBuild::ClpBuild (const ClpBuild & rhs) 
  : numberRows_(rhs.numberRows_),
    numberColumns_(rhs.numberColumns_),
    numberElements_(rhs.numberElements_),
    type_(rhs.type_)
{
  if (numberRows_) {
    firstRow_=NULL;
    buildFormat * lastRow = NULL;
    buildFormat * currentRow = (buildFormat *) rhs.firstRow_;
    for (int iRow=0;iRow<numberRows_;iRow++) {
      buildFormat * row = currentRow;
      assert (row);
      int numberElements = row->numberElements;
      int length = sizeof(buildFormat)+(numberElements-1)*(sizeof(double)+sizeof(int));
      int doubles = (length + sizeof(double)-1)/sizeof(double);
      double * copyOfRow = new double [doubles];
      memcpy(copyOfRow,row,length);
      if (!firstRow_) {
        firstRow_ = copyOfRow;
      } else {
        // update pointer
        lastRow->next = (buildFormat *) copyOfRow;
      }
      currentRow = currentRow->next; // on to next
      lastRow = (buildFormat *) copyOfRow;
    }
    currentRow_=firstRow_;
    lastRow_=(double *) lastRow;
  } else {
    currentRow_=NULL;
    firstRow_=NULL;
    lastRow_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpBuild::~ClpBuild ()
{
  buildFormat * row = (buildFormat *) firstRow_;
  for (int iRow=0;iRow<numberRows_;iRow++) {
    double * array = (double *) row;
    row = row->next;
    delete [] array;
  }
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpBuild &
ClpBuild::operator=(const ClpBuild& rhs)
{
  if (this != &rhs) {
    buildFormat * row = (buildFormat *) firstRow_;
    for (int iRow=0;iRow<numberRows_;iRow++) {
      double * array = (double *) row;
      row = row->next;
      delete [] array;
    }
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    numberElements_=rhs.numberElements_;
    type_=rhs.type_;
    if (numberRows_) {
      firstRow_=NULL;
      buildFormat * lastRow = NULL;
      buildFormat * currentRow = (buildFormat *) rhs.firstRow_;
      for (int iRow=0;iRow<numberRows_;iRow++) {
        buildFormat * row = currentRow;
        assert (row);
        int numberElements = row->numberElements;
        int length = sizeof(buildFormat)+(numberElements-1)*(sizeof(double)+sizeof(int));
        int doubles = (length + sizeof(double)-1)/sizeof(double);
        double * copyOfRow = new double [doubles];
        memcpy(copyOfRow,row,length);
        if (!firstRow_) {
          firstRow_ = copyOfRow;
        } else {
          // update pointer
          lastRow->next = (buildFormat *) copyOfRow;
        }
        currentRow = currentRow->next; // on to next
        lastRow = (buildFormat *) copyOfRow;
      }
      currentRow_=firstRow_;
      lastRow_=(double *) lastRow;
    } else {
      currentRow_=NULL;
      firstRow_=NULL;
      lastRow_=NULL;
    }
  }
  return *this;
}
// add a row
void 
ClpBuild::addRow(int numberInRow, const int * columns,
                 const double * elements, double rowLower, 
                 double rowUpper)
{
  buildFormat * lastRow = (buildFormat *) lastRow_;
  int length = sizeof(buildFormat)+(numberInRow-1)*(sizeof(double)+sizeof(int));
  int doubles = (length + sizeof(double)-1)/sizeof(double);
  double * newRow = new double [doubles];
  if (!firstRow_) {
    firstRow_ = newRow;
  } else {
    // update pointer
    lastRow->next = (buildFormat *) newRow;
  }
  lastRow_=newRow;
  currentRow_=newRow;
  // now fill in
  buildFormat * row = (buildFormat *) newRow;
  double * els = &row->restDouble[0];
  int * cols = (int *) (els+numberInRow);
  row->next=NULL;
  row->rowNumber=numberRows_;
  numberRows_++;
  row->numberElements=numberInRow;
  numberElements_ += numberInRow;
  row->lower=rowLower;
  row->upper=rowUpper;
  for (int k=0;k<numberInRow;k++) {
    int iColumn = columns[k];
    assert (iColumn>=0);
    numberColumns_ = CoinMax(numberColumns_,iColumn+1);
    els[k]=elements[k];
    cols[k]=iColumn;
  }
}
/*  Returns number of elements in a row and information in row
 */
int 
ClpBuild::row(int whichRow, double & rowLower, double & rowUpper,
              int * & indices, double * & elements)
{
  setCurrentRow(whichRow);
  return currentRow(rowLower,rowUpper,indices,elements);
}
/*  Returns number of elements in current row and information in row
    Used as rows may be stored in a chain
*/
int 
ClpBuild::currentRow(double & rowLower, double & rowUpper,
                     int * & indices, double * & elements)
{
  buildFormat * row = (buildFormat *) currentRow_;
  if (row) {
    int numberElements = row->numberElements;
    elements = &row->restDouble[0];
    indices = (int *) (elements+numberElements);
    rowLower = row->lower;
    rowUpper=row->upper;
    return numberElements;
  } else {
    return -1;
  }
}
// Set current row
void 
ClpBuild::setCurrentRow(int whichRow)
{
  if (whichRow>=0&&whichRow<numberRows_) {
    int nSkip = whichRow-1;
    buildFormat * row = (buildFormat *) firstRow_;
    // if further on then we can start from where we are
    buildFormat * current = (buildFormat *) currentRow_;
    if (current->rowNumber<=whichRow) {
      row=current;
      nSkip = whichRow-current->rowNumber;
    }
    for (int iRow=0;iRow<nSkip;iRow++) {
      row = row->next;
    }
    assert (whichRow==row->rowNumber);
    currentRow_ = (double *) row;
  }
    
}
// Returns current row number
int 
ClpBuild::currentRow() const
{
  buildFormat * row = (buildFormat *) currentRow_;
  if (row)
    return row->rowNumber;
  else
    return -1;
}
