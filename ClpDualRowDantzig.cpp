// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "ClpSimplex.hpp"
#include "ClpDualRowDantzig.hpp"
#include "OsiIndexedVector.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpDualRowDantzig::ClpDualRowDantzig () 
: ClpDualRowPivot()
{
  type_=1;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpDualRowDantzig::ClpDualRowDantzig (const ClpDualRowDantzig & source) 
: ClpDualRowPivot(source)
{  

}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpDualRowDantzig::~ClpDualRowDantzig ()
{

}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpDualRowDantzig &
ClpDualRowDantzig::operator=(const ClpDualRowDantzig& rhs)
{
  if (this != &rhs) {
    ClpDualRowPivot::operator=(rhs);
  }
  return *this;
}

// Returns pivot row, -1 if none
int 
ClpDualRowDantzig::pivotRow()
{
  assert(model_);
  int iRow;
  const int * pivotVariable = model_->pivotVariable();
  double largest=model_->currentPrimalTolerance();
  // we can't really trust infeasibilities if there is primal error
  if (model_->largestPrimalError()>1.0e-8)
    largest *= model_->largestPrimalError()/1.0e-8;
  int chosenRow=-1;
  for (iRow=0;iRow<model_->numberRows();iRow++) {
    int iPivot=pivotVariable[iRow];
    double value = model_->solution(iPivot);
    double lower = model_->lower(iPivot);
    double upper = model_->upper(iPivot);
    if (max(value-upper,lower-value)>largest) {
      int iSequence = pivotVariable[iRow];
      if (!model_->flagged(iSequence)) {
	chosenRow=iRow;
	largest=max(value-upper,lower-value);
      }
    }
  }
  return chosenRow;
}
  
/* Updates primal solution (and maybe list of candidates)
   Uses input vector which it deletes
   Computes change in objective function
*/
void 
ClpDualRowDantzig::updatePrimalSolution(
					OsiIndexedVector * primalUpdate,
					double primalRatio,
					double & objectiveChange)
{
  double * work = primalUpdate->denseVector();
  int number = primalUpdate->getNumElements();
  int * which = primalUpdate->getIndices();
  int i;
  double changeObj=0.0;
  const int * pivotVariable = model_->pivotVariable();
  for (i=0;i<number;i++) {
    int iRow=which[i];
    int iPivot=pivotVariable[iRow];
    double & value = model_->solutionAddress(iPivot);
    double cost = model_->cost(iPivot);
    double change = primalRatio*work[iRow];
    value -= change;
    changeObj -= change*cost;
    work[iRow]=0.0;
  }
  primalUpdate->setNumElements(0);
  objectiveChange += changeObj;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpDualRowPivot * ClpDualRowDantzig::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpDualRowDantzig(*this);
  } else {
    return new ClpDualRowDantzig();
  }
}

