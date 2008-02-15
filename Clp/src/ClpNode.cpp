// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpSimplex.hpp"
#include "ClpNode.hpp"
#include "ClpFactorization.hpp"
#include "ClpDualRowSteepest.hpp"

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNode::ClpNode () :
  branchingValue_(0.5),
  factorization_(NULL),
  weights_(NULL),
  status_(NULL),
  primalSolution_(NULL),
  dualSolution_(NULL),
  pivotVariables_(NULL),
  fixed_(NULL),
  sequence_(1),
  numberFixed_(0)
{
 branchState_.firstBranch=0;
 branchState_.branch=0;
}
//-------------------------------------------------------------------
// Useful Constructor from model
//-------------------------------------------------------------------
ClpNode::ClpNode (const ClpSimplex * model, const ClpNodeStuff * stuff) :
  branchingValue_(0.5),
  factorization_(NULL),
  weights_(NULL),
  status_(NULL),
  primalSolution_(NULL),
  dualSolution_(NULL),
  pivotVariables_(NULL),
  fixed_(NULL),
  sequence_(1),
  numberFixed_(0)
{
  branchState_.firstBranch=0;
  branchState_.branch=0;
  gutsOfConstructor(model,stuff);
}

//-------------------------------------------------------------------
// Most of work of constructor from model
//-------------------------------------------------------------------
void
ClpNode::gutsOfConstructor (const ClpSimplex * model, const ClpNodeStuff * stuff) 
{
  // save stuff
  factorization_ = new ClpFactorization(*model->factorization());
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  int numberTotal = numberRows+numberColumns;
  status_ = CoinCopyOfArray(model->statusArray(),numberTotal);
  primalSolution_ = CoinCopyOfArray(model->solutionRegion(),numberTotal);
  dualSolution_ = CoinCopyOfArray(model->djRegion(),numberTotal); //? has duals as well?
  pivotVariables_ = CoinCopyOfArray(model->pivotVariable(),numberRows); 
  ClpDualRowSteepest* pivot =
    dynamic_cast< ClpDualRowSteepest*>(model->dualRowPivot());
  if (pivot)
    weights_ = new ClpDualRowSteepest(*pivot);
  const double * lower = model->columnLower();
  const double * upper = model->columnUpper();
  const double * solution = model->primalColumnSolution();
  const char * integerType = model->integerInformation();
  int iColumn;
  sequence_=-1;
  double integerTolerance = stuff->integerTolerance_;
  double mostAway=integerTolerance;
  int numberAway=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (integerType[iColumn]) {
      double value = solution[iColumn];
      value = max(value,(double) lower[iColumn]);
      value = min(value,(double) upper[iColumn]);
      double nearest = floor(value+0.5);
      if (fabs(value-nearest)>integerTolerance)
	numberAway++;
      if (fabs(value-nearest)>mostAway) {
	mostAway=fabs(value-nearest);
	sequence_=iColumn;
	branchingValue_=value;
	branchState_.branch=0;
	if (value<=nearest)
	  branchState_.firstBranch=1; // up
	else
	  branchState_.firstBranch=0; // down
      }
    }
  }
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNode::ClpNode (const ClpNode & source) 
{  
  printf("ClpNode copy not implemented\n");
  abort();
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNode::~ClpNode ()
{
  delete factorization_;
  delete weights_;
  delete [] status_;
  delete [] primalSolution_;
  delete [] dualSolution_;
  delete [] pivotVariables_;
  delete [] fixed_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNode &
ClpNode::operator=(const ClpNode& rhs)
{
  if (this != &rhs) {
    printf("ClpNode = not implemented\n");
    abort();
  }
  return *this;
}
// Applies node to model
void 
ClpNode::applyNode(ClpSimplex * model, bool justBounds )
{
  // current bound
  int way=branchState_.firstBranch;
  if (branchState_.branch>0)
    way=1-way;
  if (!way) {
    // This should also do underlying internal bound
    model->setColumnUpper(sequence_,floor(branchingValue_));
  } else {
    // This should also do underlying internal bound
    model->setColumnLower(sequence_,ceil(branchingValue_));
  }
  const double * lower = model->columnLower();
  const double * upper = model->columnUpper();
  // apply dj fixings
  for (int i=0;i<numberFixed_;i++) {
    int iColumn = fixed_[i];
    if ((iColumn&0x10000000)!=0) {
      iColumn &= 0xfffffff;
      model->setColumnLower(iColumn,upper[iColumn]);
    } else {
	model->setColumnUpper(iColumn,lower[iColumn]);
    }
  }
  if (!justBounds) {
    model->setFactorization(*factorization_);
    ClpDualRowSteepest* pivot =
      dynamic_cast< ClpDualRowSteepest*>(model->dualRowPivot());
    if (pivot)
      *pivot=*weights_; // may be better to copy stuff
    int numberRows = model->numberRows();
    int numberColumns = model->numberColumns();
    int numberTotal = numberRows+numberColumns;
    CoinMemcpyN(status_,numberTotal,model->statusArray());
    CoinMemcpyN(primalSolution_,numberTotal,model->solutionRegion());
    CoinMemcpyN(dualSolution_,numberTotal,model->djRegion()); //? has duals as well?
    CoinMemcpyN(pivotVariables_,numberRows,model->pivotVariable());
  }
}
// Fix on reduced costs
int 
ClpNode::fixOnReducedCosts(ClpSimplex * model)
{
  return 0;
}
/* Way for integer variable -1 down , +1 up */
int 
ClpNode::way() const
{
  int way=branchState_.firstBranch;
  if (branchState_.branch>0)
    way=1-way;
  return way ? -1 : +1;
}
// Return true if branch exhausted
bool 
ClpNode::fathomed() const
{ 
  return branchState_.branch>=1
;
}
// Change state of variable i.e. go other way
void 
ClpNode::changeState()
{
  branchState_.branch++;
  assert (branchState_.branch<=2);
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpNodeStuff::ClpNodeStuff () :
  integerTolerance_(1.0e-7),
  integerIncrement_(1.0e-8)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpNodeStuff::ClpNodeStuff (const ClpNodeStuff & source) 
{  
  printf("ClpNodeStuff copy not implemented\n");
  abort();
}
//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpNodeStuff &
ClpNodeStuff::operator=(const ClpNodeStuff& rhs)
{
  if (this != &rhs) {
    printf("ClpNodeStuff = not implemented\n");
    abort();
  }
  return *this;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpNodeStuff::~ClpNodeStuff ()
{
}
