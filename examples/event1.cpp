
#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplex.hpp"

//#############################################################################

/************************************************************************
This keeps track of iterations using event handler
Coding was thrown together - data should be saved more elagantly
*/
#define SAVE_ITS 3
// Sequence of In variable
int sequenceIn[SAVE_ITS]={-1};
// Direction of In, 1 going up, -1 going down, 0 not a clue
int directionIn[SAVE_ITS]={-1};
// Sequence of Out variable
int sequenceOut[SAVE_ITS]={-1};
// Direction of Out, 1 to upper bound, -1 to lower bound, 0 - superbasic
int directionOut[SAVE_ITS]={-1};
// Pivot Row
int pivot[SAVE_ITS]={-1};
/** This is so user can trap events and do useful stuff.  

    ClpSimplex model_ is available as well as anything else you care 
    to pass in
*/

class MyEventHandler : public ClpEventHandler {

public:
  /**@name Overrides */
  //@{
  virtual int event(Event whichEvent);
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler(ClpSimplex *model);
  /** Destructor */
  virtual ~MyEventHandler();
  /** The copy constructor. */
  MyEventHandler(const MyEventHandler &rhs);
  /// Assignment
  MyEventHandler &operator=(const MyEventHandler &rhs);
  /// Clone
  virtual ClpEventHandler *clone() const;
  //@}

protected:
  // data goes here
};
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler()
  : ClpEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyEventHandler::MyEventHandler(const MyEventHandler &rhs)
  : ClpEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler::MyEventHandler(ClpSimplex *model)
  : ClpEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyEventHandler::~MyEventHandler()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyEventHandler &
MyEventHandler::operator=(const MyEventHandler &rhs)
{
  if (this != &rhs) {
    ClpEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpEventHandler *MyEventHandler::clone() const
{
  return new MyEventHandler(*this);
}

int MyEventHandler::event(Event whichEvent)
{
  if (whichEvent == endOfIteration) {
    // move up
    for (int i=SAVE_ITS-2;i>=0;i--) {
      sequenceIn[i+1] = sequenceIn[i];
      directionIn[i+1] = directionIn[i];
      sequenceOut[i+1] = sequenceOut[i];
      directionOut[i+1] = directionOut[i];
      pivot[i+1] = pivot[i];
    }      
    sequenceIn[0] = model_->sequenceIn();
    directionIn[0] = model_->directionIn();
    sequenceOut[0] = model_->sequenceOut();
    directionOut[0] = model_->directionOut();
    pivot[0] = model_->pivotRow();
  }
  return -1;
}

int main(int argc, const char *argv[])
{
#ifndef OSICLP
  // using ClpSimplex
  ClpSimplex model;
#else
  OsiClpSolverInterface solver1;
#endif
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
  if (argc >= 2) {
    mpsFileName = argv[1];
#ifndef OSICLP
    int numMpsReadErrors = model.readMps(mpsFileName.c_str());
#else
    int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
#endif
    if (numMpsReadErrors != 0) {
      printf("%d errors reading MPS file\n", numMpsReadErrors);
      return numMpsReadErrors;
    }
  } else {
    printf("Need mps file\n");
    return -1;
  }
  // allow Clp to track iterations
  MyEventHandler clpEventHandler;
#ifndef OSICLP
  ClpSimplex * simplex = &model;
#else
  // go over to Clp
  ClpSimplex * simplex = solver1.getModelPtr();
#endif
  simplex->passInEventHandler(&clpEventHandler);
  // If tiny problem more output
  if (simplex->numberRows()<40)
    simplex->setLogLevel(63);
  simplex->primal();
  /* print last few iterations
     then can change basis etc */
  printf("last few iterations\n");
  int numberIterations = simplex->numberIterations();
  for (int i=0;i<SAVE_ITS;i++) {
    printf("Iteration %d in %d direction %d out %d direction %d pivotrow %d\n",
	   numberIterations-i,sequenceIn[i],directionIn[i],
	   sequenceOut[i],directionOut[i],pivot[i]);
  }
  return 0;
}
