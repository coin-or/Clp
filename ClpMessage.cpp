// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpMessage.hpp"

typedef struct {
  CLP_Message internalNumber;
  int externalNumber; // or continuation
  char detail;
  const char * message;
} Clp_message;
static Clp_message us_english[]=
{
  {CLP_SIMPLEX_FINISHED,0,1,"Optimal - objective value %g"},
  {CLP_SIMPLEX_INFEASIBLE,1,1,"Primal infeasible - objective value %g"},
  {CLP_SIMPLEX_UNBOUNDED,2,1,"Dual infeasible - objective value %g"},
  {CLP_SIMPLEX_STOPPED,3,1,"Stopped - objective value %g"},
  {CLP_SIMPLEX_ERROR,4,1,"Stopped due to errors - objective value %g"},
  //{CLP_SIMPLEX_STATUS,5,1,"%d  Objective %g%? Primal infeas %g (%d)%? Dual infeas %g (%d)%? without free dual infeas (%d)"},
  {CLP_SIMPLEX_STATUS,5,1,"%d  Obj %g%? Primal inf %g (%d)%? Dual inf %g (%d)%? w.o. free dual inf (%d)"},
  {CLP_DUAL_BOUNDS,25,3,"Looking optimal checking bounds with %g"},
  {CLP_SIMPLEX_ACCURACY,6,3,"Primal error %g, dual error %g"},
  {CLP_SIMPLEX_BADFACTOR,7,2,"Singular factorization of basis - status %d"},
  {CLP_SIMPLEX_BOUNDTIGHTEN,8,3,"Bounds were tightened %d times"},
  {CLP_SIMPLEX_INFEASIBILITIES,9,1,"%d infeasibilities"},
  {CLP_SIMPLEX_FLAG,10,3,"Flagging variable %c%d"},
  {CLP_SIMPLEX_GIVINGUP,11,2,"Stopping as close enough"},
  {CLP_DUAL_CHECKB,12,2,"New dual bound of %g"},
  {CLP_DUAL_ORIGINAL,13,3,"Going back to original objective"},
  {CLP_SIMPLEX_PERTURB,14,1,"Perturbing problem by %g %% of %g - largest nonzero change %g (%% %g) - largest zero change %g"},
  {CLP_PRIMAL_ORIGINAL,15,2,"Going back to original tolerance"},
  {CLP_PRIMAL_WEIGHT,16,2,"New infeasibility weight of %g"},
  {CLP_PRIMAL_OPTIMAL,17,2,"Looking optimal with tolerance of %g"},
  {CLP_SINGULARITIES,18,2,"%d total structurals rejected in initial factorization"},
  {CLP_MODIFIEDBOUNDS,19,1,"%d variables/rows fixed as scaled bounds too close"},
  {CLP_RIMSTATISTICS1,20,2,"Absolute values of scaled objective range from %g to %g"},
  {CLP_RIMSTATISTICS2,21,2,"Absolute values of scaled bounds range from %g to %g, minimum gap %g"},
  {CLP_RIMSTATISTICS3,22,2,"Absolute values of scaled rhs range from %g to %g, minimum gap %g"},
  {CLP_POSSIBLELOOP,23,2,"Possible loop - %d matches (%x) after %d checks"},
  {CLP_SMALLELEMENTS,24,1,"Matrix will be packed to eliminate %d small elements"},
  {CLP_SIMPLEX_HOUSE1,101,32,"dirOut %d, dirIn %d, theta %g, out %g, dj %g, alpha %g"},
  {CLP_SIMPLEX_HOUSE2,102,4,"%d %g In: %c%d Out: %c%d%? dj ratio %g distance %g%? dj %g distance %g"},
  {CLP_SIMPLEX_NONLINEAR,103,4,"Primal nonlinear change %g (%d)"},
  {CLP_SIMPLEX_FREEIN,104,32,"Free column in %d"},
  {CLP_SIMPLEX_PIVOTROW,105,32,"Pivot row %d"},
  {CLP_DUAL_CHECK,106,4,"Btran alpha %g, ftran alpha %g"},
  {CLP_PRIMAL_DJ,107,4,"Btran dj %g, ftran dj %g"},
  {CLP_PACKEDSCALE_INITIAL,1001,2,"Initial range of elements is %g to %g"},
  {CLP_PACKEDSCALE_WHILE,1002,3,"Range of elements is %g to %g"},
  {CLP_PACKEDSCALE_FINAL,1003,2,"Final range of elements is %g to %g"},
  {CLP_PACKEDSCALE_FORGET,1004,2,"Not bothering to scale as good enough"},
  {CLP_INITIALIZE_STEEP,1005,3,"Initializing steepest edge weights - old %g, new %g"},
  {CLP_UNABLE_OPEN,6001,0,"Unable to open file %s for reading"},
  {CLP_BAD_BOUNDS,6002,0,"%d bad bound pairs or bad objectives were found - first at %c%d"},
  {CLP_BAD_MATRIX,6003,0,"Matrix has %d large values, first at column %d, row %d is %g"},
  {CLP_LOOP,6004,0,"Can't get out of loop - stopping (Debug)"},
  {CLP_DUPLICATEELEMENTS,26,1,"Matrix will be packed to eliminate %d duplicate elements"},
  {CLP_IMPORT_RESULT,27,1,"Model was imported from %s in %g seconds"},
  {CLP_IMPORT_ERRORS,3001,1," There were %d errors when importing model from %s"},
  {CLP_EMPTY_PROBLEM,3002,0,"Not solving empty problem - %d rows, %d columns and %d elements"},
  {CLP_CRASH,28,1,"Crash put %d variables in basis, %d dual infeasibilities"},
  {CLP_END_VALUES_PASS,29,1,"End of values pass after %d iterations"},
  {CLP_QUADRATIC_BOTH,108,32,"%s %d (%g) and %d (%g) both basic"},
  {CLP_QUADRATIC_PRIMAL_DETAILS,109,32,"coeff %g, %g, %g - dj %g - deriv zero at %g, sj at %g"},
  {CLP_IDIOT_ITERATION,30,1,"%d infeas %g, obj %g - mu %g, its %d, %d interior"}, 
  {CLP_INFEASIBLE,3003,1,"Analysis indicates model infeasible or unbounded"}, 
  {CLP_MATRIX_CHANGE,31,2,"Matrix can not be converted into %s"}, 
  {CLP_TIMING,32,1,"%s objective %.10g - %d iterations time %.2f2%?, Presolve %.2f%?, Idiot %.2f%?"}, 
  {CLP_INTERVAL_TIMING,33,2,"%s took %.2f seconds (total %.2f)"}, 
  {CLP_SPRINT,34,1,"Pass %d took %d iterations, objective %g, dual infeasibilities %g( %d)"},
  {CLP_DUMMY_END,999999,0,""}
};
static Clp_message uk_english[]=
{
  {CLP_SIMPLEX_FINISHED,0,1,"Optimal - objective value %g,\
 okay CLP can solve some LPs but you really need Xpress from Dash Associates :-)"},
  {CLP_DUMMY_END,999999,0,""}
};
/* Constructor */
ClpMessage::ClpMessage(Language language) :
  CoinMessages(sizeof(us_english)/sizeof(Clp_message))
{
  language_=language;
  strcpy(source_,"Clp");
  Clp_message * message = us_english;

  while (message->internalNumber!=CLP_DUMMY_END) {
     CoinOneMessage oneMessage(message->externalNumber,message->detail,
			       message->message);
     addMessage(message->internalNumber,oneMessage);
     message ++;
}

  // now override any language ones

  switch (language) {
  case uk_en:
    message = uk_english;
    break;

  default:
    message=NULL;
    break;
  }

  // replace if any found
  if (message) {
    while (message->internalNumber!=CLP_DUMMY_END) {
      replaceMessage(message->internalNumber,message->message);
      message ++;
    }
  }
}
