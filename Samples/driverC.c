/* Copyright (C) 2003 International Business Machines
   Corporation and others.  All Rights Reserved. */

/* This example shows the use of the "C" interface */

#include "Clp_C_Interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main (int argc, const char *argv[])
{
  /* Get default model */
  Clp_Simplex  * model = Clp_newModel();
  int status;
  // Keep names when reading an mps file
  if (argc<2)
    status=Clp_readMps(model,"../../Mps/Sample/p0033.mps",1,0);
  else
    status=Clp_readMps(model,argv[1],1,0);

  if (status) {
    fprintf(stderr,"Bad readMps %s\n",argv[1]);
    fprintf(stdout,"Bad readMps %s\n",argv[1]);
    exit(1);
  }

  if (argc<3 ||!strstr(argv[2],"primal")) {
    /* Use the dual algorithm unless user said "primal" */
    Clp_initialDualSolve(model);
  } else {
    Clp_initialPrimalSolve(model);
  }

  char modelName[80];
  Clp_problemName(model,80,modelName);
  printf("Model %s has %d rows and %d columns\n",
	 modelName,Clp_numberRows(model),Clp_numberColumns(model));

  /* remove this to print solution */

  /*exit(0); */

  /*
    Now to print out solution.  The methods used return modifiable
    arrays while the alternative names return const pointers -
    which is of course much more virtuous.

    This version just does non-zero columns
  
   */
  printf("--------------------------------------\n");

  /* Columns */

  int numberColumns = Clp_numberColumns(model);

  /* Alternatively getColSolution(model) */
  double * columnPrimal = Clp_primalColumnSolution(model);
  /* Alternatively getReducedCost(model) */
  double * columnDual = Clp_dualColumnSolution(model);
  /* Alternatively getColLower(model) */
  double * columnLower = Clp_columnLower(model);
  /* Alternatively getColUpper(model) */
  double * columnUpper = Clp_columnUpper(model);
  /* Alternatively getObjCoefficients(model) */
  double * columnObjective = Clp_objective(model);
    
  /* If we have not kept names (parameter to readMps) this will be 0 */
  assert(Clp_lengthNames(model));

  int iColumn;
  
  printf("                       Primal          Dual         Lower         Upper          Cost\n");

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double value;
    value = columnPrimal[iColumn];
    if (value>1.0e-8||value<-1.0e-8) {
      char name[20];
      Clp_columnName(model,iColumn,name);
      printf("%6d %8s",iColumn,name);
      printf(" %13g",columnPrimal[iColumn]);
      printf(" %13g",columnDual[iColumn]);
      printf(" %13g",columnLower[iColumn]);
      printf(" %13g",columnUpper[iColumn]);
      printf(" %13g",columnObjective[iColumn]);
      printf("\n");
    }
  }
  printf("--------------------------------------\n");

  return 0;
}    
