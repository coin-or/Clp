/* $Id: iis.cpp 2278 2017-10-02 09:51:14Z forrest $ */
/*
  Copyright (C) 2003, International Business Machines Corporation and others.
  All Rights Reserved.

  This sample program is designed to illustrate programming techniques
  using CoinLP, has not been thoroughly tested and comes without any
  warranty whatsoever.

  You may copy, modify and distribute this sample program without any
  restrictions whatsoever and without any payment to anyone.
*/
/* Find an irreducible infeasible subsytem.
   This is not trying to be fast - this implementation is 
   really trying to find the maximum feasible subsystem 
*/
#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include <iomanip>

int main(int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  if (argc < 2) {
    printf("please give a model\n");
    status=1;
  } else {
    status = model.readMps(argv[1], true);
  }
  if (status)
    exit(10);
  // check infeasible
  model.initialSolve();

  if (model.problemStatus()==0) {
    printf("Problem is feasible - nothing to do\n");
    exit(11);
  }
  double * rowLower = model.rowLower();
  double * rowUpper = model.rowUpper();
  int numberColumns = model.numberColumns();
  int numberRows = model.numberRows();
  if (model.numberPrimalInfeasibilities()==1) {
    const double * rowActivity = model.primalRowSolution();
    for (int iRow = 0; iRow < numberRows; iRow++) {
      if (rowActivity[iRow]<rowLower[iRow]-1.0e-5||
	  rowActivity[iRow]>rowUpper[iRow]+1.0e-5) {
	printf("The following one row is the IIS\n");
	printf("%d %s\n",iRow,model.getRowName(iRow).c_str());
	return 0;
      }
    }
  }
  // This example only looks at constraints
  // add infeasibility slacks
  ClpSimplex model2(model);
  double time1 = CoinCpuTime();

  int nAdded=0;

  CoinBigIndex * starts = new CoinBigIndex [2*numberRows+1];
  int * rows = new int [2*numberRows+1];
  double * elements = new double [2*numberRows];
  double * newObjective = new double [2*numberRows];
  memset(model2.objective(),0,numberColumns*sizeof(double));
  for (int iRow = 0; iRow < numberRows; iRow++) {
    if (rowLower[iRow]>-1.0e30) {
      rows[nAdded]=iRow;
      elements[nAdded++]=1.0;
    }
    if (rowUpper[iRow]<1.0e30) {
      rows[nAdded]=iRow;
      elements[nAdded++]=-1.0;
    }
  }
  rows[nAdded]=-1; // so can test if two slacks for a row
  starts[0]=0;
#define OBJ 10000.0
  for (int i=0;i<nAdded;i++) {
    starts[i+1]=i+1;
    newObjective[i]=OBJ;
  }
  model2.addColumns(nAdded,NULL,NULL,newObjective,starts,rows,elements);
  delete [] newObjective;
  delete [] starts;
  delete [] elements;
  int numberColumns2=numberColumns+nAdded;
  rowLower = model2.rowLower();
  rowUpper = model2.rowUpper();
  // solve
  model2.allSlackBasis();
  model2.dual();
  printf("Initial sum of infeasibilities is %g\n",model2.objectiveValue()/OBJ);
  if (model2.numberPrimalInfeasibilities()) {
    printf("ouch\n");
    model2.writeMps("bad.mps");
    abort();
  }
  model2.setLogLevel(0);
  model2.setInfeasibilityCost(1.0e15);
  //#define SWITCH_OFF_SCALING
#ifdef SWITCH_OFF_SCALING
  model2.scaling(0);
#endif
  int * coverRows = new int [numberRows];
  int * candidateRows = new int [numberRows];
  int * nextRows = new int [numberRows];
  int nCover=0;
  int nCandidate=0;
  int outRow=9999;
  const double * duals = model2.dualRowSolution();
  for (int iRow=0;iRow<numberRows;iRow++) {
    if (fabs(duals[iRow])>1.0e-5)
      nextRows[nCandidate++]=iRow;
  }
  assert (nCandidate);
  /* save basis info -
     we could reduce problem size each time but
     normally not many passes needed */
  unsigned char * statusArray = new unsigned char [numberColumns2+numberRows];
  double * solution = new double[numberColumns2];
  memcpy(statusArray,model2.statusArray(),numberColumns2+numberRows);
  memcpy(solution,model2.primalColumnSolution(),numberColumns2*sizeof(double));
  double lastObjectiveValue = model2.objectiveValue();
  double * objective = model2.objective();
  int nSolves=0;
  int nPass=0;
  while (outRow>=0) {
    memcpy(candidateRows,nextRows,nCandidate*sizeof(int));
    double sumInf=COIN_DBL_MAX;
    int nInf=numberColumns2;
    bool badAccuracy=false;
    outRow=-1;
    nPass++;
    for (int j=0;j<nCandidate;j++) {
      memcpy(model2.statusArray(),statusArray,numberColumns2+numberRows);
      memcpy(model2.primalColumnSolution(),solution,numberColumns2*sizeof(double));
      int iRow=candidateRows[j];
      /* One can free rows or zero out costs -
	 zeroing out costs is cleaner from basis */
#define ZAP_COSTS
#ifdef ZAP_COSTS
      // zero cost and solve
      int iColumn0=-1;
      int iColumn1=-1;
      for (int i=0;i<nAdded;i++) {
	if (rows[i]==iRow) {
	  iColumn0=i+numberColumns;
	  objective[iColumn0]=0.0;
	  if (rows[i+1]==iRow) {
	    iColumn1=i+1+numberColumns;
	    objective[iColumn1]=0.0;
	  }
	  break;
	}
      }
      model2.primal(1);
      if (model2.objectiveValue()>lastObjectiveValue+1.0e-1) {
	if (!badAccuracy)
	  printf("Sum infeasibilities increased from %g to %g! - problems with accuracy\n",
		 lastObjectiveValue/OBJ,model2.objectiveValue()/OBJ);
	badAccuracy=true;
      }
      nSolves++;
      objective[iColumn0]=OBJ;
      if (iColumn1>=0)
	objective[iColumn1]=OBJ;
#else
      // delete and solve
      double lower=rowLower[iRow];
      double upper=rowUpper[iRow];
      rowLower[iRow]=-COIN_DBL_MAX;
      rowUpper[iRow]=COIN_DBL_MAX;
      model2.primal(1);
      if (model2.objectiveValue()>lastObjectiveValue+1.0e-1) {
	if (!badAccuracy)
	  printf("Sum infeasibilities increased from %g to %g! - problems with accuracy\n",
		 lastObjectiveValue/OBJ,model2.objectiveValue()/OBJ);
	badAccuracy=true;
      }
      nSolves++;
      rowLower[iRow]=lower;
      rowUpper[iRow]=upper;
#endif
      double objectiveValue = model2.objectiveValue();
      if (objectiveValue<sumInf+1.0e-3) {
	int n=0;
	for (int i=numberColumns;i<numberColumns2;i++) {
	  if (solution[i]>1.0e-5)
	    n++;
	}
	if (objectiveValue>1.0e-6) {
	  if (objectiveValue<sumInf-1.0e-3||n<nInf) {
	    sumInf=objectiveValue;
	    outRow=iRow;
	    nInf=n;
	  }
	} else {
	  coverRows[nCover++]=iRow;
	  printf("Pass %d (%d solves, %.2f seconds) - candidate %d (%s) added - Final as feasible\n",nPass,nSolves,CoinCpuTime()-time1,iRow,
		 model2.getRowName(iRow).c_str());
	  outRow=-1;
	  break; // finished
	}
      }
    }
    if (outRow>=0) {
#ifdef ZAP_COSTS
      int iColumn0=-1;
      int iColumn1=-1;
      for (int i=0;i<nAdded;i++) {
	if (rows[i]==outRow) {
	  iColumn0=i+numberColumns;
	  objective[iColumn0]=0.0;
	  if (rows[i+1]==outRow) {
	    iColumn1=i+1+numberColumns;
	    objective[iColumn1]=0.0;
	  }
	  break;
	}
      }
#else
      rowLower[outRow]=-COIN_DBL_MAX;
      rowUpper[outRow]=COIN_DBL_MAX;
#endif
      coverRows[nCover++]=outRow;
      memcpy(model2.statusArray(),statusArray,numberColumns2+numberRows);
      memcpy(model2.primalColumnSolution(),solution,numberColumns2*sizeof(double));
      model2.primal(1);
      if (model2.objectiveValue()>lastObjectiveValue+1.0e-1) {
	if (!badAccuracy)
	  printf("Sum infeasibilities increased from %g to %g on cover solve! - problems with accuracy\n",
		 lastObjectiveValue/OBJ,model2.objectiveValue()/OBJ);
	badAccuracy=true;
      }
      lastObjectiveValue=model2.objectiveValue();
      memcpy(statusArray,model2.statusArray(),numberColumns2+numberRows);
      memcpy(solution,model2.primalColumnSolution(),numberColumns2*sizeof(double));
      nCandidate=0;
      for (int iRow=0;iRow<numberRows;iRow++) {
	if (fabs(duals[iRow])>1.0e-6)
	  nextRows[nCandidate++]=iRow;
      }
      printf("Pass %d (%d solves, %.2f seconds) - candidate %d (%s) added - infeasibility %g (%d)\n",
	     nPass,nSolves,CoinCpuTime()-time1,outRow,
	     model2.getRowName(outRow).c_str(),
	     model2.objectiveValue()/OBJ,nInf);
    }
  }
  model2=model;
  model2.deleteRows(nCover,coverRows);
  // make sure not unbounded
  memset(model2.objective(),0,numberColumns*sizeof(double));
  model2.dual();
  printf("The following %d rows cover the IIS\n",nCover);
  for (int i=0;i<nCover;i++) {
    int iRow=coverRows[i];
    printf("%d %s\n",iRow,model.getRowName(iRow).c_str());
  }
  if (model2.problemStatus()) {
    printf("We seem to have an accuracy problem??\n");
  } else {
    // see if we can do better
    CoinPackedMatrix * matrix = model.matrix();
    // get row copy
    CoinPackedMatrix rowCopy = *matrix;
    rowCopy.reverseOrdering();
    const int * column = rowCopy.getIndices();
    const int * rowLength = rowCopy.getVectorLengths();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    const double * element = rowCopy.getElements();
    numberRows=model2.numberRows();
    memcpy(statusArray,model2.statusArray(),numberColumns+numberRows);
    memcpy(solution,model2.primalColumnSolution(),numberColumns*sizeof(double));
    int lastRow=numberRows-1;
    for (int i=0;i<nCover;i++) {
      memcpy(model2.statusArray(),statusArray,numberColumns+numberRows);
      memcpy(model2.primalColumnSolution(),solution,numberColumns*sizeof(double));
      int iRow=coverRows[i];
      model2.addRow(rowLength[iRow],
		    column+rowStart[iRow],element+rowStart[iRow],
		    model.rowLower()[iRow],model.rowUpper()[iRow]);
      model2.dual();
      if (!model2.problemStatus()) {
	printf("%d %s should not be in cover\n",iRow,model.getRowName(iRow).c_str());
      }
      model2.deleteRows(1,&lastRow);
    }
  }
  delete [] rows;
  delete [] statusArray;
  delete [] solution;
  delete [] candidateRows;
  delete [] nextRows;
  delete [] coverRows;
  return 0;
}
