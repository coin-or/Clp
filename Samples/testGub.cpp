// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "ClpGubMatrix.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "CoinSort.hpp"
#include "CoinTime.hpp"
int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  int maxIts=0;
  int maxFactor=100;
  if (argc<2)
    status=model.readMps("../../Mps/Sample/p0033.mps");
  else
    status=model.readMps(argv[1]);
  if (status) {
    printf("errors on input\n");
    exit(77);
  }
  if (argc>2) {
    maxFactor = atoi(argv[2]);
    printf("max factor %d\n",maxFactor);
  }
  if (argc>3) {
    maxIts = atoi(argv[3]);
    printf("max its %d\n",maxIts);
  }
  // For now scaling off
  model.scaling(0);
  if (maxIts) {
    // Do partial dantzig
    ClpPrimalColumnSteepest dantzig(5);
    model.setPrimalColumnPivotAlgorithm(dantzig);
    //model.messageHandler()->setLogLevel(63);
    model.setFactorizationFrequency(maxFactor);
    model.setMaximumIterations(maxIts);
    model.primal();
    if (!model.status())
      exit(1);
  }
  // find gub
  int numberRows = model.numberRows();
  int * gubStart = new int[numberRows];
  int * gubEnd = new int[numberRows];
  int * which = new int[numberRows];
  int * whichGub = new int[numberRows];
  double * lower = new double[numberRows];
  double * upper = new double[numberRows];
  const double * rowLower = model.rowLower();
  const double * rowUpper = model.rowUpper();
  int numberColumns = model.numberColumns();
  int * mark = new int[numberColumns];
  int iRow,iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++)
    mark[iColumn]=-1;
  CoinPackedMatrix * matrix = model.matrix();
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int * column = rowCopy.getIndices();
  const int * rowLength = rowCopy.getVectorLengths();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const double * element = rowCopy.getElements();
  int putGub=numberRows;
  int putNonGub=numberRows;
  for (iRow=numberRows-1;iRow>=0;iRow--) {
    bool gubRow=true;
    int first=numberColumns+1;
    int last=-1;
    for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      if (element[j]!=1.0) {
	gubRow=false;
	break;
      } else {
	int iColumn = column[j];
	if (mark[iColumn]>=0) {
	  gubRow=false;
	  break;
	} else {
	  last=max(last,iColumn);
	  first = min(first,iColumn);
	}
      }
    }
    if (last-first+1!=rowLength[iRow]||!gubRow) {
      which[--putNonGub]=iRow;
    } else {
      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	mark[iColumn]=iRow;
      }
      putGub--;
      gubStart[putGub]=first;
      gubEnd[putGub]=last+1;
      lower[putGub]=rowLower[iRow];
      upper[putGub]=rowUpper[iRow];
      whichGub[putGub]=iRow;
    }
  }
  int numberNonGub=numberRows-putNonGub;
  int numberGub = numberRows-putGub;
  if (numberGub>0) {
    printf("** %d gub rows\n",numberGub);
    for (iColumn=0;iColumn<numberColumns;iColumn++)
      mark[iColumn]=iColumn;
    ClpSimplex model2(&model,numberNonGub,which+putNonGub,numberColumns,mark);
    ClpMatrixBase * saveMatrix = model2.clpMatrix();
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
    // sort gubs so monotonic
    int * which = new int[numberGub];
    int i;
    for (i=0;i<numberGub;i++)
      which[i]=i;
    CoinSort_2(gubStart+putGub,gubStart+putGub+numberGub,which);
    int * temp1 = new int [numberGub];
    for (i=0;i<numberGub;i++) {
      int k=which[i];
      temp1[i]=gubEnd[putGub+k];
    }
    memcpy(gubEnd+putGub,temp1,numberGub*sizeof(int));
    delete [] temp1;
    double * temp2 = new double [numberGub];
    for (i=0;i<numberGub;i++) {
      int k=which[i];
      temp2[i]=lower[putGub+k];
    }
    memcpy(lower+putGub,temp2,numberGub*sizeof(double));
    for (i=0;i<numberGub;i++) {
      int k=which[i];
      temp2[i]=upper[putGub+k];
    }
    memcpy(upper+putGub,temp2,numberGub*sizeof(double));
    delete [] temp2;
    delete [] which;
    model2.replaceMatrix(new ClpGubMatrix(clpMatrix,numberGub,
					 gubStart+putGub,gubEnd+putGub,
					 lower+putGub,upper+putGub));
    clpMatrix->setMatrixNull();
    delete clpMatrix;
#if 1
    saveMatrix = model2.clpMatrix();
    ClpGubMatrix* gubMatrix =
      dynamic_cast< ClpGubMatrix*>(saveMatrix);
    assert (gubMatrix);
    // deal with basis for gub
    for (int iSet=0;iSet<numberGub;iSet++) {
      int iRow = whichGub[putGub+iSet];
      gubMatrix->setStatus(iSet,model.getRowStatus(iRow));
    }
#endif
    // For now scaling off
    model2.scaling(0);
    // Do partial dantzig
    ClpPrimalColumnSteepest dantzig(5);
    model2.setPrimalColumnPivotAlgorithm(dantzig);
    //model2.messageHandler()->setLogLevel(63);
    model2.setFactorizationFrequency(maxFactor);
    model2.setMaximumIterations(4000000);
    double time1 = CoinCpuTime();
    FILE * fp;
    if ((fp=fopen("in.sol","r"))) {
      
      ClpGubMatrix * gubMatrix =
	dynamic_cast< ClpGubMatrix*>(model2.clpMatrix());
      assert (gubMatrix);
      double * solution = model2.primalColumnSolution();
      int numberColumns = model2.numberColumns();
      int numberRows = model2.numberRows();
      char * status = new char [numberColumns];
      int numberSets = gubMatrix->numberSets();
      char * setStatus = new char[numberSets];
      char * rowStatus = new char[numberRows];
      int i;
      int n;
      n = fread(solution,sizeof(double),numberColumns,fp);
      assert (n==numberColumns);
      n = fread(status,sizeof(char),numberColumns,fp);
      assert (n==numberColumns);
      for (i=0;i<numberColumns;i++) {
	if (status[i]==0)
	  model2.setStatus(i,ClpSimplex::basic);
	else if (status[i]==1)
	  model2.setStatus(i,ClpSimplex::atLowerBound);
	else if (status[i]==2)
	  model2.setStatus(i,ClpSimplex::atUpperBound);
	else if (status[i]==3)
	  model2.setStatus(i,ClpSimplex::isFixed);
      }
      n = fread(model2.primalRowSolution(),sizeof(double),numberRows,fp);
      assert (n==numberRows);
      n = fread(rowStatus,sizeof(char),numberRows,fp);
      assert (n==numberRows);
      for (i=0;i<numberRows;i++) {
	if (rowStatus[i]==0)
	  model2.setRowStatus(i,ClpSimplex::basic);
	else if (rowStatus[i]==1)
	  model2.setRowStatus(i,ClpSimplex::atLowerBound);
	else if (rowStatus[i]==2)
	  model2.setRowStatus(i,ClpSimplex::atUpperBound);
	else if (rowStatus[i]==3)
	  model2.setRowStatus(i,ClpSimplex::isFixed);
      }
      n = fread(setStatus,sizeof(char),numberSets,fp);
      assert (n==numberSets);
      n = fread(gubMatrix->keyVariable(),sizeof(int),numberSets,fp);
      assert (n==numberSets);
      const int * keyVariable = gubMatrix->keyVariable();
      for (i=0;i<numberSets;i++) {
	if (setStatus[i]==0)
	  gubMatrix->setStatus(i,ClpSimplex::basic);
	else if (setStatus[i]==1)
	  gubMatrix->setStatus(i,ClpSimplex::atLowerBound);
	else if (setStatus[i]==2)
	  gubMatrix->setStatus(i,ClpSimplex::atUpperBound);
	else if (setStatus[i]==3)
	  gubMatrix->setStatus(i,ClpSimplex::isFixed);
	int iKey = keyVariable[i];
	if (iKey<numberColumns)
	  model2.setStatus(iKey,ClpSimplex::superBasic); // don't want in basis as make get thrown out
      }
      fclose(fp);
      delete [] status;
      delete [] setStatus;
      delete [] rowStatus;
    }
    model2.primal(1);
    printf("Primal took %g seconds\n",CoinCpuTime()-time1);
    if (1) {
      ClpGubMatrix * gubMatrix =
	dynamic_cast< ClpGubMatrix*>(model2.clpMatrix());
      assert (gubMatrix);
      const double * solution = model2.primalColumnSolution();
      int numberColumns = model2.numberColumns();
      int numberRows = model2.numberRows();
      char * status = new char [numberColumns];
      int numberSets = gubMatrix->numberSets();
      char * setStatus = new char[numberSets];
      char * rowStatus = new char[numberRows];
      const double * lowerColumn = model2.columnLower();
      const double * upperColumn = model2.columnUpper();
      const int * backward = gubMatrix->backward();
      const int * keyVariable = gubMatrix->keyVariable();
      int i;
      FILE * fp=fopen ("xx.sol","w");
      fwrite(solution,sizeof(double),numberColumns,fp);
      for (i=0;i<numberColumns;i++) {
	ClpSimplex::Status thisStatus = model2.getStatus(i);
	if (thisStatus==ClpSimplex::basic) {
	  status[i]=0;
	} else if (thisStatus==ClpSimplex::atLowerBound) {
	  status[i]=1;
	} else if (thisStatus==ClpSimplex::atUpperBound) {
	  status[i]=2;
	} else if (thisStatus==ClpSimplex::isFixed) {
	  status[i]=3;
	} else {
	  // see if key
	  int iBack = backward[i];
	  if (iBack>=0&&i==keyVariable[iBack]) {
	    status[i]=0; // key
	    //model2.setStatus(i,ClpSimplex::basic);
	  } else {
	    if (fabs(solution[i]-lowerColumn[i])<1.0e-6)
	      status[i]=1;
	    else if (fabs(solution[i]-upperColumn[i])<1.0e-6)
	      status[i]=2;
	    else
	      abort();
	  }
	}
      }
      fwrite(status,sizeof(char),numberColumns,fp);
      fwrite(model2.primalRowSolution(),sizeof(double),numberRows,fp);
      for (i=0;i<numberRows;i++) {
	ClpSimplex::Status thisStatus = model2.getRowStatus(i);
	if (thisStatus==ClpSimplex::basic)
	  rowStatus[i]=0;
	else if (thisStatus==ClpSimplex::atLowerBound)
	  rowStatus[i]=1;
	else if (thisStatus==ClpSimplex::atUpperBound)
	  rowStatus[i]=2;
	else if (thisStatus==ClpSimplex::isFixed)
	  rowStatus[i]=3;
	else
	  abort();
      }
      fwrite(rowStatus,sizeof(char),numberRows,fp);
      for (i=0;i<numberSets;i++) {
	// set key as non-basic
	if (keyVariable[i]<numberColumns)
	  model2.setStatus(keyVariable[i],ClpSimplex::superBasic);
	ClpSimplex::Status thisStatus = gubMatrix->getStatus(i);
	if (thisStatus==ClpSimplex::basic)
	  setStatus[i]=0;
	else if (thisStatus==ClpSimplex::atLowerBound)
	  setStatus[i]=1;
	else if (thisStatus==ClpSimplex::atUpperBound)
	  setStatus[i]=2;
	else if (thisStatus==ClpSimplex::isFixed)
	  setStatus[i]=3;
	else
	  abort();
      }
      fwrite(setStatus,sizeof(char),numberSets,fp);
      fwrite(gubMatrix->keyVariable(),sizeof(int),numberSets,fp);
      fclose(fp);
      delete [] status;
      delete [] setStatus;
      delete [] rowStatus;
    }
    model2.primal();
  } else {
    // dummy gub
    ClpMatrixBase * saveMatrix = model.clpMatrix();
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
    model.replaceMatrix(new ClpGubMatrix(clpMatrix,0,NULL,NULL,NULL,NULL)); 
    // For now scaling off
    model.scaling(0);
    ClpPrimalColumnSteepest dantzig(5);
    model.setPrimalColumnPivotAlgorithm(dantzig);
    double time1 = CoinCpuTime();
    model.primal();
    printf("Primal took %g seconds\n",CoinCpuTime()-time1);
  }
  delete [] mark;
  delete [] gubStart;
  delete [] gubEnd;
  delete [] which;
  delete [] whichGub;
  delete [] lower;
  delete [] upper;
  return 0;
}    
