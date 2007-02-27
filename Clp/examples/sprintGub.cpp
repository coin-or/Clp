// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"
#include <iomanip>
static void sprint(ClpSimplex & model,
		   const int * setToRow,  const int * whichGub,
		   const int * whichColumn,  const int * startSet,
		   const char * rowType,
		   int numberGub,double smallFactor)
{
  int numberColumns = model.numberColumns();
  int numberRows = model.numberRows();

  // We will need an array to choose variables.
  int * sort = new int [numberColumns];
  float * weight = new float [numberColumns];
  const int * ordRow = setToRow+numberGub;
  int numberSort=0;
  // Set up initial list
  numberSort=numberRows;
  int i;
  for (i=0;i<numberSort;i++)
    sort[i] = i;
  // and basic
  int nOrd=numberColumns - startSet[numberGub];
  int iRow,iColumn;
  const double * columnLower = model.columnLower();
  const double * columnUpper = model.columnUpper();
  double * fullSolution = model.primalColumnSolution();
  CoinPackedMatrix * matrix = model.matrix();
  const int * columnStart = matrix->getVectorStarts();
  for (iColumn=numberRows;iColumn<numberColumns;iColumn++) {
    if(model.getColumnStatus(iColumn)==ClpSimplex::basic||
       columnStart[iColumn+1]==columnStart[iColumn]+1) {
      sort[numberSort++]=iColumn;
    } else if (fullSolution[iColumn]>columnLower[iColumn]+1.0e-7&&
	       fullSolution[iColumn]<columnUpper[iColumn]-1.0e-7) {
      sort[numberSort++]=iColumn;
    }
    if (fullSolution[iColumn]>columnLower[iColumn]-1.0e-7&&
	fullSolution[iColumn]<columnUpper[iColumn]+1.0e-7) {
    } else {
      printf("bad %d %g %g %g\n",iColumn,columnLower[iColumn],
	     fullSolution[iColumn],columnUpper[iColumn]);
    }
    if (fullSolution[iColumn]<columnLower[iColumn])
      fullSolution[iColumn]=columnLower[iColumn];
    if (fullSolution[iColumn]>columnUpper[iColumn])
      fullSolution[iColumn]=columnUpper[iColumn];
  }
  double time1 = CoinCpuTime();

  // Just do this number of passes
  int maxPass=1000;
  int iPass;
  double lastObjective=1.0e31;

  // Just take this number of columns in small problem
  int smallNumberColumns = CoinMin(((int) (smallFactor*numberRows)),numberColumns);
  smallNumberColumns = CoinMax(smallNumberColumns,3000);
  // We will not be using all rows
  int * whichRows = new int [numberRows];
  // except first time
  int nOrdRow=numberRows-numberGub;
  for (iRow=0;iRow<nOrdRow;iRow++)
    whichRows[iRow]=ordRow[iRow];
  for (iRow=0;iRow<numberGub;iRow++)
    whichRows[iRow+nOrdRow]=setToRow[iRow];
  int smallNumberRows=numberRows;
  // Basic counts for each set
  int * counts = new int[numberGub];
  memset(counts,0,numberGub*sizeof(int));
  // Negative dj counts
  int * negCounts = new int[numberGub];
  // Negative dj sum
  double * negSum = new double[numberGub];
  double originalOffset;
  model.getDblParam(ClpObjOffset,originalOffset);

  for (iPass=0;iPass<maxPass;iPass++) {
    printf("Start of pass %d\n",iPass);
    // Create small problem
    ClpSimplex small(&model,smallNumberRows,whichRows,numberSort,sort);
    small.setPerturbation(100);
    small.setInfeasibilityCost(1.0e13);
    small.defaultFactorizationFrequency();
    small.scaling(0);
    // now see what variables left out do to row solution
    double * rowSolution = model.primalRowSolution();
    memset (rowSolution,0,numberRows*sizeof(double));
    // zero out ones in small problem
    for (iColumn=0;iColumn<numberSort;iColumn++) {
      int kColumn = sort[iColumn];
      fullSolution[kColumn]=0.0;
    }
    // Get objective offset
    double offset=0.0;
    const double * objective = model.objective();
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      offset += fullSolution[iColumn]*objective[iColumn];
    small.setDblParam(ClpObjOffset,originalOffset-offset);
    model.times(1.0,fullSolution,rowSolution);

    double * lower = small.rowLower();
    double * upper = small.rowUpper();
    for (iRow=0;iRow<smallNumberRows;iRow++) {
      int jRow = whichRows[iRow];
      if (lower[iRow]>-1.0e50) 
	lower[iRow] -= rowSolution[jRow];
      if (upper[iRow]<1.0e50)
	upper[iRow] -= rowSolution[jRow];
    }
    if (iPass) {
      memset (rowSolution,0,smallNumberRows*sizeof(double));
      small.times(1.0,small.primalColumnSolution(),rowSolution);
      for (iRow=0;iRow<smallNumberRows;iRow++) {
	if (rowSolution[iRow]<lower[iRow]-1.0e-5||
	    rowSolution[iRow]>upper[iRow]+1.0e-5)
	  printf("bad %d %g %g %g\n",iRow,lower[iRow],rowSolution[iRow],
		 upper[iRow]);
      }
    }
    // Solve 
    small.setMaximumIterations(50000);
    small.primal(1);
    if (small.problemStatus()!=0) {
      //exit(88);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (fabs(fullSolution[iColumn])>1.0e-8)
	  printf("full %d %g\n",iColumn,fullSolution[iColumn]);
      }
      small.setLogLevel(63);
      small.writeMps("infeas.mps");
      small.primal();
      for (iRow=0;iRow<smallNumberRows;iRow++) {
	double value = small.primalRowSolution()[iRow];
	if (value<lower[iRow]-1.0e-7||value>upper[iRow]+1.0e-7)
	  printf("%d %g %g %g\n",iRow,lower[iRow],value,upper[iRow]);
      }
      small.dual();
      small.setLogLevel(1);
      if (small.problemStatus()!=0) {
	exit(77);
      }
    }
    // move solution back
    const double * solution = small.primalColumnSolution();
    for (iColumn=0;iColumn<numberSort;iColumn++) {
      int kColumn = sort[iColumn];
      model.setColumnStatus(kColumn,small.getColumnStatus(iColumn));
      fullSolution[kColumn]=solution[iColumn];
      if (fullSolution[iColumn]>columnLower[iColumn]-1.0e-7&&
	  fullSolution[iColumn]<columnUpper[iColumn]+1.0e-7) {
      } else {
	printf("bad2 %d\n",iColumn);
      }
    }
    model.setObjectiveValue(small.objectiveValue());
    double * rowSol1 = model.primalRowSolution();
    double * rowSol2 = small.primalRowSolution();
    double * rowPi1 = model.dualRowSolution();
    double * rowPi2 = small.dualRowSolution();
    for (iRow=0;iRow<smallNumberRows;iRow++) {
      int jRow = whichRows[iRow];
      model.setRowStatus(jRow,small.getRowStatus(iRow));
      rowSol1[jRow] = rowSol2[iRow];
      rowPi1[jRow] = rowPi2[iRow];
    }
    if ((small.objectiveValue()>lastObjective-1.0e-7&&iPass>5)||
	(!small.numberIterations()&&iPass)||
	iPass==maxPass-1) {

      break; // finished
    } else {
      lastObjective = small.objectiveValue();
      // get reduced cost for large problem
      // this assumes minimization
      // First get duals for key variables
      CoinPackedMatrix * matrix = model.matrix();
      const int * row = matrix->getIndices();
      //const int * columnLength = matrix->getVectorLengths();
      const CoinBigIndex * columnStart = matrix->getVectorStarts();
      const double * elementByColumn = matrix->getElements();
      const double * obj = model.objective();
      double * rowPi1 = model.dualRowSolution();
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if(model.getColumnStatus(iColumn)==ClpSimplex::basic) {
	  int iSet = whichGub[iColumn];
	  if (iSet>=0&&!counts[iSet]) {
	    double dj = obj[iColumn];
	    int kstart = columnStart[iColumn];
	    int kend = columnStart[iColumn+1];
	    int j;
	    int gubRow = setToRow[iSet];
	    //double x=rowPi1[gubRow];
	    rowPi1[gubRow]=0.0;
	    //bool found=false;
	    for (j=kstart;j<kend;j++) {
	      int iRow= row[j];
	      dj -= rowPi1[iRow]*elementByColumn[j];
	      //if (iRow==gubRow) {
	      //found=true;
	      //assert (elementByColumn[j]==1.0);
	      //}
	    }
	    //assert (found);
	    //assert (fabs(dj-x)<1.0e-2);
	    rowPi1[gubRow]=dj;
	  }
	}
      }
      memset(counts,0,numberGub*sizeof(int));
      memset(negCounts,0,numberGub*sizeof(int));
      memset(negSum,0,numberGub*sizeof(double));
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double dj = obj[iColumn];
	int kstart = columnStart[iColumn];
	int kend = columnStart[iColumn+1];
	int j;
	for (j=kstart;j<kend;j++) {
	  int iRow= row[j];
	  dj -= rowPi1[iRow]*elementByColumn[j];
	}
	double value = fullSolution[iColumn];
	int iSet = whichGub[iColumn];
	if (iSet<0) {
	  if (model.getColumnStatus(iColumn)==ClpSimplex::basic)
	    dj = -1.0e50;
	  else if (dj<-1.0e-7&&value<columnUpper[iColumn])
	    dj = dj;
	  else if (dj>1.0e-7&&value>columnLower[iColumn])
	    dj = -dj;
	  else if (columnUpper[iColumn]>columnLower[iColumn])
	    dj = fabs(dj);
	  else
	    dj = 1.0e50;
	  // always take
          dj=-1.0e30;
	} else {
	  if (model.getColumnStatus(iColumn)==ClpSimplex::basic) {
	    dj = -1.0e50;
	    counts[iSet]++;
	  } else if (dj<-1.0e-7&&value<columnUpper[iColumn]) {
	    dj = dj;
	    negCounts[iSet]++;
	    negSum[iSet] -= dj;
	  } else if (dj>1.0e-7&&value>columnLower[iColumn]) {
	    dj = -dj;
	    negCounts[iSet]++;
	    negSum[iSet] -= dj;
	  } else if (columnUpper[iColumn]>columnLower[iColumn]) {
	    dj = fabs(dj);
	  } else {
	    dj = 1.0e50;
	  }
	}
	weight[iColumn]=dj;
      }
      int stats[100];
      memset(stats,0,sizeof(stats));
      int iSet;
      double sumDj=0.0;
      int numberDj=0;
      for (iSet=0;iSet<numberGub;iSet++) {
	int k = counts[iSet];
	if (k<100)
	  stats[k]++;
	sumDj += negSum[iSet];
	numberDj += negCounts[iSet];
	negCounts[iSet]=iSet;
	if (k>1)
	  negSum[iSet]=-1.0e50;
	else
	  negSum[iSet]=-negSum[iSet];
      }
      double averageDj = sumDj/((double) (numberDj+1));
      double rejectDj = CoinMax(10.0*averageDj,1.0);
      int algorithm = 1;
      smallNumberRows=0;
      //for (iSet=0;iSet<100;iSet++)
      //if (stats[iSet])
      //  printf("%d sets have %d basic\n",stats[iSet],iSet);
      printf("%d negative djs summing to %g\n",numberDj,sumDj);
      for (iRow=0;iRow<nOrdRow;iRow++) {
	int kRow = ordRow[iRow];
	if (rowType[kRow]>=0)
	  whichRows[smallNumberRows++]=kRow;
      }
      memset(counts,0,numberGub*sizeof(int));
      if (!algorithm) {
	// sort
	CoinSort_2(negSum,negSum+numberGub,negCounts);
	int n=nOrd;
	int nReject=0;
	for (iSet=0;iSet<numberGub;iSet++) {
	  int jSet = negCounts[iSet];
	  int type = counts[jSet];
	  double reject = 1.0e-1;
	  if (type==2)
	    reject = rejectDj;
	  else if (type>2)
	    reject = 10.0*rejectDj;
	  counts[jSet]=1;
	  for (int j=startSet[jSet];j<startSet[jSet+1];j++) {
	    iColumn = whichColumn[j];
	    float dj = weight[iColumn];
	    if (dj<reject) {
	      weight[iColumn]=-1.0e20;
	      n++;
	    } else {
	      weight[iColumn]=0.0;
	      nReject++;
	    }
	  }
	  whichRows[smallNumberRows++]=setToRow[jSet];
	  if (n>smallNumberColumns&&negSum[iSet]>-1.0e30) {
	    break;
	  }
	}
	printf("%d variables (%d not used) in %d sets - rows %d\n",n,
	       nReject,iSet+1,
	       smallNumberRows);
	numberSort=0;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  iSet = whichGub[iColumn];
	  if (iSet<0||counts[iSet]) {
	    if (weight[iColumn]<-1.0e10)
	      sort[numberSort++] = iColumn;
	  }
	}
      } else {
	// normal way
	numberSort = 0;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  int iSet = whichGub[iColumn];
	  if (iSet<0||model.getColumnStatus(iColumn)==ClpSimplex::basic) {
	    sort[numberSort++] = iColumn;
	    if (iSet>=0)
	      counts[iSet]++;
	  }
	}
	double tolerance;
	if (numberDj+numberSort>smallNumberColumns)
	  tolerance = -model.dualTolerance();
	else 
	  tolerance = 10.0*averageDj;
	int n2=0;
	int saveN = numberSort;
	int * sort2 = sort+saveN;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  int iSet = whichGub[iColumn];
	  if (iSet>=0&&model.getColumnStatus(iColumn)!=ClpSimplex::basic) {
	    double dj = weight[iColumn];
	    if (dj<tolerance) {
	      weight[n2] = dj;
	      sort2[n2++] = iColumn;
	    }
	  }
	}
	// sort
	CoinSort_2(weight,weight+n2,sort+saveN);
	numberSort += n2;
	numberSort = CoinMin(smallNumberColumns,numberSort);
#if 1
	for (i=saveN;i<numberSort;i++) {
	  iColumn = sort[i];
	  int iSet = whichGub[iColumn];
	  assert (iSet>=0);
	  counts[iSet]++;
	}
	// take out ones with just key
	int n = numberSort;
	numberSort = 0;
	for (i=0;i<n;i++) {
	  iColumn = sort[i];
	  int iSet = whichGub[iColumn];
	  if (iSet<0||counts[iSet]>1) {
	    sort[numberSort++]=iColumn;
	  }
	}
	n=smallNumberRows;
	for (iSet=0;iSet<numberGub;iSet++) {
	  if (counts[iSet]>1) {
	    whichRows[smallNumberRows++]=setToRow[iSet];
	  } else {
	    //whichRows[smallNumberRows++]=setToRow[iSet];
	    counts[iSet]=0; // so pi will be computed
	  }
	}
	printf("%d variables in %d sets - rows %d\n",numberSort,
	       smallNumberRows-n,
	       smallNumberRows);
#else
	for (iSet=0;iSet<numberGub;iSet++) {
	  whichRows[smallNumberRows++]=setToRow[iSet];
	}
#endif
      }
    }
  }
  delete [] counts;
  delete [] negCounts;
  delete [] negSum;
  delete [] weight;
  delete [] sort;
  delete [] whichRows;
  printf("Sprint took %g seconds\n",CoinCpuTime()-time1);
}

int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  // Keep names
  if (argc<2) {
    status=model.readMps("small.mps",true);
  } else {
    status=model.readMps(argv[1],true);
  }
  if (status)
    exit(10);
  /*
    This driver implements what I called Sprint.  Cplex calls it 
    "sifting" which is just as silly.  When I thought of this trivial idea
    it reminded me of an LP code of the 60's called sprint which after
    every factorization took a subset of the matrix into memory (all
    64K words!) and then iterated very fast on that subset.  On the
    problems of those days it did not work very well, but it worked very
    well on aircrew scheduling problems where there were very large numbers
    of columns all with the same flavor.
  */

  /* The idea works best if you can get feasible easily.  To make it
     more general we can add in costed slacks */

  int originalNumberColumns = model.numberColumns();
  int numberRows = model.numberRows();

  int * sort = new int [numberRows+originalNumberColumns];
  int numberSort=0;
  // Say we are going to add slacks - if you can get a feasible
  // solution then do that at the comment - Add in your own coding here
  bool addSlacks = true;

  if (addSlacks) {
    // initial list will just be artificials
    // first we will set all variables as close to zero as possible
    int iColumn;
    const double * columnLower = model.columnLower();
    const double * columnUpper = model.columnUpper();
    double * columnSolution = model.primalColumnSolution();

    for (iColumn=0;iColumn<originalNumberColumns;iColumn++) {
      double value =0.0;
      if (columnLower[iColumn]>0.0)
	value = columnLower[iColumn];
      else if (columnUpper[iColumn]<0.0)
	value = columnUpper[iColumn];
      columnSolution[iColumn]=value;
    }
    // now see what that does to row solution
    double * rowSolution = model.primalRowSolution();
    memset (rowSolution,0,numberRows*sizeof(double));
    model.times(1.0,columnSolution,rowSolution);

    int * addStarts = new int [numberRows+1];
    int * addRow = new int[numberRows];
    double * addElement = new double[numberRows];
    const double * lower = model.rowLower();
    const double * upper = model.rowUpper();
    addStarts[0]=0;
    int numberArtificials=0;
    double * addCost = new double [numberRows];
    const double penalty=1.0e8;
    int iRow;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (lower[iRow]>rowSolution[iRow]) {
	addRow[numberArtificials]=iRow;
	addElement[numberArtificials]=1.0;
	addCost[numberArtificials]=penalty;
	numberArtificials++;
	addStarts[numberArtificials]=numberArtificials;
      } else if (upper[iRow]<rowSolution[iRow]) {
	addRow[numberArtificials]=iRow;
	addElement[numberArtificials]=-1.0;
	addCost[numberArtificials]=penalty;
	numberArtificials++;
	addStarts[numberArtificials]=numberArtificials;
      }
    }
    model.addColumns(numberArtificials,NULL,NULL,addCost,
		      addStarts,addRow,addElement);
    delete [] addStarts;
    delete [] addRow;
    delete [] addElement;
    delete [] addCost;
    // Set up initial list
    numberSort=numberArtificials;
    int i;
    for (i=0;i<numberSort;i++)
      sort[i] = i+originalNumberColumns;
  } else {
    // Get initial list in some magical way
    // Add in your own coding here
    abort();
  }
  // gub stuff
  {
    int numberColumns = model.numberColumns();
    int numberRows = model.numberRows();
    
    int * whichGub = new int[numberColumns];
    int * setToRow = new int[numberRows];
    int * ordRow = setToRow+numberRows;
    int * startSet = new int [numberRows];
    int * whichColumn = new int [numberColumns];
    int numberGub=0;
    int nOrd=0;
    int iRow,iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      whichGub[iColumn]=-1;
    }
    CoinPackedMatrix * matrix = model.matrix();
    // get row copy
    CoinPackedMatrix rowCopy = *matrix;
    rowCopy.reverseOrdering();
    const int * column = rowCopy.getIndices();
#ifndef NDEBUG
    const int * rowLength = rowCopy.getVectorLengths();
#endif
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    for (iRow=0;iRow<numberRows;iRow++) {
      // check contiguous
      assert (rowStart[iRow]+rowLength[iRow]==rowStart[iRow+1]);
    }
    const double * element = rowCopy.getElements();
    numberGub=0;
    startSet[0]=0;
    int nel=0;
    for (iRow=numberRows-1;iRow>=0;iRow--) {
      bool gubRow=true;
      for (int j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	if (element[j]!=1.0) {
	  gubRow=false;
	  break;
	} else {
	  int iColumn = column[j];
	  if (whichGub[iColumn]>=0) {
	    gubRow=false;
	    break;
	  }
	}
      }
      if (!gubRow) {
	ordRow--;
	*ordRow=iRow;
      } else {
	for (int j=rowStart[iRow];j<rowStart[iRow+1];j++) {
	  int iColumn = column[j];
	  whichGub[iColumn]=numberGub;
	  whichColumn[nel++]=iColumn;
	}
	setToRow[numberGub]=iRow;
	startSet[numberGub+1]=nel;
	numberGub++;
      }
    }
    nOrd=numberColumns-startSet[numberGub];
    printf("%d ordinary columns, %d gub sets\n",nOrd,numberGub);
    assert (numberGub<numberRows); // otherwise array bad
    char * rowType = new char[numberRows];
    memset(rowType,0,numberRows);
    sprint( model, setToRow, whichGub,
	    whichColumn,  startSet, rowType, numberGub,3);
    delete [] setToRow;
    delete [] whichGub;
    delete [] whichColumn;
    delete [] startSet;
    delete [] rowType;
    model.primal(1);
  }
  if (addSlacks) {
    int i;
    int numberColumns = model.numberColumns();
    int numberArtificials = numberColumns-originalNumberColumns;
    for (i=0;i<numberArtificials;i++)
      sort[i] = i + originalNumberColumns;
    model.deleteColumns(numberArtificials,sort);
  }
  delete [] sort;
  model.primal(1);
  return 0;
}    
