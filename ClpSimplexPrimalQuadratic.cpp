// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexPrimalQuadratic.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpPrimalColumnPivot.hpp"
#include "ClpMessage.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream>

// A sequential LP method
int 
ClpSimplexPrimalQuadratic::primalSLP(int numberPasses, double deltaTolerance)
{
  // Are we minimizing or maximizing
  double whichWay=optimizationDirection();
  // This is as a user would see

  int numberColumns = this->numberColumns();
  double * columnLower = this->columnLower();
  double * columnUpper = this->columnUpper();
  double * objective = this->objective();
  double * solution = this->primalColumnSolution();
  
  // Save objective
  
  double * saveObjective = new double [numberColumns];
  memcpy(saveObjective,objective,numberColumns*sizeof(double));

  // Get list of non linear columns
  CoinPackedMatrix * quadratic = quadraticObjective();
  if (!quadratic) {
    // no quadratic part
    return primal(0);
  }
  int numberNonLinearColumns = 0;
  int iColumn;
  int * listNonLinearColumn = new int[numberColumns];
  memset(listNonLinearColumn,0,numberColumns*sizeof(int));
  const int * columnQuadratic = quadratic->getIndices();
  const int * columnQuadraticStart = quadratic->getVectorStarts();
  const int * columnQuadraticLength = quadratic->getVectorLengths();
  const double * quadraticElement = quadratic->getElements();
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
      int jColumn = columnQuadratic[j];
      listNonLinearColumn[jColumn]=1;
      listNonLinearColumn[iColumn]=1;
    }
  }
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if(listNonLinearColumn[iColumn])
      listNonLinearColumn[numberNonLinearColumns++]=iColumn;
  }
  
  if (!numberNonLinearColumns) {
    delete [] listNonLinearColumn;
    // no quadratic part
    return primal(0);
  }

  // get feasible
  if (numberPrimalInfeasibilities())
    primal(1);
  // still infeasible
  if (numberPrimalInfeasibilities())
    return 0;

  int jNon;
  int * last[3];
  
  double * trust = new double[numberNonLinearColumns];
  double * trueLower = new double[numberNonLinearColumns];
  double * trueUpper = new double[numberNonLinearColumns];
  for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
    iColumn=listNonLinearColumn[jNon];
    trust[jNon]=0.5;
    trueLower[jNon]=columnLower[iColumn];
    trueUpper[jNon]=columnUpper[iColumn];
  }
  int iPass;
  double lastObjective=1.0e31;
  double * saveSolution = new double [numberColumns];
  double targetDrop=1.0e31;
  double objectiveOffset;
  getDblParam(ClpObjOffset,objectiveOffset);
  // 1 bound up, 2 up, -1 bound down, -2 down, 0 no change
  for (iPass=0;iPass<3;iPass++) {
    last[iPass]=new int[numberNonLinearColumns];
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) 
      last[iPass][jNon]=0;
  }
  double goodMove=false;
  for (iPass=0;iPass<numberPasses;iPass++) {
    // redo objective
    double offset=0.0;
    double objValue=-objectiveOffset;
    memcpy(objective,saveObjective,numberColumns*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      objValue += objective[iColumn]*solution[iColumn];
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      double valueI = solution[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution[jColumn];
	double elementValue = quadraticElement[j];
	objValue += 0.5*valueI*valueJ*elementValue;
	offset += 0.5*valueI*valueJ*elementValue;
	double gradientI = valueJ*elementValue;
	double gradientJ = valueI*elementValue;
	offset -= gradientI*valueI;
	objective[iColumn] += gradientI;
	offset -= gradientJ*valueJ;
	objective[jColumn] += gradientJ;
      }
    }
    printf("objective %g, objective offset %g\n",objValue,offset);
    setDblParam(ClpObjOffset,objectiveOffset-offset);
    objValue *= whichWay;
    int * temp=last[2];
    last[2]=last[1];
    last[1]=last[0];
    last[0]=temp;
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      double change = solution[iColumn]-saveSolution[iColumn];
      if (change<-1.0e-5) {
	if (fabs(change+trust[jNon])<1.0e-5) 
	  temp[jNon]=-1;
	else
	  temp[jNon]=-2;
      } else if(change>1.0e-5) {
	if (fabs(change-trust[jNon])<1.0e-5) 
	  temp[jNon]=1;
	else
	  temp[jNon]=2;
      } else {
	temp[jNon]=0;
      }
    } 
    if (objValue<=lastObjective) 
      goodMove=true;
    else
      goodMove=false;
    double maxDelta=0.0;
    double maxGap=0.0;
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      maxDelta = max(maxDelta,
		     fabs(solution[iColumn]-saveSolution[iColumn]));
      if (goodMove) {
	if (last[0][jNon]*last[1][jNon]<0) {
	  // halve
	  trust[jNon] *= 0.5;
	} else {
	  if (last[0][jNon]==last[1][jNon]&&
	      last[0][jNon]==last[2][jNon])
	    trust[jNon] *= 1.5; 
	}
      } else if (trust[jNon]>10.0*deltaTolerance) {
	trust[jNon] *= 0.5;
      }
      maxGap = max(maxGap,trust[jNon]);
    }
    std::cout<<"largest gap is "<<maxGap<<std::endl;
    if (goodMove) {
      double drop = lastObjective-objValue;
      std::cout<<"True drop was "<<drop<<std::endl;
      std::cout<<"largest delta is "<<maxDelta<<std::endl;
      if (maxDelta<deltaTolerance&&drop<1.0e-4) {
	std::cout<<"Exiting"<<std::endl;
	break;
      }
    } else {
      lastObjective += 1.0e-3; // to stop exit
    }
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      columnLower[iColumn]=max(solution[iColumn]
			       -trust[jNon],
			       trueLower[jNon]);
      columnUpper[iColumn]=min(solution[iColumn]
			       +trust[jNon],
			       trueUpper[jNon]);
    }
    if (goodMove) {
      memcpy(saveSolution,solution,numberColumns*sizeof(double));
      
      targetDrop=0.0;
      if (iPass) {
	// get reduced costs
	this->matrix()->transposeTimes(this->dualRowSolution(),
				       this->dualColumnSolution());
	double * r = this->dualColumnSolution();
	for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
	  iColumn=listNonLinearColumn[jNon];
	  double dj = objective[iColumn]-r[iColumn];
	  if (dj<0.0) 
	    targetDrop -= dj*(columnUpper[iColumn]-solution[iColumn]);
	  else
	    targetDrop -= dj*(columnLower[iColumn]-solution[iColumn]);
	}
      }
      std::cout<<"Pass - "<<iPass
	       <<", target drop is "<<targetDrop
	       <<std::endl;
      lastObjective = objValue;
      this->primal(1);
    } else {
      // bad pass - restore solution
      memcpy(solution,saveSolution,numberColumns*sizeof(double));
      iPass--;
    }
  }
  // restore solution
  memcpy(solution,saveSolution,numberColumns*sizeof(double));
  for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
    iColumn=listNonLinearColumn[jNon];
    columnLower[iColumn]=max(solution[iColumn],
			     trueLower[jNon]);
    columnUpper[iColumn]=min(solution[iColumn],
			     trueUpper[jNon]);
  }
  // redo objective
  double offset=0.0;
  double objValue=-objectiveOffset;
  memcpy(objective,saveObjective,numberColumns*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    objValue += objective[iColumn]*solution[iColumn];
  for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
    iColumn=listNonLinearColumn[jNon];
    double valueI = solution[iColumn];
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
      int jColumn = columnQuadratic[j];
      double valueJ = solution[jColumn];
      double elementValue = quadraticElement[j];
      objValue += 0.5*valueI*valueJ*elementValue;
      offset += 0.5*valueI*valueJ*elementValue;
      double gradientI = valueJ*elementValue;
      double gradientJ = valueI*elementValue;
      offset -= gradientI*valueI;
      objective[iColumn] += gradientI;
      offset -= gradientJ*valueJ;
      objective[jColumn] += gradientJ;
    }
  }
  printf("objective %g, objective offset %g\n",objValue,offset);
  setDblParam(ClpObjOffset,objectiveOffset-offset);
  this->primal(1);
  for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
    iColumn=listNonLinearColumn[jNon];
    columnLower[iColumn]= trueLower[jNon];
    columnUpper[iColumn]= trueUpper[jNon];
  }
  delete [] saveSolution;
  for (iPass=0;iPass<3;iPass++) 
    delete [] last[iPass];
  delete [] trust;
  delete [] trueUpper;
  delete [] trueLower;
  delete [] saveObjective;
  delete [] listNonLinearColumn;
  return 0;
}
// Beale's method
int 
ClpSimplexPrimalQuadratic::primalBeale()
{
  // Wolfe's method looks easier - lets try that
  // This is as a user would see

  int numberColumns = this->numberColumns();
  double * columnLower = this->columnLower();
  double * columnUpper = this->columnUpper();
  double * objective = this->objective();
  double * solution = this->primalColumnSolution();
  double * dj = this->dualColumnSolution();
  double * pi = this->dualRowSolution();

  int numberRows = this->numberRows();
  double * rowLower = this->rowLower();
  double * rowUpper = this->rowUpper();

  // and elements
  CoinPackedMatrix * matrix = this->matrix();
  const int * row = matrix->getIndices();
  const int * columnStart = matrix->getVectorStarts();
  const double * element =  matrix->getElements();
  const int * columnLength = matrix->getVectorLengths();

  // Get list of non linear columns
  CoinPackedMatrix * quadratic = quadraticObjective();
  if (!quadratic||!quadratic->getNumElements()) {
    // no quadratic part
    return primal(1);
  }

  int iColumn;
  const int * columnQuadratic = quadratic->getIndices();
  const int * columnQuadraticStart = quadratic->getVectorStarts();
  const int * columnQuadraticLength = quadratic->getVectorLengths();
  const double * quadraticElement = quadratic->getElements();
  // Get a feasible solution 
  if (numberPrimalInfeasibilities())
    primal(1);
  // still infeasible
  if (numberPrimalInfeasibilities())
    return 0;
  
  // Create larger problem
  int newNumberRows = numberRows+numberColumns;
  // See how many artificials we will need
  // For now assume worst
  int newNumberColumns = 3*numberColumns+ numberRows;
  int numberElements = 2*matrix->getNumElements()
    +quadratic->getNumElements()
    + 2*numberColumns;
  double * elements2 = new double[numberElements];
  int * start2 = new int[newNumberColumns+1];
  int * row2 = new int[numberElements];
  double * objective2 = new double[newNumberColumns];
  double * columnLower2 = new double[newNumberColumns];
  double * columnUpper2 = new double[newNumberColumns];
  double * rowLower2 = new double[newNumberRows];
  double * rowUpper2 = new double[newNumberRows];
  memset(rowLower2,0,newNumberRows*sizeof(double));
  memcpy(rowLower2,rowLower,numberRows*sizeof(double));
  memcpy(rowLower2+numberRows,objective,numberColumns*sizeof(double));
  memset(rowUpper2,0,newNumberRows*sizeof(double));
  memcpy(rowUpper2,rowUpper,numberRows*sizeof(double));
  memcpy(rowUpper2+numberRows,objective,numberColumns*sizeof(double));
  memset(objective2,0,newNumberColumns*sizeof(double));
  // Get a row copy of quadratic objective in standard format
  CoinPackedMatrix copyQ;
  copyQ.reverseOrderedCopyOf(*quadratic);
  const int * columnQ = copyQ.getIndices();
  const CoinBigIndex * rowStartQ = copyQ.getVectorStarts();
  const int * rowLengthQ = copyQ.getVectorLengths(); 
  const double * elementByRowQ = copyQ.getElements();
  newNumberColumns=0;
  numberElements=0;
  start2[0]=0;
  // x
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // Original matrix
    columnLower2[iColumn]=columnLower[iColumn];
    columnUpper2[iColumn]=columnUpper[iColumn];
    int j;
    for (j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];
	 j++) {
      elements2[numberElements]=element[j];
      row2[numberElements++]=row[j];
    }
    // Quadratic and modify djs
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];
	 j++) {
      int jColumn = columnQuadratic[j];
      double value = quadraticElement[j];
      if (iColumn!=jColumn) 
	value *= 0.5;
      dj[iColumn] += solution[jColumn]*value;
      elements2[numberElements]=-value;
      row2[numberElements++]=jColumn+numberRows;
    }
    for (j=rowStartQ[iColumn];
	 j<rowStartQ[iColumn]+rowLengthQ[iColumn];
	 j++) {
      int jColumn = columnQ[j];
      double value = elementByRowQ[j];
      if (iColumn!=jColumn) { 
	value *= 0.5;
	dj[iColumn] += solution[jColumn]*value;
	elements2[numberElements]=-value;
	row2[numberElements++]=jColumn+numberRows;
      }
    }
    start2[iColumn+1]=numberElements;
  }
  newNumberColumns=numberColumns;
  // pi
  int iRow;
  // Get a row copy in standard format
  CoinPackedMatrix copy;
  copy.reverseOrderedCopyOf(*(this->matrix()));
  // get matrix data pointers
  const int * column = copy.getIndices();
  const CoinBigIndex * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  const double * elementByRow = copy.getElements();
  for (iRow=0;iRow<numberRows;iRow++) {
    // should look at rows to get correct bounds
    columnLower2[newNumberColumns]=-COIN_DBL_MAX;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    int j;
    for (j=rowStart[iRow];
	 j<rowStart[iRow]+rowLength[iRow];
	 j++) {
      elements2[numberElements]=elementByRow[j];
      row2[numberElements++]=column[j]+numberRows;
    }
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
  }
  // djs and artificials
  double tolerance = dualTolerance();
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // dj
    columnLower2[newNumberColumns]=0.0;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    elements2[numberElements]=1.0;
    row2[numberElements++]=iColumn+numberRows;
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
    // artificial (assuming no bounds)
    if (getStatus(iColumn)==basic||solution[iColumn]>1.0e-7) {
      if (fabs(dj[iColumn])>tolerance) {
	columnLower2[newNumberColumns]=0.0;
	columnUpper2[newNumberColumns]=COIN_DBL_MAX;
	objective2[newNumberColumns]=1.0;
	if (dj[iColumn]>0.0)
	  elements2[numberElements]=1.0;
	else
	  elements2[numberElements]=-1.0;
	row2[numberElements++]=iColumn+numberRows;
	newNumberColumns++;
	start2[newNumberColumns]=numberElements;
      }
    } else if (dj[iColumn]<-tolerance) {
      columnLower2[newNumberColumns]=0.0;
      columnUpper2[newNumberColumns]=COIN_DBL_MAX;
      objective2[newNumberColumns]=1.0;
      elements2[numberElements]=-1.0;
      row2[numberElements++]=iColumn+numberRows;
      newNumberColumns++;
      start2[newNumberColumns]=numberElements;
    }
  }
  // Create model
  ClpSimplex model2;
  model2.loadProblem(newNumberColumns,newNumberRows,
		     start2,row2, elements2,
		     columnLower2,columnUpper2,
		     objective2,
		     rowLower2,rowUpper2);
  model2.allSlackBasis();
  // Move solution across
  double * solution2 = model2.primalColumnSolution();
  memcpy(solution2,solution,numberColumns*sizeof(double));
  memcpy(solution2+numberColumns,pi,numberRows*sizeof(double));
  for (iRow=0;iRow<numberRows;iRow++) 
    model2.setRowStatus(iRow,getRowStatus(iRow));
  // djs and artificials
  newNumberColumns = numberRows+numberColumns;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    model2.setStatus(iColumn,getStatus(iColumn));
    if (getStatus(iColumn)==basic) {
      solution2[newNumberColumns++]=0.0;
      if (fabs(dj[iColumn])>tolerance) {
	solution2[newNumberColumns++]=fabs(dj[iColumn]);
      }
    } else if (dj[iColumn]<-tolerance) {
      solution2[newNumberColumns++]=0.0;
      solution2[newNumberColumns++]=-dj[iColumn];
    } else {
      solution2[newNumberColumns++]=dj[iColumn];
    }
  }
  memset(model2.primalRowSolution(),0,newNumberRows*sizeof(double));
  model2.times(1.0,model2.primalColumnSolution(),
	       model2.primalRowSolution());
  // solve
  model2.primal(1);
  return 0;
}
  

