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
  abort();
  return 0;
}
  

