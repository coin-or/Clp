// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpSimplexPrimalQuadratic.hpp"
#include "ClpPrimalQuadraticDantzig.hpp"
#include "ClpFactorization.hpp"
#include "ClpNonLinearCost.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinMpsIO.hpp"
#include "ClpPrimalColumnPivot.hpp"
#include "ClpMessage.hpp"
#include <cfloat>
#include <cassert>
#include <string>
#include <stdio.h>
#include <iostream>
class tempMessage :
   public CoinMessageHandler {

public:
  virtual int print() ;
  tempMessage(ClpSimplex * model);
  ClpSimplex * model_;
};

// Constructor with pointer to model
tempMessage::tempMessage(ClpSimplex * model)
  : CoinMessageHandler(),
    model_(model)
{
}
int
tempMessage::print()
{
  static int numberFeasible=0;
  if (currentSource()=="Clp") {
    if (currentMessage().externalNumber()==5) { 
      if (!numberFeasible&&!model_->nonLinearCost()->numberInfeasibilities()) {
	numberFeasible++;
	model_->setMaximumIterations(0);
      }
    }  
  }
  return CoinMessageHandler::print();
}

// A sequential LP method
int 
ClpSimplexPrimalQuadratic::primalSLP(int numberPasses, double deltaTolerance)
{
  // Are we minimizing or maximizing
  double whichWay=optimizationDirection();
  // This is as a user would see

  int numberColumns = this->numberColumns();
  int numberRows = this->numberRows();
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
  if (!this->status()||numberPrimalInfeasibilities())
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
    if (solution[iColumn]<trueLower[jNon])
      solution[iColumn]=trueLower[jNon];
    else if (solution[iColumn]>trueUpper[jNon])
      solution[iColumn]=trueUpper[jNon];
  }
  int iPass;
  double lastObjective=1.0e31;
  double * saveSolution = new double [numberColumns];
  double * savePi = new double [numberRows];
  unsigned char * saveStatus = new unsigned char[numberRows+numberColumns];
  double targetDrop=1.0e31;
  double objectiveOffset;
  getDblParam(ClpObjOffset,objectiveOffset);
  // 1 bound up, 2 up, -1 bound down, -2 down, 0 no change
  for (iPass=0;iPass<3;iPass++) {
    last[iPass]=new int[numberNonLinearColumns];
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) 
      last[iPass][jNon]=0;
  }
  // goodMove +1 yes, 0 no, -1 last was bad - just halve gaps, -2 do nothing
  int goodMove=-2;
  char * statusCheck = new char[numberColumns];
  for (iPass=0;iPass<numberPasses;iPass++) {
    // redo objective
    double offset=0.0;
    double objValue=-objectiveOffset;
    double lambda=-1.0;
    if (goodMove>=0) {
      // get best interpolation 
      double coeff0=-objectiveOffset,coeff1=0.0,coeff2=0.0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	coeff0 += saveObjective[iColumn]*solution[iColumn];
	coeff1 += saveObjective[iColumn]*(saveSolution[iColumn]-solution[iColumn]);
      }
      for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
	iColumn=listNonLinearColumn[jNon];
	double valueI = solution[iColumn];
	double valueISave = saveSolution[iColumn];
	int j;
	for (j=columnQuadraticStart[iColumn];
	     j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	  int jColumn = columnQuadratic[j];
	  double valueJ = solution[jColumn];
	  double valueJSave = saveSolution[jColumn];
	  double elementValue = 0.5*quadraticElement[j];
	  coeff0 += valueI*valueJ*elementValue;
	  coeff1 += (valueI*valueJSave+valueISave*valueJ-2.0*valueI*valueJ)*elementValue;
	  coeff2 += (valueISave*valueJSave+valueI*valueJ-valueISave*valueJ-valueI*valueJSave)*elementValue;
	}
      }
      double lambdaValue;
      if (fabs(coeff2)<1.0e-9) {
	if (coeff1+coeff2>=0.0) 
	  lambda = 0.0;
	else
	  lambda = 1.0;
      } else {
	lambda = -(0.5*coeff1)/coeff2;
	if (lambda>1.0||lambda<0.0) {
	  if (coeff1+coeff2>=0.0) 
	    lambda = 0.0;
	  else
	    lambda = 1.0;
	}
      }
      lambdaValue = lambda*lambda*coeff2+lambda*coeff1+coeff0;
      printf("coeffs %g %g %g - lastobj %g\n",coeff0,coeff1,coeff2,lastObjective);
      printf("obj at saved %g, obj now %g zero deriv at %g - value %g\n",
	     coeff0+coeff1+coeff2,coeff0,lambda,lambdaValue);
      if (!iPass) lambda=0.0;
      if (lambda>0.0&&lambda<=1.0) {
	// update solution
	for (iColumn=0;iColumn<numberColumns;iColumn++) 
	  solution[iColumn] = lambda * saveSolution[iColumn] 
	    + (1.0-lambda) * solution[iColumn];
	if (lambda>0.999) {
	  memcpy(this->dualRowSolution(),savePi,numberRows*sizeof(double));
	  memcpy(status_,saveStatus,numberRows+numberColumns);
	}
	if (lambda>0.99999&&fabs(coeff1+coeff2)>1.0e-2) {
	  // tighten all
	  goodMove=-1;
	}
      }
    }
    memcpy(objective,saveObjective,numberColumns*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      objValue += objective[iColumn]*solution[iColumn];
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      if (getColumnStatus(iColumn)==basic) {
	if (solution[iColumn]<columnLower[iColumn]+1.0e-8)
	  statusCheck[iColumn]='l';
	else if (solution[iColumn]>columnUpper[iColumn]-1.0e-8)
	  statusCheck[iColumn]='u';
	else
	  statusCheck[iColumn]='B';
      } else {
	if (solution[iColumn]<columnLower[iColumn]+1.0e-8)
	  statusCheck[iColumn]='L';
	else
	  statusCheck[iColumn]='U';
      }
      double valueI = solution[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution[jColumn];
	double elementValue = quadraticElement[j];
	elementValue *= 0.5;
	objValue += valueI*valueJ*elementValue;
	offset += valueI*valueJ*elementValue;
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
    // goodMove +1 yes, 0 no, -1 last was bad - just halve gaps, -2 do nothing
    double maxDelta=0.0;
    if (goodMove>=0) {
      if (objValue<=lastObjective) 
	goodMove=1;
      else
	goodMove=0;
    } else {
      maxDelta=1.0e10;
    }
    double maxGap=0.0;
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      maxDelta = max(maxDelta,
		     fabs(solution[iColumn]-saveSolution[iColumn]));
      if (goodMove>0) {
	if (last[0][jNon]*last[1][jNon]<0) {
	  // halve
	  trust[jNon] *= 0.5;
	} else {
	  if (last[0][jNon]==last[1][jNon]&&
	      last[0][jNon]==last[2][jNon])
	    trust[jNon] *= 1.5; 
	}
      } else if (goodMove!=-2&&trust[jNon]>10.0*deltaTolerance) {
	trust[jNon] *= 0.2;
      }
      maxGap = max(maxGap,trust[jNon]);
    }
    std::cout<<"largest gap is "<<maxGap<<std::endl;
    if (iPass>10000) {
      for (jNon=0;jNon<numberNonLinearColumns;jNon++) 
	trust[jNon] *=0.0001;
    }
    if (goodMove>0) {
      double drop = lastObjective-objValue;
      std::cout<<"True drop was "<<drop<<std::endl;
      std::cout<<"largest delta is "<<maxDelta<<std::endl;
      if (maxDelta<deltaTolerance&&drop<1.0e-4&&goodMove&&lambda<0.99999) {
	std::cout<<"Exiting"<<std::endl;
	break;
      }
    }
    if (!iPass)
      goodMove=1;
    targetDrop=0.0;
    double * r = this->dualColumnSolution();
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      columnLower[iColumn]=max(solution[iColumn]
			       -trust[jNon],
			       trueLower[jNon]);
      columnUpper[iColumn]=min(solution[iColumn]
			       +trust[jNon],
			       trueUpper[jNon]);
    }
    if (iPass) {
      // get reduced costs
      this->matrix()->transposeTimes(savePi,
				     this->dualColumnSolution());
      for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
	iColumn=listNonLinearColumn[jNon];
	double dj = objective[iColumn]-r[iColumn];
	r[iColumn]=dj;
	if (dj<0.0) 
	  targetDrop -= dj*(columnUpper[iColumn]-solution[iColumn]);
	else
	  targetDrop -= dj*(columnLower[iColumn]-solution[iColumn]);
      }
    } else {
      memset(r,0,numberColumns*sizeof(double));
    }
#if 0
    for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
      iColumn=listNonLinearColumn[jNon];
      if (statusCheck[iColumn]=='L'&&r[iColumn]<-1.0e-4) {
	columnLower[iColumn]=max(solution[iColumn],
				 trueLower[jNon]);
	columnUpper[iColumn]=min(solution[iColumn]
				 +trust[jNon],
				 trueUpper[jNon]);
      } else if (statusCheck[iColumn]=='U'&&r[iColumn]>1.0e-4) {
	columnLower[iColumn]=max(solution[iColumn]
				 -trust[jNon],
				 trueLower[jNon]);
	columnUpper[iColumn]=min(solution[iColumn],
				 trueUpper[jNon]);
      } else {
	columnLower[iColumn]=max(solution[iColumn]
				 -trust[jNon],
				 trueLower[jNon]);
	columnUpper[iColumn]=min(solution[iColumn]
				 +trust[jNon],
				 trueUpper[jNon]);
      }
    }
#endif
    if (goodMove) {
      memcpy(saveSolution,solution,numberColumns*sizeof(double));
      memcpy(savePi,this->dualRowSolution(),numberRows*sizeof(double));
      memcpy(saveStatus,status_,numberRows+numberColumns);
      
      std::cout<<"Pass - "<<iPass
	       <<", target drop is "<<targetDrop
	       <<std::endl;
      lastObjective = objValue;
      if (targetDrop<1.0e-5&&goodMove&&iPass) {
	printf("Exiting on target drop %g\n",targetDrop);
	break;
      }
      {
	double * r = this->dualColumnSolution();
	for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
	  iColumn=listNonLinearColumn[jNon];
	  printf("Trust %d %g - solution %d %g obj %g dj %g state %c - bounds %g %g\n",
		 jNon,trust[jNon],iColumn,solution[iColumn],objective[iColumn],
		 r[iColumn],statusCheck[iColumn],columnLower[iColumn],
		 columnUpper[iColumn]);
	}
      }
      setLogLevel(63);
      this->scaling(false);
      this->primal(1);
      if (this->status()) {
	CoinMpsIO writer;
	writer.setMpsData(*matrix(), COIN_DBL_MAX,
			  getColLower(), getColUpper(),
			  getObjCoefficients(),
			  (const char*) 0 /*integrality*/,
			  getRowLower(), getRowUpper(),
			  NULL,NULL);
	writer.writeMps("xx.mps");
      }
      assert (!this->status()); // must be feasible
      goodMove=1;
    } else {
      // bad pass - restore solution
      printf("Backtracking\n");
      memcpy(solution,saveSolution,numberColumns*sizeof(double));
      memcpy(this->dualRowSolution(),savePi,numberRows*sizeof(double));
      memcpy(status_,saveStatus,numberRows+numberColumns);
      iPass--;
      goodMove=-1;
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
  delete [] statusCheck;
  delete [] savePi;
  delete [] saveStatus;
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
  // redo values
  setDblParam(ClpObjOffset,objectiveOffset);
  objectiveValue_ += optimizationDirection_*offset;
  for (jNon=0;jNon<numberNonLinearColumns;jNon++) {
    iColumn=listNonLinearColumn[jNon];
    columnLower[iColumn]= trueLower[jNon];
    columnUpper[iColumn]= trueUpper[jNon];
  }
  memcpy(objective,saveObjective,numberColumns*sizeof(double));
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
// Dantzig's method
int 
ClpSimplexPrimalQuadratic::primalQuadratic(int phase)
{
#if 0
  // Dantzig's method looks easier - lets try that

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
#if 0
  // deliberate bad solution
  // Change to use phase
  //double * saveO = new double[numberColumns];
  //memcpy(saveO,objective,numberColumns*sizeof(double));
  //memset(objective,0,numberColumns*sizeof(double));
  tempMessage messageHandler(this);;
  passInMessageHandler(&messageHandler);
  factorization()->maximumPivots(1);
  primal();
  CoinMessageHandler handler2;
  passInMessageHandler(&handler2);
  factorization()->maximumPivots(100);
  setMaximumIterations(1000);
#endif
  //memcpy(objective,saveO,numberColumns*sizeof(double));
  //printf("For testing - deliberate bad solution\n");
  //columnUpper[0]=0.0;
  //columnLower[0]=0.0;
  //quadraticSLP(50,1.0e-4);
  //primal(1);
  //columnUpper[0]=COIN_DBL_MAX;
  
  // Create larger problem
  // First workout how many rows extra
  ClpQuadraticInfo info(this);
  int numberQuadratic = info.numberQuadraticColumns();
  int newNumberRows = numberRows+numberQuadratic;
  int newNumberColumns = numberColumns + numberRows + numberQuadratic;
  int numberElements = 2*matrix->getNumElements()
    +2*quadratic->getNumElements()
    + numberQuadratic;
  double * elements2 = new double[numberElements];
  int * start2 = new int[newNumberColumns+1];
  int * row2 = new int[numberElements];
  double * objective2 = new double[newNumberColumns];
  double * columnLower2 = new double[newNumberColumns];
  double * columnUpper2 = new double[newNumberColumns];
  double * rowLower2 = new double[newNumberRows];
  double * rowUpper2 = new double[newNumberRows];
  const int * which = info.quadraticSequence();
  const int * back = info.backSequence();
  memcpy(rowLower2,rowLower,numberRows*sizeof(double));
  memcpy(rowUpper2,rowUpper,numberRows*sizeof(double));
  int iRow;
  for (iRow=0;iRow<numberQuadratic;iRow++) {
    double cost = objective[back[iRow]];
    rowLower2[iRow+numberRows]=cost;
    rowUpper2[iRow+numberRows]=cost;
  }
  memset(objective2,0,newNumberColumns*sizeof(double));
  // Get a row copy of quadratic objective in standard format
  CoinPackedMatrix copyQ;
  copyQ.reverseOrderedCopyOf(*quadratic);
  const int * columnQ = copyQ.getIndices();
  const CoinBigIndex * rowStartQ = copyQ.getVectorStarts();
  const int * rowLengthQ = copyQ.getVectorLengths(); 
  const double * elementByRowQ = copyQ.getElements();
  // Move solution across
  double * solution2 = new double[newNumberColumns];
  memset(solution2,0,newNumberColumns*sizeof(double));
  newNumberColumns=0;
  numberElements=0;
  start2[0]=0;
  // x
  memcpy(dj,objective,numberColumns*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // Original matrix
    columnLower2[iColumn]=columnLower[iColumn];
    columnUpper2[iColumn]=columnUpper[iColumn];
    solution2[iColumn]=solution[iColumn];
    int j;
    for (j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];
	 j++) {
      elements2[numberElements]=element[j];
      row2[numberElements++]=row[j];
    }
    if (which[iColumn]>=0) {
      // Quadratic and modify djs
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];
	   j++) {
	int jColumn = columnQuadratic[j];
	double value = quadraticElement[j];
	if (iColumn!=jColumn) 
	  value *= 0.5;
	dj[jColumn] += solution[iColumn]*value;
	elements2[numberElements]=-value;
	row2[numberElements++]=which[jColumn]+numberRows;
      }
      for (j=rowStartQ[iColumn];
	   j<rowStartQ[iColumn]+rowLengthQ[iColumn];
	   j++) {
	int jColumn = columnQ[j];
	double value = elementByRowQ[j];
	if (iColumn!=jColumn) { 
	  value *= 0.5;
	  dj[jColumn] += solution[iColumn]*value;
	  elements2[numberElements]=-value;
	  row2[numberElements++]=which[jColumn]+numberRows;
	}
      }
    }
    start2[iColumn+1]=numberElements;
  }
  newNumberColumns=numberColumns;
  // pi
  // Get a row copy in standard format
  CoinPackedMatrix copy;
  copy.reverseOrderedCopyOf(*(this->matrix()));
  // get matrix data pointers
  const int * column = copy.getIndices();
  const CoinBigIndex * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  const double * elementByRow = copy.getElements();
  for (iRow=0;iRow<numberRows;iRow++) {
    solution2[newNumberColumns]=pi[iRow];
    double value = pi[iRow];
    columnLower2[newNumberColumns]=-COIN_DBL_MAX;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    int j;
    for (j=rowStart[iRow];
	 j<rowStart[iRow]+rowLength[iRow];
	 j++) {
      double elementValue=elementByRow[j];
      int jColumn = column[j];
      elements2[numberElements]=elementValue;
      row2[numberElements++]=jColumn+numberRows;
      dj[jColumn]-= value*elementValue;
    }
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
  }
  // djs 
  for (iColumn=0;iColumn<numberQuadratic;iColumn++) {
    columnLower2[newNumberColumns]=-COIN_DBL_MAX;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    solution2[newNumberColumns]=dj[iColumn];
    elements2[numberElements]=1.0;
    row2[numberElements++]=back[iColumn]+numberRows;
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
  }
  // Create model 
  ClpSimplex model2(*this);
  model2.resize(0,0);
  model2.loadProblem(newNumberColumns,newNumberRows,
		     start2,row2, elements2,
		     columnLower2,columnUpper2,
		     objective2,
		     rowLower2,rowUpper2);
  delete [] objective2;
  delete [] rowLower2;
  delete [] rowUpper2;
  delete [] columnLower2;
  delete [] columnUpper2;
  // Now create expanded quadratic objective for use in primalRow
  // Later pack down in some way
  start2[0]=0;
  numberElements=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // Quadratic
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];
	 j++) {
      int jColumn = columnQuadratic[j];
      double value = quadraticElement[j];
      if (iColumn!=jColumn) 
	value *= 0.5;
      elements2[numberElements]=value;
      row2[numberElements++]=jColumn;
    }
    for (j=rowStartQ[iColumn];
	 j<rowStartQ[iColumn]+rowLengthQ[iColumn];
	 j++) {
      int jColumn = columnQ[j];
      double value = elementByRowQ[j];
      if (iColumn!=jColumn) { 
	value *= 0.5;
	elements2[numberElements]=value;
	row2[numberElements++]=jColumn;
      }
    }
    start2[iColumn+1]=numberElements;
  }
  // and pad
  for (;iColumn<newNumberColumns;iColumn++)
    start2[iColumn+1]=numberElements;
  // Load up objective
  model2.loadQuadraticObjective(newNumberColumns,start2,row2,elements2);
  delete [] start2;
  delete [] row2;
  delete [] elements2;
  model2.allSlackBasis();
  model2.scaling(false);
  model2.setLogLevel(this->logLevel());
  // Move solution across
  memcpy(model2.primalColumnSolution(),solution2,
	 newNumberColumns*sizeof(double));
  columnLower2 = model2.columnLower();
  columnUpper2 = model2.columnUpper();
  delete [] solution2;
  solution2 = model2.primalColumnSolution();
  // Compute row activities and check feasible
  double * rowSolution2 = model2.primalRowSolution();
  memset(rowSolution2,0,newNumberRows*sizeof(double));
  model2.times(1.0,solution2,rowSolution2);
  rowLower2 = model2.rowLower();
  rowUpper2 = model2.rowUpper();
#if 0
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    Status xStatus = getColumnStatus(iColumn);
    bool isSuperBasic;
    int iS = iColumn+newNumberRows;
    double value = solution2[iS];
    if (fabs(value)>dualTolerance_)
      isSuperBasic=true;
    else
      isSuperBasic=false;
    // For moment take all x out of basis
    // Does not seem right
    isSuperBasic=true;
    model2.setColumnStatus(iColumn,xStatus);
    if (xStatus==basic) {
      if (!isSuperBasic) {
	model2.setRowStatus(numberRows+iColumn,basic);
	model2.setColumnStatus(iS,superBasic);
      } else {
	model2.setRowStatus(numberRows+iColumn,isFixed);
	model2.setColumnStatus(iS,basic);
	model2.setColumnStatus(iColumn,superBasic);
      }
    } else {
      model2.setRowStatus(numberRows+iColumn,isFixed);
      model2.setColumnStatus(iS,basic);
    }
  }
  for (iRow=0;iRow<numberRows;iRow++) {
    Status rowStatus = getRowStatus(iRow);
    model2.setRowStatus(iRow,rowStatus);
    if (rowStatus!=basic) {
      model2.setColumnStatus(iRow+numberColumns,basic); // make dual basic
    }
    assert (rowSolution2[iRow]>=rowLower2[iRow]-primalTolerance_);
    assert (rowSolution2[iRow]<=rowUpper2[iRow]+primalTolerance_);
  }
  // why ?? - take duals out and adjust
  for (iRow=0;iRow<numberRows;iRow++) {
    model2.setRowStatus(iRow,basic);
    model2.setColumnStatus(iRow+numberColumns,superBasic);
    solution2[iRow+numberColumns]=0.0;
  }
#else
  for (iRow=0;iRow<numberRows;iRow++) {
    assert (rowSolution2[iRow]>=rowLower2[iRow]-primalTolerance_);
    assert (rowSolution2[iRow]<=rowUpper2[iRow]+primalTolerance_);
    model2.setRowStatus(iRow,basic);
    model2.setColumnStatus(iRow+numberColumns,superBasic);
    solution2[iRow+numberColumns]=0.0;
  }
  for (iColumn=numberRows+numberColumns;iColumn<newNumberColumns;iColumn++) {
    model2.setColumnStatus(iColumn,basic);
    model2.setRowStatus(iColumn-numberColumns,isFixed);
  }
#endif
  memset(rowSolution2,0,newNumberRows*sizeof(double));
  model2.times(1.0,solution2,rowSolution2);
  for (iColumn=0;iColumn<numberQuadratic;iColumn++) {
    int iS = back[iColumn]+newNumberRows;
    int iRow = iColumn+numberRows;
    double value = rowSolution2[iRow];
    if (value>rowUpper2[iRow]) {
      rowSolution2[iRow] = rowUpper2[iRow];
      solution2[iS]-=value-rowUpper2[iRow];
    } else {
      rowSolution2[iRow] = rowLower2[iRow];
      solution2[iS]-=value-rowLower2[iRow];
    }
  }

  
  // See if any s sub j have wrong sign and/or use djs from infeasibility objective
  double objectiveOffset;
  getDblParam(ClpObjOffset,objectiveOffset);
  double objValue = -objectiveOffset;
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    objValue += objective[iColumn]*solution2[iColumn];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double valueI = solution2[iColumn];
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
      int jColumn = columnQuadratic[j];
      double valueJ = solution2[jColumn];
      double elementValue = quadraticElement[j];
      objValue += 0.5*valueI*valueJ*elementValue;
    }
  }
  printf("Objective value %g\n",objValue);
  for (iColumn=0;iColumn<newNumberColumns;iColumn++) 
    printf("%d %g\n",iColumn,solution2[iColumn]);
#if 1
  CoinMpsIO writer;
  writer.setMpsData(*model2.matrix(), COIN_DBL_MAX,
		    model2.getColLower(), model2.getColUpper(),
		    model2.getObjCoefficients(),
		    (const char*) 0 /*integrality*/,
		    model2.getRowLower(), model2.getRowUpper(),
		    NULL,NULL);
  writer.writeMps("xx.mps");
#endif  
  // Now do quadratic
  // If we did not do Slp we should have primal feasible basic solution
  // Do safe cast as no data
  ClpSimplexPrimalQuadratic * modelPtr = 
    (ClpSimplexPrimalQuadratic *) (&model2);
  ClpPrimalQuadraticDantzig dantzigP(modelPtr,numberRows);
  modelPtr->setPrimalColumnPivotAlgorithm(dantzigP);
  modelPtr->messageHandler()->setLogLevel(63);
  modelPtr->primalQuadratic2(this,phase);
  memcpy(dualRowSolution(),model2.primalColumnSolution()+numberColumns_,numberRows_*sizeof(double));
  memcpy(primalColumnSolution(),model2.primalColumnSolution(),numberColumns_*sizeof(double));
  memset(model2.primalRowSolution(),0,newNumberRows*sizeof(double));
  model2.times(1.0,model2.primalColumnSolution(),
	       model2.primalRowSolution());
  memcpy(dualColumnSolution(),model2.primalColumnSolution()+numberRows_+numberColumns_,
	 numberColumns_*sizeof(double));
  objValue = -objectiveOffset;
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    objValue += objective[iColumn]*solution2[iColumn];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double valueI = solution2[iColumn];
    if (fabs(valueI)>1.0e-5) {
      int djColumn = iColumn+numberRows+numberColumns;
      assert(solution2[djColumn]<1.0e-7);
    }
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
      int jColumn = columnQuadratic[j];
      double valueJ = solution2[jColumn];
      double elementValue = quadraticElement[j];
      objValue += 0.5*valueI*valueJ*elementValue;
    }
  }
  printf("Objective value %g\n",objValue);
  objectiveValue_ = objValue + objectiveOffset;
  return 0;
#else
  // Get a feasible solution 
  if (numberPrimalInfeasibilities())
    primal(1);
  // still infeasible
  if (numberPrimalInfeasibilities())
    return 1;
  ClpQuadraticInfo info;
  ClpSimplexPrimalQuadratic * model2 = makeQuadratic(info);
#if 0
  CoinMpsIO writer;
  writer.setMpsData(*model2->matrix(), COIN_DBL_MAX,
		    model2->getColLower(), model2->getColUpper(),
		    model2->getObjCoefficients(),
		    (const char*) 0 /*integrality*/,
		    model2->getRowLower(), model2->getRowUpper(),
		    NULL,NULL);
  writer.writeMps("xx.mps");
#endif  
  // Now do quadratic
  ClpPrimalQuadraticDantzig dantzigP(model2,numberRows_);
  model2->setPrimalColumnPivotAlgorithm(dantzigP);
  model2->messageHandler()->setLogLevel(63);
  model2->primalQuadratic2(this,phase);
  endQuadratic(model2,info);
  return 0;
#endif
}
int ClpSimplexPrimalQuadratic::primalQuadratic2 (const ClpSimplexPrimalQuadratic * originalModel,int phase)
{

  algorithm_ = +2;

  // save data
  ClpDataSave data = saveData();
  
  // Assume problem is feasible
  // Stuff below will be moved into a class
  int numberXColumns = originalModel->numberColumns();
  int numberXRows = originalModel->numberRows();
  int baseS = numberXColumns+numberXRows;
  int * pair = new int [numberColumns_];
  int iColumn;
  for (iColumn=0;iColumn<numberXColumns;iColumn++) {
    int iS = iColumn+baseS;
    pair[iColumn]=iS;
    pair[iS]=iColumn;
  }
#if 0
  // Throw out all x variables whose dj is nonzero
  for (iColumn=0;iColumn<numberXColumns;iColumn++) {
    if (getColumnStatus(iColumn)==basic) {
      int jColumn = iColumn+numberRows_;
      if (getColumnStatus(jColumn)==basic) 
	setColumnStatus(iColumn,superBasic);
    }
  }  
  int originalNumberRows = originalModel->numberRows();
  for (int i=0;i<originalNumberRows;i++) {
    int jSequence = i+numberXColumns;
    if (getRowStatus(i)==basic) 
      setColumnStatus(jSequence,superBasic);
  }
#endif
  // Save solution
  double * saveSolution = new double [numberRows_+numberColumns_];
  assert (!scalingFlag_);
  memcpy(saveSolution,columnActivity_,numberColumns_*sizeof(double));
  memcpy(saveSolution+numberColumns_,rowActivity_,numberRows_*sizeof(double));
  // initialize - values pass coding and algorithm_ is +1
  if (!startup(1)) {

    // Setup useful stuff
    ClpQuadraticInfo info(originalModel);
    if (phase)
      memcpy(solution_,saveSolution,
	     (numberRows_+numberColumns_)*sizeof(double));
    int lastCleaned=0; // last time objective or bounds cleaned up
    int sequenceIn=-1;
    int crucialSj=-1;
    
    // special nonlinear cost
    delete nonLinearCost_;
    nonLinearCost_ = new ClpNonLinearCost(this,originalModel->numberColumns());
    // Say no pivot has occurred (for steepest edge and updates)
    pivotRow_=-2;
    
    // This says whether to restore things etc
    int factorType=0;
    /*
      Status of problem:
      0 - optimal
      1 - infeasible
      2 - unbounded
      -1 - iterating
      -2 - factorization wanted
      -3 - redo checking without factorization
      -4 - looks infeasible
      -5 - looks unbounded
    */
    while (problemStatus_<0) {
      int iRow,iColumn;
      // clear
      for (iRow=0;iRow<4;iRow++) {
	rowArray_[iRow]->clear();
      }    
      
      for (iColumn=0;iColumn<2;iColumn++) {
	columnArray_[iColumn]->clear();
      }    
      
      // give matrix (and model costs and bounds a chance to be
      // refreshed (normally null)
      matrix_->refresh(this);
      // If we have done no iterations - special
      if (lastGoodIteration_==numberIterations_)
	factorType=3;
      if (phase)
	memcpy(saveSolution,solution_,
	       (numberRows_+numberColumns_)*sizeof(double));
      if (phase&&0) {
	// Clean solution
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  if (getColumnStatus(iColumn)==isFree)
	    solution_[iColumn]=0.0;
	}
      }  
      // may factorize, checks if problem finished
      statusOfProblemInPrimal(lastCleaned,factorType,progress_);
      if (phase)
	memcpy(solution_,saveSolution,
	       (numberRows_+numberColumns_)*sizeof(double));
      
      // Compute objective function from scratch
      const CoinPackedMatrix * quadratic = originalModel->quadraticObjective();
      const int * columnQuadratic = quadratic->getIndices();
      const int * columnQuadraticStart = quadratic->getVectorStarts();
      const int * columnQuadraticLength = quadratic->getVectorLengths();
      const double * quadraticElement = quadratic->getElements();
      const double * originalCost = originalModel->objective();
      objectiveValue_=0.0;
      for (iColumn=0;iColumn<numberXColumns;iColumn++) {
	double valueI = solution_[iColumn];
	objectiveValue_ += valueI*originalCost[iColumn];
	int j;
	for (j=columnQuadraticStart[iColumn];
	     j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	  int jColumn = columnQuadratic[j];
	  double valueJ = solution_[jColumn];
	  double elementValue = 0.5*quadraticElement[j];
	  objectiveValue_ += valueI*valueJ*elementValue;
	}
      }

      // Say good factorization
      factorType=1;
      
      // Say no pivot has occurred (for steepest edge and updates)
      pivotRow_=-2;
      // Check problem phase 
      // We assume all X are feasible
      phase=0;
      if (saveSolution) {
	sequenceIn=-1;
	for (iColumn=0;iColumn<numberXColumns;iColumn++) {
	  double value = solution_[iColumn];
	  double lower = lower_[iColumn];
	  double upper = upper_[iColumn];
	  if (value>lower+primalTolerance_&&value<upper-primalTolerance_) {
	    if (getColumnStatus(iColumn)!=basic) {
	      if (phase!=2) {
		phase=2;
		sequenceIn=iColumn;
	      }
	    }
	  }
	  if (getColumnStatus(iColumn)==basic) {
	    int iS = pair[iColumn];
	    assert (iS>=0&&getColumnStatus(iS)!=basic);
	    if(fabs(solution_[iS])>dualTolerance_) {
	      if (phase==0) {
		phase=1;
		sequenceIn=iS;
	      }
	    }
	  }
	}
	int offset=numberXColumns-numberColumns_;
	for (iColumn=numberColumns_;iColumn<numberColumns_+numberXRows;iColumn++) {
	  double value = solution_[iColumn];
	  double lower = lower_[iColumn];
	  double upper = upper_[iColumn];
	  if (value>lower+primalTolerance_&&value<upper-primalTolerance_) {
	    if (getColumnStatus(iColumn)!=basic) {
	      if (phase!=2) {
		phase=2;
		sequenceIn=iColumn;
	      }
	    }
	  }
	  if (getColumnStatus(iColumn)==basic) {
	    int iS = iColumn+offset;
	    assert (getColumnStatus(iS)!=basic);
	    if(fabs(solution_[iS])>dualTolerance_) {
	      if (phase==0) {
		phase=1;
		sequenceIn=iS;
	      }
	    }
	  }
	}
      }
      if (!phase) {
	delete [] saveSolution;
	saveSolution=NULL;
      }

      // exit if victory declared
      if (!phase&&primalColumnPivot_->pivotColumn(rowArray_[1],
					  rowArray_[2],rowArray_[3],
					  columnArray_[0],
					  columnArray_[1]) < 0) {
	problemStatus_=0;
	break;
      }
      
      // Iterate
      problemStatus_=-1;
      whileIterating(originalModel,sequenceIn,&info,crucialSj,phase);
    }
  }
  // clean up
  delete [] saveSolution;
  delete [] pair;
  finish();
  restoreData(data);
  return problemStatus_;
}
/*
  Reasons to come out:
  -1 iterations etc
  -2 inaccuracy 
  -3 slight inaccuracy (and done iterations)
  -4 end of values pass and done iterations
  +0 looks optimal (might be infeasible - but we will investigate)
  +2 looks unbounded
  +3 max iterations 
*/
int
ClpSimplexPrimalQuadratic::whileIterating(
		      const ClpSimplexPrimalQuadratic * originalModel,
		      int & sequenceIn,
		      ClpQuadraticInfo * info,
		      int & crucialSj,
		      int phase)
{
  checkComplementary (info);

  int returnCode=-1;
  double saveObjective = objectiveValue_;
  int numberXColumns = originalModel->numberColumns();
  int oldSequenceIn=sequenceIn;
  // status stays at -1 while iterating, >=0 finished, -2 to invert
  // status -3 to go to top without an invert
  while (problemStatus_==-1) {
#ifdef CLP_DEBUG
    {
      int i;
      // not [1] as has information
      for (i=0;i<4;i++) {
	if (i!=1)
	  rowArray_[i]->checkClear();
      }    
      for (i=0;i<2;i++) {
	columnArray_[i]->checkClear();
      }    
    }      
#endif
    // choose column to come in
    // can use pivotRow_ to update weights
    // pass in list of cost changes so can do row updates (rowArray_[1])
    // NOTE rowArray_[0] is used by computeDuals which is a 
    // slow way of getting duals but might be used 
    // Initially Dantzig and look at s variables
    // Only do if one not already chosen
    bool cleanupIteration;
    if (phase==2) {
      // values pass
      if (sequenceIn<0) {
	// get next
	int iColumn;
	int iStart = oldSequenceIn+1;
	for (iColumn=iStart;iColumn<numberXColumns;iColumn++) {
	  double value = solution_[iColumn];
	  double lower = lower_[iColumn];
	  double upper = upper_[iColumn];
	  if (value>lower+primalTolerance_&&value<upper-primalTolerance_) {
	    if (getColumnStatus(iColumn)!=basic) {
	      sequenceIn=iColumn;
	      break;
	    }
	  }
	}
	if (sequenceIn<0) {
	  iStart=max(iStart,numberColumns_);
	  int numberXRows = originalModel->numberRows();
	  for (iColumn=numberColumns_;iColumn<numberColumns_+numberXRows;
	       iColumn++) {
	    double value = solution_[iColumn];
	    double lower = lower_[iColumn];
	    double upper = upper_[iColumn];
	    if (value>lower+primalTolerance_&&value<upper-primalTolerance_) {
	      if (getColumnStatus(iColumn)!=basic) {
		sequenceIn=iColumn;
		break;
	      }
	    }
	  }
	}
      }
      oldSequenceIn=sequenceIn;
      sequenceIn_ = sequenceIn;
      cleanupIteration=false;
      dualIn_ = solution_[sequenceIn_+numberRows_]; 
      valueIn_=solution_[sequenceIn_];
      if (dualIn_>0.0)
	directionIn_ = -1;
      else 
	directionIn_ = 1;
    } else {
      if (sequenceIn<0) {
	primalColumn(rowArray_[1],rowArray_[2],rowArray_[3],
		     columnArray_[0],columnArray_[1]);
	cleanupIteration=false;
      } else {
	sequenceIn_ = sequenceIn;
	cleanupIteration=true;
      }
    }
    pivotRow_=-1;
    sequenceOut_=-1;
    rowArray_[1]->clear();
    if (sequenceIn_>=0) {
      sequenceIn=-1;
      // we found a pivot column
      int chosen = sequenceIn_;
      // do second half of iteration
      while (chosen>=0) {
	objectiveValue_=saveObjective;
	returnCode=-1;
	pivotRow_=-1;
	sequenceOut_=-1;
	rowArray_[1]->clear();
	// we found a pivot column
	// update the incoming column
	sequenceIn_=chosen;
	chosen=-1;
	unpack(rowArray_[1]);
	factorization_->updateColumnFT(rowArray_[2],rowArray_[1]);
	if (cleanupIteration) {
	  // move back to a version of primalColumn?
	  valueIn_=solution_[sequenceIn_];
	  // should keep pivot row of crucialSj as well (speed)
	  int iSjRow=-1;
	  {
	    double * work=rowArray_[1]->denseVector();
	    int number=rowArray_[1]->getNumElements();
	    int * which=rowArray_[1]->getIndices();
	    double tj = 0.0;
	    int iIndex;
	    for (iIndex=0;iIndex<number;iIndex++) {
	      int iRow = which[iIndex];
	      double alpha = work[iRow];
	      int iPivot=pivotVariable_[iRow];
	      if (iPivot==crucialSj) {
		tj = alpha;
		iSjRow = iRow;
		double d2 = solution_[crucialSj]/tj;
		// see which way to go
		if (d2>0)
		  dj_[sequenceIn_]= -1.0;
		else
		  dj_[sequenceIn_]= 1.0;
		break;
	      }
	    }
	    if (!tj) {
	      printf("trouble\n");
	      assert (sequenceIn_>numberXColumns&&sequenceIn_<numberColumns_);
	      dj_[sequenceIn_]=solution_[sequenceIn_];
	    //assert(tj);
	    }
	  }

	  dualIn_=dj_[sequenceIn_];
	  if (sequenceIn_>numberXColumns&&
	      sequenceIn_<-numberColumns_) {
	    // We can let flip to 0.0
	    assert (cleanupIteration);
	    if (solution_[sequenceIn_]>0.0) {
	      nonLinearCost_->setBounds(sequenceIn_, 0.0,
					0.5*COIN_DBL_MAX);
	      lower_[sequenceIn_]=0.0;
	      upper_[sequenceIn_]=0.5*COIN_DBL_MAX;
	    } else {
	      nonLinearCost_->setBounds(sequenceIn_, -0.5*COIN_DBL_MAX,
					0.0);
	      lower_[sequenceIn_]=-0.5*COIN_DBL_MAX;
	      upper_[sequenceIn_]=0.0;
	    }
	  } else {
	    // Let go through bounds
	    nonLinearCost_->setBounds(sequenceIn_, -0.5*COIN_DBL_MAX,
				      0.5*COIN_DBL_MAX);
	    lower_[sequenceIn_]=-0.5*COIN_DBL_MAX;
	    upper_[sequenceIn_]=0.5*COIN_DBL_MAX;
	  }
	  lowerIn_=lower_[sequenceIn_];
	  upperIn_=upper_[sequenceIn_];
	  if (dualIn_>0.0)
	    directionIn_ = -1;
	  else 
	    directionIn_ = 1;
	} else {
	  if (sequenceIn_<numberColumns_) {
	    // Set dj as value of slack
	    crucialSj = sequenceIn_+numberRows_; // sj which should go to 0.0
	    //dualIn_=solution_[crucialSj];
	  } else {
	    // Set dj as value of pi
	    crucialSj = sequenceIn_-numberColumns_+numberXColumns; // pi which should go to 0.0
	    //dualIn_=solution_[crucialSj];
	  }
	}
	// save reduced cost
	//double saveDj = dualIn_;
	// do ratio test and re-compute dj
	// Note second parameter long enough for columns
	int result=primalRow(rowArray_[1],rowArray_[3],rowArray_[2],rowArray_[0],
			     info,
			     originalModel,crucialSj,cleanupIteration);
	saveObjective = objectiveValue_;
	if (pivotRow_>=0) {
	  // if stable replace in basis
	  int updateStatus = factorization_->replaceColumn(rowArray_[2],
							   pivotRow_,
							   alpha_);
	  // if no pivots, bad update but reasonable alpha - take and invert
	  if (updateStatus==2&&
	      lastGoodIteration_==numberIterations_&&fabs(alpha_)>1.0e-5)
	    updateStatus=4;
	  if (updateStatus==1||updateStatus==4) {
	    // slight error
	    if (factorization_->pivots()>5||updateStatus==4) {
	      returnCode=-3;
	    }
	  } else if (updateStatus==2) {
	    // major error
	    // better to have small tolerance even if slower
	    factorization_->zeroTolerance(1.0e-15);
	    int maxFactor = factorization_->maximumPivots();
	    if (maxFactor>10) {
	      if (forceFactorization_<0)
		forceFactorization_= maxFactor;
	      forceFactorization_ = max (1,(forceFactorization_>>1));
	    } 
	    // later we may need to unwind more e.g. fake bounds
	    if(lastGoodIteration_ != numberIterations_) {
	      rowArray_[1]->clear();
	      pivotRow_=-1;
	      returnCode=-4;
	      // retry on return
	      sequenceIn = sequenceIn_;
	      break;
	    } else {
	      // need to reject something
	      char x = isColumn(sequenceIn_) ? 'C' :'R';
	      handler_->message(CLP_SIMPLEX_FLAG,messages_)
		<<x<<sequenceWithin(sequenceIn_)
		<<CoinMessageEol;
	      setFlagged(sequenceIn_);
	      lastBadIteration_ = numberIterations_; // say be more cautious
	      rowArray_[1]->clear();
	      pivotRow_=-1;
	      returnCode = -5;
	      break;
	      
	    }
	  } else if (updateStatus==3) {
	    // out of memory
	    // increase space if not many iterations
	    if (factorization_->pivots()<
		0.5*factorization_->maximumPivots()&&
		factorization_->pivots()<200)
	      factorization_->areaFactor(
					 factorization_->areaFactor() * 1.1);
	    returnCode =-2; // factorize now
	  }
	  // here do part of steepest - ready for next iteration
	  primalColumnPivot_->updateWeights(rowArray_[1]);
	} else {
	  if (pivotRow_==-1) {
	    // no outgoing row is valid
	    rowArray_[0]->clear();
	    if (!factorization_->pivots()) {
	      returnCode = 2; //say looks unbounded
	      // do ray
	      primalRay(rowArray_[1]);
	    } else {
	      returnCode = 4; //say looks unbounded but has iterated
	    }
	    break;
	  } else {
	    // flipping from bound to bound
	  }
	}


	// update primal solution

	double objectiveChange=0.0;
	// Cost on pivot row may change - may need to change dualIn
	double oldCost=0.0;
	if (pivotRow_>=0)
	  oldCost = cost(pivotVariable_[pivotRow_]);
	// rowArray_[1] is not empty - used to update djs
	updatePrimalsInPrimal(rowArray_[1],theta_, objectiveChange);
	if (pivotRow_>=0)
	  dualIn_ += (oldCost-cost(pivotVariable_[pivotRow_]));
	double oldValue = valueIn_;
	if (directionIn_==-1) {
	  // as if from upper bound
	  if (sequenceIn_!=sequenceOut_) {
	    // variable becoming basic
	    valueIn_ -= fabs(theta_);
	  } else {
	    valueIn_=lowerIn_;
	  }
	} else {
	  // as if from lower bound
	  if (sequenceIn_!=sequenceOut_) {
	    // variable becoming basic
	    valueIn_ += fabs(theta_);
	  } else {
	    valueIn_=upperIn_;
	  }
	}
	objectiveChange += dualIn_*(valueIn_-oldValue);
	// outgoing
	if (sequenceIn_!=sequenceOut_) {
	  if (directionOut_>0) {
	    valueOut_ = lowerOut_;
	  } else {
	    valueOut_ = upperOut_;
	  }
	  double lowerValue = lower_[sequenceOut_];
	  double upperValue = upper_[sequenceOut_];
	  assert(valueOut_>=lowerValue-primalTolerance_&&
		 valueOut_<=upperValue+primalTolerance_);
	  // may not be exactly at bound and bounds may have changed
	  if (valueOut_<=lowerValue+primalTolerance_) {
	    directionOut_=1;
	  } else if (valueOut_>=upperValue-primalTolerance_) {
	    directionOut_=-1;
	  } else {
	    printf("*** variable wandered off bound %g %g %g!\n",
		   lowerValue,valueOut_,upperValue);
	    if (valueOut_-lowerValue<=upperValue-valueOut_) {
	      valueOut_=lowerValue;
	      directionOut_=1;
	    } else {
	      valueOut_=upperValue;
	      directionOut_=-1;
	    }
	  }
	  solution_[sequenceOut_]=valueOut_;
	  nonLinearCost_->setOne(sequenceOut_,valueOut_);
	}
	// change cost and bounds on incoming if primal
	nonLinearCost_->setOne(sequenceIn_,valueIn_); 
	int whatNext=housekeeping(objectiveChange);
	checkComplementary (info);

	if (whatNext==1) {
	  returnCode =-2; // refactorize
	} else if (whatNext==2) {
	  // maximum iterations or equivalent
	  returnCode=3;
	} else if(numberIterations_ == lastGoodIteration_
		  + 2 * factorization_->maximumPivots()) {
	  // done a lot of flips - be safe
	  returnCode =-2; // refactorize
	}
	// may need to go round again
	cleanupIteration = true;
	double newDj;
	// may not be correct on second time
	if (sequenceIn_<numberXColumns)
	  newDj = solution_[sequenceIn_+numberRows_];
	else
	  newDj = solution_[sequenceIn_-numberColumns_+numberXColumns];
	newDj=1.0; // force
	if (result&&fabs(newDj)>dualTolerance_) {
	  assert(sequenceOut_<numberXColumns||
		 sequenceOut_>=numberColumns_);
	  if (sequenceOut_<numberXColumns) {
	    chosen =sequenceOut_ + numberRows_; // sj variable
	  } else {
	    // Does this mean we can change pi
	    int iRow = sequenceOut_-numberColumns_;
	    if (iRow<numberRows_-numberXColumns) {
	      int iPi = iRow+numberXColumns;
	      printf("pi for row %d is %g\n",
		     iRow,solution_[iPi]);
	      chosen=iPi;
	    } else {
	      printf("Row %d is in column part\n",iRow);
	      abort();
	    }
	  }
	} else {
	  break;
	}
      }
      if (returnCode<-1&&returnCode>-5) {
	problemStatus_=-2; // 
      } else if (returnCode==-5) {
	// something flagged - continue;
      } else if (returnCode==2) {
	problemStatus_=-5; // looks unbounded
      } else if (returnCode==4) {
	problemStatus_=-2; // looks unbounded but has iterated
      } else if (returnCode!=-1) {
	assert(returnCode==3);
	problemStatus_=3;
      }
    } else {
      // no pivot column
#ifdef CLP_DEBUG
      if (handler_->logLevel()&32)
	printf("** no column pivot\n");
#endif
      if (nonLinearCost_->numberInfeasibilities())
	problemStatus_=-4; // might be infeasible 
      returnCode=0;
      break;
    }
  }
  return returnCode;
}
/* 
   Row array has pivot column
   This chooses pivot row.
   For speed, we may need to go to a bucket approach when many
   variables go through bounds
   On exit rhsArray will have changes in costs of basic variables
*/
int 
ClpSimplexPrimalQuadratic::primalRow(CoinIndexedVector * rowArray,
				     CoinIndexedVector * rhsArray,
				     CoinIndexedVector * spareArray,
				     CoinIndexedVector * spareArray2,
				     ClpQuadraticInfo * info,
				     const ClpSimplexPrimalQuadratic * originalModel,
				     int crucialSj,
				     bool cleanupIteration)
{
  int result=1;
  // sequence stays as row number until end
  pivotRow_=-1;
  int numberSwapped=0;
  int numberRemaining=0;

  int numberThru =0; // number gone thru a barrier
  int lastThru =0; // number gone thru a barrier on last time
  
  double totalThru=0.0; // for when variables flip
  double acceptablePivot=1.0e-7;
  if (factorization_->pivots())
    acceptablePivot=1.0e-5; // if we have iterated be more strict
  double bestEverPivot=acceptablePivot;
  int lastPivotRow = -1;
  double lastPivot=0.0;
  double lastTheta=1.0e50;
  int lastNumberSwapped=0;

  // use spareArrays to put ones looked at in
  // First one is list of candidates
  // We could compress if we really know we won't need any more
  // Second array has current set of pivot candidates
  // with a backup list saved in double * part of indexed vector

  // for zeroing out arrays after
  int maximumSwapped=0;
  // pivot elements
  double * spare;
  // indices
  int * index, * indexSwapped;
  int * saveSwapped;
  spareArray->clear();
  spareArray2->clear();
  spare = spareArray->denseVector();
  index = spareArray->getIndices();
  saveSwapped = (int *) spareArray2->denseVector();
  indexSwapped = spareArray2->getIndices();

  // we also need somewhere for effective rhs
  double * rhs=rhsArray->denseVector();

  /*
    First we get a list of possible pivots.  We can also see if the
    problem looks unbounded.

    At first we increase theta and see what happens.  We start
    theta at a reasonable guess.  If in right area then we do bit by bit.
    We save possible pivot candidates

   */

  // do first pass to get possibles 
  // We can also see if unbounded

  double * work=rowArray->denseVector();
  int number=rowArray->getNumElements();
  int * which=rowArray->getIndices();

  // we need to swap sign if coming in from ub
  double way = directionIn_;
  double maximumMovement;
  if (way>0.0) 
    maximumMovement = min(1.0e30,upperIn_-valueIn_);
  else
    maximumMovement = min(1.0e30,valueIn_-lowerIn_);

  int iIndex;

  // Work out coefficients for quadratic term
  // This is expanded one
  const CoinPackedMatrix * quadratic = quadraticObjective();
  const int * columnQuadratic = quadratic->getIndices();
  const int * columnQuadraticStart = quadratic->getVectorStarts();
  const int * columnQuadraticLength = quadratic->getVectorLengths();
  const double * quadraticElement = quadratic->getElements();
  const double * originalCost = originalModel->objective();
  // Use rhsArray
  rhsArray->clear();
  int * index2 = rhsArray->getIndices();
  int numberXColumns = originalModel->numberColumns();
  int number2=0;
  //int numberOriginalRows = originalModel->numberRows();
  // sj 
  int iSjRow=-1;
  double tj = 0.0;
  for (iIndex=0;iIndex<number;iIndex++) {
    int iRow = which[iIndex];
    double alpha = -work[iRow]*way;
    int iPivot=pivotVariable_[iRow];
    if (iPivot<numberXColumns) {
      index2[number2++]=iPivot;
      rhs[iPivot]=alpha;
      //printf("col %d alpha %g solution %g\n",iPivot,alpha,solution_[iPivot]);
    } else {
      //printf("new col %d alpha %g solution %g\n",iPivot,alpha,solution_[iPivot]);
      if (iPivot==crucialSj) {
	tj = alpha;
	iSjRow = iRow;
      }
    }
  }
  // Incoming
  if (sequenceIn_<numberXColumns) {
    index2[number2++]=sequenceIn_;
    rhs[sequenceIn_]=way;
    printf("incoming col %d alpha %g solution %g\n",sequenceIn_,way,valueIn_);
  } else {
    printf("incoming new col %d alpha %g solution %g\n",sequenceIn_,way,valueIn_);
  }
  rhsArray->setNumElements(number2);
  // Change in objective will be theta*coeff1 + theta*theta*coeff2
  double coeff1 = 0.0;
  double coeff2 = 0.0;
  if (numberIterations_>=0||cleanupIteration) {
    for (iIndex=0;iIndex<number2;iIndex++) {
      int iColumn=index2[iIndex];
      //double valueI = solution_[iColumn];
      double alphaI = rhs[iColumn];
      coeff1 += alphaI*originalCost[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution_[jColumn];
	double alphaJ = rhs[jColumn];
	double elementValue = quadraticElement[j];
	coeff1 += (valueJ*alphaI)*elementValue;
	coeff2 += (alphaI*alphaJ)*elementValue;
      }
    }
  } else {
    const int * row = matrix_->getIndices();
    const int * columnStart = matrix_->getVectorStarts();
    const int * columnLength = matrix_->getVectorLengths();
    const double * element = matrix_->getElements();
    int j;
    int jRow = sequenceIn_+originalModel->numberRows();
    printf("sequence in %d, cost %g\n",sequenceIn_,
	   upper_[jRow+numberColumns_]);
    double xx=0.0;
    for (j=0;j<numberColumns_;j++) {
      int k;
      for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	if (row[k]==jRow) {
	  printf ("col %d, el %g, sol %g, contr %g\n",
		  j,element[k],solution_[j],element[k]*solution_[j]);
	  xx+= element[k]*solution_[j];
	}
      }
    }
    printf("sum %g\n",xx);
    for (iIndex=0;iIndex<number2;iIndex++) {
      int iColumn=index2[iIndex];
      double valueI = solution_[iColumn];
      double alphaI = rhs[iColumn];
      //rhs[iColumn]=0.0;
      printf("Column %d, alpha %g sol %g cost %g, contr %g\n",
	     iColumn,alphaI,valueI,originalCost[iColumn],
	     alphaI*originalCost[iColumn]);
      coeff1 += alphaI*originalCost[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution_[jColumn];
	double alphaJ = rhs[jColumn];
	double elementValue = quadraticElement[j];
	printf("j %d alphaJ %g solJ %g el %g, contr %g\n",
	       jColumn,alphaJ,valueJ,elementValue,
	       (valueJ*alphaI)*elementValue);
	coeff1 += (valueJ*alphaI)*elementValue;
	coeff2 += (alphaI*alphaJ)*elementValue;
      }
    }
  }
  coeff2 *= 0.5;
  printf("coefficients %g %g - dualIn %g\n",coeff1,coeff2,dualIn_);
  if (!cleanupIteration) 
    assert (fabs(way*coeff1-dualIn_)<1.0e-4);
  // interesting places are where derivative zero or sj goes to zero
  double d1,d2=1.0e50;
  if (fabs(coeff2)>1.0e-9)
    d1 = - 0.5*coeff1/coeff2;
  else if (coeff1<=1.0e-9)
    d1 = maximumMovement;
  else
    d1 = 0.0;
  if (fabs(tj)<1.0e-7) {
    if (sequenceIn_<numberXColumns)
      printf("column %d is basically linear\n",sequenceIn_);
    //assert(!columnQuadraticLength[sequenceIn_]);
  } else {
    d2 = -solution_[crucialSj]/tj;
    if (d2<0.0) {
      printf("d2 would be negative at %g\n",d2);
      d2=1.0e50;
    }
  }
  printf("derivative zero at %g, sj zero at %g\n",d1,d2);
  if (d1>1.0e10&&d2>1.0e10) {
    // linear variable entering
    // maybe we should have done dual iteration to force sj to 0.0
    printf("linear variable\n");
  }
  maximumMovement = min(maximumMovement,d1);
  maximumMovement = min(maximumMovement,d2);
  d2 = min(d1,d2);
    
  rhsArray->clear();
  double tentativeTheta = maximumMovement;
  double upperTheta = maximumMovement;


  for (iIndex=0;iIndex<number;iIndex++) {

    int iRow = which[iIndex];
    double alpha = work[iRow];
    int iPivot=pivotVariable_[iRow];
    alpha *= way;
    double oldValue = solution(iPivot);
    // get where in bound sequence
    if (alpha>0.0) {
      // basic variable going towards lower bound
      double bound = lower(iPivot);
      oldValue -= bound;
    } else if (alpha<0.0) {
      // basic variable going towards upper bound
      double bound = upper(iPivot);
      oldValue = bound-oldValue;
    }
    double value = oldValue-tentativeTheta*fabs(alpha);
    assert (oldValue>=-primalTolerance_*1.0001);
    if (value<-primalTolerance_) {
      // add to list
      spare[numberRemaining]=alpha;
      rhs[iRow]=oldValue;
      index[numberRemaining++]=iRow;
      double value=oldValue-upperTheta*fabs(alpha);
      if (value<-primalTolerance_&&fabs(alpha)>=acceptablePivot)
	upperTheta = (oldValue+primalTolerance_)/fabs(alpha);
    }
  }

  // we need to keep where rhs non-zeros are
  int numberInRhs=numberRemaining;
  memcpy(rhsArray->getIndices(),index,numberInRhs*sizeof(int));
  rhsArray->setNumElements(numberInRhs);

  theta_=maximumMovement;

  bool goBackOne = false;

  if (numberRemaining) {

    // looks like pivoting
    // now try until reasonable theta
    tentativeTheta = max(10.0*upperTheta,1.0e-7);
    tentativeTheta = min(tentativeTheta,maximumMovement);
    
    // loops increasing tentative theta until can't go through
    
    while (tentativeTheta <= maximumMovement) {
      double thruThis = 0.0;
      
      double bestPivot=acceptablePivot;
      pivotRow_ = -1;
      
      numberSwapped = 0;
      
      upperTheta = maximumMovement;
      
      for (iIndex=0;iIndex<numberRemaining;iIndex++) {

	int iRow = index[iIndex];
	double alpha = spare[iIndex];
	double oldValue = rhs[iRow];
	double value = oldValue-tentativeTheta*fabs(alpha);

	if (value<-primalTolerance_) {
	  // how much would it cost to go thru
	  thruThis += alpha*
	    nonLinearCost_->changeInCost(pivotVariable_[iRow],alpha);
	  // goes on swapped list (also means candidates if too many)
	  indexSwapped[numberSwapped++]=iRow;
	  if (fabs(alpha)>bestPivot) {
	    bestPivot=fabs(alpha);
	    pivotRow_ = iRow;
	    theta_ = oldValue/bestPivot;
	  }
	} else {
	  value = oldValue-upperTheta*fabs(alpha);
	  if (value<-primalTolerance_ && fabs(alpha)>=acceptablePivot) 
	    upperTheta = (oldValue+primalTolerance_)/fabs(alpha);
	}
      }
      
      maximumSwapped = max(maximumSwapped,numberSwapped);

      double dualCheck = - 2.0*coeff2*tentativeTheta - coeff1 - 99999999;
      // but make a bit more pessimistic
      dualCheck=max(dualCheck-100.0*dualTolerance_,0.99*dualCheck);
      if (totalThru+thruThis>=dualCheck) {
	// We should be pivoting in this batch
	// so compress down to this lot

	int saveNumber = numberRemaining;
	numberRemaining=0;
	for (iIndex=0;iIndex<numberSwapped;iIndex++) {
	  int iRow = indexSwapped[iIndex];
	  spare[numberRemaining]=way*work[iRow];
	  index[numberRemaining++]=iRow;
	}
	memset(spare+numberRemaining,0,
	       (saveNumber-numberRemaining)*sizeof(double));
	int iTry;
#define MAXTRY 100
	// first get ratio with tolerance
	for (iTry=0;iTry<MAXTRY;iTry++) {
	  
	  upperTheta=maximumMovement;
	  numberSwapped = 0;
	  
	  for (iIndex=0;iIndex<numberRemaining;iIndex++) {
	    
	    int iRow = index[iIndex];
	    double alpha = fabs(spare[iIndex]);
	    double oldValue = rhs[iRow];
	    double value = oldValue-upperTheta*alpha;
	    
	    if (value<-primalTolerance_ && alpha>=acceptablePivot) 
	      upperTheta = (oldValue+primalTolerance_)/alpha;
	    
	  }
	  
	  // now look at best in this lot
	  bestPivot=acceptablePivot;
	  pivotRow_=-1;
	  for (iIndex=0;iIndex<numberRemaining;iIndex++) {
	    
	    int iRow = index[iIndex];
	    double alpha = spare[iIndex];
	    double oldValue = rhs[iRow];
	    double value = oldValue-upperTheta*fabs(alpha);
	    
	    if (value<=0) {
	      // how much would it cost to go thru
	      totalThru += alpha*
		nonLinearCost_->changeInCost(pivotVariable_[iRow],alpha);
	      // goes on swapped list (also means candidates if too many)
	      indexSwapped[numberSwapped++]=iRow;
	      if (fabs(alpha)>bestPivot) {
		bestPivot=fabs(alpha);
		theta_ = oldValue/bestPivot;
		pivotRow_=iRow;
	      }
	    } else {
	      value = oldValue-upperTheta*fabs(alpha);
	      if (value<-primalTolerance_ && fabs(alpha)>=acceptablePivot) 
		upperTheta = (oldValue+primalTolerance_)/fabs(alpha);
	    }
	  }
	  
	  maximumSwapped = max(maximumSwapped,numberSwapped);
	  if (bestPivot<0.1*bestEverPivot&&
	      bestEverPivot>1.0e-6&&bestPivot<1.0e-3) {
	    // back to previous one
	    goBackOne = true;
	    break;
	  } else if (pivotRow_==-1&&upperTheta>largeValue_) {
	    if (lastPivot>acceptablePivot) {
	      // back to previous one
	      goBackOne = true;
	      //break;
	    } else {
	      // can only get here if all pivots so far too small
	    }
	    break;
	  } else {
	    dualCheck = - 2.0*coeff2*theta_ - coeff1-9999999;
	    if (totalThru>=dualCheck) {
	      break; // no point trying another loop
	    } else {
	      // skip this lot
	      nonLinearCost_->goThru(numberSwapped,way,indexSwapped, work,rhs);
	      lastPivotRow=pivotRow_;
	      lastTheta = theta_;
	      lastThru = numberThru;
	      numberThru += numberSwapped;
	      lastNumberSwapped = numberSwapped;
	      memcpy(saveSwapped,indexSwapped,lastNumberSwapped*sizeof(int));
	      if (bestPivot>bestEverPivot)
		bestEverPivot=bestPivot;
	    }
	  }
	}
	break;
      } else {
	// skip this lot
	nonLinearCost_->goThru(numberSwapped,way,indexSwapped, work,rhs);
	lastPivotRow=pivotRow_;
	lastTheta = theta_;
	lastThru = numberThru;
	numberThru += numberSwapped;
	lastNumberSwapped = numberSwapped;
	memcpy(saveSwapped,indexSwapped,lastNumberSwapped*sizeof(int));
	if (bestPivot>bestEverPivot)
	  bestEverPivot=bestPivot;
	totalThru += thruThis;
	tentativeTheta = 2.0*upperTheta;
      }
    }
    // can get here without pivotRow_ set but with lastPivotRow
    if (goBackOne||(pivotRow_<0&&lastPivotRow>=0)) {
      // back to previous one
      pivotRow_=lastPivotRow;
      theta_ = lastTheta;
	    // undo this lot
      nonLinearCost_->goBack(lastNumberSwapped,saveSwapped,rhs);
      memcpy(indexSwapped,saveSwapped,lastNumberSwapped*sizeof(int));
      numberSwapped = lastNumberSwapped;
    }
  }

  if (pivotRow_>=0) {
    
#define MINIMUMTHETA 1.0e-12
    // will we need to increase tolerance
#ifdef CLP_DEBUG
    bool found=false;
#endif
    double largestInfeasibility = primalTolerance_;
    if (theta_<MINIMUMTHETA) {
      theta_=MINIMUMTHETA;
      for (iIndex=0;iIndex<numberSwapped;iIndex++) {
	int iRow = indexSwapped[iIndex];
#ifdef CLP_DEBUG
	if (iRow==pivotRow_)
	  found=true;
#endif
	largestInfeasibility = max (largestInfeasibility,
				    -(rhs[iRow]-fabs(work[iRow])*theta_));
      }
#ifdef CLP_DEBUG
      assert(found);
      if (largestInfeasibility>primalTolerance_&&(handler_->logLevel()&32))
	printf("Primal tolerance increased from %g to %g\n",
	       primalTolerance_,largestInfeasibility);
#endif
      primalTolerance_ = max(primalTolerance_,largestInfeasibility);
    }
    alpha_ = work[pivotRow_];
    // translate to sequence
    sequenceOut_ = pivotVariable_[pivotRow_];
    valueOut_ = solution(sequenceOut_);
    lowerOut_=lower_[sequenceOut_];
    upperOut_=upper_[sequenceOut_];

    if (way<0.0) 
      theta_ = - theta_;
    double newValue = valueOut_ - theta_*alpha_;
    if (alpha_*way<0.0) {
      directionOut_=-1;      // to upper bound
      if (fabs(theta_)>1.0e-6)
	upperOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
      else
	upperOut_ = newValue;
    } else {
      directionOut_=1;      // to lower bound
      if (fabs(theta_)>1.0e-6)
	lowerOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
      else
	lowerOut_ = newValue;
    }
    dualOut_ = reducedCost(sequenceOut_);
  } else {
    double trueMaximumMovement;
    if (way>0.0) 
      trueMaximumMovement = min(1.0e30,upperIn_-valueIn_);
    else
      trueMaximumMovement = min(1.0e30,valueIn_-lowerIn_);
    if (maximumMovement<1.0e20&&maximumMovement==trueMaximumMovement) {
      // flip
      theta_ = maximumMovement;
      pivotRow_ = -2; // so we can tell its a flip
      result=0;
      sequenceOut_ = sequenceIn_;
      valueOut_ = valueIn_;
      dualOut_ = dualIn_;
      lowerOut_ = lowerIn_;
      upperOut_ = upperIn_;
      alpha_ = 0.0;
      if (way<0.0) {
	directionOut_=1;      // to lower bound
	theta_ = lowerOut_ - valueOut_;
      } else {
	directionOut_=-1;      // to upper bound
	theta_ = upperOut_ - valueOut_;
      }
      // we may still have sj to get rid of
    } else if (fabs(maximumMovement-d2)<dualTolerance_) {
      // sj going to zero
      result=0;
      assert (pivotRow_<0);
      nonLinearCost_->setBounds(crucialSj, 0.0,0.0);
      lower_[crucialSj]=0.0;
      upper_[crucialSj]=0.0;
      setColumnStatus(crucialSj,isFixed);
      pivotRow_ = iSjRow;
      alpha_ = work[pivotRow_];
      // translate to sequence
      sequenceOut_ = pivotVariable_[pivotRow_];
      valueOut_ = solution(sequenceOut_);
      lowerOut_=lower_[sequenceOut_];
      upperOut_=upper_[sequenceOut_];
      theta_ = d2;
      if (way<0.0) 
	theta_ = - theta_;
      double newValue = valueOut_ - theta_*alpha_;
      if (alpha_*way<0.0) {
	directionOut_=-1;      // to upper bound
	if (fabs(theta_)>1.0e-6)
	  upperOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
	else
	  upperOut_ = newValue;
      } else {
	directionOut_=1;      // to lower bound
	if (fabs(theta_)>1.0e-6)
	  lowerOut_ = nonLinearCost_->nearest(sequenceOut_,newValue);
	else
	  lowerOut_ = newValue;
      }
      //????
      dualOut_ = reducedCost(sequenceOut_);
    } else {
      // need to do something
      abort();
    }
  }

  // clear arrays

  memset(spare,0,numberRemaining*sizeof(double));
  memset(saveSwapped,0,maximumSwapped*sizeof(int));

  // put back original bounds etc
  nonLinearCost_->goBackAll(rhsArray);

  rhsArray->clear();
  // Change in objective will be theta*coeff1 + theta*theta*coeff2
  objectiveValue_ += theta_*coeff1+theta_*theta_*coeff2;
  printf("New objective value %g\n",objectiveValue_);
  {
    int iColumn;
    objectiveValue_ =0.0;
    CoinPackedMatrix * quadratic = originalModel->quadraticObjective();
    const int * columnQuadratic = quadratic->getIndices();
    const int * columnQuadraticStart = quadratic->getVectorStarts();
    const int * columnQuadraticLength = quadratic->getVectorLengths();
    const double * quadraticElement = quadratic->getElements();
    int numberColumns = originalModel->numberColumns();
    const double * objective = originalModel->objective();
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      objectiveValue_ += objective[iColumn]*solution_[iColumn];
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double valueI = solution_[iColumn];
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution_[jColumn];
	double elementValue = quadraticElement[j];
	objectiveValue_ += 0.5*valueI*valueJ*elementValue;
      }
    }
    printf("Objective value %g\n",objectiveValue());
  }
  return result;
}
// Just for debug
int 
ClpSimplexPrimalQuadratic::checkComplementary (const 
					       ClpQuadraticInfo * info)
{
  const ClpSimplex * originalModel= info->originalModel();
  int numberXColumns = originalModel->numberColumns();
  int i;
  for (i=0;i<numberXColumns;i++) {
    int jSequence = i+ numberRows_;
    if (getColumnStatus(i)==basic) {
      if (getColumnStatus(jSequence)==basic)
	printf("Struct %d (%g) and %d (%g) both basic\n",
	       i,solution_[i],jSequence,solution_[jSequence]);
    }
  }
  int originalNumberRows = originalModel->numberRows();
  int offset = numberXColumns;
  for (i=0;i<originalNumberRows;i++) {
    int jSequence = i+offset;
    if (getRowStatus(i)==basic) {
      if (getColumnStatus(jSequence)==basic)
	printf("Row %d (%g) and %d (%g) both basic\n",
	       i,solution_[i],jSequence,solution_[jSequence]);
    }
  }
  return 0;
}
  
/* This creates the large version of QP and
      fills in quadratic information
*/
ClpSimplexPrimalQuadratic * 
ClpSimplexPrimalQuadratic::makeQuadratic(ClpQuadraticInfo & info)
{

  // Get list of non linear columns
  CoinPackedMatrix * quadratic = quadraticObjective();
  if (!quadratic||!quadratic->getNumElements()) {
    // no quadratic part
    return NULL;
  }

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

  int iColumn;
  const int * columnQuadratic = quadratic->getIndices();
  const int * columnQuadraticStart = quadratic->getVectorStarts();
  const int * columnQuadraticLength = quadratic->getVectorLengths();
  const double * quadraticElement = quadratic->getElements();
#if 0
  // deliberate bad solution
  // Change to use phase
  //double * saveO = new double[numberColumns];
  //memcpy(saveO,objective,numberColumns*sizeof(double));
  //memset(objective,0,numberColumns*sizeof(double));
  tempMessage messageHandler(this);;
  passInMessageHandler(&messageHandler);
  factorization()->maximumPivots(1);
  primal();
  CoinMessageHandler handler2;
  passInMessageHandler(&handler2);
  factorization()->maximumPivots(100);
  setMaximumIterations(1000);
#endif
  //memcpy(objective,saveO,numberColumns*sizeof(double));
  // Get a feasible solution 
  //printf("For testing - deliberate bad solution\n");
  //columnUpper[0]=0.0;
  //columnLower[0]=0.0;
  //quadraticSLP(50,1.0e-4);
  //primal(1);
  //columnUpper[0]=COIN_DBL_MAX;
  
  // Create larger problem
  // First workout how many rows extra
  info=ClpQuadraticInfo(this);
  int numberQuadratic = info.numberQuadraticColumns();
  int newNumberRows = numberRows+numberQuadratic;
  int newNumberColumns = numberColumns + numberRows + numberQuadratic;
  int numberElements = 2*matrix->getNumElements()
    +2*quadratic->getNumElements()
    + numberQuadratic;
  double * elements2 = new double[numberElements];
  int * start2 = new int[newNumberColumns+1];
  int * row2 = new int[numberElements];
  double * objective2 = new double[newNumberColumns];
  double * columnLower2 = new double[newNumberColumns];
  double * columnUpper2 = new double[newNumberColumns];
  double * rowLower2 = new double[newNumberRows];
  double * rowUpper2 = new double[newNumberRows];
  const int * which = info.quadraticSequence();
  const int * back = info.backSequence();
  memcpy(rowLower2,rowLower,numberRows*sizeof(double));
  memcpy(rowUpper2,rowUpper,numberRows*sizeof(double));
  int iRow;
  for (iRow=0;iRow<numberQuadratic;iRow++) {
    double cost = objective[back[iRow]];
    rowLower2[iRow+numberRows]=cost;
    rowUpper2[iRow+numberRows]=cost;
  }
  memset(objective2,0,newNumberColumns*sizeof(double));
  // Get a row copy of quadratic objective in standard format
  CoinPackedMatrix copyQ;
  copyQ.reverseOrderedCopyOf(*quadratic);
  const int * columnQ = copyQ.getIndices();
  const CoinBigIndex * rowStartQ = copyQ.getVectorStarts();
  const int * rowLengthQ = copyQ.getVectorLengths(); 
  const double * elementByRowQ = copyQ.getElements();
  // Move solution across
  double * solution2 = new double[newNumberColumns];
  memset(solution2,0,newNumberColumns*sizeof(double));
  newNumberColumns=0;
  numberElements=0;
  start2[0]=0;
  // x
  memcpy(dj,objective,numberColumns*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // Original matrix
    columnLower2[iColumn]=columnLower[iColumn];
    columnUpper2[iColumn]=columnUpper[iColumn];
    solution2[iColumn]=solution[iColumn];
    int j;
    for (j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];
	 j++) {
      elements2[numberElements]=element[j];
      row2[numberElements++]=row[j];
    }
    if (which[iColumn]>=0) {
      // Quadratic and modify djs
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];
	   j++) {
	int jColumn = columnQuadratic[j];
	double value = quadraticElement[j];
	if (iColumn!=jColumn) 
	  value *= 0.5;
	dj[jColumn] += solution[iColumn]*value;
	elements2[numberElements]=-value;
	row2[numberElements++]=which[jColumn]+numberRows;
      }
      for (j=rowStartQ[iColumn];
	   j<rowStartQ[iColumn]+rowLengthQ[iColumn];
	   j++) {
	int jColumn = columnQ[j];
	double value = elementByRowQ[j];
	if (iColumn!=jColumn) { 
	  value *= 0.5;
	  dj[jColumn] += solution[iColumn]*value;
	  elements2[numberElements]=-value;
	  row2[numberElements++]=which[jColumn]+numberRows;
	}
      }
    }
    start2[iColumn+1]=numberElements;
  }
  newNumberColumns=numberColumns;
  // pi
  // Get a row copy in standard format
  CoinPackedMatrix copy;
  copy.reverseOrderedCopyOf(*(this->matrix()));
  // get matrix data pointers
  const int * column = copy.getIndices();
  const CoinBigIndex * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  const double * elementByRow = copy.getElements();
  for (iRow=0;iRow<numberRows;iRow++) {
    solution2[newNumberColumns]=pi[iRow];
    double value = pi[iRow];
    columnLower2[newNumberColumns]=-COIN_DBL_MAX;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    int j;
    for (j=rowStart[iRow];
	 j<rowStart[iRow]+rowLength[iRow];
	 j++) {
      double elementValue=elementByRow[j];
      int jColumn = column[j];
      elements2[numberElements]=elementValue;
      row2[numberElements++]=jColumn+numberRows;
      dj[jColumn]-= value*elementValue;
    }
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
  }
  // djs 
  for (iColumn=0;iColumn<numberQuadratic;iColumn++) {
    columnLower2[newNumberColumns]=-COIN_DBL_MAX;
    columnUpper2[newNumberColumns]=COIN_DBL_MAX;
    solution2[newNumberColumns]=dj[iColumn];
    elements2[numberElements]=1.0;
    row2[numberElements++]=back[iColumn]+numberRows;
    newNumberColumns++;
    start2[newNumberColumns]=numberElements;
  }
  // Create model 
  ClpSimplex * model2 = new ClpSimplex(*this);
  model2->resize(0,0);
  model2->loadProblem(newNumberColumns,newNumberRows,
		     start2,row2, elements2,
		     columnLower2,columnUpper2,
		     objective2,
		     rowLower2,rowUpper2);
  delete [] objective2;
  delete [] rowLower2;
  delete [] rowUpper2;
  delete [] columnLower2;
  delete [] columnUpper2;
  // Now create expanded quadratic objective for use in primalRow
  // Later pack down in some way
  start2[0]=0;
  numberElements=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // Quadratic
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];
	 j++) {
      int jColumn = columnQuadratic[j];
      double value = quadraticElement[j];
      if (iColumn!=jColumn) 
	value *= 0.5;
      elements2[numberElements]=value;
      row2[numberElements++]=jColumn;
    }
    for (j=rowStartQ[iColumn];
	 j<rowStartQ[iColumn]+rowLengthQ[iColumn];
	 j++) {
      int jColumn = columnQ[j];
      double value = elementByRowQ[j];
      if (iColumn!=jColumn) { 
	value *= 0.5;
	elements2[numberElements]=value;
	row2[numberElements++]=jColumn;
      }
    }
    start2[iColumn+1]=numberElements;
  }
  // and pad
  for (;iColumn<newNumberColumns;iColumn++)
    start2[iColumn+1]=numberElements;
  // Load up objective
  model2->loadQuadraticObjective(newNumberColumns,start2,row2,elements2);
  delete [] start2;
  delete [] row2;
  delete [] elements2;
  model2->allSlackBasis();
  model2->scaling(false);
  model2->setLogLevel(this->logLevel());
  // Move solution across
  memcpy(model2->primalColumnSolution(),solution2,
	 newNumberColumns*sizeof(double));
  columnLower2 = model2->columnLower();
  columnUpper2 = model2->columnUpper();
  delete [] solution2;
  solution2 = model2->primalColumnSolution();
  // Compute row activities and check feasible
  double * rowSolution2 = model2->primalRowSolution();
  memset(rowSolution2,0,newNumberRows*sizeof(double));
  model2->times(1.0,solution2,rowSolution2);
  rowLower2 = model2->rowLower();
  rowUpper2 = model2->rowUpper();
#if 0
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    Status xStatus = getColumnStatus(iColumn);
    bool isSuperBasic;
    int iS = iColumn+newNumberRows;
    double value = solution2[iS];
    if (fabs(value)>dualTolerance_)
      isSuperBasic=true;
    else
      isSuperBasic=false;
    // For moment take all x out of basis
    // Does not seem right
    isSuperBasic=true;
    model2->setColumnStatus(iColumn,xStatus);
    if (xStatus==basic) {
      if (!isSuperBasic) {
	model2->setRowStatus(numberRows+iColumn,basic);
	model2->setColumnStatus(iS,superBasic);
      } else {
	model2->setRowStatus(numberRows+iColumn,isFixed);
	model2->setColumnStatus(iS,basic);
	model2->setColumnStatus(iColumn,superBasic);
      }
    } else {
      model2->setRowStatus(numberRows+iColumn,isFixed);
      model2->setColumnStatus(iS,basic);
    }
  }
  for (iRow=0;iRow<numberRows;iRow++) {
    Status rowStatus = getRowStatus(iRow);
    model2->setRowStatus(iRow,rowStatus);
    if (rowStatus!=basic) {
      model2->setColumnStatus(iRow+numberColumns,basic); // make dual basic
    }
    assert (rowSolution2[iRow]>=rowLower2[iRow]-primalTolerance_);
    assert (rowSolution2[iRow]<=rowUpper2[iRow]+primalTolerance_);
  }
  // why ?? - take duals out and adjust
  for (iRow=0;iRow<numberRows;iRow++) {
    model2->setRowStatus(iRow,basic);
    model2->setColumnStatus(iRow+numberColumns,superBasic);
    solution2[iRow+numberColumns]=0.0;
  }
#else
  for (iRow=0;iRow<numberRows;iRow++) {
    assert (rowSolution2[iRow]>=rowLower2[iRow]-primalTolerance_);
    assert (rowSolution2[iRow]<=rowUpper2[iRow]+primalTolerance_);
    model2->setRowStatus(iRow,basic);
    model2->setColumnStatus(iRow+numberColumns,superBasic);
    solution2[iRow+numberColumns]=0.0;
  }
  for (iColumn=numberRows+numberColumns;iColumn<newNumberColumns;iColumn++) {
    model2->setColumnStatus(iColumn,basic);
    model2->setRowStatus(iColumn-numberColumns,isFixed);
  }
#endif
  memset(rowSolution2,0,newNumberRows*sizeof(double));
  model2->times(1.0,solution2,rowSolution2);
  for (iColumn=0;iColumn<numberQuadratic;iColumn++) {
    int iS = back[iColumn]+newNumberRows;
    int iRow = iColumn+numberRows;
    double value = rowSolution2[iRow];
    if (value>rowUpper2[iRow]) {
      rowSolution2[iRow] = rowUpper2[iRow];
      solution2[iS]-=value-rowUpper2[iRow];
    } else {
      rowSolution2[iRow] = rowLower2[iRow];
      solution2[iS]-=value-rowLower2[iRow];
    }
  }

  
  // See if any s sub j have wrong sign and/or use djs from infeasibility objective
  double objectiveOffset;
  getDblParam(ClpObjOffset,objectiveOffset);
  double objValue = -objectiveOffset;
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    objValue += objective[iColumn]*solution2[iColumn];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double valueI = solution2[iColumn];
    int j;
    for (j=columnQuadraticStart[iColumn];
	 j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
      int jColumn = columnQuadratic[j];
      double valueJ = solution2[jColumn];
      double elementValue = quadraticElement[j];
      objValue += 0.5*valueI*valueJ*elementValue;
    }
  }
  printf("Objective value %g\n",objValue);
  for (iColumn=0;iColumn<newNumberColumns;iColumn++) 
    printf("%d %g\n",iColumn,solution2[iColumn]);
  return (ClpSimplexPrimalQuadratic *) model2;
}

// This moves solution back and deletes information
int 
ClpSimplexPrimalQuadratic::endQuadratic(ClpSimplexPrimalQuadratic * quadraticModel,
		   ClpQuadraticInfo & info)
{
  memcpy(dualRowSolution(),quadraticModel->primalColumnSolution()+numberColumns_,numberRows_*sizeof(double));
  const double * solution2 = quadraticModel->primalColumnSolution();
  memcpy(primalColumnSolution(),solution2,numberColumns_*sizeof(double));
  memset(quadraticModel->primalRowSolution(),0,
	 quadraticModel->numberRows()*sizeof(double));
  quadraticModel->times(1.0,quadraticModel->primalColumnSolution(),
	       quadraticModel->primalRowSolution());
  memcpy(dualColumnSolution(),
	 quadraticModel->primalColumnSolution()+numberRows_+numberColumns_,
	 numberColumns_*sizeof(double));

  int iColumn;
  double objectiveOffset;
  getDblParam(ClpObjOffset,objectiveOffset);
  double objValue = -objectiveOffset;
  double * objective = this->objective();
  for (iColumn=0;iColumn<numberColumns_;iColumn++) 
    objValue += objective[iColumn]*solution2[iColumn];
  CoinPackedMatrix * quadratic = quadraticObjective();
  if (quadratic) {
    const int * columnQuadratic = quadratic->getIndices();
    const int * columnQuadraticStart = quadratic->getVectorStarts();
    const int * columnQuadraticLength = quadratic->getVectorLengths();
    const double * quadraticElement = quadratic->getElements();
    for (iColumn=0;iColumn<numberColumns_;iColumn++) {
      double valueI = solution2[iColumn];
      if (fabs(valueI)>1.0e-5) {
	int djColumn = iColumn+numberRows_+numberColumns_;
	assert(solution2[djColumn]<1.0e-7);
      }
      int j;
      for (j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution2[jColumn];
	double elementValue = quadraticElement[j];
	objValue += 0.5*valueI*valueJ*elementValue;
      }
    }
    objectiveValue_ = objValue + objectiveOffset;
  }
  printf("Objective value %g\n",objValue);
  return 0;
}

/// Default constructor. 
ClpQuadraticInfo::ClpQuadraticInfo()
  : originalModel_(NULL),
    quadraticSequence_(NULL),
    backSequence_(NULL),
    crucialSj_(-1),
    numberXRows_(-1),
    numberXColumns_(-1),
    numberQuadraticColumns_(0)
{
}
// Constructor from original model
ClpQuadraticInfo::ClpQuadraticInfo(const ClpSimplex * model)
  : originalModel_(model),
    quadraticSequence_(NULL),
    backSequence_(NULL),
    crucialSj_(-1),
    numberXRows_(-1),
    numberXColumns_(-1),
    numberQuadraticColumns_(0)
{
  if (originalModel_) {
    numberXRows_ = originalModel_->numberRows();
    numberXColumns_ = originalModel_->numberColumns();
    quadraticSequence_ = new int[numberXColumns_];
    backSequence_ = new int[numberXColumns_];
    int i;
    numberQuadraticColumns_=numberXColumns_;
    for (i=0;i<numberXColumns_;i++) {
      quadraticSequence_[i]=i;
      backSequence_[i]=i;
    }
  }
}
// Destructor
ClpQuadraticInfo:: ~ClpQuadraticInfo()
{
  delete [] quadraticSequence_;
  delete [] backSequence_;
}
// Copy
ClpQuadraticInfo::ClpQuadraticInfo(const ClpQuadraticInfo& rhs)
  : originalModel_(rhs.originalModel_),
    quadraticSequence_(NULL),
    backSequence_(NULL),
    crucialSj_(rhs.crucialSj_),
    numberXRows_(rhs.numberXRows_),
    numberXColumns_(rhs.numberXColumns_),
    numberQuadraticColumns_(rhs.numberQuadraticColumns_)
{
  if (numberXColumns_) {
    quadraticSequence_ = new int[numberXColumns_];
    memcpy(quadraticSequence_,rhs.quadraticSequence_,
	   numberXColumns_*sizeof(int));
    backSequence_ = new int[numberXColumns_];
    memcpy(backSequence_,rhs.backSequence_,
	   numberXColumns_*sizeof(int));
  }
}
// Assignment
ClpQuadraticInfo & 
ClpQuadraticInfo::operator=(const ClpQuadraticInfo&rhs)
{
  if (this != &rhs) {
    originalModel_ = rhs.originalModel_;
    delete [] quadraticSequence_;
    quadraticSequence_ = NULL;
    delete [] backSequence_;
    backSequence_ = NULL;
    crucialSj_ = rhs.crucialSj_;
    numberXRows_ = rhs.numberXRows_;
    numberXColumns_ = rhs.numberXColumns_;
    numberQuadraticColumns_=rhs.numberQuadraticColumns_;
    if (numberXColumns_) {
      quadraticSequence_ = new int[numberXColumns_];
      memcpy(quadraticSequence_,rhs.quadraticSequence_,
	     numberXColumns_*sizeof(int));
      backSequence_ = new int[numberXColumns_];
      memcpy(backSequence_,rhs.backSequence_,
	     numberXColumns_*sizeof(int));
    }
  }
  return *this;
}
