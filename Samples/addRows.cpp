// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.


// This is a simple example to create a model by row
#include "ClpSimplex.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include <iomanip>
#include <cassert>

int main (int argc, const char *argv[])
{
  // Empty model
  ClpSimplex  model;
  
  // Objective - just nonzeros
  int objIndex[]={0,2};
  double objValue[]={1.0,4.0};
  // Upper bounds - as dense vector
  double upper[]={2.0,COIN_DBL_MAX,4.0};

  // Create space for 3 columns
  model.resize(0,3);
  // Fill in 
  int i;
  // Virtuous way
  // First objective
  for (i=0;i<2;i++)
    model.setObjectiveCoefficient(objIndex[i],objValue[i]);
  // Now bounds (lower will be zero by default but do again)
  for (i=0;i<3;i++) {
    model.setColumnLower(i,0.0);
    model.setColumnUpper(i,upper[i]);
  }
  /*
    We could also have done in non-virtuous way e.g.
    double * objective = model.objective();
    and then set directly
  */
  // Faster to add rows all at once - but this is easier to show
  // Now add row 1 as >= 2.0
  int row1Index[] = {0,2};
  double row1Value[]={1.0,1.0};
  model.addRow(2,row1Index,row1Value,
               2.0,COIN_DBL_MAX);
  // Now add row 2 as == 1.0
  int row2Index[] = {0,1,2};
  double row2Value[]={1.0,-5.0,1.0};
  model.addRow(3,row2Index,row2Value,
                 1.0,1.0);
  // solve
  model.dual();

  /*
    Adding one row at a time has a significant overhead so let's
    try a more complicated but faster way

    First time adding in 10000 rows one by one
  */
  model.allSlackBasis();
  ClpSimplex modelSave=model;
  double time1 = CoinCpuTime();
  int k;
  for ( k=0;k<10000;k++) {
    int row2Index[] = {0,1,2};
    double row2Value[]={1.0,-5.0,1.0};
    model.addRow(3,row2Index,row2Value,
                 1.0,1.0);
  }
  printf("Time for 10000 addRow is %g\n",CoinCpuTime()-time1);
  model.dual();
  model=modelSave;
  // Now use build
  CoinBuild buildObject;
  time1 = CoinCpuTime();
  for ( k=0;k<10000;k++) {
    int row2Index[] = {0,1,2};
    double row2Value[]={1.0,-5.0,1.0};
    buildObject.addRow(3,row2Index,row2Value,
                 1.0,1.0);
  }
  model.addRows(buildObject);
  printf("Time for 10000 addRow is %g\n",CoinCpuTime()-time1);
  model.dual();
  // Print column solution
  int numberColumns = model.numberColumns();

  // Alternatively getColSolution()
  double * columnPrimal = model.primalColumnSolution();
  // Alternatively getReducedCost()
  double * columnDual = model.dualColumnSolution();
  // Alternatively getColLower()
  double * columnLower = model.columnLower();
  // Alternatively getColUpper()
  double * columnUpper = model.columnUpper();
  // Alternatively getObjCoefficients()
  double * columnObjective = model.objective();
    
  int iColumn;
  
  std::cout<<"               Primal          Dual         Lower         Upper          Cost"
	   <<std::endl;

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double value;
    std::cout<<std::setw(6)<<iColumn<<" ";
    value = columnPrimal[iColumn];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = columnDual[iColumn];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = columnLower[iColumn];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = columnUpper[iColumn];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;
    value = columnObjective[iColumn];
    if (fabs(value)<1.0e5)
      std::cout<<setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14)<<value;
    else
      std::cout<<setiosflags(std::ios::scientific)<<std::setw(14)<<value;

    std::cout<<std::endl;
  }
  std::cout<<"--------------------------------------"<<std::endl;

  return 0;
}    
