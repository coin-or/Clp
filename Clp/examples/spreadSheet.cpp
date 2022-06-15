// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/* I found this useful when playing around with small problems - and trying
   to see why preprocessing is not working as well as it should.
   This just makes a simple spreadsheet for small problems. */

#include "ClpSimplex.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinMpsIO.hpp"
int main(int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find MPS file.\n");
    exit(1);
  } else {
    status = model.readMps(argv[1],true);
  }
  if (status) {
    printf("errors on input\n");
    exit(77);
  }
  model.dual(0);
  int numberRows = model.numberRows();
  int numberColumns = model.numberColumns();
  double * columnLower = model.columnLower();
  double * columnUpper = model.columnUpper();
  double * columnActivity = model.primalColumnSolution();
  double * objective = model.objective();
  double * rowLower = model.rowLower();
  double * rowUpper = model.rowUpper();
  double * rowActivity = model.primalRowSolution();
  CoinPackedMatrix * matrixByColumn = model.matrix();
  CoinPackedMatrix matrixByRow(*matrixByColumn);
  matrixByRow.reverseOrdering();
  double *elementByRow = matrixByRow.getMutableElements();
  int *column = matrixByRow.getMutableIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();
  FILE * fp = fopen("model.csv","w");
  if (!fp) {
    printf("Unable to open model.csv for writing\n");
    exit(77);
  }
  fprintf(fp,"LowerRhs,");
  int numberIntegers = 0;
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    fprintf(fp,"%s,",model.columnName(iColumn).c_str());
    if (model.isInteger(iColumn))
      numberIntegers++;
  }
  fprintf(fp,"UpperRhs,RowActivity\n");
  for (int iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]>-1.0e30)
      fprintf(fp,"%g,",rowLower[iRow]);
    else
      fprintf(fp,"-Inf,");
    CoinBigIndex start = rowStart[iRow];
    CoinBigIndex end = start+rowLength[iRow];
    CoinSort_2(column+start,column+end,elementByRow+start);
    int iLast =0;
    for (CoinBigIndex j=start;j<end;j++) {
      int iColumn = column[j];
      while (iLast<iColumn) {
	fprintf(fp,",");
	iLast++;
      }
      iLast++;
      fprintf (fp,"%g,",elementByRow[j]);
    }
    while (iLast<numberColumns) {
      fprintf(fp,",");
      iLast++;
    }
    if (rowUpper[iRow]<1.0e30)
      fprintf(fp,"%g,",rowUpper[iRow]);
    else
      fprintf(fp,"+Inf,");
    fprintf(fp,"%g\n",rowActivity[iRow]);
  }
  const double * arrays[4]={objective,columnLower,columnUpper,columnActivity};
  const char * names[4]={"Objective","Lower","Upper","Activity"};
  for (int iType=0;iType<4;iType++) {
    const double * array = arrays[iType];
    fprintf(fp,"%s",names[iType]);
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (array[iColumn]<-1.0e30)
	fprintf(fp,",-Inf");
      else if (array[iColumn]>1.0e30)
	fprintf(fp,",+Inf");
      else
	fprintf(fp,",%g",array[iColumn]);
    }
    if (iType==3) {
      if (!model.objectiveOffset())
	fprintf(fp,",,%g\n",model.objectiveValue());
      else
	fprintf(fp,",%g,%g\n",model.objectiveOffset(),model.objectiveValue());
    } else {
      fprintf(fp,"\n");
    }
  }
  if (numberIntegers) {
    fprintf(fp,"Integers");
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (model.isInteger(iColumn)) {
	if (columnUpper[iColumn]==1.0&&columnLower[iColumn]==0.0)
	  fprintf(fp,",BV");
	else
	  fprintf(fp,",UI");
      } else {
	fprintf(fp,",");
      }
    }
  }
  fclose(fp);
  return 0;
}
