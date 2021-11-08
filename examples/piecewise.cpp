// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/* This example takes a matrix (modified version of netlib/afiro.mps)
   The first and second four variables can be replaced by piecewise linear.
   Solves it as is and then reformulates  and solves the resulting problem 
   - which should start with optimal value of objective */

#include "ClpSimplex.hpp"
#include "ClpNonLinearCost.hpp"
#include "CoinMpsIO.hpp"
#include <iomanip>

int main(int argc, const char *argv[])
{
     int status;
     CoinMpsIO m;
     if (argc < 2)
          status = m.readMps("modified_afiro.mps", "");
     else
          status = m.readMps(argv[1], "");

     if (status) {
          fprintf(stdout, "Bad readMps %s\n", argv[1]);
          exit(1);
     }

     // Load up model1 - so we can use known good solution
     ClpSimplex model1;
     model1.loadProblem(*m.getMatrixByCol(),
                        m.getColLower(), m.getColUpper(),
                        m.getObjCoefficients(),
                        m.getRowLower(), m.getRowUpper());
     // put in column names
     int numberColumns1 = model1.numberColumns();
     for (int i=0;i<numberColumns1;i++) {
       std::string name = m.columnName(i);
       model1.setColumnName(i,name);
     }
     model1.dual();
     // model with four variables replaced by one (for first 8 variables)
     int delColumns[6] = {1,2,3,5,6,7};
     ClpSimplex model2 = model1;
     model2.deleteColumns(6,delColumns);
     //
     const double * lower1 = model1.columnLower();
     const double * upper1 = model1.columnUpper();
     const double * objective1 = model1.objective();
     double * lower2 = model2.columnLower();
     double * upper2 = model2.columnUpper();
     const double * objective2 = model2.objective();
     double breakpointA[5],slopeA[5];
     double breakpointB[5],slopeB[5];
     double bound;
     // First variable
     breakpointA[0] = lower1[0];
     bound = 0.0;
     for (int i=0;i<4;i++) {
       slopeA[i] = objective1[i];
       if (i)
	 assert(lower1[i]==0.0);
       bound += upper1[i];
       breakpointA[i+1] = bound;
       printf("From %g to %g, slope of A is %g\n",
	      breakpointA[i],breakpointA[i+1],slopeA[i]);
     }
     slopeA[4] = 0.0; // not used - but be virtuous
     // Second variable
     breakpointB[0] = lower1[0];
     bound = 0.0;
     for (int i=0;i<4;i++) {
       slopeB[i] = objective1[i+4];
       if (i)
	 assert(lower1[i+4]==0.0);
       bound += upper1[i+4];
       breakpointB[i+1] = bound;
       printf("From %g to %g, slope of B is %g\n",
	      breakpointB[i],breakpointB[i+1],slopeB[i]);
     }
     slopeB[4] = 0.0; // not used - but be virtuous
     // see what solution is with first try
     upper2[0] = breakpointA[4];
     upper2[1] = breakpointB[4];
     // to save a bit of coding be lazy and create new throwawy model
     ClpSimplex model3 = model2;
     const double * solution1 = model1.primalColumnSolution();
     double * solution2 = model3.primalColumnSolution();
     model3.allSlackBasis();
     solution2[0]=solution1[0]+solution1[1]+solution1[2]+solution1[3];
     solution2[1]=solution1[4]+solution1[5]+solution1[6]+solution1[7];
     model3.primal(1);
     // print solutions
     printf("Variable A %g, %g, %g, %g -> %g\n",
	    solution1[0],solution1[1],solution1[2],solution1[3],
	    solution2[0]);
     printf("Variable B %g, %g, %g, %g -> %g\n",
	    solution1[4],solution1[5],solution1[6],solution1[7],
	    solution2[1]);
     int numberColumns = model2.numberColumns();
     for (int i=2;i<numberColumns;i++)
       printf("Variable %d %s value %g -> %g\n",
	      i,model2.columnName(i).c_str(),solution1[i+6],
	      solution2[i]);
     int * segstart = new int[numberColumns+1];
     double * breakpt  = new double[2*numberColumns+8];
     double * slope  = new double[2*numberColumns+8];
     bool basic = false;
     solution2 = model2.primalColumnSolution();
     segstart[0] = 0;
     // A
     double valueA = 0.0;
     for (int i=0;i<5;i++) {
       breakpt[i] = breakpointA[i];
       slope[i] = slopeA[i];
       if (i<4) {
	 valueA += solution1[i];
	 if (model1.getColumnStatus(i)==ClpSimplex::basic)
	   basic = true;
       }
     }
     solution2[0] = valueA;
     if (basic)
       model2.setColumnStatus(0,ClpSimplex::basic);
     segstart[1] = 5;
     // B
     basic = false;
     double valueB = 0.0;
     for (int i=0;i<5;i++) {
       breakpt[i+5] = breakpointB[i];
       slope[i+5] = slopeB[i];
       if (i<4) {
	 valueB += solution1[i+4];
	 if (model1.getColumnStatus(i+4)==ClpSimplex::basic)
	   basic = true;
       }
     }
     solution2[1] = valueB;
     if (basic)
       model2.setColumnStatus(1,ClpSimplex::basic);
     int segptr = 10;
     segstart[2] = segptr;
     for (int i=2;i<numberColumns;i++) {
       breakpt[segptr] = lower2[i];
       slope[segptr++] = objective2[i];
       breakpt[segptr] = upper2[i]; 
       slope[segptr++] = 0.0;
       segstart[i+1] = segptr;
     }
     /* so at this stage -
	segstart gives where piecewise non linear information is
	for each variable. Normal variables just have two items -
	breakpt has normal lower and upper bounds and first slope item
	has normal cost.
	The piecewise linear variables have more items - the first and last
	breakpt items have lowest possible value and highest possible value
	and slope gives the increasing costin each range.
      */
     model2.scaling(0); // at present has to have scaling off - could fix
     // Create nonlinear objective
     int returnCode =
       model2.createPiecewiseLinearCosts(segstart, breakpt, slope);
     if( returnCode != 0 )  
     {
        printf("Unexpected return code %d from model.createPiecewiseLinearCosts()\n", returnCode);
        return returnCode;
     }

     // delete
     delete [] segstart;
     delete [] breakpt;
     delete [] slope; 
     
     model2.primal(1);
     printf("Variable A %g, %g, %g, %g -> %g\n",
	    solution1[0],solution1[1],solution1[2],solution1[3],
	    solution2[0]);
     printf("Variable B %g, %g, %g, %g -> %g\n",
	    solution1[4],solution1[5],solution1[6],solution1[7],
	    solution2[1]);
     for (int i=2;i<numberColumns;i++)
       printf("Variable %d %s value %g -> %g\n",
	      i,model2.columnName(i).c_str(),solution1[i+6],
	      solution2[i]);
     return 0;
}
