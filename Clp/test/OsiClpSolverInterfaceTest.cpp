// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include "OsiConfig.h"

#include <cassert>
//#include <cstdlib>
//#include <cstdio>
//#include <iostream>

#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinMessage.hpp"
#include "ClpMessage.hpp"
#include "ClpFactorization.hpp"
#include "CoinModel.hpp"
#include "CoinIndexedVector.hpp"

//#############################################################################

class OsiClpMessageTest :
   public CoinMessageHandler {

public:
  virtual int print() ;
  OsiClpMessageTest();
};

OsiClpMessageTest::OsiClpMessageTest() : CoinMessageHandler()
{
}
int
OsiClpMessageTest::print()
{
  if (currentMessage().externalNumber()==0&&currentSource()=="Clp") 
    std::cout<<"This is not actually an advertisement by Dash Associates - just my feeble attempt to test message handling and language - JJHF"<<std::endl;
  else if (currentMessage().externalNumber()==5&&currentSource()=="Osi") 
    std::cout<<"End of search trapped"<<std::endl;
  return CoinMessageHandler::print();
}

//--------------------------------------------------------------------------
// test EKKsolution methods.
int
OsiClpSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{
  
  // Test default constructor
  {
    OsiClpSolverInterface m;
    assert( m.rowsense_==NULL );
    assert( m.rhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.ws_==NULL);
    assert( m.itlimOrig_==9999999);
    assert( m.lastAlgorithm_==0);
    assert( m.integerInformation_==NULL);
  }
  
  
  {    
    CoinRelFltEq eq;
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    
    {
      OsiClpSolverInterface im;    
      
      assert( im.getNumCols() == 0 ); 
      
      assert( im.getModelPtr()!=NULL );
      // Test reset
      im.reset();
      assert( im.rowsense_==NULL );
      assert( im.rhs_==NULL );
      assert( im.rowrange_==NULL );
      assert( im.matrixByRow_==NULL );
      assert( im.ws_==NULL);
      assert( im.itlimOrig_==9999999);
      assert( im.lastAlgorithm_==0);
      assert( im.integerInformation_==NULL);
    }
    
    // Test copy constructor and assignment operator
    {
      OsiClpSolverInterface lhs;
      {      
        OsiClpSolverInterface im(m);        
        
        OsiClpSolverInterface imC1(im);
        assert( imC1.getModelPtr()!=im.getModelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        
        OsiClpSolverInterface imC2(im);
        assert( imC2.getModelPtr()!=im.getModelPtr() );
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        
        assert( imC2.getModelPtr()!=imC1.getModelPtr() );
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope
      
      assert( lhs.getModelPtr() != m.getModelPtr() );
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
    }
    // Test clone
    {
      OsiClpSolverInterface oslSi(m);
      OsiSolverInterface * siPtr = &oslSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiClpSolverInterface * oslClone = dynamic_cast<OsiClpSolverInterface*>(siClone);
      assert( oslClone != NULL );
      assert( oslClone->getModelPtr() != oslSi.getModelPtr() );
      assert( oslClone->getNumRows() == oslSi.getNumRows() );
      assert( oslClone->getNumCols() == m.getNumCols() );
      
      delete siClone;
    }
  
    // test infinity
    {
      OsiClpSolverInterface si;
      assert( eq(si.getInfinity(),OsiClpInfinity));
    }     
    
    // Test setting solution
    {
      OsiClpSolverInterface m1(m);
      int i;

      double * cs = new double[m1.getNumCols()];
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        cs[i] = i + .5;
      m1.setColSolution(cs);
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        assert(m1.getColSolution()[i] == i + .5);
      
      double * rs = new double[m1.getNumRows()];
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        rs[i] = i - .5;
      m1.setRowPrice(rs);
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        assert(m1.getRowPrice()[i] == i - .5);

      delete [] cs;
      delete [] rs;
    }
    
    
    // Test fraction Indices
    {
      OsiClpSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      assert(  fim.isContinuous(0) );
      assert(  fim.isContinuous(1) );
      assert( !fim.isContinuous(2) );
      assert( !fim.isContinuous(3) );
      assert(  fim.isContinuous(4) );
      
      assert( !fim.isInteger(0) );
      assert( !fim.isInteger(1) );
      assert(  fim.isInteger(2) );
      assert(  fim.isInteger(3) );
      assert( !fim.isInteger(4) );
      
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert(  fim.isBinary(2) );
      assert(  fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert( !fim.isIntegerNonBinary(2) );
      assert( !fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );
      // Test fractionalIndices
      {
        double sol[]={2.9,3.0};
        memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        OsiVectorInt fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        sol[0]=5 + .00001/2.;
        sol[1]=8 - .00001/2.;
        memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        sol[0]=5 + .00001*2.;
        sol[1]=8 - .00001*2.;
        memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );
      }


      
      // Change date so column 2 & 3 are integerNonBinary
      double ub[]={5.0,6.0};
      memcpy(fim.modelPtr_->columnUpper()+2,ub,2*sizeof(double));
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert( !fim.isBinary(2) );
      assert( !fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert(  fim.isIntegerNonBinary(2) );
      assert(  fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );
    }
    // Test some catches
      if (!OsiClpHasNDEBUG())
    {
      OsiClpSolverInterface solver;
#ifndef NDEBUG  /* in optimized mode, the following code crashes */
      try {
        solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
        std::cout<<"Correct throw"<<std::endl;
      }
#endif
      std::string fn = mpsDir+"exmip1";
      solver.readMps(fn.c_str(),"mps");
      try {
        solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
        std::cout<<"** Incorrect throw"<<std::endl;
        abort();
      }
#ifndef NDEBUG
      try {
        int index[]={0,20};
        double value[]={0.0,0.0,0.0,0.0};
        solver.setColSetBounds(index,index+2,value);
      }
      catch (CoinError e) {
        std::cout<<"Correct throw"<<std::endl;
      }
#endif
    }
    // Test apply cuts method
    {      
      OsiClpSolverInterface im(m);
      OsiCuts cuts;
      
      // Generate some cuts 
      {
        // Get number of rows and columns in model
        int nr=im.getNumRows();
        int nc=im.getNumCols();
        assert( nr == 5 );
        assert( nc == 8 );
        
        // Generate a valid row cut from thin air
        int c;
        {
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *el = new double[nc];
          for (c=0;c<nc;c++) el[c]=((double)c)*((double)c);
          
          OsiRowCut rc;
          rc.setRow(nc,inx,el);
          rc.setLb(-100.);
          rc.setUb(100.);
          rc.setEffectiveness(22);
          
          cuts.insert(rc);
          delete[]el;
          delete[]inx;
        }
        
        // Generate valid col cut from thin air
        {
          const double * oslColLB = im.getColLower();
          const double * oslColUB = im.getColUpper();
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *lb = new double[nc];
          double *ub = new double[nc];
          for (c=0;c<nc;c++) lb[c]=oslColLB[c]+0.001;
          for (c=0;c<nc;c++) ub[c]=oslColUB[c]-0.001;
          
          OsiColCut cc;
          cc.setLbs(nc,inx,lb);
          cc.setUbs(nc,inx,ub);
          
          cuts.insert(cc);
          delete [] ub;
          delete [] lb;
          delete [] inx;
        }
        
        {
          // Generate a row and column cut which are ineffective
          OsiRowCut * rcP= new OsiRowCut;
          rcP->setEffectiveness(-1.);
          cuts.insert(rcP);
          assert(rcP==NULL);
          
          OsiColCut * ccP= new OsiColCut;
          ccP->setEffectiveness(-12.);
          cuts.insert(ccP);
          assert(ccP==NULL);
        }
        {
          //Generate inconsistent Row cut
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          rc.setLb(3.);
          rc.setUb(4.);
          assert(!rc.consistent());
          cuts.insert(rc);
        }
        {
          //Generate inconsistent col cut
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          cc.setUbs(ne,inx,el);
          assert(!cc.consistent());
          cuts.insert(cc);
        }
        {
          // Generate row cut which is inconsistent for model m
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          assert(rc.consistent());
          assert(!rc.consistent(im));
          cuts.insert(rc);
        }
        {
          // Generate col cut which is inconsistent for model m
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={30};
          double el[ne]={2.0};
          cc.setLbs(ne,inx,el);
          assert(cc.consistent());
          assert(!cc.consistent(im));
          cuts.insert(cc);
        }
        {
          // Generate col cut which is infeasible
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={0};
          double el[ne]={2.0};
          cc.setUbs(ne,inx,el);
          cc.setEffectiveness(1000.);
          assert(cc.consistent());
          assert(cc.consistent(im));
          assert(cc.infeasible(im));
          cuts.insert(cc);
        }
      }
      assert(cuts.sizeRowCuts()==4);
      assert(cuts.sizeColCuts()==5);
      
      OsiSolverInterface::ApplyCutsReturnCode rc = im.applyCuts(cuts);
      assert( rc.getNumIneffective() == 2 );
      assert( rc.getNumApplied() == 2 );
      assert( rc.getNumInfeasible() == 1 );
      assert( rc.getNumInconsistentWrtIntegerModel() == 2 );
      assert( rc.getNumInconsistent() == 2 );
      assert( cuts.sizeCuts() == rc.getNumIneffective() +
        rc.getNumApplied() +
        rc.getNumInfeasible() +
        rc.getNumInconsistentWrtIntegerModel() +
        rc.getNumInconsistent() );
    }
    
    {    
      OsiClpSolverInterface oslSi(m);
      int nc = oslSi.getNumCols();
      int nr = oslSi.getNumRows();
      const double * cl = oslSi.getColLower();
      const double * cu = oslSi.getColUpper();
      const double * rl = oslSi.getRowLower();
      const double * ru = oslSi.getRowUpper();
      assert( nc == 8 );
      assert( nr == 5 );
      assert( eq(cl[0],2.5) );
      assert( eq(cl[1],0.0) );
      assert( eq(cu[1],4.1) );
      assert( eq(cu[2],1.0) );
      assert( eq(rl[0],2.5) );
      assert( eq(rl[4],3.0) );
      assert( eq(ru[1],2.1) );
      assert( eq(ru[4],15.0) );
      
      const double * cs = oslSi.getColSolution();
      assert( eq(cs[0],2.5) );
      assert( eq(cs[7],0.0) );
      
      assert( !eq(cl[3],1.2345) );
      oslSi.setColLower( 3, 1.2345 );
      assert( eq(oslSi.getColLower()[3],1.2345) );
      
      assert( !eq(cu[4],10.2345) );
      oslSi.setColUpper( 4, 10.2345 );
      assert( eq(oslSi.getColUpper()[4],10.2345) );

      double objValue = oslSi.getObjValue();
      assert( eq(objValue,3.5) );

      assert( eq( oslSi.getObjCoefficients()[0],  1.0) );
      assert( eq( oslSi.getObjCoefficients()[1],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[2],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[3],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[4],  2.0) );
      assert( eq( oslSi.getObjCoefficients()[5],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[6],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[7], -1.0) );
    }
    
    // Test matrixByRow method
    { 
      const OsiClpSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByRow();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    
      
      assert( smP->getMajorDim() == 5 ); 
      assert( smP->getNumElements() == 14 );
      
    }
    // Test adding several cuts
    {
      OsiClpSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      fim.initialSolve();
      OsiRowCut cuts[3];
      
      
      // Generate one ineffective cut plus two trivial cuts
      int c;
      int nc = fim.getNumCols();
      int *inx = new int[nc];
      for (c=0;c<nc;c++) inx[c]=c;
      double *el = new double[nc];
      for (c=0;c<nc;c++) el[c]=1.0e-50+((double)c)*((double)c);
      
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      el[4]=0.0; // to get inf later
      
      for (c=2;c<4;c++) {
        el[0]=1.0;
        inx[0]=c;
        cuts[c-1].setRow(1,inx,el);
        cuts[c-1].setLb(1.);
        cuts[c-1].setUb(100.);
        cuts[c-1].setEffectiveness(c);
      }
      fim.writeMps("x1.mps");
      fim.applyRowCuts(3,cuts);
      fim.writeMps("x2.mps");
      // resolve - should get message about zero elements
      fim.resolve();
      fim.writeMps("x3.mps");
      // check integer solution
      const double * cs = fim.getColSolution();
      CoinRelFltEq eq;
      assert( eq(cs[2],   1.0) );
      assert( eq(cs[3],   1.0) );
      // check will find invalid matrix
      el[0]=1.0/el[4];
      inx[0]=0;
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      fim.applyRowCut(cuts[0]);
      // resolve - should get message about zero elements
      fim.resolve();
      assert (fim.isAbandoned());
      delete[]el;
      delete[]inx;
    }
    // Test using bad basis
    {
      OsiClpSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      fim.initialSolve();
      int numberRows = fim.getModelPtr()->numberRows();
      int * rowStatus = new int[numberRows];
      int numberColumns = fim.getModelPtr()->numberColumns();
      int * columnStatus = new int[numberColumns];
      fim.getBasisStatus(columnStatus,rowStatus);
      int i;
      for (i=0;i<numberRows;i++) {
        rowStatus[i]=2;
      }        
      fim.setBasisStatus(columnStatus,rowStatus);
      fim.resolve();
      delete [] rowStatus;
      delete [] columnStatus;
    }
        // Test matrixByCol method
    {
  
      const OsiClpSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByCol();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   5.6) );
      assert( eq(ev[2],   1.0) );
      assert( eq(ev[3],   2.0) );
      assert( eq(ev[4],   1.1) );
      assert( eq(ev[5],   1.0) );
      assert( eq(ev[6],  -2.0) );
      assert( eq(ev[7],   2.8) );
      assert( eq(ev[8],  -1.0) );
      assert( eq(ev[9],   1.0) );
      assert( eq(ev[10],  1.0) );
      assert( eq(ev[11], -1.2) );
      assert( eq(ev[12], -1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==2 );
      assert( mi[2]==4 );
      assert( mi[3]==6 );
      assert( mi[4]==8 );
      assert( mi[5]==10 );
      assert( mi[6]==11 );
      assert( mi[7]==12 );
      assert( mi[8]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  4 );
      assert( ei[2]  ==  0 );
      assert( ei[3]  ==  1 );
      assert( ei[4]  ==  1 );
      assert( ei[5]  ==  2 );
      assert( ei[6]  ==  0 );
      assert( ei[7]  ==  3 );
      assert( ei[8]  ==  0 );
      assert( ei[9]  ==  4 );
      assert( ei[10] ==  2 );
      assert( ei[11] ==  3 );
      assert( ei[12] ==  0 );
      assert( ei[13] ==  4 );    
      
      assert( smP->getMajorDim() == 8 ); 
      assert( smP->getNumElements() == 14 );

      assert( smP->getSizeVectorStarts()==9 );
      assert( smP->getMinorDim() == 5 );
      
    }
    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow
    {
      OsiClpSolverInterface lhs;
      {      
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rhs_==NULL );
        assert( m.matrixByRow_==NULL );
        
        OsiClpSolverInterface siC1(m);     
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
        assert( siC1.matrixByRow_==NULL );

        const char   * siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        
        const double * siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        
        const double * siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        
        const CoinPackedMatrix * siC1mbr = siC1.getMatrixByRow();
        assert( siC1mbr != NULL );
        
        const double * ev = siC1mbr->getElements();
        assert( eq(ev[0],   3.0) );
        assert( eq(ev[1],   1.0) );
        assert( eq(ev[2],  -2.0) );
        assert( eq(ev[3],  -1.0) );
        assert( eq(ev[4],  -1.0) );
        assert( eq(ev[5],   2.0) );
        assert( eq(ev[6],   1.1) );
        assert( eq(ev[7],   1.0) );
        assert( eq(ev[8],   1.0) );
        assert( eq(ev[9],   2.8) );
        assert( eq(ev[10], -1.2) );
        assert( eq(ev[11],  5.6) );
        assert( eq(ev[12],  1.0) );
        assert( eq(ev[13],  1.9) );
        
        const CoinBigIndex * mi = siC1mbr->getVectorStarts();
        assert( mi[0]==0 );
        assert( mi[1]==5 );
        assert( mi[2]==7 );
        assert( mi[3]==9 );
        assert( mi[4]==11 );
        assert( mi[5]==14 );
        
        const int * ei = siC1mbr->getIndices();
        assert( ei[0]  ==  0 );
        assert( ei[1]  ==  1 );
        assert( ei[2]  ==  3 );
        assert( ei[3]  ==  4 );
        assert( ei[4]  ==  7 );
        assert( ei[5]  ==  1 );
        assert( ei[6]  ==  2 );
        assert( ei[7]  ==  2 );
        assert( ei[8]  ==  5 );
        assert( ei[9]  ==  3 );
        assert( ei[10] ==  6 );
        assert( ei[11] ==  0 );
        assert( ei[12] ==  4 );
        assert( ei[13] ==  7 );    
        
        assert( siC1mbr->getMajorDim() == 5 ); 
        assert( siC1mbr->getNumElements() == 14 );
        

        assert( siC1rs  == siC1.getRowSense() );
        assert( siC1rhs == siC1.getRightHandSide() );
        assert( siC1rr  == siC1.getRowRange() );

        // Change OSL Model by adding free row
        OsiRowCut rc;
        rc.setLb(-DBL_MAX);
        rc.setUb( DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);
             
        // Since model was changed, test that cached
        // data is now freed.
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
        // now updated! assert( siC1.matrixByRow_==NULL );
        
        siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        assert( siC1rs[5]=='N' );

        siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        assert( eq(siC1rhs[5],0.0 ) ); 

        siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        assert( eq(siC1rr[5],0.0) );
    
        lhs=siC1;
      }
      // Test that lhs has correct values even though siC1 has gone out of scope    
      assert( lhs.rowrange_==NULL );
      assert( lhs.rowsense_==NULL );
      assert( lhs.rhs_==NULL ); 
      assert( lhs.matrixByRow_==NULL ); 
      
      const char * lhsrs  = lhs.getRowSense();
      assert( lhsrs[0]=='G' );
      assert( lhsrs[1]=='L' );
      assert( lhsrs[2]=='E' );
      assert( lhsrs[3]=='R' );
      assert( lhsrs[4]=='R' );
      assert( lhsrs[5]=='N' );
      
      const double * lhsrhs = lhs.getRightHandSide();
      assert( eq(lhsrhs[0],2.5) );
      assert( eq(lhsrhs[1],2.1) );
      assert( eq(lhsrhs[2],4.0) );
      assert( eq(lhsrhs[3],5.0) );
      assert( eq(lhsrhs[4],15.) ); 
      assert( eq(lhsrhs[5],0.0) ); 
      
      const double *lhsrr  = lhs.getRowRange();
      assert( eq(lhsrr[0],0.0) );
      assert( eq(lhsrr[1],0.0) );
      assert( eq(lhsrr[2],0.0) );
      assert( eq(lhsrr[3],5.0-1.8) );
      assert( eq(lhsrr[4],15.0-3.0) );
      assert( eq(lhsrr[5],0.0) );      
      
      const CoinPackedMatrix * lhsmbr = lhs.getMatrixByRow();
      assert( lhsmbr != NULL );       
      const double * ev = lhsmbr->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = lhsmbr->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = lhsmbr->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    
      
      int md = lhsmbr->getMajorDim();
      assert(  md == 6 ); 
      assert( lhsmbr->getNumElements() == 14 );
    }
    
  }

  // Test add/delete columns
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    double inf = m.getInfinity();

    CoinPackedVector c0;
    c0.insert(0, 4);
    c0.insert(1, 1);
    m.addCol(c0, 0, inf, 3);
    m.initialSolve();
    double objValue = m.getObjValue();
    CoinRelFltEq eq(1.0e-2);
    assert( eq(objValue,2520.57) );
    // Try deleting first column
    int * d = new int[1];
    d[0]=0;
    m.deleteCols(1,d);
    delete [] d;
    d=NULL;
    m.resolve();
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );
    // Try deleting column we added
    int iCol = m.getNumCols()-1;
    m.deleteCols(1,&iCol);
    m.resolve();
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );

  }
  // Test branch and bound
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    // reduce printout
    m.setHintParam(OsiDoReducePrint,true,OsiHintTry);
    // test maximization
    int n  = m.getNumCols();
    int i;
    double * obj2 = new double[n];
    const double * obj = m.getObjCoefficients();
    for (i=0;i<n;i++) {
      obj2[i]=-obj[i];
    }
    m.setObjective(obj2);
    delete [] obj2;
    m.setObjSense(-1.0);
    // Save bounds
    double * saveUpper = CoinCopyOfArray(m.getColUpper(),n);
    double * saveLower = CoinCopyOfArray(m.getColLower(),n);
    m.branchAndBound();
    // reset bounds
    m.setColUpper(saveUpper);
    m.setColLower(saveLower);
    delete [] saveUpper;
    delete [] saveLower;
    // reset cutoff - otherwise we won't get solution
    m.setDblParam(OsiDualObjectiveLimit,-COIN_DBL_MAX);
    m.branchAndBound();
    double objValue = m.getObjValue();
    CoinRelFltEq eq(1.0e-2);
    assert( eq(objValue,-3089) );
    const double * cs = m.getColSolution();
    for ( i=0;i<n;i++) {
      if (cs[i]>1.0e-7)
        printf("%d has value %g\n",i,cs[i]);
    }
  }
  // Test matt
  if (fopen("../Mps/Infeas/galenet.mps","r")) {    
    OsiClpSolverInterface m;
    m.readMps("../Mps/Infeas/galenet","mps");
    m.setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);
    m.resolve();
    
    std::vector<double *> rays = m.getDualRays(1);
    std::cout << "Dual Ray: " << std::endl;
    for(int i = 0; i < m.getNumRows(); i++){
      if(fabs(rays[0][i]) > 0.00001)
        std::cout << i << " : " << rays[0][i] << std::endl;
    }
    
    std::cout << "isProvenOptimal = " << m.isProvenOptimal() << std::endl;
    std::cout << "isProvenPrimalInfeasible = " << m.isProvenPrimalInfeasible()
         << std::endl;
    
    delete [] rays[0];
    
  }
  // Test infeasible bounds
  {
    OsiClpSolverInterface solver;
    std::string fn = mpsDir+"exmip1";
    solver.readMps(fn.c_str(),"mps");
    int index[]={0};
    double value[]={1.0,0.0};
    solver.setColSetBounds(index,index+1,value);
    solver.setHintParam(OsiDoPresolveInInitial, false, OsiHintDo);
    solver.initialSolve();
    assert (!solver.isProvenOptimal());
  }

  // Build a model
  {    
    OsiClpSolverInterface model;
    std::string fn = mpsDir+"p0033";
    model.readMps(fn.c_str(),"mps");
    // Point to data
    int numberRows = model.getNumRows();
    const double * rowLower = model.getRowLower();
    const double * rowUpper = model.getRowUpper();
    int numberColumns = model.getNumCols();
    const double * columnLower = model.getColLower();
    const double * columnUpper = model.getColUpper();
    const double * columnObjective = model.getObjCoefficients();
    // get row copy
    CoinPackedMatrix rowCopy = *model.getMatrixByRow();
    const int * column = rowCopy.getIndices();
    const int * rowLength = rowCopy.getVectorLengths();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    const double * element = rowCopy.getElements();
    
    // solve
    model.initialSolve();
    // Now build new model
    CoinModel build;
    // Row bounds
    int iRow;
    for (iRow=0;iRow<numberRows;iRow++) {
      build.setRowBounds(iRow,rowLower[iRow],rowUpper[iRow]);
      build.setRowName(iRow,model.getRowName(iRow).c_str());
    }
    // Column bounds and objective
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      build.setColumnLower(iColumn,columnLower[iColumn]);
      build.setColumnUpper(iColumn,columnUpper[iColumn]);
      build.setObjective(iColumn,columnObjective[iColumn]);
      build.setColumnName(iColumn,model.getColName(iColumn).c_str());
    }
    // Adds elements one by one by row (backwards by row)
    for (iRow=numberRows-1;iRow>=0;iRow--) {
      int start = rowStart[iRow];
      for (int j=start;j<start+rowLength[iRow];j++) 
        build(iRow,column[j],element[j]);
    }
    // Now create Model
    OsiClpSolverInterface model2;
    model2.loadFromCoinModel(build);
    model2.initialSolve();
    // Save - should be continuous
    model2.writeMps("continuous");
    int * whichInteger = new int[numberColumns];
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      whichInteger[iColumn]=iColumn;
    // mark as integer
    model2.setInteger(whichInteger,numberColumns);
    delete [] whichInteger;
    // save - should be integer
    model2.writeMps("integer");
    // check names are there
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      printf("%d name %s\n",iColumn,model2.getColName(iColumn).c_str());
    }
    
    // Now do with strings attached
    // Save build to show how to go over rows
    CoinModel saveBuild = build;
    build = CoinModel();
    // Column bounds
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      build.setColumnLower(iColumn,columnLower[iColumn]);
      build.setColumnUpper(iColumn,columnUpper[iColumn]);
    }
    // Objective - half the columns as is and half with multiplier of "1.0+multiplier"
    // Pick up from saveBuild (for no reason at all)
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value = saveBuild.objective(iColumn);
      if (iColumn*2<numberColumns) {
        build.setObjective(iColumn,columnObjective[iColumn]);
      } else {
        // create as string
        char temp[100];
        sprintf(temp,"%g + abs(%g*multiplier)",value,value);
        build.setObjective(iColumn,temp);
      }
    }
    // It then adds rows one by one but for half the rows sets their values
    //      with multiplier of "1.0+1.5*multiplier"
    for (iRow=0;iRow<numberRows;iRow++) {
      if (iRow*2<numberRows) {
        // add row in simple way
        int start = rowStart[iRow];
        build.addRow(rowLength[iRow],column+start,element+start,
                     rowLower[iRow],rowUpper[iRow]);
      } else {
        // As we have to add one by one let's get from saveBuild
        CoinModelLink triple=saveBuild.firstInRow(iRow);
        while (triple.column()>=0) {
          int iColumn = triple.column();
          if (iColumn*2<numberColumns) {
            // just value as normal
            build(iRow,triple.column(),triple.value());
          } else {
            // create as string
            char temp[100];
            sprintf(temp,"%g + (1.5*%g*multiplier)",triple.value(), triple.value());
            build(iRow,iColumn,temp);
          }
          triple=saveBuild.next(triple);
        }
        // but remember to do rhs
        build.setRowLower(iRow,rowLower[iRow]);
        build.setRowUpper(iRow,rowUpper[iRow]);
      }
    }
    // If small switch on error printing
    if (numberColumns<50)
      build.setLogLevel(1);
    int numberErrors=model2.loadFromCoinModel(build);
    // should fail as we never set multiplier
    assert (numberErrors);
    build.associateElement("multiplier",0.0);
    numberErrors=model2.loadFromCoinModel(build);
    assert (!numberErrors);
    model2.initialSolve();
    // It then loops with multiplier going from 0.0 to 2.0 in increments of 0.1
    for (double multiplier=0.0;multiplier<2.0;multiplier+= 0.1) {
      build.associateElement("multiplier",multiplier);
      numberErrors=model2.loadFromCoinModel(build,true);
      assert (!numberErrors);
      model2.resolve();
    }
  }
  // Solve an lp by hand
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    m.setObjSense(-1.0);
    m.getModelPtr()->messageHandler()->setLogLevel(4);
    m.initialSolve();
    m.getModelPtr()->factorization()->maximumPivots(5);
    m.setObjSense(1.0);
    // clone to test pivot as well as primalPivot
    OsiSolverInterface * mm =m.clone();
    // enable special mode
    m.enableSimplexInterface(true);
    // need to do after clone 
    mm->enableSimplexInterface(true);
    // we happen to know that variables are 0-1 and rows are L
    int numberIterations=0;
    int numberColumns = m.getNumCols();
    int numberRows = m.getNumRows();
    double * fakeCost = new double[numberColumns];
    double * duals = new double [numberRows];
    double * djs = new double [numberColumns];
    const double * solution = m.getColSolution();
    memcpy(fakeCost,m.getObjCoefficients(),numberColumns*sizeof(double));
    while (1) {
      const double * dj;
      const double * dual;
      if ((numberIterations&1)==0) {
        // use given ones
        dj = m.getReducedCost();
        dual = m.getRowPrice();
      } else {
        // create
        dj = djs;
        dual = duals;
        m.getReducedGradient(djs,duals,fakeCost);
      }
      int i;
      int colIn=9999;
      int direction=1;
      double best=1.0e-6;
      // find most negative reduced cost
      // Should check basic - but should be okay on this problem
      for (i=0;i<numberRows;i++) {
        double value=dual[i];
        if (value>best) {
          direction=-1;
          best=value;
          colIn=-i-1;
        }
      }
      for (i=0;i<numberColumns;i++) {
        double value=dj[i];
        if (value<-best&&solution[i]<1.0e-6) {
          direction=1;
          best=-value;
          colIn=i;
        } else if (value>best&&solution[i]>1.0-1.0e-6) {
          direction=-1;
          best=value;
          colIn=i;
        }
      }
      if (colIn==9999)
        break; // should be optimal
      int colOut;
      int outStatus;
      double theta;
      int returnCode=m.primalPivotResult(colIn,direction,
					 colOut,outStatus,theta,NULL);
      assert (!returnCode);
      printf("out %d, direction %d theta %g\n",
             colOut,outStatus,theta);
      if (colIn!=colOut) {
        returnCode=mm->pivot(colIn,colOut,outStatus);
        assert (returnCode>=0);
      } else {
        // bound flip (so pivot does not make sense)
        returnCode=mm->primalPivotResult(colIn,direction,
					 colOut,outStatus,theta,NULL);
        assert (!returnCode);
      }
      numberIterations++;
    }
    delete mm;
    delete [] fakeCost;
    delete [] duals;
    delete [] djs;
    // exit special mode
    m.disableSimplexInterface();
    m.getModelPtr()->messageHandler()->setLogLevel(4);
    m.messageHandler()->setLogLevel(0);
    m.resolve();
    assert (!m.getIterationCount());
    m.setObjSense(-1.0);
    m.initialSolve();
  }
# if 0
/*
  This section stops working without setObjectiveAndRefresh. Assertion failure
  down in the guts of clp, likely due to reduced costs not properly updated.
  Leave the code in for a bit so it's easily recoverable if anyone actually
  yells about the loss. There was no response to a public announcement
  of intent to delete, but sometimes it takes a whack on the head to get
  peoples' attention. At some point, it'd be good to come back through and
  make this work again. -- lh, 100828 --
*/
  // Do parametrics on the objective by hand
  {
    std::cout << " Beginning Osi Simplex mode 2 ... " << std::endl ;
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    ClpSimplex * simplex = m.getModelPtr();
    simplex->messageHandler()->setLogLevel(4);
    m.initialSolve();
    simplex->factorization()->maximumPivots(5);
    simplex->messageHandler()->setLogLevel(63);
    m.setObjSense(1.0);
    // enable special mode
    m.enableSimplexInterface(true);
    int numberIterations=0;
    int numberColumns = m.getNumCols();
    int numberRows = m.getNumRows();
    double * changeCost = new double[numberColumns];
    double * duals = new double [numberRows];
    double * djs = new double [numberColumns];
    // As getReducedGradient mucks about with innards of Clp get arrays to save
    double * dualsNow = new double [numberRows];
    double * djsNow = new double [numberColumns];
   
    const double * solution = m.getColSolution();
    int i;
    // Set up change cost
    for (i=0;i<numberColumns;i++)
      changeCost[i]=1.0+0.1*i;;
    // Range of investigation
    double totalChange=100.0;
    double totalDone=0.0;
    while (true) {
      std::cout << " Starting iterations ... " << std::endl ;
      // Save current
      // (would be more accurate to start from scratch)
      memcpy(djsNow, m.getReducedCost(),numberColumns*sizeof(double));
      memcpy(dualsNow,m.getRowPrice(),numberRows*sizeof(double));
      // Get reduced gradient of changeCost
      m.getReducedGradient(djs,duals,changeCost);
      int colIn=9999;
      int direction=1;
      // We are going up to totalChange but we have done some
      double best=totalChange-totalDone;
      // find best ratio
      // Should check basic - but should be okay on this problem
      // We are cheating here as we know L rows
      // Really should be using L and U status but I modified from previous example
      for (i=0;i<numberRows;i++) {
        if (simplex->getRowStatus(i)==ClpSimplex::basic) {
          assert (fabs(dualsNow[i])<1.0e-4&&fabs(duals[i])<1.0e-4);
        } else {
          assert (dualsNow[i]<1.0e-4);
          if (duals[i]>1.0e-8) {
            if (dualsNow[i]+best*duals[i]>0.0) {
              best = CoinMax(-dualsNow[i]/duals[i],0.0);
              direction=-1;
              colIn=-i-1;
            }
          }
        }
      }
      for (i=0;i<numberColumns;i++) {
        if (simplex->getColumnStatus(i)==ClpSimplex::basic) {
          assert (fabs(djsNow[i])<1.0e-4&&fabs(djs[i])<1.0e-4);
        } else {
          if (solution[i]<1.0e-6) {
            assert (djsNow[i]>-1.0e-4);
            if (djs[i]<-1.0e-8) {
              if (djsNow[i]+best*djs[i]<0.0) {
                best = CoinMax(-djsNow[i]/djs[i],0.0);
                direction=1;
                colIn=i;
              }
            }
          } else if (solution[i]>1.0-1.0e-6) {
            assert (djsNow[i]<1.0e-4);
            if (djs[i]>1.0e-8) {
              if (djsNow[i]+best*djs[i]>0.0) {
                best = CoinMax(-djsNow[i]/djs[i],0.0);
                direction=-1;
                colIn=i;
              }
            }
          }
        }
      }
      if (colIn==9999)
        break; // should be optimal
      // update objective - djs is spare array
      const double * obj = m.getObjCoefficients();
      for (i=0;i<numberColumns;i++) {
        djs[i]=obj[i]+best*changeCost[i];
      }
      totalDone += best;
      printf("Best change %g, total %g\n",best,totalDone);
      m.setObjectiveAndRefresh(djs);
      int colOut;
      int outStatus;
      double theta;
      int returnCode=m.primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL);
      assert (!returnCode);
      double objValue = m.getObjValue(); // printed one may not be accurate
      printf("in %d out %d, direction %d theta %g, objvalue %g\n",
             colIn,colOut,outStatus,theta,objValue);
      numberIterations++;
    }
    // update objective to totalChange- djs is spare array
    const double * obj = m.getObjCoefficients();
    double best = totalChange-totalDone;
    for (i=0;i<numberColumns;i++) {
      djs[i]=obj[i]+best*changeCost[i];
    }
    delete [] changeCost;
    delete [] duals;
    delete [] djs;
    delete [] dualsNow;
    delete [] djsNow;
    // exit special mode
    std::cout << " Finished simplex mode 2 ; checking result." << std::endl ;
    m.disableSimplexInterface();
    simplex->messageHandler()->setLogLevel(4);
    m.resolve();
    assert (!m.getIterationCount());
  }
# endif
  // Solve an lp when interface is on
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    // enable special mode
    m.setHintParam(OsiDoScale,false,OsiHintDo);
    m.setHintParam(OsiDoPresolveInInitial,false,OsiHintDo);
    m.setHintParam(OsiDoDualInInitial,false,OsiHintDo);
    m.setHintParam(OsiDoPresolveInResolve,false,OsiHintDo);
    m.setHintParam(OsiDoDualInResolve,false,OsiHintDo);
    m.enableSimplexInterface(true);
    m.initialSolve();
  }
  // Check tableau stuff when simplex interface is on
  {    
    OsiClpSolverInterface m;
    /* 
       Wolsey : Page 130
       max 4x1 -  x2
       7x1 - 2x2    <= 14
       x2    <= 3
       2x1 - 2x2    <= 3
       x1 in Z+, x2 >= 0
    */
    
    double inf_ = m.getInfinity();
    int n_cols = 2;
    int n_rows = 3;
    
    double obj[2] = {-4.0, 1.0};
    double collb[2] = {0.0, 0.0};
    double colub[2] = {inf_, inf_};
    double rowlb[3] = {-inf_, -inf_, -inf_};
    double rowub[3] = {14.0, 3.0, 3.0};
    
    int rowIndices[5] =  {0,     2,    0,    1,    2};
    int colIndices[5] =  {0,     0,    1,    1,    1};
    double elements[5] = {7.0, 2.0, -2.0,  1.0, -2.0};
    CoinPackedMatrix M(true, rowIndices, colIndices, elements, 5);
    
    m.loadProblem(M, collb, colub, obj, rowlb, rowub);
    m.enableSimplexInterface(true);
    
    m.initialSolve();
    
    //check that the tableau matches wolsey (B-1 A)
    // slacks in second part of binvA
    double * binvA = (double*) malloc((n_cols+n_rows) * sizeof(double));
    
    printf("B-1 A");
    int i;
    for( i = 0; i < n_rows; i++){
      m.getBInvARow(i, binvA,binvA+n_cols);
      printf("\nrow: %d -> ",i);
      for(int j=0; j < n_cols+n_rows; j++){
        printf("%g, ", binvA[j]);
      }
    }
    printf("\n");
    printf("And by column");
    for( i = 0; i < n_cols+n_rows; i++){
      m.getBInvACol(i, binvA);
      printf("\ncolumn: %d -> ",i);
      for(int j=0; j < n_rows; j++){
        printf("%g, ", binvA[j]);
      }
    }
    printf("\n");
    // and when doing as expert
    m.setSpecialOptions(m.specialOptions()|512);
    ClpSimplex * clp = m.getModelPtr();
    /* Do twice -
       first time with enableSimplexInterface still set
       then without and with scaling
    */
    for (int iPass=0;iPass<2;iPass++) {
      const double * rowScale = clp->rowScale();
      const double * columnScale = clp->columnScale();
      if (!iPass)
        assert (!rowScale);
      else
        assert (rowScale); // only true for this example
      /* has to be exactly correct as in OsiClpsolverInterface.cpp
         (also redo each pass as may change
      */
      CoinIndexedVector * rowArray = clp->rowArray(1);
      CoinIndexedVector * columnArray = clp->columnArray(0);
      int n;
      int * which;
      double * array;
      printf("B-1 A");
      for( i = 0; i < n_rows; i++){
        m.getBInvARow(i, binvA,binvA+n_cols);
        printf("\nrow: %d -> ",i);
        int j;
        // First columns
        n = columnArray->getNumElements();
        which = columnArray->getIndices();
        array = columnArray->denseVector();
        for(j=0; j < n; j++){
          int k=which[j];
          if (!columnScale) {
            printf("(%d %g), ", k, array[k]);
          } else {
            printf("(%d %g), ", k, array[k]/columnScale[k]);
          }
          // zero out
          array[k]=0.0;
        }
        // say empty
        columnArray->setNumElements(0);
        // and check (would not be in any production code)
        columnArray->checkClear();
        // now rows
        n = rowArray->getNumElements();
        which = rowArray->getIndices();
        array = rowArray->denseVector();
        for(j=0; j < n; j++){
          int k=which[j];
          if (!rowScale) {
            printf("(%d %g), ", k+n_cols, array[k]);
          } else {
            printf("(%d %g), ", k+n_cols, array[k]*rowScale[k]);
          }
          // zero out
          array[k]=0.0;
        }
        // say empty
        rowArray->setNumElements(0);
        // and check (would not be in any production code)
        rowArray->checkClear();
      }
      printf("\n");
      printf("And by column (trickier)");
      const int * pivotVariable = clp->pivotVariable();
      for( i = 0; i < n_cols+n_rows; i++){
        m.getBInvACol(i, binvA);
        printf("\ncolumn: %d -> ",i);
        n = rowArray->getNumElements();
        which = rowArray->getIndices();
        array = rowArray->denseVector();
        for(int j=0; j < n; j++){
          int k=which[j];
          // need to know pivot variable for +1/-1 (slack) and row/column scaling
          int pivot = pivotVariable[k];
          if (pivot<n_cols) {
            // scaled coding is in just in case 
            if (!columnScale) {
              printf("(%d %g), ", k, array[k]);
            } else {
              printf("(%d %g), ", k, array[k]*columnScale[pivot]);
            }
          } else {
            if (!rowScale) {
              printf("(%d %g), ", k, -array[k]);
            } else {
              printf("(%d %g), ", k, -array[k]/rowScale[pivot-n_cols]);
            }
          }
          // zero out
          array[k]=0.0;
        }
        // say empty
        rowArray->setNumElements(0);
        // and check (would not be in any production code)
        rowArray->checkClear();
      }
      printf("\n");
      // now deal with next pass
      if (!iPass) {
        m.disableSimplexInterface();
        // see if we can get scaling for testing
        clp->scaling(1);
        m.enableFactorization();
      } else {
        // may not be needed - but cleaner
        m.disableFactorization();
      }
    }
    m.setSpecialOptions(m.specialOptions()&~512);
    free(binvA);
    // Do using newer interface
    m.enableFactorization();
    {
      CoinIndexedVector * rowArray = new CoinIndexedVector(n_rows);
      CoinIndexedVector * columnArray = new CoinIndexedVector(n_cols);
      int n;
      int * which;
      double * array;
      printf("B-1 A");
      for( i = 0; i < n_rows; i++){
        m.getBInvARow(i, columnArray,rowArray);
        printf("\nrow: %d -> ",i);
        int j;
        // First columns
        n = columnArray->getNumElements();
        which = columnArray->getIndices();
        array = columnArray->denseVector();
        for(j=0; j < n; j++){
          int k=which[j];
	  printf("(%d %g), ", k, array[k]);
          // zero out
          array[k]=0.0;
        }
        // say empty (if I had not zeroed array[k] I would use ->clear())
        columnArray->setNumElements(0);
        // and check (would not be in any production code)
        columnArray->checkClear();
        // now rows
        n = rowArray->getNumElements();
        which = rowArray->getIndices();
        array = rowArray->denseVector();
        for(j=0; j < n; j++){
          int k=which[j];
	  printf("(%d %g), ", k+n_cols, array[k]);
          // zero out
          array[k]=0.0;
        }
        // say empty
        rowArray->setNumElements(0);
        // and check (would not be in any production code)
        rowArray->checkClear();
      }
      printf("\n");
      printf("And by column");
      const int * pivotVariable = clp->pivotVariable();
      for( i = 0; i < n_cols+n_rows; i++){
        m.getBInvACol(i, rowArray);
        printf("\ncolumn: %d -> ",i);
        n = rowArray->getNumElements();
        which = rowArray->getIndices();
        array = rowArray->denseVector();
        for(int j=0; j < n; j++){
          int k=which[j];
	  printf("(%d %g), ", k, array[k]);
          // zero out
          array[k]=0.0;
        }
        // say empty
        rowArray->setNumElements(0);
        // and check (would not be in any production code)
        rowArray->checkClear();
      }
      printf("\n");
      delete rowArray;
      delete columnArray;
    }
    // may not be needed - but cleaner
    m.disableFactorization();
  }
  // Check tableau stuff when simplex interface is off
  {    
    OsiClpSolverInterface m;
    /* 
       Wolsey : Page 130
       max 4x1 -  x2
       7x1 - 2x2    <= 14
       x2    <= 3
       2x1 - 2x2    <= 3
       x1 in Z+, x2 >= 0
    */
    
    double inf_ = m.getInfinity();
    int n_cols = 2;
    int n_rows = 3;
    
    double obj[2] = {-4.0, 1.0};
    double collb[2] = {0.0, 0.0};
    double colub[2] = {inf_, inf_};
    double rowlb[3] = {-inf_, -inf_, -inf_};
    double rowub[3] = {14.0, 3.0, 3.0};
    
    int rowIndices[5] =  {0,     2,    0,    1,    2};
    int colIndices[5] =  {0,     0,    1,    1,    1};
    double elements[5] = {7.0, 2.0, -2.0,  1.0, -2.0};
    CoinPackedMatrix M(true, rowIndices, colIndices, elements, 5);
    
    m.loadProblem(M, collb, colub, obj, rowlb, rowub);
    // Test with scaling
    m.setHintParam(OsiDoScale,true);
    // Tell code to keep factorization
    m.setupForRepeatedUse(3,0);
    
    m.initialSolve();
    // Repeated stuff only in resolve
    // m.resolve(); (now in both (just for (3,0))
    
    //check that the tableau matches wolsey (B-1 A)
    // slacks in second part of binvA
    double * binvA = (double*) malloc((n_cols+n_rows) * sizeof(double));
    // and getBasics
    int * pivots = new int[n_rows];
    m.getBasics(pivots);
    
    printf("B-1 A");
    int i;
    for( i = 0; i < n_rows; i++){
      m.getBInvARow(i, binvA,binvA+n_cols);
      printf("\nrow: %d (pivot %d) -> ",i,pivots[i]);
      for(int j=0; j < n_cols+n_rows; j++){
        printf("%g, ", binvA[j]);
      }
    }
    printf("\n");
    printf("And by column");
    for( i = 0; i < n_cols+n_rows; i++){
      m.getBInvACol(i, binvA);
      printf("\ncolumn: %d -> ",i);
      for(int j=0; j < n_rows; j++){
        printf("%g, ", binvA[j]);
      }
    }
    printf("\n");
    free(binvA);
    delete [] pivots;
  }
  // Do common solverInterface testing 
  {
    OsiClpSolverInterface m;
    return OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
}
