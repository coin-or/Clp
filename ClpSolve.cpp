// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

// This file has higher level solve functions


#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#include "ClpInterior.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#endif
#include "ClpSolve.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpMessage.hpp"
#include "CoinTime.hpp"

#include "ClpPresolve.hpp"
#ifndef SLIM_CLP
#include "Idiot.hpp"
#ifdef WSSMP_BARRIER
#include "ClpCholeskyWssmp.hpp"
#include "ClpCholeskyWssmpKKT.hpp"
#define FAST_BARRIER
#endif
#ifdef UFL_BARRIER
#include "ClpCholeskyUfl.hpp"
#define FAST_BARRIER
#endif
#ifdef TAUCS_BARRIER
#include "ClpCholeskyTaucs.hpp"
#define FAST_BARRIER
#endif
#ifndef FAST_BARRIER
static int numberBarrier=0;
#endif
#ifdef COIN_HAS_VOL
#include "VolVolume.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinMpsIO.hpp"

//#############################################################################

class lpHook : public VOL_user_hooks {
private:
   lpHook(const lpHook&);
   lpHook& operator=(const lpHook&);
private:
   /// Pointer to dense vector of structural variable upper bounds
   double  *colupper_;
   /// Pointer to dense vector of structural variable lower bounds
   double  *collower_;
   /// Pointer to dense vector of objective coefficients
   double  *objcoeffs_;
   /// Pointer to dense vector of right hand sides
   double  *rhs_;
   /// Pointer to dense vector of senses
   char    *sense_;

   /// The problem matrix in a row ordered form
   CoinPackedMatrix rowMatrix_;
   /// The problem matrix in a column ordered form
   CoinPackedMatrix colMatrix_;

public:
   lpHook(double* clb, double* cub, double* obj,
	  double* rhs, char* sense, const CoinPackedMatrix& mat);
   virtual ~lpHook();
   
public:
   // for all hooks: return value of -1 means that volume should quit
   /** compute reduced costs    
       @param u (IN) the dual variables
       @param rc (OUT) the reduced cost with respect to the dual values
   */
   virtual int compute_rc(const VOL_dvector& u, VOL_dvector& rc);

   /** Solve the subproblem for the subgradient step.
       @param dual (IN) the dual variables
       @param rc (IN) the reduced cost with respect to the dual values
       @param lcost (OUT) the lagrangean cost with respect to the dual values
       @param x (OUT) the primal result of solving the subproblem
       @param v (OUT) b-Ax for the relaxed constraints
       @param pcost (OUT) the primal objective value of <code>x</code>
   */
   virtual int solve_subproblem(const VOL_dvector& dual, const VOL_dvector& rc,
				double& lcost, VOL_dvector& x, VOL_dvector& v,
				double& pcost);
   /** Starting from the primal vector x, run a heuristic to produce
       an integer solution  
       @param x (IN) the primal vector
       @param heur_val (OUT) the value of the integer solution (return 
       <code>DBL_MAX</code> here if no feas sol was found 
   */
   virtual int heuristics(const VOL_problem& p, 
			  const VOL_dvector& x, double& heur_val) {
      return 0;
   }
};
 
//#############################################################################

lpHook::lpHook(double* clb, double* cub, double* obj,
	       double* rhs, char* sense,
	       const CoinPackedMatrix& mat)
{
   colupper_ = cub;
   collower_ = clb;
   objcoeffs_ = obj;
   rhs_ = rhs;
   sense_ = sense;
   assert (mat.isColOrdered());
   colMatrix_.copyOf(mat);
   rowMatrix_.reverseOrderedCopyOf(mat);
}

//-----------------------------------------------------------------------------

lpHook::~lpHook()
{
}

//#############################################################################

int
lpHook::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
   rowMatrix_.transposeTimes(u.v, rc.v);
   const int psize = rowMatrix_.getNumCols();

   for (int i = 0; i < psize; ++i)
      rc[i] = objcoeffs_[i] - rc[i];
   return 0;
}

//-----------------------------------------------------------------------------

int
lpHook::solve_subproblem(const VOL_dvector& dual, const VOL_dvector& rc,
			 double& lcost, VOL_dvector& x, VOL_dvector& v,
			 double& pcost)
{
   int i;
   const int psize = x.size();
   const int dsize = v.size();

   // compute the lagrangean solution corresponding to the reduced costs
   for (i = 0; i < psize; ++i) 
      x[i] = (rc[i] >= 0.0) ? collower_[i] : colupper_[i];

   // compute the lagrangean value (rhs*dual + primal*rc)
   lcost = 0;
   for (i = 0; i < dsize; ++i)
      lcost += rhs_[i] * dual[i];
   for (i = 0; i < psize; ++i)
      lcost += x[i] * rc[i];

   // compute the rhs - lhs 
   colMatrix_.times(x.v, v.v);
   for (i = 0; i < dsize; ++i)
      v[i] = rhs_[i] - v[i];

   // compute the lagrangean primal objective
   pcost = 0;
   for (i = 0; i < psize; ++i)
      pcost += x[i] * objcoeffs_[i];

   return 0;
}

//#############################################################################
/** A quick inlined function to convert from lb/ub style constraint
    definition to sense/rhs/range style */
inline void
convertBoundToSense(const double lower, const double upper,
					char& sense, double& right,
					double& range) 
{
  range = 0.0;
  if (lower > -1.0e20) {
    if (upper < 1.0e20) {
      right = upper;
      if (upper==lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < 1.0e20) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

static int
solveWithVolume(ClpSimplex * model, int numberPasses, int doIdiot)
{
   VOL_problem volprob;
   volprob.parm.gap_rel_precision=0.00001;
   volprob.parm.maxsgriters=3000;
   if(numberPasses>3000) {
     volprob.parm.maxsgriters=numberPasses;
     volprob.parm.primal_abs_precision=0.0;
     volprob.parm.minimum_rel_ascent=0.00001;
   } else if (doIdiot>0) {
     volprob.parm.maxsgriters=doIdiot;
   }
   if (model->logLevel()<2) 
     volprob.parm.printflag=0;
   else
     volprob.parm.printflag=3;
   const CoinPackedMatrix* mat = model->matrix();
   int psize = model->numberColumns();
   int dsize = model->numberRows();
   char * sense = new char[dsize];
   double * rhs = new double [dsize];

   // Set the lb/ub on the duals
   volprob.dsize = dsize;
   volprob.psize = psize;
   volprob.dual_lb.allocate(dsize);
   volprob.dual_ub.allocate(dsize);
   int i;
   const double * rowLower = model->rowLower();
   const double * rowUpper = model->rowUpper();
   for (i = 0; i < dsize; ++i) {
     double range;
     convertBoundToSense(rowLower[i],rowUpper[i],
                         sense[i],rhs[i],range);
      switch (sense[i]) {
       case 'E':
	 volprob.dual_lb[i] = -1.0e31;
	 volprob.dual_ub[i] = 1.0e31;
	 break;
       case 'L':
	 volprob.dual_lb[i] = -1.0e31;
	 volprob.dual_ub[i] = 0.0;
	 break;
       case 'G':
	 volprob.dual_lb[i] = 0.0;
	 volprob.dual_ub[i] = 1.0e31;
	 break;
       default:
	 printf("Volume Algorithm can't work if there is a non ELG row\n");
	 return 1;
      }
   }
   // Check bounds
   double * saveLower = model->columnLower();
   double * saveUpper = model->columnUpper();
   bool good=true;
   for (i=0;i<psize;i++) {
     if (saveLower[i]<-1.0e20||saveUpper[i]>1.0e20) {
       good=false;
       break;
     }
   }
   if (!good) {
     saveLower = CoinCopyOfArray(model->columnLower(),psize);
     saveUpper = CoinCopyOfArray(model->columnUpper(),psize);
     for (i=0;i<psize;i++) {
       if (saveLower[i]<-1.0e20)
         saveLower[i]=-1.0e20;
       if(saveUpper[i]>1.0e20) 
         saveUpper[i]=1.0e20;
     }
   }
   lpHook myHook(saveLower, saveUpper,
		 model->objective(),
		 rhs, sense, *mat);

   volprob.solve(myHook, false /* no warmstart */);

   if (saveLower!=model->columnLower()) {
     delete [] saveLower;
     delete [] saveUpper;
   }
   //------------- extract the solution ---------------------------

   //printf("Best lagrangean value: %f\n", volprob.value);

   double avg = 0;
   for (i = 0; i < dsize; ++i) {
      switch (sense[i]) {
       case 'E':
	 avg += CoinAbs(volprob.viol[i]);
	 break;
       case 'L':
	 if (volprob.viol[i] < 0)
	    avg +=  (-volprob.viol[i]);
	 break;
       case 'G':
	 if (volprob.viol[i] > 0)
	    avg +=  volprob.viol[i];
	 break;
      }
   }
      
   //printf("Average primal constraint violation: %f\n", avg/dsize);

   // volprob.dsol contains the dual solution (dual feasible)
   // volprob.psol contains the primal solution
   //              (NOT necessarily primal feasible)
   CoinMemcpyN(volprob.dsol.v,dsize,model->dualRowSolution());
   CoinMemcpyN(volprob.psol.v,psize,model->primalColumnSolution());
   return 0;
}
#endif
static ClpInterior * currentModel2 = NULL;
#endif
//#############################################################################
// Allow for interrupts
// But is this threadsafe ? (so switched off by option)

#include "CoinSignal.hpp"
static ClpSimplex * currentModel = NULL;

extern "C" {
   static void signal_handler(int whichSignal)
   {
      if (currentModel!=NULL) 
	 currentModel->setMaximumIterations(0); // stop at next iterations
#ifndef SLIM_CLP
      if (currentModel2!=NULL) 
	 currentModel2->setMaximumBarrierIterations(0); // stop at next iterations
#endif
      return;
   }
}

/** General solve algorithm which can do presolve
    special options (bits)
    1 - do not perturb
    2 - do not scale
    4 - use crash (default allslack in dual, idiot in primal)
    8 - all slack basis in primal
    16 - switch off interrupt handling
    32 - do not try and make plus minus one matrix
    64 - do not use sprint even if problem looks good
 */
int 
ClpSimplex::initialSolve(ClpSolve & options)
{
  ClpSolve::SolveType method=options.getSolveType();
  //ClpSolve::SolveType originalMethod=method;
  ClpSolve::PresolveType presolve = options.getPresolveType();
  int saveMaxIterations = maximumIterations();
  int finalStatus=-1;
  int numberIterations=0;
  double time1 = CoinCpuTime();
  double timeX = time1;
  double time2;
  ClpMatrixBase * saveMatrix=NULL;
  ClpSimplex * model2 = this;
  bool interrupt = (options.getSpecialOption(2)==0);
  CoinSighandler_t saveSignal=SIG_DFL;
  if (interrupt) {
    currentModel = model2;
    // register signal handler
    saveSignal = signal(SIGINT,signal_handler);
  }
  // If no status array - set up basis
  if (!status_)
    allSlackBasis();
  ClpPresolve pinfo;
  pinfo.setSubstitution(options.substitution());
  int presolveOptions = options.presolveActions();
  bool presolveToFile = (presolveOptions&0x40000000)!=0;
  presolveOptions &= ~0x40000000;
  if ((presolveOptions&0xffff)!=0)
    pinfo.setPresolveActions(presolveOptions);
  // switch off singletons to slacks
  //pinfo.setDoSingletonColumn(false); // done by bits
  int printOptions = options.getSpecialOption(5);
  if ((printOptions&1)!=0)
    pinfo.statistics();
  double timePresolve=0.0;
  double timeIdiot=0.0;
  double timeCore=0.0;
  int savePerturbation=perturbation_;
  int saveScaling = scalingFlag_;
#ifndef SLIM_CLP
#ifndef NO_RTTI
  if (dynamic_cast< ClpNetworkMatrix*>(matrix_)) {
    // network - switch off stuff
    presolve = ClpSolve::presolveOff;
  }
#else
  if (matrix_->type()==11) {
    // network - switch off stuff
    presolve = ClpSolve::presolveOff;
  }
#endif
#endif
  // For below >0 overrides
  // 0 means no, -1 means maybe
  int doIdiot=0;
  int doCrash=0;
  int doSprint=0;
  int doSlp=0;
  int primalStartup=1;
  if (method!=ClpSolve::useDual&&method!=ClpSolve::useBarrier
      &&method!=ClpSolve::useBarrierNoCross) {
    switch (options.getSpecialOption(1)) {
    case 0:
      doIdiot=-1;
      doCrash=-1;
      doSprint=-1;
      break;
    case 1:
      doIdiot=0;
      doCrash=1;
      if (options.getExtraInfo(1)>0)
	doCrash = options.getExtraInfo(1);
      doSprint=0;
      break;
    case 2:
      doIdiot=1;
      if (options.getExtraInfo(1)>0)
	doIdiot = options.getExtraInfo(1);
      doCrash=0;
      doSprint=0;
      break;
    case 3:
      doIdiot=0;
      doCrash=0;
      doSprint=1;
      break;
    case 4:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      break;
    case 5:
      doIdiot=0;
      doCrash=-1;
      doSprint=-1;
      break;
    case 6:
      doIdiot=-1;
      doCrash=-1;
      doSprint=0;
      break;
    case 7:
      doIdiot=-1;
      doCrash=0;
      doSprint=-1;
      break;
    case 8:
      doIdiot=-1;
      doCrash=0;
      doSprint=0;
      break;
    case 9:
      doIdiot=0;
      doCrash=0;
      doSprint=-1;
      break;
    case 10:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      if (options.getExtraInfo(1)>0)
	doSlp = options.getExtraInfo(1);
      break;
    case 11:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      primalStartup=0;
      break;
    default:
      abort();
    }
  } else {
    // Dual
    switch (options.getSpecialOption(0)) {
    case 0:
      doIdiot=0;
      doCrash=0;
      doSprint=0;
      break;
    case 1:
      doIdiot=0;
      doCrash=1;
      if (options.getExtraInfo(0)>0)
	doCrash = options.getExtraInfo(0);
      doSprint=0;
      break;
    case 2:
      doIdiot=-1;
      if (options.getExtraInfo(0)>0)
	doIdiot = options.getExtraInfo(0);
      doCrash=0;
      doSprint=0;
      break;
    default:
      abort();
    }
  }
#ifdef COIN_HAS_VOL
  // Save number of idiot
  int saveDoIdiot=doIdiot;
#endif
  // Just do this number of passes in Sprint
  int maxSprintPass=100;
  bool costedSlacks=false;
  if (presolve!=ClpSolve::presolveOff) {
    int numberPasses=5;
    if (presolve==ClpSolve::presolveNumber) {
      numberPasses=options.getPresolvePasses();
      presolve = ClpSolve::presolveOn;
    } else if (presolve==ClpSolve::presolveNumberCost) {
      numberPasses=options.getPresolvePasses();
      presolve = ClpSolve::presolveOn;
      costedSlacks=true;
      // switch on singletons to slacks
      pinfo.setDoSingletonColumn(true);
    }
#ifndef CLP_NO_STD
    if (presolveToFile) {
      // PreSolve to file - not fully tested
      printf("Presolving to file - presolve.save\n");
      pinfo.presolvedModelToFile(*this,"presolve.save",1.0e-8,
			   false,numberPasses);
      model2=this;
    } else {
#endif
      model2 = pinfo.presolvedModel(*this,1.0e-8,
				    false,numberPasses,true,costedSlacks);
#ifndef CLP_NO_STD
    }
#endif
    time2 = CoinCpuTime();
    timePresolve = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Presolve"<<timePresolve<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    if (!model2) {
      handler_->message(CLP_INFEASIBLE,messages_)
	<<CoinMessageEol;
      model2 = this;
      if (options.infeasibleReturn()) {
        return -1;
      }
      presolve=ClpSolve::presolveOff;
    }
    // We may be better off using original (but if dual leave because of bounds)
    if (presolve!=ClpSolve::presolveOff&&
	numberRows_<1.01*model2->numberRows_&&numberColumns_<1.01*model2->numberColumns_
	&&model2!=this) {
      if(method!=ClpSolve::useDual||
	 (numberRows_==model2->numberRows_&&numberColumns_==model2->numberColumns_)) {
	delete model2;
	model2 = this;
	presolve=ClpSolve::presolveOff;
      }
    }
  }
  if (interrupt)
    currentModel = model2;
  // See if worth trying +- one matrix
  bool plusMinus=false;
  int numberElements=model2->getNumElements();
#ifndef SLIM_CLP
#ifndef NO_RTTI
  if (dynamic_cast< ClpNetworkMatrix*>(matrix_)) {
    // network - switch off stuff
    doIdiot=0;
    doSprint=0;
  }
#else
  if (matrix_->type()==11) {
    // network - switch off stuff
    doIdiot=0;
    doSprint=0;
  }
#endif
#endif
  int numberColumns = model2->numberColumns();
  int numberRows = model2->numberRows();
  // If not all slack basis - switch off all
  int number=0;
  int iRow;
  for (iRow=0;iRow<numberRows;iRow++)
    if (model2->getRowStatus(iRow)==basic)
      number++;
  if (number<numberRows) {
    doIdiot=0;
    doCrash=0;
    doSprint=0;
  }
  if (options.getSpecialOption(3)==0) {
    if(numberElements>100000)
      plusMinus=true;
    if(numberElements>10000&&(doIdiot||doSprint)) 
      plusMinus=true;
  }
#ifndef SLIM_CLP
  // Statistics (+1,-1, other) - used to decide on strategy if not +-1
  CoinBigIndex statistics[3]={-1,0,0};
  if (plusMinus) {
    saveMatrix = model2->clpMatrix();
#ifndef NO_RTTI
    ClpPackedMatrix* clpMatrix =
      dynamic_cast< ClpPackedMatrix*>(saveMatrix);
#else
    ClpPackedMatrix* clpMatrix = NULL;
    if (saveMatrix->type()==1)
      clpMatrix =
	static_cast< ClpPackedMatrix*>(saveMatrix);
#endif
    if (clpMatrix) {
      ClpPlusMinusOneMatrix * newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
      if (newMatrix->getIndices()) {
	model2->replaceMatrix(newMatrix);
      } else {
	handler_->message(CLP_MATRIX_CHANGE,messages_)
	  <<"+- 1"
	  <<CoinMessageEol;
        CoinMemcpyN(newMatrix->startPositive(),3,statistics);
	saveMatrix=NULL;
	plusMinus=false;
	delete newMatrix;
      }
    } else {
      saveMatrix=NULL;
      plusMinus=false;
    }
  }
#endif
  if (this->factorizationFrequency()==200) {
    // User did not touch preset
    int numberRows = model2->numberRows();
    const int cutoff1=10000;
    const int cutoff2=100000;
    const int base=75;
    const int freq0 = 50;
    const int freq1=200;
    const int freq2=400;
    const int maximum=2000;
    int frequency;
    if (numberRows<cutoff1)
      frequency=base+numberRows/freq0;
    else if (numberRows<cutoff2)
      frequency=base+cutoff1/freq0 + (numberRows-cutoff1)/freq1;
    else
      frequency=base+cutoff1/freq0 + (cutoff2-cutoff1)/freq1 + (numberRows-cutoff2)/freq2;
    model2->setFactorizationFrequency(CoinMin(maximum,frequency));
  } else if (model2!=this) {
    // make sure model2 has correct value
    model2->setFactorizationFrequency(this->factorizationFrequency());
  }
  if (method==ClpSolve::automatic) {
    if (doSprint==0&&doIdiot==0) {
      // off
      method=ClpSolve::useDual;
    } else {
      // only do primal if sprint or idiot
      if (doSprint>0) {
        method=ClpSolve::usePrimalorSprint;
      } else if (doIdiot>0) {
        method=ClpSolve::usePrimal;
      } else {
        if (numberElements<500000) {
          // Small problem
          if(numberRows*10>numberColumns||numberColumns<6000
             ||(numberRows*20>numberColumns&&!plusMinus))
            doSprint=0; // switch off sprint
        } else {
          // larger problem
          if(numberRows*8>numberColumns)
            doSprint=0; // switch off sprint
        }
        // switch off sprint or idiot if any free variable
	int iColumn;
	double * columnLower = model2->columnLower();
	double * columnUpper = model2->columnUpper();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
          if (columnLower[iColumn]<-1.0e10&&columnUpper[iColumn]>1.0e10) {
            doSprint=0;
            doIdiot=0;
            break;
          }
        }
        int nPasses=0;
        // look at rhs
        int iRow;
        double largest=0.0;
        double smallest = 1.0e30;
        double largestGap=0.0;
        int numberNotE=0;
        bool notInteger=false;
        for (iRow=0;iRow<numberRows;iRow++) {
          double value1 = model2->rowLower_[iRow];
          if (value1&&value1>-1.0e31) {
            largest = CoinMax(largest,fabs(value1));
            smallest=CoinMin(smallest,fabs(value1));
            if (fabs(value1-floor(value1+0.5))>1.0e-8) {
              notInteger=true;
              break;
            }
          }
          double value2 = model2->rowUpper_[iRow];
          if (value2&&value2<1.0e31) {
            largest = CoinMax(largest,fabs(value2));
            smallest=CoinMin(smallest,fabs(value2));
            if (fabs(value2-floor(value2+0.5))>1.0e-8) {
              notInteger=true;
              break;
            }
          }
          if (value2>value1) {
            numberNotE++;
            if (value2>1.0e31||value1<-1.0e31)
              largestGap = COIN_DBL_MAX;
            else
              largestGap = value2-value1;
          }
        }
        bool tryIt= numberRows>200&&numberColumns>2000&&numberColumns>2*numberRows;
        if (numberRows<1000&&numberColumns<3000)
          tryIt=false;
        if (notInteger)
          tryIt=false;
        if (largest/smallest>10||(largest/smallest>2.0&&largest>50))
          tryIt=false;
        if (tryIt) {
          if (largest/smallest>2.0) {
            nPasses = 10+numberColumns/100000;
            nPasses = CoinMin(nPasses,50);
            nPasses = CoinMax(nPasses,15);
            if (numberRows>20000&&nPasses>5) {
              // Might as well go for it
              nPasses = CoinMax(nPasses,71);
            } else if (numberRows>2000&&nPasses>5) {
              nPasses = CoinMax(nPasses,50);
            } else if (numberElements<3*numberColumns) {
              nPasses=CoinMin(nPasses,10); // probably not worh it
            }
          } else if (largest/smallest>1.01||numberElements<=3*numberColumns) {
            nPasses = 10+numberColumns/1000;
            nPasses = CoinMin(nPasses,100);
            nPasses = CoinMax(nPasses,30);
            if (numberRows>25000) {
              // Might as well go for it
              nPasses = CoinMax(nPasses,71);
            }
            if (!largestGap)
              nPasses *= 2;
          } else {
            nPasses = 10+numberColumns/1000;
            nPasses = CoinMin(nPasses,200);
            nPasses = CoinMax(nPasses,100);
            if (!largestGap)
              nPasses *= 2;
          }
        }
        //printf("%d rows %d cols plus %c tryIt %c largest %g smallest %g largestGap %g npasses %d sprint %c\n",
        //     numberRows,numberColumns,plusMinus ? 'Y' : 'N',
        //     tryIt ? 'Y' :'N',largest,smallest,largestGap,nPasses,doSprint ? 'Y' :'N');
        //exit(0);
        if (!tryIt||nPasses<=5)
          doIdiot=0;
        if (doSprint) {
          method = ClpSolve::usePrimalorSprint;
        } else if (doIdiot) {
          method = ClpSolve::usePrimal;
        } else {
          method = ClpSolve::useDual;
        }
      }
    }
  }
  if (method==ClpSolve::usePrimalorSprint) {
    if (doSprint<0) { 
      if (numberElements<500000) {
	// Small problem
	if(numberRows*10>numberColumns||numberColumns<6000
	   ||(numberRows*20>numberColumns&&!plusMinus))
	  method=ClpSolve::usePrimal; // switch off sprint
      } else {
	// larger problem
	if(numberRows*8>numberColumns)
	  method=ClpSolve::usePrimal; // switch off sprint
	// but make lightweight
	if(numberRows*10>numberColumns||numberColumns<6000
	   ||(numberRows*20>numberColumns&&!plusMinus))
	  maxSprintPass=5;
      }
    } else if (doSprint==0) {
      method=ClpSolve::usePrimal; // switch off sprint
    }
  }
  if (method==ClpSolve::useDual) {
    double * saveLower=NULL;
    double * saveUpper=NULL;
    if (presolve==ClpSolve::presolveOn) {
      int numberInfeasibilities = model2->tightenPrimalBounds(0.0,0);
      if (numberInfeasibilities) {
	handler_->message(CLP_INFEASIBLE,messages_)
	  <<CoinMessageEol;
	model2 = this;
	presolve=ClpSolve::presolveOff;
      }
    } else if (numberRows_+numberColumns_>5000) {
      // do anyway
      saveLower = new double[numberRows_+numberColumns_];
      CoinMemcpyN(model2->columnLower(),numberColumns_,saveLower);
      CoinMemcpyN(model2->rowLower(),numberRows_,saveLower+numberColumns_);
      saveUpper = new double[numberRows_+numberColumns_];
      CoinMemcpyN(model2->columnUpper(),numberColumns_,saveUpper);
      CoinMemcpyN(model2->rowUpper(),numberRows_,saveUpper+numberColumns_);
      int numberInfeasibilities = model2->tightenPrimalBounds();
      if (numberInfeasibilities) {
	handler_->message(CLP_INFEASIBLE,messages_)
	  <<CoinMessageEol;
        CoinMemcpyN(saveLower,numberColumns_,model2->columnLower());
        CoinMemcpyN(saveLower+numberColumns_,numberRows_,model2->rowLower());
        delete [] saveLower;
        saveLower=NULL;
        CoinMemcpyN(saveUpper,numberColumns_,model2->columnUpper());
        CoinMemcpyN(saveUpper+numberColumns_,numberRows_,model2->rowUpper());
        delete [] saveUpper;
        saveUpper=NULL;
      }
    }
#ifndef COIN_HAS_VOL
    // switch off idiot and volume for now 
    doIdiot=0;
#endif
    // pick up number passes
    int nPasses=0;
    int numberNotE=0;
#ifndef SLIM_CLP
    if ((doIdiot<0&&plusMinus)||doIdiot>0) {
      // See if candidate for idiot
      nPasses=0;
      Idiot info(*model2);
      // Get average number of elements per column
      double ratio  = ((double) numberElements/(double) numberColumns);
      // look at rhs
      int iRow;
      double largest=0.0;
      double smallest = 1.0e30;
      double largestGap=0.0;
      for (iRow=0;iRow<numberRows;iRow++) {
	double value1 = model2->rowLower_[iRow];
	if (value1&&value1>-1.0e31) {
	  largest = CoinMax(largest,fabs(value1));
	  smallest=CoinMin(smallest,fabs(value1));
	}
	double value2 = model2->rowUpper_[iRow];
	if (value2&&value2<1.0e31) {
	  largest = CoinMax(largest,fabs(value2));
	  smallest=CoinMin(smallest,fabs(value2));
	}
	if (value2>value1) {
	  numberNotE++;
	  if (value2>1.0e31||value1<-1.0e31)
	    largestGap = COIN_DBL_MAX;
	  else
	    largestGap = value2-value1;
	}
      }
      if (doIdiot<0) {
	if (numberRows>200&&numberColumns>5000&&ratio>=3.0&&
	    largest/smallest<1.1&&!numberNotE) {
	  nPasses = 71;
	}
      } 
      if (doIdiot>0) {
	nPasses=CoinMax(nPasses,doIdiot);
	if (nPasses>70) {
	  info.setStartingWeight(1.0e3);
	  info.setDropEnoughFeasibility(0.01);
	}
      }
      if (nPasses>20) {
#ifdef COIN_HAS_VOL
        int returnCode = solveWithVolume(model2,nPasses,saveDoIdiot);
        if (!returnCode) {
          time2 = CoinCpuTime();
          timeIdiot = time2-timeX;
          handler_->message(CLP_INTERVAL_TIMING,messages_)
            <<"Idiot Crash"<<timeIdiot<<time2-time1
            <<CoinMessageEol;
          timeX=time2;
        } else {
          nPasses=0;
        }
#else
        nPasses=0;
#endif
      } else {
        nPasses=0;
      }
    }
#endif
    if (doCrash) {
      switch(doCrash) {
	// standard
      case 1:
	model2->crash(1000,1);
	break;
	// As in paper by Solow and Halim (approx)
      case 2:
      case 3:
	model2->crash(model2->dualBound(),0);
	break;
        // Just put free in basis
      case 4:
        model2->crash(0.0,3);
        break;
      }
    }
    if (!nPasses) {
      int saveOptions = model2->specialOptions();
      if (model2->numberRows()>100)
	model2->setSpecialOptions(saveOptions|64); // go as far as possible
      model2->dual(0);
      model2->setSpecialOptions(saveOptions);
    } else if (!numberNotE&&0) {
      // E so we can do in another way
      double * pi = model2->dualRowSolution();
      int i;
      int numberColumns = model2->numberColumns();
      int numberRows = model2->numberRows();
      double * saveObj = new double[numberColumns];
      memcpy(saveObj,model2->objective(),numberColumns*sizeof(double));
      memcpy(model2->dualColumnSolution(),model2->objective(),
	     numberColumns*sizeof(double));
      model2->clpMatrix()->transposeTimes(-1.0,pi,model2->dualColumnSolution());
      memcpy(model2->objective(),model2->dualColumnSolution(),
	     numberColumns*sizeof(double));
      const double * rowsol = model2->primalRowSolution();
      double offset=0.0;
      for (i=0;i<numberRows;i++) {
	offset += pi[i]*rowsol[i];
      }
      double value2;
      model2->getDblParam(ClpObjOffset,value2);
      //printf("Offset %g %g\n",offset,value2);
      model2->setDblParam(ClpObjOffset,value2-offset);
      model2->setPerturbation(51);
      //model2->setRowObjective(pi);
      // zero out pi
      //memset(pi,0,numberRows*sizeof(double));
      // Could put some in basis - only partially tested
      model2->allSlackBasis(); 
      //model2->factorization()->maximumPivots(200);
      //model2->setLogLevel(63);
      // solve
      model2->dual(0);
      model2->setDblParam(ClpObjOffset,value2);
      memcpy(model2->objective(),saveObj,numberColumns*sizeof(double));
      // zero out pi
      //memset(pi,0,numberRows*sizeof(double));
      //model2->setRowObjective(pi);
      delete [] saveObj;
      //model2->dual(0);
      model2->setPerturbation(50);
      model2->primal();
    } else {
      // solve
      model2->setPerturbation(100);
      model2->dual(2);
      model2->setPerturbation(50);
      model2->dual(0);
    }
    if (saveLower) {
      CoinMemcpyN(saveLower,numberColumns_,model2->columnLower());
      CoinMemcpyN(saveLower+numberColumns_,numberRows_,model2->rowLower());
      delete [] saveLower;
      saveLower=NULL;
      CoinMemcpyN(saveUpper,numberColumns_,model2->columnUpper());
      CoinMemcpyN(saveUpper+numberColumns_,numberRows_,model2->rowUpper());
      delete [] saveUpper;
      saveUpper=NULL;
    }
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Dual"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  } else if (method==ClpSolve::usePrimal) {
#ifndef SLIM_CLP
    if (doIdiot) {
      int nPasses=0;
      Idiot info(*model2);
      // Get average number of elements per column
      double ratio  = ((double) numberElements/(double) numberColumns);
      // look at rhs
      int iRow;
      double largest=0.0;
      double smallest = 1.0e30;
      double largestGap=0.0;
      int numberNotE=0;
      for (iRow=0;iRow<numberRows;iRow++) {
	double value1 = model2->rowLower_[iRow];
	if (value1&&value1>-1.0e31) {
	  largest = CoinMax(largest,fabs(value1));
	  smallest=CoinMin(smallest,fabs(value1));
	}
	double value2 = model2->rowUpper_[iRow];
	if (value2&&value2<1.0e31) {
	  largest = CoinMax(largest,fabs(value2));
	  smallest=CoinMin(smallest,fabs(value2));
	}
	if (value2>value1) {
	  numberNotE++;
	  if (value2>1.0e31||value1<-1.0e31)
	    largestGap = COIN_DBL_MAX;
	  else
	    largestGap = value2-value1;
	}
      }
      bool increaseSprint=plusMinus;
      if (!plusMinus) {
        // If 90% +- 1 then go for sprint
        if (statistics[0]>=0&&10*statistics[2]<statistics[0]+statistics[1])
          increaseSprint=true;
      }
      bool tryIt= numberRows>200&&numberColumns>2000&&numberColumns>2*numberRows;
      if (numberRows<1000&&numberColumns<3000)
        tryIt=false;
      if (tryIt) {
	if (increaseSprint) {
	  info.setStartingWeight(1.0e3);
	  info.setReduceIterations(6);
	  // also be more lenient on infeasibilities
	  info.setDropEnoughFeasibility(0.5*info.getDropEnoughFeasibility());
	  info.setDropEnoughWeighted(-2.0);
	  if (largest/smallest>2.0) {
	    nPasses = 10+numberColumns/100000;
	    nPasses = CoinMin(nPasses,50);
	    nPasses = CoinMax(nPasses,15);
	    if (numberRows>20000&&nPasses>5) {
	      // Might as well go for it
	      nPasses = CoinMax(nPasses,71);
	    } else if (numberRows>2000&&nPasses>5) {
	      nPasses = CoinMax(nPasses,50);
	    } else if (numberElements<3*numberColumns) {
	      nPasses=CoinMin(nPasses,10); // probably not worh it
              if (doIdiot<0)
                info.setLightweight(1); // say lightweight idiot
	    } else {
              if (doIdiot<0)
                info.setLightweight(1); // say lightweight idiot
            }
	  } else if (largest/smallest>1.01||numberElements<=3*numberColumns) {
	    nPasses = 10+numberColumns/1000;
	    nPasses = CoinMin(nPasses,100);
	    nPasses = CoinMax(nPasses,30);
	    if (numberRows>25000) {
	      // Might as well go for it
	      nPasses = CoinMax(nPasses,71);
	    }
	    if (!largestGap)
	      nPasses *= 2;
	  } else {
	    nPasses = 10+numberColumns/1000;
	    nPasses = CoinMin(nPasses,200);
	    nPasses = CoinMax(nPasses,100);
	    info.setStartingWeight(1.0e-1);
	    info.setReduceIterations(6);
	    if (!largestGap)
	      nPasses *= 2;
	    //info.setFeasibilityTolerance(1.0e-7);
	  }
	  // If few passes - don't bother
	  if (nPasses<=5)
	    nPasses=0;
	} else {
	  if (doIdiot<0)
	    info.setLightweight(1); // say lightweight idiot
	  if (largest/smallest>1.01||numberNotE||statistics[2]>statistics[0]+statistics[1]) {
	    if (numberRows>25000||numberColumns>5*numberRows) {
	      nPasses = 50;
	    } else if (numberColumns>4*numberRows) {
	      nPasses = 20;
	    } else {
	      nPasses=5;
	    }
	  } else {
	    if (numberRows>25000||numberColumns>5*numberRows) {
	      nPasses = 50;
	      info.setLightweight(0); // say not lightweight idiot
	    } else if (numberColumns>4*numberRows) {
	      nPasses = 20;
	    } else {
	      nPasses=15;
	    }
	  }
	  if (ratio<3.0) { 
	    nPasses=(int) ((ratio*(double) nPasses)/4.0); // probably not worh it
	  } else {
	    nPasses = CoinMax(nPasses,5);
	  }
	  if (numberRows>25000&&nPasses>5) {
	    // Might as well go for it
	    nPasses = CoinMax(nPasses,71);
	  } else if (increaseSprint) {
	    nPasses *= 2;
	    nPasses=CoinMin(nPasses,71);
	  } else if (nPasses==5&&ratio>5.0) {
	    nPasses = (int) (((double) nPasses)*(ratio/5.0)); // increase if lots of elements per column
          }
	  if (nPasses<=5)
	    nPasses=0;
	  //info.setStartingWeight(1.0e-1);
	}
      }
      if (doIdiot>0) {
	// pick up number passes
	nPasses=options.getExtraInfo(1);
	if (nPasses>70) {
	  info.setStartingWeight(1.0e3);
	  info.setReduceIterations(6);
	  // also be more lenient on infeasibilities
	  info.setDropEnoughFeasibility(0.5*info.getDropEnoughFeasibility());
	  info.setDropEnoughWeighted(-2.0);
	} else if (nPasses>=50) {
	  info.setStartingWeight(1.0e3);
	  //info.setReduceIterations(6);
	} 
	// For experimenting
	if (nPasses<70&&(nPasses%10)>0&&(nPasses%10)<4) {
	  info.setStartingWeight(1.0e3);
	  info.setLightweight(nPasses%10); // special testing
          printf("warning - odd lightweight %d\n",nPasses%10);
	  //info.setReduceIterations(6);
	}
      }
      if (nPasses) {
	doCrash=0;
#if 0
	double * solution = model2->primalColumnSolution();
	int iColumn;
	double * saveLower = new double[numberColumns];
	memcpy(saveLower,model2->columnLower(),numberColumns*sizeof(double));
	double * saveUpper = new double[numberColumns];
	memcpy(saveUpper,model2->columnUpper(),numberColumns*sizeof(double));
	printf("doing tighten before idiot\n");
	model2->tightenPrimalBounds();
	// Move solution
	double * columnLower = model2->columnLower();
	double * columnUpper = model2->columnUpper();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (columnLower[iColumn]>0.0)
	    solution[iColumn]=columnLower[iColumn];
	  else if (columnUpper[iColumn]<0.0)
	    solution[iColumn]=columnUpper[iColumn];
	  else
	    solution[iColumn]=0.0;
	}
	memcpy(columnLower,saveLower,numberColumns*sizeof(double));
	memcpy(columnUpper,saveUpper,numberColumns*sizeof(double));
	delete [] saveLower;
	delete [] saveUpper;
#else
	// Allow for crossover
        //if (doIdiot>0)
          info.setStrategy(512|info.getStrategy());
	// Allow for scaling
	info.setStrategy(32|info.getStrategy());
	info.crash(nPasses,model2->messageHandler(),model2->messagesPointer());
#endif
	time2 = CoinCpuTime();
	timeIdiot = time2-timeX;
	handler_->message(CLP_INTERVAL_TIMING,messages_)
	  <<"Idiot Crash"<<timeIdiot<<time2-time1
	  <<CoinMessageEol;
	timeX=time2;
      }
    }
#endif
    // ?
    if (doCrash) {
      switch(doCrash) {
	// standard
      case 1:
	model2->crash(1000,1);
	break;
	// As in paper by Solow and Halim (approx)
      case 2:
	model2->crash(model2->dualBound(),0);
	break;
	// My take on it
      case 3:
	model2->crash(model2->dualBound(),-1);
	break;
        // Just put free in basis
      case 4:
        model2->crash(0.0,3);
        break;
      }
    }
#ifndef SLIM_CLP
    if (doSlp>0&&objective_->type()==2) {
      model2->nonlinearSLP(doSlp,1.0e-5);
    }
#endif
    model2->primal(primalStartup);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Primal"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
  } else if (method==ClpSolve::usePrimalorSprint) {
    // Sprint
    /*
      This driver implements what I called Sprint when I introduced the idea
      many years ago.  Cplex calls it "sifting" which I think is just as silly.
      When I thought of this trivial idea
      it reminded me of an LP code of the 60's called sprint which after
      every factorization took a subset of the matrix into memory (all
      64K words!) and then iterated very fast on that subset.  On the
      problems of those days it did not work very well, but it worked very
      well on aircrew scheduling problems where there were very large numbers
      of columns all with the same flavor.
    */
    
    /* The idea works best if you can get feasible easily.  To make it
       more general we can add in costed slacks */
    
    int originalNumberColumns = model2->numberColumns();
    int numberRows = model2->numberRows();
    
    // We will need arrays to choose variables.  These are too big but ..
    double * weight = new double [numberRows+originalNumberColumns];
    int * sort = new int [numberRows+originalNumberColumns];
    int numberSort=0;
    // We are going to add slacks to get feasible.
    // initial list will just be artificials
    // first we will set all variables as close to zero as possible
    int iColumn;
    const double * columnLower = model2->columnLower();
    const double * columnUpper = model2->columnUpper();
    double * columnSolution = model2->primalColumnSolution();
    
    for (iColumn=0;iColumn<originalNumberColumns;iColumn++) {
      double value =0.0;
      if (columnLower[iColumn]>0.0)
	value = columnLower[iColumn];
      else if (columnUpper[iColumn]<0.0)
	value = columnUpper[iColumn];
      columnSolution[iColumn]=value;
    }
    // now see what that does to row solution
    double * rowSolution = model2->primalRowSolution();
    memset (rowSolution,0,numberRows*sizeof(double));
    model2->times(1.0,columnSolution,rowSolution);
    
    int * addStarts = new int [numberRows+1];
    int * addRow = new int[numberRows];
    double * addElement = new double[numberRows];
    const double * lower = model2->rowLower();
    const double * upper = model2->rowUpper();
    addStarts[0]=0;
    int numberArtificials=0;
    double * addCost = new double [numberRows];
    const double penalty=1.0e8*optimizationDirection_;
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
    model2->addColumns(numberArtificials,NULL,NULL,addCost,
		       addStarts,addRow,addElement);
    delete [] addStarts;
    delete [] addRow;
    delete [] addElement;
    delete [] addCost;
    // look at rhs to see if to perturb
    double largest=0.0;
    double smallest = 1.0e30;
    for (iRow=0;iRow<numberRows;iRow++) {
      double value;
      value = fabs(model2->rowLower_[iRow]);
      if (value&&value<1.0e30) {
	largest = CoinMax(largest,value);
	smallest=CoinMin(smallest,value);
      }
      value = fabs(model2->rowUpper_[iRow]);
      if (value&&value<1.0e30) {
	largest = CoinMax(largest,value);
	smallest=CoinMin(smallest,value);
      }
    }
    double * saveLower = NULL;
    double * saveUpper = NULL;
    if (largest<2.01*smallest) {
      // perturb - so switch off standard
      model2->setPerturbation(100);
      saveLower = new double[numberRows];
      memcpy(saveLower,model2->rowLower_,numberRows*sizeof(double));
      saveUpper = new double[numberRows];
      memcpy(saveUpper,model2->rowUpper_,numberRows*sizeof(double));
      double * lower = model2->rowLower();
      double * upper = model2->rowUpper();
      for (iRow=0;iRow<numberRows;iRow++) {
	double lowerValue=lower[iRow], upperValue=upper[iRow];
	double value = CoinDrand48();
	if (upperValue>lowerValue+primalTolerance_) {
	  if (lowerValue>-1.0e20&&lowerValue)
	    lowerValue -= value * 1.0e-4*fabs(lowerValue); 
	  if (upperValue<1.0e20&&upperValue)
	    upperValue += value * 1.0e-4*fabs(upperValue); 
	} else if (upperValue>0.0) {
	  upperValue -= value * 1.0e-4*fabs(lowerValue); 
	  lowerValue -= value * 1.0e-4*fabs(lowerValue); 
	} else if (upperValue<0.0) {
	  upperValue += value * 1.0e-4*fabs(lowerValue); 
	  lowerValue += value * 1.0e-4*fabs(lowerValue); 
	} else {
	}
	lower[iRow]=lowerValue;
	upper[iRow]=upperValue;
      }
    }
    int i;
    // Set up initial list
    if (numberArtificials) {
      numberSort=numberArtificials;
      for (i=0;i<numberSort;i++)
	sort[i] = i+originalNumberColumns;
    } else {
      numberSort = CoinMin(numberRows_,numberColumns_);
      for (i=0;i<numberSort;i++)
	sort[i] = i;
    }
    
    // redo as will have changed
    columnLower = model2->columnLower();
    columnUpper = model2->columnUpper();
    int numberColumns = model2->numberColumns();
    double * fullSolution = model2->primalColumnSolution();
    
    // Just do this number of passes in Sprint
    if (doSprint>0)
      maxSprintPass=options.getExtraInfo(1);
    int iPass;
    double lastObjective=1.0e31;
    // It will be safe to allow dense
    model2->setInitialDenseFactorization(true);
    
    // Just take this number of columns in small problem
    int smallNumberColumns = CoinMin(3*numberRows,numberColumns);
    smallNumberColumns = CoinMax(smallNumberColumns,3000);
    smallNumberColumns = CoinMin(smallNumberColumns,numberColumns);
    //int smallNumberColumns = CoinMin(12*numberRows/10,numberColumns);
    //smallNumberColumns = CoinMax(smallNumberColumns,3000);
    //smallNumberColumns = CoinMax(smallNumberColumns,numberRows+1000);
    // We will be using all rows
    int * whichRows = new int [numberRows];
    for (iRow=0;iRow<numberRows;iRow++)
      whichRows[iRow]=iRow;
    double originalOffset;
    model2->getDblParam(ClpObjOffset,originalOffset);
    int totalIterations=0;
    for (iPass=0;iPass<maxSprintPass;iPass++) {
      //printf("Bug until submodel new version\n");
      //CoinSort_2(sort,sort+numberSort,weight);
      // Create small problem
      ClpSimplex small(model2,numberRows,whichRows,numberSort,sort);
      small.setPerturbation(model2->perturbation());
      // now see what variables left out do to row solution
      double * rowSolution = model2->primalRowSolution();
      double * sumFixed = new double[numberRows];
      memset (sumFixed,0,numberRows*sizeof(double));
      int iRow,iColumn;
      // zero out ones in small problem
      for (iColumn=0;iColumn<numberSort;iColumn++) {
	int kColumn = sort[iColumn];
	fullSolution[kColumn]=0.0;
      }
      // Get objective offset
      double offset=0.0;
      const double * objective = model2->objective();
      for (iColumn=0;iColumn<numberColumns;iColumn++) 
	offset += fullSolution[iColumn]*objective[iColumn];
      small.setDblParam(ClpObjOffset,originalOffset-offset);
      model2->times(1.0,fullSolution,sumFixed);
      
      double * lower = small.rowLower();
      double * upper = small.rowUpper();
      for (iRow=0;iRow<numberRows;iRow++) {
	if (lower[iRow]>-1.0e50) 
	  lower[iRow] -= sumFixed[iRow];
	if (upper[iRow]<1.0e50)
	  upper[iRow] -= sumFixed[iRow];
	rowSolution[iRow] -= sumFixed[iRow];
      }
      delete [] sumFixed;
      // Solve 
      if (interrupt)
	currentModel = &small;
      small.primal(1);
      totalIterations += small.numberIterations();
      // move solution back
      const double * solution = small.primalColumnSolution();
      for (iColumn=0;iColumn<numberSort;iColumn++) {
	int kColumn = sort[iColumn];
	model2->setColumnStatus(kColumn,small.getColumnStatus(iColumn));
	fullSolution[kColumn]=solution[iColumn];
      }
      for (iRow=0;iRow<numberRows;iRow++) 
	model2->setRowStatus(iRow,small.getRowStatus(iRow));
      memcpy(model2->primalRowSolution(),small.primalRowSolution(),
	     numberRows*sizeof(double));
      // get reduced cost for large problem
      memcpy(weight,model2->objective(),numberColumns*sizeof(double));
      model2->transposeTimes(-1.0,small.dualRowSolution(),weight);
      int numberNegative=0;
      double sumNegative = 0.0;
      // now massage weight so all basic in plus good djs
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double dj = weight[iColumn]*optimizationDirection_;
	double value = fullSolution[iColumn];
	if (model2->getColumnStatus(iColumn)==ClpSimplex::basic) 
	  dj = -1.0e50;
	else if (dj<0.0&&value<columnUpper[iColumn])
	  dj = dj;
	else if (dj>0.0&&value>columnLower[iColumn])
	  dj = -dj;
	else if (columnUpper[iColumn]>columnLower[iColumn])
	  dj = fabs(dj);
	else
	  dj = 1.0e50;
	weight[iColumn] = dj;
	if (dj<-dualTolerance_&&dj>-1.0e50) {
	  numberNegative++;
	  sumNegative -= dj;
	}
	sort[iColumn] = iColumn;
      }
      handler_->message(CLP_SPRINT,messages_)
	<<iPass+1<<small.numberIterations()<<small.objectiveValue()<<sumNegative
	<<numberNegative
	<<CoinMessageEol;
      if ((small.objectiveValue()*optimizationDirection_>lastObjective-1.0e-7&&iPass>5)||
	  !small.numberIterations()||
	  iPass==maxSprintPass-1||small.status()==3) {
	
	break; // finished
      } else {
	lastObjective = small.objectiveValue()*optimizationDirection_;
	// sort
	CoinSort_2(weight,weight+numberColumns,sort);
	numberSort = smallNumberColumns;
      }
    }
    if (interrupt) 
      currentModel = model2;
    for (i=0;i<numberArtificials;i++)
      sort[i] = i + originalNumberColumns;
    model2->deleteColumns(numberArtificials,sort);
    delete [] weight;
    delete [] sort;
    delete [] whichRows;
    if (saveLower) {
      // unperturb and clean
      for (iRow=0;iRow<numberRows;iRow++) {
	double diffLower = saveLower[iRow]-model2->rowLower_[iRow];
	double diffUpper = saveUpper[iRow]-model2->rowUpper_[iRow];
	model2->rowLower_[iRow]=saveLower[iRow];
	model2->rowUpper_[iRow]=saveUpper[iRow];
	if (diffLower) 
	  assert (!diffUpper||fabs(diffLower-diffUpper)<1.0e-5);
	else
	  diffLower = diffUpper;
	model2->rowActivity_[iRow] += diffLower;
      }
      delete [] saveLower;
      delete [] saveUpper;
    }
    model2->primal(0);
    model2->setPerturbation(savePerturbation);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Sprint"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    model2->setNumberIterations(model2->numberIterations()+totalIterations);
  } else if (method==ClpSolve::useBarrier||method==ClpSolve::useBarrierNoCross) {
#ifndef SLIM_CLP
    //printf("***** experimental pretty crude barrier\n");
    //#define SAVEIT 1
#ifndef SAVEIT
#define BORROW
#endif
#ifdef BORROW
    ClpInterior barrier;
    barrier.borrowModel(*model2);
#else
    ClpInterior barrier(*model2);
#endif
    if (interrupt) 
      currentModel2 = &barrier;
    int barrierOptions = options.getSpecialOption(4);
    bool aggressiveGamma=false;
    bool scale=false;
    bool doKKT=false;
    if (barrierOptions&32) {
      barrierOptions &= ~32;
      doKKT=true;
    }
    if (barrierOptions&16) {
      barrierOptions &= ~16;
      aggressiveGamma=true;
    }
    if (barrierOptions&8) {
      barrierOptions &= ~8;
      scale=true;
    }
    // If quadratic force KKT
#ifndef NO_RTTI
    ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(barrier.objectiveAsObject()));
#else
    ClpQuadraticObjective * quadraticObj = NULL;
    if (objective_->type()==2)
      quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
#ifndef FAST_BARRIER
    if (!numberBarrier)
      std::cout<<"Warning - the default ordering is just on row counts! "
	       <<"The factorization is being improved"<<std::endl;
    numberBarrier++;
#endif
    if (quadraticObj) {
      doKKT=true;
    }
    switch (barrierOptions) {
    case 0:
    default:
      if (!doKKT) {
	ClpCholeskyBase * cholesky = new ClpCholeskyBase();
	barrier.setCholesky(cholesky);
      } else {
	ClpCholeskyBase * cholesky = new ClpCholeskyBase();
	cholesky->setKKT(true);
	barrier.setCholesky(cholesky);
      }
      break;
    case 1:
      if (!doKKT) {
	ClpCholeskyDense * cholesky = new ClpCholeskyDense();
	barrier.setCholesky(cholesky);
      } else {
	ClpCholeskyDense * cholesky = new ClpCholeskyDense();
	cholesky->setKKT(true);
	barrier.setCholesky(cholesky);
      }
      break;
#ifdef WSSMP_BARRIER
    case 2:
      {
	ClpCholeskyWssmp * cholesky = new ClpCholeskyWssmp(CoinMax(100,model2->numberRows()/10));
	barrier.setCholesky(cholesky);
	assert (!doKKT);
      }
      break;
    case 3:
      if (!doKKT) {
	ClpCholeskyWssmp * cholesky = new ClpCholeskyWssmp();
	barrier.setCholesky(cholesky);
      } else {
	ClpCholeskyWssmpKKT * cholesky = new ClpCholeskyWssmpKKT(CoinMax(100,model2->numberRows()/10));
	barrier.setCholesky(cholesky);
      }
      break;
#endif
#ifdef UFL_BARRIER
    case 4:
      if (!doKKT) {
	ClpCholeskyUfl * cholesky = new ClpCholeskyUfl();
	barrier.setCholesky(cholesky);
      } else {
	ClpCholeskyUfl * cholesky = new ClpCholeskyUfl();
	cholesky->setKKT(true);
	barrier.setCholesky(cholesky);
      }
      break;
#endif
#ifdef TAUCS_BARRIER
    case 5:
      {
	ClpCholeskyTaucs * cholesky = new ClpCholeskyTaucs();
	barrier.setCholesky(cholesky);
	assert (!doKKT);
      }
      break;
#endif
    }
    int numberRows = model2->numberRows();
    int numberColumns = model2->numberColumns();
    int saveMaxIts = model2->maximumIterations();
    if (saveMaxIts<1000) {
      barrier.setMaximumBarrierIterations(saveMaxIts);
      model2->setMaximumIterations(1000000);
    }
#ifndef SAVEIT
    //barrier.setGamma(1.0e-8);
    //barrier.setDelta(1.0e-8);
    barrier.setDiagonalPerturbation(1.0e-14);
    if (aggressiveGamma) {
      barrier.setGamma(1.0e-3);
      barrier.setDelta(1.0e-3);
    }
    if (scale)
      barrier.scaling(1);
    else
      barrier.scaling(0);
    barrier.primalDual();
#elif SAVEIT==1
    barrier.primalDual();
#else
    model2->restoreModel("xx.save");
    // move solutions
    CoinMemcpyN(model2->primalRowSolution(),
		numberRows,barrier.primalRowSolution());
    CoinMemcpyN(model2->dualRowSolution(),
		numberRows,barrier.dualRowSolution());
    CoinMemcpyN(model2->primalColumnSolution(),
		numberColumns,barrier.primalColumnSolution());
    CoinMemcpyN(model2->dualColumnSolution(),
		numberColumns,barrier.dualColumnSolution());
#endif
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Barrier"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    int maxIts = barrier.maximumBarrierIterations();
    int barrierStatus=barrier.status();
    double gap = barrier.complementarityGap();
    // get which variables are fixed
    double * saveLower=NULL;
    double * saveUpper=NULL;
    ClpPresolve pinfo2;
    ClpSimplex * saveModel2=NULL;
    int numberFixed = barrier.numberFixed();
    if (numberFixed*20>barrier.numberRows()&&numberFixed>5000&&0) {
      // may as well do presolve
      int numberRows = barrier.numberRows();
      int numberColumns = barrier.numberColumns();
      int numberTotal = numberRows+numberColumns;
      saveLower = new double [numberTotal];
      saveUpper = new double [numberTotal];
      memcpy(saveLower,barrier.columnLower(),numberColumns*sizeof(double));
      memcpy(saveLower+numberColumns,barrier.rowLower(),numberRows*sizeof(double));
      memcpy(saveUpper,barrier.columnUpper(),numberColumns*sizeof(double));
      memcpy(saveUpper+numberColumns,barrier.rowUpper(),numberRows*sizeof(double));
      barrier.fixFixed();
      saveModel2=model2;
    } else if (numberFixed) {
      // Set fixed to bounds (may have restored earlier solution)
      barrier.fixFixed(false);
    }
#ifdef BORROW    
    barrier.returnModel(*model2);
    double * rowPrimal = new double [numberRows];
    double * columnPrimal = new double [numberColumns];
    double * rowDual = new double [numberRows];
    double * columnDual = new double [numberColumns];
    // move solutions other way
    CoinMemcpyN(model2->primalRowSolution(),
		numberRows,rowPrimal);
    CoinMemcpyN(model2->dualRowSolution(),
		numberRows,rowDual);
    CoinMemcpyN(model2->primalColumnSolution(),
		numberColumns,columnPrimal);
    CoinMemcpyN(model2->dualColumnSolution(),
		  numberColumns,columnDual);
#else
    double * rowPrimal = barrier.primalRowSolution();
    double * columnPrimal = barrier.primalColumnSolution();
    double * rowDual = barrier.dualRowSolution();
    double * columnDual = barrier.dualColumnSolution();
    // move solutions
    CoinMemcpyN(rowPrimal,
		numberRows,model2->primalRowSolution());
    CoinMemcpyN(rowDual,
		numberRows,model2->dualRowSolution());
    CoinMemcpyN(columnPrimal,
		numberColumns,model2->primalColumnSolution());
    CoinMemcpyN(columnDual,
		  numberColumns,model2->dualColumnSolution());
#endif
    if (saveModel2) {
      // do presolve
      model2 = pinfo2.presolvedModel(*model2,1.0e-8,
				    false,5,true);
      if (!model2) {
	model2=saveModel2;
	saveModel2=NULL;
	delete [] saveLower;
	saveLower=NULL;
	delete [] saveUpper;
	saveUpper=NULL;
      }
    }
    if (method==ClpSolve::useBarrier) {
      if (maxIts&&barrierStatus<4&&!quadraticObj) {
	//printf("***** crossover - needs more thought on difficult models\n");
#if SAVEIT==1
	model2->ClpSimplex::saveModel("xx.save");
#endif
	// make sure no status left
	model2->createStatus();
	// solve
	model2->setPerturbation(100);
#if 1
	// throw some into basis 
	{
	  int numberRows = model2->numberRows();
	  int numberColumns = model2->numberColumns();
	  double * dsort = new double[numberColumns];
	  int * sort = new int[numberColumns];
	  int n=0;
	  const double * columnLower = model2->columnLower();
	  const double * columnUpper = model2->columnUpper();
	  double * primalSolution = model2->primalColumnSolution();
	  const double * dualSolution = model2->dualColumnSolution();
	  double tolerance = 10.0*primalTolerance_;
	  int i;
	  for ( i=0;i<numberRows;i++) 
	    model2->setRowStatus(i,superBasic);
	  for ( i=0;i<numberColumns;i++) {
	    double distance = CoinMin(columnUpper[i]-primalSolution[i],
				      primalSolution[i]-columnLower[i]);
	    if (distance>tolerance) {
	      if (fabs(dualSolution[i])<1.0e-5)
		distance *= 100.0;
	      dsort[n]=-distance;
	      sort[n++]=i;
	      model2->setStatus(i,superBasic);
	    } else if (distance>primalTolerance_) {
	      model2->setStatus(i,superBasic);
	    } else if (primalSolution[i]<=columnLower[i]+primalTolerance_) {
	      model2->setStatus(i,atLowerBound);
	      primalSolution[i]=columnLower[i];
	    } else {
	      model2->setStatus(i,atUpperBound);
	      primalSolution[i]=columnUpper[i];
	    }
	  }
	  CoinSort_2(dsort,dsort+n,sort);
	  n = CoinMin(numberRows,n);
	  for ( i=0;i<n;i++) {
	    int iColumn = sort[i];
	    model2->setStatus(iColumn,basic);
	  }
	  delete [] sort;
	  delete [] dsort;
	}
	if (gap<1.0e-3*((double) (numberRows+numberColumns))) {
	  int numberRows = model2->numberRows();
	  int numberColumns = model2->numberColumns();
	  // just primal values pass
	  double saveScale = model2->objectiveScale();
	  model2->setObjectiveScale(1.0e-3);
	  model2->primal(2);
	  model2->setObjectiveScale(saveScale);
	  // save primal solution and copy back dual
	  CoinMemcpyN(model2->primalRowSolution(),
		      numberRows,rowPrimal);
	  CoinMemcpyN(rowDual,
		      numberRows,model2->dualRowSolution());
	  CoinMemcpyN(model2->primalColumnSolution(),
		      numberColumns,columnPrimal);
	  CoinMemcpyN(columnDual,
		      numberColumns,model2->dualColumnSolution());
	  //model2->primal(1);
	  // clean up reduced costs and flag variables
	  {
	    double * dj = model2->dualColumnSolution();
	    double * cost = model2->objective();
	    double * saveCost = new double[numberColumns];
	    memcpy(saveCost,cost,numberColumns*sizeof(double));
	    double * saveLower = new double[numberColumns];
	    double * lower = model2->columnLower();
	    memcpy(saveLower,lower,numberColumns*sizeof(double));
	    double * saveUpper = new double[numberColumns];
	    double * upper = model2->columnUpper();
	    memcpy(saveUpper,upper,numberColumns*sizeof(double));
	    int i;
	    double tolerance = 10.0*dualTolerance_;
	    for ( i=0;i<numberColumns;i++) {
	      if (model2->getStatus(i)==basic) {
		dj[i]=0.0;
	      } else if (model2->getStatus(i)==atLowerBound) {
		if (optimizationDirection_*dj[i]<tolerance) {
		  if (optimizationDirection_*dj[i]<0.0) {
		    //if (dj[i]<-1.0e-3)
		    //printf("bad dj at lb %d %g\n",i,dj[i]);
		    cost[i] -= dj[i];
		    dj[i]=0.0;
		  }
		} else {
		  upper[i]=lower[i];
		}
	      } else if (model2->getStatus(i)==atUpperBound) {
		if (optimizationDirection_*dj[i]>tolerance) {
		  if (optimizationDirection_*dj[i]>0.0) {
		    //if (dj[i]>1.0e-3)
		    //printf("bad dj at ub %d %g\n",i,dj[i]);
		    cost[i] -= dj[i];
		    dj[i]=0.0;
		  }
		} else {
		  lower[i]=upper[i];
		}
	      }
	    }
	    // just dual values pass
	    //model2->setLogLevel(63);
	    //model2->setFactorizationFrequency(1);
	    model2->dual(2);
	    memcpy(cost,saveCost,numberColumns*sizeof(double));
	    delete [] saveCost;
	    memcpy(lower,saveLower,numberColumns*sizeof(double));
	    delete [] saveLower;
	    memcpy(upper,saveUpper,numberColumns*sizeof(double));
	    delete [] saveUpper;
	  }
	  // and finish
	  // move solutions
	  CoinMemcpyN(rowPrimal,
		      numberRows,model2->primalRowSolution());
	  CoinMemcpyN(columnPrimal,
		      numberColumns,model2->primalColumnSolution());
	}
	double saveScale = model2->objectiveScale();
	model2->setObjectiveScale(1.0e-3);
	model2->primal(2);
	model2->setObjectiveScale(saveScale);
	model2->primal(1);
#else
	// just primal
	model2->primal(1);
#endif
      } else if (barrierStatus==4) {
	// memory problems
	model2->setPerturbation(savePerturbation);
	model2->createStatus();
	model2->dual();
      } else if (maxIts&&quadraticObj) {
	// make sure no status left
	model2->createStatus();
	// solve
	model2->setPerturbation(100);
	model2->reducedGradient(1);
      }
    }
    model2->setMaximumIterations(saveMaxIts);
#ifdef BORROW
    delete [] rowPrimal;
    delete [] columnPrimal;
    delete [] rowDual;
    delete [] columnDual;
#endif
    if (saveLower) {
      pinfo2.postsolve(true);
      delete model2;
      model2=saveModel2;
      int numberRows = model2->numberRows();
      int numberColumns = model2->numberColumns();
      memcpy(model2->columnLower(),saveLower,numberColumns*sizeof(double));
      memcpy(model2->rowLower(),saveLower+numberColumns,numberRows*sizeof(double));
      delete [] saveLower;
      memcpy(model2->columnUpper(),saveUpper,numberColumns*sizeof(double));
      memcpy(model2->rowUpper(),saveUpper+numberColumns,numberRows*sizeof(double));
      delete [] saveUpper;
      model2->primal(1);
    }
    model2->setPerturbation(savePerturbation);
    time2 = CoinCpuTime();
    timeCore = time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Crossover"<<timeCore<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
#else
    abort();
#endif
  } else {
    assert (method!=ClpSolve::automatic); // later
    time2=0.0;
  }
  if (saveMatrix) {
    // delete and replace
    delete model2->clpMatrix();
    model2->replaceMatrix(saveMatrix);
  }
  numberIterations = model2->numberIterations();
  finalStatus=model2->status();
  if (presolve==ClpSolve::presolveOn) {
    int saveLevel = logLevel();
    if ((specialOptions_&1024)==0)
      setLogLevel(CoinMin(1,saveLevel));
    else
      setLogLevel(CoinMin(0,saveLevel));
    pinfo.postsolve(true);
    factorization_->areaFactor(model2->factorization()->adjustedAreaFactor());
    time2 = CoinCpuTime();
    timePresolve += time2-timeX;
    handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Postsolve"<<time2-timeX<<time2-time1
      <<CoinMessageEol;
    timeX=time2;
    if (!presolveToFile)
      delete model2;
    if (interrupt)
      currentModel = this;
    // checkSolution(); already done by postSolve
    setLogLevel(saveLevel);
    if (finalStatus!=3&&(finalStatus||status()==-1)) {
      int savePerturbation = perturbation();
      setPerturbation(100);
      if (finalStatus==2) {
        // unbounded - get feasible first
        double save = optimizationDirection_;
        optimizationDirection_=0.0;
        primal(1);
        optimizationDirection_=save;
      }
      primal(1);
      setPerturbation(savePerturbation);
      numberIterations += this->numberIterations();
      finalStatus=status();
      time2 = CoinCpuTime();
      handler_->message(CLP_INTERVAL_TIMING,messages_)
      <<"Cleanup"<<time2-timeX<<time2-time1
      <<CoinMessageEol;
      timeX=time2;
    }
  }
  setMaximumIterations(saveMaxIterations);
  std::string statusMessage[]={"Unknown","Optimal","PrimalInfeasible","DualInfeasible","Stopped",
			       "Errors","User stopped"};
  assert (finalStatus>=-1&&finalStatus<=5);
  handler_->message(CLP_TIMING,messages_)
    <<statusMessage[finalStatus+1]<<objectiveValue()<<numberIterations<<time2-time1;
  handler_->printing(presolve==ClpSolve::presolveOn)
    <<timePresolve;
  handler_->printing(timeIdiot)
    <<timeIdiot;
  handler_->message()<<CoinMessageEol;
  if (interrupt) 
    signal(SIGINT,saveSignal);
  perturbation_=savePerturbation;
  scalingFlag_=saveScaling;
  return finalStatus;
}
// General solve
int 
ClpSimplex::initialSolve()
{
  // Default so use dual
  ClpSolve options;
  return initialSolve(options);
}
// General dual solve
int 
ClpSimplex::initialDualSolve()
{
  ClpSolve options;
  // Use dual
  options.setSolveType(ClpSolve::useDual);
  return initialSolve(options);
}
// General dual solve
int 
ClpSimplex::initialPrimalSolve()
{
  ClpSolve options;
  // Use primal
  options.setSolveType(ClpSolve::usePrimal);
  return initialSolve(options);
}
// Default constructor
ClpSolve::ClpSolve (  )
{
  method_ = automatic;
  presolveType_=presolveOn;
  numberPasses_=5;
  int i;
  for (i=0;i<6;i++)
    options_[i]=0;
  for (i=0;i<6;i++)
    extraInfo_[i]=-1;
  for (i=0;i<2;i++)
    independentOptions_[i]=0;
  // But switch off slacks
  independentOptions_[1]=512;
  // Substitute up to 3
  independentOptions_[2]=3;
  
}

// Copy constructor. 
ClpSolve::ClpSolve(const ClpSolve & rhs)
{
  method_ = rhs.method_;
  presolveType_=rhs.presolveType_;
  numberPasses_=rhs.numberPasses_;
  int i;
  for ( i=0;i<6;i++)
    options_[i]=rhs.options_[i];
  for ( i=0;i<6;i++)
    extraInfo_[i]=rhs.extraInfo_[i];
  for (i=0;i<3;i++)
    independentOptions_[i]=rhs.independentOptions_[i];
}
// Assignment operator. This copies the data
ClpSolve & 
ClpSolve::operator=(const ClpSolve & rhs)
{
  if (this != &rhs) {
    method_ = rhs.method_;
    presolveType_=rhs.presolveType_;
    numberPasses_=rhs.numberPasses_;
    int i;
    for (i=0;i<6;i++)
      options_[i]=rhs.options_[i];
    for (i=0;i<6;i++)
      extraInfo_[i]=rhs.extraInfo_[i];
    for (i=0;i<3;i++)
      independentOptions_[i]=rhs.independentOptions_[i];
  }
  return *this;

}
// Destructor
ClpSolve::~ClpSolve (  )
{
}
// See header file for deatils
void 
ClpSolve::setSpecialOption(int which,int value,int extraInfo)
{
  options_[which]=value;
  extraInfo_[which]=extraInfo;
}
int 
ClpSolve::getSpecialOption(int which) const
{
  return options_[which];
}

// Solve types
void 
ClpSolve::setSolveType(SolveType method, int extraInfo)
{
  method_=method;
}

ClpSolve::SolveType 
ClpSolve::getSolveType()
{
  return method_;
}

// Presolve types
void 
ClpSolve::setPresolveType(PresolveType amount, int extraInfo)
{
  presolveType_=amount;
  numberPasses_=extraInfo;
}
ClpSolve::PresolveType 
ClpSolve::getPresolveType()
{
  return presolveType_;
}
// Extra info for idiot (or sprint)
int 
ClpSolve::getExtraInfo(int which) const
{
  return extraInfo_[which];
}
int 
ClpSolve::getPresolvePasses() const
{
  return numberPasses_;
}
/* Say to return at once if infeasible,
   default is to solve */
void 
ClpSolve::setInfeasibleReturn(bool trueFalse)
{
  independentOptions_[0]= trueFalse ? 1 : 0;
}
