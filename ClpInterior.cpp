// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.




#include "CoinPragma.hpp"

#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "ClpInterior.hpp"
#include "ClpMatrixBase.hpp"
#ifdef PDCO
#include "ClpLsqr.hpp"
#include "ClpPdcoBase.hpp"
#endif
#include "CoinDenseVector.hpp"
#include "ClpMessage.hpp"
#include "ClpHelperFunctions.hpp"
#include "ClpCholeskyDense.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpQuadraticObjective.hpp"
#include <cfloat>

#include <string>
#include <stdio.h>
#include <iostream>
//#############################################################################

ClpInterior::ClpInterior () :

  ClpModel(),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  worstComplementarity_(0.0),
  xsize_(0.0),
  zsize_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rhs_(NULL),
   x_(NULL),
   y_(NULL),
  dj_(NULL),
  lsqrObject_(NULL),
  pdcoStuff_(NULL),
  mu_(0.0),
  objectiveNorm_(1.0e-12),
  rhsNorm_(1.0e-12),
  solutionNorm_(1.0e-12),
  dualObjective_(0.0),
  primalObjective_(0.0),
  diagonalNorm_(1.0e-12),
  stepLength_(0.995),
  linearPerturbation_(1.0e-12),
  diagonalPerturbation_(1.0e-15),
  gamma_(0.0),
  delta_(0),
  targetGap_(1.0e-12),
  projectionTolerance_(1.0e-7),
  maximumRHSError_(0.0),
  maximumBoundInfeasibility_(0.0),
  maximumDualError_(0.0),
  diagonalScaleFactor_(0.0),
  scaleFactor_(1.0),
  actualPrimalStep_(0.0),
  actualDualStep_(0.0),
  smallestInfeasibility_(0.0),
  complementarityGap_(0.0),
  baseObjectiveNorm_(0.0),
  worstDirectionAccuracy_(0.0),
  maximumRHSChange_(0.0),
  errorRegion_(NULL),
  rhsFixRegion_(NULL),
  upperSlack_(NULL),
  lowerSlack_(NULL),
  diagonal_(NULL),
  solution_(NULL),
  workArray_(NULL),
  deltaX_(NULL),
  deltaY_(NULL),
  deltaZ_(NULL),
  deltaT_(NULL),
  deltaSU_(NULL),
  deltaSL_(NULL),
  primalR_(NULL),
  dualR_(NULL),
  rhsB_(NULL),
  rhsU_(NULL),
  rhsL_(NULL),
  rhsZ_(NULL),
  rhsT_(NULL),
  rhsC_(NULL),
  zVec_(NULL),
  tVec_(NULL),
  cholesky_(NULL),
  numberComplementarityPairs_(0),
  numberComplementarityItems_(0),
  maximumBarrierIterations_(200),
  gonePrimalFeasible_(false),
  goneDualFeasible_(false),
  algorithm_(-1)
{
  memset(historyInfeasibility_,0,LENGTH_HISTORY*sizeof(double));
  solveType_=3; // say interior based life form
  cholesky_ = new ClpCholeskyDense(); // put in placeholder
}

// Subproblem constructor
ClpInterior::ClpInterior ( const ClpModel * rhs,
			   int numberRows, const int * whichRow,
			   int numberColumns, const int * whichColumn,
			   bool dropNames, bool dropIntegers)
  : ClpModel(rhs, numberRows, whichRow,
	     numberColumns,whichColumn,dropNames,dropIntegers),
    largestPrimalError_(0.0),
    largestDualError_(0.0),
    sumDualInfeasibilities_(0.0),
    sumPrimalInfeasibilities_(0.0),
    worstComplementarity_(0.0),
    xsize_(0.0),
    zsize_(0.0),
    lower_(NULL),
    rowLowerWork_(NULL),
    columnLowerWork_(NULL),
    upper_(NULL),
    rowUpperWork_(NULL),
    columnUpperWork_(NULL),
    cost_(NULL),
    rhs_(NULL),
    x_(NULL),
    y_(NULL),
    dj_(NULL),
    lsqrObject_(NULL),
    pdcoStuff_(NULL),
    mu_(0.0),
    objectiveNorm_(1.0e-12),
    rhsNorm_(1.0e-12),
    solutionNorm_(1.0e-12),
    dualObjective_(0.0),
    primalObjective_(0.0),
    diagonalNorm_(1.0e-12),
    stepLength_(0.99995),
    linearPerturbation_(1.0e-12),
    diagonalPerturbation_(1.0e-15),
    gamma_(0.0),
    delta_(0),
    targetGap_(1.0e-12),
    projectionTolerance_(1.0e-7),
    maximumRHSError_(0.0),
    maximumBoundInfeasibility_(0.0),
    maximumDualError_(0.0),
    diagonalScaleFactor_(0.0),
    scaleFactor_(0.0),
    actualPrimalStep_(0.0),
    actualDualStep_(0.0),
    smallestInfeasibility_(0.0),
    complementarityGap_(0.0),
    baseObjectiveNorm_(0.0),
    worstDirectionAccuracy_(0.0),
    maximumRHSChange_(0.0),
    errorRegion_(NULL),
    rhsFixRegion_(NULL),
    upperSlack_(NULL),
    lowerSlack_(NULL),
    diagonal_(NULL),
    solution_(NULL),
    workArray_(NULL),
    deltaX_(NULL),
    deltaY_(NULL),
    deltaZ_(NULL),
    deltaT_(NULL),
    deltaSU_(NULL),
    deltaSL_(NULL),
    primalR_(NULL),
    dualR_(NULL),
    rhsB_(NULL),
    rhsU_(NULL),
    rhsL_(NULL),
    rhsZ_(NULL),
    rhsT_(NULL),
    rhsC_(NULL),
    zVec_(NULL),
    tVec_(NULL),
    cholesky_(NULL),
    numberComplementarityPairs_(0),
    numberComplementarityItems_(0),
    maximumBarrierIterations_(200),
    gonePrimalFeasible_(false),
    goneDualFeasible_(false),
    algorithm_(-1)
{
  memset(historyInfeasibility_,0,LENGTH_HISTORY*sizeof(double));
  solveType_=3; // say interior based life form
  cholesky_= new ClpCholeskyDense();
}

//-----------------------------------------------------------------------------

ClpInterior::~ClpInterior ()
{
  gutsOfDelete();
}
//#############################################################################
/* 
   This does housekeeping
*/
int 
ClpInterior::housekeeping()
{
  numberIterations_++;
  return 0;
}
// Copy constructor. 
ClpInterior::ClpInterior(const ClpInterior &rhs) :
  ClpModel(rhs),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  worstComplementarity_(0.0),
  xsize_(0.0),
  zsize_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rhs_(NULL),
   x_(NULL),
   y_(NULL),
  dj_(NULL),
  lsqrObject_(NULL),
  pdcoStuff_(NULL),
  errorRegion_(NULL),
  rhsFixRegion_(NULL),
  upperSlack_(NULL),
  lowerSlack_(NULL),
  diagonal_(NULL),
  solution_(NULL),
  workArray_(NULL),
  deltaX_(NULL),
  deltaY_(NULL),
  deltaZ_(NULL),
  deltaT_(NULL),
  deltaSU_(NULL),
  deltaSL_(NULL),
  primalR_(NULL),
  dualR_(NULL),
  rhsB_(NULL),
  rhsU_(NULL),
  rhsL_(NULL),
  rhsZ_(NULL),
  rhsT_(NULL),
  rhsC_(NULL),
  zVec_(NULL),
  tVec_(NULL),
  cholesky_(NULL)
{
  gutsOfDelete();
  gutsOfCopy(rhs);
  solveType_=3; // say interior based life form
}
// Copy constructor from model
ClpInterior::ClpInterior(const ClpModel &rhs) :
  ClpModel(rhs),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  sumDualInfeasibilities_(0.0),
  sumPrimalInfeasibilities_(0.0),
  worstComplementarity_(0.0),
  xsize_(0.0),
  zsize_(0.0),
  lower_(NULL),
  rowLowerWork_(NULL),
  columnLowerWork_(NULL),
  upper_(NULL),
  rowUpperWork_(NULL),
  columnUpperWork_(NULL),
  cost_(NULL),
  rhs_(NULL),
   x_(NULL),
   y_(NULL),
  dj_(NULL),
  lsqrObject_(NULL),
  pdcoStuff_(NULL),
  mu_(0.0),
  objectiveNorm_(1.0e-12),
  rhsNorm_(1.0e-12),
  solutionNorm_(1.0e-12),
  dualObjective_(0.0),
  primalObjective_(0.0),
  diagonalNorm_(1.0e-12),
  stepLength_(0.99995),
  linearPerturbation_(1.0e-12),
  diagonalPerturbation_(1.0e-15),
  gamma_(0.0),
  delta_(0),
  targetGap_(1.0e-12),
  projectionTolerance_(1.0e-7),
  maximumRHSError_(0.0),
  maximumBoundInfeasibility_(0.0),
  maximumDualError_(0.0),
  diagonalScaleFactor_(0.0),
  scaleFactor_(0.0),
  actualPrimalStep_(0.0),
  actualDualStep_(0.0),
  smallestInfeasibility_(0.0),
  complementarityGap_(0.0),
  baseObjectiveNorm_(0.0),
  worstDirectionAccuracy_(0.0),
  maximumRHSChange_(0.0),
  errorRegion_(NULL),
  rhsFixRegion_(NULL),
  upperSlack_(NULL),
  lowerSlack_(NULL),
  diagonal_(NULL),
  solution_(NULL),
  workArray_(NULL),
  deltaX_(NULL),
  deltaY_(NULL),
  deltaZ_(NULL),
  deltaT_(NULL),
  deltaSU_(NULL),
  deltaSL_(NULL),
  primalR_(NULL),
  dualR_(NULL),
  rhsB_(NULL),
  rhsU_(NULL),
  rhsL_(NULL),
  rhsZ_(NULL),
  rhsT_(NULL),
  rhsC_(NULL),
  zVec_(NULL),
  tVec_(NULL),
  cholesky_(NULL),
  numberComplementarityPairs_(0),
  numberComplementarityItems_(0),
  maximumBarrierIterations_(200),
  gonePrimalFeasible_(false),
  goneDualFeasible_(false),
  algorithm_(-1)
{
  memset(historyInfeasibility_,0,LENGTH_HISTORY*sizeof(double));
  solveType_=3; // say interior based life form
  cholesky_= new ClpCholeskyDense();
}
// Assignment operator. This copies the data
ClpInterior & 
ClpInterior::operator=(const ClpInterior & rhs)
{
  if (this != &rhs) {
    gutsOfDelete();
    ClpModel::operator=(rhs);
    gutsOfCopy(rhs);
  }
  return *this;
}
void 
ClpInterior::gutsOfCopy(const ClpInterior & rhs)
{
  lower_ = ClpCopyOfArray(rhs.lower_,numberColumns_+numberRows_);
  rowLowerWork_ = lower_+numberColumns_;
  columnLowerWork_ = lower_;
  upper_ = ClpCopyOfArray(rhs.upper_,numberColumns_+numberRows_);
  rowUpperWork_ = upper_+numberColumns_;
  columnUpperWork_ = upper_;
  //cost_ = ClpCopyOfArray(rhs.cost_,2*(numberColumns_+numberRows_));
  cost_ = ClpCopyOfArray(rhs.cost_,numberColumns_);
  rhs_ = ClpCopyOfArray(rhs.rhs_,numberRows_);
   x_ = ClpCopyOfArray(rhs.x_,numberColumns_);
   y_ = ClpCopyOfArray(rhs.y_,numberRows_);
  dj_ = ClpCopyOfArray(rhs.dj_,numberColumns_+numberRows_);
#ifdef PDCO
  lsqrObject_= new ClpLsqr(*rhs.lsqrObject_);
  pdcoStuff_ = rhs.pdcoStuff_->clone();
#endif
  largestPrimalError_ = rhs.largestPrimalError_;
  largestDualError_ = rhs.largestDualError_;
  sumDualInfeasibilities_ = rhs.sumDualInfeasibilities_;
  sumPrimalInfeasibilities_ = rhs.sumPrimalInfeasibilities_;
  worstComplementarity_ = rhs.worstComplementarity_;
  xsize_ = rhs.xsize_;
  zsize_ = rhs.zsize_;
  solveType_=rhs.solveType_;
  mu_ = rhs.mu_;
  objectiveNorm_ = rhs.objectiveNorm_;
  rhsNorm_ = rhs.rhsNorm_;
  solutionNorm_ = rhs.solutionNorm_;
  dualObjective_ = rhs.dualObjective_;
  primalObjective_ = rhs.primalObjective_;
  diagonalNorm_ = rhs.diagonalNorm_;
  stepLength_ = rhs.stepLength_;
  linearPerturbation_ = rhs.linearPerturbation_;
  diagonalPerturbation_ = rhs.diagonalPerturbation_;
  gamma_=rhs.gamma_;
  delta_=rhs.delta_;
  targetGap_ = rhs.targetGap_;
  projectionTolerance_ = rhs.projectionTolerance_;
  maximumRHSError_ = rhs.maximumRHSError_;
  maximumBoundInfeasibility_ = rhs.maximumBoundInfeasibility_;
  maximumDualError_ = rhs.maximumDualError_;
  diagonalScaleFactor_ = rhs.diagonalScaleFactor_;
  scaleFactor_ = rhs.scaleFactor_;
  actualPrimalStep_ = rhs.actualPrimalStep_;
  actualDualStep_ = rhs.actualDualStep_;
  smallestInfeasibility_ = rhs.smallestInfeasibility_;
  complementarityGap_ = rhs.complementarityGap_;
  baseObjectiveNorm_ = rhs.baseObjectiveNorm_;
  worstDirectionAccuracy_ = rhs.worstDirectionAccuracy_;
  maximumRHSChange_ = rhs.maximumRHSChange_;
  errorRegion_ = ClpCopyOfArray(rhs.errorRegion_,numberRows_);
  rhsFixRegion_ = ClpCopyOfArray(rhs.rhsFixRegion_,numberRows_);
  deltaY_ = ClpCopyOfArray(rhs.deltaY_,numberRows_);
  upperSlack_ = ClpCopyOfArray(rhs.upperSlack_,numberRows_+numberColumns_);
  lowerSlack_ = ClpCopyOfArray(rhs.lowerSlack_,numberRows_+numberColumns_);
  diagonal_ = ClpCopyOfArray(rhs.diagonal_,numberRows_+numberColumns_);
  deltaX_ = ClpCopyOfArray(rhs.deltaX_,numberRows_+numberColumns_);
  deltaZ_ = ClpCopyOfArray(rhs.deltaZ_,numberRows_+numberColumns_);
  deltaT_ = ClpCopyOfArray(rhs.deltaT_,numberRows_+numberColumns_);
  deltaSU_ = ClpCopyOfArray(rhs.deltaSU_,numberRows_+numberColumns_);
  deltaSL_ = ClpCopyOfArray(rhs.deltaSL_,numberRows_+numberColumns_);
  primalR_ = ClpCopyOfArray(rhs.primalR_,numberRows_+numberColumns_);
  dualR_ = ClpCopyOfArray(rhs.dualR_,numberRows_+numberColumns_);
  rhsB_ = ClpCopyOfArray(rhs.rhsB_,numberRows_);
  rhsU_ = ClpCopyOfArray(rhs.rhsU_,numberRows_+numberColumns_);
  rhsL_ = ClpCopyOfArray(rhs.rhsL_,numberRows_+numberColumns_);
  rhsZ_ = ClpCopyOfArray(rhs.rhsZ_,numberRows_+numberColumns_);
  rhsT_ = ClpCopyOfArray(rhs.rhsT_,numberRows_+numberColumns_);
  rhsC_ = ClpCopyOfArray(rhs.rhsC_,numberRows_+numberColumns_);
  solution_ = ClpCopyOfArray(rhs.solution_,numberRows_+numberColumns_);
  workArray_ = ClpCopyOfArray(rhs.workArray_,numberRows_+numberColumns_);
  zVec_ = ClpCopyOfArray(rhs.zVec_,numberRows_+numberColumns_);
  tVec_ = ClpCopyOfArray(rhs.tVec_,numberRows_+numberColumns_);
  cholesky_ = rhs.cholesky_->clone();
  numberComplementarityPairs_ = rhs.numberComplementarityPairs_;
  numberComplementarityItems_ = rhs.numberComplementarityItems_;
  maximumBarrierIterations_ = rhs.maximumBarrierIterations_;
  gonePrimalFeasible_ = rhs.gonePrimalFeasible_;
  goneDualFeasible_ = rhs.goneDualFeasible_;
  algorithm_ = rhs.algorithm_;
}

void 
ClpInterior::gutsOfDelete()
{
  delete [] lower_;
  lower_=NULL;
  rowLowerWork_=NULL;
  columnLowerWork_=NULL;
  delete [] upper_;
  upper_=NULL;
  rowUpperWork_=NULL;
  columnUpperWork_=NULL;
  delete [] cost_;
  cost_=NULL;
  delete [] rhs_;
  rhs_ = NULL;
  delete [] x_;
  x_ = NULL;
  delete [] y_;
  y_ = NULL;
  delete [] dj_;
  dj_ = NULL;
#ifdef PDCO
  delete lsqrObject_;
  lsqrObject_ = NULL;
  //delete pdcoStuff_;
  pdcoStuff_=NULL;
#endif
  delete [] errorRegion_;
  errorRegion_ = NULL;
  delete [] rhsFixRegion_;
  rhsFixRegion_ = NULL;
  delete [] deltaY_;
  deltaY_ = NULL;
  delete [] upperSlack_;
  upperSlack_ = NULL;
  delete [] lowerSlack_;
  lowerSlack_ = NULL;
  delete [] diagonal_;
  diagonal_ = NULL;
  delete [] deltaX_;
  deltaX_ = NULL;
  delete [] deltaZ_;
  deltaZ_ = NULL;
  delete [] deltaT_;
  deltaT_ = NULL;
  delete [] deltaSU_;
  deltaSU_ = NULL;
  delete [] deltaSL_;
  deltaSL_ = NULL;
  delete [] primalR_;
  primalR_=NULL;
  delete [] dualR_;
  dualR_=NULL;
  delete [] rhsB_;
  rhsB_ = NULL;
  delete [] rhsU_;
  rhsU_ = NULL;
  delete [] rhsL_;
  rhsL_ = NULL;
  delete [] rhsZ_;
  rhsZ_ = NULL;
  delete [] rhsT_;
  rhsT_ = NULL;
  delete [] rhsC_;
  rhsC_ = NULL;
  delete [] solution_;
  solution_ = NULL;
  delete [] workArray_;
  workArray_ = NULL;
  delete [] zVec_;
  zVec_ = NULL;
  delete [] tVec_;
  tVec_ = NULL;
  delete cholesky_;
}
bool
ClpInterior::createWorkingData()
{
  bool goodMatrix=true;
  //check matrix
  if (!matrix_->allElementsInRange(this,1.0e-12,1.0e20)) {
    problemStatus_=4;
    goodMatrix= false;
  }
  int nTotal = numberRows_+numberColumns_;
  delete [] solution_;
  solution_ = new double[nTotal];
  memcpy(solution_,columnActivity_,
	 numberColumns_*sizeof(double));
  memcpy(solution_+numberColumns_,rowActivity_,
	 numberRows_*sizeof(double));
  delete [] cost_;
  cost_ = new double[nTotal];
  int i;
  double direction = optimizationDirection_*objectiveScale_;
  // direction is actually scale out not scale in
  if (direction)
    direction = 1.0/direction;
  const double * obj = objective();
  for (i=0;i<numberColumns_;i++)
    cost_[i] = direction*obj[i];
  memset(cost_+numberColumns_,0,numberRows_*sizeof(double));
  // do scaling if needed
  if (scalingFlag_>0&&!rowScale_) {
    if (matrix_->scale(this))
      scalingFlag_=-scalingFlag_; // not scaled after all
  }
  delete [] lower_;
  delete [] upper_;
  lower_ = new double[nTotal];
  upper_ = new double[nTotal];
  rowLowerWork_ = lower_+numberColumns_;
  columnLowerWork_ = lower_;
  rowUpperWork_ = upper_+numberColumns_;
  columnUpperWork_ = upper_;
  memcpy(rowLowerWork_,rowLower_,numberRows_*sizeof(double));
  memcpy(rowUpperWork_,rowUpper_,numberRows_*sizeof(double));
  memcpy(columnLowerWork_,columnLower_,numberColumns_*sizeof(double));
  memcpy(columnUpperWork_,columnUpper_,numberColumns_*sizeof(double));
  // clean up any mismatches on infinity
  for (i=0;i<numberColumns_;i++) {
    if (columnLowerWork_[i]<-1.0e30)
      columnLowerWork_[i] = -COIN_DBL_MAX;
    if (columnUpperWork_[i]>1.0e30)
      columnUpperWork_[i] = COIN_DBL_MAX;
  }
  // clean up any mismatches on infinity
  for (i=0;i<numberRows_;i++) {
    if (rowLowerWork_[i]<-1.0e30)
      rowLowerWork_[i] = -COIN_DBL_MAX;
    if (rowUpperWork_[i]>1.0e30)
      rowUpperWork_[i] = COIN_DBL_MAX;
  }
  // check rim of problem okay
  if (!sanityCheck())
    goodMatrix=false;
  if(rowScale_) {
    for (i=0;i<numberColumns_;i++) {
      double multiplier = rhsScale_/columnScale_[i];
      cost_[i] *= columnScale_[i];
      if (columnLowerWork_[i]>-1.0e50)
	columnLowerWork_[i] *= multiplier;
      if (columnUpperWork_[i]<1.0e50)
	columnUpperWork_[i] *= multiplier;
      
    }
    for (i=0;i<numberRows_;i++) {
      double multiplier = rhsScale_*rowScale_[i];
      if (rowLowerWork_[i]>-1.0e50)
	rowLowerWork_[i] *= multiplier;
      if (rowUpperWork_[i]<1.0e50)
	rowUpperWork_[i] *= multiplier;
    }
  } else if (rhsScale_!=1.0) {
    for (i=0;i<numberColumns_+numberRows_;i++) {
      if (lower_[i]>-1.0e50)
	lower_[i] *= rhsScale_;
      if (upper_[i]<1.0e50)
	upper_[i] *= rhsScale_;
      
    }
  }
  assert (!errorRegion_);
  errorRegion_ = new double [numberRows_];
  assert (!rhsFixRegion_);
  rhsFixRegion_ = new double [numberRows_];
  assert (!deltaY_);
  deltaY_ = new double [numberRows_];
  CoinZeroN(deltaY_,numberRows_);
  assert (!upperSlack_);
  upperSlack_ = new double [nTotal];
  assert (!lowerSlack_);
  lowerSlack_ = new double [nTotal];
  assert (!diagonal_);
  diagonal_ = new double [nTotal];
  assert (!deltaX_);
  deltaX_ = new double [nTotal];
  CoinZeroN(deltaX_,nTotal);
  assert (!deltaZ_);
  deltaZ_ = new double [nTotal];
  CoinZeroN(deltaZ_,nTotal);
  assert (!deltaT_);
  deltaT_ = new double [nTotal];
  CoinZeroN(deltaT_,nTotal);
  assert (!deltaSU_);
  deltaSU_ = new double [nTotal];
  CoinZeroN(deltaSU_,nTotal);
  assert (!deltaSL_);
  deltaSL_ = new double [nTotal];
  CoinZeroN(deltaSL_,nTotal);
  assert (!primalR_);
  assert (!dualR_);
  // create arrays if we are doing KKT
  if (cholesky_->type()>20) {
    primalR_ = new double [nTotal];
    CoinZeroN(primalR_,nTotal);
    dualR_ = new double [numberRows_];
    CoinZeroN(dualR_,numberRows_);
  }
  assert (!rhsB_);
  rhsB_ = new double [numberRows_];
  CoinZeroN(rhsB_,numberRows_);
  assert (!rhsU_);
  rhsU_ = new double [nTotal];
  CoinZeroN(rhsU_,nTotal);
  assert (!rhsL_);
  rhsL_ = new double [nTotal];
  CoinZeroN(rhsL_,nTotal);
  assert (!rhsZ_);
  rhsZ_ = new double [nTotal];
  CoinZeroN(rhsZ_,nTotal);
  assert (!rhsT_);
  rhsT_ = new double [nTotal];
  CoinZeroN(rhsT_,nTotal);
  assert (!rhsC_);
  rhsC_ = new double [nTotal];
  CoinZeroN(rhsC_,nTotal);
  assert (!workArray_);
  workArray_ = new double [nTotal];
  CoinZeroN(workArray_,nTotal);
  assert (!zVec_);
  zVec_ = new double [nTotal];
  CoinZeroN(zVec_,nTotal);
  assert (!tVec_);
  tVec_ = new double [nTotal];
  CoinZeroN(tVec_,nTotal);
  assert (!dj_);
  dj_ = new double [nTotal];
  if (!status_)
    status_ = new unsigned char [numberRows_+numberColumns_];
  memset(status_,0,numberRows_+numberColumns_);
  return goodMatrix;
}
void
ClpInterior::deleteWorkingData()
{
  int i;
  if (optimizationDirection_!=1.0||objectiveScale_!=1.0) {
    double scaleC = optimizationDirection_/objectiveScale_;
    // and modify all dual signs
    for (i=0;i<numberColumns_;i++) 
      reducedCost_[i] = scaleC*dj_[i];
    for (i=0;i<numberRows_;i++) 
      dual_[i] *= scaleC;
  }
  if (rowScale_) {
    double scaleR = 1.0/rhsScale_;
    for (i=0;i<numberColumns_;i++) {
      double scaleFactor = columnScale_[i];
      double valueScaled = columnActivity_[i];
      columnActivity_[i] = valueScaled*scaleFactor*scaleR;
      double valueScaledDual = reducedCost_[i];
      reducedCost_[i] = valueScaledDual/scaleFactor;
    }
    for (i=0;i<numberRows_;i++) {
      double scaleFactor = rowScale_[i];
      double valueScaled = rowActivity_[i];
      rowActivity_[i] = (valueScaled*scaleR)/scaleFactor;
      double valueScaledDual = dual_[i];
      dual_[i] = valueScaledDual*scaleFactor;
    }
  } else if (rhsScale_!=1.0) {
    double scaleR = 1.0/rhsScale_;
    for (i=0;i<numberColumns_;i++) {
      double valueScaled = columnActivity_[i];
      columnActivity_[i] = valueScaled*scaleR;
    }
    for (i=0;i<numberRows_;i++) {
      double valueScaled = rowActivity_[i];
      rowActivity_[i] = valueScaled*scaleR;
    }
  }
  delete [] cost_;
  cost_ = NULL;
  delete [] solution_;
  solution_ = NULL;
  delete [] lower_;
  lower_ = NULL;
  delete [] upper_;
  upper_ = NULL;
  delete [] errorRegion_;
  errorRegion_ = NULL;
  delete [] rhsFixRegion_;
  rhsFixRegion_ = NULL;
  delete [] deltaY_;
  deltaY_ = NULL;
  delete [] upperSlack_;
  upperSlack_ = NULL;
  delete [] lowerSlack_;
  lowerSlack_ = NULL;
  delete [] diagonal_;
  diagonal_ = NULL;
  delete [] deltaX_;
  deltaX_ = NULL;
  delete [] workArray_;
  workArray_ = NULL;
  delete [] zVec_;
  zVec_ = NULL;
  delete [] tVec_;
  tVec_ = NULL;
  delete [] dj_;
  dj_ = NULL;
}
// Sanity check on input data - returns true if okay
bool 
ClpInterior::sanityCheck()
{
  // bad if empty
  if (!numberRows_||!numberColumns_||!matrix_->getNumElements()) {
    problemStatus_=emptyProblem();
    return false;
  }
  int numberBad ;
  double largestBound, smallestBound, minimumGap;
  double smallestObj, largestObj;
  int firstBad;
  int modifiedBounds=0;
  int i;
  numberBad=0;
  firstBad=-1;
  minimumGap=1.0e100;
  smallestBound=1.0e100;
  largestBound=0.0;
  smallestObj=1.0e100;
  largestObj=0.0;
  // If bounds are too close - fix
  double fixTolerance = 1.1*primalTolerance();
  for (i=numberColumns_;i<numberColumns_+numberRows_;i++) {
    double value;
    value = fabs(cost_[i]);
    if (value>1.0e50) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value) {
      if (value>largestObj)
	largestObj=value;
      if (value<smallestObj)
	smallestObj=value;
    }
    value=upper_[i]-lower_[i];
    if (value<-primalTolerance()) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value<=fixTolerance) {
      if (value) {
	// modify
	upper_[i] = lower_[i];
	modifiedBounds++;
      }
    } else {
      if (value<minimumGap)
	minimumGap=value;
    }
    if (lower_[i]>-1.0e100&&lower_[i]) {
      value = fabs(lower_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
    if (upper_[i]<1.0e100&&upper_[i]) {
      value = fabs(upper_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
  }
  if (largestBound)
    handler_->message(CLP_RIMSTATISTICS3,messages_)
      <<smallestBound
      <<largestBound
      <<minimumGap
      <<CoinMessageEol;
  minimumGap=1.0e100;
  smallestBound=1.0e100;
  largestBound=0.0;
  for (i=0;i<numberColumns_;i++) {
    double value;
    value = fabs(cost_[i]);
    if (value>1.0e50) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value) {
      if (value>largestObj)
	largestObj=value;
      if (value<smallestObj)
	smallestObj=value;
    }
    value=upper_[i]-lower_[i];
    if (value<-primalTolerance()) {
      numberBad++;
      if (firstBad<0)
	firstBad=i;
    } else if (value<=fixTolerance) {
      if (value) {
	// modify
	upper_[i] = lower_[i];
	modifiedBounds++;
      }
    } else {
      if (value<minimumGap)
	minimumGap=value;
    }
    if (lower_[i]>-1.0e100&&lower_[i]) {
      value = fabs(lower_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
    if (upper_[i]<1.0e100&&upper_[i]) {
      value = fabs(upper_[i]);
      if (value>largestBound)
	largestBound=value;
      if (value<smallestBound)
	smallestBound=value;
    }
  }
  char rowcol[]={'R','C'};
  if (numberBad) {
    handler_->message(CLP_BAD_BOUNDS,messages_)
      <<numberBad
      <<rowcol[isColumn(firstBad)]<<sequenceWithin(firstBad)
      <<CoinMessageEol;
    problemStatus_=4;
    return false;
  }
  if (modifiedBounds)
    handler_->message(CLP_MODIFIEDBOUNDS,messages_)
      <<modifiedBounds
      <<CoinMessageEol;
  handler_->message(CLP_RIMSTATISTICS1,messages_)
    <<smallestObj
    <<largestObj
    <<CoinMessageEol;  if (largestBound)
    handler_->message(CLP_RIMSTATISTICS2,messages_)
      <<smallestBound
      <<largestBound
      <<minimumGap
      <<CoinMessageEol;
  return true;
}
/* Loads a problem (the constraints on the
   rows are given by lower and upper bounds). If a pointer is 0 then the
   following values are the default:
   <ul>
   <li> <code>colub</code>: all columns have upper bound infinity
   <li> <code>collb</code>: all columns have lower bound 0 
   <li> <code>rowub</code>: all rows have upper bound infinity
   <li> <code>rowlb</code>: all rows have lower bound -infinity
   <li> <code>obj</code>: all variables have 0 objective coefficient
   </ul>
*/
void 
ClpInterior::loadProblem (  const ClpMatrixBase& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
}
void 
ClpInterior::loadProblem (  const CoinPackedMatrix& matrix,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(matrix, collb, colub, obj, rowlb, rowub,
			rowObjective);
}

/* Just like the other loadProblem() method except that the matrix is
   given in a standard column major ordered format (without gaps). */
void 
ClpInterior::loadProblem (  const int numcols, const int numrows,
		    const CoinBigIndex* start, const int* index,
		    const double* value,
		    const double* collb, const double* colub,   
		    const double* obj,
		    const double* rowlb, const double* rowub,
		    const double * rowObjective)
{
  ClpModel::loadProblem(numcols, numrows, start, index, value,
			  collb, colub, obj, rowlb, rowub,
			  rowObjective);
}
void 
ClpInterior::loadProblem (  const int numcols, const int numrows,
			   const CoinBigIndex* start, const int* index,
			   const double* value,const int * length,
			   const double* collb, const double* colub,   
			   const double* obj,
			   const double* rowlb, const double* rowub,
			   const double * rowObjective)
{
  ClpModel::loadProblem(numcols, numrows, start, index, value, length,
			  collb, colub, obj, rowlb, rowub,
			  rowObjective);
}
// Read an mps file from the given filename
int 
ClpInterior::readMps(const char *filename,
	    bool keepNames,
	    bool ignoreErrors)
{
  int status = ClpModel::readMps(filename,keepNames,ignoreErrors);
  return status;
}
#ifdef PDCO
#include "ClpPdco.hpp"
/* Pdco algorithm - see ClpPdco.hpp for method */
int 
ClpInterior::pdco()
{
  return ((ClpPdco *) this)->pdco();
}
// ** Temporary version
int  
ClpInterior::pdco( ClpPdcoBase * stuff, Options &options, Info &info, Outfo &outfo)
{
  return ((ClpPdco *) this)->pdco(stuff,options,info,outfo);
}
#endif
#include "ClpPredictorCorrector.hpp"
// Primal-Dual Predictor-Corrector barrier
int 
ClpInterior::primalDual()
{ 
  return ((ClpPredictorCorrector *) this)->solve();
}

void 
ClpInterior::checkSolution()
{
  int iRow,iColumn;
  CoinMemcpyN(cost_,numberColumns_,reducedCost_);
  matrix_->transposeTimes(-1.0,dual_,reducedCost_);
  // Now modify reduced costs for quadratic
  double quadraticOffset=quadraticDjs(reducedCost_,
				      solution_,scaleFactor_);

  objectiveValue_ = 0.0;
  // now look at solution
  sumPrimalInfeasibilities_=0.0;
  sumDualInfeasibilities_=0.0;
  double dualTolerance =  10.0*dblParam_[ClpDualTolerance];
  double primalTolerance =  dblParam_[ClpPrimalTolerance];
  double primalTolerance2 =  10.0*dblParam_[ClpPrimalTolerance];
  worstComplementarity_=0.0;
  complementarityGap_=0.0;

  // Done scaled - use permanent regions for output
  // but internal for bounds
  const double * lower = lower_+numberColumns_;
  const double * upper = upper_+numberColumns_;
  for (iRow=0;iRow<numberRows_;iRow++) {
    double infeasibility=0.0;
    double distanceUp = min(upper[iRow]-
      rowActivity_[iRow],1.0e10);
    double distanceDown = min(rowActivity_[iRow] -
      lower[iRow],1.0e10);
    if (distanceUp>primalTolerance2) {
      double value = dual_[iRow];
      // should not be negative
      if (value<-dualTolerance) {
	sumDualInfeasibilities_ += -dualTolerance-value;
	value = - value*distanceUp;
	if (value>worstComplementarity_) 
	  worstComplementarity_=value;
	complementarityGap_ += value;
      }
    }
    if (distanceDown>primalTolerance2) {
      double value = dual_[iRow];
      // should not be positive
      if (value>dualTolerance) {
	sumDualInfeasibilities_ += value-dualTolerance;
	value =  value*distanceDown;
	if (value>worstComplementarity_) 
	  worstComplementarity_=value;
	complementarityGap_ += value;
      }
    }
    if (rowActivity_[iRow]>upper[iRow]) {
      infeasibility=rowActivity_[iRow]-upper[iRow];
    } else if (rowActivity_[iRow]<lower[iRow]) {
      infeasibility=lower[iRow]-rowActivity_[iRow];
    }
    if (infeasibility>primalTolerance) {
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance;
    }
  }
  lower = lower_;
  upper = upper_;
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double infeasibility=0.0;
    objectiveValue_ += cost_[iColumn]*columnActivity_[iColumn];
    double distanceUp = min(upper[iColumn]-
      columnActivity_[iColumn],1.0e10);
    double distanceDown = min(columnActivity_[iColumn] -
      lower[iColumn],1.0e10);
    if (distanceUp>primalTolerance2) {
      double value = reducedCost_[iColumn];
      // should not be negative
      if (value<-dualTolerance) {
	sumDualInfeasibilities_ += -dualTolerance-value;
	value = - value*distanceUp;
	if (value>worstComplementarity_) 
	  worstComplementarity_=value;
	complementarityGap_ += value;
      }
    }
    if (distanceDown>primalTolerance2) {
      double value = reducedCost_[iColumn];
      // should not be positive
      if (value>dualTolerance) {
	sumDualInfeasibilities_ += value-dualTolerance;
	value =  value*distanceDown;
	if (value>worstComplementarity_) 
	  worstComplementarity_=value;
	complementarityGap_ += value;
      }
    }
    if (columnActivity_[iColumn]>upper[iColumn]) {
      infeasibility=columnActivity_[iColumn]-upper[iColumn];
    } else if (columnActivity_[iColumn]<lower[iColumn]) {
      infeasibility=lower[iColumn]-columnActivity_[iColumn];
    }
    if (infeasibility>primalTolerance) {
      sumPrimalInfeasibilities_ += infeasibility-primalTolerance;
    }
  }
  // add in offset
  objectiveValue_ += 0.5*quadraticOffset;
}
// Set cholesky (and delete present one)
void 
ClpInterior::setCholesky(ClpCholeskyBase * cholesky)
{
  delete cholesky_;
  cholesky_= cholesky;
}
/* Borrow model.  This is so we dont have to copy large amounts
   of data around.  It assumes a derived class wants to overwrite
   an empty model with a real one - while it does an algorithm.
   This is same as ClpModel one. */
void 
ClpInterior::borrowModel(ClpModel & otherModel) 
{
  ClpModel::borrowModel(otherModel);
}
/* Return model - updates any scalars */
void 
ClpInterior::returnModel(ClpModel & otherModel)
{
  ClpModel::returnModel(otherModel);
}
// Return number fixed to see if worth presolving
int 
ClpInterior::numberFixed() const
{
  int i;
  int nFixed=0;
  for (i=0;i<numberColumns_;i++) {
    if (columnUpper_[i]<1.0e20||columnLower_[i]>-1.0e20) { 
      if (columnUpper_[i]>columnLower_[i]) { 
	if (fixedOrFree(i))
	  nFixed++;
      }
    }
  }
  for (i=0;i<numberRows_;i++) {
    if (rowUpper_[i]<1.0e20||rowLower_[i]>-1.0e20) { 
      if (rowUpper_[i]>rowLower_[i]) { 
	if (fixedOrFree(i+numberColumns_))
	  nFixed++;
      }
    }
  }
  return nFixed;
}
// fix variables interior says should be
void 
ClpInterior::fixFixed(bool reallyFix)
{
  int i;
  for (i=0;i<numberColumns_;i++) {
    if (columnUpper_[i]<1.0e20||columnLower_[i]>-1.0e20) { 
      if (columnUpper_[i]>columnLower_[i]) { 
	if (fixedOrFree(i)) {
	  if (columnActivity_[i]-columnLower_[i]<columnUpper_[i]-columnActivity_[i]) {
	    if (reallyFix)
	      columnUpper_[i]=columnLower_[i];
	    columnActivity_[i]=columnLower_[i];
	  } else {
	    if (reallyFix)
	      columnLower_[i]=columnUpper_[i];
	    columnActivity_[i]=columnUpper_[i];
	  }
	}
      }
    }
  }
  for (i=0;i<numberRows_;i++) {
    if (rowUpper_[i]<1.0e20||rowLower_[i]>-1.0e20) { 
      if (rowUpper_[i]>rowLower_[i]) { 
	if (fixedOrFree(i+numberColumns_)) {
	  if (rowActivity_[i]-rowLower_[i]<rowUpper_[i]-rowActivity_[i]) {
	    rowUpper_[i]=rowLower_[i];
	    rowActivity_[i]=rowLower_[i];
	  } else {
	    rowLower_[i]=rowUpper_[i];
	    rowActivity_[i]=rowUpper_[i];
	  }
	}
      }
    }
  }
}
/* Modifies djs to allow for quadratic.
   returns quadratic offset */
double 
ClpInterior::quadraticDjs(double * djRegion, const double * solution,
			  double scaleFactor)
{
  double quadraticOffset=0.0;
#ifndef NO_RTTI
  ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(objective_));
#else
  ClpQuadraticObjective * quadraticObj = NULL;
  if (objective_->type()==2)
    quadraticObj = (static_cast< ClpQuadraticObjective*>(objective_));
#endif
  if (quadraticObj) {
    CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
    const int * columnQuadratic = quadratic->getIndices();
    const CoinBigIndex * columnQuadraticStart = quadratic->getVectorStarts();
    const int * columnQuadraticLength = quadratic->getVectorLengths();
    double * quadraticElement = quadratic->getMutableElements();
    int numberColumns = quadratic->getNumCols();
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      double value = 0.0;
      for (CoinBigIndex j=columnQuadraticStart[iColumn];
	   j<columnQuadraticStart[iColumn]+columnQuadraticLength[iColumn];j++) {
	int jColumn = columnQuadratic[j];
	double valueJ = solution[jColumn];
	double elementValue = quadraticElement[j];
	//value += valueI*valueJ*elementValue;
	value += valueJ*elementValue;
	quadraticOffset += solution[iColumn]*valueJ*elementValue;
      }
      djRegion[iColumn] += scaleFactor*value;
    }
  }
  return quadraticOffset;
}
