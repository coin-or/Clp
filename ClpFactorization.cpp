// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpFactorization.hpp"
#include "ClpQuadraticObjective.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"
#include "ClpNetworkBasis.hpp"
#include "ClpNetworkMatrix.hpp"
//#define CHECK_NETWORK
#ifdef CHECK_NETWORK
const static bool doCheck=true;
#else
const static bool doCheck=false;
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization () :
   CoinFactorization() 
{
  networkBasis_ = NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization (const ClpFactorization & rhs) :
   CoinFactorization(rhs) 
{
  if (rhs.networkBasis_)
    networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
  else
    networkBasis_=NULL;
}

ClpFactorization::ClpFactorization (const CoinFactorization & rhs) :
   CoinFactorization(rhs) 
{
  networkBasis_=NULL;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization () 
{
  delete networkBasis_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization& rhs)
{
  if (this != &rhs) {
    CoinFactorization::operator=(rhs);
    delete networkBasis_;
    if (rhs.networkBasis_)
      networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
    else
      networkBasis_=NULL;
  }
  return *this;
}
int 
ClpFactorization::factorize ( ClpSimplex * model,
			      int solveType, bool valuesPass)
{
  ClpMatrixBase * matrix = model->clpMatrix(); 
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  // If too many compressions increase area
  if (numberPivots_>1&&numberCompressions_*10 > numberPivots_+10) {
    areaFactor_ *= 1.1;
  }
  if (!networkBasis_||doCheck) {
    status_=-99;
    int * pivotVariable = model->pivotVariable();
    //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
    while (status_<-98) {
      
      int i;
      int numberBasic=0;
      int numberRowBasic;
      // Move pivot variables across if they look good
      int * pivotTemp = model->rowArray(0)->getIndices();
      assert (!model->rowArray(0)->getNumElements());
      if (!matrix->rhsOffset(model)) {
	// Seems to prefer things in order so quickest
	// way is to go though like this
	for (i=0;i<numberRows;i++) {
	  if (model->getRowStatus(i) == ClpSimplex::basic) 
	    pivotTemp[numberBasic++]=i;
	}
	numberRowBasic=numberBasic;
	/* Put column basic variables into pivotVariable
	   This is done by ClpMatrixBase to allow override for gub
	*/
	matrix->generalExpanded(model,0,numberBasic);
      } else {
	// Long matrix - do a different way
	bool fullSearch=false;
	for (i=0;i<numberRows;i++) {
	  int iPivot = pivotVariable[i];
	  if (iPivot>=numberColumns) {
	    pivotTemp[numberBasic++]=iPivot-numberColumns;
	  }
	}
	numberRowBasic=numberBasic;
	for (i=0;i<numberRows;i++) {
	  int iPivot = pivotVariable[i];
	  if (iPivot<numberColumns) {
	    if (iPivot>=0) {
	      pivotTemp[numberBasic++]=iPivot;
	    } else {
	      // not full basis
	      fullSearch=true;
	      break;
	    }
	  }
	}
	if (fullSearch) {
	  // do slow way
	  numberBasic=0;
	  for (i=0;i<numberRows;i++) {
	    if (model->getRowStatus(i) == ClpSimplex::basic) 
	      pivotTemp[numberBasic++]=i;
	  }
	  numberRowBasic=numberBasic;
	  /* Put column basic variables into pivotVariable
	     This is done by ClpMatrixBase to allow override for gub
	  */
	  matrix->generalExpanded(model,0,numberBasic);
	}
      }
      assert (numberBasic<=model->maximumBasic());
      // see if matrix a network
#ifndef NO_RTTI
      ClpNetworkMatrix* networkMatrix =
	dynamic_cast< ClpNetworkMatrix*>(model->clpMatrix());
#else
      ClpNetworkMatrix* networkMatrix = NULL;
      if (model->clpMatrix()->type()==11)
	networkMatrix = 
	static_cast< ClpNetworkMatrix*>(model->clpMatrix());
#endif
      // If network - still allow ordinary factorization first time for laziness
      int saveMaximumPivots = maximumPivots();
      delete networkBasis_;
      networkBasis_ = NULL;
      if (networkMatrix&&!doCheck)
	maximumPivots(1);
      while (status_==-99) {
	// maybe for speed will be better to leave as many regions as possible
	gutsOfDestructor();
	gutsOfInitialize(2);
	CoinBigIndex numberElements=numberRowBasic;

	// compute how much in basis

	int i;
	// can change for gub
	int numberColumnBasic = numberBasic-numberRowBasic;

	numberElements +=matrix->fillBasis(model,
					   pivotTemp+numberRowBasic, 
					   numberRowBasic,
					   numberColumnBasic,
					   NULL,NULL,NULL);
	// and recompute as network side say different
	if (model->numberIterations())
	  numberRowBasic = numberBasic - numberColumnBasic;
	numberElements = 3 * numberBasic + 3 * numberElements + 10000;
#if 0
	// If iteration not zero then may be compressed
	getAreas ( !model->numberIterations() ? numberRows : numberBasic, 
		   numberRowBasic+numberColumnBasic, numberElements,
		   2 * numberElements );
#else
	getAreas ( numberRows,
		   numberRowBasic+numberColumnBasic, numberElements,
		   2 * numberElements );
#endif
	//fill
	//copy
	numberElements=numberRowBasic;
	for (i=0;i<numberRowBasic;i++) {
	  int iRow = pivotTemp[i];
	  indexRowU_[i]=iRow;
	  indexColumnU_[i]=i;
	  elementU_[i]=slackValue_;
	}
	// can change for gub so redo
	numberColumnBasic = numberBasic-numberRowBasic;
	numberElements +=matrix->fillBasis(model, 
					   pivotTemp+numberRowBasic, 
					   numberRowBasic, 
					   numberColumnBasic,
					   indexRowU_+numberElements, 
					   indexColumnU_+numberElements,
					   elementU_+numberElements);
#if 0
	{
	  printf("%d row basic, %d column basic\n",numberRowBasic,numberColumnBasic);
	  for (int i=0;i<numberElements;i++) 
	    printf("row %d col %d value %g\n",indexRowU_[i],indexColumnU_[i],
		   elementU_[i]);
	}
#endif
	// recompute number basic
        numberBasic = numberRowBasic+numberColumnBasic;
	lengthU_ = numberElements;

	preProcess ( 0 );
	factor (  );
	if (status_==-99) {
	  // get more memory
	  areaFactor(2.0*areaFactor());
	}
      }
      // If we get here status is 0 or -1
      if (status_ == 0) {
	// We may need to tamper with order and redo - e.g. network with side
	int useNumberRows = numberRows;
	// **** we will also need to add test in dual steepest to do
	// as we do for network
        matrix->generalExpanded(model,12,useNumberRows);
	int * permuteBack = permuteBack_;
	int * back = pivotColumnBack_;
	//int * pivotTemp = pivotColumn_;
	//ClpDisjointCopyN ( pivotVariable, numberRows , pivotTemp  );
	// Redo pivot order
	for (i=0;i<numberRowBasic;i++) {
	  int k = pivotTemp[i];
	  // so rowIsBasic[k] would be permuteBack[back[i]]
	  pivotVariable[permuteBack[back[i]]]=k+numberColumns;
	}
	for (;i<useNumberRows;i++) {
	  int k = pivotTemp[i];
	  // so rowIsBasic[k] would be permuteBack[back[i]]
	  pivotVariable[permuteBack[back[i]]]=k;
	}
	// Set up permutation vector
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	ClpDisjointCopyN ( permute_, useNumberRows , pivotColumn_  );
	ClpDisjointCopyN ( permuteBack_, useNumberRows , pivotColumnBack_  );
	if (networkMatrix) {
	  maximumPivots(saveMaximumPivots);
	  // create network factorization
	  if (doCheck)
	    delete networkBasis_; // temp
	  networkBasis_ = new ClpNetworkBasis(model,numberRows_,
					      pivotRegion_,
					      permuteBack_,
					      startColumnU_,
					      numberInColumn_,
					      indexRowU_,
					      elementU_);
	  // kill off arrays in ordinary factorization
	  if (!doCheck) {
	    gutsOfDestructor();
	    status_=0;
#if 0
	    // but put back permute arrays so odd things will work
	    int numberRows = model->numberRows();
	    pivotColumnBack_ = new int [numberRows];
	    permute_ = new int [numberRows];
	    int i;
	    for (i=0;i<numberRows;i++) {
	      pivotColumnBack_[i]=i;
	      permute_[i]=i;
	    }
#endif
	  }
	} else {
	  // See if worth going sparse and when
	  if (numberFtranCounts_>100) {
	    ftranAverageAfterL_ = CoinMax(ftranCountAfterL_/ftranCountInput_,1.0);
	    ftranAverageAfterR_ = CoinMax(ftranCountAfterR_/ftranCountAfterL_,1.0);
	    ftranAverageAfterU_ = CoinMax(ftranCountAfterU_/ftranCountAfterR_,1.0);
	    assert (ftranCountInput_&&ftranCountAfterL_&&ftranCountAfterR_);
	    if (btranCountInput_&&btranCountAfterU_&&btranCountAfterR_) {
	      btranAverageAfterU_ = CoinMax(btranCountAfterU_/btranCountInput_,1.0);
	      btranAverageAfterR_ = CoinMax(btranCountAfterR_/btranCountAfterU_,1.0);
	      btranAverageAfterL_ = CoinMax(btranCountAfterL_/btranCountAfterR_,1.0);
	    } else {
	      // we have not done any useful btrans (values pass?)
	      btranAverageAfterU_ = 1.0;
	      btranAverageAfterR_ = 1.0;
	      btranAverageAfterL_ = 1.0;
	    }
	  }
	  // scale back
	  
	  ftranCountInput_ *= 0.8;
	  ftranCountAfterL_ *= 0.8;
	  ftranCountAfterR_ *= 0.8;
	  ftranCountAfterU_ *= 0.8;
	  btranCountInput_ *= 0.8;
	  btranCountAfterU_ *= 0.8;
	  btranCountAfterR_ *= 0.8;
	  btranCountAfterL_ *= 0.8;
	}
      } else if (status_ == -1&&(solveType==0||solveType==2)) {
	// This needs redoing as it was merged coding - does not need array
	int numberTotal = numberRows+numberColumns;
	int * isBasic = new int [numberTotal];
	int * rowIsBasic = isBasic+numberColumns;
	int * columnIsBasic = isBasic;
	for (i=0;i<numberTotal;i++) 
	  isBasic[i]=-1;
	for (i=0;i<numberRowBasic;i++) {
	  int iRow = pivotTemp[i];
	  rowIsBasic[iRow]=1;
	}
	for (;i<numberBasic;i++) {
	  int iColumn = pivotTemp[i];
	  columnIsBasic[iColumn]=1;
	}
	numberBasic=0;
	for (i=0;i<numberRows;i++) 
	  pivotVariable[i]=-1;
	// mark as basic or non basic
	for (i=0;i<numberRows;i++) {
	  if (rowIsBasic[i]>=0) {
	    if (pivotColumn_[numberBasic]>=0) 
	      rowIsBasic[i]=pivotColumn_[numberBasic];
	    else
	      rowIsBasic[i]=-1;
	    numberBasic++;
	  }
	}
	for (i=0;i<numberColumns;i++) {
	  if (columnIsBasic[i]>=0) {
	    if (pivotColumn_[numberBasic]>=0) 
	      columnIsBasic[i]=pivotColumn_[numberBasic];
	    else
	      columnIsBasic[i]=-1;
	    numberBasic++;
	  }
	}
	// leave pivotVariable in useful form for cleaning basis
	int * pivotVariable = model->pivotVariable();
	for (i=0;i<numberRows;i++) {
	  pivotVariable[i]=-1;
	}
	for (i=0;i<numberRows;i++) {
	  if (model->getRowStatus(i) == ClpSimplex::basic) {
	    int iPivot = rowIsBasic[i];
	    if (iPivot>=0) 
	      pivotVariable[iPivot]=i+numberColumns;
	  }
	}
	for (i=0;i<numberColumns;i++) {
	  if (model->getColumnStatus(i) == ClpSimplex::basic) {
	    int iPivot = columnIsBasic[i];
	    if (iPivot>=0) 
	      pivotVariable[iPivot]=i;
	  }
	}
	delete [] isBasic;
	double * columnLower = model->lowerRegion();
	double * columnUpper = model->upperRegion();
	double * columnActivity = model->solutionRegion();
	double * rowLower = model->lowerRegion(0);
	double * rowUpper = model->upperRegion(0);
	double * rowActivity = model->solutionRegion(0);
	//redo basis - first take ALL columns out
	int iColumn;
	double largeValue = model->largeValue();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (model->getColumnStatus(iColumn)==ClpSimplex::basic) {
	    // take out
	    if (!valuesPass) {
	      double lower=columnLower[iColumn];
	      double upper=columnUpper[iColumn];
	      double value=columnActivity[iColumn];
	      if (lower>-largeValue||upper<largeValue) {
		if (fabs(value-lower)<fabs(value-upper)) {
		  model->setColumnStatus(iColumn,ClpSimplex::atLowerBound);
		  columnActivity[iColumn]=lower;
		} else {
		  model->setColumnStatus(iColumn,ClpSimplex::atUpperBound);
		  columnActivity[iColumn]=upper;
		}
	      } else {
		model->setColumnStatus(iColumn,ClpSimplex::isFree);
	      }
	    } else {
	      model->setColumnStatus(iColumn,ClpSimplex::superBasic);
	    }
	  }
	}
	int iRow;
	for (iRow=0;iRow<numberRows;iRow++) {
	  int iSequence=pivotVariable[iRow];
	  if (iSequence>=0) {
	    // basic
	    if (iSequence>=numberColumns) {
	      // slack in - leave
	      //assert (iSequence-numberColumns==iRow);
	    } else {
	      // put back structural
	      model->setColumnStatus(iSequence,ClpSimplex::basic);
	    }
	  } else {
	    // put in slack
	    model->setRowStatus(iRow,ClpSimplex::basic);
	  }
	}
	// Put back any key variables for gub (status_ not touched)
	matrix->generalExpanded(model,1,status_);
	// signal repeat
	status_=-99;
	// set fixed if they are
	for (iRow=0;iRow<numberRows;iRow++) {
	  if (model->getRowStatus(iRow)!=ClpSimplex::basic ) {
	    if (rowLower[iRow]==rowUpper[iRow]) {
	      rowActivity[iRow]=rowLower[iRow];
	      model->setRowStatus(iRow,ClpSimplex::isFixed);
	    }
	  }
	}
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (model->getColumnStatus(iColumn)!=ClpSimplex::basic ) {
	    if (columnLower[iColumn]==columnUpper[iColumn]) {
	      columnActivity[iColumn]=columnLower[iColumn];
	      model->setColumnStatus(iColumn,ClpSimplex::isFixed);
	    }
	  }
	}
      } 
    }
  } else {
    // network - fake factorization - do nothing
    status_=0;
  }

  if (!status_) {
    // take out part if quadratic
    if (model->algorithm()==2) {
      ClpObjective * obj = model->objectiveAsObject();
      ClpQuadraticObjective * quadraticObj = (dynamic_cast< ClpQuadraticObjective*>(obj));
      assert (quadraticObj);
      CoinPackedMatrix * quadratic = quadraticObj->quadraticObjective();
      int numberXColumns = quadratic->getNumCols();
      assert (numberXColumns<numberColumns);
      int base = numberColumns-numberXColumns;
      int * which = new int [numberXColumns];
      int * pivotVariable = model->pivotVariable();
      int * permute = pivotColumn();
      int i;
      int n=0;
      for (i=0;i<numberRows;i++) {
	int iSj = pivotVariable[i]-base;
	if (iSj>=0&&iSj<numberXColumns) 
	  which[n++]=permute[i];
      }
      if (n)
	emptyRows(n,which);
      delete [] which;
    }
  }
  return status_;
}
/* Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   If checkBeforeModifying is true will do all accuracy checks
   before modifying factorization.  Whether to set this depends on
   speed considerations.  You could just do this on first iteration
   after factorization and thereafter re-factorize
   partial update already in U */
int 
ClpFactorization::replaceColumn ( const ClpSimplex * model, 
				  CoinIndexedVector * regionSparse,
				  CoinIndexedVector * tableauColumn,
				  int pivotRow,
				  double pivotCheck ,
				  bool checkBeforeModifying)
{
  if (!networkBasis_) {
    // see if FT
    if (doForrestTomlin_)
      return CoinFactorization::replaceColumn(regionSparse,
					      pivotRow,
					      pivotCheck,
					      checkBeforeModifying);
    else
      return CoinFactorization::replaceColumnPFI(tableauColumn,
					      pivotRow,pivotCheck); // Note array

  } else {
    if (doCheck) {
      int returnCode = CoinFactorization::replaceColumn(regionSparse,
							pivotRow,
							pivotCheck,
							checkBeforeModifying);
      networkBasis_->replaceColumn(regionSparse,
				   pivotRow);
      return returnCode;
    } else {
      // increase number of pivots
      numberPivots_++;
      return networkBasis_->replaceColumn(regionSparse,
				   pivotRow);
    }
  }
}

/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnFT ( CoinIndexedVector * regionSparse,
				   CoinIndexedVector * regionSparse2)
{
#ifdef CLP_DEBUG
  regionSparse->checkClear();
#endif
  if (!networkBasis_) {
    collectStatistics_ = true;
    int returnValue= CoinFactorization::updateColumnFT(regionSparse,
						       regionSparse2);
    collectStatistics_ = false;
    return returnValue;
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    int returnCode = CoinFactorization::updateColumnFT(regionSparse,
						       regionSparse2);
    networkBasis_->updateColumn(regionSparse,save);
    int i;
    double * array = regionSparse2->denseVector();
    double * array2 = save->denseVector();
    for (i=0;i<numberRows_;i++) {
      double value1 = array[i];
      double value2 = array2[i];
      assert (value1==value2);
    }
    delete save;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse,regionSparse2,-1);
    return 1;
#endif
  }
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumn ( CoinIndexedVector * regionSparse,
				 CoinIndexedVector * regionSparse2,
				 bool noPermute) const
{
#ifdef CLP_DEBUG
  if (!noPermute)
    regionSparse->checkClear();
#endif
  if (!networkBasis_) {
    collectStatistics_ = true;
    int returnValue= CoinFactorization::updateColumn(regionSparse,
						     regionSparse2,
						     noPermute);
    collectStatistics_ = false;
    return returnValue;
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    int returnCode = CoinFactorization::updateColumn(regionSparse,
						     regionSparse2,
						     noPermute);
    networkBasis_->updateColumn(regionSparse,save);
    int i;
    double * array = regionSparse2->denseVector();
    double * array2 = save->denseVector();
    for (i=0;i<numberRows_;i++) {
      double value1 = array[i];
      double value2 = array2[i];
      assert (value1==value2);
    }
    delete save;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse,regionSparse2,-1);
    return 1;
#endif
  }
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnTranspose ( CoinIndexedVector * regionSparse,
    			  CoinIndexedVector * regionSparse2) const
{
  if (!networkBasis_) {
    collectStatistics_ = true;
    return CoinFactorization::updateColumnTranspose(regionSparse,
						    regionSparse2);
    collectStatistics_ = false;
  } else {
#ifdef CHECK_NETWORK
      CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
      int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
								regionSparse2);
      networkBasis_->updateColumnTranspose(regionSparse,save);
      int i;
      double * array = regionSparse2->denseVector();
      double * array2 = save->denseVector();
      for (i=0;i<numberRows_;i++) {
	double value1 = array[i];
	double value2 = array2[i];
	assert (value1==value2);
      }
      delete save;
      return returnCode;
#else
      return networkBasis_->updateColumnTranspose(regionSparse,regionSparse2);
#endif
  }
}
/* makes a row copy of L for speed and to allow very sparse problems */
void 
ClpFactorization::goSparse()
{
  if (!networkBasis_) 
    CoinFactorization::goSparse();
}
// Cleans up i.e. gets rid of network basis 
void 
ClpFactorization::cleanUp()
{
  delete networkBasis_;
  networkBasis_=NULL;
  resetStatistics();
}
/// Says whether to redo pivot order
bool 
ClpFactorization::needToReorder() const
{
#ifdef CHECK_NETWORK
  return true;
#endif
  if (!networkBasis_)
    return true;
  else
    return false;
}
