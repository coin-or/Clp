// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "ClpFactorization.hpp"
#ifndef SLIM_CLP
#include "ClpQuadraticObjective.hpp"
#endif
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpSimplex.hpp"
#include "ClpMatrixBase.hpp"
#ifndef SLIM_CLP
#include "ClpNetworkBasis.hpp"
#include "ClpNetworkMatrix.hpp"
//#define CHECK_NETWORK
#ifdef CHECK_NETWORK
const static bool doCheck=true;
#else
const static bool doCheck=false;
#endif
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
#ifndef SLIM_CLP
  networkBasis_ = NULL;
#endif
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpFactorization::ClpFactorization (const ClpFactorization & rhs) :
   CoinFactorization(rhs) 
{
#ifndef SLIM_CLP
  if (rhs.networkBasis_)
    networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
  else
    networkBasis_=NULL;
#endif
}

ClpFactorization::ClpFactorization (const CoinFactorization & rhs) :
   CoinFactorization(rhs) 
{
#ifndef SLIM_CLP
  networkBasis_=NULL;
#endif
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpFactorization::~ClpFactorization () 
{
#ifndef SLIM_CLP
  delete networkBasis_;
#endif
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpFactorization &
ClpFactorization::operator=(const ClpFactorization& rhs)
{
  if (this != &rhs) {
    CoinFactorization::operator=(rhs);
#ifndef SLIM_CLP
    delete networkBasis_;
    if (rhs.networkBasis_)
      networkBasis_ = new ClpNetworkBasis(*(rhs.networkBasis_));
    else
      networkBasis_=NULL;
#endif
  }
  return *this;
}
#if 0
static unsigned int saveList[10000];
int numberSave=-1;
inline bool isDense(int i) {
  return ((saveList[i>>5]>>(i&31))&1)!=0;
}
inline void setDense(int i) {
  unsigned int & value = saveList[i>>5];
  int bit = i&31;
  value |= (1<<bit);
}
#endif
int 
ClpFactorization::factorize ( ClpSimplex * model,
			      int solveType, bool valuesPass)
{
  ClpMatrixBase * matrix = model->clpMatrix(); 
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  if (!numberRows)
    return 0;
  // If too many compressions increase area
  if (numberPivots_>1&&numberCompressions_*10 > numberPivots_+10) {
    areaFactor_ *= 1.1;
  }
  //int numberPivots=numberPivots_;
#if 0
  if (model->algorithm()>0)
    numberSave=-1;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_||doCheck) {
#endif
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
#if 0
        if (numberSave>0) {
          int nStill=0;
          int nAtBound=0;
          int nZeroDual=0;
          CoinIndexedVector * array = model->rowArray(3);
          CoinIndexedVector * objArray = model->columnArray(1);
          array->clear();
          objArray->clear();
          double * cost = model->costRegion();
          double tolerance = model->primalTolerance();
          double offset=0.0;
          for (i=0;i<numberRows;i++) {
            int iPivot = pivotVariable[i];
            if (iPivot<numberColumns&&isDense(iPivot)) {
              if (model->getColumnStatus(iPivot)==ClpSimplex::basic) {
                nStill++;
                double value=model->solutionRegion()[iPivot];
                double dual = model->dualRowSolution()[i];
                double lower=model->lowerRegion()[iPivot];
                double upper=model->upperRegion()[iPivot];
                ClpSimplex::Status status;
                if (fabs(value-lower)<tolerance) {
                  status=ClpSimplex::atLowerBound;
                  nAtBound++;
                } else if (fabs(value-upper)<tolerance) {
                  nAtBound++;
                  status=ClpSimplex::atUpperBound;
                } else if (value>lower&&value<upper) {
                  status=ClpSimplex::superBasic;
                } else {
                  status=ClpSimplex::basic;
                }
                if (status!=ClpSimplex::basic) {
                  if (model->getRowStatus(i)!=ClpSimplex::basic) {
                    model->setColumnStatus(iPivot,ClpSimplex::atLowerBound);
                    model->setRowStatus(i,ClpSimplex::basic);
                    pivotVariable[i]=i+numberColumns;
                    model->dualRowSolution()[i]=0.0;
                    model->djRegion(0)[i]=0.0;
                    array->add(i,dual);
                    offset += dual*model->solutionRegion(0)[i];
                  }
                }
                if (fabs(dual)<1.0e-5)
                  nZeroDual++;
              }
            }
          }
          printf("out of %d dense, %d still in basis, %d at bound, %d with zero dual - offset %g\n",
                 numberSave,nStill,nAtBound,nZeroDual,offset);
          if (array->getNumElements()) {
            // modify costs
            model->clpMatrix()->transposeTimes(model,1.0,array,model->columnArray(0),
                                               objArray);
            array->clear();
            int n=objArray->getNumElements();
            int * indices = objArray->getIndices();
            double * elements = objArray->denseVector();
            for (i=0;i<n;i++) {
              int iColumn = indices[i];
              cost[iColumn] -= elements[iColumn];
              elements[iColumn]=0.0;
            }
            objArray->setNumElements(0);
          }
        }
#endif
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
      if (numberBasic>model->maximumBasic()) {
#if 0 // ndef NDEBUG
        printf("%d basic - should only be %d\n",
               numberBasic,numberRows);
#endif
        // Take out some
        numberBasic=numberRowBasic;
        for (int i=0;i<numberColumns;i++) {
          if (model->getColumnStatus(i) == ClpSimplex::basic) {
            if (numberBasic<numberRows)
              numberBasic++;
            else
              model->setColumnStatus(i,ClpSimplex::superBasic);
          }
        }
        numberBasic=numberRowBasic;
        matrix->generalExpanded(model,0,numberBasic);
      }
#ifndef SLIM_CLP
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
      if (networkMatrix)
	biasLU_=0; // All to U if network
      //int saveMaximumPivots = maximumPivots();
      delete networkBasis_;
      networkBasis_ = NULL;
      if (networkMatrix&&!doCheck)
	maximumPivots(1);
#endif
      //printf("L, U, R %d %d %d\n",numberElementsL(),numberElementsU(),numberElementsR());
      while (status_==-99) {
	// maybe for speed will be better to leave as many regions as possible
	gutsOfDestructor();
	gutsOfInitialize(2);
	CoinBigIndex numberElements=numberRowBasic;

	// compute how much in basis

	int i;
	// can change for gub
	int numberColumnBasic = numberBasic-numberRowBasic;

	numberElements +=matrix->countBasis(model,
					   pivotTemp+numberRowBasic, 
					   numberRowBasic,
					    numberColumnBasic);
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
	// Fill in counts so we can skip part of preProcess
	int * numberInRow = numberInRow_.array();
	int * numberInColumn = numberInColumn_.array();
	CoinZeroN ( numberInRow, numberRows_ + 1 );
	CoinZeroN ( numberInColumn, maximumColumnsExtra_ + 1 );
	double * elementU = elementU_.array();
	int * indexRowU = indexRowU_.array();
	CoinBigIndex * startColumnU = startColumnU_.array();
	for (i=0;i<numberRowBasic;i++) {
	  int iRow = pivotTemp[i];
	  indexRowU[i]=iRow;
	  startColumnU[i]=i;
	  elementU[i]=slackValue_;
	  numberInRow[iRow]=1;
	  numberInColumn[i]=1;
	}
	startColumnU[numberRowBasic]=numberRowBasic;
	// can change for gub so redo
	numberColumnBasic = numberBasic-numberRowBasic;
	matrix->fillBasis(model, 
			  pivotTemp+numberRowBasic, 
			  numberColumnBasic,
			  indexRowU_.array(), 
			  startColumnU+numberRowBasic,
			  numberInRow,
			  numberInColumn+numberRowBasic,
			  elementU_.array());
#if 0
	{
	  printf("%d row basic, %d column basic\n",numberRowBasic,numberColumnBasic);
	  for (int i=0;i<numberElements;i++) 
	    printf("row %d col %d value %g\n",indexRowU_.array()[i],indexColumnU_[i],
		   elementU_.array()[i]);
	}
#endif
	// recompute number basic
        numberBasic = numberRowBasic+numberColumnBasic;
	if (numberBasic) 
	  numberElements = startColumnU[numberBasic-1]
	    +numberInColumn[numberBasic-1];
	else
	  numberElements=0;
	lengthU_ = numberElements;
        //saveFactorization("dump.d");
	if (biasLU_>=3||numberRows_!=numberColumns_)
	  preProcess ( 2 );
	else
	  preProcess ( 3 ); // no row copy
	factor (  );
	if (status_==-99) {
	  // get more memory
	  areaFactor(2.0*areaFactor());
	} else if (status_==-1&&model->numberIterations()==0&&
                   denseThreshold_) {
          // Round again without dense
          denseThreshold_=0;
          status_=-99;
        }
      }
      // If we get here status is 0 or -1
      if (status_ == 0) {
	// We may need to tamper with order and redo - e.g. network with side
	int useNumberRows = numberRows;
	// **** we will also need to add test in dual steepest to do
	// as we do for network
        matrix->generalExpanded(model,12,useNumberRows);
	const int * permuteBack = permuteBack_.array();
	const int * back = pivotColumnBack_.array();
	//int * pivotTemp = pivotColumn_.array();
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
#if 0
        if (numberSave>=0) {
          numberSave=numberDense_;
          memset(saveList,0,((numberRows_+31)>>5)*sizeof(int));
          for (i=numberRows_-numberSave;i<numberRows_;i++) {
            int k=pivotTemp[pivotColumn_.array()[i]];
            setDense(k);
          }
        }
#endif
	// Set up permutation vector
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	ClpDisjointCopyN ( permute_.array(), useNumberRows , pivotColumn_.array()  );
	ClpDisjointCopyN ( permuteBack_.array(), useNumberRows , pivotColumnBack_.array()  );
#ifndef SLIM_CLP
	if (networkMatrix) {
	  maximumPivots(CoinMax(2000,maximumPivots()));
	  // redo arrays
	  for (int iRow=0;iRow<4;iRow++) {
	    int length =model->numberRows()+maximumPivots();
	    if (iRow==3||model->objectiveAsObject()->type()>1)
	      length += model->numberColumns();
	    model->rowArray(iRow)->reserve(length);
	  }
	  // create network factorization
	  if (doCheck)
	    delete networkBasis_; // temp
	  networkBasis_ = new ClpNetworkBasis(model,numberRows_,
					      pivotRegion_.array(),
					      permuteBack_.array(),
					      startColumnU_.array(),
					      numberInColumn_.array(),
					      indexRowU_.array(),
					      elementU_.array());
	  // kill off arrays in ordinary factorization
	  if (!doCheck) {
	    gutsOfDestructor();
	    // but make sure numberRows_ set
	    numberRows_ = model->numberRows();
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
#endif
	  // See if worth going sparse and when
	  if (numberFtranCounts_>100) {
	    ftranCountInput_= CoinMax(ftranCountInput_,1.0);
	    ftranAverageAfterL_ = CoinMax(ftranCountAfterL_/ftranCountInput_,1.0);
	    ftranAverageAfterR_ = CoinMax(ftranCountAfterR_/ftranCountAfterL_,1.0);
	    ftranAverageAfterU_ = CoinMax(ftranCountAfterU_/ftranCountAfterR_,1.0);
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
#ifndef SLIM_CLP
	}
#endif
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
	const int * pivotColumn = pivotColumn_.array();
	for (i=0;i<numberRows;i++) {
	  if (rowIsBasic[i]>=0) {
	    if (pivotColumn[numberBasic]>=0) {
	      rowIsBasic[i]=pivotColumn[numberBasic];
	    } else {
	      rowIsBasic[i]=-1;
              model->setRowStatus(i,ClpSimplex::superBasic);
            }
	    numberBasic++;
	  }
	}
	for (i=0;i<numberColumns;i++) {
	  if (columnIsBasic[i]>=0) {
	    if (pivotColumn[numberBasic]>=0) 
	      columnIsBasic[i]=pivotColumn[numberBasic];
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
              assert(model->getRowStatus(iRow)!=ClpSimplex::basic);
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
#ifndef SLIM_CLP
  } else {
    // network - fake factorization - do nothing
    status_=0;
    numberPivots_ = 0;
  }
#endif
#ifndef SLIM_CLP
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
#endif
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
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
    // see if FT
    if (doForrestTomlin_) {
      int returnCode= CoinFactorization::replaceColumn(regionSparse,
					      pivotRow,
					      pivotCheck,
					      checkBeforeModifying);
      return returnCode;
    } else {
      return CoinFactorization::replaceColumnPFI(tableauColumn,
					      pivotRow,pivotCheck); // Note array
    }

#ifndef SLIM_CLP
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
#endif
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
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
    collectStatistics_ = true;
    int returnValue= CoinFactorization::updateColumnFT(regionSparse,
						       regionSparse2);
    collectStatistics_ = false;
    return returnValue;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    double * check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumnFT(regionSparse,
						       regionSparse2);
    networkBasis_->updateColumn(regionSparse,save,-1);
    int i;
    double * array = regionSparse2->denseVector();
    int * indices = regionSparse2->getIndices();
    int n=regionSparse2->getNumElements();
    memset(check,0,numberRows_*sizeof(double));
    double * array2 = save->denseVector();
    int * indices2 = save->getIndices();
    int n2=save->getNumElements();
    assert (n==n2);
    if (save->packedMode()) {
      for (i=0;i<n;i++) {
	check[indices[i]]=array[i];
      }
      for (i=0;i<n;i++) {
	double value2 = array2[i];
	assert (check[indices2[i]]==value2);
      }
    } else {
      for (i=0;i<numberRows_;i++) {
	double value1 = array[i];
	double value2 = array2[i];
	assert (value1==value2);
      }
    }
    delete save;
    delete [] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse,regionSparse2,-1);
    return 1;
#endif
  }
#endif
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
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
    collectStatistics_ = true;
    int returnValue= CoinFactorization::updateColumn(regionSparse,
						     regionSparse2,
						     noPermute);
    collectStatistics_ = false;
    return returnValue;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    double * check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumn(regionSparse,
						     regionSparse2,
						     noPermute);
    networkBasis_->updateColumn(regionSparse,save,-1);
    int i;
    double * array = regionSparse2->denseVector();
    int * indices = regionSparse2->getIndices();
    int n=regionSparse2->getNumElements();
    memset(check,0,numberRows_*sizeof(double));
    double * array2 = save->denseVector();
    int * indices2 = save->getIndices();
    int n2=save->getNumElements();
    assert (n==n2);
    if (save->packedMode()) {
      for (i=0;i<n;i++) {
	check[indices[i]]=array[i];
      }
      for (i=0;i<n;i++) {
	double value2 = array2[i];
	assert (check[indices2[i]]==value2);
      }
    } else {
      for (i=0;i<numberRows_;i++) {
	double value1 = array[i];
	double value2 = array2[i];
	assert (value1==value2);
      }
    }
    delete save;
    delete [] check;
    return returnCode;
#else
    networkBasis_->updateColumn(regionSparse,regionSparse2,-1);
    return 1;
#endif
  }
#endif
}
/* Updates one column (FTRAN) from region2
   number returned is negative if no room
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnForDebug ( CoinIndexedVector * regionSparse,
				 CoinIndexedVector * regionSparse2,
				 bool noPermute) const
{
  if (!noPermute)
    regionSparse->checkClear();
  if (!numberRows_)
    return 0;
  collectStatistics_ = false;
  int returnValue= CoinFactorization::updateColumn(regionSparse,
                                                   regionSparse2,
                                                   noPermute);
  return returnValue;
}
/* Updates one column (BTRAN) from region2
   region1 starts as zero and is zero at end */
int 
ClpFactorization::updateColumnTranspose ( CoinIndexedVector * regionSparse,
    			  CoinIndexedVector * regionSparse2) const
{
  if (!numberRows_)
    return 0;
#ifndef SLIM_CLP
  if (!networkBasis_) {
#endif
    collectStatistics_ = true;
    return CoinFactorization::updateColumnTranspose(regionSparse,
						    regionSparse2);
    collectStatistics_ = false;
#ifndef SLIM_CLP
  } else {
#ifdef CHECK_NETWORK
    CoinIndexedVector * save = new CoinIndexedVector(*regionSparse2);
    double * check = new double[numberRows_];
    int returnCode = CoinFactorization::updateColumnTranspose(regionSparse,
							      regionSparse2);
    networkBasis_->updateColumnTranspose(regionSparse,save);
    int i;
    double * array = regionSparse2->denseVector();
    int * indices = regionSparse2->getIndices();
    int n=regionSparse2->getNumElements();
    memset(check,0,numberRows_*sizeof(double));
    double * array2 = save->denseVector();
    int * indices2 = save->getIndices();
    int n2=save->getNumElements();
    assert (n==n2);
    if (save->packedMode()) {
      for (i=0;i<n;i++) {
	check[indices[i]]=array[i];
      }
      for (i=0;i<n;i++) {
	double value2 = array2[i];
	assert (check[indices2[i]]==value2);
      }
    } else {
      for (i=0;i<numberRows_;i++) {
	double value1 = array[i];
	double value2 = array2[i];
	assert (value1==value2);
      }
    }
    delete save;
    delete [] check;
    return returnCode;
#else
    return networkBasis_->updateColumnTranspose(regionSparse,regionSparse2);
#endif
  }
#endif
}
/* makes a row copy of L for speed and to allow very sparse problems */
void 
ClpFactorization::goSparse()
{
#ifndef SLIM_CLP
  if (!networkBasis_) 
#endif
    CoinFactorization::goSparse();
}
// Cleans up i.e. gets rid of network basis 
void 
ClpFactorization::cleanUp()
{
#ifndef SLIM_CLP
  delete networkBasis_;
  networkBasis_=NULL;
#endif
  resetStatistics();
}
/// Says whether to redo pivot order
bool 
ClpFactorization::needToReorder() const
{
#ifdef CHECK_NETWORK
  return true;
#endif
#ifndef SLIM_CLP
  if (!networkBasis_)
#endif
    return true;
#ifndef SLIM_CLP
  else
    return false;
#endif
}
// Get weighted row list 
void
ClpFactorization::getWeights(int * weights) const
{
#ifndef SLIM_CLP
  if (networkBasis_) {
    // Network - just unit
    for (int i=0;i<numberRows_;i++) 
      weights[i]=1;
    return;
  }
#endif
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  int * permuteBack = pivotColumnBack_.array();
  int * indexRowU = indexRowU_.array();
  const CoinBigIndex * startColumnU = startColumnU_.array();
  const CoinBigIndex * startRowL = startRowL_.array();
  if (!startRowL||!numberInRow_.array()) {
    int * temp = new int[numberRows_];
    memset(temp,0,numberRows_*sizeof(int));
    int i;
    for (i=0;i<numberRows_;i++) {
      // one for pivot
      temp[i]++;
      CoinBigIndex j;
      for (j=startColumnU[i];j<startColumnU[i]+numberInColumn[i];j++) {
	int iRow=indexRowU[j];
	temp[iRow]++;
      }
    }
    CoinBigIndex * startColumnL = startColumnL_.array();
    int * indexRowL = indexRowL_.array();
    for (i=baseL_;i<baseL_+numberL_;i++) {
      CoinBigIndex j;
      for (j=startColumnL[i];j<startColumnL[i+1];j++) {
	int iRow = indexRowL[j];
	temp[iRow]++;
      }
    }
    for (i=0;i<numberRows_;i++) {
      int number = temp[i];
      int iPermute = permuteBack[i];
      weights[iPermute]=number;
    }
    delete [] temp;
  } else {
    int i;
    for (i=0;i<numberRows_;i++) {
      int number = startRowL[i+1]-startRowL[i]+numberInRow[i]+1;
      //number = startRowL[i+1]-startRowL[i]+1;
      //number = numberInRow[i]+1;
      int iPermute = permuteBack[i];
      weights[iPermute]=number;
    }
  }
}
