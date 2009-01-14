// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
   
#include "CoinPragma.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

int boundary_sort=1000;
int boundary_sort2=1000;
int boundary_sort3=10000;

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
// History since 1.0 at end
#define CLPVERSION "1.06.00"

#include "CoinMpsIO.hpp"
#include "CoinFileIO.hpp"

#include "ClpFactorization.hpp"
#include "CoinTime.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPresolve.hpp"
#include "CbcOrClpParam.hpp"
#include "CoinSignal.hpp"
#ifdef DMALLOC
#include "dmalloc.h"
#endif
#ifdef WSSMP_BARRIER
#define FOREIGN_BARRIER
#endif
#ifdef UFL_BARRIER
#define FOREIGN_BARRIER
#endif
#ifdef TAUCS_BARRIER
#define FOREIGN_BARRIER
#endif

static double totalTime=0.0;
static bool maskMatches(const int * starts, char ** masks,
			std::string & check);
static ClpSimplex * currentModel = NULL;

extern "C" {
   static void 
#if defined(_MSC_VER)
   __cdecl
#endif // _MSC_VER
   signal_handler(int whichSignal)
   {
      if (currentModel!=NULL) 
	 currentModel->setMaximumIterations(0); // stop at next iterations
      return;
   }
}

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

int mainTest (int argc, const char *argv[],int algorithm,
	      ClpSimplex empty, bool doPresolve,int switchOff,bool doVector);
static void statistics(ClpSimplex * originalModel, ClpSimplex * model);
static void generateCode(const char * fileName,int type);
// Returns next valid field
int CbcOrClpRead_mode=1;
FILE * CbcOrClpReadCommand=stdin;

int
#if defined(_MSC_VER)
__cdecl
#endif // _MSC_VER
main (int argc, const char *argv[])
{
  // next {} is just to make sure all memory should be freed - for debug
  {
    double time1 = CoinCpuTime(),time2;
    // Set up all non-standard stuff
    //int numberModels=1;
    ClpSimplex * models = new ClpSimplex[1];
    
    // default action on import
    int allowImportErrors=0;
    int keepImportNames=1;
    int doIdiot=-1;
    int outputFormat=2;
    int slpValue=-1;
    int cppValue=-1;
    int printOptions=0;
    int printMode=0;
    int presolveOptions=0;
    int doCrash=0;
    int doVector=0;
    int doSprint=-1;
    // set reasonable defaults
    int preSolve=5;
    bool preSolveFile=false;
    models->setPerturbation(50);
    models->messageHandler()->setPrefix(false);
    const char dirsep =  CoinFindDirSeparator();
    std::string directory;
    std::string dirSample;
    std::string dirNetlib;
    std::string dirMiplib;
    if (dirsep == '/') {
      directory = "./";
      dirSample = "../../Data/Sample/";
      dirNetlib = "../../Data/Netlib/";
      dirMiplib = "../../Data/miplib3/";
    } else {
      directory = ".\\";
      dirSample = "..\\..\\Data\\Sample\\";
      dirNetlib = "..\\..\\Data\\Netlib\\";
      dirMiplib = "..\\..\\Data\\miplib3\\";
    }
    std::string defaultDirectory = directory;
    std::string importFile ="";
    std::string exportFile ="default.mps";
    std::string importBasisFile ="";
    int basisHasValues=0;
    int substitution=3;
    int dualize=3;  // dualize if looks promising
    std::string exportBasisFile ="default.bas";
    std::string saveFile ="default.prob";
    std::string restoreFile ="default.prob";
    std::string solutionFile ="stdout";
    std::string solutionSaveFile ="solution.file";
    std::string printMask="";
    CbcOrClpParam parameters[CBCMAXPARAMETERS];
    int numberParameters ;
    establishParams(numberParameters,parameters) ;
    parameters[whichParam(BASISIN,numberParameters,parameters)].setStringValue(importBasisFile);
    parameters[whichParam(BASISOUT,numberParameters,parameters)].setStringValue(exportBasisFile);
    parameters[whichParam(PRINTMASK,numberParameters,parameters)].setStringValue(printMask);
    parameters[whichParam(DIRECTORY,numberParameters,parameters)].setStringValue(directory);
    parameters[whichParam(DIRSAMPLE,numberParameters,parameters)].setStringValue(dirSample);
    parameters[whichParam(DIRNETLIB,numberParameters,parameters)].setStringValue(dirNetlib);
    parameters[whichParam(DIRMIPLIB,numberParameters,parameters)].setStringValue(dirMiplib);
    parameters[whichParam(DUALBOUND,numberParameters,parameters)].setDoubleValue(models->dualBound());
    parameters[whichParam(DUALTOLERANCE,numberParameters,parameters)].setDoubleValue(models->dualTolerance());
    parameters[whichParam(EXPORT,numberParameters,parameters)].setStringValue(exportFile);
    parameters[whichParam(IDIOT,numberParameters,parameters)].setIntValue(doIdiot);
    parameters[whichParam(IMPORT,numberParameters,parameters)].setStringValue(importFile);
    parameters[whichParam(SOLVERLOGLEVEL,numberParameters,parameters)].setIntValue(models->logLevel());
    parameters[whichParam(MAXFACTOR,numberParameters,parameters)].setIntValue(models->factorizationFrequency());
    parameters[whichParam(MAXITERATION,numberParameters,parameters)].setIntValue(models->maximumIterations());
    parameters[whichParam(OUTPUTFORMAT,numberParameters,parameters)].setIntValue(outputFormat);
    parameters[whichParam(PRESOLVEPASS,numberParameters,parameters)].setIntValue(preSolve);
    parameters[whichParam(PERTVALUE,numberParameters,parameters)].setIntValue(models->perturbation());
    parameters[whichParam(PRIMALTOLERANCE,numberParameters,parameters)].setDoubleValue(models->primalTolerance());
    parameters[whichParam(PRIMALWEIGHT,numberParameters,parameters)].setDoubleValue(models->infeasibilityCost());
    parameters[whichParam(RESTORE,numberParameters,parameters)].setStringValue(restoreFile);
    parameters[whichParam(SAVE,numberParameters,parameters)].setStringValue(saveFile);
    parameters[whichParam(TIMELIMIT,numberParameters,parameters)].setDoubleValue(models->maximumSeconds());
    parameters[whichParam(SOLUTION,numberParameters,parameters)].setStringValue(solutionFile);
    parameters[whichParam(SAVESOL,numberParameters,parameters)].setStringValue(solutionSaveFile);
    parameters[whichParam(SPRINT,numberParameters,parameters)].setIntValue(doSprint);
    parameters[whichParam(SUBSTITUTION,numberParameters,parameters)].setIntValue(substitution);
    parameters[whichParam(DUALIZE,numberParameters,parameters)].setIntValue(dualize);
    parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].setDoubleValue(1.0e-8);
    int verbose=0;
    
    // total number of commands read
    int numberGoodCommands=0;
    bool * goodModels = new bool[1];

    // Hidden stuff for barrier
    int choleskyType = 0;
    int gamma=0;
    int scaleBarrier=0;
    int doKKT=0;
    int crossover=2; // do crossover unless quadratic
    
    int iModel=0;
    goodModels[0]=false;
    //models[0].scaling(1);
    //models[0].setDualBound(1.0e6);
    //models[0].setDualTolerance(1.0e-7);
    //ClpDualRowSteepest steep;
    //models[0].setDualRowPivotAlgorithm(steep);
    //models[0].setPrimalTolerance(1.0e-7);
    //ClpPrimalColumnSteepest steepP;
    //models[0].setPrimalColumnPivotAlgorithm(steepP);
    std::string field;
    std::cout<<"Coin LP version "<<CLPVERSION
	     <<", build "<<__DATE__<<std::endl;
    // Print command line
    if (argc>1) {
      printf("command line - ");
      for (int i=0;i<argc;i++)
        printf("%s ",argv[i]);
      printf("\n");
    }
    
    while (1) {
      // next command
      field=CoinReadGetCommand(argc,argv);
      
      // exit if null or similar
      if (!field.length()) {
	if (numberGoodCommands==1&&goodModels[0]) {
	  // we just had file name - do dual or primal
	  field="either";
	} else if (!numberGoodCommands) {
	  // let's give the sucker a hint
	  std::cout
	    <<"Clp takes input from arguments ( - switches to stdin)"
	    <<std::endl
	    <<"Enter ? for list of commands or help"<<std::endl;
	  field="-";
	} else {
	  break;
	}
      }
      
      // see if ? at end
      int numberQuery=0;
      if (field!="?"&&field!="???") {
	int length = field.length();
	int i;
	for (i=length-1;i>0;i--) {
	  if (field[i]=='?') 
	    numberQuery++;
	  else
	    break;
	}
	field=field.substr(0,length-numberQuery);
      }
      // find out if valid command
      int iParam;
      int numberMatches=0;
      int firstMatch=-1;
      for ( iParam=0; iParam<numberParameters; iParam++ ) {
	int match = parameters[iParam].matches(field);
	if (match==1) {
	  numberMatches = 1;
	  firstMatch=iParam;
	  break;
	} else {
	  if (match&&firstMatch<0)
	    firstMatch=iParam;
	  numberMatches += match>>1;
	}
      }
      if (iParam<numberParameters&&!numberQuery) {
	// found
	CbcOrClpParam found = parameters[iParam];
	CbcOrClpParameterType type = found.type();
	int valid;
	numberGoodCommands++;
	if (type==GENERALQUERY) {
	  std::cout<<"In argument list keywords have leading - "
	    ", -stdin or just - switches to stdin"<<std::endl;
	  std::cout<<"One command per line (and no -)"<<std::endl;
	  std::cout<<"abcd? gives list of possibilities, if only one + explanation"<<std::endl;
	  std::cout<<"abcd?? adds explanation, if only one fuller help"<<std::endl;
	  std::cout<<"abcd without value (where expected) gives current value"<<std::endl;
	  std::cout<<"abcd value sets value"<<std::endl;
	  std::cout<<"Commands are:"<<std::endl;
	  int maxAcross=5;
	  bool evenHidden=false;
	  if ((verbose&8)!=0) {
	    // even hidden
	    evenHidden = true;
	    verbose &= ~8;
	  }
          if (verbose)
            maxAcross=1;
	  int limits[]={1,101,201,301,401};
	  std::vector<std::string> types;
	  types.push_back("Double parameters:");
	  types.push_back("Int parameters:");
	  types.push_back("Keyword parameters:");
	  types.push_back("Actions or string parameters:");
	  int iType;
	  for (iType=0;iType<4;iType++) {
	    int across=0;
            if ((verbose%4)!=0)
              std::cout<<std::endl;
	    std::cout<<types[iType]<<std::endl;
            if ((verbose&2)!=0)
              std::cout<<std::endl;
	    for ( iParam=0; iParam<numberParameters; iParam++ ) {
	      int type = parameters[iParam].type();
	      if ((parameters[iParam].displayThis()||evenHidden)&&
		  type>=limits[iType]
		  &&type<limits[iType+1]) {
		if (!across) {
                  if ((verbose&2)==0) 
                    std::cout<<"  ";
                  else
                    std::cout<<"Command ";
                }
                std::cout<<parameters[iParam].matchName()<<"  ";
		across++;
		if (across==maxAcross) {
		  across=0;
                  if (verbose) {
                    // put out description as well
                    if ((verbose&1)!=0) 
                      std::cout<<parameters[iParam].shortHelp();
                    std::cout<<std::endl;
                    if ((verbose&2)!=0) {
                      std::cout<<"---- description"<<std::endl;
                      parameters[iParam].printLongHelp();
                      std::cout<<"----"<<std::endl<<std::endl;
                    }
                  } else {
                    std::cout<<std::endl;
                  }
		}
	      }
	    }
	    if (across)
	      std::cout<<std::endl;
	  }
	} else if (type==FULLGENERALQUERY) {
	  std::cout<<"Full list of commands is:"<<std::endl;
	  int maxAcross=5;
	  int limits[]={1,101,201,301,401};
	  std::vector<std::string> types;
	  types.push_back("Double parameters:");
	  types.push_back("Int parameters:");
	  types.push_back("Keyword parameters and others:");
	  types.push_back("Actions:");
	  int iType;
	  for (iType=0;iType<4;iType++) {
	    int across=0;
	    std::cout<<types[iType]<<std::endl;
	    for ( iParam=0; iParam<numberParameters; iParam++ ) {
	      int type = parameters[iParam].type();
	      if (type>=limits[iType]
		  &&type<limits[iType+1]) {
		if (!across)
		  std::cout<<"  ";
		std::cout<<parameters[iParam].matchName()<<"  ";
		across++;
		if (across==maxAcross) {
		  std::cout<<std::endl;
		  across=0;
		}
	      }
	    }
	    if (across)
	      std::cout<<std::endl;
	  }
	} else if (type<101) {
	  // get next field as double
	  double value = CoinReadGetDoubleField(argc,argv,&valid);
	  if (!valid) {
	    parameters[iParam].setDoubleParameter(models+iModel,value);
	  } else if (valid==1) {
	    std::cout<<" is illegal for double parameter "<<parameters[iParam].name()<<" value remains "<<
	      parameters[iParam].doubleValue()<<std::endl;
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].doubleValue()<<std::endl;
	  }
	} else if (type<201) {
	  // get next field as int
	  int value = CoinReadGetIntField(argc,argv,&valid);
	  if (!valid) {
	    if (parameters[iParam].type()==PRESOLVEPASS)
	      preSolve = value;
	    else if (parameters[iParam].type()==IDIOT)
	      doIdiot = value;
	    else if (parameters[iParam].type()==SPRINT)
	      doSprint = value;
	    else if (parameters[iParam].type()==OUTPUTFORMAT)
	      outputFormat = value;
	    else if (parameters[iParam].type()==SLPVALUE)
	      slpValue = value;
	    else if (parameters[iParam].type()==CPP)
	      cppValue = value;
	    else if (parameters[iParam].type()==PRESOLVEOPTIONS)
	      presolveOptions = value;
	    else if (parameters[iParam].type()==PRINTOPTIONS)
	      printOptions = value;
	    else if (parameters[iParam].type()==SUBSTITUTION)
	      substitution = value;
	    else if (parameters[iParam].type()==DUALIZE)
	      dualize = value;
            else if (parameters[iParam].type()==VERBOSE)
              verbose = value;
            parameters[iParam].setIntParameter(models+iModel,value);
	  } else if (valid==1) {
	    std::cout<<" is illegal for integer parameter "<<parameters[iParam].name()<<" value remains "<<
	      parameters[iParam].intValue()<<std::endl;
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].intValue()<<std::endl;
	  }
	} else if (type<301) {
	  // one of several strings
	  std::string value = CoinReadGetString(argc,argv);
	  int action = parameters[iParam].parameterOption(value);
	  if (action<0) {
	    if (value!="EOL") {
	      // no match
	      parameters[iParam].printOptions();
	    } else {
	      // print current value
	      std::cout<<parameters[iParam].name()<<" has value "<<
		parameters[iParam].currentOption()<<std::endl;
	    }
	  } else {
	    parameters[iParam].setCurrentOption(action);
	    // for now hard wired
	    switch (type) {
	    case DIRECTION:
	      if (action==0)
		models[iModel].setOptimizationDirection(1);
	      else if (action==1)
		models[iModel].setOptimizationDirection(-1);
	      else
		models[iModel].setOptimizationDirection(0);
	      break;
	    case DUALPIVOT:
	      if (action==0) {
		ClpDualRowSteepest steep(3);
		models[iModel].setDualRowPivotAlgorithm(steep);
	      } else if (action==1) {
		//ClpDualRowDantzig dantzig;
		ClpDualRowSteepest dantzig(5);
		models[iModel].setDualRowPivotAlgorithm(dantzig);
	      } else if (action==2) {
		// partial steep
		ClpDualRowSteepest steep(2);
		models[iModel].setDualRowPivotAlgorithm(steep);
	      } else {
		ClpDualRowSteepest steep;
		models[iModel].setDualRowPivotAlgorithm(steep);
	      }
	      break;
	    case PRIMALPIVOT:
	      if (action==0) {
		ClpPrimalColumnSteepest steep(3);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==1) {
		ClpPrimalColumnSteepest steep(0);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==2) {
		ClpPrimalColumnDantzig dantzig;
		models[iModel].setPrimalColumnPivotAlgorithm(dantzig);
	      } else if (action==3) {
		ClpPrimalColumnSteepest steep(4);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==4) {
		ClpPrimalColumnSteepest steep(1);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==5) {
		ClpPrimalColumnSteepest steep(2);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==6) {
		ClpPrimalColumnSteepest steep(10);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      }
	      break;
	    case SCALING:
	      models[iModel].scaling(action);
	      break;
	    case AUTOSCALE:
	      models[iModel].setAutomaticScaling(action!=0);
	      break;
	    case SPARSEFACTOR:
	      models[iModel].setSparseFactorization((1-action)!=0);
	      break;
	    case BIASLU:
	      models[iModel].factorization()->setBiasLU(action);
	      break;
	    case PERTURBATION:
	      if (action==0)
		models[iModel].setPerturbation(50);
	      else
		models[iModel].setPerturbation(100);
	      break;
	    case ERRORSALLOWED:
	      allowImportErrors = action;
	      break;
            case INTPRINT:
              printMode=action;
              break;
	    case KEEPNAMES:
	      keepImportNames = 1-action;
	      break;
	    case PRESOLVE:
	      if (action==0)
		preSolve = 5;
	      else if (action==1)
		preSolve=0;
	      else if (action==2)
		preSolve=10;
	      else
		preSolveFile=true;
	      break;
	    case PFI:
	      models[iModel].factorization()->setForrestTomlin(action==0);
	      break;
	    case CRASH:
	      doCrash=action;
	      break;
	    case VECTOR:
	      doVector=action;
	      break;
	    case MESSAGES:
	      models[iModel].messageHandler()->setPrefix(action!=0);
	      break;
	    case CHOLESKY:
	      choleskyType = action;
	      break;
	    case GAMMA:
	      gamma=action;
	      break;
	    case BARRIERSCALE:
	      scaleBarrier=action;
	      break;
	    case KKT:
	      doKKT=action;
	      break;
	    case CROSSOVER:
	      crossover=action;
	      break;
	    default:
	      abort();
	    }
	  }
	} else {
	  // action
	  if (type==EXIT)
	    break; // stop all
	  switch (type) {
	  case DUALSIMPLEX:
	  case PRIMALSIMPLEX:
	  case EITHERSIMPLEX:
	  case BARRIER:
            // synonym for dual
	  case BAB:
	    if (goodModels[iModel]) {
              double objScale = 
                parameters[whichParam(OBJSCALE2,numberParameters,parameters)].doubleValue();
              if (objScale!=1.0) {
                int iColumn;
                int numberColumns=models[iModel].numberColumns();
                double * dualColumnSolution = 
                  models[iModel].dualColumnSolution();
                ClpObjective * obj = models[iModel].objectiveAsObject();
                assert(dynamic_cast<ClpLinearObjective *> (obj));
                double offset;
                double * objective = obj->gradient(NULL,NULL,offset,true);
                for (iColumn=0;iColumn<numberColumns;iColumn++) {
                  dualColumnSolution[iColumn] *= objScale;
                  objective[iColumn] *= objScale;;
                }
                int iRow;
                int numberRows=models[iModel].numberRows();
                double * dualRowSolution = 
                  models[iModel].dualRowSolution();
                for (iRow=0;iRow<numberRows;iRow++) 
                  dualRowSolution[iRow] *= objScale;
                models[iModel].setObjectiveOffset(objScale*models[iModel].objectiveOffset());
              }
	      ClpSolve::SolveType method;
	      ClpSolve::PresolveType presolveType;
	      ClpSimplex * model2 = models+iModel;
              if (dualize) {
		bool tryIt=true;
		double fractionColumn=1.0;
		double fractionRow=1.0;
		if (dualize==3) {
		  dualize=1;
		  int numberColumns=model2->numberColumns();
		  int numberRows=model2->numberRows();
		  if (numberRows<50000||5*numberColumns>numberRows) {
		    tryIt=false;
		  } else {
		    fractionColumn=0.1;
		    fractionRow=0.1;
		  }
		}
		if (tryIt) {
		  model2 = static_cast<ClpSimplexOther *> (model2)->dualOfModel(fractionRow,fractionColumn);
		  if (model2) {
		    printf("Dual of model has %d rows and %d columns\n",
			   model2->numberRows(),model2->numberColumns());
		    model2->setOptimizationDirection(1.0);
		  } else {
		    model2 = models+iModel;
		    dualize=0;
		  }
		} else {
		  dualize=0;
		}
              }
	      ClpSolve solveOptions;
              solveOptions.setPresolveActions(presolveOptions);
              solveOptions.setSubstitution(substitution);
	      if (preSolve!=5&&preSolve) {
		presolveType=ClpSolve::presolveNumber;
                if (preSolve<0) {
                  preSolve = - preSolve;
                  if (preSolve<=100) {
                    presolveType=ClpSolve::presolveNumber;
                    printf("Doing %d presolve passes - picking up non-costed slacks\n",
                           preSolve);
                    solveOptions.setDoSingletonColumn(true);
                  } else {
                    preSolve -=100;
                    presolveType=ClpSolve::presolveNumberCost;
                    printf("Doing %d presolve passes - picking up costed slacks\n",
                           preSolve);
                  }
                } 
	      } else if (preSolve) {
		presolveType=ClpSolve::presolveOn;
	      } else {
		presolveType=ClpSolve::presolveOff;
              }
	      solveOptions.setPresolveType(presolveType,preSolve);
	      if (type==DUALSIMPLEX||type==BAB) {
		method=ClpSolve::useDual;
	      } else if (type==PRIMALSIMPLEX) {
		method=ClpSolve::usePrimalorSprint;
	      } else if (type==EITHERSIMPLEX) {
		method=ClpSolve::automatic;
	      } else {
		method = ClpSolve::useBarrier;
		if (crossover==1) {
		  method=ClpSolve::useBarrierNoCross;
		} else if (crossover==2) {
		  ClpObjective * obj = models[iModel].objectiveAsObject();
		  if (obj->type()>1) {
		    method=ClpSolve::useBarrierNoCross;
		    presolveType=ClpSolve::presolveOff;
		    solveOptions.setPresolveType(presolveType,preSolve);
		  } 
		}
	      }
	      solveOptions.setSolveType(method);
	      if(preSolveFile)
		presolveOptions |= 0x40000000;
	      solveOptions.setSpecialOption(4,presolveOptions);
	      solveOptions.setSpecialOption(5,printOptions&1);
	      if (doVector) {
		ClpMatrixBase * matrix = models[iModel].clpMatrix();
		if (dynamic_cast< ClpPackedMatrix*>(matrix)) {
		  ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
		  clpMatrix->makeSpecialColumnCopy();
		}
	      }
	      if (method==ClpSolve::useDual) {
		// dual
		if (doCrash)
		  solveOptions.setSpecialOption(0,1,doCrash); // crash
		else if (doIdiot)
		  solveOptions.setSpecialOption(0,2,doIdiot); // possible idiot
	      } else if (method==ClpSolve::usePrimalorSprint) {
		// primal
		// if slp turn everything off
		if (slpValue>0) {
		  doCrash=false;
		  doSprint=0;
		  doIdiot=-1;
		  solveOptions.setSpecialOption(1,10,slpValue); // slp
		  method=ClpSolve::usePrimal;
		}
		if (doCrash) {
		  solveOptions.setSpecialOption(1,1,doCrash); // crash
		} else if (doSprint>0) {
		  // sprint overrides idiot
		  solveOptions.setSpecialOption(1,3,doSprint); // sprint
		} else if (doIdiot>0) {
		  solveOptions.setSpecialOption(1,2,doIdiot); // idiot
		} else if (slpValue<=0) {
		  if (doIdiot==0) {
		    if (doSprint==0)
		      solveOptions.setSpecialOption(1,4); // all slack
		    else
		      solveOptions.setSpecialOption(1,9); // all slack or sprint
		  } else {
		    if (doSprint==0)
		      solveOptions.setSpecialOption(1,8); // all slack or idiot
		    else
		      solveOptions.setSpecialOption(1,7); // initiative
		  }
		}
		if (basisHasValues==-1)
		  solveOptions.setSpecialOption(1,11); // switch off values
	      } else if (method==ClpSolve::useBarrier||method==ClpSolve::useBarrierNoCross) {
		int barrierOptions = choleskyType;
		if (scaleBarrier)
		  barrierOptions |= 8;
		if (doKKT)
		  barrierOptions |= 16;
		if (gamma)
		  barrierOptions |= 32*gamma;
		if (crossover==3) 
		  barrierOptions |= 256; // try presolve in crossover
		solveOptions.setSpecialOption(4,barrierOptions);
	      }
	      int status;
	      if (cppValue>=0) {
                // generate code
                FILE * fp = fopen("user_driver.cpp","w");
	        if (fp) {
	          // generate enough to do solveOptions
                  model2->generateCpp(fp);
                  solveOptions.generateCpp(fp);
                  fclose(fp);
                  // now call generate code
                  generateCode("user_driver.cpp",cppValue);
                } else {
                  std::cout<<"Unable to open file user_driver.cpp"<<std::endl;
                }
              }
#ifdef CLP_MULTIPLE_FACTORIZATIONS   
	      int denseCode = parameters[whichParam(DENSE,numberParameters,parameters)].intValue();
	      model2->factorization()->setGoDenseThreshold(denseCode);
	      int smallCode = parameters[whichParam(SMALLFACT,numberParameters,parameters)].intValue();
	      model2->factorization()->setGoSmallThreshold(smallCode);
	      model2->factorization()->goDenseOrSmall(model2->numberRows());
#endif
              try {
                status=model2->initialSolve(solveOptions);
              }
              catch (CoinError e) {
                e.print();
                status=-1;
              }
              if (dualize) {
                int returnCode=static_cast<ClpSimplexOther *> (models+iModel)->restoreFromDual(model2);
		if (model2->status()==3)
		  returnCode=0;
                delete model2;
		if (returnCode&&dualize!=2) {
		  currentModel = models+iModel;
		  // register signal handler
		  signal(SIGINT,signal_handler);
		  models[iModel].primal(1);
		  currentModel=NULL;
		}
              }
              if (status>=0)
                basisHasValues=1;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
          case STATISTICS:
	    if (goodModels[iModel]) {
              // If presolve on look at presolved
              bool deleteModel2=false;
              ClpSimplex * model2 = models+iModel;
              if (preSolve) {
                ClpPresolve pinfo;
                int presolveOptions2 = presolveOptions&~0x40000000;
                if ((presolveOptions2&0xffff)!=0)
                  pinfo.setPresolveActions(presolveOptions2);
                pinfo.setSubstitution(substitution);
                if ((printOptions&1)!=0)
                  pinfo.statistics();
                double presolveTolerance = 
                  parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
                model2 = 
                  pinfo.presolvedModel(models[iModel],presolveTolerance,
                                       true,preSolve);
                if (model2) {
                  printf("Statistics for presolved model\n");
                  deleteModel2=true;
                } else {
                  printf("Presolved model looks infeasible - will use unpresolved\n");
                  model2 = models+iModel;
                }
              } else {
                printf("Statistics for unpresolved model\n");
                model2 =  models+iModel;
              }
              statistics(models+iModel,model2);
              if (deleteModel2)
                delete model2;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case TIGHTEN:
	    if (goodModels[iModel]) {
     	      int numberInfeasibilities = models[iModel].tightenPrimalBounds();
	      if (numberInfeasibilities)
		std::cout<<"** Analysis indicates model infeasible"<<std::endl;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case PLUSMINUS:
	    if (goodModels[iModel]) {
	      ClpMatrixBase * saveMatrix = models[iModel].clpMatrix();
	      ClpPackedMatrix* clpMatrix =
		dynamic_cast< ClpPackedMatrix*>(saveMatrix);
	      if (clpMatrix) {
		ClpPlusMinusOneMatrix * newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
		if (newMatrix->getIndices()) {
		  models[iModel].replaceMatrix(newMatrix);
		  delete saveMatrix;
		  std::cout<<"Matrix converted to +- one matrix"<<std::endl;
		} else {
		  std::cout<<"Matrix can not be converted to +- 1 matrix"<<std::endl;
		}
	      } else {
		std::cout<<"Matrix not a ClpPackedMatrix"<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case NETWORK:
	    if (goodModels[iModel]) {
	      ClpMatrixBase * saveMatrix = models[iModel].clpMatrix();
	      ClpPackedMatrix* clpMatrix =
		dynamic_cast< ClpPackedMatrix*>(saveMatrix);
	      if (clpMatrix) {
		ClpNetworkMatrix * newMatrix = new ClpNetworkMatrix(*(clpMatrix->matrix()));
		if (newMatrix->getIndices()) {
		  models[iModel].replaceMatrix(newMatrix);
		  delete saveMatrix;
		  std::cout<<"Matrix converted to network matrix"<<std::endl;
		} else {
		  std::cout<<"Matrix can not be converted to network matrix"<<std::endl;
		}
	      } else {
		std::cout<<"Matrix not a ClpPackedMatrix"<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case IMPORT:
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
              // See if gmpl file
              int gmpl=0;
              std::string gmplData;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		// See if .lp
		{
		  const char * c_name = field.c_str();
		  int length = strlen(c_name);
		  if (length>3&&!strncmp(c_name+length-3,".lp",3))
		    gmpl=-1; // .lp
		}
                bool absolutePath;
                if (dirsep=='/') {
                  // non Windows (or cygwin)
                  absolutePath=(field[0]=='/');
                } else {
                  //Windows (non cycgwin)
                  absolutePath=(field[0]=='\\');
                  // but allow for :
                  if (strchr(field.c_str(),':'))
                    absolutePath=true;
                }
		if (absolutePath) {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
                  // See if gmpl (model & data) - or even lp file
                  int length = field.size();
                  int percent = field.find('%');
                  if (percent<length&&percent>0) {
                    gmpl=1;
                    fileName = directory+field.substr(0,percent);
                    gmplData = directory+field.substr(percent+1);
                    if (percent<length-1)
                      gmpl=2; // two files
                    printf("GMPL model file %s and data file %s\n",
                           fileName.c_str(),gmplData.c_str());
		  }
		}
                std::string name=fileName;
                if (fileCoinReadable(name)) {
		  // can open - lets go for it
		  canOpen=true;
                  if (gmpl==2) {
                    FILE *fp;
                    fp=fopen(gmplData.c_str(),"r");
                    if (fp) {
                      fclose(fp);
                    } else {
                      canOpen=false;
                      std::cout<<"Unable to open file "<<gmplData<<std::endl;
                    }
                  }
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
                int status;
                if (!gmpl)
                  status =models[iModel].readMps(fileName.c_str(),
                                                 keepImportNames!=0,
                                                 allowImportErrors!=0);
                else if (gmpl>0)
                  status= models[iModel].readGMPL(fileName.c_str(),
                                                  (gmpl==2) ? gmplData.c_str() : NULL,
                                                  keepImportNames!=0);
		else
                  status= models[iModel].readLp(fileName.c_str(),1.0e-12);
		if (!status||(status>0&&allowImportErrors)) {
		  goodModels[iModel]=true;
		  // sets to all slack (not necessary?)
		  models[iModel].createStatus();
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		  // Go to canned file if just input file
		  if (CbcOrClpRead_mode==2&&argc==2) {
		    // only if ends .mps
		    char * find = reinterpret_cast<char *>(strstr(fileName.c_str(),".mps"));
		    if (find&&find[4]=='\0') {
		      find[1]='p'; find[2]='a';find[3]='r';
		      FILE *fp=fopen(fileName.c_str(),"r");
		      if (fp) {
			CbcOrClpReadCommand=fp; // Read from that file
			CbcOrClpRead_mode=-1;
		      }
		    }
		  }
		} else {
		  // errors
		  std::cout<<"There were "<<status<<
		    " errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case EXPORT:
	    if (goodModels[iModel]) {
              double objScale = 
                parameters[whichParam(OBJSCALE2,numberParameters,parameters)].doubleValue();
              if (objScale!=1.0) {
                int iColumn;
                int numberColumns=models[iModel].numberColumns();
                double * dualColumnSolution = 
                  models[iModel].dualColumnSolution();
                ClpObjective * obj = models[iModel].objectiveAsObject();
                assert(dynamic_cast<ClpLinearObjective *> (obj));
                double offset;
                double * objective = obj->gradient(NULL,NULL,offset,true);
                for (iColumn=0;iColumn<numberColumns;iColumn++) {
                  dualColumnSolution[iColumn] *= objScale;
                  objective[iColumn] *= objScale;;
                }
                int iRow;
                int numberRows=models[iModel].numberRows();
                double * dualRowSolution = 
                  models[iModel].dualRowSolution();
                for (iRow=0;iRow<numberRows;iRow++) 
                  dualRowSolution[iRow] *= objScale;
                models[iModel].setObjectiveOffset(objScale*models[iModel].objectiveOffset());
              }
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		// If presolve on then save presolved
		bool deleteModel2=false;
		ClpSimplex * model2 = models+iModel;
		if (dualize&&dualize<3) {
		  model2 = static_cast<ClpSimplexOther *> (model2)->dualOfModel();
		  printf("Dual of model has %d rows and %d columns\n",
			 model2->numberRows(),model2->numberColumns());
		  model2->setOptimizationDirection(1.0);
		  preSolve=0; // as picks up from model
		}
		if (preSolve) {
		  ClpPresolve pinfo;
		  int presolveOptions2 = presolveOptions&~0x40000000;
		  if ((presolveOptions2&0xffff)!=0)
		    pinfo.setPresolveActions(presolveOptions2);
                  pinfo.setSubstitution(substitution);
		  if ((printOptions&1)!=0)
		    pinfo.statistics();
                  double presolveTolerance = 
                    parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
		  model2 = 
		    pinfo.presolvedModel(models[iModel],presolveTolerance,
					 true,preSolve,false,false);
		  if (model2) {
		    printf("Saving presolved model on %s\n",
			   fileName.c_str());
		    deleteModel2=true;
		  } else {
		    printf("Presolved model looks infeasible - saving original on %s\n",
			   fileName.c_str());
		    deleteModel2=false;
		    model2 = models+iModel;

		  }
		} else {
		  printf("Saving model on %s\n",
			   fileName.c_str());
		}
#if 0
		// Convert names
		int iRow;
		int numberRows=model2->numberRows();
		int iColumn;
		int numberColumns=model2->numberColumns();

		char ** rowNames = NULL;
		char ** columnNames = NULL;
		if (model2->lengthNames()) {
		  rowNames = new char * [numberRows];
		  for (iRow=0;iRow<numberRows;iRow++) {
		    rowNames[iRow] = 
		      CoinStrdup(model2->rowName(iRow).c_str());
#ifdef STRIPBLANKS
		    char * xx = rowNames[iRow];
		    int i;
		    int length = strlen(xx);
		    int n=0;
		    for (i=0;i<length;i++) {
		      if (xx[i]!=' ')
			xx[n++]=xx[i];
		    }
		    xx[n]='\0';
#endif
		  }
		  
		  columnNames = new char * [numberColumns];
		  for (iColumn=0;iColumn<numberColumns;iColumn++) {
		    columnNames[iColumn] = 
		      CoinStrdup(model2->columnName(iColumn).c_str());
#ifdef STRIPBLANKS
		    char * xx = columnNames[iColumn];
		    int i;
		    int length = strlen(xx);
		    int n=0;
		    for (i=0;i<length;i++) {
		      if (xx[i]!=' ')
			xx[n++]=xx[i];
		    }
		    xx[n]='\0';
#endif
		  }
		}
		CoinMpsIO writer;
		writer.setMpsData(*model2->matrix(), COIN_DBL_MAX,
				  model2->getColLower(), model2->getColUpper(),
				  model2->getObjCoefficients(),
				  (const char*) 0 /*integrality*/,
				  model2->getRowLower(), model2->getRowUpper(),
				  columnNames, rowNames);
		// Pass in array saying if each variable integer
		writer.copyInIntegerInformation(model2->integerInformation());
		writer.setObjectiveOffset(model2->objectiveOffset());
		writer.writeMps(fileName.c_str(),0,1,1);
		if (rowNames) {
		  for (iRow=0;iRow<numberRows;iRow++) {
		    free(rowNames[iRow]);
		  }
		  delete [] rowNames;
		  for (iColumn=0;iColumn<numberColumns;iColumn++) {
		    free(columnNames[iColumn]);
		  }
		  delete [] columnNames;
		}
#else
		model2->writeMps(fileName.c_str(),(outputFormat-1)/2,1+((outputFormat-1)&1));
#endif
		if (deleteModel2)
		  delete model2;
		time2 = CoinCpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case BASISIN:
	    if (goodModels[iModel]) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
		}
		FILE *fp=fopen(fileName.c_str(),"r");
		if (fp) {
		  // can open - lets go for it
		  fclose(fp);
		  canOpen=true;
		} else {
		  std::cout<<"Unable to open file "<<fileName<<std::endl;
		}
	      }
	      if (canOpen) {
		int values = models[iModel].readBasis(fileName.c_str());
		if (values==0)
		  basisHasValues=-1;
		else
		  basisHasValues=1;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case PRINTMASK:
            // get next field
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		parameters[iParam].setStringValue(name);
                printMask = name;
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case BASISOUT:
	    if (goodModels[iModel]) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		ClpSimplex * model2 = models+iModel;
		model2->writeBasis(fileName.c_str(),outputFormat>1,outputFormat-2);
		time2 = CoinCpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
	    break;
	  case SAVE:
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"wb");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		int status;
		// If presolve on then save presolved
		bool deleteModel2=false;
		ClpSimplex * model2 = models+iModel;
		if (preSolve) {
		  ClpPresolve pinfo;
                  double presolveTolerance = 
                    parameters[whichParam(PRESOLVETOLERANCE,numberParameters,parameters)].doubleValue();
		  model2 = 
		    pinfo.presolvedModel(models[iModel],presolveTolerance,
					 false,preSolve);
		  if (model2) {
		    printf("Saving presolved model on %s\n",
			   fileName.c_str());
		    deleteModel2=true;
		  } else {
		    printf("Presolved model looks infeasible - saving original on %s\n",
			   fileName.c_str());
		    deleteModel2=false;
		    model2 = models+iModel;

		  }
		} else {
		  printf("Saving model on %s\n",
			   fileName.c_str());
		}
		status =model2->saveModel(fileName.c_str());
		if (deleteModel2)
		  delete model2;
		if (!status) {
		  goodModels[iModel]=true;
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were errors on output"<<std::endl;
		}
	      }
	    }
	    break;
	  case RESTORE:
	    {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='\\') {
		fileName = field;
	      } else if (field[0]=='~') {
		char * environVar = getenv("HOME");
		if (environVar) {
		  std::string home(environVar);
		  field=field.erase(0,1);
		  fileName = home+field;
		} else {
		  fileName=field;
		}
	      } else {
		fileName = directory+field;
	      }
	      FILE *fp=fopen(fileName.c_str(),"rb");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		int status =models[iModel].restoreModel(fileName.c_str());
		if (!status) {
		  goodModels[iModel]=true;
		  time2 = CoinCpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case MAXIMIZE:
	    models[iModel].setOptimizationDirection(-1);
	    break;
	  case MINIMIZE:
	    models[iModel].setOptimizationDirection(1);
	    break;
	  case ALLSLACK:
	    models[iModel].allSlackBasis(true);
	    break;
	  case REVERSE:
	    if (goodModels[iModel]) {
	      int iColumn;
	      int numberColumns=models[iModel].numberColumns();
	      double * dualColumnSolution = 
		models[iModel].dualColumnSolution();
	      ClpObjective * obj = models[iModel].objectiveAsObject();
	      assert(dynamic_cast<ClpLinearObjective *> (obj));
	      double offset;
	      double * objective = obj->gradient(NULL,NULL,offset,true);
	      for (iColumn=0;iColumn<numberColumns;iColumn++) {
		dualColumnSolution[iColumn] = -dualColumnSolution[iColumn];
		objective[iColumn] = -objective[iColumn];
	      }
	      int iRow;
	      int numberRows=models[iModel].numberRows();
	      double * dualRowSolution = 
		models[iModel].dualRowSolution();
	      for (iRow=0;iRow<numberRows;iRow++) {
		dualRowSolution[iRow] = -dualRowSolution[iRow];
	      }
              models[iModel].setObjectiveOffset(-models[iModel].objectiveOffset());
	    }
	    break;
	  case DIRECTORY:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]==dirsep) {
		  directory = name;
		} else {
		  directory = name+dirsep;
		}
		parameters[iParam].setStringValue(directory);
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case DIRSAMPLE:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]==dirsep) {
		  dirSample = name;
		} else {
		  dirSample = name+dirsep;
		}
		parameters[iParam].setStringValue(dirSample);
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case DIRNETLIB:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]==dirsep) {
		  dirNetlib = name;
		} else {
		  dirNetlib = name+dirsep;
		}
		parameters[iParam].setStringValue(dirNetlib);
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case DIRMIPLIB:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]==dirsep) {
		  dirMiplib = name;
		} else {
		  dirMiplib = name+dirsep;
		}
		parameters[iParam].setStringValue(dirMiplib);
	      } else {
		parameters[iParam].printString();
	      }
	    }
	    break;
	  case STDIN:
	    CbcOrClpRead_mode=-1;
	    break;
	  case NETLIB_DUAL:
	  case NETLIB_EITHER:
	  case NETLIB_BARRIER:
	  case NETLIB_PRIMAL:
	  case NETLIB_TUNE:
	    {
	      // create fields for unitTest
	      const char * fields[4];
	      int nFields=4;
	      fields[0]="fake main from unitTest";
	      std::string mpsfield = "-dirSample=";
	      mpsfield += dirSample.c_str();
	      fields[1]=mpsfield.c_str();
	      std::string netfield = "-dirNetlib=";
	      netfield += dirNetlib.c_str();
	      fields[2]=netfield.c_str();
	      fields[3]="-netlib";
	      int algorithm;
	      if (type==NETLIB_DUAL) {
		std::cerr<<"Doing netlib with dual algorithm"<<std::endl;
		algorithm =0;
	      } else if (type==NETLIB_BARRIER) {
		std::cerr<<"Doing netlib with barrier algorithm"<<std::endl;
		algorithm =2;
	      } else if (type==NETLIB_EITHER) {
		std::cerr<<"Doing netlib with dual or primal algorithm"<<std::endl;
		algorithm =3;
	      } else if (type==NETLIB_TUNE) {
		std::cerr<<"Doing netlib with best algorithm!"<<std::endl;
		algorithm =5;
                // uncomment next to get active tuning
                // algorithm=6;
	      } else {
		std::cerr<<"Doing netlib with primal agorithm"<<std::endl;
		algorithm=1;
	      }
              int specialOptions = models[iModel].specialOptions();
              models[iModel].setSpecialOptions(0);
	      mainTest(nFields,fields,algorithm,models[iModel],
		       (preSolve!=0),specialOptions,doVector!=0);
	    }
	    break;
	  case UNITTEST:
	    {
	      // create fields for unitTest
	      const char * fields[2];
	      int nFields=2;
	      fields[0]="fake main from unitTest";
	      std::string dirfield = "-dirSample=";
	      dirfield += dirSample.c_str();
	      fields[1]=dirfield.c_str();
              int specialOptions = models[iModel].specialOptions();
              models[iModel].setSpecialOptions(0);
              int algorithm=-1;
              if (models[iModel].numberRows())
                algorithm=7;
	      mainTest(nFields,fields,algorithm,models[iModel],(preSolve!=0),specialOptions,doVector!=0);
	    }
	    break;
	  case FAKEBOUND:
	    if (goodModels[iModel]) {
	      // get bound
	      double value = CoinReadGetDoubleField(argc,argv,&valid);
	      if (!valid) {
		std::cout<<"Setting "<<parameters[iParam].name()<<
		  " to DEBUG "<<value<<std::endl;
		int iRow;
		int numberRows=models[iModel].numberRows();
		double * rowLower = models[iModel].rowLower();
		double * rowUpper = models[iModel].rowUpper();
		for (iRow=0;iRow<numberRows;iRow++) {
		  // leave free ones for now
		  if (rowLower[iRow]>-1.0e20||rowUpper[iRow]<1.0e20) {
		    rowLower[iRow]=CoinMax(rowLower[iRow],-value);
		    rowUpper[iRow]=CoinMin(rowUpper[iRow],value);
		  }
		}
		int iColumn;
		int numberColumns=models[iModel].numberColumns();
		double * columnLower = models[iModel].columnLower();
		double * columnUpper = models[iModel].columnUpper();
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  // leave free ones for now
		  if (columnLower[iColumn]>-1.0e20||
		      columnUpper[iColumn]<1.0e20) {
		    columnLower[iColumn]=CoinMax(columnLower[iColumn],-value);
		    columnUpper[iColumn]=CoinMin(columnUpper[iColumn],value);
		  }
		}
	      } else if (valid==1) {
		abort();
	      } else {
		std::cout<<"enter value for "<<parameters[iParam].name()<<
		  std::endl;
	      }
	    }
	    break;
	  case REALLY_SCALE:
	    if (goodModels[iModel]) {
	      ClpSimplex newModel(models[iModel],
				  models[iModel].scalingFlag());
	      printf("model really really scaled\n");
	      models[iModel]=newModel;
	    }
	    break;
	  case USERCLP:
            // Replace the sample code by whatever you want
	    if (goodModels[iModel]) {
              ClpSimplex * thisModel = &models[iModel];
              printf("Dummy user code - model has %d rows and %d columns\n",
                     thisModel->numberRows(),thisModel->numberColumns());
	    }
	    break;
	  case HELP:
	    std::cout<<"Coin LP version "<<CLPVERSION
		     <<", build "<<__DATE__<<std::endl;
	    std::cout<<"Non default values:-"<<std::endl;
	    std::cout<<"Perturbation "<<models[0].perturbation()<<" (default 100)"
		     <<std::endl;
	    CoinReadPrintit(
		    "Presolve being done with 5 passes\n\
Dual steepest edge steep/partial on matrix shape and factorization density\n\
Clpnnnn taken out of messages\n\
If Factorization frequency default then done on size of matrix\n\n\
(-)unitTest, (-)netlib or (-)netlibp will do standard tests\n\n\
You can switch to interactive mode at any time so\n\
clp watson.mps -scaling off -primalsimplex\nis the same as\n\
clp watson.mps -\nscaling off\nprimalsimplex"
		    );
  	    break;
	  case SOLUTION:
	    if (goodModels[iModel]) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
	      FILE *fp=NULL;
	      if (field=="-"||field=="EOL"||field=="stdout") {
		// stdout
		fp=stdout;
		fprintf(fp,"\n");
	      } else if (field=="stderr") {
		// stderr
		fp=stderr;
		fprintf(fp,"\n");
	      } else {
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environVar = getenv("HOME");
		  if (environVar) {
		    std::string home(environVar);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
		}
		fp=fopen(fileName.c_str(),"w");
	      }
	      if (fp) {
		// Write solution header (suggested by Luigi Poderico)
		double objValue = models[iModel].getObjValue()*models[iModel].getObjSense();
		int iStat = models[iModel].status();
		if (iStat==0) {
		  fprintf(fp, "optimal\n" );
		} else if (iStat==1) {
		  // infeasible
		  fprintf(fp, "infeasible\n" );
		} else if (iStat==2) {
		  // unbounded
		  fprintf(fp, "unbounded\n" );
		} else if (iStat==3) {
		  fprintf(fp, "stopped on iterations or time\n" );
		} else if (iStat==4) {
		  fprintf(fp, "stopped on difficulties\n" );
		} else if (iStat==5) {
		  fprintf(fp, "stopped on ctrl-c\n" );
		} else {
		  fprintf(fp, "status unknown\n" );
		}
		fprintf(fp, "Objective value %15.8g\n", objValue);
		// make fancy later on
		int iRow;
		int numberRows=models[iModel].numberRows();
		int lengthName = models[iModel].lengthNames(); // 0 if no names
		// in general I don't want to pass around massive
		// amounts of data but seems simpler here
		std::vector<std::string> rowNames =
		  *(models[iModel].rowNames());
		std::vector<std::string> columnNames =
		  *(models[iModel].columnNames());

		double * dualRowSolution = models[iModel].dualRowSolution();
		double * primalRowSolution = 
		  models[iModel].primalRowSolution();
		double * rowLower = models[iModel].rowLower();
		double * rowUpper = models[iModel].rowUpper();
		double primalTolerance = models[iModel].primalTolerance();
		char format[6];
		sprintf(format,"%%-%ds",CoinMax(lengthName,8));
                bool doMask = (printMask!=""&&lengthName);
		int * maskStarts=NULL;
		int maxMasks=0;
		char ** masks =NULL;
		if (doMask) {
		  int nAst =0;
		  const char * pMask2 = printMask.c_str();
		  char pMask[100];
		  int iChar;
		  int lengthMask = strlen(pMask2);
		  assert (lengthMask<100);
		  if (*pMask2=='"') {
		    if (pMask2[lengthMask-1]!='"') {
		      printf("mismatched \" in mask %s\n",pMask2);
		      break;
		    } else {
		      strcpy(pMask,pMask2+1);
		      *strchr(pMask,'"')='\0';
		    }
		  } else if (*pMask2=='\'') {
		    if (pMask2[lengthMask-1]!='\'') {
		      printf("mismatched ' in mask %s\n",pMask2);
		      break;
		    } else {
		      strcpy(pMask,pMask2+1);
		      *strchr(pMask,'\'')='\0';
		    }
		  } else {
		    strcpy(pMask,pMask2);
		  }
		  if (lengthMask>lengthName) {
		    printf("mask %s too long - skipping\n",pMask);
		    break;
		  }
		  maxMasks = 1;
		  for (iChar=0;iChar<lengthMask;iChar++) {
		    if (pMask[iChar]=='*') {
		      nAst++;
		      maxMasks *= (lengthName+1);
		    }
		  }
		  int nEntries = 1;
		  maskStarts = new int[lengthName+2];
		  masks = new char * [maxMasks];
		  char ** newMasks = new char * [maxMasks];
		  int i;
		  for (i=0;i<maxMasks;i++) {
		    masks[i] = new char[lengthName+1];
		    newMasks[i] = new char[lengthName+1];
		  }
		  strcpy(masks[0],pMask);
		  for (int iAst=0;iAst<nAst;iAst++) {
		    int nOldEntries = nEntries;
		    nEntries=0;
		    for (int iEntry = 0;iEntry<nOldEntries;iEntry++) {
		      char * oldMask = masks[iEntry];
		      char * ast = strchr(oldMask,'*');
		      assert (ast);
		      int length = strlen(oldMask)-1;
		      int nBefore = ast-oldMask;
		      int nAfter = length-nBefore;
		      // and add null
		      nAfter++;
		      for (int i=0;i<=lengthName-length;i++) {
			char * maskOut = newMasks[nEntries];
   CoinMemcpyN(oldMask,nBefore,maskOut);
			for (int k=0;k<i;k++) 
			  maskOut[k+nBefore]='?';
   CoinMemcpyN(ast+1,nAfter,maskOut+nBefore+i);
			nEntries++;
			assert (nEntries<=maxMasks);
		      }
		    }
		    char ** temp = masks;
		    masks = newMasks;
		    newMasks = temp;
		  }
		  // Now extend and sort
		  int * sort = new int[nEntries];
		  for (i=0;i<nEntries;i++) {
		    char * maskThis = masks[i];
		    int length = strlen(maskThis);
		    while (maskThis[length-1]==' ')
		      length--;
		    maskThis[length]='\0';
		    sort[i]=length;
		  }
		  CoinSort_2(sort,sort+nEntries,masks);
		  int lastLength=-1;
		  for (i=0;i<nEntries;i++) {
		    int length = sort[i];
		    while (length>lastLength) 
		      maskStarts[++lastLength] = i;
		  }
		  maskStarts[++lastLength]=nEntries;
		  delete [] sort;
		  for (i=0;i<maxMasks;i++)
		    delete [] newMasks[i];
		  delete [] newMasks;
		}
                if (printMode>2) {
                  for (iRow=0;iRow<numberRows;iRow++) {
                    int type=printMode-3;
                    if (primalRowSolution[iRow]>rowUpper[iRow]+primalTolerance||
                        primalRowSolution[iRow]<rowLower[iRow]-primalTolerance) {
                      fprintf(fp,"** ");
                      type=2;
                    } else if (fabs(primalRowSolution[iRow])>1.0e-8) {
                      type=1;
                    } else if (numberRows<50) {
                      type=3;
                    } 
                    if (doMask&&!maskMatches(maskStarts,masks,rowNames[iRow]))
                      type=0;
                    if (type) {
                      fprintf(fp,"%7d ",iRow);
                      if (lengthName)
                        fprintf(fp,format,rowNames[iRow].c_str());
                      fprintf(fp,"%15.8g        %15.8g\n",primalRowSolution[iRow],
                              dualRowSolution[iRow]);
                    }
                  }
                }
		int iColumn;
		int numberColumns=models[iModel].numberColumns();
		double * dualColumnSolution = 
  models[iModel].dualColumnSolution();
		double * primalColumnSolution = 
  models[iModel].primalColumnSolution();
		double * columnLower = models[iModel].columnLower();
		double * columnUpper = models[iModel].columnUpper();
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  int type=(printMode>3) ? 1 : 0;
		  if (primalColumnSolution[iColumn]>columnUpper[iColumn]+primalTolerance||
		      primalColumnSolution[iColumn]<columnLower[iColumn]-primalTolerance) {
		    fprintf(fp,"** ");
		    type=2;
		  } else if (fabs(primalColumnSolution[iColumn])>1.0e-8) {
		    type=1;
		  } else if (numberColumns<50) {
		    type=3;
		  }
		  if (doMask&&!maskMatches(maskStarts,masks,
					   columnNames[iColumn]))
                    type =0;
		  if (type) {
		    fprintf(fp,"%7d ",iColumn);
		    if (lengthName)
		      fprintf(fp,format,columnNames[iColumn].c_str());
		    fprintf(fp,"%15.8g        %15.8g\n",
			    primalColumnSolution[iColumn],
			    dualColumnSolution[iColumn]);
		  }
		}
		if (fp!=stdout)
		  fclose(fp);
		if (masks) {
		  delete [] maskStarts;
		  for (int i=0;i<maxMasks;i++)
		    delete [] masks[i];
		  delete [] masks;
		}
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	      
	    }
	  
	    break;
	  case SAVESOL:
	    if (goodModels[iModel]) {
	      // get next field
	      field = CoinReadGetString(argc,argv);
	      if (field=="$") {
		field = parameters[iParam].stringValue();
	      } else if (field=="EOL") {
		parameters[iParam].printString();
		break;
	      } else {
		parameters[iParam].setStringValue(field);
	      }
	      std::string fileName;
              if (field[0]=='/'||field[0]=='\\') {
                fileName = field;
              } else if (field[0]=='~') {
                char * environVar = getenv("HOME");
                if (environVar) {
                  std::string home(environVar);
                  field=field.erase(0,1);
                  fileName = home+field;
                } else {
                  fileName=field;
                }
              } else {
                fileName = directory+field;
              }
              saveSolution(models+iModel,fileName);
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	      
	    }
	    break;
	  default:
	    abort();
	  }
	} 
      } else if (!numberMatches) {
	std::cout<<"No match for "<<field<<" - ? for list of commands"
		 <<std::endl;
      } else if (numberMatches==1) {
	if (!numberQuery) {
	  std::cout<<"Short match for "<<field<<" - completion: ";
	  std::cout<<parameters[firstMatch].matchName()<<std::endl;
	} else if (numberQuery) {
	  std::cout<<parameters[firstMatch].matchName()<<" : ";
	  std::cout<<parameters[firstMatch].shortHelp()<<std::endl;
	  if (numberQuery>=2) 
	    parameters[firstMatch].printLongHelp();
	}
      } else {
	if (!numberQuery) 
	  std::cout<<"Multiple matches for "<<field<<" - possible completions:"
		   <<std::endl;
	else
	  std::cout<<"Completions of "<<field<<":"<<std::endl;
	for ( iParam=0; iParam<numberParameters; iParam++ ) {
	  int match = parameters[iParam].matches(field);
	  if (match&&parameters[iParam].displayThis()) {
	    std::cout<<parameters[iParam].matchName();
	    if (numberQuery>=2) 
	      std::cout<<" : "<<parameters[iParam].shortHelp();
	    std::cout<<std::endl;
	  }
	}
      }
    }
    delete [] models;
    delete [] goodModels;
  }
  // By now all memory should be freed
#ifdef DMALLOC
  dmalloc_log_unfreed();
  dmalloc_shutdown();
#endif
  return 0;
}    
static void breakdown(const char * name, int numberLook, const double * region)
{
  double range[] = {
    -COIN_DBL_MAX,
    -1.0e15,-1.0e11,-1.0e8,-1.0e5,-1.0e4,-1.0e3,-1.0e2,-1.0e1,
    -1.0,
    -1.0e-1,-1.0e-2,-1.0e-3,-1.0e-4,-1.0e-5,-1.0e-8,-1.0e-11,-1.0e-15,
    0.0,
    1.0e-15,1.0e-11,1.0e-8,1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,
    1.0,
    1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,1.0e8,1.0e11,1.0e15,
    COIN_DBL_MAX};
  int nRanges = static_cast<int> (sizeof(range)/sizeof(double));
  int * number = new int[nRanges];
  memset(number,0,nRanges*sizeof(int));
  int * numberExact = new int[nRanges];
  memset(numberExact,0,nRanges*sizeof(int));
  int i;
  for ( i=0;i<numberLook;i++) {
    double value = region[i];
    for (int j=0;j<nRanges;j++) {
      if (value==range[j]) {
        numberExact[j]++;
        break;
      } else if (value<range[j]) {
        number[j]++;
        break;
      }
    }
  }
  printf("\n%s has %d entries\n",name,numberLook);
  for (i=0;i<nRanges;i++) {
    if (number[i]) 
      printf("%d between %g and %g",number[i],range[i-1],range[i]);
    if (numberExact[i]) {
      if (number[i])
        printf(", ");
      printf("%d exactly at %g",numberExact[i],range[i]);
    }
    if (number[i]+numberExact[i])
      printf("\n");
  }
  delete [] number;
  delete [] numberExact;
}
void sortOnOther(int * column,
		  const CoinBigIndex * rowStart,
		  int * order,
		  int * other,
		  int nRow,
		  int nInRow,
		  int where)
{
  if (nRow<2||where>=nInRow)
    return;
  // do initial sort
  int kRow;
  int iRow;
  for ( kRow=0;kRow<nRow;kRow++) {
    iRow = order[kRow];
    other[kRow]=column[rowStart[iRow]+where];
  }
  CoinSort_2(other,other+nRow,order);
  int first=0;
  iRow=order[0];
  int firstC=column[rowStart[iRow]+where];
  kRow=1;
  while (kRow<nRow) {
    int lastC=9999999;;
    for (;kRow<nRow+1;kRow++) {
      if (kRow<nRow) {
	iRow=order[kRow];
	lastC=column[rowStart[iRow]+where];
      } else {
	lastC=9999999;
      }
      if (lastC>firstC) 
	break;
    }
    // sort
    sortOnOther(column,rowStart,order+first,other,kRow-first,
		 nInRow,where+1);
    firstC=lastC;
    first=kRow;
  }
}
static void statistics(ClpSimplex * originalModel, ClpSimplex * model)
{
  int numberColumns = originalModel->numberColumns();
  const char * integerInformation  = originalModel->integerInformation(); 
  const double * columnLower = originalModel->columnLower();
  const double * columnUpper = originalModel->columnUpper();
  int numberIntegers=0;
  int numberBinary=0;
  int iRow,iColumn;
  if (integerInformation) {
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (integerInformation[iColumn]) {
        if (columnUpper[iColumn]>columnLower[iColumn]) {
          numberIntegers++;
          if (columnUpper[iColumn]==0.0&&columnLower[iColumn]==1) 
            numberBinary++;
        }
      }
    }
  }
  numberColumns = model->numberColumns();
  int numberRows = model->numberRows();
  columnLower = model->columnLower();
  columnUpper = model->columnUpper();
  const double * rowLower = model->rowLower();
  const double * rowUpper = model->rowUpper();
  const double * objective = model->objective();
  CoinPackedMatrix * matrix = model->matrix();
  CoinBigIndex numberElements = matrix->getNumElements();
  const int * columnLength = matrix->getVectorLengths();
  //const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const double * elementByColumn = matrix->getElements();
  int * number = new int[numberRows+1];
  memset(number,0,(numberRows+1)*sizeof(int));
  int numberObjSingletons=0;
  /* cType
     0 0/inf, 1 0/up, 2 lo/inf, 3 lo/up, 4 free, 5 fix, 6 -inf/0, 7 -inf/up,
     8 0/1
  */ 
  int cType[9];
  std::string cName[]={"0.0->inf,","0.0->up,","lo->inf,","lo->up,","free,","fixed,","-inf->0.0,",
                       "-inf->up,","0.0->1.0"};
  int nObjective=0;
  memset(cType,0,sizeof(cType));
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int length=columnLength[iColumn];
    if (length==1&&objective[iColumn])
      numberObjSingletons++;
    number[length]++;
    if (objective[iColumn])
      nObjective++;
    if (columnLower[iColumn]>-1.0e20) {
      if (columnLower[iColumn]==0.0) {
        if (columnUpper[iColumn]>1.0e20)
          cType[0]++;
        else if (columnUpper[iColumn]==1.0)
          cType[8]++;
        else if (columnUpper[iColumn]==0.0)
          cType[5]++;
        else
          cType[1]++;
      } else {
        if (columnUpper[iColumn]>1.0e20) 
          cType[2]++;
        else if (columnUpper[iColumn]==columnLower[iColumn])
          cType[5]++;
        else
          cType[3]++;
      }
    } else {
      if (columnUpper[iColumn]>1.0e20) 
        cType[4]++;
      else if (columnUpper[iColumn]==0.0) 
        cType[6]++;
      else
        cType[7]++;
    }
  }
  /* rType
     0 E 0, 1 E 1, 2 E -1, 3 E other, 4 G 0, 5 G 1, 6 G other, 
     7 L 0,  8 L 1, 9 L other, 10 Range 0/1, 11 Range other, 12 free 
  */ 
  int rType[13];
  std::string rName[]={"E 0.0,","E 1.0,","E -1.0,","E other,","G 0.0,","G 1.0,","G other,",
                       "L 0.0,","L 1.0,","L other,","Range 0.0->1.0,","Range other,","Free"};
  memset(rType,0,sizeof(rType));
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowLower[iRow]>-1.0e20) {
      if (rowLower[iRow]==0.0) {
        if (rowUpper[iRow]>1.0e20)
          rType[4]++;
        else if (rowUpper[iRow]==1.0)
          rType[10]++;
        else if (rowUpper[iRow]==0.0)
          rType[0]++;
        else
          rType[11]++;
      } else if (rowLower[iRow]==1.0) {
        if (rowUpper[iRow]>1.0e20) 
          rType[5]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[1]++;
        else
          rType[11]++;
      } else if (rowLower[iRow]==-1.0) {
        if (rowUpper[iRow]>1.0e20) 
          rType[6]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[2]++;
        else
          rType[11]++;
      } else {
        if (rowUpper[iRow]>1.0e20) 
          rType[6]++;
        else if (rowUpper[iRow]==rowLower[iRow])
          rType[3]++;
        else
          rType[11]++;
      }
    } else {
      if (rowUpper[iRow]>1.0e20) 
        rType[12]++;
      else if (rowUpper[iRow]==0.0) 
        rType[7]++;
      else if (rowUpper[iRow]==1.0) 
        rType[8]++;
      else
        rType[9]++;
    }
  }
  // Basic statistics
  printf("\n\nProblem has %d rows, %d columns (%d with objective) and %d elements\n",
         numberRows,numberColumns,nObjective,numberElements);
  if (number[0]+number[1]) {
    printf("There are ");
    if (numberObjSingletons)
      printf("%d singletons with objective ",numberObjSingletons);
    int numberNoObj = number[1]-numberObjSingletons;
    if (numberNoObj)
      printf("%d singletons with no objective ",numberNoObj);
    if (number[0])
      printf("** %d columns have no entries",number[0]);
    printf("\n");
  }
  printf("Column breakdown:\n");
  int k;
  for (k=0;k<static_cast<int> (sizeof(cType)/sizeof(int));k++) {
    printf("%d of type %s ",cType[k],cName[k].c_str());
    if (((k+1)%3)==0)
      printf("\n");
  }
  if ((k%3)!=0)
    printf("\n");
  printf("Row breakdown:\n");
  for (k=0;k<static_cast<int> (sizeof(rType)/sizeof(int));k++) {
    printf("%d of type %s ",rType[k],rName[k].c_str());
    if (((k+1)%3)==0)
      printf("\n");
  }
  if ((k%3)!=0)
    printf("\n");
#define SYM
#ifndef SYM
  if (model->logLevel()<2)
    return ;
#endif
  int kMax = model->logLevel()>3 ? 1000000 : 10;
  k=0;
  for (iRow=1;iRow<=numberRows;iRow++) {
    if (number[iRow]) {
      k++;
      printf("%d columns have %d entries\n",number[iRow],iRow);
      if (k==kMax)
        break;
    }
  }
  if (k<numberRows) {
    int kk=k;
    k=0;
    for (iRow=numberRows;iRow>=1;iRow--) {
      if (number[iRow]) {
        k++;
        if (k==kMax)
          break;
      }
    }
    if (k>kk) {
      printf("\n    .........\n\n");
      iRow=k;
      k=0;
      for (;iRow<numberRows;iRow++) {
        if (number[iRow]) {
          k++;
          printf("%d columns have %d entries\n",number[iRow],iRow);
          if (k==kMax)
            break;
        }
      }
    }
  }
  delete [] number;
  printf("\n\n");
  if (model->logLevel()==63
#ifdef SYM
      ||true
#endif
      ) {
    // get column copy
    CoinPackedMatrix columnCopy = *matrix;
    const int * columnLength = columnCopy.getVectorLengths();
    number = new int[numberRows+1];
    memset(number,0,(numberRows+1)*sizeof(int));
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int length=columnLength[iColumn];
      number[length]++;
    }
    k=0;
    for (iRow=1;iRow<=numberRows;iRow++) {
      if (number[iRow]) {
	k++;
      }
    }
    int * row = columnCopy.getMutableIndices();
    const CoinBigIndex * columnStart = columnCopy.getVectorStarts();
    double * element = columnCopy.getMutableElements();
    int * order = new int[numberColumns];
    int * other = new int[numberColumns];
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int length=columnLength[iColumn];
      order[iColumn]=iColumn;
      other[iColumn]=length;
      CoinBigIndex start = columnStart[iColumn];
      CoinSort_2(row+start,row+start+length,element+start);
    }
    CoinSort_2(other,other+numberColumns,order);
    int jColumn=number[0]+number[1];
    for (iRow=2;iRow<=numberRows;iRow++) {
      if (number[iRow]) {
	printf("XX %d columns have %d entries\n",number[iRow],iRow);
	int kColumn = jColumn+number[iRow];
	sortOnOther(row,columnStart,
		     order+jColumn,other,number[iRow],iRow,0);
	// Now print etc
	if (iRow<500000) {
	  for (int lColumn =jColumn;lColumn<kColumn;lColumn++) {
	    iColumn = order[lColumn];
	    CoinBigIndex start = columnStart[iColumn];
	    if (model->logLevel()==63) {
	      printf("column %d %g <= ",iColumn,columnLower[iColumn]);
	      for (CoinBigIndex i=start;i<start+iRow;i++)
		printf("( %d, %g) ",row[i],element[i]);
	      printf("<= %g\n",columnUpper[iColumn]);
	    }
	  }
	}
	jColumn =kColumn;
      }
    }
    delete [] order;
    delete [] other;
    delete [] number;
  }
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int * rowLength = rowCopy.getVectorLengths();
  number = new int[numberColumns+1];
  memset(number,0,(numberColumns+1)*sizeof(int));
  for (iRow=0;iRow<numberRows;iRow++) {
    int length=rowLength[iRow];
    number[length]++;
  }
  if (number[0])
    printf("** %d rows have no entries\n",number[0]);
  k=0;
  for (iColumn=1;iColumn<=numberColumns;iColumn++) {
    if (number[iColumn]) {
      k++;
      printf("%d rows have %d entries\n",number[iColumn],iColumn);
      if (k==kMax)
        break;
    }
  }
  if (k<numberColumns) {
    int kk=k;
    k=0;
    for (iColumn=numberColumns;iColumn>=1;iColumn--) {
      if (number[iColumn]) {
        k++;
        if (k==kMax)
          break;
      }
    }
    if (k>kk) {
      printf("\n    .........\n\n");
      iColumn=k;
      k=0;
      for (;iColumn<numberColumns;iColumn++) {
        if (number[iColumn]) {
          k++;
          printf("%d rows have %d entries\n",number[iColumn],iColumn);
          if (k==kMax)
            break;
        }
      }
    }
  }
  if (model->logLevel()==63
#ifdef SYM
      ||true
#endif
      ) {
    int * column = rowCopy.getMutableIndices();
    const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    double * element = rowCopy.getMutableElements();
    int * order = new int[numberRows];
    int * other = new int[numberRows];
    for (iRow=0;iRow<numberRows;iRow++) {
      int length=rowLength[iRow];
      order[iRow]=iRow;
      other[iRow]=length;
      CoinBigIndex start = rowStart[iRow];
      CoinSort_2(column+start,column+start+length,element+start);
    }
    CoinSort_2(other,other+numberRows,order);
    int jRow=number[0]+number[1];
    double * weight = new double[numberRows];
    double * randomColumn = new double[numberColumns+1];
    double * randomRow = new double [numberRows+1];
    int * sortRow = new int [numberRows];
    int * possibleRow = new int [numberRows];
    int * backRow = new int [numberRows];
    int * stackRow = new int [numberRows];
    int * sortColumn = new int [numberColumns];
    int * possibleColumn = new int [numberColumns];
    int * backColumn = new int [numberColumns];
    int * backColumn2 = new int [numberColumns];
    int * mapRow = new int [numberRows];
    int * mapColumn = new int [numberColumns];
    int * stackColumn = new int [numberColumns];
    double randomLower = CoinDrand48();
    double randomUpper = CoinDrand48();
    double randomInteger = CoinDrand48();
    int * startAdd = new int[numberRows+1];
    int * columnAdd = new int [2*numberElements];
    double * elementAdd = new double[2*numberElements];
    int nAddRows=0;
    startAdd[0]=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      randomColumn[iColumn] = CoinDrand48();
      backColumn2[iColumn]=-1;
    }
    for (iColumn=2;iColumn<=numberColumns;iColumn++) {
      if (number[iColumn]) {
	printf("XX %d rows have %d entries\n",number[iColumn],iColumn);
	int kRow = jRow+number[iColumn];
	sortOnOther(column,rowStart,
		     order+jRow,other,number[iColumn],iColumn,0);
	// Now print etc
	if (iColumn<500000) {
	  int nLook=0;
	  for (int lRow =jRow;lRow<kRow;lRow++) {
	    iRow = order[lRow];
	    CoinBigIndex start = rowStart[iRow];
	    if (model->logLevel()==63) {
	      printf("row %d %g <= ",iRow,rowLower[iRow]);
	      for (CoinBigIndex i=start;i<start+iColumn;i++) 
		printf("( %d, %g) ",column[i],element[i]);
	      printf("<= %g\n",rowUpper[iRow]);
	    }
	    int first = column[start];
	    double sum=0.0;
	    for (CoinBigIndex i=start;i<start+iColumn;i++) {
	      int jColumn = column[i];
	      double value = element[i];
	      jColumn -= first;
	      assert (jColumn>=0);
	      sum += value*randomColumn[jColumn];
	    }
	    if (rowLower[iRow]>-1.0e30&&rowLower[iRow])
	      sum += rowLower[iRow]*randomLower;
	    else if (!rowLower[iRow])
	      sum += 1.234567e-7*randomLower;
	    if (rowUpper[iRow]<1.0e30&&rowUpper[iRow])
	      sum += rowUpper[iRow]*randomUpper;
	    else if (!rowUpper[iRow])
	      sum += 1.234567e-7*randomUpper;
	    sortRow[nLook]=iRow;
	    randomRow[nLook++]=sum;
	    // best way is to number unique elements and bounds and use
	    if (fabs(sum)>1.0e4)
	      sum *= 1.0e-6;
	    weight[iRow]=sum;
	  }
	  assert (nLook<=numberRows);
	  CoinSort_2(randomRow,randomRow+nLook,sortRow);
	  randomRow[nLook]=COIN_DBL_MAX;
	  double last=-COIN_DBL_MAX;
	  int iLast=-1;
	  for (int iLook=0;iLook<nLook+1;iLook++) {
	    if (randomRow[iLook]>last) {
	      if (iLast>=0) {
		int n=iLook-iLast;
		if (n>1) {
		  //printf("%d rows possible?\n",n);
		}
	      }
	      iLast=iLook;
	      last=randomRow[iLook];
	    }
	  }
	}
	jRow =kRow;
      }
    }
    CoinPackedMatrix columnCopy = *matrix;
    const int * columnLength = columnCopy.getVectorLengths();
    const int * row = columnCopy.getIndices();
    const CoinBigIndex * columnStart = columnCopy.getVectorStarts();
    const double * elementByColumn = columnCopy.getElements();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int length=columnLength[iColumn];
      CoinBigIndex start = columnStart[iColumn];
      double sum = objective[iColumn];
      if (columnLower[iColumn]>-1.0e30&&columnLower[iColumn])
	sum += columnLower[iColumn]*randomLower;
      else if (!columnLower[iColumn])
	sum += 1.234567e-7*randomLower;
      if (columnUpper[iColumn]<1.0e30&&columnUpper[iColumn])
	sum += columnUpper[iColumn]*randomUpper;
      else if (!columnUpper[iColumn])
	sum += 1.234567e-7*randomUpper;
      if (model->isInteger(iColumn))
	sum += 9.87654321e-6*randomInteger;
      for (CoinBigIndex i=start;i<start+length;i++) {
	int iRow = row[i];
	sum += elementByColumn[i]*weight[iRow];
      }
      sortColumn[iColumn]=iColumn;
      randomColumn[iColumn]=sum;
    }
    {
      CoinSort_2(randomColumn,randomColumn+numberColumns,sortColumn);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	int i=sortColumn[iColumn];
	backColumn[i]=iColumn;
      }
      randomColumn[numberColumns]=COIN_DBL_MAX;
      double last=-COIN_DBL_MAX;
      int iLast=-1;
      for (int iLook=0;iLook<numberColumns+1;iLook++) {
	if (randomColumn[iLook]>last) {
	  if (iLast>=0) {
	    int n=iLook-iLast;
	    if (n>1) {
	      //printf("%d columns possible?\n",n);
	    }
	    for (int i=iLast;i<iLook;i++) {
	      possibleColumn[sortColumn[i]]=n;
	    }
	  }
	  iLast=iLook;
	  last=randomColumn[iLook];
	}
      }
      for (iRow =0;iRow<numberRows;iRow++) {
	CoinBigIndex start = rowStart[iRow];
	double sum=0.0;
	int length=rowLength[iRow];
	for (CoinBigIndex i=start;i<start+length;i++) {
	  int jColumn = column[i];
	  double value = element[i];
	  jColumn=backColumn[jColumn];
	  sum += value*randomColumn[jColumn];
	  //if (iColumn==23089||iRow==23729)
	  //printf("row %d cola %d colb %d value %g rand %g sum %g\n",
	  //   iRow,jColumn,column[i],value,randomColumn[jColumn],sum);
	}
	sortRow[iRow]=iRow;
	randomRow[iRow]=weight[iRow];
	randomRow[iRow]=sum;
      }
      CoinSort_2(randomRow,randomRow+numberRows,sortRow);
      for (iRow=0;iRow<numberRows;iRow++) {
	int i=sortRow[iRow];
	backRow[i]=iRow;
      }
      randomRow[numberRows]=COIN_DBL_MAX;
      last=-COIN_DBL_MAX;
      iLast=-1;
      // Do backward indices from order
      for (iRow=0;iRow<numberRows;iRow++) {
	other[order[iRow]]=iRow;
      }
      for (int iLook=0;iLook<numberRows+1;iLook++) {
	if (randomRow[iLook]>last) {
	  if (iLast>=0) {
	    int n=iLook-iLast;
	    if (n>1) {
	      //printf("%d rows possible?\n",n);
	      // Within group sort as for original "order"
	      for (int i=iLast;i<iLook;i++) {
		int jRow=sortRow[i];
		order[i]=other[jRow];
	      }
	      CoinSort_2(order+iLast,order+iLook,sortRow+iLast);
	    }
	    for (int i=iLast;i<iLook;i++) {
	      possibleRow[sortRow[i]]=n;
	    }
	  }
	  iLast=iLook;
	  last=randomRow[iLook];
	}
      }
      // Temp out
      for (int iLook=0;iLook<numberRows-1000000;iLook++) {
	iRow=sortRow[iLook];
	CoinBigIndex start = rowStart[iRow];
	int length=rowLength[iRow];
	int numberPossible = possibleRow[iRow];
	for (CoinBigIndex i=start;i<start+length;i++) {
	  int jColumn = column[i];
	  if (possibleColumn[jColumn]!=numberPossible)
	    numberPossible=-1;
	}
	int n=numberPossible;
	if (numberPossible>1) {
	  //printf("pppppossible %d\n",numberPossible);
	  for (int jLook=iLook+1;jLook<iLook+numberPossible;jLook++) {
	    int jRow=sortRow[jLook];
	    CoinBigIndex start2 = rowStart[jRow];
	    assert (numberPossible==possibleRow[jRow]);
	    assert(length==rowLength[jRow]);
	    for (CoinBigIndex i=start2;i<start2+length;i++) {
	      int jColumn = column[i];
	      if (possibleColumn[jColumn]!=numberPossible)
		numberPossible=-1;
	    }
	  }
	  if (numberPossible<2) {
	    // switch off
	    for (int jLook=iLook;jLook<iLook+n;jLook++) 
	      possibleRow[sortRow[jLook]]=-1;
	  }
	  // skip rest
	  iLook += n-1;
	} else {
	  possibleRow[iRow]=-1;
	}
      }
      for (int iLook=0;iLook<numberRows;iLook++) {
	iRow=sortRow[iLook];
	int numberPossible = possibleRow[iRow];
	// Only if any integers
	int numberIntegers=0;
	CoinBigIndex start = rowStart[iRow];
	int length=rowLength[iRow];
	for (CoinBigIndex i=start;i<start+length;i++) {
	  int jColumn = column[i];
	  if (model->isInteger(jColumn))
	    numberIntegers++;
	}
	if (numberPossible>1&&!numberIntegers) {
	  //printf("possible %d - but no integers\n",numberPossible);
	}
	if (numberPossible>1&&(numberIntegers||true)) {
	  // 
	  printf("possible %d - %d integers\n",numberPossible,numberIntegers);
	  int lastLook=iLook;
	  int nMapRow=-1;
	  for (int jLook=iLook+1;jLook<iLook+numberPossible;jLook++) {
	    // stop if too many failures
	    if (jLook>iLook+10&&nMapRow<0)
	      break;
	    // Create identity mapping
	    int i;
	    for (i=0;i<numberRows;i++)
	      mapRow[i]=i;
	    for (i=0;i<numberColumns;i++)
	      mapColumn[i]=i;
	    int offset=jLook-iLook;
	    int nStackC=0;
	    // build up row and column mapping
	    int nStackR=1;
	    stackRow[0]=iLook;
	    bool good=true;
	    while (nStackR) {
	      nStackR--;
	      int look1 = stackRow[nStackR];
	      int look2 = look1 + offset;
	      assert (randomRow[look1]==randomRow[look2]);
	      int row1=sortRow[look1];
	      int row2=sortRow[look2];
	      assert (mapRow[row1]==row1);
	      assert (mapRow[row2]==row2);
	      mapRow[row1]=row2;
	      mapRow[row2]=row1;
	      CoinBigIndex start1 = rowStart[row1];
	      CoinBigIndex offset2 = rowStart[row2]-start1;
	      int length=rowLength[row1];
	      assert( length==rowLength[row2]);
	      for (CoinBigIndex i=start1;i<start1+length;i++) {
		int jColumn1 = column[i];
		int jColumn2 = column[i+offset2];
		if (randomColumn[backColumn[jColumn1]]!=
		    randomColumn[backColumn[jColumn2]]) {
		  good=false;
		  break;
		}
		if (mapColumn[jColumn1]==jColumn1) {
		  // not touched
		  assert (mapColumn[jColumn2]==jColumn2);
		  if (jColumn1!=jColumn2) {
		    // Put on stack
		    mapColumn[jColumn1]=jColumn2;
		    mapColumn[jColumn2]=jColumn1;
		    stackColumn[nStackC++]=jColumn1;
		  }
		} else {
		  if(mapColumn[jColumn1]!=jColumn2||
		     mapColumn[jColumn2]!=jColumn1) {
		    // bad
		    good=false;
		    printf("bad col\n");
		    break;
		  }
		}
	      }
	      if (!good)
		break;
	      while (nStackC) {
		nStackC--;
		int iColumn = stackColumn[nStackC];
		int iColumn2 = mapColumn[iColumn];
		assert (iColumn!=iColumn2);
		int length=columnLength[iColumn];
		assert (length==columnLength[iColumn2]);
		CoinBigIndex start = columnStart[iColumn];
		CoinBigIndex offset2 = columnStart[iColumn2]-start;
		for (CoinBigIndex i=start;i<start+length;i++) {
		  int iRow = row[i];
		  int iRow2 = row[i+offset2];
		  if (mapRow[iRow]==iRow) {
		    // First (but be careful)
		    if (iRow!=iRow2) {
		      //mapRow[iRow]=iRow2;
		      //mapRow[iRow2]=iRow;
		      int iBack=backRow[iRow];
		      int iBack2=backRow[iRow2];
		      if (randomRow[iBack]==randomRow[iBack2]&&
			  iBack2-iBack==offset) {
			stackRow[nStackR++]=iBack;
		      } else {
			//printf("randomRow diff - weights %g %g\n",
			//     weight[iRow],weight[iRow2]);
			// bad
			good=false;
			break;
		      }
		    }
		  } else {
		    if(mapRow[iRow]!=iRow2||
		       mapRow[iRow2]!=iRow) {
		      // bad
		      good=false;
		      printf("bad row\n");
		      break;
		    }
		  }
		}
		if (!good)
		  break;
	      }
	    }
	    // then check OK
	    if (good) {
	      for (iRow =0;iRow<numberRows;iRow++) {
		CoinBigIndex start = rowStart[iRow];
		int length=rowLength[iRow];
		if (mapRow[iRow]==iRow) {
		  for (CoinBigIndex i=start;i<start+length;i++) {
		    int jColumn = column[i];
		    backColumn2[jColumn]=i-start;
		  }
		  for (CoinBigIndex i=start;i<start+length;i++) {
		    int jColumn = column[i];
		    if (mapColumn[jColumn]!=jColumn) {
		      int jColumn2 =mapColumn[jColumn];
		      CoinBigIndex i2 = backColumn2[jColumn2];
		      if (i2<0) {
			good=false;
		      } else if (element[i]!=element[i2+start]) {
			good=false;
		      }
		    }
		  }
		  for (CoinBigIndex i=start;i<start+length;i++) {
		    int jColumn = column[i];
		    backColumn2[jColumn]=-1;
		  }
		} else {
		  int row2 = mapRow[iRow];
		  assert (iRow = mapRow[row2]);
		  if (rowLower[iRow]!=rowLower[row2]||
		      rowLower[row2]!=rowLower[iRow])
		    good=false;
		  CoinBigIndex offset2 = rowStart[row2]-start;
		  for (CoinBigIndex i=start;i<start+length;i++) {
		    int jColumn = column[i];
		    double value = element[i];
		    int jColumn2 = column[i+offset2];
		    double value2 = element[i+offset2];
		    if (value!=value2||mapColumn[jColumn]!=jColumn2||
			mapColumn[jColumn2]!=jColumn)
		      good=false;
		  }
		}
	      }
	      if (good) {
		// check rim
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  if (mapColumn[iColumn]!=iColumn) {
		    int iColumn2 = mapColumn[iColumn];
		    if(objective[iColumn]!=objective[iColumn2])
		      good=false;
		    if (columnLower[iColumn]!=columnLower[iColumn2])
		      good=false;
		    if (columnUpper[iColumn]!=columnUpper[iColumn2])
		      good=false;
		    if (model->isInteger(iColumn)!=model->isInteger(iColumn2))
		      good=false;
		  }
		}
	      }
	      if (good) {
		// temp
		if (nMapRow<0) {
		  //const double * solution = model->primalColumnSolution();
		  // find mapped
		  int nMapColumn=0;
		  for (int i=0;i<numberColumns;i++) {
		    if (mapColumn[i]>i) 
		      nMapColumn++;
		  }
		  nMapRow=0;
		  int kRow=-1;
		  for (int i=0;i<numberRows;i++) {
		    if (mapRow[i]>i) { 
		      nMapRow++;
		      kRow=i;
		    }
		  }
		  printf("%d columns, %d rows\n",nMapColumn,nMapRow);
		  if (nMapRow==1) {
		    CoinBigIndex start = rowStart[kRow];
		    int length=rowLength[kRow];
		    printf("%g <= ",rowLower[kRow]);
		    for (CoinBigIndex i=start;i<start+length;i++) {
		      int jColumn = column[i];
		      if (mapColumn[jColumn]!=jColumn) 
			printf("* ");
		      printf("%d,%g ",jColumn,element[i]);
		    }
		    printf("<= %g\n",rowUpper[kRow]);
		  }
		}
		// temp
		int row1=sortRow[lastLook];
		int row2=sortRow[jLook];
		lastLook=jLook;
		CoinBigIndex start1 = rowStart[row1];
		CoinBigIndex offset2 = rowStart[row2]-start1;
		int length=rowLength[row1];
		assert( length==rowLength[row2]);
		CoinBigIndex put=startAdd[nAddRows];
		double multiplier=length<11 ? 2.0 : 1.125;
		double value=1.0;
		for (CoinBigIndex i=start1;i<start1+length;i++) {
		  int jColumn1 = column[i];
		  int jColumn2 = column[i+offset2];
		  columnAdd[put]=jColumn1;
		  elementAdd[put++]=value;
		  columnAdd[put]=jColumn2;
		  elementAdd[put++]=-value;
		  value *= multiplier;
		}
		nAddRows++;
		startAdd[nAddRows]=put;
	      } else {
		printf("ouch - did not check out as good\n");
	      }
	    }
	  }
	  // skip rest
	  iLook += numberPossible-1;
	}
      }
    }
    if (nAddRows) {
      double * lower = new double [nAddRows];
      double * upper = new double[nAddRows];
      int i;
      //const double * solution = model->primalColumnSolution();
      for (i=0;i<nAddRows;i++) {
	lower[i]=0.0;
	upper[i]=COIN_DBL_MAX;
      }
      printf("Adding %d rows with %d elements\n",nAddRows,
	     startAdd[nAddRows]);
      //ClpSimplex newModel(*model);
      //newModel.addRows(nAddRows,lower,upper,startAdd,columnAdd,elementAdd);
      //newModel.writeMps("modified.mps");
      delete [] lower;
      delete [] upper;
    }
    delete [] startAdd;
    delete [] columnAdd;
    delete [] elementAdd;
    delete [] order;
    delete [] other;
    delete [] randomColumn;
    delete [] weight;
    delete [] randomRow;
    delete [] sortRow;
    delete [] backRow;
    delete [] possibleRow;
    delete [] sortColumn;
    delete [] backColumn;
    delete [] backColumn2;
    delete [] possibleColumn;
    delete [] mapRow;
    delete [] mapColumn;
    delete [] stackRow;
    delete [] stackColumn;
  }
  delete [] number;
  // Now do breakdown of ranges
  breakdown("Elements",numberElements,elementByColumn);
  breakdown("RowLower",numberRows,rowLower);
  breakdown("RowUpper",numberRows,rowUpper);
  breakdown("ColumnLower",numberColumns,columnLower);
  breakdown("ColumnUpper",numberColumns,columnUpper);
  breakdown("Objective",numberColumns,objective);
}
static bool maskMatches(const int * starts, char ** masks,
			std::string & check)
{
  // back to char as I am old fashioned
  const char * checkC = check.c_str();
  int length = strlen(checkC);
  while (checkC[length-1]==' ')
    length--;
  for (int i=starts[length];i<starts[length+1];i++) {
    char * thisMask = masks[i];
    int k;
    for ( k=0;k<length;k++) {
      if (thisMask[k]!='?'&&thisMask[k]!=checkC[k]) 
	break;
    }
    if (k==length)
      return true;
  }
  return false;
}
static void clean(char * temp)
{
  char * put = temp;
  while (*put>=' ')
    put++;
  *put='\0';
}
static void generateCode(const char * fileName,int type)
{
  FILE * fp = fopen(fileName,"r");
  assert (fp);
  int numberLines=0;
#define MAXLINES 500
#define MAXONELINE 200
  char line[MAXLINES][MAXONELINE];
  while (fgets(line[numberLines],MAXONELINE,fp)) {
    assert (numberLines<MAXLINES);
    clean(line[numberLines]);
    numberLines++;
  }
  fclose(fp);
  // add in actual solve
  strcpy(line[numberLines],"5  clpModel->initialSolve(clpSolve);");
  numberLines++;
  fp = fopen(fileName,"w");
  assert (fp);
  char apo='"';	  
  char backslash = '\\';

  fprintf(fp,"#include %cClpSimplex.hpp%c\n",apo,apo);
  fprintf(fp,"#include %cClpSolve.hpp%c\n",apo,apo);

  fprintf(fp,"\nint main (int argc, const char *argv[])\n{\n");
  fprintf(fp,"  ClpSimplex  model;\n");
  fprintf(fp,"  int status=1;\n");
  fprintf(fp,"  if (argc<2)\n");
  fprintf(fp,"    fprintf(stderr,%cPlease give file name%cn%c);\n",
          apo,backslash,apo);
  fprintf(fp,"  else\n");
  fprintf(fp,"    status=model.readMps(argv[1],true);\n");
  fprintf(fp,"  if (status) {\n");
  fprintf(fp,"    fprintf(stderr,%cBad readMps %%s%cn%c,argv[1]);\n",
                apo,backslash,apo);
  fprintf(fp,"    exit(1);\n");
  fprintf(fp,"  }\n\n");
  fprintf(fp,"  // Now do requested saves and modifications\n");
  fprintf(fp,"  ClpSimplex * clpModel = & model;\n");
  int wanted[9];
  memset(wanted,0,sizeof(wanted));
  wanted[0]=wanted[3]=wanted[5]=wanted[8]=1;
  if (type>0) 
    wanted[1]=wanted[6]=1;
  if (type>1) 
    wanted[2]=wanted[4]=wanted[7]=1;
  std::string header[9]=
  { "","Save values","Redundant save of default values","Set changed values",
    "Redundant set default values","Solve","Restore values","Redundant restore values","Add to model"};
  for (int iType=0;iType<9;iType++) {
    if (!wanted[iType])
      continue;
    int n=0;
    int iLine;
    for (iLine=0;iLine<numberLines;iLine++) {
      if (line[iLine][0]=='0'+iType) {
        if (!n)
          fprintf(fp,"\n  // %s\n\n",header[iType].c_str());
        n++;
        fprintf(fp,"%s\n",line[iLine]+1);
      }
    }
  }
  fprintf(fp,"\n  // Now you would use solution etc etc\n\n");
  fprintf(fp,"  return 0;\n}\n");
  fclose(fp);
  printf("C++ file written to %s\n",fileName);
}
/*
  Version 1.00.00 October 13 2004.
  1.00.01 October 18.  Added basis handling helped/prodded by Thorsten Koch.
  Also modifications to make faster with sbb (I hope I haven't broken anything).
  1.00.02 March 21 2005.  Redid ClpNonLinearCost to save memory also redid
  createRim to try and improve cache characteristics.
  1.00.03 April 8 2005.  Added Volume algorithm as crash and made code more
  robust on testing.  Also added "either" and "tune" algorithm.
  1.01.01 April 12 2005.  Decided to go to different numbering.  Backups will
  be last 2 digits while middle 2 are for improvements.  Still take a long 
  time to get to 2.00.01
  1.01.02 May 4 2005.  Will be putting in many changes - so saving stable version
  1.02.01 May 6 2005.  Lots of changes to try and make faster and more stable in
  branch and cut.
  1.02.02 May 19 2005.  Stuff for strong branching and some improvements to simplex
  1.03.01 May 24 2006.  Lots done but I can't remember what!
  1.03.03 June 13 2006.  For clean up after dual perturbation
  1.04.01 June 26 2007.  Lots of changes but I got lazy
  1.05.00 June 27 2007.  This is trunk so when gets to stable will be 1.5
 */
