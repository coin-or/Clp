// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
   
#include "CoinPragma.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>


#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
// History since 1.0 at end
#define CLPVERSION "1.00.02"

#include "CoinMpsIO.hpp"

#include "ClpFactorization.hpp"
#include "CoinTime.hpp"
#include "ClpSimplex.hpp"
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

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

int mainTest (int argc, const char *argv[],int algorithm,
	      ClpSimplex empty, bool doPresolve,int doIdiot);
// Returns next valid field
int CbcOrClpRead_mode=1;
FILE * CbcOrClpReadCommand=stdin;
int main (int argc, const char *argv[])
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
    int printOptions=0;
    int presolveOptions=0;
    int doCrash=0;
    int doSprint=-1;
    // set reasonable defaults
    int preSolve=5;
    bool preSolveFile=false;
    models->setPerturbation(50);
    models->messageHandler()->setPrefix(false);
    std::string directory ="./";
    std::string importFile ="";
    std::string exportFile ="default.mps";
    std::string importBasisFile ="";
    int basisHasValues=0;
    std::string exportBasisFile ="default.bas";
    std::string saveFile ="default.prob";
    std::string restoreFile ="default.prob";
    std::string solutionFile ="stdout";
    CbcOrClpParam parameters[CBCMAXPARAMETERS];
    int numberParameters ;
    establishParams(numberParameters,parameters) ;
    parameters[whichParam(BASISIN,numberParameters,parameters)].setStringValue(importBasisFile);
    parameters[whichParam(BASISOUT,numberParameters,parameters)].setStringValue(exportBasisFile);
    parameters[whichParam(DIRECTORY,numberParameters,parameters)].setStringValue(directory);
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
    parameters[whichParam(SPRINT,numberParameters,parameters)].setIntValue(doSprint);
    
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
	      if (parameters[iParam].displayThis()&&type>=limits[iType]
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
	    parameters[iParam].setDoubleValue(value);
	    parameters[iParam].setDoubleParameter(models+iModel,value);
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].doubleValue()<<std::endl;
	  }
	} else if (type<201) {
	  // get next field as int
	  int value = CoinReadGetIntField(argc,argv,&valid);
	  if (!valid) {
	    parameters[iParam].setIntValue(value);
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
	    else if (parameters[iParam].type()==PRESOLVEOPTIONS)
	      presolveOptions = value;
	    else if (parameters[iParam].type()==PRINTOPTIONS)
	      printOptions = value;
	    else
	      parameters[iParam].setIntParameter(models+iModel,value);
	  } else if (valid==1) {
	    abort();
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
		ClpPrimalColumnSteepest steep(2);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==4) {
		ClpPrimalColumnSteepest steep(1);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==5) {
		ClpPrimalColumnSteepest steep(4);
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
	    if (goodModels[iModel]) {
	      ClpSolve::SolveType method;
	      ClpSolve::PresolveType presolveType;
	      ClpSimplex * model2 = models+iModel;
	      ClpSolve solveOptions;
	      if (preSolve!=5&&preSolve)
		presolveType=ClpSolve::presolveNumber;
	      else if (preSolve)
		presolveType=ClpSolve::presolveOn;
	      else
		presolveType=ClpSolve::presolveOff;
	      solveOptions.setPresolveType(presolveType,preSolve);
	      if (type==DUALSIMPLEX) {
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
		if (gamma)
		  barrierOptions |= 16;
		if (doKKT)
		  barrierOptions |= 32;
		solveOptions.setSpecialOption(4,barrierOptions);
	      }
	      model2->initialSolve(solveOptions);
	      basisHasValues=1;
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
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environ = getenv("HOME");
		  if (environ) {
		    std::string home(environ);
		    field=field.erase(0,1);
		    fileName = home+field;
		  } else {
		    fileName=field;
		  }
		} else {
		  fileName = directory+field;
                  // See if gmpl (model & data)
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
		FILE *fp=fopen(fileName.c_str(),"r");
		if (fp) {
		  // can open - lets go for it
		  fclose(fp);
		  canOpen=true;
                  if (gmpl==2) {
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
                else
                  status= models[iModel].readGMPL(fileName.c_str(),
                                                  (gmpl==2) ? gmplData.c_str() : NULL,
                                                  keepImportNames!=0);
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
		    char * find = strstr(fileName.c_str(),".mps");
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
		char * environ = getenv("HOME");
		if (environ) {
		  std::string home(environ);
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
		if (preSolve) {
		  ClpPresolve pinfo;
		  int presolveOptions2 = presolveOptions&~0x40000000;
		  if ((presolveOptions2&0xffff)!=0)
		    pinfo.setPresolveActions(presolveOptions2);
		  if ((printOptions&1)!=0)
		    pinfo.statistics();
		  model2 = 
		    pinfo.presolvedModel(models[iModel],1.0e-8,
					 true,preSolve);
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
		      strdup(model2->rowName(iRow).c_str());
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
		      strdup(model2->columnName(iColumn).c_str());
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
		  char * environ = getenv("HOME");
		  if (environ) {
		    std::string home(environ);
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
		char * environ = getenv("HOME");
		if (environ) {
		  std::string home(environ);
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
		char * environ = getenv("HOME");
		if (environ) {
		  std::string home(environ);
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
		  model2 = 
		    pinfo.presolvedModel(models[iModel],1.0e-8,
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
		char * environ = getenv("HOME");
		if (environ) {
		  std::string home(environ);
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
		dualColumnSolution[iColumn] = dualColumnSolution[iColumn];
		objective[iColumn] = -objective[iColumn];
	      }
	      int iRow;
	      int numberRows=models[iModel].numberRows();
	      double * dualRowSolution = 
		models[iModel].dualRowSolution();
	      for (iRow=0;iRow<numberRows;iRow++) 
		dualRowSolution[iRow] = dualRowSolution[iRow];
              models[iModel].setObjectiveOffset(models[iModel].objectiveOffset());
	    }
	    break;
	  case DIRECTORY:
	    {
	      std::string name = CoinReadGetString(argc,argv);
	      if (name!="EOL") {
		int length=name.length();
		if (name[length-1]=='/'||name[length-1]=='\\')
		  directory=name;
		else
		  directory = name+"/";
		parameters[iParam].setStringValue(directory);
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
	    {
	      // create fields for unitTest
	      const char * fields[4];
	      int nFields=2;
	      fields[0]="fake main from unitTest";
	      fields[1]="-netlib";
	      if (directory!="./") {
		fields[2]="-netlibDir";
		fields[3]=directory.c_str();
		nFields=4;
	      }
	      int algorithm;
	      if (type==NETLIB_DUAL) {
		std::cerr<<"Doing netlib with dual agorithm"<<std::endl;
		algorithm =0;
	      } else if (type==NETLIB_BARRIER) {
		std::cerr<<"Doing netlib with barrier agorithm"<<std::endl;
		algorithm =2;
	      } else if (type==NETLIB_EITHER) {
		std::cerr<<"Doing netlib with dual or primal agorithm"<<std::endl;
		algorithm =3;
	      } else {
		std::cerr<<"Doing netlib with primal agorithm"<<std::endl;
		algorithm=1;
	      }
	      mainTest(nFields,fields,algorithm,models[iModel],
		       (preSolve!=0),doIdiot);
	    }
	    break;
	  case UNITTEST:
	    {
	      // create fields for unitTest
	      const char * fields[3];
	      int nFields=1;
	      fields[0]="fake main from unitTest";
	      if (directory!="./") {
		fields[1]="-mpsDir";
		fields[2]=directory.c_str();
		nFields=3;
	      }
	      mainTest(nFields,fields,false,models[iModel],(preSolve!=0),
		       false);
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
	      } else if (field=="stderr") {
		// stderr
		fp=stderr;
	      } else {
		if (field[0]=='/'||field[0]=='\\') {
		  fileName = field;
		} else if (field[0]=='~') {
		  char * environ = getenv("HOME");
		  if (environ) {
		    std::string home(environ);
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
		int printAll = (printOptions>>1)&3;
		char format[6];
		sprintf(format,"%%-%ds",CoinMax(lengthName,8));
		for (iRow=0;iRow<numberRows;iRow++) {
		  int type=printAll;
		  if (primalRowSolution[iRow]>rowUpper[iRow]+primalTolerance||
		      primalRowSolution[iRow]<rowLower[iRow]-primalTolerance) {
		    fprintf(fp,"** ");
		    type=2;
		  } else if (fabs(primalRowSolution[iRow])>1.0e-8) {
		    type=1;
		  } else if (numberRows<50) {
		    type=3;
		  } 
		  if (type) {
		    fprintf(fp,"%7d ",iRow);
		    if (lengthName)
		      fprintf(fp,format,rowNames[iRow].c_str());
		    fprintf(fp,"%15.8g        %15.8g\n",primalRowSolution[iRow],
			    dualRowSolution[iRow]);
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
		  int type=printAll;
		  if (primalColumnSolution[iColumn]>columnUpper[iColumn]+primalTolerance||
		      primalColumnSolution[iColumn]<columnLower[iColumn]-primalTolerance) {
		    fprintf(fp,"** ");
		    type=2;
		  } else if (fabs(primalColumnSolution[iColumn])>1.0e-8) {
		    type=1;
		  } else if (numberColumns<50) {
		    type=3;
		  }
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
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
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
/*
  Version 1.00.00 October 13 2004.
  1.00.01 October 18.  Added basis handline helped/prodded by Thorsten Koch.
  Also modifications to make faster with sbb (I hope I haven't broken anything).
  1.00.02 March 21 2005.  Redid ClpNonLinearCost to save memory also redid
  createRim to try and improve cache characteristics.
 */
