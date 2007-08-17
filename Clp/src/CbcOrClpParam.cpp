// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcOrClpParam.hpp"

#include <string>
#include <iostream>
#include <cassert>

#ifdef COIN_HAS_CBC
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplex.hpp"
#endif
#include "CbcModel.hpp"
#endif
#ifdef COIN_HAS_CLP
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#endif
#ifdef COIN_HAS_READLINE     
#include <readline/readline.h>
#include <readline/history.h>
#endif
#ifdef COIN_HAS_CBC
// from CoinSolve
static char coin_prompt[]="Coin:";
#else
static char coin_prompt[]="Clp:";
#endif
static bool doPrinting=true;
std::string afterEquals="";
void setCbcOrClpPrinting(bool yesNo)
{
  doPrinting=yesNo;
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam () 
  : type_(INVALID),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    lengthName_(0),
    lengthMatch_(0),
    definedKeyWords_(),
    name_(),
    shortHelp_(),
    longHelp_(),
    action_(INVALID),
    currentKeyWord_(-1),
    display_(false),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    whereUsed_(7)
{
}
// Other constructors
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
	   double lower, double upper, CbcOrClpParameterType type,
		    bool display)
  : type_(type),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    whereUsed_(7)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
	   int lower, int upper, CbcOrClpParameterType type,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    whereUsed_(7)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help, 
		    std::string firstValue,
		    CbcOrClpParameterType type,int whereUsed,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(0),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_(""),
    whereUsed_(whereUsed)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
CbcOrClpParam::CbcOrClpParam (std::string name, std::string help,
		    CbcOrClpParameterType type,int whereUsed,
		    bool display)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1),
    display_(display),
    intValue_(-1),
    doubleValue_(-1.0),
    stringValue_("")
{
  whereUsed_=whereUsed;
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcOrClpParam::CbcOrClpParam (const CbcOrClpParam & rhs) 
{  
  type_ = rhs.type_;
  lowerDoubleValue_ = rhs.lowerDoubleValue_;
  upperDoubleValue_ = rhs.upperDoubleValue_;
  lowerIntValue_ = rhs.lowerIntValue_;
  upperIntValue_ = rhs.upperIntValue_;
  lengthName_ = rhs.lengthName_;
  lengthMatch_ = rhs.lengthMatch_;
  definedKeyWords_ = rhs.definedKeyWords_;
  name_ = rhs.name_;
  shortHelp_ = rhs.shortHelp_;
  longHelp_ = rhs.longHelp_;
  action_ = rhs.action_;
  currentKeyWord_ = rhs.currentKeyWord_;
  display_=rhs.display_;
  intValue_=rhs.intValue_;
  doubleValue_=rhs.doubleValue_;
  stringValue_=rhs.stringValue_;
  whereUsed_=rhs.whereUsed_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcOrClpParam::~CbcOrClpParam ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcOrClpParam &
CbcOrClpParam::operator=(const CbcOrClpParam& rhs)
{
  if (this != &rhs) {
    type_ = rhs.type_;
    lowerDoubleValue_ = rhs.lowerDoubleValue_;
    upperDoubleValue_ = rhs.upperDoubleValue_;
    lowerIntValue_ = rhs.lowerIntValue_;
    upperIntValue_ = rhs.upperIntValue_;
    lengthName_ = rhs.lengthName_;
    lengthMatch_ = rhs.lengthMatch_;
    definedKeyWords_ = rhs.definedKeyWords_;
    name_ = rhs.name_;
    shortHelp_ = rhs.shortHelp_;
    longHelp_ = rhs.longHelp_;
    action_ = rhs.action_;
    currentKeyWord_ = rhs.currentKeyWord_;
    display_=rhs.display_;
    intValue_=rhs.intValue_;
    doubleValue_=rhs.doubleValue_;
    stringValue_=rhs.stringValue_;
    whereUsed_=rhs.whereUsed_;
  }
  return *this;
}
void 
CbcOrClpParam::gutsOfConstructor()
{
  std::string::size_type  shriekPos = name_.find('!');
  lengthName_ = name_.length();
  if ( shriekPos==std::string::npos ) {
    //does not contain '!'
    lengthMatch_= lengthName_;
  } else {
    lengthMatch_=shriekPos;
    name_ = name_.substr(0,shriekPos)+name_.substr(shriekPos+1);
    lengthName_--;
  }
}
// Insert string (only valid for keywords)
void 
CbcOrClpParam::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

int 
CbcOrClpParam::matches (std::string input) const
{
  // look up strings to do more elegantly
  if (input.length()>lengthName_) {
    return 0;
  } else {
    unsigned int i;
    for (i=0;i<input.length();i++) {
      if (tolower(name_[i])!=tolower(input[i])) 
	break;
    }
    if (i<input.length()) {
      return 0;
    } else if (i>=lengthMatch_) {
      return 1;
    } else {
      // matched but too short
      return 2;
    }
  }
}
// Returns name which could match
std::string 
CbcOrClpParam::matchName (  ) const
{ 
  if (lengthMatch_==lengthName_) 
    return name_;
  else
    return name_.substr(0,lengthMatch_)+"("+name_.substr(lengthMatch_)+")";
}

// Returns parameter option which matches (-1 if none)
int 
CbcOrClpParam::parameterOption ( std::string check ) const
{
  int numberItems = definedKeyWords_.size();
  if (!numberItems) {
    return -1;
  } else {
    int whichItem=0;
    unsigned int it;
    for (it=0;it<definedKeyWords_.size();it++) {
      std::string thisOne = definedKeyWords_[it];
      std::string::size_type  shriekPos = thisOne.find('!');
      unsigned int length1 = thisOne.length();
      unsigned int length2 = length1;
      if ( shriekPos!=std::string::npos ) {
	//contains '!'
	length2 = shriekPos;
	thisOne = thisOne.substr(0,shriekPos)+
	  thisOne.substr(shriekPos+1);
	length1 = thisOne.length();
      }
      if (check.length()<=length1&&length2<=check.length()) {
	unsigned int i;
	for (i=0;i<check.length();i++) {
	  if (tolower(thisOne[i])!=tolower(check[i])) 
	    break;
	}
	if (i<check.length()) {
	  whichItem++;
	} else if (i>=length2) {
	  break;
	} 
      } else {
	whichItem++;
      }
    }
    if (whichItem<numberItems)
      return whichItem;
    else
      return -1;
  }
}
// Prints parameter options
void 
CbcOrClpParam::printOptions (  ) const
{
  std::cout<<"<Possible options for "<<name_<<" are:";
  unsigned int it;
  for (it=0;it<definedKeyWords_.size();it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type  shriekPos = thisOne.find('!');
    if ( shriekPos!=std::string::npos ) {
      //contains '!'
      thisOne = thisOne.substr(0,shriekPos)+
	"("+thisOne.substr(shriekPos+1)+")";
    }
    std::cout<<" "<<thisOne;
  }
  assert (currentKeyWord_>=0&&currentKeyWord_<(int)definedKeyWords_.size());
  std::string current = definedKeyWords_[currentKeyWord_];
  std::string::size_type  shriekPos = current.find('!');
  if ( shriekPos!=std::string::npos ) {
    //contains '!'
    current = current.substr(0,shriekPos)+
      "("+current.substr(shriekPos+1)+")";
  }
  std::cout<<";\n\tcurrent  "<<current<<">"<<std::endl;
}
// Print action and string
void 
CbcOrClpParam::printString() const
{
  if (name_=="directory")
    std::cout<<"Current working directory is "<<stringValue_<<std::endl;
  else if (name_.substr(0,6)=="printM")
    std::cout<<"Current value of printMask is "<<stringValue_<<std::endl;
  else
    std::cout<<"Current default (if $ as parameter) for "<<name_
	     <<" is "<<stringValue_<<std::endl;
}
void CoinReadPrintit(const char * input)
{
  int length =strlen(input);
  char temp[101];
  int i;
  int n=0;
  for (i=0;i<length;i++) {
    if (input[i]=='\n') {
      temp[n]='\0';
      std::cout<<temp<<std::endl;
      n=0;
    } else if (n>=65&&input[i]==' ') {
      temp[n]='\0';
      std::cout<<temp<<std::endl;
      n=0;
    } else if (n||input[i]!=' ') {
      temp[n++]=input[i];
    }
  }
  if (n) {
    temp[n]='\0';
    std::cout<<temp<<std::endl;
  }
}
// Print Long help
void 
CbcOrClpParam::printLongHelp() const
{
  if (type_>=1&&type_<400) {
    CoinReadPrintit(longHelp_.c_str());
    if (type_<SOLVERLOGLEVEL) {
      printf("<Range of values is %g to %g;\n\tcurrent %g>\n",lowerDoubleValue_,upperDoubleValue_, doubleValue_);
      assert (upperDoubleValue_>lowerDoubleValue_);
    } else if (type_<DIRECTION) {
      printf("<Range of values is %d to %d;\n\tcurrent %d>\n",lowerIntValue_,upperIntValue_,intValue_);
      assert (upperIntValue_>lowerIntValue_);
    } else if (type_<DIRECTORY) {
      printOptions();
    }
  }
}
#ifdef COIN_HAS_CBC
int
CbcOrClpParam::setDoubleParameter (OsiSolverInterface * model,double value) 
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    double oldValue=doubleValue_;
    doubleValue_=value;
    switch(type_) {
    case DUALTOLERANCE:
      model->getDblParam(OsiDualTolerance,oldValue);
      model->setDblParam(OsiDualTolerance,value);
      break;
    case PRIMALTOLERANCE:
      model->getDblParam(OsiPrimalTolerance,oldValue);
      model->setDblParam(OsiPrimalTolerance,value);
      break;
    default:
      break;
    }
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    return 0;
  }
}
#endif
#ifdef COIN_HAS_CLP
int
CbcOrClpParam::setDoubleParameter (ClpSimplex * model,double value) 
{
  double oldValue = doubleValue_;
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    doubleValue_=value;
    switch(type_) {
#ifndef COIN_HAS_CBC
    case DUALTOLERANCE:
      model->setDualTolerance(value);
      break;
    case PRIMALTOLERANCE:
      model->setPrimalTolerance(value);
      break;
#endif
    case DUALBOUND:
      model->setDualBound(value);
      break;
    case PRIMALWEIGHT:
      model->setInfeasibilityCost(value);
      break;
#ifndef COIN_HAS_CBC
    case TIMELIMIT:
      model->setMaximumSeconds(value);
      break;
#endif
    case OBJSCALE:
      model->setObjectiveScale(value);
      break;
    case RHSSCALE:
      model->setRhsScale(value);
      break;
    case PRESOLVETOLERANCE:
      model->setDblParam(ClpPresolveTolerance,value);
      break;
    default:
      break;
    }
    return 0;
  }
}
double 
CbcOrClpParam::doubleParameter (ClpSimplex * model) const
{
  double value;
  switch(type_) {
#ifndef COIN_HAS_CBC
  case DUALTOLERANCE:
    value=model->dualTolerance();
    break;
  case PRIMALTOLERANCE:
    value=model->primalTolerance();
    break;
#endif
  case DUALBOUND:
    value=model->dualBound();
    break;
  case PRIMALWEIGHT:
    value=model->infeasibilityCost();
    break;
#ifndef COIN_HAS_CBC
  case TIMELIMIT:
    value=model->maximumSeconds();
    break;
#endif
  case OBJSCALE:
    value=model->objectiveScale();
    break;
  case RHSSCALE:
    value=model->rhsScale();
    break;
  default:
    value=doubleValue_;
    break;
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (ClpSimplex * model,int value) 
{
  int oldValue = intValue_;
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    intValue_=value;
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    switch(type_) {
    case SOLVERLOGLEVEL:
      model->setLogLevel(value);
      if (value>2)
	model->factorization()->messageLevel(8);
      else
	model->factorization()->messageLevel(0);
      break;
    case MAXFACTOR:
      model->factorization()->maximumPivots(value);
      break;
    case PERTVALUE:
      model->setPerturbation(value);
      break;
    case MAXITERATION:
      model->setMaximumIterations(value);
      break;
    case SPECIALOPTIONS:
      model->setSpecialOptions(value);
#ifdef COIN_HAS_CBC
    case THREADS:
      model->setNumberThreads(value);
      break;
#endif
    default:
      break;
    }
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (ClpSimplex * model) const
{
  int value;
  switch(type_) {
#ifndef COIN_HAS_CBC
  case SOLVERLOGLEVEL:
    value=model->logLevel();
    break;
#endif
  case MAXFACTOR:
    value=model->factorization()->maximumPivots();
    break;
    break;
  case PERTVALUE:
    value=model->perturbation();
    break;
  case MAXITERATION:
    value=model->maximumIterations();
    break;
  case SPECIALOPTIONS:
    value=model->specialOptions();
    break;
#ifndef COIN_HAS_CBC
  case THREADS:
    value = model->numberThreads();
#endif
  default:
    value=intValue_;
    break;
  }
  return value;
}
#endif
int
CbcOrClpParam::checkDoubleParameter (double value) const
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    return 0;
  }
}
#ifdef COIN_HAS_CBC
double 
CbcOrClpParam::doubleParameter (OsiSolverInterface * model) const
{
  double value=0.0;
  switch(type_) {
  case DUALTOLERANCE:
    assert(model->getDblParam(OsiDualTolerance,value));
    break;
  case PRIMALTOLERANCE:
    assert(model->getDblParam(OsiPrimalTolerance,value));
    break;
  default:
    return doubleValue_;
    break;
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (OsiSolverInterface * model,int value) 
{
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    int oldValue=intValue_;
    intValue_=oldValue;
    switch(type_) {
    case SOLVERLOGLEVEL:
      model->messageHandler()->setLogLevel(value);
      break;
    default:
      break;
    }
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (OsiSolverInterface * model) const
{
  int value=0;
  switch(type_) {
  case SOLVERLOGLEVEL:
    value=model->messageHandler()->logLevel();
    break;
  default:
    value=intValue_;
    break;
  }
  return value;
}
int
CbcOrClpParam::setDoubleParameter (CbcModel &model,double value) 
{
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    double oldValue=doubleValue_;
    doubleValue_ = value;
    switch(type_) {
    case INFEASIBILITYWEIGHT:
      oldValue=model.getDblParam(CbcModel::CbcInfeasibilityWeight);
      model.setDblParam(CbcModel::CbcInfeasibilityWeight,value);
      break;
    case INTEGERTOLERANCE:
      oldValue=model.getDblParam(CbcModel::CbcIntegerTolerance);
      model.setDblParam(CbcModel::CbcIntegerTolerance,value);
      break;
    case INCREMENT:
      oldValue=model.getDblParam(CbcModel::CbcCutoffIncrement);
      model.setDblParam(CbcModel::CbcCutoffIncrement,value);
    case ALLOWABLEGAP:
      oldValue=model.getDblParam(CbcModel::CbcAllowableGap);
      model.setDblParam(CbcModel::CbcAllowableGap,value);
      break;
    case CUTOFF:
      oldValue=model.getCutoff();
      model.setCutoff(value);
      break;
    case TIMELIMIT_BAB:
      oldValue = model.getDblParam(CbcModel::CbcMaximumSeconds) ;
      {
        //OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
        //ClpSimplex * lpSolver = clpSolver->getModelPtr();
        //lpSolver->setMaximumSeconds(value);
        model.setDblParam(CbcModel::CbcMaximumSeconds,value) ;
      }
      break ;
    case DUALTOLERANCE:
    case PRIMALTOLERANCE:
      setDoubleParameter(model.solver(),value);
      return 0; // to avoid message
    default:
      break;
    }
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    return 0;
  }
}
double 
CbcOrClpParam::doubleParameter (CbcModel &model) const
{
  double value;
  switch(type_) {
  case INFEASIBILITYWEIGHT:
    value=model.getDblParam(CbcModel::CbcInfeasibilityWeight);
    break;
  case INTEGERTOLERANCE:
    value=model.getDblParam(CbcModel::CbcIntegerTolerance);
    break;
  case INCREMENT:
    value=model.getDblParam(CbcModel::CbcCutoffIncrement);
  case ALLOWABLEGAP:
    value=model.getDblParam(CbcModel::CbcAllowableGap);
    break;
  case CUTOFF:
    value=model.getCutoff();
    break;
  case TIMELIMIT_BAB:
    value = model.getDblParam(CbcModel::CbcMaximumSeconds) ;
    break ;
  case DUALTOLERANCE:
  case PRIMALTOLERANCE:
    value=doubleParameter(model.solver());
    break;
  default:
    value = doubleValue_;
    break;
  }
  return value;
}
int 
CbcOrClpParam::setIntParameter (CbcModel &model,int value) 
{
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    int oldValue=intValue_;
    intValue_ = value;
    switch(type_) {
    case LOGLEVEL:
      oldValue = model.messageHandler()->logLevel();
      model.messageHandler()->setLogLevel(CoinAbs(value));
      break;
    case SOLVERLOGLEVEL:
      oldValue = model.solver()->messageHandler()->logLevel();
      model.solver()->messageHandler()->setLogLevel(value);
      break;
    case MAXNODES:
      oldValue=model.getIntParam(CbcModel::CbcMaxNumNode);
      model.setIntParam(CbcModel::CbcMaxNumNode,value);
      break;
    case MAXSOLS:
      oldValue=model.getIntParam(CbcModel::CbcMaxNumSol);
      model.setIntParam(CbcModel::CbcMaxNumSol,value);
      break;
    case STRONGBRANCHING:
      oldValue=model.numberStrong();
      model.setNumberStrong(value);
      break;
    case NUMBERBEFORE:
      oldValue=model.numberBeforeTrust();
      model.setNumberBeforeTrust(value);
      break;
    case NUMBERANALYZE:
      oldValue=model.numberAnalyzeIterations();
      model.setNumberAnalyzeIterations(value);
      break;
    case NUMBERMINI:
      oldValue=model.sizeMiniTree();
      model.setSizeMiniTree(value);
      break;
    case CUTPASSINTREE:
      oldValue=model.getMaximumCutPasses();
      model.setMaximumCutPasses(value);
      break;
    case CUTPASS:
      oldValue=model.getMaximumCutPassesAtRoot();
      model.setMaximumCutPassesAtRoot(value);
      break;
#ifdef COIN_HAS_CBC
#ifdef CBC_THREAD
    case THREADS:
      oldValue=model.getNumberThreads();
      model.setNumberThreads(value);
      break;
#endif
#endif
    default:
      break;
    }
    if (doPrinting)
      std::cout<<name_<<" was changed from "<<oldValue<<" to "
               <<value<<std::endl;
    return 0;
  }
}
int 
CbcOrClpParam::intParameter (CbcModel &model) const
{
  int value;
  switch(type_) {
  case LOGLEVEL:
    value = model.messageHandler()->logLevel();
      break;
  case SOLVERLOGLEVEL:
    value = model.solver()->messageHandler()->logLevel();
      break;
  case MAXNODES:
    value = model.getIntParam(CbcModel::CbcMaxNumNode);
    break;
  case MAXSOLS:
    value = model.getIntParam(CbcModel::CbcMaxNumSol);
    break;
  case STRONGBRANCHING:
    value=model.numberStrong();
    break;
  case NUMBERBEFORE:
    value=model.numberBeforeTrust();
    break;
  case NUMBERANALYZE:
    value=model.numberAnalyzeIterations();
    break;
  case NUMBERMINI:
    value=model.sizeMiniTree();
    break;
  case CUTPASSINTREE:
    value=model.getMaximumCutPasses();
    break;
  case CUTPASS:
    value=model.getMaximumCutPassesAtRoot();
    break;
#ifdef COIN_HAS_CBC
#ifdef CBC_THREAD
  case THREADS:
    value = model.getNumberThreads();
#endif
#endif
  default:
    value=intValue_;
    break;
  }
  return value;
}
#endif
// Sets current parameter option using string
void 
CbcOrClpParam::setCurrentOption ( const std::string value )
{
  int action = parameterOption(value);
  if (action>=0)
    currentKeyWord_=action;
}
// Sets current parameter option
void 
CbcOrClpParam::setCurrentOption ( int value , bool printIt)
{
  if (printIt&&value!=currentKeyWord_)
    std::cout<<"Option for "<<name_<<" changed from "
             <<definedKeyWords_[currentKeyWord_]<<" to "
             <<definedKeyWords_[value]<<std::endl;

    currentKeyWord_=value;
}
void 
CbcOrClpParam::setIntValue ( int value )
{ 
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
  } else {
    intValue_=value;
  }
}
void 
CbcOrClpParam::setDoubleValue ( double value )
{ 
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
  } else {
    doubleValue_=value;
  }
}
void 
CbcOrClpParam::setStringValue ( std::string value )
{ 
  stringValue_=value;
}
static char line[1000];
static char * where=NULL;
extern int CbcOrClpRead_mode;
extern FILE * CbcOrClpReadCommand;
// Simple read stuff
std::string
CoinReadNextField()
{
  std::string field;
  if (!where) {
    // need new line
#ifdef COIN_HAS_READLINE     
    if (CbcOrClpReadCommand==stdin) {
      // Get a line from the user. 
      where = readline (coin_prompt);
      
      // If the line has any text in it, save it on the history.
      if (where) {
	if ( *where)
	  add_history (where);
	strcpy(line,where);
	free(where);
      }
    } else {
      where = fgets(line,1000,CbcOrClpReadCommand);
    }
#else
    if (CbcOrClpReadCommand==stdin) {
      fprintf(stdout,coin_prompt);
      fflush(stdout);
    }
    where = fgets(line,1000,CbcOrClpReadCommand);
#endif
    if (!where)
      return field; // EOF
    where = line;
    // clean image
    char * lastNonBlank = line-1;
    while ( *where != '\0' ) {
      if ( *where != '\t' && *where < ' ' ) {
	break;
      } else if ( *where != '\t' && *where != ' ') {
	lastNonBlank = where;
      }
      where++;
    }
    where=line;
    *(lastNonBlank+1)='\0';
  }
  // munch white space
  while(*where==' '||*where=='\t')
    where++;
  char * saveWhere = where;
  while (*where!=' '&&*where!='\t'&&*where!='\0')
    where++;
  if (where!=saveWhere) {
    char save = *where;
    *where='\0';
    //convert to string
    field=saveWhere;
    *where=save;
  } else {
    where=NULL;
    field="EOL";
  }
  return field;
}

std::string
CoinReadGetCommand(int argc, const char *argv[])
{
  std::string field="EOL";
  // say no =
  afterEquals="";
  while (field=="EOL") {
    if (CbcOrClpRead_mode>0) {
      if (CbcOrClpRead_mode<argc) {
	field = argv[CbcOrClpRead_mode++];
	if (field=="-") {
	  std::cout<<"Switching to line mode"<<std::endl;
	  CbcOrClpRead_mode=-1;
	  field=CoinReadNextField();
	} else if (field[0]!='-') {
	  if (CbcOrClpRead_mode!=2) {
	    // now allow std::cout<<"skipping non-command "<<field<<std::endl;
	    // field="EOL"; // skip
	  } else {
	    // special dispensation - taken as -import name
	    CbcOrClpRead_mode--;
	    field="import";
	  }
	} else {
	  if (field!="--") {
	    // take off -
	    field = field.substr(1);
	  } else {
	    // special dispensation - taken as -import --
	    CbcOrClpRead_mode--;
	    field="import";
	  }
	}
      } else {
	field="";
      }
    } else {
      field=CoinReadNextField();
    }
  }
  // if = then modify and save
  std::string::size_type found = field.find('=');
  if (found!=std::string::npos) {
    afterEquals = field.substr(found+1);
    field = field.substr(0,found);
  }
  //std::cout<<field<<std::endl;
  return field;
}
std::string
CoinReadGetString(int argc, const char *argv[])
{
  std::string field="EOL";
  if (afterEquals=="") {
    if (CbcOrClpRead_mode>0) {
      if (CbcOrClpRead_mode<argc) {
        if (argv[CbcOrClpRead_mode][0]!='-') { 
          field = argv[CbcOrClpRead_mode++];
        } else if (!strcmp(argv[CbcOrClpRead_mode],"--")) {
          field = argv[CbcOrClpRead_mode++];
          // -- means import from stdin
          field = "-";
        }
      }
    } else {
      field=CoinReadNextField();
    }
  } else {
    field=afterEquals;
    afterEquals = "";
  }
  //std::cout<<field<<std::endl;
  return field;
}
// valid 0 - okay, 1 bad, 2 not there
int
CoinReadGetIntField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (afterEquals=="") {
    if (CbcOrClpRead_mode>0) {
      if (CbcOrClpRead_mode<argc) {
        // may be negative value so do not check for -
        field = argv[CbcOrClpRead_mode++];
      }
    } else {
      field=CoinReadNextField();
    }
  } else {
    field=afterEquals;
    afterEquals = "";
  }
  int value=0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    const char * start = field.c_str();
    char * endPointer = NULL;
    // check valid
    value =  strtol(start,&endPointer,10);
    if (*endPointer=='\0') {
      *valid = 0;
    } else {
      *valid = 1;
      std::cout<<"String of "<<field;
    }
  } else {
    *valid=2;
  }
  return value;
}
double
CoinReadGetDoubleField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (afterEquals=="") {
    if (CbcOrClpRead_mode>0) {
      if (CbcOrClpRead_mode<argc) {
        // may be negative value so do not check for -
        field = argv[CbcOrClpRead_mode++];
      }
    } else {
      field=CoinReadNextField();
    }
  } else {
    field=afterEquals;
    afterEquals = "";
  }
  double value=0.0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    const char * start = field.c_str();
    char * endPointer = NULL;
    // check valid
    value =  strtod(start,&endPointer);
    if (*endPointer=='\0') {
      *valid = 0;
    } else {
      *valid = 1;
      std::cout<<"String of "<<field;
    }
  } else {
    *valid=2;
  }
  return value;
}
/*
  Subroutine to establish the cbc parameter array. See the description of
  class CbcOrClpParam for details. Pulled from C..Main() for clarity. 
*/
void 
establishParams (int &numberParameters, CbcOrClpParam *const parameters)
{
  numberParameters=0;
  parameters[numberParameters++]=
    CbcOrClpParam("?","For help",GENERALQUERY,7,false);
  parameters[numberParameters++]=
    CbcOrClpParam("???","For help",FULLGENERALQUERY,7,false);
  parameters[numberParameters++]=
    CbcOrClpParam("-","From stdin",
		  STDIN,3,false);
#ifdef COIN_HAS_CBC
    parameters[numberParameters++]=
      CbcOrClpParam("allow!ableGap","Stop when gap between best possible and \
best less than this",
	      0.0,1.0e20,ALLOWABLEGAP);
  parameters[numberParameters-1].setDoubleValue(0.0);
  parameters[numberParameters-1].setLonghelp
    (
     "If the gap between best solution and best possible solution is less than this \
then the search will be terminated.  Also see ratioGap."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("allS!lack","Set basis back to all slack and reset solution",
		  ALLSLACK,3);
  parameters[numberParameters-1].setLonghelp
    (
     "Mainly useful for tuning purposes.  Normally the first dual or primal will be using an all slack \
basis anyway."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("auto!Scale","Whether to scale objective, rhs and bounds of problem if they look odd",
		  "off",AUTOSCALE,7,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "If you think you may get odd objective values or large equality rows etc then\
 it may be worth setting this true.  It is still experimental and you may prefer\
 to use objective!Scale and rhs!Scale"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("barr!ier","Solve using primal dual predictor corrector algorithm",
		  BARRIER);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the  primal dual predictor \
corrector algorithm.  You will probably want to link in and choose a better \
ordering and factorization than the default ones provided.  It will also solve models \
with quadratic objectives."
     
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("basisI!n","Import basis from bas file",
		  BASISIN,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read an MPS format basis file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libz then it can read compressed\
 files 'xxxxxxxx.gz'.."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("basisO!ut","Export basis as bas file",
		  BASISOUT);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write an MPS format basis file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.bas'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("biasLU","Whether factorization biased towards U",
		  "UU",BIASLU,2,false);
  parameters[numberParameters-1].append("UX");
  parameters[numberParameters-1].append("LX");
  parameters[numberParameters-1].append("LL");
  parameters[numberParameters-1].setCurrentOption("LX");
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("branch!AndCut","Do Branch and Cut",
		  BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "This does branch and cut.  There are many parameters which can affect the performance.  \
First just try with default settings and look carefully at the log file.  Did cuts help?  Did they take too long?  \
Look at output to see which cuts were effective and then do some tuning.  You will see that the \
options for cuts are off, on, root and ifmove, forceon.  Off is \
obvious, on means that this cut generator will be tried in the branch and cut tree (you can fine tune using \
'depth').  Root means just at the root node while 'ifmove' means that cuts will be used in the tree if they \
look as if they are doing some good and moving the objective value.  Forceon is same as on but forces code to use \
cut generator at every node.  For probing forceonbut just does fixing probing in tree - not strengthening etc.  \
If pre-processing reduced the size of the \
problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics \
which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions \
were obtained fairly early in the search so the important point is to select the best variable to branch on.  \
See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching \
and trustPseudoCosts parameters."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("bscale","Whether to scale in barrier",
		  "off",BARRIERSCALE,7,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters++]=
    CbcOrClpParam("chol!esky","Which cholesky algorithm",
		  "native",CHOLESKY,7);
  parameters[numberParameters-1].append("dense");
  //#ifdef FOREIGN_BARRIER
#ifdef WSSMP_BARRIER
  parameters[numberParameters-1].append("fudge!Long");
  parameters[numberParameters-1].append("wssmp");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("fudge!Long_dummy");
  parameters[numberParameters-1].append("wssmp_dummy");
#endif
#ifdef UFL_BARRIER
  parameters[numberParameters-1].append("Uni!versityOfFlorida");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("Uni!versityOfFlorida_dummy");    
#endif
#ifdef TAUCS_BARRIER
  parameters[numberParameters-1].append("Taucs");
#define REAL_BARRIER
#else
  parameters[numberParameters-1].append("Taucs_dummy");
#endif
  parameters[numberParameters-1].setLonghelp
    (
     "For a barrier code to be effective it needs a good Cholesky ordering and factorization.  \
You will need to link in one from another source.  See Makefile.locations for some \
possibilities."
     ); 
  //#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("clique!Cuts","Whether to use Clique cuts",
		  "off",CLIQUECUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on clique cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("combine!Solutions","Whether to use combine solution heuristic",
		  "off",COMBINE);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("do");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on a heuristic which does branch and cut on the problem given by just \
using variables which have appeared in one or more solutions. \
It obviously only tries after two or more solutions. \
The Do option switches on before preprocessing."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("cost!Strategy","How to use costs as priorities",
		  "off",COSTSTRATEGY);
  parameters[numberParameters-1].append("pri!orities");
  parameters[numberParameters-1].append("column!Order?");
  parameters[numberParameters-1].append("01f!irst?");
  parameters[numberParameters-1].append("01l!ast?");
  parameters[numberParameters-1].append("length!?");
  parameters[numberParameters-1].setLonghelp
    (
     "This orders the variables in order of their absolute costs - with largest cost ones being branched on \
first.  This primitive strategy can be surprsingly effective.  The column order\
 option is obviously not on costs but easy to code here."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("cpp!Generate","Generates C++ code",
		  -1,50000,CPP);
  parameters[numberParameters-1].setLonghelp
    (
     "Once you like what the stand-alone solver does then this allows \
you to generate user_driver.cpp which approximates the code.  \
0 gives simplest driver, 1 generates saves and restores, 2 \
generates saves and restores even for variables at default value. \
4 bit in cbc generates size dependent code rather than computed values."
     );
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("crash","Whether to create basis for problem",
		  "off",CRASH);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("so!low_halim");
  parameters[numberParameters-1].append("ha!lim_solow(JJF mods)");
  //  parameters[numberParameters-1].append("4");
  //  parameters[numberParameters-1].append("5");
  parameters[numberParameters-1].setLonghelp
    (
     "If crash is set on and there is an all slack basis then Clp will flip or put structural\
 variables into basis with the aim of getting dual feasible.  On the whole dual seems to be\
 better without it and there alernative types of 'crash' for primal e.g. 'idiot' or 'sprint'. \
I have also added a variant due to Solow and Halim which is as on but just flip."); 
  parameters[numberParameters++]=
    CbcOrClpParam("cross!over","Whether to get a basic solution after barrier",
		  "on",CROSSOVER);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].append("maybe");
  parameters[numberParameters-1].append("presolve");
  parameters[numberParameters-1].setLonghelp
    (
     "Interior point algorithms do not obtain a basic solution (and \
the feasibility criterion is a bit suspect (JJF)).  This option will crossover \
to a basic solution suitable for ranging or branch and cut.  With the current state \
of quadratic it may be a good idea to switch off crossover for quadratic (and maybe \
presolve as well) - the option maybe does this."
     );
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("cutD!epth","Depth in tree at which to do cuts",
		  -1,999999,CUTDEPTH);
  parameters[numberParameters-1].setLonghelp
    (
     "Cut generators may be - off, on only at root, on if they look possible \
and on.  If they are done every node then that is that, but it may be worth doing them \
every so often.  The original method was every so many nodes but it is more logical \
to do it whenever depth in tree is a multiple of K.  This option does that and defaults \
to -1 (off)."
     );
  parameters[numberParameters-1].setIntValue(-1);
  parameters[numberParameters++]=
    CbcOrClpParam("cuto!ff","All solutions must be better than this",
		  -1.0e60,1.0e60,CUTOFF);
  parameters[numberParameters-1].setDoubleValue(1.0e50);
  parameters[numberParameters-1].setLonghelp
    (
     "All solutions must be better than this value (in a minimization sense).  \
This is also set by code whenever it obtains a solution and is set to value of \
objective for solution minus cutoff increment."
     );
  parameters[numberParameters++]=
    CbcOrClpParam("cuts!OnOff","Switches all cuts on or off",
		  "off",CUTSSTRATEGY);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "This can be used to switch on or off all cuts (apart from Reduce and Split).  Then you can do \
individual ones off or on \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("debug!In","read valid solution from file",
		  DEBUG,7,false);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read a solution file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.\n\n\
If set to create it will create a file called debug.file  after search; if set \
to createAfterPre it will create one suitable for use after preprocessing.\n\n\
The idea is that if you suspect a bad cut generator and you did not use preprocessing \
you can do a good run with debug set to 'create' and then switch on the cuts you suspect and \
re-run with debug set to 'debug.file'  Similarly if you do use preprocessing but use \
createAfterPre.  The create case has same effect as saveSolution."
     ); 
#endif 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("dextra1","Extra double parameter 1",
		  -COIN_DBL_MAX,COIN_DBL_MAX,DEXTRA1,false);
  parameters[numberParameters-1].setDoubleValue(0.0);
  parameters[numberParameters++]=
    CbcOrClpParam("dextra2","Extra double parameter 2",
		  -COIN_DBL_MAX,COIN_DBL_MAX,DEXTRA2,false);
  parameters[numberParameters-1].setDoubleValue(0.0);
  parameters[numberParameters++]=
    CbcOrClpParam("dextra3","Extra double parameter 3",
		  -COIN_DBL_MAX,COIN_DBL_MAX,DEXTRA3,false);
  parameters[numberParameters-1].setDoubleValue(0.0);
  parameters[numberParameters++]=
    CbcOrClpParam("dextra4","Extra double parameter 4",
		  -COIN_DBL_MAX,COIN_DBL_MAX,DEXTRA4,false);
  parameters[numberParameters-1].setDoubleValue(0.0);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("direction","Minimize or Maximize",
		  "min!imize",DIRECTION);
  parameters[numberParameters-1].append("max!imize");
  parameters[numberParameters-1].append("zero");
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'direction maximize' for maximization.\n\
You can also use the parameters 'maximize' or 'minimize'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("directory","Set Default directory for import etc.",
		  DIRECTORY);
  parameters[numberParameters-1].setLonghelp
    (
     "This sets the directory which import, export, saveModel, restoreModel etc will use.\
  It is initialized to './'"
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("doH!euristic","Do heuristics before any preprocessing",
		  DOHEURISTIC,3);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally heuristics are done in branch and bound.  It may be useful to do them outside. \
Doing this may also set cutoff."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("dualB!ound","Initially algorithm acts as if no \
gap between bounds exceeds this value",
		  1.0e-20,1.0e12,DUALBOUND);
  parameters[numberParameters-1].setLonghelp
    (
     "The dual algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  If a problem has both upper and\
 lower bounds then it is trivial to get dual feasible by setting non basic variables\
 to correct bound.  If the gap between the upper and lower bounds of a variable is more\
 than the value of dualBound Clp introduces fake bounds so that it can make the problem\
 dual feasible.  This has the same effect as a composite objective function in the\
 primal algorithm.  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the bounds.  OSL had a heuristic to\
 adjust bounds, maybe we need that here."
     );
  parameters[numberParameters++]=
    CbcOrClpParam("dualize","Solves dual reformulation",
		  0,2,DUALIZE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Don't even think about it."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("dualP!ivot","Dual pivot choice algorithm",
		  "auto!matic",DUALPIVOT);
  parameters[numberParameters-1].append("dant!zig");
  parameters[numberParameters-1].append("partial");
  parameters[numberParameters-1].append("steep!est");
  parameters[numberParameters-1].setLonghelp
    (
     "Clp can use any pivot selection algorithm which the user codes as long as it\
 implements the features in the abstract pivot base class.  The Dantzig method is implemented\
 to show a simple method but its use is deprecated.  Steepest is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while automatic decides at each iteration based on information\
 about the factorization."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("dualS!implex","Do dual simplex algorithm",
		  DUALSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the dual steepest edge algorithm.\
The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by dual pivot method, fake bound on variables and dual and primal tolerances."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("dualT!olerance","For an optimal solution \
no dual infeasibility may exceed this value",
		  1.0e-20,1.0e12,DUALTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally the default tolerance is fine, but you may want to increase it a\
 bit if a dual run seems to be having a hard time.  One method which can be faster is \
to use a large tolerance e.g. 1.0e-4 and dual and then clean up problem using primal and the \
correct tolerance (remembering to switch off presolve for this final short clean up phase)."
     ); 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("either!Simplex","Do dual or primal simplex algorithm",
		  EITHERSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the dual or primal algorithm,\
 based on a dubious analysis of model."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("end","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops execution ; end, exit, quit and stop are synonyms"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("error!sAllowed","Whether to allow import errors",
		  "off",ERRORSALLOWED,3);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "The default is not to use any model which had errors when reading the mps file.\
  Setting this to 'on' will allow all errors from which the code can recover\
 simply by ignoring the error.  There are some errors from which the code can not recover \
e.g. no ENDATA.  This has to be set before import i.e. -errorsAllowed on -import xxxxxx.mps."
     );
  parameters[numberParameters++]=
    CbcOrClpParam("exit","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("export","Export model as mps file",
		  EXPORT);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write an MPS format file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.mps'."
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("extra1","Extra integer parameter 1",
		  -1,COIN_INT_MAX,EXTRA1,false);
  parameters[numberParameters-1].setIntValue(-1);
  parameters[numberParameters++]=
    CbcOrClpParam("extra2","Extra integer parameter 2",
		  -1,COIN_INT_MAX,EXTRA2,false);
  parameters[numberParameters-1].setIntValue(-1);
  parameters[numberParameters++]=
    CbcOrClpParam("extra3","Extra integer parameter 3",
		  -1,COIN_INT_MAX,EXTRA3,false);
  parameters[numberParameters-1].setIntValue(-1);
  parameters[numberParameters++]=
    CbcOrClpParam("extra4","Extra integer parameter 4",
		  -1,COIN_INT_MAX,EXTRA4,false);
  parameters[numberParameters-1].setIntValue(-1);
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("fakeB!ound","All bounds <= this value - DEBUG",
		  1.0,1.0e15,FAKEBOUND,false);
#ifdef COIN_HAS_CBC
    parameters[numberParameters++]=
      CbcOrClpParam("feas!ibilityPump","Whether to try Feasibility Pump",
		    "off",FPUMP);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("do");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on feasibility pump heuristic at root. This is due to Fischetti and Lodi \
and uses a sequence of Lps to try and get an integer feasible solution. \
Some fine tuning is available by passFeasibilityPump. Do options does heuristic before preprocessing"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("fix!OnDj","Try heuristic based on fixing variables with \
reduced costs greater than this",
		  -1.0e20,1.0e20,DJFIX,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If this is set integer variables with reduced costs greater than this will be fixed \
before branch and bound - use with extreme caution!" 
     ); 
    parameters[numberParameters++]=
      CbcOrClpParam("flow!CoverCuts","Whether to use Flow Cover cuts",
		    "off",FLOWCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters-1].append("ifmove");
    parameters[numberParameters-1].append("forceOn");
    parameters[numberParameters-1].setLonghelp
    (
     "This switches on flow cover cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
    parameters[numberParameters++]=
      CbcOrClpParam("force!Solution","Whether to use given solution as crash for BAB",
		    "off",USESOLUTION);
    parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "If on then tries to branch to solution given by AMPL or priorities file."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("gamma!(Delta)","Whether to regularize barrier",
		  "off",GAMMA,7,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("gamma");
  parameters[numberParameters-1].append("delta");
  parameters[numberParameters-1].append("onstrong");
  parameters[numberParameters-1].append("gammastrong");
  parameters[numberParameters-1].append("deltastrong");
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("gomory!Cuts","Whether to use Gomory cuts",
		  "off",GOMORYCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "The original cuts - beware of imitations!  Having gone out of favor, they are now more \
fashionable as LP solvers are more robust and they interact well with other cuts.  They will almost always \
give cuts (although in this executable they are limited as to number of variables in cut).  \
However the cuts may be dense so it is worth experimenting. \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("greedy!Heuristic","Whether to use a greedy heuristic",
		  "off",GREEDY);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("do");
  //parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].setLonghelp
    (
     "Switches on a greedy heuristic which will try and obtain a solution.  It may just fix a \
percentage of variables and then try a small branch and cut run. \
The Do option switches on before preprocessing."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("heur!isticsOnOff","Switches most heuristics on or off",
		  "off",HEURISTICSTRATEGY);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "This can be used to switch on or off all heuristics.  Then you can do \
individual ones off or on.  CbcTreeLocal is not included as it dramatically \
alters search."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("help","Print out version, non-standard options and some help",
		  HELP,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This prints out some help to get user started.  If you have printed this then \
you should be past that stage:-)"
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("hot!StartMaxIts","Maximum iterations on hot start",
		  0,COIN_INT_MAX,MAXHOTITS,false);
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("idiot!Crash","Whether to try idiot crash",
		  -1,999999,IDIOT);
  parameters[numberParameters-1].setLonghelp
    (
     "This is a type of 'crash' which works well on some homogeneous problems.\
 It works best on problems with unit elements and rhs but will do something to any model.  It should only be\
 used before primal.  It can be set to -1 when the code decides for itself whether to use it,\
 0 to switch off or n > 0 to do n passes."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("import","Import model from mps file",
		  IMPORT,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read an MPS format file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz'.."
     );
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("inc!rement","A valid solution must be at least this \
much better than last integer solution",
		  -1.0e20,1.0e20,INCREMENT);
  parameters[numberParameters-1].setLonghelp
    (
     "Whenever a solution is found the bound on solutions is set to solution (in a minimization\
sense) plus this.  If it is not set then the code will try and work one out e.g. if \
all objective coefficients are multiples of 0.01 and only integer variables have entries in \
objective then this can be set to 0.01.  Be careful if you set this negative!"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("inf!easibilityWeight","Each integer infeasibility is expected \
to cost this much",
		  0.0,1.0e20,INFEASIBILITYWEIGHT);
  parameters[numberParameters-1].setLonghelp
    (
     "A primitive way of deciding which node to explore next.  Satisfying each integer infeasibility is \
expected to cost this much."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("initialS!olve","Solve to continuous",
		  SOLVECONTINUOUS);
  parameters[numberParameters-1].setLonghelp
    (
     "This just solves the problem to continuous - without adding any cuts"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("integerT!olerance","For an optimal solution \
no integer variable may be this away from an integer value",
	      1.0e-20,0.5,INTEGERTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "Beware of setting this smaller than the primal tolerance."
     ); 
#endif 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("keepN!ames","Whether to keep names from import",
		  "on",KEEPNAMES);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].setLonghelp
    (
     "It saves space to get rid of names so if you need to you can set this to off.  \
This needs to be set before the import of model - so -keepnames off -import xxxxx.mps."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("KKT","Whether to use KKT factorization",
		  "off",KKT,7,false);
  parameters[numberParameters-1].append("on");
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("knapsack!Cuts","Whether to use Knapsack cuts",
		  "off",KNAPSACKCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on knapsack cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("lift!AndProjectCuts","Whether to use Lift and Project cuts",
		  "off",LANDPCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "Lift and project cuts - may be expensive to compute. \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("local!TreeSearch","Whether to use local treesearch",
		  "off",LOCALTREE);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on a local search algorithm when a solution is found.  This is from \
Fischetti and Lodi and is not really a heuristic although it can be used as one. \
When used from Coin solve it has limited functionality.  It is not switched on when \
heuristics are switched on."
     ); 
#endif
#ifndef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("log!Level","Level of detail in Solver output",
		  -1,63,SOLVERLOGLEVEL);
#else
  parameters[numberParameters++]=
    CbcOrClpParam("log!Level","Level of detail in Coin branch and Cut output",
		  -63,63,LOGLEVEL);
  parameters[numberParameters-1].setIntValue(1);
#endif
  parameters[numberParameters-1].setLonghelp
    (
     "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("max!imize","Set optimization direction to maximize",
		  MAXIMIZE,7);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'maximize' for maximization.\n\
You can also use the parameters 'direction maximize'."
     ); 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("maxF!actor","Maximum number of iterations between \
refactorizations",
		  1,999999,MAXFACTOR);
  parameters[numberParameters-1].setLonghelp
    (
     "If this is at its initial value of 200 then in this executable clp will guess at a\
 value to use.  Otherwise the user can set a value.  The code may decide to re-factorize\
 earlier for accuracy."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("maxIt!erations","Maximum number of iterations before \
stopping",
		  0,2147483647,MAXITERATION);
  parameters[numberParameters-1].setLonghelp
    (
     "This can be used for testing purposes.  The corresponding library call\n\
      \tsetMaximumIterations(value)\n can be useful.  If the code stops on\
 seconds or by an interrupt this will be treated as stopping on maximum iterations"
     ); 
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("maxN!odes","Maximum number of nodes to do",
		  0,2147483647,MAXNODES);
  parameters[numberParameters-1].setLonghelp
    (
     "This is a repeatable way to limit search.  Normally using time is easier \
but then the results may not be repeatable."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("maxS!olutions","Maximum number of solutions to get",
		  1,2147483647,MAXSOLS);
  parameters[numberParameters-1].setLonghelp
    (
     "You may want to stop after (say) two solutions or an hour."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("min!imize","Set optimization direction to minimize",
		  MINIMIZE,7);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is minimize - use 'maximize' for maximization.\n\
This should only be necessary if you have previously set maximization \
You can also use the parameters 'direction minimize'."
     );
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("mipO!ptions","Dubious options for mip",
		  0,COIN_INT_MAX,MIPOPTIONS,false);
  parameters[numberParameters++]=
    CbcOrClpParam("more!MipOptions","More dubious options for mip",
		  -1,COIN_INT_MAX,MOREMIPOPTIONS,false);
  parameters[numberParameters++]=
    CbcOrClpParam("mixed!IntegerRoundingCuts","Whether to use Mixed Integer Rounding cuts",
		  "off",MIXEDCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on mixed integer rounding cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("mess!ages","Controls if Clpnnnn is printed",
		  "off",MESSAGES);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    ("The default behavior is to put out messages such as:\n\
   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\n\
but this program turns this off to make it look more friendly.  It can be useful\
 to turn them back on if you want to be able to 'grep' for particular messages or if\
 you intend to override the behavior of a particular message."
     );
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("miniT!ree","Size of fast mini tree",
		  0,COIN_INT_MAX,NUMBERMINI,false);
  parameters[numberParameters-1].setLonghelp
    (
     "The idea is that I can do a small tree fast. \
This is a first try and will hopefully become more sophisticated."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("miplib","Do some of miplib test set",
		  MIPLIB,3);
#endif 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("netlib","Solve entire netlib test set",
		  NETLIB_EITHER,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using dual or primal.\
The user can set options before e.g. clp -presolve off -netlib"
     ); 
#ifdef REAL_BARRIER
  parameters[numberParameters++]=
    CbcOrClpParam("netlibB!arrier","Solve entire netlib test set with barrier",
		  NETLIB_BARRIER,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using barrier.\
The user can set options before e.g. clp -kkt on -netlib"
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("netlibD!ual","Solve entire netlib test set (dual)",
		  NETLIB_DUAL,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using dual.\
The user can set options before e.g. clp -presolve off -netlib"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("netlibP!rimal","Solve entire netlib test set (primal)",
		  NETLIB_PRIMAL,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using primal.\
The user can set options before e.g. clp -presolve off -netlibp"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("netlibT!une","Solve entire netlib test set with 'best' algorithm",
		  NETLIB_TUNE,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp and then solves the netlib test set using whatever \
works best.  I know this is cheating but it also stresses the code better by doing a \
mixture of stuff.  The best algorithm was chosen on a Linux ThinkPad using native cholesky \
with University of Florida ordering."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("network","Tries to make network matrix",
		  NETWORK,7,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Clp will go faster if the matrix can be converted to a network.  The matrix\
 operations may be a bit faster with more efficient storage, but the main advantage\
 comes from using a network factorization.  It will probably not be as fast as a \
specialized network code."
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("node!Strategy","What strategy to use to select nodes",
                  "hybrid",NODESTRATEGY);
  parameters[numberParameters-1].append("fewest");
  parameters[numberParameters-1].append("depth");
  parameters[numberParameters-1].append("upfewest");
  parameters[numberParameters-1].append("downfewest");
  parameters[numberParameters-1].append("updepth");
  parameters[numberParameters-1].append("downdepth");
  parameters[numberParameters-1].setLonghelp
    (
     "Normally before a solution the code will choose node with fewest infeasibilities. \
You can choose depth as the criterion.  You can also say if up or down branch must \
be done first (the up down choice will carry on after solution). \
Default has now been changed to hybrid which is breadth first on small depth nodes then fewest."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("numberA!nalyze","Number of analysis iterations",
		  -COIN_INT_MAX,COIN_INT_MAX,NUMBERANALYZE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "This says how many iterations to spend at root node analyzing problem. \
This is a first try and will hopefully become more sophisticated."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("objective!Scale","Scale factor to apply to objective",
		  -1.0e20,1.0e20,OBJSCALE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If the objective function has some very large values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  It is applied after scaling"
     ); 
  parameters[numberParameters-1].setDoubleValue(1.0);
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("outDup!licates","takes duplicate rows etc out of integer model",
		  OUTDUPROWS,7,false);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("output!Format","Which output format to use",
		  1,6,OUTPUTFORMAT);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally export will be done using normal representation for numbers and two values\
 per line.  You may want to do just one per line (for grep or suchlike) and you may wish\
 to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal.\
 otherwise odd values gives one value per line, even two.  Values 1,2 give normal format, 3,4\
 gives greater precision, while 5,6 give IEEE values.  When used for exporting a basis 1 does not save \
values, 2 saves values, 3 with greater accuracy and 4 in IEEE."
     );
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("passC!uts","Number of cut passes at root node",
		  -999999,999999,CUTPASS);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is 100 passes if less than 500 columns, 100 passes (but \
stop if drop small if less than 5000 columns, 20 otherwise"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("passF!easibilityPump","How many passes in feasibility pump",
		  0,10000,FPUMPITS);
  parameters[numberParameters-1].setLonghelp
    (
     "This fine tunes Feasibility Pump by doing more or fewer passes."
     ); 
  parameters[numberParameters-1].setIntValue(20);
#endif 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("passP!resolve","How many passes in presolve",
		  -200,100,PRESOLVEPASS,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally Presolve does 5 passes but you may want to do less to make it\
 more lightweight or do more if improvements are still being made.  As Presolve will return\
 if nothing is being taken out, you should not normally need to use this fine tuning."
     );
#endif 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("passT!reeCuts","Number of cut passes in tree",
		  -999999,999999,CUTPASSINTREE);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is one pass"
     ); 
#endif 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("pertV!alue","Method of perturbation",
		  -5000,102,PERTVALUE,false);
  parameters[numberParameters++]=
    CbcOrClpParam("perturb!ation","Whether to perturb problem",
		  "on",PERTURBATION);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].setLonghelp
    (
     "Perturbation helps to stop cycling, but Clp uses other measures for this.\
  However large problems and especially ones with unit elements and unit rhs or costs\
 benefit from perturbation.  Normally Clp tries to be intelligent, but you can switch this off.\
  The Clp library has this off by default.  This program has it on by default."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("PFI","Whether to use Product Form of Inverse in simplex",
		  "off",PFI,7,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "By default clp uses Forrest-Tomlin L-U update.  If you are masochistic you can switch it off."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("plus!Minus","Tries to make +- 1 matrix",
		  PLUSMINUS,7,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Clp will go slightly faster if the matrix can be converted so that the elements are\
 not stored and are known to be unit.  The main advantage is memory use.  Clp may automatically\
 see if it can convert the problem so you should not need to use this."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("pO!ptions","Dubious print options",
		  0,COIN_INT_MAX,PRINTOPTIONS,false);
  parameters[numberParameters-1].setIntValue(0);
  parameters[numberParameters-1].setLonghelp
    (
     "If this is > 0 then presolve will give more information and branch and cut will give statistics"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("preO!pt","Presolve options",
		  0,COIN_INT_MAX,PRESOLVEOPTIONS,false);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("presolve","Whether to presolve problem",
		  "on",PRESOLVE);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters-1].append("more");
  parameters[numberParameters-1].append("file");
  parameters[numberParameters-1].setLonghelp
    (
     "Presolve analyzes the model to find such things as redundant equations, equations\
 which fix some variables, equations which can be transformed into bounds etc etc.  For the\
 initial solve of any problem this is worth doing unless you know that it will have no effect.  \
on will normally do 5 passes while using 'more' will do 10.  If the problem is very large you may need \
to write the original to file using 'file'."
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("preprocess","Whether to use integer preprocessing",
                  "off",PREPROCESS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("save");
  parameters[numberParameters-1].append("equal");
  parameters[numberParameters-1].append("sos");
  parameters[numberParameters-1].append("trysos");
  parameters[numberParameters-1].append("equalall");
  parameters[numberParameters-1].append("strategy");
  parameters[numberParameters-1].setLonghelp
    (
     "This tries to reduce size of model in a similar way to presolve and \
it also tries to strengthen the model - this can be very useful and is worth trying. \
 Save option saves on file presolved.mps.  equal will turn <= cliques into \
==.  sos will create sos sets if all 0-1 in sets (well one extra is allowed) \
and no overlaps.  trysos is same but allows any number extra.  equalall will turn all \
valid inequalities into equalities with integer slacks.  strategy is as \
on but uses CbcStrategy."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("preT!olerance","Tolerance to use in presolve",
		  1.0e-20,1.0e12,PRESOLVETOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "The default is 1.0e-8 - you may wish to try 1.0e-7 if presolve says the problem is \
infeasible and you have awkward numbers and you are sure the problem is really feasible."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("primalP!ivot","Primal pivot choice algorithm",
		  "auto!matic",PRIMALPIVOT);
  parameters[numberParameters-1].append("exa!ct");
  parameters[numberParameters-1].append("dant!zig");
  parameters[numberParameters-1].append("part!ial");
  parameters[numberParameters-1].append("steep!est");
  parameters[numberParameters-1].append("change");
  parameters[numberParameters-1].append("sprint");
  parameters[numberParameters-1].setLonghelp
    (
     "Clp can use any pivot selection algorithm which the user codes as long as it\
 implements the features in the abstract pivot base class.  The Dantzig method is implemented\
 to show a simple method but its use is deprecated.  Exact devex is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while change initially does dantzig until the factorization\
 becomes denser.  This is still a work in progress."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("primalS!implex","Do primal simplex algorithm",
		  PRIMALSIMPLEX);
  parameters[numberParameters-1].setLonghelp
    (
     "This command solves the current model using the primal algorithm.\
  The default is to use exact devex.\
 The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by column selection  method, infeasibility weight and dual and primal tolerances."
     );
#endif 
  parameters[numberParameters++]=
    CbcOrClpParam("primalT!olerance","For an optimal solution \
no primal infeasibility may exceed this value",
		  1.0e-20,1.0e12,PRIMALTOLERANCE);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally the default tolerance is fine, but you may want to increase it a\
 bit if a primal run seems to be having a hard time"
     ); 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("primalW!eight","Initially algorithm acts as if it \
costs this much to be infeasible",
		  1.0e-20,1.0e20,PRIMALWEIGHT);
  parameters[numberParameters-1].setLonghelp
    (
     "The primal algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  So Clp is minimizing this weight times\
 the sum of primal infeasibilities plus the true objective function (in minimization sense).\
  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the weight in order to get feasible.\
  OSL had a heuristic to\
 adjust bounds, maybe we need that here."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("printi!ngOptions","Print options",
		  "normal",INTPRINT,3);
  parameters[numberParameters-1].append("integer");
  parameters[numberParameters-1].append("special");
  parameters[numberParameters-1].append("rows");
  parameters[numberParameters-1].append("all");
  parameters[numberParameters-1].setLonghelp
    (
     "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \n\
integer - nonzero integer column variables\n\
special - in format suitable for OsiRowCutDebugger\n\
rows - nonzero column variables and row activities\n\
all - all column variables and row activities.\n\
\nFor non-integer problems 'integer' and 'special' act like 'normal'.  \
Also see printMask for controlling output."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("printM!ask","Control printing of solution on a  mask",
		  PRINTMASK,3);
  parameters[numberParameters-1].setLonghelp
    (
     "If set then only those names which match mask are printed in a solution. \
'?' matches any character and '*' matches any set of characters. \
 The default is '' i.e. unset so all variables are printed. \
This is only active if model has names."
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("prio!rityIn","Import priorities etc from file",
		  PRIORITYIN,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This will read a file with priorities from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  This can not read from compressed files. \
File is in csv format with allowed headings - name, number, priority, direction, up, down, solution.  Exactly one of\
 name and number must be given."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("probing!Cuts","Whether to use Probing cuts",
		  "off",PROBINGCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].append("forceOnBut");
  parameters[numberParameters-1].append("forceOnStrong");
  parameters[numberParameters-1].append("forceOnButStrong");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on probing cuts (either at root or in entire tree) \
See branchAndCut for information on options. \
but strong options do more probing"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("pumpT!une","Dubious ideas for feasibility pump",
		  0,100000000,FPUMPTUNE);
  parameters[numberParameters-1].setLonghelp
    (
     "This fine tunes Feasibility Pump \n\
\t>=1000000 use as accumulate switch\n\
\t>=1000 use index+1 as number of large loops\n\
\t>=100 use 0.05 objvalue as increment\n\
\t>=10 use +0.1 objvalue for cutoff (add)\n\
\t1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds"
     ); 
  parameters[numberParameters-1].setIntValue(0);
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("quit","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("ratio!Gap","Stop when gap between best possible and \
best less than this fraction of larger of two",
		  0.0,1.0e20,GAPRATIO);
  parameters[numberParameters-1].setDoubleValue(0.0);
  parameters[numberParameters-1].setLonghelp
    (
     "If the gap between best solution and best possible solution is less than this fraction \
of the objective value at the root node then the search will terminate.  See 'allowableGap' for a \
way of using absolute value rather than fraction."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("reallyO!bjectiveScale","Scale factor to apply to objective in place",
		  -1.0e20,1.0e20,OBJSCALE2,false);
  parameters[numberParameters-1].setLonghelp
    (
     "You can set this to -1.0 to test maximization or other to stress code"
     ); 
  parameters[numberParameters-1].setDoubleValue(1.0);
  parameters[numberParameters++]=
    CbcOrClpParam("reallyS!cale","Scales model in place",
		  REALLY_SCALE,7,false);
#endif
#ifdef COIN_HAS_CBC
    parameters[numberParameters++]=
      CbcOrClpParam("reduce!AndSplitCuts","Whether to use Reduce-and-Split cuts",
	      "off",REDSPLITCUTS);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters-1].append("root");
    parameters[numberParameters-1].append("ifmove");
    parameters[numberParameters-1].append("forceOn");
    parameters[numberParameters-1].setLonghelp
    (
     "This switches on reduce and split  cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("residual!CapacityCuts","Whether to use Residual Capacity cuts",
		  "off",RESIDCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "Residual capacity cuts. \
See branchAndCut for information on options."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("restore!Model","Restore model from binary file",
		  RESTORE);
  parameters[numberParameters-1].setLonghelp
    (
     "This reads data save by saveModel from the given file.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("reverse","Reverses sign of objective",
		  REVERSE,7,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Useful for testing if maximization works correctly"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("rhs!Scale","Scale factor to apply to rhs and bounds",
		  -1.0e20,1.0e20,RHSSCALE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If the rhs or bounds have some very large meaningful values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale"
     ); 
  parameters[numberParameters-1].setDoubleValue(1.0);
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
      CbcOrClpParam("Rins","Whether to try Relaxed Induced Neighborhood Search",
		    "off",RINS);
    parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on Relaxed induced neighborhood Search."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("round!ingHeuristic","Whether to use Rounding heuristic",
		  "off",ROUNDING);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("do");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on a simple (but effective) rounding heuristic at each node of tree.  \
The Do option switches on before preprocessing."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("saveM!odel","Save model to binary file",
		  SAVE);
  parameters[numberParameters-1].setLonghelp
    (
     "This will save the problem to the given file name for future use\
 by restoreModel.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("saveS!olution","saves solution to file",
		  SAVESOL);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write a binary solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows \
and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column \
activities and reduced costs - see bottom of CbcOrClpParam.cpp for code that reads or writes file."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("scal!ing","Whether to scale problem",
		  "off",SCALING);
  parameters[numberParameters-1].append("equi!librium");
  parameters[numberParameters-1].append("geo!metric");
  parameters[numberParameters-1].append("auto!matic");
  parameters[numberParameters-1].setLonghelp
    (
     "Scaling can help in solving problems which might otherwise fail because of lack of\
 accuracy.  It can also reduce the number of iterations.  It is not applied if the range\
 of elements is small.  When unscaled it is possible that there may be small primal and/or\
 infeasibilities."
     ); 
  parameters[numberParameters-1].setCurrentOption(3); // say auto
#ifndef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("sec!onds","Maximum seconds",
		  -1.0,1.0e12,TIMELIMIT);
  parameters[numberParameters-1].setLonghelp
    (
     "After this many seconds clp will act as if maximum iterations had been reached \
(if value >=0).  \
In this program it is really only useful for testing but the library function\n\
      \tsetMaximumSeconds(value)\n can be useful."
     );
#else
  parameters[numberParameters++]=
    CbcOrClpParam("sec!onds","maximum seconds",
		  -1.0,1.0e12,TIMELIMIT_BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "After this many seconds coin solver will act as if maximum nodes had been reached."
     );
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("sleep","for debug",
		  DUMMY,7,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If passed to solver fom ampl, then ampl will wait so that you can copy .nl file for debug."
     ); 
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("slp!Value","Number of slp passes before primal",
		  -1,50000,SLPVALUE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If you are solving a quadratic problem using primal then it may be helpful to do some \
sequential Lps to get a good approximate solution."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("solu!tion","Prints solution to file",
		  SOLUTION);
  parameters[numberParameters-1].setLonghelp
    (
     "This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask."
     ); 
#ifdef COIN_HAS_CLP
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("solv!e","Solve problem",
		  BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "If there are no integer variables then this just solves LP.  If there are integer variables \
this does branch and cut."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("sos!Options","Whether to use SOS from AMPL",
		  "off",SOS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setCurrentOption("on");
  parameters[numberParameters-1].setLonghelp
    (
     "Normally if AMPL says there are SOS variables they should be used, but sometime sthey should\
 be turned off - this does so."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("slog!Level","Level of detail in Solver output",
		  -1,63,SOLVERLOGLEVEL);
  parameters[numberParameters-1].setLonghelp
    (
     "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information."
     );
#else
  // allow solve as synonym for dual 
  parameters[numberParameters++]=
    CbcOrClpParam("solv!e","Solve problem using dual simplex",
		  BAB);
  parameters[numberParameters-1].setLonghelp
    (
     "Just so can use solve for clp as well as in cbc"
     ); 
#endif
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("spars!eFactor","Whether factorization treated as sparse",
		  "on",SPARSEFACTOR,7,false);
  parameters[numberParameters-1].append("off");
  parameters[numberParameters++]=
    CbcOrClpParam("special!Options","Dubious options for Simplex - see ClpSimplex.hpp",
		  0,COIN_INT_MAX,SPECIALOPTIONS,false);
  parameters[numberParameters++]=
    CbcOrClpParam("sprint!Crash","Whether to try sprint crash",
		  -1,5000000,SPRINT);
  parameters[numberParameters-1].setLonghelp
    (
     "For long and thin problems this program may solve a series of small problems\
 created by taking a subset of the columns.  I introduced the idea as 'Sprint' after\
 an LP code of that name of the 60's which tried the same tactic (not totally successfully).\
  Cplex calls it 'sifting'.  -1 is automatic choice, 0 is off, n is number of passes"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("stat!istics","Print some statistics",
		  STATISTICS);
  parameters[numberParameters-1].setLonghelp
    (
     "This command prints some statistics for the current model.\
 If log level >1 then more is printed.\
 These are for presolved model if presolve on (and unscaled)."
     );
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("stop","Stops clp execution",
		  EXIT);
  parameters[numberParameters-1].setLonghelp
    (
     "This stops the execution of Clp, end, exit, quit and stop are synonyms"
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("strengthen","Create strengthened problem",
		  STRENGTHEN,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This creates a new problem by applying the root node cuts.  All tight constraints \
will be in resulting problem"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("strong!Branching","Number of variables to look at in strong branching",
		  0,999999,STRONGBRANCHING);
  parameters[numberParameters-1].setLonghelp
    (
     "In order to decide which variable to branch on, the code will choose up to this number \
of unsatisfied variables to do mini up and down branches on.  Then the most effective one is chosen. \
If a variable is branched on many times then the previous average up and down costs may be used - \
see number before trust."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("subs!titution","How long a column to substitute for in presolve",
		  0,10000,SUBSTITUTION,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Normally Presolve gets rid of 'free' variables when there are no more than 3 \
 variables in column.  If you increase this the number of rows may decrease but number of \
 elements may increase."
     ); 
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("testO!si","Test OsiObject stuff",
		  -1,COIN_INT_MAX,TESTOSI,false);
#endif
#ifdef CBC_THREAD
  parameters[numberParameters++]=
    CbcOrClpParam("thread!s","Number of threads to try and use",
		  -100,10000,THREADS,false);
  parameters[numberParameters-1].setLonghelp
    (
     "To use multiple threads, set threads to number wanted.  It may be better \
to use one or two more than number of cpus available.  If 100+n then n threads and \
threads used in sub-trees, if 200+n use threads for root cuts, 300+n - both."
     ); 
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("tighten!Factor","Tighten bounds using this times largest \
activity at continuous solution",
		  1.0e-3,1.0e20,TIGHTENFACTOR,false);
  parameters[numberParameters-1].setLonghelp
    (
     "This sleazy trick can help on some problems."
     ); 
#endif
#ifdef COIN_HAS_CLP
  parameters[numberParameters++]=
    CbcOrClpParam("tightLP","Poor person's preSolve for now",
		  TIGHTEN,7,false);
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("trust!PseudoCosts","Number of branches before we trust pseudocosts",
		  -1,2000000,NUMBERBEFORE);
  parameters[numberParameters-1].setLonghelp
    (
     "Using strong branching computes pseudo-costs.  After this many times for a variable we just \
trust the pseudo costs and do not do any more strong branching."
     ); 
#endif
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("tune!PreProcess","Dubious tuning parameters",
		  0,20000000,PROCESSTUNE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "For making equality cliques this is minimumsize.  Also for adding \
integer slacks.  May be used for more later \
If <1000 that is what it does.  If <1000000 - numberPasses is (value/1000)-1 and tune is tune %1000. \
If >= 1000000! - numberPasses is (value/1000000)-1 and tune is tune %1000000.  In this case if tune is now still >=10000 \
numberPassesPerInnerLoop is changed from 10 to (tune-10000)-1 and tune becomes tune % 10000!!!!! - happy? - \
so to keep normal limit on cliques of 5, do 3 major passes (include presolves) but only doing one tightening pass per major pass - \
you would use 3010005 (I think)"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("two!MirCuts","Whether to use Two phase Mixed Integer Rounding cuts",
		  "off",TWOMIRCUTS);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].append("root");
  parameters[numberParameters-1].append("ifmove");
  parameters[numberParameters-1].append("forceOn");
  parameters[numberParameters-1].setLonghelp
    (
     "This switches on two phase mixed integer rounding  cuts (either at root or in entire tree) \
See branchAndCut for information on options."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("unitTest","Do unit test",
		  UNITTEST,3);
  parameters[numberParameters-1].setLonghelp
    (
     "This exercises the unit test for clp"
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("userClp","Hand coded Clp stuff",
		  USERCLP);
  parameters[numberParameters-1].setLonghelp
    (
     "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCLP in main driver and modify sample code."
     ); 
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("userCbc","Hand coded Cbc stuff",
		  USERCBC);
  parameters[numberParameters-1].setLonghelp
    (
     "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCBC in main driver and modify sample code."
     ); 
#endif
  parameters[numberParameters++]=
    CbcOrClpParam("vector","Whether to use vector? Form of matrix in simplex",
		  "off",VECTOR,7,false);
  parameters[numberParameters-1].append("on");
  parameters[numberParameters-1].setLonghelp
    (
     "If this is on and ClpPackedMatrix uses extra column copy in odd format."
     ); 
  parameters[numberParameters++]=
    CbcOrClpParam("verbose","Switches on longer help on single ?",
		  0,15,VERBOSE,false);
  parameters[numberParameters-1].setLonghelp
    (
     "Set to 1 to get short help with ? list, 2 to get long help, 3 for both.  (add 4 to just get ampl ones)."
     ); 
  parameters[numberParameters-1].setIntValue(0);
#ifdef COIN_HAS_CBC
  parameters[numberParameters++]=
    CbcOrClpParam("vub!heuristic","Type of vub heuristic",
		  -2,10,VUBTRY,false);
  parameters[numberParameters-1].setLonghelp
    (
     "If set will try and fix some integer variables"
     ); 
  parameters[numberParameters-1].setIntValue(-1);
#endif 
  assert(numberParameters<CBCMAXPARAMETERS);
}
// Given a parameter type - returns its number in list
int whichParam (CbcOrClpParameterType name, 
		int numberParameters, CbcOrClpParam *const parameters)
{
  int i;
  for (i=0;i<numberParameters;i++) {
    if (parameters[i].type()==name)
      break;
  }
  assert (i<numberParameters);
  return i;
}
#ifdef COIN_HAS_CLP
// Dump a solution to file
void saveSolution(const ClpSimplex * lpSolver,std::string fileName)
{
  FILE * fp=fopen(fileName.c_str(),"wb");
  if (fp) {
    int numberRows=lpSolver->numberRows();
    int numberColumns=lpSolver->numberColumns();
    double objectiveValue = lpSolver->objectiveValue();
    fwrite(&numberRows,sizeof(int),1,fp);
    fwrite(&numberColumns,sizeof(int),1,fp);
    fwrite(&objectiveValue,sizeof(double),1,fp);
    double * dualRowSolution = lpSolver->dualRowSolution();
    double * primalRowSolution = lpSolver->primalRowSolution();
    fwrite(primalRowSolution,sizeof(double),numberRows,fp);
    fwrite(dualRowSolution,sizeof(double),numberRows,fp);
    double * dualColumnSolution = lpSolver->dualColumnSolution();
    double * primalColumnSolution = lpSolver->primalColumnSolution();
    fwrite(primalColumnSolution,sizeof(double),numberColumns,fp);
    fwrite(dualColumnSolution,sizeof(double),numberColumns,fp);
    fclose(fp);
  } else {
    std::cout<<"Unable to open file "<<fileName<<std::endl;
  }
}
/* Restore a solution from file.
   mode 0 normal, 1 swap rows and columns and primal and dual
   if 2 set then also change signs
*/
void restoreSolution(ClpSimplex * lpSolver,std::string fileName,int mode)
{
  FILE * fp=fopen(fileName.c_str(),"rb");
  if (fp) {
    int numberRows=lpSolver->numberRows();
    int numberColumns=lpSolver->numberColumns();
    int numberRowsFile;
    int numberColumnsFile;
    double objectiveValue;
    fread(&numberRowsFile,sizeof(int),1,fp);
    fread(&numberColumnsFile,sizeof(int),1,fp);
    fread(&objectiveValue,sizeof(double),1,fp);
    double * dualRowSolution = lpSolver->dualRowSolution();
    double * primalRowSolution = lpSolver->primalRowSolution();
    double * dualColumnSolution = lpSolver->dualColumnSolution();
    double * primalColumnSolution = lpSolver->primalColumnSolution();
    if (mode) {
      // swap
      int k=numberRows;
      numberRows=numberColumns;
      numberColumns=k;
      double * temp;
      temp = dualRowSolution;
      dualRowSolution = primalColumnSolution;
      primalColumnSolution=temp;
      temp = dualColumnSolution;
      dualColumnSolution = primalRowSolution;
      primalRowSolution=temp;
    }
    if (numberRows!=numberRowsFile||numberColumns!=numberColumnsFile) {
      std::cout<<"Mismatch on rows and/or columns"<<std::endl;
    } else {
      lpSolver->setObjectiveValue(objectiveValue);
      fread(primalRowSolution,sizeof(double),numberRows,fp);
      fread(dualRowSolution,sizeof(double),numberRows,fp);
      fread(primalColumnSolution,sizeof(double),numberColumns,fp);
      fread(dualColumnSolution,sizeof(double),numberColumns,fp);
      if (mode==3) {
	int i;
	for (i=0;i<numberRows;i++) {
	  primalRowSolution[i] = -primalRowSolution[i];
	  dualRowSolution[i] = -dualRowSolution[i];
	}
	for (i=0;i<numberColumns;i++) {
	  primalColumnSolution[i] = -primalColumnSolution[i];
	  dualColumnSolution[i] = -dualColumnSolution[i];
	}
      }
    }
    fclose(fp);
  } else {
    std::cout<<"Unable to open file "<<fileName<<std::endl;
  }
}
#endif
