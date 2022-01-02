// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpConfig.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cassert>

#include "CoinParam.hpp"
#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpParam.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"

#ifdef CLP_CILK
#ifndef CLP_THREAD
#define CLP_THREAD
#endif
#endif

#if defined(CLP_HAS_WSMP) && !defined(USE_EKKWSSMP)
#ifndef CLP_THREAD
#define CLP_THREAD
#endif
#endif

#ifdef CLP_HAS_ABC
#include "AbcCommon.hpp"
#endif

#if COIN_INT_MAX == 0
#undef COIN_INT_MAX
#define COIN_INT_MAX 2147483647
#endif
#if FLUSH_PRINT_BUFFER > 2
int coinFlushBufferFlag = 0;
#endif

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpParam::ClpParam()
  : type_(CLP_PARAM_NOTUSED_INVALID)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , lengthName_(0)
  , lengthMatch_(0)
  , definedKeyWords_()
  , name_()
  , shortHelp_()
  , longHelp_()
  , action_(CLP_PARAM_NOTUSED_INVALID)
  , currentKeyWord_(-1)
  , display_(0)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
}
// Other constructors
ClpParam::ClpParam(std::string name, std::string help,
  double lower, double upper, ClpParameterType type,
  int display)
  : type_(type)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
ClpParam::ClpParam(std::string name, std::string help,
  int lower, int upper, ClpParameterType type,
  int display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(7)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
ClpParam::ClpParam(std::string name, std::string help,
  std::string firstValue,
  ClpParameterType type, int whereUsed,
  int display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(0)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , whereUsed_(whereUsed)
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(firstValue);
}
// Action
ClpParam::ClpParam(std::string name, std::string help,
  ClpParameterType type, int whereUsed,
  int display)
  : type_(type)
  , lowerDoubleValue_(0.0)
  , upperDoubleValue_(0.0)
  , lowerIntValue_(0)
  , upperIntValue_(0)
  , definedKeyWords_()
  , name_(name)
  , shortHelp_(help)
  , longHelp_()
  , action_(type)
  , currentKeyWord_(-1)
  , display_(display)
  , intValue_(-1)
  , doubleValue_(-1.0)
  , stringValue_("")
  , fakeKeyWord_(-1)
  , fakeValue_(0)
{
  whereUsed_ = whereUsed;
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpParam::ClpParam(const ClpParam &rhs)
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
  display_ = rhs.display_;
  intValue_ = rhs.intValue_;
  doubleValue_ = rhs.doubleValue_;
  stringValue_ = rhs.stringValue_;
  whereUsed_ = rhs.whereUsed_;
  fakeKeyWord_ = rhs.fakeKeyWord_;
  fakeValue_ = rhs.fakeValue_;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpParam::~ClpParam()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpParam &
ClpParam::operator=(const ClpParam &rhs)
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
    display_ = rhs.display_;
    intValue_ = rhs.intValue_;
    doubleValue_ = rhs.doubleValue_;
    stringValue_ = rhs.stringValue_;
    whereUsed_ = rhs.whereUsed_;
    fakeKeyWord_ = rhs.fakeKeyWord_;
    fakeValue_ = rhs.fakeValue_;
  }
  return *this;
}
void ClpParam::gutsOfConstructor()
{
  std::string::size_type shriekPos = name_.find('!');
  lengthName_ = static_cast< unsigned int >(name_.length());
  if (shriekPos == std::string::npos) {
    //does not contain '!'
    lengthMatch_ = lengthName_;
  } else {
    lengthMatch_ = static_cast< unsigned int >(shriekPos);
    name_ = name_.substr(0, shriekPos) + name_.substr(shriekPos + 1);
    lengthName_--;
  }
}

// Insert string (only valid for keywords)
void ClpParam::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

double
ClpParam::doubleParameter(ClpSimplex *model) const
{
  double value;
  switch (type_) {
  case CLP_PARAM_DBL_DUALTOLERANCE:
    value = model->dualTolerance();
    break;
  case CLP_PARAM_DBL_PRIMALTOLERANCE:
    value = model->primalTolerance();
    break;
  case CLP_PARAM_DBL_ZEROTOLERANCE:
    value = model->getSmallElementValue();
    break;
  case CLP_PARAM_DBL_DUALBOUND:
    value = model->dualBound();
    break;
  case CLP_PARAM_DBL_PRIMALWEIGHT:
    value = model->infeasibilityCost();
    break;
  case CLP_PARAM_DBL_TIMELIMIT:
    value = model->maximumSeconds();
    break;
  case CLP_PARAM_DBL_OBJSCALE:
    value = model->objectiveScale();
    break;
  case CLP_PARAM_DBL_RHSSCALE:
    value = model->rhsScale();
    break;
  case CLP_PARAM_DBL_PRESOLVETOLERANCE:
    value = model->presolveTolerance();
    break;
  default:
    value = doubleValue_;
    break;
  }
  return value;
}

int ClpParam::setDoubleParameter(ClpSimplex *model, double value,
                                 bool doPrinting)
{
  int returnCode;
  std::string message =
     setDoubleParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}

std::string
ClpParam::setDoubleParameterWithMessage(ClpSimplex *model, double value,
                                        int &returnCode)
{
  double oldValue = doubleValue_;
  std::ostringstream buffer;
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerDoubleValue_;
     buffer << " to " << upperDoubleValue_ << std::endl;
     returnCode = 1;
  } else {
    buffer << name_ << " was changed from ";
    buffer << oldValue << " to " << value << std::endl;
    returnCode = 0;
    doubleValue_ = value;
    switch (type_) {
    case CLP_PARAM_DBL_DUALTOLERANCE:
      model->setDualTolerance(value);
      break;
    case CLP_PARAM_DBL_PRIMALTOLERANCE:
      model->setPrimalTolerance(value);
      break;
    case CLP_PARAM_DBL_ZEROTOLERANCE:
      model->setSmallElementValue(value);
      break;
    case CLP_PARAM_DBL_DUALBOUND:
      model->setDualBound(value);
      break;
    case CLP_PARAM_DBL_PRIMALWEIGHT:
      model->setInfeasibilityCost(value);
      break;
    case CLP_PARAM_DBL_TIMELIMIT:
      model->setMaximumSeconds(value);
      break;
    case CLP_PARAM_DBL_OBJSCALE:
      model->setObjectiveScale(value);
      break;
    case CLP_PARAM_DBL_RHSSCALE:
      model->setRhsScale(value);
      break;
    case CLP_PARAM_DBL_PRESOLVETOLERANCE:
      model->setDblParam(ClpPresolveTolerance, value);
      break;
    case CLP_PARAM_DBL_PROGRESS:
      model->setMinIntervalProgressUpdate(value);
      break;
    default:
      break;
    }
    returnCode = 0;
  }
  return buffer.str();
}

void ClpParam::setDoubleValue(double value)
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
  } else {
    doubleValue_ = value;
  }
}

int ClpParam::checkDoubleParameter(double value) const
{
  if (value < lowerDoubleValue_ || value > upperDoubleValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerDoubleValue_ << " to " << upperDoubleValue_ << std::endl;
    return 1;
  } else {
    return 0;
  }
}

int ClpParam::intParameter(ClpSimplex *model) const
{
  int value;
  switch (type_) {
  case CLP_PARAM_INT_LOGLEVEL:
    value = model->logLevel();
    break;
  case CLP_PARAM_INT_MAXFACTOR:
    value = model->factorization()->maximumPivots();
    break;
    break;
  case CLP_PARAM_INT_PERTVALUE:
    value = model->perturbation();
    break;
  case CLP_PARAM_INT_MAXITERATION:
    value = model->maximumIterations();
    break;
  case CLP_PARAM_INT_SPECIALOPTIONS:
    value = model->specialOptions();
    break;
  case CLP_PARAM_INT_RANDOMSEED:
    value = model->randomNumberGenerator()->getSeed();
    break;
  case CLP_PARAM_INT_MORESPECIALOPTIONS:
    value = model->moreSpecialOptions();
    break;
  case CLP_PARAM_INT_VECTOR_MODE:
    value = model->vectorMode();
    break;
#ifdef CLP_THREAD
  case CLP_PARAM_INT_THREADS:
    value = model->numberThreads();
    break;
#endif
  default:
    value = intValue_;
    break;
  }
  return value;
}

int ClpParam::setIntParameter(ClpSimplex *model, int value,
                              bool doPrinting)
{
  int returnCode;
  std::string message = setIntParameterWithMessage(model, value, returnCode);
  if (doPrinting){
    std::cout << message << std::endl;
  }
  return returnCode;
}

std::string
ClpParam::setIntParameterWithMessage(ClpSimplex *model, int value,
                                     int &returnCode)
{
  int oldValue = intValue_;
  std::ostringstream buffer;
  if (value < lowerIntValue_ || value > upperIntValue_) {
     buffer << value << " was provided for " << name_;
     buffer << " - valid range is " << lowerIntValue_ << " to ";
     buffer << upperIntValue_ << std::endl;
     returnCode = 1;
  } else {
    intValue_ = value;
    if (type_ != CLP_PARAM_INT_RANDOMSEED){
       buffer << name_ << " was changed from " << oldValue;
       buffer << " to " << value << std::endl;
    }
    switch (type_) {
    case CLP_PARAM_INT_LOGLEVEL:
      model->setLogLevel(value);
      if (value > 2)
        model->factorization()->messageLevel(8);
      else
        model->factorization()->messageLevel(0);
      break;
    case CLP_PARAM_INT_MAXFACTOR:
      model->factorization()->maximumPivots(value);
      break;
    case CLP_PARAM_INT_PERTVALUE:
      model->setPerturbation(value);
      break;
    case CLP_PARAM_INT_MAXITERATION:
      model->setMaximumIterations(value);
      break;
    case CLP_PARAM_INT_SPECIALOPTIONS:
      model->setSpecialOptions(value);
      break;
    case CLP_PARAM_INT_RANDOMSEED: {
      if (value == 0) {
        double time = fabs(CoinGetTimeOfDay());
        while (time >= COIN_INT_MAX)
          time *= 0.5;
        value = static_cast< int >(time);
        buffer << "using time of day, " << name_ << " was changed from ";
        buffer << oldValue << " to " << value << std::endl;
      }
      model->setRandomSeed(value);
    } break;
    case CLP_PARAM_INT_MORESPECIALOPTIONS:
      model->setMoreSpecialOptions(value);
      break;
    case CLP_PARAM_INT_VECTOR_MODE:
      model->setVectorMode(value);
      break;
#ifdef CLP_THREAD
    //TODO: deal with this
    case CLP_PARAM_INT_THREADS:
      model->setNumberThreads(value);
      break;
#endif
    default:
      break;
    }
    returnCode = 0;
  }
  return buffer.str();
}

void ClpParam::setIntValue(int value)
{
  if (value < lowerIntValue_ || value > upperIntValue_) {
    std::cout << value << " was provided for " << name_ << " - valid range is " << lowerIntValue_ << " to " << upperIntValue_ << std::endl;
  } else {
    intValue_ = value;
  }
}

void ClpParam::setStringValue(std::string value)
{
  stringValue_ = value;
}

// Returns name which could match
std::string
ClpParam::matchName() const
{
  if (lengthMatch_ == lengthName_)
    return name_;
  else
    return name_.substr(0, lengthMatch_) + "(" + name_.substr(lengthMatch_) + ")";
}
// Returns length of name for printing
int ClpParam::lengthMatchName() const
{
  if (lengthName_ == lengthMatch_)
    return lengthName_;
  else
    return lengthName_ + 2;
}

int ClpParam::matches(std::string input) const
{
  // look up strings to do more elegantly
  if (input.length() > lengthName_) {
    return 0;
  } else {
    unsigned int i;
    for (i = 0; i < input.length(); i++) {
      if (tolower(name_[i]) != tolower(input[i]))
        break;
    }
    if (i < input.length()) {
      return 0;
    } else if (i >= lengthMatch_) {
      return 1;
    } else {
      // matched but too short
      return 2;
    }
  }
}
// Returns parameter option which matches (-1 if none)
int ClpParam::parameterOption(std::string check) const
{
  int numberItems = static_cast< int >(definedKeyWords_.size());
  if (!numberItems) {
    return -1;
  } else {
    int whichItem = 0;
    unsigned int it;
    for (it = 0; it < definedKeyWords_.size(); it++) {
      std::string thisOne = definedKeyWords_[it];
      std::string::size_type shriekPos = thisOne.find('!');
      size_t length1 = thisOne.length();
      size_t length2 = length1;
      if (shriekPos != std::string::npos) {
        //contains '!'
        length2 = shriekPos;
        thisOne = thisOne.substr(0, shriekPos) + thisOne.substr(shriekPos + 1);
        length1 = thisOne.length();
      }
      if (check.length() <= length1 && length2 <= check.length()) {
        unsigned int i;
        for (i = 0; i < check.length(); i++) {
          if (tolower(thisOne[i]) != tolower(check[i]))
            break;
        }
        if (i < check.length()) {
          whichItem++;
        } else if (i >= length2) {
          break;
        }
      } else {
        whichItem++;
      }
    }
    if (whichItem < numberItems) {
      return whichItem;
    } else {
      if (fakeKeyWord_ <= 0)
        return -1;
      // allow plus or minus
      int n;
      if (check.substr(0, 4) == "plus" || check.substr(0, 4) == "PLUS") {
        n = 4;
      } else if (check.substr(0, 5) == "minus" || check.substr(0, 5) == "MINUS") {
        n = 5;
      } else {
        return -1;
      }
      int value = 0;
      std::string field = check.substr(n);
      if (field != "EOL") {
        const char *start = field.c_str();
        char *endPointer = NULL;
        // check valid
        value = static_cast< int >(strtol(start, &endPointer, 10));
        if (*endPointer != '\0') {
          return -1;
        }
        if (n == 4)
          return value + 1000;
        else
          return -value - 1000;
      } else {
        return -1;
      }
    }
  }
}
// Prints parameter options
void ClpParam::printOptions() const
{
  std::cout << "<Possible options for " << name_ << " are:";
  unsigned int it;
  for (it = 0; it < definedKeyWords_.size(); it++) {
    std::string thisOne = definedKeyWords_[it];
    std::string::size_type shriekPos = thisOne.find('!');
    if (shriekPos != std::string::npos) {
      //contains '!'
      thisOne = thisOne.substr(0, shriekPos) + "(" + thisOne.substr(shriekPos + 1) + ")";
    }
    std::cout << " " << thisOne;
  }
  assert(currentKeyWord_ >= 0 && currentKeyWord_ < static_cast< int >(definedKeyWords_.size()));
  std::string current = definedKeyWords_[currentKeyWord_];
  std::string::size_type shriekPos = current.find('!');
  if (shriekPos != std::string::npos) {
    //contains '!'
    current = current.substr(0, shriekPos) + "(" + current.substr(shriekPos + 1) + ")";
  }
  std::cout << ";\n\tcurrent  " << current << ">" << std::endl;
}
void ClpParam::setCurrentOption(const std::string value)
{
  int action = parameterOption(value);
  if (action >= 0)
    currentKeyWord_ = action;
#if FLUSH_PRINT_BUFFER > 2
  if (name_ == "bufferedMode")
    coinFlushBufferFlag = action;
#endif
}
// Sets current parameter option
void ClpParam::setCurrentOption(int value, bool printIt)
{
  if (printIt && value != currentKeyWord_)
    std::cout << "Option for " << name_ << " changed from "
              << definedKeyWords_[currentKeyWord_] << " to "
              << definedKeyWords_[value] << std::endl;

#if FLUSH_PRINT_BUFFER > 2
  if (name_ == "bufferedMode")
    coinFlushBufferFlag = value;
#endif
  currentKeyWord_ = value;
}
// Sets current parameter option and returns printable string
std::string
ClpParam::setCurrentOptionWithMessage(int value)
{
  std::ostringstream buffer;
  if (value != currentKeyWord_) {
    buffer << "Option for " << name_ << " changed from ";
    char current[100];
    char newString[100];
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 ||
                                 currentKeyWord_ < fakeKeyWord_)){
       buffer << definedKeyWords_[currentKeyWord_];
    }else if (currentKeyWord_ < 0){
       buffer << "minus" << -currentKeyWord_ - 1000;
    }else{
       buffer << "plus" << currentKeyWord_ - 1000;
    }
    buffer << " to ";
    if (value >= 0 && (fakeKeyWord_ <= 0 || value < fakeKeyWord_)){
       buffer << definedKeyWords_[value];
    }else if (value < 0){
       buffer << "minus" << -value - 1000;
    }else{
       buffer << "plus" << value - 1000;
    }
    
#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode"){
      coinFlushBufferFlag = value;
    }
#endif
    currentKeyWord_ = value;
  }
  return buffer.str();
}

// Sets current parameter option using string with message
std::string
ClpParam::setCurrentOptionWithMessage(const std::string value)
{
  std::ostringstream buffer;
  int action = parameterOption(value);
  char current[100];
  if (action >= 0) {
#if FLUSH_PRINT_BUFFER > 2
    if (name_ == "bufferedMode")
      coinFlushBufferFlag = action;
#endif
    buffer << "Option for " << name_ << " changed from ";
    if (currentKeyWord_ >= 0 && (fakeKeyWord_ <= 0 ||
                                 currentKeyWord_ < fakeKeyWord_)){
       buffer << definedKeyWords_[currentKeyWord_];
    }else if (currentKeyWord_ < 0){
       buffer << "minus" << -currentKeyWord_ - 1000;
    }else{
       buffer << "plus" << currentKeyWord_ - 1000;
    }
    buffer << " to " << value;
    currentKeyWord_ = action;
  } else {
     buffer << "Option for " << name_ << " given illegal value " << value;
  }
  return buffer.str();
}

/* Returns current parameter option position
   but if fake keyword returns fakeValue_
*/
int ClpParam::currentOptionAsInteger() const
{
  int fakeInteger;
  return currentOptionAsInteger(fakeInteger);
}
/* Returns current parameter option position
   but if fake keyword returns fakeValue_ and sets
   fakeInteger to value
*/
int ClpParam::currentOptionAsInteger(int &fakeInteger) const
{
  fakeInteger = -COIN_INT_MAX;
  if (fakeKeyWord_ < 0) {
    return currentKeyWord_;
  } else if (currentKeyWord_ >= 0 && currentKeyWord_ < fakeKeyWord_) {
    return currentKeyWord_;
  } else {
    // fake
    if (currentKeyWord_ < 0)
      fakeInteger = currentKeyWord_ + 1000;
    else
      fakeInteger = currentKeyWord_ - 1000;
    return fakeValue_;
  }
}

// TODO: Fix this
// Print Long help
void ClpParam::printLongHelp() const
{
  if (type_ >= 1 && type_ < 600) {
    CoinPrintString(longHelp_);
    if (type_ < CLP_PARAM_INT_LOGLEVEL) {
      printf("<Range of values is %g to %g;\n\tcurrent %g>\n", lowerDoubleValue_, upperDoubleValue_, doubleValue_);
      assert(upperDoubleValue_ > lowerDoubleValue_);
    } else if (type_ < CLP_PARAM_STR_DIRECTION) {
      printf("<Range of values is %d to %d;\n\tcurrent %d>\n", lowerIntValue_, upperIntValue_, intValue_);
      assert(upperIntValue_ > lowerIntValue_);
    } else if (type_ < CLP_PARAM_ACTION_DIRECTORY) {
      printOptions();
    }
  }
}

// Print action and string
void ClpParam::printString() const
{
  if (name_ == "directory")
    std::cout << "Current working directory is " << stringValue_ << std::endl;
  else if (name_.substr(0, 6) == "printM")
    std::cout << "Current value of printMask is " << stringValue_ << std::endl;
  else
    std::cout << "Current default (if $ as parameter) for " << name_
              << " is " << stringValue_ << std::endl;
}

// Sets value of fake keyword to current size of keywords
void ClpParam::setFakeKeyWord(int fakeValue)
{
  fakeKeyWord_ = static_cast< int >(definedKeyWords_.size());
  assert(fakeKeyWord_ > 0);
  fakeValue_ = fakeValue;
  assert(fakeValue_ >= 0);
}

//###########################################################################
//###########################################################################

/*
  Subroutine to establish the cbc parameter array. See the description of
  class ClpParam for details. Pulled from C..Main() for clarity.
*/
void establishClpParams(std::vector< ClpParam > &parameters)
{
  parameters.clear();
  parameters.push_back(ClpParam("?", "For help", CLP_PARAM_GENERALQUERY, 7, 0));
  parameters.push_back(ClpParam("???", "For help", CLP_PARAM_FULLGENERALQUERY, 7, 0));
  parameters.push_back(ClpParam("-", "From stdin", CLP_PARAM_ACTION_STDIN, 3, 0));

// some help strings that repeat for many options
#define CUTS_LONGHELP \
  "Value 'on' enables the cut generator and CBC will try it in the branch and cut tree (see cutDepth on how to fine tune the behavior). " \
  "Value 'root' lets CBC run the cut generator generate only at the root node. " \
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as if it is doing some good and moves the objective value. " \
  "Value 'forceon' turns on the cut generator and forces CBC to use it at every node."
#define HEURISTICS_LONGHELP \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. after preprocessing. " \
  "Value 'before' means use the heuristic only if option doHeuristics is used. " \
  "Value 'both' means to use the heuristic if option doHeuristics is used and during solve."

#if CLP_HAS_ABC
  ClpParam paramAboca("abc", "Whether to visit Aboca", "off", CLP_PARAM_STR_ABCWANTED, 7, 0);
  paramAboca.append("one");
  paramAboca.append("two");
  paramAboca.append("three");
  paramAboca.append("four");
  paramAboca.append("five");
  paramAboca.append("six");
  paramAboca.append("seven");
  paramAboca.append("eight");
  paramAboca.append("on");
  paramAboca.append("decide");
  paramAboca.setFakeKeyWord(10);
  paramAboca.setLonghelp(
    "Decide whether to use A Basic Optimization Code (Accelerated?) \
and whether to try going parallel!");
  parameters.push_back(paramAboca);
#endif
  {
    ClpParam p("allC!ommands", "Whether to print less used commands",
      "no", CLP_PARAM_STR_ALLCOMMANDS);

    p.append("more");
    p.append("all");
    p.setLonghelp(
      "For the sake of your sanity, only the more useful and simple commands \
are printed out on ?.");
    parameters.push_back(p);
  }
  {
    ClpParam p("allS!lack", "Set basis back to all slack and reset solution",
      CLP_PARAM_ACTION_ALLSLACK, 3);
    p.setLonghelp(
      "Mainly useful for tuning purposes.  Normally the first dual or primal will be using an all slack \
basis anyway.");
    parameters.push_back(p);
  }
  {
    ClpParam p("auto!Scale", "Whether to scale objective, rhs and bounds of problem if they look odd",
      "off", CLP_PARAM_STR_AUTOSCALE, 7, 0);
    p.append("on");
    p.setLonghelp(
      "If you think you may get odd objective values or large equality rows etc then\
 it may be worth setting this true.  It is still experimental and you may prefer\
 to use objective!Scale and rhs!Scale.");
    parameters.push_back(p);
  }
  {
    ClpParam p("barr!ier", "Solve using primal dual predictor corrector algorithm",
      CLP_PARAM_ACTION_BARRIER);
    p.setLonghelp(
      "This command solves the current model using the  primal dual predictor \
corrector algorithm.  You may want to link in an alternative \
ordering and factorization.  It will also solve models \
with quadratic objectives.");
    parameters.push_back(p);
  }
  {
    ClpParam p("basisI!n", "Import basis from bas file",
      CLP_PARAM_ACTION_BASISIN, 3);
    p.setLonghelp(
      "This will read an MPS format basis file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libz then it can read compressed\
 files 'xxxxxxxx.gz' or xxxxxxxx.bz2.");
    parameters.push_back(p);
  }
  {
    ClpParam p("basisO!ut", "Export basis as bas file",
      CLP_PARAM_ACTION_BASISOUT);

    p.setLonghelp(
      "This will write an MPS format basis file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.bas'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("biasLU", "Whether factorization biased towards U",
      "UU", CLP_PARAM_STR_BIASLU, 2, 0);

    p.append("UX");
    p.append("LX");
    p.append("LL");
    p.setCurrentOption("LX");
    parameters.push_back(p);
  }
  {
    ClpParam p("branch!AndCut", "Do Branch and Cut",
      CLP_PARAM_ACTION_BAB);
    p.setLonghelp(
      "This does branch and cut.  There are many parameters which can affect the performance.  \
First just try with default settings and look carefully at the log file.  Did cuts help?  Did they take too long?  \
Look at output to see which cuts were effective and then do some tuning.  You will see that the \
options for cuts are off, on, root and ifmove, forceon.  Off is \
obvious. " CUTS_LONGHELP " For probing, forceonbut just does fixing probing in tree - not strengthening etc.  \
If pre-processing reduced the size of the \
problem or strengthened many coefficients then it is probably wise to leave it on.  Switch off heuristics \
which did not provide solutions.  The other major area to look at is the search.  Hopefully good solutions \
were obtained fairly early in the search so the important point is to select the best variable to branch on.  \
See whether strong branching did a good job - or did it just take a lot of iterations.  Adjust the strongBranching \
and trustPseudoCosts parameters.  If cuts did a good job, then you may wish to \
have more rounds of cuts - see passC!uts and passT!ree.");
    parameters.push_back(p);
  }
  {
    ClpParam p("bscale", "Whether to scale in barrier (and ordering speed)",
      "off", CLP_PARAM_STR_BARRIERSCALE, 7, 0);
    p.append("on");
    p.append("off1");
    p.append("on1");
    p.append("off2");
    p.append("on2");
    parameters.push_back(p);
  }
#if FLUSH_PRINT_BUFFER > 2
  {
    ClpParam p("buff!eredMode", "Whether to flush print buffer",
      "on", CLP_PARAM_STR_BUFFER_MODE);

    p.append("off");
    p.setLonghelp(
      "Default is on, off switches on unbuffered output");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#endif
  {
    ClpParam p("chol!esky", "Which cholesky algorithm",
      "native", CLP_PARAM_STR_CHOLESKY, 7);

    p.append("dense");
    //#ifdef FOREIGN_BARRIER
#ifdef CLP_HAS_WSMP
    p.append("fudge!Long");
    p.append("wssmp");
#else
    p.append("fudge!Long_dummy");
    p.append("wssmp_dummy");
#endif
#if defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD) || defined(CLP_HAS_GLPK)
    p.append("Uni!versityOfFlorida");
#else
    p.append("Uni!versityOfFlorida_dummy");
#endif
#ifdef TAUCS_BARRIER
    p.append("Taucs");
#else
    p.append("Taucs_dummy");
#endif
#ifdef CLP_HAS_MUMPS
    p.append("Mumps");
#else
    p.append("Mumps_dummy");
#endif
#ifdef PARDISO_BARRIER
    p.append("Pardiso");
#else
    p.append("Pardiso_dummy");
#endif
    p.setLonghelp(
      "For a barrier code to be effective it needs a good Cholesky ordering and factorization.  \
The native ordering and factorization is not state of the art, although acceptable.  \
You may want to link in one from another source.  See Makefile.locations for some \
possibilities.");

    parameters.push_back(p);
  }
  {
    ClpParam p("cpp!Generate", "Generates C++ code",
      -1, 50000, CLP_PARAM_INT_CPP, 1);
    p.setLonghelp(
      "Once you like what the stand-alone solver does then this allows \
you to generate user_driver.cpp which approximates the code.  \
0 gives simplest driver, 1 generates saves and restores, 2 \
generates saves and restores even for variables at default value. \
4 bit in cbc generates size dependent code rather than computed values.  \
This is now deprecated as you can call stand-alone solver - see \
Cbc/examples/driver4.cpp.");
    parameters.push_back(p);
  }
  {
    ClpParam p("crash", "Whether to create basis for problem",
      "off", CLP_PARAM_STR_CRASH);

    p.append("on");
    p.append("so!low_halim");
    p.append("lots");
    p.append("free");
    p.append("zero");
    p.append("single!ton");
#ifdef CLP_INHERIT_MODE
    p.append("dual");
    p.append("dw");
    p.append("idiot");
#else
    p.append("idiot1");
    p.append("idiot2");
    p.append("idiot3");
    p.append("idiot4");
    p.append("idiot5");
    p.append("idiot6");
    p.append("idiot7");
#endif
    p.setLonghelp(
      "If crash is set to 'on' and there is an all slack basis then Clp will flip or put structural\
     variables into the basis with the aim of getting dual feasible.  On average, dual simplex seems to perform\
     better without it and there are alternative types of 'crash' for primal simplex, e.g. 'idiot' or 'sprint'. \
    A variant due to Solow and Halim which is as 'on' but just flips is also available.");

    parameters.push_back(p);
  }
  {
    ClpParam p("cross!over", "Whether to get a basic solution with the simplex algorithm after the barrier algorithm finished",
      "on", CLP_PARAM_STR_CROSSOVER);
    p.append("off");
    p.append("maybe");
    p.append("presolve");
    p.setLonghelp(
      "Interior point algorithms do not obtain a basic solution.\
 This option will crossover to a basic solution suitable for ranging or branch and cut.  With the current state \
of the solver for quadratic programs it may be a good idea to switch off crossover for this case (and maybe \
presolve as well) - the option 'maybe' does this.");
    parameters.push_back(p);
  }
  {
    ClpParam p("decomp!ose", "Whether to try decomposition",
      -COIN_INT_MAX, COIN_INT_MAX, CLP_PARAM_INT_DECOMPOSE_BLOCKS, 1);
    p.setLonghelp(
      "0 - off, 1 choose blocks >1 use as blocks \
Dantzig Wolfe if primal, Benders if dual \
- uses sprint pass for number of passes");
    p.setIntValue(0);
    parameters.push_back(p);
  }
#if CLP_MULTIPLE_FACTORIZATIONS > 0
  {
    ClpParam p("dense!Threshold", "Threshold for using dense factorization",
      -1, 10000, CLP_PARAM_INT_DENSE, 1);
    p.setLonghelp(
      "If processed problem <= this use dense factorization");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
#endif
  {
    ClpParam p("direction", "Minimize or Maximize",
      "min!imize", CLP_PARAM_STR_DIRECTION);
    p.append("max!imize");
    p.append("zero");
    p.setLonghelp(
      "The default is minimize - use 'direction maximize' for maximization.\n\
You can also use the parameters 'maximize' or 'minimize'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("directory", "Set Default directory for import etc.",
      CLP_PARAM_ACTION_DIRECTORY);
    p.setLonghelp(
      "This sets the directory which import, export, saveModel, restoreModel etc will use.\
  It is initialized to './'");
    parameters.push_back(p);
  }
  {
    ClpParam p("dirSample", "Set directory where the COIN-OR sample problems are.",
      CLP_PARAM_ACTION_DIRSAMPLE, 7, 1);

    p.setLonghelp(
      "This sets the directory where the COIN-OR sample problems reside. It is\
 used only when -unitTest is passed to clp. clp will pick up the test problems\
 from this directory.\
 It is initialized to '../../Data/Sample'");
    parameters.push_back(p);
  }
  {
    ClpParam p("dirNetlib", "Set directory where the netlib problems are.",
      CLP_PARAM_ACTION_DIRNETLIB, 7, 1);

    p.setLonghelp(
      "This sets the directory where the netlib problems reside. One can get\
 the netlib problems from COIN-OR or from the main netlib site. This\
 parameter is used only when -netlib is passed to clp. clp will pick up the\
 netlib problems from this directory. If clp is built without zlib support\
 then the problems must be uncompressed.\
 It is initialized to '../../Data/Netlib'");
    parameters.push_back(p);
  }
  {
    ClpParam p("dirMiplib", "Set directory where the miplib 2003 problems are.",
      CLP_PARAM_ACTION_DIRMIPLIB, 7, 1);

    p.setLonghelp(
      "This sets the directory where the miplib 2003 problems reside. One can\
 get the miplib problems from COIN-OR or from the main miplib site. This\
 parameter is used only when -miplib is passed to cbc. cbc will pick up the\
 miplib problems from this directory. If cbc is built without zlib support\
 then the problems must be uncompressed.\
 It is initialized to '../../Data/miplib3'");
    parameters.push_back(p);
  }
  {
    ClpParam p("dualB!ound", "Initially algorithm acts as if no \
gap between bounds exceeds this value",
      1.0e-20, 1.0e20, CLP_PARAM_DBL_DUALBOUND);
    p.setLonghelp(
      "The dual algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  If a problem has both upper and\
 lower bounds then it is trivial to get dual feasible by setting non basic variables\
 to correct bound.  If the gap between the upper and lower bounds of a variable is more\
 than the value of dualBound Clp introduces fake bounds so that it can make the problem\
 dual feasible.  This has the same effect as a composite objective function in the\
 primal algorithm.  Too high a value may mean more iterations, while too low a bound means\
 the code may go all the way and then have to increase the bounds.  OSL had a heuristic to\
 adjust bounds, maybe we need that here.");

    parameters.push_back(p);
  }
  {
    ClpParam p("dualize", "Solves dual reformulation",
      0, 4, CLP_PARAM_INT_DUALIZE, 1);
    p.setLonghelp(
      "Don't even think about it.");

    parameters.push_back(p);
  }
  {
    ClpParam p("dualP!ivot", "Dual pivot choice algorithm",
      "auto!matic", CLP_PARAM_STR_DUALPIVOT, 7, 1);
    p.append("dant!zig");
    p.append("partial");
    p.append("steep!est");
    p.append("PEsteep!est");
    p.append("PEdantzig");
    p.setLonghelp(
      "The Dantzig method is simple but its use is deprecated.  Steepest is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while automatic decides at each iteration based on information\
 about the factorization.\
 The PE variants add the Positive Edge criterion. \
 This selects incoming variables to try to avoid degenerate moves. See also option psi.");
    parameters.push_back(p);
  }
  {
    ClpParam p("dualS!implex", "Do dual simplex algorithm",
      CLP_PARAM_ACTION_DUALSIMPLEX);
    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the dual steepest edge algorithm.\
The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by dual pivot method, fake bound on variables and dual and primal tolerances.");
    parameters.push_back(p);
  }
  {
    ClpParam p("initialS!olve", "Solve to continuous",
      CLP_PARAM_ACTION_SOLVECONTINUOUS);

    p.setLonghelp(
      "This just solves the problem to continuous - without adding any cuts");
    parameters.push_back(p);
  }
  {
    ClpParam p("dualT!olerance", "For an optimal solution \
no dual infeasibility may exceed this value",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_DUALTOLERANCE);
    p.setLonghelp(
      "Normally the default tolerance is fine, but one may want to increase it a\
 bit if the dual simplex algorithm seems to be having a hard time.  One method which can be faster is \
to use a large tolerance e.g. 1.0e-4 and the dual simplex algorithm and then to clean up the problem using the primal simplex algorithm with the \
correct tolerance (remembering to switch off presolve for this final short clean up phase).");
    parameters.push_back(p);
  }
  {
    ClpParam p("either!Simplex", "Do dual or primal simplex algorithm",
      CLP_PARAM_ACTION_EITHERSIMPLEX);
    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the dual or primal algorithm,\
 based on a dubious analysis of model.");
    parameters.push_back(p);
  }
  {
    ClpParam p("end", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops execution ; end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
  {
    ClpParam p("environ!ment", "Read commands from environment",
      CLP_PARAM_ACTION_ENVIRONMENT, 7, 0);
    p.setLonghelp(
      "This starts reading from environment variable CBC_CLP_ENVIRONMENT.");
    parameters.push_back(p);
  }
  {
    ClpParam p("error!sAllowed", "Whether to allow import errors",
      "off", CLP_PARAM_STR_ERRORSALLOWED, 3);

    p.append("on");
    p.setLonghelp(
      "The default is not to use any model which had errors when reading the mps file.\
  Setting this to 'on' will allow all errors from which the code can recover\
 simply by ignoring the error.  There are some errors from which the code can not recover \
e.g. no ENDATA.  This has to be set before import i.e. -errorsAllowed on -import xxxxxx.mps.");
    parameters.push_back(p);
  }
  {
    ClpParam p("exit", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
  {
    ClpParam p("export", "Export model as mps file",
      CLP_PARAM_ACTION_EXPORT);

    p.setLonghelp(
      "This will write an MPS format file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.mps'.  \
It can be useful to get rid of the original names and go over to using Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off before importing mps file.");
    parameters.push_back(p);
  }
  {
    ClpParam p("fact!orization", "Which factorization to use",
      "normal", CLP_PARAM_STR_FACTORIZATION);
    p.append("dense");
    p.append("simple");
    p.append("osl");
    p.setLonghelp(
#ifndef ABC_INHERIT
      "The default is to use the normal CoinFactorization, but \
other choices are a dense one, OSL's, or one designed for small problems."
#else
      "Normally the default is to use the normal CoinFactorization, but \
other choices are a dense one, OSL's, or one designed for small problems. \
However if at Aboca then the default is CoinAbcFactorization and other choices are \
a dense one, one designed for small problems or if enabled a long factorization."
#endif
    );
    parameters.push_back(p);
  }
  {
    ClpParam p("fakeB!ound", "All bounds <= this value - DEBUG",
      1.0, 1.0e15, CLP_PARAM_ACTION_FAKEBOUND, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("gamma!(Delta)", "Whether to regularize barrier",
      "off", CLP_PARAM_STR_GAMMA, 7, 1);
    p.append("on");
    p.append("gamma");
    p.append("delta");
    p.append("onstrong");
    p.append("gammastrong");
    p.append("deltastrong");
    parameters.push_back(p);
  }
  {
    ClpParam p("gsolu!tion", "Puts glpk solution to file",
      CLP_PARAM_ACTION_GMPL_SOLUTION);

    p.setLonghelp(
      "Will write a glpk solution file to the given file name.  It will use the default \
directory given by 'directory'.  A name of '$' will use the previous value for the \
name.  This is initialized to 'stdout' (this defaults to ordinary solution if stdout). \
If problem created from gmpl model - will do any reports.");
    parameters.push_back(p);
  }
  {
    ClpParam p("guess", "Guesses at good parameters", CLP_PARAM_ACTION_GUESS, 7);
    p.setLonghelp(
    "This looks at model statistics and does an initial solve \
setting some parameters which may help you to think of possibilities.");
    parameters.push_back(p);
  }
  {
    ClpParam p("idiot!Crash", "Whether to try idiot crash",
      -1, COIN_INT_MAX, CLP_PARAM_INT_IDIOT);

    p.setLonghelp(
      "This is a type of 'crash' which works well on some homogeneous problems.\
 It works best on problems with unit elements and rhs but will do something to any \
 model.  It should only be used before the primal simplex algorithm.  It can be set to -1 when the code \
 decides for itself whether to use it, 0 to switch off, or n > 0 to do n passes.");
    parameters.push_back(p);
  }
  {
    ClpParam p("import", "Import model from mps file",
      CLP_PARAM_ACTION_IMPORT, 3);
    p.setLonghelp(
      "This will read an MPS format file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  If you have libgz then it can read compressed\
 files 'xxxxxxxx.gz' or 'xxxxxxxx.bz2'.  \
If 'keepnames' is off, then names are dropped -> Rnnnnnnn and Cnnnnnnn.");
    parameters.push_back(p);
  }
  {
    ClpParam p("keepN!ames", "Whether to keep names from import",
      "on", CLP_PARAM_STR_KEEPNAMES);
    p.append("off");
    p.setLonghelp(
      "It saves space to get rid of names so if you need to you can set this to off.  \
This needs to be set before the import of model - so -keepnames off -import xxxxx.mps.");
    parameters.push_back(p);
  }
  {
    ClpParam p("KKT", "Whether to use KKT factorization in barrier",
      "off", CLP_PARAM_STR_KKT, 7, 1);
    p.append("on");
    parameters.push_back(p);
  }
  {
    ClpParam p("log!Level", "Level of detail in Coin branch and Cut output",
      -63, 63, CLP_PARAM_INT_LOGLEVEL);
    p.setIntValue(1);
    p.setLonghelp(
      "If 0 then there should be no output in normal circumstances.  1 is probably the best\
 value for most uses, while 2 and 3 give more information.");
    parameters.push_back(p);
  }
  {
    ClpParam p("max!imize", "Set optimization direction to maximize",
      CLP_PARAM_ACTION_MAXIMIZE, 7);
    p.setLonghelp(
      "The default is minimize - use 'maximize' for maximization.\n\
You can also use the parameters 'direction maximize'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("maxF!actor", "Maximum number of iterations between \
refactorizations",
      1, COIN_INT_MAX, CLP_PARAM_INT_MAXFACTOR);
    p.setLonghelp(
      "If this is left at its default value of 200 then CLP will guess a\
 value to use.  CLP may decide to re-factorize earlier for accuracy.");
    parameters.push_back(p);
  }
  {
    ClpParam p("maxIt!erations", "Maximum number of iterations before \
stopping",
      0, COIN_INT_MAX, CLP_PARAM_INT_MAXITERATION);
    p.setLonghelp(
      "This can be used for testing purposes.  The corresponding library call\n\
      \tsetMaximumIterations(value)\n can be useful.  If the code stops on\
 seconds or by an interrupt this will be treated as stopping on maximum iterations.  This is ignored in branchAndCut - use maxN!odes.");
    parameters.push_back(p);
  }
  {
    ClpParam p("min!imize", "Set optimization direction to minimize",
      CLP_PARAM_ACTION_MINIMIZE, 7);
    p.setLonghelp(
      "The default is minimize - use 'maximize' for maximization.\n\
This should only be necessary if you have previously set maximization \
You can also use the parameters 'direction minimize'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("mess!ages", "Controls if Clpnnnn is printed",
      "off", CLP_PARAM_STR_MESSAGES);

    p.append("on");
    p.setLonghelp("The default behavior is to put out messages such as:\n\
   Clp0005 2261  Objective 109.024 Primal infeas 944413 (758)\n\
but this program turns this off to make it look more friendly.  It can be useful\
 to turn them back on if you want to be able to 'grep' for particular messages or if\
 you intend to override the behavior of a particular message.  This only affects Clp not Cbc.");
    parameters.push_back(p);
  }
  {
    ClpParam p("moreS!pecialOptions", "Yet more dubious options for Simplex - see ClpSimplex.hpp",
      0, COIN_INT_MAX, CLP_PARAM_INT_MORESPECIALOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("netlib", "Solve entire netlib test set",
      CLP_PARAM_ACTION_NETLIB_EITHER, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using dual or primal.\
The user can set options before e.g. clp -presolve off -netlib");
    parameters.push_back(p);
  }
  {
    ClpParam p("netlibB!arrier", "Solve entire netlib test set with barrier",
      CLP_PARAM_ACTION_NETLIB_BARRIER, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using barrier.\
The user can set options before e.g. clp -kkt on -netlib");
    parameters.push_back(p);
  }
  {
    ClpParam p("netlibD!ual", "Solve entire netlib test set (dual)",
      CLP_PARAM_ACTION_NETLIB_DUAL, 3, 1);

    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using dual.\
The user can set options before e.g. clp -presolve off -netlib");
    parameters.push_back(p);
  }
  {
    ClpParam p("netlibP!rimal", "Solve entire netlib test set (primal)",
      CLP_PARAM_ACTION_NETLIB_PRIMAL, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using primal.\
The user can set options before e.g. clp -presolve off -netlibp");
    parameters.push_back(p);
  }
  {
    ClpParam p("netlibT!une", "Solve entire netlib test set with 'best' algorithm",
      CLP_PARAM_ACTION_NETLIB_TUNE, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp and then solves the netlib test set using whatever \
works best.  I know this is cheating but it also stresses the code better by doing a \
mixture of stuff.  The best algorithm was chosen on a Linux ThinkPad using native cholesky \
with University of Florida ordering.");
    parameters.push_back(p);
  }
  {
    ClpParam p("network", "Tries to make network matrix",
      CLP_PARAM_ACTION_NETWORK, 7, 0);
    p.setLonghelp(
      "Clp will go faster if the matrix can be converted to a network.  The matrix\
 operations may be a bit faster with more efficient storage, but the main advantage\
 comes from using a network factorization.  It will probably not be as fast as a \
specialized network code.");
    parameters.push_back(p);
  }
  {
    ClpParam p("objective!Scale", "Scale factor to apply to objective",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_OBJSCALE, 1);
    p.setLonghelp(
      "If the objective function has some very large values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  It is applied after scaling.  You are unlikely to need this.");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
  {
    ClpParam p("outDup!licates", "takes duplicate rows etc out of integer model",
      CLP_PARAM_ACTION_OUTDUPROWS, 7, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("output!Format", "Which output format to use",
      1, 6, CLP_PARAM_INT_OUTPUTFORMAT);
    p.setLonghelp(
      "Normally export will be done using normal representation for numbers and two values\
 per line.  You may want to do just one per line (for grep or suchlike) and you may wish\
 to save with absolute accuracy using a coded version of the IEEE value. A value of 2 is normal.\
 otherwise odd values gives one value per line, even two.  Values 1,2 give normal format, 3,4\
 gives greater precision, while 5,6 give IEEE values.  When used for exporting a basis 1 does not save \
values, 2 saves values, 3 with greater accuracy and 4 in IEEE.");
    parameters.push_back(p);
  }
  {
    ClpParam p("para!metrics", "Import data from file and do parametrics",
      CLP_PARAM_ACTION_PARAMETRICS, 3);
    p.setLonghelp(
      "This will read a file with parametric data from the given file name \
and then do parametrics.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to '', i.e. it must be set.  This can not read from compressed files. \
File is in modified csv format - a line ROWS will be followed by rows data \
while a line COLUMNS will be followed by column data.  The last line \
should be ENDATA. The ROWS line must exist and is in the format \
ROWS, inital theta, final theta, interval theta, n where n is 0 to get \
CLPI0062 message at interval or at each change of theta \
and 1 to get CLPI0063 message at each iteration.  If interval theta is 0.0 \
or >= final theta then no interval reporting.  n may be missed out when it is \
taken as 0.  If there is Row data then \
there is a headings line with allowed headings - name, number, \
lower(rhs change), upper(rhs change), rhs(change).  Either the lower and upper \
fields should be given or the rhs field. \
The optional COLUMNS line is followed by a headings line with allowed \
headings - name, number, objective(change), lower(change), upper(change). \
 Exactly one of name and number must be given for either section and \
missing ones have value 0.0.");
    parameters.push_back(p);
  }
  {
    ClpParam p("passP!resolve", "How many passes in presolve",
      -200, 100, CLP_PARAM_INT_PRESOLVEPASS, 1);
    p.setLonghelp(
      "Normally Presolve does 10 passes but you may want to do less to make it\
 more lightweight or do more if improvements are still being made.  As Presolve will return\
 if nothing is being taken out, you should not normally need to use this fine tuning.");
    parameters.push_back(p);
  }
  {
    ClpParam p("pertV!alue", "Method of perturbation",
      -5000, 102, CLP_PARAM_INT_PERTVALUE, 1);
    parameters.push_back(p);
  }
  {
    ClpParam p("perturb!ation", "Whether to perturb the problem",
      "on", CLP_PARAM_STR_PERTURBATION);
    p.append("off");
    p.setLonghelp(
      "Perturbation helps to stop cycling, but CLP uses other measures for this.\
  However, large problems and especially ones with unit elements and unit right hand sides or costs\
 benefit from perturbation.  Normally CLP tries to be intelligent, but one can switch this off.");
  // The Clp library has this off by default.  This program has it on by default.
    parameters.push_back(p);
  }
  {
    ClpParam p("PFI", "Whether to use Product Form of Inverse in simplex",
      "off", CLP_PARAM_STR_PFI, 7, 0);
    p.append("on");
    p.setLonghelp(
      "By default clp uses Forrest-Tomlin L-U update.  If you are masochistic you can switch it off.");
    parameters.push_back(p);
  }
  {
    ClpParam p("plus!Minus", "Tries to make +- 1 matrix",
      CLP_PARAM_ACTION_PLUSMINUS, 7, 0);
    p.setLonghelp(
      "Clp will go slightly faster if the matrix can be converted so that the elements are\
 not stored and are known to be unit.  The main advantage is memory use.  Clp may automatically\
 see if it can convert the problem so you should not need to use this.");
    parameters.push_back(p);
  }
  {
    ClpParam p("pO!ptions", "Dubious print options",
      0, COIN_INT_MAX, CLP_PARAM_INT_PRINTOPTIONS, 1);
    p.setIntValue(0);
    p.setLonghelp(
      "If this is > 0 then presolve will give more information and branch and cut will give statistics");
    parameters.push_back(p);
  }
  {
    ClpParam p("preO!pt", "Presolve options",
      0, COIN_INT_MAX, CLP_PARAM_INT_PRESOLVEOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("presolve", "Whether to presolve problem",
      "on", CLP_PARAM_STR_PRESOLVE);
    p.append("off");
    p.append("more");
    p.append("file");
    p.setLonghelp("Presolve analyzes the model to find such things as redundant equations, equations\
 which fix some variables, equations which can be transformed into bounds, etc.  For the\
 initial solve of any problem this is worth doing unless one knows that it will have no effect.  \
Option 'on' will normally do 5 passes, while using 'more' will do 10.  If the problem is very large one can \
let CLP write the original problem to file by using 'file'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("preT!olerance", "Tolerance to use in presolve",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRESOLVETOLERANCE);
    p.setLonghelp(
      "One may want to increase this tolerance if presolve says the problem is \
infeasible and one has awkward numbers and is sure that the problem is really feasible.");
    parameters.push_back(p);
  }
  {
    ClpParam p("primalP!ivot", "Primal pivot choice algorithm",
      "auto!matic", CLP_PARAM_STR_PRIMALPIVOT, 7, 1);

    p.append("exa!ct");
    p.append("dant!zig");
    p.append("part!ial");
    p.append("steep!est");
    p.append("change");
    p.append("sprint");
    p.append("PEsteep!est");
    p.append("PEdantzig");
    p.setLonghelp(
      "The Dantzig method is simple but its use is deprecated.  Exact devex is the method of choice and there\
 are two variants which keep all weights updated but only scan a subset each iteration.\
 Partial switches this on while 'change' initially does 'dantzig' until the factorization\
 becomes denser.  This is still a work in progress.\
 The PE variants add the Positive Edge criterion.\
 This selects incoming variables to try to avoid degenerate moves. \
 See also Towhidi, M., Desrosiers, J., Soumis, F., The positive edge criterion within COIN-OR's CLP;\
 Omer, J., Towhidi, M., Soumis, F., The positive edge pricing rule for the dual simplex.");

    parameters.push_back(p);
  }
  {
    ClpParam p("primalS!implex", "Do primal simplex algorithm",
      CLP_PARAM_ACTION_PRIMALSIMPLEX);

    p.setLonghelp(
      "This command solves the continuous relaxation of the current model using the primal algorithm.\
  The default is to use exact devex.\
 The time and iterations may be affected by settings such as presolve, scaling, crash\
 and also by column selection  method, infeasibility weight and dual and primal tolerances.");
    parameters.push_back(p);
  }
  {
    ClpParam p("primalT!olerance", "For a feasible solution \
no primal infeasibility, i.e., constraint violation, may exceed this value",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRIMALTOLERANCE);
    p.setLonghelp(
      "Normally the default tolerance is fine, but one may want to increase it a\
 bit if the primal simplex algorithm seems to be having a hard time.");
    parameters.push_back(p);
  }
  {
    ClpParam p("primalW!eight", "Initially algorithm acts as if it \
costs this much to be infeasible",
      1.0e-20, COIN_DBL_MAX, CLP_PARAM_DBL_PRIMALWEIGHT);
    p.setLonghelp(
      "The primal algorithm in Clp is a single phase algorithm as opposed to a two phase\
 algorithm where you first get feasible then optimal.  So Clp is minimizing this weight times\
 the sum of primal infeasibilities plus the true objective function (in minimization sense).\
  Too high a value may mean more iterations, while too low a value means\
 the algorithm may iterate into the wrong directory for long and then has to increase the weight in order to get feasible."); // OSL had a heuristic to adjust bounds, maybe we need that here.
    parameters.push_back(p);
  }
  {
    ClpParam p("psi", "Two-dimension pricing factor for Positive Edge criterion",
      -1.1, 1.1, CLP_PARAM_DBL_PSI);

    p.setDoubleValue(-0.5);
    p.setLonghelp(
      "The Positive Edge criterion has been added to \
select incoming variables to try and avoid degenerate moves. \
Variables not in the promising set have their infeasibility weight multiplied by psi, \
so 0.01 would mean that if there were any promising variables, then they would always be chosen, \
while 1.0 effectively switches the algorithm off. \
There are two ways of switching this feature on.  One way is to set psi to a positive value and then \
the Positive Edge criterion will be used for both primal and dual simplex.  The other way is to select PEsteepest \
in dualpivot choice (for example), then the absolute value of psi is used. \
Code donated by Jeremy Omer.  See \
Towhidi, M., Desrosiers, J., Soumis, F., The positive edge criterion within COIN-OR's CLP; \
Omer, J., Towhidi, M., Soumis, F., The positive edge pricing rule for the dual simplex.");  // Until this settles down it is only implemented in CLP.
    parameters.push_back(p);
  }
  {
    ClpParam p("printi!ngOptions", "Print options",
      "normal", CLP_PARAM_STR_INTPRINT, 3);
    p.append("integer");
    p.append("special");
    p.append("rows");
    p.append("all");
    p.append("csv");
    p.append("bound!ranging");
    p.append("rhs!ranging");
    p.append("objective!ranging");
    p.append("stats");
    p.append("boundsint");
    p.append("boundsall");
    p.append("fixint");
    p.append("fixall");
    p.setLonghelp(
      "This changes the amount and format of printing a solution:\nnormal - nonzero column variables \n\
integer - nonzero integer column variables\n\
special - in format suitable for OsiRowCutDebugger\n\
rows - nonzero column variables and row activities\n\
all - all column variables and row activities.\n\
\nFor non-integer problems 'integer' and 'special' act like 'normal'.  \
Also see printMask for controlling output.");
    parameters.push_back(p);
  }
  {
    ClpParam p("printM!ask", "Control printing of solution on a  mask",
      CLP_PARAM_ACTION_PRINTMASK, 3);

    p.setLonghelp(
      "If set then only those names which match mask are printed in a solution. \
'?' matches any character and '*' matches any set of characters. \
 The default is '' i.e. unset so all variables are printed. \
This is only active if model has names.");
    parameters.push_back(p);
  }

  {
    ClpParam p("progress!(Interval)", "Time interval for printing progress",
		    -COIN_DBL_MAX,COIN_DBL_MAX,
		    CLP_PARAM_DBL_PROGRESS);
    p.setLonghelp(
      "This sets a minimum interval for some printing - elapsed seconds");

    p.setDoubleValue(0.7);
    parameters.push_back(p);
  }
  {
    ClpParam p("quit", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");

    parameters.push_back(p);
  }
  {
    ClpParam p("randomS!eed", "Random seed for Clp",
      0, COIN_INT_MAX, CLP_PARAM_INT_RANDOMSEED);

    p.setLonghelp(
      "Initialization of the random seed for pseudo-random numbers used to break ties in degenerate problems. "
      "This may yield a different continuous optimum and, in the context of Cbc, different cuts and heuristic solutions. "
      "The special value of 0 lets CLP use the time of the day for the initial seed.");
    p.setIntValue(1234567);
    parameters.push_back(p);
  }
  {
    ClpParam p("restoreS!olution", "reads solution from file",
      CLP_PARAM_ACTION_RESTORESOL);

    p.setLonghelp(
      "This will read a binary solution file from the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  This reads in a file from saveSolution");
    parameters.push_back(p);
  }
  {
    ClpParam p("reallyO!bjectiveScale", "Scale factor to apply to objective in place",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_OBJSCALE2, 0);
    p.setLonghelp("You can set this to -1.0 to test maximization or other to stress code");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
  {
    ClpParam p("reallyS!cale", "Scales model in place",
      CLP_PARAM_ACTION_REALLY_SCALE, 7, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("restore!Model", "Restore model from binary file",
      CLP_PARAM_ACTION_RESTORE, 7, 1);
    p.setLonghelp(
      "This reads data save by saveModel from the given file.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");

    parameters.push_back(p);
  }
  {
    ClpParam p("reverse", "Reverses sign of objective",
      CLP_PARAM_ACTION_REVERSE, 7, 0);
    p.setLonghelp(
      "Useful for testing if maximization works correctly");
    parameters.push_back(p);
  }
  {
    ClpParam p("rhs!Scale", "Scale factor to apply to rhs and bounds",
      -COIN_DBL_MAX, COIN_DBL_MAX, CLP_PARAM_DBL_RHSSCALE, 0);
    p.setLonghelp(
      "If the rhs or bounds have some very large meaningful values, you may wish to scale them\
 internally by this amount.  It can also be set by autoscale.  This should not be needed.");
    p.setDoubleValue(1.0);
    parameters.push_back(p);
  }
  {
    ClpParam p("saveM!odel", "Save model to binary file",
      CLP_PARAM_ACTION_SAVE, 7, 1);
    p.setLonghelp(
      "This will save the problem to the given file name for future use\
 by restoreModel.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'default.prob'.");
    parameters.push_back(p);
  }
  {
    ClpParam p("saveS!olution", "saves solution to file",
      CLP_PARAM_ACTION_SAVESOL);

    p.setLonghelp(
      "This will write a binary solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'solution.file'.  To read the file use fread(int) twice to pick up number of rows \
and columns, then fread(double) to pick up objective value, then pick up row activities, row duals, column \
activities and reduced costs - see bottom of ClpParam.cpp for code that reads or writes file. \
If name contains '_fix_read_' then does not write but reads and will fix all variables");
    parameters.push_back(p);
  }
  {
    ClpParam p("scal!ing", "Whether to scale problem",
      "off", CLP_PARAM_STR_SCALING);
    p.append("equi!librium");
    p.append("geo!metric");
    p.append("auto!matic");
    p.append("dynamic");
    p.append("rows!only");
    p.setLonghelp(
      "Scaling can help in solving problems which might otherwise fail because of lack of\
 accuracy.  It can also reduce the number of iterations.  It is not applied if the range\
 of elements is small.  When the solution is evaluated in the unscaled problem, it is possible that small primal and/or\
 dual infeasibilities occur. "
 "Option 'equilibrium' uses the largest element for scaling. "
 "Option 'geometric' uses the squareroot of the product of largest and smallest element. "
 "Option 'auto' let CLP choose a method that gives the best ratio of the largest element to the smallest one.");
    p.setCurrentOption(3); // say auto
    parameters.push_back(p);
  }
  {
    ClpParam p("sec!onds", "Maximum seconds",
      -1.0, COIN_DBL_MAX, CLP_PARAM_DBL_TIMELIMIT);

    p.setLonghelp("After this many seconds clp will act as if maximum iterations had been reached \
(if value >=0).");
    parameters.push_back(p);
  }
  {
    ClpParam p("sleep", "for debug",
      CLP_PARAM_ACTION_DUMMY, 7, 0);

    p.setLonghelp(
      "If passed to solver fom ampl, then ampl will wait so that you can copy .nl file for debug.");
    parameters.push_back(p);
  }
  {
    ClpParam p("slp!Value", "Number of slp passes before primal",
      -50000, 50000, CLP_PARAM_INT_SLPVALUE, 1);
    p.setLonghelp(
      "If you are solving a quadratic problem using primal then it may be helpful to do some \
sequential Lps to get a good approximate solution.");
    parameters.push_back(p);
  }
#if CLP_MULTIPLE_FACTORIZATIONS > 0
  {
    ClpParam p("small!Factorization", "Threshold for using small factorization",
      -1, 10000, CLP_PARAM_INT_SMALLFACT, 1);
    p.setLonghelp(
      "If processed problem <= this use small factorization");
    p.setIntValue(-1);
    parameters.push_back(p);
  }
#endif
  {
    ClpParam p("solu!tion", "Prints solution to file",
      CLP_PARAM_ACTION_SOLUTION);
    p.setLonghelp(
      "This will write a primitive solution file to the given file name.  It will use the default\
 directory given by 'directory'.  A name of '$' will use the previous value for the name.  This\
 is initialized to 'stdout'.  The amount of output can be varied using printi!ngOptions or printMask.");
    parameters.push_back(p);
  }
  // allow solve as synonym for possible dual
  {
    ClpParam p("solv!e", "Solve problem using dual simplex (probably)",
      CLP_PARAM_ACTION_EITHERSIMPLEX);
    p.setLonghelp(
      "Solves ");
    parameters.push_back(p);
  }
  {
    ClpParam p("spars!eFactor", "Whether factorization treated as sparse",
      "on", CLP_PARAM_STR_SPARSEFACTOR, 7, 0);
    p.append("off");
    parameters.push_back(p);
  }
  {
    ClpParam p("special!Options", "Dubious options for Simplex - see ClpSimplex.hpp",
      0, COIN_INT_MAX, CLP_PARAM_INT_SPECIALOPTIONS, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("sprint!Crash", "Whether to try sprint crash",
      -1, COIN_INT_MAX, CLP_PARAM_INT_SPRINT);
    p.setLonghelp(
      "For long and thin problems this method may solve a series of small problems\
 created by taking a subset of the columns.  The idea as 'Sprint' was introduced by J. Forrest after\
 an LP code of that name of the 60's which tried the same tactic (not totally successfully).\
  CPLEX calls it 'sifting'.  -1 lets CLP automatically choose the number of passes, 0 is off, n is number of passes");
    parameters.push_back(p);
  }
  {
    ClpParam p("stat!istics", "Print some statistics",
      CLP_PARAM_ACTION_STATISTICS);
    p.setLonghelp(
      "This command prints some statistics for the current model.\
 If log level >1 then more is printed.\
 These are for presolved model if presolve on (and unscaled).");
    parameters.push_back(p);
  }
  {
    ClpParam p("stop", "Stops clp execution",
      CLP_PARAM_ACTION_EXIT);
    p.setLonghelp(
      "This stops the execution of Clp, end, exit, quit and stop are synonyms");
    parameters.push_back(p);
  }
  {
    ClpParam p("subs!titution", "How long a column to substitute for in presolve",
      0, 10000, CLP_PARAM_INT_SUBSTITUTION, 0);
    p.setLonghelp(
      "Normally Presolve gets rid of 'free' variables when there are no more than 3 \
 coefficients in a row.  If you increase this, the number of rows may decrease but the number of \
 coefficients may increase.");
    parameters.push_back(p);
  }
#ifdef CLP_THREAD
  {
    ClpParam p("thread!s", "Number of threads to try and use",
      -100, 100000, CBC_PARAM_INT_THREADS, 1);
    p.setIntValue(0);
    p.setLonghelp(
      "To use multiple threads, set threads to number wanted.  It may be better \
to use one or two more than number of cpus available.  If 100+n then n threads and \
search is repeatable (maybe be somewhat slower), \
if 200+n use threads for root cuts, 400+n threads used in sub-trees.");
    parameters.push_back(p);
  }
#endif
  {
    ClpParam p("tightLP", "Poor person's preSolve for now",
      CLP_PARAM_ACTION_TIGHTEN, 7, 0);
    parameters.push_back(p);
  }
  {
    ClpParam p("unitTest", "Do unit test",
      CLP_PARAM_ACTION_UNITTEST, 3, 1);
    p.setLonghelp(
      "This exercises the unit test for clp");
    parameters.push_back(p);
  }
  {
    ClpParam p("userClp", "Hand coded Clp stuff",
      CLP_PARAM_ACTION_USERCLP, 0, 0);
    p.setLonghelp(
      "There are times e.g. when using AMPL interface when you may wish to do something unusual.  \
Look for USERCLP in main driver and modify sample code.");
    parameters.push_back(p);
  }
  {
#ifndef COIN_AVX2
    ClpParam p("vector", "Whether to use vector? Form of matrix in simplex",
      "off", CLP_PARAM_STR_VECTOR, 7, 0);
    p.append("on");
    p.setLonghelp(
      "If this is on ClpPackedMatrix uses extra column copy in odd format.");
    parameters.push_back(p);
#else
    ClpParam p("vector", "Try and use vector instructions in simplex",
      "off", CLP_PARAM_STR_VECTOR, 7, 0);
    p.append("on");
    p.append("ones");
    p.setLonghelp(
      "At present only for Intel architectures - but could be extended.  \
Uses avx2 or avx512 instructions. Uses different storage for matrix - can be \
of benefit without instruction set on some problems.  \
I may add pool to switch on a pool matrix");
    parameters.push_back(p);
#endif
  }
  {
    ClpParam p("verbose", "Switches on longer help on single ?",
      0, 31, CLP_PARAM_INT_VERBOSE, 0);
    p.setLonghelp(
      "Set to 1 to get short help with ? list, 2 to get long help, 3 for both.  (add 4 to just get ampl ones).");
    p.setIntValue(0);
    parameters.push_back(p);
  }
  {
    ClpParam p("zeroT!olerance", "Kill all coefficients \
whose absolute value is less than this value",
      1.0e-100, 1.0e-5, CLP_PARAM_DBL_ZEROTOLERANCE);
    p.setLonghelp(
      "This applies to reading mps files (and also lp files \
if KILL_ZERO_READLP defined)");
    p.setDoubleValue(1.0e-20);
    parameters.push_back(p);
  }
}
// Given a parameter type - returns its number in list
int whichClpParam(const ClpParameterType &name,
  const std::vector< ClpParam > &parameters)
{
  for (int i = 0; i < (int)parameters.size(); i++) {
    if (parameters[i].type() == name)
      return i;
  }
  return std::numeric_limits< int >::max(); // should not arrive here
}
/* Restore a solution from file.
   mode 0 normal, 1 swap rows and columns and primal and dual
   if 2 set then also change signs
*/
void restoreSolution(ClpSimplex *lpSolver, std::string fileName, int mode)
{
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    int numberRowsFile;
    int numberColumnsFile;
    double objectiveValue;
    size_t nRead;
    nRead = fread(&numberRowsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&numberColumnsFile, sizeof(int), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    nRead = fread(&objectiveValue, sizeof(double), 1, fp);
    if (nRead != 1)
      throw("Error in fread");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    if (mode) {
      // swap
      int k = numberRows;
      numberRows = numberColumns;
      numberColumns = k;
      double *temp;
      temp = dualRowSolution;
      dualRowSolution = primalColumnSolution;
      primalColumnSolution = temp;
      temp = dualColumnSolution;
      dualColumnSolution = primalRowSolution;
      primalRowSolution = temp;
    }
    if (numberRows > numberRowsFile || numberColumns > numberColumnsFile) {
      std::cout << "Mismatch on rows and/or columns - giving up" << std::endl;
    } else {
      lpSolver->setObjectiveValue(objectiveValue);
      if (numberRows == numberRowsFile && numberColumns == numberColumnsFile) {
        nRead = fread(primalRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(dualRowSolution, sizeof(double), numberRows, fp);
        if (nRead != static_cast< size_t >(numberRows))
          throw("Error in fread");
        nRead = fread(primalColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
        nRead = fread(dualColumnSolution, sizeof(double), numberColumns, fp);
        if (nRead != static_cast< size_t >(numberColumns))
          throw("Error in fread");
      } else {
        std::cout << "Mismatch on rows and/or columns - truncating" << std::endl;
        double *temp = new double[CoinMax(numberRowsFile, numberColumnsFile)];
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, primalRowSolution);
        nRead = fread(temp, sizeof(double), numberRowsFile, fp);
        if (nRead != static_cast< size_t >(numberRowsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberRows, dualRowSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, primalColumnSolution);
        nRead = fread(temp, sizeof(double), numberColumnsFile, fp);
        if (nRead != static_cast< size_t >(numberColumnsFile))
          throw("Error in fread");
        CoinMemcpyN(temp, numberColumns, dualColumnSolution);
        delete[] temp;
      }
      if (mode == 3) {
        int i;
        for (i = 0; i < numberRows; i++) {
          primalRowSolution[i] = -primalRowSolution[i];
          dualRowSolution[i] = -dualRowSolution[i];
        }
        for (i = 0; i < numberColumns; i++) {
          primalColumnSolution[i] = -primalColumnSolution[i];
          dualColumnSolution[i] = -dualColumnSolution[i];
        }
      }
    }
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}
// Dump a solution to file
void saveSolution(const ClpSimplex *lpSolver, std::string fileName)
{
  if (strstr(fileName.c_str(), "_fix_read_")) {
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (fp) {
      ClpSimplex *solver = const_cast< ClpSimplex * >(lpSolver);
      restoreSolution(solver, fileName, 0);
      // fix all
      int logLevel = solver->logLevel();
      int iColumn;
      int numberColumns = solver->numberColumns();
      double *primalColumnSolution = solver->primalColumnSolution();
      double *columnLower = solver->columnLower();
      double *columnUpper = solver->columnUpper();
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        double value = primalColumnSolution[iColumn];
        if (value > columnUpper[iColumn]) {
          if (value > columnUpper[iColumn] + 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnUpper[iColumn];
        } else if (value < columnLower[iColumn]) {
          if (value < columnLower[iColumn] - 1.0e-6 && logLevel > 1)
            printf("%d value of %g - bounds %g %g\n",
              iColumn, value, columnLower[iColumn], columnUpper[iColumn]);
          value = columnLower[iColumn];
        }
        columnLower[iColumn] = value;
        columnUpper[iColumn] = value;
      }
      return;
    }
  }
  FILE *fp = fopen(fileName.c_str(), "wb");
  if (fp) {
    int numberRows = lpSolver->numberRows();
    int numberColumns = lpSolver->numberColumns();
    double objectiveValue = lpSolver->objectiveValue();
    size_t nWrite;
    nWrite = fwrite(&numberRows, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&numberColumns, sizeof(int), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    nWrite = fwrite(&objectiveValue, sizeof(double), 1, fp);
    if (nWrite != 1)
      throw("Error in fwrite");
    double *dualRowSolution = lpSolver->dualRowSolution();
    double *primalRowSolution = lpSolver->primalRowSolution();
    nWrite = fwrite(primalRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    nWrite = fwrite(dualRowSolution, sizeof(double), numberRows, fp);
    if (nWrite != static_cast< size_t >(numberRows))
      throw("Error in fwrite");
    double *dualColumnSolution = lpSolver->dualColumnSolution();
    double *primalColumnSolution = lpSolver->primalColumnSolution();
    nWrite = fwrite(primalColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    nWrite = fwrite(dualColumnSolution, sizeof(double), numberColumns, fp);
    if (nWrite != static_cast< size_t >(numberColumns))
      throw("Error in fwrite");
    fclose(fp);
  } else {
    std::cout << "Unable to open file " << fileName << std::endl;
  }
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
