// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

#include <time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

#define CLPVERSION "0.92"

//#include "CoinPackedMatrix.hpp"
//#include "CoinPackedVector.hpp"
#include "CoinMpsIO.hpp"

#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"

#include "Presolve.hpp"
#ifdef CLP_IDIOT
#include "Idiot.hpp"
#endif
// For Branch and bound
//  #include "OsiClpSolverInterface.hpp"
//  #include "OsiCuts.hpp"
//  #include "OsiRowCut.hpp"
//  #include "OsiColCut.hpp"

#ifdef DMALLOC
#include "dmalloc.h"
#endif
static double totalTime=0.0;
static double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow;        /* clock_t is same as int */
  
  ticksnow = (unsigned int)clock();
  
  cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}


//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

int mainTest (int argc, const char *argv[],bool doDual,
	      ClpSimplex empty, bool doPresolve,int doIdiot);
enum ClpParameterType {
  GENERALQUERY=-100,
  
  DUALTOLERANCE=1,PRIMALTOLERANCE,DUALBOUND,PRIMALWEIGHT,

  LOGLEVEL=101,MAXFACTOR,PERTURBATION,MAXITERATION,PRESOLVEPASS,IDIOT,
  
  DIRECTION=201,DUALPIVOT,SCALING,ERRORSALLOWED,KEEPNAMES,SPARSEFACTOR,
  PRIMALPIVOT,PRESOLVE,
  
  DIRECTORY=301,IMPORT,EXPORT,RESTORE,SAVE,DUALSIMPLEX,PRIMALSIMPLEX,BAB,
  MAXIMIZE,MINIMIZE,EXIT,STDIN,UNITTEST,NETLIB_DUAL,NETLIB_PRIMAL,SOLUTION,
  TIGHTEN,FAKEBOUND,VERSION,

  INVALID=1000
};
/// Very simple class for setting parameters
class ClpItem {

public:

  /**@name Constructor and destructor */
  //@{
  /// Constructors
  ClpItem (  );
  ClpItem (std::string name, std::string help,
	   double lower, double upper, ClpParameterType type);
  ClpItem (std::string name, std::string help,
	   int lower, int upper, ClpParameterType type);
  // Other strings will be added by insert
  ClpItem (std::string name, std::string help, std::string defaultValue,
	   ClpParameterType type);
  // Action
  ClpItem (std::string name, std::string help,
	   ClpParameterType type);
  /// Copy constructor. 
  ClpItem(const ClpItem &);
  /// Assignment operator. This copies the data
    ClpItem & operator=(const ClpItem & rhs);
  /// Destructor
  ~ClpItem (  );
  //@}

  /**@name stuff */
  //@{
  /// Insert string (only valid for keywords)
  void append(std::string keyWord);
  /// Adds one help line
  void addHelp(std::string keyWord);
  /// Returns name
  inline std::string  name(  ) const {
    return name_;
  };
  /// Returns short help
  inline std::string  shortHelp(  ) const {
    return shortHelp_;
  };
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(ClpSimplex * model, double value) const;
  /// Gets a double parameter
  double doubleParameter(ClpSimplex * model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(ClpSimplex * model, int value) const;
  /// Gets a int parameter
  int intParameter(ClpSimplex * model) const;
  /// Returns name which could match
  std::string matchName (  ) const;
  /// Returns parameter option which matches (-1 if none)
  int parameterOption ( std::string check ) const;
  /// Prints parameter options
  void printOptions (  ) const;
  /// Returns current parameter option
  inline std::string currentOption (  ) const
  { return definedKeyWords_[currentKeyWord_]; };
  /// Sets current parameter option
  inline void setCurrentOption ( int value )
  { currentKeyWord_=value; };
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches (std::string input) const;
  /// type
  inline ClpParameterType type() const
  { return type_;};
private:
  /// gutsOfConstructor
  void gutsOfConstructor();
  //@}
////////////////// data //////////////////
private:

  /**@name data 
   We might as well throw all type data in - could derive?
  */
  //@{
  // Type see ClpParameterType
  ClpParameterType type_;
  /// If double == okay
  double lowerDoubleValue_;
  double upperDoubleValue_;
  /// If int == okay
  int lowerIntValue_;
  int upperIntValue_;
  // Length of name
  unsigned int lengthName_;
  // Minimum match
  unsigned int lengthMatch_;
  /// set of valid strings
  std::vector<std::string> definedKeyWords_;
  /// Name
  std::string name_;
  /// Short help
  std::string shortHelp_;
  /// Long help
  std::vector<std::string> longHelp_;
  /// Action
  ClpParameterType action_;
  /// Current keyWord (if a keyword parameter)
  int currentKeyWord_;
  //@}
};
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpItem::ClpItem () 
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
    currentKeyWord_(-1)
{
}
// Other constructors
ClpItem::ClpItem (std::string name, std::string help,
	   double lower, double upper, ClpParameterType type)
  : type_(type),
    lowerIntValue_(0),
    upperIntValue_(0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1)
{
  lowerDoubleValue_ = lower;
  upperDoubleValue_ = upper;
  gutsOfConstructor();
}
ClpItem::ClpItem (std::string name, std::string help,
	   int lower, int upper, ClpParameterType type)
  : type_(type),
    lowerDoubleValue_(0.0),
    upperDoubleValue_(0.0),
    definedKeyWords_(),
    name_(name),
    shortHelp_(help),
    longHelp_(),
    action_(type),
    currentKeyWord_(-1)
{
  gutsOfConstructor();
  lowerIntValue_ = lower;
  upperIntValue_ = upper;
}
// Other strings will be added by append
ClpItem::ClpItem (std::string name, std::string help, 
		  std::string defaultValue,
		  ClpParameterType type)
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
    currentKeyWord_(0)
{
  gutsOfConstructor();
  definedKeyWords_.push_back(defaultValue);
}
// Action
ClpItem::ClpItem (std::string name, std::string help,
	   ClpParameterType type)
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
    currentKeyWord_(-1)
{
  gutsOfConstructor();
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpItem::ClpItem (const ClpItem & rhs) 
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
  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpItem::~ClpItem ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpItem &
ClpItem::operator=(const ClpItem& rhs)
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
  }
  return *this;
}
void 
ClpItem::gutsOfConstructor()
{
  unsigned int  shriekPos = name_.find('!');
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
ClpItem::append(std::string keyWord)
{
  definedKeyWords_.push_back(keyWord);
}

int 
ClpItem::matches (std::string input) const
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
ClpItem::matchName (  ) const
{ 
  if (lengthMatch_==lengthName_) 
    return name_;
  else
    return name_.substr(0,lengthMatch_)+"("+name_.substr(lengthMatch_)+")";
}

// Returns parameter option which matches (-1 if none)
int 
ClpItem::parameterOption ( std::string check ) const
{
  int numberItems = definedKeyWords_.size();
  if (!numberItems) {
    return -1;
  } else {
    int whichItem=0;
    unsigned int it;
    for (it=0;it<definedKeyWords_.size();it++) {
      std::string thisOne = definedKeyWords_[it];
      unsigned int  shriekPos = thisOne.find('!');
      unsigned int length1 = thisOne.length();
      unsigned int length2 = length1;
      if ( shriekPos!=std::string::npos ) {
	//contains '!'
	length2 = shriekPos;
	thisOne = thisOne.substr(0,shriekPos)+
	  thisOne.substr(shriekPos+1);
	length1 = thisOne.length();
      }
      if (check.length()<=length1) {
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
ClpItem::printOptions (  ) const
{
  std::cout<<"Possible options for "<<name_<<" are:"<<std::endl;
  unsigned int it;
  for (it=0;it<definedKeyWords_.size();it++) {
    std::string thisOne = definedKeyWords_[it];
    unsigned int  shriekPos = thisOne.find('!');
    if ( shriekPos!=std::string::npos ) {
      //contains '!'
      thisOne = thisOne.substr(0,shriekPos)+
	"("+thisOne.substr(shriekPos+1)+")";
    }
    std::cout<<thisOne<<std::endl;
  }
}
int
ClpItem::setDoubleParameter (ClpSimplex * model,double value) const
{
  double oldValue = doubleParameter(model);
  if (value<lowerDoubleValue_||value>upperDoubleValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerDoubleValue_<<" to "<<
      upperDoubleValue_<<std::endl;
    return 1;
  } else {
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    switch(type_) {
    case DUALTOLERANCE:
      model->setDualTolerance(value);
      break;
    case PRIMALTOLERANCE:
      model->setPrimalTolerance(value);
      break;
    case DUALBOUND:
      model->setDualBound(value);
      break;
    case PRIMALWEIGHT:
      model->setInfeasibilityCost(value);
      break;
    default:
      abort();
    }
    return 0;
  }
}
double 
ClpItem::doubleParameter (ClpSimplex * model) const
{
  double value;
  switch(type_) {
  case DUALTOLERANCE:
    value=model->dualTolerance();
    break;
  case PRIMALTOLERANCE:
    value=model->primalTolerance();
    break;
  case DUALBOUND:
    value=model->dualBound();
    break;
  case PRIMALWEIGHT:
    value=model->infeasibilityCost();
    break;
  default:
    abort();
  }
  return value;
}
int 
ClpItem::setIntParameter (ClpSimplex * model,int value) const
{
  int oldValue = intParameter(model);
  if (value<lowerIntValue_||value>upperIntValue_) {
    std::cout<<value<<" was provided for "<<name_<<
      " - valid range is "<<lowerIntValue_<<" to "<<
      upperIntValue_<<std::endl;
    return 1;
  } else {
    std::cout<<name_<<" was changed from "<<oldValue<<" to "
	     <<value<<std::endl;
    switch(type_) {
    case LOGLEVEL:
      model->setLogLevel(value);
      break;
    case MAXFACTOR:
      model->factorization()->maximumPivots(value);
      break;
    case DIRECTION:
      model->setOptimizationDirection(value);
      break;
    case PERTURBATION:
      model->setPerturbation(value);
      break;
    case MAXITERATION:
      model->setMaximumIterations(value);
      break;
    default:
      abort();
    }
    return 0;
  }
}
int 
ClpItem::intParameter (ClpSimplex * model) const
{
  int value;
  switch(type_) {
  case LOGLEVEL:
    value=model->logLevel();
    break;
  case MAXFACTOR:
    value=model->factorization()->maximumPivots();
    break;
  case DIRECTION:
    value=model->optimizationDirection();
    break;
  case PERTURBATION:
    value=model->perturbation();
    break;
  case MAXITERATION:
    value=model->maximumIterations();
    break;
  default:
    abort();
  }
  return value;
}
#ifdef READLINE     
#include "/usr/include/readline/readline.h"
#include "/usr/include/readline/history.h"
#include <signal.h>
static ClpSimplex * currentModel = NULL;
static void signal_handler(int whichSignal)
{
  if (currentModel!=NULL) 
    currentModel->setMaximumIterations(0); // stop at next iterations
  return;
}
#endif
// Returns next valid field
static int read_mode=1;
static char line[1000];
static char * where=NULL;

std::string
nextField()
{
  std::string field;
  if (!where) {
    // need new line
#ifdef READLINE     
    // Get a line from the user. 
    where = readline ("Clp:");
     
    // If the line has any text in it, save it on the history.
    if (where) {
      if ( *where)
	add_history (where);
      strcpy(line,where);
    }
#else
    fprintf(stdout,"Clp:");
    fflush(stdout);
    where = fgets(line,1000,stdin);
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
getCommand(int argc, const char *argv[])
{
  std::string field="EOL";
  while (field=="EOL") {
    if (read_mode>0) {
      if (read_mode<argc) {
	field = argv[read_mode++];
	if (field=="-") {
	  std::cout<<"Switching to line mode"<<std::endl;
	  read_mode=-1;
	  field=nextField();
	} else if (field[0]!='-') {
	  if (read_mode!=2) {
	    std::cout<<"skipping non-command "<<field<<std::endl;
	    field="EOL"; // skip
	  } else {
	    // special dispensation - taken as -import name
	    read_mode--;
	    field="import";
	  }
	} else {
	  if (field!="--") {
	    // take off -
	    field = field.substr(1);
	  } else {
	    // special dispensation - taken as -import --
	    read_mode--;
	    field="import";
	  }
	}
      } else {
	field="";
      }
    } else {
      field=nextField();
    }
  }
  //std::cout<<field<<std::endl;
  return field;
}
std::string
getString(int argc, const char *argv[])
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      if (argv[read_mode][0]!='-') { 
	field = argv[read_mode++];
      } else if (!strcmp(argv[read_mode],"--")) {
	field = argv[read_mode++];
	// -- means import from stdin
	field = "-";
      }
    }
  } else {
    field=nextField();
  }
  //std::cout<<field<<std::endl;
  return field;
}
// valid 0 - okay, 1 bad, 2 not there
int
getIntField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      // may be negative value so do not check for -
      field = argv[read_mode++];
    }
  } else {
    field=nextField();
  }
  int value=0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value =  atoi(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}
double
getDoubleField(int argc, const char *argv[],int * valid)
{
  std::string field="EOL";
  if (read_mode>0) {
    if (read_mode<argc) {
      // may be negative value so do not check for -
      field = argv[read_mode++];
    }
  } else {
    field=nextField();
  }
  double value=0.0;
  //std::cout<<field<<std::endl;
  if (field!="EOL") {
    // how do I check valid
    value = atof(field.c_str());
    *valid=0;
  } else {
    *valid=2;
  }
  return value;
}
int main (int argc, const char *argv[])
{
  // next {} is just to make sure all memory should be freed - for debug
  {
    double time1 = cpuTime(),time2;
#define MAXPARAMETERS 100
    ClpItem parameters[MAXPARAMETERS];
    int numberParameters=0;
    parameters[numberParameters++]=
      ClpItem("?","For help",GENERALQUERY);
    parameters[numberParameters++]=
      ClpItem("dualT!olerance","For an optimal solution \
no dual infeasibility may exceed this value",
	      1.0e-20,1.0e12,DUALTOLERANCE);
    parameters[numberParameters++]=
      ClpItem("primalT!olerance","For an optimal solution \
no primal infeasibility may exceed this value",
	      1.0e-20,1.0e12,PRIMALTOLERANCE);
    parameters[numberParameters++]=
      ClpItem("dualB!ound","Initially algorithm acts as if no \
gap between bounds exceeds this value",
	      1.0e-20,1.0e12,DUALBOUND);
    parameters[numberParameters++]=
      ClpItem("primalW!eight","Initially algorithm acts as if it \
costs this much to be infeasible",
	      1.0e-20,1.0e12,PRIMALWEIGHT);
    parameters[numberParameters++]=
      ClpItem("log!Level","Level of detail in output",
	      0,63,LOGLEVEL);
    parameters[numberParameters++]=
      ClpItem("maxF!actor","Maximum number of iterations between \
refactorizations",
	      1,999999,MAXFACTOR);
    parameters[numberParameters++]=
      ClpItem("maxIt!erations","Maximum number of iterations before \
stopping",
	      0,99999999,MAXITERATION);
    parameters[numberParameters++]=
      ClpItem("pert!urbation","Method of perturbation",
	      -50,102,PERTURBATION);
    parameters[numberParameters++]=
      ClpItem("direction","Minimize or Maximize",
	      "min!imize",DIRECTION);
    parameters[numberParameters-1].append("max!imize");
    parameters[numberParameters++]=
      ClpItem("dualP!ivot","Dual pivot choice algorithm",
	      "steep!est",DUALPIVOT);
    parameters[numberParameters-1].append("dant!zig");
    parameters[numberParameters++]=
      ClpItem("primalP!ivot","Primal pivot choice algorithm",
	      "steep!est",PRIMALPIVOT);
    parameters[numberParameters-1].append("exa!ct");
    parameters[numberParameters-1].append("dant!zig");
    parameters[numberParameters++]=
      ClpItem("scal!ing","Whether to scale problem",
	      "on",SCALING);
    parameters[numberParameters-1].append("off");
    parameters[numberParameters++]=
      ClpItem("presolve","Whether to presolve problem",
	      "off",PRESOLVE);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      ClpItem("idiot!Crash","Whether to try idiot crash",
	      0,200,IDIOT);
    parameters[numberParameters++]=
      ClpItem("passP!resolve","How many passes in presolve",
	      0,100,PRESOLVEPASS);
    parameters[numberParameters++]=
      ClpItem("spars!eFactor","Whether factorization treated as sparse",
	      "on",SPARSEFACTOR);
    parameters[numberParameters-1].append("off");
    parameters[numberParameters++]=
      ClpItem("error!sAllowed","Whether to allow import errors",
	      "off",ERRORSALLOWED);
    parameters[numberParameters-1].append("on");
    parameters[numberParameters++]=
      ClpItem("keepN!ames","Whether to keep names from import",
	      "on",KEEPNAMES);
    parameters[numberParameters-1].append("off");
    parameters[numberParameters++]=
      ClpItem("directory","Set Default import directory",
	      DIRECTORY);
    parameters[numberParameters++]=
      ClpItem("import","Import model from mps file",
	      IMPORT);
    parameters[numberParameters++]=
      ClpItem("export","Export model as mps file",
	      EXPORT);
    parameters[numberParameters++]=
      ClpItem("save!Model","Save model to binary file",
	      SAVE);
    parameters[numberParameters++]=
      ClpItem("restore!Model","Restore model from binary file",
	      RESTORE);
    parameters[numberParameters++]=
      ClpItem("dualS!implex","Do dual simplex algorithm",
	      DUALSIMPLEX);
    parameters[numberParameters++]=
      ClpItem("primalS!implex","Do primal simplex algorithm",
	      PRIMALSIMPLEX);
    parameters[numberParameters++]=
      ClpItem("branch!Andbound","Do Branch and Bound",
	      BAB);
    parameters[numberParameters++]=
      ClpItem("tight!en","Poor person's preSolve for now",
	      TIGHTEN);
    parameters[numberParameters++]=
      ClpItem("sol!ution","Prints solution to file",
	      SOLUTION);
    parameters[numberParameters++]=
      ClpItem("max!imize","Set optimization direction to maximize",
	      MAXIMIZE);
    parameters[numberParameters++]=
      ClpItem("min!imize","Set optimization direction to minimize",
	      MINIMIZE);
    parameters[numberParameters++]=
      ClpItem("exit","Stops clp execution",
	      EXIT);
    parameters[numberParameters++]=
      ClpItem("stop","Stops clp execution",
	      EXIT);
    parameters[numberParameters++]=
      ClpItem("quit","Stops clp execution",
	      EXIT);
    parameters[numberParameters++]=
      ClpItem("-","From stdin",
	      STDIN);
    parameters[numberParameters++]=
      ClpItem("stdin","From stdin",
	      STDIN);
    parameters[numberParameters++]=
      ClpItem("unitTest","Do unit test",
	      UNITTEST);
    parameters[numberParameters++]=
      ClpItem("netlib","Solev entlib test set",
	      NETLIB_DUAL);
    parameters[numberParameters++]=
      ClpItem("netlibP!rimal","Solev entlib test set",
	      NETLIB_PRIMAL);
    parameters[numberParameters++]=
      ClpItem("fakeB!ound","All bounds <= this value - DEBUG",
	      1.0,1.0e15,FAKEBOUND);
    parameters[numberParameters++]=
      ClpItem("ver!sion","Print out version",
	      VERSION);
    assert(numberParameters<MAXPARAMETERS);
    
    // total number of commands read
    int numberGoodCommands=0;
    //int numberModels=1;
    ClpSimplex * models = new ClpSimplex[1];
    bool * goodModels = new bool[1];
    
#ifdef READLINE     
    currentModel = models;
    // register signal handler
    signal(SIGINT,signal_handler);
#endif
    
    // default action on import
    int allowImportErrors=0;
    int keepImportNames=1;
    int preSolve=0;
    int doIdiot=0;
    
    int iModel=0;
    goodModels[0]=false;
    // set reasonable defaults
    models[0].scaling(1);
    models[0].setDualBound(1.0e6);
    models[0].setDualTolerance(1.0e-7);
    ClpDualRowSteepest steep;
    models[0].setDualRowPivotAlgorithm(steep);
    models[0].setPrimalTolerance(1.0e-8);
    ClpPrimalColumnSteepest steepP;
    models[0].setPrimalColumnPivotAlgorithm(steepP);
    std::string directory ="./";
    std::string field;
    
    while (1) {
      // next command
      field=getCommand(argc,argv);
      
      // exit if null or similar
      if (!field.length()) {
	if (numberGoodCommands==1&&goodModels[0]) {
	  // we just had file name
	  models[0].scaling(1);
	  int saveMaxIterations = models[0].maximumIterations();
	  models[0].dual();
	  models[0].setMaximumIterations(saveMaxIterations);
	} else if (!numberGoodCommands) {
	  // let's give the sucker a hint
	  std::cout
	    <<"Clp takes input from arguments ( - switches to stdin)"
	    <<std::endl
	    <<"Enter ? for list of commands, (-)unitTest or (-)netlib"
	    <<" for tests"<<std::endl;
	}
	break;
      }
      
      // see if ? at end
      int numberQuery=0;
      if (field!="?") {
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
      for ( iParam=0; iParam<numberParameters; iParam++ ) {
	int match = parameters[iParam].matches(field);
	if (match==1) {
	  numberMatches = 1;
	  break;
	} else {
	  numberMatches += match>>1;
	}
      }
      if (iParam<numberParameters&&!numberQuery) {
	// found
	ClpItem found = parameters[iParam];
	ClpParameterType type = found.type();
	int valid;
	numberGoodCommands++;
	if (type==GENERALQUERY) {
	  std::cout<<"In argument list keywords have leading - "
	    ", -stdin or just - switches to stdin"<<std::endl;
	  std::cout<<"One command per line (and no -)"<<std::endl;
	  std::cout<<"abcd? gives list of possibilities, if only one + explanation"<<std::endl;
	  std::cout<<"abcd?? adds explanation, if only one fuller help(LATER)"<<std::endl;
	  std::cout<<"abcd without value (where expected) gives current value"<<std::endl;
	  std::cout<<"abcd value or abcd = value sets value"<<std::endl;
	  std::cout<<"Commands are:"<<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam+=4 ) {
	    int i;
	    for (i=iParam;i<min(numberParameters,iParam+4);i++) 
	      std::cout<<parameters[i].matchName()<<"  ";
	    std::cout<<std::endl;
	  }
	} else if (type<101) {
	  // get next field as double
	  double value = getDoubleField(argc,argv,&valid);
	  if (!valid) {
	    parameters[iParam].setDoubleParameter(models+iModel,value);
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].doubleParameter(models+iModel)<<std::endl;
	  }
	} else if (type<201) {
	  // get next field as int
	  int value = getIntField(argc,argv,&valid);
	  if (!valid) {
	    if (parameters[iParam].type()==PRESOLVEPASS)
	      preSolve = value;
	    else if (parameters[iParam].type()==IDIOT)
	      doIdiot = value;
	    else
	      parameters[iParam].setIntParameter(models+iModel,value);
	  } else if (valid==1) {
	    abort();
	  } else {
	    std::cout<<parameters[iParam].name()<<" has value "<<
	      parameters[iParam].intParameter(models+iModel)<<std::endl;
	  }
	} else if (type<301) {
	  // one of several strings
	  std::string value = getString(argc,argv);
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
	      else
		models[iModel].setOptimizationDirection(-11);
	      break;
	    case DUALPIVOT:
	      if (action==0) {
		ClpDualRowSteepest steep;
		models[iModel].setDualRowPivotAlgorithm(steep);
	      } else {
		ClpDualRowDantzig dantzig;
		models[iModel].setDualRowPivotAlgorithm(dantzig);
	      }
	      break;
	    case PRIMALPIVOT:
	      if (action==0) {
		ClpPrimalColumnSteepest steep(1);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else if (action==1) {
		ClpPrimalColumnSteepest steep(0);
		models[iModel].setPrimalColumnPivotAlgorithm(steep);
	      } else {
		ClpPrimalColumnDantzig dantzig;
		models[iModel].setPrimalColumnPivotAlgorithm(dantzig);
	      }
	      break;
	    case SCALING:
	      models[iModel].scaling(1-action);
	      break;
	    case SPARSEFACTOR:
	      models[iModel].setSparseFactorization(1-action);
	      break;
	    case ERRORSALLOWED:
	      allowImportErrors = action;
	      break;
	    case KEEPNAMES:
	      keepImportNames = 1-action;
	      break;
	    case PRESOLVE:
	      preSolve = action*5;
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
	    if (goodModels[iModel]) {
	      int saveMaxIterations = models[iModel].maximumIterations();
	      ClpSimplex * model2 = models+iModel;
#ifdef USE_PRESOLVE
	      Presolve pinfo;
	      if (preSolve) {
		model2 = pinfo.presolvedModel(models[iModel],1.0e-8,false,preSolve);
		model2->checkSolution();
#ifdef CLP_DEBUG
		printf("%g %g (%d) %g (%d)\n"
		       ,model2->objectiveValue()
		       ,model2->sumDualInfeasibilities()
		       ,model2->numberDualInfeasibilities()
		       ,model2->sumPrimalInfeasibilities()
		       ,model2->numberPrimalInfeasibilities());
#endif
		if (type==DUALSIMPLEX) {
		  int numberInfeasibilities = model2->tightenPrimalBounds();
		  if (numberInfeasibilities)
		    std::cout<<"** Analysis indicates model infeasible"
			     <<std::endl;
		}
	      }
#endif
#ifdef READLINE     
	      currentModel = model2;
#endif
	      if (type==DUALSIMPLEX) {
		model2->dual();
	      } else {
#ifdef CLP_IDIOT
		if (doIdiot) {
		  Idiot info(*model2);
		  info.crash(doIdiot);
		}
#endif
		model2->primal(1);
	      }
#ifdef USE_PRESOLVE
	      if (preSolve) {
		pinfo.postsolve(true);
		
		delete model2;
		printf("Resolving from postsolved model\n");
		
#ifdef READLINE     
		currentModel = models+iModel;
#endif
		models[iModel].primal(1);
#ifdef CLP_DEBUG_not
		models[iModel].checkSolution();
		printf("%g dual %g(%d) Primal %g(%d)\n",
		       models[iModel].objectiveValue(),
		       models[iModel].sumDualInfeasibilities(),
		       models[iModel].numberDualInfeasibilities(),
		       models[iModel].sumPrimalInfeasibilities(),
		       models[iModel].numberPrimalInfeasibilities());
		{
		  Presolve pinfoA;
		  model2 = pinfoA.presolvedModel(models[iModel],1.0e-8);
		  
		  printf("Resolving from presolved optimal solution\n");
		  model2->primal(1);
		  
		  delete model2;
		}
#endif
	      }
#endif
	      models[iModel].setMaximumIterations(saveMaxIterations);
	      time2 = cpuTime();
	      totalTime += time2-time1;
	      std::cout<<"Result "<<models[iModel].status()<<
		" - "<<models[iModel].objectiveValue()<<
		" iterations "<<models[iModel].numberIterations()<<
		" took "<<time2-time1<<" seconds - total "<<totalTime<<std::endl;
	      if (models[iModel].status())
		std::cerr<<"Non zero status "<<models[iModel].status()<<
		  std::endl;
	      time1=time2;
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
	  case BAB:
#if 0
	    if (goodModels[iModel]) {
	      int saveMaxIterations = models[iModel].maximumIterations();
#ifdef READLINE     
	      currentModel = models+iModel;
#endif
	      {
		// get OsiClp stuff 
		OsiClpSolverInterface m(models+iModel);
		m.getModelPtr()->messageHandler()->setLogLevel(0);
		m.branchAndBound();
		m.resolve();
		std::cout<<"Optimal solution "<<m.getObjValue()<<std::endl;
		m.releaseClp();
	      }
	      models[iModel].setMaximumIterations(saveMaxIterations);
	      time2 = cpuTime();
	      totalTime += time2-time1;
	      std::cout<<"Result "<<models[iModel].status()<<
		" - "<<models[iModel].objectiveValue()<<
		" iterations "<<models[iModel].numberIterations()<<
		" took "<<time2-time1<<" seconds - total "<<totalTime<<std::endl;
	      if (models[iModel].status())
		std::cerr<<"Non zero status "<<models[iModel].status()<<
		  std::endl;
	      time1=time2;
	    } else {
	      std::cout<<"** Current model not valid"<<std::endl;
	    }
#endif
	    break;
	  case IMPORT:
	    {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field=="-") {
		// stdin
		canOpen=true;
		fileName = "-";
	      } else {
		if (field[0]=='/'||field[0]=='~')
		  fileName = field;
		else
		  fileName = directory+field;
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
		int status =models[iModel].readMps(fileName.c_str(),
						   keepImportNames,
						   allowImportErrors);
		if (!status||(status>0&&allowImportErrors)) {
		  goodModels[iModel]=true;
		  // sets to all slack (not necessary?)
		  models[iModel].createStatus();
		  time2 = cpuTime();
		  totalTime += time2-time1;
		  time1=time2;
		} else {
		  // errors
		  std::cout<<"There were "<<status<<
		    " errors on input"<<std::endl;
		}
	      }
	    }
	    break;
	  case EXPORT:
	    {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='~')
		fileName = field;
	      else
		fileName = directory+field;
	      FILE *fp=fopen(fileName.c_str(),"w");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		// Convert names
		int iRow;
		int numberRows=models[iModel].numberRows();

		char ** rowNames = new char * [numberRows];
		for (iRow=0;iRow<numberRows;iRow++) {
		  rowNames[iRow] = 
		    strdup(models[iModel].rowName(iRow).c_str());
		}
		int iColumn;
		int numberColumns=models[iModel].numberColumns();

		char ** columnNames = new char * [numberColumns];
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  columnNames[iColumn] = 
		    strdup(models[iModel].columnName(iColumn).c_str());
		}

		ClpSimplex& m = models[iModel];
		CoinMpsIO writer;
		writer.setMpsData(*m.matrix(), CLP_INFINITY,
				  m.getColLower(), m.getColUpper(),
				  m.getObjCoefficients(),
				  (const char*) 0 /*integrality*/,
				  m.getRowLower(), m.getRowUpper(),
				  columnNames, rowNames);
		writer.writeMps(fileName.c_str());
		for (iRow=0;iRow<numberRows;iRow++) {
		  free(rowNames[iRow]);
		}
		delete [] rowNames;
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  free(columnNames[iColumn]);
		}
		delete [] columnNames;
		time2 = cpuTime();
		totalTime += time2-time1;
		time1=time2;
	      }
	    }
	    break;
	  case SAVE:
	    {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='~')
		fileName = field;
	      else
		fileName = directory+field;
	      FILE *fp=fopen(fileName.c_str(),"wb");
	      if (fp) {
		// can open - lets go for it
		fclose(fp);
		canOpen=true;
	      } else {
		std::cout<<"Unable to open file "<<fileName<<std::endl;
	      }
	      if (canOpen) {
		int status =models[iModel].saveModel(fileName.c_str());
		if (!status) {
		  goodModels[iModel]=true;
		  time2 = cpuTime();
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
	      field = getString(argc,argv);
	      std::string fileName;
	      bool canOpen=false;
	      if (field[0]=='/'||field[0]=='~')
		fileName = field;
	      else
		fileName = directory+field;
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
		  time2 = cpuTime();
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
	  case DIRECTORY:
	    directory = getString(argc,argv);
	    break;
	  case STDIN:
	    read_mode=-1;
	    break;
	  case NETLIB_DUAL:
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
	      if (type==NETLIB_DUAL)
		std::cerr<<"Doing netlib with dual agorithm"<<std::endl;
	      else
		std::cerr<<"Doing netlib with primal agorithm"<<std::endl;
	      mainTest(nFields,fields,(type==NETLIB_DUAL),models[iModel],
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
	      double value = getDoubleField(argc,argv,&valid);
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
		    rowLower[iRow]=max(rowLower[iRow],-value);
		    rowUpper[iRow]=min(rowUpper[iRow],value);
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
		    columnLower[iColumn]=max(columnLower[iColumn],-value);
		    columnUpper[iColumn]=min(columnUpper[iColumn],value);
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
	  case VERSION:
	    std::cout<<"Coin LP version "<<CLPVERSION
		     <<", build "<<__DATE__<<std::endl;
	    break;
	  case SOLUTION:
	    if (goodModels[iModel]) {
	      // get next field
	      field = getString(argc,argv);
	      std::string fileName;
	      FILE *fp=NULL;
	      if (field=="-"||field=="EOL") {
		// stdout
		fp=stdout;
	      } else {
		if (field[0]=='/'||field[0]=='~')
		  fileName = field;
		else
		  fileName = directory+field;
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
		char format[6];
		sprintf(format,"%%-%ds",max(lengthName,8));
		for (iRow=0;iRow<numberRows;iRow++) {
		  fprintf(fp,"%7d ",iRow);
		  if (lengthName)
		    fprintf(fp,format,rowNames[iRow].c_str());
		  fprintf(fp,"%15.8g        %15.8g\n",primalRowSolution[iRow],
			  dualRowSolution[iRow]);
		}
		int iColumn;
		int numberColumns=models[iModel].numberColumns();
		double * dualColumnSolution = 
		  models[iModel].dualColumnSolution();
		double * primalColumnSolution = 
		  models[iModel].primalColumnSolution();
		for (iColumn=0;iColumn<numberColumns;iColumn++) {
		  fprintf(fp,"%7d ",iColumn);
		  if (lengthName)
		    fprintf(fp,format,columnNames[iColumn].c_str());
		  fprintf(fp,"%15.8g        %15.8g\n",
			  primalColumnSolution[iColumn],
			  dualColumnSolution[iColumn]);
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
	  std::cout<<"Short match for "<<field<<" possible completion:"
		   <<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam++ ) {
	    int match = parameters[iParam].matches(field);
	    if (match) 
	      std::cout<<parameters[iParam].matchName()<<std::endl;
	  }
	} else if (numberQuery) {
	  std::cout<<"Short match for "<<field<<" completion:"
		   <<std::endl;
	  for ( iParam=0; iParam<numberParameters; iParam++ ) {
	    int match = parameters[iParam].matches(field);
	    if (match) {
	      std::cout<<parameters[iParam].matchName()<<" : ";
	      std::cout<<parameters[iParam].shortHelp()<<std::endl;
	    }
	  }
	}
      } else {
	if (!numberQuery) 
	  std::cout<<"Multiple matches for "<<field<<" - possible completions:"
		   <<std::endl;
	else
	  std::cout<<"Completions of "<<field<<":"<<std::endl;
	for ( iParam=0; iParam<numberParameters; iParam++ ) {
	  int match = parameters[iParam].matches(field);
	  if (match) {
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
