
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpConfig.h"

#ifndef ClpParam_H
#define ClpParam_H
/**
   This has parameter handling stuff for Clp.

 */
class OsiSolverInterface;
class ClpSimplex;
/*! \brief Parameter codes

  Parameter type ranges are allocated as follows
  <ul>
    <li>   1 -- 100  double parameters
    <li> 101 -- 200  integer parameters
    <li> 201 -- 300  Clp string parameters
    <li> 301 -- 400  Cbc string parameters
    <li> 401 -- 500  (mostly) Clp actions
    <li> 501 -- 600  (mostly) Cbc actions
  </ul>

  `Actions' do not necessarily invoke an immediate action; it's just that they
  don't fit neatly into the parameters array.

  This coding scheme is in flux.
*/

enum ClpParameterType

{ //Begin used in Cbc
  CLP_PARAM_GENERALQUERY = -100,
  CLP_PARAM_FULLGENERALQUERY,
  CLP_PARAM_STR_ALLCOMMANDS,
  //End used in Cbc
  
  CLP_PARAM_DBL_PRIMALTOLERANCE = 1,
  CLP_PARAM_DBL_DUALTOLERANCE,
  CLP_PARAM_DBL_TIMELIMIT,
  CLP_PARAM_DBL_DUALBOUND,
  CLP_PARAM_DBL_PRIMALWEIGHT,
  CLP_PARAM_DBL_OBJSCALE,
  CLP_PARAM_DBL_RHSSCALE,
  CLP_PARAM_DBL_ZEROTOLERANCE,
  CLP_PARAM_DBL_PSI,
  CLP_PARAM_DBL_PROGRESS,
  CLP_PARAM_DBL_PRESOLVETOLERANCE,
  CLP_PARAM_DBL_OBJSCALE2,

  //Begin Used in Cbc
  CLP_PARAM_INT_LOGLEVEL = 101,
  CLP_PARAM_INT_OUTPUTFORMAT,
  CLP_PARAM_INT_PRINTOPTIONS,
  CLP_PARAM_INT_VERBOSE,
  //End used in Cbc
  CLP_PARAM_INT_PRESOLVEPASS,
  CLP_PARAM_INT_PRESOLVEOPTIONS,
  CLP_PARAM_INT_MAXFACTOR,
  CLP_PARAM_INT_PERTVALUE,
  CLP_PARAM_INT_MAXITERATION,
  CLP_PARAM_INT_IDIOT,
  CLP_PARAM_INT_SPRINT,
  CLP_PARAM_INT_SLPVALUE,
  CLP_PARAM_INT_SPECIALOPTIONS,
  CLP_PARAM_INT_SUBSTITUTION,
  CLP_PARAM_INT_DUALIZE,
  CLP_PARAM_INT_CPP,
  CLP_PARAM_INT_RANDOMSEED,
  CLP_PARAM_INT_MORESPECIALOPTIONS,
  CLP_PARAM_INT_DECOMPOSE_BLOCKS,
  CLP_PARAM_INT_VECTOR_MODE,

  //Begin used in Cbc
  CLP_PARAM_INT_THREADS = 151,
  //End used in Cbc
  CLP_PARAM_INT_DENSE,
  CLP_PARAM_INT_SMALLFACT,

  //Begin used in Cbc
  CLP_PARAM_STR_DIRECTION = 201,
  CLP_PARAM_STR_ERRORSALLOWED,
  CLP_PARAM_STR_INTPRINT,
  //End used in Cbc
  CLP_PARAM_STR_PRESOLVE,
  CLP_PARAM_STR_DUALPIVOT,
  CLP_PARAM_STR_SCALING,
  CLP_PARAM_STR_KEEPNAMES,
  CLP_PARAM_STR_SPARSEFACTOR,
  CLP_PARAM_STR_PRIMALPIVOT,
  CLP_PARAM_STR_CRASH,
  CLP_PARAM_STR_BIASLU,
  CLP_PARAM_STR_PERTURBATION,
  CLP_PARAM_STR_MESSAGES,
  CLP_PARAM_STR_AUTOSCALE,
  CLP_PARAM_STR_CHOLESKY,
  CLP_PARAM_STR_KKT,
  CLP_PARAM_STR_BARRIERSCALE,
  CLP_PARAM_STR_GAMMA,
  CLP_PARAM_STR_CROSSOVER,
  CLP_PARAM_STR_PFI,
  CLP_PARAM_STR_VECTOR,
  CLP_PARAM_STR_FACTORIZATION,
  CLP_PARAM_STR_ABCWANTED,
  CLP_PARAM_STR_BUFFER_MODE,

  //Begin used in Cbc
  CLP_PARAM_ACTION_DIRECTORY = 401,
  CLP_PARAM_ACTION_DIRSAMPLE,
  CLP_PARAM_ACTION_DIRNETLIB,
  CLP_PARAM_ACTION_DIRMIPLIB,
  CLP_PARAM_ACTION_IMPORT,
  CLP_PARAM_ACTION_EXPORT,
  CLP_PARAM_ACTION_RESTORE,
  CLP_PARAM_ACTION_SAVE,
  CLP_PARAM_ACTION_MAXIMIZE,
  CLP_PARAM_ACTION_MINIMIZE,
  CLP_PARAM_ACTION_EXIT,
  CLP_PARAM_ACTION_STDIN,
  CLP_PARAM_ACTION_SOLUTION,
  CLP_PARAM_ACTION_SAVESOL,
  CLP_PARAM_ACTION_UNITTEST,
  CLP_PARAM_ACTION_HELP,
  CLP_PARAM_ACTION_REVERSE,
  CLP_PARAM_ACTION_STATISTICS,
  CLP_PARAM_ACTION_PRINTMASK,
  CLP_PARAM_ACTION_ENVIRONMENT,
  CLP_PARAM_ACTION_GMPL_SOLUTION,
  CLP_PARAM_ACTION_DUMMY,
  //End used in Cbc
  
  CLP_PARAM_ACTION_DUALSIMPLEX,
  CLP_PARAM_ACTION_PRIMALSIMPLEX,
  CLP_PARAM_ACTION_EITHERSIMPLEX,
  CLP_PARAM_ACTION_NETLIB_EITHER,
  CLP_PARAM_ACTION_NETLIB_DUAL,
  CLP_PARAM_ACTION_NETLIB_PRIMAL,
  CLP_PARAM_ACTION_TIGHTEN,
  CLP_PARAM_ACTION_FAKEBOUND,
  CLP_PARAM_ACTION_PLUSMINUS,
  CLP_PARAM_ACTION_NETWORK,
  CLP_PARAM_ACTION_ALLSLACK,
  CLP_PARAM_ACTION_BARRIER,
  CLP_PARAM_ACTION_NETLIB_BARRIER,
  CLP_PARAM_ACTION_NETLIB_TUNE,
  CLP_PARAM_ACTION_REALLY_SCALE,
  CLP_PARAM_ACTION_BASISIN,
  CLP_PARAM_ACTION_BASISOUT,
  CLP_PARAM_ACTION_SOLVECONTINUOUS,
  CLP_PARAM_ACTION_CLEARCUTS,
  CLP_PARAM_ACTION_VERSION,
  CLP_PARAM_ACTION_OUTDUPROWS,
  CLP_PARAM_ACTION_USERCLP,
  CLP_PARAM_ACTION_MODELIN,
  CLP_PARAM_ACTION_PARAMETRICS,
  CLP_PARAM_ACTION_RESTORESOL,
  CLP_PARAM_ACTION_GUESS,

  //Begin used in Cbc
  CLP_PARAM_ACTION_BAB = 501,
  //End used in Cbc

  CLP_PARAM_NOTUSED_INVALID = 1000
};
#include <vector>
#include <string>

/// Very simple class for setting parameters

class CLPLIB_EXPORT ClpParam {
public:
  /**@name Constructor and destructor */
  //@{
  /// Constructors
  ClpParam();
  ClpParam(std::string name, std::string help,
    double lower, double upper, ClpParameterType type, int display = 2);
  ClpParam(std::string name, std::string help,
    int lower, int upper, ClpParameterType type, int display = 2);
  // Other strings will be added by insert
  ClpParam(std::string name, std::string help, std::string firstValue,
    ClpParameterType type, int whereUsed = 7, int display = 2);
  // Action
  ClpParam(std::string name, std::string help,
    ClpParameterType type, int whereUsed = 7, int display = 2);
  /// Copy constructor.
  ClpParam(const ClpParam &);
  /// Assignment operator. This copies the data
  ClpParam &operator=(const ClpParam &rhs);
  /// Destructor
  ~ClpParam();
  //@}

  /**@name stuff */
  //@{
  /// Insert string (only valid for keywords)
  void append(std::string keyWord);
  /// Returns name
  inline std::string name() const
  {
    return name_;
  }
  /// Returns short help
  inline std::string shortHelp() const
  {
    return shortHelp_;
  }
  /// Returns long help
  inline std::string longHelp() const
  {
    return longHelp_;
  }
  /// Returns set of valid strings
  inline const std::vector<std::string>& definedKeywords() const
  {
    return definedKeyWords_;
  }
  /// Returns the lower bound for a double-valued parameter
  inline double lowerDoubleValue() const
  {
     return lowerDoubleValue_;
  }
  /// Returns the upper bound for a double-valued parameter
  inline double upperDoubleValue() const
  {
     return upperDoubleValue_;
  }
  /// Returns the lower bound for an int-valued parameter
  inline int lowerIntValue() const
  {
     return lowerIntValue_;
  }
  /// Returns the upper bound for an int-valued parameter
  inline int upperIntValue() const
  {
     return upperIntValue_;
  }
  /// Gets a double parameter
  double doubleParameter(ClpSimplex *model) const;
  /// Sets a double parameter (nonzero code if error)
  int setDoubleParameter(ClpSimplex *model, double value,
                         bool doPrinting = false);
  /// Sets double parameter and returns printable string and error code
  std::string setDoubleParameterWithMessage(ClpSimplex *model, double value,
                                            int &returnCode);
  /// Gets double value
  inline double doubleValue() const
  {
    return doubleValue_;
  }
  /// Sets double value
  void setDoubleValue(double value);
  /// Checks a double parameter (nonzero code if error)
  int checkDoubleParameter(double value) const;
  /// Gets a int parameter
  int intParameter(ClpSimplex *model) const;
  /// Sets a int parameter (nonzero code if error)
  int setIntParameter(ClpSimplex *model, int value,
                      bool doPrinting = false);
  /// Sets int parameter and returns printable string and error code
  std::string setIntParameterWithMessage(ClpSimplex *model, int value,
                                         int &returnCode);
  // Gets int value
  inline int intValue() const
  {
    return intValue_;
  }
  /// Sets int value
  void setIntValue(int value);
  /// Gets string value
  inline std::string stringValue() const
  {
    return stringValue_;
  }
  /// Sets string value
  void setStringValue(std::string value);
  /// Returns name which could match
  std::string matchName() const;
  /// Returns length of name for ptinting
  int lengthMatchName() const;
  /// Returns 1 if matches minimum, 2 if matches less, 0 if not matched
  int matches(std::string input) const;
  /// Returns parameter option which matches (-1 if none)
  int parameterOption(std::string check) const;
  /// Prints parameter options
  void printOptions() const;
  /// Returns current parameter option
  inline std::string currentOption() const
  {
    return definedKeyWords_[currentKeyWord_];
  }
  /// Sets current parameter option
  void setCurrentOption(int value, bool printIt = false);
  /// Sets current parameter option using string
  void setCurrentOption(const std::string value);
  /// Sets current parameter option and returns printable string
  std::string setCurrentOptionWithMessage(int value);
  /// Sets current parameter option using string with message
  std::string setCurrentOptionWithMessage(const std::string value);
  /// Returns current parameter option position
  int currentOptionAsInteger() const;
  /** Returns current parameter option position
	 but if fake keyword returns a fake value and sets
	 fakeInteger to true value.  If not fake then fakeInteger is -COIN_INT_MAX
      */
  int currentOptionAsInteger(int &fakeInteger) const;
  /// type
  inline ClpParameterType type() const
  {
    return type_;
  }
  /// whether to display
  inline int displayThis() const
  {
    return display_;
  }
  /// Set Long help
  inline void setLonghelp(const std::string help)
  {
    longHelp_ = help;
  }
  /// Print Long help
  void printLongHelp() const;
  /// Print action and string
  void printString() const;
  /** 7 if used everywhere,
         1 - used by clp
         2 - used by cbc
         4 - used by ampl
     */
  inline int whereUsed() const
  {
    return whereUsed_;
  }
  /// Gets value of fake keyword
  inline int fakeKeyWord() const
  {
    return fakeKeyWord_;
  }
  /// Sets value of fake keyword
  inline void setFakeKeyWord(int value, int fakeValue)
  {
    fakeKeyWord_ = value;
    fakeValue_ = fakeValue;
  }
  /// Sets value of fake keyword to current size of keywords
  void setFakeKeyWord(int fakeValue);

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
  std::vector< std::string > definedKeyWords_;
  /// Name
  std::string name_;
  /// Short help
  std::string shortHelp_;
  /// Long help
  std::string longHelp_;
  /// Action
  ClpParameterType action_;
  /// Current keyWord (if a keyword parameter)
  int currentKeyWord_;
  /// Display on ?
  int display_;
  /// Integer parameter - current value
  int intValue_;
  /// Double parameter - current value
  double doubleValue_;
  /// String parameter - current value
  std::string stringValue_;
  /** 7 if used everywhere,
         1 - used by clp
         2 - used by cbc
         4 - used by ampl
     */
  int whereUsed_;
  /** If >=0 then integers allowed as a fake keyword
	 So minusnnnn would got to -nnnn in currentKeyword_
	 and plusnnnn would go to fakeKeyword_+nnnn
     */
  int fakeKeyWord_;
  /// Return this as main value if an integer
  int fakeValue_;
  //@}
};

CLPLIB_EXPORT
void setClpPrinting(bool yesNo);

#define CBCMAXPARAMETERS 250
/*
  Subroutine to establish the clp parameter array. See the description of
  class ClpParam for details. 
*/
CLPLIB_EXPORT
void establishClpParams(std::vector< ClpParam > &params);
// Given a parameter type - returns its number in list
CLPLIB_EXPORT
int whichClpParam(const ClpParameterType &name,
  const std::vector< ClpParam > &parameters);
// Dump/restore a solution to file
CLPLIB_EXPORT
void saveSolution(const ClpSimplex *lpSolver, std::string fileName);
CLPLIB_EXPORT
void restoreSolution(ClpSimplex *lpSolver, std::string fileName, int mode);

#endif /* ClpParam_H */

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
