/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#ifndef ClpParam_H
#define ClpParam_H

#include "ClpConfig.h"

#include "CoinParam.hpp"

class ClpParameters;

/* \file ClpParam.hpp
   \brief Declarations for parameters that control Clp.
*/

/*
 */

/*! \class ClpParam
    \brief Class for Clp control parameters

    Adds parameter type codes and push/pull functions to the generic parameter
    object.
*/

class CLPLIB_EXPORT ClpParam : public CoinParam {

public:

   /*! \name Enumeration types used to index parameters */
   //@{
   
   /*! \brief Codes to specify paramters */
   
   enum ClpParamCode {
      FIRSTPARAM = 0,
      
      // Help and Information Parameters
      FIRSTHELPPARAM,
      GENERALQUERY,
      FULLGENERALQUERY,
      LASTHELPPARAM,
      
      // ActionParameters
      FIRSTACTIONPARAM,
      // Begin used in Cbc
      MAXIMIZE,
      MINIMIZE,
      EXIT,
      END,
      QUIT,
      STOP,
      STDIN,
      UNITTEST,
      HELP,
      REVERSE,
      STATISTICS,
      ENVIRONMENT,
      GMPL_SOLUTION,
      DUMMY,
      BAB,
      // End used in Cbc
      DUALSIMPLEX,
      PRIMALSIMPLEX,
      EITHERSIMPLEX,
      NETLIB_EITHER,
      NETLIB_DUAL,
      NETLIB_PRIMAL,
      TIGHTEN,
      PLUSMINUS,
      NETWORK,
      ALLSLACK,
      BARRIER,
      NETLIB_BARRIER,
      NETLIB_TUNE,
      REALLY_SCALE,
      SOLVE,
      SOLVECONTINUOUS,
      CLEARCUTS,
      VERSION,
      OUTDUPROWS,
      USERCLP,
      MODELIN,
      PARAMETRICS,
      GUESS,
      LASTACTIONPARAM,
      
      FIRSTSTRINGPARAM,
      // Begin used in Cbc
      DIRECTORY,
      DIRSAMPLE,
      DIRNETLIB,
      DIRMIPLIB,
      IMPORT,
      EXPORT,
      PRINTMASK,
      RESTORE,
      SAVE,
      SAVESOL,
      SOLUTION,
      // End used in Cbc
      BASISIN,
      BASISOUT,
      LASTSTRINGPARAM,
      
      // On/Off Parameters
      FIRSTBOOLPARAM,
      AUTOSCALE,
      BUFFER_MODE,
      // Begin used in Cbc
      ERRORSALLOWED,
      // End used in Cbc
      KEEPNAMES,
      KKT,
      MESSAGES,
      PERTURBATION,
      PFI,
      SPARSEFACTOR,
      LASTBOOLPARAM,
      
      // Keyword Parameters
      FIRSTKWDPARAM,
      ABCWANTED,
      BARRIERSCALE,
      BIASLU,
      CHOLESKY,
      COMMANDPRINTLEVEL,
      CRASH,
      CROSSOVER,
      DIRECTION,
      DUALPIVOT,
      FACTORIZATION,
      GAMMA,
      INTPRINT,
      PRESOLVE,
      PRIMALPIVOT,
      SCALING,
      VECTOR,
      LASTKWDPARAM,
      
      // Integer Parameters
      FIRSTINTPARAM,
      // Begin Used in Cbc
      LOGLEVEL,
      OUTPUTFORMAT,
      PRINTOPTIONS,
      VERBOSE,
      THREADS,
      // End used in Cbc
      PRESOLVEPASS,
      PRESOLVEOPTIONS,
      MAXFACTOR,
      PERTVALUE,
      MAXITERATION,
      IDIOT,
      SPRINT,
      SLPVALUE,
      SPECIALOPTIONS,
      SUBSTITUTION,
      DUALIZE,
      CPP,
      RANDOMSEED,
      MORESPECIALOPTIONS,
      DECOMPOSE_BLOCKS,
      VECTOR_MODE,
      DENSE,
      SMALLFACT,
      LASTINTPARAM,
      
      // Double Paramters
      FIRSTDBLPARAM,
      PRIMALTOLERANCE,
      DUALTOLERANCE,
      FAKEBOUND,
      TIMELIMIT,
      DUALBOUND,
      PRIMALWEIGHT,
      OBJSCALE,
      RHSSCALE,
      ZEROTOLERANCE,
      PSI,
      PROGRESS,
      PRESOLVETOLERANCE,
      OBJSCALE2,
      LASTDBLPARAM,
      
      INVALID,
      
      LASTPARAM
   };

  //@}

  /*! \name Constructors and Destructors

      Be careful how you specify parameters for the constructors! There's great
      potential for confusion.
    */
  //@{
  /*! \brief Default constructor */

  ClpParam();

  /*! \brief Constructor for a parameter with a double value

      The default value is 0.0. Be careful to clearly indicate that \p lower and
      \p upper are real (double) values to distinguish this constructor from the
      constructor for an integer parameter.
    */
  ClpParam(int code, std::string name, std::string help,
           double lower = -COIN_DBL_MAX, double upper = COIN_DBL_MAX,
           double defaultValue = 0.0, std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);
  
  /*! \brief Constructor for a parameter with an integer value
    
    The default value is 0.
  */
  ClpParam(int code, std::string name, std::string help,
           int lower = -COIN_INT_MAX, int upper = COIN_INT_MAX,
           int defaultValue = 0, std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);
  
  /*! \brief Constructor for a parameter with keyword values

      The string supplied as \p firstValue becomes the first keyword.
      Additional keywords can be added using appendKwd(). Keywords are numbered
      from zero. It's necessary to specify both the first keyword (\p
      firstValue) and the default keyword index (\p dflt) in order to
      distinguish this constructor from the string and action parameter
      constructors.
    */
  ClpParam(int code, std::string name, std::string help,
           std::string defaultKwd, int defaultMode,
           std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);
  
  /*! \brief Constructor for a string parameter

      The default string value must be specified explicitly to distinguish
      a string constructor from an action parameter constructor.
    */
  ClpParam(int code, std::string name, std::string help,
           std::string defaultValue, std::string longHelp = "",
           CoinDisplayPriority displayPriority = CoinParam::displayPriorityHigh);

  /*! \brief Constructor for an action parameter */
  // No defaults to resolve ambiguity
  ClpParam(int code, std::string name, std::string help,
           std::string longHelp, CoinDisplayPriority displayPriority);

  /*! \brief Copy constructor */
  ClpParam(const ClpParam &orig);

  /*! \brief Clone */
  ClpParam *clone();

  /*! \brief Assignment */
  ClpParam &operator=(const ClpParam &rhs);

  /*! \brief  Destructor */
  ~ClpParam();

  //@}

  /*! \name Methods to query and manipulate a parameter object */
  //@{

  /*! \brief Get the parameter code  */
  inline int paramCode() const { return (paramCode_); }

  /*! \brief Set the parameter code */
  inline void setParamCode(int code) { paramCode_ = code; }

  /*! \brief Get the enclosing ClpParameters object */
  inline ClpParameters *parameters() const { return (parameters_); }

  /*! \brief Set the enclosing ClpParameters object */
  inline void setParameters(ClpParameters *p) { parameters_ = p; }

  /*! \brief A hacky function to print some information about string parameters */
   std::string printString() const;

//@}

private:
  /*! \name Data */
  //@{

  /// Parameter code
  int paramCode_;

  /// Settings object
  ClpParameters *parameters_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
