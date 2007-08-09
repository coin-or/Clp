/*
 * Include file for the configuration of Clp.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __CLPCONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_clp.h"

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#else /* HAVE_CONFIG_H */

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/***************************************************************************/
/*             HERE DEFINE THE CONFIGURATION SPECIFIC MACROS               */
/***************************************************************************/

/* This needs to be updated "by hand"!!!! */
#define CLPVERSION "1.4trunk"

/* Define to the debug sanity check level (0 is no test) */
#define COIN_CLP_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_CLP_VERBOSITY 0

/* Define to 1 if the Clp package is used */
#define COIN_HAS_CLP 1

/* Define to 1 if the CoinUtils package is used */
#define COIN_HAS_COINUTILS 1

/* Define to 1 if readline is available */
/* #undef COIN_HAS_READLINE */

#endif /* HAVE_CONFIG_H */

#endif /*__HAVE_COINUTILS_CONFIG_H__*/
