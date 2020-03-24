
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_clp_default.h */
#ifndef CLPLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define CLPLIB_EXPORT __declspec(dllexport)
#else
#define CLPLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_clp_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to 1 if the CoinUtils package is used */
#define CLP_HAS_COINUTILS 1

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
