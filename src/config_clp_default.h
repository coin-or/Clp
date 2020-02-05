
/***************************************************************************/
/*           HERE DEFINE THE PROJECT SPECIFIC PUBLIC MACROS                */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Version number of project */
#define CLP_VERSION "trunk"

/* Major Version number of project */
#define CLP_VERSION_MAJOR 9999

/* Minor Version number of project */
#define CLP_VERSION_MINOR 9999

/* Release Version number of project */
#define CLP_VERSION_RELEASE 9999

#ifndef CLP_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define CLP_EXPORT __declspec(dllimport)
#else
#define CLP_EXPORT
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
