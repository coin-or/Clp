/* Copyright (C) 2002, 2003 International Business Machines
   Corporation and others.  All Rights Reserved.*/
#ifndef ClpSimplexCDefine_H
#define ClpSimplexCDefine_H

/** This has #defines etc for the "C" interface to Clp.

*/

/* Plus infinity */
#ifndef COIN_DBL_MAX
#define COIN_DBL_MAX DBL_MAX
#endif

/* We need to allow for Microsoft */
#ifndef CLPLIBAPI

#if defined (CLPMSDLL)
#   define CLPLIBAPI __declspec(dllexport)
#   define CLPLINKAGE  __stdcall
#   define CLPLINKAGE_CB  __cdecl
#else
#   define CLPLIBAPI 
#   define CLPLINKAGE
#   define CLPLINKAGE_CB 
#endif

#endif
/** User does not need to see structure of model */
typedef void Clp_Simplex;
/** typedef for user call back.
 The cvec are constructed so don't need to be const*/
typedef  void (CLPLINKAGE_CB *clp_callback) (Clp_Simplex * model,int  msgno, int ndouble,
                            const double * dvec, int nint, const int * ivec,
                            int nchar, char ** cvec);

#if COIN_BIG_INDEX==0
typedef int CoinBigIndex;
#elif COIN_BIG_INDEX==1
typedef long CoinBigIndex;
#else
typedef long long CoinBigIndex;
#endif
#endif
