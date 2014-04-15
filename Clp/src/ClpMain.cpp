/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"

#include "AbcCommon.hpp"
#include "ClpSimplex.hpp"
#ifdef ABC_INHERIT
#include "AbcSimplex.hpp"
#endif
#ifndef ABC_INHERIT
void ClpMain0(ClpSimplex * models);
int ClpMain1(int argc, const char *argv[],ClpSimplex * model);
#else
void ClpMain0(AbcSimplex * models);
int ClpMain1(int argc, const char *argv[],AbcSimplex * model);
#endif
//#define CILK_TEST
#ifdef CILK_TEST
static void cilkTest();
#endif
int
#if defined(_MSC_VER)
__cdecl
#endif // _MSC_VER
main (int argc, const char *argv[])
{
#ifdef CILK_TEST
  cilkTest();
#endif
#ifndef ABC_INHERIT
  ClpSimplex * models = new ClpSimplex[1];
#else
  AbcSimplex * models = new AbcSimplex[1];
#endif
  std::cout << "Coin LP version " << CLP_VERSION
	    << ", build " << __DATE__ << std::endl;
  // Print command line
  if (argc > 1) {
    printf("command line - ");
    for (int i = 0; i < argc; i++)
      printf("%s ", argv[i]);
    printf("\n");
  }
  ClpMain0(models);
  int returnCode = ClpMain1(argc, argv,models);
  delete [] models;
  return returnCode;
}
/*
  Version 1.00.00 October 13 2004.
  1.00.01 October 18.  Added basis handling helped/prodded by Thorsten Koch.
  Also modifications to make faster with sbb (I hope I haven't broken anything).
  1.00.02 March 21 2005.  Redid ClpNonLinearCost to save memory also redid
  createRim to try and improve cache characteristics.
  1.00.03 April 8 2005.  Added Volume algorithm as crash and made code more
  robust on testing.  Also added "either" and "tune" algorithm.
  1.01.01 April 12 2005.  Decided to go to different numbering.  Backups will
  be last 2 digits while middle 2 are for improvements.  Still take a long
  time to get to 2.00.01
  1.01.02 May 4 2005.  Will be putting in many changes - so saving stable version
  1.02.01 May 6 2005.  Lots of changes to try and make faster and more stable in
  branch and cut.
  1.02.02 May 19 2005.  Stuff for strong branching and some improvements to simplex
  1.03.01 May 24 2006.  Lots done but I can't remember what!
  1.03.03 June 13 2006.  For clean up after dual perturbation
  1.04.01 June 26 2007.  Lots of changes but I got lazy
  1.05.00 June 27 2007.  This is trunk so when gets to stable will be 1.5
  1.11.00 November 5 2009 (Guy Fawkes) - OSL factorization and better ordering
 */
#ifdef CILK_TEST
// -*- C++ -*-

/*
 * cilk-for.cilk
 *
 * Copyright (c) 2007-2008 Cilk Arts, Inc.  55 Cambridge Street,
 * Burlington, MA 01803.  Patents pending.  All rights reserved. You may
 * freely use the sample code to guide development of your own works,
 * provided that you reproduce this notice in any works you make that
 * use the sample code.  This sample code is provided "AS IS" without
 * warranty of any kind, either express or implied, including but not
 * limited to any implied warranty of non-infringement, merchantability
 * or fitness for a particular purpose.  In no event shall Cilk Arts,
 * Inc. be liable for any direct, indirect, special, or consequential
 * damages, or any other damages whatsoever, for any use of or reliance
 * on this sample code, including, without limitation, any lost
 * opportunity, lost profits, business interruption, loss of programs or
 * data, even if expressly advised of or otherwise aware of the
 * possibility of such damages, whether in an action of contract,
 * negligence, tort, or otherwise.
 *
 * This file demonstrates a Cilk++ for loop
 */

#include <cilk/cilk.h>
//#include <cilk/cilkview.h>
#include <cilk/reducer_max.h>
#include <cstdlib>
#include <iostream>

// cilk_for granularity.
#define CILK_FOR_GRAINSIZE 128

double dowork(double i)
{
    // Waste time:
    int j;
    double k = i;
    for (j = 0; j < 50000; ++j) {
        k += k / ((j + 1) * (k + 1));
    }

    return k;
}
static void doSomeWork(double * a,int low, int high)
{
  if (high-low>300) {
    int mid=(high+low)>>1;
    cilk_spawn doSomeWork(a,low,mid);
    doSomeWork(a,mid,high);
    cilk_sync;
  } else {
    for(int i = low; i < high; ++i) {
      a[i] = dowork(a[i]);
    }
  }
}

void cilkTest()
{
    unsigned int n = 10000;
    //cilk::cilkview cv;


    double* a = new double[n];

    for(unsigned int i = 0; i < n; i++) {
        // Populate A 
        a[i] = (double) ((i * i) % 1024 + 512) / 512;
    }

    std::cout << "Iterating over " << n << " integers" << std::endl;

    //cv.start();
#if 1
    //#pragma cilk_grainsize=CILK_FOR_GRAINSIZE
    cilk_for(unsigned int i = 0; i < n; ++i) {
        a[i] = dowork(a[i]);
    }
#else
    doSomeWork(a,0,n);
#endif
    int * which =new int[n];
    unsigned int n2=n>>1;
    for (int i=0;i<n2;i++) 
      which[i]=n-2*i;
    cilk::reducer_max_index<int,double> maximumIndex(-1,0.0);
    cilk_for(unsigned int i = 0; i < n2; ++i) {
      int iWhich=which[i];
      maximumIndex.calc_max(iWhich,a[iWhich]);
    }
    int bestIndex=maximumIndex.get_index();
    int bestIndex2=-1;
    double largest=0.0;
    cilk_for(unsigned int i = 0; i < n2; ++i) {
      int iWhich=which[i];
      if (a[iWhich]>largest) {
	bestIndex2=iWhich;
	largest=a[iWhich];
      }
    }
    assert (bestIndex==bestIndex2);
    //cv.stop();
    //cv.dump("cilk-for-results", false);

    //std::cout << cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;

    exit(0);
}
#endif
