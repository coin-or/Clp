// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

#include <time.h>

#include "CoinMpsIO.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpParameters.hpp"

#include "Presolve.hpp"
#ifdef CLP_IDIOT
#include "Idiot.hpp"
#endif


#include <time.h>
#ifndef _MSC_VER
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

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

// Function Prototypes. Function definitions is in this file.
void testingMessage( const char * const msg );

//----------------------------------------------------------------
// unitTest [-mpsDir=V1] [-netlibDir=V2] [-netlib]
// 
// where:
//   -mpsDir: directory containing mps test files
//       Default value V1="../Mps/Sample"    
//   -netlibDir: directory containing netlib files
//       Default value V2="../Mps/Netlib"
//   -netlib
//       If specified, then netlib test set run
//
// All parameters are optional.
//----------------------------------------------------------------
int mainTest (int argc, const char *argv[],bool doDual,
	      ClpSimplex empty, bool doPresolve,int doIdiot)
{
  int i;

  // define valid parameter keywords
  std::set<std::string> definedKeyWords;
  definedKeyWords.insert("-mpsDir");
  definedKeyWords.insert("-netlibDir");
  definedKeyWords.insert("-netlib");

  // Create a map of parameter keys and associated data
  std::map<std::string,std::string> parms;
  for ( i=1; i<argc; i++ ) {
    std::string parm(argv[i]);
    std::string key,value;
    unsigned int  eqPos = parm.find('=');

    // Does parm contain and '='
    if ( eqPos==std::string::npos ) {
      //Parm does not contain '='
      key = parm;
    }
    else {
      key=parm.substr(0,eqPos);
      value=parm.substr(eqPos+1);
    }

    // Is specifed key valid?
    if ( definedKeyWords.find(key) == definedKeyWords.end() ) {
      // invalid key word.
      // Write help text
      std::cerr <<"Undefined parameter \"" <<key <<"\".\n";
      std::cerr <<"Correct usage: \n";
      std::cerr <<"  unitTest [-mpsDir=V1] [-netlibDir=V2] [-netlib]\n";
      std::cerr <<"  where:\n";
      std::cerr <<"    -mpsDir: directory containing mps test files\n";
      std::cerr <<"        Default value V1=\"../Mps/Sample\"\n";
      std::cerr <<"    -netlibDir: directory containing netlib files\n";
      std::cerr <<"        Default value V2=\"../Mps/Netlib\"\n";
      std::cerr <<"    -netlib\n";
      std::cerr <<"        If specified, then netlib testset run.\n";
      return 1;
    }
    parms[key]=value;
  }
  
  const char dirsep =  CoinFindDirSeparator();
  // Set directory containing mps data files.
  std::string mpsDir;
  if (parms.find("-mpsDir") != parms.end())
    mpsDir=parms["-mpsDir"] + dirsep;
  else 
    mpsDir = dirsep == '/' ? "../Mps/Sample/" : "..\\Mps\\Sample\\";
 
  // Set directory containing netlib data files.
  std::string netlibDir;
  if (parms.find("-netlibDir") != parms.end())
    netlibDir=parms["-netlibDir"] + dirsep;
  else 
    netlibDir = dirsep == '/' ? "../Mps/Netlib/" : "..\\Mps\\Netlib\\";

  testingMessage( "Testing ClpSimplex\n" );
  ClpSimplexUnitTest(mpsDir,netlibDir);
  if (parms.find("-netlib") != parms.end())
  {
    unsigned int m;
    
    // Define test problems: 
    //   mps names, 
    //   maximization or minimization, 
    //   Number of rows and columns in problem, and
    //   objective function value
    std::vector<std::string> mpsName;
    std::vector<bool> min;
    std::vector<int> nRows;
    std::vector<int> nCols;
    std::vector<double> objValue;
    std::vector<double> objValueTol;
    mpsName.push_back("25fv47");
    min.push_back(true);
    nRows.push_back(822);
    nCols.push_back(1571);
    objValueTol.push_back(1.E-10);
    objValue.push_back(5.5018458883E+03);

    mpsName.push_back("80bau3b");min.push_back(true);nRows.push_back(2263);nCols.push_back(9799);objValueTol.push_back(1.e-10);objValue.push_back(9.8722419241E+05);
    mpsName.push_back("blend");min.push_back(true);nRows.push_back(75);nCols.push_back(83);objValueTol.push_back(1.e-10);objValue.push_back(-3.0812149846e+01);
    mpsName.push_back("pilotnov");min.push_back(true);nRows.push_back(976);nCols.push_back(2172);objValueTol.push_back(1.e-10);objValue.push_back(-4.4972761882e+03);
    mpsName.push_back("maros-r7");min.push_back(true);nRows.push_back(3137);nCols.push_back(9408);objValueTol.push_back(1.e-10);objValue.push_back(1.4971851665e+06);
    
    mpsName.push_back("pilot");min.push_back(true);nRows.push_back(1442);nCols.push_back(3652);objValueTol.push_back(1.e-5);objValue.push_back(/*-5.5740430007e+02*/-557.48972927292);
    mpsName.push_back("pilot4");min.push_back(true);nRows.push_back(411);nCols.push_back(1000);objValueTol.push_back(1.e-8);objValue.push_back(-2.5811392641e+03);
    mpsName.push_back("pilot87");min.push_back(true);nRows.push_back(2031);nCols.push_back(4883);objValueTol.push_back(1.e-4);objValue.push_back(3.0171072827e+02);
    mpsName.push_back("adlittle");min.push_back(true);nRows.push_back(57);nCols.push_back(97);objValueTol.push_back(1.e-10);objValue.push_back(2.2549496316e+05);
    mpsName.push_back("afiro");min.push_back(true);nRows.push_back(28);nCols.push_back(32);objValueTol.push_back(1.e-10);objValue.push_back(-4.6475314286e+02);
    mpsName.push_back("agg");min.push_back(true);nRows.push_back(489);nCols.push_back(163);objValueTol.push_back(1.e-10);objValue.push_back(-3.5991767287e+07);
    mpsName.push_back("agg2");min.push_back(true);nRows.push_back(517);nCols.push_back(302);objValueTol.push_back(1.e-10);objValue.push_back(-2.0239252356e+07);
    mpsName.push_back("agg3");min.push_back(true);nRows.push_back(517);nCols.push_back(302);objValueTol.push_back(1.e-10);objValue.push_back(1.0312115935e+07);
    mpsName.push_back("bandm");min.push_back(true);nRows.push_back(306);nCols.push_back(472);objValueTol.push_back(1.e-10);objValue.push_back(-1.5862801845e+02);
    mpsName.push_back("beaconfd");min.push_back(true);nRows.push_back(174);nCols.push_back(262);objValueTol.push_back(1.e-10);objValue.push_back(3.3592485807e+04);
    mpsName.push_back("bnl1");min.push_back(true);nRows.push_back(644);nCols.push_back(1175);objValueTol.push_back(1.e-10);objValue.push_back(1.9776295615E+03);
    mpsName.push_back("bnl2");min.push_back(true);nRows.push_back(2325);nCols.push_back(3489);objValueTol.push_back(1.e-10);objValue.push_back(1.8112365404e+03);
    mpsName.push_back("boeing1");min.push_back(true);nRows.push_back(/*351*/352);nCols.push_back(384);objValueTol.push_back(1.e-10);objValue.push_back(-3.3521356751e+02);
    mpsName.push_back("boeing2");min.push_back(true);nRows.push_back(167);nCols.push_back(143);objValueTol.push_back(1.e-10);objValue.push_back(-3.1501872802e+02);
    mpsName.push_back("bore3d");min.push_back(true);nRows.push_back(234);nCols.push_back(315);objValueTol.push_back(1.e-10);objValue.push_back(1.3730803942e+03);
    mpsName.push_back("brandy");min.push_back(true);nRows.push_back(221);nCols.push_back(249);objValueTol.push_back(1.e-10);objValue.push_back(1.5185098965e+03);
    mpsName.push_back("capri");min.push_back(true);nRows.push_back(272);nCols.push_back(353);objValueTol.push_back(1.e-10);objValue.push_back(2.6900129138e+03);
    mpsName.push_back("cycle");min.push_back(true);nRows.push_back(1904);nCols.push_back(2857);objValueTol.push_back(1.e-9);objValue.push_back(-5.2263930249e+00);
    mpsName.push_back("czprob");min.push_back(true);nRows.push_back(930);nCols.push_back(3523);objValueTol.push_back(1.e-10);objValue.push_back(2.1851966989e+06);
    mpsName.push_back("d2q06c");min.push_back(true);nRows.push_back(2172);nCols.push_back(5167);objValueTol.push_back(1.e-7);objValue.push_back(122784.21557456);
    mpsName.push_back("d6cube");min.push_back(true);nRows.push_back(416);nCols.push_back(6184);objValueTol.push_back(1.e-7);objValue.push_back(3.1549166667e+02);
    mpsName.push_back("degen2");min.push_back(true);nRows.push_back(445);nCols.push_back(534);objValueTol.push_back(1.e-10);objValue.push_back(-1.4351780000e+03);
    mpsName.push_back("degen3");min.push_back(true);nRows.push_back(1504);nCols.push_back(1818);objValueTol.push_back(1.e-10);objValue.push_back(-9.8729400000e+02);
    mpsName.push_back("dfl001");min.push_back(true);nRows.push_back(6072);nCols.push_back(12230);objValueTol.push_back(1.e-5);objValue.push_back(1.1266396047E+07);
    mpsName.push_back("e226");min.push_back(true);nRows.push_back(224);nCols.push_back(282);objValueTol.push_back(1.e-10);objValue.push_back(-1.8751929066e+01+7.113); // The correct answer includes -7.113 term. This is a constant in the objective function. See line 1683 of the mps file.
    mpsName.push_back("etamacro");min.push_back(true);nRows.push_back(401);nCols.push_back(688);objValueTol.push_back(1.e-6);objValue.push_back(-7.5571521774e+02 );
    mpsName.push_back("fffff800");min.push_back(true);nRows.push_back(525);nCols.push_back(854);objValueTol.push_back(1.e-6);objValue.push_back(5.5567961165e+05);
    mpsName.push_back("finnis");min.push_back(true);nRows.push_back(498);nCols.push_back(614);objValueTol.push_back(1.e-6);objValue.push_back(1.7279096547e+05);
    mpsName.push_back("fit1d");min.push_back(true);nRows.push_back(25);nCols.push_back(1026);objValueTol.push_back(1.e-10);objValue.push_back(-9.1463780924e+03);
    mpsName.push_back("fit1p");min.push_back(true);nRows.push_back(628);nCols.push_back(1677);objValueTol.push_back(1.e-10);objValue.push_back(9.1463780924e+03);
    mpsName.push_back("fit2d");min.push_back(true);nRows.push_back(26);nCols.push_back(10500);objValueTol.push_back(1.e-10);objValue.push_back(-6.8464293294e+04);
    mpsName.push_back("fit2p");min.push_back(true);nRows.push_back(3001);nCols.push_back(13525);objValueTol.push_back(1.e-9);objValue.push_back(6.8464293232e+04);
    mpsName.push_back("forplan");min.push_back(true);nRows.push_back(162);nCols.push_back(421);objValueTol.push_back(1.e-6);objValue.push_back(-6.6421873953e+02);
    mpsName.push_back("ganges");min.push_back(true);nRows.push_back(1310);nCols.push_back(1681);objValueTol.push_back(1.e-5);objValue.push_back(-1.0958636356e+05);
    mpsName.push_back("gfrd-pnc");min.push_back(true);nRows.push_back(617);nCols.push_back(1092);objValueTol.push_back(1.e-10);objValue.push_back(6.9022359995e+06);
    mpsName.push_back("greenbea");min.push_back(true);nRows.push_back(2393);nCols.push_back(5405);objValueTol.push_back(1.e-10);objValue.push_back(/*-7.2462405908e+07*/-72555248.129846);
    mpsName.push_back("greenbeb");min.push_back(true);nRows.push_back(2393);nCols.push_back(5405);objValueTol.push_back(1.e-10);objValue.push_back(/*-4.3021476065e+06*/-4302260.2612066);
    mpsName.push_back("grow15");min.push_back(true);nRows.push_back(301);nCols.push_back(645);objValueTol.push_back(1.e-10);objValue.push_back(-1.0687094129e+08);
    mpsName.push_back("grow22");min.push_back(true);nRows.push_back(441);nCols.push_back(946);objValueTol.push_back(1.e-10);objValue.push_back(-1.6083433648e+08);
    mpsName.push_back("grow7");min.push_back(true);nRows.push_back(141);nCols.push_back(301);objValueTol.push_back(1.e-10);objValue.push_back(-4.7787811815e+07);
    mpsName.push_back("israel");min.push_back(true);nRows.push_back(175);nCols.push_back(142);objValueTol.push_back(1.e-10);objValue.push_back(-8.9664482186e+05);
    mpsName.push_back("kb2");min.push_back(true);nRows.push_back(44);nCols.push_back(41);objValueTol.push_back(1.e-10);objValue.push_back(-1.7499001299e+03);
    mpsName.push_back("lotfi");min.push_back(true);nRows.push_back(154);nCols.push_back(308);objValueTol.push_back(1.e-10);objValue.push_back(-2.5264706062e+01);
    mpsName.push_back("maros");min.push_back(true);nRows.push_back(847);nCols.push_back(1443);objValueTol.push_back(1.e-10);objValue.push_back(-5.8063743701e+04);
    mpsName.push_back("modszk1");min.push_back(true);nRows.push_back(688);nCols.push_back(1620);objValueTol.push_back(1.e-10);objValue.push_back(3.2061972906e+02);
    mpsName.push_back("nesm");min.push_back(true);nRows.push_back(663);nCols.push_back(2923);objValueTol.push_back(1.e-5);objValue.push_back(1.4076073035e+07);
    mpsName.push_back("perold");min.push_back(true);nRows.push_back(626);nCols.push_back(1376);objValueTol.push_back(1.e-6);objValue.push_back(-9.3807580773e+03);
    //mpsName.push_back("qap12");min.push_back(true);nRows.push_back(3193);nCols.push_back(8856);objValueTol.push_back(1.e-6);objValue.push_back(5.2289435056e+02);
    //mpsName.push_back("qap15");min.push_back(true);nRows.push_back(6331);nCols.push_back(22275);objValueTol.push_back(1.e-10);objValue.push_back(1.0409940410e+03);
    mpsName.push_back("recipe");min.push_back(true);nRows.push_back(92);nCols.push_back(180);objValueTol.push_back(1.e-10);objValue.push_back(-2.6661600000e+02);
    mpsName.push_back("sc105");min.push_back(true);nRows.push_back(106);nCols.push_back(103);objValueTol.push_back(1.e-10);objValue.push_back(-5.2202061212e+01);
    mpsName.push_back("sc205");min.push_back(true);nRows.push_back(206);nCols.push_back(203);objValueTol.push_back(1.e-10);objValue.push_back(-5.2202061212e+01);
    mpsName.push_back("sc50a");min.push_back(true);nRows.push_back(51);nCols.push_back(48);objValueTol.push_back(1.e-10);objValue.push_back(-6.4575077059e+01);
    mpsName.push_back("sc50b");min.push_back(true);nRows.push_back(51);nCols.push_back(48);objValueTol.push_back(1.e-10);objValue.push_back(-7.0000000000e+01);
    mpsName.push_back("scagr25");min.push_back(true);nRows.push_back(472);nCols.push_back(500);objValueTol.push_back(1.e-10);objValue.push_back(-1.4753433061e+07);
    mpsName.push_back("scagr7");min.push_back(true);nRows.push_back(130);nCols.push_back(140);objValueTol.push_back(1.e-6);objValue.push_back(-2.3313892548e+06);
    mpsName.push_back("scfxm1");min.push_back(true);nRows.push_back(331);nCols.push_back(457);objValueTol.push_back(1.e-10);objValue.push_back(1.8416759028e+04);
    mpsName.push_back("scfxm2");min.push_back(true);nRows.push_back(661);nCols.push_back(914);objValueTol.push_back(1.e-10);objValue.push_back(3.6660261565e+04);
    mpsName.push_back("scfxm3");min.push_back(true);nRows.push_back(991);nCols.push_back(1371);objValueTol.push_back(1.e-10);objValue.push_back(5.4901254550e+04);
    mpsName.push_back("scorpion");min.push_back(true);nRows.push_back(389);nCols.push_back(358);objValueTol.push_back(1.e-10);objValue.push_back(1.8781248227e+03);
    mpsName.push_back("scrs8");min.push_back(true);nRows.push_back(491);nCols.push_back(1169);objValueTol.push_back(1.e-5);objValue.push_back(9.0429998619e+02);
    mpsName.push_back("scsd1");min.push_back(true);nRows.push_back(78);nCols.push_back(760);objValueTol.push_back(1.e-10);objValue.push_back(8.6666666743e+00);
    mpsName.push_back("scsd6");min.push_back(true);nRows.push_back(148);nCols.push_back(1350);objValueTol.push_back(1.e-10);objValue.push_back(5.0500000078e+01);
    mpsName.push_back("scsd8");min.push_back(true);nRows.push_back(398);nCols.push_back(2750);objValueTol.push_back(1.e-10);objValue.push_back(9.0499999993e+02);
    mpsName.push_back("sctap1");min.push_back(true);nRows.push_back(301);nCols.push_back(480);objValueTol.push_back(1.e-10);objValue.push_back(1.4122500000e+03);
    mpsName.push_back("sctap2");min.push_back(true);nRows.push_back(1091);nCols.push_back(1880);objValueTol.push_back(1.e-10);objValue.push_back(1.7248071429e+03);
    mpsName.push_back("sctap3");min.push_back(true);nRows.push_back(1481);nCols.push_back(2480);objValueTol.push_back(1.e-10);objValue.push_back(1.4240000000e+03);
    mpsName.push_back("seba");min.push_back(true);nRows.push_back(516);nCols.push_back(1028);objValueTol.push_back(1.e-10);objValue.push_back(1.5711600000e+04);
    mpsName.push_back("share1b");min.push_back(true);nRows.push_back(118);nCols.push_back(225);objValueTol.push_back(1.e-10);objValue.push_back(-7.6589318579e+04);
    mpsName.push_back("share2b");min.push_back(true);nRows.push_back(97);nCols.push_back(79);objValueTol.push_back(1.e-10);objValue.push_back(-4.1573224074e+02);
    mpsName.push_back("shell");min.push_back(true);nRows.push_back(537);nCols.push_back(1775);objValueTol.push_back(1.e-10);objValue.push_back(1.2088253460e+09);
    mpsName.push_back("ship04l");min.push_back(true);nRows.push_back(403);nCols.push_back(2118);objValueTol.push_back(1.e-10);objValue.push_back(1.7933245380e+06);
    mpsName.push_back("ship04s");min.push_back(true);nRows.push_back(403);nCols.push_back(1458);objValueTol.push_back(1.e-10);objValue.push_back(1.7987147004e+06);
    mpsName.push_back("ship08l");min.push_back(true);nRows.push_back(779);nCols.push_back(4283);objValueTol.push_back(1.e-10);objValue.push_back(1.9090552114e+06);
    mpsName.push_back("ship08s");min.push_back(true);nRows.push_back(779);nCols.push_back(2387);objValueTol.push_back(1.e-10);objValue.push_back(1.9200982105e+06);
    mpsName.push_back("ship12l");min.push_back(true);nRows.push_back(1152);nCols.push_back(5427);objValueTol.push_back(1.e-10);objValue.push_back(1.4701879193e+06);
    mpsName.push_back("ship12s");min.push_back(true);nRows.push_back(1152);nCols.push_back(2763);objValueTol.push_back(1.e-10);objValue.push_back(1.4892361344e+06);
    mpsName.push_back("sierra");min.push_back(true);nRows.push_back(1228);nCols.push_back(2036);objValueTol.push_back(1.e-10);objValue.push_back(1.5394362184e+07);
    mpsName.push_back("stair");min.push_back(true);nRows.push_back(357);nCols.push_back(467);objValueTol.push_back(1.e-10);objValue.push_back(-2.5126695119e+02);
    mpsName.push_back("standata");min.push_back(true);nRows.push_back(360);nCols.push_back(1075);objValueTol.push_back(1.e-10);objValue.push_back(1.2576995000e+03);
    //mpsName.push_back("standgub");min.push_back(true);nRows.push_back(362);nCols.push_back(1184);objValueTol.push_back(1.e-10);objValue.push_back(1257.6995); 
    mpsName.push_back("standmps");min.push_back(true);nRows.push_back(468);nCols.push_back(1075);objValueTol.push_back(1.e-10);objValue.push_back(1.4060175000E+03); 
    mpsName.push_back("stocfor1");min.push_back(true);nRows.push_back(118);nCols.push_back(111);objValueTol.push_back(1.e-10);objValue.push_back(-4.1131976219E+04);
    mpsName.push_back("stocfor2");min.push_back(true);nRows.push_back(2158);nCols.push_back(2031);objValueTol.push_back(1.e-10);objValue.push_back(-3.9024408538e+04);
    //mpsName.push_back("stocfor3");min.push_back(true);nRows.push_back(16676);nCols.push_back(15695);objValueTol.push_back(1.e-10);objValue.push_back(-3.9976661576e+04);
    //mpsName.push_back("truss");min.push_back(true);nRows.push_back(1001);nCols.push_back(8806);objValueTol.push_back(1.e-10);objValue.push_back(4.5881584719e+05);
    mpsName.push_back("tuff");min.push_back(true);nRows.push_back(334);nCols.push_back(587);objValueTol.push_back(1.e-10);objValue.push_back(2.9214776509e-01);
    mpsName.push_back("vtpbase");min.push_back(true);nRows.push_back(199);nCols.push_back(203);objValueTol.push_back(1.e-10);objValue.push_back(1.2983146246e+05);
    mpsName.push_back("wood1p");min.push_back(true);nRows.push_back(245);nCols.push_back(2594);objValueTol.push_back(1.e-10);objValue.push_back(1.4429024116e+00);
    mpsName.push_back("woodw");min.push_back(true);nRows.push_back(1099);nCols.push_back(8405);objValueTol.push_back(1.e-10);objValue.push_back(1.3044763331E+00);

    double timeTaken =0.0;
  // Loop once for each Mps File
    for (m=0; m<mpsName.size(); m++ ) {
      std::cerr <<"  processing mps file: " <<mpsName[m] 
		<<" (" <<m+1 <<" out of " <<mpsName.size() <<")" <<std::endl;
    
      // Read data mps file,
      std::string fn = netlibDir+mpsName[m];
      CoinMpsIO mps;
      mps.readMps(fn.c_str(),"mps");
      double time1 = cpuTime();
      ClpSimplex solution=empty;
      solution.loadProblem(*mps.getMatrixByCol(),mps.getColLower(),
			   mps.getColUpper(),
			   mps.getObjCoefficients(),
			   mps.getRowLower(),mps.getRowUpper());

      solution.setDblParam(ClpObjOffset,mps.objectiveOffset());
#if 0
      solution.setOptimizationDirection(-1);
      {
	int j;
	double * obj = solution.objective();
	int n=solution.numberColumns();
	for (j=0;j<n;j++) 
	  obj[j] *= -1.0;
      }
#endif
      if (doPresolve) {
#ifdef USE_PRESOLVE
	Presolve pinfo;
	ClpSimplex * model2 = pinfo.presolvedModel(solution,1.0e-8);
	// change from 200
	model2->factorization()->maximumPivots(100+model2->numberRows()/50);
	if (doDual) {
	  // faster if bounds tightened
	  int numberInfeasibilities = model2->tightenPrimalBounds();
	  if (numberInfeasibilities)
	    std::cout<<"** Analysis indicates model infeasible"
		     <<std::endl;
	  if (doIdiot<0)
	    model2->crash(1000,2);
	  model2->dual();
	} else {
#ifdef CLP_IDIOT
	  if (doIdiot>0) {
	    Idiot info(*model2);
	    info.crash(doIdiot);
	  }
#endif
	  model2->primal(1);
	}
	pinfo.postsolve(true);
	
	delete model2;
#if 1
	printf("Resolving from postsolved model\n");
	// later try without (1) and check duals before solve
	solution.primal(1);
	if (solution.numberIterations())
	  printf("****** iterated %d\n",solution.numberIterations());
	/*solution.checkSolution();
	printf("%g dual %g(%d) Primal %g(%d)\n",
	       solution.objectiveValue(),
	       solution.sumDualInfeasibilities(),
	       solution.numberDualInfeasibilities(),
	       solution.sumPrimalInfeasibilities(),
	       solution.numberPrimalInfeasibilities());*/
#endif
	if (0) {
	  Presolve pinfoA;
	  model2 = pinfoA.presolvedModel(solution,1.0e-8);

	  printf("Resolving from presolved optimal solution\n");
	  model2->primal(1);
		
	  delete model2;
	}
#else
	if (doDual) {
	  if (doIdiot<0)
	    solution.crash(1000,2);
	  solution.dual();
	} else {
#ifdef CLP_IDIOT
	  if (doIdiot>0) {
	    Idiot info(solution);
	    info.crash(doIdiot);
	  }
#endif
	  solution.primal(1);
	}
#endif
      } else {
	if (doDual) {
	  if (doIdiot<0)
	    solution.crash(1000,2);
	  solution.dual();
	} else {
#ifdef CLP_IDIOT
	  if (doIdiot>0) {
	    Idiot info(solution);
	    info.crash(doIdiot);
	  }
#endif
	  solution.primal(1);
	}
      }
      double time2 = cpuTime()-time1;
      timeTaken += time2;
      printf("Took %g seconds\n",time2);
      // Test objective solution value
      {
        double soln = solution.objectiveValue();
        CoinRelFltEq eq(objValueTol[m]);
        std::cerr <<soln <<",  " <<objValue[m] <<" diff "<<
	  soln-objValue[m]<<std::endl;
        if(!eq(soln,objValue[m]))
	  printf("** difference fails\n");
      }
    }
    printf("Total time %g seconds\n",timeTaken);
  }
  else {
    testingMessage( "***Skipped Testing on netlib    ***\n" );
    testingMessage( "***use -netlib to test class***\n" );
  }
  
  testingMessage( "All tests completed successfully\n" );
  return 0;
}

 
// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
  std::cerr <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}

//--------------------------------------------------------------------------
// test factorization methods and simplex method
void
ClpSimplexUnitTest(const std::string & mpsDir,
		   const std::string & netlibDir)
{
  
  CoinRelFltEq eq(0.000001);

  {
    ClpSimplex solution;
  
    // matrix data
    //deliberate hiccup of 2 between 0 and 1
    CoinBigIndex start[5]={0,4,7,8,9};
    int length[5]={2,3,1,1,1};
    int rows[11]={0,2,-1,-1,0,1,2,0,1,2};
    double elements[11]={7.0,2.0,1.0e10,1.0e10,-2.0,1.0,-2.0,1,1,1};
    CoinPackedMatrix matrix(true,3,5,8,elements,rows,start,length);
    
    // rim data
    double objective[7]={-4.0,1.0,0.0,0.0,0.0,0.0,0.0};
    double rowLower[5]={14.0,3.0,3.0,1.0e10,1.0e10};
    double rowUpper[5]={14.0,3.0,3.0,-1.0e10,-1.0e10};
    double colLower[7]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double colUpper[7]={100.0,100.0,100.0,100.0,100.0,100.0,100.0};
    
    // basis 1
    int rowBasis1[3]={-1,-1,-1};
    int colBasis1[5]={1,1,-1,-1,1};
    solution.loadProblem(matrix,colLower,colUpper,objective,
			 rowLower,rowUpper);
    int i;
    solution.createStatus();
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	solution.setRowStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution.setRowStatus(i,ClpSimplex::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	solution.setColumnStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution.setColumnStatus(i,ClpSimplex::basic);
      }
    }
    solution.setLogLevel(3+4+8+16+32);
    solution.primal();
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	solution.setRowStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution.setRowStatus(i,ClpSimplex::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	solution.setColumnStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution.setColumnStatus(i,ClpSimplex::basic);
      }
    }
    // intricate stuff does not work with scaling
    solution.scaling(0);
    assert(!solution.factorize ( ));
    const double * colsol = solution.primalColumnSolution();
    const double * rowsol = solution.primalRowSolution();
    solution.getSolution(rowsol,colsol);
    double colsol1[5]={20.0/7.0,3.0,0.0,0.0,23.0/7.0};
    for (i=0;i<5;i++) {
      assert(eq(colsol[i],colsol1[i]));
    }
    // now feed in again without actually doing factorization
    ClpFactorization factorization2 = *solution.factorization();
    ClpSimplex solution2 = solution;
    solution2.setFactorization(factorization2);
    solution2.createStatus();
    for (i=0;i<3;i++) {
      if (rowBasis1[i]<0) {
	solution2.setRowStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution2.setRowStatus(i,ClpSimplex::basic);
      }
    }
    for (i=0;i<5;i++) {
      if (colBasis1[i]<0) {
	solution2.setColumnStatus(i,ClpSimplex::atLowerBound);
      } else {
	solution2.setColumnStatus(i,ClpSimplex::basic);
      }
    }
    // intricate stuff does not work with scaling
    solution2.scaling(0);
    solution2.getSolution(rowsol,colsol);
    colsol = solution2.primalColumnSolution();
    rowsol = solution2.primalRowSolution();
    for (i=0;i<5;i++) {
      assert(eq(colsol[i],colsol1[i]));
    }
    solution2.setDualBound(0.1);
    solution2.dual();
    objective[2]=-1.0;
    objective[3]=-0.5;
    objective[4]=10.0;
    solution.dual();
    for (i=0;i<3;i++) {
      rowLower[i]=-1.0e50;
      colUpper[i+2]=0.0;
    }
    solution.setLogLevel(3);
    solution.dual();
    double rowObjective[]={1.0,0.5,-10.0};
    solution.loadProblem(matrix,colLower,colUpper,objective,
			 rowLower,rowUpper,rowObjective);
    solution.dual();
  }
  {    
    CoinMpsIO m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    ClpSimplex solution;
    solution.loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
			 m.getObjCoefficients(),
			 m.getRowLower(),m.getRowUpper());
    solution.dual();
  }
  // test steepest edge
  {    
    CoinMpsIO m;
    std::string fn = netlibDir+"finnis";
    m.readMps(fn.c_str(),"mps");
    ClpModel model;
    model.loadProblem(*m.getMatrixByCol(),m.getColLower(),
		    m.getColUpper(),
		    m.getObjCoefficients(),
		    m.getRowLower(),m.getRowUpper());
    ClpSimplex solution(model);

    solution.scaling(1); 
    solution.setDualBound(1.0e8);
    //solution.factorization()->maximumPivots(1);
    //solution.setLogLevel(3);
    solution.setDualTolerance(1.0e-7);
    // set objective sense,
    ClpDualRowSteepest steep;
    solution.setDualRowPivotAlgorithm(steep);
    solution.setDblParam(ClpObjOffset,m.objectiveOffset());
    solution.dual();
  }
  // test normal solution
  {    
    CoinMpsIO m;
    std::string fn = netlibDir+"afiro";
    m.readMps(fn.c_str(),"mps");
    ClpSimplex solution;
    ClpModel model;
    // do twice - without and with scaling
    int iPass;
    for (iPass=0;iPass<2;iPass++) {
      // explicit row objective for testing
      int nr = m.getNumRows();
      double * rowObj = new double[nr];
      CoinFillN(rowObj,nr,0.0);
      model.loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		      m.getObjCoefficients(),
		      m.getRowLower(),m.getRowUpper(),rowObj);
      delete [] rowObj;
      solution = ClpSimplex(model);
      if (iPass) {
	solution.scaling();
      }
      solution.dual();
      solution.dual();
      // test optimal
      assert (solution.status()==0);
      int numberColumns = solution.numberColumns();
      int numberRows = solution.numberRows();
      CoinPackedVector colsol(numberColumns,solution.primalColumnSolution());
      double * objective = solution.objective();
      double objValue = colsol.dotProduct(objective);
      CoinRelFltEq eq(1.0e-8);
      assert(eq(objValue,-4.6475314286e+02));
      double * lower = solution.columnLower();
      double * upper = solution.columnUpper();
      double * sol = solution.primalColumnSolution();
      double * result = new double[numberColumns];
      CoinFillN ( result, numberColumns,0.0);
      solution.matrix()->transposeTimes(solution.dualRowSolution(), result);
      int iRow , iColumn;
      // see if feasible and dual feasible
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = sol[iColumn];
	assert(value<upper[iColumn]+1.0e-8);
	assert(value>lower[iColumn]-1.0e-8);
	value = objective[iColumn]-result[iColumn];
	assert (value>-1.0e-5);
	if (sol[iColumn]>1.0e-5)
	  assert (value<1.0e-5);
      }
      delete [] result;
      result = new double[numberRows];
      CoinFillN ( result, numberRows,0.0);
      solution.matrix()->times(colsol, result);
      lower = solution.rowLower();
      upper = solution.rowUpper();
      sol = solution.primalRowSolution();
      for (iRow=0;iRow<numberRows;iRow++) {
	double value = result[iRow];
	assert(eq(value,sol[iRow]));
	assert(value<upper[iRow]+1.0e-8);
	assert(value>lower[iRow]-1.0e-8);
      }
      delete [] result;
      // test row objective
      double * rowObjective = solution.rowObjective();
      CoinDisjointCopyN(solution.dualRowSolution(),numberRows,rowObjective);
      CoinDisjointCopyN(solution.dualColumnSolution(),numberColumns,objective);
      // this sets up all slack basis
      solution.createStatus();
      solution.dual();
      CoinFillN(rowObjective,numberRows,0.0);
      CoinDisjointCopyN(m.getObjCoefficients(),numberColumns,objective);
      solution.dual();
    }
  }
  // test unbounded
  {    
    CoinMpsIO m;
    std::string fn = netlibDir+"brandy";
    m.readMps(fn.c_str(),"mps");
    ClpSimplex solution;
    // do twice - without and with scaling
    int iPass;
    for (iPass=0;iPass<2;iPass++) {
      solution.loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		      m.getObjCoefficients(),
		      m.getRowLower(),m.getRowUpper());
      if (iPass)
	solution.scaling();
      solution.setOptimizationDirection(-1);
      // test unbounded and ray
#ifdef DUAL
      solution.setDualBound(100.0);
      solution.dual();
#else
      solution.primal();
#endif
      assert (solution.status()==2);
      int numberColumns = solution.numberColumns();
      int numberRows = solution.numberRows();
      double * lower = solution.columnLower();
      double * upper = solution.columnUpper();
      double * sol = solution.primalColumnSolution();
      double * ray = solution.unboundedRay();
      double * objective = solution.objective();
      double objChange=0.0;
      int iRow , iColumn;
      // make sure feasible and columns form ray
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = sol[iColumn];
	assert(value<upper[iColumn]+1.0e-8);
	assert(value>lower[iColumn]-1.0e-8);
	value = ray[iColumn];
	if (value>0.0)
	  assert(upper[iColumn]>1.0e30);
	else if (value<0.0)
	  assert(lower[iColumn]<-1.0e30);
	objChange += value*objective[iColumn];
      }
      // make sure increasing objective
      assert(objChange>0.0);
      double * result = new double[numberRows];
      CoinFillN ( result, numberRows,0.0);
      solution.matrix()->times(sol, result);
      lower = solution.rowLower();
      upper = solution.rowUpper();
      sol = solution.primalRowSolution();
      for (iRow=0;iRow<numberRows;iRow++) {
	double value = result[iRow];
	assert(eq(value,sol[iRow]));
	assert(value<upper[iRow]+1.0e-8);
	assert(value>lower[iRow]-1.0e-8);
      }
      CoinFillN ( result, numberRows,0.0);
      solution.matrix()->times(ray, result);
      // there may be small differences (especially if scaled)
      for (iRow=0;iRow<numberRows;iRow++) {
	double value = result[iRow];
	if (value>1.0e-8)
	  assert(upper[iRow]>1.0e30);
	else if (value<-1.0e-8)
	  assert(lower[iRow]<-1.0e30);
      }
      delete [] result;
      delete [] ray;
    }
  }
  // test infeasible
  {    
    CoinMpsIO m;
    std::string fn = netlibDir+"brandy";
    m.readMps(fn.c_str(),"mps");
    ClpSimplex solution;
    // do twice - without and with scaling
    int iPass;
    for (iPass=0;iPass<2;iPass++) {
      solution.loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		      m.getObjCoefficients(),
		      m.getRowLower(),m.getRowUpper());
      if (iPass)
	solution.scaling();
      // test infeasible and ray
      solution.columnUpper()[0]=0.0;
#ifdef DUAL
      solution.setDualBound(100.0);
      solution.dual();
#else
      solution.primal();
#endif
      assert (solution.status()==1);
      int numberColumns = solution.numberColumns();
      int numberRows = solution.numberRows();
      double * lower = solution.rowLower();
      double * upper = solution.rowUpper();
      double * ray = solution.infeasibilityRay();
      assert(ray);
      // construct proof of infeasibility
      int iRow , iColumn;
      double lo=0.0,up=0.0;
      int nl=0,nu=0;
      for (iRow=0;iRow<numberRows;iRow++) {
	if (lower[iRow]>-1.0e20) {
	  lo += ray[iRow]*lower[iRow];
	} else {
	  if (ray[iRow]>1.0e-8) 
	    nl++;
	}
	if (upper[iRow]<1.0e20) {
	  up += ray[iRow]*upper[iRow];
	} else {
	  if (ray[iRow]>1.0e-8) 
	    nu++;
	}
      }
      if (nl)
	lo=-1.0e100;
      if (nu)
	up=1.0e100;
      double * result = new double[numberColumns];
      double lo2=0.0,up2=0.0;
      CoinFillN ( result, numberColumns,0.0);
      solution.matrix()->transposeTimes(ray, result);
      lower = solution.columnLower();
      upper = solution.columnUpper();
      nl=nu=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (result[iColumn]>1.0e-8) {
	  if (lower[iColumn]>-1.0e20)
	    lo2 += result[iColumn]*lower[iColumn];
	  else
	    nl++;
	  if (upper[iColumn]<1.0e20)
	    up2 += result[iColumn]*upper[iColumn];
	  else
	    nu++;
	} else if (result[iColumn]<-1.0e-8) {
	  if (lower[iColumn]>-1.0e20)
	    up2 += result[iColumn]*lower[iColumn];
	  else
	    nu++;
	  if (upper[iColumn]<1.0e20)
	    lo2 += result[iColumn]*upper[iColumn];
	  else
	    nl++;
	}
      }
      if (nl)
	lo2=-1.0e100;
      if (nu)
	up2=1.0e100;
      // make sure inconsistency
      assert(lo2>up||up2<lo);
      delete [] result;
      delete [] ray;
    }
  }
  
}
