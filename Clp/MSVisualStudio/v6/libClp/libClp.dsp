# Microsoft Developer Studio Project File - Name="libClp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libClp - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libClp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libClp.mak" CFG="libClp - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libClp - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libClp - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libClp - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GR /GX /O2 /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libClp - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /I "..\..\..\..\CoinUtils\src" /I "..\..\..\..\BuildTools\headers" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "libClp - Win32 Release"
# Name "libClp - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\..\..\Clp\src\CbcOrClpParam.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\Clp_C_Interface.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyBase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyDense.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyTaucs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyUfl.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyWssmp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyWssmpKKT.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowDantzig.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowPivot.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowSteepest.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDummyMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDynamicExampleMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDynamicMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpEventHandler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpFactorization.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpGubDynamicMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpGubMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpHelperFunctions.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpInterior.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpLinearObjective.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMain.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMatrixBase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMessage.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpModel.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNetworkBasis.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNetworkMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNonLinearCost.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpObjective.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPackedMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPlusMinusOneMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPredictorCorrector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPresolve.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnDantzig.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnPivot.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnSteepest.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpQuadraticObjective.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplex.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexDual.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexNonlinear.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexOther.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexPrimal.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSolve.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\Idiot.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\IdiSolve.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\..\..\Clp\src\CbcOrClpParam.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\Clp_C_Interface.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpCholeskyDense.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpConfig.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowDantzig.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowPivot.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDualRowSteepest.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDummyMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDynamicExampleMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpDynamicMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpEventHandler.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpFactorization.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpGubDynamicMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpGubMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpHelperFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpInterior.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpLinearObjective.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMatrixBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpModel.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNetworkBasis.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNetworkMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpNonLinearCost.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpObjective.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPackedMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpParameters.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPlusMinusOneMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPredictorCorrector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPresolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnDantzig.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnPivot.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpPrimalColumnSteepest.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpQuadraticObjective.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplex.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexDual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexNonlinear.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexOther.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSimplexPrimal.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\ClpSolve.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\Coin_C_defines.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinBuild.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinDenseVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinDistance.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinError.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFactorization.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFileIO.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFinite.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinFloatEqual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinHelperFunctions.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinIndexedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessage.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMessageHandler.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinModel.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinModelUseful.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinMpsIO.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPackedVectorBase.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPragma.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveDoubleton.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveDual.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveDupcol.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveEmpty.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveFixed.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveForcing.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveImpliedFree.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveIsolated.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveMatrix.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolvePsdebug.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveSingleton.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveSubst.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveTighten.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveTripleton.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveUseless.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinPresolveZeros.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinShallowPackedVector.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinSignal.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinSort.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinTime.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStart.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\CoinUtils\src\CoinWarmStartBasis.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\..\BuildTools\headers\configall_system.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\BuildTools\headers\configall_system_msc.h
# End Source File
# Begin Source File

SOURCE=..\..\..\..\Clp\src\Idiot.hpp
# End Source File
# End Group
# End Target
# End Project
