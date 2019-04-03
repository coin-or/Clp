set(NAME "CLP")

# PTHREAD
if(PTHREADS_FOUND)
  set(${NAME}_PTHREADS "1" CACHE INTERNAL "Use pthread")
endif()

set(COIN_${NAME}_CHECKLEVEL "0" CACHE INTERNAL
  "${NAME} check level")
set(COIN_${NAME}_VERBOSITY "0" CACHE INTERNAL
  "${NAME} verbosity level")
configure_file(config.h.cmake.in config.h)
configure_file(config_clp.h.cmake.in config_clp.h)

set(_SRCS
  ClpCholeskyBase.cpp
  ClpCholeskyDense.cpp
  ClpConstraint.cpp
  ClpConstraintLinear.cpp
  ClpConstraintQuadratic.cpp
  Clp_C_Interface.cpp
  ClpDualRowDantzig.cpp
  ClpDualRowPivot.cpp
  ClpDualRowSteepest.cpp
  ClpDummyMatrix.cpp
  ClpDynamicExampleMatrix.cpp
  ClpDynamicMatrix.cpp
  ClpEventHandler.cpp
  ClpFactorization.cpp
  ClpGubDynamicMatrix.cpp
  ClpGubMatrix.cpp
  ClpHelperFunctions.cpp
  ClpInterior.cpp
  ClpLinearObjective.cpp
  ClpMatrixBase.cpp
  ClpMessage.cpp
  ClpModel.cpp
  ClpNetworkBasis.cpp
  ClpNetworkMatrix.cpp
  ClpNonLinearCost.cpp
  ClpNode.cpp
  ClpObjective.cpp
  ClpPackedMatrix.cpp
  ClpPlusMinusOneMatrix.cpp
  ClpPredictorCorrector.cpp
  ClpPdco.cpp
  ClpPdcoBase.cpp
  ClpLsqr.cpp
  ClpPresolve.cpp
  ClpPrimalColumnDantzig.cpp
  ClpPrimalColumnPivot.cpp
  ClpPrimalColumnSteepest.cpp
  ClpQuadraticObjective.cpp
  ClpSimplex.cpp
  ClpSimplexDual.cpp
  ClpSimplexNonlinear.cpp
  ClpSimplexOther.cpp
  ClpSimplexPrimal.cpp
  ClpSolve.cpp
  Idiot.cpp
  IdiSolve.cpp
  ClpCholeskyPardiso.cpp
  ClpPESimplex.cpp
  ClpPEDualRowDantzig.cpp
  ClpPEDualRowSteepest.cpp
  ClpPEPrimalColumnDantzig.cpp
  ClpPEPrimalColumnSteepest.cpp
  )

set(_HDRS
  ClpConfig.h
  ClpCholeskyBase.hpp
  ClpCholeskyDense.hpp
  ClpConstraint.hpp
  ClpConstraintLinear.hpp
  ClpConstraintQuadratic.hpp
  Clp_C_Interface.h
  ClpDualRowDantzig.hpp
  ClpDualRowPivot.hpp
  ClpDualRowSteepest.hpp
  ClpDummyMatrix.hpp
  ClpDynamicExampleMatrix.hpp
  ClpDynamicMatrix.hpp
  ClpEventHandler.hpp
  ClpFactorization.hpp
  ClpGubDynamicMatrix.hpp
  ClpGubMatrix.hpp
  ClpHelperFunctions.hpp
  ClpInterior.hpp
  ClpLinearObjective.hpp
  ClpMatrixBase.hpp
  ClpMessage.hpp
  ClpModel.hpp
  ClpNetworkBasis.hpp
  ClpNetworkMatrix.hpp
  ClpNonLinearCost.hpp
  ClpNode.hpp
  ClpObjective.hpp
  ClpPackedMatrix.hpp
  ClpParameters.hpp
  ClpPlusMinusOneMatrix.hpp
  ClpPredictorCorrector.hpp
  ClpPdco.hpp
  ClpPdcoBase.hpp
  ClpLsqr.hpp
  ClpPresolve.hpp
  ClpPrimalColumnDantzig.hpp
  ClpPrimalColumnPivot.hpp
  ClpPrimalColumnSteepest.hpp
  ClpQuadraticObjective.hpp
  ClpSimplex.hpp
  ClpSimplexDual.hpp
  ClpSimplexNonlinear.hpp
  ClpSimplexOther.hpp
  ClpSimplexPrimal.hpp
  ClpSolve.hpp
  Idiot.hpp
  ClpCholeskyPardiso.hpp
  ClpPESimplex.hpp
  ClpPEDualRowDantzig.hpp
  ClpPEDualRowSteepest.hpp
  ClpPEPrimalColumnDantzig.hpp
  ClpPEPrimalColumnSteepest.hpp
  )

add_library(Clp ${_SRCS} ${_HDRS})
target_include_directories(Clp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(Clp
  PUBLIC HAVE_CONFIG_H
  PRIVATE CLP_BUILD COIN_HAS_CLP)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET Clp PROPERTY CXX_STANDARD 11)
  set_property(TARGET Clp PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(Clp PUBLIC cxx_std_11)
endif()
target_link_libraries(Clp PUBLIC
  Coin::CoinUtils
  Coin::Osi)
set_target_properties(Clp PROPERTIES
  PUBLIC_HEADER "${_HDRS};CbcOrClpParam.cpp;${CMAKE_CURRENT_BINARY_DIR}/config_clp.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::Clp ALIAS Clp)

# Install
include(GNUInstallDirs)
install(TARGETS Clp
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

#################
##  ClpSolver  ##
#################
set(_SRCS
  ClpSolver.cpp
  CbcOrClpParam.cpp
  MyEventHandler.cpp
  MyMessageHandler.cpp)

set(_HDRS
  CbcOrClpParam.hpp
  MyEventHandler.hpp
  MyMessageHandler.hpp
  unitTest.cpp)

add_library(ClpSolver ${_SRCS} ${_HDRS})
target_include_directories(ClpSolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(ClpSolver
  PUBLIC HAVE_CONFIG_H
  PRIVATE CLP_BUILD COIN_HAS_CLP)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET ClpSolver PROPERTY CXX_STANDARD 11)
  set_property(TARGET ClpSolver PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(ClpSolver PUBLIC cxx_std_11)
endif()
target_link_libraries(ClpSolver PUBLIC Coin::Clp)
set_target_properties(ClpSolver PROPERTIES
  PUBLIC_HEADER "${_HDRS}"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::ClpSolver ALIAS ClpSolver)

# Install
install(TARGETS ClpSolver
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

##############
##  OsiClp  ##
##############
add_subdirectory(OsiClp)
