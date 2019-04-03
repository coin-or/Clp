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
  src/ClpCholeskyBase.cpp
  src/ClpCholeskyDense.cpp
  src/ClpConstraint.cpp
  src/ClpConstraintLinear.cpp
  src/ClpConstraintQuadratic.cpp
  src/ClpDualRowDantzig.cpp
  src/ClpDualRowPivot.cpp
  src/ClpDualRowSteepest.cpp
  src/ClpDummyMatrix.cpp
  src/ClpDynamicExampleMatrix.cpp
  src/ClpDynamicMatrix.cpp
  src/ClpEventHandler.cpp
  src/ClpFactorization.cpp
  src/ClpGubDynamicMatrix.cpp
  src/ClpGubMatrix.cpp
  src/ClpHelperFunctions.cpp
  src/ClpInterior.cpp
  src/ClpLinearObjective.cpp
  src/ClpLsqr.cpp
  src/ClpMatrixBase.cpp
  src/ClpMessage.cpp
  src/ClpModel.cpp
  src/ClpNetworkBasis.cpp
  src/ClpNetworkMatrix.cpp
  src/ClpNode.cpp
  src/ClpNonLinearCost.cpp
  src/ClpObjective.cpp
  src/ClpPEDualRowDantzig.cpp
  src/ClpPEDualRowSteepest.cpp
  src/ClpPEPrimalColumnDantzig.cpp
  src/ClpPEPrimalColumnSteepest.cpp
  src/ClpPESimplex.cpp
  src/ClpPackedMatrix.cpp
  src/ClpPdco.cpp
  src/ClpPdcoBase.cpp
  src/ClpPlusMinusOneMatrix.cpp
  src/ClpPredictorCorrector.cpp
  src/ClpPresolve.cpp
  src/ClpPrimalColumnDantzig.cpp
  src/ClpPrimalColumnPivot.cpp
  src/ClpPrimalColumnSteepest.cpp
  src/ClpQuadraticObjective.cpp
  src/ClpSimplex.cpp
  src/ClpSimplexDual.cpp
  src/ClpSimplexNonlinear.cpp
  src/ClpSimplexOther.cpp
  src/ClpSimplexPrimal.cpp
  src/ClpSolve.cpp
  src/Clp_C_Interface.cpp
  src/IdiSolve.cpp
  src/Idiot.cpp
  )

set(_HDRS
  src/ClpCholeskyBase.hpp
  src/ClpCholeskyDense.hpp
  src/ClpConfig.h
  src/ClpConstraint.hpp
  src/ClpConstraintLinear.hpp
  src/ClpConstraintQuadratic.hpp
  src/ClpDualRowDantzig.hpp
  src/ClpDualRowPivot.hpp
  src/ClpDualRowSteepest.hpp
  src/ClpDummyMatrix.hpp
  src/ClpDynamicExampleMatrix.hpp
  src/ClpDynamicMatrix.hpp
  src/ClpEventHandler.hpp
  src/ClpFactorization.hpp
  src/ClpGubDynamicMatrix.hpp
  src/ClpGubMatrix.hpp
  src/ClpInterior.hpp
  src/ClpLinearObjective.hpp
  src/ClpMatrixBase.hpp
  src/ClpMessage.hpp
  src/ClpModel.hpp
  src/ClpNetworkMatrix.hpp
  src/ClpNode.hpp
  src/ClpNonLinearCost.hpp
  src/ClpObjective.hpp
  src/ClpPackedMatrix.hpp
  src/ClpParameters.hpp
  src/ClpPdcoBase.hpp
  src/ClpPdco.hpp
  src/ClpPEDualRowDantzig.hpp
  src/ClpPEDualRowSteepest.hpp
  src/ClpPEPrimalColumnDantzig.hpp
  src/ClpPEPrimalColumnSteepest.hpp
  src/ClpPESimplex.hpp
  src/ClpPlusMinusOneMatrix.hpp
  src/ClpPresolve.hpp
  src/ClpPrimalColumnDantzig.hpp
  src/ClpPrimalColumnPivot.hpp
  src/ClpPrimalColumnSteepest.hpp
  src/ClpQuadraticObjective.hpp
  src/ClpSimplex.hpp
  src/ClpSimplexDual.hpp
  src/ClpSimplexNonlinear.hpp
  src/ClpSimplexOther.hpp
  src/ClpSimplexPrimal.hpp
  src/ClpSolve.hpp
  src/Clp_C_Interface.h
  src/Idiot.hpp
  )

include(GNUInstallDirs)

add_library(Clp ${_SRCS} ${_HDRS})
target_include_directories(Clp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
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
if(APPLE)
  set_target_properties(Clp PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(Clp PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(Clp PUBLIC
  Coin::CoinUtils
  Coin::Osi)
set_target_properties(Clp PROPERTIES
  PUBLIC_HEADER "${_HDRS};src/CbcOrClpParam.cpp;${CMAKE_CURRENT_BINARY_DIR}/config_clp.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::Clp ALIAS Clp)

# Install
install(TARGETS Clp
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

##############
##  OsiClp  ##
##############
set(_SRCS src/OsiClp/OsiClpSolverInterface.cpp)
set(_HDRS	src/OsiClp/OsiClpSolverInterface.hpp)

add_library(OsiClp ${_SRCS} ${_HDRS})
target_include_directories(OsiClp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/OsiClp>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(OsiClp
  PUBLIC HAVE_CONFIG_H
  PRIVATE CLP_BUILD)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET OsiClp PROPERTY CXX_STANDARD 11)
  set_property(TARGET OsiClp PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(OsiClp PUBLIC cxx_std_11)
endif()
if(APPLE)
  set_target_properties(OsiClp PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(OsiClp PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(OsiClp PUBLIC
  Coin::CoinUtils
  Coin::Osi
  Coin::Clp)
set_target_properties(OsiClp PROPERTIES
  PUBLIC_HEADER "${_HDRS}"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::OsiClp ALIAS OsiClp)

# Install
install(TARGETS OsiClp
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

#################
##  ClpSolver  ##
#################
set(_SRCS
  src/ClpSolver.cpp
  src/CbcOrClpParam.cpp
  src/Clp_ampl.cpp
  src/MyEventHandler.cpp
  src/MyMessageHandler.cpp)

set(_HDRS
  src/CbcOrClpParam.hpp
  src/Clp_ampl.h
  src/MyEventHandler.hpp
  src/MyMessageHandler.hpp
  src/unitTest.cpp)

add_library(ClpSolver ${_SRCS} ${_HDRS})
target_include_directories(ClpSolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
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
if(APPLE)
  set_target_properties(ClpSolver PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(ClpSolver PROPERTIES INSTALL_RPATH "$ORIGIN")
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
