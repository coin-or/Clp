set(NAME "OSICLP")

configure_file(config_osiclp.h.cmake.in config_osiclp.h)

set(_SRCS OsiClpSolverInterface.cpp)
set(_HDRS OsiClpSolverInterface.hpp)

add_library(OsiClp ${_SRCS} ${_HDRS})
target_include_directories(OsiClp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
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
target_link_libraries(OsiClp PUBLIC
  Coin::CoinUtils
  Coin::Osi
  Coin::Clp)
set_target_properties(OsiClp PROPERTIES
  PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_osiclp.h"
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
