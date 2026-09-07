
add_executable(ClpTest "ClpTest.cpp")
target_include_directories(Clp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(ClpTest
  PUBLIC HAVE_CONFIG_H
  PRIVATE CLP_BUILD COIN_HAS_CLP)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET ClpTest PROPERTY CXX_STANDARD 11)
  set_property(TARGET ClpTest PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(ClpTest PUBLIC cxx_std_11)
endif()
if(APPLE)
  set_target_properties(ClpTest PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
elseif(UNIX)
  cmake_path(RELATIVE_PATH CMAKE_INSTALL_FULL_LIBDIR
    BASE_DIRECTORY ${CMAKE_INSTALL_FULL_BINDIR}
    OUTPUT_VARIABLE libdir_relative_path)
  set_target_properties(ClpTest PROPERTIES
    INSTALL_RPATH "$ORIGIN/${libdir_relative_path}:$ORIGIN")
endif()
target_link_libraries(ClpTest PRIVATE
  Coin::CoinUtils
  Coin::Osi
  Coin::Clp)
set_target_properties(ClpTest PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_executable(Coin::ClpTest ALIAS ClpTest)

add_test(NAME ClpTest COMMAND ClpTest)
