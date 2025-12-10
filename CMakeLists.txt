# Build options
option(BUILD_TESTS "Build C++ tests with GoogleTest." ON)

# CMake configuration
cmake_minimum_required(VERSION 3.20)

# Project information
project(mdgm VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Source Library
file(GLOB 
  SRC_FILES 
  src/graph_storage.cpp
  src/undirected_graph.cpp
)
add_library(mdgm_objs OBJECT ${SRC_FILES})
target_include_directories(mdgm_objs PUBLIC 
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/inst/include
)

# Testing with GoogleTest
if(BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.17.0
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE) # for MSVC
  FetchContent_MakeAvailable(googletest)

  enable_testing()

  file(GLOB TEST_FILES tests/cpp/*.cpp)
  add_executable(mdgm_tests ${TEST_FILES})

  target_link_libraries(mdgm_tests PRIVATE
    mdgm_objs
    GTest::gtest_main
  )

  include(GoogleTest)
  gtest_discover_tests(mdgm_tests)
endif()

