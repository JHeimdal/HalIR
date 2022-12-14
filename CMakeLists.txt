cmake_minimum_required(VERSION 3.20)

project(HalIR
        VERSION 0.0.1
        DESCRIPTION "A library for cfiting synthetic spektra to IR data"
        LANGUAGES C CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_EXTENSIONS OFF)
  include(CTest)
endif()

include(FetchContent)

find_package(GSL)
add_subdirectory(src)

add_subdirectory(apps)

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
