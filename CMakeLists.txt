cmake_minimum_required(VERSION 3.15)

project(HalIR
        VERSION 0.0.1
        DESCRIPTION "A library for cfiting synthetic spektra to IR data"
        LANGUAGES C)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

add_subdirectory(src)

add_subdirectory(apps)

#if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING) AND BUILD_TESTING)
  #add_subdirectory(tests)
#endif()
