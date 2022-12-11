cmake_minimum_required(VERSION 3.8)

project(HalIR LANGUAGES C)

add_library(HalIR SHARED src/halir.c include/HalIR/halir.h)
target_include_directories(HalIR PUBLIC include)

add_executable(halir apps/halir.c)
target_link_libraries(halir PUBLIC HalIR)
