
enable_language(Fortran)
#add_definitions(-ffixed-line-length-132)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132")
set(HEADER_LIST "${HalIR_SOURCE_DIR}/include/HalIR/halir.h")
add_library(HalIR SHARED halir.c cJSON.c BD_TIPS_2021.FOR ${HEADER_LIST})

target_include_directories(HalIR PUBLIC ../include)
target_link_libraries(HalIR PUBLIC m gsl gslcblas cerf)
