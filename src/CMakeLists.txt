
enable_language(Fortran)
#add_definitions(-ffixed-line-length-132)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132")
set(HEADER_LIST "${HalIR_SOURCE_DIR}/include/HalIR/halir.h"
                "${HalIR_SOURCE_DIR}/include/HalIR/halir_cls.h")
add_library(HalIR SHARED halir_calculate.c halir_workspace.c halir_result.c halir_parser.c halir_spc.c halir_cls.c cJSON.c BD_TIPS_2021.FOR ${HEADER_LIST})

target_include_directories(HalIR PUBLIC ../include)
target_link_libraries(HalIR PUBLIC m gsl gslcblas cerf)

if(OpenMP_C_FOUND)
  target_link_libraries(HalIR PUBLIC OpenMP::OpenMP_C)
  target_compile_definitions(HalIR PRIVATE HALIR_HAVE_OPENMP)
endif()
