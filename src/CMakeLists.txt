
set(HEADER_LIST "${HalIR_SOURCE_DIR}/include/HalIR/halir.h")
add_library(HalIR SHARED halir.c cJSON.c ${HEADER_LIST})

target_include_directories(HalIR PUBLIC ../include)
target_link_libraries(HalIR PUBLIC m gsl gslcblas cerf)
