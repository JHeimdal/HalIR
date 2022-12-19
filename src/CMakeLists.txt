
set(HEADER_LIST "${HalIR_SOURCE_DIR}/include/HalIR/halir.h")
add_library(HalIR SHARED halir.c cJSON.c ${HEADER_LIST})

target_include_directories(HalIR PUBLIC ../include)
