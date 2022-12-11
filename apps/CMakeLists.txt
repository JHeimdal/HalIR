add_executable(halir halir.c cJSON.c)

target_link_libraries(halir PUBLIC HalIR)
