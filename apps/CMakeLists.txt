add_executable(halir halir.c)

#target_link_directories(halir PUBLIC HalIR)
target_link_libraries(halir PUBLIC HalIR)
