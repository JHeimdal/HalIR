add_executable(halir halir.c)

#target_link_directories(halir PUBLIC HalIR)
set(CMAKE_BUILD_TYPE Debug)
target_link_libraries(halir PUBLIC HalIR)
