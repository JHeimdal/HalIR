add_executable(halir halir.c)
add_executable(halir_cls halir_cls.c)

#target_link_directories(halir PUBLIC HalIR)
set(CMAKE_BUILD_TYPE Debug)
target_link_libraries(halir PUBLIC HalIR)
target_link_libraries(halir_cls PUBLIC HalIR)