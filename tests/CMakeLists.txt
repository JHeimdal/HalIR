
add_executable(test_input test_input.cpp)

target_compile_features(test_input PRIVATE cxx_std_17)

target_link_libraries(test_input PRIVATE HalIR)

enable_testing()
add_test(NAME test_exec COMMAND $<TARGET_FILE:halir>)
add_test(NAME input_units_test COMMAND $<TARGET_FILE:test_input>)
set_property(TEST test_exec PROPERTY WILL_FAIL TRUE)
