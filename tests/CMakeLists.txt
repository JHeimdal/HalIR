
add_executable(test_input test_input.cpp)
add_executable(test_hpar test_hpar.cpp)

target_compile_features(test_input PRIVATE cxx_std_17)

target_link_libraries(test_input PRIVATE HalIR)
target_link_libraries(test_hpar PRIVATE HalIR)

enable_testing()
add_test(NAME test_exec COMMAND $<TARGET_FILE:halir>)
set_property(TEST test_exec PROPERTY WILL_FAIL TRUE)
add_test(NAME input_units_test_F_TO_K COMMAND $<TARGET_FILE:test_input> 0)
add_test(NAME input_units_test_C_TO_K COMMAND $<TARGET_FILE:test_input> 1)
add_test(NAME input_units_test_BAR_TO_ATM COMMAND $<TARGET_FILE:test_input> 2)
add_test(NAME input_units_test_MBAR_TO_ATM COMMAND $<TARGET_FILE:test_input> 3)
add_test(NAME input_units_test_PA_TO_ATM COMMAND $<TARGET_FILE:test_input> 4)
add_test(NAME input_units_test_HPA_TO_ATM COMMAND $<TARGET_FILE:test_input> 5)
add_test(NAME input_units_test_MMHG_TO_ATM COMMAND $<TARGET_FILE:test_input> 6)

