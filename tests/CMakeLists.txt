
#add_executable(testapp testapp.cpp)

#target_compile_features(testapp PRIVATE cxx_std_17)

#target_link_libraries(testapp PRIVATE HalIR)

enable_testing()
add_test(NAME testappTest COMMAND $<TARGET_FILE:halir> input.json)
