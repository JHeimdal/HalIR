FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6
  )
FetchContent_MakeAvailable(catch)

add_executable(testapp testapp.cpp)

target_compile_features(testapp PRIVATE cxx_std_17)

target_link_libraries(testapp PRIVATE HalIR Catch2::Catch2)

enable_testing()
add_test(NAME testappTest COMMAND $<TARGET_FILE:testapp>)
