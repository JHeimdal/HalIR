#set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost 1.58.0 COMPONENTS system filesystem REQUIRED)
find_package(Armadillo)
find_package(OpenMP)
if (OPENMP_FOUND)
    #set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
include_directories( ${HalIR_SOURCE_DIR/SpecSac} ${Boost_INCLUDE_DIR} ${HalIR_SOURCE_DIR/HalIR} )
set( SRCS test.cpp )
ADD_DEFINITIONS( -std=c++11 )
link_directories( ${HalIR_BINARY_DIR/SpecSac} ${HalIR_BINARY_DIR/HalIR} ${Boost_LIBRARY_DIR} )

add_executable(mytest ${SRCS})
message(STATUS ${ARMADILLO_LIBRARIES})
target_link_libraries(mytest ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} SpecSac HalIR )


install( TARGETS mytest DESTINATION .bin )
install( DIRECTORY TestDat DESTINATION ${PROJECT_BINARY_DIR}/Test )
install( PROGRAMS genMtest.py halirtest.py run_halir_test.py spc.py plotf.py
         DESTINATION ${PROJECT_BINARY_DIR}/Test )
#message( ${PROJECT_SOURCE_DIR} )
#message( ${PROJECT_BINARY_DIR} )
