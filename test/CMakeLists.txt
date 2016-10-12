#
# mht/tests CMakeLists.txt
#

#Enable warnings
add_definitions("-Wall")

#Include the directory itself
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#The targets included
set(SOURCES run_tests.cc transforms_test.cc)

#Compile the code
add_executable(unit_tests ${SOURCES})
target_link_libraries(unit_tests emdw ${GLINEAR_LIBRARIES} ${PATRECII_LIBRARIES} ${BOOST_LIBRARIES} ${GTEST_LIBRARIES} pthread)
