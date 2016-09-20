#
# mht/tests CMakeLists.txt
#

#Enable warnings
add_definitions("-Wall")

#The targets included
set(SOURCES tests.cpp)

#Compile the code
add_executable(unitTests ${SOURCES})
target_link_libraries(unitTests ${GTEST_LIBRARIES} pthread)
