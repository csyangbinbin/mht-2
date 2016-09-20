#
# mht/src CMakeLists.txt
#

#Select C++11
set(CMAKE_CXX_COMPILER "g++-4.8")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 ")

#Enable warnings
add_definitions("-Wall")

#Include the directory itself
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#The targets included
set(SOURCES mainapp.cpp Student.cpp)

#Compile the code
add_executable(mht ${SOURCES})
target_link_libraries(mht emdw ${GLINEAR_LIBRARIES} ${PATRECII_LIBRARIES} ${BOOST_LIBRARIES})
