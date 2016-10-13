#
# mht/src CMakeLists.txt
#

#Enable warnings
add_definitions("-Wall")

#Include the headers and src folders
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.hpp")
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

#Create a library for unit testing
add_library(mht_base ${SOURCES} ${HEADERS})

#Compile the code
add_executable(mht main.cc)
target_link_libraries(mht mht_base emdw ${GLINEAR_LIBRARIES} ${PATRECII_LIBRARIES} ${BOOST_LIBRARIES})
