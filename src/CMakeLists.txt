#
# mht/src CMakeLists.txt
#

#Enable warnings
add_definitions("-Wall")

#Include the directory itself
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#The targets included
set(SOURCES main.cc)

#Compile the code
add_executable(mht ${SOURCES})
target_link_libraries(mht emdw ${GLINEAR_LIBRARIES} ${PATRECII_LIBRARIES} ${BOOST_LIBRARIES})
