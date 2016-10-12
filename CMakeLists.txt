#
# mht CMakeLists.txt
#

#Enforce CMake version
cmake_minimum_required(VERSION 2.6)

#Declare the project title
if("${PROJECT_TITLE}" STREQUAL "")
  set(PROJECT_TITLE mht)
endif()
project(${PROJECT_TITLE})

#Select C++11
set(CMAKE_CXX_COMPILER "g++-4.8")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 ")

#Status updates
message(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "Build  directory: $ENV{PWD}")
message(STATUS "Configuring mht ...")

#Locate CMake files to find directories
if("${CMAKE_MODULE_PATH}" STREQUAL "")
  set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
endif()

#The macro below forces the build directory to be different from source directory:
include(MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")

#Locate Boost
find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

#Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

#Locate gLinear
find_package(gLinear REQUIRED)
include_directories(${GLINEAR_INCLUDE_DIRS})

#Locate patrecII
find_package(PatrecII REQUIRED)
include_directories(${PATRECII_INCLUDE_DIRS})

#Locate EMDW - Cheapish relative path (Findemdw.cmake approach doesn't work yet)
set(EMDW_INCLUDE_DIRS ~/devel/emdw/src)
include_directories(${EMDW_INCLUDE_DIRS})

#Include the headers folder
include_directories(${CMAKE_SOURCE_DIR}/include)

#Add subdirectories 
add_subdirectory(src)
add_subdirectory(test)
