cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set build type: Release or Debug" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(STATUS "Configuring project in debug mode. Set -DCMAKE_BUILD_TYPE=Release for release mode")
elseif(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "Configuring project in release mode. Set -DCMAKE_BUILD_TYPE=Debug for debug mode")
else(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "Configuring project in " ${CMAKE_BUILD_TYPE} " mode")
endif()

project(MPCLib CXX)
add_subdirectory(src)
