project(MPCLib CXX)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
set(CMAKE_INCLUDE_PATH  ${CMAKE_INCLUDE_PATH} /opt/local/include)
link_directories(/opt/local/lib)

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

FIND_PACKAGE(GSL)
FIND_PACKAGE(Matheval)
FIND_PACKAGE(Config)

if (${GSL_FOUND})
    set(ENABLE_DISCRETIZATION TRUE)
    add_definitions(-DENABLE_DISCRETIZATION)
    if (${LIBCONFIGPP_FOUND} AND ${MATHEVAL_FOUND})
        set(ENABLE_CONFIG_FILE TRUE)
        add_definitions(-DENABLE_CONFIG_FILE)
    endif()
endif()

add_subdirectory(src)
