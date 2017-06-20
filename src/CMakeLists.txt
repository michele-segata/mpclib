# add additional include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# add TCLAP files to sources
file(GLOB tclap_sources tclap/*.h)

set(library_sources ${tclap_sources} Array.hh Array.cc
    QuadProg++.cc QuadProg++.hh Variable.cc Variable.h Constraint.cc
    Constraint.h matrix_utils.cc matrix_utils.h QuadraticProblem.cc
    QuadraticProblem.h MPCProblem.cc MPCProblem.h)

add_executable(test-mpc test-mpc.cc ${library_sources})