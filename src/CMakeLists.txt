# add additional include directories
include_directories(${CMAKE_SOURCE_DIR}/src)
if (${ENABLE_DISCRETIZATION})
    include_directories(${GSL_INCLUDE_DIRS})
endif()

# add TCLAP files to sources
file(GLOB tclap_sources tclap/*.h)

set(library_sources ${tclap_sources} Array.hh Array.cc
    QuadProg++.cc QuadProg++.hh Variable.cc Variable.h Constraint.cc
    Constraint.h matrix_utils.cc matrix_utils.h QuadraticProblem.cc
    QuadraticProblem.h MPCProblem.cc MPCProblem.h)

if (${ENABLE_DISCRETIZATION})
    set(LIBS ${LIBS} ${GSL_LIBRARIES})
endif()

if (${ENABLE_CONFIG_FILE})
    set(library_sources ${library_sources} ConfigLoader.cpp ConfigLoader.h
        MathEvaluator.cpp MathEvaluator.h)
    set(LIBS ${LIBS} ${MATHEVAL_LIBRARIES} ${LIBCONFIGPP_LIBRARIES})
endif()

add_executable(test-mpc test-mpc.cc ${library_sources})
target_link_libraries(test-mpc ${LIBS})