cmake_minimum_required(VERSION 3.18)

project(Pace2025)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()

set(PACE_USE_GUROBI OFF CACHE BOOL "use installed Gurobi")
set(PACE_USE_EVALMAXSAT ON CACHE BOOL "use installed or subproject-built EvalMaxSat")
set(PACE_USE_UWRMAXSAT OFF CACHE BOOL "use installed UWrMaxSat")
set(PACE_USE_ORTOOLS OFF CACHE BOOL "use installed or subproject-built ortools")

set(PACE_USE_ASAN OFF CACHE BOOL "use the address sanitizer")
set(PACE_ARCH "native" CACHE STRING "Target CPU (micro)architecture passed to the compiler via `-march`.")

set(PACE_LOG OFF CACHE BOOL "enable logging")
set(PACE_SAT_CACHE OFF CACHE BOOL "enable SAT solution caching")

if(PACE_LOG)
    add_compile_definitions(PACE_LOG)
endif()

add_compile_options("-march=${PACE_ARCH}")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in Release mode")

    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON CACHE BOOL "INTERPROCEDURAL_OPTIMIZATION is enabled by default for release builds")
    set(CMAKE_POLICY_DEFAULT_CMP0069 NEW CACHE STRING "INTERPROCEDURAL_OPTIMIZATION should no complain of OGDF's older CMake version")

    if(NOT BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Shared Libs OFF")
    endif()
else()
    message(STATUS "Building in non-Release mode ${CMAKE_BUILD_TYPE}")

    if(NOT BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS ON CACHE BOOL "Shared Libs ON")
    endif()
endif()

if(PACE_USE_ASAN)
    add_compile_options(-g -fno-omit-frame-pointer -fsanitize=address)
    add_link_options(-fno-omit-frame-pointer -fsanitize=address)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/ext)
add_executable(ogdf_dsexact
    src/ogdf_instance.cpp
    src/ogdf_io.cpp
    src/ogdf_main.cpp
    src/ogdf_solver.cpp
    src/ogdf_subsetrefine.cpp
    src/ogdf_treewidth.cpp
    src/ogdf_vcreductions.cpp
    src/ogdf_solver/sol_greedy.cpp)
add_executable(ogdf_validate src/ogdf_validator.cpp src/ogdf_io.cpp)

if(PACE_SAT_CACHE)
    add_compile_definitions(SAT_CACHE)
    target_sources(ogdf_dsexact PUBLIC src/ogdf_solver/satcache.cpp)
endif()

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIR})

# EVALMAXSAT ---------------------------------------------------------------
if(PACE_USE_EVALMAXSAT)
    add_compile_definitions(USE_EVALMAXSAT)
    target_sources(ogdf_dsexact PUBLIC src/ogdf_solver/sol_evalmaxsat.cpp)
    set(localLibs
        MaLib
        cadical
        EvalMaxSAT
        glucose
    )

    foreach(localLib ${localLibs})
        if(NOT TARGET ${localLib})
            add_subdirectory(${PROJECT_SOURCE_DIR}/ext/EvalMaxSat/lib/${localLib})
        endif()
    endforeach(localLib)

    target_link_libraries(ogdf_dsexact ${localLibs} z)
endif()

# UWrMaxSat ----------------------------------------------------------------
if(PACE_USE_UWRMAXSAT)
    add_compile_definitions(USE_UWRMAXSAT)
    target_sources(ogdf_dsexact PUBLIC src/ogdf_solver/sol_uwrmaxsat.cpp)
    find_path(UWRMAXSAT_INCLUDE_DIR "uwrmaxsat/ipamir.h" HINTS "${PROJECT_SOURCE_DIR}/ext")
    find_library(UWRMAXSAT_LIBRARY uwrmaxsat HINTS "${PROJECT_SOURCE_DIR}/ext/uwrmaxsat/build/release/lib")
    find_library(LMAXPRE_LIBRARY maxpre HINTS "${PROJECT_SOURCE_DIR}/ext/maxpre/src/lib")
    find_library(SCIP_LIBRARY scip HINTS "${PROJECT_SOURCE_DIR}/ext/scipoptsuite/build/lib")
    find_library(LSOPLEX_LIBRARY soplex-pic HINTS "${PROJECT_SOURCE_DIR}/ext/scipoptsuite/build/lib")
    find_library(LCADICAL_LIBRARY cadical HINTS "${PROJECT_SOURCE_DIR}/ext/cadical/build")

    find_package_handle_standard_args(UWrMaxSat DEFAULT_MSG UWRMAXSAT_INCLUDE_DIR UWRMAXSAT_LIBRARY
        LMAXPRE_LIBRARY SCIP_LIBRARY LSOPLEX_LIBRARY LCADICAL_LIBRARY)

    add_library(UWrMaxSat UNKNOWN IMPORTED)
    set_target_properties(UWrMaxSat PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UWRMAXSAT_INCLUDE_DIR}"
        IMPORTED_LOCATION "${UWRMAXSAT_LIBRARY}")

    add_library(LMAXPRE UNKNOWN IMPORTED)
    set_target_properties(LMAXPRE PROPERTIES IMPORTED_LOCATION "${LMAXPRE_LIBRARY}")
    add_library(LSCIP UNKNOWN IMPORTED)
    set_target_properties(LSCIP PROPERTIES IMPORTED_LOCATION "${SCIP_LIBRARY}")
    add_library(LSOPLEX UNKNOWN IMPORTED)
    set_target_properties(LSOPLEX PROPERTIES IMPORTED_LOCATION "${LSOPLEX_LIBRARY}")
    add_library(LCADICAL UNKNOWN IMPORTED)
    set_target_properties(LCADICAL PROPERTIES IMPORTED_LOCATION "${LCADICAL_LIBRARY}")

    # ar -rcs build/release/lib/libuwrmaxsat.a build/release/Debug.o build/release/Hardware_adders.o build/release/Hardware_clausify.o build/release/Hardware_reuse_sorters.o build/release/Hardware_sorters.o build/release/ipamir.o build/release/Main_utils.o build/release/MsSolver.o build/release/PbParser.o build/release/PbSolver.o build/release/PbSolver_convertAdd.o build/release/PbSolver_convertBdd.o build/release/PbSolver_convert.o build/release/PbSolver_convertSort.o build/release/Pre_separator.o build/release/SatSolver.o build/release/ScipSolver.o build/release/System.o build/release/ADTs/FEnv.o build/release/ADTs/File.o build/release/ADTs/Global.o
    # g++ build/release/Main.o build/release/lib/libuwrmaxsat.a -Wall  -L../maxpre/src/lib -lmaxpre -L../scipoptsuite/build/lib -lscip -lsoplex-pic -L/lib -L../cadical/build -lcadical -D USE_SCIP -pthread -lz -lgmp --static   -o build/release/bin/uwrmaxsat
    target_link_libraries(UWrMaxSat INTERFACE gmp z LMAXPRE LSCIP LSOPLEX LCADICAL)
    target_link_options(UWrMaxSat INTERFACE "-D USE_SCIP")

    target_link_libraries(ogdf_dsexact UWrMaxSat)
endif()

# GUROBI -------------------------------------------------------------------
if(PACE_USE_GUROBI)
    add_compile_definitions(USE_GUROBI)
    target_sources(ogdf_dsexact PUBLIC src/ogdf_solver/sol_gurobi.cpp)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
    find_package(GUROBI REQUIRED)

    include_directories(${GUROBI_INCLUDE_DIRS})
    target_link_libraries(ogdf_dsexact ${GUROBI_LIBRARY})
endif()

# OR-TOOLS -----------------------------------------------------------------
if(PACE_USE_ORTOOLS)
    find_package(ortools)
    add_compile_definitions(USE_ORTOOLS)
    target_sources(ogdf_dsexact PUBLIC src/ogdf_solver/sol_ortools.cpp)

    if(${ortools_FOUND})
        message(STATUS "ortools found")
    else()
        set(BUILD_DEPS ON CACHE BOOL "ortools-option")
        set(BUILD_DOC OFF CACHE BOOL "ortools-option")
        set(BUILD_FLATZINC OFF CACHE BOOL "ortools-option")
        set(BUILD_SAMPLES OFF CACHE BOOL "ortools-option")
        set(BUILD_EXAMPLES OFF CACHE BOOL "ortools-option")
        set(USE_HIGHS OFF CACHE BOOL "ortools-option")
        set(USE_SCIP OFF CACHE BOOL "ortools-option")
        set(USE_COINOR OFF CACHE BOOL "ortools-option")
        set(USE_GLPK OFF CACHE BOOL "ortools-option")
        set(USE_CPLEX OFF CACHE BOOL "ortools-option")
        set(USE_DOTNET_6 OFF CACHE BOOL "ortools-option")

        list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/or-tools/cmake")
        add_subdirectory("${PROJECT_SOURCE_DIR}/ext/or-tools")
        include_directories(${PROJECT_SOURCE_DIR}/ext/or-tools/ortools)
    endif()

    target_link_libraries(ogdf_dsexact ortools::ortools)
endif()

# HTD ----------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR}/ext/htd/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/htd)

target_link_libraries(ogdf_dsexact htd)

# OGDF -----------------------------------------------------------------
find_package(OGDF)
if (${OGDF_FOUND})
    message(STATUS "OGDF found")
else ()
    set(OGDF_ARCH "native" CACHE STRING "OGDF_ARCH defaults to haswell")
    set(OGDF_MEMORY_MANAGER POOL_NTS CACHE STRING "OGDF_MEMORY_MANAGER defaults to POOL_NTS")

    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        set(OGDF_USE_ASSERT_EXCEPTIONS OFF CACHE BOOL "OGDF_USE_ASSERT_EXCEPTIONS defaults to OFF")
    else ()
        set(OGDF_USE_ASSERT_EXCEPTIONS ON CACHE BOOL "OGDF_USE_ASSERT_EXCEPTIONS defaults to ON")
    endif ()

    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/ogdf/cmake")
    add_subdirectory("${PROJECT_SOURCE_DIR}/ext/ogdf")
    include_directories(${PROJECT_SOURCE_DIR}/ext/ogdf/include)
endif ()

target_link_libraries(ogdf_dsexact OGDF)
target_link_libraries(ogdf_validate OGDF)
