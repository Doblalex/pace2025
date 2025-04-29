cmake_minimum_required(VERSION 3.18)

project(Pace2025)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build Type" FORCE)
endif()

set(PACE_USE_GUROBI OFF CACHE BOOL "use installed Gurobi")
set(PACE_USE_ORTOOLS OFF CACHE BOOL "use installed or subproject-built ortools")
set(PACE_USE_OGDF ON CACHE BOOL "use installed or subproject-built OGDF")
set(PACE_USE_ASAN OFF CACHE BOOL "use the address sanitizer")
set(PACE_ARCH "native" CACHE STRING "Target CPU (micro)architecture passed to the compiler via `-march`.")

add_compile_options("-march=${PACE_ARCH}")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Building in Release mode")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON CACHE BOOL "INTERPROCEDURAL_OPTIMIZATION is enabled by default for release builds")
  set(CMAKE_POLICY_DEFAULT_CMP0069 NEW CACHE STRING "INTERPROCEDURAL_OPTIMIZATION should no complain of OGDF's older CMake version")
else()
  message(STATUS "Building in non-Release mode ${CMAKE_BUILD_TYPE}")
endif()

if(PACE_USE_ASAN)
  add_compile_options(-g -fno-omit-frame-pointer -fsanitize=address)
  add_link_options(-fno-omit-frame-pointer -fsanitize=address)
endif()

set(SOURCES
  src/main.cpp
  src/readinstance.cpp
  src/preprocessing.cpp
  src/partition.cpp
  src/heuristics.cpp
  src/exactsolver.cpp
)

set(OGDF_SOURCES
  src/ogdf_solver.cpp
  src/ogdf_greedy.cpp
  src/ogdf_instance.cpp
  src/ogdf_io.cpp
  src/ogdf_maxsat.cpp
  src/ogdf_subsetrefine.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(dsexact ${SOURCES})
add_executable(ogdf_dsexact ${OGDF_SOURCES} src/ogdf_main.cpp)
add_executable(ogdf_validate src/ogdf_validator.cpp src/ogdf_io.cpp)
add_executable(ogdf_hsexact ${OGDF_SOURCES} src/ogdf_mainhs.cpp)

# EVALMAXSAT ---------------------------------------------------------------
set(localLibs
  MaLib
  cadical
  EvalMaxSAT
)

FOREACH(localLib ${localLibs})
  if(NOT TARGET ${localLib})
    add_subdirectory(${PROJECT_SOURCE_DIR}/ext/EvalMaxSAT/lib/${localLib})
  endif()
ENDFOREACH(localLib)

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/EvalMaxSAT)
target_link_libraries(ogdf_dsexact ${localLibs} z)
target_link_libraries(ogdf_hsexact ${localLibs} z)
target_link_libraries(dsexact ${localLibs} z)

# GUROBI -------------------------------------------------------------------
if(PACE_USE_GUROBI)
  add_compile_definitions(USE_GUROBI)
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
  find_package(GUROBI REQUIRED)

  include_directories(${GUROBI_INCLUDE_DIRS})
  target_link_libraries(dsexact ${GUROBI_LIBRARY})

  # target_link_libraries(ogdf_dsexact ${GUROBI_LIBRARY})
endif()

# OR-TOOLS -----------------------------------------------------------------
if(PACE_USE_ORTOOLS)
  find_package(ortools)
  add_compile_definitions(USE_ORTOOLS)

  IF(${ortools_FOUND})
    message(STATUS "ortools found")
  ELSE()
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

  target_link_libraries(dsexact ortools::ortools)
  target_link_libraries(ogdf_hsexact ortools::ortools)
  target_link_libraries(ogdf_dsexact ortools::ortools)
endif()

# OGDF -----------------------------------------------------------------
if(PACE_USE_OGDF)
  find_package(OGDF)
  add_compile_definitions(USE_OGDF)

  IF(${OGDF_FOUND})
    message(STATUS "OGDF found")
  ELSE()
    set(OGDF_MEMORY_MANAGER POOL_NTS CACHE STRING "OGDF_MEMORY_MANAGER defaults to POOL_NTS")

    if(CMAKE_BUILD_TYPE STREQUAL "Release")
      set(OGDF_USE_ASSERT_EXCEPTIONS OFF CACHE BOOL "OGDF_USE_ASSERT_EXCEPTIONS defaults to OFF")
    else()
      set(OGDF_USE_ASSERT_EXCEPTIONS ON CACHE BOOL "OGDF_USE_ASSERT_EXCEPTIONS defaults to ON")
    endif()

    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/ogdf/cmake")
    add_subdirectory("${PROJECT_SOURCE_DIR}/ext/ogdf")
    include_directories(${PROJECT_SOURCE_DIR}/ext/ogdf/include)
  endif()

  target_link_libraries(ogdf_dsexact OGDF)
  target_link_libraries(ogdf_hsexact OGDF)
  target_link_libraries(ogdf_validate OGDF)
endif()
