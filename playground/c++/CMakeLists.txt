cmake_minimum_required(VERSION 3.18)

project(Pace2025)

option(CXX "enable C++ compilation" ON)
# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CXX)
    set (CMAKE_CXX_STANDARD 17)
    enable_language(CXX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SOURCES
  src/main.cpp
  src/readinstance.cpp
  src/preprocessing.cpp
  src/partition.cpp
  src/heuristics.cpp
  src/exactsolver.cpp
)

add_executable(dsexact ${SOURCES})
add_executable(ogdf_dsexact src/ogdf_main.cpp)

# EVALMAXSAT ---------------------------------------------------------------
set(localLibs
	MaLib
    cadical
    EvalMaxSAT
)

set(externLibs
    z
)

FOREACH(localLib ${localLibs})
    if(NOT TARGET ${localLib})
        add_subdirectory(${PROJECT_SOURCE_DIR}/ext/EvalMaxSAT/lib/${localLib})
    endif()
ENDFOREACH(localLib)


find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/EvalMaxSAT)
target_link_libraries(ogdf_dsexact ${localLibs} ${externLibs})
include_directories(${PROJECT_SOURCE_DIR}/include)
# GUROBI -------------------------------------------------------------------
set(USE_GUROBI OFF)

if(USE_GUROBI)
  add_compile_definitions(USE_GUROBI)
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
  find_package(GUROBI REQUIRED)

  include_directories(${GUROBI_INCLUDE_DIRS})
  if(CXX)
    target_link_libraries(dsexact optimized ${GUROBI_CXX_LIBRARY}
      debug ${GUROBI_CXX_DEBUG_LIBRARY})
  endif()
  target_link_libraries(dsexact ${GUROBI_LIBRARY})
endif()

# OR-TOOLS -----------------------------------------------------------------
set(USE_ORTOOLS OFF)


if(USE_ORTOOLS)
  find_package(ortools)
  add_compile_definitions(USE_ORTOOLS)
  IF(${ortools_FOUND})
    message(STATUS "ortools found")
  ELSE()
    set(BUILD_DEPS ON)
    set(BUILD_DOC OFF)
    set(BUILD_FLATZINC OFF)
    set(BUILD_SAMPLES OFF)
    set(BUILD_EXAMPLES OFF)
    set(USE_HIGHS OFF)
    set(USE_SCIP OFF)
    set(USE_COINOR OFF)
    set(USE_GLPK OFF)
    set(USE_CPLEX OFF)
    set(USE_DOTNET_6 OFF)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/or-tools/cmake")
    add_subdirectory("${PROJECT_SOURCE_DIR}/ext/or-tools")

    include_directories(${PROJECT_SOURCE_DIR}/ext/or-tools/ortools)
  endif()
  target_link_libraries(dsexact ortools::ortools)
endif()


# OGDF -----------------------------------------------------------------
set(USE_OGDF ON)


if(USE_OGDF)
  find_package(OGDF)
  add_compile_definitions(USE_OGDF)
  IF(${OGDF_FOUND})
    message(STATUS "OGDF found")
  ELSE()
    set(OGDF_MEMORY_MANAGER POOL_NTS)
    # CMAKE_BUILD_TYPE=Release
    # BUILD_SHARED_LIBS=OFF
    # CMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE
    # CMAKE_POLICY_DEFAULT_CMP0069=NEW
    # OGDF_USE_ASSERT_EXCEPTIONS=OFF
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/ext/ogdf/cmake")
    add_subdirectory("${PROJECT_SOURCE_DIR}/ext/ogdf")

    include_directories(${PROJECT_SOURCE_DIR}/ext/ogdf/include)
  endif()
  target_link_libraries(dsexact OGDF)
  target_link_libraries(ogdf_dsexact OGDF)
endif()
