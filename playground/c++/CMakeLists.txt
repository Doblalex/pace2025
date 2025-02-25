cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_COMPILER g++)

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
target_link_libraries(dsexact ${localLibs} ${externLibs})

# GUROBI -------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
find_package(GUROBI REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${GUROBI_INCLUDE_DIRS})
if(CXX)
  target_link_libraries(dsexact optimized ${GUROBI_CXX_LIBRARY}
    debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()
target_link_libraries(dsexact ${GUROBI_LIBRARY})

# OR-TOOLS -----------------------------------------------------------------
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

target_link_libraries(dsexact ortools::ortools)

# SCIP ---------------------------------------------------------------------

# set(IPOPT off)
# set(WORHP off)
# set(ZIMPL off)
# set(READLINE off)
# set(SOPLEX_DIR ${CMAKE_BINARY_DIR})
# set(SCIP_DIR ${CMAKE_BINARY_DIR})
# set(PAPILO off)
# set(SHARED off)
# add_subdirectory(${PROJECT_SOURCE_DIR}/ext/soplex)
# add_subdirectory(${PROJECT_SOURCE_DIR}/ext/scip)
# find_package(SCIP REQUIRED)

# include_directories(${SCIP_INCLUDE_DIRS})

# find_library(LIBM m)

# if(NOT LIBM)
#   set(LIBM "")
# endif()
# target_link_libraries(dsexact ${SCIP_LIBRARIES} ${LIBM})
# ---------------------------------------------------------------------











