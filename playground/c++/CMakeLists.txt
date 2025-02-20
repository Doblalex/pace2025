cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

project(Pace2025)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SOURCES
  src/main.cpp
  src/readinstance.cpp
)

add_executable(dsexact ${SOURCES})

target_include_directories(dsexact
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
)