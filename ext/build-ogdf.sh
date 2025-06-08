#!/bin/bash

git submodule update --init --recursive
cd ogdf
mkdir -p build-debug build-release

cmake -S . -B build-release \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE -DCMAKE_POLICY_DEFAULT_CMP0069=NEW \
    -DOGDF_MEMORY_MANAGER=POOL_NTS -DOGDF_USE_ASSERT_EXCEPTIONS=OFF # -DOGDF_ARCH=haswell
cmake --build build-release -j $(nproc)

cmake -S . -B build-debug \
    -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON \
    -DOGDF_MEMORY_MANAGER=MALLOC_TS -DOGDF_LEAK_CHECK=ON \
    -DOGDF_USE_ASSERT_EXCEPTIONS=ON -DOGDF_USE_ASSERT_EXCEPTIONS_WITH_STACK_TRACE=ON_LIBUNWIND
cmake --build build-debug -j $(nproc)
