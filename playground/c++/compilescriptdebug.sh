#!/usr/bin/env bash
#qsub -N test -pe pthreads 20 -e ~/error.log compilescriptdebug.sh 
source /home1/adobler/.bashrc
cd /home1/adobler/pace2025/playground/c++/build-debug
cmake -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=gold -DPACE_LOG=1 -DCMAKE_BUILD_TYPE=Debug ..
make -j 20 > makeout.out