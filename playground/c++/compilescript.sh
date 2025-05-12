#!/usr/bin/env bash
#qsub -N test -pe pthreads 20 -e ~/error.log compilescript.sh 
source /home1/adobler/.bashrc
cd /home1/adobler/pace2025/playground/c++/build-release
cmake -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=gold -DPACE_LOG=1 ..
make -j 20 > makeout.out