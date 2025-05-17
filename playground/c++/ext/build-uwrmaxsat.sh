#!/bin/bash

# get sources of scipoptsuite from https://scipopt.org/index.php#download
SCIP_VER=9.2.2

set -x
set -e

# clone the repositories into uwrmaxsat and cominisatps:
git clone https://github.com/marekpiotrow/UWrMaxSat uwrmaxsat
git clone https://github.com/marekpiotrow/cominisatps
pushd cominisatps
  rm core simp mtl utils
  ln -s minisat/core minisat/simp minisat/mtl minisat/utils .
popd

# clone and build the CaDiCaL SAT solver by Armin Biere:
git clone https://github.com/arminbiere/cadical
pushd cadical
  patch -p1 <../uwrmaxsat/cadical.patch
  ./configure --no-contracts --no-tracing
  make cadical
popd
pushd uwrmaxsat
  cp config.cadical config.mk
popd

# build the MaxPre preprocessor (if you want to use it - see Comments below):
git clone https://github.com/Laakeri/maxpre
pushd maxpre
  sed -i 's/-g/-D NDEBUG/' src/Makefile
  make lib
popd

# build the SCIP solver library (if you want to use it)
curl https://scipopt.org/download/release/scipoptsuite-${SCIP_VER}.tgz -o scipoptsuite.tgz
tar zxvf scipoptsuite.tgz
mv scipoptsuite-${SCIP_VER} scipoptsuite
pushd scipoptsuite
  mkdir build
  pushd build
    cmake -DSYM=nauty -DSHARED=off -DNO_EXTERNAL_CODE=on -DSOPLEX=on -DTPI=tny ..
    cmake --build . --config Release --target libscip libsoplex-pic
  popd
popd

# build the UWrMaxSat solver (release version, statically linked):
pushd uwrmaxsat
  grep -rl scipoptsuite-9.2.1 . | xargs -r -l1 -- sed -i 's/scipoptsuite-9.2.1/scipoptsuite/' {} ;
  make clean
  make r
popd
