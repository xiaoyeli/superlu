#!/bin/sh -e

scriptname=$(basename $0)
testname="${scriptname%.*}"

echo ""
echo "Running test: $testname ..." 

rm -r SRC TESTING CMake* || echo "no files to delete"

env CC=gcc \
  cmake -DUSE_XSDK_DEFAULTS=TRUE \
  -DBUILD_SHARED_LIBS=FALSE \
  -Denable_blaslib=TRUE \
  -DCMAKE_INSTALL_PREFIX=$PWD/../install \
  ../../superlu &> configure.$testname.out || echo "FAILED configure"

make -j12 &> make.$testname.out || echo "FAILED build"

ctest -j12 &> ctest.$testname.out || echo "FAILED tests"
