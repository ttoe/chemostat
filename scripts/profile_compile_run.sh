#!/usr/bin/env bash

echo "Profiling build:"
echo "================"
stack --work-dir .stack-work-profiling --profile build

echo ""

echo "Profiled execution:"
echo "==================="
time .stack-work-profiling/install/x86_64-osx/nightly-2018-06-29/8.4.3/bin/chemostat +RTS -p
#time stack exec -- chemostat +RTS -p
