#!/usr/bin/env bash

echo "Profiling build:"
echo "================"
stack --work-dir .stack-work-profiling --profile build

echo ""

echo "Profiled execution:"
echo "==================="
time stack exec -- chemostat +RTS -p
