#!/usr/bin/env bash

echo "Compilation:"
echo "============"
time (stack build)

echo ""

echo "Execution:"
echo "=========="
time stack exec diffeq-hs
