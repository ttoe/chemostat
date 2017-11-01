#!/usr/bin/env bash

echo "Emptying plots/"
rm -R ./plots
mkdir ./plots

echo Running "$1"

time (stack runghc "$1")

unamestr=$(uname)
if [[ "$unamestr" == 'Darwin' ]]; then
  open ./plots/*
fi
