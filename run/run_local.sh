#!/bin/bash
# Runs the vertical remapping code on a local machine.

echo '-- Setting run variables...'

# Set directory path variables.
run_dir=`pwd`
build_dir=../build

for cells in {1..10..1}
do
    $build_dir/vert_remap ncell $((2**$cells))
done

echo '-- Runs complete!'
