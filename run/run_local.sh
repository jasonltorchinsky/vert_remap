#!/bin/bash
# Runs the vertical remapping code on a local machine.

echo '-- Setting run variables...'

# Set directory path variables.
run_dir=`pwd`
build_dir=../build

for cells in {8..16..1}
do
    for ogrid in sin
    do
	for tfunc in sig exp stp
	do
	    $build_dir/vert_remap ncell $((2**$cells)) ogrid $ogrid tfunc $tfunc
	done
    done
done

echo '-- Runs complete!'
