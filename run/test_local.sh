#!/bin/bash
# Runs the vertical remapping code on a local machine.

echo '-- Setting run variables...'

# Set directory path variables.
run_dir=`pwd`
build_dir=../build

# Re-build the application
cd "$build_dir"
rm -rf *
cmake ..
cmake --build .
cd "$run_dir"

# Set run variables
cell_counts=(4 5 6 7 8)
ogrid_opts=(cub rng sin sqr uni)
tfunc_opts=(exp sig)
rngseed=10

for cells in ${cell_counts[@]}
do
    for ogrid in ${ogrid_opts[@]} 
    do
	for tfunc in ${tfunc_opts[@]}
	do
            for alg in on off new
            do
		echo "-- Run: Grid - ${ogrid}, Test Function - ${tfunc}, Algorithm - ${alg}"
		$build_dir/vert_remap ncell $((2**$cells)) ogrid $ogrid tfunc $tfunc alg $alg
            done
	done
    done
done

echo '-- Runs complete!'
