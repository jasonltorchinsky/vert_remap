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
cell_counts=(3 4 5)
ogrid_opts=(sqr)
tfunc_opts=(exp)
alg_opts=(new)
rngseed=10

for cells in ${cell_counts[@]}
do
    for ogrid in ${ogrid_opts[@]} 
    do
	for tfunc in ${tfunc_opts[@]}
	do
            for alg in ${alg_opts[@]}
            do
		echo "-- Run: Grid - ${ogrid}, Test Function - ${tfunc}, Algorithm - ${alg}"
		$build_dir/vert_remap --verbose ncell $((2**$cells)) ogrid $ogrid tfunc $tfunc alg $alg seed $rngseed
            done
	done
    done
done

echo '-- Runs complete!'
