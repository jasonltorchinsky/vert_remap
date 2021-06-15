#!/bin/bash
# Runs the vertical remapping code on a local machine.

# Set directory path variables.
echo '-- Re-building executable...'
run_dir=`pwd`
build_dir=../build

# Re-build the application
cd "$build_dir"
rm -rf *
cmake ..
cmake --build .

exec=vert_remap
if test -f "$exec"; then

    # Set run variables
    echo '-- Setting run variables...'
    cd "$run_dir"
    cell_counts=(4 5 6 7 8 9 10)
    ogrid_opts=(cub rng sqr sin)
    tfunc_opts=(exp nxp sig sqr)
    alg_opts=(on off ngh)
    rngseed=42
    
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
    
else

    echo '-- Executable build failed!'

fi
