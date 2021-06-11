#!/bin/bash
# Runs the vertical remapping code on a local machine.

# Set directory path variables.
run_dir=`pwd`
build_dir=../build

# Re-build the application
echo '-- Re-building executable...'
cd "$build_dir"
rm -rf *
cmake ..
cmake --build .


exec=vert_remap
if test -f "$exec"; then
    # Set run variables
    echo '-- Setting run variables...'
    cd "$run_dir"
    cell_counts=(3)
    ogrid_opts=(cub)
    tfunc_opts=(exp)
    rngseed=10

    for cells in ${cell_counts[@]}
    do
	for ogrid in ${ogrid_opts[@]} 
	do
	    for tfunc in ${tfunc_opts[@]}
	    do
		for alg in new
		do
		    gdb -tui --args $build_dir/vert_remap --verbose ncell $((2**$cells)) ogrid $ogrid tfunc $tfunc alg $alg
		done
	    done
	done
    done
    
    echo '-- Runs complete!'

else

    echo '-- Executable build failed!'

fi
