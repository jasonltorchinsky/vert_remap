#!/bin/bash
# Runs the vertical remapping code on a local machine.

echo '-- Setting run variables...'

# Set directory path variables.
run_dir=`pwd`
build_dir=../build

# Set run variables
ogrid_opts=(cub sig sin sqr rng)
tfunc_opts=(exp sig stp wdg)
rngseed=10

for cells in {8..16..1}
do
    for ogrid in ${ogrid_opts[@]} 
    do
	for tfunc in ${tfunc_opts[@]}
	do
            for lim in on off
            do
		$build_dir/vert_remap ncell $((2**$cells)) ogrid $ogrid tfunc $tfunc lim $lim
            done
	done
    done
done

echo '-- Runs complete!'
