#!/bin/bash
# Generates plots for all of the basic runs

echo '-- Generating plots...'

# Which plots to generate
plterr=yes
pltpwerr=no
pltapprox=no
pltmassdiff=no
pltmasserr=no
pltpwdiff=yes

# Set run variables
ogrid_opts=(cub rng sqr sin)
tfunc_opts=(exp nxp sig sqr)
alg_opts=(on off ngh)

for ogrid in ${ogrid_opts[@]}
do
    for tfunc in ${tfunc_opts[@]}
    do
	for alg in ${alg_opts[@]}
	do
	    echo "-- Run: Grid - ${ogrid}, Test Function - ${tfunc}, Algorithm - ${alg}"
	    python main.py --ogrid=$ogrid --tfunc=$tfunc \
		   --plterr=$plterr --pltpwerr=$pltpwerr \
		   --pltapprox=$pltapprox --pltmassdiff=$pltmassdiff \
		   --pltmasserr=$pltmasserr --pltpwdiff=$pltpwdiff \
		   --alg=$alg
	done
    done
done

echo '-- Plots complete!'
