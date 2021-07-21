#!/bin/bash
# Generates plots for all of the basic runs

echo '-- Generating plots...'

# Which plots to generate
plterr=no
pltpwerr=no
pltapprox=no
pltmassdiff=no
pltmasserr=no
pltpwdiff=no
plterrcomp=yes

# Set run variables
ogrid_opts=(cub)
tfunc_opts=(asr)
alg_opts=(gmb)
norm_opts=(-1 1 2)

for ogrid in ${ogrid_opts[@]}
do
    for tfunc in ${tfunc_opts[@]}
    do
	for alg in ${alg_opts[@]}
	do
            for norm in ${norm_opts[@]}
            do
	        echo "-- Run: Grid - ${ogrid}, Test Function - ${tfunc}, Algorithm - ${alg}"
        	    python main.py --ogrid=$ogrid --tfunc=$tfunc \
	               --plterr=$plterr --pltpwerr=$pltpwerr \
		       --pltapprox=$pltapprox --pltmassdiff=$pltmassdiff \
		       --pltmasserr=$pltmasserr --pltpwdiff=$pltpwdiff \
		       --plterrcomp=$plterrcomp --norm=$norm  --alg=$alg
            done
	done
    done
done

echo '-- Plots complete!'
