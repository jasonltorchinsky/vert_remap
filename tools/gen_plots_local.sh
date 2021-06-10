#!/bin/bash
# Generates plots for all of the basic runs

echo '-- Generating plots...'

# Which plots to generate
plterr=yes
pltpwerr=yes
pltapprox=no
pltmassdiff=no
pltmasserr=no

# Set run variables
ogrid_opts=(cub rng sig sin sqr uni)
tfunc_opts=(exp sig stp wdg)
alg_opts=(on off new)

for ogrid in ${ogrid_opts[@]}
do
    for tfunc in ${tfunc_opts[@]}
    do
	for alg in ${alg_opts[@]}
	do
	    python main.py --ogrid=$ogrid --tfunc=$tfunc \
		   --plterr=$plterr --pltpwerr=$pltpwerr \
		   --pltapprox=$pltapprox --pltmassdiff=$pltmassdiff \
		   --pltmasserr=$pltmasserr --alg=$alg
	done
    done
done

echo '-- Plots complete!'
