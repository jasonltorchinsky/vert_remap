#!/bin/bash
# Generates plots for all of the basic runs

echo '-- Generating plots...'

# Which plots to generate
plterr=yes
pltpwerr=yes
pltapprox=yes
pltmassdiff=yes
pltmasserr=yes

# Set run variables
ogrid_opts=(cub sig sin sqr rng)
tfunc_opts=(exp sig stp wdg)

for ogrid in ${ogrid_opts[@]}
do
    for tfunc in ${tfunc_opts[@]}
    do
	for lim in on off
	do
	    python main.py --ogrid=$ogrid --tfunc=$tfunc \
		   --plterr=$plterr --pltpwerr=$pltpwerr \
		   --pltapprox=$pltapprox --pltmassdiff=$pltmassdiff \
		   --pltmasserr=$pltmasserr --lim=$lim
	done
    done
done

echo '-- Plots complete!'
