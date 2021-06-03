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
ogrid_opts=(sqr)
tfunc_opts=(exp)

for ogrid in ${ogrid_opts[@]}
do
    for tfunc in ${tfunc_opts[@]}
    do
	for alg in on off new
	do
	    python main.py --ogrid=$ogrid --tfunc=$tfunc \
		   --plterr=$plterr --pltpwerr=$pltpwerr \
		   --pltapprox=$pltapprox --pltmassdiff=$pltmassdiff \
		   --pltmasserr=$pltmasserr --alg=$alg
	done
    done
done

echo '-- Plots complete!'
