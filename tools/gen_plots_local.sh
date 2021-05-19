#!/bin/bash
# Generates plots for all of the basic runs

echo '-- Generating plots...'

for ogrid in sqr cub sig sin
do
    for tfunc in sig exp stp
    do	
	python main.py --ogrid=$ogrid --tfunc=$tfunc
    done
done

echo '-- Plots complete!'
