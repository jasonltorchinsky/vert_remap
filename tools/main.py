"""
This script is intended to generate all desired visualizations.
Author: Jason Torchinsky
Date: 2021-05
"""

################################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path
from sys import argv
import getopt

from error_plot import gen_error_plot
from pw_error_plot import gen_pw_error_plot


################################################################################
# Directory paths

outputDir = path.join("..", "output")

################################################################################
# Read command line arguments

opts, args = getopt.getopt(argv[1:], '', {'ogrid=', 'tfunc='})

ogrid = ''
tfunc = ''
for opt, arg in opts:
    if opt == '--ogrid':
        ogrid = arg
    elif opt == '--tfunc':
        tfunc = arg


################################################################################
# Experiment parameters

ncellList = 2**np.arange(8,16,1)


if ogrid == 'sqr':
    ogridFunc = '$x^2$'
elif ogrid == 'cub':
    ogridFunc = '$x^3$'
elif ogrid == 'sig':
    ogridFunc = '$1$/$(1 + e^{-30\,(x - 0.5)})$'
elif ogrid == 'sin':
    ogridFunc = '$0.5\,(1 - cos(\pi x))$'
else:
    ogridFunc = ''

if tfunc == 'exp':
    tfuncFunc = '$e^x$'
elif tfunc == 'stp':
    tfuncFunc = '$0$ for $x \leq 0.5$, $1$ for $x > 0.5$'
elif tfunc == 'sig':
    tfuncFunc = '$1$/$(1 + e^{-30\,(x - 0.5)})$'
else:
    tfuncFunc  = ''

# Generate plots.
gen_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc)
gen_pw_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc)
