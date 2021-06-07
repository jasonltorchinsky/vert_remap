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
from approx_plot import gen_approx_plot
from massdiff_plot import gen_massdiff_plot
from masserror_plot import gen_masserror_plot


################################################################################
# Directory paths

outputDir = path.join("..", "output")

################################################################################
# Read command line arguments

opts, args = getopt.getopt(argv[1:], '', {'ogrid=', 'tfunc=', 'plterr=',
                                          'pltpwerr=', 'pltapprox=',
                                          'pltmassdiff=', 'pltmasserr=',
                                          'alg='})

ogrid = ''
tfunc = ''
plterr = ''
pltpwerr = ''
pltapprox = ''
pltmassdiff = ''
pltmasserr = ''
alg = ''
for opt, arg in opts:
    if opt == '--ogrid':
        ogrid = arg
    elif opt == '--tfunc':
        tfunc = arg
    elif opt == '--plterr':
        plterr = arg
    elif opt == '--pltpwerr':
        pltpwerr = arg
    elif opt == '--pltapprox':
        pltapprox = arg
    elif opt == '--pltmassdiff':
        pltmassdiff = arg
    elif opt == '--pltmasserr':
        pltmasserr = arg
    elif opt == '--alg':
        alg = arg


################################################################################
# Experiment parameters

ncellList = 2**np.arange(4,8,1)


if ogrid == 'sqr':
    ogridFunc = '$x^2$'
elif ogrid == 'cub':
    ogridFunc = '$x^3$'
elif ogrid == 'sig':
    ogridFunc = '$1$/$(1 + e^{-30\,(x - 0.5)})$'
elif ogrid == 'sin':
    ogridFunc = '$0.5\,(1 - cos(\pi x))$'
elif ogrid == 'rng':
    ogridFunc = 'Random, $[-h/4, h/4]$'
elif ogrid == 'uni':
    ogridFunc = 'Uniform'
else:
    ogridFunc = ''

if tfunc == 'exp':
    tfuncFunc = '$e^x$'
elif tfunc == 'stp':
    tfuncFunc = '$0$ for $x \leq 0.5$, $1$ for $x > 0.5$'
elif tfunc == 'sig':
    tfuncFunc = '$1$/$(1 + e^{-30\,(x - 0.5)})$'
elif tfunc == 'wdg':
    tfuncFunc = '$0.5\,x$ for $x \leq 0.5$, $1.5\,x - 0.5$ for $x > 0.5$'
else:
    tfuncFunc  = ''

# Generate plots.
if (plterr == 'yes'):
    gen_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
if (pltpwerr == 'yes'):
    gen_pw_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
if (pltapprox == 'yes'):
    gen_approx_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
if (pltmassdiff == 'yes'):
    gen_massdiff_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
if (pltmasserr == 'yes'):
    gen_masserror_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
