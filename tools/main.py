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
from pw_diff_plot import gen_pw_diff_plot


################################################################################
# Directory paths

outputDir = path.join("..", "output")

################################################################################
# Read command line arguments

opts, args = getopt.getopt(argv[1:], '', {'ogrid=', 'tfunc=', 'plterr=',
                                          'pltpwerr=', 'pltapprox=',
                                          'pltmassdiff=', 'pltmasserr=',
                                          'pltpwdiff=',
                                          'alg='})

ogrid = ''
tfunc = ''
plterr = ''
pltpwerr = ''
pltapprox = ''
pltmassdiff = ''
pltmasserr = ''
pltpwdiff = ''
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
    elif opt == '--pltpwdiff':
        pltpwdiff = arg
    elif opt == '--alg':
        alg = arg


################################################################################
# Experiment parameters

ncellList = 2**np.arange(4,10,1)


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
elif tfunc == 'nxp':
    tfuncFunc = '$e^{(1-x)}$'
elif tfunc == 'sqr':
    tfuncFunc = '$4\,(x-0.5)^2$'
elif tfunc == 'osc':
    tfuncFunc = '$sin(8\,\pi\,x)$'
elif tfunc == 'gau':
    tfuncFunc = '$e^{-(1/2)\,((x-0.5)/0.05)^2}$'
elif tfunc == 'asr':
    tfuncFunc = '$(16/9)\,(x - 1/4)^2$'
else:
    tfuncFunc  = ''

# Generate plots.
if (plterr == 'yes'):
    print(' ~~ Creating error plot...')
    gen_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Error plot created!')
if (pltpwerr == 'yes'):
    print(' ~~ Creating point-wise error plot...')
    gen_pw_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Point-wise error plot created!')
if (pltapprox == 'yes'):
    print(' ~~ Creating approximation plot...')
    gen_approx_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Approximation plot created!')
if (pltmassdiff == 'yes'):
    print(' ~~ Creating mass-difference plot...')
    gen_massdiff_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Mass-difference plot created!')
if (pltmasserr == 'yes'):
    print(' ~~ Creating mass-error plot...')
    gen_masserror_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Mass-error plot created~')
if (pltpwdiff == 'yes'):
    print(' ~~ Creating point-wise difference plot...')
    gen_pw_diff_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc, alg)
    print(' ~~ Point-wise difference plot created!')
