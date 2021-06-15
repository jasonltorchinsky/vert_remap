###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_pw_diff_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc,
                    alg):
    
    ###########################################################################
    # Experiment paramters
    nruns = np.size(ncellList)
    
    ###########################################################################
    # Get each approximation and each truth, although we will only use the
    # highest-resolution truth.

    gridDict = dict()
    trueDict = dict()
    approxDict = dict()

    print(" ~~ Reading output from each run and getting approximation...")
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        QNew = np.zeros(shape=(ncell))
        QTrue = np.zeros_like(QNew)
        grid2 = np.zeros_like(QNew)
        runStr = '{0:08d}'.format(ncell)
        fileName = ogrid + '_' + tfunc + '_' + runStr + '_' + alg + '.nc'
        filePath = path.join(outputDir, fileName)
        with xr.open_dataset(filePath) as ds:
            QNew[:] = np.asarray(ds.QNew)[:]
            QTrue[:] = np.asarray(ds.QTrue)[:]
            grid2[:] = np.asarray(ds.grid2_stg)[:]

        # Calculate pointwise error
        gridDict[ncell] = grid2
        trueDict[ncell] = QTrue
        approxDict[ncell] = QNew
    
    ###########################################################################
    # Plot the approximations

    print(" ~~ Plotting the pointwise differences...")
    
    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    ## Axes.
    ax = plt.axes()
    
    ## Plot the truth, approximations
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        lineLabel = 'Cell Count: ' + '{}'.format(ncell)
        plt.plot(gridDict[ncell], trueDict[ncell] - approxDict[ncell], label = lineLabel)
    
    # Add plot titles and labels
    plt.title('Original Grid: ' + ogridFunc + '\nTest Function: ' + tfuncFunc + 
              '\nAlgorithm: ' + alg.title())
    plt.xlabel('$x$')
    plt.ylabel('Difference (True - Approximation)')

    # Add legend
    plt.legend()
    
    # Save the plot to file
    fileName = 'pw_diff_' + ogrid + '_' + tfunc + '_' + alg + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
