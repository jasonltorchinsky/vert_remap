###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_approx_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc,
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

    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        Q2 = np.zeros(shape=(ncell))
        QTrue = np.zeros_like(Q2)
        grid2 = np.zeros_like(Q2)
        runStr = '{0:08d}'.format(ncell)
        fileName = ogrid + '_' + tfunc + '_' + runStr + '_' + alg + '.nc'
        filePath = path.join(outputDir, fileName)
        with xr.open_dataset(filePath) as ds:
            Q2[:] = np.asarray(ds.Q2)[:]
            QTrue[:] = np.asarray(ds.QTrue)[:]
            grid2[:] = np.asarray(ds.grid2_stg)[:]

        # Calculate pointwise error
        gridDict[ncell] = grid2
        trueDict[ncell] = QTrue
        approxDict[ncell] = Q2
    
    ###########################################################################
    # Plot the approximations

    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    ## Axes.
    ax = plt.axes()
    
    ## Plot the truth, approximations
    lineLabel = 'Exact'
    plt.plot(gridDict[ncellList[nruns-1]], trueDict[ncellList[nruns-1]], 'k--', label = lineLabel)
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        lineLabel = 'Cell Count: ' + '{}'.format(ncell)
        plt.plot(gridDict[ncell], approxDict[ncell], label = lineLabel)
    
    # Add plot titles and labels
    plt.title('Original Grid: ' + ogridFunc + '\nTest Function: ' + tfuncFunc + 
              '\nAlgorithm: ' + alg.title())
    plt.xlabel('$x$')
    plt.ylabel('$Q$')

    # Add legend
    plt.legend()
    
    # Save the plot to file
    fileName = 'approx_' + ogrid + '_' + tfunc + '_' + alg + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
