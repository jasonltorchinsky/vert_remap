###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_pw_error_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc):
    
    ###########################################################################
    # Experiment paramters
    nruns = np.size(ncellList)
    
    ###########################################################################
    # Calculate difference between each approximation and each truth.

    gridDict = dict()
    diffDict = dict()

    print(" ~~ Reading output from each run and calculating pointwise error...")
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        QNew = np.zeros(shape=(ncell))
        QTrue = np.zeros_like(QNew)
        grid2 = np.zeros_like(QNew)
        runStr = '{0:08d}'.format(ncell)
        fileName = ogrid + '_' + tfunc + '_' + runStr + '.nc'
        filePath = path.join(outputDir, fileName)
        with xr.open_dataset(filePath) as ds:
            QNew[:] = np.asarray(ds.QNew)[:]
            QTrue[:] = np.asarray(ds.QTrue)[:]
            grid2[:] = np.asarray(ds.grid2_stg)[:]

        # Calculate pointwise error
        gridDict[ncell] = grid2
        diffDict[ncell] = np.abs(QTrue - QNew)
    
    ###########################################################################
    # Plot the pointwise errors

    print(" ~~ Plotting the pointwise errors...")
    
    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    ## Axes.
    ax = plt.axes()
    
    ## Plot the errors
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        lineLabel = 'Cell Count: ' + '{}'.format(ncell)
        plt.plot(gridDict[ncell], diffDict[ncell], label = lineLabel)
    
    # Add plot titles and labels
    plt.title('Original Grid: ' + ogridFunc + ', Test Function: ' + tfuncFunc)
    plt.xlabel('$x$')
    plt.ylabel('Pointwise Error')

    # Add legend
    plt.legend()
    
    # Save the plot to file
    fileName = 'pw_err_' + ogrid + '_' + tfunc + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
