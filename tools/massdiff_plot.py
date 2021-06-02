###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_massdiff_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc,
                      alg):
    
    ###########################################################################
    # Experiment paramters
    nruns = np.size(ncellList)
    
    ###########################################################################
    # Get the difference in mass between the exact and approximate solutions.

    massDiffList = np.zeros(shape=nruns)

    print(" ~~ Reading output from each run and calculating the " +
          "total mass difference...")
    for runIdx in range(nruns):
        ncell = ncellList[runIdx]
        QNew = np.zeros(shape=(ncell))
        QTrue = np.zeros_like(QNew)
        dp2 = np.zeros_like(QNew)
        runStr = '{0:08d}'.format(ncell)
        fileName = ogrid + '_' + tfunc + '_' + runStr + '_' + alg + '.nc'
        filePath = path.join(outputDir, fileName)
        with xr.open_dataset(filePath) as ds:
            QNew[:] = np.asarray(ds.QNew)[:]
            QTrue[:] = np.asarray(ds.QTrue)[:]
            dp2[:] = np.asarray(ds.dp2)[:]

        # Calculate total mass error
        trueMass = np.sum(np.multiply(QTrue, dp2))
        approxMass = np.sum(np.multiply(QNew, dp2))
        massDiffList[runIdx] = approxMass - trueMass
        
    ###########################################################################
    # Plot the total mass differences

    print(" ~~ Plotting the total mass difference...")
    
    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    ## Axes.
    ax = plt.axes()
    ax.set_xscale("log")
    
    ## Plot the mass differences
    plt.plot(ncellList, massDiffList, 'k-')
    
    # Add plot titles and labels
    plt.title('Original Grid: ' + ogridFunc + '\nTest Function: ' + tfuncFunc + 
              '\nAlgorithm: ' + alg.title())
    plt.xlabel('Cell Count')
    plt.ylabel('Mass Difference (True - Approximation)')

    # Save the plot to file
    fileName = 'massdiff_' + ogrid + '_' + tfunc + '_' + alg + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
