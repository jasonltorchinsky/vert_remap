###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_masserror_plot(outputDir, ncellList, ogrid, ogridFunc, tfunc, tfuncFunc,
                       alg):
    
    ###########################################################################
    # Experiment paramters
    nruns = np.size(ncellList)
    
    ###########################################################################
    # Calculate global error for each run.

    # Normalized error (LeVeque, Appendeix A, Eqn. A.17, A.18), but only use
    # the number of cells.
    error1List = np.zeros(shape=nruns)
    error2List = np.zeros_like(error1List)
    errorMaxList = np.zeros_like(error1List)

    print(" ~~ Reading output from each run and calculating global " +
          "mass error...")
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

        # Calculate error in the three norms
        error1List[runIdx] = np.sum(np.multiply(dp2**2, np.abs(QNew - QTrue)))
    
        error2List[runIdx] = np.sqrt(np.sum(
            np.multiply(dp2**3, np.abs(QNew - QTrue)**2)))

        errorMaxList[runIdx] = np.max(np.multiply(dp2, np.abs(QNew - QTrue)))

    ############################################################################
    # Calculate best-fit line for the errors
    
    # Get the best-fit coefficients in log-log space, then convert back to
    # linear space for the plot
    bfLine1Coeff = np.polyfit(np.log(ncellList), np.log(error1List), deg = 1)
    bfLine1 = np.exp(bfLine1Coeff[1] + bfLine1Coeff[0]*np.log(ncellList))
    
    bfLine2Coeff = np.polyfit(np.log(ncellList), np.log(error2List), deg = 1)
    bfLine2 = np.exp(bfLine2Coeff[1] + bfLine2Coeff[0]*np.log(ncellList))
    
    bfLineMaxCoeff = np.polyfit(np.log(ncellList), np.log(errorMaxList),
                                deg = 1)
    bfLineMax = np.exp(bfLineMaxCoeff[1] + bfLineMaxCoeff[0]*np.log(ncellList))
    
    ###########################################################################
    # Plot the errors

    print(" ~~ Plotting the errors...")
    
    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    ## Axes.
    ax = plt.axes()
    ax.set_xscale("log")
    
    ax.set_yscale("log")
    
    ## Plot the errors
    plt.plot(ncellList, error1List, 'bo', label = 'One-Norm')
    plt.plot(ncellList, bfLine1, 'b--')
    plt.plot(ncellList, error2List, 'go', label = 'Two-Norm')
    plt.plot(ncellList, bfLine2, 'g--')
    plt.plot(ncellList, errorMaxList, 'ko', label = 'Max-Norm')
    plt.plot(ncellList, bfLineMax, 'k--')
    
    # Add plot titles and labels
    plt.title('Original Grid: ' + ogridFunc + '\nTest Function: ' + tfuncFunc + 
              '\nAlgorithm: ' + alg.title())
    plt.xlabel('Number of Cells')
    plt.ylabel('Mass Error')
    
    # Add legend
    plt.legend()
    
    # Add text boxes for equations of best-fit lines
    bfLine1Ord = 'Order: ' + '{:.3f}'.format(-bfLine1Coeff[0])
    ax.text(0.1, 0.05, bfLine1Ord, transform = ax.transAxes, color = 'b')
    
    bfLine2Ord = 'Order: ' + '{:.3f}'.format(-bfLine2Coeff[0])
    ax.text(0.3, 0.05, bfLine2Ord, transform = ax.transAxes, color = 'g')
    
    bfLineMaxOrd = 'Order: ' + '{:.3f}'.format(-bfLineMaxCoeff[0])
    ax.text(0.5, 0.05, bfLineMaxOrd, transform = ax.transAxes, color = 'k')
    
    # Save the plot to file
    fileName = 'masserr_' + ogrid + '_' + tfunc + '_' + alg + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
