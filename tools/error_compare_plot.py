###############################################################################
# Import packages

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
from os import mkdir, path

###############################################################################
# Plotting function

def gen_error_comp_plot(outputDir, ncellList, norm, alg):
    
    ###########################################################################
    # Experiment paramters
    nruns = np.size(ncellList)
    ogrids = ["sqr", "cub", "sin", "sig", "sng", "rng"]
    tfuncs = ["sqr", "gau", "exp", "nxp", "asr", "osc"]
    ngrids = len(ogrids)
    nfuncs = len(tfuncs)
    
    ###########################################################################
    # Calculate global error for each run.

    # Normalized error (LeVeque, Appendix A, Eqn. A.17, A.18), but only use
    # the number of cells.
    errorList = np.zeros(shape=nruns)
    accMtx = np.zeros(shape=(nfuncs,ngrids))

    for gIdx in range(len(ogrids)):
        ogrid = ogrids[gIdx]
        for fIdx in range(len(tfuncs)):
            tfunc = tfuncs[fIdx]
            for runIdx in range(nruns):
                ncell = ncellList[runIdx]
                Q2 = np.zeros(shape=(ncell))
                QTrue = np.zeros_like(Q2)
                dp2 = np.zeros_like(Q2)
                runStr = '{0:08d}'.format(ncell)
                fileName = ogrid + '_' + tfunc + '_' + runStr + '_' + alg + '.nc'
                filePath = path.join(outputDir, fileName)
                with xr.open_dataset(filePath) as ds:
                    Q2[:] = np.asarray(ds.Q2)[:]
                    QTrue[:] = np.asarray(ds.QTrue)[:]
                    dp2[:] = np.asarray(ds.dp2)[:]


                # Calculate error in desired norm
                if norm == 1:
                    errorList[runIdx] = np.sum(np.multiply(dp2, np.abs(Q2 - QTrue)))
                elif norm == 2:
                    errorList[runIdx] = np.sqrt(np.sum(np.multiply(dp2, np.abs(Q2 - QTrue)**2)))
                elif norm == -1:
                    errorList[runIdx] = np.max(np.abs(Q2 - QTrue))
                

            ############################################################################
            # Calculate best-fit line for the errors
            
            # Get the best-fit coefficients in log-log space, then convert back to
            # linear space for the plot
            bfLineCoeff = np.polyfit(np.log(ncellList), np.log(errorList), deg = 1)
            accMtx[fIdx,gIdx] = -bfLineCoeff[0]
            
    
    ###########################################################################
    # Plot the errors
    
    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)
        
    # Set up figure.
    fig = plt.figure(figsize=(10,8))
    
    if norm >= 1:
        normstr = str(norm)
    elif norm == -1:
        normstr = "Max"
        
    plt.title("Accuracy Comparison - Algorithm: " + alg + ", Norm - " + normstr)
    cmap = mpl.cm.gray
    bounds = [1, 1.5, 2.0]
    cnorm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')

    fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap),
             orientation='horizontal',
             label="Order of Accuracy")

    img = plt.imshow(accMtx, interpolation='none', cmap=cmap, norm=cnorm)

    plt.xlabel('Test Function')
    plt.xticks(np.arange(len(tfuncs)), tfuncs)
    plt.ylabel('Original Grid')
    plt.yticks(np.arange(len(ogrids)), ogrids)
    
    # Save the plot to file
    fileName = 'err_comp_' + normstr + '_' + alg + '.png'
    plt.savefig(path.join(plotsPath, fileName), dpi = 300)
