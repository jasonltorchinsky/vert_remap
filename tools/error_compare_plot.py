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
    if alg != '':
        algs = [alg]
    else:
        algs = ["on", "off", "gmb", "lco", "fff"]
    ngrids = len(ogrids)
    nfuncs = len(tfuncs)
    nalgs = len(algs)

    ###########################################################################
    # Set up the subplots

    # Make sure the output directory for plots has been created.
    plotsPath = path.join('..', 'plots')
    if not path.exists(plotsPath):
        mkdir(plotsPath)

    fig, axs = plt.subplots(1, nalgs, figsize=(10,4))

    if norm >= 1:
        normstr = str(norm) + '-Norm'
    elif norm == -1:
        normstr = "Max-Norm"

    if alg != '':
        outfileName = 'err_comp_' + normstr + '_' + alg + '.png'
    else:
        outfileName = 'err_comp_' + normstr + '.png'
        
    fig.suptitle(normstr + ' Accuracy')

    cmap = mpl.cm.gray
    bounds = [1, 1.5, 2.0]
    cnorm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='both')
    
    
    ###########################################################################
    # Create each sub-plot

    for algIdx in range(nalgs):
        alg = algs[algIdx]
        # Normalized error (LeVeque, Appendix A, Eqn. A.17, A.18), but only use
        # the number of cells.
        errorList = np.zeros(shape=nruns)
        accMtx = np.zeros(shape=(nfuncs,ngrids))
        
        for gIdx in range(ngrids):
            ogrid = ogrids[gIdx]
            for fIdx in range(nfuncs):
                tfunc = tfuncs[fIdx]
                for runIdx in range(nruns):
                    ncell = ncellList[runIdx]
                    Q2 = np.zeros(shape=(ncell))
                    QTrue = np.zeros_like(Q2)
                    dp2 = np.zeros_like(Q2)
                    runStr = '{0:08d}'.format(ncell)
                    fileName = ogrid + '_' + tfunc + '_' + runStr + '_' \
                        + alg + '.nc'
                    filePath = path.join(outputDir, fileName)
                    with xr.open_dataset(filePath) as ds:
                        Q2[:] = np.asarray(ds.Q2)[:]
                        QTrue[:] = np.asarray(ds.QTrue)[:]
                        dp2[:] = np.asarray(ds.dp2)[:]
                        
                        
                    # Calculate error in desired norm
                    if norm == 1:
                        errorList[runIdx] = np.sum(
                            np.multiply(dp2, np.abs(Q2 - QTrue)) )
                    elif norm == 2:
                        errorList[runIdx] = np.sqrt(
                            np.sum(np.multiply(dp2, np.abs(Q2 - QTrue)**2)) )
                    elif norm == -1:
                        errorList[runIdx] = np.max(np.abs(Q2 - QTrue))
                

                ###################################################################
                # Calculate the order of accuracy.
            
                bfLineCoeff = np.polyfit( np.log(ncellList), np.log(errorList),
                                          deg = 1 )
                accMtx[fIdx,gIdx] = -bfLineCoeff[0]
            
    
        ###################################################################
        # Create the subplot
        if (alg == 'on'):
            algstr = 'Strategy (0)'
        elif (alg == 'off'):
            algstr = 'Strategy (1)'
        elif (alg == 'gmb'):
            algstr = 'Strategy (2)'
        elif (alg == 'lco'):
            algstr = 'Strategy (3)'
        elif (alg == 'fff'):
            algstr = 'Strategy (4)'
        else:
            algstr = 'Strategy (' + alg + ')'
            
            
        axs[algIdx].set_title(algstr)
        axs[algIdx].imshow(accMtx, interpolation='none', cmap=cmap, norm=cnorm)

    for ax in axs.flat:
        ax.set( xlabel='Test Function',
                ylabel='Original Grid',
                xticks=np.arange(nfuncs),
                xticklabels=['$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$'],
                yticks=np.arange(ngrids),
                yticklabels=['$g_1$','$g_2$','$g_3$','$g_4$','$g_5$','$g_6$'] )

    for ax in axs.flat:
        ax.label_outer()

        
    fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap),
                 ax=axs,
                 orientation='horizontal',
                 label="Order of Accuracy")
    
    # Save the plot to file
    fig.savefig(path.join(plotsPath, outfileName), dpi = 300)
