#!/usr/bin/python
# plotbin.py
# Chris Herdman
# 12.03.2012
# 
# Plot binning analysis for raw MC data

import pyutils

import matplotlib.pyplot as plt
import numpy as np

import argparse
import pimchelp
import MCstat

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 
    
    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Plot binning analysis for \
                                            MC Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    parser.add_argument('--estimator','-e', help='A list of estimator names \
                                            that are to be plotted.', type=str)
    parser.add_argument('--skip','-s', help='Number of measurements to be \
                        skipped in the binning analysis.', type=int, default=0)
    parser.add_argument('--scale', help='Option to compare binning results \
                        for different parameters', action='store_true')
    args = parser.parse_args()
    
    fileNames = args.fileNames
    scale = args.scale
    
    if len(fileNames) < 1:
        parser.error("Need to specify at least one scalar estimator file")
    
    # We count the number of lines in the estimator file to make sure we have
    # some data and grab the headers
    headers = pimchelp.getHeadersDict(fileNames[0])
    
    # If we don't choose an estimator, provide a list of possible ones
    if not args.estimator or args.estimator not in headers:
        errorString = "Need to specify one of:\n"
        for head,index in headers.iteritems():
            errorString += "\"%s\"" % head + "   "
        parser.error(errorString)
    
    numFiles = len(fileNames)
    col = list([headers[args.estimator]])
    
    # Attempt to find a 'pretty name' for the label, otherwise just default to
    # the column heading
    label = pimchelp.Description()
    try:
        yLong = label.estimatorLongName[args.estimator]
    except:
        yLong = args.estimator
    try:
        yShort = label.estimatorShortName[args.estimator]
    except:
        yShort = args.estimator
    
    # ============================================================================
    # Figure 1 : Error vs. bin level
    # ============================================================================
    plt.figure(1)
    
    colors = ["#688EAF", "#FC991D", "#7DEB74", "#FA6781", "#8B981D", "#BB7548", 
            "#AD8FE4", "#96E4AA", "#D669B0", "#E1C947", "#A78200", "#7C9FE4", 
            "#957DA6", "#75BF38", "#C3B059", "#51C17A", "#79AEBB", "#2790AC", 
            "#688ECE", "#749DB7"]
    for i in range(3):
        colors += colors
    
    n = 0
    for fileName in fileNames:
        
        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();
        dataFile.close()
        
        if len(dataLines) > 2:
            data = np.loadtxt(fileName,usecols=col)
            if not pyutils.isList(data):
               data = list([data])
            
            delta = MCstat.bin(data[args.skip:])
            if n == 0:
                delta_ar = np.zeros((numFiles,delta.shape[0]))
            delta_ar[n,:] = delta.T
            #delta_ar[n,:len(delta)] = delta.T
            n += 1
    
    if n > 1:
        if scale:
            for m in range(n):
                plt.plot(np.arange(len(delta_ar)),delta_ar[m],marker='s',markersize=4,\
                         linestyle='-',linewidth=1.0,color=colors[m],\
                         markeredgecolor=colors[m])
        else:
            Delta = np.average(delta_ar,0)
            dDelta = np.std(delta_ar,0)/np.sqrt(n)             
            plt.errorbar(np.arange(len(Delta)),Delta,dDelta,marker='s',markersize=4,\
                     linestyle='-',linewidth=1.0,color=colors[0],\
                     markeredgecolor=colors[0])
            bin_ac = MCstat.bin_ac(Delta,dDelta)
            bin_conv = MCstat.bin_conv(Delta,dDelta)
            print 'Convergence Ratio: %1.2f+/-%1.2f'%(bin_conv['CF'],bin_conv['dCF'])
            print 'autocorrlelation time: %2.1f+/-%2.1f' % \
                                                (bin_ac['tau'],bin_ac['dtau'])
    else:
        plt.plot(delta,marker='s',markersize=4,linestyle='-',linewidth=1.0,\
             color=colors[0],markeredgecolor=colors[0])
        print 'Convergence Ratio: %1.3f' % MCstat.bin_conv(delta)['CF']
        print 'autocorrlelation time: %3.3f' % MCstat.bin_ac(delta)['tau']
    
    plt.ylabel(r"$\Delta_l$")
    plt.xlabel("$l$")
    plt.title("Bin scaling: "+yLong)
    plt.tight_layout()
    
    plt.show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
