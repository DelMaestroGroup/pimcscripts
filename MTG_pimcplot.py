#!/usr/bin/python

# pimcplot.py
# Adrian Del Maestro
# 07.20.2009
# 
# Plot rough estimators vs. MC Bins for files supplied as input
# This does the same thing as other pimcplot, but uses argparse.

import matplotlib
import os,sys
import pyutils
import loadgmt,kevent
from pylab import *
import pimchelp
import argparse
from scipy import stats
import MCstat

# -----------------------------------------------------------------------------
def parseCMD():
    ''' Setup the command line parser options. '''
    
    desc = "Performs a cumulative average plot of raw Monte Carlo data"
    parser = argparse.ArgumentParser(description=desc) 
    parser.add_argument("file", help='data file', nargs='+')
    parser.add_argument('-e', '--estimator', type=str,
            help="The estimator to be plotted.")
    parser.add_argument('-s', '--skip', type=int, default=0,
            help="The number of measurements to skip.")
    parser.add_argument('-p', '--period', type=int, default=50,
            help="The period of the average window.")
    parser.add_argument('-l', '--legend', type=str,
            help="The legend label.")
    parser.add_argument('-d', '--error', type=float,
            help="The size of the error bars.")
    parser.add_argument('--bin', action='store_true',dest='bin',
            default=True,
            help="Use the binned errorbars.")
    parser.add_argument('--ttest', action='store_true',dest='ttest',
            default=True,
            help="Perform a ttest.")

    return parser.parse_args() 

# -----------------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if ndim(data) > dim:
        numBins  = size(data,dim) 
        dataAve  = average(data,dim) 
        dataAve2 = average(data*data,dim) 
 
        try:
            bins = MCstat.bin(data) 
            dataErr = amax(bins,axis=0)
        except:
            dataErr   = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
#        bins = MCstat.bin(data) 
#        dataErr = amax(bins,axis=0)
        dataErr2 = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr

# -----------------------------------------------------------------------------
def cumulativeMovingAverage(data):
    '''Compute the cumulative mean as a function of bin index.'''
    CMA = zeros_like(data)
    CMA[0] = data[0]
    for n in range(len(data)-1):
        CMA[n+1] = (data[n+1] + n*CMA[n])/(1.0*(n+1))
    return CMA

# -----------------------------------------------------------------------------
def cumulativeMovingAverageWithError(data):
    '''Compute the cumulative mean as a function of bin index.'''
    CMA = zeros_like(data)
    CME = zeros_like(data)
    CMA[0] = data[0]
    CME[0] = 0.0;
    for n in range(len(data)-1):
        CMA[n+1] = average(data[:n+2])
        CME[n+1] = stats.sem(data[:n+2])

    return CMA,CME

# -----------------------------------------------------------------------------
def simpleMovingAverage(period,data):
    assert period == int(period) and period > 0, "Period must be an integer >0"

    import numpy as np
    weightings = np.repeat(1.0, period) / period 
    return np.convolve(data, weightings, 'valid')

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # parse the command line options
    args = parseCMD()

    fileNames = args.file
    skip = args.skip
    period = args.period
    estimator = args.estimator
    leglabel = [args.legend]
    error = args.error

    # if labels are not assigned, we default to the PIMCID
    if not leglabel:
        leglabel = []
        for n,fileName in enumerate(fileNames):
            leglabel.append(fileName[-13:-4])
            print fileName

    # We count the number of lines in the estimator file to make sure we have
    # some data and grab the headers
    headers = pimchelp.getHeadersDict(fileNames[0])

    # If we don't choose an estimator, provide a list of possible ones
    if estimator not in headers:
        errorString = "Need to specify one of:\n"
        for head,index in headers.iteritems():
            errorString += "\"%s\"" % head + "   "
        parser.error(errorString)

    numFiles = len(fileNames)
    col = list([headers[estimator]])

    # Attempt to find a 'pretty name' for the label, otherwise just default to
    # the column heading
    label = pimchelp.Description()
    try:
        yLong = label.estimatorLongName[estimator]
    except:
        yLong = estimator
    try:
        yShort = label.estimatorShortName[estimator]
    except:
        yShort = estimator

    # First we load and store all the data 
    data = []
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();
        dataFile.close()

        if len(dataLines) > 2:
            data.append(loadtxt(fileName,usecols=col,unpack=True))

    # ============================================================================
    # Figure 1 : column vs. MC Steps
    # ============================================================================
    figure(1)
    connect('key_press_event',kevent.press)

    #colors  = loadgmt.getColorList('cw/1','cw1-029',max(numFiles,2))
    #colors  = loadgmt.getColorList('cw/1','cw1-013',max(numFiles,2))
    colors  = loadgmt.getColorList('oc', 'rainbow', max(numFiles,2))
#    oc/rainbow

    for n,cdata in enumerate(data):
        plot(cdata[skip:],marker='s',color=colors[n],markeredgecolor=colors[n],\
             markersize=4,linestyle='-',linewidth=1.0)

    ylabel(yLong)
    xlabel("MC Bin Number")

    # ============================================================================
    # Figure 2 : running average of column vs. MC Bins
    # ============================================================================
    figure(2)
    connect('key_press_event',kevent.press)

    n = 0
    for n,cdata in enumerate(data):
        if size(cdata) > 1:
            
            # Get the cumulative moving average
            #if args['--error']:
            if args.error:
                cma = cumulativeMovingAverage(cdata[skip:])
                sem = error*ones_like(cma)
                #elif args['--bin']:
            elif args.bin:
                cma = cumulativeMovingAverage(cdata[skip:])
                ave,err = getStats(cdata[skip:])
                sem = err*ones_like(cma)
                print '%s:  %s = %8.4E +- %8.4E' % (leglabel[n],yShort, ave,err) 
            else:
                cma,sem = cumulativeMovingAverageWithError(cdata[skip:])

            sma = simpleMovingAverage(50,cdata[skip:])
            x = range(int(0.10*len(cma)),len(cma))
            plot(x,cma[x],color=colors[n],linewidth=1.0,marker='None',linestyle='-',
                label=leglabel[n])
            fill_between(x, cma[x]-sem[x], cma[x]+sem[x],color=colors[n], alpha=0.1)
            n += 1

    ylabel(yLong)
    xlabel("MC Bin Number")
    tight_layout()
    try:
        leg = legend(loc='best', frameon=False, prop={'size':16}, markerscale=2, ncol=2)
        for l in leg.get_lines():
            l.set_linewidth(4.0)
    except:
        pass

    # Perform a Welch's t-test
    if args.ttest:
        # We only perform the Welch's t test if we have multiple samples we are
        # comparing
        N = len(data)
        if N > 1:
            tval = zeros([N,N])
            p = zeros([N,N])

            for i in range(N):
                for j in range(i+1,N):
                    tval[i,j],p[i,j] = stats.ttest_ind(data[i][skip:], data[j][skip:], 
                                               equal_var=False)

        # ============================================================================
        # Figure 3 : plot the estimator histogram along with t-test values
        # ============================================================================
        fig = figure(3)
        connect('key_press_event',kevent.press)
        for i in range(N):
            n, bins, patches = hist(data[i], 100, normed=True, facecolor=colors[i], 
                                    alpha=0.75, label=leglabel[i],
                                    edgecolor='w')
        # Add the p-values from the t-test
        y = 0.92
        if N > 1:
            figtext(0.78, y, 't-test p values', horizontalalignment='center', 
                 verticalalignment='top', fontsize=15, backgroundcolor='white')

            for i in range(N):
                for j in range(i+1,N):
                    y -= 0.03
                    lab = 'p(' + leglabel[i] + ' - ' + leglabel[j] + ') = ' + '%4.2f'%p[i,j]
                    figtext(0.78, y, lab, horizontalalignment='center', 
                         verticalalignment='top', fontsize=12,
                       backgroundcolor='white')

        legend(loc='upper left', frameon=False)
        xlabel(yLong)
        ylabel(r'$P($' + estimator + r'$)$')
 

            
    show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
