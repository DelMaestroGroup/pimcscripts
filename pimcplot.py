#!/usr/bin/python
"""pimcplot

Description:
  Performs a cumulative average plot of raw Monte Carlo data

Usage:
    pimcplot.py [options] [--legend=<label>...] --estimator=<name> <file>...

    pimcplot.py -h | --help 

Options:
  -h, --help                    Show this screen.
  --estimator=<name>, -e <name> The estimator to be plotted.
  --skip=<n>, -s <n>            Number of measurements to be skipped [default: 0].
  --period=<m>, -p <m>          The period of the average window [default: 50].
  --legend=<label>, -l <label>  A legend label
  --error=<units>, -d           Size of the error bars
  --bin                         Use the binned errorbars
"""

# pimcplot.py
# Adrian Del Maestro
# 07.20.2009
# 
# Plot rough estimators vs. MC Bins for files supplied as input

import matplotlib
#matplotlib.use('TKAgg')

import os,sys
import pyutils
import loadgmt,kevent
from pylab import *
import pimchelp
from docopt import docopt
from scipy import stats
import MCstat

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if ndim(data) > dim:
        numBins  = size(data,dim) 
        dataAve  = average(data,dim) 
        dataAve2 = average(data*data,dim) 
        bins = MCstat.bin(data) 
        dataErr = amax(bins,axis=0)
        dataErr2 = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 

#        try:
#            bins = MCstat.bin(data) 
#            dataErr = amax(bins,axis=0)
#        except:
#            dataErr   = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
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
    args = docopt(__doc__)

    fileNames = args['<file>']
    skip = int(args['--skip'])
    period = int(args['--period'])
    estimator = args['--estimator']
    leglabel = args['--legend'] and args['--legend']
    error = args['--error'] and float(args['--error'])

    # if labels are not assigned, we default to the PIMCID
    if not leglabel:
        leglabel = []
        for n,fileName in enumerate(fileNames):
            leglabel.append(fileName[-13:-4])

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

    # ============================================================================
    # Figure 1 : column vs. MC Steps
    # ============================================================================
    figure(1)
    connect('key_press_event',kevent.press)

    colors  = loadgmt.getColorList('cw/1','cw1-029',max(numFiles,2))

    n = 0
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();
        dataFile.close()

        if len(dataLines) > 2:
            data = loadtxt(fileName,usecols=col)
            
            if not pyutils.isList(data):
               data = list([data])

            plot(data[skip:],marker='s',color=colors[n],markeredgecolor=colors[n],\
                        markersize=4,linestyle='-',linewidth=1.0)
    
            n += 1

    ylabel(yLong)
    xlabel("MC Bin Number")

    # ============================================================================
    # Figure 2 : running average of column vs. MC Bins
    # ============================================================================
    figure(2)
    connect('key_press_event',kevent.press)

    n = 0
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();

        if len(dataLines) > 2:

            data = loadtxt(fileName,usecols=col)
            if size(data) > 1:
                
                # Get the cumulative moving average
                if args['--error']:
                    cma = cumulativeMovingAverage(data[skip:])
                    sem = error*ones_like(cma)
                elif args['--bin']:
                    cma = cumulativeMovingAverage(data[skip:])
                    ave,err = getStats(data[skip:])
                    sem = err*ones_like(cma)
                    print '%s:  %s = %8.4E +- %8.4E' % (leglabel[n],yShort, ave,err) 
                else:
                    cma,sem = cumulativeMovingAverageWithError(data[skip:])

                sma = simpleMovingAverage(50,data[skip:])
                x = range(int(0.04*len(cma)),len(cma))
                plot(x,cma[x],color=colors[n],linewidth=1.0,marker='None',linestyle='-',
                    label=leglabel[n])
                fill_between(x, cma[x]-sem[x], cma[x]+sem[x],color=colors[n], alpha=0.1)
                n += 1

    ylabel(yLong)
    xlabel("MC Bin Number")
    tight_layout()
    legend(loc='best', frameon=False, prop={'size':14},ncol=2)

    show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
