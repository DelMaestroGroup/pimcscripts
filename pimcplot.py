#! /usr/bin/env python3
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
  --truncateid=<t>, -t <t>      Truncate PIMCID to last <t> characters [default: 0].
  --legend=<label>, -l <label>  A legend label
  --period=<m>, -p <m>          The period of the average window [default: 50].
  --error=<units>, -d           Size of the error bars
  --nobin                       Don't use the binned errorbars
  --nolegend                    Turn off the legend
  --ttest                       Perform a ttest
  --hline=<val>                 Include a horizontal line at <val> in the averaged plot
  --hlabel=<hl>                 A legend label for the horizontal line.
  --title=<title>               A title for the plots.
  --savefig=<figure>            A filename for saved plots (extensions supported by active matplotlib backend).
  --quiet                       Suppress output.
"""

# pimcplot.py
# Adrian Del Maestro
# 07.20.2009
# 
# Plot rough estimators vs. MC Bins for files supplied as input

from __future__ import print_function 
import numpy as np
import pimchelp
from docopt import docopt
from scipy import stats
import MCstat
import re
from os import path

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if np.ndim(data) > dim:
        numBins  = np.size(data,dim) 
        dataAve  = np.average(data,dim) 
        dataAve2 = np.average(data**2,dim) 
        bins = MCstat.bin(data) 
        dataErr = np.amax(bins,axis=0)
        dataErr2 = np.sqrt( np.abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr

# -----------------------------------------------------------------------------
def cumulativeMovingAverage(data):
    '''Compute the cumulative mean as a function of bin index.'''
    CMA = np.zeros_like(data)
    CMA[0] = data[0]
    for n in range(len(data)-1):
        CMA[n+1] = (data[n+1] + n*CMA[n])/(1.0*(n+1))
    return CMA

# -----------------------------------------------------------------------------
def cumulativeMovingAverageWithError(data):
    '''Compute the cumulative mean as a function of bin index.'''
    CMA = np.zeros_like(data)
    CME = np.zeros_like(data)
    CMA[0] = data[0]
    CME[0] = 0.0;
    for n in range(len(data)-1):
        CMA[n+1] = np.average(data[:n+2])
        CME[n+1] = stats.sem(data[:n+2])

    return CMA,CME

# -----------------------------------------------------------------------------
def simpleMovingAverage(period,data):
    assert period == int(period) and period > 0, "Period must be an integer >0"
    weightings = np.repeat(1.0, period) / period 
    return np.convolve(data, weightings, 'valid')

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # parse the command line options
    args = docopt(__doc__)

    fileNames = args['<file>']
    if args['--skip']:
        if '.' in args['--skip']:
            raw_skip = float(args['--skip'])
        else:
           raw_skip = int(args['--skip'])
           skip=raw_skip

    period = int(args['--period'])
    estimator = args['--estimator']
    leglabel = args['--legend'] and args['--legend']
    error = args['--error'] and float(args['--error'])
    val = args['--hline'] and float(args['--hline'])
    plttitle = args['--title']
    trunc = int(args['--truncateid'])
    if args['--savefig']:
        #FIXME need headless usage
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # if labels are not assigned, we default to the PIMCID
    if not leglabel:
        leglabel = []
        for n,fileName in enumerate(fileNames):
            fn1,fn2 = path.split(fileName)
            pimcid = "-".join(re.split(r'(?<!-)-',fn2)[6:]).split('.')[0] # Look back to avoid splitting negatives and take last portion to get id and strip extension (ask NSN)
            leglabel.append(pimcid[-trunc:])

    # We count the number of lines in the estimator file to make sure we have
    # some data and grab the headers
    headers = pimchelp.getHeadersDict(fileNames[0])

    # If we don't choose an estimator, provide a list of possible ones
    # if estimator not in headers:
    #     errorString = "Need to specify one of:\n"
    #     for head,index in headers.items():
    #         errorString += "\"%s\"" % head + "   "
    #     parser.error(errorString)

    numFiles = len(fileNames)
    col = list([headers[estimator]])

    # Attempt to find a 'pretty name' for the label, otherwise just default to
    # the column heading
    label = pimchelp.Description()
    try:
        yLong = label.estimatorLongName[estimator]
    except:
        yLong = estimator.replace('_', '')
    try:
        yShort = label.estimatorShortName[estimator]
    except:
        yShort = estimator.replace('_', '')

    # get a label for a possible horizontal line
    if args['--hline']:
        if args['--hlabel']:
            hlabel = args['--hlabel']
        else:
            hlabel = yShort.split()[0] + ' = ' + args['--hline']

    # First we load and store all the data 
    data = []
    for fileName in fileNames:

        dataFile = open(fileName,'r');
        dataLines = dataFile.readlines();
        dataFile.close()

        if len(dataLines) > 2:
            data.append(np.loadtxt(fileName,usecols=col,unpack=True))

    # ============================================================================
    # Figure 1 : column vs. MC Steps
    # ============================================================================
    plt.figure(1)

    # colors = ["#70D44A", "#BE5AD4", "#D04537", "#81D0D5", "#393A2E", "#C49ECA", 
    #           "#C5CB7A", "#523767", "#D39139", "#C8488C", "#CB817A", "#73D999", 
    #           "#6F2836", "#6978CF", "#588569", "#CDC3AD", "#5C7890", "#7F5327", 
    #           "#D0D33D", "#5B7D2E"]
    # colors = ["#53B0AD", "#D74C20", "#CD53DA", "#58C038", "#5F4B7A", "#49622A", 
    #           "#CE4379", "#D2912E", "#7970D2", "#749AC9", "#7D5121", "#5EAC72", 
    #           "#CB85AA", "#853B46", "#396465", "#A5A33E", "#D47F5B", "#BA4EA5", 
    #           "#C93F44", "#5D9937"]
    # colors =["#CAC5E8", "#E1D273", "#82DFCE", "#F1AF92", "#E1E7CF", "#92C798", 
    #         "#CDF197", "#DDD199", "#ECB1D1", "#91C7DE", "#E3AF6E", "#ECAAAC", 
    #         "#CEB29C", "#B2BCA9", "#B0C778", "#DBCDD7", "#B5DEE0", "#A5ECB4", 
    #         "#C5E6C0", "#E5CCC0"] 
    colors = ["#688EAF", "#FC991D", "#7DEB74", "#FA6781", "#8B981D", "#BB7548", 
            "#AD8FE4", "#96E4AA", "#D669B0", "#E1C947", "#A78200", "#7C9FE4", 
            "#957DA6", "#75BF38", "#C3B059", "#51C17A", "#79AEBB", "#2790AC", 
            "#688ECE", "#749DB7"]
    for i in range(3):
        colors += colors

    for n,cdata in enumerate(data):
        
        # get the number of lines to skip
        if not isinstance(raw_skip,int):
            skip = int(raw_skip*cdata.shape[0])

        plt.plot(cdata[skip:],marker='s',color=colors[n%len(colors)],markeredgecolor=colors[n%len(colors)],\
                         markersize=4,linestyle='-',linewidth=1.0,label=leglabel[n])

    plt.ylabel(yLong)
    plt.xlabel("MC Bin Number")
    if plttitle:
        plt.title(plttitle)
    if not args['--nolegend']:
        leg = plt.legend(loc='best', frameon=False, prop={'size':16},markerscale=2, ncol=2)
        for l in leg.get_lines():
            l.set_linewidth(4.0)
    if args['--savefig']:
        plt.savefig('1-' + args['--savefig'],bbox_inches='tight')

    # ============================================================================
    # Figure 2 : running average of column vs. MC Bins
    # ============================================================================
    plt.figure(2)

    n = 0
    for n,cdata in enumerate(data):
        if np.size(cdata) > 1:

            # get the number of lines to skip
            if not isinstance(raw_skip,int):
                skip = int(raw_skip*cdata.shape[0])

            # Get the cumulative moving average
            if args['--error']:
                cma = cumulativeMovingAverage(cdata[skip:])
                sem = error*np.ones_like(cma)
            elif args['--nobin']:
                cma,sem = cumulativeMovingAverageWithError(cdata[skip:])
            else:
                cma = cumulativeMovingAverage(cdata[skip:])
                ave,err = getStats(cdata[skip:])
                sem = err*np.ones_like(cma)
                if not args['--quiet']:
                    print('%s:  %s = %8.4E +- %8.4E' % (leglabel[n], yShort, ave, err)) 

            sma = simpleMovingAverage(50,cdata[skip:])
            x = range(int(0.10*len(cma)),len(cma))
            plt.plot(x,cma[x],color=colors[n%len(colors)],linewidth=1.0,marker='None',linestyle='-',
                label=leglabel[n])
            plt.fill_between(x, cma[x]-sem[x], cma[x]+sem[x],color=colors[n%len(colors)], alpha=0.1)
            n += 1

    # Add a possible horizontal line indicating some value
    if args['--hline']:
        plt.axhline(y=val,color='gray',linewidth=2.5, label=hlabel)

    plt.ylabel(yLong)
    plt.xlabel("MC Bin Number")
    if plttitle:
        plt.title(plttitle)

    if not args['--nolegend']:
        leg = plt.legend(loc='best', frameon=False, prop={'size':16},markerscale=2, ncol=2)
        for l in leg.get_lines():
            l.set_linewidth(4.0)

    if args['--savefig']:
        plt.savefig('2-' + args['--savefig'],bbox_inches='tight')

    # Perform a Welch's t-test
    if args['--ttest']:
        # We only perform the Welch's t test if we have multiple samples we are
        # comparing
        N = len(data)
        if N > 1:
            tval = np.zeros([N,N])
            p = np.zeros([N,N])

            for i in range(N):

                # get the number of lines to skip
                if not isinstance(raw_skip,int):
                    skip = int(raw_skip*data[i].shape[0])

                for j in range(i+1,N):
                    tval[i,j],p[i,j] = stats.ttest_ind(data[i][skip:], data[j][skip:], 
                                               equal_var=False)

        # ============================================================================
        # Figure 3 : plot the estimator histogram along with t-test values
        # ============================================================================
        fig = plt.figure(3)
        # connect('key_press_event',kevent.press)
        for i in range(N):
            n, bins, patches = plt.hist(data[i], 100, normed=True, facecolor=colors[i], 
                                    alpha=0.75, label=leglabel[i],
                                    edgecolor='w')
        # Add the p-values from the t-test
        y = 0.92
        if N > 1:
            plt.figtext(0.78, y, 't-test p values', horizontalalignment='center', 
                 verticalalignment='top', fontsize=15, backgroundcolor='white')

            for i in range(N):
                for j in range(i+1,N):
                    y -= 0.03
                    lab = 'p(' + leglabel[i] + ' - ' + leglabel[j] + ') = ' + '%4.2f'%p[i,j]
                    plt.figtext(0.78, y, lab, horizontalalignment='center', 
                         verticalalignment='top', fontsize=12,
                         backgroundcolor='white')

        if not args['--nolegend']:
            plt.legend(loc='upper left', fontsize=15, frameon=False)
        plt.xlabel(yLong)
        plt.ylabel(r'$P($' + estimator + r'$)$')
        if plttitle:
            plt.title(plttitle)
        if args['--savefig']:
            plt.savefig('3-' + args['--savefig'],bbox_inches='tight')
    if not args['--savefig']:
        plt.show()
    
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
