#!/usr/bin/python
"""pimcplot

Description:
  Performs a cumulative average plot of raw Monte Carlo data

Usage:
    pimcplot.py [options] [--legend=<label>...] --estimator=<name>... <file>...

    pimcplot.py -h | --help 

Options:
  -h, --help                    Show this screen.
  --estimator=<name>, -e <name> The estimator to be plotted.
  --skip=<n>, -s <n>            Number of measurements to be skipped [default: 0].
  --period=<m>, -p <m>          The period of the average window [default: 50].
  --legend=<label>, -l <label>  A legend label
  --error=<units>, -d           Size of the error bars
  --bin                         Use the binned errorbars
  --ttest                       Perform a ttest
  --pdf                         Save PDF file of each plot
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

from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

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
    estimators = args['--estimator']
    leglabel = args['--legend'] and args['--legend']
    error = args['--error'] and float(args['--error'])

    # number of estimators must equal number of filenames given
    if len(fileNames) != len(estimators):
        print "\nNumber of estimators must equal number of filenames supplied."
        print "Note that these must be in the same order!\n"
        sys.exit()
    
    # if labels are not assigned, we default to the PIMCID
    if not leglabel:
        leglabel = []
        for n,fileName in enumerate(fileNames):
            leglabel.append(fileName[-13:-4])
    
    numFiles = len(fileNames)
    
    # first load all of data from all files into a list
    col = []
    data = []
    for numEst in range(len(estimators)):
        headers = pimchelp.getHeadersDict(fileNames[numEst])
        col = list([headers[estimators[numEst]]])

        dataFile = open(fileNames[numEst],'r');
        dataLines = dataFile.readlines();
        dataFile.close()
    
        if len(dataLines) > 2:
            data.append(loadtxt(fileNames[numEst],usecols=col,unpack=True))

    # loop over files, loading them and plotting all data wanted.
    for numEst in range(len(estimators)):
     
        # We count the number of lines in the estimator file to make sure we have
        # some data and grab the headers
        headers = pimchelp.getHeadersDict(fileNames[numEst])
       
        # If we don't choose an estimator, provide a list of possible ones
        if estimators[numEst] not in headers:
            errorString = "Need to specify one of:\n"
            for head,index in headers.iteritems():
                errorString += "\"%s\"" % head + "   "
            parser.error(errorString)

        # Attempt to find a 'pretty name' for the label, otherwise just default to
        # the column heading
        label = pimchelp.Description()
        try:
            yLong = label.estimatorLongName[estimators[numEst]]
        except:
            yLong = estimators[numEst]
        try:
            yShort = label.estimatorShortName[estimators[numEst]]
        except:
            yShort = estimators[numEst]

        # ============================================================================
        # Figure 1 : column vs. MC Steps
        # ============================================================================
        if args['--pdf']:
            figure(1, dpi=40, figsize=(6,3.8))
        else:
            figure(1)
        connect('key_press_event',kevent.press)

        #colors  = loadgmt.getColorList('cw/1','cw1-029',max(numFiles,2))
        #colors  = loadgmt.getColorList('cw/1','cw1-013',max(numFiles,2))
        colors  = loadgmt.getColorList('oc','rainbow',max(numFiles,2))
    #    oc/rainbow

        for n,cdata in enumerate(data):
        #for cdata in data:
            #print n
            print colors[n]
            plot(cdata[skip:],marker='s',color=colors[n],markeredgecolor=colors[n],\
                 markersize=4,linestyle='-',linewidth=1.0, label=leglabel[n])

        ylabel(yLong)
        xlabel("MC Bin Number")
        if numEst == 0:
            legend(loc=3, frameon=False, ncol=2)

        if args['--pdf']:
            savefig('col_vs_MCSteps.pdf', format='pdf',
                    bbox_inches='tight')
            savefig('col_vs_MCSteps_trans.pdf', format='pdf',
                    bbox_inches='tight', transparent=True, dpi=40)
        else:
            savefig('col_vs_MCSteps_trans.png', format='png',
                    bbox_inches='tight', transparent=True)


        # ============================================================================
        # Figure 2 : running average of column vs. MC Bins
        # ============================================================================
        if args['--pdf']:
            figure(2, dpi=40, figsize=(6,3.8))
        else:
            figure(2, figsize=(6,3.8))
        connect('key_press_event',kevent.press)

        n = 0
        for n,cdata in enumerate(data):
            
            if size(cdata) > 1:
                
                # Get the cumulative moving average
                if args['--error']:
                    cma = cumulativeMovingAverage(cdata[skip:])
                    sem = error*ones_like(cma)
                elif args['--bin']:
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
                #n += 1

        ylabel(yLong)
        xlabel("MC Bin Number")
        
        if numEst == 0:
            leg = legend(loc='best', frameon=False, prop={'size':12},markerscale=2, ncol=2)
        for l in leg.get_lines():
            l.set_linewidth(4.0)

        if args['--pdf']:
            savefig('runAve_vs_MCSteps.pdf', format='pdf',
                    bbox_inches='tight')
            savefig('runAve_vs_MCSteps_trans.pdf', format='pdf',
                    bbox_inches='tight',transparent=True, dpi=40)
        else:
            savefig('runAve_vs_MCSteps_trans.png', format='png',
                    bbox_inches='tight',transparent=True)

        # Perform a Welch's t-test
        if args['--ttest']:
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
            if args['--pdf']:
                fig = figure(3, dpi=40, figsize=(6,3.8))
            else:
                fig = figure(3)
            connect('key_press_event',kevent.press)
            for i in range(N):
                n, bins, patches = hist(data[i], 100, normed=True, facecolor=colors[i], 
                        alpha=0.75, label=leglabel[i],
                        edgecolor='w')
            '''# Add the p-values from the t-test
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
                           backgroundcolor='white')'''

            #legend(loc='upper left', fontsize=15, frameon=False)
            if numEst == 0:
                legend(loc='upper left', frameon=False)
            xlabel(yLong)
            yLabb = estimators[numEst]
            if yLabb == 'Ecv/N':
                yLabb = 'E/N'
            ylabel(r'$P($' + yLabb + r'$)$')


            if args['--pdf']:
                savefig('ttest_histogram.pdf', format='pdf',
                        bbox_inches='tight')
                savefig('ttest_histogram_trans.pdf', format='pdf',
                        bbox_inches='tight', transparent=True, dpi=40)
            else:
                savefig('ttest_histogram_trans.png', format='png',
                        bbox_inches='tight', transparent=True)
    show()
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
