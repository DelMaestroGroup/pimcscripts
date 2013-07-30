#!/usr/bin/python
"""pimcplot

Description:
  Plot a scalar estimator.

Usage:
    pimcplot.py [options] [--label=<label>...] [--xLim=<x>...] [--yLim=<y>...] --estimator=<name>... <file>...

    pimcplot.py -h | --help 

Options:
  -h, --help                    Show this screen.
  --estimator=<name>, -e <name> A list of estimator names that are to be plotted.
  --plot=<p>, -p <p>            The plot type.  l=lines, p=points,e=errorbars. [default: e]
  --skip=<n>, -s <n>            Number of measurements to be skipped. Choices: l, e, p, lp, le. [default: 0].
  --period=<m>, -p <m>          The period of the average window [default: 50].
  --label=<label>, -l <label>   A legend label.
  --xLim=<x>, -x <x>            x-axis limits.
  --yLim=<y>, -y <y>            y-axis limits.
  --ndim=<d>, -n <d>            Number of spatial dimensions [default: 3].
  --same                        Plot estimators on the same plot.
"""

# plot_scalar_estimator.py
# Adrian Del Maestro
# 11.29.2011
# 
# Plot a single scalar estimator

import os,sys,loadgmt,kevent,pyutils,pimchelp,plotoptions
import numpy as np
import pylab as pl
from docopt import docopt

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # parse the command line options
    args = docopt(__doc__)

    fileNames = args['<file>']
    estimatorToPlot = args['--estimator']
    NDIM = int(args['--ndim']) 
    plotType = args['--plot']
    varLabel = args['--label']
    xLim = args['--xLim']
    yLim = args['--yLim']

    # don't want any empty lists
    #if varLabel == []:
    #    varLabel = None
    if xLim == []:
        xLim = None
    if yLim == []:
        yLim = None

    if len(fileNames) < 1:
        sys.exit("\nNeed to specify at least one scalar estimator file\n")

    # Analyze the imput files
    #reduce = pimchelp.ScalarReduce(fileNames,varLabel)
    reduce = pimchelp.ScalarReduce(fileNames,None)

    # Get a color scheme and marker list
    numColors = max(reduce.getNumVarParams(),2)
    markers,colors  = plotoptions.markersColors(numColors)

    # create a description object
    descrip = pimchelp.Description(NDIM)

    # get the plot options
    pOptions = plotoptions.plotOptions(plotType)

    # Plot each estimator
    for n,est in enumerate(estimatorToPlot):
        if args['--same']:
            pl.figure(1)
        else:
            pl.figure(n+1)
        ax = pl.subplot(111)
        for i in range(reduce.getNumVarParams()):
            if args['--same']:
                lab = varLabel[n]
            else:
                lab = reduce.getVarLabel(i)
            pOptions['color'] = colors[i]
            if 'p' in plotType:
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        markerfacecolor=colors[i],label=lab, **pOptions)
            if plotType == 'l':
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        label=lab,**pOptions)
            elif 'e' in plotType:
                if args['--same']:
                    eb = pl.errorbar(reduce.param(), reduce.estimator(est,i),
                            yerr=reduce.estimatorError(est,i), 
                            markerfacecolor=colors[i+n], ecolor=colors[i+n],
                            label=lab, **pOptions)
                else:
                    eb = pl.errorbar(reduce.param(), reduce.estimator(est,i),
                            yerr=reduce.estimatorError(est,i), 
                            markerfacecolor=colors[i], ecolor=colors[i],
                            label=lab, **pOptions)
                
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(2.0)

        #pl.tight_layout()
        pl.xlabel(descrip.paramLongName[reduce.reduceLabel])
        pl.ylabel(descrip.estimatorLongName[est])
        leg = pl.legend(frameon=False, loc='best', prop={'size':18})
        if xLim != None:
            pl.xlim(xLim[0],xLim[1])

        if yLim != None:
            pl.ylim(yLim[0],yLim[1])

    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

