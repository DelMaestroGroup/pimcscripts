"""Plot one or more reduced scalar estimators.
Author: Adrian Del Maestro
Date: 11.29.2011

Description:
    Plots a scalar estimator as a function of the reduced variable for one or
    more fixed parameters

usage: 
    plot_scalar_estimator.py --estimator=<E>... [--legend=<L>...] [--xlim=<x>...] [--ylim=<y>...] [options] <reduce-files>... 

    plot_scalar_estimator.py -h | --help

Positional arguments:
  reduce-files              Reduced scalar estimator files.

Options:
  -h, --help                Show this help message and exit

  --estimator=<E>, -e <E>  The estimators to be plotted.
  --plot_type=<P>, -p <P>  The plot type, a combination of l = lines,
                           p = points, e = errorbars, f = filled, [default: lp].
  --ndim=<d>, -d <d>       Number of spatial dimensions [default: 3].
  --legend=<L>, -l <L>     Legend key for multiple curves.
  --xlim=<x>, -x=<x>       X-axis limits.
  --ylim=<y>, -y=<y>       Y-axis limits. 
"""

from docopt import docopt
import pyutils
import pimchelp
import numpy as np
import pylab as pl
import plotoptions

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # Get the command line arguments
    args = docopt(__doc__)

    fileNames = args['<reduce-files>']
    estimatorToPlot = args['--estimator']
    NDIM = int(args['--ndim'])
    plotType = args['--plot_type']

    if args['--legend']:
        varLabel = [l for l in args['--legend']]
    else: 
        varLabel = [None]

    xLim = args['--xlim'] and [float(x) for x in args['--xlim']]
    yLim = args['--ylim'] and [float(y) for y in args['--ylim']]

    # Analyze the imput files
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
        pl.figure(n+1)
        ax = pl.subplot(111)
        for i in range(reduce.getNumVarParams()):
            lab = reduce.getVarLabel(i)
            if not lab:
                lab = varLabel[i]
            pOptions['color'] = colors[i]

            # fill between to represent error bars
            if 'f' in plotType:
                emin = reduce.estimator(est,i)-reduce.estimatorError(est,i)
                emax = reduce.estimator(est,i)+reduce.estimatorError(est,i)
                pl.fill_between(reduce.param(), emax,emin, color=colors[i],
                                alpha=0.2)
            # plot with points
            if 'p' in plotType:
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        markerfacecolor=colors[i],label=lab, **pOptions)

            # plot lines only
            if plotType == 'l':
                pl.plot(reduce.param(), reduce.estimator(est,i),
                        label=lab,**pOptions)
            # plot errorbars
            elif 'e' in plotType:
                eb = pl.errorbar(reduce.param(), reduce.estimator(est,i),
                        yerr=reduce.estimatorError(est,i), 
                        markerfacecolor=colors[i], ecolor=colors[i],
                        label=lab, **pOptions)
                # Set the width of the cap lines
                for cap in eb[1]:
                    cap.set_mew(2.0)

            
        # pl.tight_layout()
        pl.xlabel(descrip.paramLongName[reduce.reduceLabel],labelpad=16)
        pl.ylabel(descrip.estimatorLongName[est])
        leg = pl.legend(frameon=False, loc='best', prop={'size':18})
        if xLim != []:
            pl.xlim(xLim[0],xLim[1])

        if yLim != []:
            pl.ylim(yLim[0],yLim[1])

    pl.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

