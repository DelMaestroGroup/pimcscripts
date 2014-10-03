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
  --title=<t>, -t=<t>      A plot title.
  --hline=<h>              A horizontal line indicating a value.
  --hlabel=<hl>            A legend label for the horizontal line
"""

from docopt import docopt
import pyutils
import pimchelp
import numpy as np
import matplotlib.pyplot as plt
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
    title = args['--title']
    hval = args['--hline'] and float(args['--hline'])

    if args['--hline']:
        if args['--hlabel']:
            hlabel = args['--hlabel']
        else:
            hlabel = '%s = %5.2f' % (estimatorToPlot[0],hval)

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
        plt.figure(n+1)
        ax = plt.subplot(111)
        if title:
            plt.title(title)

        # Add a possible horizontal line indicating some value
        if args['--hline']:
            plt.axhline(y=hval,color='#4d4d4d',linewidth=2.0,label=hlabel)

        for i in range(reduce.getNumVarParams()):
            lab = reduce.getVarLabel(i)
            if not lab or varLabel:
                lab = varLabel[i]
            pOptions['color'] = colors[i]

            # fill between to represent error bars
            if 'f' in plotType:
                emin = reduce.estimator(est,i)-reduce.estimatorError(est,i)
                emax = reduce.estimator(est,i)+reduce.estimatorError(est,i)
                plt.fill_between(reduce.param(), emax,emin, color=colors[i],
                                alpha=0.2)
            # plot with points
            if 'p' in plotType:
                plt.plot(reduce.param(), reduce.estimator(est,i),
                        markerfacecolor=colors[i],label=lab, **pOptions)

            # plot lines only
            if plotType == 'l':
                plt.plot(reduce.param(), reduce.estimator(est,i),
                        label=lab,**pOptions)
            # plot errorbars
            elif 'e' in plotType:
                eb = plt.errorbar(reduce.param(), reduce.estimator(est,i),
                        yerr=reduce.estimatorError(est,i), 
                        markerfacecolor=colors[i], ecolor=colors[i],
                        label=lab, **pOptions)

                # Set the width of the cap lines
                # for cap in eb[1]:
                #     cap.set_mew(2.0)

            
        plt.xlabel(descrip.paramLongName[reduce.reduceLabel],labelpad=16)
        plt.ylabel(descrip.estimatorLongName[est])
        leg = plt.legend(frameon=False, loc='best', prop={'size':18})
        if xLim != []:
            plt.xlim(xLim[0],xLim[1])

        if yLim != []:
            plt.ylim(yLim[0],yLim[1])

    plt.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

