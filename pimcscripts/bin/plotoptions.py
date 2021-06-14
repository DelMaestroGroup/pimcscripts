# plotoptions.py
# Adrian Del Maestro
# 12.15.2011
# 
# Options that will be used when plotting vector and scalar estimators

from matplotlib import cm
from matplotlib import colors as mplcolors
import numpy as np

# -----------------------------------------------------------------------------
def plotOptions(plotType):
    ''' Consctruct a dictionary of plot options based on what type of plot is
    desired.'''

    # First define the common properties
    pOps = {}

    # marker properties
    pOps['markersize'] = 7
    pOps['markeredgewidth'] = 0.75
    pOps['markeredgecolor'] = '#4d4d4d'
    pOps['marker'] = 'o'

    # line properties
    pOps['color'] = 'black'
    pOps['linewidth'] = 0.75

    if 'l' in plotType:
        pOps['linestyle'] = '-' 
    else:
        pOps['linestyle'] = 'None' 

    if plotType == 'p':
        pOps['linestyle'] = 'None'

    if plotType == 'l':
        pOps['linewidth'] = 3.0
        pOps['marker'] = None
        pOps['markeredgewidth'] = 0.0
        pOps['markersize'] = 0.0
        pOps.pop('color')

    if 'e' in plotType:
        pOps['capsize'] = 4
        
    return pOps


# -----------------------------------------------------------------------------
def markersColors(numColors):
    '''return the markers and colors used for plotting.'''

    # He4 Pore PRL Colors
#    colors = ['#556270','#1C3249','#4ECDC4','#19857D','#C7F464','#FF6B6B','#A62323'] #556270

    # http://www.graphviz.org/content/color-names
    if numColors == 1:
        numColors+=1;

    color_map = 'Spectral_r'
    cmap = cm.get_cmap(color_map)

    colors = []
    for n in np.linspace(0,1.0,numColors):
        colors.append(mplcolors.to_hex(cmap(n)))

    markers = ['o','s','^','v','>','<','p','d','*','H','$\clubsuit$',
           '$\spadesuit$', 'x','$\otimes$','1','2','3','4']

    return markers,colors

