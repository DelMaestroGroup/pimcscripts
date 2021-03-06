# plotoptions.py
# Adrian Del Maestro
# 12.15.2011
# 
# Options that will be used when plotting vector and scalar estimators

import loadgmt

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
    #colors  = loadgmt.getColorList('cb/div','Spectral_08',numColors)
    #colors  = loadgmt.getColorList('cb/div','PiYG_07',numColors)
    colors  = loadgmt.getColorList('cb/qual','Set1_09',numColors)

    # colors.reverse()
    markers = loadgmt.getMarkerList()

    return markers,colors

