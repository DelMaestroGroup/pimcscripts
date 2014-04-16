# =============================================================================
# Functions used in jackknife data analysis scripts.
#
# Author:           Max Graves
# Last Revised:     18-JAN-2014
# =============================================================================

import numpy as np
import glob,argparse

def getHeadersFromFile(fileName, skipLines=0):
    ''' 
    Get the data column headers (temperatures).
    This assumes the headers are on the top line of the file.
    '''
    with open(fileName,'r') as inFile:
        inLines = inFile.readlines()
        temps = inLines[0].split()
        temps.pop(0)

    return temps


def getXlabel(reduceType):
    ''' given reduceType, return x axis label as string'''
    if reduceType == 'T':
        xLab = 'Temperature [K]'
    elif reduceType == 'u':
        xLab = 'Chemical Potential [K]'
    
    return xLab


def jackknife(data,data2=None,data3=None):
    ''' 
    Return jackknife average (accounting for bias) and error.
    This can be passed either a single numpy array or three
    numpy arrays.  If passed a single numpy array it will simply
    return mean and error for that array.  If passed three, it 
    is expected that these are the three arrays for specific heat.
    '''
    if data2!=None or data3!=None:  
        Cv=True
    else:
        Cv=False
    
    numBins = int(len(data))
    jkTerms = np.zeros(numBins)

    # compute total sums first
    totSum1 = 0.0
    if Cv:
        totSum2 = 0.0
        totSum3 = 0.0 
    for i in range(numBins):
        totSum1 += data[i]
        if Cv:
            totSum2 += data2[i]
            totSum3 += data3[i]

    # create resampled arrays
    for i in range(numBins):
        t1 = (totSum1 - data[i])/(1.0*numBins-1.0)
        if Cv:
            t2 = (totSum2 - data2[i])/(1.0*numBins-1.0)
            t3 = (totSum3 - data3[i])/(1.0*numBins-1.0)
            jkTerms[i] = t1 - t2*t2 - t3
        else:
            jkTerms[i] = t1

    # compute raw average (no resampling)
    if Cv:
        rawAve = np.mean(data) - (np.mean(data2))**2 - np.mean(data3)
    else:
        rawAve = np.mean(data)

    # compute mean and standard error
    jkAve = np.mean(jkTerms)
    ActAve = 1.0*numBins*rawAve - 1.0*(numBins-1)*jkAve
    jkVar = np.mean(jkTerms**2) - jkAve**2
    jkErr = np.sqrt((numBins-1)*jkVar)

    return ActAve, jkErr


def parseCMD():
    ''' Parse the command line. '''
    descLine = ('Script for analyzing data collected from multiple random'
            +' number seeds.')
    parser = argparse.ArgumentParser(description=descLine)
    
    # help strings.
    fNameHelp = 'Data File Names.  For checkStandardError.py, this will be \
            the directory name that holds all direcs with zAveraged files.'
    nColHelpStr = ('Column number of first occurrence of the type of data you '+
            'want to analyze from zAveraged file. Used for checkStandardError.')
    
    # parser arguments.
    parser.add_argument('fileNames', help=fNameHelp, nargs='+')
    parser.add_argument('-s', '--skip', type=int,
            default=0,
            help='Number of bins to skip')
    parser.add_argument('-c', '--crunched', action='store_true', 
            dest='Crunched', default=False,
            help='Is the data already crunched from multiple seeds?')
    parser.add_argument('-R', '--randomColors', action='store_false', 
            dest='RandomColors', default=True,
            help='Include this to NOT randomly shuffle list of plot colors.')
    parser.add_argument('-r', '--reduceType', type=str,
            default='T',
            help='What are we reducing over? [T,u]')
    parser.add_argument('-n', '--nCol', type=int,
            default=1, 
            help=nColHelpStr)
    parser.add_argument('-e', '--nEst', type=int,
            default=1, 
            help='Number of different types of estimators in data file.')   
    parser.add_argument('-S', '--bulkSeparation', type=int,
            default=1, 
            help='Enter the bulk separation.  Used for Lz scaling.')
   
    return parser.parse_args()


def scalingVariable(direcs):
    ''' 
    determine titles for plotting zAveragedNtwind.dat files
    '''
    for nd, d in enumerate(direcs):
        # determine titles for plotting
        if ('Lz' in d):
            scalingVar = 'Lz'
        else:
            extent = d[12:16]

            extentScale = False
            if nd==0:
               firstExtent = extent 
            if extent == firstExtent:
                scalingVar = 'extent'
            else:
                scalingVar = 'thickness'
    
    return scalingVar


# =============================================================================
