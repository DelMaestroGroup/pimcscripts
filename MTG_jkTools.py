import numpy as np
import glob,argparse

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='pulls down lots of files.')
    parser.add_argument('fileNames', help='Data File Name.', nargs='+')
    parser.add_argument('-t', '--typeOfAverage', type=str,
            default='jackknife',
            help='NOT WORKING YET: Do you want jackknife or bootstrap?')
    parser.add_argument('-s', '--skip', type=int,
            default=1000,
            help='Number of bins to skip')
    parser.add_argument('-c', '--crunched', action='store_true', 
            dest='Crunched', default=False,
            help='Is Cv data already crunched from multiple seeds?')

    return parser.parse_args()

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

# =============================================================================
