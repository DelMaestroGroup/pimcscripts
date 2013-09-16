import numpy as np
import glob,argparse

def jackknife(data,data2=None,data3=None):
    ''' Return jackknife average (accounting for bias) and error.'''
    numBins = int(len(data))
    jkTerms = np.zeros(numBins)

    if data2==None:     # in case of single 
        for i in range(numBins):
            jkTerms[i] = np.mean(np.delete(data,i))
        dataAve = np.mean(data)
    else:               # in case of specific heat
        for i in range(numBins):
            jkTerms[i] = ( np.mean(np.delete(data,i))
                    - np.mean(np.delete(data2,i))**2 
                    - np.mean(np.delete(data3,i)) )
        dataAve = ( np.mean(data) - np.mean(data2)**2 
                - np.mean(data3) )
    
    jkAve = np.mean(jkTerms)
    ActAve = 1.0*numBins*dataAve - 1.0*(numBins-1)*jkAve
    jkVar = np.mean(jkTerms**2) - jkAve**2
    jkErr = np.sqrt((numBins-1)*jkVar)

    return ActAve, jkErr

# =============================================================================
