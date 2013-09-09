import pylab as pl
import glob,argparse

#def parseCMD():
#    ''' parse the command line. '''
#    parser = argparse.ArgumentParser(description='pulls down lots of files.')
#    parser.add_argument('fileName', help='Data File Name.')
#    parser.add_argument('-t', '--typeOfAverage', type=str,
#            default='jackknife', 
#            help='NOT WORKING YET: Do you want jackknife or bootstrap?')
#
#    return parser.parse_args()

def jackknife(data,data2=None,data3=None):
    ''' Return jackknife average (accounting for bias) and error.'''
    numBins = len(data)
    jkTerms = pl.array([])

    if data2==None:     # in case of single 
        for i in range(numBins):
            newsample = data[:i]
            newsample = pl.append(newsample, data[i+1:])
            jkTerms = pl.append(jkTerms,pl.average(newsample))
        dataAve = pl.average(data)
    else:               # in case of specific heat
        for i in range(numBins):
            newsamp1, newsamp2, newsamp3 = data[:i], data2[:i], data3[:i]
            newsamp1 = pl.append(newsamp1,data[i+1:])
            newsamp2 = pl.append(newsamp2,data2[i+1:])
            newsamp3 = pl.append(newsamp3,data3[i+1:])
            jkTerms = pl.append(jkTerms, pl.average(newsamp1)
                    - pl.average(newsamp2)**2 - pl.average(newsamp3))
        dataAve = ( pl.average(data) - pl.average(data2)**2 
                - pl.average(data3) )
    
    jkAve = pl.average(jkTerms)
    ActAve = 1.0*numBins*dataAve - 1.0*(numBins-1)*jkAve
    jkVar = pl.average(jkTerms**2) - jkAve**2
    jkErr = pl.sqrt((numBins-1)*jkVar)

    return ActAve, jkErr

# =============================================================================
