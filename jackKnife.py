import pylab as pl
import glob,argparse
import averagingTools as aTools

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='pulls down lots of files.')
    parser.add_argument('fileNames', help='Data File Name.', nargs='+')
    parser.add_argument('-t', '--typeOfAverage', type=str,
            default='jackknife', 
            help='NOT WORKING YET: Do you want jackknife or bootstrap?')
    parser.add_argument('-s', '--skip', type=int,
            default=100,
            help='Number of bins to skip')

    return parser.parse_args()

def main():

    args = parseCMD()
    fileNames = args.fileNames
    skip = args.skip
    #data = pl.loadtxt(fileName, unpack=True, usecols=[colNum])
    temps = pl.array([])
    Cvs = pl.array([])
    CvsErr = pl.array([])
    for fileName in fileNames:
        temps = pl.append(temps, float(fileName[-40:-34]))
        EEcv, Ecv, dEdB = pl.loadtxt(fileName, unpack=True, usecols=(11,12,13))
        print len(EEcv[skip:])
        jkAve, jkErr = aTools.jackknife(EEcv[skip:],Ecv[skip:],dEdB[skip:])
        print '<est> = ',jkAve,' +/- ',jkErr
        Cvs = pl.append(Cvs,jkAve)
        CvsErr = pl.append(CvsErr,jkErr)

    # plot specific heat for QHO
    tempRange = pl.arange(0.01,1.0,0.01)
    CvAnalytic = 1.0/(4.0*(tempRange*pl.sinh(1.0/(2.0*tempRange)))**2)

    pl.plot(tempRange,CvAnalytic, label='Exact')
    pl.errorbar(temps,Cvs,CvsErr, label='PIMC',color='Violet',fmt='o')
    pl.xlabel('Temperature [K]')
    pl.ylabel('Specific Heat')
    pl.title('1D QHO -- 1 boson')
    pl.legend(loc=2)
    pl.show()




# =============================================================================
if __name__=='__main__':
    main()
