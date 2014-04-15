import pylab as pl
import glob,argparse,sys
import averagingTools as aTools

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='pulls down lots of files.')
    parser.add_argument('fileNames', help='Data File Name.', nargs='+')
    parser.add_argument('-t','--timeStepScaling',action='store_true',
            dest='timeStepScaling',default=False,
            help='Do you want to do a time step scaling?')
    parser.add_argument('-s', '--skip', type=int,
            default=1000,
            help='Number of bins to skip')

    return parser.parse_args()

def main():

    args = parseCMD()
    skip = args.skip
   
    # make array of energies
    Es, EsErr = pl.array([]),pl.array([])
    ET, ETerr = pl.array([]),pl.array([])

    # array of <x^2>
    x2s, x2sErr = pl.array([]),pl.array([])

    # temperature array
    temps = pl.array([])

    for fileName in args.fileNames:
        temps = pl.append(temps, float(fileName[-40:-34]))
        Ecv,Eth,x2 = pl.loadtxt(fileName, unpack=True, usecols=(4,-6,-5))
        Ecv,Eth,x2 = Ecv[skip:], Eth[skip:], x2[skip:]
        Es = pl.append(Es,pl.average(Ecv))
        ET = pl.append(ET, pl.average(Eth))
        x2s = pl.append(x2s,pl.average(x2))
        EsErr = pl.append(EsErr, pl.std(Ecv)/pl.sqrt(float(len(Ecv))))
        ETerr = pl.append(ETerr, pl.std(Eth)/pl.sqrt(float(len(Eth))))
        x2sErr = pl.append(x2sErr, pl.std(x2s)/pl.sqrt(float(len(x2s))))
        print 'Finished for temp: ',str(fileName[-40:-34])
    
    # plot specific heat for QHO
    tempRange = pl.arange(0.01,1.0,0.01)
    Eanalytic = 0.5/pl.tanh(1.0/(2.0*tempRange))

    pl.figure(1)
    pl.plot(tempRange,Eanalytic, label='Exact')
    pl.errorbar(temps,Es,EsErr, label='PIMC virial',color='Lime',fmt='o')
    pl.errorbar(temps,ET,ETerr, label='PIMC therm.',color='k',fmt='o')
    pl.xlabel('Temperature [K]')
    pl.ylabel('Energy [K]')
    pl.title('1D QHO -- 1 boson')
    pl.legend(loc=2)

    pl.figure(2)
    pl.plot(tempRange,Eanalytic, label='Exact')
    pl.errorbar(temps,x2s,x2sErr, label='PIMC therm.',color='k',fmt='o')
    pl.xlabel('Temperature [K]')
    pl.ylabel(r'$\langle x^2 \rangle$', fontsize=20)
    pl.title('1D QHO -- 1 boson')
    pl.legend(loc=2)   

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
