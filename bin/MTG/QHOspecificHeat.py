# =============================================================================
# 1D QHO Energy/Specific Heat Plotting Script.
# Note: This is an adapted version of jackKnife.py.
#
# Author:           Max Graves
# Last Revised:     29-JAN-2014
# =============================================================================

import pyximport; pyximport.install()
import pylab as pl
import glob,argparse,sys
import MTG_jkTools as jk
from matplotlib import rcParams

# set up latex fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

def parseCMD():
    ''' Parse the command line. '''
    desc = 'Takes (g)ce-estimator files over range of temperatures and \
            creates new file with jackknife averages for both energy and \
            specific heat.  Then plots these along with the analytic \
            results for the 1D qho with one particle.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fileNames', help='Data File Name.', nargs='+')
    parser.add_argument('-s', '--skip', type=int,
            default=1000,
            help='Number of bins to skip')

    return parser.parse_args()

def main():

    args = parseCMD()
   
    # Check if our data file exists, if not: write one.
    # Otherwise, open the file and plot.
    check = glob.glob('*JackKnifeData_Cv.dat*')
    
    if check == []:
        
        fileNames = args.fileNames
        skip = args.skip
        
        temps,Cvs,CvsErr = pl.array([]),pl.array([]),pl.array([])
        Es, EsErr = pl.array([]),pl.array([])
        ET, ETerr = pl.array([]),pl.array([])

  
        # open new data file, write headers
        fout = open('JackKnifeData_Cv.dat', 'w')
        fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n'% (
            'T', 'Ecv', 'EcvErr', 'Et', 'EtErr','Cv', 'CvErr'))

        # perform jackknife analysis of data, writing to disk
        for fileName in fileNames:
            temp = float(fileName[-40:-34])
            temps = pl.append(temps, temp)

            # grab and analyze energy data
            Ecv,Eth = pl.loadtxt(fileName, unpack=True, usecols=(4,-5))
            E = pl.average(Ecv)
            Et = pl.average(Eth)
            EErr = pl.std(Ecv)/pl.sqrt(float(len(Ecv)))
            EtErr = pl.std(Eth)/pl.sqrt(float(len(Eth)))
            Es = pl.append(Es,E)
            ET = pl.append(ET, Et)
            EsErr = pl.append(EsErr, EErr)
            ETerr = pl.append(ETerr, EtErr)
            
            # grab and analyze specific heat data
            EEcv, Ecv, dEdB = pl.loadtxt(fileName, unpack=True, usecols=(11,12,13))
            jkAve, jkErr = jk.jackknife(EEcv[skip:],Ecv[skip:],dEdB[skip:])
            Cvs = pl.append(Cvs,jkAve)
            CvsErr = pl.append(CvsErr,jkErr)
            
            fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                temp, E, EErr, Et, EtErr, jkAve, jkErr))
            print 'T = ',str(temp),' complete.'
        
        fout.close()

    else:
        print 'Found existing data file in CWD.'
        temps,Es,EsErr,ET,ETerr,Cvs,CvsErr = pl.loadtxt(
                'JackKnifeData_Cv.dat', unpack=True)

    # plot specific heat for QHO
    tempRange = pl.arange(0.01,1.0,0.01)
    Eanalytic = 0.5/pl.tanh(1.0/(2.0*tempRange))
    CvAnalytic = 1.0/(4.0*(tempRange*pl.sinh(1.0/(2.0*tempRange)))**2)

    pl.figure(1)
    ax1 = pl.subplot(211)
    pl.plot(tempRange,CvAnalytic, label='Exact')
    pl.errorbar(temps,Cvs,CvsErr, label='PIMC',color='Violet',fmt='o')
    pl.ylabel('Specific Heat',fontsize=20)
    pl.title('1D QHO -- 1 boson',fontsize=20)
    pl.legend(loc=2)

    pl.setp(ax1.get_xticklabels(), visible=False)
    ax2 = pl.subplot(212, sharex=ax1)
    pl.plot(tempRange,Eanalytic, label='Exact')
    pl.errorbar(temps,Es,EsErr, label='PIMC virial',color='Lime',fmt='o')
    pl.errorbar(temps,ET,ETerr, label='PIMC therm.',color='k',fmt='o')
    #pl.scatter(temps,Es, label='PIMC')
    pl.xlabel('Temperature [K]',fontsize=20)
    pl.ylabel('Energy [K]',fontsize=20)
    pl.legend(loc=2)

    pl.savefig('1Dqho_largerCOM_800000bins_CvANDenergy.pdf',
            format='pdf', bbox_inches='tight')

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
