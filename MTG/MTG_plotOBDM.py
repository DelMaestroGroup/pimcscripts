# =============================================================================
# Script to read in OBDM files from code and compute Fourier Transform
# in order to get back the BEC fraction.  This loops over a range of temps
# (data files) and plots condensate fraction for each temp.
#
# Author:           Max Graves
# Last Revised:     28-FEB-2013
# =============================================================================

import pylab as pl
import glob
import sys
import os

def main():

    # get all files, sorted by temp
    os.chdir('./data')
    files = sorted(glob.glob("*obdm*") and glob.glob("*dat*"))[::-1]
    if files == []:
        print 'no files located'
        sys.exit()
    os.chdir('..')

    BECfrac = pl.array([])
    Temps = pl.array([])

    # compute zero frequency FFT term for each file
    for fileName in files:

        # get header line
        estFile = open('./data/'+fileName,'r')
        estLines = estFile.readlines();
        pimcid = estLines[0]
        headers = estLines[1].split()
        estFile.close()
        
        headers.pop(0)
        headers = pl.array(headers).astype(float)

        rawData = pl.loadtxt('./data/'+fileName)

        avgs = pl.array([])
        Ncol = 1.0*rawData.size/rawData[0].size
        for j in range(int(Ncol)):
            for i in range(rawData[j].size):
                if j == 0:
                    avgs = pl.append(avgs, rawData[j][i])
                else:
                    avgs[i] += rawData[j][i]
        avgs /= Ncol

        # compute FFT
        f = pl.fft(avgs)

        BECfrac = pl.append(BECfrac, f[0].real)
        Temps = pl.append(Temps, (fileName[9:15]))

        print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
        print 'Zero Frequency FFT term: ',f[0].real
        print '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='

    # plotting
    #pl.plot(headers, avgs, marker='o', linewidth=0, markerfacecolor='None',
    #        markeredgecolor='Indigo')
    #pl.ylabel(r'$\langle \hat{\Psi}^{\dagger}(r)\ \hat{\Psi}(r\') \rangle$', 
    #        fontsize=20)
    pl.plot(Temps, BECfrac, marker='o', linewidth=0, markerfacecolor='None',
            markeredgecolor='Indigo')
    pl.ylabel('Condensate Fraction', 
            fontsize=20)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    pl.title('BEC Fraction')
    pl.grid(True)
    pl.show()

# =============================================================================            
if __name__=='__main__':
    main()
