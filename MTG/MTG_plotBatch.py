import pylab as pl
import sys,glob

def main():

    fileNames = glob.glob('*batch*')

    E,Eerr = pl.array([]),pl.array([])
    Ecv,Ecverr = pl.array([]),pl.array([])
    Cv,Cverr = pl.array([]),pl.array([])
    temps = pl.array([])

    for f in sorted(fileNames):
        EF,EerrF,EcvF,EcverrF,Cv1F,Cv1errF,Cv3F,Cv3errF = pl.loadtxt(f, 
                usecols=(6,7,16,17,32,33,34,35), unpack=True)
        temps = pl.append(temps, float(f[13:19]))
        E = pl.append(E, EF[-1])
        Ecv = pl.append(Ecv, EcvF[-1])
        Eerr = pl.append(Eerr, EerrF[-1])
        Ecverr = pl.append(Ecverr, EcverrF[-1])
        Cv = pl.append(Cv, Cv1F[-1] - EcvF[-1]**2 - Cv3F[-1])
        Cverr = pl.append(Cverr, Cv1errF[-1] - EcverrF[-1]**2 - Cv3errF[-1])

    pl.figure(1)
    pl.errorbar(temps,E,Eerr,fmt='o', label='Therm')
    pl.errorbar(temps,Ecv,Ecverr,fmt='o', label='Vir')
    pl.ylabel('Energy [K]')
    pl.xlabel('Temperature [K]')
    pl.grid(True)
    pl.legend(loc=2)
 
    pl.figure(2)
    pl.errorbar(temps,Cv,Cverr, fmt='o')
    pl.ylabel('Specific Heat [K]')
    pl.xlabel('Temperature [K]')
    pl.grid(True)
    pl.legend(loc=2)   
    pl.show()

    sys.exit()
    #headers   = pimchelp.getHeadersFromFile(fileNames[0])
    '''
    # compute single centroid virial specific heat if possible
    Cv = ave[:,headers.index('EEcv')] - ave[:,headers.index('Ecv')]**2 - ave[:,headers.index('dEdB')]
    aveNew = zeros([len(fileNames),len(headers)+1],float)
    errNew = zeros([len(fileNames),len(headers)+1],float)
    for i,a in enumerate(ave):
        a = append(a, ave[:,headers.index('EEcv')][i] \
                - ave[:,headers.index('Ecv')][i]**2 \
                - ave[:,headers.index('dEdB')][i])
        aveNew[i] = a
    for i, e in enumerate(err):
        e = append(e, err[:,headers.index('EEcv')][i] \
                - err[:,headers.index('Ecv')][i]**2 \
                - err[:,headers.index('dEdB')][i])
        errNew[i] = e
    headers.append('Cv')
    ave = aveNew
    err = errNew


    E, Eerr, Ecv, Ecverr, Cv,'''

if __name__=='__main__':
    main()
