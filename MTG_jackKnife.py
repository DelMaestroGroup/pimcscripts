import pylab as pl
import glob,argparse,sys
import MTG_jkTools as aTools

def main():

    args = aTools.parseCMD()
   
    # Check if our data file exists, if not: write one.
    # Otherwise, open the file and plot.
    check = glob.glob('*JackKnifeData_Cv.dat*')
    fileNames = args.fileNames
    skip = args.skip

    reduceType = args.reduceType

    # check which ensemble
    canonical=True
    if fileNames[0][0]=='g':
        canonical=False

    print fileNames
    
    if check == []:
        
        temps,Cvs,CvsErr = pl.array([]),pl.array([]),pl.array([])
        Es, EsErr   = pl.array([]), pl.array([])
        rhos_rhos, rhos_rhoErr  = pl.array([]), pl.array([])

        filmDenses, bulkDenses      = pl.array([]), pl.array([])
        filmDensErrs, bulkDensErrs    = pl.array([]), pl.array([])
   
        # open energy/ specific heat data file, write headers
        fout = open('JackKnifeData_Cv.dat', 'w')
        if reduceType=='T':
            fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'% (
                'T', 'E', 'Eerr', 'Cv', 'CvErr'))
        elif reduceType=='u':
            fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'% (
                'mu', 'E', 'Eerr', 'Cv', 'CvErr'))
        
        # open superfluid stiffness data file, write headers
        foutSup = open('JackKnifeData_super.dat','w')
        if reduceType == 'T':
            foutSup.write('#%15s\t%16s\t%16s\n'%(
                'T', 'rho_s/rho', 'rho_s/rhoErr'))
        elif reduceType == 'u':
            foutSup.write('#%15s\t%16s\t%16s\n'%(
                'mu', 'rho_s/rho', 'rho_s/rhoErr'))
        
         
        # open superfluid stiffness data file, write headers
        foutDens = open('JackKnifeData_bipart.dat','w')
        if reduceType == 'T':
            foutDens.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'%(
                'T', 'filmDens', 'filmDensErr', 'bulkDens', 'bulkDensErr'))
        elif reduceType == 'u':
            foutDens.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'%(
                'mu', 'filmDens', 'filmDensErr', 'bulkDens', 'bulkDensErr'))

        # perform jackknife analysis of data, writing to disk
        if args.Crunched:   # check if we have combined data
            tempList = aTools.getHeadersFromFile(fileNames[0])
            for temp in tempList:
                temps = pl.append(temps,float(temp))
            n,n2,n3 = 0,0,0
            for fileName in fileNames:
                print '\n\n---',fileName,'---\n'
                for temp in tempList:
                    print n
                    if 'Estimator' in fileName:
                        E, EEcv, Ecv, dEdB = pl.loadtxt(fileName,\
                                unpack=True, usecols=(n,n+1,n+2,n+3), delimiter=',')
                        EAve, Eerr = aTools.jackknife(E[skip:])
                        jkAve, jkErr = aTools.jackknife(
                                EEcv[skip:],Ecv[skip:],dEdB[skip:])
                        if reduceType == 'T':
                            print 'T = ',float(temp),' :'
                        elif reduceType == 'u':
                            print 'mu = ',float(temp),' :'
                        print '<E>  = ',EAve,' +/- ',Eerr
                        print '<Cv> = ',jkAve,' +/- ',jkErr
                        Es      = pl.append(Es, EAve)
                        Cvs     = pl.append(Cvs, jkAve)
                        EsErr   = pl.append(EsErr, Eerr)
                        CvsErr  = pl.append(CvsErr, jkErr)
                        fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                            float(temp), EAve, Eerr, jkAve, jkErr)) 
                        n += 4

                    elif 'Super' in fileName:
                        rhos_rho = pl.loadtxt(fileName, \
                                unpack=True, usecols=(n2,), delimiter=',')
                        superAve, superErr = aTools.jackknife(rhos_rho[skip:])
                        print 'rho_s/rho = ', superAve,' +/- ',superErr
                        rhos_rhos   = pl.append(rhos_rhos, superAve)
                        rhos_rhoErr = pl.append(rhos_rhoErr, superErr)
                        foutSup.write('%16.8E\t%16.8E\t%16.8E\n' %(
                            float(temp), superAve, superErr))
                        n2 += 1

                    elif 'BiPart' in fileName:
                        filmDens, bulkDens  = pl.loadtxt(fileName, \
                                unpack=True, usecols=(n3,n3+1), delimiter=',')
                        filmDensAve, filmDensErr = aTools.jackknife(filmDens[skip:])
                        bulkDensAve, bulkDensErr = aTools.jackknife(bulkDens[skip:])
                        print 'filmDens = ', filmDensAve,' +/- ',filmDensErr
                        print 'bulkDens = ', bulkDensAve,' +/- ',bulkDensErr
                        filmDenses      = pl.append(filmDenses, filmDensAve)
                        filmDensErrs    = pl.append(filmDensErrs, filmDensErr)
                        bulkDenses  = pl.append(bulkDenses, bulkDensAve)
                        bulkDensErrs    = pl.append(bulkDensErrs, bulkDensErr)
                        foutDens.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                            float(temp), filmDensAve, filmDensErr, 
                            bulkDensAve, bulkDensErr))
                        n3 += 2

        else:       # otherwise just read in individual (g)ce-estimator files
            # THIS IS LAGGING BEHIND THE CRUNCHED SECTION
            for fileName in fileNames:
                if canonical: 
                    temp = float(fileName[13:19])
                else:
                    temp = float(fileName[14:20])
                temps = pl.append(temps, temp)
                E, EEcv, Ecv, dEdB = pl.loadtxt(fileName, unpack=True, 
                        usecols=(4,11,12,13))
                jkAve, jkErr = aTools.jackknife(
                        EEcv[skip:],Ecv[skip:],dEdB[skip:])
                EAve, Eerr = aTools.jackknife(E[skip:])
                print 'T = ',temp
                print '<Cv> = ',jkAve,' +/- ',jkErr
                print '<E>  = ',EAve,' +/- ',Eerr
                Es      = pl.append(Es, EAve)
                Cvs     = pl.append(Cvs, jkAve)
                EsErr   = pl.append(EsErr, Eerr)
                CvsErr  = pl.append(CvsErr, jkErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), EAve, Eerr, jkAve, jkErr)) 
        
        fout.close()
        foutSup.close()
        foutDens.close()

    else:
        print 'Found existing data file in CWD.'
        temps, Es, EsErr, Cvs, CvsErr = pl.loadtxt(
                'JackKnifeData_Cv.dat', 
                unpack=True)
        temps, rhos_rhos, rhos_rhoErr = pl.loadtxt(
                'JackKnifeData_super.dat', 
                unpack=True)
        temps, filmDenses, filmDensErrs, bulkDenses, bulkDensErrs = pl.loadtxt(
                'JackKnifeData_bipart.dat',
                unpack=True)
   

    errCheck = False
    if errCheck:
        EsNorm, EsErrNorm = pl.array([]), pl.array([])
        for fileName in args.fileNames:
            #Ecv,Eth = pl.loadtxt(fileName, unpack=True, usecols=(4,-5))
            Ecv = pl.loadtxt(fileName, unpack=True, usecols=(4,))
            EsNorm = pl.append(EsNorm,pl.average(Ecv))
            #ET = pl.append(ET, pl.average(Eth))
            EsErrNorm = pl.append(EsErrNorm, pl.std(Ecv)/pl.sqrt(float(len(Ecv))))
            #ETerr = pl.append(ETerr, pl.std(Eth)/pl.sqrt(float(len(Eth))))

        pl.scatter(temps, EsErrNorm, label='Standard Error', color='Navy')
        pl.scatter(temps, EsErr, label='Jackknife Error', color='Orange')
        pl.grid()
        pl.legend()
        pl.show()

    QHO = False
    if QHO:
        # analytical solutions for 1D QHO with one particle
        tempRange = pl.arange(0.01,1.0,0.01)
        Eanalytic = 0.5/pl.tanh(1.0/(2.0*tempRange))
        CvAnalytic = 1.0/(4.0*(tempRange*pl.sinh(1.0/(2.0*tempRange)))**2)

    # some Plotting options
    ShareAxis=True      # shared x-axis for Cv and Energy
    if reduceType == 'T':
        xLab = 'Temperature [K]'
    elif reduceType == 'u':
        xLab = 'Chemical Potential [K]'

    # plot the specific heat vs. temperature
    if ShareAxis:
        ax1 = pl.subplot(211)
    else:
        pl.figure(1)
    if QHO: # plot analytic result
        pl.plot(tempRange,CvAnalytic, label='Exact')
    pl.errorbar(temps,Cvs,CvsErr, label='PIMC',color='Violet',fmt='o')
    if not ShareAxis:
        pl.xlabel(xLab, fontsize=20)
    pl.ylabel('Specific Heat', fontsize=20)
    pl.grid(True)
    pl.legend(loc=2)
    
    # plot the energy vs. temperature
    if ShareAxis:
        pl.setp(ax1.get_xticklabels(), visible=False)
        ax2 = pl.subplot(212, sharex=ax1)
    else:
        pl.figure(2)
    if QHO: # plot analytic result
        pl.plot(tempRange,Eanalytic, label='Exact')
    pl.errorbar(temps,Es,EsErr, label='PIMC virial',color='Lime',fmt='o')
    pl.xlabel(xLab, fontsize=20)
    pl.ylabel('Energy [K]', fontsize=20)
    pl.grid(True)
    pl.legend(loc=2)

    pl.savefig('Helium_critical_CVest.pdf', format='pdf',
            bbox_inches='tight')

    if ShareAxis:
        pl.figure(2)
    else:
        pl.figure(3)
    pl.errorbar(temps, rhos_rhos, rhos_rhoErr, fmt='o')
    pl.xlabel(xLab, fontsize=20)
    pl.ylabel('Superfluid Stiffness', fontsize=20)
    pl.grid(True)

    if ShareAxis:
        pl.figure(3)
    else:
        pl.figure(4)
    pl.errorbar(temps, filmDenses, filmDensErrs, label='film', fmt='o')
    pl.errorbar(temps, bulkDenses, bulkDensErrs, label='bulk', fmt='o')
    pl.xlabel(xLab, fontsize=20)
    pl.ylabel(r'$\mathrm{Density}\ [\AA^{-d}]$', fontsize=20) 
    pl.legend()
    pl.grid(True)
    
    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
