import pylab as pl
import glob,argparse,sys
import MTG_jkTools as aTools
from matplotlib import rcParams
import csv

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

def main():

    estimList = []

    QHO = False
    args = aTools.parseCMD()
    fileNames = args.fileNames
    skip = args.skip
 
    # Check if any jackknife data files exist already.
    check = glob.glob('*Jack*')
    if check != []:
        sys.exit('Check your directory for Jackknife data files')

    reduceType = args.reduceType

    # check which ensemble
    canonical=True
    if fileNames[0][0]=='g':
        canonical=False

    # create list of temperatures
    temps = pl.array([])
    tempList = aTools.getHeadersFromFile(fileNames[0])
    for temp in tempList:
        temps = pl.append(temps,float(temp))

    # perform jackknife analysis of data, writing to disk
    for fileName in fileNames:
        
        print '\n\n---',fileName,'---\n'
        
        # ========================================
        # NON TRIVIAL WINDING ESTIMATOR
        # ========================================
        if 'Ntwind' in fileName:
            estimList.append('Ntwind')
            ntWinds, ntWindErrs     = pl.array([]), pl.array([]) 
            
            fout  = open('JackKnifeData_ntWind.dat','w')
            fout.write('#%15s\t%16s\n'%(reduceType, '<W^2>'))

            # read in data and call jackknife method.
            for n,temp in enumerate(tempList):
                ntWind = pl.array([])
                with open(fileName, 'rb') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        if row[0][0] == '#':
                            continue
                        else:
                            try:
                                ntWind = pl.append(ntWind, float(row[n]))
                            except:
                                continue

                ntWindAve, ntWindErr = aTools.jackknife(ntWind[skip:])
                print reduceType,' = ', temp
                print '<W^2> = ', ntWindAve,' +/- ',ntWindErr
                print 'numBins: ',len(ntWind),'\n'
                ntWinds      = pl.append(ntWinds, ntWindAve)
                ntWindErrs    = pl.append(ntWindErrs, ntWindErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), ntWindAve, ntWindErr ))
            fout.close()
       
        # ========================================
        # ENERGY AND SPECIFIC HEAT ESTIMATOR
        # ========================================
        if 'Estimator' in fileName:
            estimList.append('Estimator')
            Cvs,CvsErr = pl.array([]),pl.array([])
            Es, EsErr   = pl.array([]), pl.array([]) 
            
            fout = open('JackKnifeData_Cv.dat', 'w')
            fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'% (
                reduceType, 'E', 'Eerr', 'Cv', 'CvErr'))
            
            # read in data and call jackknife method.
            n = 0
            for temp in tempList:
                E = pl.array([])
                Cv1 = pl.array([])
                Cv2 = pl.array([])
                Cv3 = pl.array([])
                with open(fileName, 'rb') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        if row[0][0] == '#':
                            continue
                        else:
                            try:
                                E = pl.append(E, float(row[n]))
                                Cv1 = pl.append(Cv1, float(row[n+1]))
                                Cv2 = pl.append(Cv2, float(row[n+2]))
                                Cv3 = pl.append(Cv3, float(row[n+3]))
                            except:
                                continue
                n += 4

                EAve, Eerr = aTools.jackknife(E[skip:])
                jkAve, jkErr = aTools.jackknife(
                        Cv1[skip:],Cv2[skip:],Cv3[skip:])
                print reduceType,' = ',float(temp)
                print '<E>  = ',EAve,' +/- ',Eerr
                print '<Cv> = ',jkAve,' +/- ',jkErr
                Es      = pl.append(Es, EAve)
                Cvs     = pl.append(Cvs, jkAve)
                EsErr   = pl.append(EsErr, Eerr)
                CvsErr  = pl.append(CvsErr, jkErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), EAve, Eerr, jkAve, jkErr)) 

            fout.close()

        # ========================================
        # SUPERFLUID ESTIMATORS
        # ========================================
        if 'Super' in fileName:
            estimList.append('Super')
            rhos_rhos, rhos_rhoErr  = pl.array([]), pl.array([])
            Wx2s, Wx2Err = pl.array([]), pl.array([])
            Wy2s, Wy2Err = pl.array([]), pl.array([])
            Wz2s, Wz2Err = pl.array([]), pl.array([])

            fout = open('JackKnifeData_super.dat','w')
            fout.write('#%15s,\t%16s,\t%16s,\t%16s,\t%16s,\t\
                    %16s,\t%16s,\t%16s,\t%16s\n'%(
                reduceType, 'rho_s/rho', 'rho_s/rhoErr', 'Wx^2', 'Wx2_err', 
                'Wy^2', 'Wy2_err', 'Wz^2', 'Wz2_err'))
      
            # read in data and call jackknife method.
            n = 0
            for temp in tempList:
                rhos_rho = pl.array([])
                wx2 = pl.array([])
                wy2 = pl.array([])
                wz2 = pl.array([])
                with open(fileName, 'rb') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        if row[0][0] == '#':
                            continue
                        else:
                            try:
                                rhos_rho = pl.append(rhos_rho, 
                                        float(row[n]))
                                wx2 = pl.append(wx2, float(row[n+1]))
                                wy2 = pl.append(wy2, float(row[n+2]))
                                wz2 = pl.append(wz2, float(row[n+3]))
                            except:
                                continue
                n += 4

                superAve, superErr = aTools.jackknife(rhos_rho[skip:])
                wx2Ave, wx2Err = aTools.jackknife(wx2[skip:])
                wy2Ave, wy2Err = aTools.jackknife(wy2[skip:])
                wz2Ave, wz2Err = aTools.jackknife(wz2[skip:])
                print 'rho_s/rho = ', superAve,' +/- ',superErr
                rhos_rhos   = pl.append(rhos_rhos, superAve)
                rhos_rhoErr = pl.append(rhos_rhoErr, superErr)
                Wx2s        = pl.append(Wx2s, wx2Ave)
                Wx2Err      = pl.append(Wx2Err, wx2Err)
                Wy2s        = pl.append(Wy2s, wy2Ave)
                Wy2Err      = pl.append(Wy2Err, wy2Err)
                Wz2s        = pl.append(Wz2s, wz2Ave)
                Wz2Err      = pl.append(Wz2Err, wz2Err)
                fout.write('%16.8E,\t%16.8E,\t%16.8E,\t%16.8E,\t\
                        %16.8E,\t%16.8E,\t%16.8E,\t%16.8E,\t%16.8E,\n' %(
                    float(temp), superAve, superErr, wx2Ave, wx2Err,
                    wy2Ave, wy2Err, wz2Ave, wz2Err))
            
        # ========================================
        # FILM AND BULK DENSITY ESTIMATORS
        # ========================================
        if 'Bipart' in fileName:
            estimList.append('Bipart')
            filmDenses, bulkDenses      = pl.array([]), pl.array([])
            filmDensErrs, bulkDensErrs  = pl.array([]), pl.array([])
            
            fout = open('JackKnifeData_bipart.dat','w')
            fout.write('#%15s\t%16s\t%16s\t%16s\t%16s\n'%(reduceType, 
                'filmDens', 'filmDensErr', 'bulkDens', 'bulkDensErr'))
       
            # read in data and call jackknife method.
            n = 0
            for temp in tempList:
                filmDens = pl.array([])
                bulkDens = pl.array([])
                with open(fileName, 'rb') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        if row[0][0] == '#':
                            continue
                        else:
                            try:
                                filmDens = pl.append(filmDens, 
                                        float(row[n+2]))
                                bulkDens = pl.append(bulkDens, 
                                        float(row[n+3]))
                            except:
                                continue
                n += 2

                filmDensAve, filmDensErr = aTools.jackknife(filmDens[skip:])
                bulkDensAve, bulkDensErr = aTools.jackknife(bulkDens[skip:])
                print 'filmDens = ', filmDensAve,' +/- ',filmDensErr
                print 'bulkDens = ', bulkDensAve,' +/- ',bulkDensErr
                filmDenses      = pl.append(filmDenses, filmDensAve)
                filmDensErrs    = pl.append(filmDensErrs, filmDensErr)
                bulkDenses  = pl.append(bulkDenses, bulkDensAve)
                bulkDensErrs    = pl.append(bulkDensErrs, bulkDensErr)
                fout.write('%16.8E\t%16.8E\t%16.8E\t%16.8E\t%16.8E\n' %(
                    float(temp), filmDensAve, filmDensErr, 
                    bulkDensAve, bulkDensErr))
            
            fout.close()

    # =========================================================================
    # PLOTTING SECTION
    # =========================================================================
    
    xLab = jk.getXlabel(reduceType)
    
    plotNum = 1

    # PLOT ntWinding number squared
    if 'Ntwind' in estimList:
        pl.figure(plotNum)
        pl.errorbar(temps, ntWinds, ntWindErrs, fmt='o')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel(r'$\langle \Omega^2 \rangle$', fontsize=20) 
        pl.grid(True)
        plotNum += 1

    # PLOT specific heat vs. temperature
    if 'Estimator' in estimList:
        pl.figure(2)
        ax1 = pl.subplot(211)
        #if QHO: # plot analytic result
        #    pl.plot(tempRange,CvAnalytic, label='Exact')
        pl.errorbar(temps,Cvs,CvsErr, label='PIMC',color='Violet',fmt='o')
        pl.ylabel('Specific Heat', fontsize=20)
        pl.grid(True)
        pl.legend(loc=2)
        
        # PLOT the energy vs. temperature
        pl.setp(ax1.get_xticklabels(), visible=False)
        ax2 = pl.subplot(212, sharex=ax1)
        #if QHO: # plot analytic result
        #    pl.plot(tempRange,Eanalytic, label='Exact')
        pl.errorbar(temps,Es,EsErr, label='PIMC virial',color='Lime',fmt='o')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel('Energy [K]', fontsize=20)
        pl.grid(True)
        pl.legend(loc=2)

        plotNum += 1

    # PLOT superfluid properties
    if 'Super' in estimList:
        pl.figure(plotNum)
        pl.errorbar(temps, rhos_rhos, rhos_rhoErr, fmt='o')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel('Superfluid Stiffness', fontsize=20)
        pl.grid(True)

        plotNum += 1

        pl.figure(plotNum)
        pl.errorbar(temps, Wx2s, Wx2Err, fmt='o', label=r'$\langle W_x^2 \rangle$')
        pl.errorbar(temps, Wy2s, Wy2Err, fmt='o', label=r'$\langle W_y^2 \rangle$')
        pl.errorbar(temps, Wz2s, Wz2Err, fmt='o', label=r'$\langle W_z^2 \rangle$')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel(r'$\langle W_i^2 \rangle$', fontsize=20)
        pl.legend()
        pl.grid(True)

        plotNum += 1

        pl.figure(plotNum)
        pl.errorbar(temps, Wz2s, Wz2Err, fmt='o', label=r'$\langle W_z^2 \rangle$')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel(r'$\langle W_z^2 \rangle$', fontsize=20)
        pl.legend()
        pl.grid(True)

        plotNum += 1

    # PLOT film and bulk densities
    if 'Bipart' in estimList:
        pl.figure(plotNum)
        pl.errorbar(temps, filmDenses, filmDensErrs, label='film', fmt='o')
        pl.errorbar(temps, bulkDenses, bulkDensErrs, label='bulk', fmt='o')
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel(r'$\mathrm{Density}\ [\AA^{-d}]$', fontsize=20) 
        pl.legend()
        pl.grid(True)

    pl.show()


    # BELOW IS RESIDUE FROM A PREVIOUS VERSION OF THIS SCRIPT.
    # SOME OF IT MIGHT PROVE USEFUL.
    '''else:
        print 'Found existing data file in CWD.'
        temps, Es, EsErr, Cvs, CvsErr = pl.loadtxt(
                'JackKnifeData_Cv.dat', 
                unpack=True)
        if not QHO:
            temps, rhos_rhos, rhos_rhoErr, Wx2s, Wx2Err, Wy2s, Wy2Err, Wz2s, Wz2Err = pl.loadtxt(
                    'JackKnifeData_super.dat', 
                    unpack=True)
        if excVol:
            temps, filmDenses, filmDensErrs, bulkDenses, bulkDensErrs = pl.loadtxt(
                    'JackKnifeData_bipart.dat',
                    unpack=True)'''
   

    '''
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

    if QHO:
        # analytical solutions for 1D QHO with one particle
        tempRange = pl.arange(0.01,1.0,0.01)
        Eanalytic = 0.5/pl.tanh(1.0/(2.0*tempRange))
        CvAnalytic = 1.0/(4.0*(tempRange*pl.sinh(1.0/(2.0*tempRange)))**2)
    
       pl.savefig('Helium_critical_CVest_trans.pdf', format='pdf',
            bbox_inches='tight', transparent=True)

    '''


# =============================================================================
if __name__=='__main__':
    main()
