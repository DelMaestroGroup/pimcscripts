# =============================================================================
# Pulls data for new winding number from zAveragedNtwind.dat files.
# This requires that data be arranged as:
#       windingFiles/S{}/T{}zAveragedNtwind.dat
# where the {} are film extent and temperature (in that order) and each
# NtWind averaged file has its own directory.
#
# Author:           Max Graves
# Last Revised:     28-JAN-2014
# =============================================================================

import pylab as pl
import jkTools as jk
import glob,sys,os,random,re
from matplotlib import rcParams
import clusterTools as cT
import natsort

# set up latex fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}', # load siunitx
       r'\sisetup{detect-all}'  # force siunitx to use your fonts
]


def main():

    omega = True
    energy = False
    superFrac = True

    # --- Set up all options --------------------------------------------------
    
    Ly = 12.0

    # parse command line, getting algorithmic and plotting options.
    args = jk.parseCMD()
    reduceType = args.reduceType
    direc = args.fileNames[0]

    nCol = args.nCol
    nEst = args.nEst

    # some plotting color and label options
    xLab = jk.getXlabel(reduceType)
    extent = args.bulkSeparation
    colors = ['Salmon','MediumSpringGreen','DarkViolet','Fuchsia','Blue',
            'Maroon']

    if args.RandomColors:
        random.shuffle(colors)
 
    # set up figure that displays all data
    figg = pl.figure(1)
    ax = figg.add_subplot(111)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    pl.ylabel(r'$\langle \Omega \rangle$', fontsize=20)
    pl.grid(True)
    pl.xlim([0.4,2.6])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    figg2 = pl.figure(2)
    ax = figg2.add_subplot(111)
    pl.ylabel(r'$\langle \rho_s/\rho \rangle$', fontsize=20)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    pl.grid(True)
    pl.xlim([0.4,2.6])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

 
    figg3 = pl.figure(3)
    ax = figg3.add_subplot(111)
    pl.ylabel(r'$\langle \rho_{\text{film}} \rangle$', fontsize=20)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    pl.grid(True)
    pl.xlim([0.4,2.6])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

 
    figg4 = pl.figure(4)
    ax = figg4.add_subplot(111)
    pl.ylabel(r'$\langle \rho_{\text{bulk}} \rangle$', fontsize=20)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    pl.grid(True)
    pl.xlim([0.4,2.6])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

  
    
    # --- loop over all values of S -------------------------------------------
    os.chdir(direc)
    Svals = glob.glob('S*')
    Svals = natsort.natsorted(Svals)
    print Svals

    for nS, Sval in enumerate(Svals):
        
        os.chdir(Sval)

        # store bulk separation value
        S = re.search(r'\d+',Sval).group(0)

        print S

        # set normalization
        if omega:
            normFactor = 4.0*(float(S)+Ly)**2
        else:
            normFactor = 1.0

        print normFactor

        # set label for plot
        if 'distinguishable' in Sval:
            labell = 'S = '+str(S)+', Boltzmannons'
        else:
            labell = 'S = '+str(S)+' '+r'$\si{\angstrom}$'+', Bosons'

        # get all temperature directory names
        Tdirs = sorted(glob.glob('T*'))

        # build array of norman winding values along with bins, for Lz plotting.
        windingAverages = pl.array([])
        windingErrors = pl.array([])
        filmAverages = pl.array([])
        filmErrors = pl.array([])
        bulkAverages = pl.array([])
        bulkErrors = pl.array([])
        superAverages = pl.array([])
        superErrors = pl.array([])

        Ts = pl.array([])
        #Omegas = pl.array([])
        #Errs = pl.array([])
        #Films = pl.array([])
        #Ferrs = pl.array([])
        #Bulks = pl.array([])
        #Berrs = pl.array([])
        #Supers = pl.array([])
        #Serrs = pl.array([])

        # --- loop over all temperature values --------------------------------
        for Tdir in Tdirs:
            
            os.chdir(Tdir)
            
            # build array of temperatures
            Ts = pl.append(Ts, float(Tdir[1:]))
            
            # --- angular winding ---
            f = glob.glob('*zAveragedNtwind*')[0]
            
            aCol = 3
            sCol = 4
            bCol = 5
            
            avg,stdErr = cT.crunchZfile(f,aCol,sCol,bCol,normFactor)
            windingAverages = pl.append(windingAverages, avg)
            windingErrors = pl.append(windingErrors, stdErr)
            
            # --- film densities ---
            f = glob.glob('*zAveragedBipart*')[0]
            aCol = 0
            sCol = 1
            bCol = 2

            avg,stdErr = cT.crunchZfile(f,aCol,sCol,bCol,1.0)

            filmAverages = pl.append(filmAverages, avg)
            filmErrors = pl.append(filmErrors,stdErr)

            # --- bulk densities ---
            aCol = 3
            sCol = 4
            bCol = 5

            avg,stdErr = cT.crunchZfile(f,aCol,sCol,bCol,1.0)

            bulkAverages = pl.append(bulkAverages,avg)
            bulkErrors = pl.append(bulkErrors, stdErr)

            # --- superfluid fractions ---
            f = glob.glob('*zAveragedSuper*')[0]
            aCol = 0
            sCol = 1
            bCol = 2

            avg,stdErr = cT.crunchZfile(f,aCol,sCol,bCol,1.0)

            superAverages = pl.append(superAverages,avg)
            superErrors = pl.append(superErrors, stdErr)

            # ----------------------
            os.chdir('..')
       

        # add data to plot for given S value.
        pl.figure(1)
        pl.errorbar(Ts, windingAverages, windingErrors, fmt='-o', color=colors[nS], 
                label=labell)
        
        pl.figure(2)
        pl.errorbar(Ts, superAverages, superErrors, fmt='-o', color=colors[nS], 
                label=labell)
        
        pl.figure(3)
        pl.errorbar(Ts, filmAverages, filmErrors, fmt='-o', color=colors[nS], 
                label=labell)
        
        pl.figure(4)
        pl.errorbar(Ts, bulkAverages, bulkErrors, fmt='-o', color=colors[nS], 
                label=labell)

        os.chdir('..')
       
    pl.figure(1)
    pl.legend()
    
    pl.figure(2)
    pl.legend()
    
    pl.figure(3)
    pl.legend()
    
    pl.figure(4)
    pl.legend()
    #pl.savefig('Omega_vs_T_allS.pdf', format='pdf',
    #        bbox_inches='tight')
    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
