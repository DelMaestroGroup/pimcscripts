# =============================================================================
# Pulls data for new winding number from zAveragedNtwind.dat files.
# This requires that data be arranged as:
#       windingFiles/xx.xangFilm_xx.xangExtent/zAveragedNtwind.dat
# where the xx.x are film thickness and extent (in that order) and each
# NtWind averaged file has its own directory.
#
# Author:           Max Graves
# Last Revised:     28-JAN-2014
# =============================================================================

import pylab as pl
import MTG_jkTools as jk
import glob,sys,os,random
from matplotlib import rcParams

# set up latex fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}', # load siunitx
       r'\sisetup{detect-all}'  # force siunitx to use your fonts
]

def main():

    # --- Set up all options --------------------------------------------------
    
    # determine normalization factor for NormAn winding.
    S = 5.0
    Ly = 12.0
    normFactor = 4.0*(S+Ly)**2

  
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
 
    # get names of all subdirectories LzXXX.
    os.chdir(direc)
    direcs = glob.glob('*angFilm_*')
    if (direcs == []):
        direcs = glob.glob('*Lz*')
   
    scaleVar = jk.scalingVariable(direcs)

    # --- this section is used for a single temperature -----------------------

    singleTemperature = True

    # build array of Lz values for plotting Lz vs winding.
    LzValues = pl.array([])

    # build array of norman winding values along with bins, for Lz plotting.
    allAverages = pl.array([])
    allErrors = pl.array([])

    # --- loop over and collect/analyze all data from files -------------------
    for nd, d in enumerate(sorted(direcs)):

        os.chdir('./'+d)
        f = glob.glob('*')[0]

        if (scaleVar != 'Lz'):
            thickness = d[:4]
            extent = d[12:16]
        else:
            Lz = d[2:]
            extent = args.bulkSeparation
            LzValues = pl.append(LzValues, float(Lz))
        
        # determine titles for plotting
        if (scaleVar == 'extent'):
            labell = 'thickness: '+thickness+' '+r'$[\si{\angstrom}]$'
            titlle = 'Bulk Separation: '+extent+' '+r'$[\si{\angstrom}]$'
        elif (scaleVar == 'thickness'):
            labell = 'Bulk Separation: '+extent+' '+r'$[\si{\angstrom}]$'
            titlle = 'thickness: '+str(thickness)+' '+r'$[\si{\angstrom}]$'
        elif (scaleVar == 'Lz'):
            labell = r'$L_z$: '+Lz+' '+r'$[\si{\angstrom}]$'
            titlle = 'Bulk Separation: '+str(extent)+' '+r'$[\si{\angstrom}]$'

        # get headers from file.
        headers = jk.getHeadersFromFile(f)

        # determine if this is for a single temperature
        if len(headers) > 1:
            singleTemperature = False
        
        # build array of reduced temperatures t = |1-T/T_c|
        ReducedTemps = pl.array([])
        for T in headers:
            ReducedTemps = pl.append(ReducedTemps, abs(1.0 - float(T)/2.17))
        
        AVG = pl.array([])
        STD = pl.array([])

        n = nCol-1
        for header in headers:

            avgs,stds,bins = pl.genfromtxt(f, 
                    usecols=(0+3*n, 1+3*n, 2+3*n), 
                    unpack=True, delimiter=',')

            # get rid of any items which are not numbers..
            # this is some beautiful Python juju.
            bins = bins[pl.logical_not(pl.isnan(bins))]
            stds = stds[pl.logical_not(pl.isnan(stds))]
            avgs = avgs[pl.logical_not(pl.isnan(avgs))]

            # normalize data.
            stds *= normFactor
            avgs *= normFactor
            
            weights = bins/pl.sum(bins)

            avgs *= weights
            stds *= weights

            avg = pl.sum(avgs)
            stdErr = pl.sum(stds)

            AVG = pl.append(AVG, avg)
            STD = pl.append(STD, stdErr)

            n += nEst

        # if we find a range of temperatures, plot all of that data.
        if not singleTemperature:

            pl.figure(1)
            pl.errorbar(headers, AVG, STD, fmt='o', color=colors[nd], 
                    label=labell)
            pl.xlabel(xLab, fontsize=20)
            #pl.ylabel(r'$\Omega$', fontsize=20)
            #pl.ylabel(r'$\langle \Omega^2 \rangle/2 \beta \lambda N_{\text{film}} $', fontsize=20)
            pl.ylabel(r'$\langle \Omega \rangle$', fontsize=20)
            #pl.ylabel(r'$\langle N \rangle$', fontsize=20)
            pl.title(titlle)
            pl.legend()
            pl.grid(True)
            
            pl.figure(2)
            pl.yscale('log')
            pl.xscale('log')
            pl.errorbar(ReducedTemps, AVG, STD, fmt='o', color=colors[nd], 
                    label=labell)
            pl.xlabel(xLab, fontsize=20)
            pl.ylabel(r'$\langle \Omega^2 \rangle$', fontsize=20)
            #pl.ylabel(r'$\langle N \rangle$', fontsize=20)
            pl.title(titlle)
            pl.legend()
            pl.grid(True)
        else:
            allAverages = pl.append(allAverages, avg)
            allErrors = pl.append(allErrors, stdErr)

        os.chdir('..')

    if singleTemperature:
        invLz = 1.0/LzValues
        figg = pl.figure(1)
        ax = figg.add_subplot(111)
        pl.errorbar(invLz, allAverages, allErrors, fmt='o')
        pl.xlabel(r'$1/L_z\ [\si{\angstrom}^{-1}]$', fontsize=20)
        pl.ylabel(r'$\langle \Omega \rangle$', fontsize=26)
        pl.xlim([0,0.06])
        #pl.title(r'$S =\ $'+str(extent)+r' $[\si{\angstrom}]$'+' , '+r'$T=\ $'+str(float(headers[0]))+' [K]')
        pl.grid(True)

        pl.tick_params(axis='both', which='major', labelsize=16)
        pl.tick_params(axis='both', which='minor', labelsize=16)
        xticks = ax.xaxis.get_major_ticks()
        xticks[0].set_visible(False)

        pl.savefig('Omega_vs_inverseLZ.pdf', format='pdf',
                bbox_inches='tight')

    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
