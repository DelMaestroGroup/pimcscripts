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
import glob,sys,os,random,re
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
    pl.xlim([0.5,3.0])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    
    # --- loop over all values of S -------------------------------------------
    os.chdir(direc)
    Svals = glob.glob('S*')

    for nS, Sval in enumerate(Svals):
        
        os.chdir(Sval)

        # store bulk separation value
        S = re.search(r'\d+',Sval).group(0)

        # set normalization
        normFactor = 4.0*(float(S)+Ly)**2

        # set label for plot
        if 'distinguishable' in Sval:
            labell = 'S = '+str(S)+', Boltzmannons'
        else:
            labell = 'S = '+str(S)+', Bosons'

        # get all temperature directory names
        Tdirs = sorted(glob.glob('T*'))

        Ts = pl.array([])
        Omegas = pl.array([])
        Errs = pl.array([])

        # --- loop over all temperature values --------------------------------
        for Tdir in Tdirs:
            
            os.chdir(Tdir)
            
            # build array of temperatures
            Ts = pl.append(Ts, float(Tdir[1:]))

            # get data file name
            f = glob.glob('*Ntwind*')[0]

            n = nCol-1

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

            Omegas = pl.append(Omegas, avg)
            Errs = pl.append(Errs, stdErr)

            n += nEst

            os.chdir('..')
        
        # add data to plot for given S value.
        pl.errorbar(Ts, Omegas, Errs, fmt='-o', color=colors[nS], 
                label=labell)

        os.chdir('..')

    pl.legend()
    pl.savefig('Omega_vs_T_allS.pdf', format='pdf',
            bbox_inches='tight')
    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
