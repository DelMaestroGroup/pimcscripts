# =============================================================================
# Pulls data for new winding number from zAveragedNtwind.dat files.
# This requires that data be arranged as:
#       windingFiles/xx.xangFilm_xx.xangExtent/zAveragedNtwind.dat
# where the xx.x are film thickness and extent (in that order) and each
# NtWind averaged file has its own directory.
#
# Author:           Max Graves
# Last Revised:     18-APR-2014
# =============================================================================

import pylab as pl
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
    S = 20.0
    Ly = 12.0
    normFactor = 4.0*(S+Ly)**2

    # some plotting color and label options
    colors = ['Salmon','MediumSpringGreen','DarkViolet','Fuchsia','Blue',
            'Maroon']

    RandomColors = False
    if RandomColors:
        random.shuffle(colors)
 
    # get names of all subdirectories az_XXX.
    direcs = glob.glob('az_*')

    # build array of Lz values for plotting Lz vs winding.
    azValues = pl.array([])

    # build array of norman winding values along with bins, for Lz plotting.
    allAverages = pl.array([])
    allErrors = pl.array([])

    # --- loop over and collect/analyze all data from files -------------------
    for nd, d in enumerate(sorted(direcs)):

        os.chdir('./'+d)
        f = glob.glob('*Ntwind*')[0]

        print d[3:]
        
        azValues = pl.append(azValues, float(d[3:]))
                
        avgs,stds,bins = pl.genfromtxt(f, usecols=(3, 4, 5), 
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

        allAverages = pl.append(allAverages, avg)
        allErrors = pl.append(allErrors, stdErr)

        os.chdir('..')


    invaz = 1.0/azValues
    figg = pl.figure(1)
    ax = figg.add_subplot(111)
    pl.errorbar(invaz, allAverages, allErrors, fmt='o')
    pl.xlabel(r'$1/L_z\ [\si{\angstrom}^{-1}]$', fontsize=20)
    pl.ylabel(r'$\langle \Omega \rangle$', fontsize=26)
    pl.xlim([0,0.06])
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].set_visible(False)

    #pl.savefig('Omega_vs_inverseLZ.pdf', format='pdf',
    #        bbox_inches='tight')

    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
