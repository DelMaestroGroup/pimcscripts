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

    dens = False

    # --- Set up all options --------------------------------------------------
    
    # parse command line, getting algorithmic and plotting options.
    args = jk.parseCMD()
    reduceType = args.reduceType
    direc = args.fileNames[0]

    # some plotting color and label options
    colors = ['Salmon','MediumSpringGreen','DarkViolet','Fuchsia','Blue',
            'Maroon']

    if args.RandomColors:
        random.shuffle(colors)
 
    # set up figure that displays all data
    figg = pl.figure(1)
    ax = figg.add_subplot(111)
    pl.xlabel(r'$T\ [K]$', fontsize=20)
    #pl.ylabel('Spatial Density '+r'$[\si{\angstrom}^{-d}]$', fontsize=20)
    pl.ylabel('Energy '+r'$[K]$', fontsize=20)
    pl.grid(True)
    if dens:
        pl.xlim([0.5,3.0])
        pl.ylim([0,0.07])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    if dens:
        # set up bulk SVP densities for plot
        bulkVert = -30
        minMus = 0.5
        maxMus = 2.5
        boxSubtract = 1.6
        pl.plot([minMus, maxMus], [0.02198, 0.02198], 'k-', lw=3)
        pl.annotate('3D SVP', xy=(maxMus - boxSubtract, 0.02195),  #xycoords='data',
                xytext=(-50, bulkVert), textcoords='offset points',
                bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="->",
                    connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                )

        pl.plot([minMus, maxMus], [0.0432, 0.0432], 'k-', lw=3)
        pl.annotate('2D SVP', xy=(maxMus - boxSubtract, 0.0432),  #xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="->",
                    connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                )

    # --- loop over all values of S -------------------------------------------
    os.chdir(direc)
    Svals = glob.glob('S*')

    for nS, Sval in enumerate(Svals):
        
        os.chdir(Sval)

        # store bulk separation value
        S = re.search(r'\d+',Sval).group(0)

        # set label for plot
        if 'distinguishable' in Sval:
            labell = 'S = '+str(S)+', Boltzmannons'
        else:
            labell = 'S = '+str(S)+', Bosons'

        # get all temperature directory names
        Tdirs = sorted(glob.glob('T*'))

        Ts = pl.array([])
        Films = pl.array([])
        Bulks = pl.array([])
        filmErrs = pl.array([])
        bulkErrs = pl.array([])

        # --- loop over all temperature values --------------------------------
        for Tdir in Tdirs:
            
            os.chdir(Tdir)
            
            # build array of temperatures
            Ts = pl.append(Ts, float(Tdir[1:]))

            # get data file name
            #f = glob.glob('*Bipart*')[0]
            f = glob.glob('*Estimator*')[0]

            #densData = pl.genfromtxt(f, delimiter=',')
            filmavg,filmstd,filmbins,bulkavg,bulkstd,bulkbins = pl.genfromtxt(f, 
                    unpack=True, usecols=(0,1,2,3,4,5), delimiter=',')

            # get rid of any items which are not numbers..
            # this is some beautiful Python juju.
            filmbins = filmbins[pl.logical_not(pl.isnan(filmbins))]
            filmstd = filmstd[pl.logical_not(pl.isnan(filmstd))]
            filmavg = filmavg[pl.logical_not(pl.isnan(filmavg))]
            bulkbins = bulkbins[pl.logical_not(pl.isnan(bulkbins))]
            bulkstd = bulkstd[pl.logical_not(pl.isnan(bulkstd))]
            bulkavg = bulkavg[pl.logical_not(pl.isnan(bulkavg))]


            filmweights = filmbins/pl.sum(filmbins)
            bulkweights = bulkbins/pl.sum(bulkbins)

            filmavg *= filmweights
            bulkavg *= bulkweights

            filmstd *= filmweights
            bulkstd *= bulkweights

            film = pl.sum(filmavg)
            bulk = pl.sum(bulkavg)
            filmstdErr = pl.sum(filmstd)
            bulkstdErr = pl.sum(bulkstd)

            Films = pl.append(Films, film)
            Bulks = pl.append(Bulks, bulk)
            filmErrs = pl.append(filmErrs, filmstdErr)
            bulkErrs = pl.append(bulkErrs, bulkstdErr)

            os.chdir('..')
        
        # add data to plot for given S value.
        pl.errorbar(Ts, Films, filmErrs, fmt='-o', color=colors[nS], 
                label=labell+', 2d')
        pl.errorbar(Ts, Bulks, bulkErrs, fmt = '--d', color=colors[nS],
                label=labell+', 3d')

        os.chdir('..')

    pl.legend()
    pl.savefig('Densities_vs_T_allS.pdf', format='pdf',
            bbox_inches='tight')
    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
