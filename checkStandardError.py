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

    S = 5.0
    Ly = 12.0
    normFactor = 4.0*(S+Ly)**2
   
    args = jk.parseCMD()
    reduceType = args.reduceType
    direc = args.fileNames[0]

    #dependentVar = direc[:-6]
    nCol = args.nCol
    nEst = args.nEst

    os.chdir(direc)
    direcs = glob.glob('*angFilm_*')
    if (direcs == []):
        direcs = glob.glob('*Lz*')

    # some plotting options
    xLab = jk.getXlabel(reduceType)
    colors = ['Salmon','MediumSpringGreen','DarkViolet','Fuchsia','Blue',
            'Maroon']

    # shuffle colors around randomly
    if args.RandomColors:
        random.shuffle(colors)
    
    scaleVar = jk.scalingVariable(direcs)

    for nd, d in enumerate(sorted(direcs)):

        os.chdir('./'+d)
        f = glob.glob('*')[0]

        if (scaleVar != 'Lz'):
            thickness = d[:4]
            extent = d[12:16]
        else:
            Lz = d[2:]
            extent = args.bulkSeparation
        
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

        headers = jk.getHeadersFromFile(f)
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
 
                # add current data to plot
        pl.figure(1)
        pl.errorbar(headers, AVG, STD, fmt='o', color=colors[nd], 
                label=labell)
        pl.xlabel(xLab, fontsize=20)
        pl.ylabel(r'$\Omega$', fontsize=20)
        #pl.ylabel(r'$\langle \Omega^2 \rangle/2 \beta \lambda N_{\text{film}} $', fontsize=20)
        #pl.ylabel(r'$\langle \Omega^2 \rangle$', fontsize=20)
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

        os.chdir('..')
    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
