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
   
    args = jk.parseCMD()
    reduceType = args.reduceType
    direc = args.fileNames[0]
    nEst = args.nEst

    os.chdir(direc)
    direcs = glob.glob('*angFilm_*')

    # some plotting options
    xLab = jk.getXlabel(reduceType)
    colors = ['Salmon','Blue','DarkViolet','MediumSpringGreen','Fuchsia',
            'Yellow','Maroon']

    # shuffle colors around randomly
    if args.RandomColors:
        random.shuffle(colors)
    
    extentScale = jk.thicknessORextent(direcs)

    for nd, d in enumerate(sorted(direcs)):

        os.chdir('./'+d)
        f = glob.glob('*')[0]

        thickness = d[:4]
        extent = d[12:16]

        headers = jk.getHeadersFromFile(f)
        AVG = pl.array([])
        STD = pl.array([])

        #for n in range(int(len(headers))):
        n = nEst-1
        for header in headers:

            avgs,stds,bins = pl.genfromtxt(f, 
                    usecols=(0+3*n, 1+3*n, 2+3*n), 
                    unpack=True, delimiter=',')

            # get rid of any items which are not numbers..
            # this is some beautiful Python juju.
            bins = bins[pl.logical_not(pl.isnan(bins))]
            stds = stds[pl.logical_not(pl.isnan(stds))]
            avgs = avgs[pl.logical_not(pl.isnan(avgs))]
            
            weights = bins/pl.sum(bins)

            avgs *= weights
            stds *= weights

            avg = pl.sum(avgs)
            stdErr = pl.sum(stds)

            AVG = pl.append(AVG, avg)
            STD = pl.append(STD, stdErr)

            n += nEst
 
        # determine titles for plotting
        if extentScale:
            labell = 'thickness: '+thickness+' '+r'$[\si{\angstrom}]$'
            titlle = 'extent: '+extent+' '+r'$[\si{\angstrom}]$'
        else:
            labell = 'extent: '+extent+' '+r'$[\si{\angstrom}]$'
            titlle = 'thickness: '+thickness+' '+r'$[\si{\angstrom}]$' 

        # add current data to plot
        pl.errorbar(headers, AVG, STD, fmt='o', color=colors[nd], 
                label=labell)
        pl.xlabel(xLab, fontsize=20)
        #pl.ylabel(r'$\langle \Omega^2 \rangle$', fontsize=20)
        pl.ylabel(r'$\langle N \rangle$', fontsize=20)
        pl.title(titlle)
        pl.legend()
        pl.grid(True)

        os.chdir('..')
    
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
