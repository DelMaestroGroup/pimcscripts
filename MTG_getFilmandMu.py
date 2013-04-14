# =============================================================================
# Script for plotting density of 2d and 3d regions for Gasparini project
# over a range of chemical potentials.  Does this for multiple values of
# film thickness.
#
# THIS IS QUITE FRAGILE, NEEDS TO BE MADE MORE ROBUST! 
#
# Author:           Max Graves
# Last Revised:     13-MAR-2013
# =============================================================================

import pylab as pl
import os, glob, sys, subprocess

# =============================================================================
def main():
    
    # list of directories
    dirs = (glob.glob("*film*"))

    # create array of different film thickness values
    films = pl.array([])
    for direc in dirs:
        films = pl.append(films, direc[4:8])
    films = sorted(films)

    # create reduced files in all film.. directories
    n = 1
    for num in films:  
        os.chdir('./film'+str(num)+'angs/OUTPUT/')
        alll = glob.glob("*bipart*")
        avg2D, avg3D, mus = pl.array([]), pl.array([]), pl.array([])
        for a in alll:
            try:
                two, three = pl.loadtxt(a,unpack=True)
                mus = pl.append(mus, float(a[31:39]))
                avg2D = pl.append(avg2D, pl.average(two))
                avg3D = pl.append(avg3D, pl.average(three))
            except:
                pass

        pl.figure(n)
        pl.plot(mus, avg2D, marker='o', linewidth=0,
                color='Yellow', label='film: '+str(num)+' '+r'$\AA$')
        pl.plot(mus, avg3D, marker='o', linewidth=0,
                color='ForestGreen', label='bulk')
        pl.ylabel('Density', fontsize=20)
        pl.xlabel(r'$\mu\ [K]$', fontsize=20)
        pl.legend(loc=2)
        os.chdir('../../')
        n += 1
    
    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
