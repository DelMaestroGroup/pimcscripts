# =============================================================================
# Script for plotting specific heat over a range of temperatures from
# a set of gce-virial estimator files.
#
# THIS IS QUITE FRAGILE, NEEDS TO BE MADE MORE ROBUST! 
#
# Author:           Max Graves
# Last Revised:     14-JUL-2013
# =============================================================================

import pylab as pl
import os, glob, sys, subprocess, argparse

def parseCMD():
    ''' Setup the command line parser options. '''
    
    desc = "Performs a cumulative average plot of raw Monte Carlo data"
    parser = argparse.ArgumentParser(description=desc) 
    parser.add_argument('-d', '--dir', type=str,
            help="path to the gce-virial files.")

    return parser.parse_args() 

# =============================================================================
def main():

    args = parseCMD()

    os.chdir(args.dir)

    # list of directories
    dirs = (glob.glob("*film*"))

    # create array of different film thickness values
    films = pl.array([])
    for direc in dirs:
        films = pl.append(films, direc[4:8])
    films = sorted(films)

    # create reduced files in all film.. directories
    temps = pl.array([])
    Cvs = pl.array([])
    n = 1
    for num in films:  
        os.chdir('./film'+str(num)+'angs/OUTPUT/')
        alll = glob.glob("*virial*")
        E, E2, dEdB = pl.array([]), pl.array([]), pl.array([])
        for a in alll:
            print a
            temps = pl.append(temps,float(a[11:17]))
            try:
                #E, E2, dEdB = pl.loadtxt(a, usecols=(0,4,5))
                E, E2, dEdB = pl.loadtxt(a, unpack=True, usecols=(3,7,8))
                Eavg = pl.average(E)
                E2avg = pl.average(E2)
                dEdBavg = pl.average(dEdB)
                Eavg2 = Eavg*Eavg
                Cv = Eavg2 - E2avg -dEdBavg
                print Cv
                Cvs = pl.append(Cvs, Cv)
            except:
                print 'passed'
                pass

        pl.figure(n)
        pl.plot(temps, Cvs, marker='o', linewidth=0,
                color='Violet', label='film: '+str(num)+' '+r'$\AA$')
        pl.ylabel(r'$C_v / k_B \beta^2$', fontsize=20)
        pl.xlabel(r'$\mu\ [K]$', fontsize=20)
        pl.legend(loc=2)
        os.chdir('../../')
        n += 1
    
    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
