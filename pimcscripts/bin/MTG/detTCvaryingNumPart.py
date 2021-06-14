# =============================================================================
# Script for plotting <W^2>/3 and superfluid fraction from PIMC data
#
# Author:           Max Graves
# Last Revised:     13-MAR-2013
# =============================================================================

import pylab as pl
import os, glob, sys, subprocess

def findExe(redOnePath):
    """ 
    checks for reduce-one.py
    """
    if os.path.isfile(redOnePath):
        print "reduce_one was found"
        pass
    else:
        print 'Could not find reduce-one.py\n'
        print 'This file needs to be in the same'
        print 'directory as current script.'
        sys.exit()

# =============================================================================
def main():
    
    # list of directories
    dirs = (glob.glob("N*"))

    # location of reduce-one.py
    redOnePath = '../reduce-one.py'
    findExe(redOnePath)
   
    # create array of different N values
    nums = pl.array([])
    for direc in dirs:
        nums = pl.append(nums, int(direc[1:]))
    nums = sorted(nums)

    # create reduced files in all N.. directories
    for num in nums:
        os.chdir('./N'+str(int(num))+'/OUTPUT/')
        command = ("python ../../"+redOnePath+" -r T")
        subprocess.check_call(command, shell=True)
        print "finished for num: ",num
        os.chdir('../../')

    n = 0.0
    pl.figure(1)
    ax = pl.subplot(111)
    pl.figure(2)
    bx = pl.subplot(111)
    for N in nums:
        print N
        os.chdir('N'+str(int(N))+'/OUTPUT/')

        redFile = 'super-T-reduce.dat'
        mColor = (0,0.5,1.0*n/(len(dirs)-1))
        n += 1

        Temps,Fracs,FracErr,Wx2,Wx2Err,Wy2,Wy2Err,Wz2,Wz2Err = pl.loadtxt(
                redFile, unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
            
        Winds = (Wx2 + Wy2 + Wz2)/3.0
        WindsErr = (Wx2Err + Wy2Err + Wz2Err)/3.0

        # regression
        (ar,br) = pl.polyfit(Temps,Winds,1)
        linY = pl.polyval([ar,br],Temps)
        
        # plot data and best fit line from each set of data
        ax.errorbar(Temps, Winds, yerr=WindsErr, color=mColor, fmt='o', 
                #linewidth=0,
                #markerfacecolor='None', 
                markeredgecolor=mColor)
        ax.plot(Temps,linY, linewidth=2, label='N'+str(N), color=mColor)
        bx.errorbar(Temps, Fracs, color=mColor, marker='o', linewidth=0,
                markerfacecolor='None', markeredgecolor=mColor,
                label='N'+str(N))
        
        os.chdir('../../')

    # construct universality class
    U1 = 0.516*pl.ones(len(Temps))
    ax.plot(Temps, U1, color='Black',linewidth=5,label='U(1)')
    # plot all
    ax.set_xlabel(r'$T\ (K)$', size=20)
    bx.set_xlabel(r'$T\ (K)$', size=20)
    ax.set_ylabel(r'$\frac{\langle W^2 \rangle}{3}$', rotation='horizontal',
            size=25)
    ax.set_ylim(0.3)
    bx.set_ylabel(r'$\frac{\rho_s}{\rho}$', rotation='horizontal',
            size=25)
    ax.annotate(r'$T_c^{\ (exp)}$', xy=(2.177,0.3),  xycoords='data',
            xytext=(-80, 20), textcoords='offset points', fontsize=15,
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
            connectionstyle="angle,angleA=0,angleB=120,rad=0"),
            )
    ax.grid(True)
    bx.grid(True)
    ax.legend()
    bx.legend()
    pl.show()
    

# =============================================================================
if __name__=='__main__':
    main()
