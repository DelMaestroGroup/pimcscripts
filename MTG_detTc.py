# =============================================================================
# Script for plotting <W^2>/3 and superfluid fraction from PIMC data
#
# Author:           Max Graves
# Last Revised:     13-MAR-2013
# =============================================================================

import pylab as pl
import os
import glob

def main():
    
    # list of directories
    dirs = (glob.glob("*N*"))

    nums = pl.array([])
    for direc in dirs:
        nums = pl.append(nums, int(direc[1:]))
    nums = sorted(nums)

    print nums
    n = 0.0
    pl.figure(1)
    ax = pl.subplot(111)
    pl.figure(2)
    bx = pl.subplot(111)
    for N in nums:
        print N
        os.chdir('N'+str(int(N))+'/OUTPUT/')

        Temps = pl.array([])
        Winds = pl.array([])
        Fracs = pl.array([])
        files = (glob.glob("*-super-*"))
        mColor = (0,0.5,1.0*n/(len(dirs)-1))
        n += 1

        for f in files:
            print f
            frac,Wx2,Wy2,Wz2 = pl.loadtxt(f,unpack=True, usecols=(0,1,2,3))
            
            W2 = pl.average(Wx2 + Wy2 + Wz2)
            Winds = pl.append(Winds, W2/3.0)
            Fracs = pl.append(Fracs, pl.average(frac))
            Temps = pl.append(Temps, float(f[10:15]))

        # linear regression
        (ar,br) = pl.polyfit(Temps,Winds,1)
        linY = pl.polyval([ar,br],Temps)
        
        # plot data and best fit line from each set of data
        ax.plot(Temps, Winds, color=mColor, marker='o', linewidth=0,
                markerfacecolor='None', markeredgecolor=mColor)
        ax.plot(Temps,linY, linewidth=2, label='N'+str(N), color=mColor)
        bx.plot(Temps, Fracs, color=mColor, marker='o', linewidth=0,
                markerfacecolor='None', markeredgecolor=mColor,
                label='N'+str(N))
        
        os.chdir('../../')

    # construct universality class
    print Temps
    U1 = 0.516*pl.ones(len(Temps))
    ax.plot(Temps, U1, color='Black',linewidth=5,label='U(1)')
    # plot all
    ax.set_xlabel(r'$T\ (K)$', size=20)
    bx.set_xlabel(r'$T\ (K)$', size=20)
    ax.set_ylabel(r'$\frac{\langle W^2 \rangle}{3}$', rotation='horizontal',
            size=25)
    bx.set_ylabel(r'$\frac{\rho_s}{\rho}$', rotation='horizontal',
            size=25)
    ax.grid(True)
    bx.grid(True)
    ax.legend()
    bx.legend()
    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
