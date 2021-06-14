# Calculating specific heat capacitance and entropy based on the energy and particle
# numbers from the scalar estimator file. First, a fit is performed to the data. It is 
# done in 2 steps: a polynomial is fit to the superfluid region whose upper bound is 
# controlled by the command line parameter --LT. Then, a spline is fit to the rest of 
# data. Those fits are taken to represent the data and the physical quantities are computed
# and plotted based on them. 
#
# Author:           Bohdan K.
# Last Modified:    17-OCT-2013 by Max G.
# =============================================================================

import kevent
import os,sys
import argparse
import pimcscripts.pimchelp as pimchelp
import numpy as np
import pylab as pl
from math import sqrt
from scipy import interpolate
from scipy import integrate
import matplotlib

# =============================================================================
def resizeF(x,size):
    ''' Increase a 1-d array size by adding zeros in front'''
	
    oldSize = len(x) 
    x=np.resize(x,size)
    x[size-oldSize:] = x[:oldSize]
    x[0:size-oldSize] = 0
    return x

# =============================================================================
def ConvertEnergyUnits(EoverN):
    ''' EoverN - energy over number of particles; 
        From Kelvin to Joule/gramm'''

    kb  = 1.380648*10**(-23)  # Boltzman constant
    he4 = 4.002602            # helium mass in amu
    amu = 1.660539*10**(-24)  # atomic unit mass in g
    return EoverN*kb/(he4*amu)

# =============================================================================
def EntropySpline(lowT,highT,spline):   
    '''Numerically calculates the contribution to entropy from the spline fit'''

    step = (highT-lowT)/1000.0
    T = pl.arange(lowT,highT, step)
    cv = interpolate.splev(T,spline,der=1)
    return integrate.simps(cv/T, T)

# =============================================================================
def EntropyPol(lowT,highT,Dpolyf):
    '''Numerically calculates the contribution to entropy from the polynomial fit'''

    step = (highT-lowT)/1000.0
    #print 'lowT:  ',lowT
    #print 'highT: ',highT
    T = pl.arange(lowT,highT, step)
    cv = Dpolyf(T)
    return integrate.simps(cv/T, T)

# =============================================================================
def Entropy(T,Tc,T0,Dpolyf,spline):
    '''Numerically calculates the contribution to entropy from the total fit
       from T0 to T. Tc is taken as a boundary between polynomial and spline fits'''

    if T==T0: return 0
    if T>Tc:  return EntropyPol(T0,Tc,Dpolyf) + EntropySpline(Tc,T,spline)
    if T<Tc:  return EntropyPol(T0,T,Dpolyf)

# =============================================================================
def IntersectionPoint(spline,Dpolyf,IntL,IntH):
    ''' Finds Tc where the later is the temperature where spline and polymial 
        intersect'''

    L = IntL
    H = IntH
    for i in range(1,7):
        ranT = pl.arange(L,H,(H-L)/100.0) 
        for j,T in enumerate(ranT):
            diff =  interpolate.splev(T,spline,der=1) - Dpolyf(T)
            if diff<0:
               H = T
               break
            else:
                L = T
    return (L+H)/2

# =============================================================================
def IntersectionPoint2S(LSpline,HSpline,IntL,IntH):
    ''' Find intersection of 2 spline derivatives between IntL and IntH'''

    L = IntL
    H = IntH
    for i in range(1,7):
        ranT = pl.arange(L,H,(H-L)/100.0) 
        for j,T in enumerate(ranT):
            diff =  interpolate.splev(T,HSpline,der=1) - interpolate.splev(T,LSpline,der=1) 
            if diff<0:
               H = T
               break
            else:
                L = T
    return (L+H)/2

# =============================================================================
def ZeroPoint(Dpolyf,IntPoint):
    '''Finds T0 where the later is the temperature where the polynomial fit 
       intersects the temperature axis'''
 
    L = 0
    H = IntPoint
    for i in range(1,5):
        ranT = pl.arange(L,H+(H-L)/100.0,(H-L)/100.0) 
        for j,T in enumerate(reversed(ranT)):
            diff =  Dpolyf(T)
            if diff<0:
               L = T
               break
            else:
                H = T
    return H #(L+H)/2

# =============================================================================
def PolyFitToData(x,y,ey,Orders):
    '''Fits polynomials of degree Orders to the data'''

    weights = 1.0/ey

    FitCoefs = np.zeros((max(Orders) - min(Orders) + 1,max(Orders) + 1))
    DerCoefs = np.zeros((max(Orders) - min(Orders) + 1,max(Orders))) 

    for i,order in enumerate(Orders):
       polyC  = np.copy(np.polyfit(x, y, order, w=weights))
       DpolyC = np.copy(np.polyder(polyC))

       if len(polyC) < (max(Orders) + 1):
          polyC = resizeF(polyC,max(Orders) + 1)
       FitCoefs[i,] = polyC
       if len(DpolyC) < max(Orders):
          DpolyC = resizeF(DpolyC,max(Orders))
       DerCoefs[i,] = DpolyC
    
    return FitCoefs,DerCoefs

# =============================================================================
def SplineFitToData(x,y,ey,s=-1):
    '''Fits a spline to the data of smoothness s, and extrapolates its values and 
       the values of its derivative at newx'''

    if s == -1: s = len(x)+sqrt(2*len(x))

    spline = interpolate.splrep(x, y,w=1./ey,s=s)
    newx = pl.arange(x[0]-0.25,x[-1], 0.01)
    E  = interpolate.splev(newx,spline,der=0)
    dE = interpolate.splev(newx,spline,der=1)
    
    return E, dE, spline

# =============================================================================
def parseCMD():
    ''' parse the command line.'''
    parser = argparse.ArgumentParser(description='Entropy Integration of E vs T')
    parser.add_argument('fileName', help='Reduced scalar estimator file', 
                        nargs='+')
    parser.add_argument('-t','--LT', 
            help='Estimated critical temperature. \
                    It separates data for the spline and polynomial fits', 
            type=float)
    parser.add_argument('-o','--order', 
            type=int, default=1,
            help='Order of polynomial to be fit')
    parser.add_argument('-s','--smoothness', 
            type=int, default=-1,
            help='Spline smoothness. Default sets it to the recommended value.')
    
    return parser.parse_args()
  
#------------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    args = parseCMD()
    fileNames = args.fileName
    Order = args.order
    Smoothness = args.smoothness
       
    matplotlib.rcParams['text.usetex'] = True
    linestyles = ['-','--']
    colors = [(1,0,0),(0,1,0),(1,1,0),
              (1,0,1),(0,0,0),(0,0,1),
              (1,0.5,0.5),(0.5,0.5,0.5),
              (0.5,0,0),(1,0.5,0)]
    
    # Set up plots with shared axes
    fig1 = pl.figure(1,figsize=(12,12))
    pl.connect('key_press_event',kevent.press)

    ax1 = pl.subplot(311)
    pl.setp(ax1.get_xticklabels(), visible=False)
    pl.ylabel(r'$\mathrm{E\ [K]}$',fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax1.tick_params(axis='both', which='minor', labelsize=20)
    pl.grid()
    
    ax2 = pl.subplot(312,sharex=ax1)
    pl.setp(ax2.get_xticklabels(), visible=False)
    pl.ylabel(r'$\mathrm{C_V}$',fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='minor', labelsize=20)
    pl.grid()

    ax3 = pl.subplot(313,sharex=ax1)
    pl.xlabel(r'$\mathrm{T\ [K]}$',fontsize=20)
    pl.ylabel(r'${\mathrm{Entropy}}$',fontsize=20)
    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.tick_params(axis='both', which='minor', labelsize=20)
    pl.grid()
 
    # do spline fitting for each data file
    for k,fileName in enumerate(fileNames):
        data = np.loadtxt(fileName)
        ATemps = data[:,0]
        #ATemps = data[:,0] # Bohdan's test data
        AEners = data[:,1] 
        #AEners = data[:,5] 
        AEnerE = data[:,2] 
        #AEnerE = data[:,6] 
        weights  = 1/AEnerE
       
        ax1.errorbar(ATemps, AEners, AEnerE, 
                linestyle='None', linewidth=2, marker='D', markeredgewidth=0.8,
                color=colors[k], markeredgecolor = colors[k], 
                markerfacecolor='white', markersize = 6, capsize=3)
        
        lowTRegL  = 0 
        lowTRegH  = len(ATemps[ATemps<args.LT*1.0001])-1
        highTRegL = lowTRegH+1
        highTRegH = len(ATemps)-1

        # Spline fit to T>Tc
        SEners,SdEners,HSpline = SplineFitToData(ATemps[highTRegL:highTRegH],
                AEners[highTRegL:highTRegH],AEnerE[highTRegL:highTRegH],Smoothness)

        # Spline fit to T<Tc
        #SLE,SLcv,LSpline = SplineFitToData(ATemps[lowTRegL:lowTRegH], 
        #        AEners[lowTRegL:lowTRegH],AEnerE[lowTRegL:lowTRegH],Smoothness)

        # Polynomial fit to T<Tc
        SegRange = range(lowTRegL,lowTRegH)
        FitCoefs,DerCoefs = PolyFitToData(ATemps[SegRange],AEners[SegRange],
                AEnerE[SegRange],[Order])

        polyf  = np.poly1d(FitCoefs[0])
        Dpolyf = np.poly1d(DerCoefs[0])

        interT = IntersectionPoint(HSpline,Dpolyf,ATemps[lowTRegH],
                ATemps[highTRegL])
        #interT2S = IntersectionPoint2S(LSpline,HSpline,ATemps[lowTRegH],
        #        ATemps[highTRegL])

        stepT = 0.01
        newLTs  = pl.arange(ATemps[0],interT, stepT)
        highLTs = pl.arange(interT,ATemps[highTRegH], stepT)
        lowLTs = pl.arange(ATemps[0],interT, stepT)

        ax1.plot(newLTs, polyf(newLTs),  linestyle=linestyles[0],
                linewidth=0.5, marker='None',
                color=colors[k], markerfacecolor=colors[k],
                label='Pol. %s order fit for T: %s - %s K' %(
                    Order,ATemps[0],ATemps[lowTRegH]))
        ax1.plot(highLTs, interpolate.splev(highLTs,HSpline,der=0), 
                linestyle=':', linewidth=1.5, 
                color=colors[k], markerfacecolor='white', 
                label='Spline fit for T: %s - %s K; smoothness: %s' %(
                    ATemps[highTRegL],ATemps[highTRegH],Smoothness))
        #ax1.plot(lowLTs, interpolate.splev(lowLTs,LSpline,der=0),
        #        linestyle=linestyles[0],linewidth=0.5,
        #        marker='None',color=colors[k],markerfacecolor=colors[k],
        #                label='Spline fit' )
        
        ax2.plot(newLTs, Dpolyf(newLTs), linestyle=linestyles[0],
                linewidth=2, marker='None', color=colors[k],
                markerfacecolor=colors[k])
        ax2.plot(highLTs, interpolate.splev(highLTs,HSpline,der=1), 
                linestyle='-', linewidth=2, color=colors[k], 
                markerfacecolor='white')
        #ax2.plot(lowLTs, interpolate.splev(lowLTs,LSpline,der=1), 
        #        linestyle=':', linewidth=1.5,
        #        color=colors[k],markerfacecolor='white')
        ticks = ax2.xaxis.get_major_ticks()

        np.insert(ticks,3,2)

        zeroT  = ZeroPoint(Dpolyf,interT)
        ax2.xaxis.set_ticks(np.insert(ax2.xaxis.get_majorticklocs(),0,
            round(zeroT,2)))

        ax2.text(interT,interpolate.splev(interT,HSpline,der=1),
                r'${\mathrm{T_c = %s}}  $' %(round(interT,3)),fontsize=18)
        ax2.plot([interT,interT],[0,interpolate.splev(interT,HSpline,der=1)],
                color='black',linestyle='--',linewidth=0.3)
                
        ax2.plot([zeroT,zeroT],[ax2.yaxis.get_majorticklocs()[0],
            ax2.yaxis.get_majorticklocs()[-1]],color='black',
            linewidth=0.3,linestyle = '--')
        ax2.plot([ax2.xaxis.get_majorticklocs()[0],
            ax2.xaxis.get_majorticklocs()[-1]],[0,0],
            color='black',linewidth=0.3,linestyle = '--')
        #if zeroT < 1: 
        #    newETs = pl.arange(1,ATemps[-1], stepT)
        #    Entrops= [Entropy(T,interT,1,Dpolyf,HSpline) for T in newETs] 
        #else:           
        newETs = pl.arange(zeroT,ATemps[highTRegH], stepT)
        print 'Zero     T:', zeroT
        print 'Critical T',interT
        Entrops= [Entropy(T,interT,zeroT,Dpolyf,HSpline) for T in newETs] 
 
        ax3.plot(newETs, Entrops, linestyle='-', linewidth=2, color=colors[k])
    
    lgd1 = ax1.legend(loc='best')
    lgd3 = ax3.legend(loc='best')
    
    lgd1.draggable(state=True) 
    pl.tight_layout()
    pl.savefig('Helium_Critical_E_Cv_S.pdf', format='pdf', 
            bbox_inches='tight')
    pl.show()

#==============================================================================
if __name__ == "__main__": 
    main()
