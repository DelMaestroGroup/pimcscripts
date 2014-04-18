import pylab as pl
import os, glob, sys, argparse, re
from matplotlib import rcParams
from matplotlib.ticker import FuncFormatter

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}', # load siunitx
       r'\sisetup{detect-all}'  # force siunitx to use your fonts
]

def my_formatter(x, pos):
    """
    Add a trailing zero at end of 
    """
    val_str = '{:g}'.format(x)
    
    if len(val_str) == 5:
        return val_str+'0'
    else:
        return val_str+'00'

def main():

    savePDF = False
    solidDens = False

    Lx = 12.0
    Ly = 12.0
    az = 47.5
    ay = 2.0

    bulkVol = 2.0*Lx*Ly*az

    # define some plotting colors    
    colors = ['Navy','DarkViolet','MediumSpringGreen','Fuchsia',
            'Yellow','Maroon','Salmon','Blue']
    
    # -------------------------------------------------------------------------
    # bulk and film density on same plot
    figg = pl.figure(1)
    ax = figg.add_subplot(111)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Spatial Density '+r'$[\si{\angstrom}^{-d}]$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,3.5])
    #pl.ylim([0,0.07])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # bulk SVP density line
    bulkVert = -30
    minMus = -5.5
    maxMus = -0.5
    boxSubtract = 1.6
    pl.plot([minMus, maxMus], [0.02198, 0.02198], 'k-', lw=3)
    pl.annotate('3d SVP', xy=(maxMus - boxSubtract, 0.02195),  #xycoords='data',
            xytext=(-10, bulkVert), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )

    # film SVP density line
    pl.plot([minMus, maxMus], [0.0432, 0.0432], 'k-', lw=3)
    pl.annotate('2d SVP', xy=(maxMus - boxSubtract, 0.0432),  #xycoords='data',
            xytext=(30, -30), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )
    
    if solidDens:
        pl.plot([minMus, maxMus], [0.0248, 0.0248], 'k-', lw=3)
        pl.annotate('HCP solid SVP', xy=(maxMus - boxSubtract, 0.0248),  #xycoords='data',
                xytext=(-10, 30), textcoords='offset points',
                bbox=dict(boxstyle="round", fc="0.8"),
                arrowprops=dict(arrowstyle="->",
                    connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                )
    
    # -------------------------------------------------------------------------
    # bulk density
    figg2 = pl.figure(2)
    ax2 = figg2.add_subplot(111)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Bulk Density '+r'$[\si{\angstrom}^{-3}]$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,0.5])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # set up bulk SVP densities for plot
    bulkVert = 30   # changes text box above (positive) or below (negative) line
    boxSubtract = 1.6
    pl.plot([minMus, maxMus], [0.02198, 0.02198], 'k-', lw=3)
    pl.annotate('3d SVP', xy=(maxMus - boxSubtract, 0.02198),  #xycoords='data',
            xytext=(-50, bulkVert), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )
    
    # -------------------------------------------------------------------------
    # number of particles in film region
    pl.figure(3)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$N_{\text{film}}$', fontsize=20)
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
    
    # -------------------------------------------------------------------------
    # normalized angular winding
    pl.figure(4)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\Omega$', fontsize=20)
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # -------------------------------------------------------------------------
    # film density
    figg5 = pl.figure(5)
    ax5 = figg5.add_subplot(111)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\text{Film Density}\ [\si{\angstrom}^{-2}]$', fontsize=20)
    pl.ylim([0.03,0.05])
    pl.xlim([-5.5,0.5])
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
 
    # film SVP density line
    pl.plot([minMus, maxMus], [0.0432, 0.0432], 'k-', lw=3)
    pl.annotate('2d SVP', xy=(-4, 0.0432),  #xycoords='data',
            xytext=(30, -30), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )
 
    # -------------------------------------------------------------------------
    # superfluid fraction
    pl.figure(6)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\rho_S/\rho$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,0.5])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # -------------------------------------------------------------------------
    # film/bulk densities subplot 
    pl.figure(7)
    ax7a = pl.subplot(211)
    pl.ylabel(r'$\text{Film Density}\ [\si{\angstrom}^{-2}]$', fontsize=20)
    pl.ylim([0.03,0.05])
    #pl.xlim([-5.5,0.5])
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
 
    # film SVP density line
    pl.plot([minMus, maxMus], [0.0432, 0.0432], 'k-', lw=3)
    pl.annotate('2d SVP', xy=(-4, 0.0432),  #xycoords='data',
            xytext=(30, -30), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )
    pl.setp(ax7a.get_xticklabels(), visible=False)
    
    ax7b = pl.subplot(212, sharex=ax7a)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Bulk Density '+r'$[\si{\angstrom}^{-3}]$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,-0.5])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax7b.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # set up bulk SVP densities for plot
    bulkVert = 15   # changes text box above (positive) or below (negative) line
    boxSubtract = 1.6
    pl.plot([minMus, maxMus], [0.02198, 0.02198], 'k-', lw=3)
    pl.annotate('3d SVP', xy=(-3.5, 0.02198),  #xycoords='data',
            xytext=(-50, bulkVert), textcoords='offset points',
            bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="->",
                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
            )
 
    # -------------------------------------------------------------------------
    # number of particles in bulk
    pl.figure(8)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$N_{\text{bulk}}$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,0.5])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)


    # -------------------------------------------------------------------------
    Tvals = glob.glob('*T*')
   
    Markers = ['-o','-d','-*']

    
    for nT,Tval in enumerate(sorted(Tvals)):
        
        os.chdir(Tval)

        T = Tval[1:]

        Svals = glob.glob('S*')
        
        # --- loop through known directory structure --------------------------
        for nS,Sval in enumerate(sorted(Svals)):
            
            os.chdir(Sval)

            Vdirs = glob.glob('*V*')

            print os.getcwd()

            # store bulk separation value
            S = re.search(r'\d+',Sval).group(0)
             
            # get label for plot
            if 'distinguishable' in Sval:
                labell = 'S = '+str(S)+', Boltzmannons, T='+str(T)
            else:
                labell = 'S = '+str(S)+', Bosons, T='+str(T)
                shortlabell = 'S = '+str(S)+', T='+str(T)

            # projected area of film region
            projArea = float(S)*Lx

            # accessible volume in cell
            #accessibleVol = Lx*(Ly*(float(S)+2.0*az) - float(S)*(Ly-2.0*ay))
            accessibleVol = 2.0*Lx*Ly*az + 2.0*ay*Lx*float(S)

            # multiply angular winding by norm. to be fixed in code later 
            omegaNorm = 4.0*(float(S)+1.0*Ly)**2

            # Arrays to hold all data
            Vs = pl.array([])
            Films = pl.array([])
            Bulks = pl.array([])
            filmErrs = pl.array([])
            bulkErrs = pl.array([])
            Omegas = pl.array([])
            omegaErrs = pl.array([])
            NumParts = pl.array([])
            NumPartErrs = pl.array([])
            Supers = pl.array([])
            SuperErrs = pl.array([])

            
            # pass potential shift directories for current S value
            for Vdir in Vdirs:

                os.chdir(Vdir)
                
                # get bipartition file name
                f = glob.glob('*Bipart*')[0]

                # get angular winding file name
                fw = glob.glob('*Ntwind*')[0]
                
                # get estimator file name (includes total number)
                fe = glob.glob('*Estimator*')[0]

                # get superfrac file name
                fs = glob.glob('*Super*')[0]
        
                # build array of film potential shifts from directory names
                Vs = pl.append(Vs,float(Vdir[1:])) 

                # --- Densities -----------------------------------------------
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


                # ---- angular winding ----------------------------------------
                omegaAvg,omegaStd,omegaBins = pl.genfromtxt(fw,
                        unpack=True, usecols=(3,4,5), delimiter=',')

                # get rid of any items which are not numbers..
                omegaBins = omegaBins[pl.logical_not(pl.isnan(omegaBins))]
                omegaStd = omegaStd[pl.logical_not(pl.isnan(omegaStd))]
                omegaAvg = omegaAvg[pl.logical_not(pl.isnan(omegaAvg))]

                # normalize data.
                omegaStd *= omegaNorm
                omegaAvg *= omegaNorm

                weights = omegaBins/pl.sum(omegaBins)

                omegaAvg *= weights
                omegaStd *= weights

                Omega = pl.sum(omegaAvg)
                omegaErr = pl.sum(omegaStd)

                Omegas = pl.append(Omegas, Omega)
                omegaErrs = pl.append(omegaErrs, omegaErr)
               

                # ---- total number -------------------------------------------
                numAvg,numStd,numBins = pl.genfromtxt(fe,
                        unpack=True, usecols=(12,13,14), delimiter=',')

                # get rid of any items which are not numbers..
                numBins = numBins[pl.logical_not(pl.isnan(numBins))]
                numStd = numStd[pl.logical_not(pl.isnan(numStd))]
                numAvg = numAvg[pl.logical_not(pl.isnan(numAvg))]

                weights = numBins/pl.sum(numBins)

                numAvg *= weights
                numStd *= weights

                numPart = pl.sum(numAvg)
                numPartErr = pl.sum(numStd)

                NumParts = pl.append(NumParts, numPart)
                NumPartErrs = pl.append(NumPartErrs, numPartErr)


                # ---- superfluid fraction ------------------------------------
                supAvg,supStd,supBins = pl.genfromtxt(fs,
                        unpack=True, usecols=(0,1,2), delimiter=',')

                # get rid of any items which are not numbers..
                supBins = supBins[pl.logical_not(pl.isnan(supBins))]
                supStd = supStd[pl.logical_not(pl.isnan(supStd))]
                supAvg = supAvg[pl.logical_not(pl.isnan(supAvg))]

                # normalize data.
                #supStd /= (1.0*accessibleVol)
                #supAvg /= (1.0*accessibleVol)

                weights = supBins/pl.sum(supBins)

                supAvg *= weights
                supStd *= weights

                supPart = pl.sum(supAvg)
                supPartErr = pl.sum(supStd)

                Supers = pl.append(Supers, supPart)
                SuperErrs = pl.append(SuperErrs, supPartErr)


                os.chdir('..')


            # Sort data in order of increasing chemical potential.
            # This is another bit of magical Python juju.
            Vs, Films, Bulks, filmErrs, bulkErrs, Omegas, omegaErrs,\
                    NumParts, NumPartErrs, Supers, SuperErrs = pl.asarray(
                    zip(*sorted(zip(Vs, Films, Bulks, filmErrs, bulkErrs, 
                        Omegas, omegaErrs, NumParts, NumPartErrs, Supers, 
                        SuperErrs))))


            pl.figure(1) 
            pl.errorbar(Vs, Films, filmErrs, fmt=Markers[nT], color=colors[nS],
                    label=labell+', 2d',
                    markersize=8)
            pl.errorbar(Vs, Bulks, bulkErrs, fmt=Markers[nT], mec=colors[nS],
                    label=labell+', 3d', mfc='None', 
                    markersize=8)

            pl.figure(2) 
            pl.errorbar(Vs, Bulks, bulkErrs, fmt = Markers[nT], color=colors[nS],
                    label=labell, markersize=8)

            pl.figure(3)
            pl.errorbar(Vs, Films*projArea, fmt=Markers[nT], color=colors[nS],
                    label = labell, markersize=8)

            pl.figure(4)
            pl.errorbar(Vs, Omegas, omegaErrs, fmt=Markers[nT], color=colors[nS],
                    label = labell, markersize=8)

            pl.figure(5)
            pl.errorbar(Vs, Films, filmErrs, fmt=Markers[nT], color=colors[nS],
                    label = labell, markersize=8)
     
            pl.figure(6)
            pl.errorbar(Vs, Supers, SuperErrs, fmt=Markers[nT], color=colors[nS],
                    label = labell, markersize=8)

            pl.figure(7)
            ax7a.errorbar(Vs, Films, filmErrs, fmt=Markers[nT], color=colors[nS],
                    label=shortlabell, markersize=8)
            ax7b.errorbar(Vs, Bulks, bulkErrs, fmt = Markers[nT], color=colors[nS],
                    label=shortlabell, markersize=8)
 
            pl.figure(8) 
            pl.errorbar(Vs, Bulks*bulkVol, bulkErrs*bulkVol, 
                    fmt = Markers[nT], color=colors[nS],
                    label=labell, markersize=8)


            os.chdir('..')
        
        os.chdir('..')
   
    pl.figure(1)
    pl.legend(loc=1)
    pl.tight_layout()
    if savePDF:
        pl.savefig('densities_vs_potentialShift_allS_8APR.pdf', format='pdf',
                bbox_inches='tight')

    pl.figure(2)
    pl.legend(loc=1)
    pl.tight_layout()
 
    pl.figure(3)
    pl.legend(loc=1)
    pl.tight_layout()
  
    pl.figure(4)
    pl.legend(loc=1)
    pl.tight_layout()
   
    pl.figure(5)
    pl.legend(loc=1)
    pl.tight_layout()

    pl.figure(6)
    pl.legend(loc=1)
    pl.tight_layout()
 
    pl.figure(7)
    ax7a.legend(loc=1)
    major_formatter = FuncFormatter(my_formatter)
    ax7a.yaxis.set_major_formatter(major_formatter)
    pl.tight_layout()
 
    pl.figure(8)
    pl.legend(loc=1)
    pl.tight_layout()
  

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
