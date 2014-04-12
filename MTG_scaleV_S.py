import pylab as pl
import os, glob, sys, argparse, re
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}', # load siunitx
       r'\sisetup{detect-all}'  # force siunitx to use your fonts
]

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(
            description='plots density vs Mu for lots of bulk Separations.')
    parser.add_argument('-T', '--Temperature', type=str,
            default='T',
            help='Temperature (for plot title).')


    return parser.parse_args()

def main():

    savePDF = False
    solidDens = False

    args = parseCMD()
    Temp = args.Temperature

    Lx = 12.0
    Ly = 12.0
    az = 47.5
    ay = 2.0

    Svals = glob.glob('S*')

    # define some plotting colors    
    colors = ['Salmon','Blue','DarkViolet','MediumSpringGreen','Fuchsia',
            'Yellow','Maroon']
        
    # set up figure that displays all data
    figg = pl.figure(1)
    ax = figg.add_subplot(111)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Spatial Density '+r'$[\si{\angstrom}^{-d}]$', fontsize=20)
    pl.grid(True)
    pl.xlim([-5.5,-0.5])
    #pl.ylim([0,0.07])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # set up bulk SVP densities for plot
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
    
    # plot just the bulk density as a function of chemical potential
    figg2 = pl.figure(2)
    ax2 = figg2.add_subplot(111)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Spatial Density '+r'$[\si{\angstrom}^{-d}]$', fontsize=20)
    pl.grid(True)
    #pl.xlim([0.5,3.0])
    #pl.ylim([0,0.07])
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
    
    # number of particles in film region vs chemical potential
    pl.figure(3)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel('Number of Particles', fontsize=20)
    pl.grid(True)
    #pl.xlim([0.5,3.0])
    #pl.ylim([0,2])
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
    
    # normalized angular winding vs chemical potential
    pl.figure(4)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\Omega$', fontsize=20)
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)

    # actual density in cell
    pl.figure(5)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\rho\ [\si{\angstrom}^{-3}]$', fontsize=20)
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
 
    # superfluid fraction
    pl.figure(6)
    pl.xlabel(r'$\text{Potential Shift}\ [K]$', fontsize=20)
    pl.ylabel(r'$\rho_S/\rho$', fontsize=20)
    pl.grid(True)
    pl.tick_params(axis='both', which='major', labelsize=16)
    pl.tick_params(axis='both', which='minor', labelsize=16)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].set_visible(False)
    
    # --- loop through known directory structure ------------------------------
    for nS, Sval in enumerate(sorted(Svals)):
        
        os.chdir(Sval)

        Vdirs = glob.glob('*V*')

        print os.getcwd()

        # store bulk separation value
        S = re.search(r'\d+',Sval).group(0)
         
        # get label for plot
        if 'distinguishable' in Sval:
            labell = 'S = '+str(S)+', Boltzmannons'
        else:
            labell = 'S = '+str(S)+', Bosons'

        # projected area of film region
        projArea = float(S)*Lx

        # accessible volume in cell
        accessibleVol = Lx*Ly*(2.0*az+float(S)) - Lx*float(S)*(Ly-2.0*ay)

        # multiply angular winding by normalization to be fixed in code later 
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

            # --- Densities ---------------------------------------------------
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

            
            # ---- angular winding --------------------------------------------
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
           

            # ---- total number -----------------------------------------------
            numAvg,numStd,numBins = pl.genfromtxt(fe,
                    unpack=True, usecols=(12,13,14), delimiter=',')

            # get rid of any items which are not numbers..
            numBins = numBins[pl.logical_not(pl.isnan(numBins))]
            numStd = numStd[pl.logical_not(pl.isnan(numStd))]
            numAvg = numAvg[pl.logical_not(pl.isnan(numAvg))]

            # normalize data.
            numStd /= (1.0*accessibleVol)
            numAvg /= (1.0*accessibleVol)

            weights = numBins/pl.sum(numBins)

            numAvg *= weights
            numStd *= weights

            numPart = pl.sum(numAvg)
            numPartErr = pl.sum(numStd)

            NumParts = pl.append(NumParts, numPart)
            NumPartErrs = pl.append(NumPartErrs, numPartErr)


            # ---- superfluid fraction ----------------------------------------
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


        # Sort data in order of increasing chemical potential. I love Python.
        # This is another bit of magical Python juju.
        Vs, Films, Bulks, filmErrs, bulkErrs, Omegas, omegaErrs, NumParts, NumPartErrs, Supers, SuperErrs = pl.asarray(
                zip(*sorted(zip(Vs, Films, Bulks, filmErrs, bulkErrs, 
                    Omegas, omegaErrs, NumParts, NumPartErrs, Supers, SuperErrs))))


        pl.figure(1) 
        pl.errorbar(Vs, Films, filmErrs, fmt='--o', color=colors[nS],
                label=labell+', 2d',
                markersize=8)
        pl.errorbar(Vs, Bulks, bulkErrs, fmt = '-d', color=colors[nS],
                label=labell+', 3d',
                markersize=8)

        pl.figure(2) 
        pl.errorbar(Vs, Bulks, bulkErrs, fmt = '-d', color=colors[nS],
                label=labell+', 3d',
                markersize=8)

        pl.figure(3)
        pl.errorbar(Vs, Films*projArea, fmt='-o', color=colors[nS],
                label = labell+', 2d',
                markersize=8)

        pl.figure(4)
        pl.errorbar(Vs, Omegas, omegaErrs, fmt='-o', color=colors[nS],
                label = labell, markersize=8)

        pl.figure(5)
        pl.errorbar(Vs, NumParts, NumPartErrs, fmt='-o', color=colors[nS],
                label = labell, markersize=8)
 
        pl.figure(6)
        pl.errorbar(Vs, Supers, SuperErrs, fmt='-o', color=colors[nS],
                label = labell, markersize=8)
       
        os.chdir('..')

   
    pl.figure(1)
    pl.legend(loc=1)
    if savePDF:
        pl.savefig('densities_vs_potentialShift_allS_8APR.pdf', format='pdf',
                bbox_inches='tight')

    pl.figure(2)
    pl.legend(loc=1)
 
    pl.figure(3)
    pl.legend(loc=1)
  
    pl.figure(4)
    pl.legend(loc=1)
   
    pl.figure(5)
    pl.legend(loc=1)

    pl.figure(6)
    pl.legend(loc=1)
  
    #pl.savefig('density_vs_mu_allAlphas.pdf', format='pdf',
    #        bbox_inches='tight')
    #pl.savefig('density_vs_mu_allAlphas_trans.pdf', format='pdf',
    #        bbox_inches='tight', transparent=True)
    
    # SUPERFLUID STIFFNESS
    '''fig2,ax2 = pl.subplots(1, figsize=(8,6.5))

    n = 0
    for Sval in sorted(Svals):
        os.chdir(Sval)
        a = Sval[-1]
        if Sval[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\Sval = $'+'%s: ' % a
        
        mus, stiff, stiffErr = pl.loadtxt(
                'JackKnifeData_super.dat', unpack=True,
                usecols=(0,1,2))

        pl.errorbar(mus, stiff, stiffErr, fmt='o', 
                label=lab,
                color = colors[n], markeredgecolor='DarkSlateGray',
                markersize=8)

        # determine max and min values of mu
        if pl.amax(mus) > maxMus:
            maxMus = pl.amax(mus)
        if pl.amin(mus) < minMus:
            minMus = pl.amin(mus)

        os.chdir('..')
        n += 1

    if Temp == 'T':
        pl.xlabel('Temperature [K]', fontsize=16)
    else:
        pl.xlabel('Chemical Potential [K]', fontsize=16)
    pl.ylabel('Superfluid Fraction', fontsize=16)
    if Temp != 'T':
        pl.title('T = %s K' % Temp)
    else:
        pl.title(r'$\mu\ =\ $'+ChemPot+' K')
    pl.legend()
    
    pl.savefig('superFrac_vs_mu_allAlphas.pdf', format='pdf',
            bbox_inches='tight')
  
    pl.savefig('superFrac_vs_mu_allAlphas_trans.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
 
    '''
    '''
    # WINDING NUMBER COMPONENTS
    fig3,ax3 = pl.subplots(1, figsize=(8,6.5))

    n = 0
    for Sval in sorted(Svals):
        os.chdir(Sval)
        a = Sval[-1]
        if Sval[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\Sval = $'+'%s: ' % a
        
        mus, Wx2, Wx2Err, Wy2, Wy2Err, Wz2, Wz2Err = pl.loadtxt(
                'JackKnifeData_super.dat', unpack=True,
                usecols=(0,3,4,5,6,7,8))

        pl.errorbar(mus, Wx2, Wx2Err, fmt='o', 
                label=(lab+': '+r'$\langle W_x^2 \rangle$'),
                color = colors[n], markeredgecolor='DarkSlateGray',
                markersize=8)

        pl.errorbar(mus, Wy2, Wy2Err, fmt='v', 
                label=(lab+': ' +r'$\langle W_y^2 \rangle$'),
                color = colors[n], markeredgecolor='DarkSlateGray',
                markersize=8)

        pl.errorbar(mus, Wz2, Wz2Err, fmt='s', 
                label=(lab+': '+r'$\langle W_z^2 \rangle$' ),
                color = colors[n], markeredgecolor='DarkSlateGray',
                markersize=8)

        # determine max and min values of mu
        if pl.amax(mus) > maxMus:
            maxMus = pl.amax(mus)
        if pl.amin(mus) < minMus:
            minMus = pl.amin(mus)

        os.chdir('..')
        n += 1

    if Temp == 'T':
        pl.xlabel('Temperature [K]', fontsize=16)
    else:
        pl.xlabel('Chemical Potential [K]', fontsize=16)
    pl.ylabel(r'$\langle W_i^2 \rangle$', fontsize=16)
    if Temp != 'T':
        pl.title('T = %s K' % Temp)
    else:
        pl.title(r'$\mu\ =\ $'+ChemPot+' K')
    pl.legend()
    
    pl.savefig('windingNumbers_vs_mu_allAlphas.pdf', format='pdf',
            bbox_inches='tight')
     
    pl.savefig('windingNumbers_vs_mu_allAlphas_trans.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
    
    # W_z^2
    fig4,ax4 = pl.subplots(1, figsize=(8,6.5))

    n = 0
    for Sval in sorted(Svals):
        os.chdir(Sval)
        a = Sval[-1]
        if Sval[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\Sval = $'+'%s: ' % a
        
        mus, Wz2, Wz2Err = pl.loadtxt(
                'JackKnifeData_super.dat', unpack=True,
                usecols=(0,7,8))

        pl.errorbar(mus, Wz2, Wz2Err, fmt='s', 
                #label=(r'$\Sval = $'+'%s: ' % a),
                label=lab,
                color = colors[n], markeredgecolor='DarkSlateGray',
                markersize=8)

        # determine max and min values of mu
        if pl.amax(mus) > maxMus:
            maxMus = pl.amax(mus)
        if pl.amin(mus) < minMus:
            minMus = pl.amin(mus)

        os.chdir('..')
        n += 1

    if Temp == 'T':
        pl.xlabel('Temperature [K]', fontsize=16)
    else:
        pl.xlabel('Chemical Potential [K]', fontsize=16)
    pl.ylabel(r'$\langle W_z^2 \rangle$', fontsize=16)
    if Temp != 'T':
        pl.title('T = %s K' % Temp)
    else:
        pl.title(r'$\mu\ =\ $'+ChemPot+' K')
    pl.legend()
    
    pl.savefig('windingZ_vs_mu_allAlphas.pdf', format='pdf',
            bbox_inches='tight')
 
    pl.savefig('windingZ_vs_mu_allAlphas_trans.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
    '''

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()