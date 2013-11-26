import pylab as pl
import os, glob, sys, argparse
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='plots density vs Mu for lots of alphas.')
    parser.add_argument('-T', '--Temperature', type=str,
            default='T',
            help='Temperature (for plot title).')
    parser.add_argument('-u', '--ChemPot', type=str,
            default='u',
            help='Chemical Potential (for plot title).')
    parser.add_argument('-L', '--textLocation', type=str,
            default='R', choices = 'LR',
            help='Position of text boxes for SVP densities')
    parser.add_argument('-x', '--excVol', action='store_true',
            dest='excVol', default=False,
            help='Include if all data has excluded volume.')


    return parser.parse_args()

def main():

    args = parseCMD()

    tLoc = args.textLocation
    Temp = args.Temperature
    ChemPot = args.ChemPot
    excVol = args.excVol

    if tLoc == 'R':
        boxSubtract = 1
    else:
        boxSubtract = 4

    alphas = glob.glob('*alpha*')
    if alphas==[]:
        sys.exit('Must be in CWD containing alphaN direcs')

    minMus, maxMus = 5000, -5000

    colors = ['Salmon','Blue','DarkViolet','MediumSpringGreen','Fuchsia',
            'Yellow','Maroon']
    
    if excVol:
        fig,ax = pl.subplots(1, figsize=(12,8.5))
        n = 0
        for alpha in sorted(alphas):
            os.chdir(alpha)
            a = alpha[-1]
            mus, fDens, fErr, bDens, bErr = pl.loadtxt(
                    'JackKnifeData_bipart.dat', unpack=True)

            pl.errorbar(mus, fDens, fErr, fmt='^', 
                    label=('Film: '+r'$\alpha = $'+'%s' % a),
                    color = colors[n], markeredgecolor='DarkSlateGray',
                    markersize=8)
            pl.errorbar(mus, bDens, bErr, fmt='o',
                    label=('Bulk: '+r'$\alpha = $'+'%s' % a),
                    color = colors[n], markeredgecolor='DarkSlateGray',
                    markersize=8)

            # determine max and min values of mu
            if pl.amax(mus) > maxMus:
                maxMus = pl.amax(mus)
            if pl.amin(mus) < minMus:
                minMus = pl.amin(mus)

            os.chdir('..')
            n += 1

        # set up bulk SVP densities for plot
        if Temp != 'T':
            bulkVert = -30
        else:
            bulkVert = 30
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


        if Temp == 'T':
            pl.xlabel('Temperature [K]', fontsize=16)
        else:
            pl.xlabel('Chemical Potential [K]', fontsize=16)
        pl.ylabel('Spatial Density '+r'$[\AA^{-d}]$', fontsize=16)
        if Temp != 'T':
            pl.title('T = %s K' % Temp)
        else:
            pl.title(r'$\mu\ =\ $'+ChemPot+' K')
        pl.legend(loc=2)
        
        pl.savefig('density_vs_mu_allAlphas.pdf', format='pdf',
                bbox_inches='tight')
    
    # SUPERFLUID STIFFNESS
    fig2,ax2 = pl.subplots(1, figsize=(12,8.5))

    n = 0
    for alpha in sorted(alphas):
        os.chdir(alpha)
        a = alpha[-1]
        if alpha[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\alpha = $'+'%s: ' % a
        
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
 
    # WINDING NUMBER COMPONENTS
    fig3,ax3 = pl.subplots(1, figsize=(12,8.5))

    n = 0
    for alpha in sorted(alphas):
        os.chdir(alpha)
        a = alpha[-1]
        if alpha[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\alpha = $'+'%s: ' % a
        
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
    
    # W_z^2
    fig4,ax4 = pl.subplots(1, figsize=(12,8.5))

    n = 0
    for alpha in sorted(alphas):
        os.chdir(alpha)
        a = alpha[-1]
        if alpha[-1] == '0':
            lab = 'Bulk'
        else:
            lab = r'$\alpha = $'+'%s: ' % a
        
        mus, Wz2, Wz2Err = pl.loadtxt(
                'JackKnifeData_super.dat', unpack=True,
                usecols=(0,7,8))

        pl.errorbar(mus, Wz2, Wz2Err, fmt='s', 
                #label=(r'$\alpha = $'+'%s: ' % a),
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

    pl.show()

# =============================================================================
if __name__=='__main__':
    main()
