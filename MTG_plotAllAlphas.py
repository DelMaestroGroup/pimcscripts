import pylab as pl
import os, glob, sys, argparse

def parseCMD():
    ''' Parse the command line. '''
    parser = argparse.ArgumentParser(description='plots density vs Mu for lots of alphas.')
    parser.add_argument('-T', '--Temperature', type=str,
            default='T',
            help='Temperature (for plot title).')
    parser.add_argument('-L', '--textLocation', type=str,
            default='R', choices = 'LR',
            help='Position of text boxes for SVP densities')

    return parser.parse_args()

def main():

    args = parseCMD()

    tLoc = args.textLocation
    Temp = args.Temperature
    if tLoc == 'R':
        boxSubtract = 1
    else:
        boxSubtract = 4

    alphas = glob.glob('*alpha*')
    if alphas==[]:
        sys.exit('Must be in CWD containing alphaN direcs')

    minMus, maxMus = 5000, -5000

    colors = ["#70D44A", "#BE5AD4", "#D04537", "#81D0D5", "#393A2E", "#C49ECA",
              "#C5CB7A", "#523767", "#D39139", "#C8488C", "#CB817A", "#73D999",
              "#6F2836", "#6978CF", "#588569", "#CDC3AD", "#5C7890", "#7F5327",
              "#D0D33D", "#5B7D2E"]
    colors = ["#53B0AD", "#D74C20", "#CD53DA", "#58C038", "#5F4B7A", "#49622A",
              "#CE4379", "#D2912E", "#7970D2", "#749AC9", "#7D5121", "#5EAC72",
              "#CB85AA", "#853B46", "#396465", "#A5A33E", "#D47F5B", "#BA4EA5",
              "#C93F44", "#5D9937"]
    colors =["#CAC5E8", "#E1D273", "#82DFCE", "#F1AF92", "#E1E7CF", "#92C798", 
            "#CDF197", "#DDD199", "#ECB1D1", "#91C7DE", "#E3AF6E", "#ECAAAC", 
            "#CEB29C", "#B2BCA9", "#B0C778", "#DBCDD7", "#B5DEE0", "#A5ECB4", 
            "#C5E6C0", "#E5CCC0"]
    colors = ["#688EAF", "#FC991D", "#7DEB74", "#FA6781", "#8B981D", "#BB7548", 
            "#AD8FE4", "#96E4AA", "#D669B0", "#E1C947", "#A78200", "#7C9FE4", 
            "#957DA6", "#75BF38", "#C3B059", "#51C17A", "#79AEBB", "#2790AC", 
            "#688ECE", "#749DB7"]
    

    n = 0
    for alpha in sorted(alphas):
        os.chdir(alpha)
        a = alpha[-1]
        mus, fDens, fErr, bDens, bErr = pl.loadtxt(
                'JackKnifeData_bipart.dat', unpack=True)

        pl.errorbar(mus, fDens, fErr, fmt='o', 
                label=('Film: '+r'$\alpha = $'+'%s' % a),
                color = colors[n])
        pl.errorbar(mus, bDens, bErr, fmt='o',
                label=('Bulk: '+r'$\alpha = $'+'%s' % a),
                color = colors[n+1])

        # determine max and min values of mu
        if pl.amax(mus) > maxMus:
            maxMus = pl.amax(mus)
        if pl.amin(mus) < minMus:
            minMus = pl.amin(mus)

        os.chdir('..')
        n += 2

    print pl.amax(mus)

    # set up bulk SVP densities for plot
    pl.plot([minMus, maxMus], [0.02198, 0.02198], 'k-', lw=3)
    pl.annotate('3D SVP', xy=(maxMus - boxSubtract, 0.02195),  #xycoords='data',
            xytext=(-50, -30), textcoords='offset points',
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


    pl.xlabel('Chemical Potential [K]', fontsize=20)
    pl.ylabel('Spatial Density '+r'$[\AA^{-d}]$', fontsize=20)
    pl.title('T = %s K' % Temp) 
    pl.legend(loc=2)
    pl.show()


if __name__=='__main__':
    main()
