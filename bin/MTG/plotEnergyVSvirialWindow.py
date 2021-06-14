import pylab as pl
import argparse
from matplotlib import rcParams

# set up latex fonts
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True



def parseCMD():
    ''' parse the command line. '''
    parser = argparse.ArgumentParser(description='plot E vs. W')
    parser.add_argument('fileName', help='Reduced estimator file')
    parser.add_argument('-T','--temp',type=str, help='temperature')
    parser.add_argument('-N','--num',type=str, help='num Particles')

    return parser.parse_args()

def main():

    args = parseCMD()
    temp = args.temp
    num = args.num

    W, EN, ENerr = pl.loadtxt(args.fileName, unpack=True, usecols=(0,17,18))

    pl.errorbar(W,EN,ENerr, color='Teal', fmt='o')
    pl.ylabel('Energy per Particle [K]',fontsize=20)
    pl.xlabel('Virial Window Size',fontsize=20)
    pl.grid()
    
    pl.savefig('energy_vs_virialWindow_3DHe4_N'+num+'_T'+temp+'_bin100000.pdf',
            format='pdf', bbox_inches='tight')
    
    pl.show()


if __name__=='__main__':
    main()
