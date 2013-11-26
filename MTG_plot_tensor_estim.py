"""Plot Tensor Estimator

Description:
  Performs a density plot of a tensor estimator.  Currently defaults to 3d
  Currently the extents of the simulation cell are hard coded to
  8x8x24 Angstroms.
  
Usage:
  plot_tensor_estimator.py <input-file>

  plot_tensor_estimator.py -h | --help 

Options:
  -h --help                                 Show this screen.
"""

# ===================================================================
# Plot the density of particles in a simulation cell.
#
# Author:       Max Graves
#               Adrian Del Maestro
# Date:         8-NOV-2012
# ===================================================================

import pylab as pl
import argparse
from docopt import docopt
import sys
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

# ===================================================================
def main():

    Lx  = 8.0
    Ly  = 8.0
    Lz  = 24.0

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']

    # Open up the tensor file, and determine the number of grid boxes in each
    # dimension
    inFile = open(fileName,'r')
    N = int(inFile.readline().split()[1])
    inFile.close()

    # Assuming a 3d data set, load the data
    data = pl.loadtxt(fileName).reshape([N,N,N])
        
    # plot histograms in all three projections
    pl.figure(1)
    pl.imshow(pl.sum(data,axis=2)/N, extent=[-4,4,-4,4])
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    #pl.title('Particle Density Projection (X-Y)')
    pl.colorbar(shrink=0.4)

    pl.savefig('xy_densityHistogram.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
    
    pl.figure(2)
    pl.imshow(pl.sum(data,axis=1)/N, extent=[-12,12,-4,4])
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    #pl.title('Particle Density Projection (X-Z)')
    pl.colorbar(shrink=0.4)
 
    pl.savefig('zx_densityHistogram.pdf', format='pdf',
            bbox_inches='tight', transparent=True)  
    
    pl.figure(3)
    pl.imshow(pl.sum(data,axis=0)/N, extent=[-12,12,-4,4])
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$y\  [\AA]$')
    #pl.title('Particle Density Projection (Y-Z)')
    pl.colorbar(shrink=0.4)

    pl.savefig('zy_densityHistogram.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
  
    pl.figure(4)    # all plots together
    #pl.suptitle('Particle Density Projection (X-Y) -- T = '+str(fileName[12:18]), size=20)
    gs = gridspec.GridSpec(8,8)
    xyplot = pl.subplot(gs[2:4,6:8])
    pl.imshow(pl.sum(data,axis=2)/N, extent=[-(Lx/2.0),(Lx/2.0),-(Ly/2.0),(Ly/2.0)] )
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    labels = xyplot.get_xticklabels()
    for label in labels:
        label.set_rotation(90)
    xyplot.yaxis.set_label_position('right')
    xyplot.xaxis.set_label_position('top')
    #pl.clim(0,7.5)

    zxplot = pl.subplot(gs[1:4,:6])
    pl.imshow(pl.sum(data,axis=1)/N, extent=[-(Lz/2.0),(Lz/2.0),-(Lx/2.0),(Lx/2.0)])
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    #pl.clim(0,7.5)
    pl.colorbar(shrink=0.7)

    zyplot = pl.subplot(gs[4:7,1:7])
    pl.imshow(pl.sum(data,axis=0)/N, extent=[-(Lz/2.0),(Lz/2.0),-(Ly/2.0),(Ly/2.0)])
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$y\  [\AA]$')
    #pl.clim(0,11.0)
    pl.colorbar(shrink=0.7)

    pl.savefig('all_densityHistogram.pdf', format='pdf',
            bbox_inches='tight', transparent=True)
 

    pl.show()
          
# ===================================================================
if __name__=="__main__":
    main()
