"""Plot Tensor Estimator

Description:
  Performs a density plot of a tensor estimator.  Currently defaults to 3d
  
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

# ===================================================================
def main():

    # Read in the command line arguments
    args = docopt(__doc__)
    fileName = args['<input-file>']

    # Open up the tensor file, and determine the number of grid boxes in each
    # dimension
    inFile = open(fileName,'r')
    oldFile=True
    # get size of 
    if oldFile:
        print 'Chose to work on an older position file'
        for i, line in enumerate(inFile):
            if i==1:
                N = int(line.split()[-1])+1
                break
    else:
        print 'Chose to work on a newer position file'
        N = int(inFile.readline().split()[1])
    inFile.close()

    # Assuming a 3d data set, load the data
    data = pl.loadtxt(fileName).reshape([N,N,N])
        
    # plot histograms in all three projections
    pl.figure(1)
    pl.imshow(pl.sum(data,axis=2))
    pl.xlabel(r'$y\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Y)')
    pl.colorbar(shrink=0.4)

    pl.savefig('xy_densityHistogram.pdf', format='pdf',
            bbox_inches='tight')
    
    pl.figure(2)
    pl.imshow(pl.sum(data,axis=1))
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$x\  [\AA]$')
    pl.title('Particle Density Projection (X-Z)')
    pl.colorbar(shrink=0.4)
 
    pl.savefig('zx_densityHistogram.pdf', format='pdf',
            bbox_inches='tight')  
    
    pl.figure(3)
    pl.imshow(pl.sum(data,axis=0))
    pl.xlabel(r'$z\  [\AA]$')
    pl.ylabel(r'$y\  [\AA]$')
    pl.title('Particle Density Projection (Y-Z)')
    pl.colorbar(shrink=0.4)

    pl.savefig('zy_densityHistogram.pdf', format='pdf',
            bbox_inches='tight')
   
    pl.show()
          
# ===================================================================
if __name__=="__main__":
    main()
