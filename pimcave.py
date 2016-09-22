#! /usr/bin/env python

# Adrian Del Maestro
# 07.20.2009
'''pimcave.py

Description:
    Generates averages from pimc output data. 

Usage: pimcave.py [ -s <skip>] (<file>...)


Options:
  -h, --help                Show this help message and exit
  -s <skip>, --skip=<skip>  How many input lines should we skip? [default: 0]

'''

from __future__ import print_function
from docopt import docopt
import numpy as np

# -----------------------------------------------------------------------------
def stats(data):
    '''Return the average and standard error in data. '''
    ave = np.average(data)
    ave2 = np.average(data*data)
    err = np.sqrt(np.abs(ave2-ave**2)/(1.0*data.size-1.0) ) 

    return ave,err

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():

    # parse the command line options
    args = docopt(__doc__)

    fileNames = args['<file>']
    skip = int(args['--skip'])

    for fileName in fileNames:

        # get the pimcid
        pimcid = fileName.split('-')[-1].rstrip('.dat')

        # We check to see if we are dealing with the one body density matrix
        if fileName.find('obdm') != -1:
            normalize = True

        # open the file and determine how many measurements there are
        estData = np.genfromtxt(fileName,names=True,skip_header=1, deletechars="")
        numLines = estData.size
        
        # If we have data, compute averages and error
        if numLines-skip > 0:
            print('# PIMCID %s' % pimcid)
            print('# Number Samples %6d' % (numLines-skip))
            for name in estData.dtype.names:
                ave,err = stats(estData[name][skip:])
                print('%-16s%12.5f\t%12.5f\t%5.2f' %
                      (name,ave,err,100.0*np.abs(err/ave)))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
