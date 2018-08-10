#! /usr/bin/env python

# Author:      2009-07-20 Adrian Del Maestro
# Modified:    2017-07-17 Nathan Nichols
# =============================================================================

'''pimcave.py

Description:
    Generates averages from pimc output data. 

Usage: pimcave.py [ -s <skip> -r] (<file>...)


Options:
  -h, --help                Show this help message and exit
  -s <skip>, --skip=<skip>  How many input lines should we skip? [default: 0]
  -r, --repeated_header     Deal with repeated headers
  -l, --header_lines        Number of header lines [default: 1]

'''

import argparse
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
    parser = argparse.ArgumentParser(description='Generates averages from pimc output data.')
    parser.add_argument('-s', '--skip', type=int, dest='skip', default = 0, 
                        help='How many input lines should we skip? [default: 0]')
    parser.add_argument('-l', '--header_lines', type=int, dest='header_lines', 
                        default = 1, help='How many header lines to skip? [default: 1]')
    parser.add_argument('file', type=str, nargs='+', 
                        help='File or files to average.')
    parser.add_argument('-r','--repeated_header', 
                        help='deal with duplicate headers', action='store_true')

    args = parser.parse_args()

    fileNames = args.file
    skip = args.skip
    header_lines = args.header_lines

    for fileName in fileNames:

        # get the pimcid
        pimcid = fileName[-40:].rstrip('.dat')
        
        # We check to see if we are dealing with the one body density matrix
        if fileName.find('obdm') != -1:
            normalize = True

        # open the file and determine how many measurements there are
        estData = np.genfromtxt(fileName,names=True,skip_header=header_lines, deletechars="")
        numLines = estData.size
        
        # If we have data, compute averages and error
        if numLines-skip > 0:
            print('# PIMCID %s' % pimcid)
            print('# Number Samples %6d' % (numLines-skip))
            for name in estData.dtype.names:
                ave,err = stats(estData[name][skip:])
                if args.repeated_header:
                    cname = name.partition("_")[0]
                else:
                    cname = name
                print('%-16s%12.5f\t%12.5f\t%5.2f' %
                      (cname,ave,err,100.0*np.abs(err/ave)))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
