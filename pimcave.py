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
  -l, --header_lines        Number of header lines to skip

'''

import argparse
import numpy as np
import sys
import subprocess

# -----------------------------------------------------------------------------
def stats(data):
    '''Return the average and standard error in data. '''
    ave = np.average(data)
    ave2 = np.average(data*data)
    err = np.sqrt(np.abs(ave2-ave**2)/(1.0*data.size-1.0) ) 

    return ave,err

# ----------------------------------------------------------------------
def line_counts(filename):
    '''Use wc to count the number of lines and header lines in a file. '''
    num_lines = int(subprocess.check_output(['wc', '-l', filename]).split()[0])
    num_header = str(subprocess.check_output(['grep','-o','-i','\#',filename])).count('#')
    return num_header,num_lines

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Generates averages from pimc output data.')
    parser.add_argument('-s', '--skip', 
                        help='How many input lines should we skip? [default: 0]')
    parser.add_argument('-l', '--header_lines', type=int, dest='header_lines', 
                        help='How many header lines to skip?')
    parser.add_argument('file', type=str, nargs='+', 
                        help='File or files to average.')
    parser.add_argument('-r','--repeated_header', 
                        help='deal with duplicate headers', action='store_true')

    args = parser.parse_args()

    fileNames = args.file

    if not args.skip:
        skip = 0
    else:
        if '.' in args.skip:
            skip = float(args.skip)
            if skip < 0.0 or skip >= 1.0:
                raise ValueError('skip < 0.0 or skip >= 1.0')
        else:
            skip = int(args.skip)

    for fileName in fileNames:

        # get the pimcid
        pimcid = fileName[-40:].rstrip('.dat')
        
        # We check to see if we are dealing with the one body density matrix
        if fileName.find('obdm') != -1:
            normalize = True

        try:
            # We peak into the file and determine how many header lines to skip
            header_lines = 0;
            if args.header_lines:
                header_lines = args.header_lines
            else:
                with open(fileName,'r') as inFile:
                    line = inFile.readline()
                    while line[0] == '#':
                        if ('PIMCID' in line) or ('ESTINF' in line):
                            header_lines += 1
                        line = inFile.readline()

            # open the file and determine how many measurements there are
            estData = np.genfromtxt(fileName,names=True,skip_header=header_lines, deletechars="")
            numLines = estData.size
            
            # Determine how many lines to skip if we have defined a fraction 
            if isinstance(skip, float):
                skip = int(numLines*skip)

            # If we have data, compute averages and error
            if numLines-skip > 0:
                print('# PIMCID %s' % pimcid)
                print('# Number Samples %6d' % (numLines-skip))
                for name in estData.dtype.names:
                    ave,err = stats(estData[name][skip:])

                    if err != 0.0:
                        rel_err = 100*np.abs(err/ave)
                    else:
                        rel_err = 0.0

                    if args.repeated_header:
                        cname = name.partition("_")[0]
                    else:
                        cname = name
                    print('%-16s%12.5f\t%12.5f\t%5.2f' % (cname,ave,err,rel_err))

        except:
          print(f"Couldn't Average File {fileName}")

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
