#!/usr/bin/python
# pimcbin.py
# Author:      2012-11-30 Chris Herdman
# Modified:    2017-07-17 Nathan Nichols
# =============================================================================
# 
# Reduce and do binning analysis on results for a single PIMC estimator 
# data file supplied as an input

import os,sys
import pyutils
import argparse
import numpy as np
import MCstat

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Reduce and do binning analysis on results for a single PIMC estimator file.')
    parser.add_argument('-s', '--skip', type=int, dest='skip', default = 0, help='How many input lines should we skip? [default: 0]')
    parser.add_argument('file', type=str, nargs='+',
                        help='File or files to average.')

    args = parser.parse_args()

    fileNames = args.file
    skip = args.skip

    for fileName in fileNames:
        normalize = False;

        # We check to see if we are dealing with the one body density matrix
        if fileName.find('obdm') != -1:
            normalize = True

        # We count the number of lines in the estimator file to make sure we 
        # have some data and grab the headers
        estFile = open(fileName,'r');
        estLines = estFile.readlines();
        numLines = len(estLines) - 2    # We expect two comment lines
        pimcid = estLines[0]
        print(pimcid)
        headers = estLines[1].split()
        estFile.close()

        # If we have data, compute averages and error
        if numLines-args.skip > 0:
            estData = pyutils.loadFile(fileName)

            # Now we skip data rows to test for convergence
            for n in range(args.skip):
                estData.pop(0)

            estAve = pyutils.average(estData,1)
            bins = MCstat.bin(np.array(estData))
            estErr = bins[-1,:]
            numData = len(estData)

            print pimcid, '# Number Samples %6d' %  numData
            if not normalize:
                for n,ave in enumerate(estAve):
                    if len(headers) - 1 ==  len(estAve):
                        label = headers[n+1]
                    else:
                        label = 'Col #%02d:' % n
                    print '%-16s%12.5f\t%12.5f' % (label,estAve[n],estErr[n])
            else:
                for n,ave in enumerate(estAve):
                    normAve = estAve[n]/estAve[0]
                    if abs(estAve[n]) > 1E-10:
                        normErr = (estErr[n] / estAve[n]) * normAve
                    else: 
                        normErr = 0.0;

                    if len(headers) - 1 ==  len(estAve):
                        label = headers[n+1]
                    else:
                        label = 'Col #%02d:' % n
                    print '%-16s%12.5f\t%12.5f' % (label,normAve,normErr)
    
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
