#!/usr/bin/env python

"""Fix Tiny Numbers
Author: Adrian Del Maestro
Date: 07.08.2013

Description:
  Fix an overflow for extremly small (<10-100) numbers in estimator files.
  
Usage:
  fix_tiny.py [--commit] <estimator-file>...

  fix_tiny.py -h | --help 

Options:
  --commit      Commit the changes to disk by replacing the file
  -h --help     Show this screen.
"""

from docopt import docopt
import os,sys
import re

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # read in teh command line arguments
    args = docopt(__doc__)

    fileNames =  args['<estimator-file>']

    if len(fileNames) < 1:
        print ("Need to specify at least one state file.")
        sys.exit()

    for fileName in fileNames:
        oldEstFile = file(fileName,'r')
        oldEstText = oldEstFile.read()
        oldEstFile.close()

        # now we test for any possible numbers in scientific notation jammed up
        # against each other
        def splitScientificNumber(matchObj):
            firstNum = matchObj.group('firstNum')
            secondNum = matchObj.group('secondNum')
            if secondNum[0] == '.':
                secondNum = firstNum[-1] + secondNum
                firstNum = firstNum[:-1]
            return  firstNum + ' ' + secondNum

        correctedEstText = re.sub('(?P<firstNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)(?P<secondNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)',splitScientificNumber,oldEstText)
        correctedEstText = re.sub('(?P<firstNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)(?P<secondNum>[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9][0-9])?)',splitScientificNumber,correctedEstText)

        newEst = correctedEstText.rstrip()

        
        # save the new estimator
        if args['--commit']:
            oldEstFile = file(fileName,'w')
            oldEstFile.write(newEst)
            oldEstFile.close()
        else:
            # output the new estimator
            print(newEst)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
