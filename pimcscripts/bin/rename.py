#!/usr/bin/python
# rename.py
# Adrian Del Maestro
# 09.28.2009
# 
# Rename a series of PIMC output files based on an input file

import os,sys,glob
from optparse import OptionParser

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

        # setup the command line parser options 
        parser = OptionParser() 

        # parse the command line options and get the file name
        (options, args) = parser.parse_args() 
        if len(args) != 1: 
                parser.error("need a file name")
        
        fileName = args[0]

        dirName = os.path.dirname(fileName)

        # The output file types
        fileType = ['estimator','log','obdm','pair','pcycle','state',
                    'super','worm','number','energy']

        # We parse the input file name
        fileParts = fileName.partition('log')

        # Now break up the data name into descriptive parts
        dataName = fileParts[2]
        dataName = dataName.rstrip('.dat')
        dataParts = dataName.split('-')

        # Get the ID number
        oldID = int(dataParts[-1])
        newID = oldID

        # Now we keep incrementing the ID number until we are sure it is unique
        while len(glob.glob(os.path.join(dirName,'*log*-%09d*' % newID))) > 0:
            newID += 1

        # Create the new data name
        dataName = ''
        for lab in dataParts[:-1]:
                dataName += lab + '-'
        dataName += '%09d.dat' % newID

        for ftype in fileType:
                oldName = fileParts[0] + ftype + fileParts[2]
                newName = fileParts[0] + ftype + dataName
                if os.path.exists(oldName):
                    os.popen('mv %s %s' % (oldName,newName))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
        main()
