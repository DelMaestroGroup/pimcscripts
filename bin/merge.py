#!/usr/bin/env python

#  Adrian Del Maestro
# 10.19.2009

'''merge.py

Description:
Merge the results of parallel PIMC output files, (same parameters, different
seeds) and potentially move the originals to an archive

Usage: merge.py [options] [<base_dir>]

Options:
  -h, --help                        Show this help message and exit
  -T <T>, --temperature=<T>         Simulation temperature in Kelvin
  -N <N>, --number_particles=<N>    Number of particles
  -n <n>, --density=<n>             Number density in Angstroms^{-d}
  -t <t>, --imaginary_time_step=<t> Imaginary time step
  -P <P>, --number_time_slices=<P>  Number of time slices
  -u <u>, --chemical_potential=<u>  Chemical potential in Kelvin
  -V <V>, --volume=<V>              Volume in Angstroms^d
  -L <L>, --Lz=<L>                  Length in Angstroms
  -s <skip>, --skip=<skip>          How many input lines should we skip?  [default: 0]
  -i <PIMCID>, --id=<PIMCID>        A list of PIMC ID numbers to include 
  -e <exclude> --exclude=<exclude>  A list of file types to exclude
  --seeds                           Create merge of average of individual seeds
  --canonical                       Are we in the canonical ensemble?
  --label=<label>                   A custom uuid label 
'''

from __future__ import print_function
from docopt import docopt
import os,sys
import numpy as np
import glob
import tarfile
import pimcscripts.pimchelp as pimchelp
import pimcscripts.MCstat as MCstat
import uuid
import re

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if data.ndim > dim:
        numBins  = data.shape[dim]
        dataAve  = np.average(data,dim) 
        try:
            bins = MCstat.bin(data) 
            dataErr = np.amax(bins,axis=0)
        except:
            dataAve2 = np.average(data*data,dim) 
            dataErr = np.sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 

    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr

# -----------------------------------------------------------------------------
def mergeData(pimc,etype,newID,skip,baseDir,idList=None,cyldir='',
              aveSeeds=False):
    ''' Merge the results of the PIMC data files to a single file. '''

    fileNames = pimc.getFileList(etype,idList=idList,cyldir=cyldir)
    if not fileNames:
        return

    # determine if we are trying to merge a cumulative estimator
    cumulative = etype in ['position','locsuper','planeavedensity','planeaveVext']

    # Determine if we are considering an off-diagonal estimator
    diagonalEst = not (etype == 'obdm' or etype == 'worm')

    # Open and prepare the first non-empty file
    n = 0
    empty = True
    while empty:
        with open(fileNames[n], 'r') as inFile:

            numLines = sum(1 for line in inFile) - 2
             
            if not numLines:
                n += 1
                empty = True
            else:
                empty = False

    with open(fileNames[n], 'r') as inFile:

        inLines = inFile.readlines()

        # Get the new header string
        if '#' in inLines[0]:
            header = inLines[0][2:].replace(pimc.id[n],newID)

            if aveSeeds:
                headerStats = header + f"{'Number Bins':>14}" 

            if inLines[2][0] == '#':
                header += inLines[1][2:]
                header += inLines[2][2:-1]
            else:
                header += inLines[1][2:-1]

        # get the data from the first file
        if isinstance(skip,int):
            skiprows=(1-cumulative)*(skip+2)*diagonalEst
        else:
            skiprows=(1-cumulative)*(int(skip*numLines)+2)*diagonalEst

        cdata = np.loadtxt(fileNames[n],ndmin=2,comments='#',skiprows=skiprows)

        # We either perform an average for each seed, or dump all data together
        # in one big file.
        if aveSeeds and not cumulative:
            numBins = [cdata.shape[0]]
            data = [np.average(cdata,0)]
        else:
            data = [cdata]

    # go through all other files and append data
    for i,fname in enumerate(fileNames[n+1:]):

        # Does the file exist?
        if len(glob.glob(fname)) > 0:

            # Does it contain any data?
            with open(fname, 'r') as inFile:
                numLines = sum(1 for line in inFile) - 2

            # if yes, figure out if we are skipping any rows
            if numLines:
                if isinstance(skip,int):
                    skiprows=(1-cumulative)*(skip+2)*diagonalEst
                else:
                    skiprows=(1-cumulative)*(int(skip*numLines)+2)*diagonalEst

                # load the data
                cdata = np.loadtxt(fname,ndmin=2,skiprows=skiprows)

            # if we have data, append to the array
            if cdata.size:
                if aveSeeds and not cumulative:
                    numBins.append(cdata.shape[0])
                    data.append(np.average(cdata,0))
                else:
                    if cumulative:
                        data[0] += cdata
                    else:
                        data.append(cdata)

    # Get the name of the new output file
    outName = os.path.basename(fileNames[0]).replace(pimc.id[0],newID)
    print('%-80s' % outName,end="")

    # for cumulative estimators we average, for all others we stack
    if cumulative:
        data = np.array(data[0]/(1.0*len(fileNames)))
    else:
        data = np.vstack(data)
        if aveSeeds:
            numBins = np.array(numBins)

    # open the output file for writing
    np.savetxt(baseDir + 'MERGED/' + cyldir + outName, data, fmt='%16.8E', 
               header=header, delimiter='')

    # Do the same for the summary statistics
    if aveSeeds and not cumulative:
        np.savetxt(baseDir + 'MERGED/' + cyldir + outName.replace(etype,f'bins-{etype}'), 
                   numBins, fmt='%16d', header=headerStats, delimiter='')
    print('%10d' % data.shape[0])

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # parse the command line options
    args = docopt(__doc__)

    if args['--skip']:
        if '.' in args['--skip']:
            skip = float(args['--skip'])
        else:
            skip = int(args['--skip'])

    canonical = args['--canonical']
    pimcID = args['--id']
    baseDir = args['<base_dir>']
    if not baseDir:
        baseDir = './'
    mergeDir = baseDir + 'MERGED'
    exclude_estimators = args['--exclude']

    if not exclude_estimators:
        exclude_estimators = []

    # We check if we have a MERGED directory, if not create it
    if not os.path.isdir(mergeDir):
        os.mkdir(mergeDir)

    # Create a .donotbackup file
    os.system('touch %s/.donotbackup' % mergeDir)
    
    # check if we have any cylinder estimators
    cylinder = os.path.isdir(os.path.join(baseDir,'CYLINDER'))
    
    # Create any necessary CYLINDER directories
    if cylinder and not os.path.isdir(os.path.join(mergeDir,'CYLINDER')):
        os.mkdir(os.path.join(mergeDir,'CYLINDER'))
        os.system('touch %s/CYLINDER/.donotbackup' % mergeDir)

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(canonical)
    dataName = pimchelp.getFileString_doc(args,reduce=False)

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,canonical,baseDir=baseDir)
    pimc.getSimulationParameters(idList=pimcID)

    # get a new unique uuid
    newID = str(uuid.uuid4())

    # If we have a custom label, append
    if args['--label']:
        newID = newID[:-len(args['--label'])] + args['--label']

    # Merge all the output files
    print('Merged data files:')
    for ftype in pimc.dataType:
        if ftype not in exclude_estimators:
            mergeData(pimc, ftype, newID, skip, baseDir, idList=pimcID, 
                     aveSeeds=args['--seeds']) 

    # Merge cylinder output files
    if cylinder:
        for ftype in pimc.dataType:
            if ftype not in exclude_estimators:
                mergeData(pimc, ftype, newID, skip, baseDir, idList=pimcID, 
                          cyldir='CYLINDER/', aveSeeds=args['--seeds'])

    # Prepare the new log file
    oldLogName = pimc.getFileList('log', idList=pimcID)[0]
    newLogName = os.path.basename(oldLogName).replace(pimc.id[0],newID)

    # write the new log file to disk
    with open(oldLogName,'r') as oldLogFile:
        logData = oldLogFile.read().replace(pimc.id[0],newID)
        with open(mergeDir + os.path.sep + newLogName,'w') as newLogFile:
            newLogFile.write(logData)

    # Do the same if we are merging cylinder files
    if cylinder:
        with open(oldLogName,'r') as oldLogFile:
            logData = oldLogFile.read().replace(pimc.id[0],newID)
            with open(mergeDir + '/CYLINDER/' + newLogName,'w') as newLogFile:
                newLogFile.write(logData)

    # We first create the name of the output tar file
#   mergeName = pimc.getFileList('estimator')[0].rstrip('.dat').replace('estimator','merged')
#   mergeName += '-%09d' % pimc.id[-1]
#   mergeName += '.tar.gz'
#
#   # Archive all the output files  
#   tar = tarfile.open('MERGED/'+mergeName, "w:gz")
#   for etype in pimc.outType:
#       fileNames = pimc.getFileList(etype)
#
#       # We have to exclude the merged file
#       for i,fname in enumerate(fileNames):
#           if fname.find(str(newID)) != -1:
#               fileNames.pop(i)
#           
#       # Add all the output files to the tar archive
#       for fname in fileNames:
#           tar.add(fname)
#
#           # delete the file that we have now archived
#           os.system('rm %s' % fname)
#
#   tar.close()


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__": 
    main()

