# Adrian Del Maestro
# 10.19.2009

'''merge.py

Description:
Merge the results of parallel PIMC output files, (same parameters, different
seeds) and potentially move the originals to an archive

Usage: merge.py [options]

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
  -s <skip>, --skip=<skip>          How many input lines should we skip? [default: 0]
  -i <PIMCID>, --id=<PIMCID>        A list of PIMC ID numbers to include 
  --canonical                       Are we in the canonical ensemble?
  -D <dir>, --working_directory=<dir>   The directory where we perform the merge.  [default: ]
'''

from __future__ import print_function
from optparse import OptionParser
from docopt import docopt
import os
import numpy as np
import glob
import tarfile
import pimchelp

# -----------------------------------------------------------------------------
def mergeData(pimc,etype,newID,skip,baseDir,idList=None,cyldir=''):
    ''' Merge the results of the PIMC data files to a single file. '''

    fileNames = pimc.getFileList(etype,idList=idList,cyldir=cyldir)
    if not fileNames:
        return
    
    # determine if we are trying to merge a cumulative estimator
    cumulative= etype in ['position', 'locsuper']

    # Determine if we are considering an off-diagonal estimator
    diagonalEst = not (etype == 'obdm' or etype == 'worm')

    # we skip comment rows and any data rows
    skiprows = (skip + 2)*diagonalEst

    # Open and prepare the new file
    with open(fileNames[0], 'r') as inFile:
        inLines = inFile.readlines();

        # Get the new header string
        if '#' in inLines[0]:
            header = inLines[0][2:].replace(str(pimc.id[0]),str(newID))
            header += inLines[1][2:-1]

        # get the data from the first file
        data = np.loadtxt(fileNames[0],ndmin=2,comments='#',skiprows=skiprows)

    # go through all other files and append data
    for i,fname in enumerate(fileNames[1:]):

        # Does the file exist?
        if len(glob.glob(fname)) > 0:

            cdata = np.loadtxt(fname,ndmin=2,skiprows=skiprows)
            # if we have data, append to the array
            if cdata.size:
                data = np.vstack((data,cdata))

    # Get the name of the new output file
    outName = os.path.basename(fileNames[0]).replace(str(pimc.id[0]),str(newID))
    print('%-80s' % outName,end="")

    # for cumulative estimators we average
    if cumulative:
        data = np.average(data,axis=0)

    # open the output file for writing
    np.savetxt(baseDir + 'MERGED/' + cyldir + outName, data, fmt='%16.8E', 
               header=header, delimiter='')
    print('%10d' % data.shape[0])


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # parse the command line options
    args = docopt(__doc__)

    skip = int(args['--skip'])
    canonical = args['--canonical']
    pimcID = args['--id']
    baseDir = args['--working_directory']
    mergeDir = baseDir + 'MERGED'

    # We check if we have a MERGED directory, if not create it
    if not os.path.isdir(mergeDir):
        os.mkdir(mergeDir)

    # Create a .donotbackup file
    os.system('touch %s/.donotbackup' % mergeDir)
    
    # check if we have any cylinder estimators
    cylinder = os.path.isdir(mergeDir + 'CYLINDER')
    
    # Create any necessary CYLINDER directories
    if cylinder and not os.path.isdir(mergeDir + '/CYLINDER'):
        os.mkdir(mergeDir + '/CYLINDER')
        os.system('touch %s/CYLINDER/.donotbackup' % mergeDir)

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(canonical)
    dataName = pimchelp.getFileString_doc(args,reduce=False)

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,canonical,baseDir=baseDir)
    pimc.getSimulationParameters(idList=pimcID)

    # We try to find a new PIMCID which is the average of the ones to merge, and
    # make sure it doesn't already exist
    newID = 0
    for pid in pimc.id:
        newID += int(pid)
    newID = int(newID/(1.0*len(pimc.id)))

    # Now we keep incrementing the ID number until we are sure it is unique
    while ((len(glob.glob(baseDir + '*estimator*-%09d*' % newID)) > 0) or
           (len(glob.glob(baseDir + 'MERGED/*estimator*-%09d*' % newID)) > 0)):
        newID += 1

    # Merge all the output files
    print('Merged data files:')
    for ftype in pimc.dataType:
        mergeData(pimc, ftype, newID, skip, baseDir, idList=pimcID) 

    # Merge cylinder output files
    if cylinder:
        mergeData(pimc, ftype, newID, skip, baseDir, idList=pimcID, 
                  cyldir='CYLINDER/')

    # copy over the log file
    oldLogName = pimc.getFileList('log', idList=pimcID)[0]
    newLogName = os.path.basename(oldLogName).replace(str(pimc.id[0]),str(newID))
    os.system('cp %s %s' % (oldLogName, mergeDir + os.path.sep + newLogName))

    # Do the same if we are merging cylinder files
    if cylinder:
        print("CYLINDER")
        os.system('cp %s %s' % (oldLogName, mergeDir + '/CYLINDER/' + newLogName))

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

