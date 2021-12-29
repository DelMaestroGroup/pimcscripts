#!/usr/bin/env python

# merge.py
# Adrian Del Maestro
# 10.19.2009
# 
# Merge the results of parallel PIMC output files, and move the 
# originals to an archive

import os,sys,glob,shutil
import tarfile
import pimcscripts.pimchelp as pimchelp
from optparse import OptionParser
from pylab import *


def ReadSkipping(fileName,skip=0):
    '''Read measurements from fileName while skipping first ${skip} of them'''
    
    # if fileName doesnt exist, exit
    if not(len(glob.glob(fileName)) > 0):
       return 0

    # Determine if we are considering an off-diagonal estimator
    diagonalEst = ((fileName.find('obdm') == -1) and (fileName.find('worm') == -1))   

    # Take the content from the file
    inFile = open(fileName,'r');
    inLines = inFile.readlines();
    inFile.close()

    # Strip any comment lines
    if len(inLines) != 0:
       inLines.pop(0)
       inLines.pop(0)

    #Number of measurements after skipping
    numLines = len(inLines)-diagonalEst*skip

    if diagonalEst:
       #a check in case we want to skip more lines than there are
       if numLines > 0:
          # Skip ${skip} measurements                
          for i in range(skip):
              inLines.pop(0)   
       else:
            inLines =[]
   
    return inLines, numLines      

def CreateFileSkipping(FromfileName,oldPIMCID,newPIMCID,skip=0):
    '''Create a new file whose oldPIMCID is replaced with newPIMCID and
      skip ${skip} first measurements'''

    #Read and adapt the header
    inFile = open(FromfileName,'r');
    firstLine  = inFile.readline().replace(str(oldPIMCID),str(newPIMCID));
    secondLine = inFile.readline()
    inFile.close();

    #Read the measurements we want to keep
    inList,numLines = ReadSkipping(FromfileName,skip)
    if inList == 0:
       return 0
     
    # get the output file name and open the file for writing
    outName = FromfileName.replace(str(oldPIMCID),str(newPIMCID))
    outFile = open('MERGED/' + outName,'w');
    print '%-80s' % outName,
 
    #write the header
    outFile.write(firstLine)
    outFile.write(secondLine)
    
    #if there is more measurements than we want to skip
    if numLines > 0:
       #write the content
       outFile.writelines(inList) 

    return outFile,numLines


def addFromToSkipping(FromfileName,toFile,skip=0):
    '''add PIMC measurements from FromfileName to toFile, while skipping skip first measurements'''
     
    #Read the measurements we want to keep
    inList,numLines = ReadSkipping(FromfileName,skip)
    if inList == 0:
       return -1

    #Write them
    if len(inList) != 0:
       toFile.writelines(inList)
    
    return numLines



# -----------------------------------------------------------------------------
def mergeData(pimc,type,newID,skip=0,isRestarted=False,pimcIds=None):
    ''' Merge the results of the PIMC data files to a single file. '''

    #file-names to merge
    fileNames = pimc.getFileList(type,pimcIds)
    #Total number of measurements in the merged file
    numLines = 0
    #Have we created an output file?
    fileExist = False
  
    for i,fname in enumerate(fileNames):
        # Does the file exist?
        if len(glob.glob(fname)) > 0:
           #if we havent yet the output file, create one
           if not fileExist:
               outFile,numLines = CreateFileSkipping(fname,pimcIds[0],newID,skip) 
               fileExist = (outFile != 0) 
               #taken it is a restarted job, if there is less measurements than we want to skip in 
               #the first file we will need to skip more in the files that follow
               if (isRestarted) and (numLines < 0):
                  skip = -numLines   
               #if it is not a restarted job, but the file contains less than enough measurements
               if (not isRestarted) and (numLines < 0):
                  numLines = 0 
           else:
               #if we are merging files that got started from a parent state file
               #we skip measurements only from the parent  
               if isRestarted:
                   #if the parent has more measurements that we want to skip
                   if not(numLines < 0):
                       numLines += addFromToSkipping(fname,outFile)  
                   #if not, we will have to skip some measurements in children
                   else:
                       numLines = addFromToSkipping(fname,outFile,skip)   
                       skip = -numLines     
               else:
                   numLines += addFromToSkipping(fname,outFile,skip)   
                   if numLines < 0:
                      numLines = 0         
    outFile.close() 
    print '%10d' %numLines

    # Now we check if a CYLINDER folder is present, if so, we repeat the process
    if len(glob.glob('CYLINDER')) > 0:
        numLines = 0
        fileExist = False
        # We check if we have a CYLINDER output directory, if not create it
        if len(glob.glob('MERGED/CYLINDER')) == 0:
              os.system('mkdir MERGED/CYLINDER')
        for i,fname in enumerate(fileNames):
            cylfname = 'CYLINDER/' + fname
            if len(glob.glob(cylfname)) > 0:
               if not fileExist:
                  outFile,numLines = CreateFileSkipping(cylfname,pimcIds[0],newID,skip) 
                  fileExist = (outFile != 0) 
                  if (isRestarted) and (numLines < 0):
                     skip = -numLines   
                  if (not isRestarted) and (numLines < 0):
                     numLines = 0               
               else:
                   if isRestarted:
                      if not(numLines < 0):
                          numLines += addFromToSkipping(cylfname,outFile) 
                      else:
                           numLines = addFromToSkipping(cylfname,outFile,skip)   
                           skip = -numLines     
                   else:
                        numLines += addFromToSkipping(cylfname,outFile,skip)   
                        if numLines < 0:
                           numLines = 0      
        if fileExist:
           outFile.close() 
           print '%10d' %numLines
    return

#--------------------------------------------------------------------------------------------------------------------

def getNewPIMCID(pimc):
    ''' Find a new PIMCID which is the average of the ones to merge, and make sure it doesn't already exist'''
    newID = 0
    for id in pimc.id:
        newID += int(id)
    newID = int(newID/(1.0*len(pimc.id)))

    # Now we keep incrementing the ID number until we are sure it is unique
    while ( (len(glob.glob('*log*-%09d*' % newID)) > 0) or
            (len(glob.glob('MERGED/*log*-%09d*' % newID)) > 0) or 
            (len(glob.glob('../*log*-%09d*' % newID)) > 0) ):
        newID += 1
    return newID

#--------------------------------------------------------------------------------------------------------------------

def getMergeSets(pimc,VarP):
    '''Create a dictionnary of sets of pimcids to be merged'''
    parMap = {'n':'Initial Density', 'T':'Temperature', 'N':'Initial Number Particles',
              't':'Imaginary Time Step', 'u':'Chemical Potential', 'V':'Container Volume',
              'L':'Container Length','r':'Cylindrical Container Radius'}
    
    #assign a complete name for the varying paramater
    VarP = parMap[VarP]
    
    #dictionnary key   = a value of the varying paramater
    #dicitonnary value = a list of pimcids with this value
    MergeSets= {}
    
    for PIMCid,param in pimc.params.items():
        if MergeSets.has_key(param[VarP]):
            MergeSets[param[VarP]].append(PIMCid)  
        else:
            MergeSets[param[VarP]] = [PIMCid] 
    
    #sort each subset
    for sets in MergeSets.iterkeys():
        MergeSets[sets].sort()     

    return MergeSets  

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",
                      help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",
                      help="number density in Angstroms^{-d}")
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",
                      help="imaginary time step")
    parser.add_option("-M", "--number-time-slices", dest="M", type="int",
                      help="number of time slices")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",
                      help="chemical potential in Kelvin") 
    parser.add_option("-V", "--volume", dest="V", type="float",
                      help="volume in Angstroms^d") 
    parser.add_option("-L", "--Lz", dest="L", type="float",
                      help="Length in Angstroms") 
    parser.add_option("-R", "--radius", dest="R", type="float",
                      help="radius in Angstroms") 
    parser.add_option("-v", "--varp", dest="varp",
                      choices=['T','N','n','u','t','L','V','r'], 
                      help="varying parameter, one of [T,N,n,u,t,L,V,r]") 
    parser.add_option("--canonical", action="store_true", dest="canonical",
                      help="are we in the canonical ensemble?")
    parser.add_option("--restarted", action="store_true", dest="restarted",
                      help="are we merging pimcs that got restarted from a parent state?")
    parser.add_option("-s", "--skip", dest="skip", type="int",
                      help="how many input lines should we skip?")
    parser.set_defaults(skip=0)

    parser.set_defaults(canonical=False)
    parser.set_defaults(restarted=False)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 
    if len(args) > 0: 
        parser.error("incorrect number of arguments")

    # We check if we have a MERGED directory, if not create it
    if os.path.exists('MERGED/CYLINDER') == False: 
       os.makedirs('MERGED/CYLINDER')
    
    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.canonical)

    # Form a pattern of pimc output filenames' structure satysfying our options
    dataName = pimchelp.getFileString(options,reduce=False)

    # Create the PIMC analysis helper
    pimc = pimchelp.PimcHelp(dataName,options.canonical)
 
    # Fill up the simulation parameters maps
    pimc.getSimulationParameters()

    # Delete those pimcIDs that do not satysfy parameters 
    # that are not contained in the pimc output filanames' structure
    #"implicit" parameters
    pimc.ApplyImplicitParameters()

    #if there is not need to merge with a varying parameter
    if (not options.varp):
        #Create new pimcID
        newID = getNewPIMCID(pimc)
    
        # Merge all the output files
        print 'Merged data files:'
        for type in pimc.dataType:
            mergeData(pimc,type,newID,options.skip,options.restarted)
    
        # copy over the log file
        oldLogName = pimc.getFileList('log')[0]
        newLogName = oldLogName.replace(str(pimc.id[0]),str(newID))
        os.system('cp %s %s' % (oldLogName,'MERGED/'+newLogName))

     #with a varying parameter, one needs to group corresponding pimcIds
    else:
         #group pimcIds with the same varying parameter  
         MergeSets = getMergeSets(pimc,options.varp)
         for varp in sorted(MergeSets.iterkeys()):
             mergeSet = MergeSets[varp]
             print '\nMerged data files for %s=%s:\n' %(options.varp,varp)
             print 'PIMCids to merge: %s' %mergeSet  


             #if there is only one pimcId with a varp, the just copy the files
             if (len(mergeSet) == 1):
                lsCommand = 'ls *log*%s*' %mergeSet[0]
                LogName = os.popen(lsCommand).read().split('\n')[0]
                shutil.copyfile(LogName,'MERGED/'+LogName) 
                
                for type in pimc.dataType: 
                    lsCommand = "ls *%s*%s*" %(type,mergeSet[0]) 
                    fileName = os.popen(lsCommand).read().split('\n')[0]
                    outFile,numLines = CreateFileSkipping(fileName,mergeSet[0],mergeSet[0],options.skip)
                    print '%10d' %numLines    
                    outFile.close
    
                lsCommand = "ls CYLINDER/*%s*" %mergeSet[0]
                fileNames = os.popen(lsCommand).read().split('\n')
                fileNames.pop()
                for files in fileNames:
                    outFile,numLines = CreateFileSkipping(files,mergeSet[0],mergeSet[0],options.skip)
                    print '%10d' %numLines    
                    outFile.close
             #otherwise we need to be careful what files do we merge together
             else:
                 #Create new pimcID
                 newID = getNewPIMCID(pimc)      
                 for type in pimc.dataType:
                     mergeData(pimc,type,newID,options.skip,options.restarted,mergeSet)
                 lsCommand = 'ls *log*%s*' %mergeSet[0]
                 oldLogName = os.popen(lsCommand).read().split('\n')[0]
                 newLogName = oldLogName.replace(str(mergeSet[0]),str(newID))
                 shutil.copyfile(oldLogName,'MERGED/'+newLogName)
                 #os.system('cp %s %s' % (oldLogName,'MERGED/'+newLogName))    




    # We first create the name of the output tar file
#   mergeName = pimc.getFileList('estimator')[0].rstrip('.dat').replace('estimator','merged')
#   mergeName += '-%09d' % pimc.id[-1]
#   mergeName += '.tar.gz'
#
#   # Archive all the output files  
#   tar = tarfile.open('MERGED/'+mergeName, "w:gz")
#   for type in pimc.outType:
#       fileNames = pimc.getFileList(type)
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

