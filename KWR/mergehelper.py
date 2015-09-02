import os
import glob
import sys

class PIMCHelp:
    ''' Helper methods for merging random number seeds. Note that this class assumes your 
    data files are sorted so that each directory contains only data files of identical
    simulation parameters but different process numbers'''

    # ----------------------------------------------------------------------
    def __init__(self,fileList,baseDir=''):
        
        # Upon instantiation, filter out unwanted files
        self.datfiles = self.removeFiles(fileList)
        # Determine whether we are in canonical or grand canonical ensemble
        self.prefix = self.getEnsemble()
        # Get the part of the file name containing simulation parameters.
        self.dataName = self.getDataName() 
        # The base directory for this particular instantiation of this class
        self.baseDir  = baseDir
        # Declare and fill the list of PIMCIDs
        self.IDs = []
        self.getIDs(self.datfiles)
        # Declare and fill the list of all the estimators 
        self.estimators = []
        self.getEstimators()
        # The complete set of parameters
        self.params = self.getSimulationParameters()

        self.outType  = ['estimator', 'number', 'obdm', 'pair', 'pcycle', 'super', 
                         'worm', 'radial', 'log', 'state','virial']
    
    # -----------------------------------------------------------------------------
    def removeFiles(self,fileList):
        '''Remove any files that are not data files'''
        return filter(lambda filename: (filename[-4:] == '.dat') and (filename[0] != '.'), fileList)
        
    # -----------------------------------------------------------------------------
    def getEnsemble(self):
        '''Determine the ensemble for this list of dat files'''
        if all(datFile[0:3] == 'gce' for datFile in self.datfiles):
            # We are in grand canonical
            return 'gce'
        elif all(datFile[0:2] == 'ce' for datFile in self.datfiles):
            # We are in canonical
            return 'ce'
        else:
            # We have mixed ensembles. Terminate immediately
            print 'The list of data files passed to PIMCHelp has mixed ensembles, terminating immediately'
            sys.exit(1)                
                
    # -----------------------------------------------------------------------------
    def getDataName(self):
        '''Get the part of a dat file name that has the simulation params'''
        estimator = getEstimator(self.datfiles[0])
        dataName = self.datfiles[0][len(self.prefix)+len(estimator)+2:-13]
        # If every file has the same data name, return the data name
        if all(datfile[len(self.prefix)+len(getEstimator(datfile))+2:-13] == dataName for datfile in self.datfiles):
            dataName = dataName.rstrip('-')
            return dataName
        # Else we have mixed simulation parameters, terminate immediately
        else: 
            print "The list of dat files passed to PIMCHelp has mixed simulation parameters, terminating immediately"
            sys.exit()
            
    # -----------------------------------------------------------------------------
    def getIDs(self,datFiles):
        '''Get all the PIMCIDs for a set of dat files'''
        for datfile in datFiles:   
            ID = int(datfile[-13:-4])
            if ID not in self.IDs:
                self.IDs.append(ID)    
                  
    # -----------------------------------------------------------------------------
    def getID(self,fileName): 
        ''' Return the ID number corresponding to a given filename. '''
        ID = int(fileName[-13:-4])
        return ID
        
    # -----------------------------------------------------------------------------        
    def getEstimators(self):
        '''Get all the estimators found in a given list of dat files'''
        for datfile in self.datfiles:
            parts = datfile.split('-')
            est = parts[1]
            if est not in self.estimators:
                self.estimators.append(est)
        
    # ----------------------------------------------------------------------
    def getSimulationParameters(self): 
        '''Creates a dictionary of parameter keys and their corresponding values assuming that
        the simulation parameters for every file in the directory are identical'''

        # Get the values of all simulation parameters
        logfiles=filterfilesbytype("log",self.datfiles)
        logName = logfiles[0]
        params = False
        paramsMap = {}
        with open(self.baseDir+logName, 'r') as logFile:
            for line in logFile:
                if 'Begin Simulation Parameters' in line:
                    params = True
                elif 'End Simulation Parameters' in line:
                    break

                if params and ':' in line:
                    keyVal = line.split(':')
                    paramsMap[keyVal[0].strip()] = keyVal[1].strip()

        # Add an element to the parameter map for the linear dimension (Lz) of
        # the container
        paramsMap['Container Length'] = paramsMap['Container Dimensions'].split('x')[-1]

        return paramsMap
        
# -----------------------------------------------------------------------------
def getHeadersFromFile(fileName, skipLines=0): 
    ''' Get the data column headers from a PIMC output file. '''

    with open(fileName,'r') as inFile:
        inLines = inFile.readlines();
        n = skipLines
        if 'PIMCID' in inLines[n]:
            headers = inLines[n+1].split()
        else:
            headers = inLines[n].split()
        headers.pop(0)

    return headers

# -----------------------------------------------------------------------------
def getEstimatorType(estimator):
    '''Retrieves the estimator type for a given estimator'''
    if estimator in ['estimator']:
        estType = 'Scalar'
    elif estimator in ['planedensity']:   
        estType = 'Vector'
    elif estimator in ['position','locsuper']:
        estType = 'Cumulative'
    elif estimator in ['state', 'log','worm']:
        estType = 'DoNotMerge'
    else:
        estType = None
    return estType
        
# -----------------------------------------------------------------------------        
def filterfilesbytype(filetype,datfiles):
    '''Pull out only the files of a particular type from a list of dat files'''
    return filter(lambda string: filetype in string, datfiles)       
    
# -----------------------------------------------------------------------------      
def getEstimator(datfile):
    '''Get the estimator for a given dat file'''
    parts = datfile.split('-')
    return parts[1]
               
# -------------------------------------------------------------------------------
def checkEnsemble(canonical):
    '''Make sure the correct ensemble flag is specified when processing data. Immediately
    exits if a conflict is detected.'''

    gceFiles = (len(glob.glob('gce-*'))>0)
    ceFiles = (len(glob.glob('ce-*'))>0)

    if (ceFiles and not gceFiles) and not canonical:
        sys.exit('Need to include --canonical for the canonical ensemble!')       
        
        
        
        
