# =============================================================================
# Functions used in pull/push seedDirs script.
#
# This script assumes you will not have more than 1000 different
# random number seeds.
#
# Author:           Max Graves
# Last Revised:     24-OCT-2013
# =============================================================================

import os,argparse,re,sys,glob,shutil,subprocess
import paramiko
import getpass
import pylab as pl
import numpy as np

def parseCMD():
    ''' parse the command line. '''
    # set up parser/ subparsers
    parser = argparse.ArgumentParser(description='Cluster Management Tools 1.2')
    subparsers = parser.add_subparsers(help='enter one of these before --help')
    pullParse = subparsers.add_parser('pull', help='Pull Options')
    pushParse = subparsers.add_parser('push', help='Push Options')
    
    # common options
    parser.add_argument('-U', '--UserName', type=str,
            default='mtgraves',
            help='Username for your VACC account.')
    parser.add_argument('-t', '--targetDir', type=str,
            help='Push (Pull) directories to (from) ~/(--targetDir)')

    # push options
    pushParse.add_argument('-s', '--submitJobs', action='store_true', 
            dest='submitJobs',default=False,
            help='Do you want to automatically submit the jobs?')
    pushParse.add_argument('-L', '--lowSeed', type=int,
            default=0,
            help='Low seed.')
    pushParse.add_argument('-H', '--highSeed', type=int,
            default=0,
            help='High seed. (highSeed-lowSeed = numSeeds)')
    pushParse.add_argument('-B', '--N52BulkHe', action='store_true', 
            dest='N52BulkHe',default=False,
            help='Access the submit script for N=52, 3D Bulk He-4')   
    
    # pull options
    pullParse.add_argument('-d', '--deleteDirs', action='store_true', 
            dest='delDir',default=True,
            help='Do you want to delete the individual seed direcs?')
    pullParse.add_argument('-c', '--crunch', action='store_true',
            dest='Crunch',default=False,
            help='Combine all arrays of same temperature after pulling files.')
    pullParse.add_argument('-r', '--reduceType', type=str,
            default='T',
            help='Variable to reduce over [T,u]')
    pullParse.add_argument('-T', '--trim', action='store_true',
            dest='trimData', default=True,
            help='Create new reduced (equal length column) data files.')

    return parser.parse_args()

def Credentials(userN):
    ''' get username and password '''
    print 'username:',userN
    passwd = getpass.getpass(prompt='password: ')
    return passwd

def trimData(fileNames):
    '''
    Trims data files so that the length of each column is equal
    to that of the shortest column.  This is so that loadtxt will
    work properly.  Note that this is rather naughty.
    '''
    for fileName in fileNames:
        minLen = 100000000
        numChanged = 0
        lineStop = 0
        # count number of lines that have full data 
        with open(fileName) as inFile:
            for nLine, line in enumerate(inFile):
                if line[0] != '#':
                    if int(len(line)) < minLen:
                        minLen = int(len(line))
                        numChanged += 1
                        if numChanged == 2:
                            lineStop = nLine

        with open(fileName) as inFile, open('trimmed'+fileName[7:], 'w') as outFile:
            print 'Only writing ',lineStop, ' lines.'
            for n, line in enumerate(inFile):
                if n < lineStop:
                    outFile.write( line )

def returnSeedDirName(s):
    '''
    Takes seedn (or n), returns seed00n. Takes seednn, returns seed0nn.
    '''
    if type(s) == int:
        seedNum = s
    elif type(s) == str:
        seedName = re.match(r"([a-z]+)([0-9]+)", s, re.I)
        seedNum = seedName.groups()[1]
    if len(str(seedNum))==1:
        newName = 'seed00'+str(seedNum)
    elif len(str(seedNum))==2:
        newName = 'seed0'+str(seedNum)
    elif len(str(seedNum))==3:
        newName = 'seed'+str(seedNum)
    else:
        print 'this script doesnt handle >1000 seeds'
    return newName

def repeatCheck():
    ''' check for repeated pimcIDs '''
    checkForRepeats = glob.glob('*log*')
    pimcIDs = []
    for f in checkForRepeats:
        pimcIDs.append(f[-13:-4])

    if len(pimcIDs) == len(set(pimcIDs)):
        print 'All systems are a go!'
    else:
        sys.exit('Repeats of pimcIDs exist!')

def renameFilesInDirecs(delDir):
    '''
    Rename all files within their directories and
    optionally delete all directories pulled from cluster.
    Unzip all files that are in .tar.gz format.
    '''
    dirNames = glob.glob('*seed*')
    for dirName in dirNames:
        os.chdir('./'+dirName)
        filesIn = glob.glob('*.dat*')
        for zebra in filesIn:
            if zebra[-3:]=='.gz':
                comm = ("gunzip", str(zebra))
                subprocess.check_call(comm)
                shutil.copy2(zebra[:-3], '../'+zebra[:-16]+dirName[-3:]+zebra[-13:-3])
                #print zebra[:-16]+dirName[-3:]+zebra[-13:-3]
                #print zebra[:-3]
            else:
                shutil.copy2(zebra, '../'+zebra[:-13]+dirName[-3:]+zebra[-10:])
        os.chdir('..')
        if delDir:
            shutil.rmtree('./'+dirName)
    if delDir:
        print 'Chose to delete original direcs and rename files with seed num.'

def checkIfEmpty(fName,n):
    '''
    takes the first non-header line number and returns true or false
    depending upon whether that line is blank or not.
    '''
    Empty = False
    fp = open(fName)
    numLines = 0
    for line in fp:
        numLines += 1
    fp.close()

    if numLines == n:
        Empty=True
    
    return Empty

def makeTempList(estimFiles, canonical):
    ''' make list of all temperatures, in order, based on ensemble '''
    tempList = pl.array([])
    for f in estimFiles:
        if canonical: 
            temptemp = f[13:19]
        else:
            temptemp = f[14:20]
        
        if temptemp not in tempList:
            tempList = pl.append(tempList, temptemp)
    
    return pl.sort(tempList)

def makeMuList(estimFiles, canonical):
    ''' make list of all temperatures, in order, based on ensemble '''
    muList = pl.array([])
    for f in estimFiles:
        if canonical: 
            tempMu = f[28:35]
        else:
            tempMu = f[29:36]
        
        if tempMu not in muList:
            muList = pl.append(muList, tempMu)
    
    return pl.sort(muList)

def crunchData():
    '''
    This function will grab all (g)ce-estimator files in a directory and
    if there are multiple data files for the same temperature it will write
    all of the data into a new file categorized by the temperatures.
    The actual data stored will be in columns of the three terms (in order)
    needed to compute C_v. ( E, EEcv*beta^2, Ecv*beta, dEdB*beta^2 ).
    '''

    args = parseCMD()
    reduceType = args.reduceType
    # are we reducing over Temps or Mus?
    TempRed, MuRed = False, False
    if reduceType == 'T':
        TempRed = True
    elif reduceType == 'u':
        MuRed   = True

    # glob for files
    estimFiles  = glob.glob('*estimator*')
    biPartFiles = glob.glob('*bipart_dens*')
    superFiles  = glob.glob('*super*')
    ntWindFiles = glob.glob('*ntWind*')

    # check for bipartition files
    if biPartFiles == []:
        biPart = False
    else:
        biPart = True

    # check for ntWind files
    if ntWindFiles == []:
        ntwFiles = False
    else:
        ntwFiles = True

    # check which ensemble
    canonical = True
    if estimFiles[0][0] == 'g':
        canonical = False

    # make list of all reduced variable, in order, based on ensemble
    if TempRed:
        tempList = makeTempList(estimFiles, canonical)
    elif MuRed:
        tempList = makeMuList(estimFiles, canonical)

    # lists to hold all data
    allTempsE   = []
    allTemps1   = []
    allTemps2   = []
    allTemps3   = []
    if biPart:
        allTempsfD = []
        allTempsbD = []
    allTempsSuper = []
    allTempsWx2 = []
    allTempsWy2 = []
    allTempsWz2 = []

    allTempsNTW = []

    for numTemp,temp in enumerate(tempList):
        # only grab estimator files of correct temperature
        if TempRed:
            estFiles = glob.glob('*estimator-%s*' % temp)
            superFiles = glob.glob('*super-%s*' % temp)
        elif MuRed:
            estFiles = glob.glob('*estimator*%s*' % temp)
            superFiles = glob.glob('*super*%s*' % temp)
        
        E       = pl.array([])
        EEcv    = pl.array([])
        Ecv     = pl.array([])
        dEdB    = pl.array([])
        Super   = pl.array([])
        Wx2     = pl.array([])
        Wy2     = pl.array([])
        Wz2     = pl.array([])

        NTW     = pl.array([])

        if biPart:
            if TempRed:
                bipFiles = glob.glob('*bipart_dens-%s*' % temp)
            elif MuRed:
                bipFiles = glob.glob('*bipart_dens*%s*' % temp)
            bulkDens = pl.array([])
            filmDens = pl.array([])
            for bFile in bipFiles:
                if checkIfEmpty(bFile,2):
                    pass
                else:
                    fD, bD = pl.genfromtxt(bFile,unpack=True, \
                            usecols=(0,1))
                    bulkDens = pl.append(bulkDens, bD)
                    filmDens = pl.append(filmDens, fD)
            allTempsfD += [[filmDens]]
            allTempsbD += [[bulkDens]]
        
        for tFile in estFiles:
            #print tFile
            if checkIfEmpty(tFile,2):    # check for empty file
                pass
            else:
                ET, EEcvT, EcvT, dEdBT = pl.genfromtxt(tFile, unpack=True, \
                        usecols=(4,11,12,13))
                E       = pl.append(E, ET)
                EEcv    = pl.append(EEcv, EEcvT)
                Ecv     = pl.append(Ecv, EcvT)
                dEdB    = pl.append(dEdB, dEdBT)
        
        if TempRed:
            print 'T=',temp,', bins=',len(E)
        elif MuRed:
            print 'mu=',temp,', bins=',len(E)

        for sFile in superFiles:
            if checkIfEmpty(sFile,2):
                pass
            else:
                rhos_rho, wx2, wy2, wz2 = pl.genfromtxt(sFile, unpack=True, usecols=(0,1,2,3))
                Super = pl.append(Super, rhos_rho)
                Wx2 = pl.append(Wx2, wx2)
                Wy2 = pl.append(Wy2, wy2)
                Wz2 = pl.append(Wz2, wz2)

        if ntwFiles:
            if TempRed:
                ntWindFiles = glob.glob('*ntWind-%s*' % temp)
            elif MuRed:
                ntWindFiles = glob.glob('*ntWind*%s*' % temp)
            
            for ntFile in ntWindFiles:
                if checkIfEmpty(ntFile,2):
                    pass
                else:
                    ntWsq = pl.genfromtxt(ntFile, unpack=True, usecols=(0,))
                    NTW = pl.append(NTW, ntWsq)
            allTempsNTW += [[NTW]]

        allTempsE += [[E]]
        allTemps1 += [[EEcv]]
        allTemps2 += [[Ecv]]
        allTemps3 += [[dEdB]]

        allTempsSuper += [[Super]]
        allTempsWx2 += [[Wx2]]
        allTempsWy2 += [[Wy2]]
        allTempsWz2 += [[Wz2]]

    # determine length of maximum sized array
    maxLen = 0
    for arr in range(len(allTemps1)):
        if int(len(allTemps1[arr][0])) > maxLen:
            maxLen = int(len(allTemps1[arr][0]))

    # USAGE KEY:
    # len(allTemps1) gives number of temperatures
    # allTemps[i][0] gives the array for the i^th temperature
    # allTemps[i][0][j] gives the jth element of the ith temperature array

    # create energy/ specific heat file and write the header
    fout = open('ReducedEstimatorData.dat', 'w')
    fout.write('#%15s\t%16s\t%16s\t%16s\t'% (tempList[0], '','',''))
    for temp in tempList[1:]:
        fout.write('%16s\t%16s\t%16s\t%16s\t'% (temp,'','',''))
    fout.write('\n')

    # write energy and specific heat arrays to disk
    for line in range(maxLen):
        for numT in range(len(allTemps1)):
            if int(len(allTemps1[numT][0])) <= line:
                fout.write('%16s\t%16s\t%16s\t%16s\t' % ('','','',''))
            else:
                fout.write('%16.8E,\t%16.8E,\t%16.8E,\t%16.8E,\t'% (
                    float(allTempsE[numT][0][line]),
                    float(allTemps1[numT][0][line]),
                    float(allTemps2[numT][0][line]),
                    float(allTemps3[numT][0][line]) ))
        fout.write('\n')

    fout.close()

    # create superfluid file and write the header
    fout = open('ReducedSuperData.dat', 'w')
    fout.write('#%15s\t%16s\t%16s\t%16s\t'% ( tempList[0],'','','' ))
    for temp in tempList[1:]:
        fout.write('%16s\t%16s\t%16s\t%16s\t' % ( temp,'','','' ))
    fout.write('\n')

    # write superfluid stiffness arrays to disk
    for line in range(maxLen):
        for numT in range(len(allTemps1)):
            if int(len(allTempsSuper[numT][0])) <= line:
                fout.write('%16s\t%16s\t%16s\t%16s\t' % ('','','',''))
            else:
                fout.write('%16.8E,\t%16.8E,\t%16.8E,\t%16.8E,\t'% (
                    float(allTempsSuper[numT][0][line]),
                    float(allTempsWx2[numT][0][line]),
                    float(allTempsWy2[numT][0][line]),
                    float(allTempsWz2[numT][0][line]) ))
        fout.write('\n')

    fout.close()


    if biPart:
        # create file for bipartition density
        fout = open('ReducedBiPartData.dat', 'w')
        fout.write('#%15s\t%16s\t'% (tempList[0], ''))
        for temp in tempList[1:]:
            fout.write('%16s\t%16s\t'% (temp,''))
        fout.write('\n')

        # write bipartition density arrays to disk
        for line in range(maxLen):
            for numT in range(len(allTempsfD)):
                if int(len(allTempsfD[numT][0])) <= line:
                    fout.write('%16s\t%16s\t' % ('',''))
                else:
                    fout.write('%16.8E,\t%16.8E,\t'% (
                        float(allTempsfD[numT][0][line]),
                        float(allTempsbD[numT][0][line]) ))
            fout.write('\n')

        fout.close()

    # create ntWind file and write the header
    fout = open('ReducedNTWindData.dat', 'w')
    fout.write('#%15s\t'%  tempList[0] )
    for temp in tempList[1:]:
        fout.write('%16s\t' %  temp )
    fout.write('\n')

    # write ntWind arrays to disk
    for line in range(maxLen):
        for numT in range(len(allTempsNTW)):
            if int(len(allTempsNTW[numT][0])) <= line:
                fout.write('%16s\t' % '')
            else:
                fout.write('%16.8E,\t'% (
                    float(allTempsNTW[numT][0][line]) ))
        fout.write('\n')

    fout.close()



# =============================================================================
