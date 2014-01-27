# =============================================================================
# Functions used in pull/push seedDirs script.
#
# This script assumes you will not have more than 1000 different
# random number seeds.
#
# Author:           Max Graves
# Last Revised:     7-JAN-2014
# =============================================================================

import os,argparse,re,sys,glob,shutil,subprocess
import paramiko
import getpass
import pylab as pl
import numpy as np


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


def Credentials(userN):
    ''' get username and password '''
    print 'username:',userN
    passwd = getpass.getpass(prompt='password: ')
    return passwd


def crunchData(estimTypes, colNums, observables):
    '''
    Takes a list of pimc file types ['estimator', 'super', etc..]
    and a list of list of column numbers (in order of file types) as
    [[0,1], [1,4,6], etc...] and will combine data from jobs that were 
    the same but had different random number seeds all into one
    file.
    '''

    args = parseCMD()
    reduceType = args.reduceType

    # make list of data file names
    estimFiles  = glob.glob('*%s*' % estimTypes[0])
    
    # check which ensemble
    canonical = ensembleCheck(estimFiles[0][0])

    # make list of all reduced variable, in order, based on ensemble
    indList = indepVarList(estimFiles, canonical, reduceType)
   
    for nType, estType in enumerate(estimTypes):

        # dict to hold all estType data for all independent variables (T,..)
        allTemps = {}

        # dict to hold numBins, average, stdErr for each data file
        statTemps = {}

        # loop over independent variable for given estType
        for numTemp, temp in enumerate(indList):

            # Get filenames.  The chemical potential always has a sign
            # in front and the temperature always comes after the
            # file type and a dash.
            if reduceType == 'T':
                bipFiles = glob.glob('*%s-%s*' % (estimTypes[nType], temp))
            elif reduceType == 'u':
                bipFiles = glob.glob('*%s*%s*' % (estimTypes[nType], temp))
 
            arrs    = {}
            arrs2   = {}
            for nf, f in enumerate(bipFiles):
                if checkIfEmpty(f,2):
                    print f,' is empty.'
                    pass
                else:
                    dat = pl.genfromtxt(f, unpack=True, usecols=colNums[nType])

                    # define key that can be sorted properly
                    if numTemp < 10:
                        arrName = 'arr00'+str(numTemp)
                    elif (numTemp >= 10 and numTemp < 100):
                        arrName = 'arr0'+str(numTemp)
                    else:
                        arrName = 'arr'+str(numTemp)
                    
                    # treat data from one column
                    if int(len(colNums[nType])) == 1:
                        if nf == 0:
                            arrs[arrName] = dat
                            arrs2[arrName] = [[pl.average(dat),
                                    pl.std(dat)/pl.sqrt(1.0*int(len(dat))),
                                    int(len(dat))]]
                        else:
                            arrs[arrName] = pl.append(arrs[arrName], dat)
                            arrs2[arrName] = pl.append(
                                    arrs2[arrName], [[pl.average(dat),
                                    pl.std(dat)/pl.sqrt(1.0*int(len(dat))),
                                    int(len(dat))]])

                    # treat data from multiple columns
                    else:
                        if nf == 0:
                            for n, arr in enumerate(dat):
                                arrs[arrName+'_'+str(n)] = arr
                                arrs2[arrName+'_'+str(n)] = [[pl.average(arr),
                                        pl.std(arr)/pl.sqrt(1.0*int(len(arr))),
                                        int(len(arr))]]
                        else:
                            for n, arr in enumerate(dat):
                                arrs[arrName+'_'+str(n)] = pl.append(
                                        arrs[arrName+'_'+str(n)], arr)
                                arrs2[arrName+'_'+str(n)] = pl.append(
                                        arrs2[arrName+'_'+str(n)], [[pl.average(arr),
                                            pl.std(arr)/pl.sqrt(1.0*int(len(arr))),
                                            int(len(arr))]])
           
            # construct key name.  This assumes <1000 temperatures.
            if numTemp < 10:
                allArrName = 'allTemps00'+str(numTemp)
            elif (10 <= numTemp and numTemp < 100):
                allArrName = 'allTemps0'+str(numTemp)
            else:
                allArrName = 'allTemps'+str(numTemp)

            if nType < 10:
                allArrName += '00'+str(nType)
            elif (10 <= nType and nType < 100):
                allArrName += '0'+str(nType)
            else:
                allArrName += str(nType)

               
            allTemps[allArrName] = arrs
            statTemps[allArrName] = arrs2

        # length of max. sized array for all data
        maxLen = 0
        for t in allTemps:
            for g in allTemps[t]:
                arrayLen = len(allTemps[t][g])
                if arrayLen > maxLen:
                    maxLen = arrayLen
        # open file to hold all of estType data
        newName = 'Reduced'+str(estType.capitalize())+'Data.dat'
        fout = open(newName, 'w')

        # write independent variable as top header line
        fout.write('#%15s\t' % indList[0])
        for n in range(len(colNums[nType])-1):
            fout.write('%16s\t'% '')
        for temp in indList[1:]:
            fout.write('%16s\t'% temp)
            for n in range(len(colNums[nType])-1):
                fout.write('%16s\t'% '')
        fout.write('\n')

        # write observable names as second header line
        fout.write('#%15s\t' % observables[nType][0])
        for n in range(len(observables[nType])-1):
            fout.write('%16s\t' % observables[nType][n+1])
        for temp in indList[1:]:
            for obs in observables[nType]:
                fout.write('%16s\t' % obs)
        fout.write('\n')
        
        # write data arrays to disk
        for line in range(maxLen):
            for a in sorted(allTemps.iterkeys()):
                for aa in sorted(allTemps[a]):
                    try:
                        fout.write('%16.8E,\t' % (
                            float(allTemps[a][aa][line]) ))
                    except:
                        fout.write('%16s,\t' % '')
            fout.write('\n')
        
        fout.close()

        # EXPERIMENTAL
        # length of max. sized array for averages
        maxLen2 = 0
        for t in statTemps:
            for g in statTemps[t]:
                arrayLen = len(statTemps[t][g])
                if arrayLen > maxLen2:
                    maxLen2 = arrayLen
        maxLen2 /= 3

        # open file for standard error
        newName = 'zAveraged'+str(estType.capitalize())+'Data.dat'
        fout = open(newName, 'w')

        # write independent variable as top header line
        fout.write('#%15s\t%16s\t%16s\t' % (
            indList[0], '',''))
        for n in range(len(colNums[nType])-1):
            fout.write('%16s\t%16s\t%16s\t'% ('','',''))
        for temp in indList[1:]:
            fout.write('%16s\t%16s\t%16s\t'% (temp,'',''))
            for n in range(len(colNums[nType])-1):
                fout.write('%16s\t%16s\t%16s\t'% ('','',''))
        fout.write('\n')

        # write observable names as second header line
        fout.write('#%15s\t%16s\t%16s\t' % (
            observables[nType][0],
            'std Err',
            'bins'))
        for n in range(len(observables[nType])-1):
            fout.write('%16s\t%16s\t%16s\t' % (
                observables[nType][n+1],
                'std Err',
                'bins'))
        for temp in indList[1:]:
            for obs in observables[nType]:
                fout.write('%16s\t%16s\t%16s\t' % (
                    obs, 'std Err', 'bins'))
        fout.write('\n')
        
        # write data arrays to disk
        en = 0
        for line in range(maxLen2):
            for a in sorted(statTemps.iterkeys()):
                for aa in sorted(statTemps[a]):
                    try:
                        fout.write('%15.8E,\t%15.8E,\t%15.8E,\t' % (
                            float(statTemps[a][aa][en]),
                            float(statTemps[a][aa][en+1]),
                            float(statTemps[a][aa][en+2])))
                    except:
                        fout.write('%15s,\t%15s,\t%15s,\t' % ('','',''))
                    
            fout.write('\n')
            en += 3
        
        fout.close()
       

def ensembleCheck(firstLetter):
    '''
    check the ensemble by checking the file names.
    '''
    canonical = True
    if firstLetter == 'g':
        canonical = False

    return canonical


def indepVarList(estimFiles, canonical, reduceVar):
    ''' 
    Make list of all independent variable, in order, 
    based on ensemble.  This is set up to be temperature or
    chemical potential.  If another variable is desired,
    then proper adjustments may need to be made.
    '''
    indList = pl.array([])
    for f in estimFiles:
        if reduceVar == 'T':
            if canonical:
                tempVar = f[13:19]
            else:
                tempVar = f[14:20]
        elif reduceVar == 'u':
            if canonical: 
                tempVar = f[28:35]
            else:
                tempVar = f[29:36]
        
        if tempVar not in indList:
            indList = pl.append(indList, tempVar)
    
    return pl.sort(indList)


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


def trimData(fileNames):
    '''
    Trims data files so that the length of each column is equal
    to that of the shortest column.  This is so that loadtxt will
    work properly.  This is rather naughty.
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

# =============================================================================
