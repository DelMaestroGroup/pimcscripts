# =============================================================================
# This script assumes you will not have more than 1000 different
# random number seeds.
# =============================================================================

import os,argparse,re,sys,glob,shutil
import paramiko
import getpass

def parseCMD():
    ''' parse the command line. '''
    parser = argparse.ArgumentParser(description='pulls down lots of files.')
    parser.add_argument('-U', '--UserName', type=str,
            default='mtgraves',
            help='Username for your VACC account.')
    parser.add_argument('-d', '--deleteDirs', action='store_true', 
            dest='delDir',default=True,
            help='Do you want to delete the individual seed direcs?')
    parser.add_argument('-s', '--submitJobs', action='store_true', 
            dest='submitJobs',default=False,
            help='Do you want to automatically submit the jobs?')
    parser.add_argument('-t', '--targetDir', type=str,
            help='Enter directory to build structures in: ~/(--targetDir)')
    parser.add_argument('-L', '--lowSeed', type=int,
            default=0,
            help='Low seed.')
    parser.add_argument('-H', '--highSeed', type=int,
            default=0,
            help='High seed. (highSeed-lowSeed = numSeeds)')

    return parser.parse_args()

def Credentials(userN):
    ''' get username and password '''
    print 'username:',userN
    passwd = getpass.getpass(prompt='password: ')
    return passwd

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
    rename all files within their directories and
    optionally delete all directories pulled from cluster.
    '''
    dirNames = glob.glob('*seed*')
    for dirName in dirNames:
        os.chdir('./'+dirName)
        filesIn = glob.glob('*.dat*')
        for zebra in filesIn:
            shutil.copy2(zebra, '../'+zebra[:-13]+dirName[-3:]+zebra[-10:])
        os.chdir('..')
        if delDir:
            shutil.rmtree('./'+dirName)
    if delDir:
        print 'Chose to delete original direcs'

# =============================================================================
