# =============================================================================
# Main driver for pushing multiple array jobs to the VACC.
#
# NOTE: Needs more thorough comments before being committed to TRUNK.
#
# NOTE: This script assumes you will not have more than 1000 different
# random number seeds.
#
# Author:           Max Graves
# Last Revised:     28-JAN-2014
# =============================================================================

import os,argparse,re,sys,glob,shutil
import paramiko
import getpass
import clusterTools as cT
import pylab as pl

def main():
    
    args = cT.parseCMD()
    
    passwd = cT.Credentials(args.UserName)
    
    # define all data you want to pull from pimc data files IN ORDER.
    fileTypes = ['log','estimator','super','bipart_dens','ntWind']
    colNums = [
            [4,11,12,13,-4], 
            [0,1,2,3], 
            [0,1], 
            [0,1,2] ]
    observables = [
            ['E','Cv1','Cv2','Cv3','N'],
            ['rho_s/rho','Wx^2','Wy^2','Wz^2'],
            ['filmDens','bulkDens'],
            ['W^2','W_{2d}^2','W_{3d}^2'] ]

    # grab types of estimators from above.  Remove log file.
    estimTypes = list(fileTypes)
    if 'log' in estimTypes:
        estimTypes.pop(estimTypes.index('log'))
    
    # create list of reduced file names for combined data.
    fileNames = []
    for f in estimTypes:
        fileNames.append(str('Reduced'+str(f.capitalize())+'Data.dat'))

    print '\nCreating: ',fileNames,'\n'

    # if data doesn't exist, pull it from cluster
    if not args.pulled:

        print 'Pulling data from cluster.'

        # create ssh and sftp instances
        ssh = paramiko.SSHClient() 
        ssh.load_host_keys(os.path.expanduser(
            os.path.join("~", ".ssh", "known_hosts")))
        ssh.connect('bluemoon-user2.uvm.edu', username=args.UserName, 
                password=passwd)
        sftp = ssh.open_sftp()

        # move to desired directory on cluster
        sftp.chdir(args.targetDir)

        # create list of only seed directory names, get rid of other things
        allStuffInDir = sftp.listdir()
        #for s in seedDirs: # laptop didnt like this method sometimes. wtf?
        #    if s[:4] != 'seed':
        #        seedDirs.pop(seedDirs.index(s))
        seedDirs = []
        for s in range(len(allStuffInDir)):
            if 'seed' in allStuffInDir[s]:
                seedDirs.append( allStuffInDir[s] )
        
        # pull all requested file types into organized directories
        for s in seedDirs:
            newName = cT.returnSeedDirName(s)
            if os.path.exists("./"+newName):
                sys.exit(newName+' already exists.')
            os.makedirs("./"+newName)
            os.chdir("./"+newName)
            sftp.chdir('./'+s+'/OUTPUT/')
            
            moose = sftp.listdir()
            for thing in fileTypes:
                for m in moose:
                    if thing in m:
                        sftp.get('./'+m, './'+m)
                print 'pulled ',thing,' files for ',newName

            sftp.chdir('../..')
            os.chdir('..')
        
        # Rename files to have seed number replace first three numbers of pimcID.
        # Optionally delete all seed directories that were pulled from cluster.
        cT.renameFilesInDirecs(args.delDir)
      
        # check for repeated pimcIDs -- broken.
        #cT.repeatCheck()

        # close instances of sftp and ssh
        sftp.close()
        ssh.close()

    # optionally go back and bin data from files that contained correlated
    # data.  This data was collected by accident at one point, and needed
    # to be dealt with in this manner.
    if args.binAgain:
        print 'rebinning data'
        cT.reBin()

    # optionally combine all data of the same temperature into one
    # much larger array.
    if args.Crunch:
        print 'Crunching Data'
        cT.crunchData(estimTypes,colNums,observables)

    # optionally make a trimmed version of the data files that
    # makes all arrays the length of the shortest array.
    # NOT NECESSARY BUT SOMETIMES USEFUL.
    if args.trimData:
        print 'Decided to make trimmed data files'
        cT.trimData(fileNames)
        
# =============================================================================
if __name__=='__main__':
    main()
