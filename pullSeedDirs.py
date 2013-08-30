# =============================================================================
# This script assumes you will not have more than 1000 different
# random number seeds.
# =============================================================================

import os,argparse,re,sys,glob,shutil
import paramiko
import getpass
import clusterTools as cT

def main():

    args = cT.parseCMD()

    passwd = cT.Credentials(args.UserName)

    #fileTypes = ['log','estimator']
    fileTypes = ['log']
   
    # create ssh and sftp instances
    ssh = paramiko.SSHClient() 
    ssh.load_host_keys(os.path.expanduser(
        os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect('bluemoon-user2.uvm.edu', username=args.UserName, password=passwd)
    sftp = ssh.open_sftp()

    # move to desired directory on cluster
    sftp.chdir('./projects/specificHeatTest/')

    # create list of only seed directory names, get rid of other things
    seedDirs = sftp.listdir()
    for s in seedDirs:
        if 'seed' not in s:
            seedDirs.pop(seedDirs.index(s))
     
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
  
    # check for repeated pimcIDs
    cT.repeatCheck()

    # close instances of sftp and ssh
    sftp.close()
    ssh.close() 

# =============================================================================
if __name__=='__main__':
    main()
