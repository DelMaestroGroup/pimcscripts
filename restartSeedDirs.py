# =============================================================================
# Main driver for restarting lots of jobs at once on VACC.  This works only
# for restarting single temperature, single seed jobs.  It requires a 
# directory structure like:
#      TX.X/(seedXXX)/(out,OUTPUT,submit-pimc.pbs)
# with TX.X being a single temperature with multiple seeds seedXXX and 
# each seed has only one job submitted within it such that only one
# set of gce-data files is contained in OUTPUT.
#
# NOTE: This script assumes you will not have more than 1000 different
# random number seeds (000-999).
#
# usage: python restartSeedDirs.py -t './path/to/TX.X/' push -L low -H high
#
# Author:           Max Graves
# Last Revised:     15-MAR-2014
# =============================================================================


import os,argparse,re,sys,glob,shutil,subprocess,getpass,time
import paramiko
import clusterTools as cT
import numpy as np


def main():

    # unique tag
    uTag = 'S5_T0.5'

    # set number of equilibration steps
    equilNum = 0

    # set number of bins to try
    binNum = 500
    
    # -------------------------------------------------------------------------
    # commands to add directories to path where blitz, boost, pimc sit.
    # not sure why, but this must be done every time for paramiko.
    expLibs = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/m/t/mtgraves/local/lib ; '
    expPBSstuff = 'export PATH=/opt/pbs/bin/:$PATH ; '
    # -------------------------------------------------------------------------

    # parse cmd line
    args = cT.parseCMD()

    # get username and password for VACC
    passwd = cT.Credentials(args.UserName)

    # create ssh and sftp instances
    ssh = paramiko.SSHClient() 
    ssh.load_host_keys(os.path.expanduser(
        os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect('bluemoon-user2.uvm.edu', username=args.UserName, password=passwd)
    sftp = ssh.open_sftp()

    # Move to desired directory on cluster.
    sftp.chdir(args.targetDir)

    # array of seed numbers, between -L and -H flag from cmd line.
    seedNums = np.arange(args.lowSeed,args.highSeed+1)

    # run through all seeds.
    for seedNum in seedNums:

        restartStr = ''

        # move to seedXXX direc.
        seedDirName = cT.returnSeedDirName(int(seedNum))
        sftp.chdir('./'+seedDirName)

        # if a resubmit file already exists, delete it.
        if 'resubmit-pimc.pbs' in sftp.listdir():
            sftp.remove('resubmit-pimc.pbs')

        # grab log file string from OUTPUT
        sftp.chdir('./OUTPUT')

        # get name of log file
        for f in sftp.listdir():
            if '-log-' in f:
                logFileName = f

        with sftp.open(logFileName) as inFile:
            for n,line in enumerate(inFile):
                if n == 2:
                    restartStr = line[2:-1]
        restartStr += ' >> ${PBS_O_WORKDIR}/out/pimc-0.out 2>&1'
        
        restartStr = re.sub(r'-E\s\d+',r'-E '+str(equilNum),restartStr)
        restartStr = re.sub(r'-S\s\d+',r'-S '+str(binNum),restartStr)
        
        print restartStr

        sftp.chdir('..')
        for f in sftp.listdir():
            if 'submit-pimc' in f:
                subFileName = f

        reSubName = 'resubmit-pimc.pbs'
        with sftp.open(subFileName) as inFile, sftp.open(reSubName,'w') as outFile:
            for n, line in enumerate(inFile):
                if line[:4] == 'pimc':
                    outFile.write(restartStr)
                elif r'#PBS -N' in line:
                    outFile.write(line[:-1]+uTag+'\n')
                elif r'mkdir OUTPUT' in line:
                    outFile.write(line)
                    outFile.write('gzip ${PBS_O_WORKDIR}/OUTPUT/*\n')
                    outFile.write('cp -r ${PBS_O_WORKDIR}/OUTPUT/* ./OUTPUT/\n')
                    outFile.write('gunzip OUTPUT/*\n')
                else:
                    outFile.write(line)

        # my dumb work-around for allowing writing on cluster to finish.
        time.sleep(5)

        # -----------------------------------------------------------------  
        # optionally submit jobs
        if args.submitJobs:

            # build submit command
            submitStuff = 'qsub '+reSubName
            changeDir = 'cd '+args.targetDir+'/'+seedDirName+' ; '
            subComm = changeDir+expLibs+expPBSstuff+submitStuff
           
            # submit the command
            stdin, stdout, stderr = ssh.exec_command(subComm)
            
            print 'Output: ',stdout.readlines()
            if stderr.readlines() != []:
                print 'Error: ',stderr.readlines()
   
        sftp.chdir('..')

    #sftp.close()
    #ssh.close()    # -- closing these breaks the writing for some reason.

# =============================================================================
if __name__=='__main__':
    main()
