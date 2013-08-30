# =============================================================================
# This script assumes you will not have more than 1000 different
# random number seeds.
# =============================================================================

import os,argparse,re,sys,glob,shutil,subprocess,getpass
import paramiko
import clusterTools as cT
import numpy as np

def main():

    args = cT.parseCMD()

    passwd = cT.Credentials(args.UserName)

    # create ssh and sftp instances
    ssh = paramiko.SSHClient() 
    ssh.load_host_keys(os.path.expanduser(
        os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect('bluemoon-user2.uvm.edu', username=args.UserName, password=passwd)
    sftp = ssh.open_sftp()

    # Move to desired directory on cluster from home and if it doesn't 
    # exist, then create it.  This fails if targetDir/seedXXX exists already.
    try:
        sftp.chdir(args.targetDir)
    except:
        sftp.mkdir(args.targetDir)
        sftp.chdir(args.targetDir)

    seedNums = np.arange(args.lowSeed,args.highSeed+1)
    print seedNums

    for seedNum in seedNums:
        # create seedXXX direc.
        seedDirName = cT.returnSeedDirName(int(seedNum))
        sftp.mkdir(seedDirName)
        sftp.chdir('./'+seedDirName)

        # create case structure inside of each seedXXX direc.
        sftp.mkdir('out')
        sftp.mkdir('OUTPUT')

        # path to gensubmit and submit file (must be named 'submit')
        genSubPath = '/home/max/Documents/Code/PIMC/SCRIPTS/MTG_CH_gensubmit.py'
        subFilePath = '/home/max/Documents/Code/PIMC/SCRIPTS/submitscripts/submit'

        # change random number seed in submit file
        with open(subFilePath) as inFile, open(subFilePath+'_temp', 'w') as outFile:
            for line in inFile:
                outFile.write( re.sub(r'-p 000', r'-p '+str(seedNum), line) )

        # run gensubmit
        command = ('python '+genSubPath+' '+subFilePath+'_temp --cluster bluemoon')
        subprocess.check_call(command, shell=True)

        # copy submit file over to bluemoon
        subFile = glob.glob('*submit-pimc*')
        if len(subFile) != 1:
            sys.exit('you have more than one submit file here')
        sftp.put('./'+subFile[0], './'+subFile[0])

        # optionally submit jobs
        if args.submitJobs:
            # build submit command
            expLibs = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/m/t/mtgraves/local/lib ; '
            expPBSstuff = 'export PATH=/opt/pbs/bin/:$PATH ; '
            submitStuff = 'qsub '+subFile[0]
            changeDir = 'cd '+args.targetDir+'/'+seedDirName+' ; '
            subComm = changeDir+expLibs+expPBSstuff+submitStuff
           
            # submit the command
            stdin, stdout, stderr = ssh.exec_command(subComm)
            
            print 'Output: ',stdout.readlines()
            if stderr.readlines() != []:
                print 'Error: ',stderr.readlines()
   
        sftp.chdir('..')

    sftp.close()
    ssh.close()

# =============================================================================
if __name__=='__main__':
    main()
