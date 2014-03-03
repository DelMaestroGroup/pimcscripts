# =============================================================================
# Main driver for pulling data from multiple array jobs from VACC.
#
# NOTE: This script assumes you will not have more than 1000 different
# random number seeds (000-999).
#
# Author:           Max Graves
# Last Revised:     10-FEB-2014
# =============================================================================


import os,argparse,re,sys,glob,shutil,subprocess,getpass
import paramiko
import clusterTools as cT
import numpy as np


def main():

    equilNum = 200000
    
    # -------------------------------------------------------------------------
    # NOTE: NEW USERS WILL NEED TO CHANGE THESE STRINGS!!
    # full path to gensubmit and submit file must be supplied as below.
    #
    # NOTE:  This way of doing subFilePath has been replaced by keeping submit
    # in the same directory as stateFiles on local machine in the case of
    # submissions from equilibrated states.
    genSubPath = '/home/max/Documents/Code/PIMC/SCRIPTS/MTG_CH_gensubmit.py'
    subFilePath = '/home/max/Documents/Code/PIMC/SCRIPTS/submitscripts/submit'
    
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

    # Move to desired directory on cluster from home and if it doesn't 
    # exist, then create it.  This fails if targetDir/seedXXX exists already.
    try:
        sftp.chdir(args.targetDir)
    except:
        sftp.mkdir(args.targetDir)
        sftp.chdir(args.targetDir)

    seedNums = np.arange(args.lowSeed,args.highSeed+1)

    # keep track of directory you are in in the terminal
    workingDir = os.getcwd()


    for seedNum in seedNums:

        # create seedXXX direc.
        seedDirName = cT.returnSeedDirName(int(seedNum))
        sftp.mkdir(seedDirName)
        sftp.chdir('./'+seedDirName)

        # create directory structure inside of each seedXXX direc.
        sftp.mkdir('out')
        sftp.mkdir('OUTPUT')

        # change random number seed in submit file -- this defines hacky.
        # treat new submission different than from equilibrated states.
        if args.stateFilesDir == '':
            with open(subFilePath) as inFile, open(subFilePath+'_temp', 'w') as outFile:
                for n, line in enumerate(inFile):
                    if n==0:
                        if '-p 000' not in line:
                            sys.exit('Include -p 000 in submit script!')
                    outFile.write( re.sub(r'-p 000', r'-p '+str(seedNum), line))
        else:
            subFilePath = os.getcwd()
            if subFilePath[-1]!='/':
                subFilePath +='/'
            subFilePath += args.stateFilesDir+'/../submit'
            print subFilePath
            with open(subFilePath) as inFile, open(subFilePath+'_temp', 'w') as outFile:
                for n, line in enumerate(inFile):
                    if n==0:
                        if '-p 000' not in line:
                            sys.exit('Include -p 000 in submit script!')
                    outFile.write( re.sub(r'-p 000', r'-p '+str(seedNum), line))


        # run gensubmit
        command = ('python '+genSubPath+' '+subFilePath+\
                '_temp --cluster=bluemoon '+'-m '+args.memRequest)
        subprocess.check_call(command, shell=True)

       
        # get submit file name generated from gensubmit
        subFile = glob.glob('*submit-pimc.*')
        if len(subFile) != 1:
            sys.exit('you have more than one submit file here')
        subFile = subFile[0]
        
        # optionally set up submit script to start from equilibrated state files
        if args.stateFilesDir != '':

            # get list of (g)ce-state files in directory supplied from cmdline.
            os.chdir(args.stateFilesDir)
            stateFileList = sorted(glob.glob('*state*'))
            
            # we check that there is a directory on the cluster called
            # stateFiles that contains exactly the (g)ce-state files
            # that are in the directory supplied from cmdline.
            # This double-storage is set up to make the user think about
            # what they are doing.
            try:
                sftp.chdir('../stateFiles')
                stFs = sorted(sftp.listdir())
                if stateFileList != stFs:
                    sys.exit('ERROR: State files in targeted directory on '+
                            'cluster are different than on your machine.')
                sftp.chdir('../'+seedDirName)
            except:
                sys.exit('ERROR: Must have state files in directory called '+
                        'stateFiles in the parent directory of your jobs.')

            # now go through and replace the -E #### part of the submitscript
            # -s /path/to/(g)ce-state... .  After every second occurrence of 
            # the -E flag, we skip ahead once in the list of state file names.
            sF2 = subFile[:-4]+'_temp.pbs'
            numOccur = 0
            stateNum = 0
            clusterCWDpre = '${PBS_O_WORKDIR}/../stateFiles/'
            
            os.chdir(workingDir)
            with open(subFile) as inFile, open(sF2, 'w') as outFile:
                for n, line in enumerate(inFile):
                   
                    match = re.search(r'-E\s\d+',line)
                    if match != None:
                        if args.keepEquilibrating:
                            outFile.write( re.sub(r'-E\s\d+','-E '+str(equilNum)+
                                ' -s '+clusterCWDpre+stateFileList[stateNum], line))
                        else:
                            outFile.write( re.sub(r'-E\s\d+','-s '+clusterCWDpre+stateFileList[stateNum], line))
                    else:
                        outFile.write(line)
                   
                    if numOccur%2 == 1:
                        stateNum += 1
                    
                    if match != None:
                        numOccur += 1

            os.rename(sF2,subFile)

        # copy submit file over to bluemoon
        sftp.put('./'+subFile, './'+subFile)

        # optionally submit jobs
        if args.submitJobs:

            # build submit command
            submitStuff = 'qsub '+subFile
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
