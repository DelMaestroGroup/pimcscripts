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
import natsort


def main():

    # unique tag for the name of the job to be presented in the scheduler.
    uTag = 'S05-T0.50'

    # set number of equilibration steps
    equilNum = 0

    # set number of bins to try
    binNum = 50000
    
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
        
        # get names of log files
        logFileNames = []
        for f in sftp.listdir():
            if '-log-' in f:
                logFileNames.append(f)
        
        # Naturally sort the log files (account for +/- floats in name).
        # NOTE:  This works for gce- files with all the same string
        #   but varying chemical potential.  Any other use must be 
        #   thought about first as this hasn't been tested for other uses!!!
        logFileNames = natsort.natsorted(logFileNames)

        # build array of restart strings from naturally sorted logfiles
        restartStrings = []
        
        oneLessLine = False
        for logFileName in logFileNames:
            with sftp.open(logFileName) as inFile:
                for n,line in enumerate(inFile):
                    # check for whether the restart string in gce-log has been
                    # moved up one line.  We want to write this to the resubmit
                    # script.
                    lineUp = False
                    if n == 1:
                        if len(line) < 10:
                            lineUp = True
                            continue
                        else:
                            restartStrings.append(line[2:-1])
                    if n == 1:
                        # check if log file 2nd line is blank (some are, some aren't...
                        # this is caused by bug in converting trestles (xsede) scripts
                        # over to bluemoon scripts)
                        if line == '':
                            oneLessLine = True
                            continue
                        else:
                            restartStrings.append(line[2:-1])
                    if n == 2:
                        if not lineUp:
                            restartStrings.append(line[2:-1])
                        if not oneLessLine:
                            restartStrings.append(line[2:-1])

        for nr, restartStr in enumerate(restartStrings):
            restartStrings[nr]+=' >> ${PBS_O_WORKDIR}/out/pimc-0.out 2>&1'
            restartStrings[nr] = re.sub(r'-E\s\d+',r'-E '+str(equilNum),restartStrings[nr])
            restartStrings[nr] = re.sub(r'-S\s\d+',r'-S '+str(binNum),restartStrings[nr])
            restartStrings[nr] = re.sub(r'-W\s\d+',r'-W '+str(29),restartStrings[nr])
        
        sftp.chdir('..')
        
        # check if no restart string exists, in this case just skip the seed
        # (hackyyyy -- needs to start it over in an intelligent manner.)
        if restartStrings == []:
            print 'seed number ',str(seedNum),' is empty'
            sftp.chdir('..')
            continue
        else:

            numStr = 0
            reSubName = 'resubmit-pimc.pbs'

            # determine name of the original submit script.
            for f in sftp.listdir():
                if 'submit-pimc' in f:
                    subFileName = f

            # determine if original submit script was from trestles or bluemoon
            tresTOblue = False
            with sftp.open(subFileName) as inFile:
                for n, line in enumerate(inFile):
                    if n == 6:
                        if 'uvm104' in line:
                            tresTOblue = True

            # write resubmit script for the case of original submit script being
            # from trestles.
            if tresTOblue:

                with sftp.open(reSubName,'w') as outFile:

                    # write beginning statements set up with the appropriate torque
                    # parameters.  Also, submit to the tmp directory for faster i/o.
                    outFile.write('#!/bin/bash\n\
                            \n#PBS -S /bin/bash\
                            \n#PBS -l pmem=1gb,pvmem=1gb\
                            \n#PBS -l nodes=1:ppn=1\
                            \n#PBS -l walltime=30:00:00\
                            \n#PBS -N pimc-%s\
                            \n#PBS -V\
                            \n#PBS -j oe\
                            \n#PBS -o out/pimc-${PBS_JOBID}\
                            \n#PBS -m n\n\
                            \nmkdir /tmp/${PBS_JOBID}\
                            \ncd /tmp/${PBS_JOBID}\
                            \nmkdir OUTPUT\
                            \ngzip ${PBS_O_WORKDIR}/OUTPUT/*\
                            \ncp -r ${PBS_O_WORKDIR}/OUTPUT/* ./OUTPUT/\
                            \ngunzip OUTPUT/*\
                            \necho "Starting PBS script submit-pimc.pbs at:`date`" \
                            \necho "  host:       ${PBS_O_HOST}"\
                            \necho "  node:       `cat ${PBS_NODEFILE}`"\
                            \necho "  jobid:      ${PBS_JOBID}"\n\n' % uTag)
                    
                    # write out the submit string appropriate to the job.
                    with sftp.open(subFileName) as inFile:
                        for n, line in enumerate(inFile):
                            if line[:4] == 'pimc':
                                outFile.write(restartStrings[numStr])
                                numStr += 1
                            else:
                                continue
                    
                    # write closing statements to bring data back from the node.
                    outFile.write('\n\ngzip OUTPUT/*\
                            \ncp OUTPUT/* ${PBS_O_WORKDIR}/OUTPUT\
                            \ncd ${PBS_O_WORKDIR}\
                            \ngunzip ${PBS_O_WORKDIR}/OUTPUT\
                            \nrm -r /tmp/${PBS_JOBID}')
            
            # write resubmit script for the case of original submit script being
            # from bluemoon.
            else:
                with sftp.open(subFileName) as inFile, sftp.open(reSubName,'w') as outFile:
                    for n, line in enumerate(inFile):
                        if line[:4] == 'pimc':
                            outFile.write(restartStrings[numStr])
                            numStr += 1
                        elif r'#PBS -N' in line:
                            outFile.write(line[:-2]+uTag+'\n')
                        elif r'mkdir OUTPUT' in line:
                            outFile.write(line)
                            outFile.write('gzip ${PBS_O_WORKDIR}/OUTPUT/*\n')
                            outFile.write('cp -r ${PBS_O_WORKDIR}/OUTPUT/* ./OUTPUT/\n')
                            outFile.write('gunzip OUTPUT/*\n')
                        else:
                            outFile.write(line)

            print restartStrings[0],'\n'

            # my dumb work-around for allowing writing on cluster to finish.
            time.sleep(5)

            # -----------------------------------------------------------------  
            # optionally submit jobs
            if args.submitJobs:

                #if restartStrings == []:
                #    continue
                #else:

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
