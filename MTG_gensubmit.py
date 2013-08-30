#!/usr/bin/python 
#
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

import os,sys
from argparse import ArgumentParser

# -----------------------------------------------------------------------------
def clumeq(staticPIMCOps,numOptions,optionValue,outName):
    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.sh' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''# CLUMEQ pimc submit script\n\n
#!/bin/bash
#$ -N PIMC
#$ -P zwu-174-aa
#$ -l h_rt=0:1:0 
#$ -pe default 1
#$ -q test
#$ -cwd
#$ -S /bin/bash
echo \"Starting run at: `date`\"
case $SGE_TASK_ID in\n''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)

        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n+1,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print '\nSubmit jobs with: qsub -t 1-%d %s\n' % (numOptions,fileName)

# -----------------------------------------------------------------------------
def westgrid(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a pbs submit script for westgrid. '''

    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.pbs' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=120:00:00
#PBS -N PIMC
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}\n
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)

        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)

# -----------------------------------------------------------------------------
def sharcnet(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        name = '''out/pimc-%J'''
        if optionValue.has_key('p'):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
            
        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        scriptFile.write('sleep %d\nsqsub -q serial -o %s --mpp=1G -r 7d %s\n' % (10,name,command))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)
    print '\nSubmit jobs with: ./%s\n' % (fileName)

# -----------------------------------------------------------------------------
def scinet(staticPIMCOps,numOptions,optionValue,outName):
    ''' Write a pbs submit script for scinet. '''

    if numOptions != 8:
        print 'For scinet, must submit in multiples of 8'
        return

    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.pbs' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash
# MOAB/Torque submission script for multiple serial jobs on SciNet GPC
#
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N serialx8_pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}
# Do not send email
#PBS -M agdelma@phas.ubc.ca
#PBS -m n

# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        command = './pimc.e '
        for flag,val in optionValue.iteritems():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps.rstrip('\n')
        pbsFile.write('(sleep %02d; %s) &\n' % (10,command))
    pbsFile.write('wait')
    pbsFile.close();
    print '\nSubmit job with: qsub %s\n'%fileName

# -----------------------------------------------------------------------------
def bluemoon(staticPIMCOps,numOptions,optionValue,outName,run, initCheckPt,restartJob):
    ''' Write a pbs submit script for bluemoon '''

    # Open the pbs file and write its header
    if restartJob:
        fileName = 'reSubmit-pimc%s_bluemoon.pbs' % run
    else:
        fileName = 'submit-pimc%s_bluemoon.pbs' % run
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l nodes=1:ppn=1
#PBS -l walltime=30:00:00
#PBS -N PIMC
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID}
#PBS -m n\n
# Launch job script
cd $PBS_O_WORKDIR
case ${PBS_ARRAYID} in\n''')

    # submit an array of jobs with checkpointing enabled
    if initCheckPt:
        
        # Create the command string and make the case structure
        curDir = '$PBS_O_WORKDIR/checkPt'
        for n in range(numOptions):
            if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                command = 'cr_run ./pimc.e '
            else:
                command = 'cr_run ./pimc.e -p %d ' % (n)
            
            for flag,val in optionValue.iteritems():
                command += '-%s %s ' % (flag,val[n]) 
            command += staticPIMCOps
            
            pbsFile.write('%d)\nsleep %d\n' % (n,10*n))
            pbsFile.write('echo "Starting run %s to output file at:  `date` " >> %s/pimcCP%s.out\n'% (n+1,curDir,n+1))
            pbsFile.write(r'echo "Sleeping 29 hours to allow code to run, then checkpointing... \n\n"')       
            pbsFile.write('\n\ncd $PBS_O_WORKDIR\n\n')
            pbsFile.write('%s' % command)
            pbsFile.write(' >> %s/pimcCP%s.out 2>&1 &\n\n' % (curDir,n+1))
            pbsFile.write('BGPID%s=$!\n' % (n+1))
            pbsFile.write('sleep 104400\n')
            pbsFile.write('cr_checkpoint -p $BGPID%s -f %s/job%s.checkpoint --term\n' % (n+1,curDir,n+1))
            pbsFile.write('echo "Finished run %s at: `date`"\n' % (n+1))
            pbsFile.write(';;\n##====================================================================================\n\n')

        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        pbsFile.close();
        print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)
      
    # restart an array of jobs
    elif restartJob:
        curDir = '$PBS_O_WORKDIR/checkPt'
        for n in range(numOptions):
            pbsFile.write('%d)\nsleep %d\n' % (n,10*n))
            pbsFile.write('echo "restarting run %s to output file at:  `date` " >> %s/pimcCP%s.out\n' % (n+1,curDir,n+1))
            pbsFile.write(r'echo "Sleeping 29 hours to allow code to run, then checkpointing... \n\n"')       
            pbsFile.write('\n\ncd $PBS_O_WORKDIR\n\n')
            pbsFile.write('cr_restart %s/job%s.checkpoint' % (curDir,n+1))
            pbsFile.write(' >> %s/pimcCP%s.out 2>&1 &\n\n' % (curDir,n+1))
            pbsFile.write('BGPID%s=$!\n' % (n+1))
            pbsFile.write('sleep 104400\n')
            pbsFile.write('cr_checkpoint -p $BGPID%s -f %s/job%s.checkpoint --term\n' % (n+1,curDir,n+1))
            pbsFile.write('echo "Finished run %s at: `date`"\n' % (n+1))
            pbsFile.write(';;\n##====================================================================================\n\n')

        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        pbsFile.close();
        print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)
    
    else:
        # submit an array of jobs without checkpointing enabled...not useful probably
        # Create the command string and make the case structure
        for n in range(numOptions):
            if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                command = './pimc.e '
            else:
                command =  './pimc.e -p %d ' % (n)

            for flag,val in optionValue.iteritems():
                command += '-%s %s ' % (flag,val[n])
            command += staticPIMCOps
            pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,10*n,command))
        
        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        pbsFile.close();
        
        print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)
        print 'THE SUBMIT SCRIPT DOESNT ENABLE CHECKPOINTING!!'



# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = ArgumentParser(description="Build submission scripts for various clusters") 
    parser.add_argument("file", help='configuration file')
    parser.add_argument("cluster", metavar="cluster",
            choices=['clumeq','westgrid','sharcnet','scinet','bluemoon'],\
            help="target cluster: [clumeq,westgrid,sharcnet,scinet,bluemoon]") 
    parser.add_argument('-r', type=str, dest="run", default="",
            help="optional JobId number that will be added to the scripts name")
    parser.add_argument('-c', "--CheckPtInit", action="store_true", dest="initCheckPt",
            default=False,
            help="start a job on bluemoon with TORQUE checkpointing enabled")
    parser.add_argument('-R', "--restartJob", action="store_true", dest="restartJob",
            default=False,
            help="RESTART a job on bluemoon with TORQUE checkpointing enabled")


    # parse the command line options
    args = parser.parse_args() 
    inFileName = args.file

    if (not args.cluster):
        parser.error("need to specify a cluster")

    # We open up the input file, and read in all lines.
    inFile = open(inFileName,'r')
    inLines = inFile.readlines();

    # The first line of the file contains all static pimc options
    staticPIMCOps = inLines[0].rstrip('\n')

    # The next lines contains the short-name of the pimc option
    # and a list of values.  
    optionValue = {}
    numOptions = []
    for line in inLines[1:]:
        option = line.split()
        flag = option.pop(0)
        # Now we determine if we have a range or if we have actually included the values
        if option[0].find(':') != -1:

            # determine if we have floats or ints
            if option[0].find('.') != -1:
                type_converter = lambda x: float(x)
            else:
                type_converter = lambda x: int(x)

            dataRange = option[0].split(':')

            # Parse the range string
            dataMin = type_converter(dataRange[0])
            dataMax = type_converter(dataRange[1])
            dataStep = type_converter(dataRange[2])

            # construct the option list
            vals = [dataMin]
            while vals[-1] < dataMax:
                vals.append(vals[-1]+dataStep)

            # assign the option list
            optionValue[flag] = vals
            numOptions.append(len(vals))
        else:
            # store the typed out values
            optionValue[flag] = option
            numOptions.append(len(option))
    
    # We try to extract the temperature and volume to make an output string
    outName = ''
    findInput = staticPIMCOps.split();
    n = 0
    for input in findInput:
        if input == '-T':
            outName += '-%06.3f' % float(findInput[n+1])
            break
        n += 1
    n = 0
    for input in findInput:
        if input == '-L':
            outName += '-%07.3f' % float(findInput[n+1])
            break
        n += 1

    # We make sure all option strings have the same length 
    if len(numOptions) > 0:
        for n in numOptions:
            if n != numOptions[0]:
                print 'Not all parameters have the same number of values!'
                sys.exit()
        numOptions = numOptions[0]

    if args.cluster == 'westgrid':
        westgrid(staticPIMCOps,numOptions,optionValue,outName)

    if args.cluster == 'sharcnet':
        sharcnet(staticPIMCOps,numOptions,optionValue,outName)

    if args.cluster == 'scinet':
        scinet(staticPIMCOps,numOptions,optionValue,outName)

    if args.cluster == 'clumeq':
        clumeq(staticPIMCOps,numOptions,optionValue,outName)
    
    if args.cluster == 'bluemoon':
        bluemoon(staticPIMCOps,numOptions,optionValue,outName,args.run,args.initCheckPt,args.restartJob)


# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
