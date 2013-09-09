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
            command = './pimc '
        else:
            command = './pimc -p %d ' % (n)

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
def bluemoon_blcr(staticPIMCOps,numOptions,optionValue,outName,run,initCheckPt,\
             restartJob,oneD):
    ''' Write a pbs submit script for bluemoon '''
    
    checktime=1
    maxcheck=1740
    if oneD:
        exeName = 'pimc1D'
    else:
        exeName = 'pimc'
    
    jobBase=exeName+outName
    if run:
        jobBase += '-'+run
    
    # Open the pbs file and write its header
    if restartJob:
        fileName = 'reSubmit-%s.pbs' % jobBase
    else:
        fileName = 'submit-%s.pbs' % jobBase
    pbsFile = open(fileName,'w')
    pbsFile.write('#!/bin/bash\n\n'+\
                  '#PBS -S /bin/bash\n'+\
                  '#PBS -l nodes=1:ppn=1\n'+\
                  '#PBS -l walltime=30:00:00\n'+\
                  '#PBS -N %s\n' % jobBase+\
                  '#PBS -V\n'+\
                  '#PBS -j oe\n'+\
                  '#PBS -o out/%s-${PBS_JOBID}\n' % jobBase +\
                  '#PBS -m n\n'+\
                  '#PBS -t 0-%d\n\n\n' % (numOptions-1) +\
                  'cd ${PBS_O_WORKDIR}\n'+\
                  '# Launch job script\n'+\
                  'case ${PBS_ARRAYID} in\n')
    
    # submit an array of jobs with checkpointing enabled
    if initCheckPt or restartJob:
        pbsFile.write('\n##=============================================='
                      '===================================================='
                      '===============================================\n\n')
        # Create the command string and make the case structure
        cpDir = '${PBS_O_WORKDIR}/checkPt'
        
        for n in range(numOptions):
            
            cpoutFile='%s/%s-CP%s.out'% (cpDir,jobBase,n+1)
            cpFile='%s/%s-%s.checkpoint'% (cpDir,jobBase,n+1)
            
            if initCheckPt:
                if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                    command = 'cr_run ' + exeName + ' '
                else:
                    command = 'cr_run ' + exeName + ' -p %d ' % (n)
                    
                for flag,val in optionValue.iteritems():
                    command += '-%s %s ' % (flag,val[n]) 
                command += staticPIMCOps
                
            else:
                
                command = 'cr_restart ' + cpFile
            
            pbsFile.write('\n%d)\nsleep %d\n' % (n,2*n))
            if initCheckPt:
                pbsFile.write('echo "Starting run %s to output file at:'
                              '  `date` " >> %s\n' % (n+1,cpoutFile) )
                pbsFile.write('echo "Starting run %s to output file at:'
                             '  `date` " \n' % (n+1) )
                pbsFile.write('echo "Host:" \n' )
                pbsFile.write('echo ${PBS_O_HOST} \n' )
                pbsFile.write('cat ${PBS_NODEFILE} \n' )
            else:
                pbsFile.write('if [[ ! -f %s ]]\n' % cpFile)
                pbsFile.write('then\n')
                pbsFile.write('\techo "File %s does not exist!"\n'\
                              % cpFile)
                pbsFile.write('\texit 0\n')
                pbsFile.write('fi\n')
                pbsFile.write('echo "Restarting run %s to output file at:'
                              '  `date` " >> %s\n'\
                               % (n+1,cpoutFile))
                pbsFile.write('echo "Restarting run %s to output file at:'
                              '  `date` " \n' % (n+1))
                pbsFile.write('echo "Host:" \n' )
                pbsFile.write('echo ${PBS_O_HOST} \n' )
                pbsFile.write('cat ${PBS_NODEFILE} \n' )
            pbsFile.write('echo "Will checkoint if not complete before '
                          'wallclock limit of 29 hours..."\n\n')       
            pbsFile.write('%s >> %s 2>&1 &\n\n' % (command,cpoutFile) )
            pbsFile.write('BGPID=$!\n')
            pbsFile.write('n=0\n')
            pbsFile.write('checkpointed=0\n')
            pbsFile.write('running=$(ps | grep -wce "$BGPID")\n')
            pbsFile.write('while [[ $running -eq 1 ]]\n')
            pbsFile.write('do\n')
            pbsFile.write('\tif [[ $n -eq %d ]]\n' % maxcheck)
            pbsFile.write('\tthen\n')
            pbsFile.write('\t\techo "Reached wall clock limit!"\n')
            pbsFile.write('\t\techo "checkpointing job..."\n')
            pbsFile.write('\t\tcr_checkpoint -p $BGPID -f %s --term\n'\
                          % cpFile)
            pbsFile.write('\t\tcheckpointed=1\n')
            pbsFile.write('\t\tbreak\n')
            pbsFile.write('\tfi\n')
            pbsFile.write('\tsleep %dm\n' % checktime)
            pbsFile.write('\tn=$((n + 1))\n')
            pbsFile.write('\trunning=$(ps | grep -wce "$BGPID")\n')
            pbsFile.write('done\n')
            pbsFile.write('if [[ $checkpointed -ne 1 ]]\n')
            pbsFile.write('then\n')
            pbsFile.write('\techo "Job completed!"\n')
            if restartJob:
                pbsFile.write('\techo "Moving checkpoint files from '
                              'completed job..."\n')
                pbsFile.write('\tmv %s %s/completed \n' % (cpFile,cpDir))
            pbsFile.write('fi\n')
            pbsFile.write('echo "Finished run %s at: `date`"\n' % (n+1))
            pbsFile.write('echo "Finished run at: `date` " >> %s\n' % cpoutFile )
            pbsFile.write(';;\n##=============================================='
                          '===================================================='
                          '===============================================\n\n')
        
        pbsFile.write('esac')
        pbsFile.close();
        print '\nSubmit jobs with: qsub %s\n' % fileName
        
    else:
        # submit an array of jobs without checkpointing enabled...not useful probably
        # Create the command string and make the case structure
        for n in range(numOptions):
            if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                command = 'pimc '
            else:
                command =  'pimc -p %d ' % (n)

            for flag,val in optionValue.iteritems():
                command += '-%s %s ' % (flag,val[n])
            command += staticPIMCOps
            pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
        
        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        pbsFile.close();
        
        print '\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName)
        print 'THE SUBMIT SCRIPT DOESNT ENABLE CHECKPOINTING!!'



# -----------------------------------------------------------------------------
def bluemoon(staticPIMCOps,numOptions,optionValue,outName,run,\
             resubIDX,oneD,PIGS):
    ''' Write a pbs submit script for bluemoon using internal checkpointing'''
    
    if PIGS:    
        if oneD:
            exeName = 'pimc1DOBDM'
        else:
            exeName = 'pigs'
    else:
        if oneD:
            exeName = 'pimc1D'
        else:
            exeName = 'pimc'

    jobBase=exeName+outName
    if run:
        jobBase += '-'+run
        
    jobName=jobBase +'-'+ str(resubIDX)

    # Open the pbs file and write its header
    if resubIDX > 0:
        fileName = 'reSubmit-%s.pbs' % jobName
    else:
        fileName = 'submit-%s.pbs' % jobBase
        
    ''' Local debug.... 
    'PBS_O_WORKDIR=$PWD\n'+\
    'PBS_O_HOST=A_HOST\n'+\
    'PBS_NODEFILE=nodefile.txt\n'+\
    'PBS_JOBID=my_jobid\n'+\
    'PBS_ARRAYID=2\n\n'+\    
    '''    
    pbsFile = open(fileName,'w')
    pbsFile.write('#!/bin/bash\n\n'+\
                  '#PBS -S /bin/bash\n'+\
                  '#PBS -l nodes=1:ppn=1\n'+\
                  '#PBS -l walltime=30:00:00\n'+\
                  '#PBS -N %s\n' % jobName+\
                  '#PBS -V\n'+\
                  '#PBS -j oe\n'+\
                  '#PBS -o out/%s-${PBS_ARRAYID}-${PBS_JOBID}\n' % jobName +\
                  '#PBS -m n\n'+\
                  '#PBS -t 0-%d\n\n\n' % (numOptions-1) +\
                  
                  #'cd ${PBS_O_WORKDIR}\n\n'+\
                  'mkdir /tmp/${PBS_JOBID}_${PBS_ARRAYID}\n' +\
                  'cd /tmp/${PBS_JOBID}_${PBS_ARRAYID}\n' +\
                  'mkdir OUTPUT\n' +\
                  
                  'echo "Starting PBS script %s at:`date`" \n' % fileName +\
                  'echo "\thost:\t\t${PBS_O_HOST}"\n' +\
                  'echo "\tnode:\t\t`cat ${PBS_NODEFILE}`"\n' +\
                  'echo "\tjobid:\t\t${PBS_JOBID}"\n' +\
                  'echo "\tarray job:\t${PBS_ARRAYID}" \n\n'+\
                  'case ${PBS_ARRAYID} in\n'+\
                  '\n##=============================================='
                  '===================================================='
                  '===============================================\n')
    # Create the command string and make the case structure

    for n in range(numOptions):
    
        pimcoutFile='${PBS_O_WORKDIR}/out/%s-%s.out'% (jobBase,n)
        
        if resubIDX:
            command = 'continuePIMC.sh ${PIMCID}'
        else:
            if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                command = exeName + ' '
            else:
                command = exeName + ' -p %d ' % (n)
            for flag,val in optionValue.iteritems():
                command += '-%s %s ' % (flag,val[n]) 
            command += staticPIMCOps
            
        command += ' >> ' + pimcoutFile + ' 2>&1'
                    
        if resubIDX > 0:
            pbsFile.write('\n%d)\nsleep %d\n' % (n,2*n))
            pbsFile.write('PIMCID=`getPIMCID.sh %s`\n' % pimcoutFile)
            pbsFile.write('echo "Restarting PIMCID: ${PIMCID}..."\n')
        else:
            pbsFile.write('\n%d)\nsleep %d\n' % (n,10*n))
        
        pbsFile.write('echo "%s" \n' % command )
        pbsFile.write('%s \n\n' % command )
        # was commented here down
        if (n == numOptions -1):
            pbsFile.write('resubmit=0\n')
            pbsFile.write('for i in {0..%s}\n' % (numOptions-1) )
            pbsFile.write('do\n')
            pbsFile.write('\tPIMCID=`getPIMCID.sh out/%s-${i}.out`\n' % jobBase)
            pbsFile.write('\trem=`remainingBins.sh $PIMCID`\n')
            pbsFile.write('\tif [ ${rem} -gt 0 ];then\n')
            pbsFile.write('\t\tresubmit=1\n')
            pbsFile.write('\t\tbreak\n')
            pbsFile.write('\tfi\n')
            pbsFile.write('done\n\n')
            pbsFile.write('if [ ${resubmit} -eq 1 ];then\n')
            pbsFile.write('\techo "Job array NOT complete. Resubmitting job..."\n')
            pbsFile.write('\tjobdep=`jobID2arrayID.sh ${PBS_JOBID}`\n')
            pbsFile.write('\techo "qsub -W depend=afterokarray:${jobdep} reSubmit-%s-%s.pbs"\n' % (jobBase,resubIDX+1) )
            pbsFile.write('\tqsub -W depend=afterokarray:${jobdep} reSubmit-%s-%s.pbs\n' % (jobBase,resubIDX+1) )
            pbsFile.write('\tsleep 30\n')
            pbsFile.write('else\n')
            pbsFile.write('\techo "All jobs in %s have completed!"\n' % jobBase)
            pbsFile.write('fi\n')
        # end 'was' comment
        pbsFile.write(';;\n##=============================================='
                      '===================================================='
                      '===============================================\n')
    pbsFile.write('esac\n')
    pbsFile.write('echo "Finished job %s at: `date`" \n\n' % fileName )
    pbsFile.write(  'gzip OUTPUT/*\n'                       +\
                    'cp OUTPUT/* ${PBS_O_WORKDIR}/OUTPUT\n' +\
                    'cd ${PBS_O_WORKDIR}\n'                 +\
                    'rm -r /tmp/${PBS_JOBID}_${PBS_ARRAYID}\n\n')
    pbsFile.close();
    print '\nSubmit jobs with: qsub %s\n' % fileName

    # -----------------------------------------------------------------------------
def mammouth(staticPIMCOps,numOptions,optionValue,outName,run,\
             resubIDX,oneD):
    ''' Write a pbs submit script for mammouth using internal checkpointing'''

    if oneD:
        exeName = 'pimc1D'
    else:
        exeName = 'pimc'

    jobBase=exeName+outName
    if run:
        jobBase += '-'+run

    jobName=jobBase +'-'+ str(resubIDX)

    # Open the pbs file and write its header
    if resubIDX > 0:
        fileName = 'reSubmit-%s.pbs' % jobName
    else:
        fileName = 'submit-%s.pbs' % jobBase

    ''' Local debug.... 
    'PBS_O_WORKDIR=$PWD\n'+\
    'PBS_O_HOST=A_HOST\n'+\
    'PBS_NODEFILE=nodefile.txt\n'+\
    'PBS_JOBID=my_jobid\n'+\
    'PBS_ARRAYID=2\n\n'+\    
    '''    
    pbsFile = open(fileName,'w')
    pbsFile.write('#!/bin/bash\n\n'+\
                  '#PBS -S /bin/bash\n'+\
                  '#PBS -l nodes=1:ppn=1\n'+\
                  '#PBS -l walltime=20:00:00\n'+\
                  '#PBD -q qwork@mp2\n' +\
                  '#PBS -N %s_SEGVAR\n' % jobName+\
                  '#PBS -V\n'+\
                  '#PBS -j oe\n'+\
                  '#PBS -o out/%s-SEGVAR-${PBS_JOBID}\n' % jobName +\
                  '#PBS -m n\n\n'+\
                  
                  'ARRAYID=SEGVAR\n' +\
                  'cd ${PBS_O_WORKDIR}\n\n'+\
                  'echo "Starting PBS script %s at:`date`" \n' % fileName +\
                  'echo "\thost:\t\t${PBS_O_HOST}"\n' +\
                  'echo "\tnode:\t\t`cat ${PBS_NODEFILE}`"\n' +\
                  'echo "\tjobid:\t\t${PBS_JOBID}"\n' +\
                  'echo "\tarray job:\t${ARRAYID}" \n\n'+\
                  'case ${ARRAYID} in\n'+\
                  '\n##=============================================='
                  '===================================================='
                  '===============================================\n')
    # Create the command string and make the case structure

    for n in range(numOptions):

        pimcoutFile='out/%s-%s.out'% (jobBase,n)

        if resubIDX:
            command = 'continuePIMC.sh ${PIMCID}'
        else:
            if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1):
                command = exeName + ' '
            else:
                command = exeName + ' -p %d ' % (n)
            for flag,val in optionValue.iteritems():
                command += '-%s %s ' % (flag,val[n]) 
            command += staticPIMCOps

        command += ' >> ' + pimcoutFile + ' 2>&1'

        if resubIDX > 0:
            pbsFile.write('\n%d)\nsleep %d\n' % (n,2*n))
            pbsFile.write('PIMCID=`getPIMCID.sh %s`\n' % pimcoutFile)
            pbsFile.write('echo "Restarting PIMCID: ${PIMCID}..."\n')
        else:
            pbsFile.write('\n%d)\nsleep %d\n' % (n,2*n))

        pbsFile.write('echo "%s" \n' % command )
        pbsFile.write('%s \n\n' % command )

        pbsFile.write(';;\n##=============================================='
                      '===================================================='
                      '===============================================\n')
    pbsFile.write('esac\n')
    pbsFile.write('echo "Finished job %s at: `date`" \n' % fileName )
    pbsFile.close();
    print '\nSubmit jobs with: qsub %s\n' % fileName

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = ArgumentParser(description="Build submission scripts for various clusters") 
    parser.add_argument("file", help='configuration file')
    parser.add_argument("--cluster", metavar="cluster",
            choices=['clumeq','westgrid','sharcnet','scinet','bluemoon','mammouth'],\
            help="target cluster: [clumeq,westgrid,sharcnet,scinet,bluemoon,mammouth]",
            default='bluemoon') 
    parser.add_argument('-r', type=str, dest="run", default="",
            help="optional JobId number that will be added to the scripts name")
    parser.add_argument('-c', "--CheckPtInit", action="store_true", dest="initCheckPt",
            default=False,
            help="start a job on bluemoon with TORQUE checkpointing enabled")
    parser.add_argument('-R', "--restartJob", action="store_true", dest="restartJob",
            default=False,
            help="RESTART a job on bluemoon with TORQUE checkpointing enabled")
    parser.add_argument("--oneD", action="store_true", dest="oneD",default=False,\
                        help="Run a 1D job")
    parser.add_argument("--pigs", action="store_true", dest="PIGS",default=False,\
                                            help="Run a PIGS job")
    parser.add_argument("-i","--resubIDX", dest="resubIDX",type=int,default=0,\
                            help="set resubmit index")
        
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
            epsilon = 10**(-12)
            vals = [dataMin]
            while vals[-1] < dataMax:
                if (dataMax - vals[-1]) < epsilon:
                    break
                if abs(vals[-1]+dataStep) < epsilon:
                    vals.append(0.0)
                else:
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
    '''
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
    '''
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
        bluemoon(staticPIMCOps,numOptions,optionValue,outName,args.run,\
                    args.resubIDX,args.oneD,args.PIGS)
    
    if args.cluster == 'mammouth':
        mammouth(staticPIMCOps,numOptions,optionValue,outName,args.run,\
                                args.resubIDX,args.oneD)


# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
