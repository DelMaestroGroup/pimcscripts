#! /usr/bin/env python
#
# rsubmit.py
# Adrian Del Maestro
# 09.04.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.  This restart script
# reads all log files in a current directory by default, or will
# take a list of id numbers from the command line

import sys,os
import pimchelp
import argparse
import math
import uuid

# -----------------------------------------------------------------------------
def getPIMCommand(fname):
    ''' Get the command string from the log file.'''
    # Open the log file and get the command string
    logFile = open(fname,'r')
    logLines = logFile.readlines()
    line = logLines[2]

    return line[2:]

# -----------------------------------------------------------------------------
def nasa(commands,walltime,outName,cmd,time=False,queue='broadwell'):
    ''' Write a submit script for nasa Pleides. '''

    # determine if we are adding external timing
    time_cmd = ''
    if time:
        time_cmd = '/bin/time -v '

    # determine if we can run in the normal queue
    num_hours = int(walltime.split(':')[0])
    if num_hours <= 8:
        qname = 'normal'
    else:
        qname = 'long'

    # determine how many nodes we need
    cpu_type = queue
    num_cpus = {'broadwell':28, 'skylake':40, 'cascade lake':40, 
                'haswell':24, 'ivy bridge':20, 'sandy bridge':16,
                'broadwell electra':28}
    model = {'broadwell':'bro', 'skylake':'sky_ele', 'cascade lake':'cas_ait', 
                'haswell':'has', 'ivy bridge':'ivy', 'sandy bridge':'san', 
             'broadwell electra':'bro_ele'}

    # make sure we have a valid model type
    if queue not in num_cpus:
        print(f'model {queue} does not exist!')
        sys.exit(1)

    # Open the pbs file and write its header
    if 'pigs' in commands[0]:
        lab = 'pigs'
    else:
        lab = 'pimc'

    num_nodes = int(math.ceil(len(commands)/num_cpus[cpu_type]))

    fileName = f'resubmit-{lab}{outName}.pbs'
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n\n''')
    pbsFile.write(f'#PBS -l walltime={walltime}\n')
    pbsFile.write(f'#PBS -q {qname}\n')
    pbsFile.write(f'#PBS -l select={num_nodes}:ncpus={num_cpus[cpu_type]}:model={model[cpu_type]}\n')
    pbsFile.write(f'#PBS -N {lab.upper()}{outName}\n')
    pbsFile.write('''#PBS -V
#PBS -j oe
#PBS -o ./out/

# Do not send email
#PBS -M adelmaes@uvm.edu
#PBS -m n

# Start job script

# Check if the out directory exists, if not, create it
if [ ! -d "./out" ]; then
  mkdir out
fi

# Change into submission directory
cd $PBS_O_WORKDIR

# Calling parallel jobs\n''')
    pbsFile.write(f'seq 0 {len(commands)-1} | parallel -j {num_cpus[cpu_type]} -u \
--sshloginfile $PBS_NODEFILE "cd $PWD; ./input-{lab}{outName}.sh {{}}"')
    pbsFile.close();

    print(f'\nSubmit jobs with: qsub {fileName}\n')

    # Now we create the case-file structure to submit the jobs
    fileName = f'input-{lab}{outName}.sh'
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash

echo \"Starting run at: `date`\"

# get the jobid
jobid=$1

case ${jobid} in\n''')
    # Create the command string and make the case structure
    for n,ccmd in enumerate(commands):
        command = time_cmd + ccmd.rstrip('\n')
        command += f' > out/output_{str(uuid.uuid4())[-12:]}.$jobid 2>&1'
        pbsFile.write(f'{n})\n{command}\n;;\n')
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    os.chmod(fileName,0o744)

# -----------------------------------------------------------------------------
def westgrid(logFileNames,outName):
    ''' Write a pbs submit script for westgrid. '''

    # Open the pbs file and write its header
    fileName = 'resubmit-pimc%s.pbs' % outName
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
echo \"Starting run at: `date`\"\n''')
    
    numID = len(logFileNames)

    if numID > 1:
        pbsFile.write('''\ncase ${PBS_ARRAYID} in\n''')
    
    # Create the command string and make the case structure
    for n,fname in enumerate(logFileNames):
        # Get the command string
        command = getPIMCommand(fname)
        if numID > 1:
            pbsFile.write('%d)\n%s;;\n' % (n,command))
        else:
            pbsFile.write(command)
    
    if numID > 1:
        pbsFile.write('esac\necho \"Finished run at: `date`\"')
        print('\nSubmit jobs with: qsub -t 0-%d %s\n' % (numID-1,fileName))
    else:
        pbsFile.write('echo \"Finished run at: `date`\"')
        print('\nSubmit jobs with: qsub %s\n' % fileName)
    pbsFile.close();

# -----------------------------------------------------------------------------
def sharcnet(logFileNames,outName):
    ''' Write a pbs submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'resubmit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Get the command string and output to submit file
    name = '''out/pimc-%J'''
    for n,fname in enumerate(logFileNames):
        command = getPIMCommand(fname)
        scriptFile.write('sqsub -q serial -o %s -r 6d %s' % (name,command))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)

# -----------------------------------------------------------------------------
def scinet(logFileNames,outName):
    ''' Write a pbs submit script for scinet. '''

    # Open the pbs file and write its header
    fileName = 'resubmit-pimc%s.pbs' % outName
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

    # Get the command string and output to submit file
    for n,fname in enumerate(logFileNames):
        command = getPIMCommand(fname).rstrip('\n')
        pbsFile.write('(%s) &\n' % command)
    pbsFile.write('wait')
    pbsFile.close();
    print('\nSubmit job with: qsub %s\n' % fileName)

# -----------------------------------------------------------------------------
# Begin Main Program
# -----------------------------------------------------------------------------
def main():

    # nasa processer resources
    nasa_cpus = ['broadwell', 'skylake', 'cascade lake', 'haswell', 'ivy bridge', 
                 'sandy bridge', 'broadwell electra']

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Generate a batch submission\
                                     script to restart finished runs.')
    parser.add_argument("-T", "--temperature", dest="T", type=float, \
            help="simulation temperature in Kelvin") 
    parser.add_argument("-N", "--number_particles", dest="N", type=int,\
            help="number of particles") 
    parser.add_argument("-n", "--density", dest="n", type=float,\
            help="number density in Angstroms^{-d}")
    parser.add_argument("-P", "--number_time_slices", dest="P", type=int,\
            help="number of time slices")
    parser.add_argument("-u", "--chemical_potential", dest="mu", type=float,\
            help="chemical potential in Kelvin") 
    parser.add_argument("-L", "--length", dest="L", type=float,\
            help="length in Angstroms") 
    parser.add_argument("-t", "--imaginary_time_step", dest="tau", type=float,\
            help="imaginary time step") 
    parser.add_argument("--canonical", action="store_true", dest="canonical", 
                      help="are we in the canonical ensemble?")
    parser.add_argument("-i", "--pimcid", action="append", type=str,\
            help="a list of PIMCIDs")
    parser.add_argument("-e", "--exclude", action="append", type=str,\
            help="a list of PIMC ID numbers to exclude")
    parser.add_argument("-c", "--cluster", dest="cluster",\
                        choices=['westgrid','sharcnet','scinet','vacc','local','nasa'],\
            help="target cluster: [westgrid,sharcnet,scinet,vacc,local]",required=True) 
    parser.add_argument("--time", action="store_true", dest="time", 
                      help="Prefix the binary with a timing command?") 
    parser.add_argument("-W", "--walltime", dest="walltime", type=str, 
                      help="Specify the walltime in the format HHH:MM:SS.") 
    parser.add_argument("-q", "--queue", dest="queue", choices=['torque','slurm'],
                      help="Batch system: [torque,slurm]",default='torque') 
    parser.add_argument("-m", "--model", dest="model", choices=nasa_cpus,
                      help="NASA CPU types.") 
    parser.add_argument("-l", "--label", dest="label", help="An optional label\
                       for the submission script")
    parser.add_argument("--pigs", help="Is this a pigs simulation?",
                        action="store_true")
    parser.add_argument("base_dir", help='The base directory\
                        where the data files to be reduced are located.',
                        default=os.getcwd(), nargs='?')

    # parse the command line options
    args = parser.parse_args() 

    base_dir = os.path.join(args.base_dir,'')

    # is this a pigs or pimc simulation?
    if args.pigs:
        cmd = 'pigs.e'
    else:
        cmd = 'pimc.e'

    # setup the default and user specified walltime
    wtime = {'westgrid':'00:00:00', 'sharcnet':'00:00:00', 'scinet':'00:00:00',
             'vacc':'30:00:00', 'local':'00:00:00', 'nasa':'10:00:00'}
    if args.walltime:
        wtime[args.cluster] = args.walltime

    # Get the data string and create the pimc helper object
    dataName,outName = pimchelp.getFileString(args)
    pimc = pimchelp.PimcHelp(dataName,args.canonical,baseDir=base_dir)

    # We get either all the log files in the current directory, or just the
    # requested files by their ID number
    logFileNames  = pimc.getFileList('log',idList=args.pimcid)
    if len(logFileNames) < 1:
        print(f'No log files have been found for pattern: {dataName}.')
        sys.exit()

    outName += (f'-{args.label}' if args.label else '')

    # If we have excluded any ID's we remove them from the list
    if args.exclude:
        for cid in options.exclude:
            for n,fname in enumerate(logFileNames):
                if cid in fname:
                    logFileNames.pop(n)

    restart_commands = []
    for logFile in logFileNames:
        restart_commands.append(getPIMCommand(logFile))

    # create the out directory that qsub will write status files to
    if not os.path.exists('out'):
        os.makedirs('out')

    queue = args.queue
    if args.cluster == 'nasa' and args.model:
        queue = args.model
    elif args.cluster == 'nasa' and not args.model:
        queue = 'broadwell'


    if args.cluster == 'nasa':
        nasa(restart_commands,wtime[args.cluster],outName,cmd,time=args.time,queue=queue)

    # Generate the submission script
    # gen[args.cluster](staticPIMCOps,numOptions,optionValue,wtime[args.cluster],
    #                      outName,cmd,time=args.time,queue=queue)

    # # Now create the submission files
    # if options.cluster == 'westgrid':
    #     westgrid(logFileNames,outName)

    # if options.cluster == 'sharcnet':
    #     sharcnet(logFileNames,outName)

    # if options.cluster == 'scinet':
    #     scinet(logFileNames,outName)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
