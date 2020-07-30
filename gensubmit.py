#! /usr/bin/env python
#
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

from __future__ import print_function 
import os,sys,glob,stat
import argparse

# -----------------------------------------------------------------------------
def nasa(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
         queue=None):
    ''' Write a submit script for nasa Pleides. '''

    # determine if we are adding external timing
    time_cmd = ''
    if time:
        time_cmd = '/bin/time -v '

    # determine how many nodes we need
    num_cpus = {'broadwell':28}
    num_nodes = numOptions // num_cpus['broadwell'] + 1
    cpu_type = 'broadwell'

    # Open the pbs file and write its header
    fileName = f'submit-pimc{outName}.pbs'
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n\n''')
    pbsFile.write(f'#PBS -l walltime={walltime}\n')
    pbsFile.write('#PBS -q long\n')
    pbsFile.write(f'#PBS -l select={num_nodes}:ncpus={num_cpus[cpu_type]}:model=bro\n')
    pbsFile.write('''#PBS -N PIMC
#PBS -V
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
    pbsFile.write(f'seq 0 {numOptions-1} | parallel -j {num_cpus[cpu_type]} -u \
--sshloginfile $PBS_NODEFILE "cd $PWD; ./input-pimc{outName}.sh {{}}"')
    pbsFile.close();

    # Now we create the case-file structure to submit the jobs
    fileName = f'input-pimc{outName}.sh'
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash

echo \"Starting run at: `date`\"

# get the jobid
jobid=$1

case ${jobid} in\n''')
    # Create the command string and make the case structure
    for n in range(numOptions):
        if ('p' in optionValue or staticPIMCOps.find('-p') != -1):
            command = time_cmd + './pimc.e '

        else:
            command = time_cmd + './pimc.e -p %d ' % (n)

        for flag,val in optionValue.items():
            if len(flag) == 1:
                command += f"-{flag} {val[n]} "
            else:
                command += f"--{flag}={val[n]} "
        command += staticPIMCOps 
        command += ' > out/output.$jobid 2>&1'
        pbsFile.write(f'{n})\n{command}\n;;\n')
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    os.chmod(fileName,0o744)
    
    print(f'\nSubmit jobs with: qsub -J 0-{numOptions-1:d} {fileName}\n')

# -----------------------------------------------------------------------------
def vacc(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
         queue=None):
    ''' Write a submit script for the VACC. '''

    # determine if we are adding external timing
    time_cmd = ''
    if time:
        time_cmd = '/bin/time -v '

    # Open the pbs file and write its header
    filename = ''
    if queue == 'torque':
        fileName = f'submit-pimc{outName}.pbs'
        pbsFile = open(fileName,'w')
        pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n\n''')
        pbsFile.write(f'#PBS -l walltime={walltime}\n')
        pbsFile.write('''#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_ARRAYID}-${PBS_JOBID}\n 
# Do not send email\n''')
        pbsFile.write(f'#PBS -M {os.environ["USER"]}@uvm.edu\n')
        pbsFile.write('''#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')


    elif queue=='slurm':
        fileName = f'submit-pimc{outName}.sh'
        pbsFile = open(fileName,'w')
        pbsFile.write('''#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --job-name=pimc
#SBATCH --output=out/pimc-%A-%a.log''')
        pbsFile.write(f'\n#SBATCH --time={walltime}\n')
        pbsFile.write(f'#SBATCH --mail-user={os.environ["USER"]}@uvm.edu\n')
        pbsFile.write('''#SBATCH --mail-type=FAIL,ARRAY_TASKS

# change to the directory where you submitted this script
cd ${SLURM_SUBMIT_DIR}

# Check if the out directory exists, if not, create it
if [ ! -d "./out" ]; then
  mkdir out
fi

# Executable section: echoing some Slurm data
echo "Starting sbatch script `basename $0` at:`date`"
echo "Running host:    ${SLURMD_NODENAME}"
echo "Assigned nodes:  ${SLURM_JOB_NODELIST}"
echo "Job ID:          ${SLURM_JOBID}"

''')
        # pbsFile.write(f'\n# The common PIMC options\nopts="{staticPIMCOps}"\n\n')
        pbsFile.write('''# The job array
case ${SLURM_ARRAY_TASK_ID} in
''')

    # Create the command string and make the case structure
    for n in range(numOptions):
        if ('p' in optionValue or staticPIMCOps.find('-p') != -1):
            command = time_cmd + './pimc.e '

        else:
            command = time_cmd + './pimc.e -p %d ' % (n)

        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps #'$opts'
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    if queue == 'torque':
        print(f'\nSubmit jobs with: qsub -t 0-{numOptions-1:d} {fileName}\n')
    else:
        print(f'\nSubmit jobs with: sbatch --array=0-{numOptions-1:d} {fileName}\n')

# -----------------------------------------------------------------------------
def westgrid(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
             queue=None):
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
        if ('p' in optionValue or staticPIMCOps.find('-p') != -1):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)

        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print('\nSubmit jobs with: qsub -t 0-%d %s\n' % (numOptions-1,fileName))

# -----------------------------------------------------------------------------
def sharcnet(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
            queue=None):
    ''' Write a submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        name = '''out/pimc-%J'''
        if 'p' in optionValue:
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
            
        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        scriptFile.write('sleep %d\nsqsub -q serial -o %s --mpp=1G -r 6d %s\n' % (10,name,command))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)

# -----------------------------------------------------------------------------
def local(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
          queue=None):
    ''' Write a submit script for a local machine. '''

    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# local pimc submit script\n\n''')

    # Create the command string and output to submit file
    for n in range(numOptions):
        if 'p' in optionValue:
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
            
        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        scriptFile.write('nohup nice -n 10 %s > out/out_%03d.dat &\n' % (command,n))
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)
    print('Run: ./%s'%fileName)

# -----------------------------------------------------------------------------
def scinet(staticPIMCOps,numOptions,optionValue,walltime,outName,time=False,
           queue=None):
    ''' Write a pbs submit script for scinet. '''

    if numOptions != 8:
        print('For scinet, must submit in multiples of 8')
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
        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps.rstrip('\n')
        pbsFile.write('(sleep %02d; %s) &\n' % (10,command))
    pbsFile.write('wait')
    pbsFile.close();
    print('\nSubmit job with: qsub %s\n'%fileName)

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Generate a batch submission\
                                     script from input file')
    parser.add_argument("-c", "--cluster", dest="cluster",\
                        choices=['westgrid','sharcnet','scinet','vacc','local','nasa'],\
            help="target cluster: [westgrid,sharcnet,scinet,vacc,local]",required=True) 
    parser.add_argument("-t", "--time", action="store_true", dest="time", 
                      help="Prefix the binary with a timing command?") 
    parser.add_argument("-W", "--walltime", dest="walltime", type=str, 
                      help="Specify the walltime in the format HHH:MM:SS.") 
    parser.add_argument("-q", "--queue", dest="queue", choices=['torque','slurm'],
                      help="Batch system: [torque,slurm]",default='torque') 
    parser.add_argument("-l", "--label", dest="label", help="An optional label\
                       for the submission script")
    parser.add_argument("input", help="the input file")

    # parse the command line options
    args = parser.parse_args() 
    inFileName = args.input

    # We open up the input file, and read in all lines.
    inFile = open(inFileName,'r')
    inLines = inFile.readlines();

    # Maps cluster names to functions
    gen = {'westgrid':westgrid, 'sharcnet':sharcnet, 'scinet':scinet,
           'vacc':vacc, 'local':local, 'nasa':nasa}

    # setup the default and user specified walltime
    wtime = {'westgrid':'00:00:00', 'sharcnet':'00:00:00', 'scinet':'00:00:00',
             'vacc':'30:00:00', 'local':'00:00:00', 'nasa':'10:00:00'}
    if args.walltime:
        wtime[args.cluster] = args.walltime

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
    
    # We try to extract the temperature, volume and timestep to make an output string
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
    n = 0
    for input in findInput:
        if input == '-t':
            outName += '-%07.5f' % float(findInput[n+1])
            break
        n += 1

    outName += (f'-{args.label}' if args.label else '')

    # We make sure all option strings have the same length 
    if len(numOptions) > 0:
        for n in numOptions:
            if n != numOptions[0]:
                print('Not all parameters have the same number of values!')
                sys.exit()
        numOptions = numOptions[0]

    # create the out directory that qsub will write status files to
    if not os.path.exists('out'):
        os.makedirs('out')

    # Generate the submission script
    gen[args.cluster](staticPIMCOps,numOptions,optionValue,wtime[args.cluster],
                         outName,time=args.time,queue=args.queue)


# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
