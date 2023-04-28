#! /usr/bin/env python
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

from __future__ import print_function 
import os,sys,glob,stat
import argparse
import uuid
import math
import subprocess


# -----------------------------------------------------------------------------
def validate(cmd):
    '''Validate a sample executable command for debugging purposes.'''

    # this seems to have problem on NASA machines, so we will use a shell
#    output = subprocess.run(cmd.split(' ') + ['--validate'], capture_output=True).stderr.decode("utf-8")
    output = subprocess.run(cmd + ' --validate', capture_output=True,shell=True).stderr.decode("utf-8")
    if 'SUCCESS' not in output:
        sys.exit(f'Error with command: {cmd} \n {output}')

# -----------------------------------------------------------------------------
def nasa(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
         queue='broadwell'):
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

    # num_nodes = numOptions // num_cpus[cpu_type] + 1
    num_nodes = int(math.ceil(numOptions / num_cpus[cpu_type]))

    # Open the pbs file and write its header
    if 'pigs' in cmd:
        lab = 'pigs'
    else:
        lab = 'pimc'

    fileName = f'submit-{lab}{outName}.pbs'
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
    pbsFile.write(f'seq 0 {numOptions-1} | parallel -j {num_cpus[cpu_type]} -u \
--sshloginfile $PBS_NODEFILE "cd $PWD; ./input-{lab}{outName}.sh {{}}"')
    pbsFile.close();


    # Now we create the case-file structure to submit the jobs
    ifileName = f'input-{lab}{outName}.sh'
    pbsFile = open(ifileName,'w')
    pbsFile.write('''#!/bin/bash

echo \"Starting run at: `date`\"

# get the jobid
jobid=$1

case ${jobid} in\n''')
    # Create the command string and make the case structure
    for n in range(numOptions):
        if ('p' in optionValue or staticPIMCOps.find('-p') != -1):
            command = time_cmd + f'./{cmd} '

        else:
            command = time_cmd + f'./{cmd} -p %d ' % (n)

        for flag,val in optionValue.items():
            if len(flag) == 1:
                command += f"-{flag} {val[n]} "
            else:
                command += f"--{flag}={val[n]} "
        command += staticPIMCOps 

        # test to make sure we have a valid command
        if n == 0:
            validate(command)

        command += f' > out/output_{str(uuid.uuid4())[-12:]}.$jobid 2>&1'
        pbsFile.write(f'{n})\n{command}\n;;\n')
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    os.chmod(ifileName,0o744)

    print(f'\nSubmit jobs with: qsub {fileName}\n')
    
# -----------------------------------------------------------------------------
def vacc(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
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
            command = time_cmd + f'./{cmd} '

        else:
            command = time_cmd + f'./{cmd} -p {n} '

        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps #'$opts'

        # test to make sure we have a valid command
        if n == 0:
            validate(command)

        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
    
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    if queue == 'torque':
        print(f'\nSubmit jobs with: qsub -t 0-{numOptions-1:d} {fileName}\n')
    else:
        print(f'\nSubmit jobs with: sbatch --array=0-{numOptions-1:d} {fileName}\n')

# -----------------------------------------------------------------------------
def westgrid(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
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
def sharcnet(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
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
def parallel(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
          queue=None):
    user_qos="";
    user_partition="";
    user_Node=-1;
    user_avgcores=-1;
    user_cpu=-1;
    if QUE=="slurm":
        print("On most computers that use slurm as a job allocator you must request resources, thus we will help you with that\n If you do not want that do not use -q slurm as a flag and you will get a skeleton file to print in\n")
        try:
            output = subprocess.check_output(["sacctmgr", "-p", "show", "qos"])
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {e}\nskipping step")
        else:
            output_str = output.decode("utf-8")
            output_lines = output_str.strip().split("\n")
            for line in output_lines:
                print(line)
            print("\nHere are the quality of services")
            user_qos = input("Input the qos you want to use, if none leave blank\n name:")
            while True:
                try:
                    num = input("Look at the qos of your choice and type the max number of nodes you can use, if none leave blank\n Please enter an integer greater than 1: ")
                    if num != "" and int(num)>0:
                        user_Node =int(num)
                    break
                except ValueError:
                    print("Invalid input. Please enter an integer.")

        try:
            output = subprocess.check_output(["scontrol", "show", "partitions"])
        except subprocess.CalledProcessError as e:
            print(f"Error running command: {e}\nskipping step")
        else:
            output_str = output.decode("utf-8")
            output_lines = output_str.strip().split("\n")
            for line in output_lines:
                print(line)
            print("\nHere are all the possible partitions")
            user_partition = input("Input the parition you want to use, if none leave blank\n name:")
            while True:
                try:
                    num = input("Please enter the total cpu, if none leave blank\n Please enter an integer greater than 1: ")
                    if num != "" and int(num)>0:
                        user_cpu =int(num)
                    break
                except ValueError:
                    print("Invalid input. Please enter an integer.")
            try:
                cmd = ["sinfo", "-N", "-o", "%n %C", "-p", user_partition]
                output = subprocess.check_output(cmd)
            except subprocess.CalledProcessError as e:
                print(f"Error running command: {e}\nskipping step")
            else:
                output_str = output.decode("utf-8")
                output_lines = output_str.strip().split("\n")
                for line in output_lines:
                    print(line)
                print("\nHere are all the possible nodes in the partition you chose and the amount of cores on each one")
                while True:
                    try:
                        num = input("Set the number of cores per node, I recommend a avgerage, if none leave blank\n Please enter an integer greater than 1: ")
                        if num != "" and int(num)>0:
                            user_avgcores =int(num)
                        break
                    except ValueError:
                        print("Invalid input. Please enter an integer.")
        
    ''' Write a job scirpt script for a local machine. '''
    # Open the script file and write its header
    
    fileName = 'commands-pimc%s.txt' % outName
    scriptFile = open(fileName,'w')
    # Create the command string and output to submit file
    for n in range(numOptions):
        if 'p' in optionValue:
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
        for flag,val in optionValue.items():
            command += '-%s %s ' % (flag,val[n])
        command += staticPIMCOps
        scriptFile.write('srun -n 1 -N 1 --exclusive %s\n' % (command))
    scriptFile.close();


    ''' Write a submit script for a local machine. '''
    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# local pimc submit script: write your pimc script 
#--------------------------------------------------''')
    # Create the command string and output to submit file
    if (QUE == "slurm"):
        nodes = int(math.ceil(len(range(numOptions))/user_avgcores));

    scriptFile.write('# You need to run %s jobs\n' % len(range(numOptions)))
    if (QUE =="slurm" and user_qos != ""):
         scriptFile.write('#SBATCH --qos=%s\n' % user_qos)

    if (QUE == "slurm" and user_partition != ""):
        scriptFile.write('#SBATCH --partition=%s\n' % user_partition)
    
    if(QUE == "slurm" and nodes > user_Node and user_Node != -1):
        scriptFile.write('#SBATCH -N %s\n#if you do not want all of them to run on the same cores leave blank\n' % user_Node)
    elif(QUE == "slurm" and nodes <= user_Node and user_Node !=-1):
        scriptFile.write('#SBATCH -N %s\n#if you do not want all of them to run on the same cores leave blank\n' % nodes)
    #I have it like this so that if user specificies maximum node size and cores per node then it overides maximum cpu in a partition else it is either maximum cpu or number of total tasks
    if (QUE == "slurm" and user_avgcores != -1 and user_Node != -1 and len(range(numOptions))>(user_avgcores*user_Node)):
        scriptFile.write('#SBATCH -n %s\n' % (user_avgcores*user_Node))
    elif(QUE == "slurm" and len(range(numOptions))<user_cpu):
        scriptFile.write('#SBATCH -n %s\n' % len(range(numOptions)))
    elif(QUE =="slurm" and len(range(numOptions))>user_cpu and user_cpu>1):
        scriptFile.write('#SBATCH -n %s\n' % user_cpu)
    
    if (QUE =="slurm"):
        scriptFile.write("#SBATCH --ntasks-per-core 1\n")


    scriptFile.write('''\n#--------------------------------------------------
######################
# -N     Sets the number of nodes to be allocated to job
# --ntasks-per-core     Sets the maximum number of tasks per allocated core, set this to 1.
# --ntasks-per-node      Sets the maximum number of tasks per allocated core
# --ntasks      Sets the number of tasks to be created for job, set this to the number of tasks for your job (i.e. number of lines of commands in the command file).
# Each user has a limit to number of nodes that can be requested at once use "/sacctmgr -p show qos". The number of maximum tasks that can be made at once is the number of cpu cores on every node.
#####################
\n
# write down the working directory of the pimc code to run the programs
wd= \ncd $wd\n
# Must have GNU PARALLEL to run this, if you do not have parallel then copy and paste the commands file in this file under sbatch commands and use "&"
# Ex: srun -n 1 -N 1 --exclusive ./pimc.e -p 0 -T 0.2 -t 0.002 -L 20.0 -X harmonic -I free -m 48.48 -N 1 -u 0.11 -M 256 --canonical --relaxmu --window=1 --gaussian_window_width=0.5 --relax --no_save_state --bin_size=1000 -E 100000 -S 5000 --estimator=virial --estimator="linear density rho" &
# This should set the command to the background and run them parallel, also use "wait" after all the commands. Otherwise the job will finish without running\n
# write the working directory of the parallel code
pd=
#the commands file needs to be in the working directory that the GNU parallel has access to\n
echo "Starting run at: `date`"\n''')

    scriptFile.write('''# -N 1    This is the number of arguments to pass to each job
# --delay   if running small jobs add a delay, otherwise the controlling node might be overloaded
# -j $SLURM_NTASK     The number of concurrent tasks parallel runs, otherwise the number of cpu allocated
# --joblog name     if you wish to look at the log file of tasks parallel runs add --joblog name, this creates a file in the running directory. Must also add --resume to continue the interrupted run\n\n''')
    
    scriptFile.write('$pd -N 1 -j $SLURM_NTASKS :::: commands-pimc%s.txt\n' % outName)
    
    scriptFile.write('''\n# The commands file holds a list of all the commands that will be ran. DO NOT ADD ANYTHING BESIDES EXECUTABLE LINES OF CODE
# srun is used to create job steps and use to launch processes interactively.
# --exclusive   makes srun use distinct cores for each run, if not used then the tasks may share cores
# -N 1 -n 1     Similar to the bash file headers, this allocates a single core on a single node to each task.\n\n
echo "Finished run at: `date`"\n''')


    scriptFile.close();
    os.system('chmod u+x %s'%fileName)
    print('Run: ./%s'%fileName)


#-----------------------------------------------------------------------------
def local(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
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
def scinet(staticPIMCOps,numOptions,optionValue,walltime,outName,cmd,time=False,
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

    # nasa processer resources
    nasa_cpus = ['broadwell', 'skylake', 'cascade lake', 'haswell', 'ivy bridge', 
                 'sandy bridge', 'broadwell electra']

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='Generate a batch submission\
                                     script from input file')
    parser.add_argument("-c", "--cluster", dest="cluster",\
                        choices=['westgrid','sharcnet','scinet','vacc','local','nasa', 'parallel'],\
            help="target cluster: [westgrid,sharcnet,scinet,vacc,local, parallel]",required=True) 
    parser.add_argument("-t", "--time", action="store_true", dest="time", 
                      help="Prefix the binary with a timing command?") 
    parser.add_argument("-W", "--walltime", dest="walltime", type=str, 
                      help="Specify the walltime in the format HHH:MM:SS.") 
    parser.add_argument("-q", "--queue", dest="queue", choices=['torque','slurm'],
                      help="Batch system: [torque,slurm], if need help with slurm on parallel use slurm else output is skeleton file",default='torque') 
    parser.add_argument("-m", "--model", dest="model", choices=nasa_cpus,
                      help="NASA CPU types.") 
    parser.add_argument("-l", "--label", dest="label", help="An optional label\
                       for the submission script")
    parser.add_argument("--pigs", help="Is this a pigs simulation?",
                        action="store_true")
    parser.add_argument("input", help="the input file")

    # parse the command line options
    args = parser.parse_args() 
    inFileName = args.input

    # is this a pigs or pimc simulation?
    if args.pigs:
        cmd = 'pigs.e'
    else:
        cmd = 'pimc.e'

    # We open up the input file, and read in all lines.
    inFile = open(inFileName,'r')
    inLines = inFile.readlines();

    # Maps cluster names to functions
    gen = {'westgrid':westgrid, 'sharcnet':sharcnet, 'scinet':scinet,
            'vacc':vacc, 'local':local, 'nasa':nasa, 'parallel':parallel}

    # setup the default and user specified walltime
    wtime = {'westgrid':'00:00:00', 'sharcnet':'00:00:00', 'scinet':'00:00:00',
            'vacc':'30:00:00', 'local':'00:00:00', 'nasa':'10:00:00', 'parallel':'00:00:00'}
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
        if input == '-N':
            outName += '-%04d' % float(findInput[n+1])
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

    queue = args.queue
    global QUE
    QUE = args.queue
    if args.cluster == 'nasa' and args.model:
        queue = args.model
    elif args.cluster == 'nasa' and not args.model:
        queue = 'broadwell'

    # Generate the submission script
    gen[args.cluster](staticPIMCOps,numOptions,optionValue,wtime[args.cluster],
                         outName,cmd,time=args.time,queue=queue)


# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
