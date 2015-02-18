import os 
import itertools
import argparse
import subprocess
import re
import numpy as np

def sh(cmd):
'''A useful function that allows direct access to the bash shell. Simply type what you usually
    would at the terminal into this function and you can interact with the bash shell'''
    return subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]

def makeInitialSubmit(staticOpts, opts, keys, pids):
    # Open the pbs file and write its header
    
    fileName = 'submit-pimc.pbs'
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=29:00:00
#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID} 
# Do not send email
#PBS -M kyle.robertson@uvm.edu
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
mkdir out
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')

# Create the command string and make the case structure
    n=0
    for pid in pids:
        command = './pimc.e -p %d '%pid
        cmdstring = ''
        for i in range(len(opts)):
            substring ='-'+keys[i]+' '+str(opts[i])+' '
            cmdstring+=substring
        command += cmdstring
        command += staticOpts
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
        n += 1
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    subcmd='qsub -t 0-%d %s\n' % (len(pids)-1,fileName)
    
    return fileName, subcmd

def makeJobChain(initsubcmd,numrestarts,dependentscript):
    # Open up the shell script that actually submits the initial job and
    # all dependent jobs in the chain to the queue
    fileName='jobchain.sh'
    initsubcmd=initsubcmd.rstrip('\n')
    chainfile=open(fileName,'w')
    # Write the bash header and the first submission command
    chainfile.write("#!/bin/bash\n")
    chainfile.write("one=$(%s)\n"%initsubcmd)
    chainfile.write("echo $one\n")
    # Modify the restart loop to send out the requested number of restarts
    chainfile.write("for num in `seq 1 %d`; do\n"%(numrestarts))
    chainfile.write("  two=$(qsub -W depend=afterokarray:$one %s)\n"%dependentscript)
    chainfile.write("  one=$two\ndone")
    chainfile.close()
    
    return fileName
    
def modifyRestartScript(workdir, numjobs):
    # Copy generic restart script to work dir
    sh("cp RestartPIMC.pbs %s/OldRestartPIMC.pbs"%workdir)
    # Move into work dir
    os.chdir(workdir)
    # Open up files necessary for creation of specific restart script
    oldname="OldRestartPIMC.pbs"
    oldrestartscript=open("OldRestartPIMC.pbs",'r')
    # If the restart script already exists from a previous run remove it and open a fresh copy
    sh("rm RestartPIMC.pbs")  
    newrestartscript=open("RestartPIMC.pbs",'a')
    # Write the PBS directives to the new, specific restart script
    newrestartscript.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=29:00:00
#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V 
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID} 
# Do not send email
#PBS -M kyle.robertson@uvm.edu
#PBS -m n\n''')
    # Write the array directive for the user specified number of process nums
    jobs = "-t 0-%d" % numjobs
    newrestartscript.write("#PBS %s\n\n"%jobs)
    # Add all the generic commands from the generic restart script
    for line in oldrestartscript:
        newrestartscript.write(line)
    # Remove the old restart script and move back to parent directory
    sh("rm %s"%oldname)
    os.chdir("..")

def InitiateJobChain(inFileName,numrestarts):
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
    # If the user chose to sweep through process numbers, remove this from the
    # dictionary so we don't have a separate directory for each process number
    if optionValue.has_key('p'):
        processnums=optionValue['p']
        del optionValue['p']
    # If the user did not sweep through process numbers make a dummy list to use
    # when generating the initial submit file
    else:
        processnums = [1]
        
    # Create a list of tuples. Each tuple contains a unique combination of the
    # dynamic parameters swept through in the input file. The list contains all
    # possible unique combinations. 
    valuelist = optionValue.values()
    combos=list(itertools.product(*valuelist))
    keys = optionValue.keys()
    # Loop through each unique combination of parameters
    for combo in combos:
        # Make a unique working directory based on parameter combo
        workdir=''
        for i in range(len(combo)):
            workdir += str(keys[i])+str(combo[i])+'__'
        workdir=workdir.rstrip('__')
        # If workdir is an empty string because the user only swept through 
        # process numbers, make a working directory so everything else doesn't
        # break.
        if not workdir:
            workdir="WorkingDirectory"
        os.makedirs(workdir)
        # Copy the pimc.e executable to the work dir
        sh("cp pimc.e %s"%workdir)
        # Make the initial submit file that will send out the first job array in 
        # the chain. Every job in the array will have the same physical parameters
        # but a different process number. If the user did not sweep through 
        # process numbers, the array will contain only one job. 
        initsubmitfile, qsubcmd = makeInitialSubmit(staticPIMCOps,combo,keys,processnums)
        # Move the submit file to its apppropriate work dir
        sh("mv %s %s"%(initsubmitfile,workdir))
        # Create the shell script that will send out the job chain, make it 
        # executable, and move it to its work dir
        chainFile=makeJobChain(qsubcmd,numrestarts,"RestartPIMC.pbs")
        sh("chmod u+x %s"%chainFile)
        sh("mv %s %s"%(chainFile, workdir))
        # Modify the generic restart script so it works for this particular 
        # choice of process numbers
        modifyRestartScript(workdir,len(processnums)-1)
        # Move into work dir, execute shell script that submits job chain to the queue
        # then move back to parent dir 
        os.chdir(workdir)
        os.makedirs("old_stdout_files")
        print "Queueing job chain for %s"%workdir
        sh("./%s"%chainFile)
        os.chdir("..")  
    
def ResumeJobChain(numrestarts):
    parentdir=os.getcwd()
    # For working every directory in the parent directory
    for workdir in os.listdir(parentdir):
        if os.path.isdir(workdir):
            commands=[]
            # Move intothe OUTPUT folder of each working directory and grab the 
            # restart command from the log file to add to the commands list
            os.chdir(workdir+"/OUTPUT/")
            for datfile in os.listdir(os.getcwd()):
                if "gce-log" in datfile:
                    logfile=open(datfile,'r')
                    loglines=logfile.readlines()
                    command=loglines[2]
                    command=command.lstrip('# ')
                    command=command.rstrip(' \n')
                    commands.append(command)
            os.chdir("..")
            
            # Open the resubmit pbs file and write its header
            fileName = 'resubmit-pimc.pbs'
            pbsFile = open(fileName,'w')
            pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=29:00:00
#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID} 
# Do not send email
#PBS -M kyle.robertson@uvm.edu
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
mkdir out
mv out/* old_stdout_files/
echo \"Starting run at: `date`\"
        
case ${PBS_ARRAYID} in\n''')
        
            # Create the command string and make the case structure
            n=0
            for command in commands:
                pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
                n += 1
            pbsFile.write('esac\necho \"Finished run at: `date`\"')
            pbsFile.close();
            
            subcmd='qsub -t 0-%d %s\n' % (len(commands)-1,fileName)
            os.chdir(parentdir)
            modifyRestartScript(workdir,len(commands)-1)
            os.chdir(workdir)
            chainFile=makeJobChain(subcmd,numrestarts,"RestartPIMC.pbs")
            sh("chmod u+x %s"%chainFile)
            
            print "Resuming job chain for %s"%workdir
            sh("./%s"%chainFile)
            
            os.chdir(parentdir)
    return fileName, subcmd   
        
def InitiateFromState(numrestarts,binrequest,eqrequest,pids):
    parentdir=os.getcwd()
    # For working every directory in the parent directory
    for workdir in os.listdir(parentdir):
        if os.path.isdir(workdir):
            commands=[]
            # Move into the OUTPUT directory of each working directory
            os.chdir(workdir+"/OUTPUT/")
            for datfile in os.listdir(os.getcwd()):
                if "gce-log" in datfile:
                    # Grab the command from the logfile
                    logfile=open(datfile,'r')
                    loglines=logfile.readlines()
                    command=loglines[2]
                    logfile.close()
                    command=command.lstrip('# ')
                    command=command.rstrip(' \n')
                    # Strip off the restart specific stuff
                    i=command.find("-R")
                    command=command[:i-1]
                    # Remove the process number
                    pmatch=re.search('-p ([0-9]+)',command)
                    command=command.replace(pmatch.group(0),'')
                    # Grab the statefile corresponding to the above log file
                    statefile='OUTPUT/gce-state'+datfile[7:]
                    # Add the submit from state flag
                    command=command+' -s '+statefile
                    # Modify the -S and -E flags using regular expressions to 
                    # reflect the values input by the user 
                    binmatch=re.search('-S ([0-9]+)',command)
                    bins=binmatch.group(0)
                    command=command.replace(bins,'-S %d'%binrequest)
                    eqmatch=re.search('-E ([0-9]+)',command)
                    eqsteps=eqmatch.group(0)
                    command=command.replace(eqsteps,'-E %d'%eqrequest)
                    commands.append(command)
            print commands
            # Move back to working directory
            os.chdir("..")
            
            # Open the resubmit pbs file and write its header
            fileName = 'resubmitfromstate-pimc.pbs'
            pbsFile = open(fileName,'w')
            pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=29:00:00
#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_JOBID} 
# Do not send email
#PBS -M kyle.robertson@uvm.edu
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
mkdir out
mv out/* old_stdout_files/
echo \"Starting run at: `date`\"
        
case ${PBS_ARRAYID} in\n''')
        
            # Create the command string and make the case structure
            n=0
            # For every command pulled from the log files for this working directory
            for command in commands:
                # Make a new command for each PID specified by the user
                for i in range(1, pids+1):
                    newcommand=command+" -p %d"%i
                    print newcommand
                    pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,newcommand))
                    n += 1
            pbsFile.write('esac\necho \"Finished run at: `date`\"')
            pbsFile.close();
            
            # Make qsub command
            numjobs=len(commands)*pids-1 # Subtract 1 because jobs are indexed from 0
            subcmd='qsub -t 0-%d %s\n' % (numjobs,fileName)
            os.chdir(parentdir)
            # Modify the restart script so it works for this set of jobs
            modifyRestartScript(workdir,numjobs)
            os.chdir(workdir)
            # Make the shell script responsible for submitting jobs to the queue
            chainFile=makeJobChain(subcmd,numrestarts,"RestartPIMC.pbs")
            sh("chmod u+x %s"%chainFile)
            # Send out the job chain
            print "Resuming job chain for %s"%workdir
            sh("./%s"%chainFile)
            
            os.chdir(parentdir)
    return fileName, subcmd   
                
    
def main():
    parser = argparse.ArgumentParser(description='Submits a full job chain for an arbitrary number of jobs on the VACC')
    parser.add_argument('-f', metavar='File Path',help='Path to input file')
    parser.add_argument('-R', metavar='Number of restarts',type=int,help='An integer specifying how many restarts you would like in the job chain')
    parser.add_argument('--resume',action='store_true',help='Continues a job chain from where it left off')
    parser.add_argument('--from_state',action='store_true',help="Submits a job chain from a set of state files that exist in the proper directory structure")
    parser.add_argument('-S',metavar='Number of bins desired',type=int,help='Number of bins for simulation. Only used when spawning from state')
    parser.add_argument('-E',metavar='Number of equilibration steps',type=int,help='Number of equilibration steps for simulation. Only used when spawning from state')
    parser.add_argument('-p',metavar='Process numbers',type=int,default=1,help='An integer specifying the amount of process numbers you would like to spawn from each state file')
    args = parser.parse_args()
    
    inFileName = args.f
    numrestarts=args.R
    bins=args.S
    eqsteps=args.E
    pidspawns=args.p
    
    # Catch incorrect usage 
    if args.from_state:
        if not args.E or not args.S:
            print "Must specify -E and -S flags when submitting from state"
            quit() 
            
    if not args.f and not args.resume and not args.from_state:
        print "Must specify an input file when first initiating a run" 
        quit()
   
    if numrestarts == 0:
        pass
    elif not numrestarts:
        print "Need to specify the number of restarts"
        quit()
        
    if args.resume:
        ResumeJobChain(numrestarts)
    elif args.from_state:
        InitiateFromState(numrestarts,bins,eqsteps,pidspawns)
    else:
        InitiateJobChain(inFileName,numrestarts)
        
if __name__=='__main__':       
    main()
