#!/usr/bin/python 
#
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

import os,sys,itertools
from optparse import OptionParser

# -----------------------------------------------------------------------------
def vacc(staticPIMCOps,numOptions,optionValue,outName,processbool):
    ''' Write a pbs submit script for westgrid. '''

    # Open the pbs file and write its header
    fileName = 'submit-pimc%s.pbs' % outName
    pbsFile = open(fileName,'w')
    pbsFile.write('''#!/bin/bash
#PBS -S /bin/bash\n
#PBS -l walltime=29:00:00
#PBS -l nodes=1:ppn=1,pmem=1gb
#PBS -N pimc
#PBS -V
#PBS -j oe
#PBS -o out/pimc-${PBS_ARRAYID}-${PBS_JOBID}\n 
# Do not send email
#PBS -M kyle.robertson@uvm.edu
#PBS -m n\n
# Start job script
cd $PBS_O_WORKDIR
mkdir out
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')
    valuelist = optionValue.values()
    combos=list(itertools.product(*valuelist))
    keys = optionValue.keys()
    # Create the command string and make the case structure
    n = 0
    for combo in combos:
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1 or processbool):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
        cmdstring = ''
        for i in range(len(combo)):
            substring ='-'+keys[i]+' '+str(combo[i])+' '
            cmdstring+=substring
        command += cmdstring
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
        n += 1
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print 'qsub -t 0-%d %s\n' % (len(combos)-1,fileName)

# -----------------------------------------------------------------------------
def westgrid(staticPIMCOps,numOptions,optionValue,outName,processbool):
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
mkdir out
echo \"Starting run at: `date`\"

case ${PBS_ARRAYID} in\n''')
    valuelist = optionValue.values()
    combos=list(itertools.product(*valuelist))
    keys = optionValue.keys()
    # Create the command string and make the case structure
    n = 0
    for combo in combos:
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1 or processbool):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
        cmdstring = ''
        for i in range(len(combo)):
            substring ='-'+keys[i]+' '+str(combo[i])+' '
            cmdstring+=substring
        command += cmdstring
        command += staticPIMCOps
        pbsFile.write('%d)\nsleep %d\n%s\n;;\n' % (n,2*n,command))
        n += 1
    pbsFile.write('esac\necho \"Finished run at: `date`\"')
    pbsFile.close();
    
    print 'qsub -t 0-%d %s\n' % (len(combos)-1,fileName)

# -----------------------------------------------------------------------------
def sharcnet(staticPIMCOps,numOptions,optionValue,outName,processbool):
    ''' Write a submit script for sharcnet. '''

    # Open the script file and write its header
    fileName = 'submit-pimc%s' % outName
    scriptFile = open(fileName,'w')
    scriptFile.write('''#!/bin/bash
# Sharcnet pimc submit script\n\n''')

    # Create the command string and output to submit file
    valuelist = optionValue.values()
    combos=list(itertools.product(*valuelist))
    keys = optionValue.keys()
    # Create the command string and make the case structure
    n = 0
    for combo in combos:
        name = '''out/pimc-%J'''
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1 or processbool):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
        cmdstring = ''
        for i in range(len(combo)):
            substring ='-'+keys[i]+' '+str(combo[i])+' '
            cmdstring+=substring
        command += cmdstring
        command += staticPIMCOps
        scriptFile.write('sleep %d\nsqsub -q serial -o %s --mpp=1G -r 6d %s\n' % (10,name,command))
        n += 1
    scriptFile.close();
    os.system('chmod u+x %s'%fileName)

# -----------------------------------------------------------------------------
def scinet(staticPIMCOps,numOptions,optionValue,outName,processbool):
    ''' Write a pbs submit script for scinet. '''

    valuelist = optionValue.values()
    combos=list(itertools.product(*valuelist))
    keys = optionValue.keys()
    print 'Number of parameter combinations = %d' % len(combos)
    if len(combos) % 8 != 0:
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
mkdir out
echo \"Starting run at: `date`\"\n\n''')

    n = 0
    for combo in combos:
        if (optionValue.has_key('p') or staticPIMCOps.find('-p') != -1 or processbool):
            command = './pimc.e '
        else:
            command = './pimc.e -p %d ' % (n)
        cmdstring = ''
        for i in range(len(combo)):
            substring ='-'+keys[i]+' '+str(combo[i])+' '
            cmdstring+=substring
        command += cmdstring
        command += staticPIMCOps
        pbsFile.write('(sleep %02d; %s) &\n' % (10,command))
        n += 1
    pbsFile.write('wait')
    pbsFile.close();
    print '\nSubmit job with: qsub %s\n'%fileName

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-c", "--cluster", dest="cluster", choices=['westgrid','sharcnet','scinet','vacc'],\
            help="target cluster: [westgrid,sharcnet,scinet,vacc]") 
    parser.add_option("--no_process",dest='noprocess',action='store_true',default=False,help="Do not add process number to command string")

    # parse the command line options
    (options, args) = parser.parse_args() 
    inFileName = args[0]

    if (not options.cluster):
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
    
    processbool = options.noprocess
    if options.cluster == 'westgrid':
        westgrid(staticPIMCOps,numOptions,optionValue,outName, processbool)

    if options.cluster == 'sharcnet':
        sharcnet(staticPIMCOps,numOptions,optionValue,outName, processbool)

    if options.cluster == 'scinet':
        scinet(staticPIMCOps,numOptions,optionValue,outName, processbool)

    if options.cluster == 'vacc':
        vacc(staticPIMCOps,numOptions,optionValue,outName, processbool)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
