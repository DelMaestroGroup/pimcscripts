#  Adrian Del Maestro
# 09.05.2018

'''procisf.py

Description:
Process a combined output file for the intermediate scattering function
which splits it into different files based on different q-values.

Usage: procisf.py [options]

Options:
  -h, --help                        Show this help message and exit
  -T <T>, --temperature=<T>         Simulation temperature in Kelvin
  -N <N>, --number_particles=<N>    Number of particles
  -n <n>, --density=<n>             Number density in Angstroms^{-d}
  -t <t>, --imaginary_time_step=<t> Imaginary time step
  -P <P>, --number_time_slices=<P>  Number of time slices
  -u <u>, --chemical_potential=<u>  Chemical potential in Kelvin
  -V <V>, --volume=<V>              Volume in Angstroms^d
  -L <L>, --Lz=<L>                  Length in Angstroms
  -i <PIMCID>, --id=<PIMCID>        A list of PIMC ID numbers to include 
  --canonical                       Are we in the canonical ensemble?
  -D <dir>, --working_directory=<dir>   The directory where we perform the merge.  [default: ]
'''

from __future__ import print_function
from docopt import docopt
import numpy as np
import pimchelp

# -----------------------------------------------------------------------------
def splitISF(pimc,idList,baseDir,cyldir=''):
    '''Extract the results of the intermediate scattering function from
       the combined file and write the new files.'''

    fileNames = pimc.getFileList('isf',idList=idList,cyldir=cyldir)
    if not fileNames:
        return

    n = 0
    for fileName in fileNames:

        # get the number of time slices
        pimcID = pimc.getID(fileName)
        numTimeSlices = int(pimc.params[pimcID]['Number Time Slices'])

        # get information on the q-values from the header
        with open(fileName,'r') as inFile:
            lines = inFile.readlines()
            qvals = lines[1].lstrip('#').rstrip('\n').split()
            tvals = lines[2].lstrip('#').rstrip('\n').split()

        # load all the data
        isf_data = np.loadtxt(fileName)

        # break into pieces corresonding to each q
        isf = {}
        tau = {}
        for i,cq in enumerate(qvals):
            isf[cq] = isf_data[:,numTimeSlices*i:(i+1)*numTimeSlices]
            tau[cq] = tvals[numTimeSlices*i:(i+1)*numTimeSlices]
        
        # add duplicate entry for tau = 1/T
        for i,cq in enumerate(qvals):
            isf[cq] = np.hstack((isf[cq],isf[cq][:,:1]))
            invT = float(tau[cq][1]) + float(tau[cq][-1])
            tau[cq].append('{:12.6e}'.format(invT))

        # write to separate files
        for i,cq in enumerate(qvals):
            outFile = fileName.replace('isf','isfq{:s}'.format(cq))
            header = 'PIMCID: {:d}\n'.format(pimcID)
            header += '{:>14}'.format(tau[cq][0]) 
            for ctau in tau[cq][1:]:
                header += '{:>16}'.format(ctau)
            np.savetxt(outFile, isf[cq], fmt='%16.8e',header=header, delimiter='')

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # parse the command line options
    args = docopt(__doc__)

    canonical = args['--canonical']
    pimcID = args['--id']
    baseDir = args['--working_directory']

    # get the list of files that we will process
    pimchelp.checkEnsemble(canonical)
    dataName = pimchelp.getFileString_doc(args,reduce=False)

    # setup the pimchelp object and get parameters
    pimc = pimchelp.PimcHelp(dataName,canonical,baseDir=baseDir)
    pimc.getSimulationParameters(idList=pimcID)

    # split up the ISF into q-dependent files
    splitISF(pimc,pimcID,baseDir)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

if __name__ == "__main__": 
    main()

