#!/usr/bin/env python
# reduce-one.py
# Adrian Del Maestro
# 09.03.2009
# 
# Reduce and average results for a single PIMC run based on a single parameter
# which varies.  This could be density for fixed system size, etc.

import os,sys,glob
import numpy as np
import pimcscripts.pimchelp as pimchelp
import pimcscripts.MCstat as MCstat
from collections import defaultdict
import argparse 
import subprocess

from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

# ----------------------------------------------------------------------
def line_counts(filename):
    '''Use wc to count the number of lines and header lines in a file. '''
    num_lines = int(subprocess.check_output(['wc', '-l', filename]).split()[0])
    num_header = str(subprocess.check_output(['head','-5',filename])).count('#')
    # num_header = str(subprocess.check_output(['grep','-o','-i','\#',filename])).count('#')
    return num_header,num_lines

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if data.ndim > dim:
        numBins  = data.shape[dim]
        dataAve  = np.average(data,dim) 
        try:
            bins = MCstat.bin(data) 
            dataErr = np.amax(bins,axis=0)
        except:
            dataAve2 = np.average(data*data,dim) 
            dataErr = np.sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 

#        for n,d in enumerate(dataErr):
#            if d > 2.0*dataErr2[n]:
#                dataErr[n] = 2.0*dataErr2[n]

#        try:
#            bins = MCstat.bin(data) 
#            dataErr = amax(bins,axis=0)
#        except:
#            dataErr   = sqrt( abs(dataAve2-dataAve**2)/(1.0*numBins-1.0) ) 
    else:
        dataAve = data
        dataErr = 0.0*data

    return dataAve,dataErr

# -----------------------------------------------------------------------------
def process_stats(fname,skip,get_headers=False,ave_est=False):
    '''Get the average and error for a estimator file. '''

    # determine the structure of the headers
    num_headers,num_lines = line_counts(fname)

    # Get the number of lines to skip
    if isinstance(skip, float):
        cskip = int(num_lines*skip)+num_headers
    else:
        cskip = skip+num_headers

    # Pre-averaged estimators require no processing
    if ave_est:
        x = np.arange(0,num_lines-num_headers,1)
        y = np.loadtxt(fname,ndmin=1)
        Δy = np.zeros_like(y)
        return (y,Δy),x

    # Compute the averages and error
    else:
        if get_headers:
            return getStats(np.loadtxt(fname,ndmin=2,skiprows=cskip)),pimchelp.getHeadersFromFile(fname)
        else:
            return getStats(np.loadtxt(fname,ndmin=2,skiprows=cskip))

# -----------------------------------------------------------------------------
def getScalarEst(etype,pimc,outName,reduceFlag,axis_labels,skip=0,
                 baseDir='',idList=None,ave_est=False):
    ''' Return the arrays containing the reduced averaged scalar
        estimators in question.'''

    fileNames = pimc.getFileList(etype,idList)

    # we make sure that we have a valid list of filenames
    try:
        headers = pimchelp.getHeadersFromFile(fileNames[0])

        ave = np.zeros([len(fileNames),len(headers)],float)
        err = np.zeros([len(fileNames),len(headers)],float)

        # process all files in parallel
        nc = min(len(fileNames),num_cores)
        results = Parallel(n_jobs=nc)(delayed(process_stats)(fname,skip) for fname in fileNames)
        for i,result in enumerate(results):
            ave[i,:],err[i,:] = result

        # compute single centroid virial specific heat if possible
        # if 'dEdB' in headers:
        #     Cv = ave[:,headers.index('EEcv*Beta^2')] - ave[:,headers.index('Ecv*Beta')]**2 - ave[:,headers.index('dEdB')]
        #     aveNew = np.zeros([len(fileNames),len(headers)+1],float)
        #     errNew = np.zeros([len(fileNames),len(headers)+1],float)
        #     for i,a in enumerate(ave):
        #         a = append(a, ave[:,headers.index('EEcv*Beta^2')][i] \
        #                 - ave[:,headers.index('Ecv*Beta')][i]**2 \
        #                 - ave[:,headers.index('dEdB')][i])
        #         aveNew[i] = a
        #     for i, e in enumerate(err):
        #         e = append(e, err[:,headers.index('EEcv*Beta^2')][i] \
        #                 - err[:,headers.index('Ecv*Beta')][i]**2 \
        #                 - err[:,headers.index('dEdB')][i])
        #         errNew[i] = e 
        #     headers.append('Cv')
        #     ave = aveNew
        #     err = errNew

        # output the estimator data to disk
        outFile = open(baseDir + '%s-%s' % (etype,outName),'w');

        # the param and data headers
        outFile.write('#{:>15s}'.format(reduceFlag[0]))
        for head in headers:
            outFile.write('{:>16s}{:>16s}'.format(head,'Δ{:s}'.format(head)))
        outFile.write('\n')

        # the data
        for i,f in enumerate(fileNames):
            outFile.write('%16.8E' % float(pimc.params[pimc.id[i]][reduceFlag[1]]))
            for j,h in enumerate(headers):
                outFile.write('%16.8E%16.8E' % (ave[i,j],err[i,j]))
            outFile.write('\n')
        outFile.close()

        return headers,ave,err,len(fileNames)
    except:
        print('Problem Reducing %s files' % etype)
        print(' '.join(sys.argv))
        return 0,0,0,0

# -----------------------------------------------------------------------------
def getVectorEst(etype,pimc,outName,reduceFlag,axis_labels,skip=0,baseDir='',
                 idList=None,ave_est=False):
    ''' Return the arrays consisting of the reduec averaged vector 
        estimators. '''

    xlab,ylab = axis_labels
    fileNames = pimc.getFileList(etype,idList)

    try:
        headers = pimchelp.getHeadersFromFile(fileNames[0],getEstimatorInfo=True)

        # We need to check if any information has been included in the header
        # i.e. is header a list of lists
        estInfo = ''
        if any(isinstance(h, list) for h in headers):
            estInfo = headers[0]
            headers = headers[1]

        numParams = len(fileNames)

        # Do we have a pre-averaged vector estimator? No headers are included.
        if ave_est:
            num_headers,num_lines = line_counts(fileNames[0])
            Nx = num_lines-num_headers
        else:
            Nx = len(headers)

        x   = np.zeros([numParams,Nx],float)
        ave = np.zeros([numParams,Nx],float)
        err = np.zeros([numParams,Nx],float)

        # process all files in parallel
        nc = min(len(fileNames),num_cores)
        results = Parallel(n_jobs=nc)(delayed(process_stats)(fname,skip,get_headers=True,ave_est=ave_est) 
                                   for fname in fileNames)

        # collect the results
        for i,result in enumerate(results):
            ave[i,:],err[i,:] = result[0]
            x[i,:] = result[1]

            # Compute the normalized averages and error for the OBDM
            if etype == 'obdm':
                norm = ave[i,0]
                ave[i,:] /= norm
                err[i,:] /= norm
        
        # the param and data headers
        header_p = ''
        header_d = ''
        for j in range(numParams):
            lab = '%s = %4.2f' % (reduceFlag[0],float(pimc.params[pimc.id[j]][reduceFlag[1]]))
            header_p += '{:^48s}'.format(lab)
            header_d += '{:>16s}{:>16s}{:>16s}'.format(xlab,ylab,'Δ'+ylab)
        header = estInfo + '# ' + header_p[2:] + '\n' + '# ' + header_d[2:]

        # collapse the data
        out_data = [np.vstack((x[i,:],ave[i,:],err[i,:])).T for i in range (numParams)]

        # output the vector data to disk
        outFileName = baseDir+'%s-%s' % (etype,outName)
        np.savetxt(outFileName,np.hstack(out_data),delimiter='',comments='', 
                   header=header,fmt='% 16.8E')
        return x,ave,err,len(fileNames)

    except:
        print('Problem Reducing %s files' % etype)
        print(' '.join(sys.argv))
        return 0,0,0,0

# -----------------------------------------------------------------------------
def getKappa(pimc,outName,reduceFlag,skip=0,baseDir=''):
    ''' Return the arrays containing the reduced averaged compressibility. '''

    fileNames = pimc.getFileList('estimator')
    headers   = pimchelp.getHeadersFromFile(fileNames[0])

    aveKappa = np.zeros([len(fileNames)],float)
    errKappa = np.zeros([len(fileNames)],float)

    for i,fname in enumerate(fileNames):

        # Now get the temperature, volume and linear system size
        ID = pimc.getID(fname)
        T = float(pimc.params[ID]['Temperature'])

        # We need to get the correct volume, depending on whether or not we are
        # looking at the core
        if len(glob.glob('../CYLINDER')) > 0:
            V = pi*(1.75)**2*float(pimc.params[ID]['Container Length'])
        else:
            V = float(pimc.params[ID]['Container Volume'])

        # Compute the average compressibility and its error
        estData = np.loadtxt(fname,ndmin=2)

        N     = estData[:,headers.index('N')]
        N2    = estData[:,headers.index('N^2')] 
        N3 = N*N2

        numBins = len(N)

        # Get the averages
        aveN,errN = getStats(N)
        aveN2,errN2 = getStats(N2)
        aveNN2,errNN2 = getStats(N3)


        # Get the covariance
        # This is finite, so it must be calculated!
        covNN2 = (aveNN2 - aveN*aveN2)/(1.0*numBins-1.0)

        # Get the value of rho^2 * kappa and the error
        aveKappa[i] = (aveN2-aveN**2)/(T*V)
        errKappa[i] = np.sqrt(errN2**2 + 4.0*errN**2*aveN**2 - 4.0*aveN*covNN2)/(T*V)
    
    # output the estimator data to disk
    outFile = open(baseDir+'%s-%s' % ('kappa',outName),'w');

    # the headers
    outFile.write('#%15s' % reduceFlag[0])
    outFile.write('%16s%16s' % ('kappa','+/-'))
    outFile.write('\n')

    # the data
    for i in range(len(fileNames)):
        outFile.write('%16.8E' % float(pimc.params[pimc.id[i]][reduceFlag[1]]))
        outFile.write('%16.8E%16.8E\n' % (aveKappa[i],errKappa[i]))
    outFile.close()

    return aveKappa,errKappa

# -----------------------------------------------------------------------------
def getISFEst(pimc,outName,reduceFlag,axis_labels,skip=0,baseDir='',
              idList=None,ave_est=False):
    ''' Return the arrays consisting of the reduced averaged intermediate
        scattering function. '''

    try:    
        xlab,ylab = axis_labels
        fileNames = pimc.getFileList('isf',idList)

        # get the number of time slices
        pimcID = pimc.getID(fileNames[0])
        numTimeSlices = int(pimc.params[pimcID]['Number Time Slices'])

        # get information on the q-values from the header
        with open(fileNames[0],'r') as inFile:
            lines = inFile.readlines()
            qvals = lines[1].lstrip('#').rstrip('\n').split()
            tvals = lines[2].lstrip('#').rstrip('\n').split()

        numParams = len(fileNames)

        Nq = len(qvals)
        Nx = int(len(tvals)/Nq) + 1

        ave = np.zeros([Nq,numParams,Nx],float)
        err = np.zeros([Nq,numParams,Nx],float)

        for i,fname in enumerate(fileNames):

            # load all the data
            isf_data = np.loadtxt(fname)

            # break into pieces corresonding to each q
            isf = {}
            τ = {}
            for j,cq in enumerate(qvals):
                isf[cq] = isf_data[skip:,numTimeSlices*j:(j+1)*numTimeSlices]
                τ[cq] = np.array([float(cτ) for cτ in tvals[numTimeSlices*j:(j+1)*numTimeSlices]])
            
            # add duplicate entry for tau = 1/T and get the averages and error
            for j,cq in enumerate(qvals):
                isf[cq] = np.hstack((isf[cq],isf[cq][:,:1]))
                τ[cq] = np.append(τ[cq],τ[cq][1] + τ[cq][-1])

                # Get the estimator data and compute averages
                ave[j,i,:],err[j,i,:] = getStats(isf[cq])

        # output the vector data to disk
        outFileName = baseDir+'%s-%s' % ('isf',outName)

        # the param and data headers
        header_q = ''
        header_p = ''
        header_d = ''
        for cq in qvals:
            header_q += '{:^48}'.format('q = {:s}'.format(cq))
            for j in range(numParams):
                lab = '%s = %4.2f' % (reduceFlag[0],float(pimc.params[pimc.id[j]][reduceFlag[1]]))
                header_p += '{:^48s}'.format(lab)
                header_d += '{:>16s}{:>16s}{:>16s}'.format(xlab,ylab,'Δ'+ylab)

        header = '# ' + header_q[2:] + '\n' + '# ' + header_p[2:] + '\n' + '# ' + header_d[2:]

        # collapse the data
        out_data = []
        for i,cq in enumerate(qvals):
            for j in range(numParams):
                out_data.append(np.vstack((τ[cq],ave[i,j,:],err[i,j,:])).T)

        out_data = np.hstack(out_data)
        np.savetxt(outFileName,out_data,delimiter='',comments='', header=header,fmt='% 16.8E')
        return 0,0,0,len(fileNames)

    except:
        print('Problem Reducing %s files' % 'isf')
        print(' '.join(sys.argv))
        return 0,0,0,0


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # define the mapping between short names and label names 
    shortFlags = ['n','T','N','t','u','V','L','D','q']
    parMap = {'n':'Initial Density', 'T':'Temperature', 'N':'Initial Number Particles',
              't':'Imaginary Time Step', 'u':'Chemical Potential', 
              'L':'Container Length', 'V':'Virial Window', 'M':'Update Length'}

    # setup the command line parser arguments
    parser = argparse.ArgumentParser(description='Reduce quantum Monte Carlo\
                                     output over some parameter.')
    parser.add_argument("-T", "--temperature", dest="T", type=float,
             help="simulation temperature in Kelvin") 
    parser.add_argument("-N", "--number-particles", dest="N", type=int,
               help="number of particles") 
    parser.add_argument("-n", "--density", dest="n", type=float,
               help="number density in Angstroms^{-d}")
    parser.add_argument("-t", "--imag-time-step", dest="tau", type=float,
               help="imaginary time step")
    parser.add_argument("-u", "--chemical-potential", dest="mu", type=float,
               help="chemical potential in Kelvin") 
    parser.add_argument("-L", "--Lz", dest="L", type=float,
               help="Length in Angstroms") 
    parser.add_argument("-V", "--volume", dest="V", type=float,
               help="volume in Angstroms^d") 
    parser.add_argument("-r", "--reduce", dest="reduce", required=True,
               choices=['T','N','n','u','t','L','V','M'], 
               help="variable name for reduction [T,N,n,u,t,L,V,W,M]") 
    parser.add_argument("--canonical", action="store_true", dest="canonical",
               help="are we in the canonical ensemble?")
    parser.add_argument("-R", "--radius", dest="R", type=float,
               help="radius in Angstroms") 
    parser.add_argument("-s", "--skip", help="number of measurements to skip [0]") 
    parser.add_argument("-e", "--estimator", dest="estimator", type=str,
                        action='append', help="specify a single estimator to reduce") 
    parser.add_argument("-i", "--pimcid", dest="pimcid", type=str, 
                        help="specify a single pimcid", action='append') 
    parser.add_argument("base_dir", help='The base directory\
                        where the data files to be reduced are located.',
                        default=os.getcwd(), nargs='?')
    # parser.add_argument("-p", "--plot", action="store_true", dest="plot",
    #            help="do we want to produce data plots?") 

    # parse the command line arguments
    args = parser.parse_args()

    # Determine the working directory
    baseDir = os.path.join(args.base_dir,'')

    if not args.skip:
        skip = 0
    else:
        if '.' in args.skip:
            skip = float(args.skip)
            if skip < 0.0 or skip >= 1.0:
                raise ValueError('skip < 0.0 or skip >= 1.0')
        else:
            skip = int(args.skip)

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(args.canonical)

    dataName,outName = pimchelp.getFileString(args)

    reduceFlag = []
    reduceFlag.append(args.reduce)
    reduceFlag.append(parMap[args.reduce])

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,args.canonical,baseDir=baseDir)
    pimc.getSimulationParameters()

    # Form the full output file name
    if args.R is not None:
        outName += f'-R-{args.R:04.1f}'
    # if args.pimcid is not None:
    #     outName += f'-{args.pimcid}'
    outName += '.dat'

    # possible types of estimators we may want to reduce
    est_list = pimc.dataType
    est_do = {e:False for e in est_list}

    scalar_est_list = ['estimator','energy','virial','super']
    vector_est_list = [est for est in est_list if est not in scalar_est_list]

    # setup the dispatch
    getEst = {}
    for est in est_list:
        if est in scalar_est_list:
            getEst[est] = getScalarEst
        elif est == 'isf':
            getEst[est] = getISFEst
        else:
            getEst[est] = getVectorEst

    # setup the x- and y-labels
    axis_label = {'obdm':['r [Å]','n(r)'], 'pair':['r [Å]','g(r)'], 
                   'radial':['r [Å]','ρ(r)'], 'number':['N','P(N)'],
                   'radwind':['r [Å]','ρₛ(r)'],'radarea':['r [Å]','ρₛ(r)'],
                   'planewind':['n','ρₛ(x,y)'],'planearea':['n','ρₛ(x,y)'],
                   'planedensity':['n','ρ(x,y)'], 'linedensity':['r [Å]','ρ1d(r)'],
                   'linepotential':['r [Å]','V1d(r)'], 'ssf':['q [1/Å]', 'S(q)'],
                   'ssfq':['q_index', 'S(q)'], 'planeavedensity':['n','ρ(x,y)'],
                  'lineardensity':['r [Å]','ρ(r)'], 'isf':['τ [1/K]','F(q,τ)']}

    # there are no labels for scalar estimators
    for est in scalar_est_list:
        axis_label[est] = ''

    # if we specify a single estimator, only do that one
    est_do = []
    if args.estimator is not None:
        for est in args.estimator:
            est_do.append(est)
    else:
    # otherwise do any that exist
        for est in est_list:
            if pimc.getFileList(est):
                est_do.append(est)

    # perform the reduction
    for est in est_do:
        ave_est = 'ave' in est
        est_return = getEst[est](est,pimc,outName,reduceFlag,axis_label[est],
                                 skip=skip, baseDir=baseDir,
                                 idList=args.pimcid,ave_est=ave_est)
        if est_return[-1]:
            print(f'Reduced {est} over {est_return[-1]} {parMap[args.reduce]} value(s).')

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

