#! /usr/bin/env python3
# reduce-one.py
# Adrian Del Maestro
# 09.03.2009
# 
# Reduce and average results for a single PIMC run based on a single parameter
# which varies.  This could be density for fixed system size, etc.

import os,sys,glob
import numpy as np
import pimchelp
import MCstat
from optparse import OptionParser
from matplotlib.pyplot import *
from collections import defaultdict

# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    if data.ndim > dim:
        numBins  = np.size(data,dim) 
        dataAve  = np.average(data,dim) 
        dataAve2 = np.average(data*data,dim) 
        try:
            bins = MCstat.bin(data) 
            dataErr = np.amax(bins,axis=0)
        except:
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
def getScalarEst(etype,pimc,outName,reduceFlag, skip=0, baseDir=''):
    ''' Return the arrays containing the reduced averaged scalar
        estimators in question.'''

    fileNames = pimc.getFileList(etype)
    headers   = pimchelp.getHeadersFromFile(fileNames[0])

    ave = np.zeros([len(fileNames),len(headers)],float)
    err = np.zeros([len(fileNames),len(headers)],float)
    for i,fname in enumerate(fileNames):
        # Compute the averages and error
        data = np.loadtxt(fname,ndmin=2)[skip:,:]
        ave[i,:],err[i,:] = getStats(data)
    
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

    return headers,ave,err;

# -----------------------------------------------------------------------------
def getVectorEst(etype,pimc,outName,reduceFlag,xlab,ylab, skip=0, baseDir=''):
    ''' Return the arrays consisting of the reduec averaged vector 
        estimators. '''

    fileNames = pimc.getFileList(etype)
    try:
        headers = pimchelp.getHeadersFromFile(fileNames[0],getEstimatorInfo=True)

        # We need to check if any information has been included in the header
        # i.e. is header a list of lists
        estInfo = ''
        if any(isinstance(h, list) for h in headers):
            estInfo = headers[0]
            headers = headers[1]

        numParams = len(fileNames)
        Nx = len(headers)

        x   = np.zeros([numParams,Nx],float)
        ave = np.zeros([numParams,Nx],float)
        err = np.zeros([numParams,Nx],float)
        
        for i,fname in enumerate(fileNames):
            print(fname)

            # Get the estimator data and compute averages
            data = np.loadtxt(fname,ndmin=2)[skip:,:]
            ave[i,:],err[i,:] = getStats(data)

            # get the headers
            x[i,:] = pimchelp.getHeadersFromFile(fname)

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
        return x,ave,err

    except:
        print('Problem Reducing %s files' % etype)
        return 0,0,0


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
def getISFEst(pimc,outName,reduceFlag,xlab,ylab,skip=0,baseDir=''):
    ''' Return the arrays consisting of the reduced averaged intermediate
        scattering function. '''

    try:    
        fileNames = pimc.getFileList('isf')

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
        return 0,0,0

    except:
        print('Problem Reducing %s files' % 'isf')
        return 0,0,0


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # define the mapping between short names and label names 
    shortFlags = ['n','T','N','t','u','V','L','W','D','q']
    parMap = {'n':'Initial Density', 'T':'Temperature', 'N':'Initial Number Particles',
              't':'Imaginary Time Step', 'u':'Chemical Potential', 'V':'Container Volume',
              'L':'Container Length', 'W':'Virial Window', 'M':'Update Length'}

    # setup the command line parser options 
    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-N", "--number-particles", dest="N", type="int",
                      help="number of particles") 
    parser.add_option("-n", "--density", dest="n", type="float",
                      help="number density in Angstroms^{-d}")
    parser.add_option("-t", "--imag-time-step", dest="tau", type="float",
                      help="imaginary time step")
    parser.add_option("-u", "--chemical-potential", dest="mu", type="float",
                      help="chemical potential in Kelvin") 
    parser.add_option("-L", "--Lz", dest="L", type="float",
                      help="Length in Angstroms") 
    parser.add_option("-V", "--volume", dest="V", type="float",
                      help="volume in Angstroms^d") 
    parser.add_option("-r", "--reduce", dest="reduce",
                      choices=['T','N','n','u','t','L','V','W','M'], 
                      help="variable name for reduction [T,N,n,u,t,L,V,W,M]") 
    parser.add_option("--canonical", action="store_true", dest="canonical",
                      help="are we in the canonical ensemble?")
    parser.add_option("-p", "--plot", action="store_true", dest="plot",
                      help="do we want to produce data plots?") 
    parser.add_option("-R", "--radius", dest="R", type="float",
                      help="radius in Angstroms") 
    parser.add_option("-s", "--skip", dest="skip", type="int",
                      help="number of measurements to skip") 
    parser.add_option("-e", "--estimator", dest="estimator", type="str",
                      help="specify a single estimator to reduce") 
    parser.add_option("-i", "--pimcid", dest="pimcid", type="str",
                      help="specify a single pimcid") 
    parser.set_defaults(canonical=False)
    parser.set_defaults(plot=False)
    parser.set_defaults(skip=0)

    # parse the command line options and get the reduce flag
    (options, args) = parser.parse_args() 

    # Determine the working directory
    if args:
        baseDir = args[0]
        if baseDir == '.':
            baseDir = ''
    else:
        baseDir = ''

    skip = options.skip
    
    if (not options.reduce):
        parser.error("need a correct reduce flag (-r,--reduce): [T,N,n,u,t,L,V,W,D]")

    # Check that we are in the correct ensemble
    pimchelp.checkEnsemble(options.canonical)

    dataName,outName = pimchelp.getFileString(options)
    reduceFlag = []
    reduceFlag.append(options.reduce)
    reduceFlag.append(parMap[options.reduce])

    # Create the PIMC analysis helper and fill up the simulation parameters maps
    pimc = pimchelp.PimcHelp(dataName,options.canonical,baseDir=baseDir)
    pimc.getSimulationParameters()

    # Form the full output file name
    if options.R == None:
        outName += '.dat'
    else:
        outName += '-R-%04.1f.dat' % options.R

    # possible types of estimators we may want to reduce
    estList = ['estimator', 'super', 'obdm', 'pair', 'radial', 'number', 
               'radwind', 'radarea', 'planedensity', 'planearea',
               'planewind','virial','linedensity','linepotential','energy','isf']
    estDo = {e:False for e in estList}

    # if we specify a single estimator, only do that one
    if options.estimator:
        estDo[options.estimator] = True
    # otherwise test to see if the file exists
    else:
        for e in estList:
            if pimc.getFileList(e):
                estDo[e] = True
            else:
                estDo[e] = False

    # We first reduce the scalar estimators and output them to disk
    if estDo['estimator']:
        head1,scAve1,scErr1 = getScalarEst('estimator',pimc,outName,reduceFlag,
                                           skip=skip,baseDir=baseDir)
    if estDo['energy']:
        head1,scAve1,scErr1 = getScalarEst('energy',pimc,outName,reduceFlag,
                                           skip=skip,baseDir=baseDir)

    if estDo['virial']:
        head1,scAve1,scErr1 = getScalarEst('virial',pimc,outName,reduceFlag,
                                           skip=skip,baseDir=baseDir)

    if estDo['super']:
        head2,scAve2,scErr2 = getScalarEst('super',pimc,outName,reduceFlag,
                                           skip=skip,baseDir=baseDir)

    # Now we do the normalized one body density matrix
    if estDo['obdm']:
        x1,ave1,err1 = getVectorEst('obdm',pimc,outName,reduceFlag,'r [A]','n(r)',
                                    skip=skip,baseDir=baseDir)

    # Now we do the pair correlation function
    if estDo['pair']:
        x2,ave2,err2 = getVectorEst('pair',pimc,outName,reduceFlag,'r [A]','g(r)',
                                    skip=skip,baseDir=baseDir)

    # The radial Density
    if estDo['radial']:
        x3,ave3,err3 = getVectorEst('radial',pimc,outName,reduceFlag,'r [A]','rho(r)',
                                    skip=skip,baseDir=baseDir)

    # Compute the number distribution function and compressibility if we are in
    # the grand canonical ensemble
    if estDo['number']:
        x4,ave4,err4 = getVectorEst('number',pimc,outName,reduceFlag,'N','P(N)',
                                    skip=skip,baseDir=baseDir)

# I don't know why this isn't working, MCStat is giving me an error, will
    # return to this later. AGD 
        #kappa,kappaErr = getKappa(pimc,outName,reduceFlag)

    # The radially averaged Winding superfluid density
    if estDo['radwind']:
        x5,ave5,err5 = getVectorEst('radwind',pimc,outName,reduceFlag,'r [A]','rho_s(r)',
                                    skip=skip,baseDir=baseDir)

    # The radially averaged area superfliud density
    if estDo['radarea']:
        x6,ave6,err6 = getVectorEst('radarea',pimc,outName,reduceFlag,'r [A]','rho_s(r)',
                                    skip=skip,baseDir=baseDir)

    if estDo['planewind']:
        x7,ave7,err7 = getVectorEst('planewind',pimc,outName,reduceFlag,'n','rho_s(r)',
                                    skip=skip,baseDir=baseDir)

    if estDo['planearea']:
        x8,ave8,err8 = getVectorEst('planearea',pimc,outName,reduceFlag,'n','rho_s(r)',
                                    skip=skip,baseDir=baseDir)

    if estDo['planedensity']:
        x9,ave9,err9 = getVectorEst('planedensity',pimc,outName,reduceFlag,'n','rho(r)',
                                    skip=skip,baseDir=baseDir)

    if estDo['linedensity']:
        x10,ave10,err10 = getVectorEst('linedensity',pimc,outName,reduceFlag,
                                       'r [A]','rho1d(r)',skip=skip,baseDir=baseDir)
    if estDo['linepotential']:
        x11,ave11,err11 = getVectorEst('linepotential',pimc,outName,reduceFlag,
                                       'r [A]','V1d(r)',skip=skip,baseDir=baseDir)
    if estDo['isf']:
        x11,ave11,err11 = getISFEst(pimc,outName,reduceFlag,
                                       'τ [1/K]','F(q,τ)',skip=skip,baseDir=baseDir)

    # Do we show plots?
    if options.plot:

        figNum = 1
        # Get the changing parameter that we are plotting against
        param = []
        for ID in pimc.id:
            param.append(float(pimc.params[ID][reduceFlag[1]]))
        numParams = len(param)
        markers = loadgmt.getMarkerList()
        colors  = loadgmt.getColorList('cw/1','cw1-029',10)

        # -----------------------------------------------------------------------------
        # Plot the averaged data
        # -----------------------------------------------------------------------------
        if estDo['estimator']:

            headLab = ['E/N','K/N','V/N','N', 'diagonal']
            dataCol = []
            for head in headLab:
                n = 0
                for h in head1:
                    if head == h:
                        dataCol.append(n)
                        break
                    n += 1
            yLabelCol = ['Energy / N', 'Kinetic Energy / N', 'Potential Energy / N',\
                    'Number Particles', 'Diagonal Fraction']

        
            # ============================================================================
            # Figure -- Various thermodynamic quantities
            # ============================================================================
            for n in range(len(dataCol)):
                figure(figNum)
        
                errorbar(param, scAve1[:,dataCol[n]], yerr=scErr1[:,dataCol[n]],\
                        color=colors[n],marker=markers[n],markeredgecolor=colors[n],\
                        markersize=8,linestyle='None',capsize=4)
        
                xlabel('%s'%options.reduce)
                ylabel(yLabelCol[n])
                tight_layout()
                figNum += 1
    
        # ============================================================================
        # Figure -- The superfluid density
        # ============================================================================
        if estDo['super']:
            figure(figNum)
        
            errorbar(param, scAve2[:,0], yerr=scErr2[:,0],\
                    color=colors[0],marker=markers[0],markeredgecolor=colors[0],\
                    markersize=8,linestyle='None',capsize=4)
        
            tight_layout()
            xlabel('%s'%options.reduce)
            ylabel('Superfluid Density')
    
        # ============================================================================
        # Figure -- The one body density matrix
        # ============================================================================
        if estDo['obdm']:
            figNum += 1
            figure(figNum)
            ax = subplot(111)
    
            for n in range(numParams):
                lab = '%s = %s' % (options.reduce,param[n])
                errorbar(x1[n,:], (ave1[n,:]+1.0E-15), err1[n,:],color=colors[n],marker=markers[0],\
                        markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
                #axis([0,21,1.0E-5,1.1])
            xlabel('r [Angstroms]')
            ylabel('One Body Density Matrix')
            tight_layout()
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
        # ============================================================================
        # Figure -- The pair correlation function
        # ============================================================================
        if estDo['pair']:
            figNum += 1
            figure(figNum)
        
            for n in range(numParams):
                lab = '%s = %s' % (options.reduce,param[n])
                errorbar(x2[n,:], ave2[n,:], yerr=err2[n,:],color=colors[n],marker=markers[0],\
                        markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab,capsize=6)
        
                #   axis([0,256,1.0E-5,1.2])
            xlabel('r [Angstroms]')
            ylabel('Pair Correlation Function')
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
            tight_layout()
    
        # We only plot the compressibility if we are in the grand-canonical ensemble
        if not options.canonical:
    
            # ============================================================================
            # Figure -- The Number distribution
            # ============================================================================
            if estDo['number']:
                figNum += 1
                figure(figNum)

                # Find which column contains the average number of particles
                for hn,h in enumerate(head1):
                    if h == 'N':
                        break

                for n in range(numParams): 
                    lab = '%s = %s' % (options.reduce,param[n]) 
                    aN = scAve1[n,hn] 
                    errorbar(x4[n,:]-aN, ave4[n,:], err4[n,:],color=colors[n],marker=markers[0],\
                             markeredgecolor=colors[n],\
                             markersize=8,linestyle='None',label=lab,capsize=6) 
        
                axis([-30,30,0.0,1.2])
                xlabel(r'$N-\langle N \rangle$')
                ylabel('P(N)')
                tight_layout()
                legend(loc='best', frameon=False, prop={'size':16},ncol=2)
        
                # ============================================================================
                # Figure -- The Compressibility
                # ============================================================================
                #figNum += 1
                #figure(figNum)

                #errorbar(param, kappa, yerr=kappaErr, color=colors[0],marker=markers[0],\
                #        markeredgecolor=colors[0], markersize=8,linestyle='None',capsize=6)
        
                #tight_layout()
                #xlabel('%s'%options.reduce)
                #ylabel(r'$\rho^2 \kappa$')
    
        # ============================================================================
        # Figure -- The radial density
        # ============================================================================
        if len(glob.glob('CYLINDER')) > 0:
            figNum += 1
            figure(figNum)
            ax = subplot(111)
    
            for n in range(numParams):
                lab = '%s = %s' % (options.reduce,param[n])
                errorbar(x3[n,:], (ave3[n,:]+1.0E-15), err3[n,:],color=colors[n],marker=markers[0],\
                        markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
                #axis([0,21,1.0E-5,1.1])
            tight_layout()
            xlabel('r [Angstroms]')
            ylabel('Radial Density')
            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
        show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

