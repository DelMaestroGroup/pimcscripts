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
from collections import defaultdict
import argparse 


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
def getScalarEst(etype,pimc,outName,reduceFlag,axis_labels,skip=0, baseDir='',idList=None):
    ''' Return the arrays containing the reduced averaged scalar
        estimators in question.'''

    fileNames = pimc.getFileList(etype,idList)

    # we make sure that we have a valid list of filenames
    try:
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

        return headers,ave,err,len(fileNames)
    except:
        print('Problem Reducing %s files' % etype)
        print(' '.join(sys.argv))
        return 0,0,0,0

# -----------------------------------------------------------------------------
def getVectorEst(etype,pimc,outName,reduceFlag,axis_labels,skip=0,baseDir='',
                 idList=None):
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
        Nx = len(headers)

        x   = np.zeros([numParams,Nx],float)
        ave = np.zeros([numParams,Nx],float)
        err = np.zeros([numParams,Nx],float)
        
        for i,fname in enumerate(fileNames):

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
def getISFEst(pimc,outName,reduceFlag,axis_labels,skip=0,baseDir='',idList=None):
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
    parser.add_argument("-s", "--skip", dest="skip", type=int, default=0,
               help="number of measurements to skip") 
    parser.add_argument("-e", "--estimator", dest="estimator", type=str,
                        nargs='*', help="specify a single estimator to reduce") 
    parser.add_argument("-i", "--pimcid", dest="pimcid", type=str, nargs='*',
                      help="specify a single pimcid") 
    parser.add_argument("base_dir", help='The base directory\
                        where the data files to be reduced are located.',
                        default='./', nargs='?')
    # parser.add_argument("-p", "--plot", action="store_true", dest="plot",
    #            help="do we want to produce data plots?") 

    # parse the command line arguments
    args = parser.parse_args()

    # Determine the working directory
    baseDir = os.path.join(args.base_dir,'')

    skip = args.skip
    
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
    if args.pimcid is not None:
        outName += f'-{args.pimcid[0]}'
    outName += '.dat'

    # possible types of estimators we may want to reduce
    est_list = ['estimator', 'super', 'obdm', 'pair', 'radial', 'number', 
               'radwind', 'radarea', 'planedensity', 'planearea',
               'lineardensity', 'planewind','virial','linedensity',
               'linepotential','energy','isf']
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
                   'linepotential':['r [Å]','V1d(r)'], 
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
        est_return = getEst[est](est,pimc,outName,reduceFlag,axis_label[est],
                                 skip=skip, baseDir=baseDir, idList=args.pimcid)
        if est_return[-1]:
            print(f'Reduced {est} over {est_return[-1]} {parMap[args.reduce]} value(s).')


    # Do we show plots?
    # NB: This is old code.  Plots are not taken care of with other scripts. To
    # be removed at some point.
    #if args.plot:

    #    figNum = 1
    #    # Get the changing parameter that we are plotting against
    #    param = []
    #    for ID in pimc.id:
    #        param.append(float(pimc.params[ID][reduceFlag[1]]))
    #    numParams = len(param)
    #    markers = loadgmt.getMarkerList()
    #    colors  = loadgmt.getColorList('cw/1','cw1-029',10)

    #    # -----------------------------------------------------------------------------
    #    # Plot the averaged data
    #    # -----------------------------------------------------------------------------
    #    if est_do['estimator']:

    #        headLab = ['E/N','K/N','V/N','N', 'diagonal']
    #        dataCol = []
    #        for head in headLab:
    #            n = 0
    #            for h in head1:
    #                if head == h:
    #                    dataCol.append(n)
    #                    break
    #                n += 1
    #        yLabelCol = ['Energy / N', 'Kinetic Energy / N', 'Potential Energy / N',\
    #                'Number Particles', 'Diagonal Fraction']

        
    #        # ============================================================================
    #        # Figure -- Various thermodynamic quantities
    #        # ============================================================================
    #        for n in range(len(dataCol)):
    #            figure(figNum)
        
    #            errorbar(param, scAve1[:,dataCol[n]], yerr=scErr1[:,dataCol[n]],\
    #                    color=colors[n],marker=markers[n],markeredgecolor=colors[n],\
    #                    markersize=8,linestyle='None',capsize=4)
        
    #            xlabel('%s'%args.reduce)
    #            ylabel(yLabelCol[n])
    #            tight_layout()
    #            figNum += 1
    
    #    # ============================================================================
    #    # Figure -- The superfluid density
    #    # ============================================================================
    #    if est_do['super']:
    #        figure(figNum)
        
    #        errorbar(param, scAve2[:,0], yerr=scErr2[:,0],\
    #                color=colors[0],marker=markers[0],markeredgecolor=colors[0],\
    #                markersize=8,linestyle='None',capsize=4)
        
    #        tight_layout()
    #        xlabel('%s'%args.reduce)
    #        ylabel('Superfluid Density')
    
    #    # ============================================================================
    #    # Figure -- The one body density matrix
    #    # ============================================================================
    #    if est_do['obdm']:
    #        figNum += 1
    #        figure(figNum)
    #        ax = subplot(111)
    
    #        for n in range(numParams):
    #            lab = '%s = %s' % (args.reduce,param[n])
    #            errorbar(x1[n,:], (ave1[n,:]+1.0E-15), err1[n,:],color=colors[n],marker=markers[0],\
    #                    markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
    #            #axis([0,21,1.0E-5,1.1])
    #        xlabel('r [Angstroms]')
    #        ylabel('One Body Density Matrix')
    #        tight_layout()
    #        legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
    #    # ============================================================================
    #    # Figure -- The pair correlation function
    #    # ============================================================================
    #    if est_do['pair']:
    #        figNum += 1
    #        figure(figNum)
        
    #        for n in range(numParams):
    #            lab = '%s = %s' % (args.reduce,param[n])
    #            errorbar(x2[n,:], ave2[n,:], yerr=err2[n,:],color=colors[n],marker=markers[0],\
    #                    markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab,capsize=6)
        
    #            #   axis([0,256,1.0E-5,1.2])
    #        xlabel('r [Angstroms]')
    #        ylabel('Pair Correlation Function')
    #        legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    #        tight_layout()
    
    #    # We only plot the compressibility if we are in the grand-canonical ensemble
    #    if not args.canonical:
    
    #        # ============================================================================
    #        # Figure -- The Number distribution
    #        # ============================================================================
    #        if est_do['number']:
    #            figNum += 1
    #            figure(figNum)

    #            # Find which column contains the average number of particles
    #            for hn,h in enumerate(head1):
    #                if h == 'N':
    #                    break

    #            for n in range(numParams): 
    #                lab = '%s = %s' % (args.reduce,param[n]) 
    #                aN = scAve1[n,hn] 
    #                errorbar(x4[n,:]-aN, ave4[n,:], err4[n,:],color=colors[n],marker=markers[0],\
    #                         markeredgecolor=colors[n],\
    #                         markersize=8,linestyle='None',label=lab,capsize=6) 
        
    #            axis([-30,30,0.0,1.2])
    #            xlabel(r'$N-\langle N \rangle$')
    #            ylabel('P(N)')
    #            tight_layout()
    #            legend(loc='best', frameon=False, prop={'size':16},ncol=2)
        
    #            # ============================================================================
    #            # Figure -- The Compressibility
    #            # ============================================================================
    #            #figNum += 1
    #            #figure(figNum)

    #            #errorbar(param, kappa, yerr=kappaErr, color=colors[0],marker=markers[0],\
    #            #        markeredgecolor=colors[0], markersize=8,linestyle='None',capsize=6)
        
    #            #tight_layout()
    #            #xlabel('%s'%args.reduce)
    #            #ylabel(r'$\rho^2 \kappa$')
    
    #    # ============================================================================
    #    # Figure -- The radial density
    #    # ============================================================================
    #    if len(glob.glob('CYLINDER')) > 0:
    #        figNum += 1
    #        figure(figNum)
    #        ax = subplot(111)
    
    #        for n in range(numParams):
    #            lab = '%s = %s' % (args.reduce,param[n])
    #            errorbar(x3[n,:], (ave3[n,:]+1.0E-15), err3[n,:],color=colors[n],marker=markers[0],\
    #                    markeredgecolor=colors[n], markersize=8,linestyle='None',label=lab)
    
    #            #axis([0,21,1.0E-5,1.1])
    #        tight_layout()
    #        xlabel('r [Angstroms]')
    #        ylabel('Radial Density')
    #        legend(loc='best', frameon=False, prop={'size':16},ncol=2)
    
    #    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

