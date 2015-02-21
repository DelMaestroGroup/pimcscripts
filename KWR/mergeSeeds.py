import argparse
import os
import pimchelper
import subprocess
import glob
import numpy as np

def sh(cmd):
    '''A useful function that allows direct access to the bash shell. Simply type what you usually would at 
    the terminal into this function and you can interact with the bash shell'''
    return subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]

def prepFile(mergedir,prefix,estname,dataName,newID,samplefile):
    filename = "%s-%s-%s-%d.dat"%(prefix,estname,dataName,newID)
    print "New file name is: "%filename

    if not (mergedir+filename in os.listdir(mergedir)):
        print "Merge file has not been created, prepping now"
        # Open up the merge file for appending
        mergefile = open(mergedir+filename,'a')
        # Get the necessary headers from the sample file
        print samplefile
        headers = pimchelper.getHeadersFromFile(samplefile)
        headerstring = '#'
        for label in headers:
            headerstring+='%15s'%label
        headerstring+='%15s'%"bins"
        # Write PIMCID and header to file
        mergefile.write('# PIMCID: %d\n'%newID)
        mergefile.write(headerstring+'\n')
        return mergefile
    else:
        print "Merge file has been made, do not prep file"
        mergefile=open(filename,'a')
        return mergefile
    
def averageScalarEstimator(estFile):
    '''Takes a scalar estimator file and returns the average and standard error of 
    all its measurements'''
    
    print "Averaging a scalar estimator seed ..."
    # Load in the data, defaults to floats
    data = np.loadtxt(estFile)
    # Crunch data
    aves = np.average(data,axis=0)
    stdevs = np.std(data,axis=0)
    bins=data.shape[0]
    return aves, stdevs, bins

def averageVectorEstimator(estFile):
    '''Takes a vector estimator file and returns the average and standard error of 
    all its measurements'''
    print "Averaging a vector estimator seed"
    # Load in the data, defaults to floats
    data = np.loadtxt(estFile)
    # Crunch data
    aves = np.average(data,axis=0)
    stdevs = np.std(data,axis=0)
    bins=data.shape[0]
    return aves, stdevs, bins

def main(): 
    parser = argparse.ArgumentParser(description='''Averages data from each random number seed in a set of seeds and merges the averages into one file for each estimator. The final 
    merged files exist in a directory MERGED in the same directory as the seed files.''')
    parser.add_argument('-d','--directory',dest='d',type=str,help='Parent directory where subdirectories containing data files of identical physical parameters but different seeds to be merged exist')
    #parser.add_argument('-s','--skip',dest='s',type=int,help='Number of lines to be skipped when averaging each individual seed')
    # Need to find a way of skipping a different number of bins for each subdirectory. Need to associate a number to skip with a set of parameters in a simple way. Best way
    # might be user input for each subdirectory 
    args=parser.parse_args()
    
    parentdir=args.d
    
    n=0
    for root, dirs, files in os.walk(parentdir):
        print n
        # Make the MERGED directory immediately beneath the parent directory that
        # will be the ultimate destination for all merged estimator files
        if n == 0:
            mergedir = root+'MERGED/'
            try:
                os.mkdir(mergedir)
            except OSError:
                print "MERGED directory has already been created for this parent directory ..."
                msg = "Would you like to re-merge the data? WARNING: This will delete MERGED " \
                "and everything that exists in MERGED"
                print msg
                ansr = raw_input("(y or n): ")
                while not ((ansr == 'y') or (ansr == 'n')):
                    ansr = raw_input("(y or n): ")
                if ansr == 'n':
                    print "Okay, quitting merge process now"
                    quit()
                else:
                    sh("rm -r %s"%mergedir)
                    os.mkdir(mergedir)
                
        # Begin averaging/merging process for each sub directory in the list
        # of subdirectories
        else:
            # Add the trailing slash to subdirectory name to prevent future 
            # annoyances
            root += '/'
            print "Subdirectory : %s"%root
            print "MERGED directory : %s"%mergedir
            # If there are no files at all, skip immediately
            if not files:
                print "No regular files in subdirectory, skipping it"
                n+=1
                continue
            # Create a pimchelper class that contains all the information about
            # this particular set of dat files
            pimc = pimchelper.PIMCHelp(files,root)
            # If there are not any dat files in the directory, skip it
            if not pimc.datfiles:
                print "No data files in subdirectory, skipping it"
                n+=1
                continue
            # Create a new unique PIMCID based on all the IDs in this set of dat files by averaging them.
            newID = 0
            for ID in pimc.IDs:
                newID += int(ID)
            newID = int(newID/(1.0*len(pimc.IDs)))
            # Now we keep incrementing the ID number until we are sure it is different from all surrounding IDs
            while ( (len(glob.glob(parentdir + '*/*estimator*-%09d*' % newID)) > 0 ) or
                (len(glob.glob(parentdir + '*/MERGED/*estimator*-%09d*' % newID)) > 0) ):
                newID += 1
            print "New PIMCID for MERGED set of files = %d"%newID
            # For every estimator in the list of estimators found in this sub directory
            for estimator in pimc.estimators:
                print "Assessing %s files"%estimator
                # Get all the dat files corresponding to this particular estimator
                estFiles = pimchelper.filterfilesbytype(estimator,pimc.datfiles)
                # Identify the type of estimator so that the appropriate averaging
                # procedure can be performed
                estType = pimchelper.getEstimatorType(estimator)
                if estType == 'Scalar':
                    # Prep the file in MERGED for receiving data and get the file object
                    mergefile=prepFile(mergedir,pimc.prefix,estimator,pimc.dataName,newID,root+estFiles[0])
                    # Get the average from each file
                    for seedfile in estFiles: 
                        avrgs, errs, numbins = averageScalarEstimator(root+seedfile)
                        #datarow=' '
                        #for avrg in avrgs:
                        #    datarow+='%15E'%avrg
                        #datarow+='%15d\n'%int(numbins) 
                        datarow = ' '+''.join('%15E'%avrg for avrg in avrgs)+'%15d\n'%numbins                   
                        mergefile.write(datarow)
                    mergefile.close()
                elif estType == 'Vector':
                    mergefile=prepFile(mergedir,pimc.prefix,estimator,pimc.dataName,newID,root+estFiles[0])
                    for seedfile in estFiles:
                        avrgs,errs,bins = averageVectorEstimator(root+seedfile)
                        datarow = ' '+''.join('%15E'%avrg for avrg in avrgs)+'%15d\n'%numbins 
                        #datarow=' '
                        #for avrg in avrgs:
                        #    datarow+='%15E'%avrg
                        #datarow+='%15d\n'%int(numbins)                    
                        mergefile.write(datarow)
                    mergefile.close()
                else:
                    None
        n+=1   
    print "Files merged successfully =)"
        
if __name__=='__main__':
    main()