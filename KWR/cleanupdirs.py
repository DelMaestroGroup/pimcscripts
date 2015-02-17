import os 
import subprocess
import argparse 
import distutils.core

def merge_seeds(pdir,mdir,skip,udirs):
    print "You are now merging seeds"
    data=pdir+'DATA/'
    os.chdir(data)
    # Move into each unique directory, remove the pesky pimc.out file, and call
    # merge.py skipping number of bins provided at the command line. Note that
    # the same number of measurements will be skipped across all parameter sets, 
    # which might be undesireable
    for i in udirs:
        os.chdir(i)
        os.remove('pimc.out')
        cmd='python '+mdir+' -s '+str(skip)
        subprocess.check_output(cmd,shell=True)
        os.chdir(data)
    
def main():
    parser = argparse.ArgumentParser(description='Merges the directories for differing process numbers but identical parameters on OSG')
    parser.add_argument('-p', metavar='Input Path',help='Path to parent directory for run/sweep')
    parser.add_argument('-d', metavar='Parent Directory for DATA',help='Path to the parent directory for the final DATA folder')
    parser.add_argument('--merge',help='Calls merge.py to merge .dat files for different process numbers. Must specify full path to merge.py')
    parser.add_argument('-s','--skip',default=0,help='Number of bins to be skipped by merge.py')
    args = parser.parse_args()
    
    # Move to the parent directory for the run
    os.chdir(args.p)
    # Grab all the job directories
    dirlist = []
    for jobdir in os.listdir(args.p):
        if os.path.isdir(jobdir):
            dirlist.append(jobdir)
    # Remove the process numbers
    nopdirlist = []
    for jobdir in dirlist:
        i=jobdir.find('p') 
        nopjobdir=jobdir[0:i-2]
        nopdirlist.append(nopjobdir)
    # Remove duplicate parameter sets  
    uniqueparamsets = set(nopdirlist)
    # Copy data from jobs with same parameters but different process numbers
    # to a new unifying directory.
    for uniquedir in uniqueparamsets:
        for jobdir in dirlist:
            i=jobdir.find('p') 
            nopjobdir=jobdir[0:i-2]
            if nopjobdir == uniquedir:
                distutils.dir_util.copy_tree(args.p+"/"+jobdir,args.d+"DATA/"+uniquedir)
    # Copy unifying directories to args.d/DATA/
    #for uniquedir in uniqueparamsets:
    #    distutils.dir_util.copy_tree(args.p+"/"+uniquedir,args.d+"DATA/"+uniquedir)
    
    # If specified at command line, merge seeds in DATA
    if args.merge:
        parentdatadir=args.d
        mergedir=args.merge
        skip=args.skip
        merge_seeds(parentdatadir, mergedir,skip,uniqueparamsets)

if __name__=='__main__':
    main()