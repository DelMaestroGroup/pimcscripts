import argparse
import os
import subprocess 
import numpy as np
import pylab as pl

def main():
    parser = argparse.ArgumentParser(description='Plots multiple planar density estimator files')
    parser.add_argument('-p', metavar='Input Path',help='Path to planar density estimator files')
    parser.add_argument('-o', metavar='Output Path',help='Path to desired output directory')
    parser.add_argument('-r',metavar='Reduction Variable',help='Variable to reduce planar density files over')
    args = parser.parse_args()
        
    # Reduce all plane density estimators over specified reduction variable
    print 'Reducing planar density files' 
    reductionvar = '-r '+str(args.r)
    cmd = 'python $HOME/local/pimcscripts/reduce-one.py '+reductionvar+' -e planedensity '+args.p
    subprocess.call(cmd,shell=True)
    
    #Move to directory of reduced plane density file
    os.chdir(args.p) 
    
    # Grab the reduced planar density file
    for fn in os.listdir(args.p):
        if os.path.isfile(fn):
            if ('planedensity' and 'reduce') in fn:
                name = fn 
                reduceFile = open(fn,'r')
    
    # Get values of reduced variable 
    lines = reduceFile.readlines()
    reduceVals = lines[0].split() 
    reduceVals.pop(0)
    fixreduceVals = []
    for i in range(len(reduceVals)):
        if (args.r not in reduceVals[i]) and ('=' not in reduceVals[i]):
            fixreduceVals.append(reduceVals[i])  
    print reduceVals 
    print fixreduceVals 
    reduceFile.close()
    
    # Grab data for each reduced variable value and plot
    ncount = 0
    densitycount = 1
    figcount = 1
    for i in fixreduceVals:
        print "Making density plot %i of %i" % (figcount, len(fixreduceVals))
        titlevar = i
        n, density = np.loadtxt(name, unpack=True, usecols=(ncount, densitycount))
        N = int(np.sqrt(density.shape[0]))
        density = density.reshape([N,N])
        density = density.transpose()
        pl.figure(figcount)
        pl.pcolor(density,cmap='BuGn')
        pl.xlabel("X " + r"$(\AA)$",fontsize=18)
        pl.ylabel("Y " + r"$(\AA)$",fontsize=18)
        pl.xticks(np.arange(0,110,10),range(0,11))
        pl.yticks(np.arange(0,110,10),range(0,11))
        cbar = pl.colorbar()
        cbar.set_label("Particle Density"+ r"$\left(\frac{Particles}{\AA^3}\right)$",labelpad=40,fontsize=18,rotation=270)
        pl.title( r"$%s = %s K$"%(args.r, titlevar),fontsize=24)
        pl.savefig(args.o+'planardensityplot_'+args.r+'=%s.pdf'%i,format='pdf')
        print "Finished making density plot %i of %i" % (figcount, len(fixreduceVals))
        ncount += 3
        densitycount += 3
        figcount += 1
        
if __name__=='__main__':
    main() 
