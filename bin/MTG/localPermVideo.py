# =============================================================================
# Script for creating video of the superfluid transition in Helium-4.
# Technically, what is seen is an increase in the local permutation
# number as seen in histogram form.
#
# Author:           Max Graves
# Last Revision:    4-FEB-2013
# =============================================================================

import matplotlib.gridspec as gridspec
import os,sys,glob,subprocess
import pylab as pl

# =============================================================================
def findMencoder():
    """ checks for Mencoder """
    try:
        subprocess.check_call(['mencoder'])
    except subprocess.CalledProcessError:
        print "mencoder command was found"
        pass
    except OSError:
        print 'could not find mencoder.'
        sys.exit("quitting\n")

# =============================================================================
def main():

    # holds all files, sorted by temp
    files = sorted(glob.glob("*locperm*") and glob.glob("*dat*"))[::-1]
    
    XYfileNames, YZfileNames, ZXfileNames, allfileNames = [] , [] , [] , []

    Lx = 8.0
    Ly = 8.0
    Lz = 24.0

    # create png corresponding to each file
    for fileName in files:
        # create list so mencoder makes movie in correct order
        print fileName

        # Open tensor file, determine number of grid boxes in each dimension
        inFile = open(fileName,'r')
        N = int(inFile.readline().split()[1])
        inFile.close()

        # Assuming a 3d data set, load the data
        data = pl.loadtxt(fileName).reshape([N,N,N])
            
        # plot histograms in all three projections
        pl.figure(1)    # xy plot
        pl.imshow(pl.sum(data,axis=2),extent=[-(Ly/2.0),(Ly/2.0),-(Lx/2.0),(Lx/2.0)])
        pl.xlabel(r'$y\  [\AA]$')
        pl.ylabel(r'$x\  [\AA]$')
        pl.title('Local Permutation Number (T = '+str(fileName[12:18])+')', size=20)
        pl.clim(0,7.5)
        pl.colorbar(shrink=0.4)
        XYfileNames.append(str(fileName).rstrip('.dat')+'-xy.png')
        pl.savefig(str(fileName).rstrip('.dat')+'-xy.png')
        pl.close()
        
        pl.figure(2)    # zx plot
        pl.imshow(pl.sum(data,axis=1), extent=[-(Lz/2.0),(Lz/2.0),-(Lx/2.0),(Lx/2.0)])
        pl.xlabel(r'$z\  [\AA]$')
        pl.ylabel(r'$x\  [\AA]$')
        pl.title('Local Permutation Number (T = '+str(fileName[12:18])+')', size=20)
        pl.clim(0,7.5)
        pl.colorbar(shrink=0.4)
        ZXfileNames.append(str(fileName).rstrip('.dat')+'-zx.png')
        pl.savefig(str(fileName).rstrip('.dat')+'-zx.png')
        pl.close()
       
        pl.figure(3)    # yz plot
        pl.imshow(pl.sum(data,axis=0), extent=[-(Lz/2.0),(Lz/2.0),-(Ly/2.0),(Ly/2.0)])
        pl.xlabel(r'$z\  [\AA]$')
        pl.ylabel(r'$y\  [\AA]$')
        pl.title('Local Permutation Number (T = '+str(fileName[12:18])+')', size=18)
        pl.clim(0,11.0)
        pl.colorbar(shrink=0.4)
        YZfileNames.append(str(fileName).rstrip('.dat')+'-yz.png')
        pl.savefig(str(fileName).rstrip('.dat')+'-yz.png')
        pl.close()

        pl.figure(4)    # all plots together
        pl.suptitle('Particle Density Projection (X-Y) -- T = '+str(fileName[12:18]), size=20)
        gs = gridspec.GridSpec(8,8)
        xyplot = pl.subplot(gs[2:4,6:8])
        pl.imshow(pl.sum(data,axis=2), extent=[-(Lx/2.0),(Lx/2.0),-(Ly/2.0),(Ly/2.0)] )
        pl.xlabel(r'$y\  [\AA]$')
        pl.ylabel(r'$x\  [\AA]$')
        labels = xyplot.get_xticklabels()
        for label in labels:
            label.set_rotation(90) 
        xyplot.yaxis.set_label_position('right')
        xyplot.xaxis.set_label_position('top')
        pl.clim(0,7.5)
        
        zxplot = pl.subplot(gs[1:4,:6])
        pl.imshow(pl.sum(data,axis=1), extent=[-(Lz/2.0),(Lz/2.0),-(Lx/2.0),(Lx/2.0)])
        pl.xlabel(r'$z\  [\AA]$')
        pl.ylabel(r'$x\  [\AA]$')
        pl.clim(0,7.5)
        pl.colorbar(shrink=0.7)
       
        zyplot = pl.subplot(gs[4:7,1:7])
        pl.imshow(pl.sum(data,axis=0), extent=[-(Lz/2.0),(Lz/2.0),-(Ly/2.0),(Ly/2.0)])
        pl.xlabel(r'$z\  [\AA]$')
        pl.ylabel(r'$y\  [\AA]$')
        pl.clim(0,11.0)
        pl.colorbar(shrink=0.7)

        #fig.tight_layout()
        allfileNames.append(str(fileName).rstrip('.dat')+'-all.png')
        pl.savefig(str(fileName).rstrip('.dat')+'-all.png')
        pl.close()
    
    # write fileNames in order to disk so mencoder knows the order
    fxy = open('XYfileNames.txt', 'w')
    fzx = open('ZXfileNames.txt', 'w')
    fyz = open('YZfileNames.txt', 'w')
    All = open('allfileNames.txt', 'w')
    fxy.write('\n'.join(XYfileNames))
    fzx.write('\n'.join(ZXfileNames))
    fyz.write('\n'.join(YZfileNames))
    All.write('\n'.join(allfileNames))
    fxy.close() , fyz.close() , fzx.close() , All.close()

    # create movie file from the png images just created
    findMencoder()

    # set up bash commands to run mencoder
    commandxy = ('mencoder', 'mf://@XYfileNames.txt', '-mf', 'type=png:w=800:h=600:fps=2',
           '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
           '-o', 'superFluidTransXY.avi')
    commandyz = ('mencoder', 'mf://@YZfileNames.txt', '-mf', 'type=png:w=800:h=600:fps=2',
           '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
           '-o', 'superFluidTransYZ.avi')
    commandzx = ('mencoder', 'mf://@ZXfileNames.txt', '-mf', 'type=png:w=800:h=600:fps=2',
           '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
           '-o', 'superFluidTransZX.avi')
    commandall = ('mencoder', 'mf://@allfileNames.txt', '-mf', 'type=png:w=800:h=600:fps=2',
           '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
           '-o', 'superFluidTransALL.mpg')


    print "\n\nabout to execute:\n%s\n\n" % ' '.join(commandxy)
    subprocess.check_call(commandxy)
    print "\n\nabout to execute:\n%s\n\n" % ' '.join(commandyz)
    subprocess.check_call(commandyz)
    print "\n\nabout to execute:\n%s\n\n" % ' '.join(commandzx)
    subprocess.check_call(commandzx)
    print "\n\nabout to execute:\n%s\n\n" % ' '.join(commandall)
    subprocess.check_call(commandall)



    # Clean up working directory -- delete png images and text files
    for f in glob.glob('*txt*'): os.remove(f)
    for f in glob.glob('*png*'): os.remove(f)

    print "\n Videos should now be in current working directory."
    
# =============================================================================
if __name__=='__main__':
    main()
