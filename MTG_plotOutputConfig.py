
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import argparse


# =============================================================================
def printHeader():
    '''
    print snazzy header for fun.
    '''
    print r" ____    ____ __________ ________    "
    print r"|    \  /    |          |        |   "
    print r"|     \/     |          |        |   "
    print r"|            |__      __|    ____|   "
    print r"|            |  |    |  |   |  _____ "
    print r"|   |\  /|   |  |    |  |   | |__   |"
    print r"|   | \/ |   |  |    |  |   |    |  |"
    print r"|   |    |   |  |    |  |   |____|  |"
    print r"|   |    |   |  |    |  |           |"
    print r"|___|    |___|  |____|  |___________|"


def getDims(fileName):
    """
    Reads in cell dimensions from data file.
    Taken from pimchelp.py
    """
    inFile = open(fileName,'r');
    inLines = inFile.readlines();
    cellDims = inLines[1].split()
    cellDims.pop(0)
    excDims = inLines[2].split()
    excDims.pop(0)
    inFile.close()
    return cellDims,excDims


def buildExcVol(cutX,cutY,cutZ,ax):
    """
    Builds rectangular prism (excluded volume).
    """
    r1 = Rectangle((-cutX,-cutZ),2.0*cutX,2.0*cutZ)
    r4 = Rectangle((-cutX,-cutZ),2.0*cutX,2.0*cutZ)
    r2 = Rectangle((-cutY,-cutZ),2.0*cutY,2.0*cutZ)
    r5 = Rectangle((-cutY,-cutZ),2.0*cutY,2.0*cutZ)
    r3 = Rectangle((-cutX,-cutY),2.0*cutX,2.0*cutY)
    r6 = Rectangle((-cutX,-cutY),2.0*cutX,2.0*cutY)
    r1.set_alpha(0.1)
    r2.set_alpha(0.1)
    r3.set_alpha(0.1)
    r4.set_alpha(0.1)
    r5.set_alpha(0.1)
    r6.set_alpha(0.1)
    ax.add_patch(r1), ax.add_patch(r2)
    ax.add_patch(r3), ax.add_patch(r4)
    ax.add_patch(r5), ax.add_patch(r6)

    art3d.pathpatch_2d_to_3d(r1, z=-cutY, zdir="y")
    art3d.pathpatch_2d_to_3d(r4, z=cutY, zdir="y")
    art3d.pathpatch_2d_to_3d(r2, z=cutX, zdir="x")
    art3d.pathpatch_2d_to_3d(r5, z=-cutX, zdir="x")
    art3d.pathpatch_2d_to_3d(r3, z=cutZ, zdir="z")
    art3d.pathpatch_2d_to_3d(r6, z=-cutZ, zdir="z")


def buildOuter(Lx,Ly,Lz,ax,col='Purple'):
    """
    Builds rectangular prism (outer cell walls).
    """
    r7 = Rectangle((-Lx/2,-Lz/2),Lx,Lz, color=col)
    r8 = Rectangle((-Lx/2,-Lz/2),Lx,Lz, color=col)
    r9 = Rectangle((-Ly/2,-Lz/2),Ly,Lz, color=col)
    r10 = Rectangle((-Ly/2,-Lz/2),Ly,Lz, color=col)
    r11 = Rectangle((-Lx/2,-Ly/2),Lx,Ly, color=col)
    r12 = Rectangle((-Lx/2,-Ly/2),Lx,Ly, color=col)
    r7.set_alpha(0.1)
    r8.set_alpha(0.1)
    r9.set_alpha(0.1)
    r10.set_alpha(0.1)
    r11.set_alpha(0.1)
    r12.set_alpha(0.1)
    ax.add_patch(r7), ax.add_patch(r8)
    ax.add_patch(r9), ax.add_patch(r10)
    ax.add_patch(r11), ax.add_patch(r12)

    art3d.pathpatch_2d_to_3d(r7, z=-Ly/2, zdir="y")
    art3d.pathpatch_2d_to_3d(r8, z=Ly/2, zdir="y")
    art3d.pathpatch_2d_to_3d(r9, z=Lx/2, zdir="x")
    art3d.pathpatch_2d_to_3d(r10, z=-Lx/2, zdir="x")
    art3d.pathpatch_2d_to_3d(r11, z=Lz/2, zdir="z")
    art3d.pathpatch_2d_to_3d(r12, z=-Lz/2, zdir="z")


# =============================================================================
def main():

    printHeader()

    # read in data
    fName = "gce-wl-05.000-040.441--003.000-0.00400-159288391.dat"
    data = np.loadtxt(fName)
    data = data.transpose()
    
    # determine how many beads per WL
    numBeads = 0
    for entry in data[0]:
        if entry >= numBeads:
            numBeads = entry
        else:
            break
    numBeads += 1
    numBeads = int(numBeads)
    print 'Number of beads per WL: ',numBeads

    data = data.transpose()
    WLDATA = {}

    for n, d in enumerate(data):
        num = str(int(d[1]))
        if n%50 == 0:
            WLDATA['N'+num+'_x'] = np.array([d[3]])
            WLDATA['N'+num+'_y'] = np.array([d[4]])
            WLDATA['N'+num+'_z'] = np.array([d[5]])
        WLDATA['N'+num+'_x'] = np.append(WLDATA['N'+num+'_x'], d[3])
        WLDATA['N'+num+'_y'] = np.append(WLDATA['N'+num+'_y'], d[4])
        WLDATA['N'+num+'_z'] = np.append(WLDATA['N'+num+'_z'], d[5])

    
    print WLDATA['N'+str(0)+'_x']

    #for i in range(numBeads):
    #    print WLDATA['N'+str(i)+'_x']

    # load in initial configs
    X = WLDATA['N'+str(0)+'_x']
    Y = WLDATA['N'+str(0)+'_y']
    Z = WLDATA['N'+str(0)+'_z']

    # get cell dimensions from data file
    cellDims, excDims = getDims(fileName)
    
    # assign values to variables from data file
    Lx, Ly, Lz = float(cellDims[0]), float(cellDims[1]), float(cellDims[2])
    cutX, cutY, cutZ = 0.5*Lx, float(excDims[0]), float(excDims[1])
    
    # set up plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    # set up excluded volume walls
    buildExcVol(cutX,cutY,cutZ,ax)

    # set up outer cell wells
    buildOuter(Lx,Ly,Lz,ax)

    # plot the actual configuration
    ax.scatter(X, Y, Z, s=60,color='Navy')

    #ax.set_xlim3d(-Lx/2-1, Lx/2+1)
    #ax.set_ylim3d(-Ly/2-1, Ly/2+1)
    ax.set_xlim3d(-Lz/2-1, Lz/2+1)
    ax.set_ylim3d(-Lz/2-1, Lz/2+1)
    ax.set_zlim3d(-Lz/2-1, Lz/2+1)
    ax.set_xlabel('x'+r'$[\AA]$')
    ax.set_ylabel('y'+r'$[\AA]$')
    ax.set_zlabel('z'+r'$[\AA]$')
    plt.show()






if __name__=='__main__':
    main()
