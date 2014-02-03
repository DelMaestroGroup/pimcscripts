# =============================================================================
# visPath.py
#
# Author:           Adrian Del Maestro
#                   Max Graves
# Last Revision:    14-JAN-2013
#
# Main driver for visualizing worldlines from (g)ce-wl- files.
# =============================================================================

import os,sys
from optparse import OptionParser
import visual as vis
import MTG_visTools as vt
import povexport as pov

def main():

    # parse the command line options
    parser = OptionParser()
    (options, args) = parser.parse_args()
    fileName = args[0]

    # get cell lengths
    cellDims, excDims = vt.getCellDimensions(fileName)
    L = float(cellDims[0])
    Ly = float(cellDims[1])
    Lz = float(cellDims[2])
    if len(excDims) > 0:
        ay = float(excDims[0])
        az = float(excDims[1])
    # CHECK HOW THE OUTPUT LOOKS FROM THE CODE FOR D<3 !!!!!!

    # Load the configurations from disk	
    numFrames,wl = vt.loadPIMCPaths(fileName)
    print 'Number of frames: ',numFrames
    paths = []
    print numFrames
    for t in range(numFrames):
        paths.append(vt.Path(wl[t]))

    # choose frame number
    numFrame = 0

    # time slice data
    M = paths[numFrame].numTimeSlices
    dM = 1.0*L/(1.0*(M-1))

    # set up background
    scene = vis.display(title='World Lines!!',x=0, y=0, width=800, height=844,\
            center=(0,0,0), background=(0.0,0.0,0.0))
    scene.autoscale = 0
    scene.userzoom = False
    scene.userspin = False
    #scene.zoom = 2.0
    #scene.spin = 60

    # Set up excluded volume
    if len(excDims)>0:
        excVol = vis.box(pos=(0,0,0), length=L, height=ay, width=az, opacity=0.5)

    # Set up cell walls
    excVol = vis.box(pos=(0,0,0), length=L, height=Ly, width=Lz, opacity=0.1)

    # The boundary of the simulation box in space-time -- 1D
    #ymax = -0.5*L + (M-1)*dM
    #line = [(-0.5*L,-0.5*L,0),(0.5*L,-0.5*L,0),(0.5*L,ymax,0),(-0.5*L,ymax,0),(-0.5*L,-0.5*L,0)]
    #vis.curve(pos=line,radius=0.20,color=(0.5,0.5,0.5))
    
    print "displaying scene"
    wl = vt.WLFrame(paths[numFrame], L, Ly, Lz)
    vt.getScreenShot(0)
    
    print "Creating POVray file..."
    pov.export(scene, 'kittisCattis.pov')

    print "...finished exporting."

    '''
    n = 1
    for p in paths[1:]:
        #		visual.rate(1)
        wl.update(p)
        vt.getScreenShot(n)
        n += 1
    '''


# ----------------------------------------------------------------------
if __name__ == "__main__":
	    main()
