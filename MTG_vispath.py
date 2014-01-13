'''
visPath.py

Visualize a 3D configuration from PIMC data with excluded volume
'''

import os,sys
#import loadgmt,kevent
import pimchelp
from optparse import OptionParser
from visual import *
#from pylab import *
import commands
import MTG_visTools as vt

colors1 = [(85,98,112),(78,205,196),(199,244,100),(255,107,107)]
colors = []
for color in colors1:
	c = (1.0*color[0]/255.0,1.0*color[1]/255.0,1.0*color[2]/255.0)
	colors.append(c)

XXX = -999
scale = 1.0
def main():
    parser = OptionParser()

    # parse the command line options
    (options, args) = parser.parse_args()
    fileName = args[0]

    # Get the parameter map
    L = 40

    # Load the configurations from disk	
    numFrames,wl = vt.loadPIMCPaths(fileName)
    print 'Number of frames: ',numFrames
    paths = []
    for t in range(numFrames):
        paths.append(vt.Path(wl[t]))

    # The system size and number of imaginary time steps
    #L = parMap['V']

    M = paths[0].numTimeSlices
    dM = 1.0*L/(1.0*(M-1))

    #wl = WLFrame(paths[0],L)


    scene = display(title='1D Bose Gas',x=0, y=0, width=800, height=844,\
            center=(0,0,0), background=(0.0,0.0,0.0),range=0.6*L)
    scene.autoscale = 0


    # The boundary of the simulation box in space-time
    ymax = -0.5*L + (M-1)*dM
    line = [(-0.5*L,-0.5*L,0),(0.5*L,-0.5*L,0),(0.5*L,ymax,0),(-0.5*L,ymax,0),(-0.5*L,-0.5*L,0)]
    #curve(pos=line,radius=0.20,color=(0.5,0.5,0.5))
    print 'got'

    wl = vt.WLFrame(paths[0],L)
    getScreenShot(0)

    n = 1
    for p in paths[1:]:
        #		visual.rate(1)
        wl.update(p)
        getScreenShot(n)
        n += 1


# ----------------------------------------------------------------------
if __name__ == "__main__":
	    main()
