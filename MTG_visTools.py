# =============================================================================
#   visTools.py
#
# Author:           Adrian Del Maestro
#                   Max Graves
# Date Created:     12.15.2009  (pathconfig.py)
# Last Revised:     15-JAN-2014
#
# Includes classes that store:
#       - All information about world line trajectories of particles computed 
#           from path integral monte carlo calculations.
#       - All visual data for use with vPython.
#       - Some functions used by the main driver 
#
# This was adapted from pathconfig.py and vis1D.py
# =============================================================================

import os,sys,subprocess,argparse,commands
from numpy import *
#import numpy as np
import pimchelp
from optparse import OptionParser
from visual import *
#import visual as vis


# some global variables that we need to get rid of
colors1 = [(85,98,112),(78,205,196),(199,244,100),(255,107,107)]
colors = []
for color in colors1:
    c = (1.0*color[0]/255.0,1.0*color[1]/255.0,1.0*color[2]/255.0)
    colors.append(c)


XXX = -999
scale = 1.0

# =============================================================================
# Some useful functions
# =============================================================================

def findMencoder():
    """ checks for Mencoder, returns error and exits if it doesn't find it """
    try:
        subprocess.check_call(['mencoder'])
    except subprocess.CalledProcessError:
        print "mencoder command was found"
        pass
    except OSError:
        print 'could not find mencoder.'
        sys.exit("quitting\n")


def getCellDimensions(fileName):
    ''' Get the cell dimensions. '''
    cellDims = []
    excDims = []
    with open(fileName) as inFile:
        lineCount = 0
        for line in inFile:
            lineCount += 1
            if lineCount == 3:
                cellDims = line.split()
            if (lineCount == 4 and line[0] == '#'):
                excDims = line.split()
            if (lineCount > 4):
                break
    
    cellDims.pop(0)
    if len(excDims)>0:
        excDims.pop(0)

    return cellDims, excDims


def getScreenShot(frameNum):
    ''' 
    Take a screenshot and save it to disk.
    THIS CAN PROBABLY BE DISCARDED.
    '''
    fname = 'OUTPUT/POREMOVIE/He1d-%04d.tiff' % frameNum
    commands.getoutput('/usr/sbin/screencapture -x -t tiff %s' % fname)
    crop = '/opt/local/bin/convert ' + fname + ' -crop 800x800+0+44 ' + fname
    commands.getoutput(crop)
    print "Snap %d" % frameNum


def loadPIMCPaths(fileName):
    '''
    Given a pimc worldline file, parse and extract a list, where each element 
    contains lines corresponding to a given MC step worldline configuration.
    '''

    # Open the file
    wlFile = open(fileName,'r')
    lines = wlFile.readlines()

    n = 0
    wl = []
    for line in lines:

        # Every time we hit a new configuration, reset the data array and increment
        # a counter
        if line.find('START_CONFIG') != -1:
            n += 1
            data = []
        # If we encounter an end config, append to the wl list
        elif line.find('END_CONFIG') != -1: 
            wl.append(data)
        # Otherwise, append to the data array
        elif line[0] == '#':
            a = 0
        else:
            data.append(line.split())

    return n,wl


def parseCMD():
    ''' parse the command line. '''
    parser = argparse.ArgumentParser(description='multi-purpose Python juju.')
    parser.add_argument('fileNames', help='(g)ce-wl-.. data files', nargs='+')
    parser.add_argument('-o', '--output', type=str,
            choices=['single','rotate','bins'], default='single',
            help='Enter type of output you want')

    return parser.parse_args()


def sep(p1,p2):
    ''' 
    Get the scalar separation between two points.
    THIS IS NOT BEING USED ANYMORE.  
    '''

    d = 0.0
    for i in range(3):
        d += (p1[i]-p2[i])**2

    return sqrt(d)


def writeINIfile(povFileName,HT,WD,finFrame):
    ''' writes .ini file for POV-ray animation. '''
    
    # this is what goes into the .ini file
    iniTemplate = """; POV-Ray animation ini file
Antialias=Off
Antialias_Threshold=0.1
Antialias_Depth=2

Height=%(ht)s
Width=%(wd)s

Input_File_Name=%(iFName)s

Initial_Frame=1
Final_Frame=%(finalFrame)s
Initial_Clock=0
Final_Clock=1

Cyclic_Animation=on
Pause_when_Done=off
"""
    # put quotes around pofFileName for .ini file
    povfName = '"'+povFileName+'"'
    iniText = iniTemplate % {'iFName':povfName,
            'ht':HT,
            'wd':WD,
            'finalFrame':finFrame}

    # name .ini file according to .pov file name
    iniFileName =  povFileName[:-4]+'.ini'

    # write iniText to the iniFileName
    with open(iniFileName, 'w') as outFile:
        outFile.write(iniText)


# =============================================================================
# WLFrame class - holds all display information.
# =============================================================================
class WLFrame:
    ''' 
    Holds all the display information including beads and links.
    '''
    def __init__(self, path, L, Ly=None, Lz=None):
        '''
        Parent class contains arrays for beads and links as well as
        some of the aesthetics for visual python.
        This can be passed either one, two, or three cell lengths, and
        the number passed is assumed to be the spatial dimensionality 
        of the system.
        '''
        self.L = L
        self.Ly = Ly
        self.Lz = Lz

        self.dim = 3
        # determine spatial dimension from number of arguments passed
        if (Lz == None and Ly == None):
            self.dim = 1
        elif (Lz == None and Ly != None):
            self.dim = 2

        self.dM = 1.0*self.L/(1.0*(path.numTimeSlices-1))
        self.maxNumParticles = path.numParticles + 20

        self.disBead = []
        self.disPrevBead = []
        self.disNextBead = []
        self.disPrevLink = []
        self.disNextLink = []
        beadRad = 0.08
        linkRad = 0.03

        # We first do the beads
        for m in range(path.numTimeSlices):
            bList = []
            bnList = []
            bpList = []
            nList = []
            pList = []

            # Initialize arrays for beads and links
            for n in range(self.maxNumParticles):

                # Initialize bead array
                b = sphere(pos=(0.0,0.0,0.0),radius=beadRad,color=(0.0,0.0,0.0))
                b.visible = False
                bList.append(b)

                # Initialize pbc bead arrays (next and prev beads)
                bp = sphere(pos=(0.0,0.0,0.0),radius=beadRad,color=(0.0,0.0,0.0))
                bp.visible = False
                bn = sphere(pos=(0.0,0.0,0.0),radius=beadRad,color=(0.0,0.0,0.0))
                bn.visible = False
                bpList.append(bp)
                bnList.append(bn)

                # Initialize links (prev and next curves)
                cp = curve(pos=(0.0,0.0,0.0),radius=linkRad,color=(0.0,0.0,0.0))
                cn = curve(pos=(0.0,0.0,0.0),radius=linkRad,color=(0.0,0.0,0.0))
                cp.visible = False
                cn.visible = False
                pList.append(cp)
                nList.append(cn)

            self.disBead.append(bList)
            self.disPrevBead.append(bpList)
            self.disNextBead.append(bnList)
            self.disPrevLink.append(pList)
            self.disNextLink.append(nList)

        self.update(path)


    def update(self,path):
        ''' Update all beads and links. '''

        # Update the colormap
        maxWLIndex = path.wlIndex.max()+1

        # Turn off all actors in the scene
        for m in range(path.numTimeSlices):
            for n in range(self.maxNumParticles):
                self.disBead[m][n].visible = False
                self.disPrevBead[m][n].visible = False
                self.disNextBead[m][n].visible = False
                self.disPrevLink[m][n].visible = False
                self.disNextLink[m][n].visible = False

        # update all the bead and link positions
        self.updateBeads(path)
        self.updateDirectionLinks(path,path.prev,self.disPrevLink,self.disPrevBead,0)
        self.updateDirectionLinks(path,path.next,self.disNextLink,self.disNextBead,path.numTimeSlices-1)


    def updateBeads(self,path):
        ''' Update the bead positions. '''

        for m in range(path.numTimeSlices):

            for n in range(path.numParticles):
                if path.active[m,n]:

                    if path.wlIndex[m,n] == XXX:
                        col = colors[3]
                    else:
                        col = colors[1]

                    # update position based on dimension
                    if self.dim == 1:
                        x = path.bead[m,n,0]
                        y = -0.5*self.L + m*self.dM
                        z = 0.0
                    elif self.dim == 2:
                        x = path.bead[m,n,0]
                        y = path.bead[m,n,1]
                        z = -0.5*self.L + m*self.dM
                    else:
                        x = path.bead[m,n,0]
                        y = path.bead[m,n,1]
                        z = path.bead[m,n,2]

                    p0 = (x,y,z)
                    self.disBead[m][n].pos = p0
                    self.disBead[m][n].color = col

                    if not self.disBead[m][n].visible:
                        self.disBead[m][n].visible = True;


    def updateDirectionLinks(self,path,shift,link,shiftBead,boundary):
        ''' Update the link positions. '''

        for m in range(path.numTimeSlices):
            for n in range(path.numParticles):

                if path.active[m,n]:
                    if path.wlIndex[m,n] == XXX:
                        lcol = colors[3]
                        bcol = colors[3]
                    else:
                        lcol = colors[2]
                        bcol = colors[1]


                    if m != boundary and shift[m,n,0] != XXX:
                        ms = shift[m,n,0]
                        ns = shift[m,n,1]
                        if path.active[ms][ns]:
                            p1 = copy(self.disBead[m][n].pos)
                            p2 = copy(self.disBead[ms][ns].pos)
                            
                            windBack = False
                            # enforce PBC in x-direction
                            if abs(p1[0]-p2[0]) > 0.5*self.L:
                                windBack = True
                                if p2[0] < 0:
                                    p2[0] += self.L
                                else:
                                    p2[0] -= self.L

                            # enforce PBC in y-direction
                            if (self.dim==2 or self.dim==3):
                                if abs(p1[1]-p2[1]) > 0.5*self.Ly:
                                    windBack = True
                                    if p2[1] < 0:
                                        p2[1] += self.Ly
                                    else:
                                        p2[1] -= self.Ly

                            # enforce PBC in z-direction
                            if (self.dim==3):
                                if abs(p1[2]-p2[2]) > 0.5*self.Lz:
                                    windBack = True
                                    if p2[2] < 0:
                                        p2[2] += self.Lz
                                    else:
                                        p2[2] -= self.Lz
                            
                            if windBack:
                                shiftBead[m][n].pos = p2
                                shiftBead[m][n].visible = True
                                shiftBead[m][n].color = bcol

                            link[m][n].pos = [p1,p2]
                            link[m][n].color = lcol
                            link[m][n].visible = True
                        else:
                            link[m][n].visible = False


# =============================================================================
# Path class - holds all information about system at given MC step.
# =============================================================================
class Path:
    ''' 
    Holds all the information about the worldline configuration at a given
    MC step.
    '''

    def __init__(self,wlData,num_particles=None):
        ''' 
        Given a single worldine configuration, we populate the data arrays.
        We expect data as a list containing lines of a worldine file that
        have already been split.
        '''

        # First we determine how many particles there are in this config
        self.numParticles = XXX
        self.numTimeSlices = XXX
        for line in wlData:
            if int(line[0]) > self.numTimeSlices:
                self.numTimeSlices= int(line[0])
            if int(line[1]) > self.numParticles:
                self.numParticles = int(line[1])

        self.numTimeSlices += 1
        self.numParticles  += 1

        if not num_particles or (num_particles == self.numParticles):

            # Now we initialze all data members
            self.bead    = zeros([self.numTimeSlices,self.numParticles,3],float)
            self.active  = zeros([self.numTimeSlices,self.numParticles],int)
            self.next    = zeros([self.numTimeSlices,self.numParticles,2],int)
            self.prev    = zeros([self.numTimeSlices,self.numParticles,2],int)
            self.wlIndex = zeros([self.numTimeSlices,self.numParticles],int)

            for line in wlData:
                m = int(line[0])
                n = int(line[1])
                self.wlIndex[m,n] = int(line[2])
                self.bead[m,n,0] = float(line[3])   # current x position
                self.bead[m,n,1] = float(line[4])   # current y position
                self.bead[m,n,2] = float(line[5])   # current z position
                self.prev[m,n,0] = int(line[6])
                self.prev[m,n,1] = int(line[7])
                self.next[m,n,0] = int(line[8])
                self.next[m,n,1] = int(line[9])
                self.active[m,n] = 1

            # Finally store the wlData itself
            # if any(self.wlIndex==-999):
            #     self.wlData = None
            # else:
            self.wlData = wlData

        else:
            self.wlData = None

