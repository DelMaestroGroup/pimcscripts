# =============================================================================
# partitionPOVplot.py
#
# Author:           Max Graves
# Last Revision:    25-FEB-2013
#
# This calls on the povexport library to convert a configuration of visual
# Python spherical atoms into a POV-ray image.
# =============================================================================

import random,subprocess,sys
import visual as vis
import bipartPOVexport as pov
import numpy as np

def hex_to_rgb(hexi):
    hexi = hexi.lstrip('#')
    hlen = len(hexi)
    return tuple(int(hexi[i:i+hlen/3], 16) for i in range(0, hlen, hlen/3))


def main():

    # ----- User Defined Parameters -------------------------------------------
    # POV-ray image parameters
    HT = '1600'     # height
    WD = '2400'     # width

    # cell lengths
    Lx = 3.0    # angstroms
    Ly = 3.0
    Lz = 3.0

    # define atomic sphere radius
    sphRad = 0.20   # angstroms

    # --- define colors -------------------------------
    #partOneCol = (0.1,0.1,0.1)  # black
    #partTwoCol = (3.0,0.1,0.1)  # red
    #partOneCol = (0.1,0.1,1.0)  # navy
    
    # -------------------------------------------------
    # BIPARTITION FIGURE COLORS
    partTwoCol = (2.0,2.0,2.0)  # white
    partOneCol = () # teal-ish color
    for nv,val in enumerate(hex_to_rgb('#588492')): 
        partOneCol += ((1.5/150)*val,) 
    # -------------------------------------------------
    
    # --- define partition type -----------------------
    #partitionType = 'planar'
    #partitionType = 'spherical'
    partitionType = 'particle'

    # --- delete particles intersecting big sphere ---- 
    ridIntersect = False

    # ----- Generate atomic positions -----------------------------------------

    # create cubic lattice
    lattice = np.mgrid[-Lx:Lx+1, -Ly:Ly+1, -Lz:Lz+1]

    # randomly displace atoms and store in form for visual python 
    positions = []
    for ndim,latt in enumerate(lattice): # loop over x,y,z
        atomNum = 0
        for lat in latt:
            for lai,la in enumerate(lat):
                for li, l in enumerate(la): 
                    
                    # update by random displacement
                    la[li] += sphRad*(2.0*random.random()-1.0)
                    
                    # create list of lists of atomic positions
                    if ndim==0:
                        positions.append([la[li]])
                    else:
                        positions[atomNum].append(la[li])
                    
                    atomNum += 1
    
    # ----- visual python scene -----------------------------------------------
    scene = vis.display(title='World Lines!!',x=0, y=0, 
            width=800, height=844,
            center=(0,0,0), background=vis.color.white)
    scene.autoscale = 0

    # ----- partition cell by colors ------------------------------------------
    colors = []

    # --- planar ---
    if partitionType == 'planar':
        for atomPos in positions:
            if atomPos[0] < 0:
                colors.append(partOneCol)
            else:
                colors.append(partTwoCol)
    
    # --- central sphere ---
    if partitionType == 'spherical':    
        
        # compute radius of sphere that contains half the atoms in the cell
        per=0.7
        Rcut = ((1.0/per)*4.0/(3.0*np.pi)*Lx*Ly*Lz)**(1.0/3)

        centralRad = Rcut

        if ridIntersect:
            # get rid of atoms that intersect 'bubble'
            newpos = []
            for i, atomPos in enumerate(positions):
                atomRad = np.sqrt(atomPos[0]**2 + atomPos[1]**2 + atomPos[2]**2)
                if not (atomRad > Rcut and atomRad < Rcut+2.0*sphRad):
                    newpos.append(atomPos)

            print 'Cropped out',str(len(positions)-len(newpos)),'pesky buggers.'
            positions = newpos

            centralRad = sphRad+Rcut
        
        for atomPos in positions:
            atomRad = np.sqrt(atomPos[0]**2 + atomPos[1]**2 + atomPos[2]**2)
            if atomRad < Rcut:
                colors.append(partOneCol)
            else:
                colors.append(partTwoCol)

        # add sphere around central atoms
        vis.sphere(pos=(0,0,0), 
                radius=centralRad, color=partOneCol, opacity=0.3)

    # --- random ---
    if partitionType == 'particle':
        for atomPos in positions:
            rr = random.random()
            if rr < 0.5:
                colors.append(partOneCol)
            else:
                colors.append(partTwoCol)
    

    for numAtom,atomPos in enumerate(positions):
        vis.sphere(pos=atomPos, radius=sphRad, color=colors[numAtom])

    
    # ----- povray stuff ------------------------------------------------------

    # define povray file names
    povFileName = 'bipartView.pov'
    iniFileName =  povFileName[:-4]+'.ini'
    
    # create .pov file
    pov.export(scene, povFileName)
    
    # generate single image with povray
    command = ('povray', povFileName, 'Height='+HT,'Width='+WD)
    subprocess.check_call(command)

    # make video of rotation about the sample with POV-ray 
    '''
    # generate multiple images, same image but rotating about origin!
    import MTG_visTools as vt
    vt.writeINIfile(povFileName,HT,WD,finFrame)

    command = ('povray', iniFileName)
    subprocess.check_call(command)
    
    # check for Mencoder, the video processing utility.
    vt.findMencoder()
    
    # set up bash commands to run mencoder
    command = ('mencoder', 'mf://*.png', '-mf', 'type=png:w=800:h=600:fps=1',
           '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
           '-o', 'animatedStuff.avi')

    subprocess.check_call(command)
    '''

    # try to close visual python window
    try:
        import wx
        wx.Exit()
    except:
        pass
 
# =============================================================================
if __name__ == "__main__":
    main()
