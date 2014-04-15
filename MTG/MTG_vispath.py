# =============================================================================
# visPath.py
#
# Author:           Adrian Del Maestro
#                   Max Graves
# Last Revision:    14-JAN-2013
#
# Main driver for visualizing worldlines from (g)ce-wl- files.
# =============================================================================

import os,sys,subprocess,glob
import visual as vis
import MTG_visTools as vt
import povexport as pov

def main():

    # TEST
    #HT = '300'
    #WD = '450'
    #finFrame = '30'

    # PRODUCTION Rotate Video Parameters
    HT = '1200'
    WD = '1800'
    finFrame = '120'

    # parse the command line options
    args = vt.parseCMD() 
    fileName = args.fileNames[0]
    outType = args.output
 
    # check if images already exist in current directory.  If not, create them.
    if glob.glob('*png*') == []:

        # get cell lengths
        cellDims, excDims = vt.getCellDimensions(fileName)
        L = float(cellDims[0])
        Ly = float(cellDims[1])
        Lz = float(cellDims[2])
        if len(excDims) > 0:
            ay = float(excDims[0])
            az = float(excDims[1])
        else:
            ay, az = None, None

        # Load the configurations from disk	
        numFrames,wl = vt.loadPIMCPaths(fileName)
        print 'Number of frames: ',numFrames
        
        paths = []
        for t in range(numFrames):
            paths.append(vt.Path(wl[t],Ly,Lz))

        # choose frame number from (g)ce-wl file.
        numFrame = 0

        # time slice info.
        M = paths[numFrame].numTimeSlices
        dM = 1.0*L/(1.0*(M-1))
        
        # ----- This is where you define cell characteristics -----------------
        scene = vis.display(title='World Lines!!',x=0, y=0, 
                width=800, height=844,
                center=(0,0,0), background=(1.0,1.0,1.0))
        scene.autoscale = 0
        
        # Set up excluded volume
        if len(excDims)>0:
            excVol = vis.box(pos=(0,0,0), length=L, 
                    height=ay, width=az, opacity=0.7)

        # Set up cell walls
        cellWalls = vis.box(pos=(0,0,0), length=L, 
                height=Ly, width=Lz, opacity=0.05)

        # The boundary of the simulation box in space-time -- 1D
        #ymax = -0.5*L + (M-1)*dM
        #line = [(-0.5*L,-0.5*L,0),(0.5*L,-0.5*L,0),(0.5*L,ymax,0),
        #        (-0.5*L,ymax,0),(-0.5*L,-0.5*L,0)]
        #vis.curve(pos=line,radius=0.20,color=(0.5,0.5,0.5))

        # ---------------------------------------------------------------------

        wl = vt.WLFrame(paths[numFrame], L, Ly, Lz)
        
        # define povray file names
        povFileName = 'wlview.pov'
        iniFileName =  povFileName[:-4]+'.ini'
        
        pov.export(scene, povFileName)

        # -- generate single image with povray --------------------------------
        if outType == 'single':
            command = ('povray', povFileName, 'Height='+HT,'Width='+WD)
            subprocess.check_call(command)

        # -- generate multiple images, same image but rotating about origin ---
        if outType == 'rotate':
            vt.writeINIfile(povFileName,HT,WD,finFrame)

            command = ('povray', iniFileName)
            subprocess.check_call(command)

        # -- generate images showing evolving worldlines ----------------------
        if outType == 'bins':
            # UPDATE THIS!
            '''n = 1
            for p in paths[1:]:
                #		visual.rate(1)
                wl.update(p)
                vt.getScreenShot(n)
                n += 1
            '''
  
    
    # ---- generate video from series of images -------------------------------
    if outType != 'single':
        # check for Mencoder, the video processing utility.
        vt.findMencoder()
        
        # set up bash commands to run mencoder
        command = ('mencoder', 'mf://*.png', '-mf', 'type=png:w=800:h=600:fps=1',
               '-ovc', 'lavc', '-lavcopts', 'vcodec=mpeg4', '-oac', 'copy',
               '-o', 'animatedStuff.avi')

        print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
        subprocess.check_call(command)

        print "\n The movie was written to 'animatedStuff.avi'"

 
    # try to close visual python window
    try:
        import wx
        wx.Exit()
    except:
        pass
 
# =============================================================================
if __name__ == "__main__":
    main()
