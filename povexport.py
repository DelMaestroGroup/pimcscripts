
"""
This module exports a VPython scene as POV-Ray scene description code.
Lights and camera location from the current visual scene are included.
Optionally, you may specify a list of include files,
   and pov textures for objects.
For an example of its use, see 'povexample.py'.
Currently convex, faces, and points objects are not exported.
Further documentation is found at the start of the file.
"""

# This module exports a description of a VPython scene as a .pov file that can be
# read by the multiplatform free program POV-Ray, which produces a high-quality
# ray-traced image, with shadows (by default). The output is a targa file (.tga)
# which can be further processed by many programs, including Photoshop.

# To use, import the export routine from this file ("from povexport import export").
# When you have the VPython scene you want, execute "export()", then read the resulting
# .pov file into POV-ray.

# Here are the export routine's default options:
# export(display=None, filename=None, include_list=None, xy_ratio=4./3., custom_text='', shadowless=0)
# If no display is specified, the current display is used.
# If no filename is specified, the title of your display will be used for the .pov file name.
# include_list lets you add your own POV-ray include statements to the .pov file.
# xy_ratio=4./3. matches the typical 4/3 aspect ration of POV-Ray output targa files.
#   If you specify a different resolution in POV-Ray (800x800, say, rather than 800x600),
#   you need to specify a corresponding xy_ratio value (which would be 1. for 800x800).
# custom_text lets you add any kind of POV-Ray statements to the .pov file.
# shadowless=1 makes POV-Ray not produce any shadows.
# You can give an individual object the attribute no_shadow,
#    which if True means it casts no shadow.
# You can give an individual object the attribute no_reflection,
#    which if True means it is not reflected by anything.

#--------------------------------------------------------------------------

# ruth chabay, nc state university (rwchabay@ncsu.edu)
# v1.0 2000-12-17

# Markus Gritsch (gritsch@iue.tuwien.ac.at)
# v.1.1   2001-03-09
#   *) replaced 'scene' by 'display' everywhere
#   *) added spheres at the joints of a curve
#   *) consistent pov_texture usage in process_arrow() also for the shaft
#   *) ambient light, light sources, up, and fov are now handled correctly
#   *) some cosmetic changes to the code
# v.1.1.1 2001-03-22
#   *) added 'shadowless' keyword to export()

# Ruth Chabay
# 2001-06-23
# hack to fix error in export_curve introduced by Python 2.1
# now can't assign obj.color = array(...)

# Markus Gritsch (gritsch@iue.tuwien.ac.at)
# v.1.2   2001-11-23
#   *) added 'xy_ratio' and 'custom_text' keywords to export()
#   *) the pov-strings are now directly written to a file

# Bruce Sherwood
# 2004-07-18
# add dictionary ("legal") for identifying an object so that
# povexport will continue to work with the new Boost-based VPython
# which changes the details returned by object.__class__

# Bruce Sherwood
# 2005-06-27
# in export_curve, move del's into if, otherwise empty curve gives error
# use scene.range to scale values to about 10, because
#   POV-Ray dies with very large (or very small?) numbers (e.g. 1e10)

# Bruce Sherwood
# 2005-07-22
# corrected the new scaling (see 2005-06-27) for frame handling

# Bruce Sherwood
# 2005-07-26
# corrected the new scaling (see 2005-06-27) of shaft in arrow

# Bruce Sherwood
# 2005-12-06
# corrected error in reporting the name of an unsupported object

# Ruth Chabay
# 2007-11-13
# changes to handle zero-length arrows

# Ruth Chabay
# 2008-11-05
# changes for compatibility with visual 5 (not compatible with visual 3)

# Bruce Sherwood
# 2008-11-24
# handle opacity (works with either visual 3 or visual 5, even though visual 3 won't show it)
# make povexport work with either visual 3 or visual 5
# create ellipsoids

# Bruce Sherwood
# 2008-11-25
# correct axis error in ellipsoids
# create lights with color for visual 5; create pyramids
# add more documentation

# Scott David Daniels
# 2009-05-22
# correct opacity error in creating an arrow out of a box and a pyramid

# Bruce Sherwood
# 2009-06-13
# Visual object attribute no_shadow=True means that object doesn't cast a shadow

# Bruce Sherwood
# 2009-06-14
# Handle arrow.fixedwidth correctly (make headlength 0.5*length if necessary)

# Guy Kloss
# 2009-06-23
# Further corrections to arrow, first implementation for basic label (sans serif, no box, etc.)

# Bruce Sherwood
# 2009-07-06
# Yet more corrections to arrow, to do exactly the same sizing done in the Visual C++ code

# Bruce Sherwood
# 2009-08-26
# Visual object attribute no_reflection=True means that object doesn't cast a shadow

# Max Graves        mtg6193@gmail.com
# 2014-02-04
# All sorts of stuff.  See --MTG flag.
# Added comments to functions.

## NOTE: when changing this module please change the following string:
#POVEXPORT_VERSION = "povexport 2009-08-26 for Visual 3 and Visual 5"
POVEXPORT_VERSION = "povexport 2014-01-26 for Visual 3 and Visual 5"

from visual import *
import argparse

legal = {frame:'frame', sphere:'sphere', box:'box', cylinder:'cylinder',
                   curve:'curve', ring:'ring', arrow:'arrow', label: 'label',
                   cone:'cone', ellipsoid:'ellipsoid', pyramid:'pyramid'}
ihat=vector(1, 0, 0)
jhat=vector(0, 1, 0)
khat=vector(0, 0, 1)
displayscale = 1.0 # global scale factor to adjust scene.range to 100

# DEFINE CAMERA ROTATION VARIABLES --MTG
rotateCamera1 = 0
#rotateCamera2 = '270*clock'
rotateCamera2 = 85
rotateCamera3 = 0

cameraZoom = 1.5


# DEFINE FINISH CHARACTERISTICS --MTG
amb = '0.2'
dif = '1.0'
refl = '0.8'
boxrefl = '0'
boxdif = '1.0'
boxamb = '0.2'

cylinderTrans = 1.0


def parseCMD(): # --MTG
    ''' parse the command line. '''
    parseDesc = 'turn visual python scene into .pov file'
    parser = argparse.ArgumentParser(description=parseDesc)
    parser.add_argument('fileNames', help='(g)ce-wl-.. data files', nargs='+')
    parser.add_argument('-b','--backgroundColor',type=str,
            default = 'white',
            help='Enter Background color')

    return parser.parse_args()


def version():
    ''' print version of povexport. '''
    return POVEXPORT_VERSION


def getpolar(a):
    '''
    takes a vector and find rotation angles (standard polar coord).
    '''
    xy = sqrt(a.x**2 + a.y**2)
    theta = atan2(xy, a.z)
    phi = atan2(a.y, a.x)
    
    return [theta, phi]


def find_rotations(a):
    ''' find rotations '''
    theta, phi = getpolar(a.axis)
    # find rotation around x-axis (if a.up <> jhat)
    # "undo" theta & phi rotations so can find alpha
    aup = a.up*1.0
    aup = rotate(aup, axis=khat, angle=-phi)
    aup = rotate(aup, axis=jhat, angle=pi/2.0-theta)
    a_sin = dot(cross(jhat, norm(aup)), ihat)
    a_cos = dot(norm(aup), jhat)
    alpha = atan2(a_sin, a_cos)
    
    return (alpha, theta, phi)


def process_frame(a, code):
    ''' add in frame rotations & translations (may be nested). '''
    frame_code = ''
    fr = a.frame
    while fr:
        alpha, theta, phi = find_rotations(fr)
        frx=alpha*180./pi
        fry=-90+theta*180./pi
        frz=phi*180./pi
        rrot = '    rotate <%f, %f, %f>\n'
        ttrn = '    translate <%f, %f, %f>\n'
        frame_code += rrot % (frx, fry, frz)
        frame_code += ttrn % (displayscale*fr.x, displayscale*fr.y, displayscale*fr.z)
        fr = fr.frame

    # insert frame_code at end (these rot's must be done last)
    end = code.rfind('}')
    code = code[:end] + frame_code + code[end:]

    return code


def add_texture(a, code):
    ''' add in user-specified texture (will override color) '''
    mat = None
    if hasattr(a, 'pov_texture'):
        mat = a.pov_texture
    # Maybe could interpret materials like this (need to include POV-Ray definitions):
    # elif hasattr(a, 'material'):
    # if a.material == materials.wood:
    # mat = "T_Wood20"
    if mat:
        tstring = '    texture { '+ mat + ' }\n'
        end = code.rfind('}')
        code = code[:end] + tstring + code[end:] 
    
    return code


def no_shadow(a):
    if hasattr(a,"no_shadow") and a.no_shadow:
        return "no_shadow"
    else:
        return ""


def no_reflection(a):
    if hasattr(a,"no_reflection") and a.no_reflection:
        return "no_reflection"
    else:
        return ""


def transparency(a):
    if hasattr(a,"opacity"):
        return 1.0-a.opacity
    else:
        return 0.0


def export_sphere(a):
    sphere_template = """
sphere {
    <%(posx)f, %(posy)f, %(posz)f>, %(radius)f
    pigment {color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f>}
    %(no_shadow)s
    finish {
    reflection %(reflection)s
    diffuse %(diffuse)s
    ambient %(ambient)s }
}
"""
    object_code = sphere_template % { 'posx':displayscale*a.x, 'posy':displayscale*a.y, 'posz':displayscale*a.z,
            'radius':displayscale*a.radius,
            'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
            'no_shadow':no_shadow(a),
            'reflection':refl,
            'diffuse':dif,
            'ambient':amb}
    #'no_shadow':no_shadow(a), 'no_reflection':no_reflection(a)}
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_ellipsoid(a):
    ellipsoid_template = """
sphere {
    <0, 0, 0>, %(radius)f
    pigment {color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f>}
    scale <%(sizex)f, %(sizey)f, %(sizez)f>
    rotate <%(rotx)f, %(roty)f, %(rotz)f>
    translate <%(posx)f, %(posy)f, %(posz)f>
    %(no_shadow)s
    %(no_reflection)s
}
"""
    alpha, theta, phi = find_rotations(a)
    object_code = ellipsoid_template % { 'posx':displayscale*a.x, 'posy':displayscale*a.y, 'posz':displayscale*a.z,
                    'radius':displayscale*a.size[0]/2.,
                    'sizex':1.0, 'sizey':a.size[1]/a.size[0], 'sizez':a.size[2]/a.size[0],
                    'rotx':alpha*180./pi, 'roty':-90.+theta*180./pi, 'rotz':phi*180./pi,
                    'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
                    'no_shadow':no_shadow(a), 'no_reflection':no_reflection(a) }
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_box(a):
    # create box at origin along x-axis
    # then rotate around x,y,z axes
    # then translate to final location
    box_template = """
box {
    <%(posx)f, %(posy)f, %(posz)f>, <%(pos2x)f, %(pos2y)f, %(pos2z)f>
    pigment {color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f>}
    rotate <%(rotx)f, %(roty)f, %(rotz)f>
    translate <%(transx)f, %(transy)f, %(transz)f>
    %(no_shadow)s
    finish{
    reflection %(reflection)s
    diffuse %(diffuse)s
    ambient %(ambient)s }

}
"""
    alpha, theta, phi = find_rotations(a)
    # pos of visual box is at center
    # generate two opposite corners for povray
    pos1=-displayscale*a.size / 2.0
    pos2=+displayscale*a.size / 2.0

    object_code = box_template % { 'posx':pos1.x, 'posy':pos1.y, 'posz':pos1.z,
                                   'pos2x':pos2.x, 'pos2y':pos2.y, 'pos2z':pos2.z,
                                   'rotx':alpha*180./pi, 'roty':-90.+theta*180./pi, 'rotz':phi*180./pi,
                                   'transx':displayscale*a.x, 'transy':displayscale*a.y, 'transz':displayscale*a.z,
                                   'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
                                    'no_shadow':no_shadow(a), 
                                    'reflection':boxrefl,
                                    'diffuse':boxdif,
                                    'ambient':boxamb}
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_cylinder(a):
    cylinder_template = """
cylinder {
    <%(posx)f, %(posy)f, %(posz)f>,<%(pos2x)f, %(pos2y)f, %(pos2z)f>, %(radius)f
    pigment { color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f> }
    %(no_shadow)s
    finish { 
    reflection %(reflection)s
    diffuse %(diffuse)s
    ambient %(ambient)s }
}
"""
    endpt1=displayscale*a.pos
    endpt2=displayscale*(a.pos+a.axis)
    object_code = cylinder_template % { 'posx':endpt1.x, 'posy':endpt1.y, 'posz':endpt1.z,
                                        'pos2x':endpt2.x, 'pos2y':endpt2.y, 'pos2z':endpt2.z,
                                        'red':a.red, 'green':a.green, 'blue':a.blue, 
                                        #'transparency':transparency(a),
                                        'transparency':cylinderTrans,
                                        'radius':displayscale*a.radius,
                                        'no_shadow':no_shadow(a),
                                        'reflection':refl,
                                        'diffuse':dif,
                                        'ambient':amb}
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_curve(a):
    object_code = ''
    ccyl = cylinder(visible=False) # create cylinder and sphere that can be deleted
    csph = sphere(visible=False)
    curve_no_shadow = no_shadow(a)
    curve_no_reflection = no_reflection(a)
    if len(a.pos) > 1:
        ii = 0
        while ii < len(a.pos)-1:
            endpt1 = a.pos[ii]
            endpt2 = a.pos[ii+1]
            if a.radius > 0:
                rr = a.radius
            else:
                rr = mag(endpt1-endpt2)/200.
            # create a dummy cylinder to export
            curve_color = (float(a.color[ii][0]),
                                   float(a.color[ii][1]),
                                   float(a.color[ii][2]))
            if endpt2[0] == endpt1[0] and endpt2[1] == endpt1[1] and endpt2[2] == endpt1[2]:
                ii += 1
                continue
            ccyl = cylinder(pos=endpt1, axis=(endpt2-endpt1),
                            radius=rr, color=curve_color,
                            frame=a.frame, no_shadow=curve_no_shadow,
                            no_reflection=curve_no_reflection, visible=0)
            csph = sphere(pos=endpt1, radius=rr, color=curve_color,
                          frame=a.frame, no_shadow=curve_no_shadow,
                          no_reflection=curve_no_reflection, visible=0)
            if hasattr(a, 'pov_texture'):
                ccyl.pov_texture = a.pov_texture
                csph.pov_texture = a.pov_texture
            object_code += export_cylinder(ccyl) + export_sphere(csph)
            ii += 1
        endpt1 = a.pos[ii]
        csph = sphere(pos=endpt1, radius=rr, color=curve_color,
                      frame=a.frame, no_shadow=curve_no_shadow,
                      no_reflection=curve_no_reflection, visible=0)
        if hasattr(a, 'pov_texture'):
            csph.pov_texture = a.pov_texture
        object_code += export_sphere(csph)
        del(ccyl)
        del(csph)
    
    return object_code


def export_ring(a):
    torus_template = """
torus {
    %(radius)f, %(minorradius)f
    pigment { color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f> }
    rotate <0.0, 0.0, -90.0> // align with x-axis
    rotate <%(rotx)f, %(roty)f, %(rotz)f>
    translate <%(transx)f, %(transy)f, %(transz)f>
    %(no_shadow)s
    %(no_reflection)s
}
"""
    ang = getpolar(a.axis)
    theta = ang[0]
    phi = ang[1]
    object_code = torus_template % { 'radius':displayscale*a.radius, 'minorradius':displayscale*a.thickness,
                                     'transx':displayscale*a.x, 'transy':displayscale*a.y, 'transz':displayscale*a.z,
                                     'rotx':0.0, 'roty':-90.+theta*180./pi, 'rotz':phi*180./pi,
                                     'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
                                     'no_shadow':no_shadow(a), 'no_reflection':no_reflection(a) }
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_pyramid(a):
    pyramid_template = """
object {Pyramid
    pigment { color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f> }
    scale <%(scalex)f, %(scaley)f, %(scalez)f>
    rotate <%(rotx)f, %(roty)f, %(rotz)f>
    translate <%(transx)f, %(transy)f, %(transz)f>
    %(no_shadow)s
    %(no_reflection)s
}
"""
    alpha, theta, phi = find_rotations(a)
    object_code = pyramid_template % { 'scalex':displayscale*a.size[0],
                                       'scaley':displayscale*a.size[1],
                                       'scalez':displayscale*a.size[2],
                              'rotx':0., 'roty':-90.+theta*180./pi, 'rotz':phi*180./pi,
                              'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
                              'transx':displayscale*a.x, 'transy':displayscale*a.y, 'transz':displayscale*a.z,
                              'no_shadow':no_shadow(a), 'no_reflection':no_reflection(a) }
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export_arrow(a):
    al = a.length
    if displayscale*al < 0.01:       ## omit zero-length arrows
        return
    sw = a.shaftwidth
    if sw == 0:
        sw = 0.1*al
    hl = a.headlength
    if hl == 0:
        hl = 2*sw
    hw = a.headwidth
    if hw == 0:
        hw = 3*sw
    ## The following may have been needed only because there were other problems:
    ##    hw *= 0.8  ## seems too big, decrease this 2008-11-07 rwc
    if a.fixedwidth:
        if hl > .5*al:
            hl = .5*al
    else:
        if sw < .02*al:
            swtemp = .02*al
            hw *= swtemp/sw
            hl *= swtemp/sw
            sw = swtemp
        if hl > .5*al:
            hltemp = .5*al
            hw *= hltemp/hl
            sw *= hltemp/hl
            hl = hltemp
    sl = al-hl # length of shaft
    arrow_no_shadow = no_shadow(a)
    arrow_no_reflection = no_reflection(a)

    # head is a pyramid; need to create a dummy pyramid
    apyramid = pyramid(pos=a.pos+a.axis.norm()*sl, axis=a.axis, up=a.up, frame=a.frame,
               size=(hl,hw,hw), color=a.color, no_shadow=arrow_no_shadow,
               no_reflection=arrow_no_reflection, opacity=a.opacity, visible=0)
    m1 = export_pyramid(apyramid)
    m1 = add_texture(a, m1)
    del(apyramid)

    # shaft is a box; need to create a dummy box
    abox = box(pos=(a.pos+a.axis*(sl/al)/2.0), axis=(a.axis*(sl/al)),
               up = a.up, width=sw, height=sw,
               color=a.color, frame=a.frame, opacity=a.opacity,
               no_shadow=arrow_no_shadow, no_reflection=arrow_no_reflection, visible=0)
    m2 = export_box(abox)
    m2 = add_texture(a, m2)
    del(abox)
    # concatenate pyramid & box
    object_code = m1 + m2
    
    return object_code


def export_cone(a):
    cone_template = """
cone {
    <%(posx)f, %(posy)f, %(posz)f>, %(radius)f
    <%(pos2x)f, %(pos2y)f, %(pos2z)f>, %(minorradius)f
    pigment { color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f> }
    %(no_shadow)s
    %(no_reflection)s
}
"""
    pos2 = displayscale*(a.pos+a.axis)
    object_code = cone_template % { 'radius':displayscale*a.radius, 'minorradius':0.0,
                                    'posx':displayscale*a.x, 'posy':displayscale*a.y, 'posz':displayscale*a.z,
                                    'pos2x':pos2.x, 'pos2y':pos2.y, 'pos2z':pos2.z,
                                    'red':a.red, 'green':a.green, 'blue':a.blue, 'transparency':transparency(a),
                                    'no_shadow':no_shadow(a), 'no_reflection':no_reflection(a) }
    object_code = process_frame(a,object_code)
    object_code = add_texture(a,object_code)
    
    return object_code


def export_label(a):
    label_template = """
text {
    ttf "cyrvetic.ttf" "%(text)s" 0.1, 0
    pigment {color rgbt <%(red)f, %(green)f, %(blue)f, %(transparency)f>}
    translate <%(transx)f, %(transy)f, %(transz)f>
    %(no_shadow)s
    %(no_reflection)s
}
"""
    object_code = label_template % {'transx': displayscale * a.x,
                                    'transy': displayscale * a.y,
                                    'transz': displayscale * a.z,
                                    'red': a.red, 'green': a.green, 'blue': a.blue, 'transparency': transparency(a),
                                    'text': a.text,
                                    'no_shadow': no_shadow(a), 'no_reflection': no_reflection(a)}
    object_code = process_frame(a, object_code)
    object_code = add_texture(a, object_code)
    
    return object_code


def export(display=None, filename=None, include_list=None, xy_ratio=4./3., custom_text='', shadowless=0):
    '''
    This is essentially the function which does all the grunt work
    to convert the visual python scene into a .pov file.
    '''
    global displayscale
    if display == None:         # no display specified so find active display
        b = box(visible=0)
        display = b.display
        del b

    if filename == None:
        filename = display.title + '.pov'

    if include_list == None:
        include_text = ''
        # Maybe should always include the following definitions?
        #include_text = '#include "colors.inc"\n#include "stones.inc"\n#include "woods.inc"\n#include "metals.inc"\n'
    else:
        include_text = '\n'
        for x in include_list:
            include_text = include_text + '#include "'+ x + '"\n'

    initial_comment = """// povray code generated by povexport.py
"""

    pyramid_def = """
// Four-sided pyramid from shapes2.inc, slightly modified.
#declare Pyramid = union {
    triangle { <0, -1, -1>, <0, +1, -1>, <1, 0, 0> }
    triangle { <0, +1, -1>, <0, +1, +1>, <1, 0, 0> }
    triangle { <0, -1, +1>, <0, +1, +1>, <1, 0, 0> }
    triangle { <0, -1, +1>, <0, -1, -1>, <1, 0, 0> }
    triangle { <0, -1, -1>, <0, -1, +1>, <0, 1, 1> }
    triangle { <0, -1, -1>, <0, +1, -1>, <0, 1, 1> }
    scale <1, .5, .5>        // so height = width = 1
}
"""

    ambient_template = """
global_settings { ambient_light rgb <%(red)f, %(green)f, %(blue)f> }
"""

    scalar_ambient_template = """
global_settings { ambient_light rgb <%(amb)f, %(amb)f, %(amb)f> }
"""

    background_template = """
background { color rgb <%(red)f, %(green)f, %(blue)f> }
"""

    light_template = """
light_source { <%(posx)f, %(posy)f, %(posz)f>
    color rgb <%(red)f, %(green)f, %(blue)f>
}
"""

    camera_template = """
camera {
    right <-%(xyratio)f, 0, 0>      //visual uses right-handed coord. system
    location <%(posx)f, %(posy)f, %(posz)f>
    sky <%(upx)f, %(upy)f, %(upz)f>
    look_at <%(pos2x)f, %(pos2y)f, %(pos2z)f>
    angle %(fov)f
    rotate <%(rot1)s, %(rot2)s, %(rot3)s>
}
"""

    # begin povray setup
    file = open(filename, 'wt')

    file.write( initial_comment + include_text + custom_text + pyramid_def )
    if type(display.ambient) is tuple or type(display.ambient) is list:
        file.write( ambient_template % { 
            'red':display.ambient[0]*10 ,
            'green':display.ambient[1]*10,
            'blue':display.ambient[2]*10 })
    else:
        file.write( scalar_ambient_template % { 'amb':display.ambient*10 } )
        file.write( background_template % { 
            'red':display.background[0],
            'green':display.background[1],
            'blue':display.background[2] } )

    displayscale = 10.0/(mag(scene.mouse.camera-scene.center)*tan(scene.fov/2.))

    for light in display.lights: # reproduce visual lighting (not ideal, but good approximation)
        try:
            if type(light) is distant_light:
                #pos = norm(light.direction) * 1000.0 # far away to simulate parallel light
                pos = norm(light.direction) * 100.0 # far away to simulate parallel light
            elif type(light) is local_light: 
                pos = displayscale*light.pos
            lcolor = light.color
        except: # evidently it must be Visual 3
            #pos = norm(light) * 1000.0 # far away to simulate parallel light
            pos = norm(light) * 100.0 # far away to simulate parallel light
            lcolor = (mag(light),mag(light),mag(light)) # intensity of light (Visual 3 lights are white)
        light_code = light_template % { 'posx':pos.x, 'posy':pos.y, 'posz':pos.z,
                                        'red':lcolor[0], 'green':lcolor[1], 'blue':lcolor[2] }
        if shadowless:
            # insert frame_code at end (these rot's must be done last)
            end = light_code.rfind('}')
            light_code = light_code[:end] + '    shadowless\n' + light_code[end:]
        file.write( light_code )

    cpos = displayscale*display.mouse.camera
    ctr = displayscale*display.center
    cup = display.up
    file.write( camera_template % { 'xyratio':xy_ratio,
                                    'posx':cpos.x/cameraZoom, 
                                    'posy':cpos.y/cameraZoom, 
                                    'posz':cpos.z/cameraZoom,
                                    'upx':cup.x, 'upy':cup.y, 'upz':cup.z,
                                    'pos2x':ctr.x, 'pos2y':ctr.y, 'pos2z':ctr.z,
                                    'fov':display.fov*180/pi,
                                    'rot1':rotateCamera1,
                                    'rot2':rotateCamera2,
                                    'rot3':rotateCamera3 } )

    for obj in display.objects:
        key = obj.__class__
        if key in legal:
            obj_name = legal[key]
            if obj_name != 'frame':
                function_name = 'export_' + obj_name
                function = globals().get(function_name)
                object_code = function(obj)
                if object_code is None:
                    continue
##
##                print('object_code', object_code)
##                print('obj_name', obj_name)
##                
                file.write( object_code )
        else:
            print('WARNING: export function for ' + str(obj.__class__) + ' not implemented')

    file.close()
    return 'The export() function no longer returns the scene as an\n' \
           'endless POV-Ray string, but saves it to a file instead.'

if __name__ == '__main__':
    print(__doc__)
