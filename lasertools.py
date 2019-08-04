#!/usr/bin/env python
"""
Modified by Christoph Wagner 2017

based on gcodetools, https://github.com/cnc-club/gcodetools
based on inkscape-applytransforms, https://github.com/Klowner/inkscape-applytransforms

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""
import inkex
import simplestyle
import simplepath
import cubicsuperpath
import bezmisc
import simpletransform
from simpletransform import composeTransform, fuseTransform, parseTransform, applyTransformToPath, applyTransformToPoint, formatTransform

from multiprocessing import Pool

import cProfile
import sys
sys.path.append('/usr/share/inkscape/extensions')

import os
import math
import re
import time
import datetime
import cmath
import numpy as np
import gettext
_ = gettext.gettext

if "errormsg" not in dir(inkex):
    inkex.errormsg = lambda msg: sys.stderr.write(
        (str(msg) + "\n").encode("UTF-8"))


################################################################################
# Styles and additional parameters
################################################################################

gcode = ""

noOfThreads = 4
csp = []
profiling = False    # Disable if not debugging
debug = False        # Disable if not debugging

if profiling:
    import lsprofcalltree

timestamp = datetime.datetime.now()
math.pi2 = math.pi*2
tiny_infill_factor = 2  # x times the laser beam width will be removed
straight_tolerance = 0.000001
engraving_tolerance = 0.000002
options = {}
cspm = []
offset_y = 0
defaults = {
    'header': """
M03 S1
G90
""",
    'footer': """G00 X0 Y0
M05 S0
M02
"""
}

styles = {
    "biarc_style": {
        'line':        simplestyle.formatStyle({'stroke': '#f88', 'fill': 'none', "marker-end": 'none', 'stroke-width': '0.1'}),
        'biarc1':        simplestyle.formatStyle({'stroke': '#8f8', 'fill': 'none', "marker-end": 'none', 'stroke-width': '0.5'})
    }

}

################################################################################
# Cubic Super Path additional functions
################################################################################

'''
def checkIfLineInsideShape(splitted_line):
    # check if the middle point of the first lines segment is inside the path.
    # and remove the subline if not.
    l1, l2 = splitted_line[0], splitted_line[1]
    p = [(l1[0]+l2[0])/2, (l1[1]+l2[1])/2]
    # check for tangential points
    pMod = [(l1[0]+l2[0])/2, ((l1[1]+l2[1])/2) + 0.1]

    if l1 != l2:
        if point_inside_csp(p, csp):
            if point_inside_csp(pMod, csp):
                if len(splitted_line) == 2:
                    return splitted_line
                else:
                    return [splitted_line[0], splitted_line[1]]

            else:
                print_debug("splitted_line removed: ", splitted_line)
                return [[0, 0], [0, 0]]
        else:
            return [[0, 0], [0, 0]]
    else:
        return [[0, 0], [0, 0]]


'''


def checkIfLineInsideShape(splitted_line):
    # print_("sl input", splitted_line)
    # check if the middle point of the first lines segment is inside the path.
    # and remove the subline if not.
    l1, l2 = splitted_line[0], splitted_line[1]

    if l1 == l2 and len(splitted_line) > 2:
        l2 = splitted_line[2]

    p = [(l1[0]+l2[0])/2, (l1[1]+l2[1])/2]

    if point_inside_csp(p, csp):
        if len(splitted_line) > 2:
            print_debug("len splitted line > 2", splitted_line)

        # if l1 == l2 and len(splitted_line) > 2:
        #    l2 = splitted_line[2]

        if l1 != l2:
            print_debug("splitted line ", [l1, l2])
            return [l1, l2]
        else:
            return [[0, 0], [0, 0]]
    else:
        return [[0, 0], [0, 0]]


def csp_segment_to_bez(sp1, sp2):
    return sp1[1:]+sp2[:2]


def csp_true_bounds(csp):

    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx = [float("inf"), 0, 0, 0]
    maxx = [float("-inf"), 0, 0, 0]
    miny = [float("inf"), 0, 0, 0]
    maxy = [float("-inf"), 0, 0, 0]
    for i in range(len(csp)):
        for j in range(1, len(csp[i])):
            ax, ay, bx, by, cx, cy, x0, y0 = bezmisc.bezierparameterize(
                (csp[i][j-1][1], csp[i][j-1][2], csp[i][j][0], csp[i][j][1]))
            roots = cubic_solver(0, 3*ax, 2*bx, cx) + [0, 1]
            for root in roots:
                if type(root) is complex and abs(root.imag) < 1e-10:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxx = max([x, y, i, j, root], maxx)
                    minx = min([x, y, i, j, root], minx)

            roots = cubic_solver(0, 3*ay, 2*by, cy) + [0, 1]
            for root in roots:
                if type(root) is complex and root.imag == 0:
                    root = root.real
                if type(root) is not complex and 0 <= root <= 1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxy = max([y, x, i, j, root], maxy)
                    miny = min([y, x, i, j, root], miny)
    maxy[0], maxy[1] = maxy[1], maxy[0]
    miny[0], miny[1] = miny[1], miny[0]

    return minx, miny, maxx, maxy


def csp_at_t(sp1, sp2, t):
    ax, bx, cx, dx = sp1[1][0], sp1[2][0], sp2[0][0], sp2[1][0]
    ay, by, cy, dy = sp1[1][1], sp1[2][1], sp2[0][1], sp2[1][1]

    x1, y1 = ax+(bx-ax)*t, ay+(by-ay)*t
    x2, y2 = bx+(cx-bx)*t, by+(cy-by)*t
    x3, y3 = cx+(dx-cx)*t, cy+(dy-cy)*t

    x4, y4 = x1+(x2-x1)*t, y1+(y2-y1)*t
    x5, y5 = x2+(x3-x2)*t, y2+(y3-y2)*t

    x, y = x4+(x5-x4)*t, y4+(y5-y4)*t
    return [x, y]


def cspseglength(sp1, sp2, tolerance=0.0001):
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    return bezmisc.bezierlength(bez, tolerance)


def csp_line_intersection(l1, l2, sp1, sp2):
    dd = l1[0]
    cc = l2[0]-l1[0]
    bb = l1[1]
    aa = l2[1]-l1[1]
    if aa == cc == 0:
        return []
    if aa:
        coef1 = cc/aa
        coef2 = 1
    else:
        coef1 = 1
        coef2 = aa/cc
    bez = (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:])
    ax, ay, bx, by, cx, cy, x0, y0 = bezmisc.bezierparameterize(bez)
    a = coef1*ay-coef2*ax
    b = coef1*by-coef2*bx
    c = coef1*cy-coef2*cx
    d = coef1*(y0-bb)-coef2*(x0-dd)
    roots = cubic_solver(a, b, c, d)
    retval = []
    for i in roots:
        if type(i) is complex and abs(i.imag) < 1e-7:
            i = i.real
        if type(i) is not complex and -1e-10 <= i <= 1.+1e-10:
            retval.append(i)
    return retval


# Return only points [ (x,y) ]
def line_line_intersection_points(p1, p2, p3, p4):
    if (p1[0] == p2[0] and p1[1] == p2[1]) or (p3[0] == p4[0] and p3[1] == p4[1]):
        return []
    x = (p2[0]-p1[0])*(p4[1]-p3[1]) - (p2[1]-p1[1])*(p4[0]-p3[0])
    if x == 0:  # Lines are parallel
        if (p3[0]-p1[0])*(p2[1]-p1[1]) == (p3[1]-p1[1])*(p2[0]-p1[0]):
            if p3[0] != p4[0]:
                t11 = (p1[0]-p3[0])/(p4[0]-p3[0])
                t12 = (p2[0]-p3[0])/(p4[0]-p3[0])
                t21 = (p3[0]-p1[0])/(p2[0]-p1[0])
                t22 = (p4[0]-p1[0])/(p2[0]-p1[0])
            else:
                t11 = (p1[1]-p3[1])/(p4[1]-p3[1])
                t12 = (p2[1]-p3[1])/(p4[1]-p3[1])
                t21 = (p3[1]-p1[1])/(p2[1]-p1[1])
                t22 = (p4[1]-p1[1])/(p2[1]-p1[1])
            res = []
            if (0 <= t11 <= 1 or 0 <= t12 <= 1) and (0 <= t21 <= 1 or 0 <= t22 <= 1):
                if 0 <= t11 <= 1:
                    res += [p1]
                if 0 <= t12 <= 1:
                    res += [p2]
                if 0 <= t21 <= 1:
                    res += [p3]
                if 0 <= t22 <= 1:
                    res += [p4]
            return res
        else:
            return []
    else:
        t1 = ((p4[0]-p3[0])*(p1[1]-p3[1]) - (p4[1]-p3[1])*(p1[0]-p3[0]))/x
        t2 = ((p2[0]-p1[0])*(p1[1]-p3[1]) - (p2[1]-p1[1])*(p1[0]-p3[0]))/x
        if 0 <= t1 <= 1 and 0 <= t2 <= 1:
            return [[p1[0]*(1-t1)+p2[0]*t1, p1[1]*(1-t1)+p2[1]*t1]]
        else:
            return []


def point_to_point_d2(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2


def point_to_point_d(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)


def csp_normalized_slope(sp1, sp2, t):
    ax, ay, bx, by, cx, cy = bezmisc.bezierparameterize(
        (sp1[1][:], sp1[2][:], sp2[0][:], sp2[1][:]))[0:6]
    if sp1[1] == sp2[1] == sp1[2] == sp2[0]:
        return [1., 0.]
    f1x = 3*ax*t*t+2*bx*t+cx
    f1y = 3*ay*t*t+2*by*t+cy
    if abs(f1x*f1x+f1y*f1y) > 1e-20:
        l = math.sqrt(f1x*f1x+f1y*f1y)
        return [f1x/l, f1y/l]

    if t == 0:
        f1x = sp2[0][0]-sp1[1][0]
        f1y = sp2[0][1]-sp1[1][1]
        if abs(f1x*f1x+f1y*f1y) > 1e-20:
            l = math.sqrt(f1x*f1x+f1y*f1y)
            return [f1x/l, f1y/l]
        else:
            f1x = sp2[1][0]-sp1[1][0]
            f1y = sp2[1][1]-sp1[1][1]
            if f1x*f1x+f1y*f1y != 0:
                l = math.sqrt(f1x*f1x+f1y*f1y)
                return [f1x/l, f1y/l]
    elif t == 1:
        f1x = sp2[1][0]-sp1[2][0]
        f1y = sp2[1][1]-sp1[2][1]
        if abs(f1x*f1x+f1y*f1y) > 1e-20:
            l = math.sqrt(f1x*f1x+f1y*f1y)
            return [f1x/l, f1y/l]
        else:
            f1x = sp2[1][0]-sp1[1][0]
            f1y = sp2[1][1]-sp1[1][1]
            if f1x*f1x+f1y*f1y != 0:
                l = math.sqrt(f1x*f1x+f1y*f1y)
                return [f1x/l, f1y/l]
    else:
        return [1., 0.]


def csp_normalized_normal(sp1, sp2, t):
    nx, ny = csp_normalized_slope(sp1, sp2, t)
    return [-ny, nx]


def csp_parameterize(sp1, sp2):
    return bezmisc.bezierparameterize(csp_segment_to_bez(sp1, sp2))


def csp_draw(csp, color="#05f", group=None, style="fill:none;", width=.1, comment=""):
    if csp != [] and csp != [[]]:
        if group == None:
            group = options.doc_root
        style += "stroke:"+color+";" + "stroke-width:%0.4fpx;" % width
        args = {"d": cubicsuperpath.formatPath(csp), "style": style}
        if comment != "":
            args["comment"] = str(comment)
        inkex.etree.SubElement(group, inkex.addNS('path', 'svg'), args)


def csp_subpath_line_to(subpath, points):
    # Appends subpath with line or polyline.
    if len(points) > 0:
        if len(subpath) > 0:
            subpath[-1][2] = subpath[-1][1][:]
        if type(points[0]) == type([1, 1]):
            for p in points:
                subpath += [[p[:], p[:], p[:]]]
        else:
            subpath += [[points, points, points]]
    return subpath


################################################################################
# Area Fill Stuff
################################################################################

def point_inside_csp(p, csp, on_the_path=True):

    x, y = p
    ray_intersections_count = 0

    for subpath in csp:
        for i in range(1, len(subpath)):
            sp1, sp2 = subpath[i-1], subpath[i]
            ax, bx, cx, dx = csp_parameterize(sp1, sp2)[::2]
            if ax == 0 and bx == 0 and cx == 0 and dx == x:
                # we've got a special case here
                b = csp_true_bounds([[sp1, sp2]])
                if b[1][1] <= y <= b[3][1]:
                    # points is on the path
                    return on_the_path
                else:
                    # we can skip this segment because it wont influence the answer.
                    pass
            else:
                for t in csp_line_intersection([x, y], [x, y+5], sp1, sp2):

                    if t == 0 or t == 1:
                        # we've got another special case here
                        y1 = csp_at_t(sp1, sp2, t)[1]

                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        # if t == 0 we sould have considered this case previously.

                        if t == 1:
                            # we have to check the next segmant if it is on the same side of the ray
                            st_d = csp_normalized_slope(sp1, sp2, 1)[0]

                            if st_d == 0:
                                st_d = csp_normalized_slope(sp1, sp2, 0.99)[0]

                            for j in range(1, len(subpath)+1):

                                if (i+j) % len(subpath) == 0:
                                    continue  # skip the closing segment

                                sp11, sp22 = subpath[(
                                    i-1+j) % len(subpath)], subpath[(i+j) % len(subpath)]
                                ax1, bx1, cx1, dx1 = csp_parameterize(
                                    sp1, sp2)[::2]

                                if ax1 == 0 and bx1 == 0 and cx1 == 0 and dx1 == x:
                                    continue  # this segment parallel to the ray, so skip it
                                en_d = csp_normalized_slope(sp11, sp22, 0)[0]
                                if en_d == 0:
                                    en_d = csp_normalized_slope(
                                        sp11, sp22, 0.01)[0]
                                if st_d*en_d <= 0:
                                    ray_intersections_count += 1
                                    break
                    else:
                        y1 = csp_at_t(sp1, sp2, t)[1]

                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        else:
                            if y1 > y and 3*ax*t**2 + 2*bx*t + cx != 0:  # if it's 0 the path only touches the ray
                                ray_intersections_count += 1

    return ray_intersections_count % 2 == 1


def csp_from_polyline(line):
    return [[[point[:] for _ in range(3)] for point in subline] for subline in line]


'''
def csp_close_all_subpaths(csp, tolerance=0.000001):
    for i in range(len(csp)):
        if point_to_point_d2(csp[i][0][1], csp[i][-1][1]) > tolerance**2:
            csp[i][-1][2] = csp[i][-1][1][:]
            csp[i] += [[csp[i][0][1][:] for _ in range(3)]]
        else:
            if csp[i][0][1] != csp[i][-1][1]:
                csp[i][-1][1] = csp[i][0][1][:]
    return csp
'''

################################################################################
# Some vector functions
################################################################################


def normalize(x, y):
    l = math.sqrt(x**2+y**2)
    if l == 0:
        return [0., 0.]
    else:
        return [x/l, y/l]


def cross(a, b):
    return a[1] * b[0] - a[0] * b[1]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


################################################################################
# Common functions
################################################################################

def atan2(*arg):
    if len(arg) == 1 and (type(arg[0]) == type([0., 0.]) or type(arg[0]) == type((0., 0.))):
        return (math.pi/2 - math.atan2(arg[0][0], arg[0][1])) % math.pi2
    elif len(arg) == 2:

        return (math.pi/2 - math.atan2(arg[0], arg[1])) % math.pi2
    else:
        raise ValueError("Bad argumets for atan! (%s)" % arg)


def between(c, x, y):
    return x-straight_tolerance <= c <= y+straight_tolerance or y-straight_tolerance <= c <= x+straight_tolerance


def cubic_solver(a, b, c, d):
    if a != 0:
        #    Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
        a, b, c = (b/a, c/a, d/a)
        m = 2*a**3 - 9*a*b + 27*c
        k = a**2 - 3*b
        n = m**2 - 4*k**3
        w1 = -.5 + .5*cmath.sqrt(3)*1j
        w2 = -.5 - .5*cmath.sqrt(3)*1j
        if n >= 0:
            t = m+math.sqrt(n)
            m1 = pow(t/2, 1./3) if t >= 0 else -pow(-t/2, 1./3)
            t = m-math.sqrt(n)
            n1 = pow(t/2, 1./3) if t >= 0 else -pow(-t/2, 1./3)
        else:
            m1 = pow(complex((m+cmath.sqrt(n))/2), 1./3)
            n1 = pow(complex((m-cmath.sqrt(n))/2), 1./3)
        x1 = -1./3 * (a + m1 + n1)
        x2 = -1./3 * (a + w1*m1 + w2*n1)
        x3 = -1./3 * (a + w2*m1 + w1*n1)
        return [x1, x2, x3]
    elif b != 0:
        det = c**2-4*b*d
        if det > 0:
            return [(-c+math.sqrt(det))/(2*b), (-c-math.sqrt(det))/(2*b)]
        elif d == 0:
            return [-c/(b*b)]
        else:
            return [(-c+cmath.sqrt(det))/(2*b), (-c-cmath.sqrt(det))/(2*b)]
    elif c != 0:
        return [-d/c]
    else:
        return []


################################################################################
# print_ prints any arguments into specified log file
################################################################################

def print_(*arg):
    if (os.path.isdir(options.directory)):
        f = open(options.directory+"/log.txt", "a")
        for s in arg:
            s = str(str(s).encode('unicode_escape'))+" "
            f.write(s)
        f.write("\n")
        f.close()


def print_debug(*arg):
    if debug:
        print_("DEBUG: ", *arg)


def print_time(*arg):
    global timestamp

    time = datetime.datetime.now() - timestamp
    print_(time, "   ", *arg)
    timestamp = datetime.datetime.now()

################################################################################
# Point (x,y) operations
################################################################################


class P:
    def __init__(self, x, y=None):
        if not y == None:
            self.x, self.y = float(x), float(y)
        else:
            self.x, self.y = float(x[0]), float(x[1])

    def __add__(self, other): return P(self.x + other.x, self.y + other.y)

    def __sub__(self, other): return P(self.x - other.x, self.y - other.y)

    def __neg__(self): return P(-self.x, -self.y)

    def __mul__(self, other):
        if isinstance(other, P):
            return self.x * other.x + self.y * other.y
        return P(self.x * other, self.y * other)
    __rmul__ = __mul__

    def __div__(self, other): return P(self.x / other, self.y / other)

    def mag(self): return math.hypot(self.x, self.y)

    def unit(self):
        h = self.mag()
        if h:
            return self / h
        else:
            return P(0, 0)

    def dot(self, other): return self.x * other.x + self.y * other.y

    def rot(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        return P(self.x * c - self.y * s,  self.x * s + self.y * c)

    def angle(self): return math.atan2(self.y, self.x)

    def __repr__(self): return '%f,%f' % (self.x, self.y)

    def pr(self): return "%.3f,%.3f" % (self.x, self.y)

    def to_list(self): return [self.x, self.y]

    def ccw(self): return P(-self.y, self.x)

    def l2(self): return self.x*self.x + self.y*self.y

################################################################################
# Gcodetools class
################################################################################


class laser_gcode(inkex.Effect):

    def export_gcode(self, gcode):
        gcode_pass = gcode
        for _ in range(1, self.options.passes):
            gcode += "G91\n" + "\nG90\n" + gcode_pass
        f = open(self.options.directory+self.options.file, "w")
        f.write(self.header + "\n" + gcode + self.footer)
        f.close()

    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-d", "--directory",                       action="store", type="string",
                                     dest="directory",                           default="/insert your target directory here",      help="Output directory")
        self.OptionParser.add_option("-f", "--filename",                        action="store", type="string",
                                     dest="file",                                default="output.ngc",                   help="File name")
        self.OptionParser.add_option("",   "--add-numeric-suffix-to-filename",  action="store", type="inkbool",
                                     dest="add_numeric_suffix_to_filename",      default=False,                          help="Add numeric suffix to file name")
        self.OptionParser.add_option("",   "--laser-command-perimeter",         action="store", type="string",
                                     dest="laser_command_perimeter",             default="S100",                         help="Laser gcode command Perimeter")
        self.OptionParser.add_option("",   "--laser-command",                   action="store", type="string",
                                     dest="laser_command",                       default="S100",                         help="Laser gcode command infill")
        self.OptionParser.add_option("",   "--laser-off-command",               action="store", type="string",
                                     dest="laser_off_command",                   default="S1",                           help="Laser gcode end command")
        self.OptionParser.add_option("",   "--laser-beam-with",                 action="store", type="float",
                                     dest="laser_beam_with",                     default="0.3",                          help="Laser speed (mm/min)")
        self.OptionParser.add_option("",   "--laser-speed",                     action="store", type="int",
                                     dest="laser_speed",                         default="1200",                          help="Laser speed for infill (mm/min)")
        self.OptionParser.add_option("",   "--laser-param-speed",               action="store", type="int",
                                     dest="laser_param_speed",                         default="700",                          help="Laser speed for Parameter (mm/min)")
        self.OptionParser.add_option("",   "--passes",                          action="store", type="int",
                                     dest="passes",                              default="1",                            help="Quantity of passes")
        self.OptionParser.add_option("",   "--power-delay",                     action="store", type="string",
                                     dest="power_delay",                         default="0",                            help="Laser power-on delay (ms)")
        self.OptionParser.add_option("",   "--suppress-all-messages",           action="store", type="inkbool",
                                     dest="suppress_all_messages",               default=True,                           help="Hide messages during g-code generation")
        self.OptionParser.add_option("",   "--create-log",                      action="store", type="inkbool",
                                     dest="log_create_log",                      default=True,                           help="Create log files")
        self.OptionParser.add_option("",   "--engraving-draw-calculation-paths", action="store", type="inkbool",
                                     dest="engraving_draw_calculation_paths",    default=True,                           help="Draw additional graphics to debug engraving path")
        self.OptionParser.add_option("",   "--biarc-max-split-depth",           action="store", type="int",             dest="biarc_max_split_depth",
                                     default="2",                           help="Defines maximum depth of splitting while approximating using biarcs.")
        self.OptionParser.add_option("",   "--area-fill-angle",                 action="store", type="float",           dest="area_fill_angle",
                                     default="0",                            help="Fill area with lines heading this angle")
        self.OptionParser.add_option("",   "--engraving-newton-iterations",     action="store", type="int",             dest="engraving_newton_iterations",
                                     default="20",                           help="Number of sample points used to calculate distance")
        self.OptionParser.add_option("",   "--add-contours",                    action="store", type="inkbool",
                                     dest="add_contours",                        default=True,                           help="Add contour to Gcode paths")
        self.OptionParser.add_option("",   "--add-infill",                      action="store", type="inkbool",
                                     dest="add_infill",                          default=True,                           help="Add infill to Gcode paths")
        self.OptionParser.add_option("",   "--remove-tiny-infill-paths",        action="store", type="inkbool",
                                     dest="remove_tiny_infill_paths",            default=False,                           help="Remove tiny infill paths from Gcode")
        self.OptionParser.add_option("",   "--multi_thread",                      action="store", type="inkbool",
                                     dest="multi_thread",                          default=True,                           help="Activate multithreading support")

    def parse_curve(self, p, layer):

        c = []
        if len(p) == 0:
            return []

        p = self.transform_csp(p, layer)

        # this code is intended to replace the code below.
        # At the Moment there is a problem with muliple paths, where the fist/last path will not be generated
        # TODO: Fix that
        '''
        print_("p_post_Trans ", p)
        startPoints = np.zeros(shape=[len(p), 2])
        endPoints = np.zeros(shape=[len(p), 2])

        # reduce Array size
        for i in range(0, len(p)):
            elRed = np.array(p[i])
            elRed = elRed[:, 0]
            startPoints[i] = elRed[0]
            endPoints[i] = elRed[-1]

        print_("StartPoints", startPoints)
        print_("EndPoints", endPoints)
        print_("elRed", elRed)

        sortedPoints = np.array(self.sort_points(
            startPoints[:, 0], startPoints[:, 1], endPoints[:, 0], endPoints[:, 1]))

        for point in sortedPoints:
            ind = np.argwhere(startPoints == point[0])[0, 0]
            elRed = np.array(p[ind])
            elRed = elRed[:, 0]

            c += [[[point[0], point[1]], 'move']]

            for i in range(1, len(elRed)):
                c += [[[elRed[i-1, 0], elRed[i-1, 1]], 'line', [elRed[i, 0], elRed[i, 1]]]]

            c += [[[elRed[-1, 0], elRed[-1, 1]], 'end']]

        # Sort to reduce Rapid distance
        '''

        k = range(1, len(p))
        keys = [0]
        while len(k) > 0:
            end = p[keys[-1]][-1][1]
            dist = None
            for i in range(len(k)):
                start = p[k[i]][0][1]
                dist = max((-((end[0]-start[0])**2+(end[1]-start[1])**2), i),   dist)
            keys += [k[dist[1]]]
            del k[dist[1]]
        for k in keys:
            subpath = p[k]
            c += [[[subpath[0][1][0], subpath[0][1][1]], 'move']]
            # print_([[[subpath[0][1][0], subpath[0][1][1]], 'move']])
            for i in range(1, len(subpath)):
                # print_("subpath: ", subpath[i-1])
                sp1 = [[subpath[i-1][j][0], subpath[i-1][j][1]]
                       for j in range(3)]
                sp2 = [[subpath[i][j][0], subpath[i][j][1]] for j in range(3)]
                c += [[[sp1[0][0], sp1[0][1]], 'line', [sp2[0][0], sp2[0][1]]]]
                # print_("sp1: ", sp1)
            c += [[[subpath[-1][1][0], subpath[-1][1][1]], 'end']]
            # print_([[[subpath[-1][1][0], subpath[-1][1][1]], 'end']])
        # print_("Curve: " + str(c))

        return c

    def parse_curve2d(self, p, layer):

        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)
        # print_("p: ", p)
        np_p = np.array(p)
        # print_("len p Slice: ", len(np_p[:, 0, 0, 0]))
        sortedToolpaths = self.sort_points(np_p[:, 0, 0, 0], np_p[:, 0, 0, 1], np_p[:, 1, 0, 0], np_p[:, 1, 0, 1])

        for path in sortedToolpaths:
            c += [[[path[0], path[1]], 'move']]
            c += [[[path[0], path[1]], 'line', [path[2], path[3]]]]
            c += [[[path[2], path[3]], 'end']]

        return c

    def sort_points(self, x1, y1, x2, y2):
        sortedList = np.zeros((len(x1), 4))

        xpos = np.array(x1[1:])
        xposInv = np.array(x2[1:])
        ypos = np.array(y1[1:])
        yposInv = np.array(y2[1:])

        sortedList[0] = [x1[0], y1[0], x2[0], y2[0]]
        actXPos = x2[0]
        actYPos = y2[0]

        i = 1

        while len(xpos) > 0:
            xDist = np.abs(xpos - actXPos)
            xDistInv = np.abs(xposInv - actXPos)
            yDist = np.abs(ypos - actYPos)
            yDistInv = np.abs(yposInv - actYPos)

            distances = np.array([np.add(xDist, yDist), np.add(xDistInv, yDistInv)])
            distances = np.abs(distances)
            # print_("Distances: ", distances)

            minInd = np.unravel_index(np.argmin(distances, axis=None), distances.shape)

            if minInd[0] == 0:
                sortedList[i] = [xpos[minInd[1]], ypos[minInd[1]], xposInv[minInd[1]], yposInv[minInd[1]]]
                actXPos = xposInv[minInd[1]]
                actYPos = yposInv[minInd[1]]

            else:
                sortedList[i] = [xposInv[minInd[1]], yposInv[minInd[1]], xpos[minInd[1]], ypos[minInd[1]]]
                actXPos = xpos[minInd[1]]
                actYPos = ypos[minInd[1]]

            xpos = np.delete(xpos, minInd[1])
            ypos = np.delete(ypos, minInd[1])
            xposInv = np.delete(xposInv, minInd[1])
            yposInv = np.delete(yposInv, minInd[1])

            i = i+1

        return sortedList

    def draw_curve(self, curve, layer, group=None, style=styles["biarc_style"]):

        self.get_defs()

        if group == None:
            group = inkex.etree.SubElement(self.layers[min(1, len(
                self.layers)-1)], inkex.addNS('g', 'svg'), {"gcodetools": "Preview group"})

        s = ''
        a, b, c = [0., 0.], [1., 0.], [0., 1.]
        a, b, c = self.transform(a, layer, True), self.transform(b, layer, True), self.transform(c, layer, True)

        for si in curve:

            si[0] = self.transform(si[0], layer, True)
            # print_("si ", si)
            if len(si) == 3:
                si[2] = self.transform(si[2], layer, True)

            if s != '':
                if s[1] == 'line':
                    inkex.etree.SubElement(group, inkex.addNS('path', 'svg'),
                                           {
                        'style': style['line'],
                        'd': 'M %s,%s L %s,%s' % (s[0][0], s[0][1], si[0][0], si[0][1]),
                        "gcodetools": "Preview",
                    }
                    )
            s = si

    def check_dir(self):
        if self.options.directory[-1] not in ["/", "\\"]:
            if "\\" in self.options.directory:
                self.options.directory += "\\"
            else:
                self.options.directory += "/"
        print_("Checking direcrory: '%s'" % self.options.directory)
        if (os.path.isdir(self.options.directory)):
            if (os.path.isfile(self.options.directory+'header')):
                f = open(self.options.directory+'header', 'r')
                self.header = f.read()
                f.close()
            else:
                self.header = defaults['header']
            if (os.path.isfile(self.options.directory+'footer')):
                f = open(self.options.directory+'footer', 'r')
                self.footer = f.read()
                f.close()
            else:
                self.footer = defaults['footer']

            self.header += "G21\n"

        else:
            self.error(_("Directory does not exist! Please specify existing directory at options tab!"), "error")
            return False

        if self.options.add_numeric_suffix_to_filename:
            dir_list = os.listdir(self.options.directory)
            if "." in self.options.file:
                r = re.match(r"^(.*)(\..*)$", self.options.file)
                ext = r.group(2)
                name = r.group(1)
            else:
                ext = ""
                name = self.options.file
            max_n = 0
            for s in dir_list:
                r = re.match(r"^%s_0*(\d+)%s$" %
                             (re.escape(name), re.escape(ext)), s)
                if r:
                    max_n = max(max_n, int(r.group(1)))
            filename = name + "_" + \
                ("0"*(4-len(str(max_n+1))) + str(max_n+1)) + ext
            self.options.file = filename

        print_("Testing writing rights on '%s'" %
               (self.options.directory+self.options.file))
        try:
            f = open(self.options.directory+self.options.file, "w")
            f.close()
        except:
            self.error(_("Can not write to specified file!\n%s" % (self.options.directory+self.options.file)), "error")
            return False
        return True


################################################################################
# Generate Gcode
# Generates Gcode on given curve.
# Crve defenitnion [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
################################################################################

    def generate_gcode(self, curve, layer, tool):
        global doc_height
        global offset_y

        print_debug("    Laser parameters: " + str(tool))

        def c(c):
            c = [c[i] if i < len(c) else None for i in range(6)]
            if c[5] == 0:
                c[5] = None
            s = [" X", " Y", " Z", " I", " J", " K"]
            r = ''
            for i in range(6):
                if c[i] != None:
                    r += s[i] + ("%f" % (round(c[i], 2))).rstrip('0')
            return r

        try:
            self.last_used_tool == None
        except:
            self.last_used_tool = None

        # print_("Curve: " + str(curve) + "/n")
        g = ""
        lg, f = 'G00', "F%.1f" % tool['penetration feed']

        for i in range(1, len(curve)):
            #    Creating Gcode for curve between s=curve[i-1] and si=curve[i] start at s[0] end at s[4]=si[0]
            s, si = curve[i-1], curve[i]

            si[0] = self.transform(si[0], layer, True)
            # print_("si ", si)
            if len(si) == 3:
                si[2] = self.transform(si[2], layer, True)

            s[0][1] = s[0][1] - offset_y
            si[0][1] = si[0][1] - offset_y

            feed = f if lg not in ['G01', 'G02', 'G03'] else ''
            if s[1] == 'move':
                g += "G00" + c(si[0]) + "\n" + tool['gcode before path'] + "\n"
                lg = 'G00'
            elif s[1] == 'line':
                if lg == "G00":
                    g += "G01 " + feed + "\n"
                g += "G01" + c(si[0]) + "\n"
                lg = 'G01'

            if si[1] == 'end':
                g += tool['gcode after path'] + "\n"

        return g

        #     elif s[1] == 'end':
        # g += tool['gcode after path'] + "\n"
        # lg = 'G00'

    def get_transforms(self, g):
        root = self.document.getroot()
        trans = []
        while (g != root):
            if 'transform' in g.keys():
                t = g.get('transform')
                t = simpletransform.parseTransform(t)
                trans = simpletransform.composeTransform(
                    t, trans) if trans != [] else t
            g = g.getparent()
        return trans

    def apply_transforms(self, g, csp):
        trans = self.get_transforms(g)
        if trans != []:
            trans[1][2] = 0
            print_('    applying transformation')
            # print_(trans)
            simpletransform.applyTransformToPath(trans, csp)
        return csp

    def transform(self, source_point, layer, reverse=False):
        if layer == None:
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer not in self.transform_matrix:
            for i in range(self.layers.index(layer), -1, -1):
                if self.layers[i] in self.orientation_points:
                    break

            # print_(str(self.layers))
            # print_(str("I: " + str(i)))
            if self.layers[i] not in self.orientation_points:
                self.error(_("Orientation points for '%s' layer have not been found! Please add orientation points using Orientation tab!") % layer.get(
                    inkex.addNS('label', 'inkscape')))
            elif self.layers[i] in self.transform_matrix:
                self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
            else:
                orientation_layer = self.layers[i]
                if len(self.orientation_points[orientation_layer]) > 1:
                    self.error(_("There are more than one orientation point groups in '%s' layer") % orientation_layer.get(
                        inkex.addNS('label', 'inkscape')))
                points = self.orientation_points[orientation_layer][0]
                if len(points) == 2:
                    points += [[[(points[1][0][1]-points[0][0][1])+points[0][0][0], -(points[1][0][0]-points[0][0][0])+points[0][0][1]],
                                [-(points[1][1][1]-points[0][1][1])+points[0][1][0], points[1][1][0]-points[0][1][0]+points[0][1][1]]]]
                if len(points) == 3:
                    # print_("Layer '%s' Orientation points: " % orientation_layer.get(inkex.addNS('label', 'inkscape')))
                    # for point in points:
                    #    print_(point)
                    #    Zcoordinates definition taken from Orientatnion point 1 and 2
                    self.Zcoordinates[layer] = [
                        max(points[0][1][2], points[1][1][2]), min(points[0][1][2], points[1][1][2])]
                    matrix = np.array([
                        [points[0][0][0], points[0][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[0][0][0],
                         points[0][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[0][0]
                         [0], points[0][0][1], 1],
                        [points[1][0][0], points[1][0]
                         [1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[1][0][0],
                         points[1][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[1][0]
                         [0], points[1][0][1], 1],
                        [points[2][0][0], points[2][0]
                         [1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[2][0][0],
                         points[2][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[2]
                         [0][0], points[2][0][1], 1]
                    ])

                    if np.linalg.det(matrix) != 0:
                        m = np.linalg.solve(matrix,
                                            np.array(
                                                [[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]], [
                                                    points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
                                            )
                                            ).tolist()
                        self.transform_matrix[layer] = [
                            [m[j*3+i][0] for i in range(3)] for j in range(3)]

                    else:
                        self.error(
                            _("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"))
                else:
                    self.error(_("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"))

            self.transform_matrix_reverse[layer] = np.linalg.inv(
                self.transform_matrix[layer]).tolist()
            # print_(self.transform_matrix)
            # print_(self.transform_matrix_reverse)

        x, y = source_point[0], source_point[1]
        # print_("source point x", str(x) + " " + str(y))

        if not reverse:
            t = self.transform_matrix[layer]
        else:
            t = self.transform_matrix_reverse[layer]
        return [t[0][0]*x+t[0][1]*y+t[0][2], t[1][0]*x+t[1][1]*y+t[1][2]]

    def transform_csp(self, csp_, layer, reverse=False):
        csp = [[[csp_[i][j][0][:], csp_[i][j][1][:], csp_[i][j][2][:]]
                for j in range(len(csp_[i]))] for i in range(len(csp_))]

        for i in xrange(len(csp)):

            for j in xrange(len(csp[i])):
                # print_("csp pre trans", csp[i][j])
                for k in xrange(len(csp[i][j])):
                    csp[i][j][k] = self.transform(csp[i][j][k], layer, reverse)
                # csp[i][j][0] = self.transform(csp[i][j][0], layer, reverse)
                # print_("csp post trans", csp[i][j])
        return csp

################################################################################
# Errors handling function, notes are just printed into Logfile,
# warnings are printed into log file and warning message is displayed but
# extension continues working, errors causes log and execution is halted
# Notes, warnings adn errors could be assigned to space or comma or dot
# sepparated strings (case is ignoreg).
################################################################################
    def error(self, s, type_="Warning"):

        if type_ == "Warning":
            print_(s)
            inkex.errormsg(s+"\n")
        else:
            print_(s)
            inkex.errormsg(s)
            sys.exit()


################################################################################
# Get defs from svg
################################################################################
    def get_defs(self):
        self.defs = {}

        def recursive(g):
            for i in g:
                if i.tag == inkex.addNS("defs", "svg"):
                    for j in i:
                        self.defs[j.get("id")] = i
                if i.tag == inkex.addNS("g", 'svg'):
                    recursive(i)
        recursive(self.document.getroot())


################################################################################
# Get Gcodetools info from the svg
################################################################################
    def get_info(self):
        self.selected_paths = {}
        self.paths = {}
        self.orientation_points = {}
        self.layers = [self.document.getroot()]
        self.Zcoordinates = {}
        self.transform_matrix = {}
        self.transform_matrix_reverse = {}

        def recursive_search(g, layer, selected=False):
            items = g.getchildren()
            items.reverse()
            for i in items:
                if selected:
                    self.selected[i.get("id")] = i
                if i.tag == inkex.addNS("g", 'svg') and i.get(inkex.addNS('groupmode', 'inkscape')) == 'layer':
                    self.layers += [i]
                    recursive_search(i, i)
                elif i.get('gcodetools') == "Gcodetools orientation group":
                    points = self.get_orientation_points(i)
                    if points != None:
                        self.orientation_points[layer] = self.orientation_points[layer]+[
                            points[:]] if layer in self.orientation_points else [points[:]]
                        print_("    Found orientation points in '%s' layer: %s" % (
                            layer.get(inkex.addNS('label', 'inkscape')), points))
                    else:
                        self.error(_("Warning! Found bad orientation points in '%s' layer. Resulting Gcode could be corrupt!") % layer.get(
                            inkex.addNS('label', 'inkscape')))
                elif i.tag == inkex.addNS('path', 'svg'):
                    if "gcodetools" not in i.keys():
                        self.paths[layer] = self.paths[layer] + \
                            [i] if layer in self.paths else [i]
                        if i.get("id") in self.selected:
                            self.selected_paths[layer] = self.selected_paths[layer] + [
                                i] if layer in self.selected_paths else [i]
                elif i.tag == inkex.addNS("g", 'svg'):
                    recursive_search(i, layer, (i.get("id") in self.selected))
                elif i.get("id") in self.selected:
                    self.error(_("This extension works with Paths and Dynamic Offsets and groups of them only! All other objects will be ignored!\nSolution 1: press Path->Object to path or Shift+Ctrl+C.\nSolution 2: Path->Dynamic offset or Ctrl+J.\nSolution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file."))

        recursive_search(self.document.getroot(), self.document.getroot())

    def get_orientation_points(self, g):
        items = g.getchildren()
        items.reverse()
        p2, p3 = [], []
        p = None
        for i in items:
            if i.tag == inkex.addNS("g", 'svg') and i.get("gcodetools") == "Gcodetools orientation point (2 points)":
                p2 += [i]
            if i.tag == inkex.addNS("g", 'svg') and i.get("gcodetools") == "Gcodetools orientation point (3 points)":
                p3 += [i]
        if len(p2) == 2:
            p = p2
        elif len(p3) == 3:
            p = p3
        if p == None:
            return None
        points = []
        for i in p:
            point = [[], []]
            for node in i:
                if node.get('gcodetools') == "Gcodetools orientation point arrow":
                    point[0] = self.apply_transforms(
                        node, cubicsuperpath.parsePath(node.get("d")))[0][0][1]
                if node.get('gcodetools') == "Gcodetools orientation point text":
                    r = re.match(
                        r'(?i)\s*\(\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*\)\s*', node.text)
                    point[1] = [float(r.group(1)), float(
                        r.group(2)), float(r.group(3))]
            if point[0] != [] and point[1] != []:
                points += [point]
        if len(points) == len(p2) == 2 or len(points) == len(p3) == 3:
            return points
        else:
            return None

################################################################################
# Fill area
################################################################################

    def area_fill(self):

        global gcode
        global offset_y
        global csp

        self.options.area_fill_angle = self.options.area_fill_angle * math.pi / 180

        print_("===================================================================")
        print_("Start filling area", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")

        if len(self.selected_paths) <= 0:
            self.error(_("This extension requires at least one selected path."))
            return
        if not self.check_dir():
            return

        for layer in self.layers:
            if layer in self.selected_paths:
                if self.options.laser_beam_with <= 0:
                    self.error(_("Laser beam with must be > 0!"))

                print_time("Time until path selection")

                for path in self.selected_paths[layer]:
                    lines = []

                    print_("")
                    print_("Working on path: ")
                    print_debug(path.get("style"), path.get("d"))

                    area_group = inkex.etree.SubElement(
                        path.getparent(), inkex.addNS('g', 'svg'))
                    d = path.get('d')

                    if d == None:
                        print_("omitting non-path")
                        self.error(_("Warning: omitting non-path"))
                        continue

                    print_time("Time for path selection")

                    csp = cubicsuperpath.parsePath(d)
                    csp = self.apply_transforms(path, csp)
                    # csp = csp_close_all_subpaths(csp)
                    csp = self.transform_csp(csp, layer)

                    print_debug("csp length: ", len(csp))
                    # print_("csp: ", csp)

                    print_time("Time for csp transformation")

                    # rotate the path to get bounds in defined direction.
                    a = - self.options.area_fill_angle
                    rotated_path = [[[[point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(
                        a)+point[1]*math.cos(a)] for point in sp] for sp in subpath] for subpath in csp]
                    bounds = csp_true_bounds(rotated_path)

                    # Draw the lines
                    # Get path's bounds
                    b = [0.0, 0.0, 0.0, 0.0]
                    for k in range(4):
                        i, j, t = bounds[k][2], bounds[k][3], bounds[k][4]
                        b[k] = csp_at_t(rotated_path[i][j-1],
                                        rotated_path[i][j], t)[k % 2]

                    print_time("Time for calculating bounds")

                    # Zig-zag
                    r = self.options.laser_beam_with
                    if r <= 0:
                        self.error(_("Laser diameter must be greater than 0!"), "error")
                        return

                    lines += [[]]

                    #i = b[0] - self.options.area_fill_shift*r
                    i = b[0] - r + 0.001
                    top = True
                    last_one = True
                    while (i < b[2] or last_one):
                        if i >= b[2]:
                            last_one = False
                        if lines[-1] == []:
                            lines[-1] += [[i, b[3]]]

                        if top:
                            lines[-1] += [[i, b[1]], [i+r, b[1]]]

                        else:
                            lines[-1] += [[i, b[3]], [i+r, b[3]]]

                        top = not top
                        i += r

                    print_time("Time for calculating zigzag pattern")

                    # print_("lines from infill:", lines)

                    # Rotate created paths back
                    a = self.options.area_fill_angle
                    lines = [[[point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(
                        a)+point[1]*math.cos(a)] for point in subpath] for subpath in lines]

                    # print_("lines: ", lines)
                    print_time("Time for rotating")

                    # get the intersection points
                    splitted_line = [[lines[0][0]]]

                    for l1, l2, in zip(lines[0], lines[0][1:]):
                        ints = []

                        if l1[0] == l2[0] and l1[1] == l2[1]:
                            continue
                        for i in range(len(csp)):
                            for j in range(1, len(csp[i])):
                                sp1, sp2 = csp[i][j-1], csp[i][j]
                                roots = csp_line_intersection(l1, l2, sp1, sp2)

                                for t in roots:
                                    p = tuple(csp_at_t(sp1, sp2, t))
                                    if l1[0] == l2[0]:
                                        t1 = (p[1]-l1[1])/(l2[1]-l1[1])
                                    else:
                                        t1 = (p[0]-l1[0])/(l2[0]-l1[0])

                                    if 0 <= t1 <= 1:
                                        ints += [[t1, p[0], p[1], i, j, t]]

                        ints.sort()

                        # mitigate vertical line glitch problem

                        if len(ints) % 2 != 0:
                            print_debug("removing intersection: ", ints)
                            ints = []

                        for i in ints:
                            splitted_line[-1] += [[i[1], i[2]]]
                            splitted_line += [[[i[1], i[2]]]]
                        splitted_line[-1] += [l2]
                        i = 0

                    print_time("Time for calculating intersections")
                    print_debug("number of splitted lines: ", len(splitted_line))

                    finalLines = []

                    # Mutitheading not working on Windows; temporarily disabled
                    if self.options.multi_thread and os.name != 'nt':
                        pool = Pool(processes=noOfThreads)
                        finalLines = pool.map(
                            checkIfLineInsideShape, splitted_line)

                        pool.close()
                        pool.join()

                    else:
                        while i < len(splitted_line):
                            finalLines += [checkIfLineInsideShape(splitted_line[i])]
                            i += 1

                    i = 0

                    print_debug("number of final lines before removing emptys: ", len(finalLines))
                    # remove empty elements
                    # print_("final_line: ", finalLines)
                    np_finalLines = np.array(finalLines, dtype=np.float32)
                    index_zeros = np.argwhere(np_finalLines == [[0, 0], [0, 0]])
                    np_finalLines = np.delete(np_finalLines, index_zeros, axis=0)

                    # print_("np_final_line: ", np_finalLines)
                    print_debug("number of final lines: ", len(np_finalLines))

                    if options.remove_tiny_infill_paths:
                        start_coords = np.array(np_finalLines[:, 0])
                        end_coords = np.array(np_finalLines[:, 1])

                        distances = np.array(end_coords[:, 1]-start_coords[:, 1])
                        distances = np.abs(distances)
                        np_finalLines = (np_finalLines[distances > (tiny_infill_factor * options.laser_beam_with)])
                        # print_("final_line: ", np_finalLines)

                    print_time("Time for calculating infill paths")

                    csp_line = csp_from_polyline(np_finalLines)
                    csp_line = self.transform_csp(csp_line, layer, True)

                    if self.get_transforms(layer) != []:
                        offset_y = self.get_transforms(layer)[1][2]
                    else:
                        offset_y = 0

                    print_time("Time for transforming infill paths")

                    curve = self.parse_curve2d(csp_line, layer)
                    self.draw_curve(curve, layer, area_group)

                    print_time("Time for drawing curve")

                    gcode += self.generate_gcode(curve, layer, self.tool_infill)

                    print_time("Time for generating Gcode")

        if gcode != '' and not(self.options.add_contours):
            self.export_gcode(gcode)

        print_("===================================================================")
        print_("Finished filling area", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")

################################################################################
# Engraving
################################################################################
    def engraving(self):
        global x1, y1
        global cspm
        global nlLT, i, j
        global max_dist  # minimum of tool radius and user's requested maximum distance

        def bisect(nx1, ny1, nx2, ny2):

            cosBis = math.sqrt(max(0, (1.0+nx1*nx2-ny1*ny2)/2.0))
            # We can get correct sign of the sin, assuming cos is positive
            if (abs(ny1-ny2) < engraving_tolerance) or (abs(cosBis) < engraving_tolerance):
                if (abs(nx1-nx2) < engraving_tolerance):
                    return(nx1, ny1, 0.0)
                sinBis = math.copysign(1, ny1)
            else:
                sinBis = cosBis*(nx2-nx1)/(ny1-ny2)
            # We can correct signs by noting that the dot product
            # of bisector and either normal must be >0
            costurn = cosBis*nx1+sinBis*ny1
            if costurn == 0:
                return (ny1*100, -nx1*100, 1)  # Path doubles back on itself
            sinturn = sinBis*nx1-cosBis*ny1
            if costurn < 0:
                sinturn = -sinturn
            if 0 < sinturn*114.6 < (0):
                sinturn = 0  # set to zero if less than the user wants to see.
            return (cosBis/costurn, sinBis/costurn, sinturn)
            # end bisect

        def bez_divide(a, b, c, d):

            bx = b[0]-a[0]
            by = b[1]-a[1]
            cx = c[0]-a[0]
            cy = c[1]-a[1]
            dx = d[0]-a[0]
            dy = d[1]-a[1]
            limit = 8*math.hypot(dx, dy) / \
                self.options.engraving_newton_iterations
            # LT This is the only limit we get from the user currently
            if abs(dx*by-bx*dy) < limit and abs(dx*cy-cx*dy) < limit:
                return [[a, b, c, d]]
            abx = (a[0]+b[0])/2.0
            aby = (a[1]+b[1])/2.0
            bcx = (b[0]+c[0])/2.0
            bcy = (b[1]+c[1])/2.0
            cdx = (c[0]+d[0])/2.0
            cdy = (c[1]+d[1])/2.0
            abcx = (abx+bcx)/2.0
            abcy = (aby+bcy)/2.0
            bcdx = (bcx+cdx)/2.0
            bcdy = (bcy+cdy)/2.0
            m = [(abcx+bcdx)/2.0, (abcy+bcdy)/2.0]
            return bez_divide(a, [abx, aby], [abcx, abcy], m) + bez_divide(m, [bcdx, bcdy], [cdx, cdy], d)
            # end of bez_divide

        def save_point(x, y, i, j):

            global cspm

            x = round(x, 3)  # round to 3 decimals
            y = round(y, 3)  # round to 3 decimals

            if len(cspm) > 1:
                _, xy1, _, i1, j1 = cspm[-1]
                if i == i1 and j == j1:  # one match
                    _, xy2, _, i1, j1 = cspm[-2]
                    if i == i1 and j == j1:  # two matches. Now test linearity
                        length1 = math.hypot(xy1[0]-x, xy1[1]-y)
                        length2 = math.hypot(xy2[0]-x, xy2[1]-y)
                        length12 = math.hypot(xy2[0]-xy1[0], xy2[1]-xy1[1])
                        # get the xy distance of point 1 from the line 0-2
                        if length2 > length1 and length2 > length12:  # point 1 between them
                            xydist = abs(
                                (xy2[0]-x)*(xy1[1]-y)-(xy1[0]-x)*(xy2[1]-y))/length2
                            if xydist < engraving_tolerance:  # so far so good
                                cspm.pop()
            cspm += [[[x, y], [x, y], [x, y], i, j]]

            # end of save_point

        # end of subfunction definitions. engraving() starts here:
        ###########################################################

        print_("===================================================================")
        print_("Start doing parameters", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")
        timestamp2 = time.time()

        global gcode
        r = 0,  # theoretical and tool-radius-limited radii in pixels
        x1, y1, nx, ny = 0, 0, 0, 0

        cspe = []
        global offset_y

        if len(self.selected_paths) <= 0:
            self.error(_("Please select at least one path to engrave and run again."))
            return
        if self.options.add_infill == False:
            if not self.check_dir():
                return

        # LT See if we can use this parameter for line and Bezier subdivision:
        bitlen = 20/self.options.engraving_newton_iterations

        for layer in self.layers:
            if layer in self.selected_paths:
                # Calculate scale in pixels per user unit (mm or inch)

                max_dist = self.options.laser_beam_with/2
                engraving_group = inkex.etree.SubElement(
                    self.selected_paths[layer][0].getparent(), inkex.addNS('g', 'svg'))

                for node in self.selected_paths[layer]:

                    print_("")
                    print_("Working on path: ")
                    print_debug(node.get("style"), node.get("d"))

                    if node.tag == inkex.addNS('path', 'svg'):
                        cspi = cubicsuperpath.parsePath(node.get('d'))

                        # LT: Create my own list. n1LT[j] is for subpath j
                        nlLT = []

                        '''
                        noOfElements = 0
                        delElements = 0
                                                
                        for j in xrange(len(cspi)):  # LT For each subpath...
                            # Remove zero length segments, assume closed path
                            i = 0
                            noOfElements += len(cspi[j])

                            while i < len(cspi[j]):
                                if abs(cspi[j][i-1][1][0]-cspi[j][i][1][0]) < engraving_tolerance and abs(cspi[j][i-1][1][1]-cspi[j][i][1][1]) < engraving_tolerance:
                                    cspi[j][i-1][2] = cspi[j][i][2]
                                    del cspi[j][i]
                                    delElements += 1
                                else:
                                    i += 1

                        print_("Total number of points ", noOfElements)
                        print_("Number of removed points: ", delElements)

                        '''

                        for csp in cspi:
                            # print_("csp is",csp)
                            nlLT.append([])
                            for i in range(0, len(csp)):  # LT for each point
                                sp0, sp1, sp2 = csp[i-2], csp[i-1], csp[i]
                                # LT find angle between this and previous segment
                                x0, y0 = sp1[1]
                                nx1, ny1 = csp_normalized_normal(sp1, sp2, 0)
                                nx0, ny0 = csp_normalized_normal(sp0, sp1, 1)
                                bx, by, s = bisect(nx0, ny0, nx1, ny1)
                                # record x,y,normal,ifCorner, sin(angle-turned/2)
                                nlLT[-1] += [[[x0, y0], [bx, by], True, s]]

                                # LT now do the line
                                if sp1[1] == sp1[2] and sp2[0] == sp2[1]:  # straightline
                                    nlLT[-1] += [[sp1[1], [nx1, ny1], False, i]]
                                else:  # Bezier. First, recursively cut it up:
                                    nn = bez_divide(
                                        sp1[1], sp1[2], sp2[0], sp2[1])
                                    first = True  # Flag entry to divided Bezier
                                    for bLT in nn:  # save as two line segments
                                        for seg in range(3):
                                            if seg > 0 or first:
                                                nx1 = bLT[seg][1]-bLT[seg+1][1]
                                                ny1 = bLT[seg+1][0]-bLT[seg][0]
                                                l1 = math.hypot(nx1, ny1)
                                                if l1 < engraving_tolerance:
                                                    continue
                                                nx1 = nx1/l1  # normalise them
                                                ny1 = ny1/l1
                                                nlLT[-1] += [[bLT[seg],
                                                              [nx1, ny1], False, i]]
                                                first = False
                                            if seg < 2:  # get outgoing bisector
                                                nx0 = nx1
                                                ny0 = ny1
                                                nx1 = bLT[seg+1][1] - \
                                                    bLT[seg+2][1]
                                                ny1 = bLT[seg+2][0] - \
                                                    bLT[seg+1][0]
                                                l1 = math.hypot(nx1, ny1)
                                                if l1 < engraving_tolerance:
                                                    continue
                                                nx1 = nx1/l1  # normalise them
                                                ny1 = ny1/l1
                                                # bisect
                                                bx, by, s = bisect(
                                                    nx0, ny0, nx1, ny1)
                                                nlLT[-1] += [[bLT[seg+1],
                                                              [bx, by], True, 0.]]

                        reflex = False
                        for j in xrange(len(nlLT)):  # LT6b for each subpath
                            cspm = []  # Will be my output. List of csps.
                            r = 0  # LT initial, as first point is an angle
                            for i in xrange(len(nlLT[j])):  # LT for each node
                                n0 = nlLT[j][i-2]  # previous node
                                n1 = nlLT[j][i-1]  # current node
                                n2 = nlLT[j][i]  # next node
                                # this point/start of this line
                                x1a, y1a = n1[0]
                                nx, ny = n1[1]
                                x1b, y1b = n2[0]  # next point/end of this line
                                if n1[2]:  # We're at a corner
                                    bits = 1
                                    bit0 = 0

                                else:  # line. Cut it up if long.
                                    if n0[3] > 0 and not self.options.engraving_draw_calculation_paths:
                                        bit0 = r*n0[3]  # after acute corner
                                    else:
                                        bit0 = 0.0
                                    length = math.hypot((x1b-x1a), (y1a-y1b))
                                    bit0 = (min(length, bit0))
                                    bits = int((length-bit0)/bitlen)
                                    # split excess evenly at both ends
                                    bit0 += (length-bit0-bitlen*bits)/2
                                    # print_("j,i,r,bit0,bits",j,i,w,bit0,bits)
                                for b in xrange(bits):  # divide line into bits
                                    x1 = x1a+ny*(b*bitlen+bit0)
                                    y1 = y1a-nx*(b*bitlen+bit0)

                                    if reflex:  # just after a reflex corner
                                        reflex = False
                                    if n1[2]:  # We're at a corner
                                        if n1[3] > 0:  # acute
                                            save_point(x1, y1, i, j)
                                            save_point(x1, y1, i, j)

                                    save_point(x1, y1, i, j)

                                # LT end of for each bit of this line
                                if n1[2] == True and n1[3] < 0:  # reflex angle
                                    reflex = True
                            # LT next i
                            if len(cspm) != 0:
                                cspm += [cspm[0]]

                                # for entr in cspm:
                                #    print_("cspm ", entr)
                                cspe += [cspm]

                if cspe != []:

                    # for entr in cspe:
                        #    print_("cspe ", entr)
                    curve = self.parse_curve(cspe, layer)
                    # for entr in curve:
                    #    print_("curve ", entr)
                    self.draw_curve(curve, layer, engraving_group)

                    if self.get_transforms(layer) != []:
                        offset_y = self.get_transforms(layer)[1][2]
                    else:
                        offset_y = 0

                    gcode += self.generate_gcode(curve, layer, self.tool_perimeter)

        if gcode != '':
            self.export_gcode(gcode)

        print_("===================================================================")
        print_("Finished doing parameters", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")
        print_(time.time() - timestamp2, "s for parameters")

################################################################################
# Orientation
################################################################################
    def orientation(self, layer=None):
        print_("entering orientations")
        if layer == None:
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer in self.orientation_points:
            self.error(_("Active layer already has orientation points! Remove them or select another layer!"))

        orientation_group = inkex.etree.SubElement(layer, inkex.addNS(
            'g', 'svg'), {"gcodetools": "Gcodetools orientation group"})

        if layer.get("transform") != None:
            translate = layer.get("transform").replace(
                "translate(", "").replace(")", "").split(",")
        else:
            translate = [0, 0]

        global doc_height
        doc_height = self.unittouu(self.document.getroot().xpath(
            '@height', namespaces=inkex.NSS)[0])

        if self.document.getroot().get('height') == "100%":
            doc_height = 1052.3622047
            print_("Overruding height from 100 percents to %s" % doc_height)

        print_("Document height: " + str(doc_height))

        # setup for mm as base unit
        points = [[0., 0., 0.], [100., 0., 0.], [0., 100., 0.]]
        orientation_scale = 1
        print_("orientation_scale < 0 ===> switching to mm units=%0.10f" %
               orientation_scale)

        points = points[:2]

        print_(("using orientation scale", orientation_scale, "i=", points))
        print_("translate", translate[1])
        for i in points:
            si = [i[0]*orientation_scale,
                  (i[1]*orientation_scale)+float(translate[1])]
            g = inkex.etree.SubElement(orientation_group, inkex.addNS('g', 'svg'), {
                'gcodetools': "Gcodetools orientation point (2 points)"})
            inkex.etree.SubElement(g, inkex.addNS('path', 'svg'),
                                   {
                'style':    "stroke:none;fill:#000000;",
                'd': 'm %s,%s 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000001 z z' % (si[0], -si[1]+doc_height),
                'gcodetools': "Gcodetools orientation point arrow"
            })
            t = inkex.etree.SubElement(g, inkex.addNS('text', 'svg'),
                                       {
                'style':    "font-size:10px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;fill:#000000;fill-opacity:1;stroke:none;",
                inkex.addNS("space", "xml"): "preserve",
                'x':    str(si[0]),
                'y':    str(-si[1]+doc_height),
                'gcodetools': "Gcodetools orientation point text"
            })
            t.text = "(%s; %s; %s)" % (i[0], i[1], i[2])


################################################################################
# ApplyTransform
################################################################################

    @staticmethod
    def objectToPath(node):

        # print_("convertring node to path:", node.tag)

        if node.tag == inkex.addNS('g', 'svg'):
            return node

        if node.tag == inkex.addNS('path', 'svg') or node.tag == 'path':
            for attName in node.attrib.keys():
                if ("sodipodi" in attName) or ("inkscape" in attName):
                    del node.attrib[attName]
            return node

        return node

    def recursiveFuseTransform(self, node, transf=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]):

        # print_("transforming Node:", node.tag)

        transf = composeTransform(
            transf, parseTransform(node.get("transform", None)))

        if 'transform' in node.attrib:
            del node.attrib['transform']

        if 'style' in node.attrib:
            style = node.attrib.get('style')
            style = simplestyle.parseStyle(style)
            update = False

            if 'stroke-width' in style:
                try:
                    stroke_width = self.unittouu(
                        style.get('stroke-width').strip())
                    # pixelsnap ext assumes scaling is similar in x and y
                    # and uses the x scale...
                    # let's try to be a bit smarter
                    stroke_width *= math.sqrt(transf[0]
                                              [0]**2 + transf[1][1]**2)
                    style['stroke-width'] = str(stroke_width)
                    update = True
                except AttributeError:
                    pass

            if update:
                style = simplestyle.formatStyle(style)
                node.attrib['style'] = style

        node = self.objectToPath(node)

        if 'd' in node.attrib:
            d = node.get('d')
            p = cubicsuperpath.parsePath(d)
            applyTransformToPath(transf, p)
            node.set('d', cubicsuperpath.formatPath(p))

        elif node.tag in [inkex.addNS('polygon', 'svg'),
                          inkex.addNS('polyline', 'svg')]:
            points = node.get('points')
            points = points.strip().split(' ')
            for k, p in enumerate(points):
                if ',' in p:
                    p = p.split(',')
                    p = [float(p[0]), float(p[1])]
                    applyTransformToPoint(transf, p)
                    p = [str(p[0]), str(p[1])]
                    p = ','.join(p)
                    points[k] = p
            points = ' '.join(points)
            node.set('points', points)

        elif node.tag in [inkex.addNS('rect', 'svg'),
                          inkex.addNS('text', 'svg'),
                          inkex.addNS('image', 'svg'),
                          inkex.addNS('use', 'svg')]:
            node.set('transform', formatTransform(transf))

        for child in node.getchildren():
            self.recursiveFuseTransform(child, transf)

    def applytransforms(self):
         # Apply transformations
        self.getselected()

        if self.selected:
            for _, shape in self.selected.items():
                self.recursiveFuseTransform(shape, parseTransform(None))
        else:
            self.recursiveFuseTransform(
                self.document.getroot(), parseTransform(None))
        # Transformations applied

################################################################################
# Effect
# Main function of Lasertools class
################################################################################
    def effect(self):
        global options
        options = self.options
        options.self = self
        options.doc_root = self.document.getroot()
        global print_

        if self.options.log_create_log:
            if os.path.isdir(self.options.directory):
                try:
                    if os.path.isfile(self.options.directory+"/log.txt"):
                        os.remove(self.options.directory+"/log.txt")
                    f = open(self.options.directory+"/log.txt", "a")
                    f.write(
                        "===================================================================\n")
                    f.write("Lasertools log file.\nStarted at %s.\n%s\n" % (
                        time.strftime("%d.%m.%Y %H:%M:%S"), options.directory+"/log.txt"))
                    f.write(
                        "===================================================================\n\n")
                    f.close()
                except:
                    print_ = lambda *x: None

            else:
                self.error(("Directory does not exist! Please specify existing directory at options tab!"), "error")

        else:
            print_ = lambda *x: None
        self.get_info()

        # wait to attatch debugger to process
        if debug:
            print_debug("Python version:", sys.version_info)
            print_debug("Waiting for Debugger to be attached")
            time.sleep(3)
            print_debug("Starting Program")

        if self.orientation_points == {}:
            self.orientation(self.layers[min(0, len(self.layers)-1)])
            self.get_info()

        self.tool_infill = {
            "name": "Laser Engraver Infill",
            "id": "Laser Engraver Infill",
            "penetration feed": self.options.laser_speed,
            "feed": self.options.laser_speed,
            "gcode before path": ("G04 P" + self.options.power_delay + " " + self.options.laser_command),
            "gcode after path": self.options.laser_off_command
        }

        self.tool_perimeter = {
            "name": "Laser Engraver Perimeter",
            "id": "Laser Engraver Perimeter",
            "penetration feed": self.options.laser_param_speed,
            "feed": self.options.laser_param_speed,
            "gcode before path": ("G04 P" + self.options.power_delay + " " + self.options.laser_command_perimeter),
            "gcode after path": self.options.laser_off_command
        }

        self.get_info()

        print_("Applying all transformations")
        self.applytransforms()

        if self.options.add_infill:
            self.selected_paths = self.paths
            self.area_fill()

        if self.options.add_contours:
            self.get_info()
            self.selected_paths = self.paths

            if profiling:
                if os.path.isfile(self.options.directory+"performance.prof"):
                    os.remove(self.options.directory+"/performance.prof")

                profile = cProfile.Profile()
                profile.runctx('self.engraving()', globals(),
                               locals())

                kProfile = lsprofcalltree.KCacheGrind(profile)

                kFile = open(self.options.directory+"/performance.prof", 'w+')
                kProfile.output(kFile)
                kFile.close()

            self.engraving()


e = laser_gcode()
e.affect()
