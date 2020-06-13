#!/usr/bin/env python
"""
Modified by Christoph Wagner 2020
Modified by Alain Pelletier 2020

based on gcodetools, https://gitlab.com/inkscape/extensions/-/blob/master/gcodetools.py
based on flatten, https://gitlab.com/inkscape/extensions/-/blob/master/flatten.py

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public Licensey
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""
# standard libraries
import cmath
import copy
import math
import os
import re
import sys
import time
import multiprocessing
from multiprocessing import Pool
import datetime

# 3rd party libraries
import numpy as np


# local libraries
import inkex
from inkex import CubicSuperPath, Style, bezier
from inkex.bezier import bezierparameterize
from inkex.transforms import Transform
from inkex.elements import PathElement, Group
from inkex.paths import Path

if sys.version_info[0] > 2:
    xrange = range
    unicode = str

if "errormsg" not in dir(inkex):
    inkex.errormsg = lambda msg: sys.stderr.write((str(msg) + "\n").encode("UTF-8"))


################################################################################
# Styles and additional parameters
################################################################################

PROFILING = False   # Disable if not debugging
DEBUG = False     # Disable if not debugging
TINY_INFILL_FACTOR = 2  # x times the laser beam width will be removed
ENGRAVING_TOLERANCE = 0.000002
DEFAULTS = {
    'header': """
;Inkscape Lasertools G-code

""",
    'footer': """G00 X0 Y0

"""
}

if PROFILING:
    import lsprofcalltree
    import cProfile


gcode = ""
csp = []
timestamp = datetime.datetime.now()
options = {}
cspm = []
orient_points = []


def marker_style(stroke, marker='DrawCurveMarker', width=1):
    """Set a marker style with some basic defaults"""
    return Style(stroke=stroke, fill='none', stroke_width=width,
                 marker_end='url(#{})'.format(marker))


MARKER_STYLE = {
    "biarc_style": {
        'line': marker_style('#f88', width=0.1),
        'biarc1': marker_style('#8f8', width=0.5)
    }
}


def check_if_line_inside_shape(splitted_line_csp):

    splitted_line = splitted_line_csp[0]
    csp_temp = splitted_line_csp[1]
    l1, l2 = splitted_line[0], splitted_line[1]

    p = [(l1[0]+l2[0])/2, (l1[1]+l2[1])/2]

    if point_inside_csp(p, csp_temp):
        return [l1, l2]

    else:
        return [[0, 0], [0, 0]]


def csp_bounds(csp_temp):

    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx, miny, maxx, maxy = float("inf"), float("inf"), 0, 0
    for i in range(len(csp_temp)):
        for j in range(1, len(csp_temp[i])):

            minx_temp, miny_temp, maxx_temp, maxy_temp = bounds_of_line([csp_temp[i][j-1][0][0], csp_temp[i][j-1][0][1]], [csp_temp[i][j][0][0], csp_temp[i][j][0][1]])

            maxx = max(maxx_temp, maxx)
            minx = min(minx_temp, minx)
            maxy = max(maxy_temp, maxy)
            miny = min(miny_temp, miny)

    return minx, miny, maxx, maxy


def bounds_of_line(p1, p2):

    minx = min((p1[0], p2[0]))
    maxx = max((p1[0], p2[0]))
    miny = min((p1[1], p2[1]))
    maxy = max((p1[1], p2[1]))

    return minx, miny, maxx, maxy


def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
        return []

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return [x, y]


def csp_line_intersection(l1, l2, sp1, sp2):

    precission = 10
    sp1l = sp1[0]
    sp2l = sp2[0]

    inters = line_intersection((l1, l2), (sp1l, sp2l))

    spx_min, spy_min, spx_max, spy_max = bounds_of_line(sp1l, sp2l)

    spx_min = min((sp1l[0], sp2l[0]))
    spx_max = max((sp1l[0], sp2l[0]))
    spy_min = min((sp1l[1], sp2l[1]))
    spy_max = max((sp1l[1], sp2l[1]))

    if sp1l == sp2l or l1 == l2 or inters == []:
        return []

    if round(inters[0], precission) < round(spx_min, precission) or round(inters[0], precission) > round(spx_max, precission):
        return []

    if round(inters[1], precission) < round(spy_min, precission) or round(inters[1], precission) > round(spy_max, precission):
        return []

    return inters


def csp_normalized_normal(sp1, sp2):

    dif_x = sp1[0][0] - sp2[0][0]
    dif_y = sp1[0][1] - sp2[0][1]
    len_seg = math.sqrt(dif_x**2+dif_y**2)

    if len_seg != 0:
        nx = dif_y/len_seg
        ny = dif_x/len_seg
    else:
        return[0, 1]

    return [-ny, nx]


################################################################################
# Area Fill Stuff
################################################################################

def point_inside_csp(p, csp, on_the_path=True):

    x, y = p
    ray_intersections_count = 0

    for subpath in csp:
        for i in range(1, len(subpath)):
            sp1, sp2 = subpath[i-1], subpath[i]

            p = csp_line_intersection([x, y], [x, y+5], sp1, sp2)

            if len(p) == 2:

                if p[1] == y:
                    return on_the_path

                else:
                    if p[1] > y:
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
# Common functions
################################################################################

################################################################################
# Logging
################################################################################


def print_(*arg):
    if os.path.isdir(options.directory):
        f = open(options.directory+"/log.txt", "a")
        for s in arg:
            s = str(str(s))+" "
            f.write(s)
        f.write("\n")
        f.close()


def print_debug(*arg):
    if DEBUG:
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

    def __repr__(self): return '%f,%f' % (self.x, self.y)

    def l2(self): return self.x*self.x + self.y*self.y


################################################################################
# Lasertools class
################################################################################

class laser_gcode(inkex.EffectExtension):

    def export_gcode(self, gcode):
        gcode_pass = gcode
        for _ in range(1, self.options.passes):
            if self.options.z_stepdown == 0:
                gcode += "\nG90 \n" + gcode_pass
            else:
                gcode += "\nG91 \nG0 Z%s \nG90 \n" % self.options.z_stepdown + gcode_pass

        f = open(self.options.directory+self.options.file, "w")

        if self.options.prefix_1 != "":
            self.header += self.options.prefix_1 + "\n"
        if self.options.prefix_2 != "":
            self.header += self.options.prefix_2 + "\n"
        if self.options.prefix_3 != "":
            self.header += self.options.prefix_3 + "\n"

        if self.options.suffix_1 != "":
            self.footer += self.options.suffix_1 + "\n"
        if self.options.suffix_2 != "":
            self.footer += self.options.suffix_2 + "\n"
        if self.options.suffix_3 != "":
            self.footer += self.options.suffix_3 + "\n"

        f.write(self.header + "\n" + gcode + self.footer)
        f.close()

    def add_arguments(self, pars):
        add_argument = pars.add_argument
        add_argument("-d", "--directory", dest="directory", default="/insert your target directory here", help="Output directory")
        add_argument("-f", "--filename", dest="file", default="output.ngc", help="File name")
        add_argument("--add-numeric-suffix-to-filename", type=inkex.Boolean, dest="add_numeric_suffix_to_filename", default=False, help="Add numeric suffix to file name")
        add_argument("--laser-command-perimeter", dest="laser_command_perimeter", default="S100", help="Laser gcode command Perimeter")
        add_argument("--laser-command", dest="laser_command", default="S100", help="Laser gcode command infill")
        add_argument("--laser-off-command", dest="laser_off_command", default="S1", help="Laser gcode end command")
        add_argument("--laser-beam-with", type=float, dest="laser_beam_with", default="0.3", help="Laser speed (mm/min)")
        add_argument("--infill-overshoot", type=float, dest="infill_overshoot", default="0.0", help="Overshoot to limit acceleration overburn")
        add_argument("--contour-tolerance", type=float, dest="tolerance", default="0.1", help="Tolerance for contour approximation")
        add_argument("--laser-speed", type=int, dest="laser_speed", default="1200", help="Laser speed (mm/min)")
        add_argument("--laser-param-speed", type=int, dest="laser_param_speed", default="700", help="Laser speed for Parameter (mm/min)")
        add_argument("--passes", type=int, dest="passes", default="1", help="Quantity of passes")
        add_argument("--power-delay", dest="power_delay", default="0", help="Laser power-on delay (ms)")
        add_argument("--z-stepdown", type=float, dest="z_stepdown", default="0.0", help="Z-stepdown per pass for cutting operations")
        add_argument("--linuxcnc", type=inkex.Boolean, dest="linuxcnc", default=False, help="Use G64 P0.1 trajectory planning")
        add_argument("--add-contours", type=inkex.Boolean, dest="add_contours", default=True, help="Add contour to Gcode paths")
        add_argument("--add-infill", type=inkex.Boolean, dest="add_infill", default=True, help="Add infill to Gcode paths")
        add_argument("--remove-tiny-infill-paths", type=inkex.Boolean, dest="remove_tiny_infill_paths", default=False, help="Remove tiny infill paths from Gcode")
        add_argument("--prefix1", dest="prefix_1", default="", help="First line before G-Code starts")
        add_argument("--prefix2", dest="prefix_2", default="", help="Second line before G-Code starts")
        add_argument("--prefix3", dest="prefix_3", default="", help="Third line before G-Code starts")
        add_argument("--suffix1", dest="suffix_1", default="", help="First line after G-Code ends")
        add_argument("--suffix2", dest="suffix_2", default="", help="Second line after G-Code ends")
        add_argument("--suffix3", dest="suffix_3", default="", help="Third line after G-Code ends")

        add_argument("--multi_thread", type=inkex.Boolean, dest="multi_thread", default=True, help="Activate multithreading support")

        add_argument("--suppress-all-messages", type=inkex.Boolean, dest="suppress_all_messages", default=True, help="Hide messages during g-code generation")
        add_argument("--create-log", type=inkex.Boolean, dest="log_create_log", default=True, help="Create log files")
        add_argument("--area-fill-angle", type=float, dest="area_fill_angle", default="0", help="Fill area with lines heading this angle")

        add_argument("--active-tab", dest="active_tab", default="",	help="Defines which tab is active")

    def __init__(self):
        super(laser_gcode, self).__init__()

    def parse_curve(self, p, layer, w=None, f=None):
        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)

        # Sort to reduce Rapid distance
        k = list(range(1, len(p)))
        keys = [0]
        while len(k) > 0:
            end = p[keys[-1]][-1][1]
            dist = (-10000000, -10000000)

            for i in range(len(k)):
                start = p[k[i]][0][1]
                dist = max((-((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2), i), dist)

            keys += [k[dist[1]]]
            del k[dist[1]]
        for k in keys:
            subpath = p[k]
            c += [[[subpath[0][1][0], subpath[0][1][1]], 'move', 0, 0]]
            for i in range(1, len(subpath)):
                sp1 = [[subpath[i - 1][j][0], subpath[i - 1][j][1]] for j in range(3)]
                sp2 = [[subpath[i][j][0], subpath[i][j][1]] for j in range(3)]
                c += [[[sp1[0][0], sp1[0][1]], 'line', [sp2[0][0], sp2[0][1]]]]
            c += [[[subpath[-1][1][0], subpath[-1][1][1]], 'end', 0, 0]]
        return c

    def parse_curve2d(self, p, layer):

        c = []
        if len(p) == 0:
            return []
        p = self.transform_csp(p, layer)
        np_p = np.array(p)
        sorted_tool_paths = self.sort_points(np_p[:, 0, 0, 0], np_p[:, 0, 0, 1], np_p[:, 1, 0, 0], np_p[:, 1, 0, 1])

        for path in sorted_tool_paths:
            c += [[[path[0], path[1]], 'move']]
            c += [[[path[0], path[1]], 'line', [path[2], path[3]]]]
            c += [[[path[2], path[3]], 'end']]

        return c

    def sort_points(self, x1, y1, x2, y2):
        sorted_list = np.zeros((len(x1), 4))

        xpos = np.array(x1[1:])
        xpos_inv = np.array(x2[1:])
        ypos = np.array(y1[1:])
        ypos_inv = np.array(y2[1:])

        sorted_list[0] = [x1[0], y1[0], x2[0], y2[0]]
        actXPos = x2[0]
        actYPos = y2[0]

        i = 1

        while len(xpos) > 0:
            xDist = np.abs(xpos - actXPos)
            xDistInv = np.abs(xpos_inv - actXPos)
            yDist = np.abs(ypos - actYPos)
            yDistInv = np.abs(ypos_inv - actYPos)

            distances = np.array(
                [np.add(xDist, yDist), np.add(xDistInv, yDistInv)])
            distances = np.abs(distances)

            minInd = np.unravel_index(
                np.argmin(distances, axis=None), distances.shape)

            if minInd[0] == 0:
                sorted_list[i] = [xpos[minInd[1]], ypos[minInd[1]],
                                  xpos_inv[minInd[1]], ypos_inv[minInd[1]]]
                actXPos = xpos_inv[minInd[1]]
                actYPos = ypos_inv[minInd[1]]

            else:
                sorted_list[i] = [xpos_inv[minInd[1]],
                                  ypos_inv[minInd[1]], xpos[minInd[1]], ypos[minInd[1]]]
                actXPos = xpos[minInd[1]]
                actYPos = ypos[minInd[1]]

            xpos = np.delete(xpos, minInd[1])
            ypos = np.delete(ypos, minInd[1])
            xpos_inv = np.delete(xpos_inv, minInd[1])
            ypos_inv = np.delete(ypos_inv, minInd[1])

            i = i+1

        return sorted_list

    def draw_curve(self, curve, layer, group=None, style=MARKER_STYLE["biarc_style"]):

        for i in [0, 1]:
            sid = 'biarc1_r'.format(i)
            style[sid] = style['biarc1'.format(i)].copy()

        if group is None:
            group = self.layers[min(1, len(self.layers) - 1)].add(Group(gcodetools="Preview group"))
            if not hasattr(self, "preview_groups"):
                self.preview_groups = {layer: group}
            elif layer not in self.preview_groups:
                self.preview_groups[layer] = group
            group = self.preview_groups[layer]

        s = ''

        transform = self.get_transforms(group)
        if transform:
            transform = self.reverse_transform(transform)
            transform = str(Transform(transform))

        a = [0., 0.]
        b = [1., 0.]
        c = [0., 1.]
        k = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        a = self.transform(a, layer, True)
        b = self.transform(b, layer, True)
        c = self.transform(c, layer, True)

        for sk in curve:
            si = sk[:]
            si[0] = self.transform(si[0], layer, True)
            if len(si) == 3:
                si[2] = self.transform(si[2], layer, True)

            if s != '':
                if s[1] == 'line':
                    elem = group.add(PathElement(gcodetools="Preview"))
                    elem.transform = transform
                    elem.style = style['line']
                    elem.path = 'M {},{} L {},{}'.format(s[0][0], s[0][1], si[0][0], si[0][1])

            s = si

    def check_dir(self):
        if self.options.directory[-1] not in ["/", "\\"]:
            if "\\" in self.options.directory:
                self.options.directory += "\\"
            else:
                self.options.directory += "/"
        print_("Checking direcrory: ", self.options.directory)
        if os.path.isdir(self.options.directory):
            if os.path.isfile(self.options.directory+'header'):
                f = open(self.options.directory+'header', 'r')
                self.header = f.read()
                f.close()
            else:
                self.header = DEFAULTS['header']
            if os.path.isfile(self.options.directory+'footer'):
                f = open(self.options.directory+'footer', 'r')
                self.footer = f.read()
                f.close()
            else:
                self.footer = DEFAULTS['footer']

        else:
            self.error("Directory does not exist! Please specify existing directory at options tab!", "error")
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

        print_("Testing writing rights on ", self.options.directory+self.options.file)
        try:
            f = open(self.options.directory+self.options.file, "w")
            f.close()
        except:
            self.error("Can not write to specified file!\n{}".format(self.options.directory+self.options.file), "error")
            return False
        return True

################################################################################
# Generate Gcode
# Generates Gcode on given curve.
# Crve defenitnion [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
# strategy is either infill or parameter each will calculate differently.
# infill strategy will allow for acceleration and deceleration buffer distance to prevent speed change burn.
#
################################################################################

    def generate_gcode(self, curve, layer, tool, strategy):

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

        lg, f = 'G00', "F%.1f" % tool['penetration feed']
        g = "; START " + strategy+" strategy\nG01 " + f + "\n"

        # linuxcnc trajectory planning to limit corner burns and acceleration burns
        if strategy == 'infill' and self.options.linuxcnc:
            g += "G64 P0 ;linuxcnc continuous mode trajectory planning\n"
        if strategy == 'perimeter' and self.options.linuxcnc:
            g += "G64 P0.15 ;linuxcnc blend tolerence mode trajectory planning\n"

        # set the begining past coordinates to unlikely numbers
        pastX, pastY = -10000.05, 10000.01

        for i in range(1, len(curve)):
            #    Creating Gcode for curve between s=curve[i-1] and si=curve[i] start at s[0] end at s[4]=si[0]
            s, si = curve[i-1], curve[i]
            newcoord_different = si[0][0] != pastX or si[0][1] != pastY

            #############################
            # infill strategy
            #############################
            # checks for moves and writes the G00 X0.00, Y0.00
            # move with overshoot just moves the X
            # move without overshoot moves the X and Y
            if newcoord_different and s[1] == 'move' and strategy == "infill":
                if round(self.options.infill_overshoot, 1) > 0:
                    g += "G00 X" + str(round(si[0][0], 2)) + "\n"
                else:
                    g += "G00" + c(si[0]) + "\n"
                # write past used command and coordinates
                pastX, pastY, lg = round(
                    si[0][0], 2), round(si[0][1], 2), 'G00'

                # Check if the line is going up or down.
                # sets the laser head to start moving before the laser fires
                # fires the laser arrives at destination
                # turns off the laser and overshoots the end
                # The overshoots gives a buffer for
                # accelerating and decelerating the head
                # if overshoot is selected to be 0.0, the
                # overshoot Gcode is ignored
                # if overshoot is >0.0 G00 Y to overshoot location
            elif newcoord_different and s[1] == 'line' and strategy == "infill" and lg == 'G00':
                # detect up direction
                if round(si[0][1], 2) > pastY:
                    if round(self.options.infill_overshoot, 1) > 0:
                        g += "G00 Y" + \
                            str(round(pastY-self.options.infill_overshoot, 2)) + "\n"
                        g += "G01 Y" + str(pastY) + "\n"
                    g += tool['gcode before path'] + "\n"
                    g += "G01 Y" + str(round(si[0][1], 2)) + "\n"
                    g += tool['gcode after path'] + "\n"
                    if round(self.options.infill_overshoot, 1) > 0:
                        g += "G01 Y" + \
                            str(round(
                                si[0][1]+self.options.infill_overshoot, 2)) + "\n"
                        # write past used command and coordinates
                    pastX, pastY, lg = round(
                        si[0][0], 2), round(si[0][1], 2), 'G01'
                    # detect down direction
                elif round(si[0][1], 2) < pastY:
                    if round(self.options.infill_overshoot, 1) > 0:
                        g += "G00 Y" + \
                            str(round(pastY+self.options.infill_overshoot, 2)) + "\n"
                        g += "G01 Y" + str(pastY) + "\n"
                    g += tool['gcode before path'] + "\n"
                    g += "G01 Y" + str(round(si[0][1], 2)) + "\n"
                    g += tool['gcode after path'] + "\n"
                    if round(self.options.infill_overshoot, 1) > 0:
                        g += "G01 Y" + \
                            str(round(
                                si[0][1]-self.options.infill_overshoot, 2)) + "\n"
                        # write past used command and coordinates
                    pastX, pastY, lg = round(
                        si[0][0], 2), round(si[0][1], 2), 'G01'

                    #############################
                    # perimeter strategy
                    #############################
                    # turns off laser before issuing a G00 move instruction
            elif newcoord_different and s[1] == 'move' and strategy == "perimeter":
                # turn off laser before fast move
                g += tool['gcode after path'] + "\n"
                g += "G00" + c(si[0]) + "\n"
                # write past used command and coordinates
                pastX, pastY, lg = round(
                    si[0][0], 2), round(si[0][1], 2), 'G00'
            elif newcoord_different and s[1] == 'line' and strategy == "perimeter":
                if lg == 'G00':
                    # burn laser only after a G00 move
                    g += tool['gcode before path'] + "\n"
                x, y = round(si[0][0], 2), round(si[0][1], 2)
                gx, gy = "", ""  # clear gx and gy
                if x != pastX:
                    # only include X0.00 coordinates if they are diffrent from past burn
                    gx = " X"+str(x)
                if y != pastY:
                    # only include Y0.00 coordinates if they are diffrent from past burn
                    gy = " Y"+str(y)
                g += "G01" + gx+gy + "\n"
                # write past used command and coordinates
                pastX, pastY, lg = round(si[0][0], 2), round(si[0][1], 2), 'G01'

                # Turn off laser before leaving
        g += tool['gcode after path'] + "\n;END " + strategy + "\n\n"
        return g

    def get_transforms(self, g):
        root = self.document.getroot()
        trans = []
        while g != root:
            if 'transform' in g.keys():
                t = g.get('transform')
                t = Transform(t).matrix
                trans = (Transform(t) * Transform(trans)).matrix if trans != [] else t

            g = g.getparent()
        return trans

    def reverse_transform(self, transform):
        trans = np.array(transform + ([0, 0, 1],))
        if np.linalg.det(trans) != 0:
            trans = np.linalg.inv(trans).tolist()[:2]
            return trans
        else:
            return transform

    def transform(self, source_point, layer, reverse=False):
        if layer == None:
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer not in self.transform_matrix:
            for i in range(self.layers.index(layer), -1, -1):
                if self.layers[i] in self.orientation_points:
                    break

            if self.layers[i] not in self.orientation_points:
                self.error("Orientation points for '{}' layer have not been found! Please add orientation points using Orientation tab!".format(
                    layer.get(inkex.addNS('label', 'inkscape'))), "error")
            elif self.layers[i] in self.transform_matrix:
                self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
            else:
                orientation_layer = self.layers[i]
                if len(self.orientation_points[orientation_layer]) > 1:
                    self.error("There are more than one orientation point groups in '{}' layer".format(orientation_layer.get(inkex.addNS('label', 'inkscape'))))
                points = self.orientation_points[orientation_layer][0]
                if len(points) == 2:
                    points += [[[(points[1][0][1]-points[0][0][1])+points[0][0][0], -(points[1][0][0]-points[0][0][0])+points[0][0][1]],
                                [-(points[1][1][1]-points[0][1][1])+points[0][1][0], points[1][1][0]-points[0][1][0]+points[0][1][1]]]]
                if len(points) == 3:

                    matrix = np.array([
                        [points[0][0][0], points[0][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[0][0][0], points[0][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[0][0][0], points[0][0][1], 1],
                        [points[1][0][0], points[1][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[1][0][0], points[1][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[1][0][0], points[1][0][1], 1],
                        [points[2][0][0], points[2][0][1], 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, points[2][0][0], points[2][0][1], 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, points[2][0][0], points[2][0][1], 1]
                    ])

                    if np.linalg.det(matrix) != 0:
                        m = np.linalg.solve(matrix,
                                            np.array(
                                                [[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]], [
                                                    points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
                                            )
                                            ).tolist()
                        self.transform_matrix[layer] = [[m[j*3+i][0] for i in range(3)] for j in range(3)]

                    else:
                        self.error("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)")
                else:
                    self.error("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)")

            self.transform_matrix_reverse[layer] = np.linalg.inv(self.transform_matrix[layer]).tolist()

        x, y = source_point[0], source_point[1]

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
                for k in xrange(len(csp[i][j])):
                    csp[i][j][k] = self.transform(csp[i][j][k], layer, reverse)
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
# Get Gcodetools info from the svg
################################################################################

    def get_info(self):
        self.svg.selected_paths = {}
        self.paths = {}
        self.orientation_points = {}
        self.layers = [self.document.getroot()]
        self.transform_matrix = {}
        self.transform_matrix_reverse = {}

        def recursive_search(g, layer, selected=False):
            items = g.getchildren()
            items.reverse()
            for i in items:
                if selected:
                    self.svg.selected[i.get("id")] = i
                if i.tag == inkex.addNS("g", 'svg') and i.get(inkex.addNS('groupmode', 'inkscape')) == 'layer':
                    self.layers += [i]
                    recursive_search(i, i)
                elif i.get('gcodetools') == "Gcodetools orientation group":
                    points = orient_points
                    if points != None:
                        self.orientation_points[layer] = self.orientation_points[layer]+[points[:]] if layer in self.orientation_points else [points[:]]
                        print_("    Found orientation points in '{}' layer: '{}'".format(layer.get(inkex.addNS('label', 'inkscape')), points))
                    else:
                        self.error("Warning! Found bad orientation points in '{}' layer. Resulting Gcode could be corrupt!".format(layer.get(inkex.addNS('label', 'inkscape'))))
                elif i.tag == inkex.addNS('path', 'svg'):
                    if "gcodetools" not in i.keys():
                        self.paths[layer] = self.paths[layer] + [i] if layer in self.paths else [i]
                        if i.get("id") in self.svg.selected:
                            self.svg.selected_paths[layer] = self.svg.selected_paths[layer] + [
                                i] if layer in self.svg.selected_paths else [i]
                elif i.tag == inkex.addNS("g", 'svg'):
                    recursive_search(i, layer, (i.get("id") in self.svg.selected))
                elif i.get("id") in self.svg.selected:
                    self.error("This extension works with Paths and Dynamic Offsets and groups of them only! All other objects will be ignored!\nSolution 1: press Path->Object to path or Shift+Ctrl+C.\nSolution 2: Path->Dynamic offset or Ctrl+J.\nSolution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file.")

        recursive_search(self.document.getroot(), self.document.getroot())

################################################################################
# Fill area
################################################################################

    def area_fill(self):

        global gcode
        global csp

        self.options.area_fill_angle = self.options.area_fill_angle * math.pi / 180

        print_("===================================================================")
        print_("Start filling area", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")

        if len(self.svg.selected_paths) <= 0:
            self.error("This extension requires at least one selected path.")
            return
        if not self.check_dir():
            return

        for layer in self.layers:
            if layer in self.svg.selected_paths:
                if self.options.laser_beam_with <= 0:
                    self.error("Laser beam with must be > 0!")

                print_time("Time until path selection")

                for path in self.svg.selected_paths[layer]:
                    lines = []

                    print_("")
                    print_("Working on path: ")
                    print_debug(path.get("style"), path.get("d"))

                    area_group = path.getparent().add(Group())
                    csp = path.path.to_superpath()
                    if not csp:
                        print_("omitting non-path")
                        self.error("Warning: omitting non-path")
                        continue

                    csp = self.transform_csp(csp, layer)

                    print_debug("csp length: ", len(csp))
                    print_time("Time for csp transformation")

                    # rotate the path to get bounds in defined direction.
                    a = - self.options.area_fill_angle
                    rotated_path = [[[[point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(
                        a)+point[1]*math.cos(a)] for point in sp] for sp in subpath] for subpath in csp]

                    # Draw the lines
                    # Get path's bounds
                    b = [0.0, 0.0, 0.0, 0.0]
                    b = csp_bounds(rotated_path)

                    # Zig-zag
                    r = self.options.laser_beam_with
                    if r <= 0:
                        self.error("Laser diameter must be greater than 0!", "error")
                        return

                    lines += [[]]

                    i = b[0] - r + 0.01
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

                    # Rotate created paths back
                    a = self.options.area_fill_angle
                    lines = [[[point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(
                        a)+point[1]*math.cos(a)] for point in subpath] for subpath in lines]

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
                                p = csp_line_intersection(l1, l2, sp1, sp2)

                                if len(p) > 1:

                                    if l1[0] == l2[0]:
                                        t1 = (p[1]-l1[1])/(l2[1]-l1[1])
                                    else:
                                        t1 = (p[0]-l1[0])/(l2[0]-l1[0])

                                    if 0 <= round(t1, 12) < 1:
                                        ints += [[t1, p[0], p[1]]]

                        ints.sort()

                        for i in ints:
                            splitted_line[-1] += [[i[1], i[2]]]
                            splitted_line += [[[i[1], i[2]]]]
                        splitted_line[-1] += [l2]
                        i = 0

                    print_time("Time for calculating intersections")
                    print_debug("number of splitted lines: ", len(splitted_line))

                    finalLines = []

                    # TODO: fix for Windows Systems. Causes infinite loop due to lack of Fork
                    if self.options.multi_thread and os.name != 'nt':
                        with Pool() as pool:

                            splitted_line_csp = zip(splitted_line, [csp] * len(splitted_line))
                            finalLines = pool.map(check_if_line_inside_shape, splitted_line_csp)  # 13s; s1:57

                    else:
                        # splitted_line_csp = zip(splitted_line, [csp] * len(splitted_line))
                        while i < len(splitted_line):
                            finalLines += [check_if_line_inside_shape([splitted_line[i], csp])]
                            i += 1

                    i = 0

                    print_time("Time for checking if line is insied of shape")
                    print_debug("number of final lines before removing emptys: ", len(finalLines))
                    # remove empty elements
                    # print_("final_line: ", finalLines)
                    np_finalLines = np.array(finalLines, dtype=np.float32)
                    index_zeros = np.argwhere(np_finalLines == [[0, 0], [0, 0]])
                    np_finalLines = np.delete(np_finalLines, index_zeros, axis=0)

                    print_debug("number of final lines: ", len(np_finalLines))

                    if options.remove_tiny_infill_paths:
                        start_coords = np.array(np_finalLines[:, 0])
                        end_coords = np.array(np_finalLines[:, 1])

                        distances = np.array(end_coords[:, 1]-start_coords[:, 1])
                        distances = np.abs(distances)
                        np_finalLines = (np_finalLines[distances > (TINY_INFILL_FACTOR * options.laser_beam_with)])
                        # print_("final_line: ", np_finalLines)

                    print_time("Time for calculating infill paths")

                    csp_line = csp_from_polyline(np_finalLines)

                    # TODO: Remove
                    # curve = self.parse_curve2d(csp_line, layer)
                    # self.draw_curve(curve, layer, area_group)

                    csp_line = self.transform_csp(csp_line, layer, True)

                    print_time("Time for transforming infill paths")

                    print_("csp line 1:")
                    for csp in csp_line:
                        print_(csp)

                    curve = self.parse_curve2d(csp_line, layer)
                    self.draw_curve(curve, layer, area_group)

                    print_time("Time for drawing curve")

                    gcode += self.generate_gcode(curve, layer, self.tool_infill, "infill")

                    print_time("Time for generating Gcode")

        if gcode != '' and not self.options.add_contours:
            self.export_gcode(gcode)

        print_("===================================================================")
        print_("Finished filling area", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")

################################################################################
# Engraving
################################################################################
    def engraving(self):
        global cspm

        # TODO: Candidate for refactoring
        def bisect(nx1, ny1, nx2, ny2):

            cosBis = math.sqrt(max(0, (1.0+nx1*nx2-ny1*ny2)/2.0))
            # We can get correct sign of the sin, assuming cos is positive
            if (abs(ny1-ny2) < ENGRAVING_TOLERANCE) or (abs(cosBis) < ENGRAVING_TOLERANCE):
                if abs(nx1-nx2) < ENGRAVING_TOLERANCE:
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

        def save_point_new(x, y, i, j):

            global cspm

            x = round(x, 3)  # round to 3 decimals
            y = round(y, 3)  # round to 3 decimals

            # if len(cspm) > 1:
            #     cspm += [[[x, y], [x, y], [x, y]]]

            cspm += [[[x, y], [x, y], [x, y]]]

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
                            xydist = abs((xy2[0]-x)*(xy1[1]-y)-(xy1[0]-x)*(xy2[1]-y))/length2
                            if xydist < ENGRAVING_TOLERANCE:  # so far so good
                                cspm.pop()
            cspm += [[[x, y], [x, y], [x, y], i, j]]

        # end of subfunction definitions. engraving() starts here:
        ###########################################################

        print_("===================================================================")
        print_("Start doing perimeters", time.strftime("%d.%m.%Y %H:%M:%S"))
        print_("===================================================================")
        timestamp2 = time.time()

        global gcode
        r = 0  # theoretical and tool-radius-limited radii in pixels
        x1, y1, nx, ny = 0, 0, 0, 0

        cspe = []

        if len(self.svg.selected_paths) <= 0:
            self.error("Please select at least one path to engrave and run again.")
            return
        if self.options.add_infill == False:
            if not self.check_dir():
                return

        for layer in self.layers:
            if layer in self.svg.selected_paths:
                # Calculate scale in pixels per user unit (mm or inch)

                engraving_group = self.svg.selected_paths[layer][0].getparent().add(Group())

                for node in self.svg.selected_paths[layer]:

                    print_("")
                    print_("Working on path: ")
                    print_debug(node.get("style"), node.get("d"))

                    if node.tag == inkex.addNS('path', 'svg'):
                        cspi = node.path.to_superpath()
                        nlLT = []

                        for csp in cspi:
                            nlLT.append([])

                            for i in range(0, len(csp)):  # LT for each point
                                sp0, sp1, sp2 = csp[i-2], csp[i-1], csp[i]
                                # LT find angle between this and previous segment
                                x0, y0 = sp1[1]
                                nx1, ny1 = csp_normalized_normal(sp1, sp2)
                                nx0, ny0 = csp_normalized_normal(sp0, sp1)

                                s = bisect(nx0, ny0, nx1, ny1)[2]

                                # record x,y,normal,ifCorner, sin(angle-turned/2)
                                nlLT[-1] += [[[x0, y0], [0, 0], True, s]]

                        reflex = False
                        for j in xrange(len(nlLT)):  # LT6b for each subpath
                            cspm = []  # Will be my output. List of csps.
                            for i in xrange(len(nlLT[j])):  # LT for each node
                                n1 = nlLT[j][i-1]  # current node
                                # this point/start of this line
                                x1, y1 = n1[0]

                                if n1[2]:  # We're at a corner
                                    bits = 1

                                for b in xrange(bits):  # divide line into bits

                                    if n1[2]:  # We're at a corner
                                        if n1[3] > 0:  # acute
                                            save_point(x1, y1, i, j)
                                            save_point(x1, y1, i, j)

                                    save_point(x1, y1, i, j)

                            # LT next i
                            if len(cspm) != 0:
                                cspm += [cspm[0]]
                                cspe += [cspm]

                if cspe != []:

                    # print_("cspe:")
                    # print_(cspe)

                    curve = self.parse_curve(cspe, layer)
                    self.draw_curve(curve, layer, engraving_group)

                    gcode += self.generate_gcode(curve, layer, self.tool_perimeter, "perimeter")

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
        global orient_points
        self.get_info()

        if layer is None:
            layer = self.svg.get_current_layer() if self.svg.get_current_layer() is not None else self.document.getroot()

        transform = self.get_transforms(layer)

        if transform:
            transform = self.reverse_transform(transform)
            transform = str(Transform(transform))

        print_("Inserting orientation points")

        if layer in self.orientation_points:
            self.error("Active layer already has orientation points! Remove them or select another layer!", "error")

        attr = {"gcodetools": "Gcodetools orientation group"}
        if transform:
            attr["transform"] = transform

        orientation_group = layer.add(Group(**attr))

        doc_height = self.svg.unittouu(self.document.getroot().get('height'))
        if self.document.getroot().get('height') == "100%":
            doc_height = 1052.3622047
            print_("Overriding height from 100 percents to {}".format(doc_height))

        orient_points = [[[100, doc_height], [100., 0.0, 0.0]], [[0.0, doc_height], [0.0, 0.0, 0.0]]]

    ################################################################################
    # ApplyTransform
    ################################################################################

    @staticmethod
    def objectToPath(node):

        if node.tag == inkex.addNS('g', 'svg'):
            return node

        if node.tag == inkex.addNS('path', 'svg') or node.tag == 'path':
            for attName in node.attrib.keys():
                if ("sodipodi" in attName) or ("inkscape" in attName):
                    del node.attrib[attName]
            return node

        return node

    def recursiveFuseTransform(self, node, transf=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]):

        transf = Transform(transf) * Transform(node.get("transform", None))

        if 'transform' in node.attrib:
            del node.attrib['transform']

        if 'style' in node.attrib:
            style = node.attrib.get('style')
            style = dict(Style.parse_str(style))
            update = False

            if 'stroke-width' in style:
                try:
                    stroke_width = self.unittouu(style.get('stroke-width').strip())
                    stroke_width *= math.hypot(transf[0][0], transf[1][1])
                    style['stroke-width'] = str(stroke_width)
                    update = True
                except AttributeError:
                    pass

            if update:
                node.attrib['style'] = Style(style).to_str()

        node = self.objectToPath(node)

        if 'd' in node.attrib:
            d = node.get('d')
            p = CubicSuperPath(d)
            p = Path(p).to_absolute().transform(transf, True)
            node.set('d', Path(CubicSuperPath(p).to_path()))

        elif node.tag in [inkex.addNS('polygon', 'svg'), inkex.addNS('polyline', 'svg')]:
            points = node.get('points')
            points = points.strip().split(' ')
            for k, p in enumerate(points):
                if ',' in p:
                    p = p.split(',')
                    p = [float(p[0]), float(p[1])]
                    Transform.apply_to_point(transf, p)
                    p = [str(p[0]), str(p[1])]
                    p = ','.join(p)
                    points[k] = p
            points = ' '.join(points)
            node.set('points', points)

        elif node.tag in [inkex.addNS('rect', 'svg'),
                          inkex.addNS('text', 'svg'),
                          inkex.addNS('image', 'svg'),
                          inkex.addNS('use', 'svg'),
                          inkex.addNS('circle', 'svg')]:
            node.set('transform', str(Transform(transf)))

        for child in node.getchildren():
            self.recursiveFuseTransform(child, transf)

    def applytransforms(self):
        # Apply transformations
        self.svg.get_selected()

        if self.svg.selected:
            for id, shape in self.svg.selected.items():
                self.recursiveFuseTransform(shape)
        else:
            self.recursiveFuseTransform(self.document.getroot())
        # Transformations applied

    def flatten(self, tolerance=0.1):
        print_("Elements in path: ", self.svg.get_selected(inkex.PathElement))

        for layer in self.layers:
            if layer in self.svg.selected_paths:
                for node in self.svg.selected_paths[layer]:
                    path = node.path.to_superpath()
                    bezier.cspsubdiv(path, tolerance)
                    newpath = []
                    for subpath in path:
                        first = True
                        for csp in subpath:
                            cmd = 'L'
                            if first:
                                cmd = 'M'
                            first = False
                            newpath.append([cmd, [csp[1][0], csp[1][1]]])
                    node.path = newpath

################################################################################
# Effect
# Main function of Lasertools class
################################################################################

    def effect(self):
        global options
        options = self.options
        options.self = self
        global print_

        if self.options.log_create_log:
            if os.path.isdir(self.options.directory):
                try:
                    if os.path.isfile(self.options.directory+"/log.txt"):
                        os.remove(self.options.directory+"/log.txt")
                    f = open(self.options.directory+"/log.txt", "a")
                    f.write("===================================================================\n")
                    f.write("Lasertools log file.\nStarted at %s.\n%s\n" % (time.strftime("%d.%m.%Y %H:%M:%S"), options.directory+"/log.txt"))
                    f.write("===================================================================\n\n")
                    f.close()
                except:
                    print_ = lambda *x: None

            else:
                self.error(("Directory does not exist! Please specify existing directory at options tab!"), "error")

        else:
            print_ = lambda *x: None
        self.get_info()

        # wait to attatch debugger to process
        if DEBUG:
            print_debug("Python version:", sys.version_info)
            print_debug("Waiting for Debugger to be attached")
            time.sleep(3)
            print_debug("Starting Program")

        if self.orientation_points == {}:
            self.orientation(self.layers[min(0, len(self.layers)-1)])
            print_debug("Orientation Points: ", self.orientation_points)
            self.get_info()
        # handle power on delay
        delayOn = ""
        if round(float(self.options.power_delay), 1) > 0:
            delayOn = "G04 P" + str(round(float(self.options.power_delay)/1000, 3)) + "\n"

        self.tool_infill = {
            "name": "Laser Engraver Infill",
            "id": "Laser Engraver Infill",
            "penetration feed": self.options.laser_speed,
            "feed": self.options.laser_speed,
            "gcode before path": (delayOn + self.options.laser_command),
            "gcode after path": self.options.laser_off_command
        }

        self.tool_perimeter = {
            "name": "Laser Engraver Perimeter",
            "id": "Laser Engraver Perimeter",
            "penetration feed": self.options.laser_param_speed,
            "feed": self.options.laser_param_speed,
            "gcode before path": (delayOn + self.options.laser_command_perimeter),
            "gcode after path": self.options.laser_off_command
        }

        self.get_info()

        print_("Applying all transformations")
        self.applytransforms()

        print_("Flattening beziers")
        self.svg.selected_paths = self.paths
        self.flatten(self.options.tolerance)

        if self.options.add_infill:
            self.area_fill()

        if self.options.add_contours:
            self.get_info()
            self.svg.selected_paths = self.paths

            if PROFILING:
                if os.path.isfile(self.options.directory+"performance.prof"):
                    os.remove(self.options.directory+"/performance.prof")

                profile = cProfile.Profile()
                profile.runctx('self.engraving()', globals(), locals())

                kProfile = lsprofcalltree.KCacheGrind(profile)

                kFile = open(self.options.directory+"/performance.prof", 'w+')
                kProfile.output(kFile)
                kFile.close()

            self.engraving()


if __name__ == '__main__':
    laser_gcode().run()
