#!/usr/bin/env python
"""
Modified by Christoph Wagner 2015
Modified by Jay Johnson 2015, J Tech Photonics, Inc., jtechphotonics.com
modified by Adam Polak 2014, polakiumengineering.org

based on Copyright (C) 2009 Nick Drobchenko, nick@cnc-club.ru
based on gcode.py (C) 2007 hugomatic...
based on addnodes.py (C) 2005,2007 Aaron Spike, aaron@ekips.org
based on dots.py (C) 2005 Aaron Spike, aaron@ekips.org
based on interp.py (C) 2005 Aaron Spike, aaron@ekips.org
based on bezmisc.py (C) 2005 Aaron Spike, aaron@ekips.org
based on cubicsuperpath.py (C) 2005 Aaron Spike, aaron@ekips.org

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
import inkex, simplestyle, simplepath
import cubicsuperpath, simpletransform, bezmisc

import os
import math
import bezmisc
import re
import copy
import sys
import time
import cmath
import numpy
import codecs
import random
import gettext
_ = gettext.gettext


### Check if inkex has errormsg (0.46 version doesnot have one.) Could be removed later.
if "errormsg" not in dir(inkex):
    inkex.errormsg = lambda msg: sys.stderr.write((unicode(msg) + "\n").encode("UTF-8"))


################################################################################
###
###        Styles and additional parameters
###
################################################################################

gcode = ""
math.pi2 = math.pi*2
straight_tolerance = 0.00001
straight_distance_tolerance = 0.00001
engraving_tolerance = 0.00001
loft_lengths_tolerance = 0.00000001
options = {}
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

intersection_recursion_depth = 10
intersection_tolerance = 0.00001

styles = {
        "loft_style" : {
            'main curve':    simplestyle.formatStyle({'stroke': '#88f', 'fill': 'none', 'stroke-width': '0.3', 'marker-end': 'none'}),
            },
        "biarc_style" : {
                'biarc0':    simplestyle.formatStyle({ 'stroke': '#88f', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'biarc1':    simplestyle.formatStyle({ 'stroke': '#8f8', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'line':        simplestyle.formatStyle({ 'stroke': '#f88', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'area':        simplestyle.formatStyle({ 'stroke': '#777', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.1' }),
            },
        "biarc_style_dark" : {
                'biarc0':    simplestyle.formatStyle({ 'stroke': '#33a', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'biarc1':    simplestyle.formatStyle({ 'stroke': '#3a3', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'line':        simplestyle.formatStyle({ 'stroke': '#a33', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
                'area':        simplestyle.formatStyle({ 'stroke': '#222', 'fill': 'none', "marker-end":'none' , 'stroke-width':'0.3' }),
            },
        "dxf_points":             simplestyle.formatStyle({ "stroke": "#ff0000", "fill": "#ff0000"}),
    }

################################################################################
###        Cubic Super Path additional functions
################################################################################

def csp_segment_to_bez(sp1,sp2) :
    return sp1[1:]+sp2[:2]

def csp_to_point_distance(csp, p, dist_bounds = [0,1e100], tolerance=.0001) :
    min_dist = [1e100,0,0,0]
    for j in range(len(csp)) :
        for i in range(1,len(csp[j])) :
            d = csp_seg_to_point_distance(csp[j][i-1],csp[j][i],p,sample_points = 10, tolerance = .0001)
            if d[0] < dist_bounds[0] :
                return [d[0],j,i,d[1]]
            else :
                if d[0] < min_dist[0] : min_dist = [d[0],j,i,d[1]]
    return min_dist

def csp_seg_to_point_distance(sp1,sp2,p,sample_points = 10, tolerance = .0001) :
    ax,ay,bx,by,cx,cy,dx,dy = csp_parameterize(sp1,sp2)
    dx, dy = dx-p[0], dy-p[1]
    if sample_points < 2 : sample_points = 2
    d = min( [(p[0]-sp1[1][0])**2 + (p[1]-sp1[1][1])**2,0.], [(p[0]-sp2[1][0])**2 + (p[1]-sp2[1][1])**2,1.]    )
    for k in range(sample_points) :
        t = float(k)/(sample_points-1)
        i = 0
        while i==0 or abs(f)>0.000001 and i<20 :
            t2,t3 = t**2,t**3
            f = (ax*t3+bx*t2+cx*t+dx)*(3*ax*t2+2*bx*t+cx) + (ay*t3+by*t2+cy*t+dy)*(3*ay*t2+2*by*t+cy)
            df = (6*ax*t+2*bx)*(ax*t3+bx*t2+cx*t+dx) + (3*ax*t2+2*bx*t+cx)**2 + (6*ay*t+2*by)*(ay*t3+by*t2+cy*t+dy) + (3*ay*t2+2*by*t+cy)**2
            if df!=0 :
                t = t - f/df
            else :
                break
            i += 1
        if 0<=t<=1 :
            p1 = csp_at_t(sp1,sp2,t)
            d1 = (p1[0]-p[0])**2 + (p1[1]-p[1])**2
            if d1 < d[0] :
                d = [d1,t]
    return d


def csp_seg_to_csp_seg_distance(sp1,sp2,sp3,sp4, dist_bounds = [0,1e100], sample_points = 10, tolerance=.0001) :
    # check the ending points first
    dist =    csp_seg_to_point_distance(sp1,sp2,sp3[1],sample_points, tolerance)
    dist += [0.]
    if dist[0] <= dist_bounds[0] : return dist
    d = csp_seg_to_point_distance(sp1,sp2,sp4[1],sample_points, tolerance)
    if d[0]<dist[0] :
        dist = d+[1.]
        if dist[0] <= dist_bounds[0] : return dist
    d =    csp_seg_to_point_distance(sp3,sp4,sp1[1],sample_points, tolerance)
    if d[0]<dist[0] :
        dist = [d[0],0.,d[1]]
        if dist[0] <= dist_bounds[0] : return dist
    d =    csp_seg_to_point_distance(sp3,sp4,sp2[1],sample_points, tolerance)
    if d[0]<dist[0] :
        dist = [d[0],1.,d[1]]
        if dist[0] <= dist_bounds[0] : return dist
    sample_points -= 2
    if sample_points < 1 : sample_points = 1
    ax1,ay1,bx1,by1,cx1,cy1,dx1,dy1 = csp_parameterize(sp1,sp2)
    ax2,ay2,bx2,by2,cx2,cy2,dx2,dy2 = csp_parameterize(sp3,sp4)
    #    try to find closes points using Newtons method
    for k in range(sample_points) :
        for j in range(sample_points) :
            t1,t2 = float(k+1)/(sample_points+1), float(j)/(sample_points+1)
            t12, t13, t22, t23 = t1*t1, t1*t1*t1, t2*t2, t2*t2*t2
            i = 0
            F1, F2, F = [0,0], [[0,0],[0,0]], 1e100
            x,y   = ax1*t13+bx1*t12+cx1*t1+dx1 - (ax2*t23+bx2*t22+cx2*t2+dx2), ay1*t13+by1*t12+cy1*t1+dy1 - (ay2*t23+by2*t22+cy2*t2+dy2)
            while i<2 or abs(F-Flast)>tolerance and i<30 :
                #draw_pointer(csp_at_t(sp1,sp2,t1))
                f1x = 3*ax1*t12+2*bx1*t1+cx1
                f1y = 3*ay1*t12+2*by1*t1+cy1
                f2x = 3*ax2*t22+2*bx2*t2+cx2
                f2y = 3*ay2*t22+2*by2*t2+cy2
                F1[0] = 2*f1x*x +  2*f1y*y
                F1[1] = -2*f2x*x -  2*f2y*y
                F2[0][0] =  2*(6*ax1*t1+2*bx1)*x + 2*f1x*f1x + 2*(6*ay1*t1+2*by1)*y +2*f1y*f1y
                F2[0][1] = -2*f1x*f2x - 2*f1y*f2y
                F2[1][0] = -2*f2x*f1x - 2*f2y*f1y
                F2[1][1] = -2*(6*ax2*t2+2*bx2)*x + 2*f2x*f2x - 2*(6*ay2*t2+2*by2)*y + 2*f2y*f2y
                F2 = inv_2x2(F2)
                if F2!=None :
                    t1 -= ( F2[0][0]*F1[0] + F2[0][1]*F1[1] )
                    t2 -= ( F2[1][0]*F1[0] + F2[1][1]*F1[1] )
                    t12, t13, t22, t23 = t1*t1, t1*t1*t1, t2*t2, t2*t2*t2
                    x,y   = ax1*t13+bx1*t12+cx1*t1+dx1 - (ax2*t23+bx2*t22+cx2*t2+dx2), ay1*t13+by1*t12+cy1*t1+dy1 - (ay2*t23+by2*t22+cy2*t2+dy2)
                    Flast = F
                    F = x*x+y*y
                else :
                    break
                i += 1
            if F < dist[0] and 0<=t1<=1 and 0<=t2<=1:
                dist = [F,t1,t2]
                if dist[0] <= dist_bounds[0] :
                    return dist
    return dist


def csp_to_csp_distance(csp1,csp2, dist_bounds = [0,1e100], tolerance=.0001) :
    dist = [1e100,0,0,0,0,0,0]
    for i1 in range(len(csp1)) :
        for j1 in range(1,len(csp1[i1])) :
            for i2 in range(len(csp2)) :
                for j2 in range(1,len(csp2[i2])) :
                    d = csp_seg_bound_to_csp_seg_bound_max_min_distance(csp1[i1][j1-1],csp1[i1][j1],csp2[i2][j2-1],csp2[i2][j2])
                    if d[0] >= dist_bounds[1] : continue
                    if  d[1] < dist_bounds[0] : return [d[1],i1,j1,1,i2,j2,1]
                    d = csp_seg_to_csp_seg_distance(csp1[i1][j1-1],csp1[i1][j1],csp2[i2][j2-1],csp2[i2][j2], dist_bounds, tolerance=tolerance)
                    if d[0] < dist[0] :
                        dist = [d[0], i1,j1,d[1], i2,j2,d[2]]
                    if dist[0] <= dist_bounds[0] :
                        return dist
            if dist[0] >= dist_bounds[1] :
                return dist
    return dist

def csp_split(sp1,sp2,t=.5) :
    [x1,y1],[x2,y2],[x3,y3],[x4,y4] = sp1[1], sp1[2], sp2[0], sp2[1]
    x12 = x1+(x2-x1)*t
    y12 = y1+(y2-y1)*t
    x23 = x2+(x3-x2)*t
    y23 = y2+(y3-y2)*t
    x34 = x3+(x4-x3)*t
    y34 = y3+(y4-y3)*t
    x1223 = x12+(x23-x12)*t
    y1223 = y12+(y23-y12)*t
    x2334 = x23+(x34-x23)*t
    y2334 = y23+(y34-y23)*t
    x = x1223+(x2334-x1223)*t
    y = y1223+(y2334-y1223)*t
    return [sp1[0],sp1[1],[x12,y12]], [[x1223,y1223],[x,y],[x2334,y2334]], [[x34,y34],sp2[1],sp2[2]]

def csp_true_bounds(csp) :

    # Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t)
    minx = [float("inf"), 0, 0, 0]
    maxx = [float("-inf"), 0, 0, 0]
    miny = [float("inf"), 0, 0, 0]
    maxy = [float("-inf"), 0, 0, 0]
    for i in range(len(csp)):
        for j in range(1,len(csp[i])):
            ax,ay,bx,by,cx,cy,x0,y0 = bezmisc.bezierparameterize((csp[i][j-1][1],csp[i][j-1][2],csp[i][j][0],csp[i][j][1]))
            roots = cubic_solver(0, 3*ax, 2*bx, cx)     + [0,1]
            for root in roots :
                if type(root) is complex and abs(root.imag)<1e-10:
                    root = root.real
                if type(root) is not complex and 0<=root<=1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxx = max([x,y,i,j,root],maxx)
                    minx = min([x,y,i,j,root],minx)

            roots = cubic_solver(0, 3*ay, 2*by, cy)     + [0,1]
            for root in roots :
                if type(root) is complex and root.imag==0:
                    root = root.real
                if type(root) is not complex and 0<=root<=1:
                    y = ay*(root**3)+by*(root**2)+cy*root+y0
                    x = ax*(root**3)+bx*(root**2)+cx*root+x0
                    maxy = max([y,x,i,j,root],maxy)
                    miny = min([y,x,i,j,root],miny)
    maxy[0],maxy[1] = maxy[1],maxy[0]
    miny[0],miny[1] = miny[1],miny[0]

    return minx,miny,maxx,maxy


############################################################################
### csp_segments_intersection(sp1,sp2,sp3,sp4)
###
### Returns array containig all intersections between two segmets of cubic
### super path. Results are [ta,tb], or [ta0, ta1, tb0, tb1, "Overlap"]
### where ta, tb are values of t for the intersection point.
############################################################################
def csp_segments_intersection(sp1,sp2,sp3,sp4) :
    a, b = csp_segment_to_bez(sp1,sp2), csp_segment_to_bez(sp3,sp4)

    def polish_intersection(a,b,ta,tb, tolerance = intersection_tolerance) :
        ax,ay,bx,by,cx,cy,dx,dy            = bezmisc.bezierparameterize(a)
        ax1,ay1,bx1,by1,cx1,cy1,dx1,dy1    = bezmisc.bezierparameterize(b)
        i = 0
        F, F1 =  [.0,.0], [[.0,.0],[.0,.0]]
        while i==0 or (abs(F[0])**2+abs(F[1])**2 > tolerance and i<10):
            ta3, ta2, tb3, tb2 = ta**3, ta**2, tb**3, tb**2
            F[0] = ax*ta3+bx*ta2+cx*ta+dx-ax1*tb3-bx1*tb2-cx1*tb-dx1
            F[1] = ay*ta3+by*ta2+cy*ta+dy-ay1*tb3-by1*tb2-cy1*tb-dy1
            F1[0][0] =  3*ax *ta2 + 2*bx *ta + cx
            F1[0][1] = -3*ax1*tb2 - 2*bx1*tb - cx1
            F1[1][0] =  3*ay *ta2 + 2*by *ta + cy
            F1[1][1] = -3*ay1*tb2 - 2*by1*tb - cy1
            det = F1[0][0]*F1[1][1] - F1[0][1]*F1[1][0]
            if det!=0 :
                F1 = [    [ F1[1][1]/det, -F1[0][1]/det],    [-F1[1][0]/det,  F1[0][0]/det] ]
                ta = ta - ( F1[0][0]*F[0] + F1[0][1]*F[1] )
                tb = tb - ( F1[1][0]*F[0] + F1[1][1]*F[1] )
            else: break
            i += 1

        return ta, tb


    def recursion(a,b, ta0,ta1,tb0,tb1, depth_a,depth_b) :
        global bezier_intersection_recursive_result
        if a==b :
            bezier_intersection_recursive_result += [[ta0,tb0,ta1,tb1,"Overlap"]]
            return
        tam, tbm = (ta0+ta1)/2, (tb0+tb1)/2
        if depth_a>0 and depth_b>0 :
            a1,a2 = bez_split(a,0.5)
            b1,b2 = bez_split(b,0.5)
            if bez_bounds_intersect(a1,b1) : recursion(a1,b1, ta0,tam,tb0,tbm, depth_a-1,depth_b-1)
            if bez_bounds_intersect(a2,b1) : recursion(a2,b1, tam,ta1,tb0,tbm, depth_a-1,depth_b-1)
            if bez_bounds_intersect(a1,b2) : recursion(a1,b2, ta0,tam,tbm,tb1, depth_a-1,depth_b-1)
            if bez_bounds_intersect(a2,b2) : recursion(a2,b2, tam,ta1,tbm,tb1, depth_a-1,depth_b-1)
        elif depth_a>0  :
            a1,a2 = bez_split(a,0.5)
            if bez_bounds_intersect(a1,b) : recursion(a1,b, ta0,tam,tb0,tb1, depth_a-1,depth_b)
            if bez_bounds_intersect(a2,b) : recursion(a2,b, tam,ta1,tb0,tb1, depth_a-1,depth_b)
        elif depth_b>0  :
            b1,b2 = bez_split(b,0.5)
            if bez_bounds_intersect(a,b1) : recursion(a,b1, ta0,ta1,tb0,tbm, depth_a,depth_b-1)
            if bez_bounds_intersect(a,b2) : recursion(a,b2, ta0,ta1,tbm,tb1, depth_a,depth_b-1)
        else : # Both segments have been subdevided enougth. Let's get some intersections :).
            intersection, t1, t2 =  straight_segments_intersection([a[0]]+[a[3]],[b[0]]+[b[3]])
            if intersection :
                if intersection == "Overlap" :
                    t1 = ( max(0,min(1,t1[0]))+max(0,min(1,t1[1])) )/2
                    t2 = ( max(0,min(1,t2[0]))+max(0,min(1,t2[1])) )/2
                bezier_intersection_recursive_result += [[ta0+t1*(ta1-ta0),tb0+t2*(tb1-tb0)]]

    global bezier_intersection_recursive_result
    bezier_intersection_recursive_result = []
    recursion(a,b,0.,1.,0.,1.,intersection_recursion_depth,intersection_recursion_depth)
    intersections = bezier_intersection_recursive_result
    for i in range(len(intersections)) :
        if len(intersections[i])<5 or intersections[i][4] != "Overlap" :
            intersections[i] = polish_intersection(a,b,intersections[i][0],intersections[i][1])
    return intersections


def csp_segments_true_intersection(sp1,sp2,sp3,sp4) :
    intersections = csp_segments_intersection(sp1,sp2,sp3,sp4)
    res = []
    for intersection in intersections :
        if  (
                (len(intersection)==5 and intersection[4] == "Overlap" and (0<=intersection[0]<=1 or 0<=intersection[1]<=1) and (0<=intersection[2]<=1 or 0<=intersection[3]<=1) )
             or ( 0<=intersection[0]<=1 and 0<=intersection[1]<=1 )
            ) :
            res += [intersection]
    return res


def csp_get_t_at_curvature(sp1,sp2,c, sample_points = 16):
    # returns a list containning [t1,t2,t3,...,tn],  0<=ti<=1...
    if sample_points < 2 : sample_points = 2
    tolerance = .0000000001
    res = []
    ax,ay,bx,by,cx,cy,dx,dy = csp_parameterize(sp1,sp2)
    for k in range(sample_points) :
        t = float(k)/(sample_points-1)
        i, F = 0, 1e100
        while i<2 or abs(F)>tolerance and i<17 :
            try : # some numerical calculation could exceed the limits
                t2 = t*t
                #slopes...
                f1x = 3*ax*t2+2*bx*t+cx
                f1y = 3*ay*t2+2*by*t+cy
                f2x = 6*ax*t+2*bx
                f2y = 6*ay*t+2*by
                f3x = 6*ax
                f3y = 6*ay
                d = (f1x**2+f1y**2)**1.5
                F1 = (
                         (    (f1x*f3y-f3x*f1y)*d - (f1x*f2y-f2x*f1y)*3.*(f2x*f1x+f2y*f1y)*((f1x**2+f1y**2)**.5) )    /
                                 ((f1x**2+f1y**2)**3)
                     )
                F = (f1x*f2y-f1y*f2x)/d - c
                t -= F/F1
            except:
                break
            i += 1
        if 0<=t<=1 and F<=tolerance:
            if len(res) == 0 :
                res.append(t)
            for i in res :
                if abs(t-i)<=0.001 :
                    break
            if not abs(t-i)<=0.001 :
                res.append(t)
    return res


def csp_max_curvature(sp1,sp2):
    ax,ay,bx,by,cx,cy,dx,dy = csp_parameterize(sp1,sp2)
    tolerance = .000001
    F = 0.
    i = 0
    while i<2 or F-Flast<tolerance and i<10 :
        t = .5
        f1x = 3*ax*t**2 + 2*bx*t + cx
        f1y = 3*ay*t**2 + 2*by*t + cy
        f2x = 6*ax*t + 2*bx
        f2y = 6*ay*t + 2*by
        f3x = 6*ax
        f3y = 6*ay
        d = pow(f1x**2+f1y**2,1.5)
        if d != 0 :
            Flast = F
            F = (f1x*f2y-f1y*f2x)/d
            F1 =     (
                         (    d*(f1x*f3y-f3x*f1y) - (f1x*f2y-f2x*f1y)*3.*(f2x*f1x+f2y*f1y)*pow(f1x**2+f1y**2,.5) )    /
                                 (f1x**2+f1y**2)**3
                    )
            i+=1
            if F1!=0:
                t -= F/F1
            else:
                break
        else: break
    return t


def csp_curvature_at_t(sp1,sp2,t, depth = 3) :
    ax,ay,bx,by,cx,cy,dx,dy = bezmisc.bezierparameterize(csp_segment_to_bez(sp1,sp2))

    #curvature = (x'y''-y'x'') / (x'^2+y'^2)^1.5

    f1x = 3*ax*t**2 + 2*bx*t + cx
    f1y = 3*ay*t**2 + 2*by*t + cy
    f2x = 6*ax*t + 2*bx
    f2y = 6*ay*t + 2*by
    d = (f1x**2+f1y**2)**1.5
    if d != 0 :
        return (f1x*f2y-f1y*f2x)/d
    else :
        t1 = f1x*f2y-f1y*f2x
        if t1 > 0 : return 1e100
        if t1 < 0 : return -1e100
        # Use the Lapitals rule to solve 0/0 problem for 2 times...
        t1 = 2*(bx*ay-ax*by)*t+(ay*cx-ax*cy)
        if t1 > 0 : return 1e100
        if t1 < 0 : return -1e100
        t1 = bx*ay-ax*by
        if t1 > 0 : return 1e100
        if t1 < 0 : return -1e100
        if depth>0 :
            # little hack ;^) hope it wont influence anything...
            return csp_curvature_at_t(sp1,sp2,t*1.004, depth-1)
        return 1e100


def csp_curvature_radius_at_t(sp1,sp2,t) :
    c = csp_curvature_at_t(sp1,sp2,t)
    if c == 0 : return 1e100
    else: return 1/c


def csp_special_points(sp1,sp2) :
    # special points = curvature == 0
    ax,ay,bx,by,cx,cy,dx,dy = bezmisc.bezierparameterize((sp1[1],sp1[2],sp2[0],sp2[1]))
    a = 3*ax*by-3*ay*bx
    b = 3*ax*cy-3*cx*ay
    c = bx*cy-cx*by
    roots = cubic_solver(0, a, b, c)
    res = []
    for i in roots :
        if type(i) is complex and i.imag==0:
            i = i.real
        if type(i) is not complex and 0<=i<=1:
            res.append(i)
    return res


def csp_at_t(sp1,sp2,t):
    ax,bx,cx,dx = sp1[1][0], sp1[2][0], sp2[0][0], sp2[1][0]
    ay,by,cy,dy = sp1[1][1], sp1[2][1], sp2[0][1], sp2[1][1]

    x1, y1 = ax+(bx-ax)*t, ay+(by-ay)*t
    x2, y2 = bx+(cx-bx)*t, by+(cy-by)*t
    x3, y3 = cx+(dx-cx)*t, cy+(dy-cy)*t

    x4,y4 = x1+(x2-x1)*t, y1+(y2-y1)*t
    x5,y5 = x2+(x3-x2)*t, y2+(y3-y2)*t

    x,y = x4+(x5-x4)*t, y4+(y5-y4)*t
    return [x,y]


def csp_splitatlength(sp1, sp2, l = 0.5, tolerance = 0.0001):
    bez = (sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:])
    t = bezmisc.beziertatlength(bez, l, tolerance)
    return csp_split(sp1, sp2, t)


def cspseglength(sp1,sp2, tolerance = 0.0001):
    bez = (sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:])
    return bezmisc.bezierlength(bez, tolerance)


def csplength(csp):
    total = 0
    lengths = []
    for sp in csp:
        for i in xrange(1,len(sp)):
            l = cspseglength(sp[i-1],sp[i])
            lengths.append(l)
            total += l
    return lengths, total


def csp_segments(csp):
    l, seg = 0, [0]
    for sp in csp:
        for i in xrange(1,len(sp)):
            l += cspseglength(sp[i-1],sp[i])
            seg += [ l ]

    if l>0 :
        seg = [seg[i]/l for i in xrange(len(seg))]
    return seg,l


def rebuild_csp (csp, segs, s=None):
    # rebuild_csp() adds to csp control points making it's segments looks like segs
    if s==None : s, l = csp_segments(csp)

    if len(s)>len(segs) : return None
    segs = segs[:]
    segs.sort()
    for i in xrange(len(s)):
        d = None
        for j in xrange(len(segs)):
            d = min( [abs(s[i]-segs[j]),j], d) if d!=None else [abs(s[i]-segs[j]),j]
        del segs[d[1]]
    for i in xrange(len(segs)):
        for j in xrange(0,len(s)):
            if segs[i]<s[j] : break
        if s[j]-s[j-1] != 0 :
            t = (segs[i] - s[j-1])/(s[j]-s[j-1])
            sp1,sp2,sp3 = csp_split(csp[j-1],csp[j], t)
            csp = csp[:j-1] + [sp1,sp2,sp3] + csp[j+1:]
            s = s[:j] + [ s[j-1]*(1-t)+s[j]*t   ] + s[j:]
    return csp, s

def csp_line_intersection(l1,l2,sp1,sp2):
    dd=l1[0]
    cc=l2[0]-l1[0]
    bb=l1[1]
    aa=l2[1]-l1[1]
    if aa==cc==0 : return []
    if aa:
        coef1=cc/aa
        coef2=1
    else:
        coef1=1
        coef2=aa/cc
    bez = (sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:])
    ax,ay,bx,by,cx,cy,x0,y0=bezmisc.bezierparameterize(bez)
    a=coef1*ay-coef2*ax
    b=coef1*by-coef2*bx
    c=coef1*cy-coef2*cx
    d=coef1*(y0-bb)-coef2*(x0-dd)
    roots = cubic_solver(a,b,c,d)
    retval = []
    for i in roots :
        if type(i) is complex and abs(i.imag)<1e-7:
            i = i.real
        if type(i) is not complex and -1e-10<=i<=1.+1e-10:
            retval.append(i)
    return retval


def csp_split_by_two_points(sp1,sp2,t1,t2) :
    if t1>t2 : t1, t2 = t2, t1
    if t1 == t2 :
        sp1,sp2,sp3 =  csp_split(sp1,sp2,t)
        return [sp1,sp2,sp2,sp3]
    elif t1 <= 1e-10 and t2 >= 1.-1e-10 :
        return [sp1,sp1,sp2,sp2]
    elif t1 <= 1e-10:
        sp1,sp2,sp3 = csp_split(sp1,sp2,t2)
        return [sp1,sp1,sp2,sp3]
    elif t2 >= 1.-1e-10 :
        sp1,sp2,sp3 = csp_split(sp1,sp2,t1)
        return [sp1,sp2,sp3,sp3]
    else:
        sp1,sp2,sp3 = csp_split(sp1,sp2,t1)
        sp2,sp3,sp4 = csp_split(sp2,sp3,(t2-t1)/(1-t1) )
        return [sp1,sp2,sp3,sp4]


def csp_subpath_split_by_points(subpath, points) :
    # points are [[i,t]...] where i-segment's number
    points.sort()
    points = [[1,0.]] + points + [[len(subpath)-1,1.]]
    parts = []
    for int1,int2 in zip(points,points[1:]) :
        if int1==int2 :
            continue
        if int1[1] == 1. :
            int1[0] += 1
            int1[1] = 0.
        if int1==int2 :
            continue
        if int2[1] == 0. :
            int2[0] -= 1
            int2[1] = 1.
        if int1[0] == 0 and int2[0]==len(subpath)-1:# and small(int1[1]) and small(int2[1]-1) :
            continue
        if int1[0]==int2[0] :    # same segment
            sp = csp_split_by_two_points(subpath[int1[0]-1],subpath[int1[0]],int1[1], int2[1])
            if sp[1]!=sp[2] :
                parts += [   [sp[1],sp[2]]     ]
        else :
            sp5,sp1,sp2 = csp_split(subpath[int1[0]-1],subpath[int1[0]],int1[1])
            sp3,sp4,sp5 = csp_split(subpath[int2[0]-1],subpath[int2[0]],int2[1])
            if int1[0]==int2[0]-1 :
                parts += [    [sp1, [sp2[0],sp2[1],sp3[2]], sp4]  ]
            else :
                parts += [  [sp1,sp2]+subpath[int1[0]+1:int2[0]-1]+[sp3,sp4]  ]
    return parts


def csp_from_arc(start, end, center, r, slope_st) :
    # Creates csp that approximise specified arc
    r = abs(r)
    alpha = (atan2(end[0]-center[0],end[1]-center[1]) - atan2(start[0]-center[0],start[1]-center[1])) % math.pi2

    sectors = int(abs(alpha)*2/math.pi)+1
    alpha_start = atan2(start[0]-center[0],start[1]-center[1])
    cos_,sin_ = math.cos(alpha_start), math.sin(alpha_start)
    k = (4.*math.tan(alpha/sectors/4.)/3.)
    if dot(slope_st , [- sin_*k*r, cos_*k*r]) < 0 :
        if alpha>0 : alpha -= math.pi2
        else: alpha += math.pi2
    if abs(alpha*r)<0.001 :
        return []

    sectors = int(abs(alpha)*2/math.pi)+1
    k = (4.*math.tan(alpha/sectors/4.)/3.)
    result = []
    for i in range(sectors+1) :
        cos_,sin_ = math.cos(alpha_start + alpha*i/sectors), math.sin(alpha_start + alpha*i/sectors)
        sp = [ [], [center[0] + cos_*r, center[1] + sin_*r], [] ]
        sp[0] = [sp[1][0] + sin_*k*r, sp[1][1] - cos_*k*r ]
        sp[2] = [sp[1][0] - sin_*k*r, sp[1][1] + cos_*k*r ]
        result += [sp]
    result[0][0] = result[0][1][:]
    result[-1][2] = result[-1][1]

    return result


def point_to_arc_distance(p, arc):
    ###        Distance calculattion from point to arc
    P0,P2,c,a = arc
    dist = None
    p = P(p)
    r = (P0-c).mag()
    if r>0 :
        i = c + (p-c).unit()*r
        alpha = ((i-c).angle() - (P0-c).angle())
        if a*alpha<0:
            if alpha>0:    alpha = alpha-math.pi2
            else: alpha = math.pi2+alpha
        if between(alpha,0,a) or min(abs(alpha),abs(alpha-a))<straight_tolerance :
            return (p-i).mag(), [i.x, i.y]
        else :
            d1, d2 = (p-P0).mag(), (p-P2).mag()
            if d1<d2 :
                return (d1, [P0.x,P0.y])
            else :
                return (d2, [P2.x,P2.y])


def csp_to_arc_distance(sp1,sp2, arc1, arc2, tolerance = 0.0001 ): # arc = [start,end,center,alpha]
    n, i = 10, 0
    d, d1, dl = (0,(0,0)), (0,(0,0)), 0
    while i<1 or (abs(d1[0]-dl[0])>tolerance and i<4):
        i += 1
        dl = d1*1
        for j in range(n+1):
            t = float(j)/n
            p = csp_at_t(sp1,sp2,t)
            d = min(point_to_arc_distance(p,arc1), point_to_arc_distance(p,arc2))
            d1 = max(d1,d)
        n=n*2
    return d1[0]

def csp_point_inside_bound(sp1, sp2, p):
    bez = [sp1[1],sp1[2],sp2[0],sp2[1]]
    x,y = p
    c = 0
    for i in range(4):
        [x0,y0], [x1,y1] = bez[i-1], bez[i]
        if x0-x1!=0 and (y-y0)*(x1-x0)>=(x-x0)*(y1-y0) and x>min(x0,x1) and x<=max(x0,x1) :
            c +=1
    return c%2==0


def csp_bound_to_point_distance(sp1, sp2, p):
    if csp_point_inside_bound(sp1, sp2, p) :
        return 0.
    bez = csp_segment_to_bez(sp1,sp2)
    min_dist = 1e100
    for i in range(0,4):
        d = point_to_line_segment_distance_2(p, bez[i-1],bez[i])
        if d <= min_dist : min_dist = d
    return min_dist


def line_line_intersect(p1,p2,p3,p4) : # Return only true intersection.
    if (p1[0]==p2[0] and p1[1]==p2[1]) or (p3[0]==p4[0] and p3[1]==p4[1]) : return False
    x = (p2[0]-p1[0])*(p4[1]-p3[1]) - (p2[1]-p1[1])*(p4[0]-p3[0])
    if x==0 : # Lines are parallel
        if (p3[0]-p1[0])*(p2[1]-p1[1]) == (p3[1]-p1[1])*(p2[0]-p1[0]) :
            if p3[0]!=p4[0] :
                t11 = (p1[0]-p3[0])/(p4[0]-p3[0])
                t12 = (p2[0]-p3[0])/(p4[0]-p3[0])
                t21 = (p3[0]-p1[0])/(p2[0]-p1[0])
                t22 = (p4[0]-p1[0])/(p2[0]-p1[0])
            else:
                t11 = (p1[1]-p3[1])/(p4[1]-p3[1])
                t12 = (p2[1]-p3[1])/(p4[1]-p3[1])
                t21 = (p3[1]-p1[1])/(p2[1]-p1[1])
                t22 = (p4[1]-p1[1])/(p2[1]-p1[1])
            return ("Overlap" if (0<=t11<=1 or 0<=t12<=1) and (0<=t21<=1 or  0<=t22<=1) else False)
        else: return False
    else :
        return (
                    0<=((p4[0]-p3[0])*(p1[1]-p3[1]) - (p4[1]-p3[1])*(p1[0]-p3[0]))/x<=1 and
                    0<=((p2[0]-p1[0])*(p1[1]-p3[1]) - (p2[1]-p1[1])*(p1[0]-p3[0]))/x<=1 )


def line_line_intersection_points(p1,p2,p3,p4) : # Return only points [ (x,y) ]
    if (p1[0]==p2[0] and p1[1]==p2[1]) or (p3[0]==p4[0] and p3[1]==p4[1]) : return []
    x = (p2[0]-p1[0])*(p4[1]-p3[1]) - (p2[1]-p1[1])*(p4[0]-p3[0])
    if x==0 : # Lines are parallel
        if (p3[0]-p1[0])*(p2[1]-p1[1]) == (p3[1]-p1[1])*(p2[0]-p1[0]) :
            if p3[0]!=p4[0] :
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
            if (0<=t11<=1 or 0<=t12<=1) and (0<=t21<=1 or  0<=t22<=1) :
                if 0<=t11<=1 : res += [p1]
                if 0<=t12<=1 : res += [p2]
                if 0<=t21<=1 : res += [p3]
                if 0<=t22<=1 : res += [p4]
            return res
        else: return []
    else :
        t1 = ((p4[0]-p3[0])*(p1[1]-p3[1]) - (p4[1]-p3[1])*(p1[0]-p3[0]))/x
        t2 = ((p2[0]-p1[0])*(p1[1]-p3[1]) - (p2[1]-p1[1])*(p1[0]-p3[0]))/x
        if 0<=t1<=1 and 0<=t2<=1 : return [ [p1[0]*(1-t1)+p2[0]*t1, p1[1]*(1-t1)+p2[1]*t1] ]
        else : return []


def point_to_point_d2(a,b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2


def point_to_point_d(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)


def point_to_line_segment_distance_2(p1, p2,p3) :
    # p1 - point, p2,p3 - line segment
    #draw_pointer(p1)
    w0 = [p1[0]-p2[0], p1[1]-p2[1]]
    v = [p3[0]-p2[0], p3[1]-p2[1]]
    c1 = w0[0]*v[0] + w0[1]*v[1]
    if c1 <= 0 :
        return w0[0]*w0[0]+w0[1]*w0[1]
    c2 = v[0]*v[0] + v[1]*v[1]
    if c2 <= c1 :
        return  (p1[0]-p3[0])**2 + (p1[1]-p3[1])**2
    return (p1[0]- p2[0]-v[0]*c1/c2)**2 + (p1[1]- p2[1]-v[1]*c1/c2)


def line_to_line_distance_2(p1,p2,p3,p4):
    if line_line_intersect(p1,p2,p3,p4) : return 0
    return min(
            point_to_line_segment_distance_2(p1,p3,p4),
            point_to_line_segment_distance_2(p2,p3,p4),
            point_to_line_segment_distance_2(p3,p1,p2),
            point_to_line_segment_distance_2(p4,p1,p2))


def csp_seg_bound_to_csp_seg_bound_max_min_distance(sp1,sp2,sp3,sp4) :
    bez1 = csp_segment_to_bez(sp1,sp2)
    bez2 = csp_segment_to_bez(sp3,sp4)
    min_dist = 1e100
    max_dist = 0.
    for i in range(4) :
        if csp_point_inside_bound(sp1, sp2, bez2[i]) or csp_point_inside_bound(sp3, sp4, bez1[i]) :
            min_dist = 0.
            break
    for i in range(4) :
        for j in range(4) :
            d = line_to_line_distance_2(bez1[i-1],bez1[i],bez2[j-1],bez2[j])
            if d < min_dist : min_dist = d
            d = (bez2[j][0]-bez1[i][0])**2 + (bez2[j][1]-bez1[i][1])**2
            if max_dist < d  : max_dist = d
    return min_dist, max_dist


def csp_reverse(csp) :
    for i in range(len(csp)) :
        n = []
        for j in csp[i] :
            n = [  [j[2][:],j[1][:],j[0][:]]  ] + n
        csp[i] = n[:]
    return csp


def csp_normalized_slope(sp1,sp2,t) :
    ax,ay,bx,by,cx,cy,dx,dy=bezmisc.bezierparameterize((sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:]))
    if sp1[1]==sp2[1]==sp1[2]==sp2[0] : return [1.,0.]
    f1x = 3*ax*t*t+2*bx*t+cx
    f1y = 3*ay*t*t+2*by*t+cy
    if abs(f1x*f1x+f1y*f1y) > 1e-20 :
        l = math.sqrt(f1x*f1x+f1y*f1y)
        return [f1x/l, f1y/l]

    if t == 0 :
        f1x = sp2[0][0]-sp1[1][0]
        f1y = sp2[0][1]-sp1[1][1]
        if abs(f1x*f1x+f1y*f1y) > 1e-20 :
            l = math.sqrt(f1x*f1x+f1y*f1y)
            return [f1x/l, f1y/l]
        else :
            f1x = sp2[1][0]-sp1[1][0]
            f1y = sp2[1][1]-sp1[1][1]
            if f1x*f1x+f1y*f1y != 0 :
                l = math.sqrt(f1x*f1x+f1y*f1y)
                return [f1x/l, f1y/l]
    elif t == 1 :
        f1x = sp2[1][0]-sp1[2][0]
        f1y = sp2[1][1]-sp1[2][1]
        if abs(f1x*f1x+f1y*f1y) > 1e-20 :
            l = math.sqrt(f1x*f1x+f1y*f1y)
            return [f1x/l, f1y/l]
        else :
            f1x = sp2[1][0]-sp1[1][0]
            f1y = sp2[1][1]-sp1[1][1]
            if f1x*f1x+f1y*f1y != 0 :
                l = math.sqrt(f1x*f1x+f1y*f1y)
                return [f1x/l, f1y/l]
    else :
        return [1.,0.]


def csp_normalized_normal(sp1,sp2,t) :
    nx,ny = csp_normalized_slope(sp1,sp2,t)
    return [-ny, nx]


def csp_parameterize(sp1,sp2):
    return bezmisc.bezierparameterize(csp_segment_to_bez(sp1,sp2))


def csp_concat_subpaths(*s):

    def concat(s1,s2) :
        if s1 == [] : return s2
        if s2 == [] : return s1
        if (s1[-1][1][0]-s2[0][1][0])**2 + (s1[-1][1][1]-s2[0][1][1])**2 > 0.00001 :
            return s1[:-1]+[ [s1[-1][0],s1[-1][1],s1[-1][1]],  [s2[0][1],s2[0][1],s2[0][2]] ] + s2[1:]
        else :
            return s1[:-1]+[ [s1[-1][0],s2[0][1],s2[0][2]] ] + s2[1:]

    if len(s) == 0 : return []
    if len(s) ==1 : return s[0]
    result = s[0]
    for s1 in s[1:]:
        result = concat(result,s1)
    return result


def csp_draw(csp, color="#05f", group = None, style="fill:none;", width = .1, comment = "") :
    if csp!=[] and csp!=[[]] :
        if group == None : group = options.doc_root
        style += "stroke:"+color+";"+ "stroke-width:%0.4fpx;"%width
        args = {"d": cubicsuperpath.formatPath(csp), "style":style}
        if comment!="" : args["comment"] = str(comment)
        inkex.etree.SubElement( group, inkex.addNS('path','svg'), args )


def csp_subpaths_end_to_start_distance2(s1,s2):
    return (s1[-1][1][0]-s2[0][1][0])**2 + (s1[-1][1][1]-s2[0][1][1])**2


def csp_clip_by_line(csp,l1,l2) :
    result = []
    for i in range(len(csp)):
        s = csp[i]
        intersections = []
        for j in range(1,len(s)) :
            intersections += [  [j,int_] for int_ in csp_line_intersection(l1,l2,s[j-1],s[j])]
        splitted_s = csp_subpath_split_by_points(s, intersections)
        for s in splitted_s[:] :
            clip = False
            for p in csp_true_bounds([s]) :
                if (l1[1]-l2[1])*p[0] + (l2[0]-l1[0])*p[1] + (l1[0]*l2[1]-l2[0]*l1[1])<-0.01 :
                    clip = True
                    break
            if clip :
                splitted_s.remove(s)
        result += splitted_s
    return result


def csp_subpath_line_to(subpath, points) :
    # Appends subpath with line or polyline.
    if len(points)>0 :
        if len(subpath)>0:
            subpath[-1][2] = subpath[-1][1][:]
        if type(points[0]) == type([1,1]) :
            for p in points :
                subpath += [ [p[:],p[:],p[:]] ]
        else:
            subpath += [ [points,points,points] ]
    return subpath


def csp_join_subpaths(csp) :
    result = csp[:]
    done_smf = True
    joined_result = []
    while done_smf :
        done_smf = False
        while len(result)>0:
            s1 = result[-1][:]
            del(result[-1])
            j = 0
            joined_smf = False
            while j<len(joined_result) :
                if csp_subpaths_end_to_start_distance2(joined_result[j],s1) <0.000001 :
                    joined_result[j] = csp_concat_subpaths(joined_result[j],s1)
                    done_smf = True
                    joined_smf = True
                    break
                if csp_subpaths_end_to_start_distance2(s1,joined_result[j]) <0.000001 :
                    joined_result[j] = csp_concat_subpaths(s1,joined_result[j])
                    done_smf = True
                    joined_smf = True
                    break
                j += 1
            if not joined_smf : joined_result += [s1[:]]
        if done_smf :
            result = joined_result[:]
            joined_result = []
    return joined_result


def triangle_cross(a,b,c):
    return (a[0]-b[0])*(c[1]-b[1]) - (c[0]-b[0])*(a[1]-b[1])

################################################################################
###        Area Fill stuff
################################################################################

def point_inside_csp(p,csp, on_the_path = True) :
   
    x,y = p
    ray_intersections_count = 0
    for subpath in csp :
        
        for i in range(1, len(subpath)) :
            sp1, sp2 = subpath[i-1], subpath[i]
            ax,ay,bx,by,cx,cy,dx,dy = csp_parameterize(sp1,sp2)
            if  ax==0 and bx==0 and cx==0 and dx==x : 
                #we've got a special case here
                b = csp_true_bounds( [[sp1,sp2]])
                if  b[1][1]<=y<=b[3][1] :
                    # points is on the path 
                    return on_the_path
                else :
                    # we can skip this segment because it wont influence the answer.
                    pass    
            else: 
                for t in csp_line_intersection([x,y],[x,y+5],sp1,sp2) :
                    if t == 0 or t == 1 :
                        #we've got another special case here
                        x1,y1 = csp_at_t(sp1,sp2,t)
                        if y1==y : 
                            # the point is on the path 
                            return on_the_path
                        # if t == 0 we sould have considered this case previously. 
                        if t == 1 :
                            # we have to check the next segmant if it is on the same side of the ray
                            st_d = csp_normalized_slope(sp1,sp2,1)[0]
                            if st_d == 0 : st_d = csp_normalized_slope(sp1,sp2,0.99)[0]
                            
                            for j in range(1, len(subpath)+1):
                                if (i+j) % len(subpath) == 0  : continue # skip the closing segment 
                                sp11,sp22 = subpath[(i-1+j) % len(subpath)], subpath[(i+j) % len(subpath)]
                                ax1,ay1,bx1,by1,cx1,cy1,dx1,dy1 = csp_parameterize(sp1,sp2)
                                if  ax1==0 and bx1==0 and cx1==0 and dx1==x : continue # this segment parallel to the ray, so skip it 
                                en_d = csp_normalized_slope(sp11,sp22,0)[0]
                                if en_d == 0 : en_d = csp_normalized_slope(sp11,sp22,0.01)[0]
                                if st_d*en_d <=0 : 
                                    ray_intersections_count += 1
                                    break 
                    else :    
                        x1,y1 = csp_at_t(sp1,sp2,t)
                        if y1==y : 
                             # the point is on the path 
                            return on_the_path
                        else :
                            if y1>y and 3*ax*t**2 + 2*bx*t + cx !=0 : # if it's 0 the path only touches the ray
                                ray_intersections_count += 1     
    return ray_intersections_count%2 == 1

def csp_from_polyline(line) :
    return [ [ [point[:] for k in range(3) ] for point in subline ]  for subline in line ]

def csp_close_all_subpaths(csp, tolerance = 0.000001):
    for i in range(len(csp)):
        if point_to_point_d2(csp[i][0][1] , csp[i][-1][1])> tolerance**2 :
            csp[i][-1][2] = csp[i][-1][1][:] 
            csp[i] += [ [csp[i][0][1][:] for j in range(3)] ]
        else: 
            if csp[i][0][1] != csp[i][-1][1] : 
                csp[i][-1][1] = csp[i][0][1][:]
    return csp

################################################################################
###        Bezier additional functions
################################################################################

def bez_bounds_intersect(bez1, bez2) :
    return bounds_intersect(bez_bound(bez2), bez_bound(bez1))


def bez_bound(bez) :
    return [
                min(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
                min(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
                max(bez[0][0], bez[1][0], bez[2][0], bez[3][0]),
                max(bez[0][1], bez[1][1], bez[2][1], bez[3][1]),
            ]


def bounds_intersect(a, b) :
    return not ( (a[0]>b[2]) or (b[0]>a[2]) or (a[1]>b[3]) or (b[1]>a[3]) )


def tpoint((x1,y1),(x2,y2),t):
    return [x1+t*(x2-x1),y1+t*(y2-y1)]


def bez_to_csp_segment(bez) :
    return [bez[0],bez[0],bez[1]], [bez[2],bez[3],bez[3]]


def bez_split(a,t=0.5) :
    a1 = tpoint(a[0],a[1],t)
    at = tpoint(a[1],a[2],t)
    b2 = tpoint(a[2],a[3],t)
    a2 = tpoint(a1,at,t)
    b1 = tpoint(b2,at,t)
    a3 = tpoint(a2,b1,t)
    return [a[0],a1,a2,a3], [a3,b1,b2,a[3]]


def bez_at_t(bez,t) :
    return csp_at_t([bez[0],bez[0],bez[1]],[bez[2],bez[3],bez[3]],t)


def bez_to_point_distance(bez,p,needed_dist=[0.,1e100]):
    # returns [d^2,t]
    return csp_seg_to_point_distance(bez_to_csp_segment(bez),p,needed_dist)


def bez_normalized_slope(bez,t):
    return csp_normalized_slope([bez[0],bez[0],bez[1]], [bez[2],bez[3],bez[3]],t)

################################################################################
###    Some vector functions
################################################################################

def normalize((x,y)) :
    l = math.sqrt(x**2+y**2)
    if l == 0 : return [0.,0.]
    else :         return [x/l, y/l]


def cross(a,b) :
    return a[1] * b[0] - a[0] * b[1]


def dot(a,b) :
    return a[0] * b[0] + a[1] * b[1]


def rotate_ccw(d) :
    return [-d[1],d[0]]


def vectors_ccw(a,b):
    return a[0]*b[1]-b[0]*a[1] < 0


def vector_from_to_length(a,b):
    return math.sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))

################################################################################
###    Common functions
################################################################################

def matrix_mul(a,b) :
    return [ [ sum([a[i][k]*b[k][j] for k in range(len(a[0])) ])   for j in range(len(b[0]))]   for i in range(len(a))]
    try :
        return [ [ sum([a[i][k]*b[k][j] for k in range(len(a[0])) ])   for j in range(len(b[0]))]   for i in range(len(a))]
    except :
        return None


def transpose(a) :
    try :
        return [ [ a[i][j] for i in range(len(a)) ] for j in range(len(a[0])) ]
    except :
        return None


def det_3x3(a):
    return  float(
        a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[1][0]*a[2][1]*a[0][2]
        - a[0][2]*a[1][1]*a[2][0] - a[0][0]*a[2][1]*a[1][2] - a[0][1]*a[2][2]*a[1][0]
        )


def inv_3x3(a): # invert matrix 3x3
    det = det_3x3(a)
    if det==0: return None
    return    [
        [  (a[1][1]*a[2][2] - a[2][1]*a[1][2])/det,  -(a[0][1]*a[2][2] - a[2][1]*a[0][2])/det,  (a[0][1]*a[1][2] - a[1][1]*a[0][2])/det ],
        [ -(a[1][0]*a[2][2] - a[2][0]*a[1][2])/det,   (a[0][0]*a[2][2] - a[2][0]*a[0][2])/det, -(a[0][0]*a[1][2] - a[1][0]*a[0][2])/det ],
        [  (a[1][0]*a[2][1] - a[2][0]*a[1][1])/det,  -(a[0][0]*a[2][1] - a[2][0]*a[0][1])/det,  (a[0][0]*a[1][1] - a[1][0]*a[0][1])/det ]
    ]


def inv_2x2(a): # invert matrix 2x2
    det = a[0][0]*a[1][1] - a[1][0]*a[0][1]
    if det==0: return None
    return [
            [a[1][1]/det, -a[0][1]/det],
            [-a[1][0]/det, a[0][0]/det]
            ]


def small(a) :
    global small_tolerance
    return abs(a)<small_tolerance


def atan2(*arg):
    if len(arg)==1 and ( type(arg[0]) == type([0.,0.]) or type(arg[0])==type((0.,0.)) ) :
        return (math.pi/2 - math.atan2(arg[0][0], arg[0][1]) ) % math.pi2
    elif len(arg)==2 :

        return (math.pi/2 - math.atan2(arg[0],arg[1]) ) % math.pi2
    else :
        raise ValueError, "Bad argumets for atan! (%s)" % arg


def draw_text(text,x,y,style = None, font_size = 20) :
    if style == None :
        style = "font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;fill:#000000;fill-opacity:1;stroke:none;"
    style += "font-size:%fpx;"%font_size
    t = inkex.etree.SubElement(    options.doc_root, inkex.addNS('text','svg'), {
                            'x':    str(x),
                            inkex.addNS("space","xml"):"preserve",
                            'y':    str(y)
                        })
    text = str(text).split("\n")
    for s in text :
        span = inkex.etree.SubElement( t, inkex.addNS('tspan','svg'),
                        {
                            'x':    str(x),
                            'y':    str(+y),
                            inkex.addNS("role","sodipodi"):"line",
                        })
        y += font_size
        span.text = s


def draw_pointer(x,color = "#f00", figure = "cross", comment = "", width = .1) :
    if figure ==  "line" :
        s = ""
        for i in range(1,len(x)/2) :
            s+= " %s, %s " %(x[i*2],x[i*2+1])
        inkex.etree.SubElement( options.doc_root, inkex.addNS('path','svg'), {"d": "M %s,%s L %s"%(x[0],x[1],s), "style":"fill:none;stroke:%s;stroke-width:%f;"%(color,width),"comment":str(comment)} )
    else :
        inkex.etree.SubElement( options.doc_root, inkex.addNS('path','svg'), {"d": "m %s,%s l 10,10 -20,-20 10,10 -10,10, 20,-20"%(x[0],x[1]), "style":"fill:none;stroke:%s;stroke-width:%f;"%(color,width),"comment":str(comment)} )


def straight_segments_intersection(a,b, true_intersection = True) : # (True intersection means check ta and tb are in [0,1])
    ax,bx,cx,dx, ay,by,cy,dy = a[0][0],a[1][0],b[0][0],b[1][0], a[0][1],a[1][1],b[0][1],b[1][1]
    if (ax==bx and ay==by) or (cx==dx and cy==dy) : return False, 0, 0
    if (bx-ax)*(dy-cy)-(by-ay)*(dx-cx)==0 :    # Lines are parallel
        ta = (ax-cx)/(dx-cx) if cx!=dx else (ay-cy)/(dy-cy)
        tb = (bx-cx)/(dx-cx) if cx!=dx else (by-cy)/(dy-cy)
        tc = (cx-ax)/(bx-ax) if ax!=bx else (cy-ay)/(by-ay)
        td = (dx-ax)/(bx-ax) if ax!=bx else (dy-ay)/(by-ay)
        return ("Overlap" if 0<=ta<=1 or 0<=tb<=1 or  0<=tc<=1 or  0<=td<=1 or not true_intersection else False), (ta,tb), (tc,td)
    else :
        ta = ( (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy) ) / ( (bx-ax)*(dy-cy)-(by-ay)*(dx-cx) )
        tb = ( ax-cx+ta*(bx-ax) ) / (dx-cx) if dx!=cx else ( ay-cy+ta*(by-ay) ) / (dy-cy)
        return (0<=ta<=1 and 0<=tb<=1 or not true_intersection), ta, tb



def isnan(x): return type(x) is float and x != x

def isinf(x): inf = 1e5000; return x == inf or x == -inf

def between(c,x,y):
    return x-straight_tolerance<=c<=y+straight_tolerance or y-straight_tolerance<=c<=x+straight_tolerance


def cubic_solver(a,b,c,d):
    if a!=0:
    #    Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
        a,b,c = (b/a, c/a, d/a)
        m = 2*a**3 - 9*a*b + 27*c
        k = a**2 - 3*b
        n = m**2 - 4*k**3
        w1 = -.5 + .5*cmath.sqrt(3)*1j
        w2 = -.5 - .5*cmath.sqrt(3)*1j
        if n>=0 :
            t = m+math.sqrt(n)
            m1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
            t = m-math.sqrt(n)
            n1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
        else :
            m1 = pow(complex((m+cmath.sqrt(n))/2),1./3)
            n1 = pow(complex((m-cmath.sqrt(n))/2),1./3)
        x1 = -1./3 * (a + m1 + n1)
        x2 = -1./3 * (a + w1*m1 + w2*n1)
        x3 = -1./3 * (a + w2*m1 + w1*n1)
        return [x1,x2,x3]
    elif b!=0:
        det = c**2-4*b*d
        if det>0 :
            return [(-c+math.sqrt(det))/(2*b),(-c-math.sqrt(det))/(2*b)]
        elif d == 0 :
            return [-c/(b*b)]
        else :
            return [(-c+cmath.sqrt(det))/(2*b),(-c-cmath.sqrt(det))/(2*b)]
    elif c!=0 :
        return [-d/c]
    else : return []


################################################################################
###        print_ prints any arguments into specified log file
################################################################################

def print_(*arg):
    f = open(options.log_filename,"a")
    for s in arg :
        s = str(unicode(s).encode('unicode_escape'))+" "
        f.write( s )
    f.write("\n")
    f.close()


################################################################################
###        Point (x,y) operations
################################################################################
class P:
    def __init__(self, x, y=None):
        if not y==None:
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
        if h: return self / h
        else: return P(0,0)
    def dot(self, other): return self.x * other.x + self.y * other.y
    def rot(self, theta):
        c = math.cos(theta)
        s = math.sin(theta)
        return P(self.x * c - self.y * s,  self.x * s + self.y * c)
    def angle(self): return math.atan2(self.y, self.x)
    def __repr__(self): return '%f,%f' % (self.x, self.y)
    def pr(self): return "%.3f,%.3f" % (self.x, self.y)
    def to_list(self): return [self.x, self.y]
    def ccw(self): return P(-self.y,self.x)
    def l2(self): return self.x*self.x + self.y*self.y


################################################################################
###
###        Biarc function
###
###        Calculates biarc approximation of cubic super path segment
###        splits segment if needed or approximates it with straight line
###
################################################################################
def biarc(sp1, sp2, z1, z2, depth=0):
    def biarc_split(sp1,sp2, z1, z2, depth):
        if depth<options.biarc_max_split_depth:
            sp1,sp2,sp3 = csp_split(sp1,sp2)
            l1, l2 = cspseglength(sp1,sp2), cspseglength(sp2,sp3)
            if l1+l2 == 0 : zm = z1
            else : zm = z1+(z2-z1)*l1/(l1+l2)
            return biarc(sp1,sp2,z1,zm,depth+1)+biarc(sp2,sp3,zm,z2,depth+1)
        else: return [ [sp1[1],'line', 0, 0, sp2[1], [z1,z2]] ]

    P0, P4 = P(sp1[1]), P(sp2[1])
    TS, TE, v = (P(sp1[2])-P0), -(P(sp2[0])-P4), P0 - P4
    tsa, tea, va = TS.angle(), TE.angle(), v.angle()
    if TE.mag()<straight_distance_tolerance and TS.mag()<straight_distance_tolerance:
        # Both tangents are zerro - line straight
        return [ [sp1[1],'line', 0, 0, sp2[1], [z1,z2]] ]
    if TE.mag() < straight_distance_tolerance:
        TE = -(TS+v).unit()
        r = TS.mag()/v.mag()*2
    elif TS.mag() < straight_distance_tolerance:
        TS = -(TE+v).unit()
        r = 1/( TE.mag()/v.mag()*2 )
    else:
        r=TS.mag()/TE.mag()
    TS, TE = TS.unit(), TE.unit()
    tang_are_parallel = ((tsa-tea)%math.pi<straight_tolerance or math.pi-(tsa-tea)%math.pi<straight_tolerance )
    if ( tang_are_parallel  and
                ((v.mag()<straight_distance_tolerance or TE.mag()<straight_distance_tolerance or TS.mag()<straight_distance_tolerance) or
                    1-abs(TS*v/(TS.mag()*v.mag()))<straight_tolerance)    ):
                # Both tangents are parallel and start and end are the same - line straight
                # or one of tangents still smaller then tollerance

                # Both tangents and v are parallel - line straight
        return [ [sp1[1],'line', 0, 0, sp2[1], [z1,z2]] ]

    c,b,a = v*v, 2*v*(r*TS+TE), 2*r*(TS*TE-1)
    if v.mag()==0:
        return biarc_split(sp1, sp2, z1, z2, depth)
    asmall, bsmall, csmall = abs(a)<10**-10,abs(b)<10**-10,abs(c)<10**-10
    if         asmall and b!=0:    beta = -c/b
    elif     csmall and a!=0:    beta = -b/a
    elif not asmall:
        discr = b*b-4*a*c
        if discr < 0:    raise ValueError, (a,b,c,discr)
        disq = discr**.5
        beta1 = (-b - disq) / 2 / a
        beta2 = (-b + disq) / 2 / a
        if beta1*beta2 > 0 :    raise ValueError, (a,b,c,disq,beta1,beta2)
        beta = max(beta1, beta2)
    elif    asmall and bsmall:
        return biarc_split(sp1, sp2, z1, z2, depth)
    alpha = beta * r
    ab = alpha + beta
    P1 = P0 + alpha * TS
    P3 = P4 - beta * TE
    P2 = (beta / ab)  * P1 + (alpha / ab) * P3


    def calculate_arc_params(P0,P1,P2):
        D = (P0+P2)/2
        if (D-P1).mag()==0: return None, None
        R = D - ( (D-P0).mag()**2/(D-P1).mag() )*(P1-D).unit()
        p0a, p1a, p2a = (P0-R).angle()%(2*math.pi), (P1-R).angle()%(2*math.pi), (P2-R).angle()%(2*math.pi)
        alpha =  (p2a - p0a) % (2*math.pi)
        if (p0a<p2a and  (p1a<p0a or p2a<p1a))    or    (p2a<p1a<p0a) :
            alpha = -2*math.pi+alpha
        if abs(R.x)>1000000 or abs(R.y)>1000000  or (R-P0).mag<.1 :
            return None, None
        else :
            return  R, alpha
    R1,a1 = calculate_arc_params(P0,P1,P2)
    R2,a2 = calculate_arc_params(P2,P3,P4)
    if R1==None or R2==None or (R1-P0).mag()<straight_tolerance or (R2-P2).mag()<straight_tolerance    : return [ [sp1[1],'line', 0, 0, sp2[1], [z1,z2]] ]

    d = csp_to_arc_distance(sp1,sp2, [P0,P2,R1,a1],[P2,P4,R2,a2])
    if d > 1 and depth<options.biarc_max_split_depth     : return biarc_split(sp1, sp2, z1, z2, depth)
    else:
        if R2.mag()*a2 == 0 : zm = z2
        else : zm  = z1 + (z2-z1)*(abs(R1.mag()*a1))/(abs(R2.mag()*a2)+abs(R1.mag()*a1))
        return [    [ sp1[1], 'arc', [R1.x,R1.y], a1, [P2.x,P2.y], [z1,zm] ], [ [P2.x,P2.y], 'arc', [R2.x,R2.y], a2, [P4.x,P4.y], [zm,z2] ]        ]


def biarc_curve_segment_length(seg):
    if seg[1] == "arc" :
        return math.sqrt((seg[0][0]-seg[2][0])**2+(seg[0][1]-seg[2][1])**2)*seg[3]
    elif seg[1] == "line" :
        return math.sqrt((seg[0][0]-seg[4][0])**2+(seg[0][1]-seg[4][1])**2)
    else:
        return 0


def biarc_curve_clip_at_l(curve, l, clip_type = "strict") :
    # get first subcurve and ceck it's length
    subcurve, subcurve_l, moved = [], 0, False
    for seg in curve:
        if seg[1] == "move" and moved or seg[1] == "end" :
            break
        if seg[1] == "move" : moved = True
        subcurve_l += biarc_curve_segment_length(seg)
        if seg[1] == "arc" or seg[1] == "line" :
            subcurve += [seg]

    if subcurve_l < l and clip_type == "strict" : return []
    lc = 0
    if (subcurve[-1][4][0]-subcurve[0][0][0])**2 + (subcurve[-1][4][1]-subcurve[0][0][1])**2 < 10**-7 : subcurve_closed = True
    i = 0
    reverse = False
    while lc<l :
        seg = subcurve[i]
        if reverse :
            if seg[1] == "line" :
                seg = [seg[4], "line", 0 , 0, seg[0], seg[5]] # Hmmm... Do we have to swap seg[5][0] and seg[5][1] (zstart and zend) or not?
            elif seg[1] == "arc" :
                seg = [seg[4], "arc", seg[2] , -seg[3], seg[0], seg[5]] # Hmmm... Do we have to swap seg[5][0] and seg[5][1] (zstart and zend) or not?
        ls = biarc_curve_segment_length(seg)
        if ls != 0 :
            if l-lc>ls :
                res += [seg]
            else :
                if seg[1] == "arc" :
                    r  = math.sqrt((seg[0][0]-seg[2][0])**2+(seg[0][1]-seg[2][1])**2)
                    x,y = seg[0][0]-seg[2][0], seg[0][1]-seg[2][1]
                    a = seg[3]/ls*(l-lc)
                    x,y = x*math.cos(a) - y*math.sin(a),  x*math.sin(a) + y*math.cos(a)
                    x,y = x+seg[2][0], y+seg[2][1]
                    res += [[ seg[0], "arc",  seg[2], a, [x,y], [seg[5][0],seg[5][1]/ls*(l-lc)]  ]]
                if seg[1] == "line" :
                    res += [[ seg[0], "line",  0, 0, [(seg[4][0]-seg[0][0])/ls*(l-lc),(seg[4][1]-seg[0][1])/ls*(l-lc)], [seg[5][0],seg[5][1]/ls*(l-lc)]  ]]
        i += 1
        if i >= len(subcurve) and not subcurve_closed:
            reverse = not reverse
        i = i%len(subcurve)
    return res

################################################################################
###        Polygon class
################################################################################
class Polygon:
    def __init__(self, polygon=None):
        self.polygon = [] if polygon==None else polygon[:]


    def move(self, x, y) :
        for i in range(len(self.polygon)) :
            for j in range(len(self.polygon[i])) :
                self.polygon[i][j][0] += x
                self.polygon[i][j][1] += y


    def bounds(self) :
        minx,miny,maxx,maxy = 1e400, 1e400, -1e400, -1e400
        for poly in self.polygon :
            for p in poly :
                if minx > p[0] : minx = p[0]
                if miny > p[1] : miny = p[1]
                if maxx < p[0] : maxx = p[0]
                if maxy < p[1] : maxy = p[1]
        return minx*1,miny*1,maxx*1,maxy*1


    def width(self):
        b = self.bounds()
        return b[2]-b[0]


    def rotate_(self,sin,cos) :
        for i in range(len(self.polygon)) :
            for j in range(len(self.polygon[i])) :
                x,y = self.polygon[i][j][0], self.polygon[i][j][1]
                self.polygon[i][j][0] = x*cos - y*sin
                self.polygon[i][j][1] = x*sin + y*cos


    def rotate(self, a):
        cos, sin = math.cos(a), math.sin(a)
        self.rotate_(sin,cos)


    def drop_into_direction(self, direction, surface) :
        # Polygon is a list of simple polygons
        # Surface is a polygon + line y = 0
        # Direction is [dx,dy]
        if len(self.polygon) == 0 or len(self.polygon[0])==0 : return
        if direction[0]**2 + direction[1]**2 <1e-10 : return
        direction = normalize(direction)
        sin,cos = direction[0], -direction[1]
        self.rotate_(-sin,cos)
        surface.rotate_(-sin,cos)
        self.drop_down(surface, zerro_plane = False)
        self.rotate_(sin,cos)
        surface.rotate_(sin,cos)


    def centroid(self):
        centroids = []
        sa = 0
        for poly in self.polygon:
            cx,cy,a = 0,0,0
            for i in range(len(poly)):
                [x1,y1],[x2,y2] = poly[i-1],poly[i]
                cx += (x1+x2)*(x1*y2-x2*y1)
                cy += (y1+y2)*(x1*y2-x2*y1)
                a  += (x1*y2-x2*y1)
            a *= 3.
            if abs(a)>0 :
                cx /= a
                cy /= a
                sa += abs(a)
                centroids += [ [cx,cy,a] ]
        if sa == 0 : return    [0.,0.]
        cx,cy = 0.,0.
        for c in centroids :
            cx += c[0]*c[2]
            cy += c[1]*c[2]
        cx /= sa
        cy /= sa
        return [cx,cy]


    def drop_down(self, surface, zerro_plane = True) :
        # Polygon is a list of simple polygons
        # Surface is a polygon + line y = 0
        # Down means min y (0,-1)
        if len(self.polygon) == 0 or len(self.polygon[0])==0 : return
        # Get surface top point
        top = surface.bounds()[3]
        if zerro_plane : top = max(0, top)
        # Get polygon bottom point
        bottom = self.bounds()[1]
        self.move(0, top - bottom + 10)
        # Now get shortest distance from surface to polygon in positive x=0 direction
        # Such distance = min(distance(vertex, edge)...)  where edge from surface and
        # vertex from polygon and vice versa...
        dist = 1e300
        for poly in surface.polygon :
            for i in range(len(poly)) :
                for poly1 in self.polygon :
                    for i1 in range(len(poly1)) :
                        st,end = poly[i-1], poly[i]
                        vertex = poly1[i1]
                        if st[0]<=vertex[0]<= end[0] or end[0]<=vertex[0]<=st[0] :
                            if st[0]==end[0] : d = min(vertex[1]-st[1],vertex[1]-end[1])
                            else : d = vertex[1] - st[1] - (end[1]-st[1])*(vertex[0]-st[0])/(end[0]-st[0])
                            if dist > d  : dist = d
                        # and vice versa just change the sign because vertex now under the edge
                        st,end = poly1[i1-1], poly1[i1]
                        vertex = poly[i]
                        if st[0]<=vertex[0]<=end[0] or end[0]<=vertex[0]<=st[0] :
                            if st[0]==end[0] : d = min(- vertex[1]+st[1],-vertex[1]+end[1])
                            else : d =  - vertex[1] + st[1] + (end[1]-st[1])*(vertex[0]-st[0])/(end[0]-st[0])
                            if dist > d  : dist = d

        if zerro_plane and dist > 10 + top : dist = 10 + top
        #print_(dist, top, bottom)
        #self.draw()
        self.move(0, -dist)


    def draw(self,color="#075",width=.1) :
        for poly in self.polygon :
            csp_draw( [csp_subpath_line_to([],poly+[poly[0]])], color=color,width=width )


    def add(self, add) :
        if type(add) == type([]) :
            self.polygon += add[:]
        else :
            self.polygon += add.polygon[:]


    def point_inside(self,p) :
        inside = False
        for poly in self.polygon :
            for i in range(len(poly)):
                st,end = poly[i-1], poly[i]
                if p==st or p==end : return True # point is a vertex = point is on the edge
                if st[0]>end[0] : st, end = end, st # This will be needed to check that edge if open only at rigth end
                c = (p[1]-st[1])*(end[0]-st[0])-(end[1]-st[1])*(p[0]-st[0])
                #print_(c)
                if st[0]<=p[0]<end[0] :
                    if c<0 :
                        inside = not inside
                    elif c == 0 : return True # point is on the edge
                elif st[0]==end[0]==p[0] and (st[1]<=p[1]<=end[1] or end[1]<=p[1]<=st[1]) : # point is on the edge
                    return True
        return inside


    def hull(self) :
        # Add vertices at all self intersection points.
        hull = []
        for i1 in range(len(self.polygon)):
            poly1 = self.polygon[i1]
            poly_ = []
            for j1 in range(len(poly1)):
                s, e = poly1[j1-1],poly1[j1]
                poly_ += [s]

                # Check self intersections
                for j2 in range(j1+1,len(poly1)):
                    s1, e1 = poly1[j2-1],poly1[j2]
                    int_ = line_line_intersection_points(s,e,s1,e1)
                    for p in int_ :
                        if point_to_point_d2(p,s)>0.000001 and point_to_point_d2(p,e)>0.000001 :
                            poly_ += [p]
                # Check self intersections with other polys
                for i2 in range(len(self.polygon)):
                    if i1==i2 : continue
                    poly2 = self.polygon[i2]
                    for j2 in range(len(poly2)):
                        s1, e1 = poly2[j2-1],poly2[j2]
                        int_ = line_line_intersection_points(s,e,s1,e1)
                        for p in int_ :
                            if point_to_point_d2(p,s)>0.000001 and point_to_point_d2(p,e)>0.000001 :
                                poly_ += [p]
            hull += [poly_]
        # Create the dictionary containing all edges in both directions
        edges = {}
        for poly in self.polygon :
            for i in range(len(poly)):
                s,e = tuple(poly[i-1]), tuple(poly[i])
                if (point_to_point_d2(e,s)<0.000001) : continue
                break_s, break_e = False, False
                for p in edges :
                    if point_to_point_d2(p,s)<0.000001 :
                        break_s = True
                        s = p
                    if point_to_point_d2(p,e)<0.000001 :
                        break_e = True
                        e = p
                    if break_s and break_e : break
                l = point_to_point_d(s,e)
                if not break_s and not break_e :
                    edges[s] = [ [s,e,l] ]
                    edges[e] = [ [e,s,l] ]
                    #draw_pointer(s+e,"red","line")
                    #draw_pointer(s+e,"red","line")
                else :
                    if e in edges :
                        for edge in edges[e] :
                            if point_to_point_d2(edge[1],s)<0.000001 :
                                break
                        if point_to_point_d2(edge[1],s)>0.000001 :
                            edges[e] += [ [e,s,l] ]
                            #draw_pointer(s+e,"red","line")

                    else :
                        edges[e] = [ [e,s,l] ]
                        #draw_pointer(s+e,"green","line")
                    if s in edges :
                        for edge in edges[s] :
                            if  point_to_point_d2(edge[1],e)<0.000001 :
                                break
                        if point_to_point_d2(edge[1],e)>0.000001 :
                            edges[s] += [ [s,e, l] ]
                            #draw_pointer(s+e,"red","line")
                    else :
                        edges[s] = [ [s,e,l] ]
                        #draw_pointer(s+e,"green","line")


        def angle_quadrant(sin,cos):
            # quadrants are (0,pi/2], (pi/2,pi], (pi,3*pi/2], (3*pi/2, 2*pi], i.e. 0 is in the 4-th quadrant
            if sin>0 and cos>=0 : return 1
            if sin>=0 and cos<0 : return 2
            if sin<0 and cos<=0 : return 3
            if sin<=0 and cos>0 : return 4


        def angle_is_less(sin,cos,sin1,cos1):
            # 0 = 2*pi is the largest angle
            if [sin1, cos1] == [0,1] : return True
            if [sin, cos] == [0,1] : return False
            if angle_quadrant(sin,cos)>angle_quadrant(sin1,cos1) :
                return False
            if angle_quadrant(sin,cos)<angle_quadrant(sin1,cos1) :
                return True
            if sin>=0 and cos>0 : return sin<sin1
            if sin>0 and cos<=0 : return sin>sin1
            if sin<=0 and cos<0 : return sin>sin1
            if sin<0 and cos>=0 : return sin<sin1


        def get_closes_edge_by_angle(edges, last):
            # Last edge is normalized vector of the last edge.
            min_angle = [0,1]
            next = last
            last_edge = [(last[0][0]-last[1][0])/last[2], (last[0][1]-last[1][1])/last[2]]
            for p in edges:
                #draw_pointer(list(p[0])+[p[0][0]+last_edge[0]*40,p[0][1]+last_edge[1]*40], "Red", "line", width=1)
                #print_("len(edges)=",len(edges))
                cur = [(p[1][0]-p[0][0])/p[2],(p[1][1]-p[0][1])/p[2]]
                cos, sin = dot(cur,last_edge),  cross(cur,last_edge)
                #draw_pointer(list(p[0])+[p[0][0]+cur[0]*40,p[0][1]+cur[1]*40], "Orange", "line", width=1, comment = [sin,cos])
                #print_("cos, sin=",cos,sin)
                #print_("min_angle_before=",min_angle)

                if     angle_is_less(sin,cos,min_angle[0],min_angle[1]) :
                    min_angle = [sin,cos]
                    next = p
                #print_("min_angle=",min_angle)

            return next

        # Join edges together into new polygon cutting the vertexes inside new polygon
        self.polygon = []
        len_edges = sum([len(edges[p]) for p in edges])
        loops = 0

        while len(edges)>0 :
            poly = []
            if loops > len_edges  : raise ValueError, "Hull error"
            loops+=1
            # Find left most vertex.
            start = (1e100,1)
            for edge in edges :
                start = min(start, min(edges[edge]))
            last = [(start[0][0]-1,start[0][1]),start[0],1]
            first_run = True
            loops1 = 0
            while (last[1]!=start[0] or first_run) :
                first_run = False
                if loops1 > len_edges  : raise ValueError, "Hull error"
                loops1 += 1
                next = get_closes_edge_by_angle(edges[last[1]],last)
                #draw_pointer(next[0]+next[1],"Green","line", comment=i, width= 1)
                #print_(next[0],"-",next[1])

                last = next
                poly += [ list(last[0]) ]
            self.polygon += [ poly ]
            # Remove all edges that are intersects new poly (any vertex inside new poly)
            poly_ = Polygon([poly])
            for p in edges.keys()[:] :
                if poly_.point_inside(list(p)) : del edges[p]
        self.draw(color="Green", width=.1)




################################################################################
###
###        Gcodetools class
###
################################################################################

class laser_gcode(inkex.Effect):

    def export_gcode(self,gcode):
        gcode_pass = gcode
        for x in range(1,self.options.passes):
            gcode += "G91\n" + "\nG90\n" + gcode_pass
        f = open(self.options.directory+self.options.file, "w")
        f.write(self.header + "\n" + gcode + self.footer)
        f.close()

    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("-d", "--directory",                       action="store", type="string",          dest="directory",                           default="/home/chris/GCodeExport",      help="Output directory")
        self.OptionParser.add_option("-f", "--filename",                        action="store", type="string",          dest="file",                                default="output.ngc",                   help="File name")
        self.OptionParser.add_option("",   "--add-numeric-suffix-to-filename",  action="store", type="inkbool",         dest="add_numeric_suffix_to_filename",      default=False,                          help="Add numeric suffix to file name")
        self.OptionParser.add_option("",   "--laser-command",                   action="store", type="string",          dest="laser_command",                       default="S100",                         help="Laser gcode command")
        self.OptionParser.add_option("",   "--laser-off-command",               action="store", type="string",          dest="laser_off_command",                   default="S1",                           help="Laser gcode end command")
        self.OptionParser.add_option("",   "--laser-beam-with",                 action="store", type="float",           dest="laser_beam_with",                     default="0.3",                          help="Laser speed (mm/min)")
        self.OptionParser.add_option("",   "--laser-speed",                     action="store", type="int",             dest="laser_speed",                         default="600",                          help="Laser speed (mm/min)")
        self.OptionParser.add_option("",   "--passes",                          action="store", type="int",             dest="passes",                              default="1",                            help="Quantity of passes")
        self.OptionParser.add_option("",   "--power-delay",                     action="store", type="string",          dest="power_delay",                         default="0",                            help="Laser power-on delay (ms)")
        self.OptionParser.add_option("",   "--suppress-all-messages",           action="store", type="inkbool",         dest="suppress_all_messages",               default=True,                           help="Hide messages during g-code generation")
        self.OptionParser.add_option("",   "--create-log",                      action="store", type="inkbool",         dest="log_create_log",                      default=True,                           help="Create log files")
        self.OptionParser.add_option("",   "--log-filename",                    action="store", type="string",          dest="log_filename",                        default='/home/chris/GCodeExport/log.txt',                             help="Create log files")
        self.OptionParser.add_option("",   "--engraving-draw-calculation-paths",action="store", type="inkbool",         dest="engraving_draw_calculation_paths",    default=True,                           help="Draw additional graphics to debug engraving path")
        self.OptionParser.add_option("",   "--unit",                            action="store", type="string",          dest="unit",                                default="G21 (All units in mm)",        help="Units either mm or inches")
        self.OptionParser.add_option("",   "--active-tab",                      action="store", type="string",          dest="active_tab",                          default="",                             help="Defines which tab is active")
        self.OptionParser.add_option("",   "--biarc-max-split-depth",           action="store", type="int",             dest="biarc_max_split_depth",               default="10",                           help="Defines maximum depth of splitting while approximating using biarcs.")
        self.OptionParser.add_option("",   "--area-fill-angle",                 action="store", type="float",           dest="area_fill_angle",                     default="0",                            help="Fill area with lines heading this angle")
        self.OptionParser.add_option("",   "--area-fill-shift",                 action="store", type="float",           dest="area_fill_shift",                     default="0",                            help="Shift the lines by tool d * shift")
        self.OptionParser.add_option("",   "--area-fill-method",                action="store", type="string",          dest="area_fill_method",                    default="zig-zag",                      help="Filling method either zig-zag or spiral")
        self.OptionParser.add_option("",   "--engraving-newton-iterations",     action="store", type="int",             dest="engraving_newton_iterations",         default="10",                           help="Number of sample points used to calculate distance")    
        self.OptionParser.add_option("",   "--add-contours",                    action="store", type="inkbool",         dest="add_contours",                        default=True,                           help="Add contour to Gcode paths")
        self.OptionParser.add_option("",   "--add-infill",                      action="store", type="inkbool",         dest="add_infill",                          default=True,                           help="Add infill to Gcode paths")
        
    def parse_curve(self, p, layer, w = None, f = None):
              
        c = []
        if len(p)==0 :
            return []
        p = self.transform_csp(p, layer)

        ### Sort to reduce Rapid distance
        k = range(1,len(p))
        keys = [0]
        while len(k)>0:
            end = p[keys[-1]][-1][1]
            dist = None
            for i in range(len(k)):
                start = p[k[i]][0][1]
                dist = max(   ( -( ( end[0]-start[0])**2+(end[1]-start[1])**2 ) ,i)    ,   dist )
            keys += [k[dist[1]]]
            del k[dist[1]]
        for k in keys:
            subpath = p[k]
            c += [ [    [subpath[0][1][0],subpath[0][1][1]]   , 'move', 0, 0] ]
            for i in range(1,len(subpath)):
                sp1 = [  [subpath[i-1][j][0], subpath[i-1][j][1]] for j in range(3)]
                sp2 = [  [subpath[i  ][j][0], subpath[i  ][j][1]] for j in range(3)]
                c += biarc(sp1,sp2,0,0) if w==None else biarc(sp1,sp2,-f(w[k][i-1]),-f(w[k][i]))
            c += [ [ [subpath[-1][1][0],subpath[-1][1][1]]  ,'end',0,0] ]
            print_("Curve: " + str(c) + "/n")
        return c


    def draw_curve(self, curve, layer, group=None, style=styles["biarc_style"]):

        self.get_defs()
        # Add marker to defs if it doesnot exists

        for i in [0,1]:
            style['biarc%s_r'%i] = simplestyle.parseStyle(style['biarc%s'%i])
            style['biarc%s_r'%i]["marker-start"] = "none"
            del(style['biarc%s_r'%i]["marker-end"])
            style['biarc%s_r'%i] = simplestyle.formatStyle(style['biarc%s_r'%i])

        if group==None:
            group = inkex.etree.SubElement( self.layers[min(1,len(self.layers)-1)], inkex.addNS('g','svg'), {"gcodetools": "Preview group"} )
        s, arcn = '', 0

        a,b,c = [0.,0.], [1.,0.], [0.,1.]
        k = (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1])
        a,b,c = self.transform(a, layer, True), self.transform(b, layer, True), self.transform(c, layer, True)
        if ((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))*k > 0 : reverse_angle = 1
        else : reverse_angle = -1
        for sk in curve:
            si = sk[:]
            si[0], si[2] = self.transform(si[0], layer, True), (self.transform(si[2], layer, True) if type(si[2])==type([]) and len(si[2])==2 else si[2])

            if s!='':
                if s[1] == 'line':
                    inkex.etree.SubElement(    group, inkex.addNS('path','svg'),
                            {
                                'style': style['line'],
                                'd':'M %s,%s L %s,%s' % (s[0][0], s[0][1], si[0][0], si[0][1]),
                                "gcodetools": "Preview",
                            }
                        )
                elif s[1] == 'arc':
                    arcn += 1
                    sp = s[0]
                    c = s[2]
                    s[3] = s[3]*reverse_angle

                    a =  ( (P(si[0])-P(c)).angle() - (P(s[0])-P(c)).angle() )%math.pi2 #s[3]
                    if s[3]*a<0:
                        if a>0:    a = a-math.pi2
                        else: a = math.pi2+a
                    r = math.sqrt( (sp[0]-c[0])**2 + (sp[1]-c[1])**2 )
                    a_st = ( math.atan2(sp[0]-c[0],- (sp[1]-c[1])) - math.pi/2 ) % (math.pi*2)
                    st = style['biarc%s' % (arcn%2)][:]
                    if a>0:
                        a_end = a_st+a
                        st = style['biarc%s'%(arcn%2)]
                    else:
                        a_end = a_st*1
                        a_st = a_st+a
                        st = style['biarc%s_r'%(arcn%2)]
                    inkex.etree.SubElement(    group, inkex.addNS('path','svg'),
                         {
                            'style': st,
                             inkex.addNS('cx','sodipodi'):        str(c[0]),
                             inkex.addNS('cy','sodipodi'):        str(c[1]),
                             inkex.addNS('rx','sodipodi'):        str(r),
                             inkex.addNS('ry','sodipodi'):        str(r),
                             inkex.addNS('start','sodipodi'):     str(a_st),
                             inkex.addNS('end','sodipodi'):       str(a_end),
                             inkex.addNS('open','sodipodi'):    'true',
                             inkex.addNS('type','sodipodi'):    'arc',
                            "gcodetools": "Preview",
                        })
            s = si


    def check_dir(self):
        if self.options.directory[-1] not in ["/","\\"]:
            if "\\" in self.options.directory :
                self.options.directory += "\\"
            else :
                self.options.directory += "/"
        print_("Checking direcrory: '%s'"%self.options.directory)
        if (os.path.isdir(self.options.directory)):
            if (os.path.isfile(self.options.directory+'header')):
                f = open(self.options.directory+'header', 'r')
                self.header = f.read()
                f.close()
            else:
                self.header = defaults['header']
            if (os.path.isfile(self.options.directory+'footer')):
                f = open(self.options.directory+'footer','r')
                self.footer = f.read()
                f.close()
            else:
                self.footer = defaults['footer']

            if self.options.unit == "G21 (All units in mm)" :
                self.header += "G21\n"
            elif self.options.unit == "G20 (All units in inches)" :
                self.header += "G20\n"
        else:
            self.error(_("Directory does not exist! Please specify existing directory at options tab!"),"error")
            return False

        if self.options.add_numeric_suffix_to_filename :
            dir_list = os.listdir(self.options.directory)
            if "." in self.options.file :
                r = re.match(r"^(.*)(\..*)$",self.options.file)
                ext = r.group(2)
                name = r.group(1)
            else:
                ext = ""
                name = self.options.file
            max_n = 0
            for s in dir_list :
                r = re.match(r"^%s_0*(\d+)%s$"%(re.escape(name),re.escape(ext) ), s)
                if r :
                    max_n = max(max_n,int(r.group(1)))
            filename = name + "_" + ( "0"*(4-len(str(max_n+1))) + str(max_n+1) ) + ext
            self.options.file = filename

        print_("Testing writing rights on '%s'"%(self.options.directory+self.options.file))
        try:
            f = open(self.options.directory+self.options.file, "w")
            f.close()
        except:
            self.error(_("Can not write to specified file!\n%s"%(self.options.directory+self.options.file)),"error")
            return False
        return True


################################################################################
###
###        Generate Gcode
###        Generates Gcode on given curve.
###
###        Crve defenitnion [start point, type = {'arc','line','move','end'}, arc center, arc angle, end point, [zstart, zend]]
###
################################################################################
    def generate_gcode(self, curve, layer, depth):
        tool = self.tools
        print_("Tool in g-code generator: " + str(tool))
        def c(c):
            c = [c[i] if i<len(c) else None for i in range(6)]
            if c[5] == 0 : c[5]=None
            s = [" X", " Y", " Z", " I", " J", " K"]
            r = ''
            for i in range(6):
                if c[i]!=None:
                    r += s[i] + ("%f" % (round(c[i],2))).rstrip('0')
            return r

        def calculate_angle(a, current_a):
            return  min(
                        [abs(a-current_a%math.pi2+math.pi2), a+current_a-current_a%math.pi2+math.pi2],
                        [abs(a-current_a%math.pi2-math.pi2), a+current_a-current_a%math.pi2-math.pi2],
                        [abs(a-current_a%math.pi2),             a+current_a-current_a%math.pi2])[1]
        if len(curve)==0 : return ""

        try :
            self.last_used_tool == None
        except :
            self.last_used_tool = None
        print_("working on curve")
        print_("Curve: " + str(curve) + "/n")
        g = ""

        lg, f =  'G00', "F%.1f"%tool['penetration feed']
        penetration_feed = "F%s"%tool['penetration feed']
        current_a = 0
        for i in range(1,len(curve)):
        #    Creating Gcode for curve between s=curve[i-1] and si=curve[i] start at s[0] end at s[4]=si[0]
            s, si = curve[i-1], curve[i]
            feed = f if lg not in ['G01','G02','G03'] else ''
            if s[1]    == 'move':
                g += "G00" + c(si[0]) + "\n" + tool['gcode before path'] + "\n"
                lg = 'G00'
            elif s[1] == 'end':
                g += tool['gcode after path'] + "\n"
                lg = 'G00'
            elif s[1] == 'line':
                if lg=="G00": g += "G01 " + feed + "\n"
                g += "G01" + c(si[0]) + "\n"
                lg = 'G01'
            elif s[1] == 'arc':
                r = [(s[2][0]-s[0][0]), (s[2][1]-s[0][1])]
                if lg=="G00": g += "G01 " + feed + "\n"
                if (r[0]**2 + r[1]**2)>.1:
                    r1, r2 = (P(s[0])-P(s[2])), (P(si[0])-P(s[2]))
                    if abs(r1.mag()-r2.mag()) < 0.001 :
                        g += ("G02" if s[3]<0 else "G03") + c(si[0]+[ None, (s[2][0]-s[0][0]),(s[2][1]-s[0][1])  ]) + "\n"
                    else:
                        r = (r1.mag()+r2.mag())/2
                        g += ("G02" if s[3]<0 else "G03") + c(si[0]) + " R%f" % (r) + "\n"
                    lg = 'G02'
                else:
                    g += "G01" + c(si[0]) + " " + feed + "\n"
                    lg = 'G01'
        if si[1] == 'end':
            g += tool['gcode after path'] + "\n"
        return g


    def get_transforms(self,g):
        root = self.document.getroot()
        trans = []
        while (g!=root):
            if 'transform' in g.keys():
                t = g.get('transform')
                t = simpletransform.parseTransform(t)
                trans = simpletransform.composeTransform(t,trans) if trans != [] else t               
            g=g.getparent()
        return trans


    def apply_transforms(self,g,csp):
        trans = self.get_transforms(g)           
        if trans != []:
            trans[1][2] = 0
            print_('appliing Transformation', time.time())
            print_(trans)
            simpletransform.applyTransformToPath(trans, csp)
        return csp


    def transform(self, source_point, layer, reverse=False):
        if layer == None :
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer not in self.transform_matrix:
            for i in range(self.layers.index(layer),-1,-1):
                if self.layers[i] in self.orientation_points :
                    break

            #print_(str(self.layers))
            #print_(str("I: " + str(i)))
            print_("Transform: " + str(self.layers[i]) , time.time())
            if self.layers[i] not in self.orientation_points:
                self.error(_("Orientation points for '%s' layer have not been found! Please add orientation points using Orientation tab!") % layer.get(inkex.addNS('label','inkscape')),"no_orientation_points")
            elif self.layers[i] in self.transform_matrix :
                self.transform_matrix[layer] = self.transform_matrix[self.layers[i]]
            else :
                orientation_layer = self.layers[i]
                if len(self.orientation_points[orientation_layer])>1 :
                    self.error(_("There are more than one orientation point groups in '%s' layer") % orientation_layer.get(inkex.addNS('label','inkscape')),"more_than_one_orientation_point_groups")
                points = self.orientation_points[orientation_layer][0]
                if len(points)==2:
                    points += [ [ [(points[1][0][1]-points[0][0][1])+points[0][0][0], -(points[1][0][0]-points[0][0][0])+points[0][0][1]], [-(points[1][1][1]-points[0][1][1])+points[0][1][0], points[1][1][0]-points[0][1][0]+points[0][1][1]] ] ]
                if len(points)==3:
                    print_("Layer '%s' Orientation points: " % orientation_layer.get(inkex.addNS('label','inkscape')), time.time())
                    for point in points:
                        print_(point)
                    #    Zcoordinates definition taken from Orientatnion point 1 and 2
                    self.Zcoordinates[layer] = [max(points[0][1][2],points[1][1][2]), min(points[0][1][2],points[1][1][2])]
                    matrix = numpy.array([
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

                    if numpy.linalg.det(matrix)!=0 :
                        m = numpy.linalg.solve(matrix,
                            numpy.array(
                                [[points[0][1][0]], [points[0][1][1]], [1], [points[1][1][0]], [points[1][1][1]], [1], [points[2][1][0]], [points[2][1][1]], [1]]
                                        )
                            ).tolist()
                        self.transform_matrix[layer] = [[m[j*3+i][0] for i in range(3)] for j in range(3)]

                    else :
                        self.error(_("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"),"wrong_orientation_points")
                else :
                    self.error(_("Orientation points are wrong! (if there are two orientation points they sould not be the same. If there are three orientation points they should not be in a straight line.)"),"wrong_orientation_points")

            self.transform_matrix_reverse[layer] = numpy.linalg.inv(self.transform_matrix[layer]).tolist()
            print_("\n Layer '%s' transformation matrixes:" % layer.get(inkex.addNS('label','inkscape')), time.time() )
            print_(self.transform_matrix)
            print_(self.transform_matrix_reverse)
            
        x,y = source_point[0], source_point[1]
        if not reverse :
            t = self.transform_matrix[layer]
        else :
            t = self.transform_matrix_reverse[layer]
        return [t[0][0]*x+t[0][1]*y+t[0][2], t[1][0]*x+t[1][1]*y+t[1][2]]


    def transform_csp(self, csp_, layer, reverse = False):
        csp = [  [ [csp_[i][j][0][:],csp_[i][j][1][:],csp_[i][j][2][:]]  for j in range(len(csp_[i])) ]   for i in range(len(csp_)) ]
        
        for i in xrange(len(csp)):
            for j in xrange(len(csp[i])):
                for k in xrange(len(csp[i][j])):
                    csp[i][j][k] = self.transform(csp[i][j][k],layer, reverse)
        return csp


################################################################################
###        Errors handling function, notes are just printed into Logfile,
###        warnings are printed into log file and warning message is displayed but
###        extension continues working, errors causes log and execution is halted
###        Notes, warnings adn errors could be assigned to space or comma or dot
###        sepparated strings (case is ignoreg).
################################################################################
    def error(self, s, type_= "Warning"):
        notes = "Note "
        warnings = """
                        Warning tools_warning
                        bad_orientation_points_in_some_layers
                        more_than_one_orientation_point_groups
                        more_than_one_tool
                        orientation_have_not_been_defined
                        tool_have_not_been_defined
                        selection_does_not_contain_paths
                        selection_does_not_contain_paths_will_take_all
                        selection_is_empty_will_comupe_drawing
                        selection_contains_objects_that_are_not_paths
                        """
        errors = """
                        Error
                        wrong_orientation_points
                        area_tools_diameter_error
                        no_tool_error
                        active_layer_already_has_tool
                        active_layer_already_has_orientation_points
                    """
        if type_.lower() in re.split("[\s\n,\.]+", errors.lower()) :
            print_(s)
            inkex.errormsg(s+"\n")
            sys.exit()
        elif type_.lower() in re.split("[\s\n,\.]+", warnings.lower()) :
            print_(s)
            if not self.options.suppress_all_messages :
                inkex.errormsg(s+"\n")
        elif type_.lower() in re.split("[\s\n,\.]+", notes.lower()) :
            print_(s)
        else :
            print_(s)
            inkex.errormsg(s)
            sys.exit()


################################################################################
###        Get defs from svg
################################################################################
    def get_defs(self) :
        self.defs = {}
        def recursive(g) :
            for i in g:
                if i.tag == inkex.addNS("defs","svg") :
                    for j in i:
                        self.defs[j.get("id")] = i
                if i.tag ==inkex.addNS("g",'svg') :
                    recursive(i)
        recursive(self.document.getroot())


################################################################################
###
###        Get Gcodetools info from the svg
###
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
                if i.tag == inkex.addNS("g",'svg') and i.get(inkex.addNS('groupmode','inkscape')) == 'layer':
                    self.layers += [i]
                    recursive_search(i,i)
                elif i.get('gcodetools') == "Gcodetools orientation group" :
                    points = self.get_orientation_points(i)
                    if points != None :
                        self.orientation_points[layer] = self.orientation_points[layer]+[points[:]] if layer in self.orientation_points else [points[:]]
                        print_("Found orientation points in '%s' layer: %s" % (layer.get(inkex.addNS('label','inkscape')), points))
                    else :
                        self.error(_("Warning! Found bad orientation points in '%s' layer. Resulting Gcode could be corrupt!") % layer.get(inkex.addNS('label','inkscape')), "bad_orientation_points_in_some_layers")
                elif i.tag == inkex.addNS('path','svg'):
                    if "gcodetools"  not in i.keys() :
                        self.paths[layer] = self.paths[layer] + [i] if layer in self.paths else [i]
                        if i.get("id") in self.selected :
                            self.selected_paths[layer] = self.selected_paths[layer] + [i] if layer in self.selected_paths else [i]
                elif i.tag == inkex.addNS("g",'svg'):
                    recursive_search(i,layer, (i.get("id") in self.selected) )
                elif i.get("id") in self.selected :
                    self.error(_("This extension works with Paths and Dynamic Offsets and groups of them only! All other objects will be ignored!\nSolution 1: press Path->Object to path or Shift+Ctrl+C.\nSolution 2: Path->Dynamic offset or Ctrl+J.\nSolution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file."),"selection_contains_objects_that_are_not_paths")
 

        recursive_search(self.document.getroot(),self.document.getroot())

    def get_orientation_points(self,g):
        items = g.getchildren()
        items.reverse()
        p2, p3 = [], []
        p = None
        for i in items:
            if i.tag == inkex.addNS("g",'svg') and i.get("gcodetools") == "Gcodetools orientation point (2 points)":
                p2 += [i]
            if i.tag == inkex.addNS("g",'svg') and i.get("gcodetools") == "Gcodetools orientation point (3 points)":
                p3 += [i]
        if len(p2)==2 : p=p2
        elif len(p3)==3 : p=p3
        if p==None : return None
        points = []
        for i in p :
            point = [[],[]]
            for  node in i :
                if node.get('gcodetools') == "Gcodetools orientation point arrow":
                    point[0] = self.apply_transforms(node,cubicsuperpath.parsePath(node.get("d")))[0][0][1]
                if node.get('gcodetools') == "Gcodetools orientation point text":
                    r = re.match(r'(?i)\s*\(\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*;\s*(-?\s*\d*(?:,|\.)*\d*)\s*\)\s*',node.text)
                    point[1] = [float(r.group(1)),float(r.group(2)),float(r.group(3))]
            if point[0]!=[] and point[1]!=[]:    points += [point]
        if len(points)==len(p2)==2 or len(points)==len(p3)==3 : return points
        else : return None

################################################################################
###
###        dxfpoints
###
################################################################################
    def dxfpoints(self):
        if self.selected_paths == {}:
            self.error(_("Noting is selected. Please select something to convert to drill point (dxfpoint) or clear point sign."),"warning")
        for layer in self.layers :
            if layer in self.selected_paths :
                for path in self.selected_paths[layer]:
                    if self.options.dxfpoints_action == 'replace':
                        path.set("dxfpoint","1")
                        r = re.match("^\s*.\s*(\S+)",path.get("d"))
                        if r!=None:
                            print_(("got path=",r.group(1)))
                            path.set("d","m %s 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000001 z" % r.group(1))
                            path.set("style",styles["dxf_points"])

                    if self.options.dxfpoints_action == 'save':
                        path.set("dxfpoint","1")

                    if self.options.dxfpoints_action == 'clear' and path.get("dxfpoint") == "1":
                        path.set("dxfpoint","0")

        
################################################################################
###
###        Engraving
###
################################################################################
    def engraving(self) :
        global x1,y1,rx,ry
        global cspm, wl
        global nlLT, i, j        
        global max_dist #minimum of tool radius and user's requested maximum distance
        
        def bisect((nx1,ny1),(nx2,ny2)) :
          
            cosBis = math.sqrt(max(0,(1.0+nx1*nx2-ny1*ny2)/2.0))
            #We can get correct sign of the sin, assuming cos is positive
            if (abs(ny1-ny2)< engraving_tolerance)  or (abs(cosBis) < engraving_tolerance) :
                if (abs(nx1-nx2)< engraving_tolerance): return(nx1,ny1,0.0)
                sinBis = math.copysign(1,ny1)
            else :
                sinBis = cosBis*(nx2-nx1)/(ny1-ny2)
            #We can correct signs by noting that the dot product
            # of bisector and either normal must be >0
            costurn=cosBis*nx1+sinBis*ny1
            if costurn == 0 : return (ny1*100,-nx1*100,1) #Path doubles back on itself
            sinturn=sinBis*nx1-cosBis*ny1
            if costurn<0 :  sinturn=-sinturn
            if 0 < sinturn*114.6 < (0) :
                sinturn=0 #set to zero if less than the user wants to see.
            return (cosBis/costurn,sinBis/costurn, sinturn)
            #end bisect

        def get_radius_to_line((x1,y1),(nx1,ny1), (nx2,ny2),(x2,y2),(nx23,ny23),(x3,y3),(nx3,ny3)):
           
            global max_dist

            #Start by converting coordinates to be relative to x1,y1
            x2,y2= x2-x1, y2-y1
            x3,y3= x3-x1, y3-y1
            
            #Make sure the line in question is facing x1,y1 and vice versa
            dist=-x2*nx23-y2*ny23
            if dist<0 : return max_dist
            denom=1.-nx23*nx1-ny23*ny1
            if denom < engraving_tolerance : return max_dist

            #radius and centre are:
            r=dist/denom
            cx=r*nx1
            cy=r*ny1
            #if c is not between the angle bisectors at the ends of the line, ignore
            #Use vector cross products. Not sure if I need the .0001 safety margins:
            if (x2-cx)*ny2 > (y2-cy)*nx2 +0.0001 :
                return max_dist 
            if (x3-cx)*ny3 < (y3-cy)*nx3  -0.0001 :
                return max_dist 
            return min(r, max_dist)
            #end of get_radius_to_line

        def get_radius_to_point((x1,y1),(nx,ny), (x2,y2)):
         
            global max_dist

            #Start by converting coordinates to be relative to x1,y1
            x2,y2= x2-x1, y2-y1
            denom=nx**2+ny**2-1
            if denom<=engraving_tolerance : #Not a corner bisector
                if denom==-1 : #Find circle centre x1,y1
                    return math.sqrt(x2**2+y2**2)
                #if x2,y2 not in front of the normal...
                if x2*nx+y2*ny <=0 : return max_dist
                #print_("Straight",x1,y1,nx,ny,x2,y2)
                return (x2**2+y2**2)/(2*(x2*nx+y2*ny) )
            #It is a corner bisector, so..
            discriminator = (x2*nx+y2*ny)**2 - denom*(x2**2+y2**2)
            if discriminator < 0 :
                return max_dist #this part irrelevant
            r=(x2*nx+y2*ny -math.sqrt(discriminator))/denom
            #print_("Corner",x1,y1,nx,ny,x1+x2,y1+y2,discriminator,r)
            return min(r, max_dist)
            #end of get_radius_to_point

        def bez_divide(a,b,c,d):
          
            bx=b[0]-a[0]
            by=b[1]-a[1]
            cx=c[0]-a[0]
            cy=c[1]-a[1]
            dx=d[0]-a[0]
            dy=d[1]-a[1]
            limit=8*math.hypot(dx,dy)/self.options.engraving_newton_iterations
            #LT This is the only limit we get from the user currently
            if abs(dx*by-bx*dy)<limit and abs(dx*cy-cx*dy)<limit :
                return [[a,b,c,d]]
            abx=(a[0]+b[0])/2.0
            aby=(a[1]+b[1])/2.0
            bcx=(b[0]+c[0])/2.0
            bcy=(b[1]+c[1])/2.0
            cdx=(c[0]+d[0])/2.0
            cdy=(c[1]+d[1])/2.0
            abcx=(abx+bcx)/2.0
            abcy=(aby+bcy)/2.0
            bcdx=(bcx+cdx)/2.0
            bcdy=(bcy+cdy)/2.0
            m=[(abcx+bcdx)/2.0,(abcy+bcdy)/2.0]
            return bez_divide(a,[abx,aby],[abcx,abcy],m) + bez_divide(m,[bcdx,bcdy],[cdx,cdy],d)
            #end of bez_divide

        def get_biggest((x1,y1),(nx,ny)):
           
            global max_dist, nlLT, i, j
            n1 = nlLT[j][i-1] #current node
            jjmin = -1
            iimin = -1
            r = max_dist
            # set limits within which to look for lines
            xmin, xmax = x1+r*nx-r, x1+r*nx+r
            ymin, ymax = y1+r*ny-r, y1+r*ny+r
            for jj in xrange(0,len(nlLT)) : #for every subpath of this object
                for ii in xrange(0,len(nlLT[jj])) : #for every point and line
                    if nlLT[jj][ii-1][2] : #if a point
                        if jj==j : #except this one
                            if abs(ii-i)<3 or abs(ii-i)>len(nlLT[j])-3 : continue
                        t1=get_radius_to_point((x1,y1),(nx,ny),nlLT[jj][ii-1][0] )
                        #print_("Try pt   i,ii,t1,x1,y1",i,ii,t1,x1,y1)
                    else: #doing a line
                        if jj==j : #except this one
                            if abs(ii-i)<2 or abs(ii-i)==len(nlLT[j])-1 : continue
                            if abs(ii-i)==2  and nlLT[j][(ii+i)/2-1][3]<=0 : continue
                            if (abs(ii-i)==len(nlLT[j])-2) and nlLT[j][-1][3]<=0 : continue
                        nx2,ny2 = nlLT[jj][ii-2][1]
                        x2,y2 = nlLT[jj][ii-1][0]
                        nx23,ny23 = nlLT[jj][ii-1][1]
                        x3,y3 = nlLT[jj][ii][0]
                        nx3,ny3 = nlLT[jj][ii][1]
                        if nlLT[jj][ii-2][3]>0 : #acute, so use normal, not bisector
                            nx2=nx23
                            ny2=ny23
                        if nlLT[jj][ii][3]>0 : #acute, so use normal, not bisector
                            nx3=nx23
                            ny3=ny23
                        x23min,x23max=min(x2,x3),max(x2,x3)
                        y23min,y23max=min(y2,y3),max(y2,y3)
                        #see if line in range
                        if n1[2]==False and (x23max<xmin or x23min>xmax or y23max<ymin or y23min>ymax) : continue
                        t1=get_radius_to_line((x1,y1),(nx,ny), (nx2,ny2),(x2,y2),(nx23,ny23), (x3,y3),(nx3,ny3))
                        #print_("Try line i,ii,t1,x1,y1",i,ii,t1,x1,y1)
                    if 0<=t1<r :
                        r = t1
                        iimin = ii
                        jjmin = jj
                        xmin, xmax = x1+r*nx-r, x1+r*nx+r
                        ymin, ymax = y1+r*ny-r, y1+r*ny+r
                #next ii
            #next jj
            return (jjmin,iimin,r)
            #end of get_biggest

        
        def save_point((x,y),w,i,j,ii,jj):
           
            global wl, cspm
            x=round(x,3) #round to 3 decimals
            y=round(y,3) #round to 3 decimals
            w=round(w,3) #round to 3 decimals
            if len(cspm)>1 :
                xy1a,xy1,xy1b,i1,j1,ii1,jj1=cspm[-1]
                w1=wl[-1]
                if i==i1 and j==j1 and ii==ii1 and jj==jj1 : #one match
                    xy1a,xy2,xy1b,i1,j1,ii1,jj1=cspm[-2]
                    w2=wl[-2]
                    if i==i1 and j==j1 and ii==ii1 and jj==jj1 : #two matches. Now test linearity
                        length1=math.hypot(xy1[0]-x,xy1[1]-y)
                        length2=math.hypot(xy2[0]-x,xy2[1]-y)
                        length12=math.hypot(xy2[0]-xy1[0],xy2[1]-xy1[1])
                        #get the xy distance of point 1 from the line 0-2
                        if length2>length1 and length2>length12 : #point 1 between them
                            xydist=abs( (xy2[0]-x)*(xy1[1]-y)-(xy1[0]-x)*(xy2[1]-y) )/length2
                            if xydist<engraving_tolerance : #so far so good
                                wdist=w2+(w-w2)*length1/length2 -w1
                                if abs(wdist)<engraving_tolerance :
                                    #print_("pop",j,i,xy1)
                                    cspm.pop()
                                    wl.pop()
            cspm+=[ [ [x,y],[x,y],[x,y],i,j,ii,jj ] ]
            wl+=[w]
            #end of save_point

        #end of subfunction definitions. engraving() starts here:
        ###########################################################

        global gcode
        r,w, wmax = 0,0,0 #theoretical and tool-radius-limited radii in pixels
        x1,y1,nx,ny =0,0,0,0
        cspe =[]
        we = []

        if len(self.selected_paths)<=0:
            self.error(_("Please select at least one path to engrave and run again."),"warning")
            return
        if not self.check_dir() : return
        #Find what units the user uses
        unit=" mm"
        if self.options.unit == "G20 (All units in inches)" :
            unit=" inches"
        elif self.options.unit != "G21 (All units in mm)" :
            self.error(_("Unknown unit selected. mm assumed"),"warning")
        
        #LT See if we can use this parameter for line and Bezier subdivision:
        bitlen=20/self.options.engraving_newton_iterations

        for layer in self.layers :
            if layer in self.selected_paths :
                #Calculate scale in pixels per user unit (mm or inch)
               
                #TODO
                orientation_scale = 1
              
                max_dist = self.options.laser_beam_with/2
                                
                engraving_group = inkex.etree.SubElement( self.selected_paths[layer][0].getparent(), inkex.addNS('g','svg') )
            
                for node in self.selected_paths[layer] :    
                    if node.tag == inkex.addNS('path','svg'):
                        cspi = cubicsuperpath.parsePath(node.get('d'))
                        #LT: Create my own list. n1LT[j] is for subpath j
                        nlLT = []
                        for j in xrange(len(cspi)): #LT For each subpath...
                            # Remove zero length segments, assume closed path
                            i = 0 #LT was from i=1
                            while i<len(cspi[j]):
                                if abs(cspi[j][i-1][1][0]-cspi[j][i][1][0])<engraving_tolerance and abs(cspi[j][i-1][1][1]-cspi[j][i][1][1])<engraving_tolerance:
                                    cspi[j][i-1][2] = cspi[j][i][2]
                                    del cspi[j][i]
                                else:
                                    i += 1
                        for csp in cspi:                            
                                           
                            #print_("csp is",csp)
                            nlLT.append ([])
                            for i in range(0,len(csp)): #LT for each point
                                #n = []
                                sp0, sp1, sp2 = csp[i-2], csp[i-1], csp[i]                                
                                #LT find angle between this and previous segment
                                x0,y0 = sp1[1]
                                nx1,ny1 = csp_normalized_normal(sp1,sp2,0)
                                #I don't trust this function, so test result
                                if abs(1-math.hypot(nx1,ny1))> 0.00001 :
                                    print_("csp_normalised_normal error t=0",nx1,ny1,sp1,sp2)
                                    self.error(_("csp_normalised_normal error. See log."),"warning")

                                nx0, ny0 = csp_normalized_normal(sp0,sp1,1)
                                if abs(1-math.hypot(nx0,ny0))> 0.00001 :
                                    print_("csp_normalised_normal error t=1",nx0,ny0,sp1,sp2)
                                    self.error(_("csp_normalised_normal error. See log."),"warning")
                                bx,by,s=bisect((nx0,ny0),(nx1,ny1))
                                #record x,y,normal,ifCorner, sin(angle-turned/2)
                                nlLT[-1] += [[ [x0,y0],[bx,by], True, s]]

                                #LT now do the line
                                if sp1[1]==sp1[2] and sp2[0]==sp2[1] : #straightline
                                    nlLT[-1]+=[[sp1[1],[nx1,ny1],False,i]]
                                else : #Bezier. First, recursively cut it up:
                                    nn=bez_divide(sp1[1],sp1[2],sp2[0],sp2[1])
                                    first=True #Flag entry to divided Bezier
                                    for bLT in nn : #save as two line segments
                                        for seg in range(3) :
                                            if seg>0 or first :
                                                nx1=bLT[seg][1]-bLT[seg+1][1]
                                                ny1=bLT[seg+1][0]-bLT[seg][0]
                                                l1=math.hypot(nx1,ny1)
                                                if l1<engraving_tolerance :
                                                    continue
                                                nx1=nx1/l1 #normalise them
                                                ny1=ny1/l1
                                                nlLT[-1]+=[[bLT[seg],[nx1,ny1], False,i]]
                                                first=False
                                            if seg<2 : #get outgoing bisector
                                                nx0=nx1
                                                ny0=ny1
                                                nx1=bLT[seg+1][1]-bLT[seg+2][1]
                                                ny1=bLT[seg+2][0]-bLT[seg+1][0]
                                                l1=math.hypot(nx1,ny1)
                                                if l1<engraving_tolerance :
                                                    continue
                                                nx1=nx1/l1 #normalise them
                                                ny1=ny1/l1
                                                #bisect
                                                bx,by,s=bisect((nx0,ny0),(nx1,ny1))
                                                nlLT[-1] += [[bLT[seg+1],[bx,by], True, 0.]]
                            
                       
                        reflex=False
                        for j in xrange(len(nlLT)): #LT6b for each subpath
                            cspm=[] #Will be my output. List of csps.
                            wl=[] #Will be my w output list
                            w = r = 0 #LT initial, as first point is an angle
                            for i in xrange(len(nlLT[j])) : #LT for each node
                                #LT Note: Python enables wrapping of array indices
                                # backwards to -1, -2, but not forwards. Hence:
                                n0 = nlLT[j][i-2] #previous node
                                n1 = nlLT[j][i-1] #current node
                                n2 = nlLT[j][i] #next node
                                #if n1[2] == True and n1[3]==0 : # A straight angle
                                    #continue
                                x1a,y1a = n1[0] #this point/start of this line
                                nx,ny = n1[1]
                                x1b,y1b = n2[0] #next point/end of this line
                                if n1[2] == True : # We're at a corner
                                    bits=1
                                    bit0=0
                                    #lastr=r #Remember r from last line
                                    lastw=w #Remember w from last line
                                    w = max_dist
                                    if n1[3]>0 : #acute. Limit radius
                                        len1=math.hypot( (n0[0][0]-n1[0][0]),( n0[0][1]-n1[0][1]) )
                                        if i<(len(nlLT[j])-1) :
                                            len2=math.hypot( (nlLT[j][i+1][0][0]-n1[0][0]),(nlLT[j][i+1][0][1]-n1[0][1]) )
                                        else:
                                            len2=math.hypot( (nlLT[j][0][0][0]-n1[0][0]),(nlLT[j][0][0][1]-n1[0][1]) )
                                        #set initial r value, not to be exceeded
                                        w = math.sqrt(min(len1,len2))/n1[3]
                                else: #line. Cut it up if long.
                                    if n0[3]>0 and not self.options.engraving_draw_calculation_paths :
                                        bit0=r*n0[3] #after acute corner
                                    else : bit0=0.0
                                    length=math.hypot((x1b-x1a),(y1a-y1b))
                                    bit0=(min(length,bit0))
                                    bits=int((length-bit0)/bitlen)
                                    #split excess evenly at both ends
                                    bit0+=(length-bit0-bitlen*bits)/2
                                    #print_("j,i,r,bit0,bits",j,i,w,bit0,bits)
                                for b in xrange(bits) : #divide line into bits
                                    x1=x1a+ny*(b*bitlen+bit0)
                                    y1=y1a-nx*(b*bitlen+bit0)
                                    jjmin,iimin,w=get_biggest( (x1,y1), (nx,ny))
                                    #print_("i,j,jjmin,iimin,w",i,j,jjmin,iimin,w)
                                    w = min(r, self.options.laser_beam_with)
                                    wmax=max(wmax,w)
                                    if reflex : #just after a reflex corner
                                        reflex = False
                                        if w<lastw :                                             
                                            save_point((n0[0][0]+n0[1][0]*w,n0[0][1]+n0[1][1]*w),w,i,j,iimin,jjmin)
                                    if n1[2] == True : # We're at a corner
                                        if n1[3]>0 : #acute
                                            save_point((x1+nx*w,y1+ny*w),w,i,j,iimin,jjmin)                                            
                                            save_point((x1,y1),0,i,j,iimin,jjmin)
                                        elif n1[3]<0  : #reflex
                                            if w>lastw :                                                
                                                wmax=max(wmax,w)
                                                save_point((x1+nx*w,y1+ny*w),w,i,j,iimin,jjmin)
                                    elif b>0 and n2[3]>0 and not self.options.engraving_draw_calculation_paths : #acute corner coming up
                                        if jjmin==j and iimin==i+2 : break                                
                                    save_point((x1+nx*w,y1+ny*w),w,i,j,iimin,jjmin)

                                #LT end of for each bit of this line
                                if n1[2] == True and n1[3]<0 : #reflex angle
                                    reflex=True
                                lastw = w #remember this w
                            #LT next i

                            cspm+=[cspm[0]]
                            print_("cspm/n",cspm, "/n")
                            #wl+=[wl[0]]
                            #print_("wl",wl)
                                            
                            if self.options.engraving_draw_calculation_paths==True:
                                node =  inkex.etree.SubElement(    engraving_group, inkex.addNS('path','svg'),                                         {
                                                         "d":     cubicsuperpath.formatPath([cspm]),
                                                        'style':    styles["biarc_style"]['biarc1'],
                                                        "lasertools": "Engraving calculation paths",
                                                    })
                                for i in xrange(len(cspm)):
                                    inkex.etree.SubElement(    engraving_group, inkex.addNS('path','svg'), 
                                            {"lasertools": "Engraving calculation paths", 'style':    "fill:none; fill-opacity:0.46; stroke:#000000; stroke-width:0.1;", inkex.addNS('cx','sodipodi'):        str(cspm[i][1][0]),  inkex.addNS('cy','sodipodi'):        str(cspm[i][1][1]),inkex.addNS('rx','sodipodi'):        str(wl[i]), inkex.addNS('ry','sodipodi'):        str(wl[i]), inkex.addNS('type','sodipodi'):    'arc'})
                            cspe += [cspm]
                            wluu = [] #width list in user units: mm/inches
                            for w in wl :
                                wluu+=[ w ]
                            #print_("wl in pixels",wl)
                            
                            we   +=    [wluu]
                       
                if cspe!=[]:                    
                    curve = self.parse_curve(cspe, layer, None, None) #convert to lines    
                    #print_("cspe= /n", cspe)                
                    self.draw_curve(curve, layer, engraving_group)                    
                    gcode += self.generate_gcode(curve, layer, 0)
                    
        if gcode!='':            
            self.export_gcode(gcode)

################################################################################
###
###        Orientation
###
################################################################################
    def orientation(self, layer=None) :
        print_("entering orientations")
        if layer == None :
            layer = self.current_layer if self.current_layer is not None else self.document.getroot()
        if layer in self.orientation_points:
            self.error(_("Active layer already has orientation points! Remove them or select another layer!"),"active_layer_already_has_orientation_points")

        orientation_group = inkex.etree.SubElement(layer, inkex.addNS('g','svg'), {"gcodetools":"Gcodetools orientation group"})

        # translate == ['0', '-917.7043']
        if layer.get("transform") != None :
            translate = layer.get("transform").replace("translate(", "").replace(")", "").split(",")
        else :
            translate = [0,0]

        # doc height in pixels (38 mm == 143.62204724px)
        doc_height = self.unittouu(self.document.getroot().xpath('@height', namespaces=inkex.NSS)[0])

        if self.document.getroot().get('height') == "100%" :
            doc_height = 1052.3622047
            print_("Overruding height from 100 percents to %s" % doc_height)

        print_("Document height: " + str(doc_height));

        if self.options.unit == "G21 (All units in mm)" :
            points = [[0.,0.,0.],[100.,0.,0.],[0.,100.,0.]]
            orientation_scale = 1
            print_("orientation_scale < 0 ===> switching to mm units=%0.10f"%orientation_scale )
        elif self.options.unit == "G20 (All units in inches)" :
            points = [[0.,0.,0.],[5.,0.,0.],[0.,5.,0.]]
            orientation_scale = 90
            print_("orientation_scale < 0 ===> switching to inches units=%0.10f"%orientation_scale )

        points = points[:2]

        print_(("using orientation scale",orientation_scale,"i=",points))
        for i in points :            
            si = [i[0]*orientation_scale, (i[1]*orientation_scale)+float(translate[1])]
            g = inkex.etree.SubElement(orientation_group, inkex.addNS('g','svg'), {'gcodetools': "Gcodetools orientation point (2 points)"})
            inkex.etree.SubElement(    g, inkex.addNS('path','svg'),
                {
                    'style':    "stroke:none;fill:#000000;",
                    'd':'m %s,%s 2.9375,-6.343750000001 0.8125,1.90625 6.843748640396,-6.84374864039 0,0 0.6875,0.6875 -6.84375,6.84375 1.90625,0.812500000001 z z' % (si[0], -si[1]+doc_height),
                    'gcodetools': "Gcodetools orientation point arrow"
                })
            t = inkex.etree.SubElement(    g, inkex.addNS('text','svg'),
                {
                    'style':    "font-size:10px;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;fill:#000000;fill-opacity:1;stroke:none;",
                    inkex.addNS("space","xml"):"preserve",
                    'x':    str(si[0]+10),
                    'y':    str(-si[1]-10+doc_height),
                    'gcodetools': "Gcodetools orientation point text"
                })
            t.text = "(%s; %s; %s)" % (i[0],i[1],i[2])

    ################################################################################
    ###
    ###        Fill area
    ###
    ################################################################################

    def area_fill(self):
        
        global gcode
        self.options.area_fill_angle = self.options.area_fill_angle * math.pi / 180
        if len(self.selected_paths)<=0:
            self.error(_("This extension requires at least one selected path."),"warning")
            return
        for layer in self.layers :
            if layer in self.selected_paths :                
                if self.options.laser_beam_with<=0 : 
                    self.error(_("Laser beam with must be > 0!"))
                
                for path in self.selected_paths[layer]:
                    lines = []
                    print_(("doing path",    path.get("style"), path.get("d")))
                    area_group = inkex.etree.SubElement( path.getparent(), inkex.addNS('g','svg') )
                    d = path.get('d')
                    if d==None:
                        print_("omitting non-path")
                        self.error(_("Warning: omitting non-path"),"selection_contains_objects_that_are_not_paths")
                        continue
                    csp = cubicsuperpath.parsePath(d)
                    csp = self.apply_transforms(path, csp)
                    csp = csp_close_all_subpaths(csp)
                    csp = self.transform_csp(csp, layer)
                    #maxx = max([x,y,i,j,root],maxx)
                    
                    # rotate the path to get bounds in defined direction.
                    a = - self.options.area_fill_angle
                    rotated_path = [   [ [ [point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(a)+point[1]*math.cos(a)]  for point in sp] for sp in subpath] for subpath in csp  ]
                    bounds =  csp_true_bounds(rotated_path)
                    
                    # Draw the lines 
                    # Get path's bounds
                    b = [0.0, 0.0, 0.0, 0.0] 
                    for k in range(4):
                        i, j, t = bounds[k][2], bounds[k][3], bounds[k][4]
                        b[k] = csp_at_t(rotated_path[i][j-1],rotated_path[i][j],t)[k%2]

                        
                    # Zig-zag
                    r = self.options.laser_beam_with
                    if r<=0 :
                        self.error('Tools diameter must be greater than 0!', 'error')
                        return

                    lines += [ [] ]

                    if self.options.area_fill_method == 'zig-zag' :
                        i = b[0] - self.options.area_fill_shift*r
                        top = True
                        last_one = True
                        while (i<b[2] or last_one) : 
                            if i>=b[2] : last_one = False
                            if lines[-1] == [] :
                                lines[-1] += [  [i,b[3]]  ]

                            if top :
                                lines[-1] += [ [i,b[1]],[i+r,b[1]] ]

                            else :
                                    lines[-1] += [ [i,b[3]], [i+r,b[3]] ]

                            top = not top
                            i += r
                    else :
                    
                        w, h  = b[2]-b[0] + self.options.area_fill_shift*r , b[3]-b[1] +  self.options.area_fill_shift*r
                        x,y = b[0] - self.options.area_fill_shift*r, b[1] - self.options.area_fill_shift*r
                        lines[-1] += [  [x,y] ]
                        stage = 0
                        start = True
                        while w>0 and h>0 :
                            stage = (stage+1)%4
                            if   stage == 0 :
                                y -= h
                                h -= r
                            elif stage == 1:
                                x += w
                                if not start:
                                    w -= r
                                start = False    
                            elif stage == 2 :
                                y += h
                                h -= r
                            elif stage == 3:
                                x -= w
                                w -=r
                            
                            lines[-1] += [ [x,y] ]

                        stage = (stage+1)%4                            
                        if w <= 0 and h>0 :
                            y = y-h if stage == 0 else y+h  
                        if h <= 0  and w>0 :
                            x = x-w if stage == 3 else x+w  
                        lines[-1] += [ [x,y] ]
                    # Rotate created paths back
                    a =  self.options.area_fill_angle
                    lines = [ [ [point[0]*math.cos(a) - point[1]*math.sin(a), point[0]*math.sin(a)+point[1]*math.cos(a)] for point in subpath] for subpath in lines  ]

                    # get the intersection points
                    
                    splitted_line = [ [lines[0][0]] ]
                    intersections = {}
                    for l1,l2, in zip(lines[0],lines[0][1:]): 
                        ints = []
                        
                        if l1[0]==l2[0] and l1[1]==l2[1] : continue
                        for i in range(len(csp)) :
                            for j in range(1,len(csp[i]))  :
                                sp1,sp2 = csp[i][j-1], csp[i][j]
                                roots = csp_line_intersection(l1,l2,sp1,sp2)
                                for t in roots :
                                    p = tuple(csp_at_t(sp1,sp2,t))
                                    if l1[0]==l2[0] :
                                        t1 = (p[1]-l1[1])/(l2[1]-l1[1])
                                    else :
                                        t1 = (p[0]-l1[0])/(l2[0]-l1[0])
                                    if 0<=t1<=1    :
                                        ints += [[t1, p[0],p[1], i,j,t]]
                                        if p in intersections :
                                            intersections[p]  += [ [i,j,t] ]  
                                        else :     
                                            intersections[p]  = [ [i,j,t] ]  
                                        #p = self.transform(p,layer,True)
                                        #draw_pointer(p)
                        ints.sort()                                                
                        for i in ints:
                            splitted_line[-1] +=[ [ i[1], i[2]] ]
                            splitted_line += [ [ [ i[1], i[2]] ] ]
                        splitted_line[-1] += [  l2  ]
                        i = 0
                                                           
                    while i < len(splitted_line) :
                        # check if the middle point of the first lines segment is inside the path.
                        # and remove the subline if not.
                        l1,l2 = splitted_line[i][0],splitted_line[i][1]
                        p = [(l1[0]+l2[0])/2, (l1[1]+l2[1])/2]
                        if not point_inside_csp(p, csp):
                            #i +=1                         
                            del splitted_line[i]
                        else :
                            i += 1                    
                    
                    # if we've used spiral method we'll try to save the order of cutting
                    do_not_change_order = self.options.area_fill_method == 'spiral' 
                                       
                    # and apply back transrormations to draw them
                    csp_line = csp_from_polyline(splitted_line)
                    csp_line = self.transform_csp(csp_line, layer, True)
                                
                    curve = self.parse_curve(csp_line, layer)
                    self.draw_curve(curve, layer, area_group)
                    gcode += self.generate_gcode(curve, layer, 0)
                    
                    
        if gcode!='' and not(self.options.add_contours):            
            self.export_gcode(gcode)                    

################################################################################
###
###        Effect
###
###        Main function of Gcodetools class
###
################################################################################
    def effect(self) :
        global options
        options = self.options
        options.self = self
        options.doc_root = self.document.getroot()
        # define print_ function
        global print_
        if self.options.log_create_log :
            try :
                if os.path.isfile(self.options.log_filename) : os.remove(self.options.log_filename)
                f = open(self.options.log_filename,"a")
                f.write("Gcodetools log file.\nStarted at %s.\n%s\n" % (time.strftime("%d.%m.%Y %H:%M:%S"),options.log_filename))
                f.write("%s tab is active.\n" % self.options.active_tab)
                f.close()
            except :
                print_  = lambda *x : None
        else : print_  = lambda *x : None
        self.get_info()        
        
        if self.orientation_points == {} :
            self.error(_("Orientation points have not been defined! A default set of orientation points has been automatically added."),"warning")
            self.orientation( self.layers[min(0,len(self.layers)-1)] )
            self.get_info()
        

        self.tools = {
            "name": "Laser Engraver",
            "id": "Laser Engraver",
            "penetration feed": self.options.laser_speed,
            "feed": self.options.laser_speed,
            "gcode before path": ("G04 P" + self.options.power_delay + " " + self.options.laser_command ),
            "gcode after path": (self.options.laser_off_command  + "\n" + "G00")
        }       
              
        self.get_info()
        if self.options.add_infill:        
            self.selected_paths = self.paths
            print_("Start filling area at", time.time())
            self.area_fill()
            
        if self.options.add_contours: 
            self.selected_paths = self.paths
            print_("Start building contour at", time.time())
            self.get_info() 
            self.engraving()
               

e = laser_gcode()
e.affect()
