#/usr/bin/env python3

import sys
import math
import numpy as np
from numpy import linalg
from scipy.special import comb

def circ_rad(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    c = (x1-x2)**2 + (y1-y2)**2
    a = (x2-x3)**2 + (y2-y3)**2
    b = (x3-x1)**2 + (y3-y1)**2
    s = 2*(a*b + b*c + c*a) - (a*a + b*b + c*c)
    px = (a*(b+c-a)*x1 + b*(c+a-b)*x2 + c*(a+b-c)*x3) / s
    py = (a*(b+c-a)*y1 + b*(c+a-b)*y2 + c*(a+b-c)*y3) / s
    ar = math.sqrt(a)
    br = math.sqrt(b)
    cr = math.sqrt(c)
    r = ar*br*cr / math.sqrt((ar+br+cr)*(-ar+br+cr)*(ar-br+cr)*(ar+br-cr))
    return (px, py), r


def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """

    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i


def bezier_curve(points, nTimes=1000):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.

       points should be a list of lists, or list of tuples
       such as [ [1,1],
                 [2,3],
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000

        See http://processingjs.nihongoresources.com/bezierinfo/
    """

    nPoints = len(points)
    xPoints = np.array([p[0] for p in points])
    yPoints = np.array([p[1] for p in points])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)

    return xvals, yvals

def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

class RoutPlanner:
    def __init__(self, rad):
        self.stations = []
        self.round = 0
        self.collet_mode = True
        self.rad = rad
        self.thrust = 100
        self.rad2dev = [(45000, 5),
                        (23000, 10),
                        (16000, 15),
                        (11500, 20),
                        (8600, 25),
                        (6700, 30),
                        (5200, 35),
                        (4350, 40),
                        (3500, 45),
                        (2800, 50),
                        (2200, 55),
                        (1760, 60),
                        (1390, 65),
                        (1100, 70),
                        (980, 75),
                        (970, 80)]
        self.rad2dev.reverse()

    def add_station(self, p):
        if len(self.stations) and p == self.stations[0]:
            self.collet_mode = False
        else:
            self.stations.append(p)

    def plan(self, pos, angle, target, last_pos):
        self.rad = 500
        i = self.stations.index(target)
        if i != 0:
            self.stations = self.stations[i:] + self.stations[:i]
        p1 =(self.rad, 0.)
        #h = np.angle(pos[0] + pos[1]*1j, deg=True) + angle
        traj = np.array(pos) - np.array(last_pos)
        h = np.angle(traj[0] + traj[1]*1j, deg=True)
        print("traj={} h={}".format(traj, h), file=sys.stderr)
        p1 = rotate(p1, degrees=h)
        p1 = np.array(p1) + np.array(pos)
        p21 = np.array(pos) - np.array(target)
        p21 = p21/linalg.norm(p21)
        p23 = np.array(self.stations[1]) - np.array(target)
        p23 = p23/linalg.norm(p23)
        ph2 = p21 + p23
        ph2 = ph2/linalg.norm(ph2)
        p2 = rotate(ph2, degrees=90)
        p2 = p2*self.rad
        p2 = p2 + target
        verts = [pos, p1, p2, target]
        return bezier_curve(verts, nTimes=30)

    def calc_dev(self, r, fac=1.0):
        for (calr1, dev1),(calr2, dev2) in zip(self.rad2dev[:-1], self.rad2dev[1:]):
            calr1 *= fac
            calr2 *= fac
            if r >= calr1 and r < calr2:
                f = float(r - calr1)/float(calr2 - calr1)
                ang = f * dev1 + (1.0 - f) * dev2
                print ("dev = " + str(ang), file=sys.stderr)

        if ang > 30 and self.thrust == 100:
            self.thrust = 50
            return self.calc_dev(r, 0.5)
        else:
            self.thrust = 100
        return ang

    def calc_direction(self, pos, angle, target, last_pos):
        self.thrust = 100
        xvals, yvals = self.plan(pos, angle, target, last_pos)
        print ("pos={}, x0={} y0={}".format(pos, xvals[0], yvals[0]), file=sys.stderr)
        c, r = circ_rad(pos,
                (xvals[0], yvals[0]),
                (xvals[1], yvals[1]))

        pc = np.array([xvals[0] - pos[0], yvals[0] - pos[1]])
        crel = np.array(c) - np.array(pos)
        d = np.cross(pc, crel)
        dev = self.calc_dev(r)

        if d < 0:
            dev = -dev
        print("c = {} r = {} dev={} d={}".format(c,r,dev,d ), file=sys.stderr)
        print (pos, file=sys.stderr)
        pc = pc / linalg.norm(pc)
        print ("pc0={}".format(pc), file=sys.stderr)
        pc = rotate(pc, degrees=dev)
        print ("pc1={}".format(pc), file=sys.stderr)
        pc *= 300
        pc = pos + pc
        print ("using bezier", file=sys.stderr)
        return pc, self.thrust, False

class Naive:
    def __init__(self, rad):
        self.ln = (-1, -1)
        self.rp = RoutPlanner(rad)
        self.boosted = False
        self.lastInitiated = False
        self.last_pos = (0, 0)

    def calc_direction(self, pos, dist, angle, target):
        opt = False
        vel = -1
        if self.lastInitiated:
            vel = math.sqrt((pos[0] - self.last_pos[0])**2 + (pos[1] - self.last_pos[1])**2)
        self.lastInitiated = True
        print ("vel=" + str(vel), file=sys.stderr)
        if self.rp.collet_mode == False:
            opt = True
            if vel > 0:
                self.rp.rad = vel * 1.5
            target, thrust, boost = self.rp.calc_direction(pos, angle, target, self.last_pos)
        elif self.ln != target:
            self.rp.add_station(target)
            self.ln = target

        self.last_pos = pos
        ax = math.atan(300.0/float(dist)) * 180 / 3.141

        boost = False
        if not opt:
            thrust = 100
        aangle = math.fabs(angle)
        eangle = aangle - ax
        print ("--" + str(eangle) + " " + str(dist), file=sys.stderr)

        #if eangle > 0 and not opt:
        #    thrust = int(100 * ( 1.0 -  eangle/60.))
        if not opt:
            if dist < 2800:
                thrust = 50
            if dist > 8000 and eangle < 5 and not self.boosted:
                boost = True
                self.boosted = True
            if eangle > 3 and dist < 2000:
                thrust = 10
        return target, thrust, boost

class calib_circ:

    def __init__(self, dev):
        self.hist = []
        self.r = []
        self.cent = []
        self.dev = dev

    def calc_direction(self, pos, dist, angle, n):
        if not len(self.hist):
            target = n
        else:
            dr = np.array(pos) - np.array(self.hist[-1])
            dr = dr / linalg.norm(dr)
            dr = rotate(dr, degrees=self.dev)
            dr = dr * 100
            dr = np.array(pos) + dr
            target = (dr[0], dr[1])
        if len(self.hist) > 2:
            c, r = circ_rad(self.hist[-2], self.hist[-1], pos)
            print ("r = {}".format(r), file=sys.stderr)
        self.hist.append(pos)
        return target, 50, False

if __name__ == "__main__":
    calib = True
    lnx, lny = -1, -1
    n = 0
    hist = []
    opt = False
    algo = calib_circ(45) if calib else Naive(700.)

    # game loop
    while True:
        # x: x position of your pod
        # y: y position of your pod
        # next_checkpoint_x: x position of the next check point
        # next_checkpoint_y: y position of the next check point

        inp = [int(i) for i in input().split()]

        i = -1

        x, y, nx, ny, dist, angle = inp

        [int(i) for i in input().split()]
        target, thrust, boost = algo.calc_direction((x, y), dist, angle, (nx, ny))

        if boost:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " BOOST")
        else:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " " + str(thrust))




