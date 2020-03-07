#/usr/bin/env python3
#/usr/bin/env python3

import sys
import math
import numpy as np
from numpy import linalg
from scipy.special import comb

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

    def add_station(self, p):
        if len(self.stations) and p == self.stations[0]:
            self.collet_mode = False
        else:
            self.stations.append(p)

    def plan(self, pos, angle, target):
        i = self.stations.index(target)
        if i != 0:
            self.stations = self.stations[i:] + self.stations[:i]
        p1 =(self.rad, 0.)
        h = np.angle(pos[0] + pos[1]*1j, deg=True) + angle
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

    def calc_direction(self, pos, angle, target):
        xvals, yvals = self.plan(pos, angle, target)
        pc = np.array([xvals[1] - pos[0], yvals[1] - pos[1]])
        pc = pc / linalg.norm(pc)
        pc *= 300
        pc = pos + pc
        print ("using bezier", file=sys.stderr)
        return pc, 100, False

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
        self.last_pos = pos
        self.lastInitiated = True
        print ("vel=" + str(vel), file=sys.stderr)
        if self.rp.collet_mode == False:
            opt = True
            if vel > 0:
                self.rp.rad = vel * 1.5
            target, thrust, boost = self.rp.calc_direction(pos, angle, target)
        elif self.ln != target:
            self.rp.add_station(target)
            self.ln = target

        ax = math.atan(300.0/float(dist)) * 180 / 3.141

        boost = False
        thrust = 100
        aangle = math.fabs(angle)
        eangle = aangle - ax
        print ("--" + str(eangle) + " " + str(dist), file=sys.stderr)

        #if eangle > 0 and not opt:
        #    thrust = int(100 * ( 1.0 -  eangle/60.))
        if dist < 2800:
            thrust = 50
        if dist > 8000 and eangle < 5 and not self.boosted:
            boost = True
            self.boosted = True
        if eangle > 3 and dist < 2000:
            thrust = 10

        return target, thrust, boost

if __name__ == "__main__":
    lnx, lny = -1, -1
    n = 0
    hist = []
    opt = False
    algo = Naive(700.)
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


