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

def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

def calc_node_approach(p1, p2, p3, rad):
    p21 = np.array(p1) - np.array(p2)
    p21 = p21 / linalg.norm(p21)
    p23 = np.array(p3) - np.array(p2)
    p23 = p23 / linalg.norm(p23)
    angle = np.math.atan2(np.linalg.det([p21,p23]),np.dot(p21,p23))
    angle_deg = np.degrees(angle)
    pivot = p21 + p23
    pivot = (pivot / linalg.norm(pivot)) * rad
    par1 = pivot * ( 1 - np.cos(angle / 2.0))
    cross = np.cross(p23, p21)
    rot90 = 90 if cross > 0 else -90
    ppepend = rotate(pivot, degrees=rot90)
    ppepend = ppepend / linalg.norm(ppepend)
    ppepend = ppepend * rad * np.sin(angle / 2.0)
    target = ppepend + par1
    target0 = ppepend + ((-1)*par1)
    target += np.array(p2)
    target0 += np.array(p1)
    pivot = pivot + np.array(p2)
    return (pivot[0], pivot[1]), ((target0[0], target0[1]), (target[0], target[1])), angle

class BreakBeforeTarget:
    def __init__(self):
        pass

    def get_thrust(self, dist, angle, target, rad):
        return 50 if dist < 2800 else 100

class CurvatureThustStrategy:
    def __init__(self):
        self.thrust = 100
        pass

    def get_thrust(self, dist, angle, target, rad):
        self.thrust = int(round(r/3300. * 100))
        if self.thrust > 100:
            self.thrust = 100
        elif self.thrust < 20:
            self.thrust = 20
        return self.thrust

def closest_point_to_segment(p, q, x):
    """
    return the closest point (as coefficient) from point x to line p-q
    """
    x = np.array(x)
    q = np.array(q)
    p = np.array(p)
    v = q - p
    s = x - p
    alpha = np.dot(s, v)/np.dot(v, v)
    return alpha


class Planner:
    def __init__(self, fac):
        self.fac = fac
        self.pivot = (0, 0)
        self.curve_st = [(0, 0), (0, 0)]
        self.r100 = 1200
        self.state = 0
        self.stations = []
        self.thrust_stg = CurvatureThustStrategy()
        self.thrust, self.thrust2 = 0, 0

    def plan(self, post, dist, angle, target, stations):
        x1 = stations[-1][0]
        y1 = stations[-1][1]
        x2 = stations[0][0]
        y2 = stations[0][1]
        d = math.sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))
        if d > (2 * self.r100 + 600):
            rad = self.r100
        else:
            rad = (d - 600) / 2.0
        self.fac = rad / d
        self.pivot, self.curve_st, self.angle = calc_node_approach(stations[-1], stations[0], stations[1], rad)
        self.state = 0
        self.stations = stations
        self.thrust2 = self.curve_st(rad)
        self.thrust = 100
        self.max_alpha = 1.0 + 300.0 / d

    def act(self, post, dist, angle, target, last_pos, p_center, rad):
        alpha = closest_point_to_segment(self.stations[-1], self.stations[0])
        if alpha < self.fac:
            return self.stations[0], self.thrust2, False
        elif alpha < (1.0 - self.fac):
            return self.stations[1], self.thrust, False
        elif alpha <= self.max_alpha:
            return target, self.thrust2, False
        else:
            return target, 30, False

class Collector:
    def __init__(self):
        self.stations = []
        self.collecting = True
        self.last_target = (-1, -1)
        self.last_pos = []
        self.algo = Basic()
        self.thrustStrategies = [CurvatureThustStrategy()]

    def act(self, pos, dist, angle, target):
        if self.collecting and len(self.stations) and self.stations[0] == target:
            self.collecting = False
            self.algo = Planner()

        if self.collecting and self.last_target != target:
            self.stations.append(target)
        elif self.last_target != target:
            self.stations = self.stations[1:].append(self.stations[0])

        if self.last_target != target:
            self.algo.plan(pos, dist, angle, target, self.stations)
            self.last_target = target

        rad = -1
        if len(self.last_pos) >= 3:
            p_cneter, rad = circ_rad(self.last_pos[-2], self.last_pos[-1], pos)
        new_target, thrust, boost = self.algo.act(pos, dist, angle, target, self.last_pos, p_cneter, rad)
        self.last_pos.append(pos)
        if len(self.last_pos) > 10:
            self.last_pos = self.last_pos[1:]
        for strategy in self.thrustStrategies:
            thr2 = strategy.get_thrust(pos, dist, angle, target, rad)
            if thr2 < thrust:
                thr2 = thrust
        return new_target, thrust, boost

if __name__ == "__main__":
    calib = False
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

