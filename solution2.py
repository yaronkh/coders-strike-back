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

def half_angle(p1, p2, p3):
    p21 = np.array(p1) - np.array(p2)
    p21 = p21 / linalg.norm(p21)
    p23 = np.array(p3) - np.array(p2)
    p23 = p23 / linalg.norm(p23)
    angle = np.math.atan2(np.linalg.det([p21,p23]),np.dot(p21,p23))
    angle = np.degrees(angle)
    ph = p21 + p23
    ph = ph / linalg.norm(ph)
    return (ph[0], ph[1]), angle

def calc_pivot_point(p1, p2, p3, r):
    ph, angle = half_angle(p1, p2, p3)

class BrakeBeforeTarget:
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

