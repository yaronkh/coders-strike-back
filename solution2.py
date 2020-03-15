import sys
import math
import numpy as np
from numpy import linalg
from scipy.special import comb

class Params:
    r100 = 1800.
    thrust_rad_100 = 3300.
    dist_break_distance = 1800.
    minimal_straight_dist = 1000.
    break_dist = 0
    break_fac = 0.25

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
    #calculate relative unit vectors
    p21 = np.array(p1) - np.array(p2)
    p21 = p21 / linalg.norm(p21)
    p23 = np.array(p3) - np.array(p2)
    p23 = p23 / linalg.norm(p23)
    #calculate the angle of stations2
    angle = np.math.atan2(np.linalg.det([p21,p23]),np.dot(p21,p23))
    angle_deg = np.degrees(angle)
    if angle_deg > 175.0 or angle_deg < -175:
        t_pivot = p21 * rad + np.array(p2)
        target = p2
        target0 = p21 * 300 + np.array(p2)
    else:
        #pivot is the curve central point
        pivot = p21 + p23
        pivot = (pivot / linalg.norm(pivot)) * rad
        t_pivot = pivot + p2

        cross = np.cross(p21, p23)
        rot90 = 90 if cross < 0 else -90
        #direct is the vector from the curve central point to the rotation start point
        direct = rotate(p21, degrees=rot90) * rad
        break_vec = p21 * Params.break_dist
        target = t_pivot + direct + break_vec
        target0 = target + p21 * Params.minimal_straight_dist
    fac1 = closest_point_to_segment(p1, p2, target) - Params.break_fac
    fac0 = closest_point_to_segment(p1, p2, target0) - Params.break_fac

    return (t_pivot[0], t_pivot[1]), ((target0[0], target0[1]), (target[0], target[1])), (fac0, fac1), angle_deg

class BreakBeforeTarget:
    def __init__(self):
        pass

    def get_thrust(self, dist, angle, target, rad):
        return 50 if dist < Params.dist_break_distance else 100

class CurvatureThustStrategy:
    def __init__(self):
        self.thrust = 100
        pass

    def get_thrust(self, dist, angle, target, rad):
        self.thrust = int(round(rad/Params.thrust_rad_100 * 100))
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
    def __init__(self):
        self.pivot = (0, 0)
        self.curve_st = [(0, 0), (0, 0)]
        self.r100 = Params.r100
        self.state = 0
        self.stations = []
        self.thrust_stg = CurvatureThustStrategy()
        self.thrust, self.thrust2 = 0, 0
        self.rfacs = (0, 0)
        self.rad = 0

    def plan(self, post, dist, angle, target, stations):
        x1 = stations[-1][0]
        y1 = stations[-1][1]
        x2 = stations[0][0]
        y2 = stations[0][1]
        d = math.sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))
        if d > (2 * self.r100 + Params.minimal_straight_dist):
            rad = self.r100
        else:
            rad = (d - 600) / 2.0
        self.rad = rad
        print ("PLAN: {},{},{}".format(stations[-1], stations[0], stations[1]), file=sys.stderr)
        self.pivot, self.curve_st, self.rfacs, self.angle = calc_node_approach(stations[-1], stations[0], stations[1], rad)
        print ("RES: pivpt={},st={},facts={},angle={}".format(self.pivot, self.curve_st, self.rfacs, self.angle), file=sys.stderr)
        self.state = 0
        self.stations = stations
        self.thrust2 = 100
        self.thrust = 100
        self.max_alpha = 1.0 + 300.0 / d

    def act(self, pos, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        alpha = closest_point_to_segment(self.stations[-1], self.stations[0], pos)
        if alpha < self.rfacs[0]:
            return self.curve_st[0], self.thrust2, False
        elif alpha < self.rfacs[1]:
            return self.curve_st[1], 60, False
        elif alpha <= self.max_alpha:
            thrust = 50 if angle > 1 else 90
            return target, thrust, False
        else:
            thrust = 30 if angle > 3 else 100
            return target, thrust, False

class Clash:
    def __init__(self):
        self.dist = -1

    def plan(self, post, dist, angle, target, stations):
        pass

    def act(self, pos, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        ox , oy = opponent_pos
        mx, my = pos
        self.dist = math.sqrt((ox - mx)*(ox - mx) + (oy - my)*(oy - my))
        print ("CLASH DIST={}".format(self.dist), file=sys.stderr)
        return opponent_pos, 100, True

class BlindPlanner:
    def __init__(self):
        pass

    def plan(self, post, dist, angle, target, stations):
        pass

    def act(self, post, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        boost = False
        if dist > 7000 and angle < 5:
            boost = True
            thrust = 100
        elif dist < 2800:
            thrust = 80 if angle < 1 else 50
        else:
            thrust = 100
        return target, thrust, boost

class Collector:
    def __init__(self):
        self.stations = []
        self.collecting = True
        self.last_target = (-1, -1)
        self.last_pos = []
        self.algo = Clash()
        self.clash = True
        self.thrustStrategies = []
        self.lap = 1
        self.stat_in_lap = 0

    def act(self, pos, dist, angle, target, opponent_pos):
        if self.last_target != target:
            self.stat_in_lap += 1
            print ("STATION={}".format(self.stat_in_lap), file=sys.stderr)

        if self.collecting and len(self.stations) and self.last_target != target and self.stations[0] == target:
            self.collecting = False
            self.lap = 1
            print ("COLLECTED: {}".format(self.stations), file=sys.stderr)
            self.algo = Planner()

        if self.collecting and self.last_target != target:
            self.stations.append(target)

        if not self.collecting and self.last_target != target and ((self.stat_in_lap % len(self.stations)) == 1):
            print ("LAP: {}".format(self.lap), file=sys.stderr)
            self.lap +=  1
        if self.lap == 3 and (self.stat_in_lap % len(self.stations) == 0) and self.last_target != target:
            print ("PHOTOFINISH..........", file=sys.stderr)
            self.algo = BlindPlanner()

        if self.last_target != target:
            self.algo.plan(pos, dist, angle, target, self.stations)

        rad = -1
        if 0 and len(self.last_pos) >= 3:
            p_center, rad = circ_rad(self.last_pos[-2], self.last_pos[-1], pos)
        else:
            p_center, rad = (0., 0.), 23000
        new_target, thrust, boost = self.algo.act(pos, dist, angle, target, self.last_pos, p_center, rad, opponent_pos)
        if self.clash and self.collecting and self.algo.dist > 2000:
            self.clash = False
            self.algo = BlindPlanner()
        if not self.collecting and self.last_target != target:
            print ("ROTATING", file=sys.stderr)
            self.stations = self.stations[1:] +[self.stations[0]]

        else:
            print ("NOT ROTATING", file=sys.stderr)
        self.last_target = target
        self.last_pos.append(pos)
        if len(self.last_pos) > 10:
            self.last_pos = self.last_pos[1:]
        for strategy in self.thrustStrategies:
            thr2 = strategy.get_thrust(dist, angle, target, rad)
            if thr2 < thrust:
                thr2 = thrust
        return new_target, thrust, boost

if __name__ == "__main__":
    calib = False
    lnx, lny = -1, -1
    n = 0
    hist = []
    opt = False
    algo = calib_circ(45) if calib else Collector()

    # game loop
    while True:
        # x: x position of your pod
        # y: y position of your pod
        # next_checkpoint_x: x position of the next check point
        # next_checkpoint_y: y position of the next check point

        inp = [int(i) for i in input().split()]

        i = -1

        x, y, nx, ny, dist, angle = inp

        ox, oy = [int(i) for i in input().split()]
        target, thrust, boost = algo.act((x, y), dist, angle, (nx, ny), (ox, oy))

        if boost:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " BOOST")
        else:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " " + str(thrust))

