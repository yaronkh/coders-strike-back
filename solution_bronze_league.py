import sys
import time
import math
import numpy as np
from numpy import linalg
from scipy.special import comb

class Params:
    def __init__(self):
        self.r1planner = Planner()
        #self.r1planner = BlindPlanner()
        self.r100 = 1800.
        self.thrust_rad_100 = 3300.
        self.dist_break_distance = 500
        self.minimal_straight_dist = 1000.
        self.break_dist = 0
        self.break_fac = 0.30
        self.chase_max = 2000
        self.chase_max_angle = 5
        self.side_puch_distance = 1000
        self.side_punch_max_angle = 45
        self.num_punchs = 1
        self.Kp = 0.1
        self.Ki = 0
        self.Kd = 0
        self.vel_const = 8.2
        self.gtKp = 1.
        self.gtKi = 0
        self.gtKd = 0


class ShoshParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.dist_break_distance = 0

class TrampolineParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.dist_break_distance = 0

class MacbilitParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.r100 = 1800.
        self.dist_break_distance = 0

class MehumashParams(Params):
    def __init__(self):
        Params.__init__(self)
        #self.vel_const = 2.
        self.r100 = 3000.
        #self.break_fac = 0.0
        self.break_fac = 0.14
        self.vel_const = 18
        #self.r100 = 6000
        self.dist_break_distance = 0

class HostileParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.break_fac = 0.3
        self.r100 = 2200.
        self.vel_const = 8.2
        self.dist_break_distance = 0


class ZigzagParams(Params):
    def __init__(self):
        Params.__init__(self)
        #self.vel_const = 2.
        self.r100 = 2000.
        #self.break_fac = 0.0
        self.break_fac = 0.2
        self.vel_const = 8.2
        #self.r100 = 6000
        self.dist_break_distance = 0

class ArrowParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.dist_break_distance = 0
        self.vel_const = 6.0

class DaltonParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.r1planner = BlindPlanner()

class TriangleParams(Params):
    def __init__(self):
        Params.__init__(self)
        self.r1planner = BlindPlanner()
        self.dist_break_distance = 300

def to_array(p):
    return np.array((p[0], p[1]))

def direct(p1, p2):
    p = p2 - p1
    return p / linalg.norm(p)

def rel_angle(p1, p2):
    return np.math.atan2(np.linalg.det([p1, p2]),np.dot(p1, p2))

def point_to_line_dist(line, p):
    (x1, y1), (x2, y2) = line
    x0, y0 = p
    d = -((y2 - y1)*x0 - (x2 - x1) *y0 +x2*y1 - y2*x1)/math.sqrt((y2 - y1)*(y2 - y1) + (x2 - x1)*(x2 - x1))
    return d

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    ret = np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    cross = np.cross(v1_u, v2_u)
    return ret if cross >= 0 else -ret


def get_intersect(a1, a2, b1, b2):
    """
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """
    s = np.vstack([a1,a2,b1,b2])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection
    if z == 0:                          # lines are parallel
        return (float('inf'), float('inf'))
    return (x/z, y/z)

def find_rad_from_two_points_and_tangent(p1, tang, p2):
    """
    p1 - point on the circle
    tang - direction vector relative to p1
    p2 - point on the circle
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = (p1 + p2) / 2.0
    p31 = (p1 - p3)
    a = np.angle(rel_angle(tang, p31))
    print ("rel_angle={}".format(a), file=sys.stderr)
    if a < 3.0 or a > 177.0 or a < -177.0:
        return 23000, (-1, -1)
    p4 = rotate(p31, degrees=90) + p3
    pend = rotate(tang, degrees=90)
    p5 = p1 + pend
    c = get_intersect((p1[0], p1[1]), (p5[0], p5[1]), (p3[0], p3[1]), (p4[0], p4[1]))
    dx = p1[0] - c[0]
    dy = p1[1] - c[1]
    return math.sqrt(dx*dx + dy*dy), c

def circ_rad(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    c = (x1-x2)**2 + (y1-y2)**2
    a = (x2-x3)**2 + (y2-y3)**2
    b = (x3-x1)**2 + (y3-y1)**2
    s = 2*(a*b + b*c + c*a) - (a*a + b*b + c*c)
    if s == 0.0:
        return (0, 0), 100000.
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
    print("ANGLE={}".format(angle_deg), file=sys.stderr)
    extra_fac = 0
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
        break_vec = p21 * Hub.params.break_dist
        target = t_pivot + direct + break_vec
        target0 = target + p21 * Hub.params.minimal_straight_dist
    #if math.fabs(angle_deg) < 51:
    #    extra_fac = 0.3
    fac1 = closest_point_to_segment(p1, p2, target) - (Hub.params.break_fac + extra_fac)
    fac0 = closest_point_to_segment(p1, p2, target0) - (Hub.params.break_fac + extra_fac)

    return (t_pivot[0], t_pivot[1]), ((target0[0], target0[1]), (target[0], target[1])), (fac0, fac1), angle_deg

class BreakBeforeTarget:
    def __init__(self):
        pass

    def get_thrust(self, thrust, dist, angle, target, rad):
        if dist < Hub.params.dist_break_distance:
            thrust = 50
        return thrust

class CurvatureThustStrategy:
    def __init__(self):
        self.thrust = 100
        pass

    def get_thrust(self, dist, angle, target, rad):
        self.thrust = int(round(rad/Hub.params.thrust_rad_100 * 100))
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
        self.state = 0
        self.stations = []
        self.thrust_stg = CurvatureThustStrategy()
        self.thrust, self.thrust2 = 0, 0
        self.rfacs = (0, 0)
        self.rad = 0
        self.regulator = PIDThrustRegulator()
        self.gt_regulator = GoToTargetRegulator()
        self.punch_mode = False
        self.breaker = BreakBeforeTarget()
        self.state = 0


    def plan(self, pos, dist, angle, target, stations):
        self.state = 0
        self.gt_regulator.reset()
        self.num_punch = Hub.params.num_punchs
        self.r100 = Hub.params.r100
        self.regulator.reset()
        self.punch_mode = False
        self.num_punch = Hub.params.num_punchs
        x1 = stations[-1][0]
        y1 = stations[-1][1]
        x2 = stations[0][0]
        y2 = stations[0][1]
        d = math.sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))
        if d > (2 * self.r100 + Hub.params.minimal_straight_dist):
            rad = self.r100
            print ("regular: dist={},rad={}".format(d, rad), file=sys.stderr)
        else:
            rad = (d - 600) / 2.0
            print ("special: dist={},rad={}".format(d, rad), file=sys.stderr)
        self.rad = rad
        print ("PLAN: {},{},{}".format(stations[-1], stations[0], stations[1]), file=sys.stderr)
        self.pivot, self.curve_st, self.rfacs, self.angle = calc_node_approach(stations[-1], stations[0], stations[1], rad)
        print ("RES: pivpt={},st={},facts={},angle={}".format(self.pivot, self.curve_st, self.rfacs, self.angle), file=sys.stderr)
        self.state = 0
        self.stations = stations
        #self.thrust2 = 100
        #self.thrust = 100
        self.max_alpha = 1.0 + 300.0 / d

    def act(self, pos, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        alpha = closest_point_to_segment(self.stations[-1], self.stations[0], pos)
        print ("alpha={}".format(alpha), file=sys.stderr)
        if dist > 1200 and self.num_punch > 0 and not self.punch_mode and Opponent.me.can_deliver_puch(Opponent.other, target, angle):
            print ("PUNCH!!!!!!", file=sys.stderr)
            self.num_punch -= 1
            if self.num_punch == 0:
                self.punch_mode = True
            self.gt_regulator.reset()
            npos, thrust, boost =  Opponent.other.next_pos(), 100, True
        else:
            if alpha < self.rfacs[0]:
                if self.state != 0:
                    self.gt_regulator.reset()
                self.stat = 0
                thrust = self.regulator.act(self.curve_st[0], angle)
                #return self.curve_st[0], self.thrust2, False
                print ("state0: moving to {}".format(self.curve_st[0]), file=sys.stderr)
                npos, thrust, boost = self.curve_st[0], thrust, False
                npos = self.gt_regulator.act(npos, pos)
            elif alpha < self.rfacs[1]:
                if self.state != 1:
                    self.gt_regulator.reset()
                self.state = 1
                thrust = self.regulator.act(self.curve_st[1], angle)
                print ("state1: moving to {}".format(self.curve_st[1]), file=sys.stderr)
                npos, thrust, boost = self.curve_st[1], thrust, False
                npos = self.gt_regulator.act(npos, pos)
            elif alpha <= self.max_alpha:
                if self.state != 2:
                    self.gt_regulator.reset()
                self.state = 2
                thrust = self.regulator.act(target, angle)
                print ("state2: moving to {}".format(target), file=sys.stderr)
                npos, thrust, boost = target, thrust, False
                npos = self.gt_regulator.act(npos, pos)
            else:
                #thrust = 30 if angle > 3 else 100
                thrust = self.regulator.act(target, angle)
                print ("over: moving to {}".format(target), file=sys.stderr)

                npos, thrust, boost = target, thrust, False
                npos = self.gt_regulator.act(npos, pos)
        thrust = self.breaker.get_thrust(thrust, dist, angle, target, rad)
        return npos, thrust, boost

class Clash:
    def __init__(self):
        self.dist = -1

    def plan(self, pos, dist, angle, target, stations):
        pass

    def act(self, pos, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        ox , oy = opponent_pos
        mx, my = pos
        self.dist = math.sqrt((ox - mx)*(ox - mx) + (oy - my)*(oy - my))
        print ("CLASH DIST={}".format(self.dist), file=sys.stderr)
        return opponent_pos, 100, True

class BlindPlanner:
    def __init__(self, br=True):
        self.is_chasing = False
        self.punch_mode = True
        self.regulator = PIDThrustRegulator()
        self.gt_regulator = GoToTargetRegulator()
        self.breaker = BreakBeforeTarget() if br else None
        self.aim = True
        pass

    def plan(self, pos, dist, angle, target, stations):
        self.aim = True
        self.regulator.reset()
        self.punch_mode = True
        self.gt_regulator.reset()
        pass

    def act(self, pos, dist, angle, target, last_pos, p_center, rad, opponent_pos):
        boost = False
        if math.fabs(angle) >= 90 and self.aim:
            return target, 0, False
        else:
            self.aim = False
        angle_ = math.fabs(angle)
        thrust = self.regulator.act(target, angle)
        tc = self.gt_regulator.act(target, pos)
        if not self.punch_mode and Opponent.me.can_deliver_puch(Opponent.other, target, angle):
            print ("PUNCH!!!!!!", file=sys.stderr)
            self.punch_mode = False
            return Opponent.other.next_pos(), 100, True
        print ("NOT PUNCHING {}".format(self.punch_mode), file=sys.stderr)
        if self.breaker != None:
            thrust = self.breaker.get_thrust(thrust, dist, angle, target, rad)
        return tc, thrust, boost

class GoToTargetRegulator:
    def __init__(self):
        self.reset()

    def reset(self):
        self.e = 0
        self.ie = 0
        self.last_e = 0
        self.dedt = 0

    def act(self, target, pos):
        dir_ = Opponent.me.dir
        if not Opponent.me.has_dir:
            return target
        pt = np.array(target) - np.array(pos)
        self.e = angle_between(pt, dir_)
        if math.fabs(self.e) > 85:
            self.reset()
            return target
        self.ie += self.e
        self.dedt = self.last_e - self.e
        Kp = Hub.params.gtKp
        Ki = Hub.params.gtKi
        Kd = Hub.params.gtKd
        c_angle = -(Kp*self.e + Ki*self.ie + Kd*self.dedt)
        self.last_e = self.e
        pt = rotate(pt, degrees=c_angle)
        res = np.array(pos) + pt
        return (res[0], res[1])


class PIDThrustRegulator:
    def __init__(self):
        self.reset()

    def reset(self):
        self.e = 0
        self.last_e = 0
        self.ie = 0
        self.dedt = 0
        self.thrust = 100

    def act(self, target, angle):
        self.e, r_a = self.error(target, angle)
        self.ie += self.e
        self.dedt = self.last_e - self.e
        thrust = self.thrust + math.fabs(math.sin(r_a))*(Hub.params.Kp*self.e + Hub.params.Ki*self.ie + Hub.params.Kd*self.dedt)
        self.last_e = self.e
        print("PID thrust={} e={}".format(thrust, self.e), file=sys.stderr)
        if thrust > 100.0:
            thrust = 100.0
        elif thrust < 0:
            thrust = 0
        self.thrust = int(thrust)
        return self.thrust

    def calc_radial_thrust(self, pos, c, target, angle):
        c = to_array(c)
        pos = to_array(pos)
        target = to_array(target)
        p_pos_target = direct(pos, target)
        p_head = rotate(p_pos_target, degrees=angle)
        p_pos_c = direct(pos, c)
        print ("p_head={},p_pos_c={}".format(p_head, p_pos_c), file=sys.stderr)
        rel_ang = rel_angle(p_head, p_pos_c)
        rel_ang_deg = np.degrees(rel_ang)
        radial_thrust = self.thrust * math.cos(rel_ang)
        print ("rel_ang={}, thrust={}, radial_thrust={}".format(rel_ang_deg, self.thrust, radial_thrust), file=sys.stderr)
        return rel_ang

    def error(self, target, angle):
        target = np.array(target)
        pos = np.array(Opponent.me.pos[-1])
        targ_dir = target - pos
        targ_dir = targ_dir / linalg.norm(targ_dir)
        my_dir = Opponent.me.dir
        if not Opponent.me.has_dir or len(Opponent.me.pos) < 3:
            error = 0
            r_a = 0
        else:
            pos0 = Opponent.me.pos[-3]
            pos1 = Opponent.me.pos[-2]
            pos2 = Opponent.me.pos[-1]
            #_, cur_rad = circ_rad( pos0, pos1, pos2)
            vel = linalg.norm(Opponent.me.vel)
            rad, c = find_rad_from_two_points_and_tangent((pos[0], pos[1]), (my_dir[0], my_dir[1]), target)
            vel_targ = Hub.params.vel_const * math.sqrt(rad)
            r_a = self.calc_radial_thrust(pos, c, target, angle)
            #dist = linalg.norm(pos - target)
            print("vel={} vel_targ={} rad={}".format(vel, vel_targ, rad), file=sys.stderr)
            #error = np.math.atan2(np.linalg.det([my_dir, targ_dir]),np.dot(my_dir, targ_dir))
            #if rad > 23000:
            #    rad = 23000
            #if cur_rad > 23000:
            #    cur_rad = 23000
            error = vel_targ - vel
        return error, r_a

class Opponent:
    me = None
    other = None

    def __init__(self):
        self.pos = []
        self.dir = None
        self.has_dir = False
        self.vel = None

    def report_pos(self, pos):
        self.pos.append(pos)
        if len(self.pos) > 10:
            self.pos = self.pos[1:]
        if len(self.pos) >= 2:
            self.vel = np.array(self.pos[-1]) - np.array(self.pos[-2])
            self.dir = self.vel / linalg.norm(self.vel)
            self.has_dir = True

    def dist(self, target):
        if len(self.pos) == 0:
            return 9999999999

        mx, my = self.pos[-1]
        tx, ty = target
        return math.sqrt((mx - tx)**2 + (my - ty)**2)

    def face_dir(self, angle, target):
        r_target = np.array(target) - np.array(self.pos[-1])
        r_target = r_target / linalg.norm(r_target)
        x_dir = np.array([1, 0])
        targ_ang = np.math.atan2(np.linalg.det([x_dir, r_target]),np.dot(x_dir, r_target))
        targ_ang = np.degrees(targ_ang)
        print ("targ_ang={}".format(targ_ang), file=sys.stderr)
        res = targ_ang - angle
        if res > 180:
            res = res - 360
        elif res < -180:
            res = 360 + res
        return res

    def can_deliver_puch(self, other, target, angle):
        d_to_other = self.dist(other.pos[-1])
        print ("DISTANCE={}".format(d_to_other), file=sys.stderr)
        if d_to_other > Hub.params.side_puch_distance:
                return False
        fc_dir = self.face_dir(angle, target)
        fdir = self.direction(other)
        rel_angle = math.fabs(fc_dir - fdir)
        print ("ABD ANGLE={} angle={} fc_dir={} fdir={}".format(rel_angle, angle, fc_dir, fdir), file=sys.stderr)
        if rel_angle <= Hub.params.side_punch_max_angle:
            print ("MAY PUNCH", file=sys.stderr)
            return True

        return False

    def next_pos(self):
        npos = np.array(self.pos[-1]) + np.array(self.vel)
        return npos[0], npos[1]

    def direction(self, other):
        p1 = np.array(self.pos[-1]) + np.array(self.vel)
        p2 = np.array(other.pos[-1]) + np.array(other.vel)
        p = p2 - p1
        x = np.array([1, 0])
        angle = np.math.atan2(np.linalg.det([x, p]),np.dot(x, p))
        angle_deg = np.degrees(angle)
        return angle_deg

    def is_chasing(self, other, target):
        if len(self.pos) < 2 or len(other.pos) < 2:
            return False

        if self.dist(other.pos[-1]) > Hub.params.chase_max:
            return False

        if self.dist(target) < other.dist(target):
            return False

        p1 = np.array(self.dir)
        p2 = np.array(other.dir)
        angle = np.math.atan2(np.linalg.det([p1, p2]),np.dot(p1, p2))
        angle_deg = np.degrees(angle)
        if math.fabs(angle_deg) > Hub.params.chase_max_angle:
            return False
        return True

class Arena:
    def __init__(self, name, stations, params=None):
        self.name = name
        self.stations = stations
        self.params = params

    def point_in_track(self, p):
        for ps in self.stations:
            if math.fabs(p[0] - ps[0]) < 70 and math.fabs(ps[1] - ps[1]) < 70:
                return True
        return False

    def opt_params(self):
        if self.params != None:
            Hub.params = self.params

class ArenaDetector:
    tracks = [Arena("hostile", [(13890, 1958), (8009, 3263), (2653, 7002), (10035, 5969)], HostileParams()),
              Arena("hostile2", [(9409, 7247), (5984, 4264), (14637, 1420), (3470, 7203)], HostileParams()),
              Arena("pyramid", [(7637, 5988), (3133, 7544), (9544, 4408), (14535, 7770), (6312, 4294), (7782, 851)]),
              Arena("triangle",[(6012, 5376), (11308, 2847), (7482, 6917)], TriangleParams()),
              Arena("dalton", [(9589, 1395), (3618, 4434), (8011, 7920), (13272, 5530)], DaltonParams()),
              Arena("makbilit", [(12942, 7222), (5655, 2587), (4107, 7431), (13475, 2323)], MacbilitParams()),
              Arena("arrow", [(10255, 4946), (6114, 2172), (3048, 5194), (6276, 7762), (14119, 7768), (13897, 1216)], ArrowParams()),
              Arena("Shosh",  [(9116, 1857), (5007, 5288), (11505, 6074)], ShoshParams()),
              Arena("Til",  [(10558, 5973), (3565, 5194), (13578, 7574), (12430, 1357)]),
              Arena("trapez",  [(11201, 5443), (7257, 6657), (5452, 2829), (10294, 3341)]),
              Arena("Mehumash", [(4049, 4630), (13054, 1928), (6582, 7823), (7494, 1330), (12701, 7080)], MehumashParams()),
              Arena("Trampoline",  [(3307, 7251), (14572, 7677), (10588, 5072), (13100, 2343), (4536, 2191), (7359, 4930)], TrampolineParams()),
              Arena("Zigzag", [(10660, 2295), (8695, 7469), (7225, 2174), (3596, 5288), (13862, 5092)],ZigzagParams())
            ]
    def __init__(self):
        self.detected_track = None
        pass

    def try_detect(self, stations):
        num_tracks = 0

        for track in ArenaDetector.tracks:
            detected = True
            for s in stations:
                if not track.point_in_track(s):
                    detected = False
                    break
            if detected:
                print ("Arena suspected: {}".format(track.name), file=sys.stderr)
                num_tracks += 1
                self.detected_track = track
        if num_tracks == 1:
            print ("Single Arena detected: {}".format(self.detected_track.name), file=sys.stderr)
            return self.detected_track

class Collector:
    def __init__(self):
        self.stations = []
        self.arena = None
        self.collecting = True
        self.last_target = (-1, -1)
        self.last_pos = []
        self.algo = Clash()
        self.clash = True
        self.thrustStrategies = []
        self.lap = 1
        self.stat_in_lap = 0
        self.opponent_pos = []
        self.arena_detector = ArenaDetector()

    def act(self, pos, dist, angle, target, opponent_pos):
        if self.clash and self.collecting:
            #if dist > 7000 and self.clash and self.collecting:
            self.clash = False
            self.algo = BlindPlanner()
        if self.last_target != target:
            self.stat_in_lap += 1
            print ("STATION={}".format(self.stat_in_lap), file=sys.stderr)

        if self.collecting and len(self.stations) and self.last_target != target and self.stations[0] == target:
            self.collecting = False
            self.lap = 1
            self.algo = Hub.params.r1planner
            print ("COLLECTED: {} {}".format(self.stations, self.algo), file=sys.stderr)

        if self.collecting and self.last_target != target:
            self.stations.append(target)
            if self.arena == None:
                self.arena = self.arena_detector.try_detect(self.stations)
                if self.arena is not None:
                    print ("Arena detected: {}".format(self.arena.name), file=sys.stderr)
                    self.arena.opt_params()
                    #self.algo = Hub.params.r1planner
                    #self.collecting = False
                    #self.stations = self.arena.stations
                else:
                    print("ARENA NOT DETECTED", file=sys.stderr)

        if not self.collecting and self.last_target != target and ((self.stat_in_lap % len(self.stations)) == 1):
            print ("LAP: {}".format(self.lap), file=sys.stderr)
            self.lap +=  1
        if self.lap == 3 and (self.stat_in_lap % len(self.stations) == 0) and self.last_target != target:
            print ("PHOTOFINISH..........", file=sys.stderr)
            self.algo = BlindPlanner(br=False)

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

class Hub:
    params = Params()

if __name__ == "__main__":
    calib = False
    lnx, lny = -1, -1
    n = 0
    hist = []
    opt = False
    algo = calib_circ(45) if calib else Collector()
    Opponent.me = Opponent()
    Opponent.other = Opponent()

    # game loop
    while True:
        # x: x position of your pod
        # y: y position of your pod
        # next_checkpoint_x: x position of the next check point
        # next_checkpoint_y: y position of the next check point

        inp = [int(i) for i in input().split()]
        i = -1

        x, y, nx, ny, dist, angle = inp
        Opponent.me.report_pos((x, y))

        ox, oy = [int(i) for i in input().split()]
        Opponent.other.report_pos((ox, oy))
        target, thrust, boost = algo.act((x, y), dist, angle, (nx, ny), (ox, oy))

        if boost:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " BOOST")
        else:
            print(str(int(target[0])) + " " + str(int(target[1]))+ " " + str(thrust))


