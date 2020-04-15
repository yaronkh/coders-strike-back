import sys
import time
import math
import numpy as np
from numpy import linalg
from scipy.special import comb
import copy


class PodParams:
    max_thrust = 200
    rot_vel = 18
    fric_thruts_to_thrust_ratio = 0.85
    def __init__(self):
        self.planner = Planner3()
        self.side_punch_max_angle = 30
        self.shell_punch_dist = 0
        self.num_punchs = 5
        self.rel_vel_for_crash_with_shield = 200
        self.num_of_turns_to_impact = 3
        self.pod_rad = 500
        self.brake_with_shield_vel = 300

class OpponentPodParams(PodParams):
    def __init__(self):
        self.planner = Planner3()
        self.side_punch_max_angle = 30
        self.collision_min_vel = 100
        self.num_punchs = 0

class OffencePodParams(PodParams):
    def __init__(self):
        PodParams.__init__(self)

class DefencePodParams(PodParams):
    def __init__(self):
        PodParams.__init__(self)
        self.side_punch_max_angle = 60
        self.shell_punch_dist = 0
        self.num_punchs = 0

class ArenaParams:
    station_rad = 600
    friction_fac = 0.85 # the current speed vector of each pod is multiplied by that factor
    station_tolerance = 450
    def __init__(self):
        self.r100 = 3600.
        self.minimal_straight_dist = 0.
        self.hard_drift_turns = 4
        self.super_hard_drift_turns = 0.0
        self.start_with_boost = False
        self.defender_dist_spare = 0
        self.break_dist = 0
        self.break_fac = 0.0
        self.gtKp = 0.9
        self.gtKi = 0
        self.gtKd = 0
        self.max_face_dir_error = 45

class ShoshParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)

class TrampolineParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)

class MacbilitParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r100 = 1800.

class TilParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)

class MehumashParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r100 = 3000.

class HostileParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r100 = 2200.


class ZigzagParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r100 = 2000.

class ArrowParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)

class PyramidParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.gtKp = 0.9

class TrapezParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.vel_const = 10.0
        self.super_hard_drift_turns = 1.5

class DaltonParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r1planner = Planner3()
        self.gtKp = 1.0
        self.start_with_boost = True

class TriangleParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)

def to_array(p):
    return np.array((p[0], p[1]))

def direct(p1, p2):
    p = p2 - p1
    return p / linalg.norm(p)

def dist_pnts(m, t):
    mx, my = m
    tx, ty = t
    return math.sqrt((mx - tx)**2 + (my - ty)**2)

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

def angleabs2angle(angle_abs, target, pos):
    """
    convert absolute heading angle to angle relative to target
    """
    angle_abs_rad = np.radians(angle_abs)
    angle_uv = (math.cos(angle_abs_rad), math.sin(angle_abs_rad))
    t = np.array(target) - np.array(pos)
    angle = angle_between(angle_uv, t)
    return angle


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
    if cross >= 0:
        return ret
    else:
        return -ret

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
    print ("rel angle {} {}".format(tang, p31), file=sys.stderr)
    a = 180. - np.degrees(rel_angle(tang, p31))
    print ("rel_angle={}".format(a), file=sys.stderr)
    if a < 3.0 or a > 177.0 or a < -177.0:
        return 23000, (-1, -1)
    p4 = rotate(p31, degrees=90) + p3
    pend = rotate(tang, degrees=90)
    p5 = p1 + pend
    c = get_intersect((p1[0], p1[1]), (p5[0], p5[1]), (p3[0], p3[1]), (p4[0], p4[1]))
    print ("center={}".format(c), file=sys.stderr)
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

def calc_node_approach(p1, p2, p3, rad, params):
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
        break_vec = p21 * params.break_dist
        target = t_pivot + direct + break_vec
        target0 = target + p21 * params.minimal_straight_dist
    #if math.fabs(angle_deg) < 51:
    #    extra_fac = 0.3
    fac1 = location_along_segment(p1, p2, target) - (params.break_fac + extra_fac)
    fac0 = location_along_segment(p1, p2, target0) - (params.break_fac + extra_fac)

    return (t_pivot[0], t_pivot[1]), ((target0[0], target0[1]), (target[0], target[1])), (fac0, fac1), angle_deg

def location_along_segment(p, q, x):
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

class PodKinematics:
    braking_dist_t100 = 0
    series_converge_fac = 1.0 / (1.0 - ArenaParams.friction_fac)

    @staticmethod
    def braking_dist(tangent_vel):
        return tangent_vel * PodKinematics.series_converge_fac

    @staticmethod
    def max_vel(thrust):
        return thrust * PodKinematics.series_converge_fac

PodKinematics.braking_dist_t100 = PodKinematics.braking_dist(PodKinematics.max_vel(PodParams.max_thrust))

class Simulator:
    def __init__(self):
        pass


    def calc_next_turn(self, tracker, target, thrust, boost, shield):
        pface = self.calc_next_rotation(tracker, target)
        thrust_vec = pface * thrust
        new_vel = np.array(tracker.vel) * ArenaParams.friction_fac + thrust_vec
        npos = np.array(tracker.pos[-1]) + new_vel
        new_vel = np.trunc(new_vel)
        new_angle = np.math.atan2(pface[1], pface[1])
        return (npos[0], npos[1]), (new_vel[0], new_vel[1]), new_angle

    def calc_next_rotation(self, tracker, target):
        p12 = unit_vector(np.array(target) - np.array(tracker.pos[-1]))
        angle = np.radians(tracker.angle_abs)
        pface = (math.cos(angle), math.sin(angle))
        angle = angle_between(pface, (p12[0], p12[1]))
        if math.fabs(angle) <=PodParams.rot_vel:
            return  (p12[0], p12[1])
        d =  PodParams.rot_vel if angle > 0 else -PodParams.rot_vel
        return rotate(np.array(pface), d)

class Planner3:
    def __init__(self):
        self.tracker = None
        self.gt_regulator = GoToTargetRegulator()

    def plan(self, pos, angle_abs, target, stations, tracker, params):
        angle = angleabs2angle(angle_abs, target, pos)
        self.tracker = tracker
        self.gt_regulator.reset(ArenaParams.station_tolerance, tracker, params)

    def act(self, pos,  angle_abs, target):
        angle = angleabs2angle(angle_abs, target, pos)
        tc, thrust = self.gt_regulator.act(target, pos, angle)
        if self.gt_regulator.is_pointing:
            v1 = linalg.norm(self.tracker.vel)
            if v1 > 10:
                s = v1 /(1 - ArenaParams.friction_fac)
                dist = dist_pnts(pos, target)
                if dist < (s / 2.0):
                    return self.tracker.stations[1], 0, False, False
        return tc, thrust, False, False

class GuardPostStrategy:
    """
    A strategy where the pod rush to a given position
    and does not let the opponent leader to reach it.
    The strategy ends when the opponent reach that target
    """
    def __init__(self, traker):
        self.tracker = traker
        self.target = (0, 0)
        self.algo = BlindPlanner()
        self.station = []
        self.in_position = False
        self.rot_to_new_pos = False

    def plan(self, pos, angle, target, stations, tracker, params
            ):
        """
        stations - the list of stations, the post to defend is the first station
        """
        self.stations = stations[:]
        self.tracker = tracker
        self.thrust = PodParams.max_thrust
        self.target = target
        self.in_position = False
        self.rot_to_new_pos = False
        self.algo.plan(pos, angle, target, stations, tracker, params)
        self.tolerance = ArenaParams.station_rad - self.tracker.pod_params.pod_rad

    def act(self, pos, angle_abs, target):
        rel = np.array(target) - np.array(pos)
        dist = linalg.norm(rel)
        targ_dir = rel / dist
        vel_targ = np.dot(self.tracker.vel, targ_dir)
        brake_dist = PodKinematics.braking_dist(vel_targ)
        if dist < self.tolerance or brake_dist >= dist:
            self.in_position = True
            return opponent_leader.pos[-1], 0, False, False
        elif self.in_position:
            angle_abs_rad = np.radians(angle_abs)
            h = (math.cos(angle_abs_rad), math.sin(angle_abs_rad))
            angle_to_new = math.fabs(angle_between(h, targ_dir))
            print ("ANGLE ADJUST = {}".format(angle_to_new), file=sys.stderr)
            if angle_to_new > 45:
                return target, 0, False, False
            return self.algo.act(pos, angle_abs, target)
        else:
            return self.algo.act(pos, angle_abs, target)

class BlindPlanner:
    def __init__(self, br=True):
        self.is_chasing = False
        self.gt_regulator = GoToTargetRegulator()
        self.aim = True
        self.tracker = None
        pass

    def plan(self, pos, angle, target, stations, tracker, params):
        self.tracker = tracker
        self.aim = True
        self.gt_regulator.reset(ArenaParams.station_tolerance, tracker, params)

    def act(self, pos, angle_abs, target):
        angle = angleabs2angle(angle_abs, target, pos)
        boost = False
        tc, thrust = self.gt_regulator.act(target, pos, angle)
        print ("START WITH BOOST {} {}".format(self.tracker.arena_params.start_with_boost, self.tracker.boost), file=sys.stderr)
        if self.tracker.arena_params.start_with_boost and  not self.tracker.boost:
            boost = True
        return tc, thrust, boost, False

class GoToTargetRegulator:
    def __init__(self):
        self.reset()

    def reset(self, tolerance=0, tracker=None, params=None):
        self.e = 0
        self.ie = 0
        self.last_e = 0
        self.dedt = 0
        self.params = params
        self.tracker = tracker
        self.tolerance = tolerance
        self.is_pointing = False

    def act(self, target, pos, angle):
        if self.is_pointing_target(target):
            self.is_pointing = True
            ret = (self.t[0], self.t[1]), PodParams.max_thrust
            print ("POINTED TO TARGET", file=sys.stderr)
            return ret
        self.is_pointing = False
        dir_ = self.tracker.direction()
        pt = np.array(target) - np.array(pos)
        if linalg.norm(self.tracker.vel) < 1.0:
            thrust = PodParams.max_thrust
        else:
            thrust = self.regulate_thrust(target, angle)
        if dir_[0] != 0. and dir_[1] != 0.:
            self.e = angle_between(pt, dir_)
            print ("angle between is {}".format(self.e), file=sys.stderr)
        else:
            self.e = 0
        if math.fabs(self.e) > 85:
            v = linalg.norm(self.tracker.vel)
            if v > 200:
                self.reset(self.tolerance, self.tracker, self.params)
                #thrust = 11
                return target, thrust
        self.ie += self.e
        self.dedt = self.last_e - self.e
        Kp = self.params.gtKp
        Ki = self.params.gtKi
        Kd = self.params.gtKd
        c_angle = -(Kp*self.e + Ki*self.ie + Kd*self.dedt)
        print ("c_angle={}".format(c_angle), file=sys.stderr)
        #if c_angle > 89:
        #    c_angle = 89
        #    thrust = 20
        #elif c_angle < -89:
        #    c_angle = 89
        #    thrust = 20
        self.last_e = self.e
        pt = rotate(pt, degrees=c_angle)
        res = np.array(pos) + pt

        return (res[0], res[1]), thrust

    def calc_target_vel(self, target):
        p1 = self.tracker.pos[-1]
        tg = unit_vector(self.tracker.vel)
        rad, _ = find_rad_from_two_points_and_tangent(p1, tg, target)
        rad += 600
        alphad = math.fabs(self.tracker.pod_deflection())
        alpha = math.fabs(np.radians(alphad))
        vt = rad * ( 1 - ArenaParams.friction_fac) * math.tan(alpha)
        #print ("rad={} alpha={},VT={}".format(rad, alphad, vt), file=sys.stderr)
        # verify that the pod is fast enough to handle that radius and velocity
        omega = np.radians(PodParams.rot_vel)
        vt2 = rad * omega
        if vt2 < vt:
            print ("reducing velocity to meet max radial vel", file=sys.stderr)
            vt = vt2
        return vt, alpha

    def is_pointing_target(self, target):
        if linalg.norm(self.tracker.vel) < 2:
            print ("VELOCITY TOO SMALL TO TELL POINTING", file=sys.stderr)
            return False
        direc = unit_vector(np.array(self.tracker.vel))
        p1 = np.array(self.tracker.pos[-1])
        p2 = p1 + direc
        perp = np.array(rotate(direc, degrees=90))
        p3 = np.array(target)
        p4 = p3 + perp
        self.t = get_intersect(p1, p2, p3, p4)
        fac = location_along_segment(p1, p2, self.t)
        d = dist_pnts(self.t, target)
        print ("distance from target={} targ={} pinpoint={} fac={}".format(d, target, self.t, fac), file=sys.stderr)
        if d < self.tolerance and fac >= 0:
            return True
        return False

    def regulate_thrust(self, target, angle):
        pod_def = self.tracker.pod_deflection()
        print ("pod deflection={}".format(pod_def), file=sys.stderr)
        if angle < 1.0 and math.fabs(self.tracker.pod_deflection()) < 1:
            v_tar, alpha = 99999999, 0
        else:
            v_tar, alpha = self.calc_target_vel(target)
        thrust = v_tar*(1 - ArenaParams.friction_fac) / math.cos(alpha)
        print ("v_tar={} alpha={} thrust={}".format(v_tar, alpha, thrust), file=sys.stderr)
        thrust /= PodParams.fric_thruts_to_thrust_ratio
        if thrust > PodParams.max_thrust:
            thrust = PodParams.max_thrust
        elif thrust < 0:
            thrust = 0
        self.thrust = int(thrust)
        return self.thrust

class Tracker:
    me = None
    other = None

    def __init__(self, num_laps, me):
        self.num_laps = 0
        self.arena_params = None
        self.pod_params = None
        self.pos = []
        self.vel = None
        self.angle = 0
        self.next_cp = 0
        self.prev_cp = 0
        self.stations = []
        self.time_left_to_punch = 0
        self.num_laps = num_laps
        self.passed_stations = -1
        self.me = me # type: boolean
        self.boost = False
        self.boost_turns = 0

    def act(self):
        print ("VEL {}".format(linalg.norm(self.vel)), file=sys.stderr)
        if self.prev_cp != self.next_cp:
            print ("TRACHER:ACT now planning", file=sys.stderr)
            self.pod_params.planner.plan(self.pos[-1], self.angle, self.stations[0], self.stations, self, self.arena_params)

        shield = self.need_protection()
        if self.time_left_to_punch > 0:
            self.time_left_to_punch -= 1
        elif not shield and self.attempt_punch():
            return
        tc, thrust, boost, shield2 = self.pod_params.planner.act(self.pos[-1], self.angle, self.stations[0])
        shield |= shield2
        self.write(tc, thrust, boost, shield)

    def write(self, tc, thrust, boost, shield):
        if not self.boost and boost:
            self.boost_turns = 4
            self.boost = boost
        if self.boost_turns > 0:
            self.boost_turns -= 1
        if shield:
            print(str(int(tc[0])) + " " + str(int(tc[1]))+ " SHIELD")
        elif boost:
            print(str(int(tc[0])) + " " + str(int(tc[1]))+ " BOOST")
        else:
            print(str(int(tc[0])) + " " + str(int(tc[1]))+ " " + str(thrust))

    def direction(self):
        vel_abs = linalg.norm(self.vel)
        if vel_abs == 0:
            vel_abs = 1.
        return self.vel[0] / vel_abs, self.vel[1] / vel_abs

    def report_pos(self, pos, vel, angle, next_cp):
        self.pos.append(pos)
        if len(self.pos) > 10:
            self.pos = self.pos[1:]
        self.angle = angle
        self.vel = vel
        if next_cp != self.next_cp:
            self.passed_stations += 1
            self.stations = arena_detector.stations[next_cp:] + arena_detector.stations[:next_cp]
        self.prev_cp = self.next_cp
        self.next_cp = next_cp

    @staticmethod
    def leader(p1, p2):
        if p1.passed_stations < p2.passed_stations:
            return p2, p1
        elif p1.passed_stations > p2.passed_stations:
            return p1, p2
        else:
            fac_p1 = math.fabs(1.0 - location_along_segment(p1.stations[-1], p1.stations[0], p1.pos[-1]))
            fac_p2 = math.fabs(1.0 - location_along_segment(p2.stations[-1], p2.stations[0], p2.pos[-1]))
            if fac_p1 < fac_p2:
                return p1, p2
            else:
                return p2, p1

    def configure(self, arena_detector, num_laps):
        self.num_laps  = num_laps
        self.arena_params = copy.deepcopy(arena_detector.detected_track.params) if arena_detector.detected_track != None  and arena_detector.detected_track.params != None else ArenaParams()
        print ("ARENA PARAMS={}".format(self.arena_params), file=sys.stderr)
        self.stations = arena_detector.stations[:]

    def dist(self, target):
        if len(self.pos) == 0:
            return 9999999999

        mx, my = self.pos[-1]
        tx, ty = target
        return math.sqrt((mx - tx)**2 + (my - ty)**2)

    def direction_rel_to(self, other):
        p1 = np.array(self.pos[-1]) + np.array(self.vel)
        p2 = np.array(other.pos[-1]) + np.array(other.vel)
        p = p2 - p1
        x = np.array([1, 0])
        angle = np.math.atan2(np.linalg.det([x, p]),np.dot(x, p))
        angle_deg = np.degrees(angle)
        return angle_deg

    def rel_vel(self, other):
        v = np.array(self.vel) - np.array(other.vel)
        return linalg.norm(v)

    def pod_deflection(self):
        al = np.radians(self.angle)
        face = (math.cos(al), math.sin(al))

        print ("ANGLE={} face={} vel={}".format(self.angle, face, self.vel), file=sys.stderr)
        return angle_between(self.vel, face)

    def going_to_collide(self, other, num_turns):
        rvel = np.array(other.vel) - np.array(self.vel)
        rpos = np.array(self.pos[-1]) - np.array(other.pos[-1])
        rpos_u = unit_vector(rpos)
        angle = angle_between(rvel, rpos)
        angle = math.fabs(angle)
        if angle >= 90:
            return False
        vc = np.dot(rvel, rpos_u)
        d = linalg.norm(rpos) - 2 * self.pod_params.pod_rad
        print ("GOING_TO_COLLIDE d={} vc={} num_turns={}".format(d, vc, num_turns), file=sys.stderr)
        return (d / vc) <= num_turns

    def attempt_punch(self):
        if len(self.pos) < 2:
            return False
        for other in Tracker.other:
            if self.can_deliver_puch(other):
                rvel = self.rel_vel(other)
                shell = self.going_to_collide(other, 1.0) and (rvel > self.pod_params.rel_vel_for_crash_with_shield or self.boost_turns > 0)
                print ("PUNCH {} shell={} boost_turns={}!!!!!!".format(rvel, shell, self.boost_turns), file=sys.stderr)
                self.write(other.next_pos(), PodParams.max_thrust, True, shell)
                self.time_left_to_punch = self.pod_params.num_punchs
                return True
        return False

    def need_protection(self):
        for other in Tracker.other:
            if len(other.pos) < 2:
                continue
            rvel = self.rel_vel(other)
            if rvel <= self.pod_params.rel_vel_for_crash_with_shield:
                continue
            if self.going_to_collide(other, 1.0):
                return True
        return False


    def can_deliver_puch(self, other):
        if not self.going_to_collide(other, self.pod_params.num_of_turns_to_impact):
            return False
        fc_dir = np.degrees(math.atan2(self.vel[1], self.vel[0]))
        fdir = self.direction_rel_to(other)
        rel_angle = math.fabs(fc_dir - fdir)
        if rel_angle <= self.pod_params.side_punch_max_angle:
            print ("MAY PUNCH", file=sys.stderr)
            return True

        return False

    def next_pos(self):
        """
        try to estimate the position in the next loop
        """
        if len(self.pos) < 2:
            return None
        dir_ = unit_vector(np.array(self.vel))
        if linalg.norm(self.vel) < 1.0:
            return self.pos[-1]
        rel = np.array(self.pos[-1]) - np.array(self.pos[-2])
        dir_ *= linalg.norm(rel)
        ret = np.array(self.pos[-1]) + dir_
        return ret[0], ret[1]

class Defender(Tracker):
    def __init__(self, num_laps, me):
        Tracker.__init__(self, num_laps, me)
        self.algo = GuardPostStrategy(num_laps)
        self.stations = []
        self.station_id = -1
        self.leader_origin = (0, 0)
        self.leader_target = (0, 0)
        self.leader_next = (0, 0)

    def calc_target_full_block(self, leader):
        p23 = np.array(leader.pos[-1]) - np.array(self.leader_target)
        p23 = p23 / linalg.norm(p23) * 800
        target = np.array(self.leader_target) + p23
        return target

    def calc_safe_target_pos(self):
        p21 = unit_vector(np.array(self.leader_origin) - np.array(self.leader_target))
        p23 = unit_vector(np.array(self.leader_next) - np.array(self.leader_target))
        p_dir = (p21 + p23) * 600 + np.array(self.leader_target)
        return (p_dir[0], p_dir[1])

    def is_self_blocking(self, leader):
        if all_leader.me == False:
            return False

        print ("ME STATIONS:{} {}".format(all_leader.passed_stations, self.station_id), file=sys.stderr)
        if all_leader.passed_stations >= self.station_id:
            return False

        return True

    def can_reach_target_before_leader(self, leader):
        leader_num_stations = self.station_id - leader.passed_stations
        if leader_num_stations <= 0:
            return False
        if self.algo.in_position:
            return True
        pnts0 = [leader.pos[-1]] + leader.stations[:leader_num_stations - 1]
        pnts1 = leader.stations[:leader_num_stations]
        leader_dist = 0.0
        for p1, p2 in zip(pnts0, pnts1):
            leader_dist += dist_pnts(p1, p2)

        my_dist = self.dist(self.leader_target)
        print ("MY_DIST = {}, leader_DIST={}".format(my_dist, leader_dist), file=sys.stderr)
        return (my_dist + self.arena_params.defender_dist_spare) <= leader_dist

    def act(self):
        leader = opponent_leader
        #print ("CAN_REACH_LEADER={}".format(self.can_reach_target_before_leader(leader)), file=sys.stderr)
        jump = 0
        while not self.can_reach_target_before_leader(leader):
            self.leader_origin = leader.stations[jump - 1]
            p23 = np.array(self.leader_origin) - np.array(leader.stations[jump])
            self.leader_target = leader.stations[jump]
            print ("SETTING TARGET TO {}".format(self.leader_target), file=sys.stderr)
            jn = (jump + 1) % len(leader.stations)
            self.leader_next = leader.stations[jn]
            p23 = p23 * 0.2
            target = leader.stations[1] + p23
            self.station_id = leader.passed_stations + jump + 1
            self.stations = leader.stations[jn:] +leader.stations[: jn]
            self.algo.plan(self.pos[-1], self.angle, target, self.stations, self, self.arena_params)
            jump += 1
        #print ("CAN_REACH_LEADER={}".format(self.can_reach_target_before_leader(leader)), file=sys.stderr)

        if self.is_self_blocking(leader):
            target = self.calc_safe_target_pos()
        else:
            target = self.calc_target_full_block(leader)

        shield = self.need_protection()
        print ("defender shield={}".format(shield), file=sys.stderr)

        if self.time_left_to_punch > 0:
            self.time_left_to_punch -= 1
        elif not shield and self.attempt_punch():
            return

        tc, thrust, boost, shield2 = self.algo.act(self.pos[-1], self.angle, target)
        shield |= shield2
        self.write(tc, thrust, boost, shield)

    def report_pos(self, pos, vel, angle, next_cp):
        self.pos.append(pos)
        if len(self.pos) > 10:
            self.pos = self.pos[1:]
        self.angle = angle
        self.vel = vel

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

class ArenaDetector:
    tracks = [Arena("hostile", [(13890, 1958), (8009, 3263), (2653, 7002), (10035, 5969)], HostileParams()),
              Arena("hostile2", [(9409, 7247), (5984, 4264), (14637, 1420), (3470, 7203)], HostileParams()),
              Arena("pyramid", [(7637, 5988), (3133, 7544), (9544, 4408), (14535, 7770), (6312, 4294), (7782, 851)], PyramidParams()),
              Arena("triangle",[(6012, 5376), (11308, 2847), (7482, 6917)], TriangleParams()),
              Arena("dalton", [(9589, 1395), (3618, 4434), (8011, 7920), (13272, 5530)], DaltonParams()),
              Arena("makbilit", [(12942, 7222), (5655, 2587), (4107, 7431), (13475, 2323)], MacbilitParams()),
              Arena("arrow", [(10255, 4946), (6114, 2172), (3048, 5194), (6276, 7762), (14119, 7768), (13897, 1216)], ArrowParams()),
              Arena("Shosh",  [(9116, 1857), (5007, 5288), (11505, 6074)], ShoshParams()),
              Arena("Til",  [(10558, 5973), (3565, 5194), (13578, 7574), (12430, 1357)], TilParams()),
              Arena("trapez",  [(11201, 5443), (7257, 6657), (5452, 2829), (10294, 3341)], TrapezParams()),
              Arena("Mehumash", [(4049, 4630), (13054, 1928), (6582, 7823), (7494, 1330), (12701, 7080)], MehumashParams()),
              Arena("Trampoline",  [(3307, 7251), (14572, 7677), (10588, 5072), (13100, 2343), (4536, 2191), (7359, 4930)], TrampolineParams()),
              Arena("Zigzag", [(10660, 2295), (8695, 7469), (7225, 2174), (3596, 5288), (13862, 5092)],ZigzagParams())
            ]
    def __init__(self):
        self.detected_track = None
        self.stations = []
        pass

    def try_detect(self, stations):
        self.stations = stations
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

arena_detector = ArenaDetector()
defence_params = DefencePodParams()
offence_params = OffencePodParams()

opponent_leader = None
opponent_follower = None
all_leader = None
all_second = None

if __name__ == "__main__":
    num_laps = int(input())
    Tracker.me = Tracker(num_laps, me=True), Defender(num_laps, me=True)
    Tracker.other = Tracker(num_laps, me=False), Tracker(num_laps, me=False)
    Tracker.me[0].pod_params = offence_params
    Tracker.me[1].pod_params = defence_params
    Tracker.other[0].pod_params = OpponentPodParams()
    Tracker.other[1].pod_params = OpponentPodParams()

    check_point_count = int(input())
    stations = []
    print ("num_check_points={}".format(check_point_count), file=sys.stderr)
    for _ in range(check_point_count):
        xs, ys = input().split()
        stations.append((int(xs), int(ys)))

    arena_detector.try_detect(stations)
    all_trackers = (Tracker.me[0], Tracker.me[1], Tracker.other[0], Tracker.other[1])

    Tracker.me[0].configure(arena_detector, num_laps)
    Tracker.me[1].configure(arena_detector, num_laps)

    # game loop
    while True:
        for tracker in all_trackers:
            s = input()
            print ("========={}".format(s), file=sys.stderr)
            x, y, vx, vy, angle, next_cp = [int(i) for i in s.split()]
            tracker.report_pos((x, y), (vx, vy), angle, next_cp)

        opponent_leader, opponent_follower = Tracker.leader(Tracker.other[0], Tracker.other[1])
        all_leader, all_second = Tracker.leader(Tracker.me[0], opponent_leader)

        print ("OFFENCE DATA", file=sys.stderr)
        Tracker.me[0].act()
        print ("DEFENCE DATA", file=sys.stderr)
        Tracker.me[1].act()
