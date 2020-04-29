#!/usr/bin/python3
import sys
import time
import math
import numpy as np
from numpy import linalg
from scipy.special import comb
import copy
import random
import struct


class PodParams:
    max_thrust = 200
    rot_vel = 18
    fric_thruts_to_thrust_ratio = 0.85
    def __init__(self):
        #self.planner                       = Genetic()
        self.planner                       = Planner3()
        self.side_punch_max_angle          = 30
        self.shell_punch_dist              = 0
        self.num_punchs                    = 5
        self.rel_vel_for_crash_with_shield = 300
        self.num_of_turns_to_impact        = 3
        self.pod_rad                       = 350
        self.brake_with_shield_vel         = 300

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
        self.pursuit_look_ahead_turns = 3

class ArenaParams:
    station_rad = 600
    friction_fac = 0.85 # the current speed vector of each pod is multiplied by that factor
    station_tolerance = 450
    def __init__(self):
        self.r100                   = 3600.
        self.minimal_straight_dist  = 0.
        self.hard_drift_turns       = 4
        self.super_hard_drift_turns = 0.0
        self.start_with_boost       = False
        self.defender_dist_spare    = 1000
        self.break_dist             = 0
        self.break_fac              = 0.0
        self.gtKp                   = 1.0
        self.gtKi                   = 0
        self.gtKd                   = 0
        self.max_face_dir_error     = 45
        self.mitbader_vel           = 130

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
        self.gtKp = 1.0

class TrapezParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.vel_const = 10.0
        self.super_hard_drift_turns = 1.5
        self.mitbader_vel           = 100

class DaltonParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.r1planner = Planner3()
        self.gtKp = 0.5
        self.start_with_boost = True

class TriangleParams(ArenaParams):
    def __init__(self):
        ArenaParams.__init__(self)
        self.mitbader_vel           = 100
        #self.planner = BlindPlanner()

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
    #print ("rel angle {} {}".format(tang, p31), file=sys.stderr)
    a = 180. - np.degrees(rel_angle(tang, p31))
    #print ("rel_angle={}".format(a), file=sys.stderr)
    if a < 3.0 or a > 177.0 or a < -177.0:
        return 23000, (-1, -1)
    p4 = rotate(p31, degrees=90) + p3
    pend = rotate(tang, degrees=90)
    p5 = p1 + pend
    c = get_intersect((p1[0], p1[1]), (p5[0], p5[1]), (p3[0], p3[1]), (p4[0], p4[1]))
    #print ("center={}".format(c), file=sys.stderr)
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
    #print("ANGLE={}".format(angle_deg), file=sys.stderr)
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

    def calc_next_turn(self, target, pos, vel, angle, thrust, boost, shield):
        #Rotation: the pod rotates to face the target point, with a maximum of 18 degrees (except for the 1rst round).
        pface = self.calc_next_rotation(pos, angle, target)

        #Acceleration: the pod's facing vector is multiplied by the given thrust value. The result is added to the current speed vector.
        thrust_vec = pface * thrust
        new_vel = (np.array(vel) + thrust_vec)#* ArenaParams.friction_fac

        #Movement: The speed vector is added to the position of the pod. If a collision would occur at this point, the pods rebound off each other.
        npos = np.array(pos) + new_vel

        #Friction: the current speed vector of each pod is multiplied by 0.85
        new_vel *= ArenaParams.friction_fac
        new_vel = np.trunc(new_vel)

        #The speed's values are truncated and the position's values are rounded to the nearest integer.
        new_angle = np.degrees(np.math.atan2(pface[1], pface[0]))
        if new_angle < 0:
            new_angle += 360
        return (int(npos[0]), int(npos[1])), (int(new_vel[0]), int(new_vel[1])), int(new_angle)

    def calc_next_rotation(self, pos, angle, target):
        p12 = unit_vector(np.array(target) - np.array(pos))
        angle = np.radians(angle)
        pface = (math.cos(angle), math.sin(angle))
        angle = angle_between(pface, (p12[0], p12[1]))
        if math.fabs(angle) <= PodParams.rot_vel:
            return  np.array((p12[0], p12[1]))
        d =  PodParams.rot_vel if angle > 0 else -PodParams.rot_vel
        r =rotate(np.array(pface), degrees=d)
        return rotate(np.array(pface), degrees=d)

class Simulator_:
    def __init__(self):
        pass

    def calc_next_turn(self, target, pos, vel, angle, thrust, boost, shield):
        #Rotation: the pod rotates to face the target point, with a maximum of 18 degrees (except for the 1rst round).
        pface = self.calc_next_rotation(pos, angle, target)

        #Acceleration: the pod's facing vector is multiplied by the given thrust value. The result is added to the current speed vector.
        thrust_vec = pface * thrust
        new_vel = (np.array(vel) + thrust_vec)#* ArenaParams.friction_fac

        #Movement: The speed vector is added to the position of the pod. If a collision would occur at this point, the pods rebound off each other.
        npos = np.array(pos) + new_vel

        #Friction: the current speed vector of each pod is multiplied by 0.85
        new_vel *= ArenaParams.friction_fac
        new_vel = np.trunc(new_vel)

        #The speed's values are truncated and the position's values are rounded to the nearest integer.
        new_angle = np.degrees(np.math.atan2(pface[1], pface[0]))
        if new_angle < 0:
            new_angle += 360
        return (npos[0], npos[1]), (new_vel[0], new_vel[1]), new_angle

    def calc_next_rotation(self, pos, angle, target):
        p12 = unit_vector(np.array(target) - np.array(pos))
        angle = np.radians(angle)
        pface = (math.cos(angle), math.sin(angle))
        angle = angle_between(pface, (p12[0], p12[1]))
        if math.fabs(angle) <= PodParams.rot_vel:
            return  np.array((p12[0], p12[1]))
        d =  PodParams.rot_vel if angle > 0 else -PodParams.rot_vel
        return rotate(np.array(pface), degrees=d)

#Genetic algoriths stuff starts here
class Genetic:
    NUM_OPT_COMMANDS = 6
    NUM_BITS_PER_COMMAND = 8
    CMD_RES_DIVIDER = (1 << NUM_BITS_PER_COMMAND)
    SINGLE_CMD_BITS = NUM_OPT_COMMANDS * NUM_BITS_PER_COMMAND
    THRUST_RESOLUTION = PodParams.max_thrust / (CMD_RES_DIVIDER - 1)
    ANGLE_RES = 360 / (CMD_RES_DIVIDER - 1)
    FlT = [False, True]
    GENERATION_SIZE = 8

    def __init__(self):
        self.simulator = Simulator()
        self.tracker = None
        self.hit_dist = 0.

    @staticmethod
    def chromo_encode(thrust_cmds, angle_cmds):
        res = [False] * Genetic.SINGLE_CMD_BITS * 2
        for i in range(Genetic.NUM_OPT_COMMANDS):
            t = int(np.around(thrust_cmds[i] / Genetic.THRUST_RESOLUTION))
            ang = int(np.around(angle_cmds[i] / Genetic.ANGLE_RES))
            for b in range(Genetic.NUM_BITS_PER_COMMAND):
                offs = (i + 1) * Genetic.NUM_BITS_PER_COMMAND - b - 1
                res[offs] = Genetic.FlT[t & 0x1]
                res[Genetic.SINGLE_CMD_BITS + offs] = Genetic.FlT[ang & 0x1]
                t = (t >> 1)
                ang = (ang >> 1)
            assert t == 0
            assert ang == 0
        return res

    @staticmethod
    def chromo_decode(chromo):
        thrust_cmds = [False] * Genetic.NUM_OPT_COMMANDS
        angle_cmds = thrust_cmds[:]
        for i in range(Genetic.NUM_OPT_COMMANDS):
            t = 0
            ang = 0
            for b in range(Genetic.NUM_BITS_PER_COMMAND):
                t = (t << 1)
                ang = (ang << 1)
                offs = i * Genetic.NUM_BITS_PER_COMMAND + b
                if chromo[offs]:
                    t |= 1
                if chromo[Genetic.SINGLE_CMD_BITS + offs]:
                    ang |= 1
            thrust_cmds[i] = int(np.around(t * Genetic.THRUST_RESOLUTION))
            angle_cmds[i] = int(np.around(ang * Genetic.ANGLE_RES))
        return thrust_cmds, angle_cmds

    def simu(self, thurst_v, angle_v):
        pos = []
        n_pos = self.tracker.pos[-1]
        n_vel = self.tracker.vel
        n_angle = self.tracker.angle
        for t, a in zip(thurst_v, angle_v):
            a_rad = np.radians(a)
            target = (math.cos(a_rad) * 30000, math.sin(a_rad) * 30000)
            n_pos, n_vel, n_angle = self.simulator.calc_next_turn(target, n_pos, n_vel, n_angle, t, False, False)
            pos.append(n_pos)
        return pos

    def target(self, chromo):
        thrust, angle = Genetic.chromo_decode(chromo)
        #print ("THRUST={},ANG={}".format(thrust, angle), file=sys.stderr)
        pos_v = self.simu(thrust, angle)
        grade0 = 0.
        grade1 = 0.

        for p in pos_v:
            d = dist_pnts(p, self.tracker.stations[0])
            if d <= self.hit_dist:
                grade0 = 1.0
        if grade0 > 0:
            d = dist_pnts(pos_v[-1], self.tracker.stations[1])
            if d <= self.hit_dist:
                grade1 = 1.0
            else:
                grade1 = self.hit_dist / d
        else:
            d = dist_pnts(pos_v[-1], self.tracker.stations[0])
            if d <= self.hit_dist:
                grade0 = 1.0
            else:
                grade0 = self.hit_dist / d
        return grade0 * 0.7 + grade1 * 0.3

    def guess_initial_set(self, num_guesses):
        res = []
        tot_guesses = []
        for _ in range(num_guesses):
            thrust_v = np.random.randint(0, high=PodParams.max_thrust, size=Genetic.NUM_OPT_COMMANDS)
            ang_v = np.random.randint(0, 360, size=Genetic.NUM_OPT_COMMANDS)
            tot_guesses.append((thrust_v, ang_v))
            res.append(Genetic.chromo_encode(thrust_v, ang_v))
        #print ("Guesses={}".format(tot_guesses), file=sys.stderr)
        return res

    def fitness(self, data):
        #for g in data:
        #    t = self.target(g)
        #    #print ("g={} fit={}".format(g, t), file=sys.stderr)
        return [self.target(g) for g in data]

    def natural_selection(self, FP):
        s = 0
        FPT = [0.0]
        for p in FP:
            FPT.append(p + s)
            s += p
        FPT.append(1.0)
        Pr = sorted(np.random.uniform(0.0, 1.0, len(FP)))
        iFP = 0
        while len(Pr):
            if (Pr[0] >= FPT[iFP] and Pr[0] < FPT[iFP + 1]):
                yield(iFP)
                Pr = Pr[1:]
            else:
                iFP += 1
        return

    @staticmethod
    def breed(parent1, parent2):
        child = []
        childP1 = []
        childP2 = []

        geneA = int(random.random() * len(parent1))
        geneB = int(random.random() * len(parent1))

        startGene = min(geneA, geneB)
        endGene = max(geneA, geneB)

        return parent1[:startGene] + parent2[startGene:endGene] + parent1[endGene:]

    @staticmethod
    def breedPopulation(matingpool, eliteSize):
        children = []
        length = len(matingpool) - eliteSize
        pool = random.sample(matingpool, len(matingpool))

        for i in range(0,eliteSize):
            children.append(matingpool[i])

        for i in range(0, length):
            child = Genetic.breed(pool[i], pool[len(matingpool)-i-1])
            children.append(child)
        return children

    @staticmethod
    def sort_parents(parents):
        return sorted(parents, key=lambda x:x[1], reverse=True)

    @staticmethod
    def mutate(children, rate):
        l = len(children[0])
        n_children = len(children)
        tot_num_genes = n_children * l
        num_mutations = int(np.rint(rate * tot_num_genes))
        num_indxs = np.random.randint(0, high=tot_num_genes - 1, size=num_mutations)
        for i in num_indxs:
            #print (i)
            i_child = int(i / l)
            i_gene = i % l
            #print ("mutating {},{}".format(i_child, i_gene))
            children[i_child][i_gene] = not children[i_child][i_gene]

    def plan(self, pos, angle_abs, target, stations, tracker, params):
        self.tracker = tracker
        self.hit_dist = ArenaParams.station_rad - self.tracker.pod_params.pod_rad
        self.seg = Planner3()
        self.seg.plan(pos,angle_abs, target, stations, tracker, params)

    def get_planner_init(self, pos, angle_abs, target):
        thrust_cmds = []
        ang_cmds = []
        pos = self.tracker.pos[-1]
        vel = self.tracker.vels[-1]
        angle = self.tracker.angle

        targs = []

        for _ in range(Genetic.NUM_OPT_COMMANDS):
            target, thrust, boost, shield = self.seg.act(pos, angle_abs, target)
            targs.append(target)
            p12 = np.array(target) - np.array(pos)
            ang_rad = math.atan2(p12[1], p12[0])
            ang_deg = np.degrees(ang_rad)
            if ang_deg < 0:
                ang_deg += 360
            thrust_cmds.append(thrust)
            ang_cmds.append(ang_deg)
            pos, vel, angle = self.simulator.calc_next_turn(target, pos, vel, angle, thrust, boost, shield)
            d = dist_pnts(pos, self.tracker.stations[0])
            #print ("PD={}".format(d), file=sys.stderr)
        #print ("REPLAY to see if we get the same distance", file=sys.stderr)

        #pos = self.tracker.pos[-1]
        #vel = self.tracker.vels[-1]
        #angle = self.tracker.angle
        #for targ, t, a in zip(targs, thrust_cmds, ang_cmds):
        #    p12 = np.array(targ) - np.array(pos)
        #    ang_targ = np.degrees(np.math.atan2(p12[1], p12[0]))
        #    targ = (0.1 * targ[0], 0.1 * targ[1])
        #    print ("ang_targ={} a={}".format(ang_targ, a), file=sys.stderr)
        #    pos, vel, angle = self.simulator.calc_next_turn(targ, pos, vel, angle, t, False, False)
        #ppp = self.simu(thrust_cmds, ang_cmds)
        #d = dist_pnts(ppp[-1], self.tracker.stations[0])
        ##d = dist_pnts(pos, self.tracker.stations[0])
        #print ("after2 replay the distance is {}".format(d), file=sys.stderr)
        return Genetic.chromo_encode(thrust_cmds, ang_cmds)


    def act(self, pos,  angle_abs, target):
        start = time.time()
        planner_chromo = self.get_planner_init(pos, angle_abs, target)
        crms = self.guess_initial_set(Genetic.GENERATION_SIZE)
        crms = [planner_chromo] + crms
        #sys.exit(255)
        #print ("fit={}".format(gg, file=sys.stderr))

        best_sol = []
        best_fit = 0.
        while (time.time() - start) < 0.06:
            fit = self.fitness(crms)

            for f, c in zip(fit, crms):
                if f > best_fit:
                    best_fit = f
                    best_sol = c[:]
            tot = sum(fit)
            FP = [ g / tot for g in fit]
            #print ("FP={}".format(FP), file=sys.stderr)
            #print (FP)
            ns = list(self.natural_selection(FP))
            #print ("ns={}".format(ns), file=sys.stderr)
            parents = [(crms[i], fit[i]) for i in ns]
            #print (parents)
            parents = Genetic.sort_parents(parents)
            par_chromo = [p[0] for p in parents]
            eliteSize = 2
            crms = Genetic.breedPopulation(par_chromo, eliteSize)
            Genetic.mutate(crms, 0.001)
        #print("best_fit={}".format(best_fit), file=sys.stderr)
        t_v, a_v = Genetic.chromo_decode(best_sol)
        thrust = t_v[0]
        ang = a_v[0]
        ang_r = np.radians(ang)
        n_target = (math.cos(ang_r) * 30000, math.sin(ang_r) * 30000)
        end = time.time()
        print ("time passed={}".format(end - start), file=sys.stderr)
        return n_target, thrust, False, False

#Genetic algoriths stuff ends here

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
                num_turns = s / v1 * 0.5
                n_p = np.array(self.tracker.stations[1]) - np.array(target)
                ang_abs_rad = np.radians(angle_abs)
                dir_ = np.array((math.cos(ang_abs_rad), math.sin(ang_abs_rad)))
                dang = math.fabs(angle_between(n_p, dir_))
                num_turns_to_rot = dang / PodParams.rot_vel
                #print ("num_turns_to_rot={} num_turns={}".format(num_turns_to_rot, num_turns), file=sys.stderr)
                #if num_turns > num_turns_to_rot:
                #    #print ("NO NEED TO NEUTRAL SO LONG", file=sys.stderr)
                #    num_turns = num_turns_to_rot
                dist = dist_pnts(pos, target)
                turns_to_targ = dist / v1
                if turns_to_targ <= num_turns:
                    return self.tracker.stations[1], 0, False, False
        return tc, thrust, False, False

class Pursuit:
    """
    a strategy to pursuit the opponent leader and deflect him
    from getting to the next target
    """
    def __init__(self, tracker):
        self.tracker = tracker
        self.algo = BlindPlanner()

    def plan(self, pos, angle, target, stations, tracker, params
            ):
        """
        stations - the list of stations, the post to defend is the first station
        """
        self.stations       = stations[:]
        self.tracker        = tracker
        self.algo.plan(pos, angle, target, stations, tracker, params)
        self.tolerance      = ArenaParams.station_rad - self.tracker.pod_params.pod_rad

    def act(self, pos, angle_abs, target):
        opp = opponent_leader
        target, vel, angle = opp.predictor.predict(self.tracker.pod_params.pursuit_look_ahead_turns)
        print ("PREDICTION={} mypos={} otherpos={}".format(target, self.tracker.pos[-1], opp.pos[-1]), file=sys.stderr)
        f = ForceBasedCollisionAvoidance(self.tracker)
        t = f.get_evading_force(Tracker.me[0], 3200)
        #if t[0] != 0 or t[1] != 0:
        #    targ, thrust = self.tracker.vec_to_thrust_target(t)
        #    return targ, thrust, False, False
        return self.algo.act(pos, angle_abs, target)

class GuardPostStrategy:
    """
    A strategy where the pod rush to a given position
    and does not let the opponent leader to reach it.
    The strategy ends when the opponent reach that target
    """
    def __init__(self, tracker):
        self.tracker        = tracker
        self.target         = (0, 0)
        self.algo           = BlindPlanner()
        self.pursuit        = Pursuit(tracker)
        self.station        = []
        self.in_position    = False
        self.rot_to_new_pos = False

    def plan(self, pos, angle, target, stations, tracker, params
            ):
        """
        stations - the list of stations, the post to defend is the first station
        """
        self.stations       = stations[:]
        self.tracker        = tracker
        self.target         = target
        self.in_position    = False
        self.rot_to_new_pos = False
        self.algo.plan(pos, angle, target, stations, tracker, params)
        self.pursuit.plan(pos, angle, target, stations, tracker, params)
        self.tolerance      = ArenaParams.station_rad - self.tracker.pod_params.pod_rad
        self.pursuing       = False
        print ("PPPPPPPPPPPLAN", file=sys.stderr)

    def act(self, pos, angle_abs, target):
        d = dist_pnts(self.tracker.pos[-1], opponent_leader.pos[-1])
        f = ForceBasedCollisionAvoidance(opponent_leader)
        in_way, fac = f.in_the_way(self.tracker)
        is_in_yeshoret = self.tracker.is_in_yeshoret(opponent_leader)
        print ("pusu={},is_way={},yesho={}".format(self.pursuing, in_way, is_in_yeshoret), file=sys.stderr)
        face = self.tracker.face_firection()
        p_rel = np.array(opponent_leader.pos[-1]) - np.array(self.tracker.pos[-1])
        f_alfa = angle_between(face, p_rel)
        if not self.tracker.is_in_yeshoret(opponent_leader) or f_alfa > 15:
            in_way = False

        if (d < 4000 and fac >= 1) or in_way or self.pursuing:
            self.pursuing = True
            return self.pursuit.act(pos, angle_abs, target)
        else:
            return self.act_prepare(pos, angle_abs, target)

    def act_prepare(self, pos, angle_abs, target):
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
            #print ("ANGLE ADJUST = {}".format(angle_to_new), file=sys.stderr)
            if angle_to_new > 30:
                return target, 0, False, False
            return self.algo.act(pos, angle_abs, target)
        else:
            ar = self.algo.act(pos, angle_abs, target)
            f = ForceBasedCollisionAvoidance(self.tracker)
            t = f.get_evading_force(Tracker.me[0], 3200)
            if t[0] != 0 or t[1] != 0:
                targ, thrust = self.tracker.vec_to_thrust_target(t)
                return targ, thrust, False, False
            else:
                return ar

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
        #print ("START WITH BOOST {} {}".format(self.tracker.arena_params.start_with_boost, self.tracker.boost), file=sys.stderr)
        if self.tracker.arena_params.start_with_boost and  not self.tracker.boost:
            boost = True
        return tc, thrust, boost, False

class GoToTargetRegulator:
    def __init__(self):
        self.reset()

    def reset(self, tolerance=0, tracker=None, params=None):
        self.e               = 0
        self.ie              = 0
        self.last_e          = 0
        self.dedt            = 0
        self.params          = params
        self.tracker         = tracker
        self.tolerance       = tolerance
        self.is_pointing     = False
        self.thrust_reducing = PodParams.max_thrust
        self.is_reducing     = False

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
        if dir_[0] != 0. or dir_[1] != 0.:
            self.e = angle_between(pt, dir_)
            #print ("angle between is {}".format(self.e), file=sys.stderr)
        else:
            self.e = 0
        #print ("RESULTED ERROR ANGLE={}".format(self.e), file=sys.stderr)
        if math.fabs(self.e) > 89:
            v = linalg.norm(self.tracker.vel)
            if True or v > 200:
                self.reset(self.tolerance, self.tracker, self.params)
                #thrust = 11
                if not self.is_reducing:
                    self.is_reducing = True
                    self.thrust_reducing = thrust - 30
                else:
                    self.thrust_reducing -= 30
                if self.thrust_reducing < 80:
                    self.thrust_reducing = 80
                return target, self.thrust_reducing
                #return target, thrust
        self.ie += self.e
        self.dedt = self.last_e - self.e
        Kp = self.params.gtKp
        Ki = self.params.gtKi
        Kd = self.params.gtKd
        c_angle = -(Kp*self.e + Ki*self.ie + Kd*self.dedt)
        #print ("c_angle={},angle={}".format(c_angle, angle), file=sys.stderr)
        #if c_angle > 89:
        #    c_angle = 89
        #    thrust = 20
        #elif c_angle < -89:
        #    c_angle = 89
        #    thrust = 20
        self.last_e = self.e
        pt1 = rotate(pt, degrees=c_angle)
        a_pt = rotate(pt, degrees=angle)
        diff = math.fabs(angle_between(pt1, a_pt))
        #print ("DIFF={}".format(diff), file=sys.stderr)
        res = np.array(pos) + pt1
        d = dist_pnts(pos, target)
        print ("DDD={}".format(d), file=sys.stderr)
        if diff > 18 and d < 3500:
            if not self.is_reducing:
                self.is_reducing = True
                self.thrust_reducing = thrust - 30
            else:
                self.thrust_reducing -= 30
            if self.thrust_reducing < 80:
                self.thrust_reducing = 80
            return (res[0], res[1]), self.thrust_reducing
        self.is_reducing = False
        return (res[0], res[1]), thrust

    def calc_target_vel(self, target):
        p1 = self.tracker.pos[-1]
        tg = unit_vector(self.tracker.vel)
        rad, _ = find_rad_from_two_points_and_tangent(p1, tg, target)
        rad += 600
        alphad = math.fabs(self.tracker.pod_deflection())
        alpha = math.fabs(np.radians(alphad))
        vt = rad * ( 1 - ArenaParams.friction_fac) * math.tan(alpha)
        #print ("vt={}".format(vt), file=sys.stderr)
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
        #print ("distance from target={} targ={} pinpoint={} fac={}".format(d, target, self.t, fac), file=sys.stderr)
        if d < self.tolerance and fac >= 0:
            return True
        return False

    def regulate_thrust(self, target, angle):
        pod_def = self.tracker.pod_deflection()
        #print ("pod deflection={}".format(pod_def), file=sys.stderr)
        if angle < 1.0 and math.fabs(self.tracker.pod_deflection()) < 1:
            v_tar, alpha = 99999999, 0
        else:
            v_tar, alpha = self.calc_target_vel(target)
        thrust = v_tar*(1 - ArenaParams.friction_fac) / math.cos(alpha)
        #print ("v_tar={} alpha={} thrust={}".format(v_tar, alpha, thrust), file=sys.stderr)
        thrust /= PodParams.fric_thruts_to_thrust_ratio
        if thrust > PodParams.max_thrust:
            thrust = PodParams.max_thrust
        elif thrust < 0:
            thrust = 0
        self.thrust = int(thrust)
        return self.thrust

class ForceBasedCollisionAvoidance:
    def __init__(self, tracker):
        self.tracker = tracker

    def calc_tot_evading_force(self):
        d2next = dist_pnts(self.tracker.pos[-1], self.tracker.stations[0])
        #print ("dnext={}".format(dnext), file=sys.stderr)
        if d2next < 1500:
            return np.array((0, 0))

        t1 = self.get_evading_force(opponent_leader, 1600)
        t2 = self.get_evading_force(opponent_follower, 1600)
        #print ("t1={},t2={}".format(t1, t2), file=sys.stderr)
        return t1 + t2

    def in_the_way(self, other):
        p_t = np.array(self.tracker.stations[0]) - np.array(self.tracker.pos[-1])
        p_other = np.array(other.pos[-1]) - np.array(self.tracker.pos[-1])
        res = math.fabs(angle_between(p_t, p_other))
        f = location_along_segment(self.tracker.pos[-1], self.tracker.stations[0], other.pos[-1])
        print ("F={}, al={}".format(f, res), file=sys.stderr)
        if f <= 1 and f > 0:
            return (res < 15), f
        return False, f

    def get_evading_force(self, other, rad):
        d1 = dist_pnts(other.pos[-1], self.tracker.pos[-1])
        #print ("d1={}".format(d1), file=sys.stderr)
        if d1 > rad:
            return np.array((0, 0))
        vrel = np.array(other.vel) - np.array(self.tracker.vel)
        #if linalg.norm(vrel) < 300:
        #    return np.array((0, 0))
        p1 = other.pos[-1]
        p2 = (p1[0] + vrel[0], p1[1] + vrel[1])
        coll = point_to_line_dist((p1, p2), (0., 0.))
        if coll < 250:
            fac = location_along_segment(p1, p2, self.tracker.pos[-1])
            print ("DETECTED COLLISION ac={}".format(fac), file=sys.stderr)
            if fac <= 0: #check if our pod is behind
                return ((0., 0.))
            t = vrel * ( 1 - ArenaParams.friction_fac) * 0.5
            #if self.in_the_way(other):
            #    t = rotate(t, degrees=90)
            return t
        return np.array((0, 0))

class Tracker:
    me = None
    other = None

    def __init__(self, num_laps, me, id=-9999):
        self.num_laps           = 0
        self.arena_params       = None
        self.pod_params         = None
        self.pos                = []
        self.angles             = []
        self.vels               = []
        self.thrusts            = []
        self.vel                = None
        self.angle              = 0
        self.next_cp            = 0
        self.prev_cp            = 0
        self.stations           = []
        self.time_left_to_punch = 0
        self.num_laps           = num_laps
        self.passed_stations    = -1
        self.me                 = me # type: boolean
        self.boost              = False
        self.boost_turns        = 0
        self.predictor          = MotionPredictor(self)
        self.predictions        = []
        self.id = id

    def thrust_vec(self, thrust, target):
        return unit_vector(np.array(target) - np.array(self.pos[-1])) * thrust

    def vec_to_thrust_target(self, thrust_vec):
        t = linalg.norm(thrust_vec)
        pos = np.array(self.pos[-1]) + thrust_vec
        if t > PodParams.max_thrust:
            t = PodParams.max_thrust
        return (pos[0], pos[1]), int(np.trunc(t))

    def act(self):
        #print ("VEL {}".format(linalg.norm(self.vel)), file=sys.stderr)
        if self.prev_cp != self.next_cp:
            print ("TRACKER:ACT now planning", file=sys.stderr)
            self.pod_params.planner.plan(self.pos[-1], self.angle, self.stations[0], self.stations, self, self.arena_params)

        shield = self.need_protection()
        #if self.time_left_to_punch > 0:
        #    self.time_left_to_punch -= 1
        #elif not shield and self.attempt_punch():
        #    return
        tc, thrust, boost, shield2 = self.pod_params.planner.act(self.pos[-1], self.angle, self.stations[0])
        f = ForceBasedCollisionAvoidance(self)
        t = f.calc_tot_evading_force()
        #print ("EVADE={}".format(t), file=sys.stderr)
        if t[0] == 0.0 and t[1] == 0.0:
            shield |= shield2
            self.write(tc, thrust, boost, shield)
        else:
            t0 = self.thrust_vec(thrust, tc)
            #print ("TARGET THRUSR={}".format(t0), file=sys.stderr)
            tot = t0 + 3 * t
            tc2, thrust = self.vec_to_thrust_target(tot)
            self.write(tc2, thrust, False, False)

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

    def face_firection(self):
        angle_rad = np.radians(self.angles[-1])
        return math.cos(angle_rad), math.sin(angle_rad)

    def calc_last_thrust(self):
        if len(self.vels) < 2:
            return 0
        vn1 = np.array(self.vels[-1]) / ArenaParams.friction_fac
        vn  = np.array(self.vels[-2])
        t   = vn1 - vn
        return linalg.norm(t)

    @staticmethod
    def record_data(data, storage):
        storage.append(data)
        if len(storage) > 10:
            storage = storage[1:]
        return storage

    def store_data(self, pos, vel, angle):
        thrust = self.calc_last_thrust()
        #print ("LAST THRUST={}".format(thrust), file=sys.stderr)
        self.pos     = Tracker.record_data(pos, self.pos)
        self.angles  = Tracker.record_data(angle, self.angles)
        self.vels    = Tracker.record_data(vel, self.vels)
        self.thrusts = Tracker.record_data(thrust, self.thrusts)

    def report_pos(self, pos, vel, angle, next_cp):
        self.store_data(pos, vel, angle)
        self.angle = angle
        self.vel = vel
        if next_cp != self.next_cp:
            self.passed_stations += 1
            self.stations = arena_detector.stations[next_cp:] + arena_detector.stations[:next_cp]
        self.prev_cp = self.next_cp
        self.next_cp = next_cp
        if len(self.predictions) >= 1:
            d = dist_pnts(pos, self.predictions[-1])
            print ("prediction dist={}".format(d), file=sys.stderr)
        ppos, vel_, ang_ = self.predictor.predict(1)
        self.predictions = Tracker.record_data(ppos, self.predictions)

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

        #print ("ANGLE={} face={} vel={}".format(self.angle, face, self.vel), file=sys.stderr)
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
        #print ("GOING_TO_COLLIDE d={} vc={} num_turns={}".format(d, vc, num_turns), file=sys.stderr)
        return (d / vc) <= num_turns

    def attempt_punch(self):
        if len(self.pos) < 2:
            return False
        for other in Tracker.other:
            if self.can_deliver_puch(other):
                rvel = self.rel_vel(other)
                shell = self.going_to_collide(other, 1.0) and (rvel > self.pod_params.rel_vel_for_crash_with_shield or self.boost_turns > 0)
                #print ("PUNCH {} shell={} boost_turns={}!!!!!!".format(rvel, shell, self.boost_turns), file=sys.stderr)
                self.write(other.next_pos(), PodParams.max_thrust, True, shell)
                self.time_left_to_punch = self.pod_params.num_punchs
                return True
        return False

    def need_protection(self):
        for other in Tracker.other:
            if len(other.pos) < 2:
                continue
            rvel = self.rel_vel(other)
            print ("rvel={}".format(rvel), file=sys.stderr)
            if rvel <= self.pod_params.rel_vel_for_crash_with_shield:
                continue
            if self.going_to_collide(other, 1.2):
                return True
        return False


    def can_deliver_puch(self, other):
        if not self.going_to_collide(other, self.pod_params.num_of_turns_to_impact):
            return False
        fc_dir = np.degrees(math.atan2(self.vel[1], self.vel[0]))
        fdir = self.direction_rel_to(other)
        rel_angle = math.fabs(fc_dir - fdir)
        if rel_angle <= self.pod_params.side_punch_max_angle:
            #print ("MAY PUNCH", file=sys.stderr)
            return True

        return False

    def next_pos(self):
        """
        try to estimate the position in the next loop
        """
        ppos, vel, angle = self.predictor.predict(1)
        return ppos

class Defender(Tracker):
    def __init__(self, num_laps, me, id=-9999):
        Tracker.__init__(self, num_laps, me, id)
        self.algo               = GuardPostStrategy(num_laps)
        self.pursuit            = False
        self.pusuit_has_planned = False
        self.stations           = []
        self.station_id         = -1
        self.leader_origin      = (0, 0)
        self.leader_target      = (0, 0)
        self.leader_next        = (0, 0)

    def calc_target_full_block(self, leader):
        p23 = np.array(leader.pos[-1]) - np.array(self.leader_target)
        p23 = p23 / linalg.norm(p23) * 1200
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

        #print ("ME STATIONS:{} {}".format(all_leader.passed_stations, self.station_id), file=sys.stderr)
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
        #print ("MY_DIST = {}, leader_DIST={}".format(my_dist, leader_dist), file=sys.stderr)
        return (my_dist + self.arena_params.defender_dist_spare) <= leader_dist

    def act(self):
        if self.pursuit:
            if not self.pusuit_has_planned:
                self.algo.plan(self.pos[-1], self.angle, self.stations[0], self.stations, self, self.arena_params)
                self.pusuit_has_planned = True
            target = opponent_leader.stations[0]
            shield = self.need_protection()
            #print ("defender shield={}".format(shield), file=sys.stderr)

            if self.time_left_to_punch > 0:
                self.time_left_to_punch -= 1
            elif not shield and self.attempt_punch():
                return

            tc, thrust, boost, shield2 = self.algo.act(self.pos[-1], self.angle, target)
            shield |= shield2
            self.write(tc, thrust, boost, shield)
        else:
            self.act_guard()

    def is_in_yeshoret(self, other):
        return self.station_id == (other.passed_stations + 1)

    def act_guard(self):
        leader = opponent_leader
        #print ("CAN_REACH_LEADER={}".format(self.can_reach_target_before_leader(leader)), file=sys.stderr)
        jump = 0
        while not self.can_reach_target_before_leader(leader):
            self.leader_origin = leader.stations[jump - 1]
            p23 = np.array(self.leader_origin) - np.array(leader.stations[jump])
            self.leader_target = leader.stations[jump]
            #print ("SETTING TARGET TO {}".format(self.leader_target), file=sys.stderr)
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
        #print ("defender shield={}".format(shield), file=sys.stderr)

        if self.time_left_to_punch > 0:
            self.time_left_to_punch -= 1
        elif not shield and self.attempt_punch():
            return

        tc, thrust, boost, shield2 = self.algo.act(self.pos[-1], self.angle, target)
        shield |= shield2
        self.write(tc, thrust, boost, shield)

    def report_pos(self, pos, vel, angle, next_cp):
        self.store_data(pos, vel, angle)
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
    tracks = [
              Arena("hostile",    [(13890, 1958), (8009,  3263), (2653,  7002),  (10035,            5969)], HostileParams()),
              Arena("hostile2",   [(9409,  7247), (5984,  4264), (14637, 1420),  (3470,             7203)], HostileParams()),
              Arena("pyramid",    [(7637,  5988), (3133,  7544), (9544,  4408),  (14535,            7770),  (6312,             4294),  (7782,             851)],  PyramidParams()),
              Arena("triangle",   [(6012,  5376), (11308, 2847), (7482,  6917)], TriangleParams()),
              Arena("dalton",     [(9589,  1395), (3618,  4434), (8011,  7920),  (13272,            5530)], DaltonParams()),
              Arena("makbilit",   [(12942, 7222), (5655,  2587), (4107,  7431),  (13475,            2323)], MacbilitParams()),
              Arena("arrow",      [(10255, 4946), (6114,  2172), (3048,  5194),  (6276,             7762),  (14119,            7768),  (13897,            1216)], ArrowParams()),
              Arena("Shosh",      [(9116,  1857), (5007,  5288), (11505, 6074)], ShoshParams()),
              Arena("Til",        [(10558, 5973), (3565,  5194), (13578, 7574),  (12430,            1357)], TilParams()),
              Arena("trapez",     [(11201, 5443), (7257,  6657), (5452,  2829),  (10294,            3341)], TrapezParams()),
              Arena("Mehumash",   [(4049,  4630), (13054, 1928), (6582,  7823),  (7494,             1330),  (12701,            7080)], MehumashParams()),
              Arena("Trampoline", [(3307,  7251), (14572, 7677), (10588, 5072),  (13100,            2343),  (4536,             2191),  (7359,             4930)], TrampolineParams()),
              Arena("Zigzag",     [(10660, 2295), (8695,  7469), (7225,  2174),  (3596,             5288),  (13862,            5092)], ZigzagParams())
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

class MotionPredictor:
     def __init__(self, tracker):
         self.steady_state = False
         self.tracker      = tracker
         self.simulator    = Simulator()

     def predict_vals(self, t0 ,t1):
         return 2 * t1 - t0

     def predict(self, turns):
         if not self.steady_state:
             l = len(self.tracker.pos)
             if l == 1:
                 return self.tracker.pos[-1], self.tracker.vels[-1], self.tracker.angles[-1]
             elif l == 2:
                 p21 = (np.array(self.tracker.pos[-1]) - np.array(self.tracker.pos[-2])) * turns
                 res = p21 + np.array(self.tracker.pos[-2])
                 return (res[0], res[1]), self.tracker.vels[-1], self.tracker.angles[-1]
             else:
                 self.steady_state = True
         t1 = self.predict_vals(self.tracker.thrusts[-2], self.tracker.thrusts[-1])
         pos, vel, angle = self.next_turn(turns, self.tracker.stations[0], self.tracker.pos[-1], self.tracker.vels[-1], self.tracker.angles[-1], self.tracker.thrusts[-1], t1)
         return pos, vel, angle

     def next_turn(self, num_turns, target, pos, vel, angle, t0, t1):
         pos, vel, angle = self.simulator.calc_next_turn(target, pos, vel, angle, t1, False, False)
         num_turns -= 1
         if num_turns == 0:
             return pos, vel, angle
         t = self.predict_vals(t0, t1)
         return self.next_turn(num_turns, target, pos, vel, angle, t1, t)

     #def predict_target(self, pos, angle, thrust):
     #    ang_rad = np.radians(angle)
     #    if angle <= 1.0:
     #        pface = np.array((math.cos(ang_rad), math.sin(ang_rad))) * 3000.
     #        targ = np.array(pos) + pface
     #        return (targ[0], targ[1])



arena_detector = ArenaDetector()
defence_params = DefencePodParams()
offence_params = OffencePodParams()

opponent_leader   = None
opponent_follower = None
all_leader        = None
all_second        = None

if __name__ == "__main__":
    num_laps = int(input())
    Tracker.me = [Tracker(num_laps, me=True, id=0), Tracker(num_laps, me=True, id=1)]
    Tracker.other = [Tracker(num_laps, me=False), Tracker(num_laps, me=False)]
    Tracker.me[0].pod_params = offence_params
    Tracker.me[1].pod_params = offence_params#defence_params
    Tracker.other[0].pod_params = OpponentPodParams()
    Tracker.other[1].pod_params = OpponentPodParams()

    check_point_count = int(input())
    stations = []
    #print ("num_check_points={}".format(check_point_count), file=sys.stderr)
    for _ in range(check_point_count):
        xs, ys = input().split()
        stations.append((int(xs), int(ys)))

    arena_detector.try_detect(stations)
    all_trackers = [Tracker.me[0], Tracker.me[1], Tracker.other[0], Tracker.other[1]]

    Tracker.me[0].configure(arena_detector, num_laps)
    Tracker.me[1].configure(arena_detector, num_laps)

    # game loop
    i = 0
    assigned = False
    while True:
        for tracker in all_trackers:
            s = input()
            print ("========={}".format(s), file=sys.stderr)
            x, y, vx, vy, angle, next_cp = [int(i) for i in s.split()]
            tracker.report_pos((x, y), (vx, vy), angle, next_cp)

        my_leader, my_follower = Tracker.leader(Tracker.me[0], Tracker.me[1])
        opponent_leader, opponent_follower = Tracker.leader(Tracker.other[0], Tracker.other[1])
        all_leader, all_second = Tracker.leader(my_leader, opponent_leader)

        print ("OFFENCE DATA", file=sys.stderr)
        Tracker.me[0].act()
        print ("DEFENCE DATA", file=sys.stderr)
        Tracker.me[1].act()
        if assigned == False and i >= 20 or (my_leader.passed_stations - my_follower.passed_stations) > 1:
            assigned = True
            id = my_follower.id
            Tracker.me[id] = Defender(num_laps, me=True, id=id)
            Tracker.me[id].passed_stations = my_follower.passed_stations
            Tracker.me[id].pod_params = defence_params
            Tracker.me[id].configure(arena_detector, num_laps)
            Tracker.me[id].stations = my_follower.stations[:]
            Tracker.me[id].pos = my_follower.pos[:]
            Tracker.me[id].vels = my_follower.vels[:]
            Tracker.me[id].angles = my_follower.angles[:]
            all_trackers[id] = Tracker.me[id]
            my_follower = Tracker.me[id]
        i += 1

