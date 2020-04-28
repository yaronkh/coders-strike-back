#!/usr/bin/python3
import fcntl
import sys
import os
import time
import math
import subprocess
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

class ArenaParams:
    station_rad = 600
    friction_fac = 0.85 # the current speed vector of each pod is multiplied by that factor
    station_tolerance = 450

class PodParams:
    max_thrust = 200
    rot_vel = 18
    pod_rad = 200

def dist_pnts(m, t):
    mx, my = m
    tx, ty = t
    return math.sqrt((mx - tx)**2 + (my - ty)**2)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def rotate(p, origin=(0, 0), degrees=0):
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(p)
    return np.squeeze((R @ (p.T-o.T) + o.T).T)

def vec_angle(p):
    ret = math.atan2(p[1], p[0])
    ret = np.degrees(ret)
    if ret < 0:
        ret += 360
    return ret

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

def dist_pnts(m, t):
    mx, my = m
    tx, ty = t
    return math.sqrt((mx - tx)**2 + (my - ty)**2)

class Pod:
    def __init__(self, pos, angle, stations):
        self.pos       = pos
        self.angle     = angle
        self.vel       = (0, 0)
        self.target    = (0, 0)
        self.stations  = stations[:]
        self.thrust    = 0
        self.shield    = False
        self.boost     = 0
        self.boosted   = False
        self.sh_turns  = 0
        self.cp_id     = 0
        self.simulator = Simulator()

    def write(self, pipe):
        next_cp = (self.cp_id + 1)% len(self.stations)
        write_to_pipe("{} {} {} {} {} {}\n".format(self.pos[0], self.pos[1], self.vel[0], self.vel[1], self.angle, next_cp), pipe)

    def read(self, r):
        if self.sh_turns > 0:
            self.sh_turns -= 1
        px, py, thrust = r.split(" ")
        if thrust == "BOOST":
            if self.boosted:
                self.thrust = PodParams.max_thrust
            else:
                self.boosted = True
                self.thrust = 650
        elif thrust == "SHIELD":
            self.sh_turns = 3
            self.shield = True
        else:
            self.thrust = int(thrust, 10)
        if self.sh_turns > 0:
            self.thrust = 0
        self.target = (int(px, 10), int(py, 10))

    def next(self):
       ret = self.simulator.calc_next_turn(self.target, self.pos, self.vel, self.angle, self.thrust, self.boost, self.shield)
       self.pos, self.vel, self.angle = ret
       d = dist_pnts(self.pos, self.stations[(self.cp_id + 1) % len(self.stations)])
       print ("d={}".format(d))
       if d <= (ArenaParams.station_rad - PodParams.pod_rad):
           self.cp_id += 1


f = open("/tmp/input", 'w')
def write_to_pipe(s, pipe):
    pipe.write(s.encode())
    f.write(s)
    pipe.flush()

class Arena:
    def __init__(self, name, stations):
        self.name = name
        self.stations = stations

    def point_in_track(self, p):
        for ps in self.stations:
            if math.fabs(p[0] - ps[0]) < 70 and math.fabs(ps[1] - ps[1]) < 70:
                return True
        return False


tracks = {
      Arena("hostile",    [(13890, 1958), (8009,  3263), (2653,  7002),   (10035, 5969)]),
      Arena("hostile2",   [(9409,  7247), (5984,  4264), (14637, 1420),   (3470,  7203)]),
      Arena("pyramid",    [(7637,  5988), (3133,  7544), (9544,  4408),   (14535, 7770),   (6312,  4294),   (7782,  851)]),
      Arena("triangle",   [(6012,  5376), (11308, 2847), (7482,  6917)]),
      Arena("dalton",     [(9589,  1395), (3618,  4434), (8011,  7920),   (13272, 5530)]),
      Arena("makbilit",   [(12942, 7222), (5655,  2587), (4107,  7431),   (13475, 2323)]),
      Arena("arrow",      [(10255, 4946), (6114,  2172), (3048,  5194),   (6276,  7762),   (14119, 7768),   (13897, 1216)]),
      Arena("Shosh",      [(9116,  1857), (5007,  5288), (11505, 6074)]),
      Arena("Til",        [(10558, 5973), (3565,  5194), (13578, 7574),   (12430, 1357)]),
      Arena("trapez",     [(11201, 5443), (7257,  6657), (5452,  2829),   (10294, 3341)]),
      Arena("Mehumash",   [(4049,  4630), (13054, 1928), (6582,  7823),   (7494,  1330),   (12701, 7080)]),
      Arena("Trampoline", [(3307,  7251), (14572, 7677), (10588, 5072),   (13100, 2343),   (4536,  2191),   (7359,  4930)]),
      Arena("Zigzag",     [(10660, 2295), (8695,  7469), (7225,  2174),   (3596,  5288),   (13862, 5092)])
}

TRACKS = {}

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

if __name__ == "__main__":
    os.environ["PYTHONUNBUFFERED"] = "1"
    for t in tracks:
        TRACKS[t.name] = t
    sel_track = TRACKS["pyramid"]
    num_laps = 1
    code = sys.argv[1]
    p = subprocess.Popen([code],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #don't block on subprocess stderr
    fd = p.stderr.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)

    write_to_pipe("{}\n".format(num_laps), p.stdin)
    write_to_pipe("{}\n".format(len(sel_track.stations)), p.stdin)
    sx=[]
    sy=[]
    for s in sel_track.stations:
        write_to_pipe("{} {}\n".format(s[0], s[1]), p.stdin)
        sx.append(s[0])
        sy.append(s[1])

    p01 = np.array(sel_track.stations[0]) - np.array(sel_track.stations[-1])
    init_ang = int(vec_angle(p01))
    runner = Pod((4000, 4000), 0, sel_track.stations)
    follower = Pod((50, 0), 0, sel_track.stations)
    opp1 = Pod((0, 0), 0, sel_track.stations)
    opp2 = Pod((0, 0), 0, sel_track.stations)
    all_pods = [runner, follower, opp1, opp2]
    x=[]
    y=[]
    for _ in range(320):
        for pod in all_pods:
            pod.write(p.stdin)
        r = p.stdout.readline().decode()
        runner.read(r)
        r = p.stdout.readline().decode()
        follower.read(r)
        runner.next()
        x.append(runner.pos[0])
        y.append(runner.pos[1])
        while True:
            try:
                inp = p.stderr.readline().decode()
                if inp == '':
                    break
                #print ("DBG: {}".format(inp[:-1]))
            except:
                break
    fig, ax = plt.subplots()
    ax.scatter(sx, sy,s=[1200]*len(sx),alpha=0.5)
    ax.scatter(x, y, s=[50]*len(x), alpha=1)
    fig.tight_layout()
    plt.show()

