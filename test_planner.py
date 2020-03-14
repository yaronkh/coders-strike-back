import math
import matplotlib.pyplot as plt
import solution2
from solution2 import *

planner = Planner()
#stations = [(0, 0), (1000, 0), (1000, 1000), (0, 1000)]
#stations = [(0, 0), (10000, 0), (5000, 5000), (2500, 7500)]
#stations = [(0, 0), (10000, 0), (0, 10000), (10000, 10000)]
stations =  [(12911, 7240), (5657, 2595), (4084, 7448), (13523, 2310)]
#stations.reverse()
xvals = [p[0] for p in stations] + [stations[0][0]]
yvals = [p[1] for p in stations] + [stations[0][1]]
plt.plot(xvals, yvals)
for i in range(len(stations)):
    planner.plan((0, 0), 100.0, 0.0, (100, 0), stations)
    print("pivot={},curve_st={}, facs={},angle={}".format(planner.pivot, planner.curve_st, planner.rfacs, planner.angle))
    st_x = [p[0] for p in planner.curve_st]
    st_y = [p[1] for p in planner.curve_st]
    ax = plt.gca()
    plt.scatter(st_x, st_y)
    plt.scatter(planner.pivot[0], planner.pivot[1])
    #c = plt.Circle(planner.pivot, planner.rad, color='red')
    #ax.add_artist(c)
    stations = stations[1:] + [stations[0]]
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
