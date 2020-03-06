#/usr/bin/env python3

import math
import matplotlib.pyplot as plt
from solution import *

rp = RoutPlanner(500.0)
rp.add_station((0, 0))
rp.add_station((1200, 0))
rp.add_station((600, 600))

xvals, yvals = rp.plan((0, 0), 90., (1200, 0))
plt.xlim(-2300, 1500)
plt.ylim(-2300, 1500)
plt.plot(xvals, yvals)
plt.show()
