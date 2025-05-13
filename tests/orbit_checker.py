import numpy as np
import matplotlib.pyplot as plt

# import sys
# import time as t
# sys.path.append('..')

from IRSMicroLensing import OrbitCalc_dv
from IRSMicroLensing import CelestialBodyData as cb

t_span = 1*(86400*365)
dt = t_span / 100000

bodies = np.array([
    [0, 0, 0, 0, 0, 0, cb.sun['mass'], -1],
    [cb.earth['a']*cb.AU2KM, cb.earth['e'], 0, 0, 0, 0, cb.earth['mass'], 0],
    [cb.moon['a'], cb.moon['e'], 0, 0, 0, 0, cb.moon['mass'], 1]
])

plotter = OrbitCalc_dv.OrbitPropagator(bodies=bodies, dt=dt, t_span=t_span, anim=True, limits=[1, 1e6, 'f'])
