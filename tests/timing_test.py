import time as t

t_i = t.time()
from IRSMicroLensing import IRSMain
print('Main:', t.time() - t_i)

t_i = t.time()
from IRSMicroLensing import IRSPlotter
print('Plotter:', t.time() - t_i)

t_i = t.time()
from IRSMicroLensing import IRSLightCurves
print('Light curves:', t.time() - t_i)

t_i = t.time()
from IRSMicroLensing import IRSCaustics
print('Caustics:', t.time() - t_i)

t_i = t.time()
from IRSMicroLensing import IRSOrbiting
print('Orbiting:', t.time() - t_i)

t_i = t.time()
from IRSMicroLensing import IRSFunctions
print('Functions:', t.time() - t_i)