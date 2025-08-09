import pickle
import time as t
import sys
sys.path.append('.')

import numpy as np

init_time = t.time()
from IRSMicroLensing import IRSCaustics as IRSC
from IRSMicroLensing import IRSFunctions as IRSF
print(f'Custom library import time: {(t.time() - init_time):.3} seconds')

''' Preparing lens parameters '''
# Big planet parameters
s1 = 0.8
alpha1 = 0
q1 = 1e-3

# Annulus parameters
num_theta = 5000
num_r = 4 * num_theta

lens_attributes = [
    [0, 0, 0.001, 1],
    [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), 0.001, q1]
]

ang_width, thickness, (y_plus, y_minus), cusp_points = IRSC.IRSCaustics.ang_width_thickness_calculator(lens_attributes)

print(f'Angular width: {ang_width}')
print(f'Thickness: {thickness}')
print(f'Annulus lower bound: {y_plus}')
print(f'Annulus upper bound: {y_minus}')

lens_parameters = {
    'pixels': 1000,
    'ang_width': ang_width,
    'thickness': thickness,
    'y_plus': y_plus,
    'y_minus': y_minus,
    'lens_att': lens_attributes,
    'num_theta': num_theta,
    'num_r': num_r
}

print(f'Number of rays: {(num_r * num_theta):.4e}')
print('=========================================================')

''' Simulating magnification map '''
binary_lens = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)

# num = int(num_theta / 2)
num = 100

subdivisions = np.arange(0, num+1, 1).astype(int)
print(f'Number of iterations: {len(subdivisions)}')

all_times = np.zeros_like(subdivisions, dtype=np.float64)

init_time = t.time()
binary_lens.calculate()
all_times[0] = t.time() - init_time

for i, subdivision in enumerate(subdivisions[1:]):
    print(f'Subdivisions: {subdivision}')
    init_time = t.time()
    binary_lens.series_calculate(subdivisions=subdivision)
    all_times[i+1] = t.time() - init_time

    dat = np.vstack((np.float64(subdivisions), all_times))

    np.save('./Unity/Data Files/timing_data_2.npy', dat)