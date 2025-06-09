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

# Small planet parameters
s2 = 0.8
alpha2 = 30
q2 = 1e-3 * q1

# Map parameters
pixels = 2000
delta = 0.01

triple_lens_attributes = [
    [0, 0, 0.001, 1],
    [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), 0.001, q1],
    [s2*np.cos(np.deg2rad(alpha2)), s2*np.sin(np.deg2rad(alpha2)), 0.001, q2]
]

triple_lens_arr = np.array(triple_lens_attributes)

CM_sum = 0

for i in range(3):
    CM_sum += triple_lens_arr[i, :2] * triple_lens_arr[i, 3]

total_M = np.sum(triple_lens_arr[:, 3])

lens_CM = CM_sum / total_M

single_lens_attributes = [
    [lens_CM[0], lens_CM[1], 0.001, 1]
]

ang_width, thickness, (y_plus, y_minus), cusp_points = IRSC.IRSCaustics.ang_width_thickness_calculator(triple_lens_attributes)
num_r, num_theta = IRSC.IRSCaustics.num_ray_calculator(pixels, ang_width, delta, y_plus, y_minus)

print(f'Angular width: {ang_width}')
print(f'Thickness: {thickness}')
print(f'Annulus lower bound: {y_plus}')
print(f'Annulus upper bound: {y_minus}')
print(f'Number of rays in theta: {num_theta}')
print(f'Number of rays in r: {num_r}')

single_lens_parameters = {
    'pixels': pixels,
    'ang_width': ang_width,
    'thickness': thickness,
    'y_plus': y_plus,
    'y_minus': y_minus,
    'lens_att': single_lens_attributes,
    'num_theta': num_theta,
    'num_r': num_r
}

triple_lens_parameters = single_lens_parameters.copy()
triple_lens_parameters.update({
    'lens_att': triple_lens_attributes
})

print(f'Number of rays: {(num_r * num_theta):.4e}')
print('=========================================================')

''' Simulating single lens magnification map '''
single_lens = IRSC.IRSCaustics(annulus_param_dict=single_lens_parameters)
single_lens_magnifications = single_lens.series_calculate(cm_offset='auto', subdivisions=100)

print('=========================================================')

''' Saving class data to file '''
init_time = t.time()
with open('./Unity/Simulations/single_lens_1e11_cm.pkl', 'wb') as single_lens_file:
    pickle.dump(single_lens, single_lens_file)

print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

print('Done')