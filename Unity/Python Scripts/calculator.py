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

# Annulus parameters
num_theta = 5000
num_r = 4 * num_theta

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
    [0.0002470252450260668, 1.3212232075066055e-07, 0.001, 1]
]

ang_width, thickness, (y_plus, y_minus), cusp_points = IRSC.IRSCaustics.ang_width_thickness_calculator(triple_lens_attributes)

print(f'Angular width: {ang_width}')
print(f'Thickness: {thickness}')
print(f'Annulus lower bound: {y_plus}')
print(f'Annulus upper bound: {y_minus}')

single_lens_parameters = {
    'pixels': 1000,
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
single_lens_magnifications = single_lens.series_calculate(cm_offset='auto')

print('=========================================================')

# ''' Simulating triple lens magnification map '''
# triple_lens = IRSC.IRSCaustics(annulus_param_dict=triple_lens_parameters)
# triple_lens_magnifications = triple_lens.series_calculate(cm_offset='auto')

# print('=========================================================')

''' Saving class data to file '''
init_time = t.time()
with open('./Unity/Simulations/single_lens_cusp.pkl', 'wb') as single_lens_file:
    pickle.dump(single_lens, single_lens_file)

# with open('./Unity/Analysis 6-4/triple_lens_small.pkl', 'wb') as triple_lens_file:
#     pickle.dump(triple_lens, triple_lens_file)

print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')