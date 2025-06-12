import pickle
import time as t
import sys
import argparse
sys.path.append('.')

import platform
import psutil

print(platform.machine())
print(platform.version())
print(platform.system())
print(platform.processor())

memory = psutil.virtual_memory()

print(memory.total / (1024 ** 3))
print(memory.used / (1024 ** 3))
print(memory.available / (1024 ** 3))

import numpy as np

init_time = t.time()
from IRSMicroLensing import IRSCaustics as IRSC
from IRSMicroLensing import IRSFunctions as IRSF
print(f'Custom library import time: {(t.time() - init_time):.3} seconds')

# Initialize parser
parser = argparse.ArgumentParser(description='Compute single lens magnification map offsetted by triple lens effective center.')

# Adding arguments
parser.add_argument('-s2', '--sep2', help='Seperation of secondary small planet')
parser.add_argument('-a2', '--angle2', help='Angle of secondary small planet')
parser.add_argument('-pmr', '--planet_mass_ratio', help='Small planet / big planet mass ratio')
parser.add_argument('-l', '--lenses', help='Shoot single, binary, or triple lens')
parser.add_argument('-o', '--origin', help='Coordinate origin')

args = vars(parser.parse_args())
print(args)

print()
print('=========================================================')
print('Lens parameters:')

''' Preparing lens parameters '''
# Big planet parameters
s1 = 0.8
alpha1 = 0
q1 = 1e-3

# Small planet parameters
s2 = np.float64(args['sep2'])
alpha2 = np.float64(args['angle2'])
q2 = np.float64(args['planet_mass_ratio']) * q1

print(f'q1 = {q1}')
print(f'q2 = {q2}')
print(f's1 = {s1}')
print(f's2 = {s2}')
print(f'alpha1 = {alpha1}')
print(f'alpha2 = {alpha2}')

# Defining single lens attributes
single_lens_attributes = np.array([
    [0, 0, 1]
])

# Defining binary lens attributes
binary_lens_attributes = np.array([
    [0, 0, 1],
    [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), q1],
])

# Triple lens attributes
triple_lens_attributes = np.array([
    [0, 0, 1],
    [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), q1],
    [s2*np.cos(np.deg2rad(alpha2)), s2*np.sin(np.deg2rad(alpha2)), q2]
])

# Rotation matrix for first planet
first_planet_DCM = np.array([
    [np.cos(np.deg2rad(alpha1)), -np.sin(np.deg2rad(alpha1))],
    [np.sin(np.deg2rad(alpha1)), np.cos(np.deg2rad(alpha1))]
])

# Rotation matrix for second planet
second_planet_DCM = np.array([
    [np.cos(np.deg2rad(alpha2)), -np.sin(np.deg2rad(alpha2))],
    [np.sin(np.deg2rad(alpha2)), np.cos(np.deg2rad(alpha2))]
])

# Calculating translation of coordinate origin for each planet
q1 = triple_lens_attributes[1, 2]
q2 = triple_lens_attributes[2, 2]

delta1_rot = np.array([q1 / ((1 + q1) * (s1 + 1/s1)), 0]).reshape(-1, 1)
delta2_rot = np.array([q2 / ((1 + q2) * (s2 + 1/s2)), 0]).reshape(-1, 1)

delta1 = np.dot(first_planet_DCM, delta1_rot).reshape(2)
delta2 = np.dot(second_planet_DCM, delta2_rot).reshape(2)

total_offset = (q1*delta1 + q2*delta2)/(q1 + q2)
print(f'Total offset: {total_offset}')

# Correcting binary lens attributes
binary_lens_attributes[:, :2] -= delta1

# Correcting triple lens attributes
if args['origin'] == 'triple_offset':
    triple_lens_attributes[:, :2] -= total_offset
elif args['origin'] == 'binary_offset':
    triple_lens_attributes[:, :2] -= delta1
elif args['origin'] == 'cm':
    CM_sum = 0

    for i in range(3):
        CM_sum += triple_lens_attributes[i, :2] * triple_lens_attributes[i, 2]

    total_M = np.sum(triple_lens_attributes[:, 2])

    lens_CM = CM_sum / total_M

    triple_lens_attributes[:, :2] -= lens_CM

# Map parameters
pixels = 2000
delta = 0.01

ang_width, thickness, (y_plus, y_minus), cusp_points = IRSC.IRSCaustics.ang_width_thickness_calculator(triple_lens_attributes)
num_r, num_theta = IRSC.IRSCaustics.num_ray_calculator(pixels, ang_width, delta, y_plus, y_minus)

print(f'Angular width: {ang_width}')
print(f'Thickness: {thickness}')
print(f'Annulus lower bound: {y_minus}')
print(f'Annulus upper bound: {y_plus}')
print(f'Number of rays in theta: {num_theta}')
print(f'Number of rays in r: {num_r}')

single_lens_parameters = {
    'pixels': pixels,
    'ang_width': ang_width,
    'thickness': thickness,
    'y_plus': y_plus,
    'y_minus': y_minus,
    'lens_att': single_lens_attributes.tolist(),
    'num_theta': num_theta,
    'num_r': num_r
}

binary_lens_parameters = single_lens_parameters.copy()
binary_lens_parameters.update({
    'lens_att': binary_lens_attributes.tolist()
})

triple_lens_parameters = single_lens_parameters.copy()
triple_lens_parameters.update({
    'lens_att': triple_lens_attributes.tolist()
})

print(f'Number of rays: {(num_r * num_theta):.4e}')
print('=========================================================')

''' Simulating L lens magnification map '''
if args['lenses'] == 'single':
    param_dict = single_lens_parameters
    file_path = f'./Unity/Simulations/Collection_pmr0.001/single_1e11.pkl'

elif args['lenses'] == 'binary':
    param_dict = binary_lens_parameters
    file_path = f'./Unity/Simulations/Collection_pmr0.001/binary_1e11.pkl'

elif args['lenses'] == 'triple':
    param_dict = triple_lens_parameters
    file_path = f'./Unity/Simulations/Collection_pmr0.001/triple_1e11_{int(alpha2)}_{q2:.0e}_{args["origin"]}.pkl'

else:
    raise ValueError(f'Wrong lens configuration passed in. Got {args["lenses"]}.')

print(f'Shooting {args["lenses"]}:')

calculator = IRSC.IRSCaustics(annulus_param_dict=param_dict)
magnifications = calculator.series_calculate(cm_offset='auto')

print('=========================================================')

''' Saving class data to file '''
init_time = t.time()
with open(file_path, 'wb') as calculator_file:
    pickle.dump(calculator, calculator_file)

print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

print('Done')