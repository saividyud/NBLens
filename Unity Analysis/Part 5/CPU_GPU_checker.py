# Importing crap
import sys
sys.path.append('.')

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# import IRSMicroLensing.IRSFunctions as IRSF
# import IRSMicroLensing.IRSCaustics2 as IRSC
import IRSMicroLensing.IRSCausticsGPU as IRSCGPU

import numexpr

import pickle

import csv

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.titlesize'] = 20
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['figure.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['mathtext.fontset'] = 'cm'

q1 = 1e-3
s1 = 0.2
q2 = 3e-4
s2 = 0.6
alpha = 45

triple_lens_att_1 = np.array([
    [0, 0, 1], # Star at origin with mass 1
    [s1, 0, q1], # Big planet
    [s2*np.cos(np.radians(alpha)), s2*np.sin(np.radians(alpha)), q2] # Small planet at angle alpha
])

center_of_mass_big = np.array([q1*s1 / (1 + q1), 0])
center_of_mass_small = np.array([q2*s2 / (1 + q2)*np.cos(np.radians(alpha)), q2*s2 / (1 + q2)*np.sin(np.radians(alpha))])

center_of_mass = np.average(triple_lens_att_1[:, :2], axis=0, weights=triple_lens_att_1[:, 2])
center_of_magnification = np.array([q2 / ((1 + q2) * (s1 + 1/s1)), 0])

triple_lens_att_1[:, :2] -= center_of_magnification

pixels = 1000
ang_width = 0.2
min_mag = 6.754609827054506
y_plus = 1.0769985822006392
y_minus = 0.9285063760849294
num_r = 185900
num_theta = 46475
max_source_radius = 0.005273437499999999
min_source_radius = 5.2734374999999986e-05

print(f'Pixels: {pixels}')
print(f'Angular width: {ang_width}')
print(f'Minimum magnification: {min_mag}')
print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, thickness = {y_plus - y_minus}')
print(f'Average radius of annulus: {(y_plus + y_minus)/2}')
print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}')
print(f'Total number of rays: {(num_r * num_theta):.3e}')
print(f'Maximum source radius: {max_source_radius}')
print(f'Minimum source radius: {min_source_radius}')

triple_lens_parameters_1 = {
    'lens_att': triple_lens_att_1.tolist(),
    'pixels': pixels,
    'ang_width': ang_width,
    'num_r': num_r,
    'num_theta': num_theta,
    'y_plus': y_plus,
    'y_minus': y_minus,
    'thickness': y_plus - y_minus
}

print('Shooting GPU simulation...')
triple_sim_gpu = IRSCGPU.IRSCaustics(annulus_param_dict=triple_lens_parameters_1)
triple_sim_gpu.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

file_path = './Unity Analysis/Part 5/Data Files/triple_sim_gpu3.pkl'
with open(file_path, 'wb') as triple_sim_gpu_file:
    pickle.dump(triple_sim_gpu, triple_sim_gpu_file)

print('Done GPU simulation')

# print('Shooting CPU simulation...')
# triple_sim_cpu = IRSC.IRSCaustics(annulus_param_dict=triple_lens_parameters_1)
# triple_sim_cpu.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

# file_path = './Unity Analysis/Part 5/Data Files/triple_sim_cpu.pkl'
# with open(file_path, 'wb') as triple_sim_cpu_file:
#     pickle.dump(triple_sim_cpu, triple_sim_cpu_file)

# print('Done CPU simulation')