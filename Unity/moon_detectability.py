import psutil
memory = psutil.virtual_memory()
print(f'Total memory: {memory.total / (1024**3)}')
print(f'Available memory: {memory.available / (1024**3)}')
print(f'Free memory: {memory.free / (1024**3)}')
print(f'Inactive memory: {memory.inactive / (1024**3)}')
print(f'Used memory: {memory.used / (1024**3)}')
print()

import sys
sys.path.append('.')

import time as t

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')

from IRSMicroLensing import IRSCaustics as IRSC
from IRSMicroLensing import IRSFunctions as IRSF

def get_class_size(obj: object):
    total_memory = sys.getsizeof(obj)

    for var, val in vars(obj).items():
        # print(f'{var}: {sys.getsizeof(val) / (1024 ** 3)}')
        total_memory += sys.getsizeof(val)

    return total_memory / (1024 ** 3)

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
num_theta = 20000
num_r = 4 * num_theta

single_lens_attributes = [
    [0, 0, 0.001, 1]
]

triple_lens_attributes = [
    [0, 0, 0.001, 1],
    [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), 0.001, q1],
    [s2*np.cos(np.deg2rad(alpha2)), s2*np.sin(np.deg2rad(alpha2)), 0.001, q2]
]

ang_width, thickness, (y_plus, y_minus), cusp_points = IRSC.IRSCaustics.ang_width_thickness_calculator(triple_lens_attributes)

print(f'Angular width: {ang_width}')
print(f'Thickness: {thickness}')
print(f'Annulus lower bound: {y_plus}')
print(f'Annulus upper bound: {y_minus}')
print()

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

''' Simulating single lens magnification map '''
single_lens = IRSC.IRSCaustics(annulus_param_dict=single_lens_parameters)
single_lens_magnifications = single_lens.plot(cm_offset='auto')
single_lens.fig_c.savefig('./Unity/Single Lens Magnifications.png', dpi=300)

print()
print(f'Single lens class memory: {get_class_size(single_lens)}')
memory = psutil.virtual_memory()
print(f'Available memory: {memory.available / (1024**3)}')
print(f'Free memory: {memory.free / (1024**3)}')
print(f'Inactive memory: {memory.inactive / (1024**3)}')
print(f'Used memory: {memory.used / (1024**3)}')
print()

''' Similating triple lens magnification map '''
triple_lens = IRSC.IRSCaustics(annulus_param_dict=triple_lens_parameters)
triple_lens_magnifications = triple_lens.plot(cm_offset='auto')
triple_lens.fig_c.savefig('./Unity/Triple Lens Magnifications.png', dpi=300)

print()
print(f'Triple lens class memory: {get_class_size(triple_lens)}')
memory = psutil.virtual_memory()
print(f'Available memory: {memory.available / (1024**3)}')
print(f'Free memory: {memory.free / (1024**3)}')
print(f'Inactive memory: {memory.inactive / (1024**3)}')
print(f'Used memory: {memory.used / (1024**3)}')
print()

''' Defining source profile '''
radius = 3e-3
LD = 0.5

source_profile = IRSF.IRSFunctions.source_profile(ang_res=single_lens.param_dict['ang_res'], rad=radius, profile_type='LD', LD=LD)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot()

img = ax.imshow(source_profile, cmap='afmhot', extent=[-radius, radius, -radius, radius])
plt.colorbar(img, ax=ax, label='Source Brightness')

ax.set_xlabel('X [$\\theta_E$]')
ax.set_ylabel('Y [$\\theta_E$]')
ax.set_title('Source Profile')

ax.set_aspect('equal')

fig.savefig('./Unity/Source Profile.png', dpi=300)

''' Convolving source profile with single lens magnification maps '''
convolved_single_lens = single_lens.convolve(source_profile=source_profile)

''' Convolving source profile with triple lens magnification maps '''
convolved_triple_lens = triple_lens.convolve(source_profile=source_profile)

''' Plotting convolved brightnesses '''
fig = plt.figure(figsize=(14, 8))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

ax1.imshow(convolved_single_lens, cmap='gray', extent=[-ang_width/2, ang_width/2, -ang_width/2, ang_width/2])
ax2.imshow(convolved_triple_lens, cmap='gray', extent=[-ang_width/2, ang_width/2, -ang_width/2, ang_width/2])

ax1.set_xlabel('X [$\\theta_E$]')
ax1.set_ylabel('Y [$\\theta_E$]')
ax1.set_title('Single Lens Convolved Brightnesses')
ax1.set_aspect('equal')

ax2.set_xlabel('X [$\\theta_E$]')
ax2.set_ylabel('Y [$\\theta_E$]')
ax2.set_title('Triple Lens Convolved Brightnesses')
ax2.set_aspect('equal')

fig.savefig('./Unity/Convolved Brightnesses.png', dpi=300)

''' Fractional deviations '''
fractional_deviations = (convolved_triple_lens - convolved_single_lens) / convolved_single_lens

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot()

# img = ax.imshow(fractional_deviations, extent=[-ang_width/2, ang_width/2, -ang_width/2, ang_width/2], cmap='grey')

img = plt.contour(single_lens.X_pix, single_lens.Y_pix, fractional_deviations,
                   levels=[-0.30, -0.10, -0.03, -0.01, 0, 0.01, 0.03, 0.10, 0.30],
                   colors=['blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red']
)

plt.colorbar(img)

ax.set_aspect('equal')
ax.set_xlim(-ang_width/2, ang_width/2)
ax.set_ylim(-ang_width/2, ang_width/2)

fig.savefig('./Unity/Contours.png')

init_time = t.time()
np.save('./Unity/Analysis 6-4/single_lens_magnifications.npy', single_lens_magnifications)
np.save('./Unity/Analysis 6-4/triple_lens_magnifications.npy', triple_lens_magnifications)
np.save('./Unity/Analysis 6-4/convolved_single_lens.npy', convolved_single_lens)
np.save('./Unity/Analysis 6-4/convolved_triple_lens.npy', convolved_triple_lens)
np.save('./Unity/Analysis 6-4/fractional_deviations.npy', fractional_deviations)
print(f'Saving arrays to file: {(t.time() - init_time):.2}')