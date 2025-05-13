import time as t

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as patches

# Timing imports
# start = t.time()
from IRSMicroLensing import IRSCaustics as IRSC
# print(f"Import time: {t.time() - start:.2f} seconds")

# start = t.time()
from IRSMicroLensing import IRSFunctions as IRSF
# print(f"Import time: {t.time() - start:.2f} seconds")

# Colormap
cmap = 'gray'

# Defining binary lens parameters
binary_lens_attributes = [
    [0, 0, 0.01, 1],
    [0.8, 0, 0.01, 1e-3]
]

# Creating a list of all lenses
complete_lens_attributes = binary_lens_attributes.copy()

for i in range(3, 10):
    complete_lens_attributes.append(
        [0.8, 0.001, 0.01, 10**(-i)] # Small mass
    )

# Finding best angular width and thickness
ang_width, thickness, caustic_cusps = IRSC.IRSCaustics.ang_width_thickness_calculator(lens_att=complete_lens_attributes)
print()
print(f"Angular width: {ang_width}")
print(f"Thickness: {thickness}")

# Probing mass space to see what the smallest mass detectable is
trinary_lens_attributes_list = []

for i in range(3, 10):
    trinary = binary_lens_attributes.copy()
    trinary.append(
        [0.8, 0.001, 0.01, 10**(-i)] # Small mass
    )

    trinary_lens_attributes_list.append(trinary)

# Creating instance of IRSCaustics for binary lens
binary_param_dict = {'pixels': 5000, 'ang_width': ang_width, 'lens_att': binary_lens_attributes, 'thickness': thickness, 'num_r': 12000, 'num_theta': 3000}
binary_mag_map = IRSC.IRSCaustics(annulus_param_dict=binary_param_dict)
print()
print('Calculating binary lens magnification map')
binary_magnifications = binary_mag_map.plot(show_lenses=False, cmap=cmap, show_plot=False)

# Defining source profile for convolution
source_profile = IRSF.IRSFunctions.source_profile(ang_res=binary_mag_map.param_dict['ang_res'], rad=1e-3, profile_type='LD', LD=0)
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot()
# ax.imshow(source_profile, cmap='afmhot')

# Convolving binary magnification map with the source profile
convolved_binary_brightnesses = binary_mag_map.convolve(source_profile=source_profile)

# Defining list of instances of IRSCaustics for trinary lens
trinary_mag_map_list = []

# Creating a list of convolved brightnesses for each trinary lens
convolved_trinary_brightnesses_list = []

# List of convolved brightness difference for each trinary lens
convolved_brightness_differences_list = []

# Main for loop for each trinary lens configuration
for i in range(len(trinary_lens_attributes_list)):
    # Creating instance of IRSCaustics for trinary lens
    trinary_param_dict = binary_param_dict.copy()
    trinary_param_dict['lens_att'] = trinary_lens_attributes_list[i]
    trinary_mag_map_list.append(IRSC.IRSCaustics(annulus_param_dict=trinary_param_dict))

    print()
    print(f"Calculating trinary lens magnification map {i+1}/{len(trinary_lens_attributes_list)}")
    trinary_magnifications = trinary_mag_map_list[i].plot(show_lenses=False, cmap=cmap, show_plot=False)

    # Convolving each magnification map with the source profile
    convolved_trinary_brightnesses_list.append(trinary_mag_map_list[i].convolve(source_profile=source_profile))

    # Plotting difference between convolved magnification maps
    convolved_brightness_differences_list.append(np.abs((convolved_trinary_brightnesses_list[i] - convolved_binary_brightnesses) / convolved_binary_brightnesses))


fig = plt.figure(figsize=(16, 8))
axes = fig.subplots(1, 8)

for i in range(len(trinary_lens_attributes_list)):
    ax = axes[i]
    ax.set_xlabel('x [theta_E]')
    ax.set_ylabel('y [theta_E]')
    ax.set_xlim(-binary_mag_map.ang_width/2, binary_mag_map.ang_width/2)
    ax.set_ylim(-binary_mag_map.ang_width/2, binary_mag_map.ang_width/2)

    if i == 0:
        ax.set_title('Binary Lens')
        img = ax.imshow(convolved_binary_brightnesses, cmap=cmap, extent=[-binary_mag_map.ang_width/2, binary_mag_map.ang_width/2, -binary_mag_map.ang_width/2, binary_mag_map.ang_width/2])

    else:
        ax.set_title(f"Trinary Lens {i}")

        img = ax.imshow(convolved_brightness_differences_list[i-1], cmap=cmap, extent=[-binary_mag_map.ang_width/2, binary_mag_map.ang_width/2, -binary_mag_map.ang_width/2, binary_mag_map.ang_width/2])

plt.show()