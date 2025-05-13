import numpy as np
import numba as nb
import scipy as sci
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

import time as t

from IRSMicroLensing import IRSCaustics as IRSC

''' Computing magnification map '''
lens_att = [
    [0, 0, 0.01, 1],
    [0.3, 0.1, 0.01, 0.01],
    [0.1, -0.2, 0.01, 0.01]
]

param_dict = {'pixels': 5000, 'ang_width': 'auto', 'lens_att': lens_att, 'thickness': 'auto', 'num_r': 12000, 'num_theta': 3000}

caus = IRSC.IRSCaustics(annulus_param_dict=param_dict)
magnifications = caus.plot(show_lenses=False)

ax = caus.fig_c.gca()

for group in caus.points:
    ax.scatter(group[:, 0], group[:, 1])

print(caus.param_dict)

''' Computing Gaussian profile '''
ang_res = caus.param_dict['ang_res']

R = 1e-3

pix_width = 1*R / ang_res

x = np.linspace(-5*R, 5*R, int(pix_width))
y = np.linspace(-5*R, 5*R, int(pix_width))

X, Y = np.meshgrid(x, y)

r2 = X**2 + Y**2

source_profile_gauss = 1.5 * np.exp(-r2*0.5 / R**2)

''' Computing limb-darkened profiles '''
ang_res = caus.param_dict['ang_res']

R = 1e-3

pix_width = 1*R / ang_res

ang_width = 1.5*R

x = np.linspace(-ang_width, ang_width, int(pix_width))
y = np.linspace(-ang_width, ang_width, int(pix_width))

X, Y = np.meshgrid(x, y)

r2 = X**2 + Y**2

# Discretizing limb-darkening coefficient
LD_coeffs = np.linspace(0, 1, 9)

# Defining source profile brightnesses
cos_theta = np.nan_to_num(np.sqrt(1 - (r2 / R**2)), nan=0)

source_profiles = np.zeros(shape=(9, X.shape[0], X.shape[1]))

for i, LD_coeff in enumerate(LD_coeffs):
    source_profiles[i, :, :] = (1 - (LD_coeff * (1 - (3/2 * cos_theta)))) * (r2 <= R**2)

''' Plotting all limb-darkened source profiles '''
# fig = plt.figure(figsize=(10, 8))
# gs = gridspec.GridSpec(3, 3, figure=fig, wspace=0.01)
# axes = fig.subplots(4, 3)

# for i, source_profile in enumerate(source_profiles):
#     if i < 3:
#         ax = axes[0, i]
#     elif i < 6:
#         ax = axes[1, i-3]
#     else:
#         ax = axes[2, i-6]
    
#     plot = ax.imshow(source_profile, cmap='afmhot', extent=[-ang_width, ang_width, -ang_width, ang_width], vmin=0, vmax=1.5)
#     ax.axis('off')

#     props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
#     ax.text(0.02, 0.98, f'$\Gamma_\lambda = {LD_coeffs[i]}$', va='top', zorder=10, bbox=props, transform=ax.transAxes)

# for i in range(0, 3):
#     axes[3, i].axis('off')

# cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
# cbar_ax.axis('off')

# axes[3, 1].imshow(source_profile_gauss, cmap='afmhot', extent=[-ang_width, ang_width, -ang_width, ang_width], vmin=0, vmax=1.5)
# props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
# axes[3, 1].text(0.02, 0.98, f'Gaussian Profile', va='top', zorder=10, bbox=props, transform=axes[3, 1].transAxes)

# bar = fig.colorbar(plot, ax=cbar_ax, location='right', fraction=1.5, aspect=20)
# bar.set_label('Normalized Brightness')

# fig.suptitle('Source Profile Brightnesses with Different LD Coefficients', y=0.92)

# plt.subplots_adjust(wspace=0.01, hspace=0.05)

''' Convolving all source profiles with magnification map '''
convolved_magnifications_log = []

for i, source_profile in enumerate(source_profiles):
    i_time = t.time()
    ans = sci.signal.fftconvolve(magnifications, source_profile, 'same')
    print(i, t.time() - i_time)

    magnifications_log = np.where(ans == 0, 0.1, ans)
    magnifications_log = np.log10(magnifications_log)
    magnifications_log = np.where(magnifications_log == -1, 0, magnifications_log)

    convolved_magnifications_log.append(magnifications_log)

''' Plotting all convolved magnification maps '''
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(3, 3, figure=fig, wspace=0.01)
axes = fig.subplots(4, 3)

for i, mag_map in enumerate(convolved_magnifications_log):
    if i < 3:
        ax = axes[0, i]
    elif i < 6:
        ax = axes[1, i-3]
    else:
        ax = axes[2, i-6]
    
    plot = ax.imshow(mag_map, cmap='gray', extent=[-caus.ang_width/2, caus.ang_width/2, -caus.ang_width/2, caus.ang_width/2], vmin=0, vmax=10)
    ax.axis('off')

    props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
    ax.text(0.02, 0.98, f'$\Gamma_\lambda = {LD_coeffs[i]}$', va='top', zorder=10, bbox=props, transform=ax.transAxes)

for i in range(0, 3):
    axes[3, i].axis('off')

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
cbar_ax.axis('off')

axes[3, 1].imshow(caus.magnifications_log, cmap='gray', extent=[-caus.ang_width/2, caus.ang_width/2, -caus.ang_width/2, caus.ang_width/2], vmin=0, vmax=10)
props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
axes[3, 1].text(0.02, 0.98, f'Unconvolved Mag Map', va='top', zorder=10, bbox=props, transform=axes[3, 1].transAxes)

bar = fig.colorbar(plot, ax=cbar_ax, location='right', fraction=1.5, aspect=20)
bar.set_label('$log_{10}$ Magnifications')

fig.suptitle('Convolved Magnification Maps', y=0.92)

# fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, hspace=0.02)
plt.subplots_adjust(wspace=0.01, hspace=0.05)

plt.show()