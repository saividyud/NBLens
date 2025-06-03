from IRSMicroLensing import IRSCaustics as IRSC

lens_att = [
    [0, 0, 0.01, 1],
    [0.8, 0, 0.01, 0.001]
]

num_theta = 5000
num_r = 2 * num_theta

binary_param_dict = {'pixels': 1000, 'ang_width': 'auto', 'lens_att': lens_att, 'thickness': 'auto', 'num_r': num_r, 'num_theta': num_theta}

binary_mag_map = IRSC.IRSCaustics(annulus_param_dict=binary_param_dict)

# print(binary_mag_map.param_dict)
# for i in range(10):
binary_magnifications = binary_mag_map.plot(show_lenses=False, show_plot=False, cm_offset=[0, 0], show_dev=False, save_plot=False, print_stats=False)

# print(binary_mag_map.y_plus, binary_mag_map.y_minus)

# for caustic in binary_mag_map.caustic_cusps:
#     binary_mag_map.ax_c.scatter(caustic[:, 0], caustic[:, 1], s=10, c='red')

# plt.show()