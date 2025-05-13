from IRSMicroLensing import IRSCaustics_annulus as IRSC
import matplotlib.pyplot as plt

# lens_att = [
#     [-0.5, 0, 0.01, 1],
#     [0.5, 0, 0.01, 1]
# ]

# param_dict = {'pixels': 1000, 'rays_per_pixel': 4, 'ang_width': 5, 'lens_att': lens_att}

# caus = IRSC.IRSCaustics(whole_plane_param_dict=param_dict)
# mag_plot = caus.plot()

# print(caus.param_dict)

# plt.show()

# lens_att = [
#     [0, 0, 0.01, 1],
#     [0.3, 0, 0.01, 0.1]
# ]

# param_dict = {'pixels': 2000, 'rays_per_pixel': 4, 'ang_width': 3, 'lens_att': lens_att}

# caus = IRSC.IRSCaustics(whole_plane_param_dict=param_dict)
# mag_plot = caus.plot()

# print(caus.param_dict)

# plt.show()

''' Earth Example '''
lens_att = [
    [0, 0, 0.01, 1],
    [0.25, 0, 0.01, 1e-6],
    [0.25, 0.01, 0.01, 1e-8]
]

param_dict = {'pixels': 5000, 'ang_width': 1e-6, 'lens_att': lens_att, 'thickness': 1e-6, 'num_r': 12000, 'num_theta': 3000}

caus = IRSC.IRSCaustics(annulus_param_dict=param_dict)
# mag_plot = caus.plot(zoom=(1e-6, 1e-6))
mag_plot = caus.plot(file_save=False)

print(caus.param_dict)

plt.show()

''' Jupiter Example '''
# lens_att = [
#     [0, 0, 0.01, 1],
#     [0.6, 0, 0.01, 0.001]
# ]

# param_dict = {'pixels': 3000, 'ang_width': 0.01, 'lens_att': lens_att, 'thickness': 0.01, 'num_r': 12000, 'num_theta': 3000}
# # param_dict = {'pixels': 1000, 'ang_width': 4, 'rays_per_pixel': 4, 'lens_att': lens_att}

# caus = IRSC.IRSCaustics(annulus_param_dict=param_dict)
# # caus = IRSC.IRSCaustics(whole_plane_param_dict=param_dict)
# mag_plot = caus.plot(show_lenses=False)

# print(caus.param_dict)

# plt.show()