from IRSMicroLensing import IRSCaustics as IRSC
import matplotlib.pyplot as plt

lens_attributes = [
    [0, 0, 0.1, 1], # Sun
    [0.25, 0, 0.01, 1e-6] # Earth
]

param_dict = {'lens_att': lens_attributes, 'pixels': 5000, 'num_r': 12000, 'num_theta': 3000, 'ang_width': 'auto', 'thickness': 'auto'}

simulation = IRSC.IRSCaustics(annulus_param_dict=param_dict)
magnifications = simulation.plot(show_lenses=False, file_save=False)

plt.show()