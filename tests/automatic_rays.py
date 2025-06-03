import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('.')

from IRSMicroLensing import IRSCaustics as IRSC

def get_class_size(obj):
    total_size = sys.getsizeof(obj)
    for attribute_name in vars(obj):
        attribute = getattr(obj, attribute_name)
        total_size += sys.getsizeof(attribute)
    return total_size

plt.switch_backend('agg')  # Use a non-interactive backend for testing

lens_att = [
    [0, 0, 0.01, 1],
    [0.8, 0, 0.01, 0.001]
]

num_theta = 3990
num_r = 4 * num_theta

param_dict_1 = {'pixels': 1000, 'ang_width': 'auto', 'lens_att': lens_att, 'thickness': 'auto', 'num_r': num_r, 'num_theta': num_theta}

mag_map_1 = IRSC.IRSCaustics(annulus_param_dict=param_dict_1)
print(mag_map_1.param_dict)
print(mag_map_1.y_plus, mag_map_1.y_minus)

magnifications = mag_map_1.plot(save_plot=True)

print(np.isinf(mag_map_1.magnifications_log).any())

print(get_class_size(mag_map_1) / 1e6)
print(mag_map_1.X.nbytes / 1e6)