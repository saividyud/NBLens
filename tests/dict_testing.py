import numpy as np

lens_att = [
    [0, 0, 0.01, 1],
    [0.4, 0.2, 0.01, 0.001],
    [0, -0.5, 0.01, 0.001]
]

param_dict = {'pixels': 10000, 'ang_width': 'auto', 'lens_att': lens_att, 'thickness': 'auto', 'num_r': 12000, 'num_theta': 3000}

thickness = 0.005

if 'thickness' in param_dict.keys():
    if param_dict['thickness'] == 'auto':
        param_dict['thickness'] = thickness
    else:
        pass

print(param_dict)