# Checking what simulations did not get run
import os
import numexpr

ss = ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1/0.9', '1/0.8', '1/0.7', '1/0.6', '1/0.5', '1/0.4', '1/0.3', '1/0.2']
qs = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
ss = [numexpr.evaluate(s).item() for s in ss]

sim_dir = './Unity/Simulations/Part_4_Binary_Collection/'

for s in ss:
    for q in qs:
        file_name = f'binary_{q:.0e}_{s:.2e}.pkl'
        file_path = sim_dir + file_name
        if os.path.exists(file_path):
            print(f'File {file_name} exists.')
        else:
            print(f'File {file_name} does NOT exist.')