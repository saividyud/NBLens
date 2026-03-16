# Checking what simulations did not get run
import os
import numexpr

# ── Parameter grid (5×5×5 = 125 simulations) ──
SS1 = ['0.2', '0.6', '1.0', '1/0.6', '1/0.2']
SS2 = ['0.2', '0.6', '1.0', '1/0.6', '1/0.2']
ALPHAS = [0.0, 45.0, 90.0, 135.0, 180.0]

sim_dir = './Unity/Simulations/Part_5_Triple_Collection/'

for s1 in SS1:
    for s2 in SS2:
        for alpha in ALPHAS:
            file_name = f'triple_{numexpr.evaluate(s1).item():.2e}_{numexpr.evaluate(s2).item():.2e}_{alpha:.0f}.pkl'
            file_path = sim_dir + file_name
            if os.path.exists(file_path):
                print(f'File {file_name} exists.')
            else:
                print(f'File {file_name} does NOT exist.')