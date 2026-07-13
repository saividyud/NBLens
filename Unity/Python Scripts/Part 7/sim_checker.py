# Checking what simulations did not get run
import os
import numexpr

# ── Parameter grid (1×6×2×5 = 60 simulations) ──
# s1 is fixed at 1.2; s2, q2, and alpha vary.
S1 = '1.2'
SS2 = ['0.2', '0.5', '0.8', '1/0.8', '1/0.5', '1/0.2']
Q2S = [1e-4, 3e-4]
ALPHAS = [0.0, 45.0, 90.0, 135.0, 180.0]

sim_dir = './Unity/Simulations/Part_7_Triple_Collection/'

for s2 in SS2:
    for q2 in Q2S:
        for alpha in ALPHAS:
            file_name = f'triple_{numexpr.evaluate(s2).item():.2e}_{q2:.2e}_{alpha:.0f}.pkl'
            running_file_name = file_name + '.running'
            file_path = sim_dir + file_name
            running_file_path = sim_dir + running_file_name
            if os.path.exists(file_path):
                print(f'File {file_name} exists.')
            elif os.path.exists(running_file_path):
                print(f'File {file_name} is running.')
            else:
                print(f'File {file_name} does NOT exist.')
