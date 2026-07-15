if __name__ == '__main__':

    #%% Importing modules
    import sys
    sys.path.append('.')

    import numpy as np
    import pandas as pd
    import numexpr

    import IRSMicroLensing.IRSFunctions as IRSF
    import IRSMicroLensing.IRSCaustics as IRSC

    import csv
    import os

    import argparse

    #%% Setting plot parameters
    # plt.rcParams['font.family'] = 'Times New Roman'
    # plt.rcParams['figure.titlesize'] = 20
    # plt.rcParams['figure.titleweight'] = 'bold'
    # plt.rcParams['figure.figsize'] = (10, 8)
    # plt.rcParams['axes.titlesize'] = 16
    # plt.rcParams['axes.labelsize'] = 14
    # plt.rcParams['figure.labelsize'] = 14
    # plt.rcParams['xtick.labelsize'] = 12
    # plt.rcParams['ytick.labelsize'] = 12
    # plt.rcParams['legend.fontsize'] = 12
    # plt.rcParams['mathtext.fontset'] = 'cm'

    #%% Defining some useful functions
    def log_array(min_pow, max_pow):
        powers = np.arange(min_pow, max_pow + 1)
        base_values = np.array([1, 3])
        
        # Use broadcasting to create the sequence
        pos_arr = np.concatenate([base_values * (10.0**p) for p in powers])

        neg_arr = -pos_arr[::-1]

        return np.concatenate([neg_arr, np.array([0]), pos_arr])

    def colors_by_log(arr):
        colors = []
        for val in arr:
            if val < 0:
                colors.append('blue')
            elif val == 0:
                colors.append('green')
            elif val > 0:
                colors.append('red')
        return colors

    #%% Defining arguments
    parser = argparse.ArgumentParser(description='Analyze microlensing simulations for fractional area above thresholds.')
    parser.add_argument('--multiplier', type=float, default=1.0, help='Multiplier for output filename.')
    args = parser.parse_args()
    multiplier = args.multiplier

    #%% Defining parameter space
    # ── Parameter grid (1×6×2×5 = 60 simulations) ──
    # s1 is fixed at 1.2; s2, q2, and alpha vary.
    S1_str = '1.2'
    SS2_str = ['0.2', '0.5', '0.8', '1/0.8', '1/0.5', '1/0.2']
    Q2S = [1e-4, 3e-4]
    ALPHAS = [0.0, 45.0, 90.0, 135.0, 180.0]
    TOTAL_SIMS = len(SS2_str) * len(Q2S) * len(ALPHAS)  # 60

    s1 = np.float64(numexpr.evaluate(S1_str).item())
    ss2 = [numexpr.evaluate(s) for s in SS2_str]
    alphas = [np.float64(numexpr.evaluate(str(alpha)).item()) for alpha in ALPHAS]

    q1 = 1e-3

    triple_directory = './Unity/Simulations/Part_7_Triple_Collection'
    binary_directory = './Unity/Simulations/Part_7_Binary_Collection'

    #%% Defining variables to be stored
    thresholds = np.logspace(-8, 0, num=100)
    headers = ['s1', 's2', 'q2', 'alpha', 'source_radius', 'LD_coeff'] + [f'{thresh:.1e}' for thresh in thresholds]
    output = []

    #%% Reading in binary lens simulation (s1 and q1 are fixed, so this reference is constant)
    binary_file = f'{binary_directory}/binary_{q1:.0e}_{s1:.2e}.pkl'
    binary_sim = IRSC.caustic_reader(binary_file)

    #%% Looping through parameter space
    output_directory = './Unity Analysis/Part 7/Data Files'
    os.makedirs(output_directory, exist_ok=True)

    with open(f'{output_directory}/analysis_output_{multiplier:.0e}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)

        for j, s2 in enumerate(ss2):
            for l, q2 in enumerate(Q2S):
                for k, alpha in enumerate(alphas):
                    print('=' * 50)
                    print(f'Processing s1={s1}, s2={s2}, q2={q2}, alpha={alpha}...')

                    #%% Reading in triple lens simulation
                    triple_file = f'{triple_directory}/triple_{s2:.2e}_{q2:.2e}_{alpha:.0f}.pkl'
                    triple_sim = IRSC.caustic_reader(triple_file)

                    #%% Extracting important parameters from simulations
                    pixels = binary_sim.pixels
                    ang_width = binary_sim.ang_width
                    ang_res = binary_sim.ang_res

                    binary_magnifications = binary_sim.magnifications
                    triple_magnifications = triple_sim.magnifications

                    min_source_radius = 0.001/2
                    max_source_radius = 0.1/2

                    X_pix, Y_pix = np.meshgrid(np.linspace(-ang_width/2, ang_width/2, pixels), np.linspace(-ang_width/2, ang_width/2, pixels))

                    #%% Defining source profile
                    print('  Defining source profile...')
                    radius = max_source_radius / 10 * multiplier
                    LD_coeff = 0.5
                    print('-' * 50)
                    source_profile = IRSF.IRSFunctions.source_profile(ang_res=binary_sim.ang_res, profile_type='LD', rad=radius, LD=LD_coeff)
                    print('-' * 50)

                    #%% Convolving magnification maps with source profile
                    print('  Convolving single lens magnification map...')
                    binary_conv_mags = binary_sim.convolve(source_profile=source_profile)
                    print('  Convolving binary lens magnification map...')
                    triple_conv_mags = triple_sim.convolve(source_profile=source_profile)

                    #%% Taking fractional difference between maps
                    print('  Calculating fractional difference map...')
                    frac_diff_map = (triple_conv_mags - binary_conv_mags) / binary_conv_mags

                    # print('  Saving fractional difference map...')
                    # np.save(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/Frac_Maps_{q:.0e}/frac_diff_map_q{q:.0e}_s{s:.2e}.npy', frac_diff_map)

                    #%% Plotting fractional difference map
                    # print('  Plotting fractional difference map...')
                    # fig = plt.figure()
                    # ax = fig.add_subplot()
                    # img = ax.contour(
                    #     X_pix, Y_pix, np.flip(frac_diff_map, axis=0),
                    #     levels=log_array(-3, -1).tolist(),
                    #     colors=colors_by_log(log_array(-3, -1)),
                    #     linewidths=0.5
                    # )
                    
                    # fig.colorbar(img, label='Fractional Difference')
                    # ax.set_title(f'Fractional Difference Map ($s={ss_str[i]}, q={q:.0e}, r={radius:.2e} \\theta_E$)')
                    # ax.set_xlabel('X [$\\theta_E$]')
                    # ax.set_ylabel('Y [$\\theta_E$]')

                    # plt.savefig(f'./Unity Analysis/Spring 2026/Figures/frac_diff_map_q{q:.0e}_s{s:.2e}_r{radius:.2e}.png', dpi=100)

                    #%% Calculating fractional area of fractional map above threshold
                    print('  Calculating area above threshold...')
                    fractional_area_above_thresholds = []
                    # actual_area_above_thresholds = []

                    for threshold in thresholds:
                        above_threshold = np.abs(frac_diff_map) > threshold
                        total_area_above_threshold = np.sum(above_threshold)
                        fractional_area_above_threshold = total_area_above_threshold / (pixels**2)
                        fractional_area_above_thresholds.append(fractional_area_above_threshold)
                        # actual_area_above_thresholds.append(fractional_area_above_threshold * (ang_width**2))

                    # actual_area_above_thresholds = np.array(actual_area_above_thresholds)
                    fractional_area_above_thresholds = np.array(fractional_area_above_thresholds)
                    
                    #%% Saving fractional deviations to csv
                    output.append([S1_str, SS2_str[j], Q2S[l], ALPHAS[k], radius, LD_coeff] + fractional_area_above_thresholds.tolist())
                    writer.writerow(output[-1])
