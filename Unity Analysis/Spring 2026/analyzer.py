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
    ss_str = ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1/0.9', '1/0.8', '1/0.7', '1/0.6', '1/0.5', '1/0.4', '1/0.3', '1/0.2']
    ss = [numexpr.evaluate(s) for s in ss_str]
    qs = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3]

    binary_directory = './Unity/Simulations/Binary_Collection'
    single_directory = './Unity/Simulations/Single_Collection'

    #%% Reading in attributes from csv file
    single_attributes = pd.read_csv('./Unity Analysis/Spring 2026/single_attributes.csv')
    binary_attributes = pd.read_csv('./Unity Analysis/Spring 2026/binary_attributes.csv')

    #%% Defining variables to be stored
    thresholds = np.logspace(-8, 0, num=100)
    headers = ['s', 'q', 'source_radius', 'LD_coeff'] + [f'{thresh:.1e}' for thresh in thresholds]
    output = []

    #%% Ensuring folders are created for output
    if not os.path.exists(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/'):
        os.makedirs(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/')
    
    for q in qs:
        if not os.path.exists(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/Frac_Maps_{q:.0e}/'):
            os.makedirs(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/Frac_Maps_{q:.0e}/')

    #%% Looping through parameter space
    with open(f'./Unity Analysis/Spring 2026/Data Files/analysis_output_{multiplier:.0e}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)

        for i, q in enumerate(qs):
            for j, s in enumerate(ss):
                print('=' * 50)
                print(f'Processing s={ss_str[j]}, q={q:.0e}...')

                #%% Reading in both single and binary lens simulations
                if s < 1.0:
                    single_file = f'{single_directory}/Single_Collection_{q:.0e}/single_{q:.0e}_{s:.2e}.pkl'
                elif s == 1.0:
                    single_file = f'{single_directory}/Single_Collection_{q:.0e}/single_{q:.0e}_9.00e-01.pkl'
                else:
                    single_file = f'{single_directory}/Single_Collection_{q:.0e}/single_{q:.0e}_{1/(s):.2e}.pkl'
                
                single_sim = IRSC.caustic_reader(single_file)

                binary_file = f'{binary_directory}/Binary_Collection_{q:.0e}/binary_{q:.0e}_{s:.2e}.pkl'
                binary_sim = IRSC.caustic_reader(binary_file)

                if binary_sim.pixels != single_sim.pixels:
                    print('  Pixel mismatch between single and binary simulations. Skipping...')
                    continue

                #%% Extracting important parameters from simulations
                pixels = single_sim.pixels
                ang_width = single_sim.ang_width
                ang_res = single_sim.ang_res

                single_magnifications = single_sim.magnifications
                binary_magnifications = binary_sim.magnifications

                min_source_radius = binary_attributes.loc[(binary_attributes['s'] == ss_str[j]) & (binary_attributes['q'] == q), 'min_source_radius'].values[0]
                max_source_radius = binary_attributes.loc[(binary_attributes['s'] == ss_str[j]) & (binary_attributes['q'] == q), 'max_source_radius'].values[0]

                X_pix, Y_pix = np.meshgrid(np.linspace(-ang_width/2, ang_width/2, pixels), np.linspace(-ang_width/2, ang_width/2, pixels))

                #%% Defining source profile
                print('  Defining source profile...')
                radius = max_source_radius / 10 * multiplier
                LD_coeff = 0.5
                print('-' * 50)
                source_profile = IRSF.IRSFunctions.source_profile(ang_res=single_sim.ang_res, profile_type='LD', rad=radius, LD=LD_coeff)
                print('-' * 50)

                #%% Convolving magnification maps with source profile
                print('  Convolving single lens magnification map...')
                single_conv_mags = single_sim.convolve(source_profile=source_profile)
                print('  Convolving binary lens magnification map...')
                binary_conv_mags = binary_sim.convolve(source_profile=source_profile)

                #%% Taking fractional difference between maps
                print('  Calculating fractional difference map...')
                frac_diff_map = (binary_conv_mags - single_conv_mags) / single_conv_mags

                print('  Saving fractional difference map...')
                np.save(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/Frac_Maps_{q:.0e}/frac_diff_map_q{q:.0e}_s{s:.2e}.npy', frac_diff_map)

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
                # print('  Calculating area above threshold...')
                # fractional_areas_above_threshold = []

                # for threshold in thresholds:
                #     above_threshold = np.abs(frac_diff_map) > threshold
                #     total_area_above_threshold = np.sum(above_threshold)
                #     fractional_area_above_threshold = total_area_above_threshold / (pixels**2)
                #     fractional_areas_above_threshold.append(fractional_area_above_threshold)

                # fractional_areas_above_threshold = np.array(fractional_areas_above_threshold)
                
                # #%% Saving fractional deviations to csv
                # output.append([ss_str[i], q, radius, LD_coeff] + fractional_areas_above_threshold.tolist())
                # writer.writerow(output[-1])