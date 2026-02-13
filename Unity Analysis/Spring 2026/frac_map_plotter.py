if __name__ == '__main__':

    #%% Importing modules
    import sys
    sys.path.append('.')

    import numpy as np
    import numexpr
    import pandas as pd

    import IRSMicroLensing.IRSCaustics as IRSC
    import IRSMicroLensing.IRSFunctions as IRSF

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    import argparse
    import sys
    import os

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['figure.titlesize'] = 20
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 20
    plt.rcParams['figure.labelsize'] = 20
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.title_fontsize'] = 16
    plt.rcParams['legend.fontsize'] = 16
    plt.rcParams['mathtext.fontset'] = 'cm'

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
    # parser = argparse.ArgumentParser(description='Plot fractional area maps for microlensing simulations.')
    # parser.add_argument('--mass_ratio', type=str, help='Mass ratio for output filename.')
    # parser.add_argument('--separation', type=str, help='Separation value for filtering parameter space.')
    # args = parser.parse_args()
    # mass_ratio = args.mass_ratio
    # separation = args.separation
    # s = numexpr.evaluate(str(separation))
    # q = numexpr.evaluate(str(mass_ratio))

    # print(f'Arguments parsed successfully. Starting analysis with mass_ratio = {mass_ratio} and separation = {separation}...')
    
    #%% Defining parameter space
    ss_str = ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1/0.9', '1/0.8', '1/0.7', '1/0.6', '1/0.5', '1/0.4', '1/0.3', '1/0.2']
    ss = np.array([numexpr.evaluate(s) for s in ss_str])
    qs = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
    multipliers = [1e-1, 3e-1, 1e0, 3e0, 1e1]
    
    binary_directory = './Unity/Simulations/Binary_Collection'
    single_directory = './Unity/Simulations/Single_Collection'

    #%% Reading in attributes from csv file
    single_attributes = pd.read_csv('./Unity Analysis/Spring 2026/single_attributes.csv')
    binary_attributes = pd.read_csv('./Unity Analysis/Spring 2026/binary_attributes.csv')

    print('Parameter space defined successfully. Starting to read fractional area maps...')

    #%% Ensuring folders are created for output
    # if not os.path.exists(f'./Unity Analysis/Spring 2026/Figures/Mag_Maps/'):
    #     os.makedirs(f'./Unity Analysis/Spring 2026/Figures/Mag_Maps/')

    # if not os.path.exists(f'./Unity Analysis/Spring 2026/Figures/Mag_Maps/Mag_Maps_{q:.0e}/'):
    #     os.makedirs(f'./Unity Analysis/Spring 2026/Figures/Mag_Maps/Mag_Maps_{q:.0e}/')

    # for s_iter in ss:
    #     if not os.path.exists(f'./Unity Analysis/Spring 2026/Figures/Frac_Maps/Frac_Maps_{q:.0e}/Frac_Maps_{s_iter:.2e}/'):
    #         os.makedirs(f'./Unity Analysis/Spring 2026/Figures/Frac_Maps/Frac_Maps_{q:.0e}/Frac_Maps_{s_iter:.2e}/')

    #%% Plotting
    fig = plt.figure(figsize=(20, 8))
    fig.set_tight_layout(True)
    axes = fig.subplots(2, 9)

    print('Starting to loop through parameter space and plot fractional area maps...')

    # #%% Looping through parameter space
    # def update(frame):
    #     multiplier = multipliers[frame]
    #     print(f'Processing multiplier: {multiplier:.0e}')

    # def update(frame):
        # s = ss[frame]
        # separation = ss_str[frame]

    ss_str = ['0.2', '0.4', '0.6', '0.8', '1.0', '1/0.8', '1/0.6', '1/0.4', '1/0.2']
    ss = np.array([numexpr.evaluate(s) for s in ss_str])
    multipliers = [1e-1, 1e1]
    q = 1e-4

    for i, multiplier in enumerate(multipliers):
        for j, s in enumerate(ss):
            separation = ss_str[j]
            
            print('=' * 50)
            print(f'Processing s={s:.2e}, q={q:.0e}, multiplier={multiplier:.0e}...')

            fig.suptitle(f'Fractional Deviation Maps: $q = {q:.0e}$', y=1.0)

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

            #%% Extracting important parameters from simulations
            pixels = single_sim.pixels
            ang_width = single_sim.ang_width
            ang_res = single_sim.ang_res

            single_magnifications = single_sim.magnifications
            binary_magnifications = binary_sim.magnifications

            min_source_radius = binary_attributes.loc[(binary_attributes['s'] == f'{separation}') & (binary_attributes['q'] == q), 'min_source_radius'].values[0]
            max_source_radius = binary_attributes.loc[(binary_attributes['s'] == f'{separation}') & (binary_attributes['q'] == q), 'max_source_radius'].values[0]

            X_pix, Y_pix = np.meshgrid(np.linspace(-ang_width/2, ang_width/2, pixels), np.linspace(-ang_width/2, ang_width/2, pixels))

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

            print('  Calculating fractional difference map...')
            frac_diff_map = (binary_conv_mags - single_conv_mags) / single_conv_mags

            ax = axes[i, j]
            ax.clear()
            ax.set_title(f'$\\rho/\\xi={multiplier:.0e}$\n$s={separation}$')

            img = ax.imshow(frac_diff_map, extent=(-ang_width/2, ang_width/2, -ang_width/2, ang_width/2), cmap='Greys', vmin=np.min(frac_diff_map), vmax=np.max(frac_diff_map))

            img = ax.contour(
                X_pix, Y_pix, np.flip(frac_diff_map, axis=0),
                levels=[-1e-2, 1e-2],
                colors=['blue', 'red'],
                linewidths=0.5
            )

            ax.set_aspect('equal')
            
            # Ensuring x ticks are written in scientific notation
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

            # Deleting variables to free up memory
            del img
            del frac_diff_map
            del single_conv_mags
            del binary_conv_mags
            del source_profile
            del single_sim
            del binary_sim

            print(f'Contour plot created successfully for s={s:.2f}, q={q:.0e}, multiplier={multiplier:.0e}.')

    axes[1, 4].set_xlabel('X [$\\theta_E$]')
    axes[0, 0].set_ylabel('Y [$\\theta_E$]')
    axes[1, 0].set_ylabel('Y [$\\theta_E$]')
    
    # ani = animation.FuncAnimation(fig, update, frames=len(ss), repeat=False)
    # ani.save(f'./Unity Analysis/Spring 2026/Figures/mag_map_animation.mp4', writer='ffmpeg', dpi=300)

    fig.savefig(f'./Unity Analysis/Spring 2026/Figures/contours.png', dpi=300)
    # fig.savefig(f'./Unity Analysis/Spring 2026/Figures/Mag_Maps/Mag_Maps_{q:.0e}/mag_map_plot_{s:.2e}.png', dpi=300)



