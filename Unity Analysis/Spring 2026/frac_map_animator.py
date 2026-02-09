if __name__ == '__main__':

    #%% Importing modules
    import sys
    sys.path.append('.')

    import numpy as np
    import numexpr
    import pandas as pd

    import IRSMicroLensing.IRSCaustics as IRSC

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    import argparse
    import sys

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['figure.titlesize'] = 20
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['figure.figsize'] = (14, 8)
    plt.rcParams['axes.titlesize'] = 20
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
    parser = argparse.ArgumentParser(description='Plot fractional area maps for microlensing simulations.')
    parser.add_argument('--multiplier', type=float, default=1.0, help='Multiplier for output filename.')
    args = parser.parse_args()
    multiplier = args.multiplier
    
    print(f'Arguments parsed successfully. Starting analysis with multiplier = {multiplier}...')
    
    #%% Defining parameter space
    ss_str = ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1/0.9', '1/0.8', '1/0.7', '1/0.6', '1/0.5', '1/0.4', '1/0.3', '1/0.2']
    ss = np.array([numexpr.evaluate(s) for s in ss_str])
    qs = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
    multipliers = [1e-1, 3e-1, 1e0, 3e0, 1e1]
    binary_attributes = pd.read_csv('./Unity Analysis/Spring 2026/binary_attributes.csv')

    print('Parameter space defined successfully. Starting to read fractional area maps...')

    #%% Plotting
    fig = plt.figure(figsize=(22, 26))
    fig.set_tight_layout(True)
    axes = fig.subplots(9, 7)

    fig.suptitle('Fractional Area v Source Star Radius ($s \leq 1$)', y=1.0)

    print('Starting to loop through parameter space and plot fractional area maps...')

    # #%% Looping through parameter space
    # def update(frame):
    #     multiplier = multipliers[frame]
    #     print(f'Processing multiplier: {multiplier:.0e}')

    for i, s in enumerate(ss[ss <= 1]):
        for j, q in enumerate(qs):
            print(f'Processing s={s:.2f}, q={q:.0e} with multiplier={multiplier:.0e}...')

            ax = axes[i, j]
            ax.set_title(f'$s={ss_str[i]}$, $q={q:.0e}$')

            frac_map = np.load(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/Frac_Maps_{q:.0e}/frac_diff_map_q{q:.0e}_s{s:.2e}.npy')
            print(f'Fractional area map loaded successfully for s={s:.2f}, q={q:.0e}, multiplier={multiplier:.0e}.')

            ang_width = binary_attributes.loc[(binary_attributes['q'] == q) & (binary_attributes['s'] == ss_str[i]), 'ang_width'].values[0]
            pixels = binary_attributes.loc[(binary_attributes['q'] == q) & (binary_attributes['s'] == ss_str[i]), 'pixels'].values[0]
            print(f'Attributes retrieved successfully for s={s:.2f}, q={q:.0e}: ang_width={ang_width}, pixels={pixels}.')

            # X_pix, Y_pix = np.meshgrid(np.linspace(-ang_width/2, ang_width/2, pixels), np.linspace(-ang_width/2, ang_width/2, pixels))
            X_pix = np.linspace(-ang_width/2, ang_width/2, pixels)
            Y_pix = np.linspace(-ang_width/2, ang_width/2, pixels)
            print(f'Pixel grid created successfully for s={s:.2f}, q={q:.0e}.')

            img = ax.contour(
                X_pix, Y_pix, np.flip(frac_map, axis=0),
                levels=log_array(-3, -1).tolist(),
                colors=colors_by_log(log_array(-3, -1)),
                linewidths=0.5
            )

            ax.set_aspect('equal')

            # Deleting variables to free up memory
            del img
            del frac_map

            print(f'Contour plot created successfully for s={s:.2f}, q={q:.0e}.')
    # ani = animation.FuncAnimation(fig, update, frames=len(multipliers), repeat=False)
    fig.savefig(f'./Unity Analysis/Spring 2026/Figures/frac_map_plot_{multiplier:.0e}_s<1.png', dpi=300)



