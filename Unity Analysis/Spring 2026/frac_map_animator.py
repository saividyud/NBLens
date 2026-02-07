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
    
    #%% Defining parameter space
    ss_str = ['0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', '1/0.9', '1/0.8', '1/0.7', '1/0.6', '1/0.5', '1/0.4', '1/0.3', '1/0.2']
    ss = np.array([numexpr.evaluate(s) for s in ss_str])
    qs = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
    multipliers = [1e-1, 3e-1, 1e0, 3e0, 1e1]
    binary_attributes = pd.read_csv('./Unity Analysis/Spring 2026/binary_attributes.csv')

    #%% Plotting
    fig = plt.figure(figsize=(22, 26))
    fig.set_tight_layout(True)
    axes = fig.subplots(9, 7)

    fig.suptitle('Fractional Area v Source Star Radius ($s \leq 1$)', y=1.0)

    #%% Looping through parameter space
    def update(frame):
        multiplier = multipliers[frame]

        for i, s in enumerate(ss[ss <= 1]):
            for j, q in enumerate(qs):
                ax = axes[i, j]
                ax.set_title(f'$s={s}$, $q={q:.0e}$')

                # binary_sim = IRSC.caustic_reader(f'./Unity/Simulations/Binary_Collection/Binary_Collection_{q:.0e}/binary_{q:.0e}_{s:.2e}.pkl')

                frac_map = np.load(f'./Unity/Simulations/Frac_Maps/Frac_Maps_{multiplier:.0e}/frac_diff_map_q{q:.0e}_s{s:.2e}.npy')

                ang_width = binary_attributes.loc[(binary_attributes['q'] == q) & (binary_attributes['s'] == s), 'ang_width'].values[0]
                pixels = binary_attributes.loc[(binary_attributes['q'] == q) & (binary_attributes['s'] == s), 'pixels'].values[0]

                X_pix, Y_pix = np.meshgrid(np.linspace(-ang_width/2, ang_width/2, pixels), np.linspace(-ang_width/2, ang_width/2, pixels))

                img = ax.contour(
                    X_pix, Y_pix, np.flip(frac_map, axis=0),
                    levels=log_array(-3, -1).tolist(),
                    colors=colors_by_log(log_array(-3, -1)),
                    linewidths=0.5
                )

    ani = animation.FuncAnimation(fig, update, frames=len(multipliers), repeat=False)
    ani.save('./Unity Analysis/Spring 2026/Figures/frac_map_animation.gif', writer='pillow', fps=1)



