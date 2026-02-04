#%% Imports

import sys
sys.path.append('.')

import numpy as np
import pandas as pd
import numexpr

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

import IRSMicroLensing.IRSFunctions as IRSF
import IRSMicroLensing.IRSCaustics as IRSC

import VBMicrolensing
import MulensModel as mm
import os

import csv

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['figure.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.title_fontsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['mathtext.fontset'] = 'cm'

#%% Reading data files
df_list = []
thresholds_list = []
radius_list = []

one_e_4_list = []
three_e_4_list = []
one_e_3_list = []
three_e_3_list = []
one_e_2_list = []

multipliers = [1e-1, 3e-1, 1e0, 3e0, 1e1]

for multiplier in multipliers:
    df = pd.read_csv(f'./Unity Analysis/Spring 2026/Data Files/analysis_output_{multiplier:.0e}.csv')
    df_list.append(df)

#%% Defining s and q values
ss_str_1 = np.array(df_list[0]['s'].to_list(), dtype=str)
ss_1 = np.array([numexpr.evaluate(str(s)) for s in df_list[0]['s'].to_list()])
qs_1 = np.array(df_list[0]['q'].to_list())

ss_str_unique_1 = np.unique(ss_str_1)
ss_unique_1 = np.unique(ss_1)
qs_unique_1 = np.unique(qs_1)

print(ss_str_unique_1)
print(qs_unique_1)

# #%% Plotting q and s with fractional area as color scale for different thresholds
# fig = plt.figure(figsize=(14, 8))
# fig.set_tight_layout(True)
# axes = fig.subplots(5, 5)

# fig.suptitle('Detectibility as a Function of q and s', y=0.97)

# for i, threshold in enumerate([1.1e-04, 3.4e-04, 1.0e-03, 3.1e-03, 1.1e-02]):
#     for j, multiplier in enumerate(multipliers):
#         ax = axes[i, j]

#         ax.set_title(f'$x={threshold:.0e}$, $\\rho/\\xi={multiplier}$')

#         df_subset = df_list[j]

#         # Fixing holes in data
#         empty_arr = ['1/0.2', 3e-4] + [np.nan] * 102
#         df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[116.5])], ignore_index=False)

#         empty_arr = ['1/0.2', 1e-3] + [np.nan] * 102
#         df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[117.5])], ignore_index=False)
#         df_subset = df_subset.sort_index().reset_index(drop=True)

#         empty_arr = ['1/0.4', 1e-3] + [np.nan] * 102
#         df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[103.5])], ignore_index=False)

#         df_subset = df_subset.sort_index().reset_index(drop=True)

#         # ax.imshow(np.array(df_subset[f'{threshold:.1e}']).reshape((17, 7)).T, extent=[
#         #     ss_unique_1.min(), ss_unique_1.max(),
#         #     qs_unique_1.min(), qs_unique_1.max()
#         # ], aspect='auto', origin='lower', norm=plt.Normalize(vmin=0, vmax=1), cmap='viridis')

#         ax.contourf(
#             np.linspace(0, 1, 17),
#             qs_unique_1,
#             np.array(df_subset[f'{threshold:.1e}']).reshape((17, 7)).T,
#             levels=np.linspace(0, 1, 11),
#             cmap='viridis',
#             vmin=0, vmax=1
#         )

#         ax.set_yscale('log')

#         # ax.grid(True, which='both', linestyle='--', linewidth=0.5)
#         ax.set_xticks(np.linspace(0, 1, 17))
#         ax.set_xticklabels(ss_str_unique_1)

# axes[-1, 2].set_xlabel('Separation [$\\theta_E$]')
# axes[2, 0].set_ylabel('Mass Ratio')


# plt.show()
# %% Animating q and s with fractional area as color scale for different thresholds
fig = plt.figure(figsize=(14, 4))
# fig.set_tight_layout(True)
axes = fig.subplots(2, 5)

fig.suptitle('Detectibility as a Function of q and s', y=0.97)

# Adding colorbar
bar = plt.colorbar(
    plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1)),
    ax=axes[1, 2], orientation='horizontal',
    fraction=0.05, pad=0.1,
    label='Fractional Area Above Threshold'
)

def update(frame):
    fig.suptitle('Detectibility as a Function of q and s\nThreshold: {:.1e}'.format(frame), y=0.97)

    for ax in axes[0, :]:
        ax.cla()

    threshold = frame

    for j, multiplier in enumerate(multipliers):
        ax = axes[0, j]

        ax.set_title(f'$\\rho/\\xi={multiplier}$')

        df_subset = df_list[j]

        # Fixing holes in data
        empty_arr = ['1/0.2', 3e-4] + [np.nan] * 102
        df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[116.5])], ignore_index=False)

        empty_arr = ['1/0.2', 1e-3] + [np.nan] * 102
        df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[117.5])], ignore_index=False)
        df_subset = df_subset.sort_index().reset_index(drop=True)

        empty_arr = ['1/0.4', 1e-3] + [np.nan] * 102
        df_subset = pd.concat([df_subset, pd.DataFrame([empty_arr], columns=df_subset.columns, index=[103.5])], ignore_index=False)

        df_subset = df_subset.sort_index().reset_index(drop=True)

        ax.contourf(
            np.linspace(0, 1, 17),
            qs_unique_1,
            np.array(df_subset[f'{threshold:.1e}']).reshape((17, 7)).T,
            levels=np.linspace(0, 1, 11),
            cmap='viridis',
            vmin=0, vmax=1
        )

        ax.set_yscale('log')

        # ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.set_xticks(np.linspace(0, 1, 17))
        ax.set_xticklabels(ss_str_unique_1)

    axes[2].set_xlabel('Separation [$\\theta_E$]')
    axes[0].set_ylabel('Mass Ratio')

ani = animation.FuncAnimation(
    fig, update,
    frames=np.array(df_list[0].columns[4:], dtype=float).tolist(),
    interval=50
)

# ani.save('./Unity Analysis/Spring 2026/Figures/detectibility_animation.mp4', writer='ffmpeg', dpi=300)
plt.show()