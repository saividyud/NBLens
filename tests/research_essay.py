from IRSMicroLensing import IRSMicroLensSimulation as IRSMLS
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

pixels = 2500
ang_width = 7

r = np.sqrt(2.3)
theta = np.deg2rad(30)

lens_att1 = [
    [0, 0, 0.01, 1],
    [r*np.cos(theta), -r*np.sin(theta), 0.01, 0.001]
]

lens_att2 = [
    [0, 0, 0.01, 1]
]

lc1 = IRSMLS.IRSLightCurves(pixels=pixels, ang_width=ang_width, lens_att=lens_att1, source_r=0.05, u0=0.4, source_type='uniform', steps=300)
lc1.plot(show_plot=False)
print()
lc2 = IRSMLS.IRSLightCurves(pixels=pixels, ang_width=ang_width, lens_att=lens_att2, source_r=0.05, u0=0.4, source_type='uniform', steps=300)
lc2.plot(show_plot=False)

fig = plt.figure('Light Curves', figsize=(10, 8))
ax = fig.add_subplot()

ax.plot(lc1.time, lc1.image_fluxes, label='Star with Planet')
ax.plot(lc2.time, lc2.image_fluxes, label='Star without Planet', linestyle='--')

# \n$u_0$ = {lc1.source_att[0][1]} $\\theta_E$\n$\\rho$ = {lc1.source_att[0][2]} $\\theta_E$
ax.set_title(f'Light Curve for Gravitational Microlensing\nMax Flux = {round(np.nanmax(lc1.image_fluxes), 4)}')
ax.set_xlabel(f't [$t_E$]')
ax.set_ylabel('Normalized Flux')
# ax.set_xlim(min(lc1.time), max(lc1.time))
ax.set_xlim(-3, 3)
ax.set_ylim(0.75, 2.75)
ax.grid()

ax.legend()

fig.savefig('./figures/Comparison Light Curve.png', dpi=500)

# ax.plot(lc.time, lc.analytic_mags)

plt.show()