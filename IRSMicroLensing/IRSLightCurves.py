from .imports import *
from .IRSMain import IRSMain
from .IRSCaustics import IRSCaustics

class IRSLightCurves(IRSMain):
    '''
    Plots the light curve for a microlensing event from the entrance into the detector to exit.

    Parameters
    ----------
    pixels : int
        Number of pixels along a side of the detector
    ang_width : float
        Angular width of lens region in terms of Einstein ring radius (theta_e)
    lens_att : Lx4 array, with L number of lenses
        Lens attributes in format:
            [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]
    source_type : string
        Type of source brightness distribution
        Supported types:
            uniform: homogenous circular distribution
            Gauss: Gaussian circular distribution
    u0 : float
        Impact parameter in units of theta_e
     steps : int
        Number of steps for light curve
    source_r : float
        Radius of source in units of theta_e
    save_plot : bool
        If the 
    animate : bool
        If the transit should be animated
    save_anim : bool
        If the animation should be saved
    show_mm : bool
        If the MulensModel light curve should be plotted

    Returns
    -------
    IRSLightCurves object
    '''
    def __init__(self, pixels: int, ang_width: int | float, lens_att: list, source_r: int | float, u0: int | float, steps: int, source_type='uniform', *args, **kwargs):
        self.steps = steps

        if isinstance(ang_width, (int, float)):
            pass
        else:
            raise TypeError(f'Attribute "ang_width" must be an integer or float. Got {type(ang_width)}.')

        # Initializing list of source positions (varying x positions and keeping y position as the impact parameter + center of mass)
        source_att = np.zeros((self.steps+1, 3))
        source_att[:, 0] = np.linspace(-ang_width/2 - source_r, ang_width/2 + source_r, steps+1)
        source_att[:, 1] = u0
        source_att[:, 2] = source_r
        
        source_att = source_att.tolist()

        for i in range(self.steps+1):
            source_att[i].append(source_type)

        super().__init__(pixels=pixels, ang_width=ang_width, lens_att=lens_att, source_att=source_att)
    
    def plot(self, print_stats=True, animate=False, show_mm=False, show_plot=True, save=False, show_caustics=False, cm_offset: tuple | list = [0, 0]):
        self.save = save
        self.print_stats = print_stats
        self.show_mm = show_mm
        self.show_plot = show_plot
        self.animate = animate
        self.show_caustics = show_caustics
        self.cm_offset = cm_offset

        # Calcualte lens center of mass (1x2 NDArray)
        self.lens_CM = self.calc_CM()

        # Translating lens positions so the center of mass is at (0, 0)
        self.lens_att[:, :2] = self.lens_att[:, :2] - self.lens_CM + self.cm_offset

        # Time array in units of t_e, from entrance into the detector to exit
        self.time = np.linspace(-self.ang_width/2, self.ang_width/2, self.steps+1) - self.cm_offset[0]

        # Initializing arrays of image flux values and list of IRSMicroLensing objects
        init_time = t.time()
        self.image_fluxes, self.image_brightnesses = super().calculate(light_curves=True, print_stats=self.print_stats)
        final_time = t.time() - init_time
        if self.print_stats:
            print('---------------------')
            print(f'Computational time: {round(final_time, 3)} seconds')
            print(f'Calculation time per step: {round(final_time/self.steps, 3)} seconds')

        # Calculating analytic light curve
        if self.show_mm:
            self.analytic_mags = self.calc_mm_mags()

        if self.animate:
            self.fig = self.plot_animation()

            # Animation object
            self.anim = animation.FuncAnimation(fig=self.fig, func=self.animate_func, interval=10, frames=(self.steps+1), cache_frame_data=False)

            if self.save:
                # Saving the Animation
                f = f"{self.import_file}.gif"
                print('Created PillowWriter')
                writergif = animation.PillowWriter(fps=self.steps/6)
                print('Saving animation')
                self.anim.save(f, writer=writergif)

        elif self.show_plot:
            # Plotting light curve
            self.fig = self.plot_lightcurve()
            
            if self.save:
                # Saving light curve
                self.fig.savefig(f'./figures/{self.import_file}.png', dpi=500)

        return self.image_fluxes

    def calc_mm_mags(self):
        '''
        Calculating magnifications for MulensModel.

        Parameters
        ----------
        None

        Returns
        -------
        analytic_mags : 1x(num+1) list
            Magnifications returned by running MulensModel on given parameters
        '''
        if self.L == 1:
            # Defining single lens model using given t_0, t_E, u_0, and rho
            self.param_dict = {'t_0': 0, 't_E': 1, 'u_0': self.source_att[0][1], 'rho': self.source_att[0][2]} 
            model = mm.Model(self.param_dict)

            # Setting magnification method (finite source)
            model.set_default_magnification_method('finite_source_uniform_WittMao94')

            # Retrieving magnifications with given time array
            analytic_mags = model.get_magnification(self.time)

        elif self.L == 2:
            # Finding which index the bigger mass was passed in (primary lens)
            big_mass = np.where(self.lens_att[:, 3] == np.max(self.lens_att[:, 3]))[0][0]

            # Secondary lens
            small_mass = int(not big_mass)

            # Defining unit source vector
            uhat = [1, 0]

            # Defining unit binary axis vector (from primary to secondary lens)
            v = [self.lens_att[small_mass, 0] - self.lens_att[big_mass, 0], self.lens_att[small_mass, 1] - self.lens_att[big_mass, 1]]
            vhat = v / np.linalg.norm(v)

            # Finding counterclockwise angle between binary axis and source trajectory (alpha)
            # First quadrant
            if vhat[0] > 0 and vhat[1] > 0:
                alpha = np.arctan(vhat[1]/vhat[0])
            
            # Second quadrant
            elif vhat[0] < 0 and vhat[1] > 0:
                alpha = np.pi + np.arctan(vhat[1]/vhat[0])

            # Third quadrant
            elif vhat[0] < 0 and vhat[1] < 0:
                alpha = np.pi + np.arctan(vhat[1]/vhat[0])

            # Fourth quadrant
            elif vhat[0] > 0 and vhat[1] < 0:
                alpha = np.arctan(vhat[1]/vhat[0])

            # Positive x direction
            elif vhat[0] > 0 and vhat[1] == 0:
                alpha = 0

            # Positive y direction
            elif vhat[0] == 0 and vhat[1] > 0:
                alpha = np.pi/2

            # Negative x direction
            elif vhat[0] < 0 and vhat[1] == 0:
                alpha = np.pi

            # Negative y direction
            elif vhat[0] == 0 and vhat[1] < 0:
                alpha = -np.pi/2

            # Calculating distance between lenses (s)
            s = np.linalg.norm(v)

            # Calculating mass ratio between lenses (q)
            q = self.lens_att[small_mass, 3] / self.lens_att[big_mass, 3]

            # Defining binary lens model
            self.param_dict = {'t_0': 0, 't_E': 1, 'u_0': self.source_att[0][1], 'rho': self.source_att[0][2], 'q': q, 'alpha': np.rad2deg(alpha), 's': s}
            model = mm.Model(self.param_dict)

            # Setting magnification method (finite source with limb darkening)
            model.set_default_magnification_method('VBBL')

            # Retrieving magnifications with given time array
            analytic_mags = model.get_magnification(self.time)
        
        else:
            self.show_mm = False
            analytic_mags = 0
            warnings.warn(f'MulensModel only valid for 1 or 2 lenses. Got {self.L} lenses.', SyntaxWarning)

        '''
        Finite source methods offered by MulensModel (single lens):
            finite_source_uniform_Gould94
            finite_source_uniform_Gould94_direct
            finite_source_uniform_WittMao94
            finite_source_LD_WittMao94
            finite_source_LD_Yoo04
            finite_source_LD_Yoo04_direct
            finite_source_uniform_Lee09
            finite_source_LD_Lee09
        '''
        '''
        Point source methods offered by MulensModel (binary lens):
            point_source
            quadrapole
            hexadecapole
            point_source_point_lens

        Finite source methods offered by MulensModel (binary lens):
            VBBL (limb darkening)
            Adaptive_Contouring (limb darkening)
        '''

        return analytic_mags

    def plot_lightcurve(self):
        '''
        Plots the light curve with given trajectories.

        Parameters
        ----------
        show_mm : bool
            If the MulensModel magnifications should be plotted or not
        
        Returns
        -------
        matplotlib.figure.Figure
        '''
        # Readying figure
        fig = plt.figure('Light Curves', figsize=(10, 8))
        ax = fig.add_subplot()

        # Plotting analytic light curve
        ax.plot(self.time, self.image_fluxes, label='Simulated Light Curve')

        # Plotting MulensModel analytic light curve if enabled
        if self.show_mm:
            model_param_str = ''
            for key, value in self.param_dict.items():
                if key != list(self.param_dict)[-1]:
                    if key == 'rho':
                        model_param_str += f'$\\{key}$' + f' = {round(value, 4)}\n'
                    elif key == 'alpha':
                        model_param_str += f'$\\{key}$' + f' = {round(value, 4)}°\n'
                    else:
                        model_param_str += f'${key}$' + f' = {round(value, 4)}\n'
                else:
                    model_param_str += f'${key}$' + f' = {round(value, 4)}'

            ax.plot(self.time, self.analytic_mags, label='MulensModel Light Curve')
            props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
            ax.text(0.02, 0.98, model_param_str, va='top', zorder=10, bbox=props, transform=ax.transAxes)
            ax.legend()

        # Formatting plot
        ax.set_title(f'Light Curve for Gravitational Microlensing\n$u_0$ = {self.source_att[0][1]} $\\theta_E$\n$\\rho$ = {self.source_att[0][2]} $\\theta_E$\nMax Flux = {round(np.nanmax(self.image_fluxes), 4)}')
        ax.set_xlabel(f't [$t_E$]')
        ax.set_ylabel('Normalized Flux')
        ax.set_xlim(min(self.time), max(self.time))
        ax.grid()

        return fig

    def plot_animation(self):
        '''
        Initializes the animation figure and plots and formats them.

        Parameters
        ----------
        None

        Returns
        -------
        matplotlib.figure.Figure
        '''
        # Initializing figure
        fig = plt.figure('Microlensing Animation', figsize=(12, 8))
        fig.tight_layout()

        # Creating image plane axis
        self.image_ax = fig.add_subplot(1, 2, 1) # Contains image contour and source patch
        
        # Creating light curve axis
        self.lc_ax = fig.add_subplot(1, 2, 2) # Contains light curve line and MulensModel line (if enabled)

        # Creating static axis (does not change while animating, on top of image axis)
        self.static_ax = fig.add_subplot(1, 2, 1) # Contains lens patches, Einstein ring patch, source trajectory line, and center of mass marker

        # Initializing lens patch list
        self.lenses = []
        for i in range(self.L):
            # Plotting lensing star
            if self.L < 11:
                self.lenses.append(patches.Circle((self.lens_att[i, 0], self.lens_att[i, 1]), self.lens_att[i, 2], color=list(colors.TABLEAU_COLORS.values())[i], label=f'r = {round(self.lens_att[i, 2], 3)}, $\\epsilon$ = {round(self.lens_att[i, 3]/self.total_M, 3)}' + '$M_{total}$'))
            else:
                self.lenses.append(patches.Circle((self.lens_att[i, 0], self.lens_att[i, 1]), self.lens_att[i, 2], color='cyan'))

        # Initializing source patch list
        self.sources = []
        for i in range(self.steps+1):
            self.sources.append(patches.Circle((self.source_att[i][0], self.source_att[i][1]), self.source_att[i][2], color='yellow', fill=False))

        # Initializing Einstein ring patch
        self.e_ring_patch = patches.Circle(self.cm_offset, 1, color='cyan', fill=False, linestyle='dashed', linewidth=0.5)

        # Displaying model parameters (if MulensModel is enabled)
        if self.show_mm:
            model_param_str = ''
            for key, value in self.param_dict.items():
                if key == 'alpha' or key == 'rho':
                    model_param_str += f'$\\{key}$' + f' = {round(value, 4)}\n'
                else:
                    model_param_str += f'${key}$' + f' = {round(value, 4)}\n'
            fig.text(0.125, 0.16, model_param_str, va='top')

        '''Formatting axes'''
        # Image axis on top of static axis
        self.image_ax.set_zorder(5)
        self.static_ax.set_zorder(10)

        # Setting static axes limits and aspect ratio
        self.static_ax.set_xlim(-self.ang_width/2, self.ang_width/2)
        self.static_ax.set_ylim(-self.ang_width/2, self.ang_width/2)
        self.static_ax.set_aspect('equal')

        # Turning off axes
        self.static_ax.axis('off')

        # Setting title of image
        self.image_ax.set_title(f'Image Brightness')

        # Setting image labels
        self.image_ax.set_xlabel('X [$\\theta_E$]')
        self.image_ax.set_ylabel('Y [$\\theta_E$]')

        # Setting image limits and axes aspect ratio
        self.image_ax.set_xlim(-self.ang_width/2, self.ang_width/2)
        self.image_ax.set_ylim(-self.ang_width/2, self.ang_width/2)
        self.image_ax.set_aspect('equal')

        # Setting title to light curve axis
        self.lc_ax.set_title(f'Light Curve')
        
        # Setting light curve labels
        self.lc_ax.set_xlabel(f't [$t_E$]')
        self.lc_ax.set_ylabel('Normalized Flux')
        
        # Seting light curve limits and adding grid
        self.lc_ax.set_xlim(min(self.time), max(self.time))
        self.lc_ax.set_ylim(min(self.image_fluxes), max(self.image_fluxes))
        self.lc_ax.grid()

        # Calculating caustics
        if self.show_caustics:
            self.caustic_ax = fig.add_subplot(1, 2, 1)
            
            self.caustic_ax.set_zorder(6)

            # Setting caustic axes limits and aspect ratio
            self.caustic_ax.set_xlim(-max(self.time), max(self.time))
            self.caustic_ax.set_ylim(-max(self.time), max(self.time))
            self.caustic_ax.set_aspect('equal')

            # Turning off axes
            self.caustic_ax.axis('off')

            mag_plot = IRSCaustics(pixels=self.pixels, ang_width=self.ang_width, lens_att=self.lens_att_list)
            self.magnifications = mag_plot.plot(print_stats=False, show_plot=False)

            # self.magnifications = np.flip(self.magnifications, axis=0)

            # Replacing all 0 values with 0.1 to plot in log10 space
            magnifications_log = np.where(self.magnifications == 0, 0.1, self.magnifications)
            magnifications_log = np.log10(magnifications_log)
            magnifications_log = np.where(magnifications_log == -1, 0, magnifications_log)

            caustics_countour = self.caustic_ax.imshow(magnifications_log, cmap='plasma', vmin=0, alpha=0.5, extent=[-self.ang_width/2, self.ang_width/2, -self.ang_width/2, self.ang_width/2])

        '''Plotting on static axes'''
        # Plotting lens center of mass
        lens_CM_point = self.static_ax.scatter(self.cm_offset[0], self.cm_offset[1], s=50, marker='+', c='cyan')

        # Plotting source trajectory
        source_traj = self.static_ax.axhline(y=self.source_att[0][1], xmin=-self.ang_width/2, xmax=self.ang_width/2, color='yellow', linestyle='dashed', linewidth=0.5)

        # Plotting lens patches
        for lens in self.lenses:
            self.static_ax.add_patch(lens)

        # Plotting Einstein ring patch
        self.static_ax.add_patch(self.e_ring_patch)

        # Printing lens attributes
        if self.L <= 2:
            if self.show_mm: location = 'upper right'; pos = (1, -0.1)
            else: location = 'upper center'; pos = (0.5, -0.1)
            self.static_ax.legend(loc=location, bbox_to_anchor=pos, ncol=1, frameon=False)
        elif self.L > 2 and self.L < 11:
            self.static_ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False)

        '''Plotting on image axes'''
        # Plotting contour of image
        # self.image_cont = self.image_ax.contourf(self.X, self.Y, self.image_brightnesses[0], levels=50, cmap='afmhot', vmin=0.0, vmax=1.0)
        self.image_brightnesses[0] = np.flip(self.image_brightnesses[0], axis=0)
        self.image_cont = self.image_ax.imshow(self.image_brightnesses[0], cmap='afmhot', vmin=0.0, vmax=1.0, extent=[-self.ang_width/2, self.ang_width/2, -self.ang_width/2, self.ang_width/2])

        # Plotting source star moving across director
        self.image_patches = self.image_ax.add_patch(self.sources[0])

        '''Plotting on light curve axes'''
        # Plotting simulated curve
        self.lc_line = self.lc_ax.plot(self.time[0], self.image_fluxes[0], zorder=10, label='Simulated Light Curve', color='blue')

        # Plotting analytic light curve if MulensModel was enabled
        if self.show_mm:
            self.mm_line = self.lc_ax.plot(self.time[0], self.analytic_mags[0], label='MulensModel Light Curve', color='red')
            self.lc_ax.set_ylim(min(self.image_fluxes), max(max(self.image_fluxes), max(self.analytic_mags)))

            # Adding legend to differentiate between MulensModel line and simulation line
            self.lc_ax.legend(loc='upper right', bbox_to_anchor=(1.14, -0.026), frameon=False)

        # Creating title for whole figure
        title = f'Flux = {round(self.image_fluxes[0], 4)}\nx = {round(self.source_att[0][0], 4)}\n{self.L} lenses'
        fig.suptitle(title)

        return fig

    def animate_func(self, frame):
        '''
        Iterable function to draw animated artists.

        Parameters
        ----------
        frame : int
            Identifying which frame to draw

        Returns
        -------
        None
        '''
        # Clearing axes

        # Removing all lines in light curve
        for line in self.lc_ax.get_lines():
            line.remove()

        # Removing image contour map
        # for coll in self.image_cont.collections:
        #     coll.remove()
        self.image_cont.remove()

        # Removing source patch
        self.image_patches.remove()

        # Plotting contour of image
        # self.image_cont = self.image_ax.contourf(self.X, self.Y, self.image_brightnesses[frame], levels=50, cmap='afmhot', vmin=0.0, vmax=1.0)
        mags = np.flip(self.image_brightnesses[frame], axis=0)
        self.image_cont = self.image_ax.imshow(mags, cmap='afmhot', vmin=0.0, vmax=1.0, extent=[-self.ang_width/2, self.ang_width/2, -self.ang_width/2, self.ang_width/2])

        # Plotting source star moving across director
        self.image_patches = self.image_ax.add_patch(self.sources[frame])

        # Plotting simulated curve
        self.lc_line = self.lc_ax.plot(self.time[:frame], self.image_fluxes[:frame], zorder=10, label='Simulated Light Curve', color='blue')

        # Plotting analytic light curve if MulensModel was enabled
        if self.show_mm:
            self.mm_line = self.lc_ax.plot(self.time[:frame], self.analytic_mags[:frame], label='MulensModel Light Curve', color='red')
        
        # Creating title for whole figure
        title = f'Flux = {round(self.image_fluxes[frame], 4)}\nx = {round(self.source_att[frame][0], 4)}\n{self.L} lenses'
        self.fig.suptitle(title)

    @property
    def steps(self):
        '''
        Type : int

        Number of steps for light curve.
        '''
        return self._steps
    
    @steps.setter
    def steps(self, val):
        if isinstance(val, int):
            # Checking if steps is greater than 0
            if val > 0:
                self._steps = val
            else:
                raise ValueError(f'Attribute "steps" must be positive. Got {val}.')
        else:
            raise TypeError(f'Attribute "steps" must be an integer. Got {type(val)}.')
