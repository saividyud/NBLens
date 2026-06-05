from .imports import *
from .IRSCaustics import IRSCaustics
from .IRSFunctions import IRSFunctions

class IRSOrbiting(OrbitPropagator):
    '''
    Generates magnification maps for each point in an orbit.

    Parameters
    ----------
    bodies : Lx8 NDArray
        An Lx8 dimensional numpy array with L bodies and 8 terms for each body: either classical orbital elements or state vector + mass + central body indicator
            bodes = np.array([[rx1, ry1, rz1, vx1, vy1, vz1, m1, -1],            -1: No central body / is the central body (will read in Cartesian state vector)
                            [a, e, i, ta, aop, raan, m2, x],                      x: Central body is index x
                            ...                                                   NOTE: 7th index term is inputted ONLY if "coes" is TRUE
                            [rxN, ryN, rzN, vxN, vyN, vzN, mN, -1]])

    dt : float
        A small change in time to numerically integrate with (can greatly vary the speed of animation)

    t_span : float
        The full timespan of the animation (not in real time) or 0 if only one period is desired

    sim_dist : 1x2 list or ndarray
        Relevant distances in kiloparsecs in format:
            [Ds, Dl] -> Units: [kpc, kpc]

    dv : list (optional)
        Adds a change in velocity that the user specifies for a certain body at a certain time. 'r' is for rectangular, while 's' is for spherical. Format detailed below:
            dv = [magnitude, yaw angle, pitch angle, body index, time, time units, 's']
            dv = [dvx, dvy, dvz, body index, time, time units, 'r']

    coes : bool (optional)
        Boolean value for if the input contains Keplerian orbital elements

    data : bool (optional)
        Boolean for if a txt file of the states is needed

    Returns
    -------
    magnifications : SxNxNteps NDArray
        List of magnification maps for each position of lenses
    '''
    def __init__(self, bodies, dt, t_span, sim_dist, dv=0, coes=True, data=False):
        # Saving passed parameters
        self.sim_dist = sim_dist
        
        # Calculating states in 3D for each time step
        super().__init__(bodies=bodies, dt=dt, t_span=t_span, dv=dv, coes=coes, data=data, spheres=0, 
                                          limits=0, deg=True, anim=False, save_plot=False, dark_plot=False, title='Title', label=0)
        
        # Extracting positions and velocities in 3D
        self.positions_3D = self.states[:, :, :3]
        self.velocities_3D = self.states[:, :, 3:]

    
    '''
    Calculates the magnification maps for each time step.

    Parameters
    ----------
    pixels : int
        Number of pixels along a side of the detector

    ang_width : float
        Angular width of lens region in terms of Einstein ring radius (theta_e)

    theta : float
        asdf

    phi : float
        asdf

    '''
    def calculate(self, pixels, ang_width, rays_per_pixel, theta, phi, annulus=0, zoom: tuple | list = None):
        # Saving passed parameters
        self.pixels = pixels
        self.ang_width = ang_width
        self.rays_per_pixel = rays_per_pixel
        self.annulus = annulus
        self.zoom = zoom
        self.theta = np.deg2rad(theta)
        self.phi = np.deg2rad(phi)

        # Derived constants
        self.ang_res = self.ang_width / self.pixels # [theta_e/pixel]

        # Computing direction cosine matrix (DCM) for projection onto sky plane
        DCM = np.array([[np.cos(theta), np.sin(theta) * np.sin(phi), np.sin(theta) * np.cos(phi)],
                        [0, np.cos(phi), -np.sin(phi)],
                        [-np.sin(theta), np.cos(theta) * np.sin(phi), np.cos(theta) * np.cos(phi)]])
        
        # Calculating projected positions in 2D and setting third parameter to masses in Msun
        self.positions_2D = self.positions_3D

        for i, slice in enumerate(self.positions_3D):
            for j, row in enumerate(slice):
                self.positions_2D[i, j, :] = np.dot(DCM, row)
        
        self.positions_2D[:, :, 2] = self.mass/cb.sun['mass']

        # Translating to keep center of mass at origin
        for i, slice in enumerate(self.positions_2D):
            CM_sum = 0

            # Iterating through each body and multiplying the positions by the mass
            for j, body in enumerate(slice):
                CM_sum += body[:2] * self.mass[j]
            
            # Calculating distance to center of mass
            r_CM = CM_sum / np.sum(self.mass)

            # Translating center of mass to be at (0, 0) at all times
            slice[:, :2] = slice[:, :2] - r_CM

        # Calculating theta_E and r_E
        self.theta_E, self.r_E = IRSFunctions.e_ring(sim_dist=self.sim_dist, M=self.mass/cb.sun['mass'])

        # Rescaling all positions to be in terms of theta_E
        self.positions_2D[:, :, :2] = self.positions_2D[:, :, :2] / self.r_E.to(u.km).value

        # Initializing magnifications tensor
        self.magnifications = np.zeros(shape=(self.steps, self.pixels, self.pixels))
        self.magnifications_log = np.zeros(shape=(self.steps, self.pixels, self.pixels))

        # Calculating magnification map for each slice
        for step, slice in tqdm(enumerate(self.positions_2D), total=self.steps):
            # IRSFunctions.progress_bar(step, self.steps)

            # Inserting radius for each lens for this time step
            slice = np.insert(slice, 2, slice[:, 2]/100, 1)

            # Forming a list of lens attributes to pass into IRSCaustics
            lens_att = slice.tolist()

            # Calculating magnification map for this time step
            mag_map = IRSCaustics(pixels=self.pixels, lens_att=lens_att, ang_width=self.ang_width, rays_per_pixel=self.rays_per_pixel, annulus=self.annulus)
            self.magnifications[step] = mag_map.plot(show_plot=False, print_stats=False)
            self.magnifications_log[step] = mag_map.magnifications_log
            
            # Finding translated x and y coordinates after zoom
            x_lower_bound = int(self.pixels/2) - m.ceil(self.zoom[0]/(2*self.ang_res))
            x_upper_bound = int(self.pixels/2) + m.ceil(self.zoom[0]/(2*self.ang_res))

            y_lower_bound = int(self.pixels/2) - m.ceil(self.zoom[1]/(2*self.ang_res))
            y_upper_bound = int(self.pixels/2) + m.ceil(self.zoom[1]/(2*self.ang_res))

            magnifications_log_zoomed = self.magnifications_log[step, y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]

            if step == 0:
                self.magnifications_log_zoomed = np.zeros(shape=(self.steps, np.shape(magnifications_log_zoomed)[0], np.shape(magnifications_log_zoomed)[1]))

            self.magnifications_log_zoomed[step] = magnifications_log_zoomed


    def plot(self, save_plot=False, show_mm=False, show_lenses=False, show_dev=False, show_axes=True, print_stats=True, cmap='gray'):
        '''
        Plots an animation of teh magnification maps for an orbiting lens system.

        Parameters
        ----------
        save_plot : bool, optional
            If plots should be saved
        show_mm : bool, optional
            If MulensModel caustics should be plotted (only supported for two lenses)
        show_lenses : bool, optional
            If the lenses and Einstein ring should be plotted on top of caustics
        show_dev : bool, optional
            Calculates the deviation from single lens if a binary lens is passed
        show_axes : bool, optional
            If the plot should have formatting or take up the entire figure
        print_stats : bool, optional
            If the program outputs computation time and steps
        show_plot : bool, optional
            If the program creates the plot
        file_save : bool, optional
            Saves magnification map data in CSV file in ./Mag Map Data/{filename}.csv
        cmap: : string, optional
            Colormap of images

        Returns
        -------
        magnifications : NxN NDArray
            Matrix of magnification values for each pixel
        '''
        # Saving passed parameters
        self.save_plot = save_plot
        self.show_mm = show_mm
        self.show_lenses = show_lenses
        self.show_dev = show_dev
        self.show_axes = show_axes
        self.print_stats = print_stats
        self.cmap = cmap

        # Initializing figure
        self.fig = plt.figure(figsize=(10, 8))
        self.ax = self.fig.add_subplot()

        self.ax.set_xlabel('X [$\\theta_E$]')
        self.ax.set_ylabel('Y [$\\theta_E$]')

        # Plotting first magnification map
        self.mag_map = self.ax.imshow(self.magnifications_log_zoomed[0], cmap=self.cmap, extent=[-self.zoom[0]/2, self.zoom[0]/2, -self.zoom[1]/2, self.zoom[1]/2], animated=True)

        if self.show_lenses:
            self.lenses = []
            self.lens_patches = []
            for lens in range(self.count):
                self.lenses.append(patches.Circle((self.positions_2D[0, lens, 0], self.positions_2D[0, lens, 1]), self.positions_2D[0, lens, 2]/100))
                self.lens_patches.append(self.ax.add_patch(self.lenses[lens]))
            
        # self.ax_scatter_plots = []
        # if self.show_lenses:
        #     for lens in range(self.count):
        #         self.ax_scatter_plots.append(self.ax.scatter(self.positions_2D[0, lens, 0], self.positions_2D[0, lens, 1]))

        # Animating magnification map
        self.anim = animation.FuncAnimation(fig=self.fig, func=self.animate_func, frames=self.steps)

        if self.save_plot:
            f = '../movies/anim.mp4'
            writergif = animation.FFMpegWriter(fps=self.steps/6)

            # writergif.setup(fig=fig, outfile=f, dpi=1200)

            print('Saving animation')
            self.anim.save(f, writer=writergif)

    def animate_func(self, frame):
        # Code for calculating timing and scaling appropriately
        time_f = np.max(self.ets)
        time = self.ets[frame]

        if time_f >= (86400 * 365.25):
            time_f = time_f / (86400 * 365.25)
            time = self.ets[frame] / (86400 * 365.25)
            units = 'years'
        elif time_f >= 86400:
            time_f = time_f / 86400
            time = self.ets[frame] / 86400
            units = 'days'
        elif time_f >= 3600:
            time_f = time_f / 3600
            time = self.ets[frame] / 3600
            units = 'hrs'
        else:
            units = 's'

        # Setting the title
        self.fig.suptitle('Time = ' + str(np.round(time, 2)) + ' ' + units)
        
        # Removing magnification map
        self.mag_map.remove()

        if self.show_lenses:
            for lens in self.lens_patches:
                lens.remove()

        self.ax.set_title('Magnification Map')

        # Drawing next frame on magnification map
        self.mag_map = self.ax.imshow(self.magnifications_log_zoomed[frame], cmap='gray', extent=[-self.zoom[0]/2, self.zoom[0]/2, -self.zoom[1]/2, self.zoom[1]/2], animated=True)

        if self.show_lenses:
            self.lenses = []
            self.lens_patches = []
            for lens in range(self.count):
                self.lenses.append(patches.Circle((self.positions_2D[frame, lens, 0], self.positions_2D[frame, lens, 1]), 0.01, color='red'))
                self.lens_patches.append(self.ax.add_patch(self.lenses[lens]))

    @property
    def zoom(self):
        '''
        Type : tuple or list

        Zoom of axes when viewing caustics.
        '''
        return self._zoom
    
    @zoom.setter
    def zoom(self, val):
        # Checking if zoom is a tuple or list
        if isinstance(val, (tuple, list)):
            # Checking if length of zoom is 2
            if len(val) == 2:
                for element in val:
                    # Checking if each element in zoom is an integer or float
                    if isinstance(element, (int, float)):
                        # Checking if each element in zoom is greater than 0
                        if element > 0:
                            self._zoom = val
                        else:
                            raise ValueError(f'Each element in attribute "zoom" must be positive. Got {element}.')
                    else:
                        raise TypeError(f'Each element in attribute "zoom" must be an integer or float. Got {type(element)}.')
            else:
                raise ValueError(f'Attribute "zoom" must be of length 2. got {len(val)}.')
        elif isinstance(val, type(None)):
            self._zoom = [self.ang_width, self.ang_width]
        else:
            raise TypeError(f'Attribute "zoom" must be an tuple or list. Got {type(val)}.')
