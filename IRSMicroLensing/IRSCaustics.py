from .imports import *
from .IRSMain import IRSMain

class IRSCaustics(IRSMain):
    '''
    Plots the magnification map for L number of lenses.

    Parameters
    ----------
    pixels : int
        Number of pixels along a side of the detector
    ang_width : float
        Angular width of lens region in terms of Einstein ring radius (theta_e)
    lens_att : 1x4 array
        Lens attributes in format:
            [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]

    Returns
    -------
    IRSCaustics object
    '''
    def __init__(self, pixels: int, ang_width: int | float, lens_att: list, rays_per_pixel: int = 1, annulus: float = 0, *args, **kwargs):
        self.rays_per_pixel = rays_per_pixel
        self.annulus = annulus
        super().__init__(pixels=pixels, ang_width=ang_width, source_att=None, lens_att=lens_att, *args, **kwargs)

    def plot(self, zoom: tuple | list = None, cm_offset: tuple | list = [0, 0], save_plot=False, show_mm=False, show_lenses=False, show_dev=False, show_axes=True, print_stats=True, show_plot=True, file_save=False, cmap='gray'):
        '''
        Creates magnification map for lens system.

        Parameters
        ----------
        zoom : tuple or int
            How much the plot should be zoomed in by
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
        self.zoom = zoom
        self.cm_offset = cm_offset
        self.save_plot = save_plot
        self.show_mm = show_mm
        self.show_lenses = show_lenses
        self.show_dev = show_dev
        self.show_axes = show_axes
        self.print_stats = print_stats
        self.show_plot = show_plot
        self.file_save = file_save
        self.cmap = cmap

        # Compiling Numba jit functions
        init_rand_arr = np.random.randint(0, 10, size=(2, 2))
        IRSCaustics.calc_uniques(init_rand_arr)

        count_rand_arr = np.random.randint(0, 10, size=2)
        magnifications = np.zeros(shape=(100, 100), dtype=np.int64)
        IRSCaustics.calc_mags(100, magnifications, init_rand_arr, count_rand_arr)

        X_comp = np.random.random((3, 3))
        Y_comp = np.random.random((3, 3))
        ann_comp = 0.01
        IRSCaustics.create_annulus(ann_comp, X_comp, Y_comp)

        begin_time = t.time()

        # Calculatiing lens center of mass
        self.lens_CM = self.calc_CM()

        # Translating lens positions so the center of mass is at offsetted center of mass
        self.lens_att[:, :2] = self.lens_att[:, :2] - self.lens_CM + self.cm_offset

        # Creating mesh grid
        init_time = t.time()

        if self.annulus != 0:
            x_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)
            y_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)

            # Finding translated x and y coordinates after zoom
            x_lower_bound = int(self.pixels*self.rays_per_pixel/2) - m.ceil((2+self.annulus)/(2*self.ang_res))*self.rays_per_pixel
            x_upper_bound = int(self.pixels*self.rays_per_pixel/2) + m.ceil((2+self.annulus)/(2*self.ang_res))*self.rays_per_pixel

            y_lower_bound = int(self.pixels*self.rays_per_pixel/2) - m.ceil((2+self.annulus)/(2*self.ang_res))*self.rays_per_pixel
            y_upper_bound = int(self.pixels*self.rays_per_pixel/2) + m.ceil((2+self.annulus)/(2*self.ang_res))*self.rays_per_pixel

            x_rays = x_rays[x_lower_bound:x_upper_bound]
            y_rays = y_rays[y_lower_bound:y_upper_bound]

            # Creating meshgrid of ray coordinates
            X_R, Y_R = np.meshgrid(x_rays, y_rays)

            ray_valid_coords = IRSCaustics.create_annulus(self.annulus, X_R, Y_R)
            self.X = X_R[ray_valid_coords].reshape((-1, 1))
            self.Y = Y_R[ray_valid_coords].reshape((-1, 1))
        else:
            x_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)
            y_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)
            self.X, self.Y = np.meshgrid(x_rays, y_rays)

        final_time = t.time() - init_time
        if self.print_stats: print(f'Creating mesh grid: {round(final_time, 3)} seconds')

        # Calculating source pixels
        init_time = t.time()
        self.xs, self.ys = self.calc_source_pixels()
        final_time = t.time() - init_time
        if self.print_stats: print(f'Calculating source pixels: {round(final_time, 3)} seconds')

        # Calculating indices of translated pixel after deflection
        init_time = t.time()
        self.indx, self.indy = self.trans_ind()
        final_time = t.time() - init_time
        if self.print_stats: print(f'Calculating indices of translated pixel after deflection: {round(final_time, 3)} seconds')

        # Calculating translated pixels
        init_time = t.time()

        # Finding wherever indx or indy is nan
        bool_arr_x = np.isnan(self.indx)
        bool_arr_y = np.isnan(self.indy)

        # Replacing nan values with some out of bounds value and making the numbers into integers
        indx_nonan = np.where(bool_arr_x, self.pixels*2, self.indx)
        indy_nonan = np.where(bool_arr_y, self.pixels*2, self.indy)

        indx = np.where(bool_arr_x, int(self.pixels*2), indx_nonan.astype(int))
        indy = np.where(bool_arr_y, int(self.pixels*2), indy_nonan.astype(int))

        # Combining indx and indy into matrix of 2-element arrays (x and y coordinates)
        comb_mat = np.stack((indx, indy), axis=2)

        final_time = t.time() - init_time
        if self.print_stats: print(f'Calculating translated pixels: {round(final_time, 3)} seconds')

        # Calculating repeated coordinates and their counts in comb_mat
        init_time = t.time()
        # repetitions, counts = np.unique(ar=comb_mat.reshape(-1, 2), axis=0, return_counts=True)
        stacked_mat = comb_mat.reshape(-1, 2)
        repetitions, counts = IRSCaustics.calc_uniques(stacked_mat)
        final_time = t.time() - init_time
        if self.print_stats: print(f'Finding pixel repetitions and counts: {round(final_time, 3)} seconds')

        # Iterating through the array of counts to find the number of times each coordinate was repeated and increment that coordinate magnification by 1
        init_time = t.time()
        # shape = np.shape(self.X)
        magnifications = np.zeros(shape=(self.pixels, self.pixels), dtype=np.int64)
        self.magnifications = IRSCaustics.calc_mags(self.pixels, magnifications, repetitions, counts) / self.rays_per_pixel**2
        final_time = t.time() - init_time
        if self.print_stats: print(f'Incrementing pixel magnifications based on counts and repetitions: {round(final_time, 3)} seconds')

        self.magnifications = np.flip(self.magnifications, axis=0)

        # Calculating MulensModel analytic caustic curves
        init_time = t.time()
        if self.show_mm:
            self.mm_x, self.mm_y = self.calc_mm_caustics()
            final_time = t.time() - init_time
            if self.print_stats: print(f'Calculating MulensModel analytic caustic curves: {round(final_time, 3)} seconds')

        # Calculating analytic magnification maps for single lens
        init_time = t.time()
        if self.show_dev:
            self.a_mags = self.calc_a_mags()

            # Calculating deviation from magnification of biggest mass
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                self.delta = (self.magnifications - self.a_mags) / self.a_mags
            final_time = t.time() - init_time
            if self.print_stats: print(f'Calculating analytic magnification map: {round(final_time, 3)} seconds')

        # Replacing all 0 values with 0.1 to plot in log10 space
        magnifications_log = np.where(self.magnifications == 0, 0.1, self.magnifications)
        magnifications_log = np.log10(magnifications_log)
        self.magnifications_log = np.where(magnifications_log == -1, 0, magnifications_log)

        # Plotting magnifications
        if self.show_plot:
            init_time = t.time()
            self.fig_c = self.plot_mags_map()
            final_time = t.time() - init_time
            if self.print_stats: print(f'Plotting magnification map: {round(final_time, 3)} seconds')

        if self.save_plot:
            init_time = t.time()
            self.fig_c.savefig(f'../figures/{self.import_file}.png', dpi=500)
            final_time = t.time() - init_time
            if self.print_stats: print(f'Saving magnification map: {round(final_time, 3)} seconds')

        if self.file_save:
            init_time = t.time()
            self.write_to_file()
            final_time = t.time() - init_time
            if self.print_stats: print(f'Saving magnification data in file ./Mag Map Data/{self.import_file}.txt: {round(final_time, 3)} seconds')
        
        end_time = t.time() - begin_time
        if self.print_stats:
            print('---------------------')
            print(f'Total time: {round(end_time, 3)} seconds')

        return self.magnifications

    @staticmethod
    @nb.jit(nb.int64[:, :](nb.int32, nb.int64[:, :], nb.int64[:, :], nb.int64[:]), nopython=True, fastmath=True)
    def calc_mags(pixels, magnifications, repetitions, counts):
        '''
        Calculates magnifications using Numba's jit method with C-like compiling for faster computing.

        Parameters
        ----------
        pixels : int32
        magnifications : 2D int64 Numpy array
        repetitions : 2D int64 Numpy array
        counts : 1D int64 Numpy array

        Returns
        -------
        magnifications : 2D int64 Numpy array
        '''
        for i, count in enumerate(counts):
            if pixels*2 not in repetitions[i]:
                magnifications[repetitions[i, 0], repetitions[i, 1]] += count

        return magnifications

    @staticmethod
    @nb.jit(nb.types.Tuple((nb.int64[:, :], nb.int64[:]))(nb.int64[:, :]), parallel=True, nopython=True, fastmath=True)
    def calc_uniques(sorted_mat):
        n, m = sorted_mat.shape
        assert m >= 0

        isUnique = np.zeros(n, np.bool_)
        uniqueCount = 1
        if n > 0:
            isUnique[0] = True
        for i in nb.prange(1, n):
            isUniqueVal = False
            for j in range(m):
                isUniqueVal |= sorted_mat[i, j] != sorted_mat[i-1, j]
            isUnique[i] = isUniqueVal
            uniqueCount += isUniqueVal

        uniqueValues = np.empty((uniqueCount, m), np.int64)
        duplicateCounts = np.zeros(len(uniqueValues), np.int64)

        cursor = 0
        for i in range(n):
            cursor += isUnique[i]
            for j in range(m):
                uniqueValues[cursor-1, j] = sorted_mat[i, j]
            duplicateCounts[cursor-1] += 1

        return uniqueValues, duplicateCounts

    @staticmethod
    def sort_mat(mat):
        n, m = mat.shape

        for i in range(m):
            kind = 'stable' if i > 0 else None
            mat = mat[np.argsort(mat[:,m-1-i], kind=kind)]

        return mat
    
    @staticmethod
    @nb.jit(nb.bool_[:, :](nb.float64, nb.float64[:, :], nb.float64[:, :]), parallel=True, nopython=True, fastmath=True)
    def create_annulus(annulus, X_R, Y_R):
        '''
        Calculates the coordinates of the rays within the annulus.

        Parameters
        ----------
        annulus : float
            Angular width of annulus in theta_E
        X_R : NxN NDArray
            X meshgrid of all the rays in the image plane
        Y_R : NxN NDArray
            Y meshgrid of all the rays in the image plane

        Returns
        -------
        ray_valid_coords : NxN NDArray
            Boolean array of all the coordinates within the annulus
        '''
        # Defining inner and outer radii of the annulus
        lower_bound = 1.0-annulus/2
        upper_bound = 1.0+annulus/2

        # Calculating the distance from the origin squared of each ray
        ray_distance_squared = X_R**2 + Y_R**2

        # Finding which ray indices are within the upper and lower bound of the annulus
        ray_valid_coords = np.logical_and(ray_distance_squared >= lower_bound**2, ray_distance_squared <= upper_bound**2)

        return ray_valid_coords

    def calc_mm_caustics(self):
        '''
        Calculates MulensModel caustics (only for two lenses).

        Parmameters
        -----------
        None

        Returns
        -------
        mm_x : 1x5000 Numpy array
            Array of x coordinates of calculated caustics
        mm_y : 1x5000 Numpy array
            Array of y coordinates of calculated caustics
        '''
        if self.L == 2:
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
            if vhat[0] > 0 and vhat[1] > 0:
                alpha = np.arctan(vhat[1]/vhat[0])
            elif vhat[0] < 0 and vhat[1] > 0:
                alpha = np.pi + np.arctan(vhat[1]/vhat[0])
            elif vhat[0] < 0 and vhat[1] < 0:
                alpha = np.pi + np.arctan(vhat[1]/vhat[0])
            elif vhat[0] > 0 and vhat[1] < 0:
                alpha = np.arctan(vhat[1]/vhat[0])
            elif vhat[0] > 0 and vhat[1] == 0:
                alpha = 0
            elif vhat[0] == 0 and vhat[1] > 0:
                alpha = np.pi/2
            elif vhat[0] < 0 and vhat[1] == 0:
                alpha = np.pi
            elif vhat[0] == 0 and vhat[1] < 0:
                alpha = -np.pi/2

            # Calculating distance between lenses (s)
            s = np.linalg.norm(v)

            # Calculating mass ratio between lenses (q)
            q = self.lens_att[small_mass, 3] / self.lens_att[big_mass, 3]

            # Initializing MulensModel Caustics class
            model = mm.Caustics(q=q, s=s)

            # Initializing parameter dictionary
            self.param_dict = {'q': q, 's': s, 'alpha': alpha}
            
            # Retrieving x and y points of caustics
            caustic_points = np.array(model.get_caustics(n_points=5000)).transpose()

            # Creating rotation matrix (to rotate caustic points into correct binary axis)
            cos, sin = np.cos(alpha), np.sin(alpha)
            Rot = np.array([[cos, -sin], [sin, cos]])

            # Initializing array of rotated caustic points (2x5000)
            rotated_caustic_points = np.zeros(shape=np.shape(caustic_points))

            # Calculating dot product of rotation matrix and each position calculated by MulensModel
            for i, pos in enumerate(caustic_points):
                rotated_caustic_points[i] = np.dot(Rot, pos)

            # Extracting array of x and array of y values
            mm_x = rotated_caustic_points.transpose()[0]
            mm_y = rotated_caustic_points.transpose()[1]
        
        else:
            self.show_mm = False
            mm_x, mm_y = 0, 0
            warnings.warn(f'MulensModel Caustics only valid for 2 lenses. Got {self.L} lenses.', SyntaxWarning)

        return mm_x, mm_y

    def calc_a_mags(self):
        '''
        Calculates analytic magnification map for a binary lens:
            A(u) = (u^2 + 2)/(u * sqrt(u^2 + 4))

        Parameters
        ----------
        None

        Returns
        -------
        a_mags : NxN NDArray
            Matrix of magnifications for each pixel due to analytic function
        '''
         # Which lens is the biggest
        center_lens = np.argmax(self.lens_att[:, 3])

        # Position of biggest lens
        pos = self.lens_att[center_lens, :2]

        # Calculating seperations from biggest lens for each pixel
        u = np.sqrt((self.X - pos[0])**2 + (self.Y - pos[1])**2)

        # Calculating analytic magnifications for each pixel
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            a_mags = (u**2 + 2) / (u * np.sqrt(u**2 + 4))

        return a_mags

    def plot_mags_map(self):
        '''
        Plots the magnification map as a filled contour plot.

        Parameters
        ----------
        None

        Returns
        -------
        matplotlib.figure.Figure
        '''
        # Initializing figure
        if self.show_axes:
            fig = plt.figure('Magnification Map and Caustics', figsize=(10, 8))
        else:
            fig = plt.figure('Magnification Map and Caustics', figsize=(10, 10))
        ax = fig.add_subplot()

        # Finding translated x and y coordinates after zoom
        x_lower_bound = int(self.pixels/2) - m.ceil(self.zoom[0]/(2*self.ang_res))
        x_upper_bound = int(self.pixels/2) + m.ceil(self.zoom[0]/(2*self.ang_res))

        y_lower_bound = int(self.pixels/2) - m.ceil(self.zoom[1]/(2*self.ang_res))
        y_upper_bound = int(self.pixels/2) + m.ceil(self.zoom[1]/(2*self.ang_res))

        x_zoomed = self.X[y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]
        y_zoomed = self.Y[y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]

        if self.show_dev:
            # Zooming into seperations based on zoom
            delta_zoomed = self.delta[y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]

            # Plotting deviation map with set view max
            # plot = ax.contourf(x_zoomed, y_zoomed, delta_zoomed, cmap='bone', levels=100)
            pic = ax.imshow(delta_zoomed, cmap=self.cmap, vmin=0, extent=[-self.zoom[0]/2, self.zoom[0]/2, -self.zoom[1]/2, self.zoom[1]/2])

        else:
            magnifications_log_zoomed = self.magnifications_log[y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]
            # magnifications_zoomed = self.magnifications[y_lower_bound:y_upper_bound, x_lower_bound:x_upper_bound]

            # Plotting magnification map with set view max
            # plot = ax.contourf(x_zoomed, y_zoomed, magnifications_log_zoomed, cmap='viridis', levels=100, vmin=0)
            # mag_log_zoom_smooth = ndi.gaussian_filter(magnifications_log_zoomed, [1, 1], mode='constant')
            plot = ax.imshow(magnifications_log_zoomed, cmap=self.cmap, vmin=0, extent=[-self.zoom[0]/2, self.zoom[0]/2, -self.zoom[1]/2, self.zoom[1]/2])

            # ax.scatter(self.X, self.Y, s=1, c='red')

        if self.show_lenses:
            for i in range(self.L):
                # Plotting lensing objects
                ax.add_patch(patches.Circle((self.lens_att[i, 0], self.lens_att[i, 1]), self.lens_att[i, 2], color='cyan'))
            
            ax.add_patch(patches.Circle(self.cm_offset, 1, color='cyan', fill=False, linestyle='dashed'))

            # Plotting center of mass
            ax.scatter(self.cm_offset[0], self.cm_offset[1], s=50, marker='+', c='cyan')

        if self.show_mm:
            mm_scat = ax.scatter(self.mm_x, self.mm_y, zorder=10, c='r', s=0.01, alpha=1)
        
        if self.show_axes:
            bar = plt.colorbar(plot)
            if self.show_dev:
                bar.set_label('Deviations')
            else:
                bar.set_label('$log_{10}$ Magnification')

            ax.set_title(f'Magnification Map\n{self.L} lenses')
            ax.set_xlabel('X [$\\theta_E$]')
            ax.set_ylabel('Y [$\\theta_E$]')
       
        else:
            ax.set_position([0, 0, 1, 1])
            ax.axis('off')
        
        ax.axis('scaled')
        ax.set_xlim(-self.zoom[0]/2, self.zoom[0]/2)
        ax.set_ylim(-self.zoom[1]/2, self.zoom[1]/2)

        # Displaying model parameters (if MulensModel is enabled)
        if self.show_mm:
            model_param_str = ''
            for key, value in self.param_dict.items():
                if key == 'alpha':
                    model_param_str += f'$\\{key}$' + f' = {round(np.rad2deg(value), 4)}°'
                else:
                    model_param_str += f'${key}$' + f' = {round(value, 4)}\n'
            
            props = dict(boxstyle='round', facecolor='white', edgecolor='lightgray', alpha=0.8)
            ax.text(0.02, 0.98, model_param_str, va='top', zorder=10, bbox=props, transform=ax.transAxes)

        return fig

    def write_to_file(self):
        '''
        Writes magnifications to text file.
        '''
        np.savetxt(f'../datafiles/{self.import_file}.txt', self.magnifications)

    def calc_caustic_points(self):
        '''
        Calculates caustic points by finding maxima of each row of image.
        '''
        neighborhood = 5

        mag_max = ndi.maximum_filter(self.magnifications_log, neighborhood)
        maxima = self.magnifications_log == mag_max
        mag_min = ndi.minimum_filter(self.magnifications_log, neighborhood)

        threshold = np.amax(self.magnifications_log) * 0.5

        diff = (mag_max - mag_min) > threshold
        maxima[diff == 0] = 0

        labeled, num_objects = ndi.label(maxima)
        slices = ndi.find_objects(labeled)

        x, y = [], []

        for dy, dx in slices:
            x_center = (dx.start + dx.stop - 1)/2
            x.append(x_center)
            y_center = (dy.start + dy.stop - 1)/2    
            y.append(y_center)

        plt.figure(figsize=(10, 8))
        plt.scatter(x, y)

        # points = np.zeros(shape=np.shape(self.magnifications_log), dtype=np.bool_)

        # for i, row in enumerate(self.delta):
        #     points[i, sig.argrelextrema(row, np.greater)[0]] = True
    
        # labeled, num_objects = ndi.label(points)
        # slices = ndi.find_objects(labeled)

        # x, y = [], []

        # for dy, dx in slices:
        #     x_center = (dx.start + dx.stop - 1)/2
        #     x.append(x_center)
        #     y_center = (dy.start + dy.stop - 1)/2    
        #     y.append(y_center)

        # plt.figure(figsize=(10, 8))
        # plt.scatter(x, y)

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
