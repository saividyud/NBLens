from .imports import *
from .IRSMain import IRSMain

class IRSCaustics(IRSMain):
    '''
    Plots the magnification map for L number of lenses.

    Parameters
    ----------
    **ONLY CHOOSE ONE PARAMETER TO PASS IN AT A TIME**

    whole_plane_param_dict : dict
        Parameter dictionary for sampling the whole source plane
        Acceptable parameters below
            **CHOOSE 2 OF THE FOLLOWING 4 PARAMETERS TO FULLY DEFINE SOURCE AND IMAGE PLANES**
            pixels : int
                Number of pixels along a side of the map
            ang_width : float, str = 'auto'
                Angular width of map in terms of Einstein ring radius [theta_e]
            ang_res : float
                Angular resolution of magnification map [theta_e / pixel]
            pixel_density : float
                Density of pixels in the magnification map [pixel / theta_e]
            
            **REQUIRED PARAMETERS**
            rays_per_pixel : int
                Number of rays to shoot along a side of a pixel (total number of rays is squared) [ray / pixel]
            lens_att : 1x4 array
                Lens attributes in format:
                    [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]

    annulus_param_dict : dict
        Parameter dictionary for only sampling within the Einstein ring radius
        Acceptable parameters below
            **CHOOSE 2 OF THE FOLLOWING 4 PARAMETERS TO FULLY DEFINE SOURCE PLANE**
            pixels : int
                Number of pixels along a side of the map
            ang_width : float
                Angular width of map in terms of Einstein ring radius [theta_e]
            ang_res : float
                Angular resolution of magnification map [theta_e / pixel]
            pixel_density : float
                Density of pixels in the magnification map [pixel / theta_e]

            **CHOOSE 1 OF THE FOLLOWING 2 PARAMETERS TO FULLY DEFINE RADIAL SAMPLING**
            dr : float
                Angular seperation of rays radially [theta_e]
            num_r : int
                Number of rays radially [rays]
            
            **CHOOSE 1 OF THE FOLLOWING 2 PARAMETERS TO DEFINE TANGENTIAL SAMPLING**
            dtheta : float
                Angular seperation of rays tangentially [degrees]
            num_theta : int
                Number of rays tangentially [rays]

            **REQUIRED PARAMETERS**
            thickness : float, str = 'auto'
                Thickness of the annulus [theta_e]
            lens_att : 1x4 array
                Lens attributes in format:
                    [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]

    Returns
    -------
    IRSCaustics object
    '''
    def __init__(self, whole_plane_param_dict: dict = None, annulus_param_dict: dict = None):
        if whole_plane_param_dict != None:
            # Checking if rays_per_pixel exists
            if whole_plane_param_dict.get('rays_per_pixel') == None:
                raise AttributeError('Parameter dict must contain rays_per_pixel.')
            
            # Making sure source and image planes are fully defined
            self.param_dict = self.whole_param_checker(whole_plane_param_dict)

            # Adding some more important parameters to dictionary
            self.param_dict.update({'num_rays': int((self.param_dict['pixels']*self.param_dict['rays_per_pixel'])**2)})

            # Deleting thickness if auto was passed into angular width
            self.param_dict.update({'thickness': 0.0})
            del self.param_dict['thickness']
            
            # Changing mode for later use
            self.mode = 'whole'

            # Extracting and naming important parameters
            self.rays_per_pixel = self.param_dict['rays_per_pixel']
            self.num_rays = self.param_dict['num_rays']

        elif annulus_param_dict != None:
            # Making sure source plane and annulus shooting region are fully defined
            self.param_dict = self.annulus_param_checker(self.whole_param_checker(annulus_param_dict))

            # Adding some more important parameters to dictionary
            self.param_dict.update({'num_rays': int(self.param_dict['num_r']*self.param_dict['num_theta'])})

            # Changing mode for later use
            self.mode = 'annulus'

            # Extracting and naming important parameters
            self.thickness = self.param_dict['thickness']
            self.num_r = self.param_dict['num_r']
            self.num_theta = self.param_dict['num_theta']
            self.num_rays = self.param_dict['num_rays']

        # Initializing parent class
        # Creates some class variables: self.pixels, self.ang_width, self.lens_att, self.ang_res, self.L, self.total_M, self.import_file
        super().__init__(pixels=self.param_dict['pixels'], ang_width=self.param_dict['ang_width'], source_att=None, lens_att=self.param_dict['lens_att'])

    def whole_param_checker(self, passed_params: dict):
        # Checking if lens_att exists
        if passed_params.get('lens_att') == None:
            raise AttributeError('Parameter dict must contain lens_att.')
        
        # Extracting pertinent arguments from dictionary
        pixels, ang_width, ang_res, pixel_density = passed_params.get('pixels'), passed_params.get('ang_width'), passed_params.get('ang_res'), passed_params.get('pixel_density')

        # Checking if ang_width should be calculated
        if ang_width == 'auto':
            # If the angular width is automatic, need to calculate the angular width using analytic equations...
            ang_width, thickness = self.ang_width_thickness_calculator(passed_params.get('lens_att'))

        # Need to check if either the pixels or the angular resolution was passed
        if pixels != None and ang_width != None:
            ang_res = ang_width / pixels
            pixel_density = 1 / ang_res
        elif pixels != None and ang_res != None:
            ang_width = pixels * ang_res
            pixel_density = 1 / ang_res
        elif pixels != None and pixel_density != None:
            ang_res = 1 / pixel_density
            ang_width = pixels * ang_res
        elif ang_width != None and ang_res != None:
            pixels = int(ang_width / ang_res)
            pixel_density = 1 / ang_res
        elif ang_width != None and pixel_density != None:
            ang_res = 1 / pixel_density
            pixels = int(ang_width / ang_res)

        full_params = {}
        full_params.update(passed_params)
        full_params.update({'pixels': pixels, 'ang_width': ang_width, 'ang_res': ang_res, 'pixel_density': pixel_density})
        
        if 'thickness' in full_params.keys():
            if full_params['thickness'] == 'auto':
                full_params['thickness'] = thickness
            else:
                pass

        return full_params
    
    def annulus_param_checker(self, passed_params: dict):
        # Checking if thickness exists
        if passed_params.get('thickness') == None:
            raise AttributeError('Parameter dict must contain thickness.')
        
        # Testing dr and num_r
        dr, num_r = passed_params.get('dr'), passed_params.get('num_r')

        if dr != None:
            num_r = passed_params['thickness'] / dr
        elif num_r != None:
            dr = passed_params['thickness'] / num_r

        # Testing dtheta and num_theta
        dtheta, num_theta = passed_params.get('dtheta'), passed_params.get('num_theta')

        if dtheta != None:
            num_theta = 360 / dtheta
        elif num_theta != None:
            dtheta = 360 / num_theta

        full_params = {}
        full_params.update(passed_params)
        full_params.update({'dr': dr, 'num_r': num_r, 'dtheta': dtheta, 'num_theta': num_theta})

        return full_params

    def ang_width_thickness_calculator(self, lens_att: list):
        '''
        Calculates the angular width of sampling region and thickness of shooting region using Equations 8-10 in https://arxiv.org/pdf/astro-ph/0505363.

        Parameters
        ----------
        lens_att : list
            Passed through via passed_params dict

        Returns
        -------
        ang_width : float
            Angular width of ray sampling region
        thickness : float
            Thickness of annulus shooting region
        '''
        lens_att = np.array(lens_att)
        
        # Number of lens objects
        L = np.shape(lens_att)[0]

        # Initializing maximum distance from origin list
        max_dist_rot = []
        max_dist = []
        self.points = np.zeros(shape=(L-1, 4, 2))

        # Initializing padding
        padding = 2

        # Calculating the center of masses of each individual star-planet combination
        lens_CMs_rot = np.zeros(shape=(L-1, 2))
        
        for i, lens in enumerate(lens_att[1:]):
            two_lenses = np.array([lens_att[0], lens])

            # Finding which index the bigger mass was passed in (primary lens)
            big_mass = np.where(two_lenses[:, 3] == np.max(two_lenses[:, 3]))[0][0]

            # Secondary lens
            small_mass = int(not big_mass)

            # Defining unit source vector
            uhat = [1, 0]

            # Defining unit binary axis vector (from primary to secondary lens)
            v = [two_lenses[small_mass, 0] - two_lenses[big_mass, 0], two_lenses[small_mass, 1] - two_lenses[big_mass, 1]]
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
            q = two_lenses[small_mass, 3] / two_lenses[big_mass, 3]

            # Calculating center of mass in binary axis
            lens_CMs_rot[i, 0] = q*s / (1 + q)

            # Creating rotation matrix
            cos, sin = np.cos(alpha), np.sin(alpha)
            Rot = np.array([[cos, -sin], [sin, cos]])

            # Finding cos(phi_c) - Equation 10
            cos_phi_c = 3/4 * (s + s**-1) * (1 - np.sqrt(1 - 32/9*(s + s**-1)**-2))

            # Finding phi_c through arccos
            phi_c = np.arccos(cos_phi_c)

            # Finding positive x and negative x positions - Equation 8
            x_pos_rot = q / ((1 + s)*(1 + s**-1)) - lens_CMs_rot[i, 0]
            x_neg_rot = -q / ((1 - s)*(1 - s**-1)) - lens_CMs_rot[i, 0]

            # Finding positive y and negative y positions - Equation 9
            y_pos_rot = (2*q * np.abs(np.sin(phi_c)**3)) / (s + s**-1 - 2*cos_phi_c)**2 - lens_CMs_rot[i, 1]
            y_neg_rot = -(2*q * np.abs(np.sin(phi_c)**3)) / (s + s**-1 - 2*cos_phi_c)**2 - lens_CMs_rot[i, 1]

            # Storing the maximum distance from the origin of the caustic
            max_dist_rot.append(np.max(np.abs(np.array([x_pos_rot, x_neg_rot, y_pos_rot, y_neg_rot]))))

            # Rotating bounding points back from binary axis to inertial axis
            x_pos = np.dot(Rot, np.array([x_pos_rot, 0]))
            x_neg = np.dot(Rot, np.array([x_neg_rot, 0]))

            y_pos = np.dot(Rot, np.array([0, y_pos_rot]))
            y_neg = np.dot(Rot, np.array([0, y_neg_rot]))

            # Finding maximum distance from the origin and storing it
            self.points[i, :, :] = np.array([x_pos, x_neg, y_pos, y_neg])
            max_dist.append(np.max(np.abs(self.points[i, :, :])))

        # Now finding the maximum of the list of maximum distances from the origin of each caustic with a padding
        ang_width = 2*max(max_dist) * padding

        # Calculating thickness from maximum distance from origin
        thickness = 2*max(max_dist_rot) * padding

        return ang_width, thickness

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
            Saves magnification map data in CSV file in ../datafiles/{filename}.csv
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
        IRSCaustics.calc_mags_whole(100, magnifications, init_rand_arr, count_rand_arr)
        IRSCaustics.calc_mags_annulus(100, magnifications, init_rand_arr, count_rand_arr, 0.1, 1.0)

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

        if self.mode == 'whole':
            # Defining vector of x coordinates for each ray
            x_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)
            
            # Defining vector of y coordinates for each ray
            y_rays = np.linspace(-self.ang_width/2 - self.ang_res/2 + self.ang_res/(2*self.rays_per_pixel), self.ang_width/2 + self.ang_res/2 - self.ang_res/(2*self.rays_per_pixel), self.rays_per_pixel*self.pixels)
            
            # Calculating meshgrid of X and Y coordinates
            self.X, self.Y = np.meshgrid(x_rays, y_rays)
            
        elif self.mode == 'annulus':
            # Defining vector of radial coordinates for each ray
            rs = np.linspace(1-(self.thickness/2), 1+(self.thickness/2), self.num_r).reshape(-1, 1)

            # Defining vector of tangential coordinates for each ray
            thetas = np.linspace(0, 2*np.pi - (2*np.pi/self.num_theta), self.num_theta).reshape(1, -1)

            # Calculating meshgrid of X and Y coordinates
            self.X = np.dot(rs, np.cos(thetas))
            self.Y = np.dot(rs, np.sin(thetas))

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
        stacked_mat = comb_mat.reshape(-1, 2)
        repetitions, counts = IRSCaustics.calc_uniques(stacked_mat)
        final_time = t.time() - init_time
        if self.print_stats: print(f'Finding pixel repetitions and counts: {round(final_time, 3)} seconds')

        # Iterating through the array of counts to find the number of times each coordinate was repeated and increment that coordinate magnification by 1
        init_time = t.time()
        
        magnifications = np.zeros(shape=(self.pixels, self.pixels), dtype=np.int64)
        
        if self.mode == 'whole':
            self.magnifications = IRSCaustics.calc_mags_whole(self.pixels, magnifications, repetitions, counts) / self.rays_per_pixel**2
        elif self.mode == 'annulus':
            # Calculating area of annulus
            A_ann = 2 * np.pi * 1.0 * self.thickness
            
            # Calculating ray density within annulus
            sigma_ann = self.num_rays / A_ann

            # Calculating area of pixel
            A_pix = (self.ang_res * 1.0)**2

            # Calculating magnifications
            self.magnifications = IRSCaustics.calc_mags_annulus(self.pixels, magnifications, repetitions, counts, sigma_ann, A_pix)
        
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
        # self.magnifications_log = self.magnifications

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
    def calc_mags_whole(pixels, magnifications, repetitions, counts):
        '''
        Calculates magnifications for the whole plane of rays using Numba's jit method with C-like compiling for faster computing.

        Parameters
        ----------
        pixels : int32
        magnifications : 2D int64 Numpy array
        repetitions : 2D int64 Numpy array
            Repeating coordinates
        counts : 1D int64 Numpy array
            Number of repetitions for each corresponding coordinate

        Returns
        -------
        magnifications : 2D int64 Numpy array
        '''
        for i, count in enumerate(counts):
            if pixels*2 not in repetitions[i]:
                magnifications[repetitions[i, 0], repetitions[i, 1]] += count

        return magnifications
    
    @staticmethod
    @nb.jit(nb.int64[:, :](nb.int32, nb.int64[:, :], nb.int64[:, :], nb.int64[:], nb.float64, nb.float64), nopython=True, fastmath=True)
    def calc_mags_annulus(pixels, magnifications, repetitions, counts, sigma_ann, A_pix):
        '''
        Calculates magnifications for an annulus of rays using Numba's jit method with C-like compiling for faster computing.

        Parameters
        ----------
        pixels : int32
        magnifications : 2D int64 Numpy array
        repetitions : 2D int64 Numpy array
            Repeating coordinates
        counts : 1D int64 Numpy array
            Number of repetitions for each corresponding coordinate
        sigma_ann : float64
            Ray density within annulus [rays / theta_e^2]
        A_pix : float64
            Area of each pixel in source plane [theta_e^2]

        Returns
        -------
        magnifications : 2D int64 Numpy array
        '''
        # Iterating through all the repitions of pixel coordinates and counts
        for i, count in enumerate(counts):
            if pixels*2 not in repetitions[i]:
                # Calculating local ray density of pixel
                sigma_pix = count / A_pix

                # Incrementing magnification of that pixel as the ratio of ray densities in source and image planes
                magnifications[repetitions[i, 0], repetitions[i, 1]] += sigma_pix / sigma_ann

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
        np.savetxt(f'./datafiles/{self.import_file}.txt', self.magnifications)

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
