from .imports import *

class IRSMain(object):
    '''
    Calculates image fluxes and magnifications for gravitational microlensing

    Parameters
    ----------
    pixels : int
        Number of pixels along a side of the detector
    ang_width : float
        Angular width of lens region in terms of Einstein ring radius (theta_e)
    source_att : Sx4 array, with S number of sources
        Each source's attributes in format:
            [x, y, r, type] -> Units: [theta_e, theta_e, theta_e, string]
        Supported types:
            uniform: homogenous circular distribution
            Gauss: Gaussian circular distribution
    lens_att : Lx4 array
        Each lens' attributes in format:
            [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]

    Returns
    -------
    IRSMain object

    Methods
    -------
    .e_ring()
        Calculates Einstein ring radius of lens.
    .calc_CM()
        Calculates the center of mass of lenses.
    .init_grid()
        Creates meshgrid.
    .calc_sources()
        Calculates source plane brightnesses.
    .circle()
        Initializes brightness levels due to sources.
    .calc_source_pixels()
        Calculates the source pixels that each ray lands on after being deflected by lens.
    .trans_ind()
        Calculates the translated indices for each pixel (which source pixel each image pixel lands on).
    .calc_image()
        Maps source brightnesses per pixel to image using translated image pixels.
    '''
    def __init__(self, pixels: int, ang_width: int | float, source_att: list, lens_att: list, *args, **kwargs):
        self.types = ['Gauss', 'uniform']

        self.pixels = pixels # [pixels]
        self.ang_width = ang_width # [theta_e]
        self.source_att = source_att # [theta_e, theta_e, theta_e, string]
        self.lens_att_list = lens_att # [theta_e, theta_e, theta_e, Msun]

        # Calculating angular resolution of detector
        self.ang_res = self.ang_width / self.pixels # [theta_e/pixel]

        # Giving all theta_e entries units
        self.lens_att = np.array(self.lens_att_list)

        # Number of lens objects
        self.L = np.shape(self.lens_att)[0]

        # Total mass of lens objects
        self.total_M = np.sum(self.lens_att[:, 3])

        # Defining import file
        for frame in inspect.stack()[1:]:
            if frame.filename[0] != '<':
                self.import_file = frame.filename.split('/')[-1].split('.')[0]

    def calculate(self, light_curves=False, print_stats=True):
        '''
        Calculates brightnesses and fluxes for source and image.

        Parameters
        ----------
        light_curves : bool
            If image fluxes should be calculated for each passed source seperately

        Returns
        -------
        image_flux : float
            Flux of image
        image_brightnesses : NxN NDArray
            Brightnesses for each pixel

        -or-

        image_flux : 1xS NDArray
            Fluxes of each seperate source
        image_brightnesses : NxNxS NDArray
            Brightness for each pixel for each seperate source
        '''
        self.light_curves = light_curves
        self.print_stats = print_stats
        
        # Creating mesh grid
        init_time = t.time()
        # if self.print_stats: print('1. Creating mesh grid')
        self.X, self.Y = self.init_grid() # theta_e
        print('Initializing grid:', round(t.time() - init_time, 3), 'seconds')

        # Calculating source plane brightnesses per pixel
        # if self.print_stats: print('2. Calculating source plane brightnesses per pixel')
        self.source_brightnesses = self.calc_source()
        print('Calculating source plane brightnesses:', round(t.time() - init_time, 3), 'seconds')
        
        if self.light_curves:
            # Calculating magnification map for lenses
            init_time = t.time()
            caustic = IRSCaustics(pixels=self.pixels, ang_width=self.ang_width, lens_att=self.lens_att_list)
            magnifications = caustic.plot(print_stats=False, show_plot=False)
            print('Calculating magnifications:', round(t.time() - init_time, 3), 'seconds')

            # Calculating image brightnesses by convoluting the magnification map with the source brightness profile
            init_time = t.time()
            self.image_brightnesses = magnifications * self.source_brightnesses
            print('Convoluting magnification map and source profiles:', round(t.time() - init_time, 3), 'seconds')

            # Finding the sum of brightnesses for each different source
            init_time = t.time()
            self.source_bright_sums = np.nansum(np.nansum(self.source_brightnesses, axis=1), axis=1)
            print('Summing source brightnesses:', round(t.time() - init_time, 3), 'seconds')
            init_time = t.time()
            self.image_bright_sums = np.nansum(np.nansum(self.image_brightnesses, axis=1), axis=1)
            print('Summing image brightnesses:', round(t.time() - init_time, 3), 'seconds')

            # Replacing where the sum is 0 with -1
            init_time = t.time()
            source_bright_sums = np.where(self.source_bright_sums == 0, -1, self.source_bright_sums)
            print('Searching for 0 brightness:', round(t.time() - init_time, 3), 'seconds')

            # Calculating image fluxes by replacing wherever the sum is -1 with 0
            init_time = t.time()
            self.image_flux = np.where(self.source_bright_sums == 0, 0, self.image_bright_sums / source_bright_sums)
            print('Calculating image fluxes by replacing 0 brightnesses:', round(t.time() - init_time, 3), 'seconds')
        
        else:
            # Calculating source pixels
            if self.print_stats: print('3. Calculating source pixels')
            self.xs, self.ys = self.calc_source_pixels()

            # Calculating indices of translated pixel after deflection
            if self.print_stats: print('4. Calculating indices of translated pixel after deflection')
            self.indx, self.indy = self.trans_ind()
            
            # Calculating image brightnesses from translated source pixels
            if self.print_stats: print('5. Calculating image brightnesses from translated source pixels')
            self.image_brightnesses = self.calc_image()

            # Calculating source and image flux
            if self.print_stats: print('6. Calculating source and image flux')
            if np.nansum(self.source_brightnesses) == 0:
                # Returning 0 if there are no sources in the image
                self.image_flux = 0
            else:
                # Normalizing the image flux
                self.image_flux = np.nansum(self.image_brightnesses) / np.nansum(self.source_brightnesses)
        
        return self.image_flux, self.image_brightnesses

    def calc_CM(self):
        '''
        Calculates center of mass of lenses.

        Parameters
        ----------
        None

        Returns
        -------
        lens_CM : 1x2 NDArray
            [x, y] positions of lens center of mass
        '''
        CM_sum = 0

        for i in range(self.L):
            CM_sum += self.lens_att[i, :2] * self.lens_att[i, 3]
        
        return CM_sum / self.total_M

    def init_grid(self):
        '''
        Creates meshgrid.

        Parameters
        ----------
        None

        Returns
        -------
        X : NxN Numpy array
            X-values of meshgrid
        Y : NxN Numpy array
            Y-values of meshgrid
        '''
        # Creating mesh grid
        x = np.linspace(-self.ang_width/2, self.ang_width/2, self.pixels)
        y = np.linspace(-self.ang_width/2, self.ang_width/2, self.pixels)
        X, Y = np.meshgrid(x, y)

        return X, Y # theta_e

    def calc_source(self):
        '''
        Calculates source plane brightnesses.

        Parameters
        ----------
        None

        Returns
        -------
        ZS : NxN Numpy array
            Brightnesses for each point in grid
        '''
        if not self.light_curves:
            # Initializing brightness levels (for mesh grid)
            shape = np.shape(self.X)
            ZS = np.zeros(shape)

            # Calculating brightness levels for each source and summing them up
            for pos in self.source_att:
                ZS += self.circle(mesh=[self.X, self.Y], rad=pos[2], xc=pos[0], yc=pos[1], type=pos[3]) # NxN matrix

        else:
            # Initializing brightness levels (for each source)
            shape = (len(self.source_att), np.shape(self.X)[0], np.shape(self.X)[1])
            ZS = np.zeros(shape)

            # Calculating brightness levels for each source and storing them seperately in each matrix
            for i, pos in enumerate(self.source_att):
                ZS[i] = self.circle(mesh=[self.X, self.Y], rad=pos[2], xc=pos[0], yc=pos[1], type=pos[3]) # NxNxS tensor

        return ZS

    def circle(self, mesh, rad, xc=0.0, yc=0.0, type='uniform'):
        '''
        Initializes brightness levels due to sources.

        Parameters
        ----------
        mesh : 1x2 array of NxN Numpy arrays
            Mesh of X and Y values
        rad : float
            Radius of source
        xc : float, optional
            X-position of source
        yc : float, optional
            Y-position of source
        type : str, optional
            Type of distrubution for source (defaulted to Gaussian distribution)
            Supported types:
                uniform: homogenous circular distribution
                Gauss: Gaussian circular distribution

        Returns
        -------
        a : array-like
            Brightness values for each x and y coordinate in mesh
        '''
        X = mesh[0]
        Y = mesh[1]

        r2 = (X - xc)**2 + (Y - yc)**2
        
        if type == 'uniform':
            a = r2 <= rad**2
        elif type == 'Gauss':
            a = np.exp(-r2*0.5 / rad**2)
        
        return a

    def calc_source_pixels(self, X, Y):
        '''
        Calculates the source pixels that each ray lands on after being deflected by lens.

        Parameters
        ----------
        None

        Returns
        -------
        xs : NxN Numpy array
            X-coordinates of translated pixel
        ys : NxN Numpy array
            Y-coordinates of translated pixel
        '''
        # Translating into complex coordinates (source pixels)
        init_time = t.time()
        # z = self.X + self.Y*1j # [theta_e] NxN
        z = np.empty(X.shape, dtype=np.complex128)
        z.real = X
        z.imag = Y
        # print(1, t.time() - init_time)

        init_time = t.time()
        zbar = np.conj(z) # [theta_e] NxN
        # print(2, t.time() - init_time)

        # Translating into complex coordinates (lens coordinates)
        init_time = t.time()
        zm = self.lens_att[:, 0] + self.lens_att[:, 1]*1j # [theta_e] 1xL
        # print(3, t.time() - init_time)

        init_time = t.time()
        epsilon = self.lens_att[:, 3] / self.total_M # [dimensionless] 1xL
        # print(4, t.time() - init_time)
        
        init_time = t.time()
        zmbar = np.conj(zm) # [theta_e] 1xL
        # print(5, t.time() - init_time)
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            init_time = t.time()
            sums = np.zeros(shape=np.shape(X), dtype=np.complex128)
            # print(6, t.time() - init_time)

            init_time = t.time()

            # sum = np.sum(epsilon[:, np.newaxis, np.newaxis] / (zbar - zmbar[:, np.newaxis, np.newaxis]), axis=0)

            # sum = IRSMain.lens_eq(self.L, sums, zbar, zmbar, epsilon)
            zeta = z - IRSMain.lens_eq(self.L, sums, zbar, zmbar, epsilon)
            # print(7, t.time() - init_time)

        # Extracting positions from complex number
        init_time = t.time()
        xs = np.real(zeta) # [theta_e]
        # print(9, t.time() - init_time)

        init_time = t.time()
        ys = np.imag(zeta) # [theta_e]
        # print(10, t.time() - init_time)

        return xs, ys # [theta_e]
    
    @staticmethod 
    @nb.jit(nb.complex128[:, :](nb.int64, nb.complex128[:, :], nb.complex128[:, :], nb.complex128[:], nb.float64[:]), nopython=True, parallel=True, fastmath=True, cache=True)
    def lens_eq(L, sums, zbar, zmbar, epsilon):
        '''
        Calculates the sum in the lens equation using Numba's jit compilation for faster computation.

        Parameters
        ----------
        L : int64
            Number of lenses
        sums : 2D array of complex128
            Zeros array for Numba to compute sums
        zbar : 2D array of complex128
            Complex conjugate of mesh grid
        zmbar : 1D array of complex128
            Complex conjugate of locations of lens masses
        epsilon : 1D array of float64
            Mass fraction of total mass of each lens
        '''
        # sums = np.sum(epsilon / (zbar - zmbar), axis=0)
        for i in range(L):
            sums += epsilon[i] / (zbar - zmbar[i]) # [theta_e]

        return sums

    def trans_ind(self):
        '''
        Calculates the translated indices for each pixel (which source pixel each image pixel lands on).

        Parameters
        ----------
        None

        Returns
        -------
        indx : NxN Numpy array
            Translated source X value for each image pixel
        indy : NxN Numpy array
            Translated source Y value for each image pixel
        '''
        bool_x = (self.xs < -self.ang_width/2) | (self.xs > self.ang_width/2)
        indx = np.where(bool_x, np.nan, self.xs)

        bool_y = (self.ys < -self.ang_width/2) | (self.ys > self.ang_width/2)
        indy = np.where(bool_y, np.nan, self.ys)

        indx = indx/self.ang_res + self.pixels/2
        indy = indy/self.ang_res + self.pixels/2

        # indx = indx/self.ang_res + indx.shape[0]/2
        # indy = indy/self.ang_res + indy.shape[0]/2

        return indx, indy

    def calc_image(self):
        '''
        Maps source brightnesses per pixel to image using translated image pixels.

        Parameters
        ----------
        None

        Returns
        -------
        image_brightnesses : NxN Numpy array
            Brightness per pixel of image
        '''
        # Boolean array on which elements are nan and which are not
        bool_arr = np.isnan(self.indx) | np.isnan(self.indy)

        # Replace these nan values with some arbitrary value to avoid error later (these values are not accessed due to bool_arr)
        indx_nonan = np.where(bool_arr, 0, self.indx)
        indy_nonan = np.where(bool_arr, 0, self.indy)

        # The elements that are nan result in zero, if not nan, then are assigned brightnesses from source_brightnesses matrix
        image_brightnesses = np.where(bool_arr, 0, self.source_brightnesses[indy_nonan.astype(int), indx_nonan.astype(int)])
        
        return image_brightnesses

# Error handling of attributes
    @property
    def pixels(self):
        '''
        Type : int

        Number of pixels per side of square detector.
        '''
        return self._pixels
    
    @pixels.setter
    def pixels(self, val):
        # Checking if pixels is an integer
        if isinstance(val, int):
            # Checking if pixels is greater than 0
            if val > 0:
                self._pixels = val
            else:
                raise ValueError(f'Attribute "pixels" must be positive. Got {val}.')
        else:
            raise TypeError(f'Attribute "pixels" must be an integer. Got {type(val)}.')

    @property
    def ang_width(self):
        '''
        Type : int or float

        Angular width of square detector in terms of $\\theta_E$.
        '''
        return self._ang_width

    @ang_width.setter
    def ang_width(self, val):
        # Checking if ang_width is a float or an int
        if isinstance(val, (float, int)):
            # Checking if ang_width is greater than 0
            if val > 0:
                self._ang_width = val
            else:
                raise ValueError(f'Attribute "ang_width" must be positive. Got {val}.')
        else:
            raise TypeError(f'Attribute "ang_width" must be an integer or float. Got {type(val)}.')

    @property
    def lens_att_list(self):
        '''
        Type : list

        List of lens attributes in format:
            [x, y, r, M] -> Units: [theta_e, theta_e, theta_e, Msun]
        '''
        return self._lens_att_list
    
    @lens_att_list.setter
    def lens_att_list(self, val):
        # Checking if lens_att is a list
        if isinstance(val, list):
            # Checking if lens_att is a 2 dimensional list
            if len(np.shape(val)) == 2:
                # Now checking each element in lens_att
                for lens in val:
                    # Checking if each element in lens_att is a list
                    if isinstance(lens, list):
                        if np.shape(lens)[0] == 4:
                            # Checking if each element inside each lens is a float or an int
                            for n in lens:
                                if isinstance(n, (int, float)):
                                    self._lens_att_list = val
                                else:
                                    raise TypeError(f'Each element in list in attribute "lens_att" must be a float or integer. Got {type(n)}.')
                        else:
                            raise ValueError(f'Each list in attribute "lens_att" must be of length 4. Got {np.shape(lens)}.')
                    else:
                        raise TypeError(f'Each element in attribute "lens_att" must be a list. Got {type(lens)}.')
            else:
                raise ValueError(f'Attribute "lens_att" must 2 dimensional. Got {len(np.shape(val))}.')
        else:
            raise TypeError(f'Attribute "lens_att" must be a list. Got {type(val)}.')
        
    @property
    def source_att(self):
        '''
        Type : list

        List of lens attributes in format:
            [x, y, r, type] -> Units: [theta_e, theta_e, theta_e, string]
        '''
        return self._source_att
    
    @source_att.setter
    def source_att(self, val):
        # Checking if source_att is a list
        if isinstance(val, list):
            # Checking if source_att is a 2 dimensional list
            if len(np.shape(val)) == 2:
                # Now checking each element in source_att
                for source in val:
                    # Checking if each element in source_att is a list
                    if isinstance(source, list):
                        if np.shape(source)[0] == 4:
                            # Checking if each element inside each source is a float or an int
                            for i, n in enumerate(source):
                                # Checking if first three elements each source is a float or an int
                                if i <= 2:
                                    if isinstance(n, (int, float)):
                                        # Checking if third element (radius) is greater than 0
                                        if i == 2:
                                            if n != 0:
                                                pass
                                            else:
                                                raise ValueError(f'Source radius must be nonzero. Got {n}.')
                                    else:
                                        raise TypeError(f'The first three elements in list in attribute "source_att" must be a float or integer. Got {type(n)}.')
                                
                                # Checking if fourth element each source is a string
                                elif i == 3:
                                    if isinstance(n, str):
                                        # Checking if the string is either 'uniform' or 'Gauss
                                        if n in self.types:
                                            self._source_att = val
                                        else:
                                            raise ValueError(f'Passed type is not in supported types: {self.types}. Got {n}.')
                                    else:
                                        raise TypeError(f'The fourth element in list in attribute "source_att" must be a string. Got {type(n)}.')
                        else:
                            raise ValueError(f'Each list in attribute "source_att" must be 4 dimensional. Got {np.shape(source)}.')
                    else:
                        raise TypeError(f'Each element in attribute "source_att" must be a list. Got {type(source)}.')
            else:
                raise ValueError(f'Attribute "source_att" must 2 dimensional. Got {len(np.shape(val))}.')
        
        # Checking if source_att is None type (only passed through IRSCaustics child class)
        elif isinstance(val, type(None)):
            self._source_att = val
        else:
            raise TypeError(f'Attribute "source_att" must be a list. Got {type(val)}.')
