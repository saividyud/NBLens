from .imports import *

class IRSFunctions:
    '''
    Contains useful functions regarding microlensing.

    Methods
    -------
    e_ring : Calculates Einstein ring radii with simulation distances and lens masses.
    plot_system_proj : Given orbital parameters, plot 3D system and its 2D coordinates projected onto the plane of the sky (x-y plane)
    get_proj_coords : Given orbital parameters of a 3D system, provide 2D coordinates projected onto the plane of the sky (x-y plane)
    '''

    def progress_bar(iteration, total):
        '''
        Displays a progress bar for iteration information.

        Parameters
        ----------
        iteration : int
            Step of loop
        total : int
            Total number of loops

        Returns
        -------
        None
        '''
        percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(100 * iteration // total)
        bar = '█' * filledLength + '-' * (100 - filledLength)
        print(f'\r |{bar}| {percent}%', end = '\r')
        if iteration == total:
            print(f' |{bar}| {percent}%', end = '\n')

    def e_ring(sim_dist : list | np.ndarray | int | float | astropy.units.quantity.Quantity, M : list | np.ndarray | int | float | astropy.units.quantity.Quantity):
        '''
        Calculates Einstein ring radii with simulation distances and lens masses.

        Parameters
        ----------
        sim_dist : 1x2 list or ndarray
            Relevant distances in kiloparsecs in format:
                [Ds, Dl] -> Units: [kpc, kpc]
        M : 1xN list or ndarray
            Masses of lenses in units of Msun

        Returns
        -------
        theta_E : Astropy.units.quantity.Quantity
            Angular Einstein ring radius in mas
        r_E : Astropy.units.quantity.Quantity
            Einstein ring radius in au
        '''
        # Error handling in inputs
        if isinstance(sim_dist, (list, np.ndarray)):
            if len(sim_dist) == 2:
                for elem in sim_dist:
                    if isinstance(elem, (int, float)):
                        dist_units = False
                    elif isinstance(elem, astropy.units.quantity.Quantity):
                        dist_units = True
                    else:
                        raise TypeError(f'Element in attribute sim_dist must be int, float, or astropy.units.quantity.Quantity. Got {type(elem)}.')
            else:
                raise AttributeError(f'Attribute sim_dist must be of length 2. Got {len(sim_dist)}.')
        else:
            raise TypeError(f'Attribute sim_dist must be a list or numpy.ndarray. Got {type(sim_dist)}.')
        
        if isinstance(M, (list, np.ndarray)):
            for elem in M:
                if isinstance(elem, (int, float)):
                    M_units = False
                elif isinstance(elem, astropy.units.quantity.Quantity):
                    M_units = True
                else:
                    raise TypeError(f'Element in attribute M must be int, float, or astropy.units.quantity.Quantity. Got {type(elem)}.')
                
        elif isinstance(M, (int, float)):
            M_units = False

        elif isinstance(M, astropy.units.quantity.Quantity):
            M_units = True

        else:
            raise TypeError(f'Attribute M must be a list, numpy.ndarray, int, float, or astropy.units.quantity.Quantity. Got {type(M)}.')

        # Define important constants
        G = const.G
        c = const.c

        # If units were not passed in, give units to them
        if not dist_units:
            sim_dist *= u.kpc # [kpc, kpc]
        
        if not M_units:
            M *= u.M_sun # [M_sun]

        # Giving SI units to basic quantities
        D_rel = (1/(1/sim_dist[1] - 1/sim_dist[0])).to(u.m) # [m]
        M = M.to(u.kg) # [kg]
        
        # Defining and calculating angular Einstein ring radius of lens
        theta_e = np.sqrt(4*G*np.sum(M) / (D_rel*c**2))*u.rad

        # Calculating Einstein ring radius of lens
        r_E = sim_dist[1] * theta_e.value

        return theta_e.to(u.mas), r_E.to(u.au)

    def source_profile(ang_res, rad, xc=0.0, yc=0.0, profile_type='uniform', LD=0, supersample=10, scale_factor=1):
        '''
        Initializes brightness levels due to sources using row-wise supersampling.

        Parameters
        ----------
        ang_res : float
            Angular resolution of grid
        rad : float
            Radius of source
        xc : float, optional
            X-position of source
        yc : float, optional
            Y-position of source
        profile_type : str, optional
            Type of distribution for source ('uniform', 'Gauss', 'LD')
        LD : float, optional
            Limb-darkening coefficient (used if profile_type is 'LD')
        supersample : int, optional
            Number of subpixels per pixel dimension (e.g., 10 means 10x10 per pixel)
        scale_factor : int, optional
            Scale factor for the source profile (default is 1, no scaling)

        Returns
        -------
        a : 2D np.ndarray
            Brightness values for each pixel in the grid
        '''
        # Pixel grid size
        pix_width = scale_factor * int(2 * rad / ang_res)
        ang_width = rad

        # Define pixel center coordinates
        x = np.linspace(-ang_width, ang_width, pix_width)
        y = np.linspace(-ang_width, ang_width, pix_width)
        X, Y = np.meshgrid(x, y)

        a = np.zeros_like(X, dtype=np.float64)

        # Supersampling grid offsets
        offset = np.linspace(-0.5, 0.5, supersample, endpoint=False) + 0.5 / supersample
        dx, dy = np.meshgrid(offset, offset)
        dx = dx.flatten()
        dy = dy.flatten()
        subpixel_count = len(dx)

        for i in tqdm(range(X.shape[0])):
            y_center = Y[i, 0]
            x_centers = X[i, :]

            # Shape: (num_pixels, num_subpixels)
            sub_x = x_centers[:, None] + ang_res * dx[None, :]
            sub_y = y_center + ang_res * dy[None, :]

            r2 = (sub_x - xc)**2 + (sub_y - yc)**2
            mask = r2 <= rad**2

            if profile_type == 'uniform':
                a[i, :] = np.sum(mask, axis=1) / subpixel_count

            elif profile_type == 'Gauss':
                brightness = np.exp(-0.5 * r2 / rad**2)
                brightness[~mask] = 0  # zero outside circle
                a[i, :] = np.mean(brightness, axis=1)

            elif profile_type == 'LD':
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    cos_theta = np.sqrt(1 - r2 / rad**2)
                    cos_theta = np.nan_to_num(cos_theta, nan=0.0)
                LD_profile = (1 - LD * (1 - 1.5 * cos_theta)) / (2 * np.pi)
                LD_profile[~mask] = 0
                a[i, :] = np.mean(LD_profile, axis=1)

        # Normalize total brightness
        a = a / np.sum(a)

        return a

    
    def find_cusp_points(caustic_points):
        points = caustic_points.T  # Shape: (400, 2)

        vectors = np.diff(points, axis=0)  # Shape: (399, 2)

        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        unit_vectors = vectors / norms

        dot_products = np.sum(unit_vectors[:-1] * unit_vectors[1:], axis=1)
        dot_products = np.clip(dot_products, -1.0, 1.0)  # Clip to avoid numerical errors
        angles = np.rad2deg(np.arccos(dot_products))  # In degrees

        cusp_indices = np.where(angles > 50)[0] + 1  # +1 to shift to the index of the corner point

        return caustic_points[0, cusp_indices], caustic_points[1, cusp_indices]

    def plot_system_proj(a_pl,phi,i,r_E=None,OMEGA=0,omega=0,theta=0,ax=None,kwargs={},proj_kwargs={'color':'crimson'}):
        """
        Given orbital parameters, plot 3D system and its 2D coordinates projected onto the plane of the sky (x-y plane)
        
        Parameters
        ----------
        a_pl : semimajor axis in AU
        phi : planet orbital inclination with respect to the system's disk in deg
        i : inclination of system disk
        r_E : einstein ring radius in AU
            --if given, all distance values will be given in units of Einstein ring radius
        OMEGA : longitude of ascending node in degrees
            --changes line of apses
        omega : argument of periastron
        theta : phase, the current angular position in the plane of the planet's orbit
            
        ax : ploting axes
        kwargs : planet orbit curve keyword arguments
        proj_kwargs : projection orbit curve keyword arguments 
        
        Returns
        -------
        3D plot of system with projected orbits, the x-y plane is the plane of the sky
        """
        

        x1,x2=ax.get_xlim();y1,y2=ax.get_ylim();z1,z2=ax.get_zlim()
        ax.autoscale()

        # set axis parameters
        ax.set_aspect('equal')
        ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z')
        
        # radian conversions
        phi = np.radians(phi)
        i = np.radians(i)
        incl=i+phi
        
        OMEGA = np.radians(OMEGA) 
        omega = np.radians(omega)
        
        theta=np.radians(theta)
        theta_range = np.linspace(0,2*np.pi,100)
        
        # gives positions in units of the Eintein ring radius
        if r_E != None:
            a_pl = a_pl/r_E 
        
        # plot system  
        if i > np.pi or i < 0 or phi > np.pi or phi < 0:
            return(print('i and phi must be in range (0,180) degrees'))
            
        else:
            x = a_pl*np.cos(theta_range)
            y = a_pl*np.sin(theta_range)
            
            # plot full orbit projected onto xy plane
            planet_proj = ax.plot((x*np.cos(omega)-y*np.sin(omega))*np.cos(OMEGA)-(x*np.sin(omega)+y*np.cos(omega))*np.sin(OMEGA)*np.cos(incl),
                                (x*np.cos(omega)-y*np.sin(omega))*np.sin(OMEGA)+(x*np.sin(omega)+y*np.cos(omega))*np.cos(OMEGA)*np.cos(incl),
                                lw=1,linestyle='--',**proj_kwargs)

            # plot full 3D orbit
            planet_3D = ax.plot((x*np.cos(omega)-y*np.sin(omega))*np.cos(OMEGA)-(x*np.sin(omega)+y*np.cos(omega))*np.sin(OMEGA)*np.cos(incl),
                                (x*np.cos(omega)-y*np.sin(omega))*np.sin(OMEGA)+(x*np.sin(omega)+y*np.cos(omega))*np.cos(OMEGA)*np.cos(incl),
                                (x*np.sin(omega)+y*np.cos(omega))*np.sin(incl),**kwargs)
            
            
            x = a_pl*np.cos(theta)
            y = a_pl*np.sin(theta)
            
            # plot planet's 3D position at given theta
            planet_point_3D = ax.plot((x*np.cos(omega)-y*np.sin(omega))*np.cos(OMEGA)-(x*np.sin(omega)+y*np.cos(omega))*np.sin(OMEGA)*np.cos(incl),
                                    (x*np.cos(omega)-y*np.sin(omega))*np.sin(OMEGA)+(x*np.sin(omega)+y*np.cos(omega))*np.cos(OMEGA)*np.cos(incl),
                                    (x*np.sin(omega)+y*np.cos(omega))*np.sin(incl),marker='o',color='black',**kwargs)
            
            # set axes limits
            x3,x4=ax.get_xlim();y3,y4=ax.get_ylim();z3,z4=ax.get_zlim()
            arr = np.array([x3,x4,y3,y4,z3,z4])
            lim = np.max(np.absolute(arr))
            ax.set_xlim(-lim,lim),ax.set_ylim(-lim,lim),ax.set_zlim(-lim,lim)

    def get_proj_coords(a_pl,phi,i,r_E=None,OMEGA=0,omega=0,theta=0):
        """
        Given orbital parameters of a 3D system, provide 2D coordinates projected onto the plane of the sky (x-y plane)
        
        Parameters
        ----------
        a_pl : semimajor axis in AU
        phi : planet orbital inclination with respect to the system's disk in deg
        i : inclination of system disk
        r_E : einstein ring radius in AU
            --if given, all distance values will be given in units of Einstein ring radius
        OMEGA : longitude of ascending node in degrees
            --changes line of apses
        omega : argument of periastron
        theta : phase, the current angular position in the plane of the planet's orbit
        
        Returns
        -------
        2D positions of planets projected onto the plane of the sky (x-y plane)
        """
        
        
        # radian conversions
        phi = np.radians(phi)
        i = np.radians(i)
        incl=i+phi
        
        OMEGA = np.radians(OMEGA) 
        omega = np.radians(omega)
        
        theta=np.radians(theta)
        
        # gives positions in units of the Eintein ring radius
        if r_E != None:
            a_pl = a_pl/r_E
            
        # calculate coordinates
        if i > np.pi or i < 0 or phi > np.pi or phi < 0:
            return(print('i and phi must be in range (0,180) degrees'))  

        else:
            x = a_pl*np.cos(theta)
            y = a_pl*np.sin(theta)
            
            # planet position at given theta projected onto xy plane
            x_proj = (x*np.cos(omega)-y*np.sin(omega))*np.cos(OMEGA)-(x*np.sin(omega)+y*np.cos(omega))*np.sin(OMEGA)*np.cos(incl)
            y_proj = (x*np.cos(omega)-y*np.sin(omega))*np.sin(OMEGA)+(x*np.sin(omega)+y*np.cos(omega))*np.cos(OMEGA)*np.cos(incl)

            return (x_proj,y_proj)

    try: 
        class _MyBinaryLens(mm.BinaryLensPointSourceWM95Magnification):
            """Custom BinaryLens class potentially for modifications or extensions."""
            def __init__(self, q, s):
                # Initialize the parent MulensModel.BinaryLens
                self.s = s

                parameters = {'q': q, 's': s, 't_0': 0, 't_E': 1, 'u_0': 0.006, 'alpha': 0}
                model_params = mm.ModelParameters(parameters=parameters)

                traj = mm.Trajectory(times=0, parameters=model_params)

                super().__init__(trajectory=traj)

            def get_n_images(self, x, y):
                """
                Calculates the positions of images for a source at position (x, y).
                Relies on the _verify_polynomial_roots method from the parent class.
                """
                # The length of the list of roots corresponds to the number of images
                magnification = self._get_1_magnification(x, y, separation=self.s)

                # Getting image positions with the origin at the planet (0, 0)
                results = np.array(self._verify_polynomial_roots())

                # Shifting image positions to have origin at star
                for i in range(len(results)):
                    results[i] += self.s + 0.0j

                return magnification, results
    except:
        warnings.warn('Test version of MulensModel not found. Some functions may not work.')

    def _separate_circles(points, q, s):
        """
        Separate points into three circles.
        
        Parameters
        ----------
        points : Numpy NDArray
            Array of shape (N, 2) with all points
        q : float
            Planet to star mass ratio
        s : float
            x-coordinate of the off-center circle centroid (y = 0)
        
        Returns
        -------
        off_center : Numpy NDArray
            Points belonging to the off-center circle
        small_origin : Numpy NDArray
            Points belonging to the smaller origin-centered circle (avg radius < 1)
        large_origin : Numpy NDArray
            Points belonging to the larger origin-centered circle (avg radius > 1)
        """
        # Off-center circle parameters
        off_center = np.array([s, 0])
        off_radius = np.sqrt(q)

        # Distance of each point to off-center centroid
        d_off = np.linalg.norm(points - off_center, axis=1)

        # Assign off-center points if their radius is close to expected
        mask_off = np.abs(d_off - off_radius) < off_radius
        off_circle = points[mask_off]

        # Remaining points are origin-centered
        origin_points = points[~mask_off]

        # Compute radii of all origin-centered points
        r_origin = np.linalg.norm(origin_points, axis=1)

        # Split into two groups: below 1 and above 1
        small_mask = r_origin < 1
        small_origin = origin_points[small_mask]
        large_origin = origin_points[~small_mask]

        # (Optional) sanity check: swap if mislabeled
        if np.mean(np.linalg.norm(small_origin, axis=1)) > np.mean(np.linalg.norm(large_origin, axis=1)):
            small_origin, large_origin = large_origin, small_origin

        return off_circle, small_origin, large_origin

    def _image_calculator(u, q, s):
        '''
        Calculates the major and minor radii for an annulus for a given source ring of radius u and binary lens configuration.

        Parameters
        ----------
        u : float
            Radius of source ring in units of Einstein radius
        q : float
            Planet to star mass ratio
        s : float
            Distance between star and planet in units of Einstein radius
        
        Returns
        -------
        y_plus : float
            Major image radius in units of Einstein radius
        y_minus : float
            Minor image radius in units of Einstein radius
        min_mag : float
            Minimum magnification in the source ring
        '''
        # Defining center of mass
        center_of_mass = np.array([q*s / (1 + q), 0])

        # Defining ring slightly larger than magnification map extent
        thetas = np.linspace(0, 2*np.pi, 1000, endpoint=False)
        xs = u * np.cos(thetas)
        ys = u * np.sin(thetas)

        # Initializing MulensModel binary lens
        model = IRSFunctions._MyBinaryLens(q, s)

        # Defining lists to hold all image positions and magnifications
        all_points = []
        ring_magnifications = []

        # Iterating through all source positions in ring
        for x, y in zip(xs, ys):
            # Calculating corresponding magnifications and image positions
            magnification, image_positions = model.get_n_images(x - center_of_mass[0], y - center_of_mass[1])

            # Storing all image positions
            for image_position in image_positions:
                all_points.append([image_position.real, image_position.imag])

            # Storing magnifications
            ring_magnifications.append(magnification)

        # Converting to Numpy arrays
        all_points = np.array(all_points)
        ring_magnifications = np.array(ring_magnifications)

        # Separating points into three critical curves: off-centered circle around planet, smaller origin-centered circle, larger origin-centered circle
        offset_circle, small_circle, large_circle = IRSFunctions._separate_circles(all_points, q, s)

        # Calculating polar coordinates of the points of each circle
        small_rs = np.sqrt(small_circle[:, 0]**2 + small_circle[:, 1]**2)
        small_thetas = np.arctan2(small_circle[:, 1], small_circle[:, 0]) * 180 / np.pi

        big_rs = np.sqrt(large_circle[:, 0]**2 + large_circle[:, 1]**2)
        big_thetas = np.arctan2(large_circle[:, 1], large_circle[:, 0]) * 180 / np.pi

        # Averaging radii of each circle to get y+ and y-
        # y+ is the average radius of the larger origin-centered circle
        # y- is the average radius of the smaller origin-centered circle
        y_plus = np.mean(big_rs)
        y_minus = np.mean(small_rs)

        return y_plus, y_minus, np.min(ring_magnifications)

    def _ang_width_calculator(lens_att: list | np.ndarray, pixels_in_small_source: int = 10, cm_offset: str | tuple | list = [0, 0], cm_translation: tuple | list = [0, 0]):
        '''
        Calculates the angular width of the magnification map using Equations 8-10 in https://arxiv.org/pdf/astro-ph/0505363.

        Parameters
        ----------
        lens_att : list
            Passed through via passed_params dict
        pixels_in_small_source : int, optional
            Number of pixels to cover the diameter of the smallest source. Default is 10
        cm_offset : str or tuple or list, optional
            Center of mass offset to be added to caustic positions. If 'auto', then the center of mass of each star-planet combination is used. Default is [0, 0]
        cm_translation : tuple or list, optional
            Center of mass translation to be applied to caustic positions. Default is [0, 0]

        Returns
        -------
        pixels : int
            Number of pixels along one side of the square magnification map
        ang_width : float
            Angular width of ray sampling region
        (qs, ss) : tuple of lists
            Mass ratios and separations of each star-planet combination
        points : Lx4x2 NDArray
            Locations of caustic cusps in the inertial frame
        (max_source_radius, min_source_radius) : tuple of floats
            Maximum and minimum source radii to be used for the simulation
        '''
        lens_att = np.array(lens_att)
        
        # Number of lens objects
        L = np.shape(lens_att)[0]

        # If radii
        if np.shape(lens_att)[1] == 4:
            mass_ind = 3
        else:
            mass_ind = 2

        # Initializing maximum distance from origin list
        max_dist_rot = []
        max_dist = []
        points = np.zeros(shape=(L-1, 4, 2))

        caustic_sizes = []

        qs = []
        ss = []

        # Calculating the center of masses of each individual star-planet combination
        lens_CMs_rot = np.zeros(shape=(L-1, 2))
        
        if L == 1:
            raise ValueError('Only one lens object detected. Cannot calculate angular width for single lens using this function.')

        else:
            for i, lens in enumerate(lens_att[1:]):
                two_lenses = np.array([lens_att[0], lens])

                # Finding which index the bigger mass was passed in (primary lens)
                big_mass = np.where(two_lenses[:, mass_ind] == np.max(two_lenses[:, mass_ind]))[0][0]

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
                ss.append(s)

                # Calculating mass ratio between lenses (q)
                q = two_lenses[small_mass, mass_ind] / two_lenses[big_mass, mass_ind]
                qs.append(q)

                # Calculating center of mass in binary axis
                lens_CMs_rot[i, 0] = q*s / (1 + q)

                # Creating rotation matrix
                cos, sin = np.cos(alpha), np.sin(alpha)
                Rot = np.array([[cos, -sin], [sin, cos]])

                # Finding cos(phi_c) - Equation 10
                cos_phi_c = 3/4 * (s + s**-1) * (1 - np.sqrt(1 - 32/9*(s + s**-1)**-2))

                # Finding phi_c through arccos
                phi_c = np.arccos(cos_phi_c)

                if cm_offset == 'auto':
                    lens_CM = lens_CMs_rot[i]

                # Finding positive x and negative x positions - Equation 8
                x_pos_rot = q / ((1 + s)*(1 + s**-1)) - lens_CMs_rot[i, 0]
                x_neg_rot = -q / ((1 - s)*(1 - s**-1)) - lens_CMs_rot[i, 0]

                # Finding positive y and negative y positions - Equation 9
                y_pos_rot = (2*q * np.abs(np.sin(phi_c)**3)) / (s + s**-1 - 2*cos_phi_c)**2 - lens_CMs_rot[i, 1]
                y_neg_rot = -(2*q * np.abs(np.sin(phi_c)**3)) / (s + s**-1 - 2*cos_phi_c)**2 - lens_CMs_rot[i, 1]

                # Storing the maximum distance from the origin of the caustic
                max_dist_rot.append(np.max(np.abs(np.array([x_pos_rot, x_neg_rot, y_pos_rot, y_neg_rot]))))

                # Rotating bounding points back from binary axis to inertial axis
                x_pos = np.dot(Rot, np.array([x_pos_rot, 0])) + lens_CM - cm_translation
                x_neg = np.dot(Rot, np.array([x_neg_rot, 0])) + lens_CM - cm_translation

                y_pos = np.dot(Rot, np.array([0, y_pos_rot])) # + lens_CM - cm_translation
                y_neg = np.dot(Rot, np.array([0, y_neg_rot])) # + lens_CM - cm_translation

                # Finding maximum distance from the origin and storing it
                points[i, :, :] = np.array([x_pos, x_neg, y_pos, y_neg])
                max_dist.append(np.max(np.abs(points[i, :, :])))

                width = abs(x_pos_rot - x_neg_rot)
                height = abs(y_pos_rot - y_neg_rot)

                caustic_sizes.append(max(width, height))

        # Finding the maximum distance of a cusp from the origin among all caustics
        max_cusp_distance = np.max(max_dist)

        # Finding maximum caustic size
        max_caustic_size = np.max(caustic_sizes)

        # Finding maximum source radius
        max_source_radius = (10 * max_caustic_size) / 2
        min_source_radius = max_source_radius / 100

        # Adding source diameter as padding
        ang_width_source = 2 * (max_cusp_distance + 2.1*max_source_radius)

        # Angular width regardless of source size
        ang_width_nosource = 4 * max_cusp_distance

        # Choose largest width:
        ang_width = max(ang_width_source, ang_width_nosource)

        # Making map twice as wide as the minimum required
        ang_width *= 2

        # Determining the number of pixels such that the source radius is covered by at least a certain number of pixels
        ang_res = min_source_radius*2 / pixels_in_small_source
        pixels = int(ang_width / ang_res) + 1

        return pixels, ang_width, (qs, ss), points, (max_source_radius, min_source_radius)

    def _annulus_bounds_calculator(ang_width: float, qs: list, ss: list, annulus_center: list | np.ndarray = None):
        '''
        Calculates the annulus bounds and minimum magnification for a given angular width and lens configuration.

        Parameters
        ----------
        ang_width : float
            Angular width of ray sampling region
        qs : list
            Mass ratios of each star-planet combination. For a single lens, pass None
        ss : list
            Separations of each star-planet combination. For a single lens, pass None
        annulus_center : list, optional
            Center of annulus. Default is None, which means no offset is applied

        Returns
        -------
        (y_plus, y_minus) : tuple
            Major and minor image positions
        min_mag : float
            Minimum magnification in the source ring
        '''
        # Number of lens objects
        if qs is None or ss is None:
            L = 1
        else:
            L = len(qs) + 1

        # Finding the maximum distance from the origin of magnification map
        u = 1.05 * ang_width / np.sqrt(2)
        
        if annulus_center is not None:
            annulus_center = -np.array(annulus_center)
            
            if annulus_center[0] > 0 and annulus_center[1] > 0:
                angle = np.arctan(annulus_center[1]/annulus_center[0])
            elif annulus_center[0] < 0 and annulus_center[1] > 0:
                angle = np.pi + np.arctan(annulus_center[1]/annulus_center[0])
            elif annulus_center[0] < 0 and annulus_center[1] < 0:
                angle = np.pi + np.arctan(annulus_center[1]/annulus_center[0])
            elif annulus_center[0] > 0 and annulus_center[1] < 0:
                angle = np.arctan(annulus_center[1]/annulus_center[0])
            elif annulus_center[0] > 0 and annulus_center[1] == 0:
                angle = 0
            elif annulus_center[0] == 0 and annulus_center[1] > 0:
                angle = np.pi/2
            elif annulus_center[0] < 0 and annulus_center[1] == 0:
                angle = np.pi
            elif annulus_center[0] == 0 and annulus_center[1] < 0:
                angle = -np.pi/2

            u = np.sqrt(u**2 + np.linalg.norm(annulus_center)**2 - 2*u*np.linalg.norm(annulus_center)*np.cos(angle - np.pi/4))
        
        # Finding major and minor image positions (y+ > 0 and y- < 0)
        if L == 1:
            # For a single lens
            y_plus = np.abs(0.5 * (u + np.sqrt(u**2 + 4)))
            y_minus = np.abs(0.5 * (u - np.sqrt(u**2 + 4)))

            # Calculate analytic magnification for a single lens
            min_mag = (u**2 + 2) / (u * np.sqrt(u**2 + 4))

        else:
            # For multiple lenses, find the annulus bounds for each star-planet combination and take the bounds that give the largest annulus
            y_pluses = []
            y_minuses = []
            min_mags = []

            for q, s in zip(qs, ss):
                y_plus, y_minus, min_mag = IRSFunctions._image_calculator(u, q, s)
                y_pluses.append(y_plus)
                y_minuses.append(y_minus)
                min_mags.append(min_mag)

            y_plus = np.max(y_pluses)
            y_minus = np.min(y_minuses)

            # Taking the minimum magnification among all star-planet combinations
            min_mag = np.min(min_mags)

        return (y_plus, y_minus), min_mag

    def _num_ray_calculator(pixels: int, ang_width: float, y_plus: float, y_minus: float, min_mag: float, delta: float = 0.01, r_theta_ratio: int | float = 4):
        '''
        Calculates the number of rays to shoot in r and theta directions to achieve the desired error in the image plane.

        Parameters
        ----------
        pixels : int
            Number of pixels along one side of the square magnification map
        ang_width : float
            Angular width of ray sampling region
        y_plus : float
            Major image position
        y_minus : float
            Minor image position
        min_mag : float
            Minimum magnification in the source ring
        delta : float, optional
            Desired error in the source plane in units of Einstein radius. Default is 0.01
        r_theta_ratio : int or float, optional
            Ratio of number of rays in r to number of rays in theta. Default is 4

        Returns
        -------
        num_r : int
            Number of rays in radial direction
        num_theta : int
            Number of rays in angular direction
        '''
        # Width of pixel
        dx = ang_width / (pixels - 1)

        # Area of annulus
        a_ann = (np.pi * (y_plus**2 - y_minus**2))

        # Ray density in annulus
        n_image = 1 / (min_mag * dx**2 * delta**2)

        # Total number of rays
        N_rays = n_image * a_ann

        # Number of rays in theta
        num_theta = int(np.ceil(np.sqrt(N_rays / r_theta_ratio)))

        # Number of rays in r
        num_r = int(np.ceil(r_theta_ratio * num_theta))

        return num_r, num_theta