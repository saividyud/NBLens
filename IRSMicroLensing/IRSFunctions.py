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