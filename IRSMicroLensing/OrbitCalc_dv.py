import math as m
import time as t

# Importing numpy Python library
import numpy as np

# Importing matplotlib libraries
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from functools import partial
import mpl_toolkits.axes_grid1
import matplotlib.widgets as w

from tqdm import tqdm

# Custom library for constants of different bodies
from . import CelestialBodyData as cb

d2r = m.pi / 180 # Converting degrees into radians
r2d = 180 / m.pi # Converting radians into degrees

class OrbitPropagator:
    def __init__(self, bodies, dt, t_span, spheres=0, limits=0, dv=0, coes=True, data=False, deg=True, anim=True, save_plot=False, dark_plot=True, title="Plot", label=0):
        '''
        The class OrbitPropagator takes the following inputs:
       
        Parameters
        ----------
        bodies : Numpy NDArray
            An Nx8 dimensional numpy array with N bodies and 8 terms for each body: either classical orbital elements or state vector + mass + central body indicator
                bodes = np.array([[rx1, ry1, rz1, vx1, vy1, vz1, m1, -1],            -1: No central body / is the central body (will read in Cartesian state vector)
                                [a, e, i, ta, aop, raan, m2, x],                      x: Central body is index x
                                ...                                                   NOTE: 7th index term is inputted ONLY if "coes" is TRUE
                                [rxN, ryN, rzN, vxN, vyN, vzN, mN, -1]])

        dt : float
            A small change in time to numerically integrate with (can greatly vary the speed of animation)

        t_span : float
            The full timespan of the animation (not in real time) or 0 if only one period is desired

        spheres : list
            Array of dictionaries from Celestial Body Data to plot spheres instead of points, defaulted to 0
                spheres = [cb.earth, cb.sun, ... cb.pluto] (for example)

        limits : list
            This variable passes in the user's choice of view for the plot. It is defaulted to 0 if there are no inputs. Options detailed below:
                Square view fixed: This option allows the user to stay at one point in space centered around the origin with a specified maximum viewing distance. The format is:
                    [maximum distance from the origin, 's']
                Square view max for one body: This option allows the user to stay at one point in space with the viewing distance set by the maximum value of the selected body. The format is:
                    [body index in 'bodies' matrix, 'f']
                Square view following: This option is the most useful and allows the user to follow a body around at a specified distance away from the body. The format is:
                    [body index in 'bodies' matrix, maximum distance away from the body, 'f']

        dv : list
            Adds a change in velocity that the user specifies for a certain body at a certain time. 'r' is for rectangular, while 's' is for spherical. Format detailed below:
                dv = [magnitude, theta, phi, body index, time, time units, 's']
                dv = [dvx, dvy, dvz, body index, time, time units, 'r']

        coes : bool
            Boolean value for if the input contains Keplerian orbital elements

        data : bool
            Boolean for if a txt file of the states is needed

        deg : bool
            Boolean for if the angles passed in via coes are degrees or radians

        anim : bool
            Never chage; can only animate (more funcionality coming later)

        dark_plot : bool
            Changes plot to dark mode (defaulted to True)

        title : bool
            Pass in the title to the plot that is going to show

        label : list
            1xN dimensional string array with the labels of each body in bodies; defaulted to 0 if no input
        '''
    
        self.bodies = bodies
        self.deg = deg
        self.count = 0
        self.save = save_plot
        self.dark_plot = dark_plot

        for i in self.bodies:                       # Counting how many bodies were passed in
            self.count += 1

        self.mass = np.zeros(self.count)            # Initialize 1D mass array of bodies
        a = 0                                       # Initialize iterating variable
        for i in self.bodies:
            self.mass[a] = i[6]                     # Saving masses into seperate array
            a += 1

        self.mu = self.mass * cb.G                  # Translates masses into mass parameters mu

        self.bodies = np.delete(self.bodies, 6, 1)              # Deletes the 6th index column in each row in bodies: mass

        if t_span != 0:
            self.t_span = t_span                        # Initialize constants from function parameters
        else:
            self.t_span = (((4 * np.pi**2) / self.mu[0]) * (self.bodies[1, 0]**3)) ** 0.5
        
        if coes:                                    # Checking if inputs are Keplerian orbital elements
            self.central_coes()
            self.coes2rv()

        self.bodies = np.delete(self.bodies, 6, 1)              # Deletes the 6th index column in each row in bodies: central body designation

        self.label = label                          # Initialize label variable

        self.dt = dt
        self.steps = int(self.t_span / self.dt)               # Calculate number of steps needed
        self.ets = np.arange(0, self.t_span, self.dt)         # Set ets array
        self.title = title                          # Sets the title
        self.spheres = spheres
        self.limits = limits

        try:                                        # Checking if dv was passed in
            self.dv = iter(dv)
        except TypeError:
            self.dv_bool = False
        else:
            self.dv_bool = True

            self.dv_body = int(dv[3])                    # Extracting body index that the dv is being applied to
            self.dv_time = float(dv[4])                    # Extracting time that the dv is being applied to
            self.dv_time_units = dv[5]              # Extracting time units that the dv is being applied to
            self.dv_type = dv[6]                    # Extracting type of dv passed in, either spherical or rectangular

            if self.dv_time_units == 'minutes':             # Converting passed time when dv is being applied to seconds
                self.dv_time *= 60
            elif self.dv_time_units == 'hours':
                self.dv_time *= 3600
            elif self.dv_time_units == 'days':
                self.dv_time *= 86400
            elif self.dv_time_units == 'years':
                self.dv_time *= 86400 * 365.25

            if self.dv_type == 's':
                self.dv = self.convert_dv(dv[:3])       # Converting spherical coordinates into rectangular in flight frame
            elif self.dv_type == 'r':
                self.dv = np.array([dv[0], dv[1], dv[2]])       # Extracting rectangular coordinates in inertial frame

            self.dv_step = np.abs(self.ets - self.dv_time).argmin()

        self.states = np.zeros((self.steps, self.count, 6))         # Initializing state vector for each body: steps slices, count rows, 6 columns
        self.states[0] = self.bodies                                      # First slice of state vector will be the initial conditions of each body
        self.propagate_orbit()                      # Calling the propagate orbit function to form the rest of the states

        if data:
            self.file_write()
        
        if anim:
            if self.limits == 0:
                self.max_val = np.max(np.abs(self.states[:, :, :3]))                            # Finds the view that the user wants and sets the maximum and minimum values
            elif self.limits[1] == 's':
                self.max_val = self.limits[0]
                self.limits = 0
            elif self.limits[1] == 'f':
                self.max_val = np.max(np.abs(self.states[:, self.limits[0], :3]))
                self.limits = 0
            else:
                self.max_val = 0

            OrbitPlotter(self.count, self.states, self.steps, self.ets, self.max_val, self.limits, self.label, self.dark_plot, self.title, self.save, self.spheres) # Runs the animation

    # Extracts dvs for different bodies and converts them
    def extract_dv(self):
        dv_count = np.ndim(self.dv)

    # Converts given dv in spherical coordinates with respect to flight path into Cartesian coordinates with respect to flight path
    def convert_dv(self, dv):
        dv_float = []
        for val in dv:
            dv_float.append(float(val))

        dv_xyz = np.array([
            dv_float[0] * m.cos(dv_float[2] * d2r) * m.cos(dv_float[1] * d2r),
            dv_float[0] * m.cos(dv_float[2] * d2r) * m.sin(dv_float[1] * d2r),
            dv_float[0] * m.sin(dv_float[2] * d2r)
        ])
        return dv_xyz

    # Now calculating acclerations of each body
    '''
    The function accelerations takes the following inputs:
        self: To be a part of the whole OrbitPropagator class
        slice: Taking one state matrix at a specific moment in time, it is an Nx6 dimensional matrix for N bodies and three positions, three velocities
        t: Defaulted to 0 and not used; only for satisfying Runge-Kutta numerical calculations
    
    Returns: 
        A single Nx6 dimensional matrix: [[vx_1, vy_1, vz_1, ax_1, ay_1, az_1],
                                          [vx_2, vy_2, vz_2, ax_2, ay_2, az_2],
                                          ...
                                          [vx_N, vy_N, vz_N, ax_N, ay_N, az_N]]
    '''
    def accelerations(self, slice, t=0):
        self.rs = np.delete(slice, [3, 4, 5], 1)                   # Deletes the velocity values from the slice matrix and results in an Nx3 matrix of only positions
        self.a = np.zeros((self.count, 3))                          # Initialize Nx3 array of accelerations, with N rows for bodies and 3 columns for acceleration in each direction
        self.slice_dot = np.zeros((self.count, 6))                  # Initializing Nx6 array of velocities and acclerations

        for i in range(self.count):                             # Calculating the sum of acclerations on each body
            self.a[i] = self.sum_accel(i)

        self.slice_dot[:, :3] = slice[:, 3:]                       # Extracts all rows and the last three columns of slice and puts it into first three columns of state dot
        self.slice_dot[:, 3:] = self.a                              # Extracts all rows and columns of a and puts into last three columns of slice dot
        return self.slice_dot

    # Function for summation in acceleration calculation
    def sum_accel(self, i):
        self.sums = 0                                               # Initializing sum variable
        for k in range(self.count):                                 # Iterating through the bodies
            if k != i:                                              # Checking if body is not itself
                self.partial = -self.mu[k] * (self.rs[i] - self.rs[k]) / (np.linalg.norm(self.rs[i] - self.rs[k]))**3     # Calculating one term of summation
                self.sums = self.sums + self.partial                # Adding one term to the total sum
        return self.sums
    
    # Function for numerical integration using Runge-Kutta 4-step method
    '''
    The function rk4_step takes in the following parameters:
        self: To be a part of the OrbitPropagator class
        f: The function to integrate (acclerations)
        y: A starting state (Nx6 dimension matrix with N rows for bodies and 6 columns for state)
        t: A specific time for the function (not used in this function)

    Returns:
        Weighted average of values for the true value via integration, an Nx6 dimension matrix with first three columns positions and last three columns velocities
    '''
    def rk4_step(self, f, y, t, h):
        k1 = f(y, t)                            # Returns the derivative of the state vector at that time
        k2 = f(y + 0.5*k1*h, t + 0.5*h)         # Returns a small step forward in time
        k3 = f(y + 0.5*k2*h, t + 0.5*h)         # Returns a small step forward in time
        k4 = f(y + k3*h, t + h)                 # Returns the last part needed for Runge-Kutta

        return y + (h / 6.0) * (k1 + (2*k2) + (2*k3) + k4)
    
    # Propagating the actual orbit using iterative integration
    '''
    The function propagate_orbit takes in the following parameters:
        self: To be a part of the OrbitPropagator class
    '''
    def propagate_orbit(self):
        print('PROPAGATING ORBIT...')
        time0 = t.time()
        for step in tqdm(range(self.steps - 1)):                  # Finding the state in the next point in time using the previous slice    
            if self.dv_bool:             
                if (step == self.dv_step) & (self.dv_type == 'r'):
                    self.states[step, self.dv_body, 3:] += self.dv
            
            self.states[step + 1] = self.rk4_step(self.accelerations, self.states[step], 0, self.dt)
            
            if self.dv_bool: 
                if (step == self.dv_step) & (self.dv_type == 's'):
                    self.dv_inert = self.calc_dv(step)
                    self.states[step, self.dv_body, 3:] += self.dv_inert
                    self.states[step + 1] = self.rk4_step(self.accelerations, self.states[step], 0, self.dt)

            # self.printProgressBar(step, self.steps - 2, prefix='Propagation:', suffix='Complete')
        print('Computational time: ' + str(round(t.time() - time0, 4)) + ' seconds') 
        print('Time steps: ' + str(self.steps))
        print('Integrations: ' + str(self.count * self.steps))

    # Finding vectors from point before the dv is applied to point when dv must be applied to potential point after the time that dv should be applied
    def calc_dv(self, step):
        p1 = self.states[step - 1, self.dv_body, :3]            # Extracting relevant positions from calculated states
        p2 = self.states[step, self.dv_body, :3]
        p3 = self.states[step + 1, self.dv_body, :3]

        v12 = p2 - p1                                           # Calculating vector from point 1 to point 2
        xhat = p3 - p2                                          # Calculating vector from point 2 to point 3

        zhat = np.cross(v12, xhat)                              # Finding vector perpendicular to orbit plane
        if zhat[2] < 0:
            zhat *= -1                                          # Making perpendicular point in positive z direction

        yhat = np.cross(zhat, xhat)                             # Finding orthogonal vector to flight direction within orbit plane

        xhat = xhat / np.linalg.norm(xhat)                      # Normalizing three direction unit vectors
        yhat = yhat / np.linalg.norm(yhat)
        zhat = zhat / np.linalg.norm(zhat)

        inert_frame_translation = np.array([xhat, yhat, zhat]).transpose()          # Forming translational matrix to go from flight frame to inertial frame

        return np.dot(inert_frame_translation, self.dv)

    # Writes all the states to a file for exporting
    def file_write(self):
        # Path for my files to go... comment out if sharing
        path = '/Users/saividyud/Documents/VSCode/N-Body Project/Data Files/' + str(self.title) + '.csv'
        
        # Path for sharing... creates a txt file with the same name as the graph with all the data in the same folder as your code
        # path = str(self.title) + '.txt'
        
        self.f = open(path, 'w')          # Opens a file for writing with the same name as the plot
        for slice in self.states:
            for row in slice:
                for i in range(6):
                    if i != 5:
                        self.f.write(str(row[i]) + ',')         # Prints data in Nx6 matrices, with N bodies and 6 elements in each (state vectors) seperated by a space
                    else:
                        self.f.write(str(row[i]))
                self.f.write('\n')                          # Each individual matrix represents a small step in time, seperated by two new line characters
            self.f.write('0\n')
        self.f.close()
    
    # Translating given Keplerian Orbital elements into Cartesian state vector
    '''
    The function coes2rv takes in the following parameters:
        self: To be a part of the OrbitPropagator class

    Returns:
        Nx6 dimensional array of initial states for each body
    '''
    def coes2rv(self):
        state0 = self.bodies                                    # Input beginning states into own variable
    
        sets = max(state0[:, 6]) + 1                           # Finding how many sets of central and small bodies there are

        # Finding the center of mass of each body pairs

        
        kepCount = 0
        for i in state0:
            if i[6] != -1:                                      # Excluding the rows with Cartesian coordinate inputs
                kepCount += 1
        
        # Convert all angles that are degrees into radians
        if self.deg:
            for i in range(kepCount + 1):
                if state0[i, 6] != -1:                          # Excluding the rows with Cartesian coordinate inputs
                    a = 2
                    while a < 6:
                        state0[i, a] *= d2r
                        a += 1

        E = np.zeros(self.count)                                  # Initialize 1xM array of eccentric anomalies, one for each body
        r_norm = np.zeros(self.count)                             # Initialize 1xM array of magnitudes of posititions, one for each body
        peri_state = np.zeros((self.count, 6))                    # Initialize Mx6 matrix of perifocal positions and velocities

        a = 0
        for i in state0:
            if i[6] != -1:
                # r_norm = a * (1 - e**2) / (1 + e*m.cos(ta))
                r_norm[a] = i[0] * (1 - i[1]**2) / (1 + i[1]*m.cos(i[3]))                           # Calulating the norm of the position vector for each body
            a += 1

        a = 0
        for i in state0:
            if i[6] != -1:
                # E = 2*m.atan(m.sqrt((1 - e) / (1 + e)) * m.tan(ta / 2.0))
                E[a] = 2 * m.atan(m.sqrt((1 - i[1]) / (1 + i[1])) * m.tan(i[3] / 2.0))              # Calculate Eccentric Anomaly of each body:
            a += 1

        a = 0
        for i in state0:
            if i[6] != -1:
                # r_perif = r_norm * np.array([m.cos(ta), m.sin(ta), 0])
                peri_state[a, :3] = r_norm[a] * np.array([m.cos(i[3]), m.sin(i[3]), 0])                                 # Calculate perifocal positions for each body
                # v_perif = m.sqrt(mu * a) / r_norm * np.array([-m.sin(E), m.cos(E) * m.sqrt(1 - e**2), 0])
                peri_state[a, 3:] = m.sqrt(self.mu[int(i[6])] * i[0]) / r_norm[a] * np.array([-m.sin(E[a]), m.cos(E[a]) * m.sqrt(1 - i[1]**2), 0])              # Calculate perifocal velocities for each body
            a += 1

        a = 0
        for i in state0:
            if i[6] != -1:
                perif2eci = np.transpose(self.eci2perif(i[5], i[4], i[2]))                          # Calculating positions and velocities by using rotation matrices
                state0[a, :3] = np.dot(perif2eci, peri_state[a, :3])
                state0[a, 3:6] = np.dot(perif2eci, peri_state[a, 3:])
            a += 1

        a = 0
        for i in state0:                                    # Now translating solved states to be around their central body
            if i[6] != -1:
                i += state0[int(i[6])]
            a += 1

    # Inerital to perifocal rotation matrix
    '''
    The function eci2perif takes in the following parameters:
        self: To be a part of the OrbitPropagator class
        raan: Takes the right ascension of ascending node for rotation calculation
        aop: Takes argument of periapse for rotation calculation
        i: Takes inclination for rotation calculation

    Returns:
        3x3 dimensional matrix for rotation
    '''
    def eci2perif(self, raan, aop, i):
        row0 = [-m.sin(raan) * m.cos(i) * m.sin(aop) + m.cos(raan) * m.cos(aop), m.cos(raan) * m.cos(i) * m.sin(aop) + m.sin(raan) * m.cos(aop), m.sin(i) * m.sin(aop)]
        row1 = [-m.sin(raan) * m.cos(i) * m.cos(aop) - m.cos(raan) * m.sin(aop), m.cos(raan) * m.cos(i) * m.cos(aop) - m.sin(raan) * m.sin(aop), m.sin(i) * m.cos(aop)]
        row2 = [m.sin(raan) * m.sin(i), -m.cos(raan) * m.sin(i), m.cos(i)]

        return np.array([row0, row1, row2])    

    # Calculates the Keplerian orbital elements of the central bodies and replaces the state vectors with translated Keplerian orbital elements
    '''
    The function coes2rv takes in the following parameters:
        self: To be a part of the OrbitPropagator class

    Returns:
        Nx6 dimensional array of initial states for each body (all coes)
    '''
    def central_coes(self):
        pass

    # Progress bar function: not written by me all credit goes to Greenstick on Quora: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    def printProgressBar (self, iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        # time.sleep(1)
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        # print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total:
            print(f'{prefix} |{bar}| {percent}% {suffix}', end = '\n')          # This line was edited by me

class OrbitPlotter:
    def __init__(self, count, states, steps, ets, max_val, limits, label, dark_plot, title, save, spheres):
        self.count = count
        self.states = states
        self.steps = steps
        self.ets = ets
        self.max_val = max_val
        self.limits = limits
        self.label = label
        self.title = title
        self.save = save
        self.spheres = spheres
        self.dark_plot = dark_plot

        if self.dark_plot:
            plt.style.use('dark_background')

        self.lines = []

        self.plot_orbits()

        self.comp_times[self.comp_times == 0] = np.nan
        print('Average time to plot each frame: ' + str(round(np.nanmean(self.comp_times), 4)) + ' seconds')
        print('Total time of animation: ' + str(round(np.nansum(self.comp_times), 4)) + ' seconds')

    # Plotting orbits with different settings
    '''
    The function plot_orbits takes in the following parameters:
        self: To be a part of the OrbitPropagator class

    Returns:
        An animation of the 3D plot of orbits or choatic body motion
    '''
    def plot_orbits(self):
        self.fig = plt.figure(figsize=(10, 8))               # Readying the plot
        self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')
        self.ax.ticklabel_format(style='sci')

        self.runs = True
        self.forwards = True
        self.frame = 0
        self.fast = 1

        self.widget_setup()

        self.begin_time = t.time()
        self.comp_times = np.zeros(self.steps)
        self.plot_animation = animation.FuncAnimation(fig=self.fig, func=self.animate_func, interval=1, frames=self.frame_func(), cache_frame_data=False)
        self.ax.set_box_aspect([1, 1, 1])                # Setting aspect ratios for the plot
        self.ax.set_aspect('equal')

        if self.save:
            # Saving the Animation
            f = r"/Users/saividyud/Desktop/animate_func.gif"
            writergif = animation.PillowWriter(fps=self.steps/6)
            self.plot_animation.save(f, writer=writergif)
        else:
            plt.show()

    # Setting up the repeated animation function for continuous plotting
    def animate_func(self, frame):
        self.slider.set_val(frame)
        self.ax.clear()                 # Clears the plot to plot new updated lines

        for i in range(self.count):                   # Iterates through each body and plots all the positions until the frame
            if self.spheres == 0 or self.spheres[i] == 0:
                self.ax.scatter(self.states[frame, i, 0], self.states[frame, i, 1], self.states[frame, i, 2], zorder=30)       # Plotting a point for the position
            else:
                _u, _v = np.mgrid[0:2*np.pi:30j, 0:np.pi:30j]                       # Plots the current position as a sphere
                _x = self.spheres[i]['radius'] * np.cos(_u) * np.sin(_v)
                _y = self.spheres[i]['radius'] * np.sin(_u) * np.sin(_v)
                _z = self.spheres[i]['radius'] * np.cos(_v)
                self.ax.plot_surface(_x + self.states[frame, i, 0], _y + self.states[frame, i, 1], _z + self.states[frame, i, 2], cmap=self.spheres[i]['colors'], zorder=-50)
            if self.label != 0:
                self.lines.append(self.ax.plot(self.states[:frame+1, i, 0], self.states[:frame+1, i, 1], self.states[:frame+1, i, 2], zorder=5, label=self.label[i]))
            else:
                self.lines.append(self.ax.plot(self.states[:frame+1, i, 0], self.states[:frame+1, i, 1], self.states[:frame+1, i, 2], zorder=5))
        
        self.ax.set_xlabel(['X (km)'])              # Setting the labels for the axes
        self.ax.set_ylabel(['Y (km)'])
        self.ax.set_zlabel(['Z (km)'])

        # Checking for different limit inputs
        if self.limits != 0:
            self.ax.set_xlim([(self.states[frame, self.limits[0], 0] - self.limits[1]), (self.states[frame, self.limits[0], 0] + self.limits[1])])
            self.ax.set_ylim([(self.states[frame, self.limits[0], 1] - self.limits[1]), (self.states[frame, self.limits[0], 1] + self.limits[1])])
            self.ax.set_zlim([(self.states[frame, self.limits[0], 2] - self.limits[1]), (self.states[frame, self.limits[0], 2] + self.limits[1])])
        else:
            self.ax.set_xlim([-self.max_val, self.max_val])
            self.ax.set_ylim([-self.max_val, self.max_val])
            self.ax.set_zlim([-self.max_val, self.max_val])

        if self.label != 0:
            self.ax.legend()

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
        self.fig.suptitle(self.title + '\nTime = ' + str(np.round(time, 2)) + ' ' + units)

        self.comp_times[frame] = t.time() - self.begin_time
        self.begin_time = t.time()
        # self.ax.view_init(30, -60)

    def frame_func(self):
        while self.runs:
            self.frame += self.fast * (self.forwards - (not self.forwards))
            increment = self.frame + self.fast * (self.forwards - (not self.forwards))
            if increment > 0 and increment < self.steps:
                yield self.frame
            else:
                self.stop()
                yield self.frame

    # Setting up widgets to be used by GUI
    def widget_setup(self):
        # Adding GUI axes above 3D plot axes
        gui_ax = self.fig.add_axes([0.203, 0.875, 0.62, 0.04])

        # Dividing GUI axes
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(gui_ax)

        # Adding 5 additional axes into GUI axes for 5 buttons + slider
        step_back_ax = divider.append_axes('right', size='80%', pad=0.05, facecolor='green')
        back_ax = divider.append_axes('right', size='80%', pad=0.05, facecolor='blue')
        stop_ax = divider.append_axes('right', size='80%', pad=0.05, facecolor='red')
        forward_ax = divider.append_axes('right', size='100%', pad=0.05, facecolor='pink')
        slider_ax = divider.append_axes('right', size='500%', pad=0.1, facecolor='yellow')

        # Initializing buttons and slider with colors and parent axes
        self.step_back_btn = w.Button(gui_ax, label='$\u29CF$', color='dimgray', hovercolor='darkgray')
        self.rewind_btn = w.Button(step_back_ax, label='$\u25C0$', color='dimgray', hovercolor='darkgray')
        self.stop_btn = w.Button(back_ax, label='$\u25A0$', color='dimgray', hovercolor='darkgray')
        self.forward_btn = w.Button(stop_ax, label='$\u25B6$', color='dimgray', hovercolor='darkgray')
        self.step_forward_btn = w.Button(forward_ax, label='$\u29D0$', color='dimgray', hovercolor='darkgray')
        self.slider = w.Slider(slider_ax, '', valmin=0, valmax=self.steps-2, valinit=0, initcolor='none', track_color='dimgray', valstep=1)

        # Mapping buttons and sliders to functions for functionality
        self.forward_btn.on_clicked(self.forward)
        self.rewind_btn.on_clicked(self.backward)
        self.stop_btn.on_clicked(self.stop)
        self.step_back_btn.on_clicked(self.step_back)
        self.step_forward_btn.on_clicked(self.step_forward)
        self.slider.on_changed(self.slider_scroll)

    # Stops the animation
    def stop(self, event=0):
        self.runs = False
        self.plot_animation.pause()

    # Resumes the animation and moves it forward
    def forward(self, event=0):
        if not self.runs:
            self.forwards = True
            self.fast = 1
        elif self.forwards:
            self.fast = self.fast * 2
        elif self.fast > 1:
            self.fast = int(self.fast / 2)
        elif self.fast == 1:
            self.forwards = True

        self.start()

    # Resumes the animation and moves it backward
    def backward(self, event=0):
        if not self.runs:
            self.forwards = False
            self.fast = 1
        elif not self.forwards:
            self.fast = self.fast * 2
        elif self.fast > 1:
            self.fast = int(self.fast / 2)
        elif self.fast == 1:
            self.forwards = False

        self.start()

    # Step time step forward
    def step_forward(self, event=0):
        self.forwards = True
        self.one_step()

    # Step time step backward
    def step_back(self, event=0):
        self.forwards = False
        self.one_step()

    # Increment step by 1 and redraw plot
    def one_step(self):
        increment = self.frame + self.fast * (self.forwards - (not self.forwards))
        if increment < self.steps-1 and increment > 0:
            self.frame = increment
        elif self.frame == 0 and self.forwards:
            self.frame += 1
        elif self.frame == self.steps and not self.forwards:
            self.frame -= 1

        self.animate_func(self.frame)

    # Resumes the animation
    def start(self):
        self.runs = True
        self.plot_animation.resume()

    # Defines the slider functionality by setting the frames value to the slider value
    def slider_scroll(self, event=0):
        self.frame = self.slider.val

