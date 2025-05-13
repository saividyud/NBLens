from .imports import *
from .IRSMain import IRSMain

class IRSPlotter(IRSMain):
    '''
    Plots image and source planes

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
    IRSPlotter object
    '''
    def __init__(self, pixels: int, ang_width: int | float, source_att: list, lens_att: list, *args, **kwargs):
        super().__init__(pixels, ang_width, source_att, lens_att, *args, **kwargs)
        super().calculate(light_curves=False)

    def plot(self, save=False):
        '''
        Plots contour maps of image and source planes.

        Parameters
        ----------
        save : bool
            If the plot should be saved in the current working directory

        Returns
        -------
        None
        '''
        # Calculating lens center of mass
        self.lens_CM = self.calc_CM()

        # Plotting both source and image color maps
        print('7. Plotting both source and image color maps')
        self.ticks = np.linspace(0, 1, 6)

        self.source_lenses = []
        self.image_lenses = []
        
        for i in range(self.L):
            # Plotting lensing objects
            self.source_lenses.append(patches.Circle((self.lens_att[i, 0], self.lens_att[i, 1]), self.lens_att[i, 2], color='cyan'))
            self.image_lenses.append(patches.Circle((self.lens_att[i, 0], self.lens_att[i, 1]), self.lens_att[i, 2], color='cyan'))

        # Plot lens Einstein ring radius
        self.source_e_ring = patches.Circle(self.lens_CM, 1, color='cyan', fill=False, linestyle='dashed')
        self.image_e_ring = patches.Circle(self.lens_CM, 1, color='cyan', fill=False, linestyle='dashed')

        self.image_cmap = self.plot_image()
        self.source_cmap = self.plot_sources()

        if save:
            print('8. Saving figures')
            self.source_cmap.savefig('../figures/Source Brightnesses.png', dpi=500)
            self.image_cmap.savefig('../figures/Image Brightnesses.png', dpi=500)

    def plot_sources(self):
        '''
        Plots source plane with brightnesses.

        Parameters
        ----------
        mesh : 1x2 array of NxN Numpy arrays
            Mesh of X and Y values

        Returns
        -------
        matplotlib.figure.Figure object
        '''
        # Initializing figure
        fig = plt.figure('Source Plane', figsize=(10, 8))
        ax = fig.add_subplot()

        # Plotting sources with brightnesses
        plot = ax.imshow(np.flip(self.source_brightnesses, 0), cmap='afmhot', vmin=0, extent=[-self.ang_width/2, self.ang_width/2, -self.ang_width/2, self.ang_width/2])
        bar = fig.colorbar(plot)
        bar.set_ticks(self.ticks)
        
        # Labeling and formatting
        for lens in self.source_lenses:
            ax.add_patch(lens)
        ax.add_patch(self.source_e_ring)
        ax.scatter(self.lens_CM[0], self.lens_CM[1], s=50, marker='+', c='cyan')
        bar.set_label('Normalized Flux')
        ax.axis('scaled')
        ax.set_title(f'Source Brightness\nResolution = {round(self.ang_res, 5)} $\\theta_E$/pixel\nTotal Flux = 1.0')
        ax.set_xlabel('X [$\\theta_E$]')
        ax.set_ylabel('Y [$\\theta_E$]')

        return fig

    def plot_image(self):
        '''
        Plots image plane with brightnesses.

        Parameters
        ----------
        None

        Returns
        -------
        matplotlib.figure.Figure object
        '''
        # Initializing figure
        fig = plt.figure('Image Plane', figsize=(10, 8))
        ax = fig.add_subplot()

        # Plotting sources with brightnesses
        plot = ax.imshow(np.flip(self.image_brightnesses, 0), cmap='afmhot', vmin=0, extent=[-self.ang_width/2, self.ang_width/2, -self.ang_width/2, self.ang_width/2])
        bar = fig.colorbar(plot)
        bar.set_ticks(self.ticks)

        # Labeling and formatting
        for lens in self.image_lenses:
            ax.add_patch(lens)
        ax.add_patch(self.image_e_ring)
        ax.scatter(self.lens_CM[0], self.lens_CM[1], s=50, marker='+', c='cyan')
        bar.set_label('Normalized Flux')
        ax.axis('scaled')
        ax.set_title(f'Image Brightness\nResolution = {round(self.ang_res, 5)} $\\theta_E$/pixel\nTotal Flux = {round(self.image_flux, 4)}')
        ax.set_xlabel('X [$\\theta_E$]')
        ax.set_ylabel('Y [$\\theta_E$]')

        return fig