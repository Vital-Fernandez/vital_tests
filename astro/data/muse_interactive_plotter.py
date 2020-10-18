import numpy as np
from pathlib import Path
from mpdaf.obj import Cube
from astropy.wcs import WCS
from matplotlib import pyplot as plt, rcParams, gridspec

# -------------- Script functions and variables

DEFAULT_FIG = {'figure.figsize': (14, 7),
               'axes.titlesize': 14,
               'axes.labelsize': 14,
               'legend.fontsize': 12,
               'xtick.labelsize': 12,
               'ytick.labelsize': 12}

DEFAULT_SPEC = {'xlabel': r'Wavelength $(\AA)$',
                'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'}

DEFAULT_IMAGE = {'xlabel': r'RA',
                 'ylabel': r'DEC',
                 'title': f'Cube flux slice'}

PERCENT_ARRAY = np.array([0, 80, 90, 95, 99.5, 99.90, 99.99, 100])


def load_muse_cube(file_address):

    cube_array = Cube(filename=str(file_address))
    header_obj = cube_array.data_header

    cube_array.wave.info()
    dw = header_obj['CD3_3']
    w_min = header_obj['CRVAL3']
    nPixels = header_obj['NAXIS3']
    w_max = w_min + dw * nPixels
    wave_array = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave_array, cube_array, header_obj


def compute_muse_flux_images(lineArea_dict, flux_cube, redshift, percent_array):

    # Containers for output data
    lineFlux_dict = {}
    levelFlux_dict = {}
    levelText_dict = {}

    # Loop through the input wavelength regions
    for lineLabel, lineLimits in lineArea_dict.items():

        # Extract cube slice using mpdaf defult tools. This requires the input wavelengths to be
        # on the same scale as in the cube
        line_image = flux_cube.get_image(np.array(lineLimits) * (1 + redshift), subtract_off=True)
        flux_image = line_image.data.data

        # Personal scale for the image flux
        log_flux = np.log10(flux_image)
        flux_contours_i = np.zeros(flux_image.shape)
        idcs_removed = np.logical_or(log_flux < 0.0, np.isnan(log_flux))
        flux_contours_i[~idcs_removed] = log_flux[~idcs_removed]

        # Define image countours based on the flux percentiles
        levelFlux_i = np.percentile(flux_contours_i[flux_contours_i > 0], percent_array)
        levels_text_i = ['None'] * len(levelFlux_i)
        for idx, per in enumerate(percent_array):
            levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$'

        # Store the data
        lineFlux_dict[lineLabel] = flux_contours_i
        levelFlux_dict[lineLabel] = levelFlux_i
        levelText_dict[lineLabel] = levels_text_i

    return lineFlux_dict, levelFlux_dict, levelText_dict


class VoxelPlotter(object):

    """
    This class produces an interative matplotlib window for the muse data cubes. On the left an
    """

    def __init__(self, obj_wave, obj_cube, image_bg, voxel_coord=None, image_fg=None, flux_levels=None,
                 fig_user_conf={}, ax_user_conf={}):

        self.fig = None
        self.ax0, self.ax1, self.in_ax = None, None, None
        self.axlim_dict = {}
        self.grid_mesh = None
        self.cube_data = obj_cube
        self.wave = obj_wave
        self.image_bg = image_bg
        self.image_fg = image_fg
        self.flux_levels = flux_levels
        self.axConf = dict(image={}, spectrum={})

        # Plot Configuration
        defaultConf = DEFAULT_FIG.copy()
        defaultConf.update(fig_user_conf)
        rcParams.update(defaultConf)

        # Figure structure
        self.fig = plt.figure(figsize=(18, 5))
        gs = gridspec.GridSpec(nrows=1, ncols=2, figure=self.fig, width_ratios=[1, 2], height_ratios=[1])
        cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        aee = self.fig.canvas.mpl_connect('axes_enter_event', self.on_enter_axes)

        # Axes configuration
        sky_wcs = WCS(self.cube_data.data_header)
        self.ax0 = self.fig.add_subplot(gs[0], projection=sky_wcs, slices=('x', 'y', 1))
        self.ax1 = self.fig.add_subplot(gs[1])

        imgAxConf = DEFAULT_IMAGE.copy()
        if 'image' in ax_user_conf:
            imgAxConf.update(ax_user_conf['image'])
        self.axConf['image'].update(imgAxConf)

        specAxConf = DEFAULT_SPEC.copy()
        if 'spectrum' in ax_user_conf:
            specAxConf.update(ax_user_conf['spectrum'])
        self.axConf['spectrum'].update(specAxConf)

        # Image mesh grid
        frame_size = self.cube_data.shape
        y, x = np.arange(0, frame_size[1]), np.arange(0, frame_size[2])
        self.grid_mesh = np.meshgrid(x, y)

        # Generate the plot
        self.plot_map_voxel(self.image_bg, voxel_coord, self.image_fg, self.flux_levels, conf_dict=self.axConf)
        plt.show()

        return

    def plot_map_voxel(self, image_bg, voxel_coord=None, image_fg=None, flux_levels=None, conf_dict={}):

        # Image color format
        cmap = plt.cm.gray
        norm = plt.Normalize(image_bg.min(), image_bg.max())
        image_bg_rgb = cmap(norm(image_bg))

        # Emphasize input coordinate
        if voxel_coord is not None:
            idx_j, idx_i = voxel_coord
            image_bg_rgb[idx_j, idx_i, :3] = [1, 0, 0]
            self.ax0.plot(idx_i, idx_j, '+', color='red')

        # Plot background image
        self.ax0.imshow(image_bg_rgb, cmap='gray', vmin=0.0, aspect='auto')

        # Plot contours image
        if image_fg is not None:
            CS3 = self.ax0.contour(self.grid_mesh[0], self.grid_mesh[1], image_fg, levels=flux_levels, alpha=0.5)

        # Voxel spectrum
        if voxel_coord is not None:
            idx_j, idx_i = voxel_coord
            flux_voxel = self.cube_data[:, idx_j, idx_i].data.data
            self.ax1.step(self.wave, flux_voxel)

        # Setting voxel spectrum title as the current coordinate
        conf_dict['spectrum']['title'] = f'Voxel {idx_j} - {idx_i}'

        # Update the axis
        self.ax0.update(conf_dict['image'])
        self.ax1.update(conf_dict['spectrum'])

        return

    def on_click(self, event, mouse_trigger_buttton=3):

        if self.in_ax == self.ax0:

            if event.button == mouse_trigger_buttton:
                # Save axes zoom
                self.save_zoom()

                # Save clicked coordinates for next plot
                idx_j, idx_i = np.rint(event.ydata).astype(int), np.rint(event.xdata).astype(int)
                print(f'Current voxel: {idx_j}-{idx_i} (mouse button {event.button})')

                # Remake the drawing
                self.ax0.clear()
                self.ax1.clear()
                self.plot_map_voxel(self.image_bg, (idx_j, idx_i), self.image_fg, self.flux_levels,
                                    conf_dict=self.axConf)

                # Reset the image
                self.reset_zoom()
                self.fig.canvas.draw()

    def on_enter_axes(self, event):
        self.in_ax = event.inaxes

    def save_zoom(self):
        self.axlim_dict['image_xlim'] = self.ax0.get_xlim()
        self.axlim_dict['image_ylim'] = self.ax0.get_ylim()
        self.axlim_dict['spec_xlim'] = self.ax1.get_xlim()
        self.axlim_dict['spec_ylim'] = self.ax1.get_ylim()

    def reset_zoom(self):
        self.ax0.set_xlim(self.axlim_dict['image_xlim'])
        self.ax0.set_ylim(self.axlim_dict['image_ylim'])
        self.ax1.set_xlim(self.axlim_dict['spec_xlim'])
        self.ax1.set_ylim(self.axlim_dict['spec_ylim'])


# -------------- An example of these functions (if this script is run, this is the starting point)
if __name__ == '__main__':

    # Declare the fits file location (you can use a standard string too)
    cube_address = Path('D:\Google drive\Astrophysics\Datos\MUSE - Amorin\CGCG007.fits')

    # Load the data
    z_obj = 0.004691
    wave_obj, cube_obj, header_obj = load_muse_cube(cube_address)
    wave_rest = wave_obj / (1 + z_obj)

    # # Declare cube image slices for the analysis
    lineAreas = {'H1_6563A': (6558.0, 6568.0),
                 'S3_6312A': (6308.15, 6317.25),
                 'O3_5007A': (5002.0, 5013.0),
                 'S2_6717A': (6717.0, 6734.0)}
    lineFlux_dict, levelFlux_dict, levelText_dict = compute_muse_flux_images(lineAreas, cube_obj, z_obj,
                                                                            percent_array=PERCENT_ARRAY)

    # Define data for the plot
    line = 'S3_6312A'
    initial_coord = (170, 170)
    background_image = lineFlux_dict['H1_6563A']
    contours_image = lineFlux_dict['S3_6312A']
    fluxLevels = levelFlux_dict['S3_6312A']

    # Labels for the plot overwritting the default ones
    plotConf = {'image': {'xlabel': r'RA',
                          'ylabel': r'DEC',
                          'title': r'Muse slice with H$\alpha$ background'
                                   f'\n{line} foreground' }}

    # Run the plotter
    vp = VoxelPlotter(wave_rest, cube_obj, background_image, initial_coord, image_fg=contours_image,
                           flux_levels=fluxLevels[2:], ax_user_conf=plotConf)
