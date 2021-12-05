import numpy as np
from astropy.wcs import WCS
from matplotlib import pyplot as plt, rcParams, gridspec, cm, colors
from src.specsiser.tools.line_measure import STANDARD_PLOT, STANDARD_AXES
from astropy.io import fits

STANDARD_PLOT = {'figure.figsize': (20, 14), 'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 12,
                 'xtick.labelsize': 12, 'ytick.labelsize': 12}

background_color = np.array((43, 43, 43))/255.0
foreground_color = np.array((179, 199, 216))/255.0

DARK_PLOT = {'figure.figsize': (14, 7),
             'axes.titlesize': 14,
             'axes.labelsize': 14,
             'legend.fontsize': 12,
             'xtick.labelsize': 12,
             'ytick.labelsize': 12,
             'text.color': foreground_color,
             'figure.facecolor': background_color,
             'axes.facecolor': background_color,
             'axes.edgecolor': foreground_color,
             'axes.labelcolor': foreground_color,
             'xtick.color': foreground_color,
             'ytick.color': foreground_color,
             'legend.edgecolor': 'inherit',
             'legend.facecolor': 'inherit'}

line_regions = {'H1_6563A': np.array([6850, 6910]),
                'S2_6731A_b': np.array([7027.5, 7057.5]),
                'O3_4363A': np.array([4565.0, 4575.0]),
                'S3_6312A': np.array([6606.5, 6617.0]),
                'O3_5007A': np.array([5232.0, 5260.0]),
                'S3_9069A': np.array([9492.5, 9506.5]),
                'S3_9531A': np.array([9975.5, 9995.0])}

megaradrp_modes = {'bias'       : 'MegaraBiasImage',
                   'arc'        : 'MegaraArcCalibration',
                   'trace_map'  : 'MegaraTraceMap',
                   'slit_flat'  : 'MegaraSlitFlat',
                   'model_map'  : 'MegaraModelMap',
                   'fiber_flat' : 'MegaraFiberFlatImage',
                   'lcb_acq'    : 'MegaraLcbAcquisition',
                   'lcb_std'    : 'MegaraLcbStdStar',
                   'lcb_image'  : 'MegaraLcbImage'}

def open_manga_cubes(file_address):

    with fits.open(file_address) as hdul:
        hdr = hdul['FLUX'].header
        wave = hdul['WAVE'].data
        flux = hdul['FLUX'].data
        err = 1/np.sqrt(hdul['IVAR'].data)
        pixMask = hdul['MASK'].data

        return wave, flux, err, hdr

class VoxelPlotter(object):

    """
    This class produces an interative matplotlib window for the muse data cubes. On the left axis with the cube slice
    image you can right click a voxel for its corresponding spectrum to be plotted on the right axis.
    """

    def __init__(self, obj_wave, obj_cube, image_bg, voxel_coord=None, image_fg=None, flux_levels=None,
                fig_user_conf={}, ax_user_conf={}, header=None):

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
        self.header = header

        # Plot Configuration
        defaultConf = STANDARD_PLOT.copy()
        defaultConf.update(fig_user_conf)
        rcParams.update(defaultConf)

        # Figure structure
        self.fig = plt.figure(figsize=(18, 5))
        gs = gridspec.GridSpec(nrows=1, ncols=2, figure=self.fig, width_ratios=[1, 2], height_ratios=[1])
        cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        aee = self.fig.canvas.mpl_connect('axes_enter_event', self.on_enter_axes)

        # Axes configuration
        sky_wcs = WCS(self.header)
        self.ax0 = self.fig.add_subplot(gs[0], projection=sky_wcs, slices=('x', 'y', 1))
        self.ax1 = self.fig.add_subplot(gs[1])

        imgAxConf = {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Cube flux slice'}
        if 'image' in ax_user_conf:
            imgAxConf.update(ax_user_conf['image'])
        self.axConf['image'].update(imgAxConf)

        specAxConf = STANDARD_AXES.copy()
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
        # cmap = plt.cm.gray
        # min_flux = np.nanpercentile(image_bg, 30)
        #
        # norm = plt.Normalize(min_flux, image_bg.max())
        # image_bg_rgb = cmap(norm(image_bg))

        # # Plot backgrou nd image
        # self.ax0.imshow(image_bg_rgb, cmap='gray', vmin=0.0, aspect='auto')

        min_flux = np.nanpercentile(image_bg, 60)
        norm_color_bg = colors.SymLogNorm(linthresh=min_flux,
                                          vmin=min_flux,
                                          base=10)
        self.ax0.imshow(image_bg, cmap=cm.gray, norm=norm_color_bg)

        # Emphasize input coordinate
        if voxel_coord is not None:
            idx_j, idx_i = voxel_coord
            # image_bg_rgb[idx_j, idx_i, :3] = [1, 0, 0]
            self.ax0.plot(idx_i, idx_j, '+', color='red')

        # Plot contours image
        if image_fg is not None:
            CS3 = self.ax0.contour(self.grid_mesh[0], self.grid_mesh[1], image_fg, cmap='viridis', levels=flux_levels,
                                   norm=colors.LogNorm())

        # Voxel spectrum
        if voxel_coord is not None:
            idx_j, idx_i = voxel_coord
            flux_voxel = self.cube_data[:, idx_j, idx_i]
            self.ax1.step(self.wave, flux_voxel)

        conf_dict['spectrum']['title'] = f'Voxel {idx_j} - {idx_i}'

        # Update the axis
        self.ax0.update(conf_dict['image'])
        self.ax1.update(conf_dict['spectrum'])

        return

    def on_click(self, event, mouse_trigger_buttton=3):

        """
        This method defines launches the new plot selection once the user clicks on an image voxel. By default this is a
        a right click on a minimum three button mouse
        :param event: This variable represents the user action on the plot
        :param mouse_trigger_buttton: Number-coded mouse button which defines the button launching the voxel selection
        :return:
        """

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
                self.plot_map_voxel(self.image_bg, (idx_j, idx_i), self.image_fg, self.flux_levels, conf_dict=self.axConf)

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