import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube
from astropy.wcs import WCS
from matplotlib import pyplot as plt, rcParams, gridspec
from src.specsiser.tools.line_measure import STANDARD_PLOT, STANDARD_AXES
import pyneb as pn

lineAreas = {'H1_6563A': (6558.0, 6568.0),
             'S3_6312A': (6308.15, 6317.25),
             'O3_5007A': (5002.0, 5013.0),
             'S2_6717A': (6717.0, 6734.0)}


def store_frame_to_fits(fits_address, fits_hdu, ext_name):

    if fits_address.is_file():
        try:
            fits.update(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name, verify=True)
        except KeyError:
            fits.append(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name)
    else:
        fits_hdu.writeto(fits_address, overwrite=True, output_verify='fix')

    return


def import_fits_data(file_address, instrument, frame_idx=None):

    cube = Cube(filename=str(file_address))
    header = cube.data_header

    cube.wave.info()
    dw = header['CD3_3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, cube, header


# STANDARD_IMAGE = {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Cube flux slice'}
#
# default_linelog_types = {'index': '<U50',
#                  'wavelength': '<f8',
#                  'intg_flux': '<f8',
#                  'intg_err': '<f8',
#                  'gauss_flux': '<f8',
#                  'gauss_err': '<f8',
#                  'eqw': '<f8',
#                  'eqw_err': '<f8',
#                  'ion': '<U50',
#                  'pynebCode': '<f8',
#                  'pynebLabel': '<f8',
#                  'lineType': '<f8',
#                  'latexLabel': '<U50',
#                  'blended_label': '<U50',
#                  'w1': '<f8',
#                  'w2': '<f8',
#                  'w3': '<f8',
#                  'w4': '<f8',
#                  'w5': '<f8',
#                  'w6': '<f8',
#                  'm_cont': '<f8',
#                  'n_cont': '<f8',
#                  'cont': '<f8',
#                  'std_cont': '<f8',
#                  'peak_flux': '<f8',
#                  'peak_wave': '<f8',
#                  'snr_line': '<f8',
#                  'snr_cont': '<f8',
#                  'amp': '<f8',
#                  'mu': '<f8',
#                  'sigma': '<f8',
#                  'amp_err': '<f8',
#                  'mu_err': '<f8',
#                  'sigma_err': '<f8',
#                  'v_r': '<f8',
#                  'v_r_err': '<f8',
#                  'sigma_vel': '<f8',
#                  'sigma_err_vel': '<f8',
#                  'observation': '<U50',
#                  'comments': '<U50',
#                  'obsFlux': '<f8',
#                  'obsFluxErr': '<f8',
#                  'f_lambda': '<f8',
#                  'obsInt': '<f8',
#                  'obsIntErr': '<f8'}
#
# background_color = np.array((43, 43, 43))/255.0
# foreground_color = np.array((179, 199, 216))/255.0
# red_color = np.array((43, 43, 43))/255.0
# yellow_color = np.array((191, 144, 0))/255.0
#
# DARK_PLOT = {'figure.figsize': (14, 7),
#              'axes.titlesize': 14,
#              'axes.labelsize': 14,
#              'legend.fontsize': 12,
#              'xtick.labelsize': 12,
#              'ytick.labelsize': 12,
#              'text.color': foreground_color,
#              'figure.facecolor': background_color,
#              'axes.facecolor': background_color,
#              'axes.edgecolor': foreground_color,
#              'axes.labelcolor': foreground_color,
#              'xtick.color': foreground_color,
#              'ytick.color': foreground_color,
#              'legend.edgecolor': 'inherit',
#              'legend.facecolor': 'inherit'}
#
# label_Conver = {'H1_6563A': 'Halpha',
#                'H1_9229A': 'HPas9',
#                'H1_9015A': 'HPas10',
#                'H1_8863A': 'HPas11',
#                'H1_8750A': 'HPas12'}
#
#
# latex_Conver = {'H1_6563A': r'H\alpha',
#                 'H1_4861A': r'H\beta',
#                 'H1_9229A': r'H_{Pas,\,9}',
#                 'H1_9015A': r'H_{Pas,\,10}',
#                 'H1_8863A': r'H_{Pas,\,11}',
#                 'H1_8750A': r'H_{Pas,\,12}',
#                 'v_r': r'$v_{r}\,(km/s)$',
#                 'sigma_vel': r'$\sigma_{int}\,(km/s)$'}
#
# dinamicLines = {'H1_6563A': r'$H\alpha_{Narrow}$',
#               'H1_6563A_w1': r'$H\alpha_{Broad\,1}$',
#               'H1_6563A_w2': r'$H\alpha_{Broad\,2}$',
#               'H1_4861A': r'$H\beta_{Narrow}$',
#               'H1_4861A_w1': r'$H\beta_{Broad\,1}$',
#               'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
#               'O3_5007A_w1': r'$[OIII]5007\AA_{Broad\,1}$',
#               'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
#               'O3_5007A_w1': r'$[OIII]5007\AA_{Broad\,1}$',
#               'S3_9069A': r'$[SIII]9069\AA_{Narrow}$',
#               'He1_5876A': r'$HeI\,5876\AA_{Narrow}$'}
#
grid_columns = {'logZ': 'logOH',
                'logU': 'logU',
                'logNO': 'logNO',
                'carbon': 'carbon',
                'o3726': 'O2_3726A',
                'o3729': 'O2_3729A',
                'ne3869': 'Ne3_3869A',
                'ne3968': 'Ne3_3968A',
                'h3970': 'H1_3970A',
                's4070': 'S2_4069A',
                's4078': 'S2_4078A',
                'h4102': 'H1_4102A',
                'c4267': 'C2_4267A',
                'h4341': 'H1_4341A',
                'o4363': 'O3_4363A',
                'he4471': 'He1_4471A',
                'o4651': 'O1_4651A',
                'c4659': 'C2_4659A',
                'fe4668': 'Fe3_4668A',
                'he4686': 'He2_4686A',
                'ar4711': 'Ar4_4711A',
                'ar4740': 'Ar4_4740A',
                'h4861': 'H1_4861A',
                'o4959': 'O3_4959A',
                'o5007': 'O3_5007A',
                'ar5192': 'Ar3_5192A',
                'n5198': 'N1_5198A',
                'n5200': 'N1_5200A',
                'cl5518': 'Cl3_5518A',
                'cl5538': 'Cl3_5538A',
                'n5755': 'N2_5755A',
                'he5876': 'He1_5876A',
                'o6300': 'O1_6300A',
                's6312': 'S3_6312A',
                'n6548': 'N2_6548A',
                'h6563': 'H1_6563A',
                'n6584': 'N2_6584A',
                'he6678': 'He1_6678A',
                's6716': 'S2_6716A',
                's6731': 'S2_6731A',
                'he7065': 'He1_7065A',
                'ar7135': 'Ar3_7136A',
                'o7323': 'O2_7319A',
                'o7332': 'O2_7330A',
                'ar7751': 'Ar3_7751A',
                's9069': 'S3_9069A',
                's9532': 'S3_9531A'}

# grid_HII_CHI_mistry_conversion = {'logOH': '12+log(O/H)',
#                                 'logNO': 'log(N/O)',
#                                 'logU': 'log(U)',
#                                 'O2_3726A_b': 'OII_3727',
#                                 'Ne3_3869A': 'NeIII_3868',
#                                 'O3_4363A': 'OIII_4363',
#                                 'O3_4959A': 'OIII_4959',
#                                 'O3_5007A': 'OIII_5007',
#                                 'N2_6584A': 'NII_6584',
#                                 'S2_6716A_b': 'SII_6725',
#                                 'S2_6716A': 'SII_6716',
#                                 'S2_6731A': 'SII_6731'}
#
#
# def reconstruct_wavelength(header):
#     dw = header['CDELT3']
#     w_min = header['CRVAL3']
#     nPixels = header['NAXIS3']
#     w_max = w_min + dw * nPixels
#     return np.linspace(w_min, w_max, nPixels, endpoint=False)
#
# def red_corr_HalphaHbeta_ratio(lines_df, default_cHbeta):
#
#     # Normalizing flux
#     if ('H1_6563A' in lines_df.index) and ('H1_4861A' in lines_df.index):
#         flux_Halpha = lines_df.loc['H1_6563A', 'intg_flux']
#         flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
#         halpha_norm = flux_Halpha / flux_Hbeta
#
#         rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
#         rc.setCorr(obs_over_theo=halpha_norm / 2.86, wave1=6563., wave2=4861.)
#         cHbeta = rc.cHbeta
#     else:
#         rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=default_cHbeta)
#         cHbeta = 'none'
#
#     return cHbeta, rc
#
#
# def image_array_binning(flux_image, percentil_scale, bad_value=np.nan, below_value=None, above_value=None):
#
#     '''
#     This function returns the bin index for an 2D array according to a percentil_scale which is supposed to be an ordered array.
#     The user can introduce a default value for nan entries in the 2D array to be assigned in the output index array
#     The user can assign a defaul value to be considered for entries in the 2D array which lie outside the bin scale
#     scales. If not specified these entries are assigned to the lower and higher bins respectively
#     :param percentil_scale:
#     :param flux_image:
#     :param bins:
#     :param bad_value:
#     :param below_value:
#     :param above_value:
#     :return:
#     '''
#
#     bins = np.percentile(flux_image, percentil_scale)
#
#     binned_array = np.digitize(flux_image, bins, right=False) - 1.0
#
#     # Bad entries assignment
#     idcs_bad = np.isnan(flux_image)
#     if idcs_bad.any():
#         bad_value = np.nan if bad_value is None else bad_value
#         if np.isnan(bad_value):
#             bad_bin = np.nan
#         else:
#             bad_bin = np.argmin(np.logical_xor(bad_value > bins, bins[-1] < bins[0]))
#         binned_array[idcs_bad] = bad_bin
#
#     # Values below binning
#     below_value = bins[0] if below_value is None else below_value
#     idcs_below = flux_image < below_value
#     if idcs_below.any():
#         below_bin = np.argmin(np.logical_xor(below_value > bins, bins[-1] < bins[0]))
#         binned_array[idcs_below] = below_bin
#
#     # Values above binning
#     above_value = bins[-1] if above_value is None else above_value
#     idcs_above = flux_image > above_value
#     if idcs_above.any():
#         above_bin = np.argmin(np.logical_xor(above_value > bins, bins[-1] < bins[0]))
#         binned_array[idcs_above] = above_bin
#
#     return binned_array
#
#
# def compute_line_flux_image(lineArea_dict, flux_cube, redshift, percent_array):
#
#     # Containers for output data
#     lineFlux_dict = {}
#     levelFlux_dict = {}
#     levelText_dict = {}
#
#     # Loop through the input wavelength regions
#     for lineLabel, lineLimits in lineArea_dict.items():
#
#         # Extract cube slice using mpdaf defult tools. This requires the input wavelengths to be
#         # on the same scale as in the cube
#         line_image = flux_cube.get_image(np.array(lineLimits) * (1 + redshift), subtract_off=True)
#         flux_image = line_image.data.data
#
#         # Personal scale for the image flux
#         log_flux = np.log10(flux_image)
#         flux_contours_i = np.zeros(flux_image.shape)
#         idcs_removed = np.logical_or(log_flux < 0.0, np.isnan(log_flux))
#         flux_contours_i[~idcs_removed] = log_flux[~idcs_removed]
#
#         # Define image countours based on the flux percentiles
#         levelFlux_i = np.percentile(flux_contours_i[flux_contours_i>0], percent_array)
#         levels_text_i = ['None'] * len(levelFlux_i)
#         for idx, per in enumerate(percent_array):
#             levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$'
#
#         # Store the data
#         lineFlux_dict[lineLabel] = flux_contours_i
#         levelFlux_dict[lineLabel] = levelFlux_i
#         levelText_dict[lineLabel] = levels_text_i
#
#     return lineFlux_dict, levelFlux_dict, levelText_dict
#
#
# def plot_voxel_flux(image_bg, voxel_coord, wave_voxel, flux_voxel, image_fg=None, flux_levels=None, sky_wcs=None,
#                     image_conf={}, spec_conf={}):
#
#     # Plot Configuration
#     defaultConf = STANDARD_PLOT.copy()
#     rcParams.update(defaultConf)
#
#     # Selecting plotting value pixels
#     frame_size = image_bg.shape
#     x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
#     X, Y = np.meshgrid(x, y)
#
#     fig = plt.figure(figsize=(18, 5))
#     gs = gridspec.GridSpec(nrows=1, ncols=2, figure=fig, width_ratios=[1, 2], height_ratios=[1])
#     ax0 = fig.add_subplot(gs[0], projection=sky_wcs, slices=('x', 'y', 1))
#     ax1 = fig.add_subplot(gs[1])
#
#     # Cube image
#     cmap = plt.cm.gray
#     norm = plt.Normalize(image_bg.min(), image_bg.max())
#     image_bg_rgb = cmap(norm(image_bg))
#
#     # Mark the pixels
#     idx_j, idx_i = voxel_coord
#     image_bg_rgb[idx_j, idx_i, :3] = [1, 0, 0]
#     ax0.plot(idx_i, idx_j, '+', color='red')
#
#     ax0.imshow(image_bg_rgb,  cmap='gray', vmin=0.0, aspect='auto')
#
#     if image_fg is not None:
#         CS3 = ax0.contour(X, Y, image_fg, levels=flux_levels, alpha=0.5)
#
#     # Voxel spectrum
#     ax1.step(wave_voxel, flux_voxel)
#
#     # Plot format
#     image_conf_default = {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Cube slice'}
#     spec_conf_default = STANDARD_AXES.copy()
#     spec_conf_default.update({'title': f'Voxel {idx_j} - {idx_i} slice'})
#
#     image_conf_default.update(image_conf)
#     spec_conf_default.update(spec_conf)
#
#     ax0.update(image_conf_default)
#     ax1.update(spec_conf_default)
#     plt.show()
#
#     # plotAddress = dataFolder / fileList[i].replace('.fits', '_armFluxComparison.png')
#     # plt.savefig(plotAddress, dpi=200, bbox_inches='tight')
#
#     return
#
# def voxel_security_check(linesDF):
#
#     check = False
#
#     if 'H1_4861A' in linesDF.index:
#         if 'S3_6312A' in linesDF.index:
#             check = True
#
#     return check
#
# def fits_db(fits_address, model_db, ext_name):
#
#     total_params_list = np.array(list(model_db['Fitting_results'].keys()))
#     inputLabels = model_db['Input_data']['lineLabels_list']
#     inputFlux = model_db['Input_data']['inputFlux_array']
#     inputErr = model_db['Input_data']['inputErr_array']
#
#     # Save variable traces
#     list_columns = []
#     for i, param in enumerate(model_db['trace'].varnames):
#         trace = model_db['trace'][param]
#         if param in total_params_list:
#             col_param = fits.Column(name=param, format='E', array=trace)
#             list_columns.append(col_param)
#         else:
#             if ('_Op' not in param) and ('_log__' not in param):
#                 col_param = fits.Column(name=param, format='E', array=trace)
#                 list_columns.append(col_param)
#
#     # Save flux traces
#     if 'calcFluxes_Op' in model_db['trace'].varnames:
#         trace = model_db['trace']['calcFluxes_Op']
#         for i in range(trace.shape[1]):
#             flux_trace = trace[:, i]
#             col_param = fits.Column(name=inputLabels[i], format='E', array=flux_trace)
#             list_columns.append(col_param)
#
#     cols = fits.ColDefs(list_columns)
#     hdu = fits.BinTableHDU.from_columns(cols, name=ext_name)
#
#     for i, label in enumerate(inputLabels):
#         hdu.header[f'hierarch flux_{label}'] = inputFlux[i]
#         hdu.header[f'hierarch err_{label}'] = inputErr[i]
#
#     if fits_address.is_file():
#         try:
#             fits.update(fits_address, data=hdu.data, header=hdu.header, extname=ext_name, verify=True)
#         except KeyError:
#             fits.append(fits_address, data=hdu.data, header=hdu.header, extname=ext_name)
#     else:
#         hdu.writeto(fits_address, overwrite=True, output_verify='fix')
#
#     return
#
# def store_frame_to_fits(fits_address, fits_hdu, ext_name):
#
#     if fits_address.is_file():
#         try:
#             fits.update(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name, verify=True)
#         except KeyError:
#             fits.append(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name)
#     else:
#         fits_hdu.writeto(fits_address, overwrite=True, output_verify='fix')
#
#     return
#
#
# class IFU_Cube_Plotter(object):
#
#     """
#     This class produces an interative matplotlib window for the muse data cubes. On the left axis with the cube slice
#     image you can right click a voxel for its corresponding spectrum to be plotted on the right axis.
#     """
#
#     def __init__(self, obj_wave, obj_data, image_bg, voxel_coord=None, image_fg=None, flux_levels=None,
#                 fig_user_conf={}, ax_user_conf={}, header=None):
#
#         self.fig = None
#         self.ax0, self.ax1, self.in_ax = None, None, None
#         self.axlim_dict = {}
#         self.grid_mesh = None
#         self.cube_data = obj_data
#         self.wave = obj_wave
#         self.image_bg = image_bg
#         self.image_fg = image_fg
#         self.flux_levels = flux_levels
#         self.axConf = dict(image={}, spectrum={})
#
#         # Plot Configuration
#         defaultConf = STANDARD_PLOT.copy()
#         defaultConf.update(fig_user_conf)
#         rcParams.update(defaultConf)
#
#         # Figure structure
#         self.fig = plt.figure(figsize=(18, 5))
#         gs = gridspec.GridSpec(nrows=1, ncols=2, figure=self.fig, width_ratios=[1, 2], height_ratios=[1])
#         cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
#         aee = self.fig.canvas.mpl_connect('axes_enter_event', self.on_enter_axes)
#
#         # Axes configuration
#         if header is not None:
#             sky_wcs = WCS(header)
#             self.ax0 = self.fig.add_subplot(gs[0], projection=sky_wcs, slices=('x', 'y', 1))
#         else:
#             self.ax0 = self.fig.add_subplot()
#         self.ax1 = self.fig.add_subplot(gs[1])
#
#         imgAxConf = STANDARD_IMAGE.copy()
#         if 'image' in ax_user_conf:
#             imgAxConf.update(ax_user_conf['image'])
#         self.axConf['image'].update(imgAxConf)
#
#         specAxConf = STANDARD_AXES.copy()
#         if 'spectrum' in ax_user_conf:
#             specAxConf.update(ax_user_conf['spectrum'])
#         self.axConf['spectrum'].update(specAxConf)
#
#         # Image mesh grid
#         frame_size = image_bg.shape
#         y, x = np.arange(0, frame_size[0]), np.arange(0, frame_size[1])
#         self.grid_mesh = np.meshgrid(x, y)
#
#         # Generate the plot
#         self.plot_map_voxel(self.image_bg, voxel_coord, self.image_fg, self.flux_levels, conf_dict=self.axConf)
#         plt.show()
#
#         return
#
#     def plot_map_voxel(self, image_bg, voxel_coord=None, image_fg=None, flux_levels=None, conf_dict={}):
#
#         # Image color format
#         cmap = plt.cm.gray
#         norm = plt.Normalize(image_bg.min(), image_bg.max())
#         image_bg_rgb = cmap(norm(image_bg))
#
#         # Plot background image
#         self.ax0.imshow(image_bg_rgb, cmap='gray', vmin=0.0, aspect='auto')
#
#         # Emphasize input coordinate
#         if voxel_coord is not None:
#             idx_j, idx_i = voxel_coord
#             image_bg_rgb[idx_j, idx_i, :3] = [1, 0, 0]
#             self.ax0.plot(idx_i, idx_j, '+', color='red')
#         else:
#             idx_j, idx_i = int(image_bg_rgb.shape[0]/2), int(image_bg_rgb.shape[1]/2)
#
#         # Plot contours image
#         if image_fg is not None:
#             CS3 = self.ax0.contour(self.grid_mesh[0], self.grid_mesh[1], image_fg, levels=flux_levels, alpha=0.5)
#
#         # Voxel spectrum
#         flux_voxel = self.cube_data[:, idx_j, idx_i]
#         self.ax1.step(self.wave, flux_voxel)
#
#         conf_dict['spectrum']['title'] = f'Voxel {idx_j} - {idx_i}'
#
#         # Update the axis
#         self.ax0.update(conf_dict['image'])
#         self.ax1.update(conf_dict['spectrum'])
#         # self.ax1.set_yscale('log')
#         return
#
#     def on_click(self, event, mouse_trigger_buttton=3):
#
#         """
#         This method defines launches the new plot selection once the user clicks on an image voxel. By default this is a
#         a right click on a minimum three button mouse
#         :param event: This variable represents the user action on the plot
#         :param mouse_trigger_buttton: Number-coded mouse button which defines the button launching the voxel selection
#         :return:
#         """
#
#         if self.in_ax == self.ax0:
#
#             if event.button == mouse_trigger_buttton:
#
#                 # Save axes zoom
#                 self.save_zoom()
#
#                 # Save clicked coordinates for next plot
#                 idx_j, idx_i = np.rint(event.ydata).astype(int), np.rint(event.xdata).astype(int)
#                 print(f'Current voxel: {idx_j}-{idx_i} (mouse button {event.button})')
#
#                 # Remake the drawing
#                 self.ax0.clear()
#                 self.ax1.clear()
#                 self.plot_map_voxel(self.image_bg, (idx_j, idx_i), self.image_fg, self.flux_levels, conf_dict=self.axConf)
#
#                 # Reset the image
#                 self.reset_zoom()
#                 self.fig.canvas.draw()
#
#     def on_enter_axes(self, event):
#         self.in_ax = event.inaxes
#
#     def save_zoom(self):
#         self.axlim_dict['image_xlim'] = self.ax0.get_xlim()
#         self.axlim_dict['image_ylim'] = self.ax0.get_ylim()
#         self.axlim_dict['spec_xlim'] = self.ax1.get_xlim()
#         self.axlim_dict['spec_ylim'] = self.ax1.get_ylim()
#
#     def reset_zoom(self):
#         self.ax0.set_xlim(self.axlim_dict['image_xlim'])
#         self.ax0.set_ylim(self.axlim_dict['image_ylim'])
#         self.ax1.set_xlim(self.axlim_dict['spec_xlim'])
#         self.ax1.set_ylim(self.axlim_dict['spec_ylim'])
#
#
# class VoxelPlotter(object):
#
#     """
#     This class produces an interative matplotlib window for the muse data cubes. On the left axis with the cube slice
#     image you can right click a voxel for its corresponding spectrum to be plotted on the right axis.
#     """
#
#     def __init__(self, obj_wave, obj_cube, image_bg, voxel_coord=None, image_fg=None, flux_levels=None,
#                 fig_user_conf={}, ax_user_conf={}, header=None):
#
#         self.fig = None
#         self.ax0, self.ax1, self.in_ax = None, None, None
#         self.axlim_dict = {}
#         self.grid_mesh = None
#         self.cube_data = obj_cube
#         self.wave = obj_wave
#         self.image_bg = image_bg
#         self.image_fg = image_fg
#         self.flux_levels = flux_levels
#         self.axConf = dict(image={}, spectrum={})
#         self.header = header
#
#         # Plot Configuration
#         defaultConf = STANDARD_PLOT.copy()
#         defaultConf.update(fig_user_conf)
#         rcParams.update(defaultConf)
#
#         # Figure structure
#         self.fig = plt.figure(figsize=(18, 5))
#         gs = gridspec.GridSpec(nrows=1, ncols=2, figure=self.fig, width_ratios=[1, 2], height_ratios=[1])
#         cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
#         aee = self.fig.canvas.mpl_connect('axes_enter_event', self.on_enter_axes)
#
#         # Axes configuration
#         sky_wcs = WCS(self.header)
#         self.ax0 = self.fig.add_subplot(gs[0], projection=sky_wcs, slices=('x', 'y', 1))
#         self.ax1 = self.fig.add_subplot(gs[1])
#
#         imgAxConf = STANDARD_IMAGE.copy()
#         if 'image' in ax_user_conf:
#             imgAxConf.update(ax_user_conf['image'])
#         self.axConf['image'].update(imgAxConf)
#
#         specAxConf = STANDARD_AXES.copy()
#         if 'spectrum' in ax_user_conf:
#             specAxConf.update(ax_user_conf['spectrum'])
#         self.axConf['spectrum'].update(specAxConf)
#
#         # Image mesh grid
#         frame_size = self.cube_data.shape
#         y, x = np.arange(0, frame_size[1]), np.arange(0, frame_size[2])
#         self.grid_mesh = np.meshgrid(x, y)
#
#         # Generate the plot
#         self.plot_map_voxel(self.image_bg, voxel_coord, self.image_fg, self.flux_levels, conf_dict=self.axConf)
#         plt.show()
#
#         return
#
#     def plot_map_voxel(self, image_bg, voxel_coord=None, image_fg=None, flux_levels=None, conf_dict={}):
#
#         # Image color format
#         cmap = plt.cm.gray
#         norm = plt.Normalize(image_bg.min(), image_bg.max())
#         image_bg_rgb = cmap(norm(image_bg))
#
#         # Plot background image
#         self.ax0.imshow(image_bg_rgb, cmap='gray', vmin=0.0, aspect='auto')
#
#         # Emphasize input coordinate
#         if voxel_coord is not None:
#             idx_j, idx_i = voxel_coord
#             image_bg_rgb[idx_j, idx_i, :3] = [1, 0, 0]
#             self.ax0.plot(idx_i, idx_j, '+', color='red')
#
#         # Plot contours image
#         if image_fg is not None:
#             CS3 = self.ax0.contour(self.grid_mesh[0], self.grid_mesh[1], image_fg, levels=flux_levels, alpha=0.5)
#
#         # Voxel spectrum
#         if voxel_coord is not None:
#             idx_j, idx_i = voxel_coord
#             flux_voxel = self.cube_data[:, idx_j, idx_i]
#             self.ax1.step(self.wave, flux_voxel)
#
#         conf_dict['spectrum']['title'] = f'Voxel {idx_j} - {idx_i}'
#
#         # Update the axis
#         self.ax0.update(conf_dict['image'])
#         self.ax1.update(conf_dict['spectrum'])
#
#         return
#
#     def on_click(self, event, mouse_trigger_buttton=3):
#
#         """
#         This method defines launches the new plot selection once the user clicks on an image voxel. By default this is a
#         a right click on a minimum three button mouse
#         :param event: This variable represents the user action on the plot
#         :param mouse_trigger_buttton: Number-coded mouse button which defines the button launching the voxel selection
#         :return:
#         """
#
#         if self.in_ax == self.ax0:
#
#             if event.button == mouse_trigger_buttton:
#
#                 # Save axes zoom
#                 self.save_zoom()
#
#                 # Save clicked coordinates for next plot
#                 idx_j, idx_i = np.rint(event.ydata).astype(int), np.rint(event.xdata).astype(int)
#                 print(f'Current voxel: {idx_j}-{idx_i} (mouse button {event.button})')
#
#                 # Remake the drawing
#                 self.ax0.clear()
#                 self.ax1.clear()
#                 self.plot_map_voxel(self.image_bg, (idx_j, idx_i), self.image_fg, self.flux_levels, conf_dict=self.axConf)
#
#                 # Reset the image
#                 self.reset_zoom()
#                 self.fig.canvas.draw()
#
#     def on_enter_axes(self, event):
#         self.in_ax = event.inaxes
#
#     def save_zoom(self):
#         self.axlim_dict['image_xlim'] = self.ax0.get_xlim()
#         self.axlim_dict['image_ylim'] = self.ax0.get_ylim()
#         self.axlim_dict['spec_xlim'] = self.ax1.get_xlim()
#         self.axlim_dict['spec_ylim'] = self.ax1.get_ylim()
#
#     def reset_zoom(self):
#         self.ax0.set_xlim(self.axlim_dict['image_xlim'])
#         self.ax0.set_ylim(self.axlim_dict['image_ylim'])
#         self.ax1.set_xlim(self.axlim_dict['spec_xlim'])
#         self.ax1.set_ylim(self.axlim_dict['spec_ylim'])