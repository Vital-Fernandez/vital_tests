import numpy as np
import lime
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.wcs import WCS
from pathlib import Path

# Input files
fname = '/home/vital/Astrodata/J0823_p2806_v2.fits'
cfgname = './example_cfg.toml'
combine_mask_file = './J0823_p2806_SN_mask.fits'
output_log = './J0823_p2806_Temp_O3_99p_log.fits'

# Load the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726)
cube.check.cube('O2_3726A', rest_frame=True, fname=output_log)

# Export the line measurements as spatial maps:
# param_list = ['z_line', 'peak_wave', 'peak_flux', 'intg_flux', 'v_r', 'v_r_err', 'chisqr', 'redchi']
# param_list = ['intg_flux']
# lines_list = ['H1_4861A', 'O3_4959A', 'O3_5007A', 'O2_3726A', 'O2_3729A']
# output_log = './J0823_p2806_FG_O2_60p_single_comp_log.fits'
# lime.save_parameter_maps(output_log, './maps', param_list, lines_list, mask_file=combine_mask_file,
#                          mask_list=['FG_O2_60p'], output_file_prefix='J0823_')


# # Export the line measurements as spatial maps:
# param_list = ['z_line', 'peak_wave', 'peak_flux', 'intg_flux', 'profile_flux', 'v_r', 'v_r_err', 'chisqr', 'redchi']
# # param_list = ['intg_flux']
# lines_list = ['H1_4861A', 'O3_4959A', 'O3_5007A', 'O2_3726A', 'O2_3729A']
# output_log = './J0823_p2806_Temp_O3_99p_log.fits'
# lime.save_parameter_maps(output_log, './maps', param_list, lines_list, mask_file=combine_mask_file,
#                          mask_list=['Temp_O3_99p'], output_file_prefix='J0823_TO3_99p_')






# Loop through the line ratios
lines_list = ['H1_4861A', 'O3_5007A']
fits_file = './maps/J0823_intg_flux.fits'
for line in lines_list:

    # Get the astronomical coordinates from one of the headers of the lines log
    hdr = fits.getheader(fits_file, line)
    ratio_map = fits.getdata(fits_file, line)

    wcs_maps = WCS(hdr)

    # Create the plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=wcs_maps, slices=('x', 'y'))
    im = ax.imshow(ratio_map, vmin=np.nanpercentile(ratio_map, 16), vmax=np.nanpercentile(ratio_map, 84))
    cbar = fig.colorbar(im, ax=ax)
    ax.update({'title': f'{line} {Path(fits_file).stem}', 'xlabel': r'RA', 'ylabel': r'DEC'})
    plt.show()


# State line ratios for the plots
fits_file = './maps/J0823_intg_flux.fits'
lines_ratio = {'O3': ['O3_5007A', 'O3_4959A'],
               'O3O2': ['O3_5007A', 'O2_3726A']}

# Loop through the line ratios
for ion, lines in lines_ratio.items():

    # Recover the parameter measurements
    latex_array = lime.label_decomposition(lines, params_list=['latex_label'])[0]
    ratio_map = fits.getdata(fits_file, lines[0]) / fits.getdata(fits_file, lines[1])
    Halpha = fits.getdata(fits_file, lines[0])
    Hbeta = fits.getdata(fits_file, lines[1])

    # Get the astronomical coordinates from one of the headers of the lines log
    hdr = fits.getheader(fits_file, lines[0])
    wcs_maps = WCS(hdr)

    # Create the plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=wcs_maps, slices=('x', 'y'))
    im = ax.imshow(ratio_map, vmin=np.nanpercentile(ratio_map, 16), vmax=np.nanpercentile(ratio_map, 84))
    cbar = fig.colorbar(im, ax=ax)
    ax.update({'title': f'J0823 flux ratio: {latex_array[0]} / {latex_array[1]}', 'xlabel': r'RA', 'ylabel': r'DEC'})
    plt.show()


# Loop through the line ratios
lines_list = ['H1_4861A', 'O3_5007A']
fits_file = './maps/J0823_redchi.fits'
for line in lines_list:

    # Get the astronomical coordinates from one of the headers of the lines log
    hdr = fits.getheader(fits_file, line)
    ratio_map = fits.getdata(fits_file, line)

    wcs_maps = WCS(hdr)

    # Create the plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=wcs_maps, slices=('x', 'y'))
    im = ax.imshow(ratio_map, vmin=np.nanpercentile(ratio_map, 16), vmax=np.nanpercentile(ratio_map, 84))
    cbar = fig.colorbar(im, ax=ax)
    ax.update({'title': f'{line} {Path(fits_file).stem}', 'xlabel': r'RA', 'ylabel': r'DEC'})
    plt.show()


# Loop through the line ratios
lines_list = ['H1_4861A', 'O3_5007A']
fits_file = './maps/J0823_peak_wave.fits'
for line in lines_list:

    # Get the astronomical coordinates from one of the headers of the lines log
    hdr = fits.getheader(fits_file, line)
    ratio_map = fits.getdata(fits_file, line)

    wcs_maps = WCS(hdr)

    # Create the plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=wcs_maps, slices=('x', 'y'))
    im = ax.imshow(ratio_map, vmin=np.nanpercentile(ratio_map, 16), vmax=np.nanpercentile(ratio_map, 84))
    cbar = fig.colorbar(im, ax=ax)
    ax.update({'title': f'{line} {Path(fits_file).stem}', 'xlabel': r'RA', 'ylabel': r'DEC'})
    plt.show()