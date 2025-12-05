from astropy.io import fits
from pathlib import Path
import re
import lime
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astro.stsci.plots import cos_image_plotter
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astro.stsci.tools import generate_apperture_mask, image_gauss_fitting


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]

# Get list of objects
target_list = list(cfg_sample['aper_mean_coord'].keys())

# Loop throught he objects and measure the FWHM
counter = 0
fwhm_dict = {}
for j, object in enumerate(target_list):

    df_obj = sample_df.loc[sample_df.object == object]
    idx_image = sample_df.filepath.str.contains(cfg_sample['favoured_images'][object][0])
    image_path = obs_folder / sample_df.loc[idx_image].filepath.values[0]
    instr_im = sample_df.loc[idx_image].instr.values[0]
    grating = sample_df.loc[idx_image].grating.values[0]

    # Get figure data
    hdr = fits.getheader(image_path, extname='SCI')
    imdata = fits.getdata(image_path, extname='SCI')
    wcs = WCS(hdr)

    groups_dict = cfg_sample['aper_mean_coord'].get(object)
    group_members = list(groups_dict.keys()) if isinstance(groups_dict, dict) else [object]

    for subObject in group_members:

        # Get the object image
        if counter >= 0:

            # Get the pointing aperture
            if not isinstance(groups_dict, dict):
                ra, dec = cfg_sample['aper_mean_coord'][object]
            else:
                ra, dec = cfg_sample['aper_mean_coord'][object][subObject]
            coord_dict = {'Mean apperture': {'ra': ra * u.deg,
                                             'dec': dec * u.deg,
                                             'radius': 1.25 * u.arcsec}}

            # Generate the apperture mask
            mask_points = cfg_sample['apperture_masks'].get(object)
            mask = generate_apperture_mask(imdata.shape, wcs, **coord_dict['Mean apperture'], mask_points=mask_points)
            imdata_masked = imdata * mask

            # Cut to the image region
            section = (2*1.25 * u.arcsec, 2*1.25 * u.arcsec)  # (ny, nx) or single quantity for square
            c_image = SkyCoord(ra=coord_dict['Mean apperture']['ra'],
                               dec=coord_dict['Mean apperture']['dec'],
                               frame="icrs")
            cutout = Cutout2D(imdata_masked, position=c_image, size=section, wcs=WCS(hdr), mode="trim", fill_value=np.nan)

            # Compute the object FWHM
            dec_sum = np.nansum(cutout.data, axis=0)
            comps_dict = image_gauss_fitting(dec_sum)

            # Plot the image
            print(f'{counter}) {object} {instr_im}, {grating}')
            cos_image_plotter(cutout.data, cutout.wcs, object, coord_dict, instr_im='COS', mask_points=mask_points,
                              dec_sum=dec_sum, comps_dict=comps_dict, title=f'{subObject}')

            # Store the results
            fwhm_dict[subObject] = float(np.round(comps_dict['g1_fwhm'], 1))

        counter += 1

# Save measurements
out_folder = project_folder/'LyC_leakers_COS/acq_image_fwhm.toml'
lime.save_cfg(out_folder, fwhm_dict, section_name='LyC_acq_image_fwhm_pixels')

#     # Get the object pointing data
#     idcs_spectra = df_obj.object.index.get_level_values('state') == 'x1dsum'
#     spec_list = df_obj.loc[idcs_spectra, 'filepath'].values
#     coord_dict = {}
#
#     # Sub-identifier for group
#     sub_group = None if object not in cfg_sample['multi_target_labels'] else cfg_sample['multi_target_labels'][object]
#
#     for spec_fname in spec_list:
#         fname = obs_folder/spec_fname
#         spec_hdr = fits.getheader(fname, ext=1)
#         spec_hdr0 = fits.getheader(fname, ext=0)
#         identifier = spec_hdr0['TARGNAME']
#         coord_dict[identifier] =  {'pa': spec_hdr["PA_APER"] * u.deg,
#                                    'orien': spec_hdr["ORIENTAT"] * u.deg,
#                                    'disp': spec_hdr["DISPAXIS"],
#                                    'ra': spec_hdr["RA_APER"] * u.deg,
#                                    'dec': spec_hdr["DEC_APER"] * u.deg,
#                                    'radius': float(re.findall(r"[-+]?\d*\.\d+|\d+", spec_hdr["S_REGION"])[2]) * u.deg}
#
#     # Get figure data
#     hdr = fits.getheader(image_path, extname='SCI')
#     imdata = fits.getdata(image_path, extname='SCI')
#
#     # Plot the image
#     # cos_image_plotter(imdata, WCS(hdr), object, coord_dict, subgroup=sub_group)
#
#     # Cut to the image region
#     section = (1.25 * u.arcsec, 1.25 * u.arcsec)  # (ny, nx) or single quantity for square
#     c_image = SkyCoord(ra=hdr['RA_APER']*u.deg, dec=hdr['DEC_APER']*u.deg, frame="icrs")
#     cutout = Cutout2D(imdata, position=c_image, size=section, wcs=WCS(hdr), mode="trim", fill_value=np.nan)
#
#     # Generate pixel sum and gaussian fitting
#     dec_sum = np.nansum(cutout.data, axis=0)
#     comps_dict = image_gauss_fitting(dec_sum)
#
#     # Plot the images
#     cos_image_plotter(cutout.data, WCS(hdr), object, coord_dict, dec_sum=dec_sum, comps_dict=comps_dict)
#
#     # Store the results
#     fwhm_dict[object] = float(np.round(comps_dict['g1_fwhm'], 1))
#
# # Save measurements
# out_folder = project_folder/'LyC_leakers_COS/acq_image_fwhm.toml'
# lime.save_cfg(out_folder, fwhm_dict, section_name='LyC_acq_image_fwhm_pixels')

