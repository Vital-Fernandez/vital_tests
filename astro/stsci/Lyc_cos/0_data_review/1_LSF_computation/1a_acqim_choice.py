from pathlib import Path
import numpy as np
import re
import lime
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astro.stsci.plots import cos_image_plotter
from astropy.coordinates import SkyCoord


# Data location
project_folder = Path('/')
obs_folder = Path('/home/vital/Astrodata/STScI')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

for object in sample_df.object.unique():
    print(object, sample_df.loc[sample_df.object == object].index.get_level_values('id').unique())

object_images = ['IZw18_SE']

for j, object in enumerate(object_images):

    labels = sample_df.loc[sample_df.object == object].index.get_level_values('id').unique()
    df_obj = sample_df.loc[sample_df.object == object]
    print(object, labels)

    # Get the object pointing data
    idcs_spectra = df_obj.object.index.get_level_values('state') == 'x1dsum'
    spec_list = df_obj.loc[idcs_spectra, 'filepath'].values
    coord_dict = {}

    # Sub-identifier for group
    sub_group = None if object not in cfg_sample['multi_target_labels'] else cfg_sample['multi_target_labels'][object]

    for fname in spec_list:
        fname = obs_folder/fname
        spec_hdr = fits.getheader(fname, ext=1)
        spec_hdr0 = fits.getheader(fname, ext=0)
        identifier = spec_hdr0['TARGNAME']#f'{spec_hdr0['TARGNAME']}_{fname.stem}'
        coord_dict[identifier] =  {'pa': spec_hdr["PA_APER"] * u.deg,
                                   'orien': spec_hdr["ORIENTAT"] * u.deg,
                                   'disp': spec_hdr["DISPAXIS"],
                                   'ra': spec_hdr["RA_APER"] * u.deg,
                                   'dec': spec_hdr["DEC_APER"] * u.deg,
                                   'PID': spec_hdr0['PROPOSID'],
                                   'radius': float(re.findall(r"[-+]?\d*\.\d+|\d+", spec_hdr["S_REGION"])[2]) * u.deg}

        print('x1dsum file name:', fname, identifier, spec_hdr["RA_APER"], spec_hdr["DEC_APER"], spec_hdr['ROOTNAME'])

    # Show all the images available
    idcs_image = (df_obj.index.get_level_values('state').isin(['flt', 'cal', 'mos', 'drz']) ) & (df_obj['grating'] != 'G185M')
    for i, idx_image in enumerate(df_obj.loc[idcs_image].index):
        if i > 0:
            image_path = obs_folder/df_obj.loc[idx_image].filepath
            print(f'{i}) {image_path.name}')

            # Get figure data
            hdr = fits.getheader(image_path, extname='SCI')
            hdr0 = fits.getheader(image_path, ext=0)
            imdata = fits.getdata(image_path, extname='SCI')

            try:
                wcs = WCS(hdr)
            except:
                wcs = WCS(fits.getheader(image_path, ext=2))

            # print('ACQ IMAGE PID:', hdr0['PROPOSID'])
            #
            # ny, nx = imdata.data.shape
            # center = WCS(hdr).pixel_to_world(nx / 2, ny / 2)
            #
            # wcs_out = WCS(naxis=2)
            # wcs_out.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            # wcs_out.wcs.crval = [center.ra.deg, center.dec.deg]
            # wcs_out.wcs.crpix = [nx / 2, ny / 2]
            # pixscale = 0.023  # arcsec/pix, example — use your COS value
            # wcs_out.wcs.cdelt = np.array([-pixscale / 3600., pixscale / 3600.])

            # wcs_in = WCS(hdr)
            # ny, nx = imdata.data.shape
            #
            # # Center in sky coords from original WCS
            # center = SkyCoord.from_pixel(nx / 2, ny / 2, wcs_in)
            #
            # wcs_out = WCS(naxis=2)
            # wcs_out.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            # wcs_out.wcs.crval = [center.ra.deg, center.dec.deg]  # sky center (deg)
            # wcs_out.wcs.crpix = [nx / 2, ny / 2]  # image center (pix)
            #
            # # Pixel scales from original WCS (deg/pix)
            # pixscale_deg = wcs_in.proj_plane_pixel_scales()  # array [dx_deg, dy_deg]
            #
            # # Option A: assume square pixels
            # cdelt1 = pixscale_deg[0].to(u.deg).value
            # cdelt2 = pixscale_deg[1].to(u.deg).value
            # wcs_out.wcs.cdelt = np.array([-cdelt1, +cdelt2])

            print('  ', df_obj.loc[idx_image].name[1], df_obj.loc[idx_image].filepath, df_obj.loc[idx_image].grating, imdata.shape)

            # Plot the image
            cos_image_plotter(imdata, wcs, object, coord_dict, subgroup=sub_group)
