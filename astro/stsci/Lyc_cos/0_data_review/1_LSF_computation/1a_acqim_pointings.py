from pathlib import Path
import numpy as np
import re
import lime
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astro.stsci.plots import cos_image_plotter


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Remove excluded objects
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]

# Get the images we want
object_images = list(cfg_sample['favoured_images'].keys())

counter = 0
mean_cord_dict = {}
for j, object in enumerate(object_images):

    df_obj = sample_df.loc[sample_df.object == object]
    idx_image = sample_df.filepath.str.contains(cfg_sample['favoured_images'][object][0])
    image_path = obs_folder / sample_df.loc[idx_image].filepath.values[0]
    instr_im = sample_df.loc[idx_image].instr.values[0]
    grating = sample_df.loc[idx_image].grating.values[0]

    groups_dict = cfg_sample['multi_target_labels'].get(object)
    group_members = [object] if groups_dict is None else list(groups_dict.keys())

    for subObject in group_members:

        # Get the object image
        if counter >= 0:

            if groups_dict is None:
                sub_labels = df_obj.index.get_level_values('id').unique()
            else:
                sub_labels = cfg_sample['multi_target_labels'][object][subObject]

            idcs_spectra = (df_obj.object.index.get_level_values('state') == 'x1dsum')
            idcs_spectra = idcs_spectra & df_obj.index.get_level_values('id').isin(sub_labels)
            spec_list = df_obj.loc[idcs_spectra, 'filepath'].values
            coord_dict = {}

            for spec_fname in spec_list:
                fname = obs_folder/spec_fname
                spec_hdr = fits.getheader(fname, ext=1)
                spec_hdr0 = fits.getheader(fname, ext=0)
                identifier = spec_hdr0['ROOTNAME']
                coord_dict[identifier] =  {'pa': spec_hdr["PA_APER"] * u.deg,
                                           'orien': spec_hdr["ORIENTAT"] * u.deg,
                                           'disp': spec_hdr["DISPAXIS"],
                                           'ra': spec_hdr["RA_APER"] * u.deg,
                                           'dec': spec_hdr["DEC_APER"] * u.deg,
                                           'radius': float(re.findall(r"[-+]?\d*\.\d+|\d+", spec_hdr["S_REGION"])[2]) * u.deg}

            # Get figure data
            hdr = fits.getheader(image_path, extname='SCI')
            imdata = fits.getdata(image_path, extname='SCI')

            # Generate dictionary with coordinates
            mean_cord_dict[subObject] = [float(np.mean([v["ra"].value for v in coord_dict.values()])),
                                         float(np.mean([v["dec"].value for v in coord_dict.values()]))]

            # Get mask line coords
            mask_points = cfg_sample['apperture_masks'].get(object)

            # Plot the image  similar [ 'lc3401m4q', 'lc3401lwq', 'le2401j2q']
            print(f'{counter}) {object} {instr_im}, {grating}')
            cos_image_plotter(imdata, WCS(hdr), object, coord_dict, instr_im=instr_im, mask_points=mask_points,
                              title=f'{subObject} ({len(spec_list)} x1dsum) ')

        counter += 1

if len(mean_cord_dict)> 0:
    lime.save_cfg(project_folder/'/LyC_leakers_COS/aper_cords.toml', mean_cord_dict, section_name='aper_mean_coord')












# for j, object in enumerate(object_images):
#
#     if j >= 0:
#
#         # Get the object image
#         df_obj = sample_df.loc[sample_df.object == object]
#         idx_image = df_obj.filepath.str.contains(cfg_sample['favoured_images'][object][0])
#         image_path = obs_folder/df_obj.loc[idx_image].filepath.values[0]
#         instr_im = df_obj.loc[idx_image].instr.values[0]
#         grating = df_obj.loc[idx_image].grating.values[0]
#
#         # Get the object pointing data
#         idcs_spectra = df_obj.object.index.get_level_values('state') == 'x1dsum'
#         spec_list = df_obj.loc[idcs_spectra, 'filepath'].values
#         coord_dict = {}
#
#         # Sub-identifier for group
#         sub_group = None if object not in cfg_sample['multi_target_labels'] else cfg_sample['multi_target_labels'][object]
#
#         for spec_fname in spec_list:
#             fname = obs_folder/spec_fname
#             spec_hdr = fits.getheader(fname, ext=1)
#             spec_hdr0 = fits.getheader(fname, ext=0)
#             identifier = spec_hdr0['TARGNAME']
#             coord_dict[identifier] =  {'pa': spec_hdr["PA_APER"] * u.deg,
#                                        'orien': spec_hdr["ORIENTAT"] * u.deg,
#                                        'disp': spec_hdr["DISPAXIS"],
#                                        'ra': spec_hdr["RA_APER"] * u.deg,
#                                        'dec': spec_hdr["DEC_APER"] * u.deg,
#                                        'radius': float(re.findall(r"[-+]?\d*\.\d+|\d+", spec_hdr["S_REGION"])[2]) * u.deg}
#
#         # Get figure data
#         hdr = fits.getheader(image_path, extname='SCI')
#         imdata = fits.getdata(image_path, extname='SCI')
#
#         # Get mask line coords
#         mask_points = cfg_sample['apperture_masks'].get(object)
#
#         # Plot the image
#         print(f'{j}) {object} {instr_im}, {grating}')
#         cos_image_plotter(imdata, WCS(hdr), object, coord_dict, subgroup=sub_group, instr_im=instr_im, mask_points=mask_points)