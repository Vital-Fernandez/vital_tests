from pathlib import Path
from matplotlib import pyplot as plt
import lime
from astropy.io import fits
import numpy as np
from lime.transitions import _DATABASE_FILE

# linedb = lime.load_frame(_DATABASE_FILE)
# lime.save_frame(Path(_DATABASE_FILE).parent/'lines_database_v2.0.3.xlsx', linedb)
#
# lime.lineDB
# new_database = '/home/vital/PycharmProjects/lime/src/lime/resources/lines_database_v2.0.4.txt'
# lime.lineDB.set_database(ref_bands=new_database, vacuum_waves=True)
#
# # Data location
# project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
# obs_folder = Path('/home/vital/Astrodata/STScI')
#
# # Configuration
# sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
# cfg_sample = lime.load_cfg(project_folder/'samples.toml')
#
# target_single = ['UGC4483', 'VIIZw403', 'NGC2366', 'UGCA281', 'Pox186']
# target_single = ['UGCA281']
# idcs_manual = sample_df.object.isin(target_single) & (sample_df.index.get_level_values('state') == 'hasp_mast')
# obj_list = sample_df.loc[idcs_manual, 'object'].unique()
#
# for i, obj_name in enumerate(obj_list):
#     n_objects = sample_df.loc[(sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'hasp_mast')].index.size
#     idcs = (sample_df.object == obj_name) & idcs_manual
#     sub_sample = sample_df.loc[idcs]
#     n_cspec_manual = sub_sample.index.size
#
#     # Manual cspec
#     if sub_sample.index.size == 0:
#         print(f'-{obj_name} NO SPECTRA')
#
#     else:
#         redshift = cfg_sample['Galaxy_redshifts'][obj_name]
#         print(f'- {obj_name} at z = {redshift}')
#
#         for j, idx in enumerate(sub_sample.index):
#             file_path = obs_folder/sub_sample.loc[idx, 'filepath']
#
#             if j == 0:
#                 print(fits.info(file_path))
#                 spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
#                 label = f'{obj_name} - {idx[2]} - {idx[3]}'
#                 bands = spec.retrieve.lines_frame()
#                 spec.plot.spectrum(bands=bands, label=label, maximize=True)



# lime.lineDB
new_database = '/home/vital/PycharmProjects/lime/src/lime/resources/lines_database_v2.0.4.txt'
lime.lineDB.set_database(ref_bands=new_database, vacuum_waves=True)

file_path = '/home/vital/Astrodata/STScI/LyC_leakers_COS/Direct_downloads/LF9G01010/hst_17515_cos_mrk-209_g130m_lf9g01_cspec.fits'
spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=0.000932, norm_flux=1e-17)
bands = spec.retrieve.lines_frame()
spec.plot.spectrum(bands=bands)

