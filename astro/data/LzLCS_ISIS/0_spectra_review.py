import numpy
import lime
import numpy as np
from lime.io import load_fits
from pathlib import Path
from shutil import copy as shu_copy

conf_file = 'LzLCS_ISIS_cfg.ini'
obsCfg = lime.load_cfg(conf_file)

dataFolder = Path(obsCfg['data_location']['data_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
arm_list = obsCfg['sample_data']['arm_list']
objList = obsCfg['sample_data']['objName_list']
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_ISIS/data/reference_mask.txt'
norm_flux = obsCfg['sample_data']['norm_flux']

for i, specName in enumerate(specNameList):

    # Join the arms spectra
    for arm in arm_list:
        file_name = f'{specName}_{arm}_f_w_e_flux_nearest.fits'
        wave, data, hdr = load_fits(dataFolder/objList[i]/file_name, instrument='ISIS', frame_idx=0)
        flux = data[0][0]

        if arm == 'Blue':
            wave_joined = wave
            flux_joined = flux
        else:
            wave_joined = np.concatenate([wave_joined, wave])
            flux_joined = np.concatenate([flux_joined, flux])

    # Treat the object
    obj_mask = dataFolder/objList[i]/f'{objList[i]}_Blue_mask.txt'
    # spec = lime.Spectrum(wave_joined, flux_joined, redshift=zList[i], norm_flux=norm_flux)
    # spec.plot_spectrum(spec_label=objList[i], frame='rest')
    lime.MaskInspector(obj_mask, wave_joined, flux_joined, redshift=zList[i], norm_flux=norm_flux, y_scale='natural',
                       n_cols=3, n_rows=2)

# for i, specName in enumerate(specNameList):
#     for arm in arm_list:
#
#         if i > 6:
#
#             # Load the spectra
#             file_name = f'{specName}_{arm}_f_w_e_flux_nearest.fits'
#             wave, data, hdr = load_fits(dataFolder/objList[i]/file_name, instrument='ISIS', frame_idx=0)
#             flux = data[0][0]
#
#             # Lime spectrum object
#             print(f'- ({i}) {objList[i]}: {arm}')
#             spec = lime.Spectrum(wave, flux, redshift=zList[i], norm_flux=norm_flux)
#             spec.plot_spectrum(spec_label=objList[i], frame='rest')
#
#             # Adjust mask to object
#             obj_mask = dataFolder/objList[i]/f'{objList[i]}_{arm}_mask.txt'
#             # shu_copy(refMask, obj_mask)
#             lime.MaskInspector(obj_mask, wave, flux, redshift=zList[i], norm_flux=norm_flux, y_scale='natural')
#
