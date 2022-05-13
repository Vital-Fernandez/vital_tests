import numpy as np
import pandas as pd
import lime
from lime.io import load_fits
from pathlib import Path
import os

conf_file = 'LzLCS_ISIS_cfg.ini'
obsCfg = lime.load_cfg(conf_file, obj_section={'sample_data': 'objName_list'}, def_cfg_sec='default_line_fitting')

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_folder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
arm_list = obsCfg['sample_data']['arm_list']
objList = obsCfg['sample_data']['objName_list']
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_ISIS/data/reference_mask.txt'
norm_flux = obsCfg['sample_data']['norm_flux']

for i, specName in enumerate(specNameList):

    if not os.path.exists(results_folder/objList[i]):
        os.mkdir(results_folder/objList[i])

    for arm in arm_list:

        # Input files
        file_spec = dataFolder/objList[i]/f'{specName}_{arm}_f_w_e_flux_nearest.fits'
        file_mask = dataFolder/objList[i]/f'{objList[i]}_{arm}_mask.txt'

        # Load the data
        wave, data, hdr = load_fits(file_spec, instrument='ISIS', frame_idx=0)
        flux = data[0][0]
        mask = lime.load_lines_log(file_mask)
        obj_cfg = obsCfg[f'{objList[i]}_line_fitting']

        # Lime spectrum object
        print(f'- ({i}) {objList[i]}: {arm}')
        spec = lime.Spectrum(wave, flux, redshift=zList[i], norm_flux=norm_flux)
        # spec.plot_spectrum(spec_label=objList[i], frame='rest')

        # Loop throught the lines
        for line in mask.index:
            mask_waves = mask.loc[line, 'w1':'w6'].values
            spec.fit_from_wavelengths(line, mask_waves, obj_cfg)
            spec.display_results(output_address=results_folder/objList[i]/f'{line}_components.png')
            spec.plot_line_velocity(output_address=results_folder/objList[i]/f'{line}_kinematics.png')

        # Join the logs
        if arm == 'Blue':
            log = spec.log.copy()
            mask_arm = mask.copy()
        else:
            log = pd.concat([log, spec.log.copy()])
            mask_arm = pd.concat([mask_arm, mask.copy()])

    # Compute the w80
    idx = log.columns.get_loc('v_95')
    log.insert(idx + 1, 'w80', np.nan)
    log['w80'] = log['v_90'] - log['v_10']

    # Save joined dataframe
    lime.save_line_log(log, results_folder/objList[i]/f'{objList[i]}_log.txt')
    lime.save_line_log(mask_arm, results_folder/objList[i]/f'{objList[i]}_mask.txt')
    lime.save_line_log(log, results_folder/f'sample_lines_log.xlsx', ext=objList[i])

