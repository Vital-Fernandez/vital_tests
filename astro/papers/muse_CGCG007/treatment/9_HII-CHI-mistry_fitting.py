import numpy as np
import pandas as pd
import lime as lm
import src.specsiser as sr
from pathlib import Path
import time
from astropy.io import fits
from astropy.table import Table
from astro.papers.gtc_greenpeas.common_methods import epm_HII_CHI_mistry
from astro.papers.muse_CGCG007.muse_CGCG007_methods import grid_HII_CHI_mistry_conversion as labelConver


# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
hcm_folder = Path(obsData['data_location']['hcm_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

# Load grid interpolation functions
# grid_file = Path(r'D:\Dropbox\Astrophysics\Tools\HCm-Teff_v5.01\C17_bb_Teff_30-90_pp.dat')
# grid_columns = ['logOH', 'Teff', 'logU', 'OII_3727', 'OIII_5007', 'SII_6725', 'SIII_9069', 'He1_4471', 'He1_5876', 'He2_4686']
# grid_DF = pd.read_csv(grid_file, index_col=False, delim_whitespace=True, comment='#', names=grid_columns)
# model_variables = ['logOH', 'Teff', 'logU']


grid_file = hcm_folder/'C17_popstar_logU-NO_adapted_emp_v5.0.dat'
grid_columns = ['logOH', 'logNO', 'logU', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6725']
grid_DF = pd.read_csv(grid_file, index_col=False, delim_whitespace=True, comment='#', names=grid_columns)
model_variables = ['logOH', 'logU', 'logNO']
# grid_DF.query("7.5 <= logOH <= 7.6 and -3 >= logU <= -2.75")
# grid_DF.query("7.5 <= logOH <= 7.6")

axes_cords = {}
for i, ax_name in enumerate(model_variables):
    axes_cords[ax_name] = np.unique(grid_DF[ax_name].values)

result = 7.3, -2.55, -1.500
for i, value in enumerate(result):
    idx_i = np.searchsorted(axes_cords, value)

# Interpolator functions
gw = sr.ModelGridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_DF, axes_columns=model_variables)
grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')


for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fitsLog_address = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Output data
    HIIchimistry_fits = objFolder / f'{obj}_HIIchimistry.fits'
    hdul_lineslog = fits.HDUList()

    # Loop throught the line regions
    start = time.time()
    for idx_region in [0]:
    # for idx_region in [0, 1, 2]:

        region_label = f'region_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Loop through the region voxels
        n_voxels = idcs_voxels.shape[0]
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            print(f'-- Treating voxel {idx_voxel}/{n_voxels} ({idx_pair})')
            idx_j, idx_i = idx_pair
            logLabel = f'{idx_j}-{idx_i}_linelog'
            # output_fit_file = objFolder/'voxel_treatments'/f'{idx_j}-{idx_i}_HII_CHI_mistry_fit.txt'

            # Load voxel lines log
            linesLog_BinTable = fits.getdata(fitsLog_address, logLabel, ver=1)
            linesDF = Table(linesLog_BinTable).to_pandas()
            linesDF.set_index('index', inplace=True)

            # Prepare data for HII-CHI-mistry
            idcs_inputLines = linesDF.index.isin(labelConver.keys())
            input_lines = linesDF.loc[idcs_inputLines].index.values

            HII_CHI_mistry_DF = pd.DataFrame()
            HII_CHI_mistry_DF.loc[0, 'ID'] = logLabel
            flux_Hbeta = linesDF.loc['H1_4861A', 'intg_flux']
            for lineLabel in input_lines:
                HII_CHI_mistry_DF.loc[0, labelConver[lineLabel]] = linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta
                HII_CHI_mistry_DF.loc[0, f'e{labelConver[lineLabel]}'] = linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta
            lineSA = HII_CHI_mistry_DF.to_records(index=False) #column_dtypes=default_linelog_types, index_dtypes='<U50')

            # Run HII-CHI-mistry
            outputSA = epm_HII_CHI_mistry(lineSA, output_file='None', n=200, sed=1, inter=1, HCm_folder=hcm_folder)
            linesCol = fits.ColDefs(outputSA)
            linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_HIIchimistry')
            hdul_lineslog.append(linesHDU)

        # Store the drive
        hdul_lineslog.writeto(HIIchimistry_fits, overwrite=True, output_verify='fix')
    end = time.time()

print(f'- Execution time {(end - start)/60:.3f} min')