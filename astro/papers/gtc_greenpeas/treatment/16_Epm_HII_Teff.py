from pathlib import Path
import numpy as np
import pandas as pd
import src.specsiser as sr
from astro.papers.gtc_greenpeas.common_methods import epm_fitting

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'
cycle = 'it3'

HCm_folder = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/'
grid_file = 'C17_bb_Teff_30-90_pp.dat'

conversion_dict = dict(O2_3726A_m='OII_3727',
                       Ne3_3869A='NeIII_3868',
                       O3_4363A='OIII_4363',
                       O3_4959A='OIII_4959',
                       O3_5007A='OIII_5007',
                       He1_4471A='HeI_4471',
                       He1_5876A='HeI_5876',
                       He2_4686A='HeII_4686',
                       N2_6584A='NII_6584',
                       S2_6716A='SII_6716',
                       S2_6731A='SII_6731')

HCm_conf = dict(n=200,
            sed=2,
            geo=1,
            inter=1,
            grid_file=f'{HCm_folder}/{grid_file}')

# Analyse the spectrum
for i, obj in enumerate(objList):

        if i < 3:

            print(f'- Treating {obj}')

        # --------------------------- Convert linelogs to EPM format ----------------------------
            # Declare input files
            objFolder = resultsFolder/f'{obj}'
            results_file = objFolder/f'{obj}_{ext}_measurements.txt'
            lineLog_file = objFolder/f'{obj}_{ext}_linesLog_{cycle}.txt'
            fitResults_file = objFolder/f'{obj}_{ext}_fitting_{cycle}.txt'

            # Declare output files
            epmLineLog_file = objFolder/f'{obj}_{ext}_epmLinesLog_{cycle}.txt'
            HII_Tef_fit_file = objFolder/f'{obj}_{ext}_HII_Tef_fit_{cycle}.txt'

            # Load the data
            linesDF = sr.lineslogFile_to_DF(lineLog_file)
            results_dict = sr.loadConfData(results_file, group_variables=False)
            fit_results_dict = sr.loadConfData(fitResults_file, group_variables=False)

            # Dataframe for epm lines log
            epm_log = pd.DataFrame()

            #  Compute oxygen abundance
            abund_dict = fit_results_dict['Fitting_results']
            O2_abund = np.power(10, abund_dict['O2'] - 12)
            O3_abund = np.power(10, abund_dict['O3'] - 12)
            OH = np.log10(O2_abund + O3_abund) + 12

            epm_log.loc['0', '12logOH'] = OH[0]
            epm_log.loc['0', 'e12logOH'] = OH[1]

            # Recover the lines observed and label for HII-Chemistry
            for sr_label, epm_label in conversion_dict.items():
                    if sr_label in linesDF.index:
                            epm_log.loc['0', epm_label] = linesDF.loc[sr_label, 'obsInt']
                            epm_log.loc['0', f'e{epm_label}'] = linesDF.loc[sr_label, 'obsIntErr']

            # Safe the dataframe
            with open(epmLineLog_file, 'wb') as output_file:
                    string_DF = epm_log.to_string(index=False)
                    output_file.write(string_DF.encode('UTF-8'))

            # --------------------------- Run Epm fitting ----------------------------

            epm_fitting(str(epmLineLog_file), HII_Tef_fit_file, **HCm_conf, HCm_folder=HCm_folder)

            # Move data to folder
            headers = ['ID', 'O2Hb', 'eO2Hb', 'O3Hb', 'eO3Hb', '4471Hb', 'e4471Hb', '5876Hb', 'e5876Hb', 'He2Hb', 'eHe2Hb', 'S2Hb', 'eS2Hb', 'S3Hb', 'eS3Hb', 'O/H','eO/H','Teff','eTeff','logU','elogU']
            HII_chem_array = np.loadtxt(HII_Tef_fit_file)
            HII_chemDF = pd.DataFrame(columns=headers)
            HII_chemDF.loc[0] = HII_chem_array

            HII_Chem_dict = {'logU': [HII_chemDF.loc[0, 'logU'], HII_chemDF.loc[0, 'elogU']],
                            'Teff': [HII_chemDF.loc[0, 'Teff'], HII_chemDF.loc[0, 'eTeff']],
                            'O/H':  [HII_chemDF.loc[0, 'O/H'], HII_chemDF.loc[0, 'eO/H']]}

            sr.parseConfDict(results_file, HII_Chem_dict, f'HII_Tef_fit_{grid_file}_{cycle}', clear_section=True)

