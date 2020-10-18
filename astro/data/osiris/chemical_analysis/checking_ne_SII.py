import numpy as np
from pathlib import Path
import src.specsiser as sr
import pyneb as pn

# Import the observation data
obsData = sr.loadConfData('../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)

S2 = pn.Atom('S', 2)
O3 = pn.Atom('O', 3)

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
    simFolder, simConf = fitsFolder / 'chemical_analysis',  fitsFile.replace('_BR.fits', '_config.txt')
    inputLinesLog = f'{objName}_inputLinesLog.txt'
    outputDb = f'{objName}_fitting.db'
    outputTxt = f'{objName}_fitting.txt'
    print(f'-{objName}')

    obj1_model = sr.SpectraSynthesizer()
    lm = sr.LineMesurer(linesDF_address=lineLogFolder / lineLogFile)

    blended_dict = obsData['blended_groups']
    blended_list = []
    for group in blended_dict:
        blended_list += blended_dict[group].split('-')

    idcs_blended = lm.linesDF.index.isin(blended_list)

    # Asociate the corresponding flux to the appropiate line
    lm.linesDF.insert(loc=1, column='obsFlux', value=np.nan)
    lm.linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)
    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
    lm.linesDF.loc[idcs_blended, 'obsFlux'] = lm.linesDF.loc[idcs_blended, 'gauss_flux']/flux_Hbeta
    lm.linesDF.loc[idcs_blended, 'obsFluxErr'] = lm.linesDF.loc[idcs_blended, 'gauss_err']/flux_Hbeta
    lm.linesDF.loc[~idcs_blended, 'obsFlux'] = lm.linesDF.loc[~idcs_blended, 'intg_flux']/flux_Hbeta
    lm.linesDF.loc[~idcs_blended, 'obsFluxErr'] = lm.linesDF.loc[~idcs_blended, 'intg_err']/flux_Hbeta
    lm.linesDF.rename(columns={'wavelength': 'obsWave'}, inplace=True)
    # lm.save_lineslog(lm.linesDF, simFolder/inputLinesLog)

    S2_6716_flux, S2_6716_err = lm.linesDF.loc['S2_6716A', 'obsFlux'], lm.linesDF.loc['S2_6731A', 'obsFluxErr']
    S2_6731_flux, S2_6731_err = lm.linesDF.loc['S2_6731A', 'obsFlux'], lm.linesDF.loc['S2_6731A', 'obsFluxErr']

    O3_4363_flux, O3_4363_err = lm.linesDF.loc['O3_4363A', 'obsFlux'], lm.linesDF.loc['O3_4363A', 'obsFluxErr']
    O3_5007_flux, O3_5007_err = lm.linesDF.loc['O3_5007A', 'obsFlux'], lm.linesDF.loc['O3_5007A', 'obsFluxErr']

    neSII = S2.getTemDen(S2_6716_flux/S2_6731_flux, tem=16000, wave1=6716, wave2=6731)
    TeOIII = O3.getTemDen(O3_4363_flux/O3_5007_flux, den=neSII, wave1=4363, wave2=5007)

    print('--', S2_6716_flux/S2_6731_flux, neSII)
    print('--', O3_4363_flux/O3_5007_flux, TeOIII)
