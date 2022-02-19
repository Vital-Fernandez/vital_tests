import pandas as pd
import pyneb as pn
import src.specsiser as sr
from pathlib import Path
from delete.data_printing import PdfPrinter

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)

dataFolder = Path(obsData['file_information']['data_folder'])
fileList = obsData['file_information']['files_list']

objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])

# Table settings
tableHeaders = [r'$Transition$', '$EW(\AA)$', '$F(\lambda)$', '$I(\lambda)$']
scaleTable = 1000.0
dict_failures = {}

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ('_BR', '_B', '_R'):

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        lineLog_file = dataFolder/'flux_analysis'/f'{obj}{ext}_linesLog.txt'
        pdfTableFile = dataFolder/'flux_analysis'/f'{obj}{ext}_linesTable'
        txtTableFile = dataFolder/'flux_analysis'/f'{obj}{ext}_linesTable.txt'
        print(f'\n- {i}: {lineLog_file}')

        # Set wavelength and flux
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        wave_rest = wave / (1 + z)
        idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

        # Load line measurer object
        lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], linesDF_address=lineLog_file, normFlux=flux_norm)
        pdf = PdfPrinter()
        tableDF = pd.DataFrame(columns=tableHeaders[1:])

        # Measure line fluxes
        idcs_lines = ~lm.linesDF.index.str.contains('_b')
        obsLines = lm.linesDF.loc[idcs_lines].index.values

        # Measure line fluxes
        pdf.create_pdfDoc(pdfTableFile, pdf_type='table')
        pdf.pdf_insert_table(tableHeaders)

        # Normalizing flux
        if 'H1_6563A' in lm.linesDF.index:
            flux_Halpha = lm.linesDF.loc['H1_6563A', 'gauss_flux']
            flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
            halpha_norm = flux_Halpha / flux_Hbeta

            rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
            rc.setCorr(obs_over_theo=halpha_norm / 2.86, wave1=6563., wave2=4861.)
            cHbeta = rc.cHbeta
        else:
            flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
            rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=obsData[obj]['cHbeta'])
            cHbeta = float(rc.cHbeta)

        for lineLabel in obsLines:

            label_entry = lm.linesDF.loc[lineLabel, 'latexLabel']
            wavelength = lm.linesDF.loc[lineLabel, 'wavelength']
            eqw, eqwErr = lm.linesDF.loc[lineLabel, 'eqw'], lm.linesDF.loc[lineLabel, 'eqw_err']

            flux_intg = lm.linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta * scaleTable
            flux_intgErr = lm.linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta * scaleTable
            flux_gauss = lm.linesDF.loc[lineLabel, 'gauss_flux'] / flux_Hbeta * scaleTable
            flux_gaussErr = lm.linesDF.loc[lineLabel, 'gauss_err'] / flux_Hbeta * scaleTable

            if (lm.linesDF.loc[lineLabel, 'blended_label'] != 'None') and ('_m' not in lineLabel):
                flux, fluxErr = flux_gauss, flux_gaussErr
                label_entry = label_entry + '$_{gauss}$'
            else:
                flux, fluxErr = flux_intg, flux_intgErr

            # Correct the flux
            corr = rc.getCorrHb(wavelength)
            intensity, intensityErr = flux * corr, fluxErr * corr

            eqw_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(eqw, eqwErr)
            flux_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(flux, fluxErr)
            intensity_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(intensity, intensityErr)

            # Add row of data
            row_i = [label_entry, eqw_entry, flux_entry, intensity_entry]
            lastRow_check = True if lineLabel == obsLines[-1] else False

            pdf.addTableRow([label_entry, eqw_entry, flux_entry, intensity_entry], last_row=lastRow_check)
            tableDF.loc[lineLabel] = row_i[1:]

        # Data last rows
        row_Hbetaflux = [r'$H\beta$ $(erg\,cm^{-2} s^{-1} \AA^{-1})$', '', flux_Hbeta * flux_norm,
                         flux_Hbeta * flux_norm * rc.getCorr(4861)]
        row_cHbeta = [r'$c(H\beta)$', '', cHbeta, '']

        pdf.addTableRow(row_Hbetaflux, last_row=False)
        pdf.addTableRow(row_cHbeta, last_row=False)
        tableDF.loc[row_Hbetaflux[0]] = row_Hbetaflux[1:]
        tableDF.loc[row_cHbeta[0]] = row_cHbeta[1:]

        # Format last rows
        pdf.table.add_hline()
        pdf.table.add_hline()

        # Save the pdf table
        try:
            pdf.generate_pdf(clean_tex=True)
        except:
            print('-- PDF compilation failure')
            dict_failures[i] = fits_file

        # Save the txt table
        with open(txtTableFile, 'wb') as output_file:
            string_DF = tableDF.to_string()
            string_DF = string_DF.replace('$', '')
            output_file.write(string_DF.encode('UTF-8'))

        # Increase counter for obj number
        counter += 1

print('Total errors', len(dict_failures.keys()))
print(dict_failures.keys())
print(dict_failures.values())


