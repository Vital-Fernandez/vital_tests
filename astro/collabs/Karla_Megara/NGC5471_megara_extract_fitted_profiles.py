import numpy as np
import lime
from astropy.io import fits
from pathlib import Path

folder_name = Path('/home/vital/Downloads')
cfg_name = './NGC5471_cfg.toml'
lime.show_instrument_cfg()
file_list = ['NGC5471_datacube_LR-B_900_scale03_drp_extracted_spectrum_knotA_center_9_9.fits',
             'NGC5471_datacube_LR-R_900_scale03_drp_extracted_spectrum_knotA_center_9_9.fits',
             'NGC5471_datacube_LR-U_900_scale03_drp_nosky_extracted_spectrum_knotA_center_9_9.fits',
             'NGC5471_datacube_LR-V_3600_scale03_drp_extracted_spectrum_knotA_center_9_9.fits']

redshift = 0.00097817
for i, fname in enumerate(file_list):
    print(fname)
    if i == 1:
        with fits.open(folder_name/fname) as hdul:
            data = hdul[1].data
            header = hdul[1].header

            wave = np.asarray(data["Wavelength"], dtype=float)
            flux_jy = np.asarray(data["Flux"], dtype=float)

            spec = lime.Spectrum(wave, flux_jy, redshift=redshift, units_flux='Jy')
            spec.unit_conversion(flux_units_out='FLAM')

            # Generate bands
            ext_root = f'NGC5471_{fname.split('_')[2]}_{fname.split('_')[3]}'
            bands_fname = f'{ext_root}_bands.txt'
            # spec.check.bands(fname=bands_fname, show_continua=True)

            # Run fits
            spec.fit.frame(bands_fname, cfg_name, cont_source='adjacent')

            #Save measurements
            fname = f'{ext_root}_plot.png'
            # spec.plot.grid(fname=fname, rest_frame=True)
            # spec.plot.grid(rest_frame=True)

            fname = f'{ext_root}_line_fluxes.txt'
            spec.save_frame(fname=fname)

            fname = f'{ext_root}_line_profiles.txt'
            spec.save_spectrum(fname=fname, line_label='H1_6563A')




