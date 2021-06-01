import numpy as np
from pathlib import Path
import src.specsiser as sr

spec_address = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/gp030321_B.fits')
wave, flux_array, header = sr.import_fits_data(spec_address, instrument='OSIRIS')
flux, err, normFlux = flux_array[0][0], flux_array[3][0], 1e-17
lm_osiris = sr.LineMesurer(wave, flux, input_err=err, normFlux=normFlux)
lm_osiris.plot_spectrum(continuumFlux=lm_osiris.errFlux, log_scale=True, axConf={'title': 'OSIRIS spectrum'})
lineWaves = np.array([5600.0, 5635.0, 5651.0, 5675.0, 5697.0, 5729.0])
lm_osiris.fit_from_wavelengths('H1_4861A', lineWaves)
lm_osiris.print_results(show_plot=True)
#
spec_address = Path('/home/vital/Astro-data/Observations/MUSE - Amorin/CGCG007.fits') #D:/Google drive/Astrophysics/Datos/MUSE - Amorin/CGCG007.fits')
wave, cube, header = sr.import_fits_data(spec_address, instrument='MUSE')
idx_j, idx_i = 171, 171
flux_voxel = cube[:, idx_j, idx_i].data.data
flux_err = cube[:, idx_j, idx_i].var.data
lm_muse = sr.LineMesurer(wave, flux_voxel, input_err=flux_err)
lm_muse.plot_spectrum(continuumFlux=lm_muse.errFlux, log_scale=True, axConf={'title': 'MUSE spectrum'})
lineWaves = np.array([4835.0, 4868.0, 4877.0, 4891.0, 4913.0, 4941.0])
lm_muse.fit_from_wavelengths('H1_4861A', lineWaves)
lm_muse.print_results(show_plot=True)
#
spec_address = Path('/home/vital/Dropbox/Astrophysics/Data/xshooter-multicomponent/fJ0925sum.05_s.fits')
err_address = Path('/home/vital/Dropbox/Astrophysics/Data/xshooter-multicomponent/fJ0925sum.05_e.fits')
wave, flux_array, header = sr.import_fits_data(spec_address, instrument='OSIRIS')
wave_e, err_array, header_e = sr.import_fits_data(err_address, instrument='OSIRIS')
normFlux = 1e-18
lm_xshooter = sr.LineMesurer(wave, flux_array[0], input_err=err_array[0], normFlux=normFlux)
lm_xshooter.plot_spectrum(continuumFlux=lm_xshooter.errFlux, log_scale=True, axConf={'title': 'XSHOOTER spectrum'})
lineWaves = np.array([6055.0, 6074.0, 6086.0, 6115.0, 6122.0, 6141.0])
lm_xshooter.fit_from_wavelengths('H1_4861A', lineWaves)
lm_xshooter.print_results(show_plot=True)

spec_address = Path('/home/vital/Astro-data/Observations/broad_multiComponent/B6479s.fits')
wave_array, flux_array, header = sr.import_fits_data(spec_address, instrument='ISIS', frame_idx=0)
wave, flux, err = wave_array, flux_array[0][0], flux_array[1][0]
normFlux = 1e-18
lm_isis = sr.LineMesurer(wave, flux, input_err=err, normFlux=normFlux)
lm_isis.plot_spectrum(continuumFlux=lm_isis.errFlux, log_scale=False, axConf={'title': 'ISIS spectrum'})
lineWaves = np.array([6290.0, 6310.0, 6318.0, 6332.0, 6340.0, 6360.0])
lm_isis.fit_from_wavelengths('H1_4861A', lineWaves)
lm_isis.print_results(show_plot=True)

# /home/vital/Astro-data/broad_multiComponent

spec_address = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/gp030321_SDSS.fits')
wave, data_dict, header = sr.import_fits_data(spec_address, instrument='SDSS')
flux, err, normFlux = data_dict['flux'], data_dict['ivar'], 1
lm_osiris = sr.LineMesurer(wave, flux, input_err=err, normFlux=normFlux)
lm_osiris.plot_spectrum(continuumFlux=lm_osiris.errFlux, log_scale=True, axConf={'title': 'SDSS spectrum'})
lineWaves = np.array([5600.0, 5635.0, 5651.0, 5675.0, 5697.0, 5729.0])
lm_osiris.fit_from_wavelengths('H1_4861A', lineWaves)
lm_osiris.print_results(show_plot=True)


# import numpy as np
# from pathlib import Path
# import src.specsiser as sr
#
# spec_address = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/gp030321_B.fits')
# wave, flux_array, header = sr.import_fits_data(spec_address, instrument='OSIRIS')
# flux, err = flux_array[0][0], flux_array[3][0]
# lm_osiris = sr.LineMesurer(wave, flux, input_err=err)
# lm_osiris.plot_spectrum_components(continuumFlux=err, log_scale=True, axConf={'title': 'OSIRIS spectrum'})
#
# spec_address = Path('D:/Google drive/Astrophysics/Datos/MUSE - Amorin/CGCG007.fits')
# wave, cube, header = sr.import_fits_data(spec_address, instrument='MUSE')
# idx_j, idx_i = 171, 171
# flux_voxel = cube[:, idx_j, idx_i].data.data
# flux_err = cube[:, idx_j, idx_i].var.data
# lm_muse = sr.LineMesurer(wave, flux_voxel, input_err=flux_err)
# lm_muse.plot_spectrum_components(continuumFlux=flux_err, log_scale=True, axConf={'title': 'MUSE spectrum'})
#
# spec_address = Path('D:/Dropbox/Astrophysics/Data/xshooter-multicomponent/fJ0925sum.05_s.fits')
# err_address = Path('D:/Dropbox/Astrophysics/Data/xshooter-multicomponent/fJ0925sum.05_e.fits')
# wave, flux_array, header = sr.import_fits_data(spec_address, instrument='OSIRIS')
# wave_e, err_array, header_e = sr.import_fits_data(err_address, instrument='OSIRIS')
# lm_osiris = sr.LineMesurer(wave, flux_array[0])
# lm_osiris.plot_spectrum_components(continuumFlux=err_array[0], log_scale=True, axConf={'title': 'XSHOOTER spectrum'})
#
# spec_address = Path('D:/Google drive/Astrophysics/Datos/broad_multiComponent/B6479s.fits')
# wave_array, flux_array, header = sr.import_fits_data(spec_address, instrument='ISIS')
# wave, flux, err = wave_array, flux_array[0][0], flux_array[1][0]
# lm_isis = sr.LineMesurer(wave, flux, input_err=err)
# lm_isis.plot_spectrum_components(continuumFlux=err, log_scale=False, axConf={'title': 'ISIS spectrum'})