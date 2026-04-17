import numpy as np
import lime as lime

lime.lineDB.reset('./lines_database_v2.0.5.txt')
lime.lineDB.set_database(vacuum_waves=True)

fits_path = 'Spectra_R10000.txt'
fit_cfg = lime.load_cfg("LIME_R10000_cfg_v2.toml")

wave_arr, flux_arr = np.loadtxt(fits_path, unpack=True, skiprows=1)

spec = lime.Spectrum.from_file(fits_path, instrument='text', redshift=0, skiprows=1, units_flux="1e-17*FLAM")

bands_df = spec.retrieve.lines_frame(band_vsigma=70,
                                     grouped_lines=list(fit_cfg['default_line_fitting'].keys()),
                                     fit_cfg=fit_cfg)

spec.fit.frame(bands_df, fit_cfg)

spec.plot.spectrum()