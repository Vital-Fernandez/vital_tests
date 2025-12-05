import numpy as np
import pandas as pd
from astropy.io import fits
import lime

lime.theme.set_style('dark')

new_database = '/home/vital/PycharmProjects/lime/src/lime/resources/lines_database_v2.0.4.txt'
lime.lineDB.set_database(ref_bands=new_database, vacuum_waves=True)

# # Load the NGC 628 fits HDU table
# fname = '/home/vital/Dropbox/Astrophysics/Data/NGC628_NIRSPECMSA/allslits_spectra_extcorr.fits'
# data = fits.getdata(fname, ext=1)
#
# # Convert RecArray to dataframe
# data = np.asarray(data)
# if not data.dtype.isnative:
#     data = data.byteswap().view(data.dtype.newbyteorder('='))
# df = pd.DataFrame.from_records(data)
# df["source"] = df["source"].astype(str)
# df["grating"] = df["grating"].astype(str)
#
# # Extract object spectrum G140M, G235M, and G395M.
# idcs = (df.source == '100307') & (df.slit == 8) & (df.grating == 'G395M')
# df_obj = df.loc[idcs]
#
# spec = lime.Spectrum(input_wave=df_obj['wavelength'].to_numpy(),
#                      input_flux=df_obj['flux_sub'].to_numpy(),
#                      input_err=df_obj['unc_sub'].to_numpy(),
#                      units_wave='um',
#                      units_flux='FLAM',
#                      redshift=0.002192)
#
# # Unit conversion
# spec.unit_conversion('AA')
#
# # Plot the bands
# bands = spec.retrieve.lines_frame()
# spec.plot.spectrum(bands=bands, log_scale=True, rest_frame=True)


file_path = '/home/vital/Astrodata/STScI/LyC_leakers_COS/Direct_downloads/LF9G01010/hst_17515_cos_mrk-209_g130m_lf9g01_cspec.fits'
spec = lime.Spectrum.from_file(file_path,  instrument='cos', redshift=0.000932, norm_flux=1e-17)
bands = spec.retrieve.lines_frame(band_vsigma=20)
spec.plot.spectrum(bands=bands, rest_frame=True)



# # Object bands
# bands = spec.retrieve.lines_frame()
#
# # Spectrum bands
# spec.plot.spectrum(bands=bands, log_scale=True, rest_frame=True)

# # Save to a text file
# spec.retrieve.spectrum('/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/Online_example_data/NGC628_G140M_NIRSPEC.txt')
