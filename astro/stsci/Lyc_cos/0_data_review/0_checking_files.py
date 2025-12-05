# import lime
# from astropy.io import fits
from matplotlib import pyplot as plt
#
# fname = '/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC leakers COS/observations/COS_proposal/LF9G01010/lf9g01010_x1dsum.fits'
# spec = lime.Spectrum.from_file(fname, instrument='cos, redshift=0)
# spec.plot.spectrum()
import lime

from astropy.io import fits

# Observations files (X1D)
# file1 = '/home/vital/Astro-data/STScI/LyC_leakers_COS/Calibrations/lf9g01wpq_x1d.fits'
# file2 = '/home/vital/Astro-data/STScI/LyC_leakers_COS/Calibrations/lf9g01wrq_x1d.fits'

# MAST results (X1D)
# file1 = '/home/vital/Astrodata/STScI/LyC_leakers_COS/Calibrations/hasp_products/hst_17515_cos_mrk-209_g130m_lf9g01_cspec.fits'
# file2 = '/home/vital/Astrodata/STScI/LyC_leakers_COS/Calibrations/hasp_products/hst_17515_cos_mrk-209_g130m_lf9g_cspec.fits'

# file1 ='/home/vital/Astrodata/STScI/LyC_leakers_COS/Direct_downloads/LEW206010/lew206biq_x1d.fits'
# file2 ='/home/vital/Astrodata/STScI/LyC_leakers_COS/Direct_downloads/LEW206010/lew206beq_x1d.fits'

file1 ='/home/vital/Astrodata/STScI/LyC_leakers_COS/Calibrations/hasp_products/hst_15352_cos_knot-a_g130m_ldi707_cspec.fits'
file2 ='/home/vital/Astrodata/STScI/LyC_leakers_COS/Calibrations/hasp_products/hst_15352_cos_knot-b_g130m_ldi7_cspec.fits'

# open files
with fits.open(file1) as hdul1, fits.open(file2) as hdul2:
    hdr1 = hdul1[0].header
    hdr2 = hdul2[0].header

    # compare keys in both
    for key in set(hdr1.keys()) | set(hdr2.keys()):  # union of keys
        val1 = hdr1.get(key, None)
        val2 = hdr2.get(key, None)
        if val1 != val2:
            print(f"{key}: {val1}  !=  {val2}")


spec1 = lime.Spectrum.from_file(file1, instrument='cos', redshift=0.000932)
spec2 = lime.Spectrum.from_file(file2, instrument='cos', redshift=0.000932)
spec1.plot.spectrum(label='file1', in_fig=None)
spec1.plot.ax.step(spec2.wave, spec2.flux, where='mid', label='file2', linestyle=':')
spec1.plot.ax.legend()
plt.tight_layout()
plt.show()


# # fits.info(fname)
# #
# # Open the FITS file
# with fits.open(fname) as hdul:
#     hdul.info()  # Show the structure of the file
#
#     # The science data are usually in extension 1
#     data = hdul[1].data
#     header0 = hdul[0].header
#     header1 = hdul[1].header
#
# # Example: inspect columns
# print(data.columns)
#
# # Access wavelength, flux, and error arrays
# wavelength = data['WAVELENGTH'][0]  # [0] because rows contain arrays
# flux = data['FLUX'][0]
# error = data['ERROR'][0]
# print(header0['TARGNAME'])
# print("First 10 wavelength points:", wavelength[:10])
# print("First 10 flux points:", flux[:10])
# fig, ax = plt.subplots()
# ax.step(wavelength, flux)
# plt.show()
