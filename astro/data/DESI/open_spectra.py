import numpy as np
from astropy.io import fits
from pathlib import Path
import lime
import requests


def open_desi_spectra(file_path, obj_idtarget=None, obj_idrows=None):

    # Confirm the user provides and identity
    if (obj_idtarget is None) and (obj_idrows is None):
        assert (obj_idtarget is None) and (obj_idrows is None), 'You need to apply introduce and IDTARGET or .fits rows'

    obj_idtarget = np.atleast_1d(obj_idtarget)

    # Container for output spectra data:
    spectra_dict = {}

    # Confirm file location
    file_path = Path(file_path)
    if not file_path.is_file():
        assert file_path, f'Input file not found at "{file_path}"'

    # Open the file
    with fits.open(fits_path) as hdulist:

        # Get target ID rows if none are provided  by TARGETIDs
        if obj_idrows is None:
            file_idtargets = hdulist["FIBERMAP"].data['TARGETID']

            obj_idrows = np.where(np.isin(file_idtargets, obj_idtarget))[0]

        # Check the target is on file
        assert obj_idrows.size > 0, f'Input TARGETID(s): {obj_idtarget},\nnot found in input file: {file_path}'

        for hdu in hdulist:
            print(hdu.header)

        # Load the EXP_FIBERMAP
        try:
            exp_targetid_table = hdulist["EXP_FIBERMAP"].data['TARGETID']
            exp_rows = np.where(np.isin(exp_targetid_table, obj_idtarget))[0]
        except KeyError:
            exp_targetid_table = None

        for i, id_obj in enumerate(obj_idtarget):
            spectra_dict[id_obj] = {}
            for band in ['B', 'R', 'Z']:
                spectra_dict[id_obj][band] = {'wave': hdulist[f'{band}_WAVELENGTH'].data,
                                              'flux': hdulist[f'{band}_FLUX'].data[obj_idrows[i], :],
                                              'ivar': hdulist[f'{band}_IVAR'].data[obj_idrows[i], :]}

    return spectra_dict


def open_desi_spectra_online(file_path, obj_idtarget=None, obj_idrows=None):

    r = requests.head(file_path)

    if r.status_code == 200:
        obj_idtarget = np.atleast_1d(obj_idtarget)

        with fits.open(file_path, use_fsspec=True, fsspec_kwargs={"anon": True}) as hdulist:

            file_idtargets = hdulist["FIBERMAP"].data['TARGETID']
            obj_idrows = np.where(np.isin(file_idtargets, obj_idtarget))[0]

            spectra_dict = {}
            for i, id_obj in enumerate(obj_idtarget):
                spectra_dict[id_obj] = {}
                for band in ['B', 'R', 'Z']:
                    spectra_dict[id_obj][band] = {'wave': hdulist[f'{band}_WAVELENGTH'].data,
                                                  'flux': hdulist[f'{band}_FLUX'].section[obj_idrows[i], :],
                                                  'ivar': hdulist[f'{band}_IVAR'].section[obj_idrows[i], :]}

            # cutout = hdul[1].section[10:20, 30:50]

    return spectra_dict


specprod_dir = '.'
special_ID = 39627835576420141
fits_path = f'{specprod_dir}/coadd-sv1-other-27256.fits'

# spectra_dict = open_desi_spectra(fits_path, obj_idtarget=special_ID, obj_idrows=None)
# wave = spectra_dict[special_ID]['B']['wave']
# flux = spectra_dict[special_ID]['B']['flux']
# spec = lime.Spectrum(wave, flux)
# spec.plot.spectrum()

fits_url = 'https://data.desi.lbl.gov/public/edr/spectro/redux/fuji/healpix/sv1/other/272/27256/coadd-sv1-other-27256.fits'
spectra_dict = open_desi_spectra_online(fits_url, obj_idtarget=special_ID, obj_idrows=None)
wave = spectra_dict[special_ID]['B']['wave']
flux = spectra_dict[special_ID]['B']['flux']
spec = lime.Spectrum(wave, flux)
spec.plot.spectrum()

