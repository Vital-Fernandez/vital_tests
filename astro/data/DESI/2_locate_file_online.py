from pathlib import Path
import numpy as np
import pandas as pd
import lime
from lime.read_fits import check_url_status

from astropy.io import fits
from astropy.table import Table
from desi_functions import desi_mask
from matplotlib import pyplot as plt

# Open DESI healpix master file

root_url = 'https://data.desi.lbl.gov/public/edr/spectro/redux/'
id_target = 39627835576420141
release = "fuji"
program = 'dark'

ref_fits = f'zall-pix-{release}.fits'
root_url = f'https://data.desi.lbl.gov/public/edr/spectro/redux/{release}'
ref_fits_url = f'{root_url}/zcatalog/{ref_fits}'

with fits.open(ref_fits_url, use_fsspec=True) as hdul:

    # Index the object
    idx_target = np.where((hdul['ZCATALOG'].data['TARGETID'] == id_target) &
                          (hdul['ZCATALOG'].data['PROGRAM'] == program))[0]
    n_targets = idx_target.sum()

    # Get healpix, survey and redshift
    hpx = hdul['ZCATALOG'].data['HEALPIX'][idx_target]
    survey = hdul['ZCATALOG'].data['SURVEY'][idx_target]
    redshift = hdul["ZCATALOG"].data['Z'][idx_target]

    # Compute the url address
    url_list = []
    for i, idx in enumerate(idx_target):
        hpx_number = hpx[i]
        hpx_ref = f'{hpx_number}'[:-2]
        target_dir = f"/healpix/{survey[i]}/{program}/{hpx_ref}/{hpx_number}"
        coadd_fname = f"coadd-{survey[i]}-{program}-{hpx_number}.fits"
        url_target = f'{root_url}/{target_dir}/{coadd_fname}'

        # check_url_status(url_target)
        url_list.append(url_target)

print(url_list)
print(redshift)

# fujidata = Table(fitsio.read(os.path.join(specprod_dir, "zcatalog", "zall-pix-{}.fits".format(specprod))))

# # -- get all targets with NSPEC=5
# t_fivespec = fujidata[fujidata["ZCAT_NSPEC"] == 5]
#
# # -- unique TARGETID of each object with five spectra
# targids = np.unique(t_fivespec["TARGETID"])

# # -- get the data for all observations of this TARGETID
# special_ID = 39627835576420141
# these_spec = t_fivespec[t_fivespec["TARGETID"] == special_ID]

# # -- get the data for all observations of this TARGETID
# special_ID = 39627835576420141
# these_spec = t_fivespec[t_fivespec["TARGETID"] == special_ID]
#
# # -- get the SURVEY, PROGRAM, SPECTYPE, and redshift values for each of the five spectra of this object
# survey = these_spec["SURVEY"].data.astype(str)
# program = these_spec["PROGRAM"].data.astype(str)
# redshift = np.round(these_spec["Z"].data, 5)
# spectype = these_spec["SPECTYPE"].data.astype(str)
#
# print("\tSURVEY  PROGRAM  SPECTYPE  REDSHIFT")

# idx = np.where((fujidata["TARGETID"] == tid) & (fujidata["SURVEY"] == survey) & (fujidata["PROGRAM"] == program))[0][0]

# # -- healpix values are integers but are converted here into a string for easier access to the file path
# hpx = fujidata["HEALPIX"].astype(str)
#
# specprod_dir = f"/global/cfs/cdirs/desi/spectro/redux/{specprod}"
# target_dir = f"{specprod_dir}/healpix/{survey}/{program}/{hpx[idx][:-2]}/{hpx[idx]}"
# coadd_fname = f"coadd-{survey}-{program}-{hpx[idx]}.fits"
#
# # -- read in the spectra with desispec
# coadd_obj = desispec.io.read_spectra(f"{target_dir}/{coadd_fname}")
# coadd_tgts = coadd_obj.target_ids().data
#
# # -- select the spectrum of  targetid
# row = (coadd_tgts == fujidata["TARGETID"][idx])
# coadd_spec = coadd_obj[row]