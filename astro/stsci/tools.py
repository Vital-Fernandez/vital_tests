import urllib
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from hasp import wrapper
from lmfit.models import GaussianModel, LinearModel
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion
import shutil

### trying scipy instead ####
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(((x - mean) / stddev) ** 2) / 2)


### trying the function implemented in IDL gaussfit ####
def fit_func(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z ** 2 / 2) + a3 + a4 * x + a5 * x ** 2
    return y


def run_hasp_wrapper(in_folder, out_folder, cross_program=True):

    # Run HASP wrapper
    if out_folder.is_dir():
        print(f'\n2) Cross program ', out_folder)
        wrapper.main(in_folder, outdir=out_folder, clobber=True, cross_program=cross_program)
    else:
        raise BrokenPipeError(f'Output dir not found {out_folder}')

    return


def list_files_with_extension(base_dir, extension):
    """
    Return a list of files with a given extension under base_dir, recursively.

    Args:
        base_dir (str or Path): The directory to search.
        extension (str): File extension (e.g. ".txt", ".fits").

    Returns:
        list[Path]: List of Path objects pointing to the files.
    """
    base_path = Path(base_dir)
    return [f for f in base_path.rglob(f"*{extension}") if f.is_file()]


def add_cos_obs(file_list, df, sample_name, state, ref_dict):

    for spec_file in file_list:
        print(spec_file)
        hdr0 = fits.getheader(spec_file, ext=0)
        hdr1 = fits.getheader(spec_file, ext=1)

        row_data = {}
        row_data['sample'] = f'{sample_name}_{hdr0['PROPOSID']}'
        row_data['id'] = f"{hdr0['TARGNAME']}"

        # Extra ID:
        if 'ASN_ID' in hdr0:
            offset_id = f"{hdr0['ASN_ID']}"
        elif 'IPPPSSOO' in hdr0:
            offset_id = f"{hdr0['IPPPSSOO']}_{hdr0['GRATING']}"
        else:
            offset_id = f"acqim_{hdr0['ROOTNAME']}"

        if state == 'x1d':
            offset_id = f"{hdr0['ROOTNAME']}"

        if state == 'counts':
            offset_id = f"{spec_file.parts[-2]}_{hdr0['ROOTNAME']}"

        if state in ['flt', 'cal']:
            offset_id = f"{spec_file.parts[-2]}_{hdr0['ROOTNAME']}"

        if 'MULTI' in row_data['id']:
            offset_id = f"NUM_EXP{hdr0['NUM_EXP']}_PINAME{hdr0['PINAME']}"

        row_data['offset_id'] = offset_id

        row_data['state'] = state

        # row_data['object'] = hdr0['TARGNAME'].lower()
        if hdr0['TARGNAME'] in ref_dict['Galaxy_alias']:
            row_data['object'] = hdr0['TARGNAME']
        else:
            row_data['object'] = next((k for k, v in ref_dict['Galaxy_alias'].items() if hdr0['TARGNAME'] in v), None)

        if row_data['object'] is None:
            print(hdr0['TARGNAME'])

        if 'MULTI' in row_data['id']:
            row_data['object'] = row_data['id'].split('_')[0]

        if row_data['object'] is None and state == 'counts':
            print('MISSINGGGG ', row_data['id'])

        row_data['RA'] = f"{hdr0['RA_TARG']}" if 'RA_TARG' in hdr0 else None
        row_data['DEC'] = f"{hdr0['DEC_TARG']}" if 'RA_TARG' in hdr0 else None
        row_data['PID'] = hdr0['PROPOSID']
        row_data['redshift'] = ref_dict['Galaxy_redshifts'].get(row_data['object'])

        row_data['instr'] =  hdr0['instrume']

        # Grating
        if 'OPT_ELEM' in hdr0:
            row_data['grating'] = f"{hdr0['OPT_ELEM']}"
        elif 'GRATING' in hdr0:
            row_data['grating'] = hdr0['GRATING']
        else:
            row_data['grating'] = f"{hdr0['INSTRUME']}"

        row_data['filecode'] =  hdr0['ASN_ID'] if 'ASN_ID' in hdr0  else spec_file.parts[-2]
        row_data['filepath'] =  Path(*spec_file.parts[5:])

        row_data['cenwave'] = hdr0.get("CENWAVE")
        row_data['life_adj'] = hdr0.get("LIFE_ADJ")
        row_data['disptab'] = hdr0.get("DISPTAB")
        row_data['detector'] = hdr0.get("DETECTOR")

        # Append
        if row_data['grating'] != 'G140L':
            if not isinstance(df.index, pd.MultiIndex):   # True
                df.loc[len(df)] = row_data
            else:
                idx = (row_data["sample"], row_data["id"], row_data["offset_id"], row_data["state"])
                if idx not in df.index:
                    df.loc[idx] = {k: v for k, v in row_data.items() if k not in df.index.names}

    return


def fetch_files(det, grating, lpPos, cenwave, disptab, datadir):


    """
    Given all the inputs, this will download both
    the LSF and Disptab files to use in the convolution and return their paths.

    Input:
    det (str): The detector used
    grating (str): Type of grating used
    lpPos (str): Lifetime position used
    cenwave (str): Central wavelength used
    disptab (str): DISPTAB used (will get the path in the function)

    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file
    """

    # Link to where all the files live
    COS_site_rootname = ("https://www.stsci.edu/files/live/sites/www/files/"
                        "home/hst/instrumentation/cos/"
                        "performance/spectral-resolution/_documents/")

    # Only one file for NUV
    if det == "NUV":
        LSF_file_name = "nuv_model_lsf.dat"

    # FUV files follow a naming pattern
    elif det == "FUV":
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    # Where to find file online
    LSF_file_webpath = COS_site_rootname + LSF_file_name
    urllib.request.urlretrieve(LSF_file_webpath, str(datadir / LSF_file_name))

    # Where to save file to locally
    print(f"Downloaded LSF file to {str(datadir / LSF_file_name)}")

    # # And we'll need to get the DISPTAB file as well
    # disptab_path = str(datadir / disptab)
    # urllib.request.urlretrieve(f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}", disptab_path)
    #
    # print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name


def image_gauss_fitting(y):

    # x grid (replace with your wavelength array if you have one)
    x = np.arange(len(y), dtype=float)

    # ------------------------------
    # Build model: two Gaussians + linear continuum (m*x + b)
    # ------------------------------
    g1 = GaussianModel(prefix='g1_')
    cont = LinearModel(prefix='c_')  # parameters: c_slope, c_intercept
    model = g1 + cont
    params = model.make_params()

    # Rough linear baseline from endpoints (robust to peaks)
    m0 = (np.median(y[-2:]) - np.median(y[:2])) / (x[-1] - x[0] + 1e-12)
    b0 = np.median(y[:2]) - m0 * x[0]
    params['c_slope'].set(m0, vary=False)
    params['c_intercept'].set(b0, vary=False)

    # Centers (bounded to data range)
    idx_max = np.argmax(y)
    span = (x.max() - x.min()) / 8.0
    params['g1_amplitude'].set(value=y[idx_max], min=0)
    params['g1_center'].set(value=idx_max, min=x.min(), max=x.max())
    params['g1_sigma'].set(value=span, min=1)

    # Fit
    result = model.fit(y, params, x=x)

    # Plot
    xx = np.linspace(x.min(), x.max(), 2000)
    y_fit = result.eval(x=xx)
    comp = model.eval_components(params=result.params, x=xx)

    # Results
    arr_dict = {'c_': comp['c_'], 'g1_': comp['g1_'], 'g1_fwhm': result.params['g1_fwhm'].value}

    return arr_dict


def generate_apperture_mask(imshape, wcs, ra, dec, radius, mask_points):

    # Generate the circular mask
    sky_region = CircleSkyRegion(center=SkyCoord(ra=ra, dec=dec, frame="icrs"), radius=radius)
    circle_pix = sky_region.to_pixel(wcs=wcs)
    circle_mask = circle_pix.to_mask(mode="center").to_image(imshape).astype(bool)

    if mask_points is not None:

        # Convert points to coordinates
        c1 = SkyCoord(mask_points[0][0], mask_points[0][1])
        c2 = SkyCoord(mask_points[1][0], mask_points[1][1])

        # Generate the line mask
        yy, xx = np.mgrid[:imshape[0], :imshape[1]]
        x1, y1 = wcs.world_to_pixel(c1)
        x2, y2 = wcs.world_to_pixel(c2)
        sign = (x2 - x1) * (yy - y1) - (y2 - y1) * (xx - x1)

        # Combine them and create an overlay
        circle_mask = circle_mask & (sign > 0)
        # overlay = np.zeros((*side_mask.shape, 4))
        # overlay[..., 3] = 0.0  # fully transparent everywhere
        # overlay[side_mask] = [0, 0, 0, 1]  # black & opaque where mask is True

    return circle_mask


def move_files(file_list, src_root_path, dest_path):

    # Clear the folder if it already exists:
    if dest_path.exists():
        shutil.rmtree(dest_path)

    # Recreate the folder
    dest_path.mkdir(parents=True, exist_ok=True)

    for src in file_list:
        src_path = src_root_path / src
        if src_path.exists():
            print(src_path, '->', dest_path / src_path.name)
            shutil.copy(src_path, dest_path / src_path.name)
        else:
            print(f"-------- File not found: {src_path}")

    return