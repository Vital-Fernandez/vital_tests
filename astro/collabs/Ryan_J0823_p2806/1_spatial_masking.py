from pathlib import Path
from astropy.io import fits
import lime


def combine_fits(file1, file2, output_file, rename1=None, rename2=None, overwrite=True):

    """
    Combine HDUs from two FITS files into a single output file.

    Parameters
    ----------
    file1, file2 : str or Path
        Input FITS filenames.
    output_file : str or Path
        Output FITS filename.
    rename1, rename2 : dict or None
        Optional dictionaries to rename extensions from file1 and file2.
        Keys can be:
          - integer HDU index (0, 1, 2, ...)
          - current EXTNAME string
        Values are the new EXTNAME strings.

        Example:
            rename1 = {
                1: "SCI_A",           # rename HDU 1 of file1
                "ERR": "ERR_A"        # rename any HDU named "ERR" in file1
            }
    overwrite : bool
        If True, overwrite existing output_file.

    Notes
    -----
    - PRIMARY HDU comes from file1.
    - All non-primary HDUs from file1 and file2 are appended, in that order.
    """

    file1, file2, output_file = Path(file1), Path(file2), Path(output_file)

    rename1 = rename1 or {}
    rename2 = rename2 or {}

    def maybe_rename(hdu, idx, rename_dict):

        """Rename EXTNAME of an HDU in-place if a rule is found."""
        if idx == 0:
            return

        old_name = hdu.header.get("EXTNAME", None)

        new_name = None

        # Priority: index-based rename, then name-based rename
        if idx in rename_dict:
            new_name = rename_dict[idx]
        elif old_name in rename_dict:
            new_name = rename_dict[old_name]

        if new_name is not None:
            hdu.header["EXTNAME"] = (new_name, "Renamed by combine_fits")

    # Open both files in context managers
    with fits.open(file1, mode="readonly") as hdul1, \
         fits.open(file2, mode="readonly") as hdul2:

        # Start new HDUList with PRIMARY from file1
        new_hdus = [hdul1[0].copy()]

        # Add non-primary HDUs from file1
        for i in range(1, len(hdul1)):
            hdu = hdul1[i].copy()
            maybe_rename(hdu, i, rename1)
            new_hdus.append(hdu)

        # Add non-primary HDUs from file2
        for j in range(1, len(hdul2)):
            hdu = hdul2[j].copy()
            maybe_rename(hdu, j, rename2)
            new_hdus.append(hdu)

        # Build output HDUList and write
        hdul_out = fits.HDUList(new_hdus)
        hdul_out.writeto(output_file, overwrite=overwrite)

    return



fname = '/home/vital/Astrodata/J0823_p2806_v2.fits'
spatial_O2nebular_mask_file = './J0823_p2806_SN_O2_3726A_mask.fits'
spatial_O3auroral_mask_file = './J0823_p2806_SN_O3_4363A_mask.fits'
combine_mask_file = './J0823_p2806_SN_mask.fits'
spatial_He2_mask_file = './J0823_p2806_SN_He2_4686A_mask.fits'

# Extra maps
peak_wave_file = './maps/J0823_peak_wave.fits'
wave_map = fits.getdata(peak_wave_file, 'O3_5007A')
mask_name = 'SOUTH'
mask_file = f'./J0823_{mask_name}_mask.fits'
wave_mask = wave_map < 5241
wave_mask = wave_mask.astype(int)

# Empty primary HDU (required by FITS standard)
primary_hdu = fits.PrimaryHDU()
hdr = fits.Header({'PARAM': 'wave'})
image_hdu = fits.ImageHDU(data=wave_mask, header=hdr, name='mask_name')
hdul = fits.HDUList([primary_hdu, image_hdu])
hdul.writeto(mask_file, overwrite=True)

# Open the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726)
cube.check.cube('O2_3726A', masks_file=mask_file, rest_frame=True)

# cube.check.cube('O2_3726A', line_fg='O3_5007A', cont_pctls_fg=[98.5, 99.5], rest_frame=True)

# # SN masks
# cube.spatial_masking('O2_3726A', param='SN_line', contour_pctls=[64], fname=spatial_O2nebular_mask_file)
# cube.spatial_masking('O3_4363A', param='SN_line', contour_pctls=[98, 99], fname=spatial_O3auroral_mask_file)

# cube.check.cube('O3_5007A', masks_file=spatial_O3auroral_mask_file, min_pctl_bg=60, rest_frame=True, masks_cmap='viridis')
# cube.plot.cube('O2_3726A', masks_file=spatial_O3auroral_mask_file)





# # Combine the files
# mask1_rename = {'MASK_0': 'FG_O2_60p'}
# mask2_rename = {'MASK_0': 'Temp_O3_99p', 'MASK_1': 'Temp_O3_98p'}
# combine_fits(spatial_O2nebular_mask_file, spatial_O3auroral_mask_file, combine_mask_file, mask1_rename, mask2_rename)
# cube.plot.cube('O2_3726A', masks_file=combine_mask_file)
#
# cube.spatial_masking('H1_4861A', param='SN_cont', contour_pctls=[55, 60, 97.5], fname=spatial_O2nebular_mask_file)
# cube.check.cube('O2_3726A', masks_file=spatial_O2nebular_mask_file, min_pctl_bg=55, rest_frame=True, masks_cmap='viridis')