from astropy.io import fits


def compare_header_values(file1, file2, ext=0):
    """
    Compare header values of two FITS files for a given HDU.

    Parameters:
        file1 (str): path to the first FITS file
        file2 (str): path to the second FITS file
        ext (int): HDU extension (default 0 = primary)

    Returns:
        dict: keys with different values and their values in each file
    """
    hdr1 = fits.getheader(file1, ext)
    hdr2 = fits.getheader(file2, ext)

    # Find keys that exist in both headers
    common_keys = set(hdr1.keys()) & set(hdr2.keys())

    # Keep only those with different values
    different_values = {key: (hdr1[key], hdr2[key])
                        for key in common_keys
                        if hdr1[key] != hdr2[key]}

    return different_values


# Example usage
diffs = compare_header_values("/home/vital/Astrodata/STScI/LyC_leakers_COS/hasp_multi/hst_cos_haro11-multi_g130m-lp04_cspec.fits",
                              "/home/vital/Astrodata/STScI/LyC_leakers_COS/hasp_multi/hst_cos_haro11-multi_g130m-lp02_cspec.fits",
                              ext=0)

print("Header keywords with different values:")
for key, (val1, val2) in diffs.items():
    print(f"{key}: file1={val1} | file2={val2}")