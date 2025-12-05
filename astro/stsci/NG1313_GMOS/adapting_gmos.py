from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import lime

fname = '/home/vital/Downloads/S20240221S0084.fits'

fits.info(fname)

# Open the FITS file
with fits.open(fname) as hdul:
    # Print a summary of the HDUs in the file
    hdul.info()

    # Access the primary HDU (usually the first one)
    primary_header = hdul[0].header
    data = hdul[1].data
    hdr = hdul[1].header

    print("\n--- Primary Header ---")
    print(repr(primary_header))

    # Access the science data. For GMOS, science data is often in an extension.
    # The extension number can vary, but [1] is a common starting point for imaging.
    try:
        science_data = hdul['SCI'].data
        science_header = hdul['SCI'].header
        print("\n--- Science Extension Header ---")
        print(repr(science_header))

        # Display the image data
        plt.imshow(science_data, cmap='viridis')
        plt.title('GMOS Science Image Data')
        plt.xlabel('X-pixels')
        plt.ylabel('Y-pixels')
        plt.colorbar(label='Pixel Value')
        plt.show()

        # You can also access data by its extension name, if it has one (e.g., 'SCI').
    except KeyError:
        print("\nNo 'SCI' extension found. Checking for the first image extension.")
        # If 'SCI' is not found, loop through and look for the first image data.
        for hdu in hdul:
            if hdu.data is not None and isinstance(hdu.data, np.ndarray):
                print(f"Found image data in HDU {hdul.index(hdu)}")
                plt.imshow(hdu.data, cmap='viridis')
                plt.title(f'GMOS Image Data (HDU {hdul.index(hdu)})')
                plt.colorbar()
                plt.show()
                break
