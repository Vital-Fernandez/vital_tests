import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zipfile
from pathlib import Path
import lime
from matplotlib import pyplot, rc_context
from scipy.ndimage import gaussian_filter1d


lime.theme.set_style(fig_cfg={'figure.figsize': (10, 10), 'figure.dpi': 400,
                              'xtick.color': 'none', 'ytick.color': 'none',  # Hides ticks
                              'xtick.labelsize': 0, 'ytick.labelsize': 0,  # Hides tick labels
                              'axes.labelsize': 0})

fig_conf = lime.theme.fig_defaults()

list_spectra = ['/home/vital/Downloads/Z004/SSP-CHA-stellar_Z0.004_logt5.00.dat',
                #'/home/vital/Downloads/Z004/SSP-CHA-stellar_Z0.004_logt7.00.dat',
                '/home/vital/Downloads/Z004/SSP-CHA-stellar_Z0.004_logt8.00.dat',
                #'/home/vital/Downloads/Z004/SSP-CHA-stellar_Z0.004_logt9.00.dat',
                '/home/vital/Downloads/Z004/SSP-CHA-stellar_Z0.004_logt10.00.dat',
                ]

# Access the Viridis colormap
cmap = plt.get_cmap('jet')
num_lines = len(list_spectra)
mag_color = [0.28, 0.7, 0.95]
with rc_context(fig_conf):

    fig, ax = plt.subplots()

    for i, address in enumerate(list_spectra):
        SSP_i = Path(address)
        wave, flux = np.loadtxt(SSP_i, unpack=True, skiprows=1)

        # Gaussian smoothing
        sigma = 40  # Standard deviation for Gaussian kernel
        flux = gaussian_filter1d(flux, sigma)

        idx0, idxf = np.searchsorted(wave, (3680, 7000))
        idx_Halpha = np.searchsorted(wave, (6000, 6100))
        flux_norm = flux[idx_Halpha[0]:idx_Halpha[1]].sum()

        ax.step(wave[idx0:idxf], flux[idx0:idxf]/flux_norm, color=cmap(mag_color[i]))

    ax.set_yscale('log')
    ax.axis('off')
    plt.savefig('/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/stellar_transparent.png',
                bbox_inches='tight',
                transparent=True)
    plt.savefig('/home/vital/Dropbox/Astrophysics/Tools/SpectralSynthesis/stellar_white.png', bbox_inches='tight',
                transparent=False)
    plt.show()


# root_zip = zipfile.Path(models_path)
#
# print(list(root_zip.iterdir()))
#
# directory = pathlib.Path("source_dir/")
#
# with zipfile.ZipFile(SSP_zip, 'r') as zip_ref:
#     for file in zip_ref.namelist():
#         print(file)
# print(SSP_zip.is_file())
# with zipfile.ZipFile(SSP_zip.as_posix(), 'r') as zip_ref:
#     directories = set()
#     for zip_info in zip_ref.infolist():
#         if zip_info.filename.endswith('/'):  # check if it is a directory
#             directories.add(zip_info.filename)

   # for file_path in directory.iterdir():
   #     archive.write(file_path, arcname=file_path.name)