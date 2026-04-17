from pathlib import Path
import numpy as np
import lime
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib import pyplot as plt, rc_context
from astropy.visualization import LogStretch


def latex_sci(x, precision=2):
    s = f"{x:.{precision}e}"
    mantissa, exp = s.split("e")
    return f"${mantissa} \\times 10^{{{int(exp)}}}$"

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')
den_dict = {}

# Open density file
density_map_fname = project_folder/'LyC_leakers_COS/Column_density/MolecularHydrogen_DensityMap_AIT_-150_150.fits'
data =fits.getdata(density_map_fname, ext=0)
hdr = fits.getheader(density_map_fname, ext=0)
wcs = WCS(hdr)

lime.theme.set_style('dark')
with rc_context(lime.theme.fig_defaults()):

    # Create teh figure
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection=wcs)

    # Image plot
    norm = ImageNormalize(data, stretch=LogStretch())
    im = ax.imshow(data, origin='lower', cmap='inferno', norm=norm)

    # Loop through the objects and plot them on the plot
    for obj, coord_deg in cfg_sample['aper_mean_coord'].items():

        if obj not in ['Haro11_C', 'Haro11_C', 'IZw18_SE']:
            coord = SkyCoord(ra=coord_deg[0] * u.deg, dec=coord_deg[1] * u.deg, frame='icrs')
            x, y = wcs.world_to_pixel(coord)
            x, y = int(np.round(x)), int(np.round(y))
            label = f'{obj} ' + r'$log(N_{HI}) = $' + latex_sci(data[y, x], 2)
            den_dict[obj] = float(data[y, x])

            ax.annotate(f'{label}', xy=(x, y), xytext=(20, 20), textcoords='offset pixels', color='white', fontsize=5,
                bbox=dict(boxstyle="round,pad=0.3", fc="black", alpha=0.6), arrowprops=dict(arrowstyle='->', color='white'))

    # Grid
    ax.grid(color='white', ls=':', alpha=0.5)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("H I Column Density  N(H)  [cm⁻²]", fontsize=10)

    # Axis formatting
    ax.coords[0].set_axislabel('Galactic Longitude (l)', fontsize=20)
    ax.coords[1].set_axislabel('Galactic Latitude (b)', fontsize=20)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.set_title("Milky Way HI Column Density Map", fontsize=25)

    # Display show
    # plt.show()
    plt.savefig(project_folder/'LyC_leakers_COS/Column_density/HI_density_map_dark.png')

    # Save measurements
    out_folder = project_folder/'LyC_leakers_COS/Column_density/LyC_leakers_hydrogen_column_density.toml'
    lime.save_cfg(out_folder, den_dict, section_name='H_column_density')

