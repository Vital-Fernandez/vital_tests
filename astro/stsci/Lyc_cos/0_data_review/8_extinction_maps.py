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
import dustmaps
from dustmaps.planck import PlanckQuery


def latex_sci(x, precision=2):
    s = f"{x:.{precision}e}"
    mantissa, exp = s.split("e")
    return f"${mantissa} \\times 10^{{{int(exp)}}}$"


# Extinction functions
# import dustmaps.planck
# dustmaps.planck.fetch()
planck = PlanckQuery()

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')

# Open density file
density_map_fname = project_folder/'LyC_leakers_COS/Column_density/MolecularHydrogen_DensityMap_AIT_-150_150.fits'
data =fits.getdata(density_map_fname, ext=0)
hdr = fits.getheader(density_map_fname, ext=0)
wcs = WCS(hdr)

ext_dict = {}
for obj, coord_deg in cfg_sample['aper_mean_coord'].items():
    coord = SkyCoord(ra=coord_deg[0] * u.deg, dec=coord_deg[1] * u.deg, frame='icrs')
    ebv = planck(coord)
    ext_dict[obj] = float(ebv)
    print(obj, ebv)

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
            label = f'{obj} ' + r'$(E(B-V) = $' + f'{ext_dict[obj]:0.3f}'

            ax.annotate(f'{label}', xy=(x, y), xytext=(20, 20), textcoords='offset pixels', color='white', fontsize=5,
                bbox=dict(boxstyle="round,pad=0.3", fc="black", alpha=0.6), arrowprops=dict(arrowstyle='->', color='white'))

    # Grid
    ax.grid(color='white', ls=':', alpha=0.5)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("H I Column Density  N(H)  [cm⁻²]", fontsize=10)

    # Axis formatting
    ax.coords[0].set_axislabel('Galactic Longitude (l)', fontsize=15)
    ax.coords[1].set_axislabel('Galactic Latitude (b)', fontsize=15)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.set_title("Milky Way color excess in log(N(H)) map", fontsize=18)

    # Display show
    # plt.show()
    plt.savefig(project_folder/'LyC_leakers_COS/Column_density/planck_ebv_map_over_logN(H).png')

# Save measurements
out_folder = project_folder/'LyC_leakers_COS/planck_ebv.toml'
lime.save_cfg(out_folder, ext_dict, section_name='planck_extinction_values')

# Linear scatter plot
col_density = cfg_sample["hydrogen_column_density_MW"]
extinction  = cfg_sample["planck_extinction_values"]

# Keys to exclude
exclude = {"IZw18_SE", "Haro11_B", "Haro11_C"}

# Keep only keys present in both sections and not excluded
common_keys = (set(col_density) & set(extinction)) - exclude

# Sort by column density (ascending)
sorted_keys = sorted(common_keys, key=lambda k: col_density[k])

x = np.array([col_density[k] for k in sorted_keys])
y = np.array([extinction[k]  for k in sorted_keys])

# --- Plot ---
fig, ax = plt.subplots(figsize=(9, 6))

ax.scatter(x, y, color="steelblue", s=70, zorder=3)

for k, xi, yi in zip(sorted_keys, x, y):
    ax.annotate(k, (xi, yi), textcoords="offset points",
                xytext=(5, 4), fontsize=7.5, color="0.3")

ax.set_xlabel(r"$\log\,N(\mathrm{H})_\mathrm{MW}$  [cm$^{-2}$]", fontsize=12)
ax.set_ylabel(r"Planck $E(B-V)$  [mag]", fontsize=12)
ax.set_title("MW H column density vs. Planck extinction", fontsize=13)
ax.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(project_folder / 'LyC_leakers_COS/Column_density/planck_ebv_versus_logN(H).png')
# plt.show()