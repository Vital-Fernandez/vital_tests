from pathlib import Path
import lime
import pandas as pd
from astro.stsci.plots import cos_image_plotter
from astropy.io import fits
from astropy.wcs import WCS

pd.set_option('display.max_colwidth', None)


# Data location
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
obs_folder = Path('/home/vital/Astrodata/STScI')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv',
                            levels=['sample', 'id', 'offset_id', 'state'])
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

idcs_leda = sample_df['object'] == 'UM461'
obj_df = sample_df.loc[idcs_leda]

idcs = obj_df.index.get_level_values('state').isin(['flt', 'drz', 'mos'])

print(obj_df.loc[idcs, ['PID', 'filepath']])

# Show all the images available
for i, idx_image in enumerate(obj_df.loc[idcs].index):
    if i > 0:
        image_path = obs_folder / obj_df.loc[idx_image].filepath
        print(f'{i}) {image_path.name}')

        # Get figure data
        hdr = fits.getheader(image_path, extname='SCI')
        hdr0 = fits.getheader(image_path, ext=0)

        try:
            imdata = fits.getdata(image_path, extname='SCI')

            try:
                wcs = WCS(hdr)
            except:
                wcs = WCS(fits.getheader(image_path, ext=2))

            print('  ', obj_df.loc[idx_image].name[1], obj_df.loc[idx_image].filepath, obj_df.loc[idx_image].grating,
                  imdata.shape)

            # Plot the image
            cos_image_plotter(imdata, wcs, 'UM461')
        except:
            print('NO DATA')

