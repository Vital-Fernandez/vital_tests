import numpy as np
import pandas as pd
from astropy.io import fits
from pathlib import Path
import lime

# Load the new database
new_database = '/home/vital/PycharmProjects/lime/src/lime/resources/lines_database_v2.0.4.txt'
lime.lineDB.set_database(ref_bands=new_database, vacuum_waves=True)

# Specify the data
fname = '/home/vital/Dropbox/Astrophysics/Data/STScI_projects/NGC_628_MSA_nirspec/allslits_spectra_extcorr.fits'
bands_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/NGC_628_MSA_nirspec/bands')
output_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/NGC_628_MSA_nirspec/results')
cfg = lime.load_cfg('ngc628_cfg.toml')

# Convert RecArray to dataframe
data = fits.getdata(fname, ext=1)
data = np.asarray(data)
if not data.dtype.isnative:
    data = data.byteswap().view(data.dtype.newbyteorder('='))
df = pd.DataFrame.from_records(data)
df["source"] = df["source"].astype(str)
df["grating"] = df["grating"].astype(str)
df = df.loc[df.source == '100310']

# Loop through the objects (21 cosmic rays, 25 absorptions)
counter = 0
obj_list =  np.sort(df.source.unique())
for i, obj_name in enumerate(obj_list):
    slit_list = np.sort(df.loc[df.source == obj_name, 'slit'].unique())
    for j, slit in enumerate(slit_list):
        grating_list = np.sort(df.loc[(df.source == obj_name) & (df.slit == slit), 'grating'].unique())
        for k, grat in enumerate(grating_list):
            if counter >= 0:
                obj_ref = f'obj{obj_name}_slit{slit}_grat{grat}'
                print(f'{counter}), {obj_ref}')
                df_obj = df.loc[(df.source == obj_name) & (df.slit == slit) & (df.grating == grat)]
                bands_file = bands_folder/f'{obj_ref}.txt'

                # Create the spectrum
                spec = lime.Spectrum(input_wave=df_obj['wavelength'].to_numpy(),
                                     input_flux=df_obj['flux_sub'].to_numpy(),
                                     input_err=df_obj['unc_sub'].to_numpy(),
                                     units_wave='um',
                                     units_flux='FLAM',
                                     redshift=0.002192)

                # Unit conversion
                spec.unit_conversion('AA')

                # Components detection
                spec.infer.components()

                # Get object bands
                bands = spec.retrieve.lines_frame(automatic_grouping=True, fit_cfg=cfg, obj_cfg_prefix=obj_ref,
                                                  components=['emission', 'doublet-em'])
                # spec.plot.spectrum(bands=bands, log_scale=True, rest_frame=True, maximize=True, show_components=True)

                if not bands_file.is_file():
                    lime.save_frame(bands_file, bands)

                # Manual Review
                # spec.check.bands(bands_file, fit_cfg=cfg, obj_cfg_prefix=obj_ref)

                # Fit the lines
                spec.fit.frame(bands_file,  fit_cfg=cfg, obj_cfg_prefix=obj_ref)
                # spec.fit.bands('He1_10832A', bands_file,  fit_cfg=cfg, obj_cfg_prefix=obj_ref)
                # bands_all = spec.retrieve.lines_frame()
                # spec.plot.spectrum(bands=bands_all, rest_frame=True)
                fstem = f'obj-{obj_name}_slit-{slit}_grat-{grat}'
                spec.plot.grid()
                # spec.save_frame(fname=output_folder/f'{fstem}_measurements.txt')
                # spec.save_frame(fname=output_folder/f'{fstem}_measurements.csv')

            # Add a counter
            counter += 1