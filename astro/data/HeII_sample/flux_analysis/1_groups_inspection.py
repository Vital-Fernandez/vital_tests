import pathlib
import numpy as np
import pandas as pd
import src.specsiser as sr

from matplotlib import pyplot as plt, rcParams


class KeyPushPlots(object):

    background = np.array((43, 43, 43)) / 255.0
    foreground = np.array((179, 199, 216)) / 255.0
    red = np.array((43, 43, 43)) / 255.0
    yellow = np.array((191, 144, 0)) / 255.0

    STANDARD_PLOT = {'figure.figsize': (14, 7),
                     'axes.titlesize': 14,
                     'axes.labelsize': 14,
                     'legend.fontsize': 12,
                     'xtick.labelsize': 12,
                     'ytick.labelsize': 12,
                     'figure.facecolor': background,
                     'axes.facecolor': background,
                     'axes.edgecolor': foreground,
                     'axes.labelcolor': foreground,
                     'axes.titlecolor': foreground,
                     'xtick.color': foreground,
                     'ytick.color': foreground
                     }
    STANDARD_AXES = {'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'}

    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    def __init__(self, input_list, dataframe_address):
        # self.fig = plt.figure()
        # self.ax = self.fig.gca()
        self.fig, self.ax = plt.subplots()
        self.bound_keys, self.bound_cid = [], {}
        self.catalogue_DF = input_list
        self.current_obj = None
        self.output_address = dataframe_address

    def add_step_through(self, gen, key):
        key = key[0]# make a single char
        if key in self.bound_keys:
            raise RuntimeError("key %s already bound"%key)

        idx_data, objName, wave, flux = next(gen)
        self.idx, self.current_obj = idx_data, objName

        self.ax.plot(wave, flux)
        self.ax.axvline(4686, linestyle='--', linewidth=0.5, color=self.foreground)
        self.ax.set_yscale('log')

        self.STANDARD_AXES['title'] = self.current_obj
        self.ax.update(self.STANDARD_AXES)

        self.fig.canvas.draw()
        self.bound_keys.append(key)

        def ontype(event):
            if event.key in ('enter', 'down'):
                if event.key == 'enter':
                    self.catalogue_DF.loc[self.current_obj, 'valid'] = True
                    print(f'- Saving {self.current_obj} ({self.idx})')
                if event.key == 'down':
                    self.catalogue_DF.loc[self.current_obj, 'valid'] = False
                    print(f'- Excluding {self.current_obj} ({self.idx})')
                self.store_dataframe()

                try:
                    plt.cla()

                    idx_data, objName, wave, flux = next(gen)
                    self.idx, self.current_obj = idx_data, objName

                    self.ax.plot(wave, flux)
                    self.ax.axvline(4686, linestyle='--', linewidth=0.5, color=self.foreground)
                    self.ax.set_yscale('log')

                    self.STANDARD_AXES['title'] = self.current_obj
                    self.ax.update(self.STANDARD_AXES)

                    self.fig.canvas.draw()

                except StopIteration:
                    self.fig.canvas.mpl_disconnect(self.bound_cid[key])
                    del self.bound_cid[key]
                    self.bound_keys.remove(key)

        self.bound_cid[key] = self.fig.canvas.mpl_connect('key_press_event', ontype)

    def store_dataframe(self):

        with open(self.output_address, 'wb') as output_file:
            string_DF = self.catalogue_DF.to_string()
            output_file.write(string_DF.encode('UTF-8'))

        return

def spec_generator(folder_data, obj_list, flux_norm, idx_start=0):

    for idx_obj, obj_name in enumerate(obj_list[idx_start:]):
        wave, data, hdrs = sr.import_fits_data(folder_data / f'{obj_name}.fits', instrument='SDSS')
        flux = data['flux'] * flux_norm
        z_i = hdrs[1]["z"][0]
        lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux, crop_waves=(1+z_i) * np.array([4685-100, 5100]))
        yield idx_obj, obj_name, lm.wave, lm.flux


conf_file_address = '../sampleHeII.ini'
obsData = sr.loadConfData(conf_file_address)

fits_folder = pathlib.Path(obsData['data_location']['fits_folder'])
data_folder = pathlib.Path(obsData['data_location']['treatment_folder'])
catalogueDataframe = data_folder/f'AVO_catalogue_dataframe.txt'
normFlux = obsData['sample_data']['norm_flux']

# Load the list of objects
logDF = sr.lineslogFile_to_DF(catalogueDataframe)
# logDF['valid'] = False
# if 'valid' not in logDF:
#     logDF['valid'] = False
print('empezamos', (logDF['valid']==True).sum())

# for obj in ['1734-53034-56', '859-52317-170', '616-52442-101', '525-52295-193', '2956-54525-19', '2173-53874-491']:
#     logDF.loc[obj, 'valid'] = True
#
# with open(catalogueDataframe, 'wb') as output_file:
#     string_DF = logDF.to_string()
#     output_file.write(string_DF.encode('UTF-8'))

first_obj = 0
list_groups = np.arange(8)
for idx_group in list_groups:

    idcs_obj = logDF['intensity_Group'] == idx_group
    objList = logDF.loc[idcs_obj].index.values

    print(f'Inspecting group {idx_group} with {np.sum(idcs_obj)} objects')
    # spec_group = spec_generator(fits_folder, objList, normFlux, idx_start=first_obj)
    #
    # pta = KeyPushPlots(logDF, catalogueDataframe)
    # pta.add_step_through(spec_group, 'enter')
    # plt.show()
# print('empezamos', (logDF['valid']==True).sum())
#
# for idx_group in [6]:
#
#     idcs_obj = logDF['intensity_Group'] == idx_group
#     subGroupDF = logDF.loc[idcs_obj]
#     subGroupDF.sort_values(by='He2_4685A', ascending=False, inplace=True)
#     objList = subGroupDF.index.values
#
#     print(f'Inspecting group {idx_group} with {np.sum(idcs_obj)} objects')
#     spec_group = spec_generator(fits_folder, objList, normFlux)
#
#     pta = KeyPushPlots(logDF, catalogueDataframe)
#     pta.add_step_through(spec_group, 'enter')
#     plt.show()


# pa incluir 1734-53034-56, 859-52317-170,616-52442-101, 525-52295-193, 2956-54525-19
    # print(f'Inspecting group {idx_group} with {np.sum(idcs_obj)} objects')
    # for i, obj in enumerate(objList):
    #
    #     spec_address = fits_folder/f'{obj}.fits'
    #     wave, data, hdrs = sr.import_fits_data(spec_address, instrument='SDSS')
    #     flux = data['flux'] * normFlux
    #     z_i = hdrs[1]["z"][0]
    #
    #     lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux)
    #     lm.plot_spectrum_components(specLabel=f'Galaxy {obj}',
    #                                 axConf={'ylabel': r'Flux $(10^{17}\,erg\,cm^{-2} s^{-1} \AA^{-1})$'})