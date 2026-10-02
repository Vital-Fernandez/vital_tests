import numpy as np
import pandas as pd
from pathlib import Path
from specsy.models.ssp import StellarBinaries

source_folder = Path('/home/vital/Astrodata/VMS_BPASS/bursts_MP25')
path_list = [f.name for f in source_folder.iterdir() if f.is_file()]

vms_ssp = StellarBinaries.from_source('VMS_BPASS', source_folder, file_list=path_list, load_spectra=True)
# vms_ssp.plot_age_metallicity()
print(np.power(10, vms_ssp.ages))
print(vms_ssp)


high_res = {"figure.dpi": 300,
    "figure.figsize": [11, 6],
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}

ref = pd.read_csv('./martins25_fig11_digitized.csv')
vms_ssp.plot_ionizing_hardness(ages_myr=[1, 2], libraries=['VMS_gr', 'VMS_grzsc'], z_sun=0.01, ylim=(1e-9, 1.5),
                               energies=(13.6, 100.01),
                               lib_labels={'VMS_gr': r'$\dot{M}$ loss not scaled with $Z$',
                                           'VMS_grzsc': r'$\dot{M}$ loss scaled with $Z$'},
                               fig_cfg=high_res, fname=None, #f'/home/vital/IonizationPhotons_VMS-Martins_AgeMetallicity_v4.png',
                               # ref_frame=ref, ref_label='Martins+25 Fig. 11'
                               )