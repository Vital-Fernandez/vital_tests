import re
import lime
import numpy as np
from pathlib import Path
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
from astro.stsci.tools import get_masked_bands, CORNER_CFG, plot_chains, plot_sesamme_samples, plot_ssp_params





# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
ssp_wave_arr = np.loadtxt('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/pyStarburst99_wave_arr.txt')
stellar_folder = project_folder / 'LyC_leakers_COS' / 'SSPs'
SESAMME_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/SESAMME_fittings')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C']


# Loop through the targets and generate the HASP files
results_dict = {}
for i, obj in enumerate(target_list):

    # Declare inputs
    print(f'{i}) Galaxy {obj}')
    input_pname = stellar_folder / f'{obj}_pyStarburst99_spectrum_Lsun.txt'
    spec = lime.Spectrum.from_file(input_pname, instrument='text')

    wave_arr = spec.wave_rest.data
    flux_arr = spec.flux.data * spec.norm_flux
    err_arr = spec.err_flux.data * spec.norm_flux
    mask = ~spec.wave_rest.mask

    # Prepare stellar inputs
    runname = f'v3_default_pystarburst99'
    fname_sesamme = SESAMME_folder / runname / f'{obj}_SESAMME.h5'
    fname_sesamme.parent.mkdir(parents=True, exist_ok=True)

    fname = SESAMME_folder / runname/f'{obj}_{runname}_stats.txt'
    df = pd.read_csv(fname, sep="\t")
    df = df.set_index("parameter")
    results_dict[obj] = df

plot_ssp_params(results_dict, target_list, metallicity_grid=[0.0, 1e-05, 0.0004, 0.002, 0.006, 0.014, 0.02],\
                fname=f'SESAMME_{runname}_results_grid.png')

    #     '/home/vital/Astrodata/STScI/LyC_leakers_COS/SESAMME_fittings/v3_default_pystarburst99/Haro11_A_v3_default_pystarburst99_stats.txt'
    # df = pd.read_csv(path, sep="\t")
    # return df.set_index("parameter").to_dict("index")
#
# """
# Plot pyStarburst99 stellar-continuum fit results across the galaxy sample.
#
# Reads one *_stats.txt file per galaxy (tab-separated: parameter, 16th, 50th,
# 84th percentile), of the form produced by the SSP fitting pipeline, e.g.:
#
#     1789139204346_Haro11_C_v3_default_pystarburst99_stats.txt
#
# Each filename is <timestamp>_<galaxy>_<version>_<tag>_pystarburst99_stats.txt.
# The leading timestamp changes per run, so files are located with a glob
# that matches on galaxy name + version/tag and ignores the timestamp. If a
# galaxy has more than one matching file (re-runs), the most recently
# modified one is used by default.
#
# Edit STATS_FOLDER, VERSION, TAG below to match your directory layout before
# running -- these three are the only path-dependent settings.
# """
#
# from pathlib import Path
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # ---------------------------------------------------------------- settings
# STATS_FOLDER = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/SESAMME_fittings/v3_default_pystarburst99/')   # <-- adjust to where *_stats.txt live
# VERSION = "v3"
# TAG = "default"
# SUFFIX = "pystarburst99_stats.txt"
#
# OUT = "./"
#
# # same target order used in the distance/extinction plots, for consistency
# # across the paper's figures
# TARGET_LIST = [
#     'Haro11_A', 'Haro11_B', 'Haro11_C',
#     'SBS0335052', 'SBS1159545', 'Haro2', 'Pox186', 'UM461', 'MRK1450',
#     'He2-10', 'IZw18', 'NGC4861', 'IZw18_SE', 'NGC1705', 'SBS1415437',
#     'UGCA281', 'UGC4483', 'VIIZw403', 'NGC2366',
# ]
#
# # parameter row label -> (y-axis label, log-scale?)
# PARAMS = {
#     "log(age/yr)":        (r"$\log(\mathrm{age}/\mathrm{yr})$", False),
#     "log(Z/Z$_{\\odot}$)": (r"$\log(Z/Z_\odot)$", False),
#     "E(B-V)":              (r"internal $E(B-V)$  [mag]", False),
#     "log(A)":               (r"$\log(A)$", False),
# }
#
# C_FIT = "#2b6cb0"
#
#
# # ------------------------------------------------------------------ loading
# def find_stats_file(obj, folder=STATS_FOLDER, version=VERSION, tag=TAG,
#                     suffix=SUFFIX):
#     """Locate a galaxy's stats file regardless of its leading timestamp."""
#     pattern = f"{obj}_{version}_{tag}_{suffix}"
#     matches = sorted(folder.glob(pattern), key=lambda p: p.stat().st_mtime)
#     if not matches:
#         return None
#     if len(matches) > 1:
#         print(f"  note: {len(matches)} matches for {obj}, using most recent "
#               f"({matches[-1].name})")
#     return matches[-1]
#
#
# def read_stats(path):
#     """parameter -> {16th, 50th, 84th}"""
#     df = pd.read_csv(path, sep="\t")
#     return df.set_index("parameter").to_dict("index")
#
#
# def build_results(names=TARGET_LIST, **kwargs):
#     results = {}
#     for obj in names:
#         path = find_stats_file(obj, **kwargs)
#         if path is None:
#             print(f"  no stats file for {obj}")
#             results[obj] = None
#             continue
#         results[obj] = read_stats(path)
#     return results
#
#
# # ------------------------------------------------------------------ plotting
# def plot_ssp_params(results, names=TARGET_LIST, params=PARAMS, out_prefix="ssp_fit"):
#     xi = {n: i for i, n in enumerate(names)}
#     x = np.arange(len(names))
#
#     n_p = len(params)
#     ncols = 2
#     nrows = int(np.ceil(n_p / ncols))
#     fig, axes = plt.subplots(nrows, ncols, figsize=(13.5, 4.4 * nrows), squeeze=False)
#     axes = axes.flatten()
#
#     for ax, (pname, (ylabel, logscale)) in zip(axes, params.items()):
#         vals, err_lo, err_hi, xs = [], [], [], []
#         for n in names:
#             row = results.get(n)
#             if row is None or pname not in row:
#                 continue
#             p16, p50, p84 = row[pname]["16th"], row[pname]["50th"], row[pname]["84th"]
#             xs.append(xi[n])
#             vals.append(p50)
#             err_lo.append(max(p50 - p16, 0.0))
#             err_hi.append(max(p84 - p50, 0.0))
#
#         ax.errorbar(xs, vals, yerr=[err_lo, err_hi], fmt="o", ms=7,
#                     color=C_FIT, mec="white", mew=1.1, ecolor=C_FIT,
#                     elinewidth=1.6, capsize=3, zorder=4)
#
#         for n in names:
#             if results.get(n) is None or pname not in (results.get(n) or {}):
#                 y0, y1 = ax.get_ylim()
#                 ax.annotate("no fit", xy=(xi[n], (y0 + y1) / 2),
#                             ha="center", va="center", fontsize=7.5, color="0.5",
#                             bbox=dict(boxstyle="round,pad=0.25", fc="white",
#                                       ec="0.7", lw=0.8, alpha=0.9))
#
#         if logscale:
#             ax.set_yscale("log")
#         ax.set_xlim(-0.7, len(names) - 0.3)
#         ax.set_xticks(x)
#         ax.set_xticklabels(names, rotation=55, ha="right", fontsize=8.5)
#         ax.set_ylabel(ylabel, fontsize=11)
#         ax.grid(axis="y", ls=":", lw=0.7, color="0.8", zorder=0)
#         for i in x[1::2]:
#             ax.axvspan(i - 0.5, i + 0.5, color="0.5", alpha=0.05, zorder=0, lw=0)
#         ax.set_title(pname, fontsize=11, loc="left")
#
#     for ax in axes[n_p:]:
#         ax.axis("off")
#
#     fig.suptitle("pyStarburst99 fit results across the sample", fontsize=13)
#     fig.tight_layout()
#     fig.savefig(f"{OUT}{out_prefix}.png", dpi=200)
#     fig.savefig(f"{OUT}{out_prefix}.pdf")
#     print(f"saved {OUT}{out_prefix}.png")
#
#
# def print_table(results, names=TARGET_LIST, params=PARAMS):
#     rows = {}
#     for n in names:
#         row = results.get(n)
#         rows[n] = {
#             pname: (row[pname]["50th"] if row and pname in row else np.nan)
#             for pname in params
#         }
#     table = pd.DataFrame(rows).T
#     table.index.name = "galaxy"
#     print("\nMedian values")
#     print(table.round(3).to_string(na_rep="--"))
#     return table
#
#
# if __name__ == "__main__":
#     results = build_results()
#     print_table(results)
#     plot_ssp_params(results)