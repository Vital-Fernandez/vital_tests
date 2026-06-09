import re
import lime
import pandas as pd
from pathlib import Path
from astro.stsci.tools import unpack_lines, absorption_spectrum, add_opacity_profile

import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def plot_logN_grid(df: pd.DataFrame) -> plt.Figure:
    """
    df : MultiIndex DataFrame with levels ["object", "line"].
         Must have columns: "particle", "logN", "logN_err".

    Color rules per line suffix:
        ends with "_o-MW"  → blue
        ends with "_Second component"   → red
        anything else      → orange
    """
    COLORS = {"Milky way": "#2196F3", "Second component": "#E53935", "First component": "#FB8C00"}
    MARKERS = {"Milky way": "o",       "Second component": "s",       "First component": "^"}

    def _category(line_label: str) -> str:
        if str(line_label).endswith("_o-MW"):
            return "Milky way"
        if str(line_label).endswith("_k-1"):
            return "Second component"
        return "First component"

    particles = sorted(df["particle"].unique())
    n = len(particles)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 4, nrows * 3.5),
        constrained_layout=True,
    )
    axes_flat = [axes] if n == 1 else list(
        axes.flat if hasattr(axes, "flat") else [axes]
    )

    objects = df.index.get_level_values("object").unique()
    obj_positions = {obj: i for i, obj in enumerate(objects)}

    for ax, particle in zip(axes_flat, particles):
        sub = df[df["particle"] == particle]

        for obj, grp in sub.groupby(level="object"):
            x = obj_positions[obj]
            lines = grp.index.get_level_values("line")

            for cat in ("Milky way", "Second component", "First component"):
                mask = [_category(l) == cat for l in lines]
                if not any(mask):
                    continue

                # keep unique (logN, logN_err) pairs for this category
                unique_rows = (
                    grp.loc[mask]
                    .drop_duplicates(subset=["logN", "logN_err"])
                )
                for _, row in unique_rows.iterrows():
                    ax.errorbar(
                        x, row["logN"],
                        yerr=row["logN_err"],
                        fmt=MARKERS[cat],
                        color=COLORS[cat],
                        capsize=3,
                        markersize=6,
                        elinewidth=1,
                        alpha=0.85,
                    )

        ax.set_title(particle, fontsize=11, pad=4)
        ax.set_xticks(range(len(objects)))
        ax.set_xticklabels(list(objects), rotation=80, ha="center", fontsize=8)
        ax.set_ylabel("log N", fontsize=9)
        ax.set_xlim(-0.6, len(objects) - 0.4)
        ax.grid(axis="y", linewidth=0.4, alpha=0.5)

    # hide unused axes
    for ax in axes_flat[n:]:
        ax.set_visible(False)

    # global legend
    legend_handles = [
        mlines.Line2D([], [], color=COLORS[cat], marker=MARKERS[cat],
                      linestyle="None", markersize=7, label=cat)
        for cat in ("Milky way", "Second component", "First component")
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        fontsize=10,
        frameon=True,
        bbox_to_anchor=(0.5, -0.02),
    )

    return fig

def plot_param_grid(
    df: pd.DataFrame,
    parameter: str = "b",
    include_ow_mw: bool = True,
    ncols: int | None = None,
) -> plt.Figure:
    """
    One subplot per unique object; x-axis = all unique particles in df.

    Parameters
    ----------
    df : MultiIndex DataFrame with levels ["object", "line"].
         Must have columns: "particle", , _err.
    parameter : column to plot on the y-axis (default "b").
                The error column is assumed to be f"{parameter}_err".
    include_ow_mw : if False, rows whose line label ends with "_o-MW"
                   are excluded from the plot and the legend.
    ncols : fixed number of grid columns; if None, uses ceil(sqrt(n)).
    """
    PALETTE = {
        "Milky way":    {"color": "#2196F3", "marker": "o"},
        "Second component":     {"color": "#E53935", "marker": "s"},
        "First component": {"color": "#FB8C00", "marker": "^"},
    }
    err_col = f"{parameter}_err"

    def _category(line_label: str) -> str:
        if str(line_label).endswith("_o-MW"):
            return "o-MW"
        if str(line_label).endswith("_k-1"):
            return "Second component"
        return "First component"

    # all particles present in the full dataframe (fixed x-axis)
    all_particles = sorted(df["particle"].unique())
    part_positions = {p: i for i, p in enumerate(all_particles)}

    objects = df.index.get_level_values("object").unique()
    n = len(objects)
    _ncols = ncols if ncols is not None else math.ceil(math.sqrt(n))
    _nrows = math.ceil(n / _ncols)

    active_cats = ["Second component", "First component"] + (["o-MW"] if include_ow_mw else [])

    fig, axes = plt.subplots(
        _nrows, _ncols,
        figsize=(_ncols * 4, _nrows * 3.5),
        constrained_layout=True,
    )
    axes_flat = (
        [axes] if n == 1
        else list(axes.flat if hasattr(axes, "flat") else [axes])
    )

    for ax, obj in zip(axes_flat, objects):
        sub = df.xs(obj, level="object")
        if not include_ow_mw:
            lines = sub.index.get_level_values("line")
            sub = sub[[not str(l).endswith("_o-MW") for l in lines]]

        for particle, grp in sub.groupby("particle"):
            x = part_positions[particle]
            lines = grp.index.get_level_values("line")

            for cat in active_cats:
                mask = [_category(l) == cat for l in lines]
                if not any(mask):
                    continue

                unique_rows = (
                    grp.loc[mask]
                    .drop_duplicates(subset=[parameter, err_col])
                )
                for _, row in unique_rows.iterrows():
                    ax.errorbar(
                        x, row[parameter],
                        yerr=row[err_col],
                        fmt=PALETTE[cat]["marker"],
                        color=PALETTE[cat]["color"],
                        capsize=3,
                        markersize=6,
                        elinewidth=1,
                        alpha=0.85,
                    )

        ax.set_title(obj, fontsize=11, pad=4)
        ax.set_xticks(range(len(all_particles)))
        ax.set_xticklabels(all_particles, rotation=80, ha="right", fontsize=8)
        ax.set_ylabel(parameter, fontsize=9)
        ax.set_xlim(-0.6, len(all_particles) - 0.4)
        ax.grid(axis="y", linewidth=0.4, alpha=0.5)

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    legend_handles = [
        mlines.Line2D(
            [], [],
            color=PALETTE[cat]["color"],
            marker=PALETTE[cat]["marker"],
            linestyle="None", markersize=7, label=cat,
        )
        for cat in active_cats
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(active_cats),
        fontsize=10,
        frameon=True,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(parameter, fontsize=13, fontweight="normal", y=1.01)

    return fig


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
metals_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'
lyman_alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
lsf_obj_folder = project_folder/'LyC_leakers_COS'/'lsf_objects'
results_folder =  project_folder/'LyC_leakers_COS'/'line_measurements'

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']
lime_voigtfit_conv = cfg_sample['LiMe_voigtfit_conversion']
lyAlpha_data = lime.load_cfg(project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'/'LyAlpha_results.toml')
LyC_list = list(cfg_sample['LyC_target_lines'].keys())

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


version = 'v1'
data_df = lime.load_frame(results_folder/f'LyC_voigtprofiles_{version}.txt', levels=['object', 'line'])

fig = plot_logN_grid(data_df)
fig.savefig(results_folder/f"species_logN_grid_{version}.pdf", bbox_inches="tight")

fig = plot_param_grid(data_df, parameter='v_r', include_ow_mw=False)
fig.savefig(results_folder/f"object_vel_grid_{version}.pdf", bbox_inches="tight")

fig = plot_param_grid(data_df, parameter='b', include_ow_mw=False)
fig.savefig(results_folder/f"object_b_grid_{version}.pdf", bbox_inches="tight")

fig = plot_param_grid(data_df, parameter='z_line', include_ow_mw=False)
fig.savefig(results_folder/f"object_z_grid_{version}.pdf", bbox_inches="tight")

#
# fig = plot_param_grid(data_df, parameter='v_r', include_ow_mw=False)

# plt.show()

# survey_dict = {}
# for i, obj in enumerate(target_list):
#
#     if i >= 0:
#
#         # Object folders the files
#         print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
#         input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
#         output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
#         opacityLymanAlpha_df_path = lyman_alpha_folder/f'{obj}_LyAlpha_lines_frame.txt'
#         opacity_df_path = metals_folder/f'{obj}_metals_lines_frame.txt'
#         voigfitz_reg = metals_folder/f'{obj}_metals_best_fit.reg'
#         obj_lsf_file = f"{lsf_obj_folder}/{obj}_hasp_lsf.txt"
#
#         HI_cfg = lyAlpha_data[f'{obj}_LyAlpha_results']
#         metals_cfg = cfg_sample['voigtfit_metals_params'][obj]
#         alpha_cfg = cfg_sample['voigtfit_lyAlpha_params'][obj]
#
#         # Index the files
#         idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')
#
#         # Load the files
#         pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
#         spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj],
#                                        **metals_cfg['from_file'])
#         # Rebinned
#         spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)
#
#         # Divide by the HI profile
#         norm_Lyman_alpha = absorption_spectrum(opacityLymanAlpha_df_path, spec)
#         norm_Lyman_alpha[norm_Lyman_alpha < 0.01] = 1
#         spec.flux = spec.flux /norm_Lyman_alpha
#         spec.err_flux = spec.err_flux /norm_Lyman_alpha
#
#         # Normalization
#         spec = spec.retrieve.normalization(**metals_cfg['normalization'])
#
#         # Spectra masking
#         if 'Metals_retrieve' in metals_cfg:
#             spec = spec.retrieve.spectrum(**metals_cfg['Metals_retrieve'])
#
#         # Save the spectrum
#         spec.save_spectrum(fname=metals_folder/f'{obj}_metals_spec.txt')
#
#         # Sort the line selection line lists
#         list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
#                                                  mask_groups=['airglow', 'MW', 'nebular', 'stellar'])


