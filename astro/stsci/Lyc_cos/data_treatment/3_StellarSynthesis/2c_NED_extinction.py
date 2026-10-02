"""
Compare Planck E(B-V) against NED's Galactic extinction for the sample.

The user picks R_V at the top. Both catalogues are first reduced to a colour
excess E(B-V) using the *native* A_V/E(B-V) coefficient of the curve NED used
(validated against the Landolt B/V ratios in the NED table itself, see below),
and only then converted to A_V with the selected R_V. That keeps the
Planck/NED ratio independent of R_V.

Native coefficients recovered from the data:
  SF11  (2011ApJ...737..103S)  A_B/A_V = 1.3198  ->  A_V/E(B-V) = 3.102
  SFD98 (1998ApJ...500..525S)  A_B/A_V = 1.3002  ->  A_V/E(B-V) = 3.315
matching the published Landolt coefficients (3.626/2.742 and 4.315/3.315).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import lime
# ---------------------------------------------------------------- user choice
R_V = 3.1          # <-- select R_V here; rescales both catalogues identically
SHOW_SFD98 = True  # also show NED's legacy SFD98 column

UP = "./"
# OUT = "/mnt/user-data/outputs/"

# native A_V/E(B-V) of each NED extinction column
RV_SF11 = 3.102
RV_SFD98 = 3.315

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# ------------------------------------------------------------- Planck E(B-V)
planck_ebv = cfg_sample['planck_extinction']
# planck_ebv = {
#     'Haro11_A': 0.016438599675893784, 'Haro11_B': 0.016438599675893784,
#     'Haro11_C': 0.016438599675893784, 'SBS0335052': 0.05214526131749153,
#     'IZw18': 0.04414275661110878,     'IZw18_SE': 0.04414275661110878,
#     'Haro2': 0.018410447984933853,    'SBS1159545': 0.017896920442581177,
#     'NGC1705': 0.013517812825739384,  'SBS1415437': 0.01302222441881895,
#     'Pox186': 0.07350756227970123,    'UGC4483': 0.038351867347955704,
#     'VIIZw403': 0.045312829315662384, 'NGC2366': 0.041468776762485504,
#     'MRK1450': 0.015245169401168823,  'UM461': 0.02103571593761444,
#     'NGC4861': 0.015333977527916431,  'UGCA281': 0.020988432690501213,
#     'He2-10': 0.12479636818170547,
# }

# ------------------------------------------------------------------ NED data
summary = pd.read_csv(UP + "ned_output.csv").set_index("input_name")
ext = pd.read_csv(UP + "ned_output_extinction.csv")

# SF11: full-precision A_V straight from OverviewOfObject
ned_sf11 = {k: v / RV_SF11 for k, v in summary["A_V"].items() if np.isfinite(v)}

# SFD98: Landolt V row of the legacy reference in the per-band table
sfd_rows = ext[(ext["Bandpass"] == "Landolt V")
               & (ext["Reference"] == "1998ApJ...500..525S")]
ned_sfd98 = {r["input_name"]: r["Galactic Extinction"] / RV_SFD98
             for _, r in sfd_rows.iterrows()}

# -------------------------------------------------------------------- layout




names = sorted(planck_ebv, key=lambda k: -planck_ebv[k])
xi = {n: i for i, n in enumerate(names)}
x = np.arange(len(names))

C_PLANCK = "#2b6cb0"
C_SF11 = "#c05621"
C_SFD = "#718096"

rx = [xi[n] for n in names if n in ned_sf11]
ry = [planck_ebv[n] / ned_sf11[n] for n in names if n in ned_sf11]
med = np.median(ry)

print(f"\n{'galaxy':<12} {'E(B-V) Planck':>14} {'SF11':>9} {'SFD98':>9} "
      f"{'P/SF11':>8} {'A_V Planck':>11} {'A_V SF11':>9}")
for n in names:
    s = ned_sf11.get(n, np.nan)
    d = ned_sfd98.get(n, np.nan)
    print(f"{n:<12} {planck_ebv[n]:>14.5f} {s:>9.5f} {d:>9.5f} "
          f"{planck_ebv[n]/s if np.isfinite(s) else np.nan:>8.2f} "
          f"{R_V*planck_ebv[n]:>11.4f} "
          f"{R_V*s if np.isfinite(s) else np.nan:>9.4f}")
print(f"\nmedian Planck/SF11 = {med:.3f}   "
      f"mean = {np.mean(ry):.3f}   scatter = {np.std(ry):.3f}")


fig, (ax, axr) = plt.subplots(
    2, 1, figsize=(13.5, 8.0), sharex=True,
    gridspec_kw=dict(height_ratios=[3.0, 1.0], hspace=0.06))

# ---- main panel: A_V ----
ax.plot(x - 0.20, [R_V * planck_ebv[n] for n in names], "o", ms=8,
        color=C_PLANCK, mec="white", mew=1.1, zorder=5)

sf_x = [xi[n] for n in names if n in ned_sf11]
sf_y = [R_V * ned_sf11[n] for n in names if n in ned_sf11]
ax.plot(sf_x, sf_y, "o", ms=8, color=C_SF11, mec="white", mew=1.1, zorder=5)

if SHOW_SFD98:
    sd_x = [xi[n] + 0.20 for n in names if n in ned_sfd98]
    sd_y = [R_V * ned_sfd98[n] for n in names if n in ned_sfd98]
    ax.plot(sd_x, sd_y, "s", ms=7, mfc="none", mec=C_SFD, mew=1.6, zorder=4)

# # connect Planck <-> SF11 to make the offset legible
# for n in names:
#     if n in ned_sf11:
#         ax.plot([xi[n] - 0.20, xi[n]],
#                 [R_V * planck_ebv[n], R_V * ned_sf11[n]],
#                 "-", color="0.6", lw=0.9, zorder=2)

# for n in names:
#     if n not in ned_sf11:
#         ax.annotate("no NED\nmatch", xy=(xi[n], R_V * planck_ebv[n] * 1.30),
#                     ha="center", va="bottom", fontsize=7.5, color="0.45",
#                     bbox=dict(boxstyle="round,pad=0.25", fc="white",
#                               ec="0.7", lw=0.8, alpha=0.9))

ax.set_yscale("log")
ax.set_ylim(0.016, 0.75)
ax.set_xlim(-0.8, len(names) - 0.2)
ax.set_ylabel(r"$A_V = R_V\,E(B-V)$   [mag]", fontsize=11)
ax.grid(axis="y", which="major", ls=":", lw=0.7, color="0.75", zorder=0)
ax.grid(axis="y", which="minor", ls=":", lw=0.4, color="0.9", zorder=0)
for i in x[1::2]:
    ax.axvspan(i - 0.5, i + 0.5, color="0.5", alpha=0.055, zorder=0, lw=0)

secax = ax.secondary_yaxis("right", functions=(lambda a: a / R_V,
                                               lambda e: e * R_V))
secax.set_ylabel(r"$E(B-V)$   [mag]", fontsize=11)

handles = [
    Line2D([], [], marker="o", ls="none", ms=8, color=C_PLANCK, mec="white",
           label="Planck  $E(B-V)$"),
    Line2D([], [], marker="o", ls="none", ms=8, color=C_SF11, mec="white",
           label="NED  Schlafly & Finkbeiner (2011)"),
]
if SHOW_SFD98:
    handles.append(Line2D([], [], marker="s", ls="none", ms=7, mfc="none",
                          mec=C_SFD, mew=1.6, label="NED  SFD (1998), legacy"))
ax.legend(handles=handles, loc="upper right", fontsize=9.5, framealpha=0.95,
          borderpad=0.7, labelspacing=0.6)
ax.set_title(rf"Galactic foreground extinction:  Planck vs NED     "
             rf"($R_V$ = {R_V:g})", fontsize=13, pad=12, loc="left")

# ---- ratio panel ----

axr.plot(rx, ry, "o", ms=7, color=C_PLANCK, mec="white", mew=1.0, zorder=5)
axr.axhline(1.0, color="0.3", lw=1.1, zorder=3)
axr.axhline(med, color=C_SF11, lw=1.3, ls="--", zorder=3)
axr.annotate(f"median = {med:.2f}", xy=(len(names) - 0.5, med), xytext=(-6, 5),
             textcoords="offset points", ha="right", va="bottom",
             fontsize=9, color=C_SF11)

axr.set_ylabel(r"$\dfrac{E(B-V)_{\rm Planck}}{E(B-V)_{\rm NED,SF11}}$",
               fontsize=11)
axr.set_ylim(0.8, 2.1)
axr.grid(axis="y", ls=":", lw=0.7, color="0.8", zorder=0)
for i in x[1::2]:
    axr.axvspan(i - 0.5, i + 0.5, color="0.5", alpha=0.055, zorder=0, lw=0)
axr.set_xticks(x)
axr.set_xticklabels(names, rotation=55, ha="right", fontsize=9.5)

fig.subplots_adjust(left=0.085, right=0.925, top=0.93, bottom=0.20)
fig.savefig("./extinction_comparison.png", dpi=200)
# fig.savefig(OUT + "extinction_comparison.pdf")
plt.show()
print("saved")

