import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import quad
from scipy.interpolate import interp1d
from pathlib import Path
import lime

C = 299792.458
H0, OM, W0 = 72.7, 0.266, -0.96


def d_lum(z, h=H0 / 100, om=OM, w0=W0):
    dh = C / (100.0 * h)
    dc = dh * quad(lambda x: 1.0 / np.sqrt(om * (1 + x) ** 3
                                           + (1 - om) * (1 + x) ** (3 * (1 + w0))),
                   0.0, z, epsabs=1e-12, epsrel=1e-12)[0]
    return (1 + z) * dc


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')
adopted = cfg_sample['Galaxy_redshifts']

dl = {k: d_lum(z) for k, z in adopted.items()}
# asymmetric cosmological-parameter errors, recomputed as before
err_p, err_m = {}, {}
for k, z in adopted.items():
    d0 = dl[k]
    dp, dm = [], []
    for idx, lo, hi in [(0, 0.687, 0.767), (1, 0.187, 0.386), (2, -1.17, -0.43)]:
        for val in (lo, hi):
            p = [H0 / 100, OM, W0]
            p[idx] = val
            d = d_lum(z, *p)
            (dp if d > d0 else dm).append(abs(d - d0))
    err_p[k] = np.sqrt(sum(x ** 2 for x in dp))
    err_m[k] = np.sqrt(sum(x ** 2 for x in dm))

# ---- NED tables ----
up = "./"
summary = pd.read_csv(up + "ned_output.csv")
dists = pd.read_csv(up + "ned_output_distances.csv")

summary = summary.set_index("input_name")
ned_z = summary["z"].to_dict()
ned_dmean = summary["D_mean"].to_dict()
ned_dsem = summary["D_mean_sem"].to_dict()

# ---- ordering: by adopted redshift, descending ----
names = sorted(adopted, key=lambda k: -adopted[k])
x = np.arange(len(names))
xi = {n: i for i, n in enumerate(names)}

# ---- secondary-axis transform: distance <-> redshift under the same cosmology ----
zgrid = np.logspace(np.log10(5e-5), np.log10(0.06), 400)
dgrid = np.array([d_lum(z) for z in zgrid])
_z_of_d = interp1d(dgrid, zgrid, bounds_error=False, fill_value="extrapolate")
_d_of_z = interp1d(zgrid, dgrid, bounds_error=False, fill_value="extrapolate")


def d2z(d):
    d = np.asarray(d, dtype=float)
    return np.where(d > 0, _z_of_d(np.clip(d, dgrid[0], dgrid[-1])), np.nan)


def z2d(z):
    z = np.asarray(z, dtype=float)
    return np.where(z > 0, _d_of_z(np.clip(z, zgrid[0], zgrid[-1])), np.nan)


# ---- final table: the three series plotted above ----
def _d_from_ned_z(name):
    z = ned_z.get(name, np.nan)
    return float(z2d(z)) if np.isfinite(z) and z > 0 else np.nan

table = pd.DataFrame(
    {
        "D_wCDM":     [dl[n] for n in names],                    # blue
        "D_NED_z":    [_d_from_ned_z(n) for n in names],         # green
        "D_NED_mean": [ned_dmean.get(n, np.nan) for n in names], # orange
    },
    index=pd.Index(names, name="galaxy"),
)

print("\nDistances [Mpc]")
print(table.round(2).to_string(na_rep="--"))


C_WCDM = "#2b6cb0"
C_NED = "#c05621"
C_Z = "#2f855a"

fig, ax = plt.subplots(figsize=(13.5, 7.2))

# --- series 1: wCDM luminosity distance from the adopted redshifts ---
d_arr = np.array([dl[n] for n in names])
ax.errorbar(x - 0.28, d_arr,
            yerr=[[err_m[n] for n in names], [err_p[n] for n in names]],
            fmt="o", ms=8, color=C_WCDM, mfc=C_WCDM, mec="white", mew=1.2,
            ecolor=C_WCDM, elinewidth=1.6, capsize=3, zorder=4)

# --- series 2: every NED redshift-independent estimate ---
rng = np.random.default_rng(7)
for name, grp in dists.groupby("input_name"):
    if name not in xi:
        continue
    d = grp["Distance"].astype(float).values
    jitter = np.linspace(-0.13, 0.13, len(d)) if len(d) > 1 else np.array([0.0])
    ax.plot(xi[name] + jitter, d, "o", ms=4.5, color=C_NED, alpha=0.75,
            mec="none", zorder=3)

# mean +/- SEM of the redshift-independent estimates
for name in names:
    dm, ds = ned_dmean.get(name, np.nan), ned_dsem.get(name, np.nan)
    if np.isfinite(dm):
        ax.errorbar(xi[name], dm, yerr=(ds if np.isfinite(ds) else None),
                    fmt="_", ms=22, mew=2.4, color=C_NED,
                    ecolor=C_NED, elinewidth=2.4, capsize=0, zorder=5)

# --- series 3: NED redshift, on the right-hand axis ---
zx, zy = [], []
for name in names:
    z = ned_z.get(name, np.nan)
    if np.isfinite(z) and z > 0:
        zx.append(xi[name] + 0.28)
        zy.append(z2d(z))          # placed via the shared distance<->z mapping
ax.plot(zx, zy, "*", ms=13, color=C_Z, mec="white", mew=0.7, zorder=6)


# VII Zw 403 has a negative NED heliocentric redshift -> cannot be placed
if np.isfinite(ned_z.get("VIIZw403", np.nan)) and ned_z["VIIZw403"] < 0:
    ax.annotate("NED z < 0\n(z = %.6f)" % ned_z["VIIZw403"],
                xy=(xi["VIIZw403"], 0.62), ha="center", va="bottom",
                fontsize=7.5, color=C_Z,
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec=C_Z, lw=0.8, alpha=0.9))

# galaxies NED failed to resolve
for name in names:
    if not np.isfinite(ned_z.get(name, np.nan)) and \
       not np.isfinite(ned_dmean.get(name, np.nan)):
        ax.annotate("no NED\nmatch", xy=(xi[name], 0.62), ha="center", va="bottom",
                    fontsize=7.5, color="0.45",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec="0.7", lw=0.8, alpha=0.9))

ax.plot(zx, zy, "*", ms=13, color=C_Z, mec="white", mew=0.7, zorder=6)

# --- 1000 km/s peculiar-velocity boundary ---
cz = np.array([C * adopted[n] for n in names])
below = cz < 1000.0

if below.any() and not below.all():
    b = np.argmax(below)
    x_boundary = b - 0.5
    ax.axvline(x_boundary, color="0.35", ls="--", lw=1.4, zorder=2)
    ax.annotate(
        r"$cz = 1000\ \mathrm{km\,s^{-1}}$",
        xy=(x_boundary, ax.get_ylim()[1]), xytext=(6, -6),
        textcoords="offset points", ha="left", va="top",
        fontsize=8.5, color="0.3",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.5", lw=0.8, alpha=0.9),
    )

# VII Zw 403 has a negative NED heliocentric redshift -> cannot be placed
...

ax.set_yscale("log")
ax.set_ylim(0.5, 200)
ax.set_xlim(-0.7, len(names) - 0.3)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=55, ha="right", fontsize=9.5)
ax.set_ylabel("Distance  [Mpc]", fontsize=11)
ax.grid(axis="y", which="major", ls=":", lw=0.7, color="0.75", zorder=0)
ax.grid(axis="y", which="minor", ls=":", lw=0.4, color="0.88", zorder=0)
for i in x[1::2]:
    ax.axvspan(i - 0.5, i + 0.5, color="0.5", alpha=0.055, zorder=0, lw=0)

secax = ax.secondary_yaxis("right", functions=(d2z, z2d))
secax.set_ylabel(r"redshift  $z$   (flat $w$CDM mapping)", fontsize=11)

handles = [
    Line2D([], [], marker="o", ls="none", ms=8, color=C_WCDM, mec="white",
           label=r"$D_L$ from adopted $z$  (flat $w$CDM: $h$=0.727, "
                 r"$\Omega_m$=0.266, $w_0$=$-$0.96)"),
    Line2D([], [], marker="o", ls="none", ms=5, color=C_NED, alpha=0.8,
           label="NED redshift-independent distances (individual estimates)"),
    Line2D([], [], marker="_", ls="none", ms=16, mew=2.4, color=C_NED,
           label=r"NED mean redshift-independent distance $\pm$ SEM"),
    Line2D([], [], marker="*", ls="none", ms=13, color=C_Z, mec="white",
           label=r"NED heliocentric $z$  (read on right axis)"),
]
ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.95,
          borderpad=0.7, labelspacing=0.7)

ax.set_title("Redshift-derived vs. redshift-independent distances",
             fontsize=13, pad=12, loc="left")
fig.tight_layout()
plt.show()
fig.savefig("./distance_comparison.png", dpi=200)
# fig.savefig("/mnt/user-data/outputs/distance_comparison.pdf")
print("saved")

