import os
import numpy as np

from astropy.io import fits
from astropy.table import Table, vstack, join
from astropy.convolution import convolve, Gaussian1DKernel
from desi_functions import desi_mask, bgs_mask, scnd_mask
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
# specprod = "fuji"    # Internal name for the EDR
# specprod_dir = "/global/cfs/cdirs/desi/public/edr/spectro/redux/fuji/"
# print(specprod_dir)


def open_desi_spectra(file_path, obj_idtarget=None, obj_idrows=None):

    # Confirm the user provides and identity
    if (obj_idtarget is None) and (obj_idrows is None):
        assert (obj_idtarget is None) and (obj_idrows is None), 'You need to apply introduce and IDTARGET or .fits rows'

    # Confirm file location
    file_path = Path(file_path)
    if not file_path.is_file():
        assert file_path, f'Input file not found at "{file_path}"'

    # Open the file
    with fits.open(fits_path) as hdulist:

        # Get target ID rows if none are provided  by TARGETIDs
        if obj_idrows is None:
            file_idtargets = hdulist["FIBERMAP"].data['TARGETID']
            obj_idrows = np.where(np.isin(file_idtargets, np.atleast_1d(obj_idtarget)))[0]
        else:
            obj_idrows = np.asarray(obj_idrows)

        # Check the target is on file
        assert obj_idrows.size > 0, f'Input TARGETID(s): {obj_idtarget},\nnot found in input file: {file_path}'

        # Load the EXP_FIBERMAP



        # elif
        #     # Indexes of file on spectra
        #     rows = np.where(np.isin(file_targetids, targetids))[0]
        #
        #     if 'coadd' in fits_path:
        #         exp_targetids = hdulist["EXP_FIBERMAP"].data['TARGETID']

    return

specprod_dir = '.'
zpix_file = f'{specprod_dir}/zall-pix-fuji.fits'
fujidata = Table.read(zpix_file, hdu=1)

#-- SV1/2/3
is_sv1 = (fujidata["SURVEY"].astype(str).data == "sv1")
is_sv2 = (fujidata["SURVEY"].astype(str).data == "sv2")
is_sv3 = (fujidata["SURVEY"].astype(str).data == "sv3")

#-- all SV data
is_sv = (is_sv1 | is_sv2 | is_sv3)

#-- commissioning data
is_cmx = (fujidata["SURVEY"].astype(str).data == "cmx")

#-- special tiles
is_special = (fujidata["SURVEY"].astype(str).data == "special")

#-- all unique targets in SV
targids_sv = np.unique(fujidata["TARGETID"][is_sv])

print("Unique targets in SV:     {0:8}".format( len(targids_sv) ))

print("\nTotal unique targets in combined catalog: {}\n".format( len(np.unique(fujidata["TARGETID"])) ))

print("Unique targets in SV1:    {0:8}".format( len(np.unique(fujidata["TARGETID"][is_sv1])) ))
print("Unique targets in SV2:    {0:8}".format( len(np.unique(fujidata["TARGETID"][is_sv2])) ))
print("Unique targets in SV3:    {0:8}".format( len(np.unique(fujidata["TARGETID"][is_sv3])) ))

print("\nUnique targets in CMX:    {0:8}".format( len(np.unique(fujidata["TARGETID"][is_cmx])) ))

print("\nUnique targets in Special:{0:8}".format( len(np.unique(fujidata["TARGETID"][is_special])) ))

print("Primary spectra in...\n")
print("  SV:              {0:8}".format( np.count_nonzero(is_sv & fujidata["SV_PRIMARY"])))
print("  combined catalog:{0:8}\n".format( np.count_nonzero(fujidata["ZCAT_PRIMARY"])))

print("  SV1:    {0:8}".format( np.count_nonzero(is_sv1 & fujidata["SV_PRIMARY"]) ))
print("  SV2:    {0:8}".format( np.count_nonzero(is_sv2 & fujidata["SV_PRIMARY"]) ))
print("  SV3:    {0:8}\n".format( np.count_nonzero(is_sv3 & fujidata["SV_PRIMARY"]) ))

print("  CMX:    {0:8}\n".format( np.count_nonzero(is_cmx & fujidata["ZCAT_PRIMARY"]) ))

print("  Special:{0:8}".format( np.count_nonzero(is_special & fujidata["ZCAT_PRIMARY"]) ))

#-- select the target mask for each type of target (BGS, LRG, ELG, QSO, etc.)

bgs_tgtmask  = desi_mask["BGS_ANY"]
lrg_tgtmask  = desi_mask["LRG"]
elg_tgtmask  = desi_mask["ELG"]
qso_tgtmask  = desi_mask["QSO"]
mws_tgtmask  = desi_mask["MWS_ANY"]
scnd_tgtmask = desi_mask["SCND_ANY"]

desi_target = fujidata["DESI_TARGET"]

#-- all BGS targets
is_bgs = (desi_target & bgs_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & bgs_tgtmask != 0) | \
         (fujidata["SV2_DESI_TARGET"] & bgs_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & bgs_tgtmask != 0)

#-- all LRG targets
is_lrg = (desi_target & lrg_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & lrg_tgtmask != 0) | \
         (fujidata["SV2_DESI_TARGET"] & lrg_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & lrg_tgtmask != 0)

#-- all ELG targets
is_elg = (desi_target & elg_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & elg_tgtmask != 0) | \
         (fujidata["SV2_DESI_TARGET"] & elg_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & elg_tgtmask != 0)

#-- all QSO targets
is_qso = (desi_target & qso_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & qso_tgtmask != 0) | \
         (fujidata["SV2_DESI_TARGET"] & qso_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & qso_tgtmask != 0)

#-- all Milky Way targets
is_mws = (desi_target & mws_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & mws_tgtmask != 0) | \
         (fujidata["SV2_DESI_TARGET"] & mws_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & mws_tgtmask != 0)

#-- all secondary targets
is_scnd = (desi_target & scnd_tgtmask != 0) | (fujidata["SV1_DESI_TARGET"] & scnd_tgtmask != 0) | \
          (fujidata["SV2_DESI_TARGET"] & scnd_tgtmask != 0) | (fujidata["SV3_DESI_TARGET"] & scnd_tgtmask != 0)

#-- total number of spectra by target type
n_bgs  = np.count_nonzero(is_bgs)
n_lrg  = np.count_nonzero(is_lrg)
n_elg  = np.count_nonzero(is_elg)
n_qso  = np.count_nonzero(is_qso)
n_mws  = np.count_nonzero(is_mws)
n_scnd = np.count_nonzero(is_scnd)

counts = [n_bgs, n_lrg, n_elg, n_qso, n_mws, n_scnd]

#-- number of primary spectra by target type
is_primary = fujidata["ZCAT_PRIMARY"]

n_bgs_prim  = np.count_nonzero(is_bgs & is_primary)
n_lrg_prim  = np.count_nonzero(is_lrg & is_primary)
n_elg_prim  = np.count_nonzero(is_elg & is_primary)
n_qso_prim  = np.count_nonzero(is_qso & is_primary)
n_mws_prim  = np.count_nonzero(is_mws & is_primary)
n_scnd_prim = np.count_nonzero(is_scnd & is_primary)

# counts_prim = [n_bgs_prim, n_lrg_prim, n_elg_prim, n_qso_prim, n_mws_prim, n_scnd_prim]
#
# fig, ax = plt.subplots(1, 1, figsize=(8,6))
#
# targets = ["BGS", "LRG", "ELG", "QSO", "MWS", "SCND"]
#
# ax.bar(targets, counts, color="purple", alpha=0.4, label="All spectra (includes duplicate targets)")
# ax.bar(targets, counts_prim, color="purple", alpha=1, label="Primary spectra (unique targets)")
# ax.set_ylabel("Number of spectra")
# ax.semilogy()
# ax.set_ylim(1e5, 2e6)
#
# for i in range(len(targets)):
#     ax.text(targets[i], counts[i], counts[i], ha="center", va="bottom", fontsize=16)
#     ax.text(targets[i], counts_prim[i], counts_prim[i], ha="center", va="top", fontsize=16, color="white")
#
# plt.legend(fontsize=18, frameon=False, markerfirst=False)
# plt.tight_layout()
# plt.show()
#
# fig, axes = plt.subplots(2, 2, figsize=(16, 10))
#
# bins = np.arange(0, 4, 0.1)
#
# # -------------------------------------------
#
# type_masks = (is_bgs, is_lrg, is_elg, is_qso)
# colors = ("C0", "C1", "C2", "C3")
#
# for ax, type_mask, name, color in zip(np.concatenate(axes), type_masks, targets[:-2], colors):
#     kwargs = dict(color=color, bins=bins)
#
#     ax.hist(fujidata["Z"][type_mask], **kwargs, alpha=0.3, label=f"Total: {len(fujidata[type_mask])}")
#     ax.hist(fujidata["Z"][type_mask], bins=bins, histtype="step", color="black")
#
#     mask = type_mask & is_sv
#     ax.hist(fujidata["Z"][mask], **kwargs, label=f"SV only: {len(fujidata[mask])}")
#     ax.hist(fujidata["Z"][mask], bins=bins, histtype="step", color="black")
#
#     ax.legend(fontsize=20, markerfirst=False, handletextpad=0.25, frameon=False)
#     ax.set_title(f"{name} targets", fontsize=22)
#
# for ax in axes[-1]:
#     ax.set_xlabel("Redshift")
#
# for ax in np.concatenate(axes):
#     ax.set_xlim(-0.1, 4.0)
#     ax.set_xticks(np.arange(0, 4.1, 0.5))
#     ax.set_ylabel("N(z)")
#     ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
#
# plt.tight_layout()
# plt.subplots_adjust(hspace=0.2, wspace=0.15)
# plt.show()

#-- get all targets with NSPEC=5
t_fivespec = fujidata[fujidata["ZCAT_NSPEC"]==5]

#-- unique TARGETID of each object with five spectra
targids = np.unique(t_fivespec["TARGETID"])

#-- get the data for all observations of this TARGETID
special_ID = 39627835576420141
these_spec = t_fivespec[t_fivespec["TARGETID"] == special_ID]

#-- get the SURVEY, PROGRAM, SPECTYPE, and redshift values for each of the five spectra of this object
survey   = these_spec["SURVEY"].data.astype(str)
program  = these_spec["PROGRAM"].data.astype(str)
redshift = np.round(these_spec["Z"].data,5)
spectype = these_spec["SPECTYPE"].data.astype(str)

print("\tSURVEY  PROGRAM  SPECTYPE  REDSHIFT")
for i in range(5):
    print("{0:1}/5\t{1:7} {2:8} {3:8} {4:8}".format(i+1, survey[i], program[i], spectype[i], redshift[i]))

# Reading a spectrum
tid = special_ID
survey = "sv1"
program = "other"

idx = np.where((fujidata["TARGETID"] == tid) & (fujidata["SURVEY"] == survey) & (fujidata["PROGRAM"] == program))[0][0]

# -- healpix values are integers but are converted here into a string for easier access to the file path
hpx = fujidata["HEALPIX"].astype(str)

if "sv" in survey:
    specprod = "fuji"

specprod_dir = f"/global/cfs/cdirs/desi/spectro/redux/{specprod}"
target_dir = f"{specprod_dir}/healpix/{survey}/{program}/{hpx[idx][:-2]}/{hpx[idx]}"
coadd_fname = f"coadd-{survey}-{program}-{hpx[idx]}.fits"
fits_path = f"{target_dir}/{coadd_fname}"
fits_path = f"./{coadd_fname}"

targetids = special_ID    # list of target ids to read from file
rows = None               # list of rows to read from file (index of target ID on file)

open_desi_spectra(fits_path, targetids, rows)




