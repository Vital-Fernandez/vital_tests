from pathlib import Path
from matplotlib import pyplot as plt
import lime
from astropy.io import fits
import numpy as np

# Data location
project_folder = Path('/')
obs_folder = Path('/home/vital/Astrodata/STScI')

# Configuration
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

idcs_manual = sample_df.index.get_level_values('state') == 'hasp_aspec'
obj_list = sample_df.loc[idcs_manual, 'object'].unique()

for i, obj_name in enumerate(obj_list):
    n_objects = sample_df.loc[(sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'x1d')].index.size
    idcs = (sample_df.object == obj_name) & idcs_manual
    sub_sample = sample_df.loc[idcs]
    n_cspec_manual = sub_sample.index.size

    # Manual cspec
    if sub_sample.index.size == 0:
        print(f'-{obj_name} NO SPECTRA')

    else:
        redshift = cfg_sample['Galaxy_redshifts'][obj_name]
        print(f'- {obj_name} at z = {redshift}')

        for j, idx in enumerate(sub_sample.index):
            file_path = obs_folder/sub_sample.loc[idx, 'filepath']

            if j == 0:
                print(fits.info(file_path))
                spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
                label = f'{obj_name} - {idx[2]} - {idx[3]}'
                spec.plot.spectrum(label=label, maximize=True, in_fig=None)
            else:
                spec_j = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
                label = f'{obj_name} - {idx[2]} - {idx[3]}'
                spec.plot.ax.step(spec_j.wave, spec_j.flux, where='mid', linestyle=':', color='black', label=label, zorder=5,)
            print('--', file_path, idx)

        # MAST cspec
        idcs = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'hasp_mult')
        sub_sample = sample_df.loc[idcs]
        n_cspec_mast = sub_sample.index.size
        if sub_sample.index.size == 0:
            print(f'-{obj_name} NO SPECTRA')
        else:
            for j, idx in enumerate(sub_sample.index):
                file_path = obs_folder / sub_sample.loc[idx, 'filepath']
                spec_j = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
                label = f'{obj_name} - {idx[2]} - {idx[3]}'
                spec.plot.ax.step(spec_j.wave, spec_j.flux, where='mid', linestyle='--', label=label, zorder=1)

        spec.plot.ax.legend()
        spec.plot.ax.set_title(f'{obj_name}, n_x1d = {n_objects}, n_cspec_manual={n_cspec_manual}')#, n_cspec_mast={n_cspec_mast}')
        spec.plot.ax.set_yscale('symlog', linthresh=10)
        plt.tight_layout()
        plt.show()

# # Check all ids are in samples
# full_ids = np.sort(sample_df.index.get_level_values('id').unique())
# hasp_ids = np.sort(hasp_sample.index.get_level_values('id').unique())
# print("✅ All IDs are in HASP sample" if set(hasp_ids).issubset(full_ids) else "❌ IDs from sample missing in HASP")
#
# # A values mapping to multiple Bs
# group_PID = sample_df.groupby("object")["PID"].nunique()
# obj_list_singleID = group_PID[~(group_PID > 1)].index.tolist()
# obj_list_multID = group_PID[group_PID > 1].index.tolist()
# output_folder_single = obs_folder/'LyC_leakers_COS/hasp_single'
# output_folder_mult = obs_folder/'LyC_leakers_COS/hasp_multi'





# for i, obj_name in enumerate(sample_df.object.unique()):
#     idcs = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'hasp_manual')
#     sub_sample = sample_df.loc[idcs]
#
#     if sub_sample.index.size == 0:
#         print(f'-{obj_name} NO SPECTRA')
#
#     else:
#         redshift = cfg_sample['Galaxy_redshifts'][obj_name]
#         print(f'- {obj_name} at z = {redshift}')
#
#         for j, idx in enumerate(sub_sample.index):
#             file_path = obs_folder/sub_sample.loc[idx, 'filepath']
#
#             if j == 0:
#                 print(fits.info(file_path))
#                 spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
#                 label = f'{obj_name} - {idx[2]} - {idx[3]}'
#                 spec.plot.spectrum(label=label, maximize=True, in_fig=None)
#             else:
#                 spec_j = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
#                 label = f'{obj_name} - {idx[2]} - {idx[3]}'
#                 spec.plot.ax.step(spec_j.wave, spec_j.flux, where='mid', linestyle=':', label=label)
#             print('--', file_path, idx)
#
#         spec.plot.ax.legend()
#         spec.plot.ax.set_yscale('symlog', linthresh=10)
#         plt.tight_layout()
#         plt.show()

        # idcs_mast = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'hasp_mast')
        # sub_sample_mast = sample_df.loc[idcs_mast]
        # spec_mast = lime.Spectrum.from_file(obs_folder/sub_sample_mast['filepath'].values[0], instrument='cos', redshift=redshift)
        # ref_mast = f'{obj_name}_{sub_sample_mast.index.values[0][3]}{sub_sample_mast.index.values[0][2]}'
        #
        # self.ax.set_yscale('symlog', linthresh=10)
        # spec.plot.spectrum(label=ref, in_fig=None)
        # spec.plot.ax.step(spec_mast.wave, spec_mast.flux, where='mid', label=ref_mast, linestyle=':')
        # spec.plot.ax.legend()
        # plt.tight_layout()
        # plt.show()

        # print(f'{obj_name} ({len(sub_sample.index)})')


# input_folder = obs_folder/'LyC_leakers_COS'/'Calibrations'
# output_folder = input_folder/'hasp_products'
# wrapper.main(input_folder, outdir=output_folder)

# list_spec = ['hst_17515_cos_mrk-209_g130m_lf9g01_cspec.fits', 'hst_17515_cos_mrk-209_g130m_lf9g_cspec.fits']
# list_spec2 = ['/home/vital/Astro-data/STScI/LyC_leakers_COS/Proposal_galaxies/lf9g01wpq_x1d.fits',
#               '/home/vital/Astro-data/STScI/LyC_leakers_COS/Proposal_galaxies/lf9g01wrq_x1d.fits']
#
# obj_name = 'UGC4483'
# i, redshift = 0, cfg_sample['Galaxy_redshifts']['UGC4483']
# spec = lime.Spectrum.from_file(output_folder/list_spec[i], instrument='cos', redshift=redshift)
# spec.plot.spectrum(label='lf9g01_cspec', in_fig=None)
#
# for spec_extra in list_spec[1:]:
#     ref = "_".join(spec_extra.rsplit(".", 1)[0].rsplit("_", 2)[1:])
#     spec2 = lime.Spectrum.from_file(output_folder/spec_extra, instrument='cos', redshift=redshift)
#     spec.plot.ax.step(spec2.wave, spec2.flux, label=ref, linestyle='--', where='mid')
#
# for spec_extra in list_spec2:
#     ref = Path(spec_extra).stem
#     spec2 = lime.Spectrum.from_file(spec_extra, instrument='cos', redshift=redshift)
#     spec.plot.ax.step(spec2.wave, spec2.flux, label=ref, linestyle=':', where='mid')
#
# spec.plot.ax.legend()
# plt.tight_layout()
# plt.show()


# Program_10033 = {'name': '10033',
#                  'url': 'https://stsci.box.com/shared/static/wnmma331eixblioo5jfnmjnzz2fjt01r.gz'}
# Program_11839 = {'name': '11839',
#                  'url': 'https://stsci.box.com/shared/static/a05mixac3hreg6g07f1kkogcpphdo26a.gz'}
# Program_13471 = {'name': '13471',
#                  'url': 'https://stsci.box.com/shared/static/jac75olate8hjalvc3tgor1fuvxupiqm.gz'}
#
# HD104237E = {'name': 'hd-104237e',
#              'url': 'https://stsci.box.com/shared/static/irelgbjs7zhdiljksao4su2c6tfhyntx.gz'}
# V_HK_ORI = {'name': 'v-hk-ori',
#             'url': 'https://stsci.box.com/shared/static/az3ytnwohj0t4wnc4m8oqbw9ziozwqx1.gz'}
#
# program = Program_10033
#
# # Set the tree
# program_name = program['name']
# url = program['url']
# if os.path.isdir(program_name):
#     shutil.rmtree(program_name)
# filename = program_name + '.tar.gz'
# _ = urllib.request.urlretrieve(url, filename)
# data_tarfile = tarfile.open(filename, mode='r|gz')
# data_tarfile.extractall(filter='data')
#
# # Run wrapper
# indir = program_name + '/input/'
# wrapper.main(indir, outdir=indir)

# class TestWrapper():
#
#     def test_10033(self):
#         program = Program_10033
#         self.setup_tree(program)
#         self.run_wrapper(program['name'])
#         report = self.compare_outputs(program['name'])
#         self.cleanup(program['name'])
#         if report is not None:
#             raise AssertionError(report)
#         return
#
#     def test_11839(self):
#         program = Program_11839
#         self.setup_tree(program)
#         self.run_wrapper(program['name'])
#         report = self.compare_outputs(program['name'])
#         self.cleanup(program['name'])
#         if report is not None:
#             raise AssertionError(report)
#         return
#
#     def test_13471(self):
#         program = Program_13471
#         self.setup_tree(program)
#         self.run_wrapper(program['name'])
#         report = self.compare_outputs(program['name'])
#         self.cleanup(program['name'])
#         if report is not None:
#             raise AssertionError(report)
#         return
#
#     def setup_tree(self, program):
#         program_name = program['name']
#         url = program['url']
#         if os.path.isdir(program_name):
#             shutil.rmtree(program_name)
#         filename = program_name + '.tar.gz'
#         _ = urllib.request.urlretrieve(url, filename)
#         data_tarfile = tarfile.open(filename, mode='r|gz')
#         data_tarfile.extractall(filter='data')
#         return
#
#     def run_wrapper(self, program):
#         indir = program + '/input/'
#         wrapper.main(indir, outdir=indir)
#         return
#
#     def compare_outputs(self, program):
#         report = None
#         # Outputs from current run are in ./, truth files to compare
#         # with are in ./truth
#         all_ok = True
#         fitsdiff_report = ''
#         keywords_to_ignore = ['DATE', 'FITS_SW', 'FILENAME',
#                               'HLSP_VER', 'S_REGION', 'CAL_VER']
#         new_hlsps = glob.glob(program + '/input/hst_*')
#         for new_product in new_hlsps:
#             truth_filename = self.get_truth_filename(program, new_product)
#             fdiff = FITSDiff(new_product, truth_filename,
#                              ignore_keywords=keywords_to_ignore,
#                              rtol=3.0e-7)
#             fitsdiff_report += fdiff.report()
#             if not fdiff.identical and all_ok:
#                 all_ok = False
#         if not all_ok:
#             report = os.linesep + fitsdiff_report
#             return report
#         print(fitsdiff_report)
#         return None
#
#     def get_truth_filename(self, program, product):
#         # Get the truth filename.  The data release might be different
#         filename = os.path.basename(product)
#         truth_filename = program + '/truth/' + filename
#         return truth_filename
#
#     def cleanup(self, program):
#         shutil.rmtree(program)
#         os.remove(program + '.tar.gz')
#         return


