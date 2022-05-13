import numpy as np
from astropy.wcs import WCS
from matplotlib import pyplot as plt, rcParams, gridspec, cm, colors
# from src.specsiser.tools.line_measure import STANDARD_PLOT, STANDARD_AXES
from astropy.io import fits

STANDARD_PLOT = {'figure.figsize': (20, 14), 'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 12,
                 'xtick.labelsize': 12, 'ytick.labelsize': 12}

background_color = np.array((43, 43, 43))/255.0
foreground_color = np.array((179, 199, 216))/255.0

line_regions = {'H1_6563A': np.array([6850, 6910]),
                'S2_6731A_b': np.array([7027.5, 7057.5]),
                'O3_4363A': np.array([4565.0, 4575.0]),
                'S3_6312A': np.array([6606.5, 6617.0]),
                'O3_5007A': np.array([5232.0, 5260.0]),
                'S3_9069A': np.array([9492.5, 9506.5]),
                'S3_9531A': np.array([9975.5, 9995.0])}

megaradrp_modes = {'bias'       : 'MegaraBiasImage',
                   'arc'        : 'MegaraArcCalibration',
                   'trace_map'  : 'MegaraTraceMap',
                   'slit_flat'  : 'MegaraSlitFlat',
                   'model_map'  : 'MegaraModelMap',
                   'fiber_flat' : 'MegaraFiberFlatImage',
                   'lcb_acq'    : 'MegaraLcbAcquisition',
                   'lcb_std'    : 'MegaraLcbStdStar',
                   'lcb_image'  : 'MegaraLcbImage'}

def open_manga_cubes(file_address):

    with fits.open(file_address) as hdul:
        hdr = hdul['FLUX'].header
        wave = hdul['WAVE'].data
        flux = hdul['FLUX'].data
        err = 1/np.sqrt(hdul['IVAR'].data)
        pixMask = hdul['MASK'].data

        return wave, flux, err, hdr