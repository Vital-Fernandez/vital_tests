import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import plot_parameter_image, compute_parameter_distributions

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Extensions for fitting files
regions_list = ['MASK_0', 'MASK_1', 'MASK_2']
param_list = np.array(['logOH', 'logNO', 'logU'])

# Measurement files and reference files
ref_simulations = ['HIICHImistry-PopStar']
tech_label = 'HII-CHI-mistry'
params_conv_dict = {'logOH': 'OH', 'logNO': 'NO'}

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = objFolder/'chemistry'
    db_address = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Header with astronomical data
    hdr = fits.getheader(maskFits_address, extname='MASK_0')

    # Loop through the measurement fits
    for j, ref_fit in enumerate(ref_simulations):

        # Plot the parameter value against the image background
        plot_parameter_image(db_address, param_list, chemFolder, ref_fit, tech_label)

        # Generate the parameter maps and store the distributions
        compute_parameter_distributions(param_list, chemFolder, ref_fit, maskFits_address, regions_list, tech_label,
                                        params_conv_dict)