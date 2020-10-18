import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import fitelp.constants as constants
from fitelp.make_latex_tables import average_velocities_table_to_latex, halpha_regions_table_to_latex
from fitelp.bpt_plotting import bpt_plot
from fitelp.kinematics_calculations import RegionCalculations
from fitelp.fit_line_profiles import plot_profiles
from fitelp.line_profile_info import RegionParameters
import src.specsiser as sr


# Path to the directory you wish to save the ouput plots, tables and results.
constants.OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Output_Files')

# Path to the directory containing your input data files (Optional).
constants.DATA_FILES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Input_Data_Files')

# Set up example region to simultaneously fit multiple emission lines
example_object = RegionParameters(region_name='example-object',
                                blue_spec_file='exampleB.fc.fits', #fits file path of the blue spectrum
                                  red_spec_file='exampleR.fc.fits', #fits file path of the red spectrum
                                  blue_spec_error_file='exampleB_ErrorFlux.fc.fits', #fits file path of the blue spectrum error
                                  red_spec_error_file='exampleR_ErrorFlux.fc.fits', #fits file path of the red spectrum error
                                  scale_flux=1e15, # Scales the fluxes
                                  center_list={'low': [3649.11211, 3661.84195, 3648.06497],
                                               'high': [3653.62683, 3664.44579, 3659.42848]}, # Center values of the Gaussian components for each zone
                                  sigma_list={'low': [37.8404349, 17.4726107, 73.8997483],
                                              'high': [30.1209633, 14.8025090, 57.1031261]}, # Sigma values of the Gaussian components for each zone
                                  lin_slope={'low': -5.8183e-07, 'high': -4.5958e-07}, # Linear slope values representing the continuum
                                  lin_int={'low': 0.00890187, 'high': 0.00789864}, # Linear intercept values representing the continuum
                                  num_comps={'low': 3, 'high': 3}, #Number of Gaussian components for each zone
                                  component_labels=['Narrow 1', 'Narrow 2', 'Broad'], #Labels for each of the gaussian components
                                  component_colors=['b', 'r', 'g'], #Colour to plot each of the gaussian components
                                  plotting_x_range=[3400, 4000], # xrange of velocities or delta velocities
                                  sigma_instr_blue=4.8, # Instrumental profile on blue arm
                                  sigma_inst_red=5.9, # Instrumental profile on red arm
                                  distance=52.8,  #Distance to object in Mpc
                                  em_lines_for_avg_vel_calc=['H-Alpha', 'OIII-5007A', 'NII-6584A', 'SII-6717A'], # List of emission-lines to use to calculate the average radial velocity
                                  plot_residuals=True,
                                  show_systemic_velocity=False, # Assumed False if not defined
                                  systemic_velocity=3650  # Center of most importan emission-line required only if showSystemicVelocity is True
                                  )

# Add emission-lines to fit
example_object.add_em_line(name='H-Alpha', plot_color= 'y', order=  20, filter= 'red', min_idx=  2800, max_idx=3140, rest_wavelength=  6562.82, num_comps=3, amp_list=  [44.999084, 18.236959, 9.312178], zone= 'low', sigma_tsquared=  164.96, comp_limits= {'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from= None)
example_object.add_em_line(name='OIII-5007A', plot_color= 'c', order=  4, filter= 'red', min_idx=  2535, max_idx=2870, rest_wavelength=  5006.84, num_comps=3, amp_list=  [15.056061, 4.566674, 5.261243], zone= 'high', sigma_tsquared=  10.39,  comp_limits= {'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from= None)
example_object.add_em_line(name='SII-6717A', plot_color= 'r', order=  22, filter= 'red', min_idx=  530, max_idx=850, rest_wavelength=  6716.47, num_comps=3, amp_list=  [4.68714, 2.259787, 1.839585], zone= 'low', sigma_tsquared=  5.19, comp_limits= {'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from= None)
example_object.add_em_line(name='SII-6731A', plot_color= '#58D68D', order=  22, filter= 'red', min_idx=  836, max_idx=1150, rest_wavelength=  6730.85, num_comps=3, amp_list=  [3.34683, 1.878706, 1.238023], zone= 'low', sigma_tsquared=  5.19,comp_limits= {'a': np.inf, 'c': False, 's': False}, copy_from= 'SII-6717A')

# You may fit multiple objects by adding extra objects to this list
regions_parameters = [example_object,]

region_array = []
for rp in regions_parameters:
    region = RegionCalculations(rp, xAxis='vel', initVals='vel')
    region_array.append(region.lineInArray)
    plot_profiles(['H-Alpha', 'OIII-5007A', 'NII-6584A', 'SII-6717A'], rp, nameForComps='SII-6717A',
                  title=rp.regionName + ' Strongest Emission Lines', sortedIndex=[0, 1, 2, 3], logscale=True,
                  ymin=None)



