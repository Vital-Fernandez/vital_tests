import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm, norm
import pyneb as pn


# Scipy formula for truncation coefficient
def truncation_limits(mu, sigma, lower_limit, upper_limit):
    return (lower_limit - mu) / sigma, (upper_limit - mu) / sigma


# Function to generate a truncated normal function
def truncated_gaussian(diag_int, diag_err, n_steps, low_limit=-np.infty, up_limit=np.infty):
    a, b = truncation_limits(diag_int, diag_err, low_limit, up_limit)
    output_dist = truncnorm.rvs(a, b, loc=diag_int, scale=diag_err, size=n_steps)
    return output_dist


# Division with error propagation
def div_with_err(a, a_err, b, b_err, sig=3):
    c = a/b
    c_err = c * np.sqrt(np.power(a_err/a, 2) + np.power(b_err/b, 2))
    return np.round(c, sig), np.round(c_err, sig)


# Declare observational data
objList = ['J1152', 'J0925'] # List of objects

compLabels = {'J1152': ['NARROW'], # Number of components per line
              'J0925': ['NARROW']} # Number of components per line

# Fluxes of the emission lines normalized by Hbeta
data_dict = {'J1152': {'O2_3726_flux': np.array([0.000]), # One value per component
                       'O2_3726_errF': np.array([0.000]), # One error per component
                       'O2_3729_flux': np.array([0.000]), # One value per component
                       'O2_3729_errF': np.array([0.000])}, # One error per component

             'J0925': {'O2_3726_flux': np.array([0.000]), # One value per component
                       'O2_3726_errF': np.array([0.000]), # One error per component
                       'O2_3729_flux': np.array([0.000]), # One value per component
                       'O2_3729_errF': np.array([0.000])}} # One error per component

O3_temp = {'J1152': np.array([10000.0, 900]),  # Magnitude and uncertainty
           'J0925': np.array([10000.0, 410])}  # Magnitude and uncertainty

O2_diagnostic_label = 'L(3729)/L(3726)'

# Plots configuration
colors = {'J1152': ['tab:blue', 'tab:red', 'tab:purple'],
          'J0925': ['tab:blue', 'tab:green', 'tab:red', 'tab:purple']}

compLabels = {'J1152': ['NARROW'],
              'J0925': ['NARROW']}

# Number of chains for Monte Carlo error propagation
n_steps = 5000
percentiles_array = np.array([0.0013, 0.0227, 0.1587, 0.5, 0.8413, 0.9772, 0.9987])

# Declaring pyneb objects for ion emissivity computations (using default atomic data)
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)

# Limits for the truncated ROII ratio for densities below and above ~75 and ~10000 cm^-3 respectively
ROII_limits = np.array([1.48, 0.35])

# Loop through the galaxies
neOII_dict, neOIIdist_dict = {}, {}
for i, galaxy in enumerate(objList):

    # Normal distribution for the temperature
    temp, temp_err = O3_temp[galaxy]
    temp_dist = norm.rvs(loc=temp, scale=temp_err, size=n_steps)

    # Dictionary to store the results
    neOII_dict[galaxy], neOIIdist_dict[galaxy] = {}, {}

    # Loop through the emission components
    for j, k_comp in enumerate(compLabels[galaxy]):

        # Compute intensity ratio
        linesDict = data_dict[galaxy]
        ROII, ROII_err = div_with_err(linesDict['O2_3729_flux'][j], linesDict['O2_3729_errF'][j],
                                      linesDict['O2_3726_flux'][j], linesDict['O2_3726_errF'][j])

        # Truncated normal distribution for the emission fluxes
        RSII_dist = truncated_gaussian(ROII, ROII_err, n_steps, low_limit=ROII_limits[1], up_limit=ROII_limits[0])

        # Compute density mean value
        neOII = O2.getTemDen(ROII, tem=temp, to_eval=O2_diagnostic_label)
        neOII_dict[galaxy][k_comp] = neOII

        # Compute density distribution
        neOII_dist = O2.getTemDen(RSII_dist, tem=temp_dist, to_eval=O2_diagnostic_label)
        neOIIdist_dict[galaxy][k_comp] = neOII_dist

        # Compute the percentiles
        comp_den_percentiles = np.percentile(neOII_dist, percentiles_array * 100)
        # comp_den_percentiles = np.empty(len(percentiles_array))
        # for i_per, percen in enumerate(percentiles_array):
        #     comp_den_percentiles[i_per] = np.percentile(neOII_dist, percentiles_array)

        # Print results
        results_line = f'{galaxy}-{k_comp}: mean density {neOII:.2f}, ' \
                       f'distribution {np.nanmean(neOII_dist):.2f} +/- {np.nanstd(neOII_dist):.2f}, ' \
                       f'emissivity errors {np.isnan(neOII_dist).sum()}'
        percentiles_line = f'percentiles: {[np.round(per, 2) for per in comp_den_percentiles]}'
        print(results_line)
        print(percentiles_line)
        print()



for i, galaxy in enumerate(objList):

    line_components = compLabels[galaxy]
    fig, axs = plt.subplots(1, len(line_components), figsize=(len(line_components)*3, 3), sharey=True)

    for j, comp in enumerate(line_components):
        axs[j].hist(neOIIdist_dict[galaxy][comp], histtype='stepfilled', bins=30, alpha=.7, density=False, color=colors[galaxy][j])
        axs[j].axvline(neOII_dict[galaxy][comp], color=colors[galaxy][j], linestyle='--', label=f'ne[SII] mean: {neOII_dict[galaxy][comp]:.0f} ' + r'$cm^{-3}$')

        # Plot distribution uncertainty region
        dist_mean, dist_std = neOIIdist_dict[galaxy][comp].mean(), neOIIdist_dict[galaxy][comp].std()
        data_label = f'ne[SII] distribution: ' + r'${:.0f}\pm{:.0f}\,cm^{{-3}}$'.format(dist_mean, dist_std)
        axs[j].axvspan(dist_mean - dist_std, dist_mean + dist_std, alpha=0.2, color='grey', label=data_label)

        # Plot labels
        axs[j].set_title(comp)
        axs[j].legend()

    fig.suptitle(f'Galaxy {galaxy}')
    plt.tight_layout()
    plt.show()
