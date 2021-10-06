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
objList = ['J1152', 'J0925']

compLabels = {'J1152': ['NARROW', 'BROAD', 'GLOBAL'],
              'J0925': ['NARROW1', 'NARROW2', 'BROAD', 'GLOBAL']}

data_dict = {'J1152': {'S2_6717_flux': np.array([15.948, 11.279, 27.227]),
                       'S2_6717_errF': np.array([1.634, 2.502, 2.988]),
                       'S2_6731_flux': np.array([13.706, 13.69, 27.396]),
                       'S2_6731_errF': np.array([1.67, 2.549, 3.047])},

             'J0925': {'S2_6717_flux': np.array([16.978, 16.176, 2.876, 36.030]),
                       'S2_6717_errF': np.array([0.826, 0.633, 1.609, 1.916]),
                       'S2_6731_flux': np.array([12.477, 14.289, 4.29, 31.056]),
                       'S2_6731_errF': np.array([0.852, 0.723, 1.743, 2.07])}}

O3_temp = {'J1152': np.array([13430.0, 900]),
           'J0925': np.array([15010.0, 410])}

Hbeta_flux = {'J1152': np.array([127.579, 130.119, 257.698]),
              'J0925': np.array([73.504, 96.340, 85.013, 254.856])}

# Number of chains for Monte Carlo error propagation
n_steps = 5000
percentiles_array = np.array([0.0013, 0.0227, 0.1587, 0.5, 0.8413, 0.9772, 0.9987])

# Declaring pyneb objects for ion emissivity computations (using default atomic data)
S2 = pn.Atom('S', 2)
O3 = pn.Atom('O', 3)
H1 = pn.RecAtom('H', 1)

# Limits for the truncated RSII ratio for densities below and above ~75 and ~10000 cm^-3 respectively
RSII_limits = np.array([1.35, 0.50])

# Loop through the galaxies
neSII_dict, neSIIdist_dict = {}, {}
for i, galaxy in enumerate(objList):

    # Normal distribution for the temperature
    temp, temp_err = O3_temp[galaxy]
    temp_dist = norm.rvs(loc=temp, scale=temp_err, size=n_steps)

    # Dictionary to store the results
    neSII_dict[galaxy], neSIIdist_dict[galaxy] = {}, {}

    # Loop through the emission components
    for j, k_comp in enumerate(compLabels[galaxy]):

        # Compute intensity ratio
        linesDict = data_dict[galaxy]
        RSII, RSII_err = div_with_err(linesDict['S2_6717_flux'][j], linesDict['S2_6717_errF'][j],
                                      linesDict['S2_6731_flux'][j], linesDict['S2_6731_errF'][j])

        # Truncated normal distribution for the emission fluxes
        RSII_dist = truncated_gaussian(RSII, RSII_err, n_steps, low_limit=RSII_limits[1], up_limit=RSII_limits[0])

        # Compute density mean value
        neSII = S2.getTemDen(RSII, tem=temp, to_eval='L(6717)/L(6731)')
        neSII_dict[galaxy][k_comp] = neSII

        # Compute density distribution
        neSII_dist = S2.getTemDen(RSII_dist, tem=temp_dist, to_eval='L(6717)/L(6731)')
        neSIIdist_dict[galaxy][k_comp] = neSII_dist

        # Compute the percentiles
        comp_den_percentiles = np.percentile(neSII_dist, percentiles_array*100)
        # comp_den_percentiles = np.empty(len(percentiles_array))
        # for i_per, percen in enumerate(percentiles_array):
        #     comp_den_percentiles[i_per] = np.percentile(neSII_dist, percentiles_array)

        # Print results
        results_line = f'{galaxy}-{k_comp}: mean density {neSII:.2f}, ' \
                       f'distribution {np.nanmean(neSII_dist):.2f} +/- {np.nanstd(neSII_dist):.2f}, ' \
                       f'emissivity errors {np.isnan(neSII_dist).sum()}'
        percentiles_line = f'percentiles: {[np.round(per, 2) for per in comp_den_percentiles]}'
        print(results_line)
        print(percentiles_line)
        print()


# Plot the results
colors = {'J1152': ['tab:blue', 'tab:red', 'tab:purple'],
          'J0925': ['tab:blue', 'tab:green', 'tab:red', 'tab:purple']}

compLabels = {'J1152': ['NARROW', 'BROAD', 'GLOBAL'],
              'J0925': ['NARROW1', 'NARROW2', 'BROAD', 'GLOBAL']}

for i, galaxy in enumerate(objList):

    line_components = compLabels[galaxy]
    fig, axs = plt.subplots(1, len(line_components), figsize=(len(line_components)*3, 3), sharey=True)

    for j, comp in enumerate(line_components):
        axs[j].hist(neSIIdist_dict[galaxy][comp], histtype='stepfilled', bins=30, alpha=.7, density=False, color=colors[galaxy][j])
        axs[j].axvline(neSII_dict[galaxy][comp], color=colors[galaxy][j], linestyle='--', label=f'ne[SII] mean: {neSII_dict[galaxy][comp]:.0f} ' + r'$cm^{-3}$')

        # Plot distribution uncertainty region
        dist_mean, dist_std = neSIIdist_dict[galaxy][comp].mean(), neSIIdist_dict[galaxy][comp].std()
        data_label = f'ne[SII] distribution: ' + r'${:.0f}\pm{:.0f}\,cm^{{-3}}$'.format(dist_mean, dist_std)
        axs[j].axvspan(dist_mean - dist_std, dist_mean + dist_std, alpha=0.2, color='grey', label=data_label)

        # Plot labels
        axs[j].set_title(comp)
        axs[j].legend()

    fig.suptitle(f'Galaxy {galaxy}')
    plt.tight_layout()
    plt.show()
