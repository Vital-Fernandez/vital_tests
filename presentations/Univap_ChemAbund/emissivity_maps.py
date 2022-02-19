import numpy as np
import pyneb as pn
from matplotlib import pyplot as plt, rcParams, cm
from astro.data.muse.common_methods import DARK_PLOT
from scipy.optimize import curve_fit
from delete.data_printing import label_decomposition

defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

def fitEmis(func_emis, xy_space, line_emis, p0=None):
    p1, p1_cov = curve_fit(func_emis, xy_space, line_emis, p0)
    return p1, p1_cov

def emisEquation_HeI(xy_space, a, b, c, d):
    temp_range, den_range = xy_space
    return (a + b * den_range) * np.power(temp_range/10000.0, c + d * den_range)

def emisTeNe_2DPlot(coeffs, func_emis, te_ne_grid, emis_grid, line_label):

    # Generate fitted surface points
    matrix_edge = int(np.sqrt(te_ne_grid[0].shape[0]))
    surface_points = func_emis(te_ne_grid, *coeffs)

    # Compare pyneb values with values from fitting
    percentage_difference = (1 - surface_points/emis_grid) * 100

    ion, wave, latex_label = label_decomposition(line_label, scalar_output=True)

    # Generate figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot plane
    ax.imshow(surface_points.reshape((matrix_edge, matrix_edge)),
               aspect=0.03,
               extent=(te_ne_grid[1].min(), te_ne_grid[1].max(), te_ne_grid[0].min(), te_ne_grid[0].max()))

    # Points with error below 1.0 are transparent:
    idx_interest = percentage_difference < 0.5
    ax.scatter(te_ne_grid[1][idx_interest], te_ne_grid[0][idx_interest], c="None", edgecolors='black', linewidths=0.35, label='Error below 1%')

    if idx_interest.sum() < emis_grid.size:

        # Plot grid points
        sc_im = ax.scatter(te_ne_grid[1][~idx_interest], te_ne_grid[0][~idx_interest], c=percentage_difference[~idx_interest],
                    edgecolors='black', linewidths=0.1, cmap=cm.OrRd, label='Error above 1%')

        # Color bar
        cbar = fig.colorbar(sc_im)
        cbar.ax.set_ylabel('Discrepancy (%)', rotation=270, labelpad=20)

    # Add labels
    ax.set_ylabel('Temperature $(K)$', fontsize=15)
    ax.set_xlabel('Density ($cm^{-3}$)', fontsize=15)
    title_label = f'{latex_label}, emissivity grid versus parametrisation'
    ax.set_title(title_label, fontsize=15)

    ax.set_ylim(te_ne_grid[0].min(), te_ne_grid[0].max())
    ax.set_xlim(te_ne_grid[1].min(), te_ne_grid[1].max())

    # Display the plot
    ax.legend(loc='lower right', framealpha=1)
    plt.show()

    return




# Range of temperatures and densities
Te_range = np.linspace(7500, 25000, 20)
ne_array = np.linspace(1, 500, 20)

# Param grids
X, Y = np.meshgrid(Te_range, ne_array)
XX = X.flatten()
YY = Y.flatten()

# Pyneb ion
H1 = pn.RecAtom('H', 1)
He1 = pn.RecAtom('He', 1)

HBeta = H1.getEmissivity(XX, YY, wave=4861, product=False)

lineLabel, wave = 'He1_7065A', 7065.0
lineEmis = He1.getEmissivity(XX, YY, wave=wave, product=False) / HBeta

synth_coefs = {'He1_4471A': np.array([2.0301, 1.5e-5, 0.1463, -0.0005]),
               'He1_5876A': np.array([0.745, -5.1e-5, 0.226, -0.0011]),
               'He1_6678A': np.array([2.612, -0.000146, 0.2355, -0.0016]),
               'He1_7065A': np.array([4.329, -0.0024, -0.368, -0.0017])}

# New emissivities
emiss_coeffs, cov1 = fitEmis(emisEquation_HeI, (XX, YY), lineEmis, p0=synth_coefs[lineLabel])
print(emiss_coeffs)
# emiss_coeffs = synth_coefs[lineLabel]
# print(emiss_coeffs)

emisTeNe_2DPlot(emiss_coeffs, emisEquation_HeI, (XX, YY), lineEmis, lineLabel)
