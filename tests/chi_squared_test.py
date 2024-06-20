import numpy as np
import matplotlib.pyplot as plt
from lime import theme
from matplotlib import pyplot as plt, rc_context
theme.set_style('dark', fig_cfg={"axes.labelsize": 14, "legend.fontsize": 12})

# Example data
true_values = np.array([6, 1000, 3000])
observed_values = np.array([5, 1004, 3005])
uncertainties = np.array([2, 2, 2])  # Original uncertainties

# Define a range for the observed value of the first data point (-5% to +5% of its true value)
percent_range = np.linspace(0, 0.05, 20)

line_list = [r'$[OIII]4363\AA$', r'$[OIII]4959\AA$', r'$[OIII]5007\AA$']

fig_cfg = theme.fig_defaults()
with rc_context(fig_cfg):

    # Plotting the results
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, theo_value in enumerate(true_values):
        chi_squared_values = []
        theo_value = true_values[i]
        for percent in percent_range:
            varied_true = theo_value * (1 + percent)  # Varying the first observed value
            current_true_values = true_values.copy()
            current_true_values[i] = varied_true  # Updating the first observed value
            chi_squared = np.sum(((observed_values - current_true_values) ** 2) / (uncertainties ** 2))
            chi_squared_values.append(chi_squared)

        chi_squared_values = np.array(chi_squared_values)/3
        ax.plot(percent_range * 100, chi_squared_values, linestyle='-', label=f'Line {line_list[i]}')

    ax.update({'xlabel': 'Line flux variation in percentage',
               'ylabel': r'$\chi^{2}$ variation'})

    ax.legend()
    plt.show()