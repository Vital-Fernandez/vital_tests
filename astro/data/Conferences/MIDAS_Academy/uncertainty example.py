from matplotlib import pyplot as plt, rc_context
import numpy as np
import lime

lime.theme.set_style('dark')
fig_conf = lime.theme.fig_defaults()

from matplotlib import patheffects

with rc_context(fig_conf):
    # Plot a straight diagonal line with ticked style path
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot a curved line with ticked style path
    nx = 101
    x = np.linspace(0.0, 1.0, nx)
    y = 0.3*np.sin(x*8) + 0.4 + (1 * x + 2)
    ax.plot(x, y, label="Curve")

    ax.legend()

    plt.show()